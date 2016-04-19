#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

struct PCD {
	int numPoints;
	float* float_data;
};
enum PCD_data_storage {
	ASCII,
	BINARY,
	NONE
};
struct Point {
	float x,y,z;
};
struct Plane {
	float a,b,c,d;
};
struct Line {
	float a,b,c;
	float left,right;
	float x0,y0;
};
struct Color {
	unsigned char r,g,b;
};
struct Quadrilateral {
	float x0,y0,z0;
	float x1,y1,z1;
	float x2,y2,z2;
	float x3,y3,z3;
};

std::vector<Plane> planes;
std::vector<Line> lines;
std::vector<Color> colors;
std::vector<int> matchIndex;
std::vector<Quadrilateral> quads;
std::vector<Point> pointcloud;
double cameraX=0,cameraY=-5,cameraZ=3;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.1;
int numSegments = 2;

PCD* NewPCD(const char* fileName) {
	PCD* pcd = new PCD();
	PCD_data_storage data_storage = NONE;
	FILE* f = fopen(fileName, "r");
	if (!f) {
		printf("File not found: %s\n", fileName);
		return NULL;
	}
	char buf[256];
	int pointsParsed = 0;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf, "POINTS %d", &pcd->numPoints) == 1) {
			pcd->float_data = (float*)malloc(3 * pcd->numPoints * sizeof(float));
		} else if (strncmp(buf,"DATA ascii",10)==0) {
			data_storage = ASCII;
		} else if (strncmp(buf,"DATA binary_compressed",23)==0) {
			data_storage = BINARY;
			fread(pcd->float_data,sizeof(float),pcd->numPoints*3,f);
			break;
		}
		else if (data_storage == ASCII) {
			if (sscanf(buf, "%f %f %f", pcd->float_data+pointsParsed * 3, pcd->float_data+pointsParsed * 3 + 1,
				pcd->float_data+pointsParsed * 3 + 2) >= 3) {
				pointsParsed++;
			}
		}
	}
	fclose(f);
	return pcd;
}

Plane segmentPlane(PCD* cloud,int iter,float threshold) {
	int numPoints = cloud->numPoints;
	float* float_data = cloud->float_data;
	Plane bestPlane = {0,0,0,0};
	int maxInliers = 0;
	int i;
	for (i=0;i<iter;i++) {
		//Pick 3 points
		int p0 = rand() % numPoints;
		int p1 = rand() % numPoints;
		int p2 = rand() % numPoints;
		float x[3] = {float_data[p0*3], float_data[p1*3], float_data[p2*3]};
		float y[3] = {float_data[p0*3+1], float_data[p1*3+1], float_data[p2*3+1]};
		float z[3] = {float_data[p0*3+2], float_data[p1*3+2], float_data[p2*3+2]};
		Plane currentPlane = {
			(y[1]-y[0])*(z[2]-z[0]) - (y[2]-y[0])*(z[1]-z[0]),
			(z[1]-z[0])*(x[2]-x[0]) - (z[2]-z[0])*(x[1]-x[0]),
			(x[1]-x[0])*(y[2]-y[0]) - (x[2]-x[0])*(y[1]-y[0]),
			0
		};
		currentPlane.d = -(currentPlane.a * x[0] + currentPlane.b * y[0] + currentPlane.c * z[0]);
		if (currentPlane.a == 0 && currentPlane.b == 0 && currentPlane.c ==0 )
			continue; //picked collinear points
		float distanceThreshold = threshold * sqrt(
			currentPlane.a * currentPlane.a + 
			currentPlane.b * currentPlane.b + 
			currentPlane.c * currentPlane.c
			);
		int numInliers = 0;
		for (int j=0;j<numPoints;j++) {
			if ( fabs( currentPlane.a * float_data[j * 3] +
				 currentPlane.b * float_data[j * 3 + 1] +
				 currentPlane.c * float_data[j * 3 + 2] +
				 currentPlane.d )
				 < distanceThreshold)	
				numInliers++;
		}
		if (numInliers > maxInliers) {
			maxInliers = numInliers;
			bestPlane = currentPlane;
		}
	}
	printf("RANSAC: found plane with %d inliers\n",maxInliers);
	//filter inliers
	float* buffer = new float[numPoints * 3];
	memcpy(buffer,float_data,numPoints*3*sizeof(float));
	free(cloud->float_data);
	cloud->numPoints = numPoints - maxInliers;
	cloud->float_data = (float*) malloc(cloud->numPoints * 3 * sizeof(float));
	float norm = sqrt(bestPlane.a*bestPlane.a + bestPlane.b*bestPlane.b + bestPlane.c*bestPlane.c);
	bestPlane.a /= norm;
	bestPlane.b /= norm;
	bestPlane.c /= norm;
	bestPlane.d /= norm;
	bool valid = fabs(bestPlane.c) < 0.5;
	int k=0;
	for (int j=0;j<numPoints;j++) {
		if ( fabs( bestPlane.a * buffer[j * 3] +
			 bestPlane.b * buffer[j * 3 + 1] +
			 bestPlane.c * buffer[j * 3 + 2] +
			 bestPlane.d )
			 >= threshold) {
			cloud->float_data[k++] = buffer[j*3];
			cloud->float_data[k++] = buffer[j*3+1];
			cloud->float_data[k++] = buffer[j*3+2];
		} else if (valid) {
			Point p = {buffer[j*3],buffer[j*3+1],buffer[j*3+2]};
			pointcloud.push_back(p);
		}
	}
	delete[] buffer;
	if (valid)
		return bestPlane;
	else
		return segmentPlane(cloud,iter,threshold);
}

Line segmentLine(PCD* cloud,int iter,float threshold) {
	int numPoints = cloud->numPoints;
	float* float_data = cloud->float_data;
	Line bestLine;
	int maxInliers = 0;
	int i;
	for (i=0;i<iter;i++) {
		//Pick 3 points
		int p0 = rand() % numPoints;
		int p1 = rand() % numPoints;
		float x[2] = {float_data[p0*3], float_data[p1*3]};
		float y[2] = {float_data[p0*3+1], float_data[p1*3+1]};
		Line currentLine = {
			y[1] - y[0],
			x[0] - x[1],
			0,0,0,0,0
		};
		currentLine.c = -(currentLine.a * x[0] + currentLine.b * y[0]);
		if (currentLine.a == 0 && currentLine.b == 0 )
			continue;
		float norm = sqrt(currentLine.a*currentLine.a + currentLine.b*currentLine.b);
		currentLine.a /= norm;
		currentLine.b /= norm;
		currentLine.c /= norm;
		int numInliers = 0;
		for (int j=0;j<numPoints;j++) {
			if ( fabs( currentLine.a * float_data[j * 3] +
				 currentLine.b * float_data[j * 3 + 1] +
				 currentLine.c)
				 < threshold) {	
				numInliers++;
			}
		}
		if (numInliers > maxInliers) {
			maxInliers = numInliers;
			bestLine = currentLine;
		}
	}
	printf("RANSAC: found line with %d inliers\n",maxInliers);
	//filter inliers
	float* buffer = new float[numPoints * 3];
	memcpy(buffer,float_data,numPoints*3*sizeof(float));
	free(cloud->float_data);
	cloud->numPoints = numPoints - maxInliers;
	cloud->float_data = (float*) malloc(cloud->numPoints * 3 * sizeof(float));
	int k=0;
	bool setInitial = true;
	for (int j=0;j<numPoints;j++) {
		if ( fabs( bestLine.a * buffer[j * 3] +
			 bestLine.b * buffer[j * 3 + 1] +
			 bestLine.c)
			 >= threshold) {
			cloud->float_data[k++] = buffer[j*3];
			cloud->float_data[k++] = buffer[j*3+1];
			cloud->float_data[k++] = buffer[j*3+2];
		} else {
			if (setInitial) {
				bestLine.x0 = buffer[j*3];
				bestLine.y0 = buffer[j*3+1];
				setInitial = false;
			} else {
				float t = bestLine.b * (buffer[j*3]-bestLine.x0) - bestLine.a * (buffer[j*3+1]-bestLine.y0);
				if (t < bestLine.left)
					bestLine.left = t;
				if (t > bestLine.right)
					bestLine.right = t;
			}
		}
	}
	delete[] buffer;
	return bestLine;
}

float getAngle(Line l, Plane p) {
	float cosx = p.a*l.b - p.b*l.a;
	float magnitude_plane = sqrt(p.a*p.a+p.b*p.b+p.c*p.c);
	float magnitude_line = sqrt(l.a*l.a+l.b*l.b);
	cosx /= magnitude_plane * magnitude_line;
	float theta = acos(fabs(cosx));
	return M_PI / 2 - theta;
}

float norm3d(float x,float y,float z) {
	return sqrt(x*x+y*y+z*z);
}

float heron(float a,float b,float c) {
	float s = (a+b+c) / 2;
	return sqrt(s*(s-a)*(s-b)*(s-c));
}

float getQuadArea(Quadrilateral q) {
	float a = norm3d(q.x1-q.x0, q.y1-q.y0, q.z1-q.z0);
	float b = norm3d(q.x2-q.x1, q.y2-q.y1, q.z2-q.z1);
	float c = norm3d(q.x2-q.x0, q.y2-q.y0, q.z2-q.z0);
	float d = norm3d(q.x3-q.x2, q.y3-q.y2, q.z3-q.z2);
	float e = norm3d(q.x3-q.x0, q.y3-q.y0, q.z3-q.z0);
	return heron(a,b,c) + heron(c,d,e);
}

float getArea(Line l,Plane p) {
	float x0 = l.b * l.left + l.x0;
	float y0 = -l.a * l.left + l.y0;
	float z0 = 0;
	float x1 = l.b * l.right + l.x0;
	float y1 = -l.a * l.right + l.y0;
	float z1 = 0;
	float p_norm = sqrt(p.a*p.a + p.b*p.b + p.c*p.c);
	float d2 = (p.a*x0+p.b*y0+p.d)/p_norm;
	float x2 = x0 - d2*p.a;
	float y2 = y0 - d2*p.b;
	float z2 = -d2*p.c;
	float d3 = (p.a*x1+p.b*y1+p.d)/p_norm;
	float x3 = x1 - d3*p.a;
	float y3 = y1 - d3*p.b;
	float z3 = -d3*p.c;
	Quadrilateral q = {
		x0,y0,z0,
		x2,y2,z2,
		x1,y1,z1,
		x3,y3,z3
	};
	quads.push_back(q);
	float A = getQuadArea(q);
	return A;
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glLineWidth(5.0);
	glBegin(GL_LINES);
	for (size_t i=0;i<lines.size();i++) {
		glColor3ub(colors[i].r,colors[i].g,colors[i].b);
		float x0 = lines[i].b * lines[i].left + lines[i].x0;
		float y0 = -lines[i].a * lines[i].left + lines[i].y0;
		float x1 = lines[i].b * lines[i].right + lines[i].x0;
		float y1 = -lines[i].a * lines[i].right + lines[i].y0;
		glVertex3d(x0,y0,0);
		glVertex3d(x1,y1,0);
	}
	glEnd();
	glPointSize(3.0);
	glBegin(GL_POINTS);
	glColor3ub(255,255,255);
	for (size_t n = 0; n < pointcloud.size(); n++){
		glVertex3d(pointcloud[n].x, pointcloud[n].y, pointcloud[n].z);
	}
	glEnd();
	glBegin(GL_QUADS);
	for (size_t i=0;i<matchIndex.size();i++) {
		Color c = colors[matchIndex.size()+i];
		Quadrilateral q = quads[i];
		glColor3ub(c.r,c.g,c.b);
		glVertex3d(q.x0,q.y0,q.z0);
		glVertex3d(q.x1,q.y1,q.z1);
		glVertex3d(q.x2,q.y2,q.z2);
		glVertex3d(q.x3,q.y3,q.z3);
	}
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();
	glPopMatrix();
	glPopAttrib();
}

int main(int argc, char* argv[]) {
	if (argc < 4) {
		printf("./map_err refCloud.pcd cameraCloud.pcd laserCloud.pcd\n");
		return 1;
	}

	srand(0);
	PCD* refCloud = NewPCD(argv[1]);
	PCD* cameraCloud = NewPCD(argv[2]);
	PCD* laserCloud = NewPCD(argv[3]);
	int numMapPoints = refCloud->numPoints + cameraCloud->numPoints;
	PCD map = {numMapPoints, new float[3*numMapPoints]};
	memcpy(map.float_data,refCloud->float_data,refCloud->numPoints * 3 * sizeof(float));
	memcpy(map.float_data + refCloud->numPoints * 3,cameraCloud->float_data,cameraCloud->numPoints * 3 * sizeof(float));

	for (int i=0;i<numSegments;i++) {
		planes.push_back(segmentPlane(&map,1000,0.1));
		lines.push_back(segmentLine(laserCloud,1000,0.1));
		Color c = {(unsigned char)(rand()%255),(unsigned char)(rand()%255),(unsigned char)(rand()%255)};
		Color d = {(unsigned char)(rand()%255),(unsigned char)(rand()%255),(unsigned char)(rand()%255)};
		colors.push_back(c);
		colors.push_back(d);
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Point Cloud", NULL);
	SDL_SetVideoMode(1100,750, 32, SDL_OPENGL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)640/480,1,1000);

	for (int i=0;i<numSegments;i++) {
		Line l = lines[i];
		float minAngle = getAngle(l,planes[0]);
		int minIndex = 0;
		for (int j=1;j<numSegments;j++) {
			float ar = getAngle(l,planes[j]);
			if (ar < minAngle) {
				minAngle = ar;
				minIndex = j;
			}
		}
		matchIndex.push_back(minIndex);
		float minArea = getArea(l,planes[minIndex]);
		printf("Line %d (length %f): matched plane %d area=%f angle=%f\n",i,l.right-l.left,minIndex,minArea,minAngle);
	}

	SDL_Event event;
	double hz = 60;
	bool updated = false;
	draw();
	while (SDL_PollEvent(&event)); //clear event buffer
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case SDLK_LEFT:
						cameraX -= 1;
						break;
						case SDLK_RIGHT:
						cameraX += 1;
						break;
						case SDLK_UP:
						cameraZ += 1;
						break;
						case SDLK_DOWN:
						cameraZ -= 1;
						break;
						default:
						break;
					}
					updated = true;
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_WHEELUP) {
						cameraX /= scrollSpeed;
						cameraY /= scrollSpeed;
						cameraZ /= scrollSpeed;
						updated = true;
					} else if (event.button.button == SDL_BUTTON_WHEELDOWN) {
						cameraX *= scrollSpeed;
						cameraY *= scrollSpeed;
						cameraZ *= scrollSpeed;
						updated = true;
					} else {
						mouseIndex = event.button.button == SDL_BUTTON_LEFT ? 1 : 2;
						previousX = event.button.x;
						previousY = event.button.y;
					}
					break;
				case SDL_MOUSEBUTTONUP:
					mouseIndex = 0;
					break;
				case SDL_MOUSEMOTION:
					if (mouseIndex == 1) {
						double rho = sqrt(cameraX*cameraX+cameraY*cameraY);
						double xstep = cameraY / rho;
						double ystep = -cameraX / rho;
						cameraX += 0.05 * (event.motion.x-previousX) * xstep;
						cameraY += 0.05 * (event.motion.x-previousX) * ystep;
						cameraZ += 0.05 * (event.motion.y-previousY);
						previousX = event.motion.x;
						previousY = event.motion.y;
						updated = true;
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
			if (updated) {
				draw();
				updated = false;
			}
		}
		usleep(1000000/hz);
	}
	delete[] map.float_data;
}
