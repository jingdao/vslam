#include <stdio.h>
#include <string.h>
#include <vector>
#include <math.h>
#include <float.h>
#include <stdlib.h>
#include <time.h>
#include <SDL/SDL.h>
#include "geometry.h"
#define ZERO_THRESHOLD 0.01
#define SAVE_IMAGES 0
#define NOISE 20

std::vector<Point> vertices;
std::vector<Triangle> faces;
std::vector<Plane> planes;
int imageWidth=640,imageHeight=480;
double sensorWidth=0.004,sensorHeight=0.003,focalLength=0.005;
double distortion_k1=0;
int numSteps=20;
double stepSize=1,turnRate=0.1;
SDL_Surface* screen = NULL;
SDL_Surface* screen_left = NULL;
SDL_Surface* screen_right = NULL;
int count = 0;
//unsigned char faceColor[] = {30,50,70,70,170,170,170,170,70,70,230,230};
//unsigned char faceColor[] = {50,50,100,100,50,50,100,100,150,150,150,150,150,150,150,150};
std::vector<int> faceColor;

void pointDetection(int* src,int width,int height);

double magnitude(Vector v) {
	return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vector normalize(Vector v) {
	double m = magnitude(v);
	Vector w = {v.x/m, v.y/m, v.z/m};
	return w;
}

bool readPLY(char* filename, std::vector<Point> *vertices, std::vector<Triangle> *faces) {
	FILE* f = fopen(filename, "r");
	if (!f) {
		printf("File not found: %s\n", filename);
		return false;
	}
	char buf[256];
	int numVertex,numFace;
	while (fgets(buf, 256, f)) {
		if (sscanf(buf,"element vertex %d",&numVertex)==1) {
		} else if (sscanf(buf,"element face %d",&numFace)==1) {
		} else if (strncmp(buf,"end_header",10)==0) {
			for (int i=0;i<numVertex;i++) {
				if (!fgets(buf,256,f))
					break;
				Point p;
				if (sscanf(buf, "%lf %lf %lf",&(p.x),&(p.y),&(p.z)) == 3) {
					vertices->push_back(p);
				} else {
					printf("Error parsing %s\n",filename);
					printf("Line %d: %s\n",i,buf);
					break;
				}
			}
			for (int i=0;i<numFace;i++) {
				if (!fgets(buf,256,f))
					break;
				Triangle t;
				int r,g,b;
				int scanned = sscanf(buf, "3 %lu %lu %lu %d %d %d",&(t.id1),&(t.id2),&(t.id3),&r,&g,&b);
				if (scanned == 6) {
					faces->push_back(t);
					faceColor.push_back(r<<16|g<<8|b);
				} else if (scanned == 3) {
					faces->push_back(t);
					faceColor.push_back(rand()%(1<<24));
				} else {
					printf("Error parsing %s\n",filename);
					printf("Line %d: %s\n",i,buf);
					break;
				}
			}
			break;
		}
	}
	fclose(f);
	return true;
}

void saveImage(char* imageName,unsigned char* imageBuffer,int width,int height,bool grayscale) {
	FILE* f = fopen(imageName,"w");
	if (!f)
		return;
	if (grayscale) {
		fprintf(f,"P5\n%d %d\n255\n",width,height);
		fwrite(imageBuffer,1,width*height,f);
	} else {
		fprintf(f,"P6\n%d %d\n255\n",width,height);
//		fwrite(imageBuffer,1,width*height*3,f);
		int* c = (int*)imageBuffer;
		unsigned char r,g,b;
		for (int i=0;i<height;i++) {
			for (int j=0;j<width;j++) {
				b = *c & 0xFF;
				g = (*c >> 8) & 0xFF;
				r = (*c >> 16) & 0xFF;
				fwrite(&r,1,1,f);
				fwrite(&g,1,1,f);
				fwrite(&b,1,1,f);
				c++;
			}
		}
	}
	printf("Saved to %s\n",imageName);
	fclose(f);
}

Plane getPlane(std::vector<Point> *vertices, Triangle tri) {
	Point p1 = (*vertices)[tri.id1];
	Point p2 = (*vertices)[tri.id2];
	Point p3 = (*vertices)[tri.id3];
	Vector v1 = {p1.x - p2.x, p1.y - p2.y, p1.z - p2.z};
	Vector v2 = {p1.x - p3.x, p1.y - p3.y, p1.z - p3.z};
	Vector crossProduct = {
		v1.y * v2.z - v1.z * v2.y,
		v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x
	};
	crossProduct = normalize(crossProduct);
	Plane plane = {
		crossProduct.x,
		crossProduct.y,
		crossProduct.z,
		p1.x*crossProduct.x + p1.y*crossProduct.y + p1.z*crossProduct.z
	};
	return plane;
}

bool intersects(Point rayOrigin,Vector rayDirection,Plane plane,double *distance) {
	//if ray perpendicular to triangle
	if (fabs(rayDirection.x*plane.a + rayDirection.y*plane.b + rayDirection.z*plane.c) < ZERO_THRESHOLD)
		return false;
	//let intersection P = rayOrigin + distance * rayDirection
	//					P lies in plane described by ax+by+cz=d
	//calculate distance
	*distance = -(plane.a*rayOrigin.x + plane.b*rayOrigin.y + plane.c*rayOrigin.z - plane.d) / 
				(plane.a*rayDirection.x + plane.b*rayDirection.y + plane.c*rayDirection.z);
	if (*distance <= 0)
		return false;
	return true;
}

bool triangleContains(std::vector<Point> *vertices,Triangle tri,Plane plane,Point target) {
	Point p1 = (*vertices)[tri.id1];
	Point p2 = (*vertices)[tri.id2];
	Point p3 = (*vertices)[tri.id3];
	Vector v1 = {p2.x - p1.x, p2.y - p1.y, p2.z - p1.z};
	Vector v2 = {p3.x - p2.x, p3.y - p2.y, p3.z - p2.z};
	Vector v3 = {p1.x - p3.x, p1.y - p3.y, p1.z - p3.z};
	Vector c1 = {target.x - p1.x, target.y - p1.y, target.z - p1.z};
	Vector c2 = {target.x - p2.x, target.y - p2.y, target.z - p2.z};
	Vector c3 = {target.x - p3.x, target.y - p3.y, target.z - p3.z};
	double tripleproduct1 =
		plane.a * (v1.y * c1.z - v1.z * c1.y) + 
		plane.b * (v1.z * c1.x - v1.x * c1.z) + 
		plane.c * (v1.x * c1.y - v1.y * c1.x);
	double tripleproduct2 =
		plane.a * (v2.y * c2.z - v2.z * c2.y) + 
		plane.b * (v2.z * c2.x - v2.x * c2.z) + 
		plane.c * (v2.x * c2.y - v2.y * c2.x);
	double tripleproduct3 =
		plane.a * (v3.y * c3.z - v3.z * c3.y) + 
		plane.b * (v3.z * c3.x - v3.x * c3.z) + 
		plane.c * (v3.x * c3.y - v3.y * c3.x);

	return tripleproduct1 >= 0 && tripleproduct2 >= 0 && tripleproduct3 >= 0;
}

void drawPoint(int* imageBuffer,int width,int height,Coord point,int pointSize,int color) {
	int xmin = point.x - pointSize/2 < 0 ? 0 : point.x - pointSize/2;
	int xmax = point.x + pointSize/2 >= width ? width-1 : point.x + pointSize/2;
	int ymin = point.y - pointSize/2 < 0 ? 0 : point.y - pointSize/2;
	int ymax = point.y + pointSize/2 >= height ? height-1 : point.y + pointSize/2;
	int* p = imageBuffer + ymin*width + xmin;
	for (int i=ymin;i<=ymax;i++) {
		for (int j=xmin;j<=xmax;j++) {
			*p++ = color;
		}
		p += width - 1 - xmax + xmin;
	}
}

void draw(unsigned char* imageBuffer,Point rayOrigin,double yaw) {
	memset(imageBuffer,0,imageWidth*imageHeight*4);
	printf("%.1f %.1f %.2f\n",rayOrigin.x,rayOrigin.y,yaw);
	for (int i=0;i<imageHeight;i++) {
		for (int j=0;j<imageWidth;j++) {
			Vector principalDirection = {
				1.0 * (j-imageWidth/2) / imageWidth * sensorWidth,
				focalLength,
				1.0 * (i-imageHeight/2) / imageHeight * sensorHeight
			};
			double distortion = 1 + distortion_k1 * (principalDirection.x*principalDirection.x + principalDirection.z*principalDirection.z);
			principalDirection.x *= distortion;
			principalDirection.z *= distortion;
			principalDirection = normalize(principalDirection);
			Vector rayDirection = {
				principalDirection.x * cos(yaw) - principalDirection.y * sin(yaw),
				principalDirection.x * sin(yaw) + principalDirection.y * cos(yaw),
				principalDirection.z
			};
			bool isValid = false;
			int faceIndex;
			double minDistance = DBL_MAX;
			//find intersection for each triangle
			for (size_t k=0;k<faces.size();k++) {
				double distance;
				if (intersects(rayOrigin,rayDirection,planes[k],&distance)) {
					if (distance < minDistance) {
						Point intersection = {
							rayOrigin.x + rayDirection.x * distance,
							rayOrigin.y + rayDirection.y * distance,
							rayOrigin.z + rayDirection.z * distance
						};
						if (triangleContains(&vertices,faces[k],planes[k],intersection)) {
							isValid = true;
							faceIndex = k;
							minDistance = distance;
						}
					}
				}
			}
			if (isValid) {
				int* imagePointer = (int*) (imageBuffer + (i*imageWidth+j)*4);
				*imagePointer = faceColor[faceIndex];
#if NOISE
				*imagePointer += (rand()%NOISE - NOISE/2)<<16;
				*imagePointer += (rand()%NOISE - NOISE/2)<<8;
				*imagePointer += (rand()%NOISE - NOISE/2);
#endif
			}
		}
	}
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("Usage: ./pixel_tracing input.ply outputFolder\n");
		return 1;
	}

	srand(time(NULL));
	if (!readPLY(argv[1],&vertices,&faces))
		return 1;
	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Pixel tracing", NULL);
	screen = SDL_SetVideoMode(imageWidth*2,imageHeight, 32, SDL_SWSURFACE);
	screen_left = SDL_CreateRGBSurface(SDL_SWSURFACE,imageWidth,imageHeight, 32,0xFF0000,0xFF00,0xFF,0);
	screen_right = SDL_CreateRGBSurface(SDL_SWSURFACE,imageWidth,imageHeight, 32,0xFF0000,0xFF00,0xFF,0);

	//get bounding box
	double minX=vertices[0].x,maxX=vertices[0].x;
	double minY=vertices[0].y,maxY=vertices[0].y;
	double minZ=vertices[0].z,maxZ=vertices[0].z;
	for (size_t i=1;i<vertices.size();i++) {
		if (vertices[i].x < minX) minX = vertices[i].x;
		else if (vertices[i].x > maxX) maxX = vertices[i].x;
		if (vertices[i].y < minY) minY = vertices[i].y;
		else if (vertices[i].y > maxY) maxY = vertices[i].y;
		if (vertices[i].z < minZ) minZ = vertices[i].z;
		else if (vertices[i].z > maxZ) maxZ = vertices[i].z;
	}
	printf("Bounding box: x:(%.2f %.2f) y:(%.2f %.2f) z:(%.2f %.2f)\n",minX,maxX,minY,maxY,minZ,maxZ);
	Point centroid = {
		(minX + maxX) / 2,
		(minY + maxY) / 2,
		(minZ + maxZ) / 2
	};

	//get normals
	for (size_t i=0;i<faces.size();i++) {
		Plane v = getPlane(&vertices,faces[i]);
		planes.push_back(v);
	}

	Point cameraOrigin = {0,0,1};
	Vector principalDirection = {0,1,0};
	SDL_Rect dstrect = {imageWidth,0,0,0};
	double yaw = 0;
	draw((unsigned char*)screen_left->pixels,cameraOrigin,yaw);
	char buf[128];
#if SAVE_IMAGES
	sprintf(buf,"%s/%d.ppm",argv[2],count++);
	saveImage(buf,(unsigned char*)screen_left->pixels,imageWidth,imageHeight,false);
#endif
	SDL_BlitSurface(screen_left,NULL,screen,NULL);
	memcpy(screen_right->pixels,screen_left->pixels,imageWidth*imageHeight*4);
	pointDetection((int*)screen_right->pixels,imageWidth,imageHeight);
	SDL_BlitSurface(screen_right,NULL,screen,&dstrect);
	SDL_Flip(screen);
	SDL_Event event;
	while (true) {
		while (SDL_PollEvent(&event)) {
			switch(event.type){
				case SDL_KEYDOWN:
					switch( event.key.keysym.sym ){
						case SDLK_LEFT:
						yaw += turnRate;
						principalDirection.x = -sin(yaw);
						principalDirection.y = cos(yaw);
						break;
						case SDLK_RIGHT:
						yaw -= turnRate;
						principalDirection.x = -sin(yaw);
						principalDirection.y = cos(yaw);
						break;
						case SDLK_UP:
						cameraOrigin.x += principalDirection.x * stepSize;
						cameraOrigin.y += principalDirection.y * stepSize;
						break;
						case SDLK_DOWN:
						cameraOrigin.x -= principalDirection.x * stepSize;
						cameraOrigin.y -= principalDirection.y * stepSize;
						break;
						case 'q':
						exit(0);
						default:
						break;
					}
					draw((unsigned char*)screen_left->pixels,cameraOrigin,yaw);
#if SAVE_IMAGES
					sprintf(buf,"%s/%d.ppm",argv[2],count++);
					saveImage(buf,(unsigned char*)screen_left->pixels,imageWidth,imageHeight,true);
#endif
					SDL_BlitSurface(screen_left,NULL,screen,NULL);
					memcpy(screen_right->pixels,screen_left->pixels,imageWidth*imageHeight*4);
					pointDetection((int*)screen_right->pixels,imageWidth,imageHeight);
					SDL_BlitSurface(screen_right,NULL,screen,&dstrect);
					SDL_Flip(screen);
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
	}
	return 0;
}
