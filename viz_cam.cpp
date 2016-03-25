#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <vector>
#include <SDL/SDL.h>
#include <GL/gl.h>
#include <GL/glu.h>

struct Point {
	float x,y,z;
};

struct CamModel {
	Point center,ul,ur,bl,br;
};

double cameraX=0,cameraY=-5,cameraZ=5;
double centerX=0,centerY=0,centerZ=0;
double upX=0,upY=0,upZ=1;
int screenWidth = 1100, screenHeight = 750;
float cameraSize = 0.1;
std::vector<CamModel> camera_location;
std::vector<Point> pointcloud;
std::vector<Point> lidarcloud;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.1;
float Twc[] = {0,-0.25,-0.18};
float cloudMin=0,cloudMax=0;

void colormap (float f,unsigned char *r,unsigned char *g,unsigned char *b) {
	*r=0;
	*g=0;
	*b=0;
	if (f<=0) {
		*b = 128;
	} else if (f <= 0.25) {
		*g = (unsigned char) f / 0.25 * 255;
		*b = (unsigned char) 128 * (1 - f / 0.25);
	} else if (f <= 0.5) {
		*g = 255;
		*r = (unsigned char) (f - 0.25) / 0.25 * 255;
	} else if (f <= 0.75) {
		*r = 255;
		*g = (unsigned char) 255 + (0.5 - f) / 0.25 * 127;
	} else if (f <= 1) {
		*r = 255;
		*g = (unsigned char) 128 * (1 - f) / 0.25;
	} else {
		*r = 255;
	}
}

void quaternionToRotation(float qx,float qy,float qz,float qw,float* R) {
	R[0] = 1 - 2*qy*qy - 2 * qz*qz;
	R[1] = 2*qx*qy - 2 * qz*qw;
	R[2] = 2*qx*qz + 2 * qy*qw;
	R[3] = 2*qx*qy + 2 * qz*qw;
	R[4] = 1 - 2*qx*qx - 2 * qz*qz;
	R[5] = 2*qy*qz - 2 * qx*qw;
	R[6] = 2*qx*qz - 2 * qy*qw;
	R[7] = 2*qy*qz + 2 * qx*qw;
	R[8] = 1 - 2*qx*qx - 2 * qy*qy;
}

Point transformPoint(Point p,float* R,float* T) {
	Point res;
	res.x = R[0] * p.x  + R[1] * p.y + R[2] * p.z + T[0];
	res.y = R[3] * p.x  + R[4] * p.y + R[5] * p.z + T[1];
	res.z = R[6] * p.x  + R[7] * p.y + R[8] * p.z + T[2];
	return res;
}

void drawLine(Point p1,Point p2) {
	glBegin(GL_LINES);
	glVertex3d(p1.x,p1.y,p1.z);
	glVertex3d(p2.x,p2.y,p2.z);
	glEnd();
}

void draw() {
	glPushAttrib(GL_ALL_ATTRIB_BITS);
	glPushMatrix();
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	gluLookAt(cameraX,cameraY,cameraZ,centerX,centerY,centerZ,upX,upY,upZ);

	glLineWidth(1.0);
	for (size_t n = 0; n < camera_location.size(); n++){
		glColor3ub(0,255,0);
		drawLine(camera_location[n].center,camera_location[n].ul);
		drawLine(camera_location[n].center,camera_location[n].ur);
		drawLine(camera_location[n].center,camera_location[n].bl);
		drawLine(camera_location[n].center,camera_location[n].br);
		drawLine(camera_location[n].ul,camera_location[n].bl);
		drawLine(camera_location[n].ur,camera_location[n].br);
		drawLine(camera_location[n].ul,camera_location[n].ur);
		drawLine(camera_location[n].bl,camera_location[n].br);
	}

	glPointSize(5.0);
	glBegin(GL_POINTS);
	unsigned char r,g,b;
	for (size_t n = 0; n < pointcloud.size(); n++){
//		colormap((pointcloud[n].z - cloudMin)/(cloudMax - cloudMin),&r,&g,&b);
		colormap(0.5 + pointcloud[n].z / 4,&r,&g,&b);
		glColor3ub(r,g,b);
		glVertex3d(pointcloud[n].x,pointcloud[n].y,pointcloud[n].z);
	}
	glColor3ub(255,255,255);
	for (size_t n = 0; n < lidarcloud.size(); n++){
		glVertex3d(lidarcloud[n].x,lidarcloud[n].y,lidarcloud[n].z);
	}
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

int main(int argc,char* argv[]) {
	if (argc < 3) {
		printf("./viz_cam pose_stamped.txt map_point.txt [lidar_map.txt]\n");
		return 1;
	}

	FILE* f = fopen(argv[1],"r");
	char buffer[128];
	float R[9];
	float Q[4];
	float T[3];
	float t;
	Point center = {0,0,0};
	Point ul = {-cameraSize,cameraSize,cameraSize};
	Point ur = {cameraSize,cameraSize,cameraSize};
	Point bl = {-cameraSize,cameraSize,-cameraSize};
	Point br = {cameraSize,cameraSize,-cameraSize};
	while (fgets(buffer,128,f)) {
		if (sscanf(buffer,"%f %f %f %f %f %f %f %f\n",&t,T,T+1,T+2,Q,Q+1,Q+2,Q+3) == 8) {
			T[0] += Twc[0];
			T[1] += Twc[1];
			T[2] += Twc[2];
			CamModel cam;
			quaternionToRotation(Q[1],Q[2],Q[3],Q[0],R);
			cam.center = transformPoint(center,R,T);
			cam.ul = transformPoint(ul,R,T);
			cam.ur = transformPoint(ur,R,T);
			cam.bl = transformPoint(bl,R,T);
			cam.br = transformPoint(br,R,T);
			camera_location.push_back(cam);
		}
	}
	fclose(f);
	f = fopen(argv[2],"r");
	while (fgets(buffer,128,f)) {
		Point p;
		if (sscanf(buffer,"%f %f %f",&p.x,&p.y,&p.z)==3) {
			if (pointcloud.size()==0 || p.z < cloudMin)
				cloudMin = p.z;
			if (pointcloud.size()==0 || p.z > cloudMax)
				cloudMax = p.z;
			pointcloud.push_back(p);
		}
	}
	fclose(f);
	if (argc > 3) {
		f = fopen(argv[3],"r");
		if (f) {
			while (fgets(buffer,128,f)) {
				Point p;
				if (sscanf(buffer,"%f %f %f",&p.x,&p.y,&p.z)==3)
					lidarcloud.push_back(p);
			}
			fclose(f);
		}
	}

	SDL_Init(SDL_INIT_VIDEO);
	SDL_WM_SetCaption("Viz Cam", NULL);
	SDL_SetVideoMode(screenWidth,screenHeight, 32, SDL_OPENGL);
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	gluPerspective(70,(double)screenWidth/screenHeight,1,1000);

	int interval = 10000;
	SDL_Event event;
	while (SDL_PollEvent(&event)); //clear event buffer
	draw();
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
					draw();
					break;
				case SDL_MOUSEBUTTONDOWN:
					if (event.button.button == SDL_BUTTON_WHEELUP) {
						cameraX /= scrollSpeed;
						cameraY /= scrollSpeed;
						cameraZ /= scrollSpeed;
						draw();
					} else if (event.button.button == SDL_BUTTON_WHEELDOWN) {
						cameraX *= scrollSpeed;
						cameraY *= scrollSpeed;
						cameraZ *= scrollSpeed;
						draw();
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
						draw();
					}
					break;
				case SDL_QUIT:
					exit(0);
					break;
			}
		}
		usleep(interval);
	}

//	atexit(SQL_Quit);

	return 0;
}
