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
bool* camera_valid;
FILE* keymatch_file,*mappoint_file;
std::vector<CamModel> camera_location;
std::vector<float> rotations;
std::vector<float> translations;
std::vector<Point> lines;
Point target;
int mouseIndex = 0;
int previousX,previousY;
double scrollSpeed = 1.1;
double lineLength = 5;
float fx = 971.760406;
float fy = 971.138862;
float cx = 319.500000;
float cy = 239.500000;
char buffer[2048];

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

Point scalePoint(Point p,float length) {
	Point q;
	double magnitude = sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
	q.x = p.x / magnitude * length;
	q.y = p.y / magnitude * length;
	q.z = p.z / magnitude * length;
	return q;
}

void getKeymatch() {
	memset(camera_valid,0,camera_location.size()*sizeof(bool));
	lines.clear();
	if (fgets(buffer,2048,keymatch_file)) {
		char* tok = strtok(buffer," ");
		while (tok) {
			int id = atoi(tok);
			camera_valid[id] = true;
			tok = strtok(NULL," \n");
			float u = atof(tok);
			tok = strtok(NULL," \n");
			float v = atof(tok);
			tok = strtok(NULL," \n");
			float* Rwc = &rotations[id*9];
			float* Twc = &translations[id*3];
			Point cam_center = {
				Twc[0],
				Twc[1],
				Twc[2]
			};
			Point map_point = {
				(u-cx) * fy,
				fx * fy,
				(cy-v) * fx
			};
			lines.push_back(cam_center);
			lines.push_back(transformPoint(scalePoint(map_point,lineLength),Rwc,Twc));
		}
	}
}

void getMappoint() {
	if (fgets(buffer,2048,mappoint_file)) {
		sscanf(buffer,"%f %f %f",&target.x,&target.y,&target.z);
	}
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

	glLineWidth(2.0);
	glColor3ub(0,255,0);
	for (size_t n = 0; n < camera_location.size(); n++){
		if (camera_valid[n]) {
			drawLine(camera_location[n].center,camera_location[n].ul);
			drawLine(camera_location[n].center,camera_location[n].ur);
			drawLine(camera_location[n].center,camera_location[n].bl);
			drawLine(camera_location[n].center,camera_location[n].br);
			drawLine(camera_location[n].ul,camera_location[n].bl);
			drawLine(camera_location[n].ur,camera_location[n].br);
			drawLine(camera_location[n].ul,camera_location[n].ur);
			drawLine(camera_location[n].bl,camera_location[n].br);
		}
	}

	glLineWidth(1.0);
	glColor3ub(255,255,255);
	for (size_t n = 0; n < lines.size() / 2; n++) {
		drawLine(lines[2*n],lines[2*n+1]);
	}

	glPointSize(5.0);
	glBegin(GL_POINTS);
	glColor3ub(255,0,0);
	glVertex3d(target.x,target.y,target.z);
	glEnd();

	glFlush();
	SDL_GL_SwapBuffers();

	glPopMatrix();
	glPopAttrib();
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./plot_epipole pose_stamped.txt key.match map_point.txt [point_index]\n");
		return 1;
	}
	int point_index = 0;
	if (argc > 4)
		point_index = atoi(argv[4]);

	FILE* f = fopen(argv[1],"r");
	float R[9];
	float Q[4];
	float T[3];
	float t;
	Point center = {0,0,0};
	Point ul = {-cameraSize,cameraSize,cameraSize};
	Point ur = {cameraSize,cameraSize,cameraSize};
	Point bl = {-cameraSize,cameraSize,-cameraSize};
	Point br = {cameraSize,cameraSize,-cameraSize};
	while (fgets(buffer,2048,f)) {
		if (sscanf(buffer,"%f %f %f %f %f %f %f %f\n",&t,T,T+1,T+2,Q,Q+1,Q+2,Q+3) == 8) {
			CamModel cam;
			quaternionToRotation(Q[1],Q[2],Q[3],Q[0],R);
			cam.center = transformPoint(center,R,T);
			cam.ul = transformPoint(ul,R,T);
			cam.ur = transformPoint(ur,R,T);
			cam.bl = transformPoint(bl,R,T);
			cam.br = transformPoint(br,R,T);
			camera_location.push_back(cam);
			for (int i=0;i<9;i++)
				rotations.push_back(R[i]);
			for (int i=0;i<3;i++)
				translations.push_back(T[i]);
		}
	}
	fclose(f);
	camera_valid = (bool*) malloc(camera_location.size()*sizeof(bool));
	keymatch_file = fopen(argv[2],"r");
	for (int i=0;i<point_index;i++)
		fgets(buffer,2048,keymatch_file);
	getKeymatch();
	mappoint_file = fopen(argv[3],"r");
	for (int i=0;i<point_index;i++)
		fgets(buffer,2048,mappoint_file);
	getMappoint();

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
						case 'n':
						getKeymatch();
						getMappoint();
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
	free(camera_valid);
	fclose(keymatch_file);
	fclose(mappoint_file);

	return 0;
}
