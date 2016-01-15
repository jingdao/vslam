#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <unistd.h>
#include <time.h>
#include <vector>

struct Point {
	float x,y,z;
};

float fx = 971.760406;
float fy = 971.138862;
float cx = 319.500000;
float cy = 239.500000;
float angle_threshold = 1.0 / 360 * M_PI;
std::vector<float> rotations;
std::vector<float> translations;
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

Point getPlane(Point p1,Point p2,Point v1) {
	Point v2 = {
		p2.x - p1.x,
		p2.y - p1.y,
		p2.z - p1.z
	};
	Point normal = {
		v1.y * v2.z - v1.z * v2.y,
		v1.z * v2.x - v1.x * v2.z,
		v1.x * v2.y - v1.y * v2.x
	};
	return normal;
}

float magnitude(Point p) {
	return sqrt(p.x*p.x+p.y*p.y+p.z*p.z);
}

float getAngle(Point plane,Point v3) {
	float cosx = plane.x*v3.x+plane.y*v3.y+plane.z*v3.z;
	cosx /= magnitude(plane) * magnitude(v3);
	return M_PI / 2 - acos(fabs(cosx));
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./check_epipole pose_stamped.txt key.match valid.match\n");
		return 1;
	}
	srand(time(NULL));
	FILE* f = fopen(argv[1],"r");
	float R[9];
	float Q[4];
	float T[3];
	float t;
	while (fgets(buffer,2048,f)) {
		if (sscanf(buffer,"%f %f %f %f %f %f %f %f\n",&t,T,T+1,T+2,Q,Q+1,Q+2,Q+3) == 8) {
			quaternionToRotation(Q[1],Q[2],Q[3],Q[0],R);
			for (int i=0;i<9;i++)
				rotations.push_back(R[i]);
			for (int i=0;i<3;i++)
				translations.push_back(T[i]);
		}
	}
	fclose(f);
	
	FILE* keymatch_file = fopen(argv[2],"r");	
	FILE* validmatch_file = fopen(argv[3],"w");
	while (fgets(buffer,2048,keymatch_file)) {
		std::vector<Point> projections;
		std::vector<Point> displacements;
		std::vector<int> id_list;
		std::vector<float> u_list;
		std::vector<float> v_list;
		char* tok = strtok(buffer," ");
		while (tok) {
			int id = atoi(tok);
			tok = strtok(NULL," \n");
			float u = atof(tok);
			tok = strtok(NULL," \n");
			float v = atof(tok);
			tok = strtok(NULL," \n");
			float* Rwc = &rotations[id*9];
			float* Twc = &translations[id*3];
			float zero[] = {0,0,0};
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
			projections.push_back(transformPoint(map_point,Rwc,zero));
			displacements.push_back(cam_center);
			id_list.push_back(id);
			u_list.push_back(u);
			v_list.push_back(v);
		}
		//RANSAC to get best plane
		Point best_plane;
		float best_error;
		for (unsigned int i=0;i<(projections.size()-1)*2;i++) {
			unsigned int id1 = rand() % projections.size(); 
			unsigned int id2;
			do {
				id2 = rand() % projections.size(); 
			} while (id2 == id1);
			Point plane = getPlane(displacements[id1],displacements[id2],projections[id1]);
			float error = 0;
			for (unsigned int j=0;j<projections.size();j++)
				if (j!=id1)
					error += getAngle(plane,projections[j]);
			if (i==0 || error < best_error) {
				best_plane = plane;
				best_error = error;
			}
		}
		//write valid matches to file
		int numValid = 0;
		char* buffer_c = buffer;
		for (unsigned int i=0;i<projections.size();i++) {
			float theta = getAngle(best_plane,projections[i]);
			if (theta < angle_threshold) {
				if (buffer_c != buffer)
					buffer_c += sprintf(buffer_c," ");
				buffer_c += sprintf(buffer_c,"%d %f %f",id_list[i],u_list[i],v_list[i]);
				numValid++;	
			}
		}
		if (numValid >= 2)
			fprintf(validmatch_file,"%s\n",buffer);
	}
	fclose(keymatch_file);
	fclose(validmatch_file);
}
