#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <Eigen/Core>
//https://en.wikipedia.org/wiki/Essential_matrix#3D_points_from_corresponding_image_points

#define DEBUG_SINGLE 0

typedef Eigen::Matrix<double,3,3,Eigen::RowMajor> Mat3;
double fx = 567;
double fy = 567;
double cx = 319.500000;
double cy = 239.500000;
//lidar to camera transformation
Mat3 Rcl;
Eigen::Vector3d Tcl;

void quaternionToRotation(double qx,double qy,double qz,double qw,double* R) {
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

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./match_g2o pose_stamped.txt key.match map_point.txt\n");
		return 1;
	}
	srand(time(NULL));
	Rcl << 1,0,0,0,0,1,0,-1,0;
	Tcl << 0,0.06,0;

	//read pose file
	FILE* pose_stamped = fopen(argv[1],"r");
	if (!pose_stamped)
		return 1;
	char buffer[2048];
	std::vector<Mat3> rotations;
	std::vector<Eigen::Vector3d> translations;
	while (fgets(buffer,2048,pose_stamped)) {
		double t,x,y,z,qx,qy,qz,qw;
		if (sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",&t,&x,&y,&z,&qw,&qx,&qy,&qz)==8) {
			double r[9];
			quaternionToRotation(qx,qy,qz,qw,r);
			Mat3 Rwl;
			memcpy(Rwl.data(),r,9*sizeof(double));
			Eigen::Vector3d Twl(x,y,z);
			rotations.push_back(Rcl * Rwl.transpose());
			translations.push_back(- Rcl * Rwl.transpose() * Twl + Tcl);
		} else {
			printf("Error parsing: %s\n",buffer);
		}
	}
	fclose(pose_stamped);

	struct timespec start,end;
	clock_gettime(CLOCK_MONOTONIC,&start);
	int count_points = 0;
	double RMSE = 0;
	FILE* key_match = fopen(argv[2],"r");
	FILE* map_point = fopen(argv[3],"w");
	if (!(key_match && map_point))
		return 1;
	while (fgets(buffer,2048,key_match)) {
#if DEBUG_SINGLE
		printf("key.match: %s",buffer);
#endif
		int id;
		char* tok = strtok(buffer," ");
		std::vector<double> uc,vc;
		std::vector<int> index;
		while (tok) {
			id = atoi(tok);
			index.push_back(id);
			tok = strtok(NULL," \n");
			double u = atof(tok);
			tok = strtok(NULL," \n");
			double v = atof(tok);
			tok = strtok(NULL," \n");
			uc.push_back(u - cx);
			vc.push_back(cy - v);
		}

		//optimize
		Eigen::Vector3d bestEstimate, x1, x2;
		double leastError = -1;
		for (unsigned int i=0;i<index.size()-1;i++) {
			for (unsigned int j=i+1;j<index.size();j++) {
				double u1 = uc[i], v1 = vc[i];
				double u2 = uc[j], v2 = vc[j];
				Mat3 R1 = rotations[index[i]], R2 = rotations[index[j]];
				Eigen::Vector3d T1 = translations[index[i]], T2 = translations[index[j]];
				Mat3 Rcc = R2 * R1.transpose();
				Eigen::Vector3d Tcc = T1 - R1 * R2.transpose() * T2;
				Eigen::Vector3d r1 = Rcc.row(0), r2 = Rcc.row(1), r3 = Rcc.row(2);
				Eigen::Vector3d uv(-u1/fx,-v1/fy,1);
				Eigen::Vector3d mult_u = r1 + u2/fx * r3;
				Eigen::Vector3d mult_v = r2 + v2/fy * r3;
				double z_est[2] = {mult_u.dot(Tcc) / mult_u.dot(uv),
								mult_v.dot(Tcc) / mult_v.dot(uv) } ;
				for (int k=0;k<1;k++) {
					x1 << -u1*z_est[k]/fx, -v1*z_est[k]/fy, z_est[k];
					x2 = Rcc * (x1 - Tcc);
					if (x1(2) >= 0 || x2(2) >= 0)
						break;
					double u_est = -fx * x2(0) / x2(2);
					double v_est = -fy * x2(1) / x2(2);
					double error = (u_est-u2) * (u_est-u2) + (v_est-v2) * (v_est-v2);
					if (leastError < 0 || error < leastError) {
						leastError = error;
						bestEstimate = R1.transpose() * (x1 - T1);
					}
				}
			}
		}
		
		//record result
		RMSE += leastError;
#if DEBUG_SINGLE
		printf("reprojection: ");
		char* c = buffer;
		c += sprintf(c,"transformation:\n");
		for (unsigned int i=0;i<index.size();i++) {
			id = index[i];
			Eigen::Vector3d xc = rotations[id] * bestEstimate + translations[id];
			double u = - fx * xc(0) / xc(2);
			double v = - fy * xc(1) / xc(2);
			printf("%d %f %f ",id,u+cx,cy-v);
			c += sprintf(c,"%4.2f %4.2f %4.2f\n%4.2f %4.2f %4.2f\n%4.2f %4.2f %4.2f\n",
						rotations[id](0,0),rotations[id](0,1),rotations[id](0,2),
						rotations[id](1,0),rotations[id](1,1),rotations[id](1,2),
						rotations[id](2,0),rotations[id](2,1),rotations[id](2,2));
			c += sprintf(c,"[%4.2f %4.2f %4.2f]\n",translations[id](0),translations[id](1),translations[id](2));
		}
		printf("\n%s",buffer);
		printf("estimate: %f %f %f %lu %f\n",bestEstimate(0),bestEstimate(1),bestEstimate(2),index.size(),leastError);
#endif
		fprintf(map_point,"%f %f %f %lu %f\n",bestEstimate(0),bestEstimate(1),bestEstimate(2),index.size(),leastError);
		count_points++;
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double dt = end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec);
	RMSE = sqrt(RMSE / count_points);
	printf("Optimized %d map points (%fs, RMSE = %f)\n",count_points, dt, RMSE);
	fclose(key_match);
	fclose(map_point);

	return 0;
}
