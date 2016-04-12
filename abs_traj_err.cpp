#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <Eigen/Geometry>

std::vector<Eigen::Vector3d>  hector_pose;
std::vector<double> hector_time;
std::vector<Eigen::Vector3d>  camera_pose;
std::vector<double> camera_time;

Eigen::Vector3d findClosest(double t) {
	int i1 = 0;
	int i2 = hector_time.size()-1;
	while (true) {
		if (i2 == i1 + 1)
			break;
		int mid = (i1+i2)/2;
		if (t == hector_time[mid])
			return hector_pose[mid];
		else if (t > hector_time[mid])
			i1 = mid;
		else
			i2 = mid;
	}
	if (t - hector_time[i1] > hector_time[i2] - t)
		return hector_pose[i2];
	else
		return hector_pose[i1];
}

int main(int argc, char* argv[]) {
	if (argc < 3) {
		printf("./abs_traj_err hector_pose.txt camera_pose.txt [aligned_pose.txt]\n");
		return 1;
	}

	char buffer[256];
	FILE* f1 = fopen(argv[1],"r");
	if (!f1)
		return 1;
	while (fgets(buffer,256,f1)) {
		float t,x,y,z,qw,qx,qy,qz;
		if (sscanf(buffer,"%f %f %f %f %f %f %f %f",&t,&x,&y,&z,&qw,&qx,&qy,&qz)==8) {
			Eigen::Vector3d v(x,y,z);
			hector_pose.push_back(v);
			hector_time.push_back(t);
		}
	}
	FILE* f2 = fopen(argv[2],"r");
	if (!f2)
		return 1;
	while (fgets(buffer,256,f2)) {
		float t,x,y,z;
		if (sscanf(buffer,"%f %f %f %f",&t,&x,&y,&z)==4) {
			Eigen::Vector3d v(x,y,z);
			camera_pose.push_back(v);
			camera_time.push_back(t);
		}
	}
	printf("Read %lu hector %lu camera pose\n",hector_pose.size(),camera_pose.size());

	Eigen::Matrix<double,3,Eigen::Dynamic> src;
	Eigen::Matrix<double,3,Eigen::Dynamic> dst;
	src.resize(3,camera_pose.size());
	dst.resize(3,camera_pose.size());

	for (size_t i=0;i<camera_pose.size();i++) {
		src.col(i) = camera_pose[i];
		dst.col(i) = findClosest(camera_time[i]);
	}

	Eigen::Matrix<double,4,4> T = Eigen::umeyama(src,dst,true);
	Eigen::Matrix<double,3,3> R = T.block<3,3>(0,0);
	Eigen::Vector3d t = T.block<3,1>(0,3);
	double c = R.col(0).norm();
	R /= c;
	std::cout << "T=" << T << "\n";
	std::cout << "c=" << c << "\n";
	std::cout << "R=" << R << "\n";
	std::cout << "t=" << t << "\n";

	if (argc >= 4) {
		double RMSE = 0;
		FILE* aligned_pose = fopen(argv[3],"w");
		for (int i=0;i<src.cols();i++) {
			Eigen::Vector3d y = dst.col(i);
			Eigen::Vector3d y_est = c * R * src.col(i) + t;
			Eigen::Vector3d dy = y - y_est;
			RMSE += dy.dot(dy);
			fprintf(aligned_pose,"%f %f %f %f\n",camera_time[i],y_est(0),y_est(1),y_est(2));
		}
		RMSE  = sqrt(RMSE / src.size());
		printf("RMSE: %f\n",RMSE);
		fclose(aligned_pose);
	}
}
