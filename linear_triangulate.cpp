#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <opencv/cv.h>
float fx = 567;
float fy = 567;
float cx = 319.500000;
float cy = 239.500000;
float invfx = 1.0/fx;
float invfy = 1.0/fy;
int ORB_levels = 8;
float ORB_scale_factor = 1.2;
float ratioFactor = 1.5*ORB_scale_factor;
cv::Mat Rcl;
cv::Mat Tcl;
std::vector<float> vScaleFactor;
std::vector<float> vSigma2;
#define DEBUG_OUTLIER 0

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

float linear_triangulate(cv::KeyPoint kp1, cv::KeyPoint kp2, cv::Mat Rcw1, cv::Mat tcw1, cv::Mat Rcw2, cv::Mat tcw2, cv::Mat *x3D, float* depth) {
    cv::Mat Ow1=-Rcw1.t()*tcw1;
	cv::Mat Rwc1 = Rcw1.t();
	cv::Mat Tcw1(3,4,CV_32F);
	Rcw1.copyTo(Tcw1.colRange(0,3));
	tcw1.copyTo(Tcw1.col(3));
    cv::Mat Ow2=-Rcw2.t()*tcw2;
	cv::Mat Rwc2 = Rcw2.t();
	cv::Mat Tcw2(3,4,CV_32F);
	Rcw2.copyTo(Tcw2.colRange(0,3));
	tcw2.copyTo(Tcw2.col(3));
	cv::Mat xn1 = (cv::Mat_<float>(3,1) << (kp1.pt.x-cx)*invfx, (kp1.pt.y-cy)*invfy, 1.0 );
	cv::Mat ray1 = Rwc1*xn1;
	cv::Mat xn2 = (cv::Mat_<float>(3,1) << (kp2.pt.x-cx)*invfx, (kp2.pt.y-cy)*invfy, 1.0 );
	cv::Mat ray2 = Rwc2*xn2;
	const float cosParallaxRays = ray1.dot(ray2)/(cv::norm(ray1)*cv::norm(ray2));

	if(cosParallaxRays<0 || cosParallaxRays>0.9998) {
#if DEBUG_OUTLIER
		printf("cosParallaxRays: %f\n",cosParallaxRays);
#endif
		return -1;
	}

	// Linear Triangulation Method
	cv::Mat A(4,4,CV_32F);
	A.row(0) = xn1.at<float>(0)*Tcw1.row(2)-Tcw1.row(0);
	A.row(1) = xn1.at<float>(1)*Tcw1.row(2)-Tcw1.row(1);
	A.row(2) = xn2.at<float>(0)*Tcw2.row(2)-Tcw2.row(0);
	A.row(3) = xn2.at<float>(1)*Tcw2.row(2)-Tcw2.row(1);

	cv::Mat w,u,vt;
	cv::SVD::compute(A,w,u,vt,cv::SVD::MODIFY_A| cv::SVD::FULL_UV);

	*x3D = vt.row(3).t();

	if(x3D->at<float>(3)==0) {
#if DEBUG_OUTLIER
		printf("scale: %f\n",x3D->at<float>(3));
#endif
		return -1;
	}

	// Euclidean coordinates
	*x3D = x3D->rowRange(0,3)/x3D->at<float>(3);
	cv::Mat x3Dt = x3D->t();

	//Check triangulation in front of cameras
	float z1 = Rcw1.row(2).dot(x3Dt)+tcw1.at<float>(2);
	if(z1<=0) {
#if DEBUG_OUTLIER
		printf("z1: %f\n",z1);
#endif
		return -1;
	}

	float z2 = Rcw2.row(2).dot(x3Dt)+tcw2.at<float>(2);
	if(z2<=0) {
#if DEBUG_OUTLIER
		printf("z2: %f\n",z2);
#endif
		return -1;
	}

	//Check reprojection error in first keyframe
	float sigmaSquare1 = vSigma2[kp1.octave];
	float x1 = Rcw1.row(0).dot(x3Dt)+tcw1.at<float>(0);
	float y1 = Rcw1.row(1).dot(x3Dt)+tcw1.at<float>(1);
	float invz1 = 1.0/z1;
	float u1 = fx*x1*invz1+cx;
	float v1 = fy*y1*invz1+cy;
	float errX1 = u1 - kp1.pt.x;
	float errY1 = v1 - kp1.pt.y;
	float err1 = (errX1*errX1+errY1*errY1);
	if(err1 >5*sigmaSquare1) {
#if DEBUG_OUTLIER
		printf("err1: %f\n",err1);
#endif
		return -1;
	}

	//Check reprojection error in second keyframe
	float sigmaSquare2 = vSigma2[kp2.octave];
	float x2 = Rcw2.row(0).dot(x3Dt)+tcw2.at<float>(0);
	float y2 = Rcw2.row(1).dot(x3Dt)+tcw2.at<float>(1);
	float invz2 = 1.0/z2;
	float u2 = fx*x2*invz2+cx;
	float v2 = fy*y2*invz2+cy;
	float errX2 = u2 - kp2.pt.x;
	float errY2 = v2 - kp2.pt.y;
	float err2 = (errX2*errX2+errY2*errY2);
	if(err2>5*sigmaSquare2) {
#if DEBUG_OUTLIER
		printf("err2: %f\n",err2);
#endif
		return -1;
	}

	//Check scale consistency
	cv::Mat normal1 = *x3D-Ow1;
	float dist1 = cv::norm(normal1);

	cv::Mat normal2 = *x3D-Ow2;
	float dist2 = cv::norm(normal2);

	if(dist1==0 || dist2==0) {
#if DEBUG_OUTLIER
		printf("dist1: %f dist2: %f\n",dist1,dist2);
#endif
		return -1;
	}

	float ratioDist = dist1/dist2;
	float ratioOctave = vScaleFactor[kp1.octave]/vScaleFactor[kp2.octave];
	if(ratioDist*ratioFactor<ratioOctave || ratioDist>ratioOctave*ratioFactor) {
#if DEBUG_OUTLIER
		printf("ratioDist: %f ratioOctave: %f\n",ratioDist,ratioOctave);
#endif
		return -1;
	}

	*depth = dist1;
	return err1 + err2;
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./linear_triangulate pose_stamped.txt cvkey.match map_point.txt\n");
		return 1;
	}
	srand(time(NULL));
	Rcl = (cv::Mat_<float>(3,3) << 1,0,0,0,0,-1,0,1,0);
	Tcl = (cv::Mat_<float>(3,1) << 0,0.06,0);
	for (int i=0;i<ORB_levels;i++) {
		if (i==0) {
			vScaleFactor.push_back(1);
			vSigma2.push_back(1);
		} else {
			vScaleFactor.push_back(vScaleFactor[i-1] * ORB_scale_factor);
			vSigma2.push_back(vScaleFactor[i] * vScaleFactor[i]);
		}
	}

	//read pose file
	FILE* pose_stamped = fopen(argv[1],"r");
	if (!pose_stamped)
		return 1;
	char buffer[1024];
	std::vector<cv::Mat> rotations;
	std::vector<cv::Mat> translations;
	while (fgets(buffer,1024,pose_stamped)) {
		double t,x,y,z,qx,qy,qz,qw;
		if (sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",&t,&x,&y,&z,&qw,&qx,&qy,&qz)==8) {
			double r[9];
			quaternionToRotation(qx,qy,qz,qw,r);
			cv::Mat Rwl = (cv::Mat_<float>(3,3) << r[0],r[1],r[2],r[3],r[4],r[5],r[6],r[7],r[8]);
			cv::Mat Twl = (cv::Mat_<float>(3,1) << x,y,z);
			rotations.push_back(Rcl * Rwl.t());
			translations.push_back(- Rcl * Rwl.t() * Twl + Tcl);
		} else {
			printf("Error parsing: %s\n",buffer);
		}
	}
	fclose(pose_stamped);

	struct timespec start,end;
	clock_gettime(CLOCK_MONOTONIC,&start);
	int count_points = 0;
	int cjIterations = 0;
	float RMSE = 0;
	int numInliers = 0;
	int numOutliers = 0;
	std::vector<cv::Mat> map_position;
	std::vector<float> map_error;
	std::vector<float> depth;
	float lidar_depth = -1;
	FILE* key_match = fopen(argv[2],"r");
	FILE* map_point = fopen(argv[3],"w");
	if (!(key_match && map_point))
		return 1;
	while (fgets(buffer,1024,key_match)) {
		float u1,v1,u2,v2;
		int o1,o2;
		if (sscanf(buffer,"%f %f %d %f %f %d",&u1,&v1,&o1,&u2,&v2,&o2)==6) {
			cv::KeyPoint kp1(u1,v1,1,-1,0,o1,-1);
			cv::KeyPoint kp2(u2,v2,1,-1,0,o2,-1);
			cv::Mat x3D;
			float d;
			float err = linear_triangulate(kp1,kp2,rotations[0],translations[0],rotations[1],translations[1],&x3D,&d);
			if (err < 0) {
				numOutliers++;
			} else {
				numInliers++;
				map_position.push_back(x3D);
				map_error.push_back(err);
				depth.push_back(d);
				RMSE += err;
			}
		}
		count_points++;
	}
	for (size_t i=0;i<map_position.size();i++) {
		fprintf(map_point,"%f %f %f %d %f\n",map_position[i].at<float>(0,0),map_position[i].at<float>(1,0),map_position[i].at<float>(2,0),2,map_error[i]);
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double dt = end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec);
	RMSE = sqrt(RMSE / map_error.size());
	std::sort(depth.begin(),depth.end());
	printf("Optimized %d map points (%d iter, %fs, RMSE = %f)\n",count_points, cjIterations, dt, RMSE);
	printf("triangulation: %d inliers %d outliers\n",numInliers,numOutliers);
	printf("depth max: %f median: %f lidar: %f (f = %f)\n",*std::max_element(depth.begin(),depth.end()),depth[depth.size()/2],lidar_depth,fx);
	fclose(key_match);
	fclose(map_point);
}
