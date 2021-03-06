#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include "g2o/core/eigen_types.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/core/jacobian_workspace.h"
#include "g2o/core/robust_kernel_impl.h"

#define DEBUG_SINGLE 0

typedef Eigen::Matrix<double,3,3,Eigen::RowMajor> Mat3;
double fx = 567;
double fy = 567;
double cx = 319.500000;
double cy = 239.500000;
int numIterations = 10;
int numSeeds = 20;
double huber_threshold = -1;
//double error_threshold = 5;
double error_threshold = 10;
//double huber_threshold = sqrt(5.991);
//lidar to camera transformation
Mat3 Rcl;
Eigen::Vector3d Tcl;
std::vector<Eigen::Vector3d> lidarcloud;

class VertexMapPoint : public g2o::BaseVertex<3, Eigen::Vector3d> {
	public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	VertexMapPoint() {}

	virtual bool read(std::istream& /*is*/) {return false;}
	virtual bool write(std::ostream& /*os*/) const {return false;}
	virtual void setToOriginImpl() {}
	virtual void oplusImpl(const double* update) {
		Eigen::Vector3d::ConstMapType v(update);
		_estimate += v;
	}
};

class EdgeProjection : public g2o::BaseUnaryEdge<2, Eigen::Vector2d, VertexMapPoint> {
	public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	EdgeProjection() {}
	Mat3 Rcw;
	Eigen::Vector3d Tcw;
	virtual bool read(std::istream& /*is*/) {return false;}
	virtual bool write(std::ostream& /*os*/) const {return false;}

	void computeError() {
		const VertexMapPoint* mp = static_cast<const VertexMapPoint*>(vertex(0));
		Eigen::Vector3d xc = Rcw * mp->estimate() + Tcw;
		double u = fx * xc(0) / xc(2);
		double v = fy * xc(1) / xc(2);
		_error(0) = u - measurement()(0);
		_error(1) = v - measurement()(1);
	}

	void linearizeOplus() {
		const VertexMapPoint* mp = static_cast<const VertexMapPoint*>(vertex(0));
		Eigen::Vector3d xc = Rcw * mp->estimate() + Tcw;
		_jacobianOplusXi(0,0) = - fx / xc(2) / xc(2) * (Rcw(2,0) * xc(0) - Rcw(0,0) * xc(2));
		_jacobianOplusXi(0,1) = - fx / xc(2) / xc(2) * (Rcw(2,1) * xc(0) - Rcw(0,1) * xc(2));
		_jacobianOplusXi(0,2) = - fx / xc(2) / xc(2) * (Rcw(2,2) * xc(0) - Rcw(0,2) * xc(2));
		_jacobianOplusXi(1,0) = - fy / xc(2) / xc(2) * (Rcw(2,0) * xc(1) - Rcw(1,0) * xc(2));
		_jacobianOplusXi(1,1) = - fy / xc(2) / xc(2) * (Rcw(2,1) * xc(1) - Rcw(1,1) * xc(2));
		_jacobianOplusXi(1,2) = - fy / xc(2) / xc(2) * (Rcw(2,2) * xc(1) - Rcw(1,2) * xc(2));
	}
};

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

void fixScale(std::vector<Eigen::Vector3d> *scan, std::vector<Eigen::Vector3d> *map, Mat3 Rcw, Eigen::Vector3d Tcw) {
	Eigen::Vector3d midpoint = (*map)[0];
	for (size_t i=1;i<map->size();i++) {
		if (fabs((*map)[i](2)) < fabs(midpoint(2)))
			midpoint = (*map)[i];
	}
	Eigen::Vector3d mc = Rcw * midpoint + Tcw;
	Eigen::Vector3d ref = Rcw * (*scan)[0] + Tcw;
	double minAngle = mc.cross(ref).norm() / ref.norm();
	for (size_t j=1;j<scan->size();j++) {
		Eigen::Vector3d rc = Rcw * (*scan)[j] + Tcw;
		double angle = mc.cross(rc).norm() / rc.norm();
		if (angle < minAngle) {
			ref = rc;
			minAngle = angle;
		}
	}
	double scale = ref.norm() / mc.norm();
	for (size_t i=0;i<map->size();i++) {
		Eigen::Vector3d xc = Rcw * (*map)[i] + Tcw;
		xc *= scale;
		(*map)[i] = Rcw.transpose() * (xc - Tcw); 
	}
	printf("Fixed scale by %f\n",scale);
}

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./match_g2o pose_stamped.txt key.match map_point.txt [-m lidar_map.txt -f fx -i nIter]\n");
		return 1;
	}
	srand(time(NULL));
	Rcl << 1,0,0,0,0,-1,0,1,0;
//	Tcl << 0,-0.25,-0.18;
	Tcl << 0,0.06,0;

	//read pose file
	FILE* pose_stamped = fopen(argv[1],"r");
	if (!pose_stamped)
		return 1;
	char buffer[1024];
	std::vector<Mat3> rotations;
	std::vector<Eigen::Vector3d> translations;
	while (fgets(buffer,1024,pose_stamped)) {
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

	for (int i=4;i<argc-1;i++) {
		if (strncmp(argv[i],"-m",2)==0) {
			FILE *f = fopen(argv[++i],"r");
			if (f) {
				while (fgets(buffer,128,f)) {
					float x,y,z;
					if (sscanf(buffer,"%f %f %f",&x,&y,&z)==3) {
						Eigen::Vector3d v(x,y,z);
						lidarcloud.push_back(v);
					}
				}
				fclose(f);
			}
		} else if (strncmp(argv[i],"-f",2)==0) {
			fx = atof(argv[++i]);
			fy = fx;
		} else if (strncmp(argv[i],"-i",2)==0) {
			numIterations = atoi(argv[++i]);
		}
	}

	struct timespec start,end;
	clock_gettime(CLOCK_MONOTONIC,&start);
	int count_points = 0;
	int cjIterations = 0;
	double RMSE = 0;
	int numInliers = 0;
	int numOutliers = 0;
	std::vector<Eigen::Vector3d> map_position;
	std::vector< std::vector<int> > map_index;
	std::vector<double> map_error;
	std::vector<double> depth;
	double lidar_depth = -1;
	FILE* key_match = fopen(argv[2],"r");
	FILE* map_point = fopen(argv[3],"w");
	if (!(key_match && map_point))
		return 1;
	while (fgets(buffer,1024,key_match)) {
#if DEBUG_SINGLE
		printf("key.match: %s",buffer);
#endif
		//initialize solver
		g2o::SparseOptimizer optimizer;
		g2o::BlockSolverX::LinearSolverType *linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
		g2o::BlockSolverX *blocksolver = new g2o::BlockSolverX(linearSolver);
		g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blocksolver);
		optimizer.setAlgorithm(solver);
		optimizer.setVerbose(false);

		//add vertex
		VertexMapPoint *vt = new VertexMapPoint();
		vt->setId(0);
		optimizer.addVertex(vt);

		//add edges
		int id;
		double u,v;
		char* tok = strtok(buffer," ");
		std::vector<int> index;
		while (tok) {
			id = atoi(tok);
			index.push_back(id);
			tok = strtok(NULL," \n");
			u = atof(tok);
			tok = strtok(NULL," \n");
			v = atof(tok);
			tok = strtok(NULL," \n");
			EdgeProjection* e = new EdgeProjection();
			e->setInformation(Eigen::Matrix<double,2,2>::Identity());
			e->setVertex(0,vt);
			e->setMeasurement(Eigen::Vector2d(u-cx,v-cy));
			e->Rcw = rotations[id];
			e->Tcw = translations[id];
			if (huber_threshold > 0) {
				g2o::RobustKernelHuber* rk = new g2o::RobustKernelHuber;
				e->setRobustKernel(rk);
				rk->setDelta(huber_threshold);
			}
			optimizer.addEdge(e);
		}

		//optimize
		Eigen::Vector3d bestEstimate;
		double leastError;
		for (int i=0;i<numSeeds;i++) {
			vt->setEstimate(Eigen::Vector3d(1.0 / RAND_MAX * rand(), 1.0 / RAND_MAX * rand(), 1.0 / RAND_MAX * rand()));
			optimizer.initializeOptimization();
			cjIterations += optimizer.optimize(numIterations);
			if (i==0 || (optimizer.activeChi2() < leastError && vt->estimate()(1) > 0)) {
				bestEstimate = vt->estimate();
				leastError = optimizer.activeChi2();
			}
		}
		
		if (leastError < error_threshold && bestEstimate(1) > 0) {
			numInliers++;
		} else {
			numOutliers++;
			continue;
		}

		//record result
		RMSE += leastError / index.size();
		depth.push_back(bestEstimate(1));
#if DEBUG_SINGLE
		printf("reprojection: ");
		char* c = buffer;
		c += sprintf(c,"transformation:\n");
		for (unsigned int i=0;i<index.size();i++) {
			id = index[i];
			Eigen::Vector3d xc = rotations[id] * bestEstimate + translations[id];
			double u = fx * xc(0) / xc(2);
			double v = fy * xc(1) / xc(2);
			printf("%d %f %f ",id,u+cx,v+cy);
			c += sprintf(c,"%4.2f %4.2f %4.2f\n%4.2f %4.2f %4.2f\n%4.2f %4.2f %4.2f\n",
						rotations[id](0,0),rotations[id](0,1),rotations[id](0,2),
						rotations[id](1,0),rotations[id](1,1),rotations[id](1,2),
						rotations[id](2,0),rotations[id](2,1),rotations[id](2,2));
			c += sprintf(c,"[%4.2f %4.2f %4.2f]\n",translations[id](0),translations[id](1),translations[id](2));
		}
		printf("\n%s",buffer);
		printf("estimate: %f %f %f %lu %f\n",bestEstimate(0),bestEstimate(1),bestEstimate(2),index.size(),leastError);
#else
//		fprintf(map_point,"%f %f %f %lu %f\n",bestEstimate(0),bestEstimate(1),bestEstimate(2),index.size(),leastError);
		map_position.push_back(bestEstimate);
		map_index.push_back(index);
		map_error.push_back(leastError);

#endif
		count_points++;
	}
	if (lidarcloud.size() > 0) {
//		fixScale(&lidarcloud,&map_position,rotations[0],translations[0]);
		std::vector<double> dv;
		for (size_t i=0;i<lidarcloud.size();i++) {
			if (lidarcloud[i](1) > 0)
				dv.push_back(lidarcloud[i](1));
		}
		std::sort(dv.begin(),dv.end());
		lidar_depth = dv[dv.size()/2];
	}
	for (size_t i=0;i<map_position.size();i++) {
		fprintf(map_point,"%f %f %f %lu %f\n",map_position[i](0),map_position[i](1),map_position[i](2),map_index[i].size(),map_error[i]);
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double dt = end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec);
	RMSE = sqrt(RMSE / map_position.size());
	std::sort(depth.begin(),depth.end());
	printf("Optimized %d map points (%d iter, %fs, RMSE = %f)\n",count_points, cjIterations, dt, RMSE);
	printf("triangulation: %d inliers %d outliers\n",numInliers,numOutliers);
	printf("depth max: %f median: %f lidar: %f (f = %f)\n",*std::max_element(depth.begin(),depth.end()),depth[depth.size()/2],lidar_depth,fx);
	fclose(key_match);
	fclose(map_point);

	return 0;
}
