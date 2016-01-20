#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include <iostream>
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
double fx = 971.760406;
double fy = 971.138862;
double cx = 319.500000;
double cy = 239.500000;
int numIterations = 10;
int numSeeds = 10;
double huber_threshold = -1;
//double huber_threshold = sqrt(5.991);
//lidar to camera transformation
Mat3 Rcl;
Eigen::Vector3d Tcl;


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
		double u = - fx * xc(0) / xc(2);
		double v = - fy * xc(1) / xc(2);
		_error(0) = u - measurement()(0);
		_error(1) = v - measurement()(1);
	}

	void linearizeOplus() {
		const VertexMapPoint* mp = static_cast<const VertexMapPoint*>(vertex(0));
		Eigen::Vector3d xc = Rcw * mp->estimate() + Tcw;
		_jacobianOplusXi(0,0) = fx / xc(2) / xc(2) * (Rcw(2,0) * xc(0) - Rcw(0,0) * xc(2));
		_jacobianOplusXi(0,1) = fx / xc(2) / xc(2) * (Rcw(2,1) * xc(0) - Rcw(0,1) * xc(2));
		_jacobianOplusXi(0,2) = fx / xc(2) / xc(2) * (Rcw(2,2) * xc(0) - Rcw(0,2) * xc(2));
		_jacobianOplusXi(1,0) = fy / xc(2) / xc(2) * (Rcw(2,0) * xc(1) - Rcw(1,0) * xc(2));
		_jacobianOplusXi(1,1) = fy / xc(2) / xc(2) * (Rcw(2,1) * xc(1) - Rcw(1,1) * xc(2));
		_jacobianOplusXi(1,2) = fy / xc(2) / xc(2) * (Rcw(2,2) * xc(1) - Rcw(1,2) * xc(2));
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

int main(int argc,char* argv[]) {
	if (argc < 4) {
		printf("./match_g2o pose_stamped.txt key.match map_point.txt\n");
		return 1;
	}
	srand(time(NULL));
	Rcl << 1,0,0,0,0,1,0,-1,0;
	Tcl << 0,0.25,0.18;

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

	struct timespec start,end;
	clock_gettime(CLOCK_MONOTONIC,&start);
	int count_points = 0;
	int cjIterations = 0;
	double RMSE = 0;
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
			e->setMeasurement(Eigen::Vector2d(u-cx,cy-v));
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
			if (i==0 || optimizer.activeChi2() < leastError) {
				bestEstimate = vt->estimate();
				leastError = optimizer.activeChi2();
			}
		}
		

		//record result
		RMSE += leastError / index.size();
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
#else
		fprintf(map_point,"%f %f %f %lu %f\n",bestEstimate(0),bestEstimate(1),bestEstimate(2),index.size(),leastError);
#endif
		count_points++;
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double dt = end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec);
	RMSE = sqrt(RMSE / count_points);
	printf("Optimized %d map points (%d iter, %fs, RMSE = %f)\n",count_points, cjIterations, dt, RMSE);
	fclose(key_match);
	fclose(map_point);

	return 0;
}
