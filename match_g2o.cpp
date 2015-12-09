#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <vector>
#include "g2o/core/eigen_types.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
#include "g2o/core/jacobian_workspace.h"

typedef Eigen::Matrix<double,3,3,Eigen::RowMajor> Mat3;
double fx = 971.760406;
double fy = 971.138862;
double cx = 319.500000;
double cy = 239.500000;
//lidar to camera transformation
Mat3 Rcl;
Eigen::Vector3d Tcl(0,0,0);


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

class EdgeProjection : public g2o::BaseUnaryEdge<1, Eigen::Vector2d, VertexMapPoint> {
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
		_error(0) = 0;
		_error(0) += (u - measurement()(0)) * (u - measurement()(0));
		_error(0) += (v - measurement()(1)) * (v - measurement()(1));
	}

//	void linearizeOplus() {
//		_jacobianOplusXi(0,0) = - measurement()(0);
//		_jacobianOplusXi(0,1) = -1;
//	}
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

	//read pose file
	FILE* pose_stamped = fopen(argv[1],"r");
	if (!pose_stamped)
		return 1;
	char buffer[1024];
	std::vector<Mat3> rotations;
	std::vector<Eigen::Vector3d> translations;
	while (fgets(buffer,1024,pose_stamped)) {
		double t,x,y,z,qx,qy,qz,qw;
		if (sscanf(buffer,"%lf %lf %lf %lf %lf %lf %lf %lf",&t,&x,&y,&z,&qx,&qy,&qz,&qw)==8) {
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
	FILE* key_match = fopen(argv[2],"r");
	FILE* map_point = fopen(argv[3],"w");
	if (!(key_match && map_point))
		return 1;
	while (fgets(buffer,1024,key_match)) {
		//initialize solver
		g2o::SparseOptimizer optimizer;
		g2o::BlockSolverX::LinearSolverType *linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
		g2o::BlockSolverX *blocksolver = new g2o::BlockSolverX(linearSolver);
		g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blocksolver);
		optimizer.setAlgorithm(solver);

		//add vertex
		VertexMapPoint *vt = new VertexMapPoint();
		vt->setId(0);
		vt->setEstimate(Eigen::Vector3d(1.0 / RAND_MAX * rand(), 1.0 / RAND_MAX * rand(), 1.0 / RAND_MAX * rand()));
		optimizer.addVertex(vt);

		//add edges
		int id;
		double u,v;
		char* tok = strtok(buffer," ");
		while (tok) {
			id = atoi(tok);
			tok = strtok(NULL," \n");
			u = atof(tok);
			tok = strtok(NULL," \n");
			v = atof(tok);
			tok = strtok(NULL," \n");
			EdgeProjection* e = new EdgeProjection();
			e->setInformation(Eigen::Matrix<double,1,1>::Identity());
			e->setVertex(0,vt);
			e->setMeasurement(Eigen::Vector2d(u-cx,v-cy));
			e->Rcw = rotations[id];
			e->Tcw = translations[id];
			optimizer.addEdge(e);
		}

		//optimize
		int nIter = 10;
		optimizer.initializeOptimization();
		optimizer.setVerbose(false);
		optimizer.optimize(nIter);

		//record result
		fprintf(map_point,"%f %f %f\n",vt->estimate()(0),vt->estimate()(1),vt->estimate()(2));
		count_points++;
	}
	clock_gettime(CLOCK_MONOTONIC,&end);
	double dt = end.tv_sec - start.tv_sec + 0.000000001 * (end.tv_nsec - start.tv_nsec);
	printf("Optimized %d map points (%fs)\n",count_points, dt);
	fclose(key_match);
	fclose(map_point);

	return 0;
}
