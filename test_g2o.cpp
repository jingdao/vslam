#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <iostream>

#include "g2o/core/eigen_types.h"
//#include "g2o/stuff/sampler.h"
//#include "g2o/stuff/command_args.h"
#include "g2o/core/sparse_optimizer.h"
#include "g2o/core/block_solver.h"
#include "g2o/core/solver.h"
#include "g2o/core/optimization_algorithm_levenberg.h"
//#include "g2o/core/optimization_algorithm_gauss_newton.h"
#include "g2o/core/base_vertex.h"
#include "g2o/core/base_unary_edge.h"
#include "g2o/solvers/dense/linear_solver_dense.h"
//#include "g2o/solvers/csparse/linear_solver_csparse.h"

//http://dspguru.com/dsp/howtos/how-to-generate-white-gaussian-noise
//Box-Mueller transform
double randn() {
	static double z0,z1;
	static bool generate = false;
	generate = !generate;
	if (generate) {
		double U1,U2,V1,V2,S;
		do {
			U1 = rand() * (1.0 / RAND_MAX);
			U2 = rand() * (1.0 / RAND_MAX);
			V1 = 2 * U1 - 1;
			V2 = 2 * U2 - 1;
			S = V1 * V1 + V2 * V2;
		} while (S >=1);

		z0 = sqrt(-2 * log(S) / S) * V1;
		z1 = sqrt(-2 * log(S) / S) * V2;
		return z0;
	} else {
		return z1;
	}
}

class VertexLine : public g2o::BaseVertex<2, Eigen::Vector2d> {
	public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	VertexLine() {}

	virtual bool read(std::istream& /*is*/) {return false;}
	virtual bool write(std::ostream& /*os*/) const {return false;}
	virtual void setToOriginImpl() {}
	virtual void oplusImpl(const double* update) {
		Eigen::Vector2d::ConstMapType v(update);
		_estimate += v;
	}
};

class EdgeLine : public g2o::BaseUnaryEdge<1, Eigen::Vector2d, VertexLine> {
	public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
	EdgeLine() {}
	virtual bool read(std::istream& /*is*/) {return false;}
	virtual bool write(std::ostream& /*os*/) const {return false;}

	void computeError() {
		const VertexLine* line = static_cast<const VertexLine*>(vertex(0));
		const double& m = line->estimate()(0);
		const double& c = line->estimate()(1);
		const double& x = measurement()(0);
		const double& y = measurement()(1);
		_error(0) = y - (m*x + c);
	}
};

int main() {
	//initialize solver
	g2o::SparseOptimizer optimizer;
	g2o::BlockSolverX::LinearSolverType *linearSolver = new g2o::LinearSolverDense<g2o::BlockSolverX::PoseMatrixType>();
	g2o::BlockSolverX *blocksolver = new g2o::BlockSolverX(linearSolver);
	g2o::OptimizationAlgorithmLevenberg* solver = new g2o::OptimizationAlgorithmLevenberg(blocksolver);
	optimizer.setAlgorithm(solver);
	srand(time(NULL));

	//add vertex
	VertexLine *v = new VertexLine();
	v->setId(0);
	v->setEstimate(Eigen::Vector2d(0.0,2.0));
	optimizer.addVertex(v);
	//add edges
	for (int i=0;i<100;i++) {
		double x = i;
		double y = 2 * x + 3 + randn();
		EdgeLine* e = new EdgeLine();
		e->setInformation(Eigen::Matrix<double,1,1>::Identity());
		e->setVertex(0,v);
		e->setMeasurement(Eigen::Vector2d(x,y));
		optimizer.addEdge(e);
	}

	//optimize
	int nIter = 10;
	optimizer.initializeOptimization();
	optimizer.setVerbose(true);
	optimizer.optimize(nIter);
	printf("Optimizing for %d iterations\n",nIter);

	//get result
	printf("m: %f c:%f\n",v->estimate()(0),v->estimate()(1));

	return 0;
}
