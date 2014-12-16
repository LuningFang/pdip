//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010 Alessandro Tasora
// Copyright (c) 2013 Project Chrono
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

#ifndef CHLCPINTERIORPOINT_H
#define CHLCPINTERIORPOINT_H

//////////////////////////////////////////////////
//
//   ChLcpInteriorPoint.h
//
//   Lulu is a cutie hahaha
// ------------------------------------------------
//             www.deltaknowledge.com
// ------------------------------------------------
///////////////////////////////////////////////////



#include "ChLcpIterativeSolver.h"


namespace chrono
{

/// An iterative LCP solver based on projective
/// fixed point method, with overrelaxation
/// and immediate variable update as in SOR methods.
/// The problem is described by a variational inequality VI(Z*x-d,K):
///
///  | M -Cq'|*|q|- | f|= |0| , l \in Y, C \in Ny, normal cone to Y
///  | Cq -E | |l|  |-b|  |c|
///
/// Also Z symmetric by flipping sign of l_i: |M  Cq'|*| q|-| f|=|0|
///                                           |Cq  E | |-l| |-b| |c|
/// * case linear problem:  all Y_i = R, Ny=0, ex. all bilaterals
/// * case LCP: all Y_i = R+:  c>=0, l>=0, l*c=0
/// * case CCP: Y_i are friction cones

class ChApi ChLcpInteriorPoint : public ChLcpIterativeSolver
{
protected:
			//
			// DATA
			//


public:
			//
			// CONSTRUCTORS
			//

	ChLcpInteriorPoint(
				int mmax_iters=50,      ///< max.number of iterations
				bool mwarm_start=false,	///< uses warm start?
				double mtolerance=0.0,  ///< tolerance for termination criterion
				double momega=1.0       ///< overrelaxation criterion
				)
			: ChLcpIterativeSolver(mmax_iters,mwarm_start, mtolerance,momega)
			{};

	virtual ~ChLcpInteriorPoint() {};

			//
			// FUNCTIONS
			//

				/// Performs the solution of the LCP.
				/// return  the maximum constraint violation after termination.
	virtual double Solve(
				ChLcpSystemDescriptor& sysd		///< system description with constraints and variables
				);



	ChMatrixDynamic<double> opt_r;
	ChMatrixDynamic<double> ff;
	ChMatrixDynamic<double> lambda_k;
	ChMatrixDynamic<double> xk;
	ChMatrixDynamic<double> r_dual;
	ChMatrixDynamic<double> r_cent;
	ChMatrixDynamic<double> d_x;
	ChMatrixDynamic<double> d_lambda;
	ChMatrixDynamic<double> Schur_rhs;
	ChMatrixDynamic<double> grad_f;

	// tmp vectors for line search
	ChMatrixDynamic<double> ff_tmp;
	ChMatrixDynamic<double> lambda_k_tmp;
	ChMatrixDynamic<double> xk_tmp;
	ChMatrixDynamic<double> x_old_tmp;
	ChMatrixDynamic<double> r_dual_tmp;
	ChMatrixDynamic<double> r_cent_tmp;




	// evaluate contraints
private:
	void evaluateConstraints(double* fric, int nc, bool tmp);
	void computeSchurRHS(double* grad_f, double *fric, int nc, double t);
	void computeSchurKKT(double* grad_f, double* fric, int nc, double t, bool tmp);
	void computeBlockDiagonal(double* val_diagonal, double *fric, int nc, double t);
	void computeSpikeRHS(double *rhs, double *fric, int nc, double t);


//	std::vector<int> findModifiedIndex(int nc, double res);


};



} // END_OF_NAMESPACE____




#endif  // END of ChLcpIterativeSOR.h
