//
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2010-2012 Alessandro Tasora
// All rights reserved.
//
// Use of this source code is governed by a BSD-style license that can be
// found in the LICENSE file at the top level of the distribution
// and at http://projectchrono.org/license-chrono.txt.
//

///////////////////////////////////////////////////
//
//   ChLcpInteriorPoint.cpp
//
//
//    file for CHRONO HYPEROCTANT LCP solver
//
// ------------------------------------------------
//             www.deltaknowledge.com
// ------------------------------------------------
///////////////////////////////////////////////////


#include "ChLcpInteriorPoint.h"
#include "ChLcpConstraintTwoFrictionT.h"
#include "pdip_solve.h"

namespace chrono
{

double ChLcpInteriorPoint::Solve(
		ChLcpSystemDescriptor& sysd		///< system description with constraints and variables
)
{
	std::cout << "-------using interior point solver!!------" << std::endl;
	std::vector<ChLcpConstraint*>& mconstraints = sysd.GetConstraintsList();
	std::vector<ChLcpVariables*>&  mvariables	= sysd.GetVariablesList();



	ChMatrixDynamic <double> mv0;
	ChSparseMatrix mM;
	ChSparseMatrix mCq;
	ChSparseMatrix mE;
	ChMatrixDynamic <double> mf;
	ChMatrixDynamic <double> mb;
	ChMatrixDynamic <double> mfric;

	sysd.ConvertToMatrixForm(&mCq, &mM, &mE, &mf, &mb, &mfric);
	sysd.FromVariablesToVector(mv0);
	
/*		ChStreamOutAsciiFile file_V0( "dump_V_old.dat" ) ;
		mv0.StreamOUTdenseMatlabFormat(file_V0) ;
		ChStreamOutAsciiFile file_M ( "dump_M.dat" ) ;
		mM.StreamOUTsparseMatlabFormat ( file_M ) ;

		ChStreamOutAsciiFile file_Cq ( "dump_Cq.dat" ) ;
		mCq.StreamOUTsparseMatlabFormat ( file_Cq ) ;

		ChStreamOutAsciiFile file_E ( "dump_E.dat" ) ;
		mE.StreamOUTsparseMatlabFormat ( file_E ) ;

		ChStreamOutAsciiFile file_f ( "dump_f.dat" ) ;
		mf.StreamOUTdenseMatlabFormat ( file_f ) ;

		ChStreamOutAsciiFile file_b ( "dump_b.dat" ) ;
		mb.StreamOUTdenseMatlabFormat ( file_b ) ;

		ChStreamOutAsciiFile file_fric ( "dump_fric.dat" ) ;
		mfric.StreamOUTdenseMatlabFormat ( file_fric ) ;

		printf("Successfully writing files!\n");
*/
/*
			file_f.GetFstream().close();
			file_fric.GetFstream().close();
			file_V0.GetFstream().close();
			file_M.GetFstream().close();
			file_Cq.GetFstream().close();
			file_b.GetFstream().close();

*/	
	int nBodies = mM.GetColumns()/6;
	size_t nVariables = mvariables.size();
	size_t nConstraints = sysd.CountActiveConstraints();
	int numContacts = nConstraints/3;
	size_t nOfConstraints = mconstraints.size();

	/* ALWYAS DO THIS IN THE LCP SOLVER!!!*/
	for (unsigned int ic = 0; ic < nConstraints; ic++)
		mconstraints[ic]->Update_auxiliary();
	//Get sparse info for contact jacobian and Minv matrix to pass on to Ang's solver
	std::vector<int> index_i_Cq;
	std::vector<int> index_j_Cq;
	std::vector<double> val_Cq;
	double val;
//	fprintf(stderr, "------------Cq(from C::E)----------\n");
	for (int ii = 0; ii < mCq.GetRows(); ii++){
		for (int jj = 0; jj < mCq.GetColumns(); jj++){
			val = mCq.GetElement(ii,jj);
			if (val){
				index_i_Cq.push_back(jj);
				index_j_Cq.push_back(ii);
				val_Cq.push_back(val);
//				fprintf(stderr, "%d %d %.20g\n", ii, jj, val);
			}
		}
	}

	

/*	for (int iv = 0; iv < mvariables.size(); iv++)
		if (mvariables[iv]->IsActive())
			mvariables[iv]->Compute_invMb_v(mvariables[iv]->Get_qb(), mvariables[iv]->Get_fb());
*/

//	int count = 0;
//	for (std::vector<int>::iterator it = index_i_Cq.begin(); it != index_i_Cq.end(); it ++){
//		std::cout << "(" << index_i_Cq[count] <<"," << index_j_Cq[count] <<"):" << val_Cq[count] << std::endl;
//		count ++;
//	}

	// Minv matrix
	std::vector<int> index_i_Minv;
	std::vector<int> index_j_Minv;
	std::vector<double> val_Minv;
	for (int i = 0; i < nBodies*6; i++){
		index_i_Minv.push_back(i);
		index_j_Minv.push_back(i);
		val_Minv.push_back(1.0/mM.GetElement(i,i));
	}

	// create reference to pass on to SPIKE
	int *Cq_i = &index_i_Cq[0];
	int *Cq_j = &index_j_Cq[0];
	int Cq_nnz = val_Cq.size();
	double *Cq_val = &val_Cq[0];

	int *Minv_i = &index_i_Minv[0];
	int *Minv_j = &index_j_Minv[0];
	double *Minv_val = &val_Minv[0];

	// formulate rhs of optimization problem f(x) = 1/2 *x'*N*x + r'*x
	ChMatrixDynamic<double> opt_r_tmp(nConstraints,1);
	ChMatrixDynamic<double> k(6*nBodies,1);
	ChMatrixDynamic<double> Mk(6*nBodies,1);
	
	int lu = 0;
	// assemble r vector
	/** 1. put [M^-1]*k in q sparse vector of each variable **/
	for (unsigned int iv = 0; iv < nVariables; iv ++)
		if (mvariables[iv]->IsActive()){
			mvariables[iv]->Compute_invMb_v(mvariables[iv]->Get_qb(), mvariables[iv]->Get_fb());
			ChMatrixDynamic <double> tmp1 = mvariables[iv]->Get_fb();
			ChMatrixDynamic <double> tmp2 = mvariables[iv]->Get_qb();
			k(6*lu + 0,0) = tmp1(0, 0);
			k(6*lu + 1,0) = tmp1(1, 0);
			k(6*lu + 2,0) = tmp1(2, 0);
			k(6*lu + 3,0) = tmp1(3, 0);
			k(6*lu + 4,0) = tmp1(4, 0);
			k(6*lu + 5,0) = tmp1(5, 0);
			Mk(6*lu + 0, 0) = tmp2(0, 0);
			Mk(6*lu + 1, 0) = tmp2(1, 0);
			Mk(6*lu + 2, 0) = tmp2(2, 0);
			Mk(6*lu + 3, 0) = tmp2(3, 0);
			Mk(6*lu + 4, 0) = tmp2(4, 0);
			Mk(6*lu + 5, 0) = tmp2(5, 0);
			lu ++;
//			fprintf(stderr, "Body %d k: %.12f %.12f %.12f\n", iv,  k(0,0), k(1,0), k(2,0));
//			fprintf(stderr, "Body %d M^[-1]*k: %.12f %.12f %.12f\n", iv,  Mk(0,0), Mk(1,0), Mk(2,0));
	}
	/** 2. now do rhs = D'*q = D'*(M^-1)*k **/
	int s_i = 0;
	opt_r.Resize(nConstraints,1);
	for (unsigned int ic = 0; ic < nConstraints; ic ++)
		if (mconstraints[ic]->IsActive()){
			opt_r(s_i,0) = mconstraints[ic]->Compute_Cq_q();
			++s_i;	
		}

	/** 3.  rhs = rhs + c **/
	sysd.BuildBiVector(opt_r_tmp);
	opt_r.MatrInc(opt_r_tmp);
	


	//printChMatrix(opt_r);
	///////////////////
	//velocity update//
	///////////////////
	ChMatrixDynamic<> mq;
	sysd.FromVariablesToVector(mq, true);

	

	////////////////////////////
	//assign solver parameters//
	////////////////////////////
	double barrier_t = 1;
	double eta_hat;
	int numStages = 1000;
	int mu1 = 10;
	double b1 = 0.5;
	double a1 = 0.01;

	// assign vectors here
	ff.Resize(numContacts*2,1);
	lambda_k.Resize(numContacts*2,1); /*initialize lambda_k*/
	xk.Resize(numContacts*3,1);
	r_dual.Resize(numContacts*3,1);
	r_cent.Resize(numContacts*2,1);
	d_x.Resize(numContacts*3,1);
	d_lambda.Resize(numContacts*2,1);
	Schur_rhs.Resize(3*numContacts,1);
	grad_f.Resize(3*numContacts,1);


	if (mconstraints.size() == 0){
			sysd.FromVectorToConstraints(xk);
			sysd.FromVectorToVariables(mq);

			for (size_t ic = 0; ic < mconstraints.size(); ic ++){
				if (mconstraints[ic]->IsActive())
					mconstraints[ic]->Increment_q(mconstraints[ic]->Get_l_i());
			}

		return 1e-8;
	}



	double *BlockDiagonal_val = new double[9*numContacts];
	int *BlockDiagonal_i = new int[9*numContacts];
	int *BlockDiagonal_j = new int[9*numContacts];
	double *spike_rhs = new double[3*numContacts];

	int tmp0, tmp1, tmp2;
	for (int i = 0; i < numContacts; i ++){
		tmp0 = 3*i;
		tmp1 = 3*i+1;
		tmp2 = 3*i+2;
		*(BlockDiagonal_i + 9*i) = tmp0;
		*(BlockDiagonal_i + 9*i+1) = tmp0;
		*(BlockDiagonal_i + 9*i+2) = tmp0;
		*(BlockDiagonal_i + 9*i+3) = tmp1;
		*(BlockDiagonal_i + 9*i+4) = tmp1;
		*(BlockDiagonal_i + 9*i+5) = tmp1;
		*(BlockDiagonal_i + 9*i+6) = tmp2;
		*(BlockDiagonal_i + 9*i+7) = tmp2;
		*(BlockDiagonal_i + 9*i+8) = tmp2;

		*(BlockDiagonal_j + 9*i) = tmp0;
		*(BlockDiagonal_j + 9*i+1) = tmp1;
		*(BlockDiagonal_j + 9*i+2) = tmp2;
		*(BlockDiagonal_j + 9*i+3) = tmp0;
		*(BlockDiagonal_j + 9*i+4) = tmp1;
		*(BlockDiagonal_j + 9*i+5) = tmp2;
		*(BlockDiagonal_j + 9*i+6) = tmp0;
		*(BlockDiagonal_j + 9*i+7) = tmp1;
		*(BlockDiagonal_j + 9*i+8) = tmp2;
	}



	// initialize xk
	for (int i = 0; i < numContacts; i ++){
		xk(3*i, 0) = 1;
		xk(3*i+1, 0) = 0;
		xk(3*i+2, 0) = 0;
	}

	evaluateConstraints(mfric.GetAddress(), numContacts, false);



	//initialize lambda
	for (int i = 0; i < lambda_k.GetRows(); i++)
		lambda_k(i,0) = -1/(barrier_t * ff(i,0));

	/////////////////////////////
	////GO THROUGH EACH STAGE////
	/////////////////////////////
	for (int stage = 0; stage < numStages; stage++){
		eta_hat = - lambda_k.MatrDot(&lambda_k, &ff);
		barrier_t = mu1 * (2*numContacts)/eta_hat;
		// assemble grad_f = N*x + r
		sysd.ShurComplementProduct(grad_f, &xk, 0);
//		fprintf(stderr, "----------N*x----------\n");
//		for (int i = 0; i < grad_f.GetRows(); i++)
//		fprintf(stderr, "%.20f\n", grad_f.GetElementN(i));
		grad_f.MatrInc(opt_r);
//		fprintf(stderr, "----------grad_f----------\n");
//		for (int i = 0; i < grad_f.GetRows(); i++)
//			fprintf(stderr, "%.20f\n", grad_f.GetElementN(i));
	
		// compute r_d and r_c for schur implementation
		computeSchurRHS(grad_f.GetAddress(), mfric.GetAddress(), numContacts, barrier_t);
//		fprintf(stderr, "----------r_dual----------\n");
//		for (int i = 0; i < r_dual.GetRows(); i++)
//			fprintf(stderr, "%.16f\n", r_dual.GetElementN(i));
//		fprintf(stderr, "----------r_cent----------\n");
//		for (int i = 0; i < r_cent.GetRows(); i++)
//			fprintf(stderr, "%.16f\n", r_cent.GetElementN(i));



		// assemble block diagonal matrix
		computeBlockDiagonal(BlockDiagonal_val, mfric.GetAddress(), numContacts, barrier_t);

		// assemble rhs vector for spike solver
		computeSpikeRHS(spike_rhs, mfric.GetAddress(), numContacts, barrier_t);
//		fprintf(stderr, "----------spike_rhs----------\n");
//		for (int i = 0; i < 3*numContacts; i++){
//			fprintf(stderr, "%.16f\n", *(spike_rhs + i));
//		}

		double *spike_dx = new double [3*numContacts];

		//call ang's solver here....
		bool solveSuc = solveSPIKE(nBodies, numContacts, Cq_i, Cq_j, Cq_nnz, Cq_val, Minv_i, Minv_j, Minv_val, BlockDiagonal_i, BlockDiagonal_j, BlockDiagonal_val, spike_dx, spike_rhs, false);
		
		if (solveSuc == false){
			std::cerr << "Solve Failed!" << std::endl;
			return 0;
		}
		// assume d_x is calculated perfectly!
		for (int i = 0; i < numContacts; i++){
			d_x(3*i,0) = *(spike_dx + 3*i);
			d_x(3*i+1,0) = *(spike_dx + 3*i + 1);
			d_x(3*i+2,0) = *(spike_dx + 3*i + 2);
		}
	

//		fprintf(stderr, "----d_x-------\n");
//		printChMatrix(d_x);
/*		fprintf(stderr, "-------d_x---------\n");
		for (int i = 0; i < d_x.GetRows(); i++){
			fprintf(stderr, "%.20f\n", d_x(i,0));
		}
*/
		// free the heap!
		//delete [] spike_dx;

		// evaluate d_lambda
		for (int i = 0; i < numContacts; i++){
			d_lambda(i) = lambda_k(i,0)/ff(i,0) * (pow(mfric(3*i,0),2)*xk(3*i,0)*d_x(3*i,0) - xk(3*i+1,0)*d_x(3*i+1,0) -xk(3*i+2,0)*d_x(3*i+2,0) - r_cent(i,0) );
			d_lambda(i + numContacts) = lambda_k(i+numContacts,0)/ff(i+numContacts)*(d_x(3*i) - r_cent(i + numContacts));
		}
/*		fprintf(stderr, "----------d_lambda----------\n");
		for (int i = 0; i < 2*numContacts; i++){
			fprintf(stderr, "%.16f\n", d_lambda(i,0));
		}
*/
		///////////////
		//LINE SEARCH//
		///////////////
		double s_max = 1;
		double tmp;
		for (int i = 0; i < 2*numContacts; i ++){
			if (d_lambda(i,0) < 0){
				tmp = -lambda_k(i,0)/d_lambda(i,0);
//				fprintf(stderr, "i = %d, tmp = %.20f\n", i, tmp);
				if (tmp < s_max){
					s_max = tmp;
				}

			}
		}
		double bla = 0.99;
		double ss = bla * s_max;
//		fprintf(stderr, "s_max = %.20g\n", s_max);

		ff_tmp.Resize(2*numContacts,1);
		lambda_k_tmp.Resize(2*numContacts,1);
		xk_tmp.Resize(3*numContacts,1);;
		r_dual_tmp.Resize(3*numContacts,1);;
		r_cent_tmp.Resize(3*numContacts,1);;

		bool DO = true;
		int count = 0;
//		fprintf(stderr, "----line search----\n");
		while (DO){
			xk_tmp = d_x;
//			fprintf(stderr, "-----d_x----\n");
//			for (int i = 0; i < 3*numContacts; i ++)
//				fprintf(stderr, "%.20g\n", xk_tmp(i,0));
			xk_tmp.MatrScale(ss);
//			fprintf(stderr, "-----ss*d_x----\n");
//			for (int i = 0; i < 3*numContacts; i ++)
//				fprintf(stderr, "%.20g\n", xk_tmp(i,0));
			xk_tmp.MatrAdd(xk,xk_tmp);
//			fprintf(stderr, "-----xk+ss*d_x----\n");
//			for (int i = 0; i < 3*numContacts; i ++)
//				fprintf(stderr, "%.20g\n", xk_tmp(i,0));
			evaluateConstraints(mfric.GetAddress(), numContacts, true);
//			fprintf(stderr, "-----tmp_ff----\n");
//			for (int i = 0; i < 2*numContacts; i ++)
//				fprintf(stderr, "%.20g\n", ff_tmp(i,0));
//			fprintf(stderr, "max_ff = %.20g\n", ff_tmp.Max());
			if (ff_tmp.Max()<0){
				DO = false;
			}
			else{
				count++;
				ss = b1 * ss;
//			fprintf(stderr,"ss[%d] = %.20g\n", count, ss); 
			}
		}

		DO = true;
		double norm_r_t = sqrt(pow(r_dual.NormTwo(),2) + pow(r_cent.NormTwo(),2));
		double norm_r_t_ss;
		count = 0;
		while (DO){
			xk_tmp = d_x;
			xk_tmp.MatrScale(ss);
			xk_tmp.MatrAdd(xk,xk_tmp);

			lambda_k_tmp = d_lambda;
			lambda_k_tmp.MatrScale(ss);
			lambda_k_tmp.MatrAdd(lambda_k, lambda_k_tmp);
			evaluateConstraints(mfric.GetAddress(),numContacts,true);
			sysd.ShurComplementProduct(grad_f, &xk_tmp, 0);
			grad_f.MatrInc(opt_r);
			computeSchurKKT(grad_f.GetAddress(), mfric.GetAddress(), numContacts, barrier_t, true);
			norm_r_t_ss = sqrt(pow(r_dual_tmp.NormTwo(),2) + pow(r_cent_tmp.NormTwo(),2));
			if (norm_r_t_ss < (1 - a1*ss)*norm_r_t)
				DO = false;
			else{
				count ++;
				ss = b1*ss;
//				fprintf(stderr,"ss[%d] = %.20g\n", count, ss); 
		}

		}
		// upadate xk and lambda_k
		d_x.MatrScale(ss);
//		fprintf(stderr, "-------ss*d_x---------\n");
//		for (int i = 0; i < d_x.GetRows(); i++)
//			fprintf(stderr, "%.20f\n", d_x(i,0));
		
//		fprintf(stderr, "----------xk = xk + ss*d_x--------\n");
		xk.MatrInc(d_x);
//		for (int i = 0; i < xk.GetRows(); i++)
//			fprintf(stderr, "%.20f\n", xk(i,0));
		
		d_lambda.MatrScale(ss);
		lambda_k.MatrInc(d_lambda);
//		fprintf(stderr, "-------lambda_k------\n");
//		for (int i = 0; i < lambda_k.GetRows(); i++)
//			fprintf(stderr, "%.20f\n", lambda_k(i,0));

		sysd.ShurComplementProduct(grad_f, &xk, 0);
		grad_f.MatrInc(opt_r);
		evaluateConstraints(mfric.GetAddress(), numContacts, false);
		computeSchurKKT(grad_f.GetAddress(), mfric.GetAddress(), numContacts, barrier_t, false);
//		std::cout << "----r_dual-----" << std::endl;
//		for (int i = 0; i < r_dual.GetRows(); i++)
//			std::cout << r_dual(i,0) << std::endl;
//		std::cout << "-----r_cent-----" << std::endl;
//		for (int i = 0; i < r_cent.GetRows(); i++)
//			std::cout << r_cent(i,0) << std::endl;
fprintf(stderr, "stage[%d], rd = %e, rg = %e, s = %f, t = %f\n", stage+1, r_dual.NormInf(), r_cent.NormInf(), ss, barrier_t);


		if (r_cent.NormInf() < 1e-8 ||stage == (numStages - 1)){
			fprintf(stderr, "solution found after %d stages!\n", stage+1);
//			fprintf(stderr, "stage[%d], rd = %e, rg = %e, s = %f, t = %f\n", stage+1, r_dual.NormInf(), r_cent.NormInf(), ss, barrier_t);
			delete [] BlockDiagonal_val;
			delete [] BlockDiagonal_i;
			delete [] BlockDiagonal_j;
			delete [] spike_rhs;

			/////////////////////////////////////////////
			//set small-magnitude contact force to zero//
			/////////////////////////////////////////////
			
//			for (int i = 0; i < numContacts; i++){
//				if (sqrt(pow(xk(3*i,0),2) + pow(xk(3*i+1,0),2) + pow(xk(3*i+2,0),2)) < 1e-6){
///					xk(3*i,0) = 0;
//					xk(3*i+1,0) = 0;
//					xk(3*i+2, 0) = 0;
//				}
//			}

			sysd.FromVectorToConstraints(xk);
			sysd.FromVectorToVariables(mq);

			for (size_t ic = 0; ic < mconstraints.size(); ic ++){
				if (mconstraints[ic]->IsActive())
					mconstraints[ic]->Increment_q(mconstraints[ic]->Get_l_i());
			}

//			return r_cent.NormInf();
			return 1e-8;

		}
		evaluateConstraints(mfric.GetAddress(), numContacts, false);

	}

}


void ChLcpInteriorPoint::computeSpikeRHS(double *rhs, double *fric, int nc, double t){
	double gamma_n, gamma_u, gamma_v;
	double l1, l2;
	double c1, c2;
	double one = 1;

	for (int i = 0; i < nc; i++){
		gamma_n = xk.GetElementN(3*i);
		gamma_u = xk.GetElementN(3*i + 1);
		gamma_v = xk.GetElementN(3*i + 2);
		l1 = lambda_k.GetElementN(i);
		l2 = lambda_k.GetElementN(i + nc);
		c1 = ff.GetElementN(i);
		c2 = ff.GetElementN(i + nc);
		*(rhs + 3*i) = - pow(*(fric+3*i),2) * gamma_n *(l1 + one/(c1*t)) - (l2 + one/(c2*t)) - r_dual.GetElementN(3*i);
		*(rhs + 3*i+1) = gamma_u * (l1 + one/(c1*t))- r_dual.GetElementN(3*i+1);
		*(rhs + 3*i+2) = gamma_v * (l1 + one/(c1*t))- r_dual.GetElementN(3*i+2);
	}
}


void ChLcpInteriorPoint::evaluateConstraints(double* fric, int nc, bool tmp){
	double *px;
	double *pff;
	double half = 0.5;
	if (tmp == true){
		px = xk_tmp.GetAddress();
		pff = ff_tmp.GetAddress();
	}
	else {
		px = xk.GetAddress();
		pff = ff.GetAddress();
	}
	double gamma_n, gamma_u, gamma_v, mu;
	for (int i = 0; i < nc; i ++){
		gamma_n = *(px + 3*i);
		gamma_u = *(px + 3*i+1);
		gamma_v = *(px + 3*i+2);
		mu = *(fric + 3*i);
		*(pff + i) = half*(pow(gamma_u,2) + pow(gamma_v,2) - pow(mu,2) * pow(gamma_n,2));
		*(pff + i + nc) = -gamma_n;
	}
}


/* computing the rhs for the case of schur */
void ChLcpInteriorPoint::computeSchurRHS(double* grad_f, double* fric, int nc, double t){
	double *prd = r_dual.GetAddress();
	double *prc = r_cent.GetAddress();
	double *plambda = lambda_k.GetAddress();
	double *px = xk.GetAddress();
	double *pff = ff.GetAddress();
	double one = 1;
	for (int i = 0; i < nc; i ++){
		*(prd + 3*i) = (*(grad_f + 3*i)- *(plambda + i)*pow(*(fric+3*i),2)*(*(px + 3*i)) - *(plambda + i +nc));
		*(prd + 3*i + 1) = (*(grad_f + 3*i + 1)+ *(plambda + i)*(*(px + 3*i + 1)));
		*(prd + 3*i + 2) = (*(grad_f + 3*i + 2)+ *(plambda + i)*(*(px + 3*i + 2)));
		*(prc + i) = *(pff + i) + one/(*(plambda + i)*t);
		*(prc + nc +i) = *(pff + nc +i) + one/(*(plambda + nc + i)*t);
	}
}

/* computing the rhs for the case of KKT */
void ChLcpInteriorPoint::computeSchurKKT(double* grad_f, double* fric, int nc, double t, bool tmp){
	double *prd, *prc, *plambda, *px, *pff, one;
	one = 1;
	if (tmp == true){
		prd = r_dual_tmp.GetAddress();
		prc = r_cent_tmp.GetAddress();
		plambda = lambda_k_tmp.GetAddress();
		px = xk_tmp.GetAddress();
		pff = ff_tmp.GetAddress();
	}
	else{
		prd = r_dual.GetAddress();
		prc = r_cent.GetAddress();
		plambda = lambda_k.GetAddress();
		px = xk.GetAddress();
		pff = ff.GetAddress();
	}

	for (int i = 0; i < nc; i ++){
		*(prd + 3*i) = (*(grad_f + 3*i)- *(plambda + i)*pow(*(fric+3*i),2)*(*(px + 3*i)) - *(plambda + i +nc));
		*(prd + 3*i + 1) = (*(grad_f + 3*i + 1)+ *(plambda + i)*(*(px + 3*i + 1)));
		*(prd + 3*i + 2) = (*(grad_f + 3*i + 2)+ *(plambda + i)*(*(px + 3*i + 2)));
		*(prc + i) = -*(pff + i)*(*(plambda+i)) - one/t;
		*(prc + nc +i) = -*(pff + nc +i)*(*(plambda+nc+i)) - one/t;
	}
}

/* formulate block diagonal matrix in row fashion*/
void ChLcpInteriorPoint::computeBlockDiagonal(double* val_diagonal, double *fric, int nc, double t){
	double *plambda = lambda_k.GetAddress();
	double *px = xk.GetAddress();
	double *pff = ff.GetAddress();
	double b11, b12, b13, b22, b23, b33;
	double tmp1, tmp2, tmp3;
	int count = 0;
	double mu;
	double gamma_n, gamma_u, gamma_v;
	for (int i = 0; i < nc; i++){
		mu = *(fric + 3*i);
		tmp1 = *(plambda + i)/(*(pff + i));
		tmp2 = tmp1 * pow(mu,2);
		tmp3 = *(plambda + i + nc)/(*(pff + i + nc));
		gamma_n = *(px + 3*i);
		gamma_u = *(px + 3*i + 1);
		gamma_v = *(px + 3*i + 2);

		b11 = -pow(mu,2)*pow(gamma_n,2)*tmp2 - tmp3 - pow(mu,2)*(*(plambda + i));
		b12 = gamma_n * gamma_u * tmp2;
		b13 = gamma_n * gamma_v * tmp2;
		b22 = -pow(gamma_u,2) * tmp1 + (*(plambda+i));
		b23 = -gamma_u * gamma_v * tmp1;
		b33 = -pow(gamma_v,2) * tmp1 + (*(plambda+i));

		/* assign value to sparse block! */
		*(val_diagonal + 9*i) = b11;
		*(val_diagonal + 9*i + 1) = b12;
		*(val_diagonal + 9*i + 2) = b13;
		*(val_diagonal + 9*i + 3) = b12;
		*(val_diagonal + 9*i + 4) = b22;
		*(val_diagonal + 9*i + 5) = b23;
		*(val_diagonal + 9*i + 6) = b13;
		*(val_diagonal + 9*i + 7) = b23;
		*(val_diagonal + 9*i + 8) = b33;
//		std::cout << "Block " << i << std::endl;
//		std::cout << b11 << " " << b12 << " " << b13 << std::endl;
//		std::cout << b12 << " " << b22 << " " << b23 << std::endl;
//		std::cout << b13 << " " << b23 << " " << b33 << std::endl;

	}
}



} // END_OF_NAMESPACE____


