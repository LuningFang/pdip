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


#include "ChLcpInteriorPointActiveSets.h"
#include "ChLcpConstraintTwoFrictionT.h"
#include "pdip_solve.h"

namespace chrono
{

double ChLcpInteriorPointActiveSets::Solve(
		ChLcpSystemDescriptor& sysd		///< system description with constraints and variables
)
{
	std::cout << "-------using interior point active sets solver!!------" << std::endl;
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

		ChStreamOutAsciiFile file_V0( "dump_V_old.dat" ) ;
		mv0.StreamOUTdenseMatlabFormat(file_V0) ;
		ChStreamOutAsciiFile file_M ( "dump_M.dat" ) ;
		mM.StreamOUTsparseMatlabFormat ( file_M ) ;

		ChStreamOutAsciiFile file_Cq ( "dump_Cq.dat" ) ;
		mCq.StreamOUTsparseMatlabFormat ( file_Cq ) ;

		ChStreamOutAsciiFile file_E ( "dump_E.dat" ) ;
		mE.StreamOUTsparseMatlabFormat ( file_E ) ;
		ChStreamOutAsciiFile file_fric ( "dump_fric.dat" ) ;
		mfric.StreamOUTdenseMatlabFormat ( file_fric ) ;
		ChStreamOutAsciiFile file_f ( "dump_f.dat" ) ;
		mf.StreamOUTdenseMatlabFormat ( file_f ) ;


		ChStreamOutAsciiFile file_b ( "dump_b.dat" ) ;
		mb.StreamOUTdenseMatlabFormat ( file_b ) ;

	
		printf("Successfully writing files!\n");



/*	file_f.GetFstream().close();
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
	int nc_original = numContacts;

	/* ALWYAS DO THIS IN THE LCP SOLVER!!!*/
	for (unsigned int ic = 0; ic < nConstraints; ic++)
		mconstraints[ic]->Update_auxiliary();





	//Get sparse info of contact Jacobian Cq
	std::vector<int> index_i_Cq;
	std::vector<int> index_j_Cq;
	std::vector<double> val_Cq;
	double val;
	for (int ii = 0; ii < mCq.GetRows(); ii++){
		for (int jj = 0; jj < mCq.GetColumns(); jj++){
			val = mCq.GetElement(ii,jj);
			if (val){
				index_i_Cq.push_back(jj);
				index_j_Cq.push_back(ii);
				val_Cq.push_back(val);
			}
		}
	}

	///////////////////////////////
	//Parameters for reduced PDIP//
	///////////////////////////////
	bool modified = false;
	double res = 5e-4;
	std::vector<int> throwaway;


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


	/////////////////////////////////////////
	//Formulate rhs of optimization problem//
	/////////////////////////////////////////
	ChMatrixDynamic <double> opt_r_tmp(nConstraints,1);
	for (unsigned int iv = 0; iv < nVariables; iv ++)  // M^[-1] * k
		if (mvariables[iv]->IsActive()){
			mvariables[iv]->Compute_invMb_v(mvariables[iv]->Get_qb(), mvariables[iv]->Get_fb());
			ChMatrix<double> k = mvariables[iv]->Get_fb();
			ChMatrix<double> Mk = mvariables[iv]->Get_qb();
		}

	int s_i = 0; //
	opt_r.Resize(nConstraints,1);
	for (unsigned int ic = 0; ic < nConstraints; ic ++)
		if (mconstraints[ic]->IsActive()){
			opt_r(s_i,0) = mconstraints[ic]->Compute_Cq_q();
			++s_i;	
		}
	sysd.BuildBiVector(opt_r_tmp);
	opt_r.MatrInc(opt_r_tmp);

	///////////////////
	//velocity update//
	///////////////////
	ChMatrixDynamic<double> mq;
	sysd.FromVariablesToVector(mq, true);


//	fprintf(stderr, "----mq----\n");
//	printChMatrix(mq);


	////////////////////////////
	//assign solver parameters//
	////////////////////////////
	double barrier_t = 1;
	double eta_hat;
	int numStages = 725;
	int mu1 = 10;
	double b1 = 0.8;
	double a1 = 0.01;
	
	double res1 = 1e-10;
	double res2 = 1e-8;

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


	int throwaway_size = 0;
	int erase_start;
	int erase_end;


	bool DO = true;
	ChMatrixDynamic<double> xk_padded(3*numContacts,1); //xk padded with zeros
	bool contactTable[numContacts];
	for (int i = 0; i < numContacts; i++)
		contactTable[i] = true; // set true for non-zero items


	/////////////////////////////
	////GO THROUGH EACH STAGE////
	/////////////////////////////
	for (int stage = 0; stage < numStages; stage++){

		///////////////////
		// Modified PDIP //
		///////////////////





		throwaway_size = throwaway.size();
		for (int i = 0; i < throwaway.size(); i++)
			fprintf(stderr, "%d ", throwaway[i]);

		fprintf(stderr, "thrown away\n");
		if (throwaway_size != 0)
			modified = true;
		// shrink problem size!
		if (modified == true){
			int throw_nc;
			std::vector <int> eraseCqIndex;  // vector stores the index


			eraseCqIndex.resize(2*throwaway_size);

			// get the index Cq_i, Cq_j, Cq_val to be thrown away
			int i = 0;
			int j = 0;

/*			fprintf(stderr, "-------index_j_Cq before shrinked-----\n");
			for (int i = 0; i < index_j_Cq.size(); i ++)
				fprintf(stderr, "loc[%d]: %d\n", i , index_j_Cq[i]);

			fprintf(stderr, "-------Cq before shrinked-----\n");
			for (int i = 0; i < index_j_Cq.size(); i ++)
				fprintf(stderr, "%d, %d, %.20f\n", index_i_Cq[i], index_j_Cq[i], val_Cq[i] );
*/
			////////////////////////////////
			//Need to rewrite this part...//
			////////////////////////////////
			DO = true;
			while (DO){
				
				throw_nc = throwaway[j];
				
				while (index_j_Cq[i] < 3*throw_nc && i < index_j_Cq.size()){
					index_j_Cq[i] = index_j_Cq[i] - 3*j;
					i ++;
					
				}
				
				eraseCqIndex[2*j] = i;
				
				while (index_j_Cq[i] <= 3*throw_nc+2 && i < index_j_Cq.size()){
						index_j_Cq[i] = index_j_Cq[i] - 3*j;	
						i ++;
				}
				
				eraseCqIndex[2*j+1] = i;
				
				j = j + 1;
				
				if ( i >= index_j_Cq.size())
					i--;
				while ( j >= throwaway_size && i < index_j_Cq.size()){
					index_j_Cq[i] = index_j_Cq[i] - 3*j;
					i ++;
					if (i >= index_j_Cq.size()){
						DO = false;
						break;
					}
						
				}


		}
		

/*			fprintf(stderr, "-------index_j_Cq after shrinked-----\n");
			for (int i = 0; i < index_j_Cq.size(); i ++)
				fprintf(stderr, "loc[%d]: %d\n", i , index_j_Cq[i]);
*/

/*		fprintf(stderr, "----eraseCqIndex-----\n");
		for (int i = 0; i < 2*throwaway_size; i ++)
			fprintf(stderr, "%d \n", eraseCqIndex[i]);
*/
			int numDeleted = 0;
			for (i = 0; i < throwaway_size; i ++){
				erase_start = eraseCqIndex[2*i];
				erase_end = eraseCqIndex[2*i+1];

				index_i_Cq.erase(index_i_Cq.begin()+erase_start - numDeleted,index_i_Cq.begin()+erase_end - numDeleted);

				index_j_Cq.erase(index_j_Cq.begin()+erase_start - numDeleted,index_j_Cq.begin()+erase_end - numDeleted);
/*				fprintf(stderr, "erase #%d from location[%d] to [%d]\n", i, erase_start - numDeleted, erase_end - numDeleted);
				for (int i = 0; i < index_j_Cq.size(); i++)
					fprintf(stderr, "loc[%d]: %d\n", i , index_j_Cq[i]);
*/				

				val_Cq.erase(val_Cq.begin()+erase_start - numDeleted,val_Cq.begin()+erase_end - numDeleted);
				
				numDeleted += (erase_end - erase_start);
			}

/*			fprintf(stderr, "-------Cq_shrinked-----\n");
			for (int i = 0; i < index_j_Cq.size(); i ++)
				fprintf(stderr, "%d, %d, %.20f\n", index_i_Cq[i], index_j_Cq[i], val_Cq[i] );
*/
			Cq_i = &index_i_Cq[0];
			Cq_j = &index_j_Cq[0];
			Cq_nnz = val_Cq.size();
			Cq_val = &val_Cq[0];
			

			/////////////////////////////////////
			//Shrink x_k, opt_r, m_fric, lambda//
			/////////////////////////////////////
			{
				ChMatrixDynamic<double> xk_shrinked(3*(numContacts - throwaway_size),1);
				ChMatrixDynamic<double> opt_r_shrinked(3*(numContacts - throwaway_size),1);
				ChMatrixDynamic<double> m_fric_shrinked(3*(numContacts - throwaway_size),1);
				ChMatrixDynamic<double> lambda_shrinked(2*(numContacts - throwaway_size),1);
				ChMatrixDynamic<double> ff_shrinked(2*(numContacts - throwaway_size),1);

				int index = 0; // index going through throwaway vector
				int j=0; // index going through the shrinked variables
				for (int i = 0; i < numContacts; i ++){
					if ( throwaway[index] == i){
						index = index + 1;
						if (index >= throwaway.size())
							index = throwaway.size()-1; //check index goes out of bounds.
					}
					else{
						xk_shrinked.SetElementN(3*j  , xk.GetElementN(3*i  ));
						xk_shrinked.SetElementN(3*j+1, xk.GetElementN(3*i+1));
						xk_shrinked.SetElementN(3*j+2, xk.GetElementN(3*i+2));

						opt_r_shrinked.SetElementN(3*j  , opt_r.GetElementN(3*i  ));
						opt_r_shrinked.SetElementN(3*j+1, opt_r.GetElementN(3*i+1));
						opt_r_shrinked.SetElementN(3*j+2, opt_r.GetElementN(3*i+2));

						m_fric_shrinked.SetElementN(3*j  , mfric.GetElementN(3*i  ));
						m_fric_shrinked.SetElementN(3*j+1, mfric.GetElementN(3*i+1));
						m_fric_shrinked.SetElementN(3*j+2, mfric.GetElementN(3*i+2));

						lambda_shrinked.SetElementN(j  , lambda_k.GetElementN(i  ));
						lambda_shrinked.SetElementN(j + numContacts - throwaway_size, lambda_k.GetElementN(i + numContacts));

						ff_shrinked.SetElementN(j  , ff.GetElementN(i  ));
						ff_shrinked.SetElementN(j + numContacts - throwaway_size, ff.GetElementN(i + numContacts));


						j = j + 1;

					}
				}


				xk.Resize(3*(numContacts - throwaway_size),1);
				xk = xk_shrinked;

				opt_r.Resize(3*(numContacts - throwaway_size),1);
				opt_r = opt_r_shrinked;

				mfric.Resize(3*(numContacts - throwaway_size),1);
				mfric = m_fric_shrinked;

				lambda_k.Resize(2*(numContacts - throwaway_size),1);
				lambda_k = lambda_shrinked;
	
				ff.Resize(2*(numContacts - throwaway_size),1);
				ff = ff_shrinked;
	


			}
/*			
				fprintf(stderr, "-----opt_r_shrinked----\n");
				for (int i = 0; i < opt_r.GetRows(); i ++)
					fprintf(stderr, "%.20f\n", opt_r.GetElementN(i));
			
				fprintf(stderr, "-----lambda_k_shrinked----\n");
				for (int i = 0; i < lambda_k.GetRows(); i ++)
					fprintf(stderr, "%.20f\n", lambda_k.GetElementN(i));
*/

			numContacts = numContacts - throwaway_size;
			//////////////////////////////
			//Resize remaining variables//
			//////////////////////////////
			grad_f.Resize(3*numContacts,1);
			ff.Resize(2*numContacts,1);
			r_dual.Resize(3*numContacts,1);
			r_cent.Resize(2*numContacts,1);
			d_x.Resize(3*numContacts,1);
			d_lambda.Resize(2*numContacts,1);
			Schur_rhs.Resize(3*numContacts,1);
		}
		////////////////////////////////
		//generate xk_padded (always!)//
		////////////////////////////////
		if (xk.GetRows() != 3*nc_original){
			int j = 0; // index goes through shrinked xk
			for (int i = 0; i < nc_original; i ++){
				if (contactTable[i] == true){
					xk_padded.SetElementN(3*i  , xk.GetElementN(3*j  ));
					xk_padded.SetElementN(3*i+1, xk.GetElementN(3*j+1));
					xk_padded.SetElementN(3*i+2, xk.GetElementN(3*j+2));
					j ++;
				}
				else{
					xk_padded.SetElementN(3*i, 0);
					xk_padded.SetElementN(3*i+1,0);
					xk_padded.SetElementN(3*i+2,0);
				}
			}
		}
		else
			xk_padded = xk;
/*				fprintf(stderr, "-----ff_shrinked----\n");
				for (int i = 0; i < ff.GetRows(); i ++)
					fprintf(stderr, "%.20f\n", ff.GetElementN(i));
*/

		eta_hat = - lambda_k.MatrDot(&lambda_k, &ff);
		barrier_t = mu1 * (2*numContacts)/eta_hat;
		// assemble grad_f = N*x + r
		//
		if (grad_f.GetRows() != 3*nc_original){
			ChMatrixDynamic<double> grad_f_padded(3*nc_original,1);
			sysd.ShurComplementProduct(grad_f_padded, &xk_padded, 0);

			// map grad_f_padded to grad_f_shrinked
			int j = 0; // index goes through shrinked grad_f
			for (int i = 0; i < nc_original; i ++){
				if (contactTable[i] == true){
					grad_f.SetElementN(3*j, grad_f_padded.GetElementN(3*i));
					grad_f.SetElementN(3*j+1, grad_f_padded.GetElementN(3*i+1));
					grad_f.SetElementN(3*j+2, grad_f_padded.GetElementN(3*i+2));
					j ++;
				}
			}



	}
		else {
			sysd.ShurComplementProduct(grad_f, &xk_padded, 0);
}
				grad_f.MatrInc(opt_r);

/*				fprintf(stderr, "-----grad_f_shrinked----\n");
				for (int i = 0; i < grad_f.GetRows(); i ++)
					fprintf(stderr, "%.20f\n", grad_f.GetElementN(i));
*/
			// compute r_d and r_c for schur implementation
		computeSchurRHS(grad_f.GetAddress(), mfric.GetAddress(), numContacts, barrier_t);


		// assemble block diagonal matrix
		computeBlockDiagonal(BlockDiagonal_val, mfric.GetAddress(), numContacts, barrier_t);

		// assemble rhs vector for spike solver
		computeSpikeRHS(spike_rhs, mfric.GetAddress(), numContacts, barrier_t);

		double *spike_dx = new double [3*numContacts];

		//call ang's solver here....
		bool solveSuc = solveSPIKE(nBodies, numContacts, Cq_i, Cq_j, Cq_nnz, Cq_val, Minv_i, Minv_j, Minv_val, BlockDiagonal_i, BlockDiagonal_j, BlockDiagonal_val, spike_dx, spike_rhs, false);

		if (solveSuc == false)
			std::cerr << "Solve Failed!" << std::endl;
		for (int i = 0; i < numContacts; i++){
			d_x(3*i,0) = *(spike_dx + 3*i);
			d_x(3*i+1,0) = *(spike_dx + 3*i + 1);
			d_x(3*i+2,0) = *(spike_dx + 3*i + 2);
		}

/*				fprintf(stderr, "-----spike_rhs_shrinked----\n");
				for (int i = 0; i < 3*numContacts; i ++)
					fprintf(stderr, "%.20f\n", *(spike_rhs+i));
			
				fprintf(stderr, "-----spike_dx_shrinked----\n");
				for (int i = 0; i < d_x.GetRows(); i ++)
					fprintf(stderr, "%.20f\n", d_x.GetElementN(i));
*/
		// free the heap!
		delete [] spike_dx;
		
		// evaluate d_lambda
		for (int i = 0; i < numContacts; i++){
			d_lambda(i) = lambda_k(i,0)/ff(i,0) * (pow(mfric(3*i,0),2)*xk(3*i,0)*d_x(3*i,0) - xk(3*i+1,0)*d_x(3*i+1,0) -xk(3*i+2,0)*d_x(3*i+2,0) - r_cent(i,0) );
			d_lambda(i + numContacts) = lambda_k(i+numContacts,0)/ff(i+numContacts)*(d_x(3*i) - r_cent(i + numContacts));
		}

/*		fprintf(stderr, "------d_x------\n");
		for (int i = 0; i < 3*numContacts; i++)
			fprintf(stderr, "%.20f\n", d_x(i));
*/


		///////////////
		//LINE SEARCH//
		///////////////
		double s_max = 1;
		double tmp;
		for (int i = 0; i < 2*numContacts; i ++){
			if (d_lambda(i,0) < 0){
				tmp = -lambda_k(i,0)/d_lambda(i,0);
				if (tmp < s_max){
					s_max = tmp;
				}

			}
		}
		double bla = 0.99;
		double ss = bla * s_max;

		ff_tmp.Resize(2*numContacts,1);
		lambda_k_tmp.Resize(2*numContacts,1);
		xk_tmp.Resize(3*numContacts,1);;
		x_old_tmp.Resize(3*numContacts,1);;
		r_dual_tmp.Resize(3*numContacts,1);;
		r_cent_tmp.Resize(3*numContacts,1);;

		DO = true;
		int count = 0;
		while (DO){
			xk_tmp = d_x;
			xk_tmp.MatrScale(ss);
			xk_tmp.MatrAdd(xk,xk_tmp);
			evaluateConstraints(mfric.GetAddress(), numContacts, true);
			if (ff_tmp.Max()<0){
				DO = false;
			}
			else{
				count++;
				ss = b1 * ss;
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


			{

				ChMatrixDynamic<double> xk_tmp_padded(3*nc_original,1);
				int j = 0; // index goes through shrinked xk_tmp
				for (int i = 0; i < nc_original; i ++){
					if (contactTable[i] == true){
						xk_tmp_padded.SetElementN(3*i  , xk_tmp.GetElementN(3*j  ));
						xk_tmp_padded.SetElementN(3*i+1, xk_tmp.GetElementN(3*j+1));
						xk_tmp_padded.SetElementN(3*i+2, xk_tmp.GetElementN(3*j+2));
						j ++;
					}
					else{
						xk_tmp_padded.SetElementN(3*i, 0);
						xk_tmp_padded.SetElementN(3*i+1,0);
						xk_tmp_padded.SetElementN(3*i+2,0);
					}
				}

/*		fprintf(stderr, "------xk_tmp - xk_tmp_padded------\n");
		for (int i = 0; i < 3*numContacts; i++)
			fprintf(stderr, "%.20f\n", xk_tmp(i)- xk_tmp_padded(i));
*/

				ChMatrixDynamic<double> grad_f_padded(3*nc_original,1);
				sysd.ShurComplementProduct(grad_f_padded, &xk_padded, 0);

				// map grad_f_padded to grad_f_shrinked
				j = 0; // index goes through shrinked grad_f
				for (int i = 0; i < nc_original; i ++){
					if (contactTable[i] == true){
						grad_f.SetElementN(3*j, grad_f_padded.GetElementN(3*i));
						grad_f.SetElementN(3*j+1, grad_f_padded.GetElementN(3*i+1));
						grad_f.SetElementN(3*j+2, grad_f_padded.GetElementN(3*i+2));
						j ++;
					}
				}
			}
			grad_f.MatrInc(opt_r);

/*			fprintf(stderr, "------grad_f = N*xk +r ----with ss = %f---\n", ss);
			for (int i = 0; i < 3*numContacts; i++)
				fprintf(stderr, "%.20f\n", grad_f(i));
*/

			computeSchurKKT(grad_f.GetAddress(), mfric.GetAddress(), numContacts, barrier_t, true);
			norm_r_t_ss = sqrt(pow(r_dual_tmp.NormTwo(),2) + pow(r_cent_tmp.NormTwo(),2));
			if (norm_r_t_ss < (1 - a1*ss)*norm_r_t)
				DO = false;
			else{
				count ++;
				ss = b1*ss;
			}

		}
		// upadate xk and lambda_k
		d_x.MatrScale(ss);

		xk.MatrInc(d_x);

		d_lambda.MatrScale(ss);
		lambda_k.MatrInc(d_lambda);

		{



			ChMatrixDynamic<double> grad_f_padded(3*nc_original,1);
			sysd.ShurComplementProduct(grad_f_padded, &xk_padded, 0);

			// map grad_f_padded to grad_f_shrinked
			int j = 0; // index goes through shrinked grad_f
			for (int i = 0; i < nc_original; i ++){
				if (contactTable[i] == true){
					grad_f.SetElementN(3*j, grad_f_padded.GetElementN(3*i));
					grad_f.SetElementN(3*j+1, grad_f_padded.GetElementN(3*i+1));
					grad_f.SetElementN(3*j+2, grad_f_padded.GetElementN(3*i+2));
					j ++;
				}
			}
		}
		grad_f.MatrInc(opt_r);

		evaluateConstraints(mfric.GetAddress(), numContacts, false);
		computeSchurKKT(grad_f.GetAddress(), mfric.GetAddress(), numContacts, barrier_t, false);
		fprintf(stderr, "stage[%d], rd = %e, rg = %e, s = %f, t = %f\n", stage+1, r_dual.NormInf(), r_cent.NormInf(), ss, barrier_t);


/*		fprintf(stderr, "-------xk-------\n");
		for (int i = 0; i < xk.GetRows(); i++){
			fprintf(stderr, "%.20f\n", xk.GetElementN(i));
		}
		fprintf(stderr, "-------contactTable-------\n");	
		for (int i = 0; i < nc_original; i++){
			fprintf(stderr, "%d\n", contactTable[i]);
		}
*/

		throwaway = findModifiedIndex(numContacts, res);
		modified = false;
		if (throwaway.size() != 0){
			modified = true;
/*			for (int i = 0; i < throwaway.size(); i++)
				fprintf(stderr, "%d\n",throwaway[i]);
*/		}


		// update contact table
		if (modified == true){
			int i = 0;
			int nc_zero = throwaway.at(i);
			int count = 0;
			std::vector<int> tmp;

			for (int countTable = 0; countTable < 3*nc_original; countTable ++){
				if (count == nc_zero && contactTable[countTable] == true){

					tmp.push_back(countTable);
					i++;

					if ( i >= throwaway.size())
						break;
					else
						nc_zero = throwaway.at(i);
				}

				if (contactTable[countTable] == true)
					count ++;
			}

			for (i = 0; i < tmp.size(); i ++){
				contactTable[tmp.at(i)] = false;
			}

/*		    fprintf(stderr, "-----countTable----\n");
			for (int i = 0; i < nc_original; i ++)
				fprintf(stderr, "contactTable[%d] = %d\n", i, contactTable[i]);
*/		}



		if (r_cent.NormInf() < res1 || stage+1 >= numStages){
						//			fprintf(stderr, "stage[%d], rd = %e, rg = %e, s = %f, t = %f\n", stage+1, r_dual.NormInf(), r_cent.NormInf(), ss, barrier_t);
			delete [] BlockDiagonal_val;
			delete [] BlockDiagonal_i;
			delete [] BlockDiagonal_j;
			delete [] spike_rhs;

			int j = 0; // index goes through shrinked xk
			for (int i = 0; i < nc_original; i ++){
				if (contactTable[i] == true){
					xk_padded.SetElementN(3*i  , xk.GetElementN(3*j  ));
					xk_padded.SetElementN(3*i+1, xk.GetElementN(3*j+1));
					xk_padded.SetElementN(3*i+2, xk.GetElementN(3*j+2));
					j ++;
				}
				else{
					xk_padded.SetElementN(3*i, 0);
					xk_padded.SetElementN(3*i+1,0);
					xk_padded.SetElementN(3*i+2,0);
				}
			}


/*			fprintf(stderr, "----------xk--------");
			for (int i = 0; i < xk_padded.GetRows(); i++){
				fprintf(stderr, "%.20f\n", xk_padded.ElementN(i));
			}
*/

			sysd.FromVectorToConstraints(xk_padded);
			sysd.FromVectorToVariables(mq);

			for (size_t ic = 0; ic < nConstraints; ic ++){
				if (mconstraints[ic]->IsActive())
					mconstraints[ic]->Increment_q(mconstraints[ic]->Get_l_i());
			}
			fprintf(stderr, "solution found after %d stages!\n", stage+1);

			return r_cent.NormInf();


		}

	}

}

std::vector<int> ChLcpInteriorPointActiveSets::findModifiedIndex(int nc, double res){
	std::vector<int> index;
	for (int i = 0; i < nc; i++){
		if (ff.GetElementN(nc + i) < res && ff.GetElementN(nc + i) > -res){
			index.push_back(i);
		}

	}
	return index;
}

void ChLcpInteriorPointActiveSets::computeSpikeRHS(double *rhs, double *fric, int nc, double t){
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


void ChLcpInteriorPointActiveSets::evaluateConstraints(double* fric, int nc, bool tmp){
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
void ChLcpInteriorPointActiveSets::computeSchurRHS(double* grad_f, double* fric, int nc, double t){
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
void ChLcpInteriorPointActiveSets::computeSchurKKT(double* grad_f, double* fric, int nc, double t, bool tmp){
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
void ChLcpInteriorPointActiveSets::computeBlockDiagonal(double* val_diagonal, double *fric, int nc, double t){
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

	}
}



} // END_OF_NAMESPACE____


