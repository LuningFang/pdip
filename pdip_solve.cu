#include <iostream>
#include <stdlib.h>

#include <cusp/array1d.h>
#include <cusp/csr_matrix.h>
#include <cusp/coo_matrix.h>
#include <cusp/blas.h>
#include <cusp/multiply.h>
#include <cusp/elementwise.h>
#include <cusp/transpose.h>

#include <thrust/copy.h>

#include <spike/solver.h>
#include <spike/spmv.h>
#include <spike/exception.h>

#include "pdip_solve.h"

typedef typename cusp::csr_matrix<int, double, cusp::device_memory> Matrix_csr;
typedef typename cusp::coo_matrix<int, double, cusp::device_memory> Matrix_coo;
typedef typename cusp::array1d<double, cusp::device_memory>         Vector;
typedef typename spike::SpmvCusp<Matrix_csr>                        SpmvFunctor;
typedef typename spike::Solver<Vector, double>                      SpikeSolver;

bool
solveSPIKE(const int nb, const int nc, const int* Cq_i, const int* Cq_j, const int Cq_nnz, const double *Cq_val, const int* Minv_i, const int* Minv_j, const double * Minv_val, const int* BD_i, const int* BD_j, const double * BD_val, double* x, const double* rhs, bool print)
{
	spike::GPUTimer gputimer;
	
	gputimer.Start();
	
	// assemble coo matrix for Cq, Minv and BlockDiagonal
	cusp::coo_matrix<int, double, cusp::host_memory> cooCq(6*nb, 3*nc, Cq_nnz);
	thrust::copy(Cq_i, Cq_i + Cq_nnz, cooCq.row_indices.begin());
	thrust::copy(Cq_j, Cq_j + Cq_nnz, cooCq.column_indices.begin());
	thrust::copy(Cq_val, Cq_val + Cq_nnz, cooCq.values.begin());

	cusp::coo_matrix<int, double, cusp::host_memory> cooMinv(6*nb, 6*nb, 6*nb);
	thrust::copy(Minv_i, Minv_i + 6*nb, cooMinv.row_indices.begin());
	thrust::copy(Minv_j, Minv_j + 6*nb, cooMinv.column_indices.begin());
	thrust::copy(Minv_val, Minv_val + 6*nb, cooMinv.values.begin());
//		fprintf(stderr, "---------Minv----------\n");
//		for (int i = 0; i < cooMinv.num_entries;i++){
//		fprintf(stderr, "%d		%d		%.20f\n", cooMinv.row_indices[i], cooMinv.column_indices[i], cooMinv.values[i]);
//		}

	cusp::coo_matrix<int, double, cusp::host_memory> cooBD(3*nc, 3*nc, 9*nc);
	thrust::copy(BD_i, BD_i + 9*nc, cooBD.row_indices.begin());
	thrust::copy(BD_j, BD_j + 9*nc, cooBD.column_indices.begin());
	thrust::copy(BD_val, BD_val + 9*nc, cooBD.values.begin());
	
	// move sparse matrices onto device in coo format
	Matrix_coo cooCq_d = cooCq;
	Matrix_coo cooMinv_d = cooMinv;
	Matrix_coo cooBD_d = cooBD;
	Matrix_coo cooA_d;
	Matrix_csr csrA;

	// assemble right hand side
	Vector cusp_rhs(rhs, rhs + 3*nc);
	Vector cusp_x(3*nc, 0);


	// assemble coefficient matrix A = Cq^T*Minv*Cq + BD
	{
		Matrix_coo CqT;
		Matrix_coo C;
		Matrix_coo N;
		cusp::transpose(cooCq_d, CqT);
		cooCq_d.sort_by_row_and_column();
		cooMinv_d.sort_by_row_and_column();
		cusp::multiply(cooMinv_d, cooCq_d, C);
		C.sort_by_row_and_column();
		CqT.sort_by_row_and_column();
		
//		cusp::coo_matrix<int, double, cusp::host_memory> CqTh = CqT;
//		fprintf(stderr, "---------Cq^T----------\n");
//		for (int i = 0; i < CqTh.num_entries;i++){
//		fprintf(stderr, "%d		%d		%.20f\n", CqTh.row_indices[i], CqTh.column_indices[i], CqTh.values[i]);
//		}
//		cusp::coo_matrix<int, double, cusp::host_memory> Ch = C;
//		fprintf(stderr, "---------Minv*Cq----------\n");
//		for (int i = 0; i < Ch.num_entries;i++){
//		fprintf(stderr, "%d		%d		%.20f\n", Ch.row_indices[i], Ch.column_indices[i], Ch.values[i]);
//		}
		cusp::multiply(CqT, C, N);
		cusp::add(N, cooBD_d, cooA_d);
		csrA = cooA_d;
		
//		cusp::coo_matrix<int, double, cusp::host_memory> Nh = N;
//		fprintf(stderr, "--------N(should be the same)!----------\n");
//		for (int i = 0; i < Nh.num_entries;i++){
//		fprintf(stderr, "%d		%d		%.20f\n", Nh.row_indices[i], Nh.column_indices[i], Nh.values[i]);
//		}
//		cusp::coo_matrix<int, double, cusp::host_memory> cooBDh = cooBD_d;
//		fprintf(stderr, "--------Block Diagonal!----------\n");
//		for (int i = 0; i < cooBDh.num_entries;i++){
//		fprintf(stderr, "%d		%d		%.20f\n", cooBDh.row_indices[i], cooBDh.column_indices[i], cooBDh.values[i]);
//		}
		cusp::coo_matrix<int, double, cusp::host_memory> cooA_h = cooA_d;
		if ( print == true){
			fprintf(stderr, "--------coefficient A!----------\n");
			for (int i = 0; i < cooA_h.num_entries;i++){
			fprintf(stderr, "%d		%d		%.20f\n", cooA_h.row_indices[i], cooA_h.column_indices[i], cooA_h.values[i]);
		}
	}


	}
	
	// setup spike solver
	spike::Options opts;
	opts.relTol = 1e-10;
	opts.safeFactorization = true;
	int numPart = 1;	

	SpmvFunctor mySpmv(csrA);
	SpikeSolver mySolver(numPart, opts);
	
	bool success = false;

	try {
		mySolver.setup(csrA);
		success = mySolver.solve(mySpmv, cusp_rhs, cusp_x);
	} catch (const std::bad_alloc &e) {
		std::cerr << "Exception (bad_alloc): " << e.what() << std::endl;
		return false;
	} catch (const spike::system_error& e) {
		std::cerr << "Exception (system_error): " << e.what() << " Error code: " << e.reason() << std::endl;
		return false;
	}

	if (success){
		cusp::array1d<double, cusp::host_memory> cusp_x_h = cusp_x;
		thrust::copy(cusp_x_h.begin(), cusp_x_h.end(), x);
//		cusp::array1d<double, cusp::host_memory> cusp_rhs_h = cusp_rhs;
//		fprintf(stderr, "--------cusp_rhs---------size of %d\n", cusp_rhs_h.size());
//		for (int i = 0; i < cusp_rhs_h.size(); i++){
//			fprintf(stderr, "%.20f\n", cusp_rhs_h[i]);
//		}	
	
//		fprintf(stderr, "--------cusp_x---------\n");
//		for (int i = 0; i < cusp_x_h.size(); i++){
//			fprintf(stderr, "%.20f\n", cusp_x_h[i]);
//		}


	}
	else {
		std::cerr << "The SPIKE solver fails to solve the equation!" << std::endl;
		return false;
	}
	gputimer.Stop();
	fprintf(stderr, "Solve time: %g ", gputimer.getElapsed());

	return true;
}
