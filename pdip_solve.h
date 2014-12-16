#ifndef CH_PDIP_SOLVE_H
#define CH_PDIP_SOLVE_H

bool solveSPIKE(const int nb, const int nc, const int* Cq_i, const int* Cq_j, const int Cq_nnz, const double *Cq_val, const int* Minv_i, const int* Minv_j, const double * Minv_val, const int* BD_i, const int* BD_j, const double * BD_val, double* x, const double* rhs, bool print);

#endif
