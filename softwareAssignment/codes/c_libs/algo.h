#ifndef ALGO_H
#define ALGO_H


void jacobi_algo(double *A, double *eigenvalues, double *V, int n);

void sort_desc(double *sigma, double *V, int n);

void compute_svd(const double *A, int m, int n,
                 double **U_out, double **S_out, double **V_out, int *r_out);

void reconstruct_image(const double *U, const double *S, const double *V,
                       int m, int n, int r, int k, double *Ak);

#endif
