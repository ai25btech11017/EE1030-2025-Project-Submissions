#ifndef matrixop_h
#define matrixop_h

void mat_mul(const double *A, const double *B, double *C, int m, int p, int n);

void trap(const double *A, double *B, int m, int n);

double frobenius_norm(const double *A, const double *B, int m, int n);

double frobenius_norm_single(const double *A, int m, int n);

void ata(const double *A, double *AtA, int m, int n);


#endif
