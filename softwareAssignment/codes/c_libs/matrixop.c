#include <stdio.h>
#include "matrixop.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>


void mat_mul(const double *A, const double *B, double *C, int m, int p, int n) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double sum = 0;
            for (int z = 0; z < p; z++){
                sum += A[indx(i, z, p)] * B[indx(z, j, n)];
            }
            C[indx(i, j, n)] = sum;
        }
}

void(const double *A, double *B, int m, int n) {
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++)
            B[indx(j, i, m)] = A[indx(i, j, n)];
}


void ata(const double *A, double *AtA, int m, int n) {
    for (int i = 0; i < n * n; i++) AtA[i] = 0.0;
    for (int i = 0; i < m; i++)
        for (int j = 0; j < n; j++) {
            double aij = A[indx(i, j, n)];
            for (int k = 0; k < n; k++)
                AtA[indx(j, k, n)] += aij * A[indx(i, k, n)];
        }
}


double frobenius_norm(const double *A, const double *B, int m, int n) {
    double sum = 0.0;
    for (int i = 0; i < m * n; i++) {
        double diff = A[i] - B[i];
        sum += diff * diff;
    }
    return sqrt(sum);
}
double frobenius_norm_single(const double *A, int m, int n) {
    double sum = 0.0;
    for (int i = 0; i < m * n; i++)
        sum += A[i] * A[i];
    return sqrt(sum);
}
