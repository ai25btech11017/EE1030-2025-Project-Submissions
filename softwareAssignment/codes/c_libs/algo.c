#include <stdio.h>
#include "algo.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include"matrixop.h"
#include"pgmfunctions.h"
void jacobi_algo(double *A, double *eigenvalues, double *V, int n) {
    for (int i = 0; i < n * n; i++) V[i] = 0.0;
    for (int i = 0; i < n; i++) V[indx(i, i, n)] = 1.0;

    for (int iter = 0; iter < MAX_ITER; iter++) {
        double max_val = 0.0;
        int p = 0, q = 1;

        for (int i = 0; i < n; i++)
            for (int j = i + 1; j < n; j++)
                if (fabs(A[indx(i, j, n)]) > max_val) {
                    max_val = fabs(A[indx(i, j, n)]);
                    p = i; q = j;
                }

        if (max_val < EPSILON) break;

        double App = A[indx(p, p, n)];
        double Aqq = A[indx(q, q, n)];
        double Apq = A[indx(p, q, n)];

        double phi = 0.5 * atan2(2 * Apq, Aqq - App);
        double c = cos(phi);
        double s = sin(phi);

        for (int i = 0; i < n; i++) {
            double Aip = A[indx(i, p, n)];
            double Aiq = A[indx(i, q, n)];
            A[indx(i, p, n)] = c * Aip - s * Aiq;
            A[indx(i, q, n)] = s * Aip + c * Aiq;
        }

        for (int j = 0; j < n; j++) {
            double Apj = A[indx(p, j, n)];
            double Aqj = A[indx(q, j, n)];
            A[indx(p, j, n)] = c * Apj - s * Aqj;
            A[indx(q, j, n)] = s * Apj + c * Aqj;
        }

        A[indx(p, q, n)] = A[indx(q, p, n)] = 0.0;
        A[indx(p, p, n)] = c * c * App - 2 * s * c * Apq + s * s * Aqq;
        A[indx(q, q, n)] = s * s * App + 2 * s * c * Apq + c * c * Aqq;

        for (int i = 0; i < n; i++) {
            double Vip = V[indx(i, p, n)];
            double Viq = V[indx(i, q, n)];
            V[indx(i, p, n)] = c * Vip - s * Viq;
            V[indx(i, q, n)] = s * Vip + c * Viq;
        }
    }

    for (int i = 0; i < n; i++)
        eigenvalues[i] = A[indx(i, i, n)];
}


void sort_desc(double *sigma, double *V, int n) {
    for (int i = 0; i < n - 1; i++) {
        int max_idx = i;
        for (int j = i + 1; j < n; j++)
            if (sigma[j] > sigma[max_idx])
                max_idx = j;

        if (max_idx != i) {
            double temp = sigma[i];
            sigma[i] = sigma[max_idx];
            sigma[max_idx] = temp;

            for (int k = 0; k < n; k++) {
                double tmp = V[indx(k, i, n)];
                V[indx(k, i, n)] = V[indx(k, max_idx, n)];
                V[indx(k, max_idx, n)] = tmp;
            }
        }
    }
}


void compute_svd(const double *A, int m, int n,
                 double **U_out, double **S_out, double **V_out, int *r_out) {

    int r = (m < n) ? m : n;
    *r_out = r;

    double *AtA = (double *)malloc(n * n * sizeof(double));
    ata(A, AtA, m, n);

    double *eigen = (double *)malloc(n * sizeof(double));
    double *Vfull = (double *)malloc(n * n * sizeof(double));
    jacobi_algo(AtA, eigen, Vfull, n);

    double *sigma_all = (double *)malloc(n * sizeof(double));
    for (int i = 0; i < n; i++)
        sigma_all[i] = (eigen[i] < 0) ? 0.0 : sqrt(eigen[i]);

    sort_desc(sigma_all, Vfull, n);

    double *S = (double *)malloc(r * sizeof(double));
    for (int i = 0; i < r; i++) S[i] = sigma_all[i];

    double *V = (double *)malloc(n * r * sizeof(double));
    for (int j = 0; j < r; j++)
        for (int i = 0; i < n; i++)
            V[indx(i, j, r)] = Vfull[indx(i, j, n)];

    double *U = (double *)malloc(m * r * sizeof(double));
    for (int j = 0; j < r; j++) {
        double sigma = S[j];
        for (int i = 0; i < m; i++) {
            double sum = 0.0;
            for (int k = 0; k < n; k++)
                sum += A[indx(i, k, n)] * V[indx(k, j, r)];
            U[indx(i, j, r)] = (sigma > 1e-14) ? sum / sigma : 0.0;
        }
    }

    free(AtA); free(eigen); free(Vfull); free(sigma_all);

    *U_out = U;
    *S_out = S;
    *V_out = V;
}


void reconstruct_image(const double *U, const double *S, const double *V,
                       int m, int n, int r, int k, double *Ak) {

    for (int i = 0; i < m * n; i++) Ak[i] = 0.0;

    for (int t = 0; t < k; t++)
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                Ak[indx(i, j, n)] += U[indx(i, t, r)] * S[t] * V[indx(j, t, r)];
}



int main() {
    char filename[100];
    printf("Enter input PGM filename: ");
    scanf("%s", filename);

    double *A;
    int m, n;
    read_pgm(filename, &A, &m, &n);
    

    double *U, *S, *V;
    int r;
    compute_svd(A, m, n, &U, &S, &V, &r);

    

    int k;
    printf("\nEnter truncation rank k (1_%d): ", r);
    scanf("%d", &k);
    

    double *Ak = (double *)malloc(m * n * sizeof(double));
    reconstruct_image(U, S, V, m, n, r, k, Ak);

    char out_name[128];
    sprintf(out_name, "reconstructed_%d.pgm", k);
    write_pgm(out_name, Ak, m, n);

    printf("\nReconstructed image saved as: %s\n", out_name);
    printf("Frobenius error = %.6f\n", frobenius_norm(A, Ak, m, n));
    printf("Relative error= %.4f\n", frobenius_norm(A,Ak,m,n)/frobenius_norm_single(A,m,n));

    free(A);
    free(U);
    free(S);
    free(V);
    free(Ak);

    return 0;
}
