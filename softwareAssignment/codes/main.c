#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include"pgmfinctions.h"
#include"matrix.h"
#include"algo.h"
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
