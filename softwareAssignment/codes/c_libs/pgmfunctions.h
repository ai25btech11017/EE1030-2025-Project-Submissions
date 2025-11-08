#ifndef pgmfunctions_h
#define pgmfunctions_h

void read_pgm(const char *filename, double **image, int *rows, int *cols);

void write_pgm(const char *filename, const double *image, int rows, int cols);

#endif
