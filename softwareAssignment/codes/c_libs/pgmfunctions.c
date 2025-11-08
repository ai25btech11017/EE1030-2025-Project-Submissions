#include<stdio.h>
#include"pgmfunctions.h"
void read_pgm(const char *filename, double **image, int *rows, int *cols) {
    FILE *file = fopen(filename, "rb");

    char magic[3];
    fscanf(file, "%2s", magic);
    

    int isP2 = !strcmp(magic, "P2");
    int isP5 = !strcmp(magic, "P5");
    

   
    int ch = fgetc(file);
    while (ch == '#') {
        while (ch != '\n' && ch != EOF) ch = fgetc(file);
        ch = fgetc(file);
    }
    ungetc(ch, file);

    int width, height, maxval;
    fscanf(file, "%d %d %d", &width, &height, &maxval);
    fgetc(file); 

    *rows = height;
    *cols = width;
    *image = (double *)malloc(height * width * sizeof(double));

    if (isP2) {
        for (int i = 0; i < height * width; i++) {
            int val;
            fscanf(file, "%d", &val);
            (*image)[i] = val;
        }
    } else {  
        for (int i = 0; i < height * width; i++) {
            unsigned char val;
            fread(&val, 1, 1, file);
            (*image)[i] = val;
        }
    }

    fclose(file);
}


void write_pgm(const char *filename, const double *image, int rows, int cols) {
    FILE *file = fopen(filename, "w");
    

    fprintf(file, "P2\n%d %d\n255\n", cols, rows);

    for (int i = 0; i < rows * cols; i++) {
        int val = (int)round(image[i]);
        if (val < 0) val = 0;
        if (val > 255) val = 255;

        fprintf(file, "%d ", val);
        if ((i + 1) % cols == 0) fprintf(file, "\n");
    }

    fclose(file);
}
