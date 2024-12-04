#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define STEPMAX 10000
#define INTV 100
#define L 10.0
#define D 1.0
#define DT 0.002

double u[N+1][N+1] = {0.0};

void writeVtk(int count) {
    FILE *fp;
    char fname[256];
    int i, j;

    sprintf(fname, "data_%08d.vtk", count);
    fp = fopen(fname, "w");
    fprintf(fp, "# vtk DataFile Version 4.1 \n");
    fprintf(fp, "DIFFUSION_EQUATION_2D\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS \n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", N+1, N+1);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %f %f 0 \n", L/N, L/N);
    fprintf(fp, "POINT_DATA %d\n", (N+1)*(N+1));
    fprintf(fp, "SCALARS U double 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (j = 0; j < N+1; j++) {
        for (i = 0; i < N+1; i++) {
            fprintf(fp, "%f\n", u[i][j]);
        }
    }

    fclose(fp);
}

void init() {
    int i, j;

    for (i = 0; i < N+1; i++) {
        for (j = 0; j < N+1; j++) {
            u[i][j] = 0.0;
        }
    }

    for (j = 0; j < N+1; j++)
        u[0][j] = 1.0;
}

int main() {
    int i, j;
    int step;
    double u0[N+1][N+1];
    double dx = L / N;

    init();

    step = 0;
    printf("%d\n", step);
    writeVtk(step);

    for (step = 1; step <= STEPMAX; step++) {
        for (i = 0; i < N+1; i++) {
            for (j = 0; j < N+1; j++)
                u0[i][j] = u[i][j];
        }

        for (i = 1; i < N; i++) {
            for (j = 1; j < N; j++) {
                u[i][j] = u0[i][j] + D / (dx * dx) * 
                          (u0[i-1][j] + u0[i+1][j] + u0[i][j-1] + u0[i][j+1] - 4.0 * u0[i][j]) * DT;
            }
        }

        for (j = 0; j < N+1; j++) {
            u[0][j] = 1.0;
            u[N][j] = 0.0;
        }
        for (i = 0; i < N+1; i++) {
            u[i][0] = 0.0;
            u[i][N] = 0.0;
        }
        if (step % INTV == 0) {
            printf("%d\n", step);
            writeVtk(step);
        }
    }

    return 0;
}