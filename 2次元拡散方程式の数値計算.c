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
    sprintf(fname, "data_%08d.vtk", count);
    fp = fopen(fname, "w");
    fprintf(fp, "# vtk DataFile Version 3.0\n");
    fprintf(fp, "2D Diffusion Equation\n");
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET STRUCTURED_POINTS\n");
    fprintf(fp, "DIMENSIONS %d %d 1\n", N+1, N+1);
    fprintf(fp, "ORIGIN 0 0 0\n");
    fprintf(fp, "SPACING %f %f 1\n", L / N, L / N);
    fprintf(fp, "POINT_DATA %d\n", (N+1) * (N+1));
    fprintf(fp, "SCALARS concentration float 1\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (int j = 0; j <= N; j++) {
        for (int i = 0; i <= N; i++) {
            fprintf(fp, "%f\n", u[i][j]);
        }
    }
    fclose(fp);
}

void init() {
    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            u[i][j] = 0.0;
        }
    }

    for (int i = N/2 - 5; i <= N/2 + 5; i++) {
        for (int j = N/2 - 5; j <= N/2 + 5; j++) {
            u[i][j] = 1.0;
        }
    }
}

void update() {
    double temp[N+1][N+1];

    for (int i = 0; i <= N; i++) {
        for (int j = 0; j <= N; j++) {
            temp[i][j] = u[i][j];
        }
    }

    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            u[i][j] = temp[i][j] + D * DT / ((L / N) * (L / N)) * 
                       (temp[i-1][j] + temp[i+1][j] + temp[i][j-1] + temp[i][j+1] - 4.0 * temp[i][j]);
        }
    }

    for (int i = 0; i <= N; i++) {
        u[0][i] = 0.0;
        u[N][i] = 0.0;
        u[i][0] = 0.0;
        u[i][N] = 0.0;
    }
}

int main() {
    int step = 0;

    init();
    writeVtk(step);

    for (step = 1; step <= STEPMAX; step++) {
        update();

        if (step % INTV == 0) {
            printf("%d\n", step);
            writeVtk(step);
        }
    }

    return 0;
}
