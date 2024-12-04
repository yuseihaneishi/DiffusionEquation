#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 100
#define STEPMAX 10000
#define INTV 100
#define L 10.0
#define D 1.0
#define DT 0.002

double u[N+1] = {0.0};
double u_new[N+1] = {0.0};

void writeVtk(int count) {

	FILE *fp;
	char fname[256];
	double uu[N+1][150] = {0};
	int i, j;

	for (i = 0; i < N + 1; i++) {
		j = (int)(u[i] / 0.01) + 10;
		uu[i][j] = 1.0;
	}

	sprintf(fname, "data_%08d.vtk", count);
	fp = fopen(fname, "w");
	fprintf(fp, "# vtk DataFile Version 4.1 \n");
	fprintf(fp, "DIFFUSION_EQUATION_1D\n");
	fprintf(fp, "ASCII\n");
	fprintf(fp, "DATASET STRUCTURED_POINTS \n");
	fprintf(fp, "DIMENSIONS %d %d 1\n", N + 1, 150);
	fprintf(fp, "ORIGIN 0 0 0\n");
	fprintf(fp, "SPACING %f %f 0 \n", 1.0, 1.0);
	fprintf(fp, "POINT_DATA %d\n", (N + 1) * 150);
	fprintf(fp, "SCALARS U double 1\n");
	fprintf(fp, "LOOKUP_TABLE default\n");
	for (j = 0; j < 150; j++) {
		for (i = 0; i < N + 1; i++) {
			fprintf(fp, "%f\n", uu[i][j]);
		}
	}

	fclose(fp);
}

int main() {

	int step;
	double dx = L / N;
	double alpha = D * DT / (dx * dx);

	if (alpha > 0.5) {
		printf("安定条件× alpha = %f\n", alpha);
		return -1;
	}

	for (int i = 0; i <= N; i++)
		u[i] = 0.0; 
	int index_4_9 = (int)(4.9 / dx + 0.5);
	int index_5_0 = (int)(5.0 / dx + 0.5);
	int index_5_1 = (int)(5.1 / dx + 0.5);
	u[index_4_9] = 1.0;
	u[index_5_0] = 1.0;
	u[index_5_1] = 1.0;

	step = 0;
	printf("Step: %d\n", step);
	writeVtk(step);

	for (step = 1; step <= STEPMAX; step++) {
		for (int i = 1; i < N; i++)
			u_new[i] = u[i] + alpha * (u[i + 1] - 2 * u[i] + u[i - 1]);

		u_new[0] = 0.0;
		u_new[N] = 0.0;

		for (int i = 0; i <= N; i++)
			u[i] = u_new[i];

		if (step % INTV == 0) {
			printf("Step: %d\n", step);
			writeVtk(step);
		}
	}

	return 0;
}
