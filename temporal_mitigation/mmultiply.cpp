#include "inverse.h"

void mmultiply(double A[N_TM * N_TM], double B[N_TM * N_TM], double C[N_TM * N_TM]) {
	int i, j;
	double Abuf[N_TM][N_TM], Bbuf[N_TM][N_TM];
#pragma HLS array_partition variable = Abuf block factor = 4 dim = 2
#pragma HLS array_partition variable = Bbuf block factor = 4 dim = 1

	for (i = 0; i < N_TM; i++) {
		for (j = 0; j < N_TM; j++) {
#pragma HLS PIPELINE
			Abuf[i][j] = A[i * N_TM + j];
			Bbuf[i][j] = B[i * N_TM + j];
		}
	}

	for (i = 0; i < N_TM; i++) {
		for (j = 0; j < N_TM; j++) {
#pragma HLS PIPELINE
			double result1 = 0.0;
			for (int k = 0; k < N_TM; k++) {
				result1 += Abuf[i][k] * Bbuf[k][j];
			}
			C[i * N_TM + j] = result1;
		}
	}
}

void mmultiply2(double A[N_TM * N_TM], double B[N_TM * N_TM], float C[N_TM * N_TM]) {
	int i, j;
	double Abuf[N_TM][N_TM], Bbuf[N_TM][N_TM];
#pragma HLS array_partition variable = Abuf block factor = 4 dim = 2
#pragma HLS array_partition variable = Bbuf block factor = 4 dim = 1

	for (i = 0; i < N_TM; i++) {
		for (j = 0; j < N_TM; j++) {
#pragma HLS PIPELINE
			Abuf[i][j] = A[i * N_TM + j];
			Bbuf[i][j] = B[i * N_TM + j];
		}
	}

	for (i = 0; i < N_TM; i++) {
		for (j = 0; j < N_TM; j++) {
#pragma HLS PIPELINE
			float result1 = 0.0;
			for (int k = 0; k < N_TM; k++) {
				result1 += Abuf[i][k] * Bbuf[k][j];
			}
			C[i * N_TM + j] = result1;
		}
	}
}
