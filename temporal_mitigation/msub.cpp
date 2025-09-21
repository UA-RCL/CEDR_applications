#include "inverse.h"

void msub(double A[N_TM * M_TM], double Ai[N_TM * M_TM], double B[N_TM * M_TM], double Bi[N_TM * M_TM], double C[N_TM * M_TM], double Ci[N_TM * M_TM]) {
	double Abuf[N_TM][M_TM], Aibuf[N_TM][M_TM], Bbuf[N_TM][M_TM], Bibuf[N_TM][M_TM];

	for (int i = 0; i < N_TM; i++) {
		for (int j = 0; j < M_TM; j++) {
#pragma HLS PIPELINE
			Abuf[i][j] = A[i * M_TM + j];
			Aibuf[i][j] = Ai[i * M_TM + j];
			Bbuf[i][j] = B[i * M_TM + j];
			Bibuf[i][j] = Bi[i * M_TM + j];
		}
	}

	for (int i = 0; i < N_TM; i++) {
		for (int j = 0; j < M_TM; j++) {
#pragma HLS PIPELINE
			C[i * M_TM + j] = double(Abuf[i][j] - Bbuf[i][j]);
			Ci[i * M_TM + j] = double(Aibuf[i][j] - Bibuf[i][j]);
		}
	}
}
