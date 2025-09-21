#include "inverse.h"

void mmadd(double A[N_TM * N_TM], double B[N_TM * N_TM], double C[N_TM * N_TM]) {
	double Abuf[N_TM][N_TM], Bbuf[N_TM][N_TM];

	for (int i = 0; i < N_TM; i++) {
		for (int j = 0; j < N_TM; j++) {
#pragma HLS PIPELINE
			Abuf[i][j] = A[i * N_TM + j];
			//	               Aibuf[i][j] = Ai[i*M + j];
			Bbuf[i][j] = B[i * N_TM + j];
			//	               Bibuf[i][j] = Bi[i*M + j];
		}
	}

	for (int i = 0; i < N_TM; i++) {
		for (int j = 0; j < N_TM; j++) {
#pragma HLS PIPELINE
			C[i * N_TM + j] = double(Abuf[i][j] + Bbuf[i][j]);
			//  C[i*N_TM+j]*=-1;
			//	 	        	  Ci[i*M+j]=double(Aibuf[i][j] - Bibuf[i][j]);
		}
	}
}
