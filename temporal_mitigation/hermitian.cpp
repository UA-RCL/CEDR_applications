#include "inverse.h"

void hermitian(double S[N_TM * M_TM], double Si[N_TM * M_TM], double Shermitian[M_TM * N_TM], double Shermitianimag[M_TM * N_TM]) {
	int i, j;

	for (i = 0; i < N_TM; i++) {
		for (j = 0; j < M_TM; j++) {
			Shermitian[j * N_TM + i] = S[i * M_TM + j];
			//		  cout<<Shermitian[j*N_TM+i]<<" ";
			Shermitianimag[j * N_TM + i] = -Si[i * M_TM + j];
		}
		//		cout<<"\n";
	}

}
