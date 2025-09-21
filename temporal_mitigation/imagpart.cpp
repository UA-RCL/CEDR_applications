#include "inverse.h"

void imagpart(double A[N_TM * N_TM], double B[N_TM * N_TM], double inv1[N_TM * N_TM], double inv2[N_TM * N_TM], double resultimag[N_TM * N_TM]) {
#pragma HLS inline region recursive
	double buffer2[N_TM * N_TM];
	double intmedbuff1[N_TM * N_TM], intmedbuff2[N_TM * N_TM];

	mmultiply(A, inv2, intmedbuff1);
	mmultiply(intmedbuff1, A, intmedbuff2);

	for (int i = 0; i < N_TM; i++) {
		for (int j = 0; j < N_TM; j++) {
			buffer2[i * N_TM + j] = B[i * N_TM + j] + intmedbuff2[i * N_TM + j];
			buffer2[i * N_TM + j] *= -1;
			resultimag[i * N_TM + j] = buffer2[i * N_TM + j];
			cout << buffer2[i * N_TM + j] << " ";
		}
		cout << endl;
	}
	// display(buffer1);

	//	alternateinverse(buffer2,intmedt2);
	//	display(intmedt2);
}
