#include "inverse.h"

void realpart(double A[N_TM * N_TM], double B[N_TM * N_TM], double inv1[N_TM * N_TM], double inv2[N_TM * N_TM], double resultreal[N_TM * N_TM]) {
#pragma HLS inline region recursive
	double buffer1[N_TM * N_TM];
	double intmedb1[N_TM * N_TM], intmedb2[N_TM * N_TM];

	mmultiply(B, inv1, intmedb1);
	mmultiply(intmedb1, B, intmedb2);

	for (int i = 0; i < N_TM; i++) {
#pragma HLS loop_tripcount min = 4 max = 10
		for (int j = 0; j < N_TM; j++) {
#pragma HLS loop_tripcount min = 4 max = 10
			buffer1[i * N_TM + j] = A[i * N_TM + j] + intmedb2[i * N_TM + j];
			cout << buffer1[i * N_TM + j] << " ";
			resultreal[i * N_TM + j] = buffer1[i * N_TM + j];
		}
		cout << endl;
	}
	// display(buffer1);

	//	alternateinverse(buffer1,intmedt1);
	//	display(intmedt1,4,4);
}
