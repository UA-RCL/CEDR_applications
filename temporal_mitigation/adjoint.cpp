#include "inverse.h"

void adjoint(double A[N_TM * N_TM], double adj[N_TM][N_TM]) {
	if (N_TM == 1) {
		adj[0][0] = 1;
		return;
	}

	// temp is used to store cofactors of A[][]
	int sign = 1;
	double temp[N_TM * N_TM];

	for (int i = 0; i < N_TM; i++) {
		for (int j = 0; j < N_TM; j++) {
			// Get cofactor of A[i][j]
			getCofactor(A, temp, i, j, N_TM);

			// sign of adj[j][i] positive if sum of row
			// and column indexes is even.
			sign = ((i + j) % 2 == 0) ? 1 : -1;

			// Interchanging rows and columns to get the
			// transpose of the cofactor matrix
			adj[j][i] = (sign) * (determinant(temp, N_TM - 1));
		}
	}
}
