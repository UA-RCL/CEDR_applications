// C++ program to find adjoint and inverse of a matrix
//#include<bits/stdc++.h>
#include "inverse.h"

// Function to calculate and store inverse, returns false if
// matrix is singular
void inverse(double A[N_TM * N_TM], double inverse[N_TM * N_TM]) {
	//#pragma HLS inline region recursive
	// Find determinant of A[][]
	double det = determinant(A, N_TM);
	cout << "Value of determinant is " << det;
	if (det == 0) {
		cout << "Singular matrix, can't find its inverse";
		return;
	}

	// Find adjoint
	double adj[N_TM][N_TM], result = 0.0;
	adjoint(A, adj);

	// Find Inverse using formula "inverse(A) = adj(A)/det(A)"
	for (int i = 0; i < N_TM; i++) {
#pragma HLS loop_tripcount min = 4 max = 10
		for (int j = 0; j < N_TM; j++) {
#pragma HLS loop_tripcount min = 4 max = 10
			result = divide(adj[i][j], det);
			//        	cout<<"Result from division of "<<adj[i][j]<<"with"<<det<<"is"<<result;
			//        	g=adj[i][j]/det;
			//#pragma HLS RESOURCE variable=g core=FDiv
			inverse[i * N_TM + j] = result;
		}
	}
	//    return true;
	cout << "Returning from the inverse function";
}
