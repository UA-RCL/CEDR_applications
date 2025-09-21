#ifndef _INVERSE_H_
#define _INVERSE_H_

#include <bits/stdc++.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include <complex>
#include <iostream>

using namespace std;
//#ifndef _INVERSE_H_
//#define _INVERSE_H_

#define N_TM  4
#define M_TM 64

void getCofactor(double A[N_TM * N_TM], double temp[N_TM * N_TM], int p, int q, int n);
#pragma SDS data data_mover(m : AXIDMA_SIMPLE, invOut : AXIDMA_SIMPLE)
void alternateinverse(double m[16], double invOut[16]);
void alternateinverse2(float m[16], double invOut[16]);
//#pragma SDS data data_mover(m:AXIDMA_SIMPLE,invOut:AXIDMA_SIMPLE)
void scalableinverse(double m[16], double invOut[16]);
double determinant(double A[N_TM * N_TM], int n);
void adjoint(double A[N_TM * N_TM], double adj[N_TM][N_TM]);
void inverse(double A[N_TM * N_TM], double inverse[N_TM * N_TM]);
void display(double *A, int, int);
// void realpart(double A[N_TM *N_TM], double B[N_TM *N_TM] , double inv1[N_TM *N_TM] ,double inv2[N_TM *N_TM] , double intmedt1[N_TM *N_TM]);
// void imagpart(double A[N_TM *N_TM],double B[N_TM *N_TM],double inv1[N_TM *N_TM], double inv2[N_TM *N_TM],double intmedt2[N_TM *N_TM]);
#pragma SDS data access_pattern(A : SEQUENTIAL, B : SEQUENTIAL, C : SEQUENTIAL)
void mmultiply(double A[N_TM * N_TM], double B[N_TM * N_TM], double C[N_TM * N_TM]);
void mmultiply2(double A[N_TM * N_TM], double B[N_TM * N_TM], float C[N_TM * N_TM]);
//#pragma SDS data access_pattern(S:SEQUENTIAL ,Si:SEQUENTIAL,Shermitian:SEQUENTIAL , Shermitianimag:SEQUENTIAL)
void hermitian(double S[N_TM * M_TM], double Si[N_TM * M_TM], double Shermitian[N_TM * M_TM], double Shermitianimag[N_TM * M_TM]);
double divide(double, double);

//#pragma SDS data access_pattern(A:SEQUENTIAL,B:SEQUENTIAL,C:SEQUENTIAL)
void mmadd(double A[N_TM * N_TM], double B[N_TM * N_TM], double C[N_TM * N_TM]);

#pragma SDS data access_pattern(A                \
                                : SEQUENTIAL, Ai \
                                : SEQUENTIAL, B  \
                                : SEQUENTIAL, Bi \
                                : SEQUENTIAL, C  \
                                : SEQUENTIAL, Ci \
                                : SEQUENTIAL)
void mmult_f(double A[N_TM * N_TM], double Ai[N_TM * N_TM], double B[N_TM * N_TM], double Bi[N_TM * N_TM], double C[N_TM * N_TM], double Ci[N_TM * N_TM]);

#pragma SDS data access_pattern(A                \
                                : SEQUENTIAL, Ai \
                                : SEQUENTIAL, B  \
                                : SEQUENTIAL, Bi \
                                : SEQUENTIAL, C  \
                                : SEQUENTIAL, Ci \
                                : SEQUENTIAL)
void mmult4(double A[N_TM * N_TM], double Ai[N_TM * N_TM], double B[N_TM * N_TM], double Bi[N_TM * N_TM], double C[N_TM * N_TM], double Ci[N_TM * N_TM]);

#pragma SDS data access_pattern(A                \
                                : SEQUENTIAL, Ai \
                                : SEQUENTIAL, B  \
                                : SEQUENTIAL, Bi \
                                : SEQUENTIAL, C  \
                                : SEQUENTIAL, Ci \
                                : SEQUENTIAL)
void mmult64(double A[N_TM * N_TM], double Ai[N_TM * N_TM], double B[N_TM * M_TM], double Bi[N_TM * M_TM], double C[N_TM * M_TM], double Ci[N_TM * M_TM]);

#pragma SDS data access_pattern(A                \
                                : SEQUENTIAL, Ai \
                                : SEQUENTIAL, B  \
                                : SEQUENTIAL, Bi \
                                : SEQUENTIAL, C  \
                                : SEQUENTIAL, Ci \
                                : SEQUENTIAL)
void msub(double A[N_TM * M_TM], double Ai[N_TM * M_TM], double B[N_TM * M_TM], double Bi[N_TM * M_TM], double C[N_TM * M_TM], double Ci[N_TM * M_TM]);

// void inverse()

#endif
