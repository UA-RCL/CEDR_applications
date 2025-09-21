#include "inverse.h"

#include "libcedr.h"

#define PROGPATH "./input/"
#define ZIN PROGPATH "z.txt"
#define ZIMAGIN PROGPATH "zimag.txt"
#define SIN PROGPATH "s.txt"
#define SIMAGIN PROGPATH "simag.txt"

#define KERN_ENTER(STR)
#define KERN_EXIT(STR)

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <numeric>

int main(int argc, char **argv) {
  int N = N_TM;
  int M = M_TM;
  // Initializing the main variables that will be required in the overall
  // computation of the Z response

  // Initializing the Z signal which will have 4*64 dimension
  double *Z, *Zi;
  Z = (double *)malloc(N * M * sizeof(double));
  Zi = (double *)malloc(N * M * sizeof(double));

  // Now defining the jammer signal which will have the same dimensions as the
  // message signal , The jammer is denoted by S
  double *S, *Si;
  S = (double *)malloc(N * M * sizeof(double));
  Si = (double *)malloc(N * M * sizeof(double));

  // now defining the argument files which will contain the corresponding values
  // of Z and S
  FILE *Zreal, *Zimag, *Sreal, *Simag, *output;
  Zreal = fopen(ZIN, "r");
  Zimag = fopen(ZIMAGIN, "r");
  Sreal = fopen(SIN, "r");
  Simag = fopen(SIMAGIN, "r");
  
  if (Zreal == NULL) {
      printf("[ERROR] Error in opening Zreal file. Exiting...\n");
      exit(1);
  }
  if (Zimag == NULL) {
      printf("[ERROR] Error in opening Zimag file. Exiting...\n");
      exit(1);
  }
  if (Sreal == NULL) {
      printf("[ERROR] Error in opening Sreal file. Exiting...\n");
      exit(1);
  }
  if (Simag == NULL) {
      printf("[ERROR] Error in opening Simag file. Exiting...\n");
      exit(1);
  }

  // now copying the contents of the files into the arrays that have been
  // assigned for the signal and the jammer
  int scan_check;
  for (int i = 0; i < N; i++) {
    for (int j = 0; j < M; j++) {
      scan_check = fscanf(Zreal, "%lf", &Z[i * M + j]);
      Z[i * M + j] /= 10.0f;
      scan_check = fscanf(Zimag, "%lf", &Zi[i * M + j]);
      Zi[i * M + j] /= 10.0f;
      scan_check = fscanf(Sreal, "%lf", &S[i * M + j]);
      S[i * M + j] /= 10.0f;
      scan_check = fscanf(Simag, "%lf", &Si[i * M + j]);
      Si[i * M + j] /= 10.0f;
    }
  }
  fclose(Zreal);
  fclose(Zimag);
  fclose(Sreal);
  fclose(Simag);

  // Computing the hermitian of S
  double *Shermitian, *Shermitianimag;
  Shermitian = (double *)malloc(M * N * sizeof(double));
  Shermitianimag = (double *)malloc(M * N * sizeof(double));
  // Now computing the result from the first multiplication (Z*S^H)--> Result 1
  double *result1, *result1imag;
  result1 = (double *)malloc(N * N * sizeof(double));
  result1imag = (double *)malloc(N * N * sizeof(double));
  // Now computing the second matrix multiplication (S*S^H) ---> Result2
//  double *result2, *result2imag;
  double *result2imag;
  cedr_re_flt_type *result2;
  result2 = (cedr_re_flt_type *)malloc(N * N * sizeof(cedr_re_flt_type));
  result2imag = (double *)malloc(N * N * sizeof(double));
  double *inv1, *inv2, *intmedt1, *intmedt2, *intmedb1, *intmedb2, *intmedb3,
      *intmedb4;
  double *buffer1, *buffer2, *buffer3, *buffer4;
//  double *bufferinv1, *bufferinv2, *bufferinv3, *bufferinv4, *bufferinv5;
  double *bufferinv1, *bufferinv4, *bufferinv5;
  cedr_re_flt_type *bufferinv2, *bufferinv3;
  // To store inverse of A[][]
  inv1 = (double *)malloc(N * N * sizeof(double));
  inv2 = (double *)malloc(N * N * sizeof(double));
  intmedt1 = (double *)malloc(N * N * sizeof(double));
  intmedt2 = (double *)malloc(N * N * sizeof(double));
  intmedb1 = (double *)malloc(N * N * sizeof(double));
  intmedb2 = (double *)malloc(N * N * sizeof(double));
  intmedb3 = (double *)malloc(N * N * sizeof(double));
  intmedb4 = (double *)malloc(N * N * sizeof(double));
  buffer1 = (double *)malloc(N * N * sizeof(double));
  buffer2 = (double *)malloc(N * N * sizeof(double));
  buffer3 = (double *)malloc(N * N * sizeof(double));
  buffer4 = (double *)malloc(N * N * sizeof(double));
  // The following arrays are used for the inverse computation
  bufferinv1 = (double *)malloc(N * N * sizeof(double));
  bufferinv2 = (cedr_re_flt_type *)malloc(N * N * sizeof(cedr_re_flt_type));
  bufferinv3 = (cedr_re_flt_type *)malloc(N * N * sizeof(cedr_re_flt_type));
  bufferinv4 = (double *)malloc(N * N * sizeof(double));
  bufferinv5 = (double *)malloc(N * N * sizeof(double));
  // Now computing the result of (Z*S^H)*(S.S^H)^-1  ---> result3 which is a 4*4
  // and 4*4 multiplication
  double *result3, *result3imag;
  result3 = (double *)malloc(N * N * sizeof(double));
  result3imag = (double *)malloc(N * N * sizeof(double));
  // Now computing the final matrix multiplication which is result3*S
  // == result4, this is 4*4 and 4*64 multiplication
  double *result4, *result4imag;
  result4 = (double *)malloc(N * M * sizeof(double));
  result4imag = (double *)malloc(N * M * sizeof(double));
  // Now we have to compute the final operation which is matrix subtraction : (Z
  // - result4) ---> Zr  4*64 - 4*64
  double *zres, *zresimag;
  zres = (double *)malloc(N * M * sizeof(double));
  zresimag = (double *)malloc(N * M * sizeof(double));

  struct timespec start, end;
  std::vector<long long> elapsed_times;
  long long elapsed_time;
  int num_measurements = 1000;
  for(int iii=0;iii<1000;iii++){
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  hermitian(S, Si, Shermitian, Shermitianimag);

  /******** API Change ********/
  // mmult(Z, Zi, Shermitian, Shermitianimag, result1, result1imag);
  // lol scope
  {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N * M * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(M * N * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N * N * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N * M; i++) {
      A_api[i].re = (cedr_re_flt_type) Z[i];
      A_api[i].im = (cedr_re_flt_type) Zi[i];
      B_api[i].re = (cedr_re_flt_type) Shermitian[i];
      B_api[i].im = (cedr_re_flt_type) Shermitianimag[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N, M, N, 0,0,0,0, .argc=argc, .argv=argv);

    for (size_t i = 0; i < N * N; i++) {
      result1[i] = (double) C_api[i].re;
      result1imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }
  /******** API Change ********/

  /******** API Change ********/
  // mmult(S, Si, Shermitian, Shermitianimag, result2, result2imag);
  // lol scope
  {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N * M * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(M * N * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N * N * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N * M; i++) {
      A_api[i].re = (cedr_re_flt_type) S[i];
      A_api[i].im = (cedr_re_flt_type) Si[i];
      B_api[i].re = (cedr_re_flt_type) Shermitian[i];
      B_api[i].im = (cedr_re_flt_type) Shermitianimag[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N, M, N, 0,0,0,0, .argc=argc, .argv=argv);

    for (size_t i = 0; i < N * N; i++) {
      result2[i] = C_api[i].re;
      result2imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }
  /******** API Change ********/


  // Computing the inverse of a == Real part
  alternateinverse2((float *)result2, inv1);

  // compute res1 =inv(inv1)*result2imag == > bufferinv1
  mmultiply(inv1, result2imag, bufferinv1);

  // compute res2 = inv(result2imag*res1+result2)
  // result2imag*bufferinv1 == > bufferinv2
  mmultiply2(result2imag, bufferinv1, (float *) bufferinv2);

  // bufferinv2+result2 ==> bufferinv3
  CEDR_ZIP_flt((cedr_cmplx_flt_type*)bufferinv2, (cedr_cmplx_flt_type*)result2, (cedr_cmplx_flt_type*)bufferinv3, N*N/2, ZIP_ADD, .argc=argc, .argv=argv);

  // Now computing inv(bufferinv3) ==> bufferinv4
  // The following 'bufferinv4' is real part of the inverse
  alternateinverse2((float *)bufferinv3, bufferinv4);

  // Now computing the imaginary part of the inverse
  mmultiply(bufferinv1, bufferinv4, bufferinv5);
  for (int k = 0; k < 4; k++) {
    for (int m = 0; m < 4; m++) {
      bufferinv5[k * 4 + m] *= -1;
    }
  }
  // Final result is -bufferinv5 == > imag part of inverse


  /******** API Change ********/
  // mmult4(result1, result1imag, bufferinv4, bufferinv5, result3, result3imag);
  // lol scope
  {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N * N * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(N * N * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N * N * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N * N; i++) {
      A_api[i].re = (cedr_re_flt_type) result1[i];
      A_api[i].im = (cedr_re_flt_type) result1imag[i];
      B_api[i].re = (cedr_re_flt_type) bufferinv4[i];
      B_api[i].im = (cedr_re_flt_type) bufferinv5[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N, N, N, 0,0,0,0, .argc=argc, .argv=argv);

    for (size_t i = 0; i < N * N; i++) {
      result3[i] = (double) C_api[i].re;
      result3imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }
  /******** API Change ********/

  /******** API Change ********/
  // mmult64(result3, result3imag, S, Si, result4, result4imag);
  // lol scope
  {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N * N * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(N * M * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N * M * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N * N; i++) {
      A_api[i].re = (cedr_re_flt_type) result3[i];
      A_api[i].im = (cedr_re_flt_type) result3imag[i];
    }
    for (size_t i = 0; i < N * M; i++) {
      B_api[i].re = (cedr_re_flt_type) S[i];
      B_api[i].im = (cedr_re_flt_type) Si[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N, N, M, 0,0,0,0, .argc=argc, .argv=argv);

    for (size_t i = 0; i < N * M; i++) {
      result4[i] = (double) C_api[i].re;
      result4imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }
  /******** API Change ********/


  msub(Z, Zi, result4, result4imag, zres, zresimag);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec);
    elapsed_times.push_back(elapsed_time);
  }
  long long min = *std::min_element(elapsed_times.begin(), elapsed_times.end());
  long long max = *std::max_element(elapsed_times.begin(), elapsed_times.end());
  double average = std::accumulate(elapsed_times.begin(), elapsed_times.end(), 0.0) / num_measurements;
  std::sort(elapsed_times.begin(), elapsed_times.end());
  double median = num_measurements % 2 == 0
                  ? (elapsed_times[num_measurements / 2 - 1] + elapsed_times[num_measurements / 2]) / 2.0
                  : elapsed_times[num_measurements / 2];
  printf("Average: %.0f ns\n", average);
  printf("Median: %.0f ns\n", median);
  printf("Min: %lld ns\n", min);
  printf("Max: %lld ns\n", max);
  printf("%.0f\n%.0f\n%lld\n%lld\n", average, median, min, max);
  
  // Printing in output file
  output = fopen("./output/temporal_mitigation_output.txt", "w");
  if (output != NULL) {
    fprintf(output, "Real part: \n");
    for (int r = 0; r < N; r++) {
      for (int c = 0; c < M; c++)
        fprintf(output, "%lf, ", zres[(r * M) + c]);
      fprintf(output, "\n");
    }
    fprintf(output, "\nImag part: \n");
    for (int r = 0; r < N; r++) {
      for (int c = 0; c < M; c++)
        fprintf(output, "%lf, ", zresimag[(r * M) + c]);
      fprintf(output, "\n");
    }
    fclose(output);
  }

  // Now freeing up the variables used
  free(Z);
  free(Zi);
  free(S);
  free(Si);
  free(Shermitian);
  free(Shermitianimag);
  free(result1);
  free(result1imag);
  free(result2);
  free(result2imag);
  free(inv1);
  free(inv2);
  free(bufferinv1);
  free(bufferinv2);
  free(bufferinv3);
  free(bufferinv4);
  free(bufferinv5);
  free(intmedt1);
  free(intmedt2);
  free(result3);
  free(result3imag);
  free(result4);
  free(result4imag);
  free(zres);
  free(zresimag);

  // cout << "[Temporal Mitigation] Execution is complete..." << std::endl;

  return 0;
}
