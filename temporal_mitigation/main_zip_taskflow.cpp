#include "inverse.h"

#include "libcedr.h"
#include "libtaskflow.h"

#define PROGPATH "./input/"
#define ZIN PROGPATH "z.txt"
#define ZIMAGIN PROGPATH "zimag.txt"
#define SIN PROGPATH "s.txt"
#define SIMAGIN PROGPATH "simag.txt"

#define KERN_ENTER(STR)
#define KERN_EXIT(STR)

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <numeric>

int main(int argc, char **argv) {
  int N = N_TM;
  int M = M_TM;
  int N_temp = N_TM;
  int M_temp = M_TM;
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

  tf::Taskflow taskflow;
  std::map<std::string, cedr_task_config_t*> task_configs;
  
  taskflow.name("Temporal Mitigation");

  // Computing the hermitian of S
  double *Shermitian, *Shermitianimag;
  Shermitian = (double *)malloc(M * N * sizeof(double));
  Shermitianimag = (double *)malloc(M * N * sizeof(double));

  auto t_hermitian = taskflow.emplace([S, Si, &Shermitian, &Shermitianimag]() {
    hermitian(S, Si, Shermitian, Shermitianimag);
  }).name("hermitian-non_kernel");

  // Now computing the result from the first multiplication (Z*S^H)--> Result 1
  double *result1, *result1imag;
  result1 = (double *)malloc(N * N * sizeof(double));
  result1imag = (double *)malloc(N * N * sizeof(double));

  /******** API Change ********/
  // mmult(Z, Zi, Shermitian, Shermitianimag, result1, result1imag);
  // lol scope
  task_configs["gemm0-CEDR_GEMM_flt-1"] = NULL;
  cedr_task_config_t *&confg0 = task_configs["gemm0-CEDR_GEMM_flt-1"];
  auto t_gemm0 = taskflow.emplace([Z, Zi, Shermitian, Shermitianimag, &result1, &result1imag, N_temp, M_temp, argc, argv, &confg0]() {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N_temp * M_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(M_temp * N_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N_temp * N_temp * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N_temp * M_temp; i++) {
      A_api[i].re = (cedr_re_flt_type) Z[i];
      A_api[i].im = (cedr_re_flt_type) Zi[i];
      B_api[i].re = (cedr_re_flt_type) Shermitian[i];
      B_api[i].im = (cedr_re_flt_type) Shermitianimag[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, M_temp, N_temp, 0,0,0,0, .argc=argc, .argv=argv);
    //CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, M_temp, N_temp, 0,0,0,0, .argc=argc, .argv=argv, .config=confg0, .task_id=0);

    for (size_t i = 0; i < N_temp * N_temp; i++) {
      result1[i] = (double) C_api[i].re;
      result1imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }).name("gemm0-CEDR_GEMM_flt-1");
  /******** API Change ********/

  // Now computing the second matrix multiplication (S*S^H) ---> Result2
//  double *result2, *result2imag;
  double *result2imag;
  cedr_re_flt_type *result2;
  result2 = (cedr_re_flt_type *)malloc(N * N * sizeof(cedr_re_flt_type));
  result2imag = (double *)malloc(N * N * sizeof(double));

  /******** API Change ********/
  // mmult(S, Si, Shermitian, Shermitianimag, result2, result2imag);
  // lol scope
  task_configs["gemm1-CEDR_GEMM_flt-1"] = NULL;
  cedr_task_config_t *&confg1 = task_configs["gemm1-CEDR_GEMM_flt-1"];
  auto t_gemm1 = taskflow.emplace([S, Si, Shermitian, Shermitianimag, &result2, &result2imag, N_temp, M_temp, argc, argv, &confg1]() {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N_temp * M_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(M_temp * N_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N_temp * N_temp * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N_temp * M_temp; i++) {
      A_api[i].re = (cedr_re_flt_type) S[i];
      A_api[i].im = (cedr_re_flt_type) Si[i];
      B_api[i].re = (cedr_re_flt_type) Shermitian[i];
      B_api[i].im = (cedr_re_flt_type) Shermitianimag[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, M_temp, N_temp, 0,0,0,0, .argc=argc, .argv=argv);
    //CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, M_temp, N_temp, 0,0,0,0, .argc=argc, .argv=argv, .config=confg1, .task_id=0);

    for (size_t i = 0; i < N_temp * N_temp; i++) {
      result2[i] = C_api[i].re;
      result2imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }).name("gemm1-CEDR_GEMM_flt-1");
  /******** API Change ********/

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

  // Computing the inverse of a == Real part
  auto t_alternateinverse = taskflow.emplace([result2, &inv1]() {
    alternateinverse2((float *)result2, inv1);
  }).name("alternateinverse-non_kernel");

  // compute res1 =inv(inv1)*result2imag == > bufferinv1
  auto t_mmultiply = taskflow.emplace([inv1, result2imag, &bufferinv1]() {
    mmultiply(inv1, result2imag, bufferinv1);
  }).name("mmultiply-non_kernel");

  // compute res2 = inv(result2imag*res1+result2)
  // result2imag*bufferinv1 == > bufferinv2
  auto t_mmultiply2 = taskflow.emplace([result2imag, bufferinv1, &bufferinv2]() {
    mmultiply2(result2imag, bufferinv1, (float *) bufferinv2);
  }).name("mmultiply2-non_kernel");

  // bufferinv2+result2 ==> bufferinv3
  task_configs["zip1-CEDR_GEMM_flt-1"] = NULL;
  cedr_task_config_t *&confz1 = task_configs["zip1-CEDR_GEMM_flt-1"];
  auto t_zip1 = taskflow.emplace([bufferinv2, result2, &bufferinv3, N_temp, argc, argv, &confz1]() {
    CEDR_ZIP_flt((cedr_cmplx_flt_type*)bufferinv2, (cedr_cmplx_flt_type*)result2, (cedr_cmplx_flt_type*)bufferinv3, N_temp*N_temp/2, ZIP_ADD, .argc=argc, .argv=argv);
    //CEDR_ZIP_flt((cedr_cmplx_flt_type*)bufferinv2, (cedr_cmplx_flt_type*)result2, (cedr_cmplx_flt_type*)bufferinv3, N_temp*N_temp/2, ZIP_ADD, .argc=argc, .argv=argv, .config=confz1, .task_id=0);
  }).name("zip1-CEDR_ZIP_flt-1");

  // Now computing inv(bufferinv3) ==> bufferinv4
  // The following 'bufferinv4' is real part of the inverse
  auto t_alternateinverse2 = taskflow.emplace([bufferinv3, &bufferinv4]() {
    alternateinverse2((float *)bufferinv3, bufferinv4);
  }).name("alternateinverse2-non_kernel");

  // Now computing the imaginary part of the inverse
  auto t_mmultiply0 = taskflow.emplace([bufferinv1, bufferinv4, &bufferinv5]() {
    mmultiply(bufferinv1, bufferinv4, bufferinv5);
    for (int k = 0; k < 4; k++) {
      for (int m = 0; m < 4; m++) {
        bufferinv5[k * 4 + m] *= -1;
      }
    }
  }).name("mmultiply0-non_kernel");
  // Final result is -bufferinv5 == > imag part of inverse

  // Now computing the result of (Z*S^H)*(S.S^H)^-1  ---> result3 which is a 4*4
  // and 4*4 multiplication
  double *result3, *result3imag;
  result3 = (double *)malloc(N * N * sizeof(double));
  result3imag = (double *)malloc(N * N * sizeof(double));

  /******** API Change ********/
  // mmult4(result1, result1imag, bufferinv4, bufferinv5, result3, result3imag);
  // lol scope
  task_configs["gemm2-CEDR_GEMM_flt-1"] = NULL;
  cedr_task_config_t *&confg2 = task_configs["gemm2-CEDR_GEMM_flt-1"];
  auto t_gemm2 = taskflow.emplace([result1, result1imag, bufferinv4, bufferinv5, &result3, &result3imag, N_temp, argc, argv, &confg2]() {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N_temp * N_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(N_temp * N_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N_temp * N_temp * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N_temp * N_temp; i++) {
      A_api[i].re = (cedr_re_flt_type) result1[i];
      A_api[i].im = (cedr_re_flt_type) result1imag[i];
      B_api[i].re = (cedr_re_flt_type) bufferinv4[i];
      B_api[i].im = (cedr_re_flt_type) bufferinv5[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, N_temp, N_temp, 0,0,0,0, .argc=argc, .argv=argv);
    //CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, N_temp, N_temp, 0,0,0,0, .argc=argc, .argv=argv, .config=confg2, .task_id=0);

    for (size_t i = 0; i < N_temp * N_temp; i++) {
      result3[i] = (double) C_api[i].re;
      result3imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }).name("gemm2-CEDR_GEMM_flt-1");
  /******** API Change ********/

  // Now computing the final matrix multiplication which is result3*S
  // == result4, this is 4*4 and 4*64 multiplication
  double *result4, *result4imag;
  result4 = (double *)malloc(N * M * sizeof(double));
  result4imag = (double *)malloc(N * M * sizeof(double));

  /******** API Change ********/
  // mmult64(result3, result3imag, S, Si, result4, result4imag);
  // lol scope
  task_configs["gemm3-CEDR_GEMM_flt-1"] = NULL;
  cedr_task_config_t *&confg3 = task_configs["gemm3-CEDR_GEMM_flt-1"];
  auto t_gemm3 = taskflow.emplace([result3, result3imag, S, Si, &result4, &result4imag, N_temp, M_temp, argc, argv, &confg3]() {
    cedr_cmplx_flt_type* A_api = (cedr_cmplx_flt_type*) malloc(N_temp * N_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* B_api = (cedr_cmplx_flt_type*) malloc(N_temp * M_temp * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type* C_api = (cedr_cmplx_flt_type*) malloc(N_temp * M_temp * sizeof(cedr_cmplx_flt_type));

    for (size_t i = 0; i < N_temp * N_temp; i++) {
      A_api[i].re = (cedr_re_flt_type) result3[i];
      A_api[i].im = (cedr_re_flt_type) result3imag[i];
    }
    for (size_t i = 0; i < N_temp * M_temp; i++) {
      B_api[i].re = (cedr_re_flt_type) S[i];
      B_api[i].im = (cedr_re_flt_type) Si[i];
    }

    CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, N_temp, M_temp, 0,0,0,0, .argc=argc, .argv=argv);
    //CEDR_GEMM_flt(A_api, B_api, C_api, N_temp, N_temp, M_temp, 0,0,0,0, .argc=argc, .argv=argv, .config=confg3, .task_id=0);

    for (size_t i = 0; i < N_temp * M_temp; i++) {
      result4[i] = (double) C_api[i].re;
      result4imag[i] = (double) C_api[i].im;
    }

    free(A_api);
    free(B_api);
    free(C_api);
  }).name("gemm3-CEDR_GEMM_flt-1");
  /******** API Change ********/

  // Now we have to compute the final operation which is matrix subtraction : (Z
  // - result4) ---> Zr  4*64 - 4*64
  double *zres, *zresimag;
  zres = (double *)malloc(N * M * sizeof(double));
  zresimag = (double *)malloc(N * M * sizeof(double));

  auto t_msub = taskflow.emplace([Z, Zi, result4, result4imag, &zres, &zresimag,N_temp]() {
    msub(Z, Zi, result4, result4imag, zres, zresimag);
  }).name("msub-non_kernel");

  t_gemm0.succeed(t_hermitian); //result1
  t_gemm1.succeed(t_hermitian); //result2
  t_alternateinverse.succeed(t_gemm1); //inv1
  t_mmultiply.succeed(t_alternateinverse); //bufferinv1
  t_mmultiply2.succeed(t_mmultiply); //bufferinv2
  t_zip1.succeed(t_mmultiply2); //bufferinv3
  t_alternateinverse2.succeed(t_zip1); //bufferinv4
  t_mmultiply0.succeed(t_alternateinverse2); //bufferinv5
  t_gemm2.succeed(t_gemm0, t_mmultiply0); //result3
  t_gemm3.succeed(t_gemm2); //result4
  t_msub.succeed(t_gemm3); //zres
  #ifdef DAG_PARSE
    taskflow.dump(std::cout);
  #endif
  //CEDR_DAG_EXTRACT(&taskflow, &task_configs);
  //CEDR_RUN_DAG(&taskflow, 1);
  tf::Executor executor;
  struct timespec start, end;
  std::vector<long long> elapsed_times;
  long long elapsed_time;
  int num_measurements = 1;
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  CEDR_DAG_EXTRACT(&taskflow, &task_configs);
  clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  elapsed_time = (end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec);
  printf("Scheduling: %lld ns\n", elapsed_time);
  for(int iii=0;iii<1;iii++){
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    //executor.run(taskflow).wait();
    CEDR_RUN_DAG(&taskflow, 1000);
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
