#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <numeric>

#include "libcedr.h"
#include "libtaskflow.h"
#define PROGPATH "./input/"
#define RAWDATA PROGPATH "rawdata_rda.txt"
#define OUTPUT "SAR_RDA-output.txt"

// Define KERN_ENTER and KERN_EXIT as NOPs
#define KERN_ENTER(KERN_STR)
#define KERN_EXIT(KERN_STR)

#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>

void swap(cedr_re_flt_type *, cedr_re_flt_type *);
void fftshift(cedr_cmplx_flt_type *, double);

void swap(cedr_re_flt_type *v1, cedr_re_flt_type *v2) {
  cedr_re_flt_type tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

void fftshift(cedr_cmplx_flt_type *data, double count) {
  int k = 0;
  int c = (double)floor(count / 2);
  // For odd and for even numbers of element use different algorithm
  if ((int)count % 2 == 0) {
    for (k = 0; k < c; k += 1) {
      swap(&data[k].re, &data[k + c].re);
      swap(&data[k].im, &data[k + c].im);
    }
  } else {
    cedr_cmplx_flt_type tmp;
    tmp.re = data[0].re;
    tmp.im = data[0].im;
    for (k = 0; k < c; k += 1) {
      data[k].re = data[c + k + 2].re;
      data[k].im = data[c + k + 2].im;
      data[c + k + 2].re = data[k + 2].re;
      data[c + k + 2].im = data[k + 2].im;
    }
    data[c].re = tmp.re;
    data[c].im = tmp.im;
  }
}

int main(int argc, char** argv) {
  double c = 3e8;

  int i, j;

  int Nslow;
  int Nfast;
  double v;
  double Xmin;
  double Xmax;
  double Yc;
  double Y0;
  double Tr;
  double Kr;
  double h;
  double lambda;

  double R0;
  double Ka;
  cedr_cmplx_flt_type *fft_out_0;
  cedr_cmplx_flt_type *fft_inp_1;
  cedr_cmplx_flt_type *fft_out_1;
  cedr_cmplx_flt_type *fft_inp_2;
  cedr_cmplx_flt_type *fft_out_2;
  cedr_cmplx_flt_type *fft_inp_3;

  FILE *fp;

  double *ta;

  double Rmin, Rmax;
  double *tr;

  cedr_cmplx_flt_type *g;
  cedr_cmplx_flt_type *g2;
  cedr_cmplx_flt_type *H;

  cedr_cmplx_flt_type *s0;
  cedr_cmplx_flt_type *S1;

  // Azimuth Compression
  double *sac;
  
  // printf("[SAR] Starting execution of the non-kernel thread\n");

  Nslow = 256;
  Nfast = 512;
  v = 150;
  Xmin = 0;
  Xmax = 50;
  Yc = 10000;
  Y0 = 500;
  Tr = 2.5e-6;
  Kr = 2e13;
  h = 5000;
  lambda = 0.0566;

  R0 = sqrt(Yc * Yc + h * h);
  Ka = 2 * v * v / lambda / R0;
  s0 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  fft_out_0 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  fft_inp_1 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  fft_out_1 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  fft_inp_2 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  fft_out_2 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  fft_inp_3 = (cedr_cmplx_flt_type*) malloc(Nslow * Nfast * sizeof(cedr_cmplx_flt_type));
  sac = (double*) malloc(Nslow * Nfast * sizeof(double));

  fp = fopen(RAWDATA, "r");
  int scan_result;
  for (i = 0; i < Nslow * Nfast; i++) {
    scan_result = fscanf(fp, "%f", &s0[i].re);
    scan_result = fscanf(fp, "%f", &s0[i].im);
  }
  fclose(fp);

  ta = (double*) malloc(Nslow * sizeof(double));
  ta[0] = 0;
  for (i = 1; i < Nslow; i++) {
    ta[i] = ta[i - 1] + (Xmax - Xmin) / v / (Nslow - 1);
  }

  Rmin = sqrt((Yc - Y0) * (Yc - Y0) + h * h);
  Rmax = sqrt((Yc + Y0) * (Yc + Y0) + h * h);
  tr = (double*) malloc(Nfast * sizeof(double));
  tr[0] = 0;
  for (i = 1; i < Nfast; i++) {
    tr[i] = tr[i - 1] + (2 * Rmax / c + Tr - 2 * Rmin / c) / (Nfast - 1);
  }
  
  //Set up g and H
  g = (cedr_cmplx_flt_type*) malloc(Nfast * sizeof(cedr_cmplx_flt_type));
  g2 = (cedr_cmplx_flt_type*) malloc(Nfast * sizeof(cedr_cmplx_flt_type));
  for (i = 0; i < Nfast; i += 1) {
    if (tr[i] > -Tr / 2 && tr[i] < Tr / 2) {
      g[i].re = cos(M_PI * Kr * tr[i] * tr[i]);
      g[i].im = -sin(M_PI * Kr * tr[i] * tr[i]);
    } else {
      g[i].re = 0;
      g[i].im = 0;
    }
  }

  H = (cedr_cmplx_flt_type*) malloc(Nslow * sizeof(cedr_cmplx_flt_type));
  for (i = 0; i < Nslow; i += 1) {
    if (ta[i] > -Tr / 2 * (Xmax - Xmin) / v / (2 * Rmax / c + Tr - 2 * Rmin / c) &&
        ta[i] < Tr / 2 * (Xmax - Xmin) / v / (2 * Rmax / c + Tr - 2 * Rmin / c)) {
      H[i].re = cos(M_PI * Ka * ta[i] * ta[i]);
      H[i].im = sin(M_PI * Ka * ta[i] * ta[i]);
    } else {
      H[i].re = 0;
      H[i].im = 0;
    }
  }

  bool forwardTrans = true;

  tf::Taskflow taskflow;
  std::map<std::string, cedr_task_config_t*> task_configs;
  
  taskflow.name("SAR");

  task_configs["task0-CEDR_FFT_flt-1"] = NULL;
  cedr_task_config_t *&conf0 = task_configs["task0-CEDR_FFT_flt-1"];
  auto t_fft0 = taskflow.emplace([g, &g2, Nfast, forwardTrans, argc, argv, &conf0]() {
    CEDR_FFT_flt(g, g2, Nfast, forwardTrans, .argc=argc, .argv=argv);
    //CEDR_FFT_flt(g, g2, Nfast, forwardTrans, .argc=argc, .argv=argv, .config=conf0, .task_id=0);
  }).name("task0-CEDR_FFT_flt-1");
  
  int first_fft1 = 0;
  int last_fft1 = Nslow;
  task_configs["phase1-CEDR_FFT_flt-" + std::to_string(Nslow)] = NULL;
  cedr_task_config_t *&confp1 = task_configs["phase1-CEDR_FFT_flt-" + std::to_string(Nslow)];
  auto fft_phase1 = taskflow.for_each_index(std::ref(first_fft1), std::ref(last_fft1), 1, [s0, &fft_out_0, Nfast, argc, argv, &confp1](int i){
    CEDR_FFT_flt((s0+i*Nfast), (fft_out_0+i*Nfast), Nfast, true, .argc=argc, .argv=argv);
    //CEDR_FFT_flt((s0+i*Nfast), (fft_out_0+i*Nfast), Nfast, true, .argc=argc, .argv=argv, .config=confp1, .task_id=i);
  }).name("phase1-CEDR_FFT_flt-" + std::to_string(Nslow));

  int first_shift1 = 0;
  int last_shift1 = Nslow;
  auto shift_phase1 = taskflow.for_each_index(std::ref(first_shift1), std::ref(last_shift1), 1, [&fft_out_0, Nfast](int i){
    fftshift((fft_out_0+i*Nfast), Nfast);
  }).name("fftshift1-non_kernel");

  int first_zip1 = 0;
  int last_zip1 = Nslow;
  task_configs["phase1-CEDR_ZIP_flt-" + std::to_string(Nslow)] = NULL;
  cedr_task_config_t *&confp_zip1 = task_configs["phase1-CEDR_ZIP_flt-" + std::to_string(Nslow)];
  auto zip_phase1 = taskflow.for_each_index(std::ref(first_zip1), std::ref(last_zip1), 1, [fft_out_0, g2, &fft_inp_1, Nfast, argc, argv, &confp_zip1](int i){
    CEDR_ZIP_flt(&(fft_out_0[i*Nfast]), g2, &(fft_inp_1[i*Nfast]), Nfast, ZIP_MULT, .argc=argc, .argv=argv);
    //CEDR_ZIP_flt(&(fft_out_0[i*Nfast]), g2, &(fft_inp_1[i*Nfast]), Nfast, ZIP_MULT, .argc=argc, .argv=argv, .config=confp_zip1, .task_id=i);
/*
    for (j = 0; j < Nfast; j += 1) {
      fft_inp_1[j].re = (fft_out_0+i*Nfast)[j].re * g2[j].re - (fft_out_0+i*Nfast)[j].im * g2[j].im;
      fft_inp_1[j].im = (fft_out_0+i*Nfast)[j].im * g2[j].re + (fft_out_0+i*Nfast)[j].re * g2[j].im;
    }
*/
  }).name("phase1-CEDR_ZIP_flt-" + std::to_string(Nslow));

  int first_fft2 = 0;
  int last_fft2 = Nslow;
  task_configs["phase2-CEDR_FFT_flt-" + std::to_string(Nslow)] = NULL;
  cedr_task_config_t *&confp2 = task_configs["phase2-CEDR_FFT_flt-" + std::to_string(Nslow)];
  auto fft_phase2 = taskflow.for_each_index(std::ref(first_fft2), std::ref(last_fft2), 1, [fft_inp_1, &fft_out_1, Nfast, argc, argv, &confp2](int i){
    CEDR_FFT_flt(&(fft_inp_1[i*Nfast]), &(fft_out_1[i*Nfast]), Nfast, false, .argc=argc, .argv=argv);
    //CEDR_FFT_flt(&(fft_inp_1[i*Nfast]), &(fft_out_1[i*Nfast]), Nfast, false, .argc=argc, .argv=argv, .config=confp2, .task_id=i);
  }).name("phase2-CEDR_FFT_flt-" + std::to_string(Nslow));

  int first_trans1 = 0;
  int last_trans1 = Nslow;
  auto trans_phase1 = taskflow.for_each_index(std::ref(first_trans1), std::ref(last_trans1), 1, [fft_out_1, &fft_inp_3, Nslow, Nfast](int i){
    for (int j = 0; j < Nfast; j += 1) {
      fft_inp_3[j * Nslow + i].re = fft_out_1[i*Nfast + j].re;
      fft_inp_3[j * Nslow + i].im = fft_out_1[i*Nfast + j].im;
    }
  }).name("trans-non_kernel");

  // Azimuth FFT
  S1 = (cedr_cmplx_flt_type*) malloc(Nfast * Nslow * sizeof(cedr_cmplx_flt_type));
  int first_fft3 = 0;
  int last_fft3 = Nfast;
  task_configs["phase3-CEDR_FFT_flt-" + std::to_string(Nfast)] = NULL;
  cedr_task_config_t *&confp3 = task_configs["phase3-CEDR_FFT_flt-" + std::to_string(Nfast)];
  auto fft_phase3 = taskflow.for_each_index(std::ref(first_fft3), std::ref(last_fft3), 1, [fft_inp_3, &S1, Nslow, argc, argv, &confp3](int i){
    CEDR_FFT_flt(&(fft_inp_3[i*Nslow]), &(S1[i*Nslow]), Nslow, true, .argc=argc, .argv=argv);
    //CEDR_FFT_flt(&(fft_inp_3[i*Nslow]), &(S1[i*Nslow]), Nslow, true, .argc=argc, .argv=argv, .config=confp3, .task_id=i);
  }).name("phase3-CEDR_FFT_flt-" + std::to_string(Nfast));

  int first_shift2 = 0;
  int last_shift2 = Nfast;
  auto shift_phase2 = taskflow.for_each_index(std::ref(first_shift2), std::ref(last_shift2), 1, [&S1, Nslow](int i){
    fftshift((S1+i*Nslow), Nslow);
  }).name("fftshift2-non_kernel");

  // Azimuth Compression
  int first_zip2 = 0;
  int last_zip2 = Nfast;
  task_configs["phase2-CEDR_ZIP_flt-" + std::to_string(Nfast)] = NULL;
  cedr_task_config_t *&confp_zip2 = task_configs["phase2-CEDR_ZIP_flt-" + std::to_string(Nfast)];
  auto zip_phase2 = taskflow.for_each_index(std::ref(first_zip2), std::ref(last_zip2), 1, [S1, H, &fft_inp_2, Nslow, argc, argv, &confp_zip2](int i){
    CEDR_ZIP_flt((S1+i*Nslow), H, (fft_inp_2+i*Nslow), Nslow, ZIP_MULT, .argc=argc, .argv=argv);
    //CEDR_ZIP_flt((S1+i*Nslow), H, (fft_inp_2+i*Nslow), Nslow, ZIP_MULT, .argc=argc, .argv=argv, .config=confp_zip2, .task_id=i);
/*
    for (j = 0; j < Nfast; j += 1) {
      fft_inp_2[j].re = (S1+i*Nslow)[j].re * H[j].re - (S1+i*Nslow)[j].im * H[j].im;
      fft_inp_2[j].im = (S1+i*Nslow)[j].im * H[j].re + (S1+i*Nslow)[j].re * H[j].im;
    }
*/
  }).name("phase2-CEDR_ZIP_flt-" + std::to_string(Nfast));

  int first_fft4 = 0;
  int last_fft4 = Nfast;
  task_configs["phase4-CEDR_FFT_flt-" + std::to_string(Nfast)] = NULL;
  cedr_task_config_t *&confp4 = task_configs["phase4-CEDR_FFT_flt-" + std::to_string(Nfast)];
  auto fft_phase4 = taskflow.for_each_index(std::ref(first_fft4), std::ref(last_fft4), 1, [fft_inp_2, &fft_out_2, Nslow, argc, argv, &confp4](int i){
    CEDR_FFT_flt((fft_inp_2+i*Nslow), (fft_out_2+i*Nslow), Nslow, false, .argc=argc, .argv=argv);
    //CEDR_FFT_flt((fft_inp_2+i*Nslow), (fft_out_2+i*Nslow), Nslow, false, .argc=argc, .argv=argv, .config=confp4, .task_id=i);
  }).name("phase4-CEDR_FFT_flt-" + std::to_string(Nfast));

  int first_shift3 = 0;
  int last_shift3 = Nfast;
  auto shift_phase3 = taskflow.for_each_index(std::ref(first_shift3), std::ref(last_shift3), 1, [&fft_out_2, Nslow](int i){
    fftshift((fft_out_2+i*Nslow), Nslow);
  }).name("fftshift3-non_kernel");

  int first_sqrt1 = 0;
  int last_sqrt1 = Nfast;
  auto sqrt_phase1 = taskflow.for_each_index(std::ref(first_sqrt1), std::ref(last_sqrt1), 1, [fft_out_2, &sac, Nslow, Nfast](int i){
    for (int j = 0; j < Nslow; j++) {
      sac[i+j*Nfast] = sqrt(fft_out_2[i*Nslow+j].re * fft_out_2[i*Nslow+j].re + fft_out_2[i*Nslow+j].im * fft_out_2[i*Nslow+j].im);
    }
  }).name("sqrt-non_kernel");

  fft_phase1.succeed(t_fft0);
  fft_phase1.precede(shift_phase1);
  shift_phase1.precede(zip_phase1);
  zip_phase1.precede(fft_phase2);
  fft_phase2.precede(trans_phase1);
  trans_phase1.precede(fft_phase3);
  fft_phase3.precede(shift_phase2);
  shift_phase2.precede(zip_phase2);
  zip_phase2.precede(fft_phase4);
  fft_phase4.precede(shift_phase3);
  shift_phase3.precede(sqrt_phase1);
  #if defined(DAG_PARSE)
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
  printf("SAR Average: %.0f ns\n", average);
  printf("SAR Median: %.0f ns\n", median);
  printf("SAR Min: %lld ns\n", min);
  printf("SAR Max: %lld ns\n", max);
  printf("%.0f\n%.0f\n%lld\n%lld\n", average, median, min, max);

  fp = fopen("./output/SAR_output.txt", "w");
  if (fp != NULL) {
    for (i = 0; i < Nslow; i++) {
      for (j = 0; j < Nfast; j++) {
        fprintf(fp, "%lf ", sac[j + i * Nfast]);
      }
      fprintf(fp, "\n");
      fflush(fp);
    }
    fclose(fp);
  }
  //printf("[SAR] Execution is complete...\n");
  return 0;
}
