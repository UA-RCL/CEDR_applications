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
#include <omp.h>
#define PROGPATH "./input/"
#define RAWDATA PROGPATH "rawdata_rda.txt"
#define OUTPUT "SAR_RDA-output.txt"

// Define KERN_ENTER and KERN_EXIT as NOPs
#define KERN_ENTER(KERN_STR)
#define KERN_EXIT(KERN_STR)

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
  for (i = 0; i < Nslow * Nfast; i++) {
    fscanf(fp, "%f", &s0[i].re);
    fscanf(fp, "%f", &s0[i].im);
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

  size_t fast = 512;
  size_t slow = 256;
  bool forwardTrans = true;

  struct timespec start, end;
  std::vector<long long> elapsed_times;
  long long elapsed_time;
  int num_measurements = 100;
  for(int iii=0;iii<100;iii++){
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  CEDR_FFT_flt(g, g2, fast, forwardTrans, .argc=argc, .argv=argv);
  // printf("[SAR] Kernel execution is complete!\n");
  cedr_task_config_t config;
  //CEDR_FFT_flt_sched(Nslow, &config);
	#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nslow; i++) {
    CEDR_FFT_flt((s0+i*Nfast), (fft_out_0+i*Nfast), fast, true, .argc=argc, .argv=argv);//, .config=&config, .task_id=i);
  }
  //#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nslow; i++) {
    fftshift((fft_out_0+i*Nfast), Nfast);
  }
  //CEDR_ZIP_flt_sched(Nslow, &config);
  #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nslow; i++) {
    CEDR_ZIP_flt(&(fft_out_0[i*Nfast]), g2, &(fft_inp_1[i*Nfast]), Nfast, ZIP_MULT, .argc=argc, .argv=argv);//, .config=&config, .task_id=i);
/*
    for (j = 0; j < Nfast; j += 1) {
      fft_inp_1[j].re = (fft_out_0+i*Nfast)[j].re * g2[j].re - (fft_out_0+i*Nfast)[j].im * g2[j].im;
      fft_inp_1[j].im = (fft_out_0+i*Nfast)[j].im * g2[j].re + (fft_out_0+i*Nfast)[j].re * g2[j].im;
    }
*/
  }
  //CEDR_FFT_flt_sched(Nslow, &config);
  #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nslow; i++) {
    CEDR_FFT_flt(&(fft_inp_1[i*Nfast]), &(fft_out_1[i*Nfast]), fast, false,argc=argc, .argv=argv);//, .config=&config, .task_id=i);
  }
  
  for (i = 0; i < Nslow; i++) {
    for (j = 0; j < Nfast; j += 1) {
      fft_inp_3[j * Nslow + i].re = fft_out_1[i*Nfast + j].re;
      fft_inp_3[j * Nslow + i].im = fft_out_1[i*Nfast + j].im;
    }
  }

  // Azimuth FFT
  S1 = (cedr_cmplx_flt_type*) malloc(Nfast * Nslow * sizeof(cedr_cmplx_flt_type));
  //CEDR_FFT_flt_sched(Nfast, &config);
  #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nfast; i++) {
    CEDR_FFT_flt(&(fft_inp_3[i*Nslow]), &(S1[i*Nslow]), slow, true,argc=argc, .argv=argv);//, .config=&config, .task_id=i);
  }
  //#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nfast; i++) {
    fftshift((S1+i*Nslow), Nslow);
  }
  // Azimuth Compression
  //CEDR_ZIP_flt_sched(Nfast, &config);
  #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nfast; i++) {
    //KERN_ENTER(make_label("ZIP[multiply][%d][float64][complex]", Nslow));
    CEDR_ZIP_flt((S1+i*Nslow), H, (fft_inp_2+i*Nslow), Nslow, ZIP_MULT, .argc=argc, .argv=argv);//, .config=&config, .task_id=i);
/*
    for (j = 0; j < Nfast; j += 1) {
      fft_inp_2[j].re = (S1+i*Nslow)[j].re * H[j].re - (S1+i*Nslow)[j].im * H[j].im;
      fft_inp_2[j].im = (S1+i*Nslow)[j].im * H[j].re + (S1+i*Nslow)[j].re * H[j].im;
    }
*/
  }
  //CEDR_FFT_flt_sched(Nfast, &config);
  #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nfast; i++) {
    CEDR_FFT_flt((fft_inp_2+i*Nslow), (fft_out_2+i*Nslow), slow, false,argc=argc, .argv=argv);//, .config=&config, .task_id=i);
  }

  //#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)
  for (i = 0; i < Nfast; i++) {
    fftshift((fft_out_2+i*Nslow), Nslow);
  }
  for (i = 0; i < Nfast; i++) {
    for (j = 0; j < Nslow; j++) {
      sac[i+j*Nfast] = sqrt(fft_out_2[i*Nslow+j].re * fft_out_2[i*Nslow+j].re + fft_out_2[i*Nslow+j].im * fft_out_2[i*Nslow+j].im);
    }
  }

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
