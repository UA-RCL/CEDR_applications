#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>
#include "libcedr.h"
#include <omp.h>                                                    
#define PROGPATH "./input/"               
#define PDPULSE PROGPATH "input_pd_pulse.txt"
#define PDPS PROGPATH "input_pd_ps.txt"   
#define OUTPUT PROGPATH "output_pd_f.txt" 
                                          
#define KERN_ENTER(str)
#define KERN_EXIT(str) 

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <numeric>


int ARGC;
char** ARGV;

/* Function Declarations */               
void general_func(cedr_re_flt_type *, cedr_re_flt_type *, size_t, cedr_re_flt_type **, size_t, size_t);
void swap(cedr_re_flt_type *, cedr_re_flt_type *);  
void fftshift(cedr_re_flt_type *, cedr_re_flt_type);
                       
/* Function Definitions */                
void general_func(cedr_re_flt_type *x, cedr_re_flt_type *y, size_t n_samp, cedr_re_flt_type **q, size_t m, size_t n_samples) {
  int argc = ARGC;
  char** argv = ARGV;
  size_t len;    
  len = 2 * n_samp;  
  cedr_re_flt_type *c = (cedr_re_flt_type*) malloc(m * 2 * len * sizeof(cedr_re_flt_type));
  cedr_re_flt_type *d = (cedr_re_flt_type*) malloc(m * 2 * len * sizeof(cedr_re_flt_type));
                                          
  size_t x_count = 0;
  size_t y_count = 0;

  cedr_re_flt_type *z = (cedr_re_flt_type*) malloc(m * 2 * (n_samp) * sizeof(cedr_re_flt_type));    
  for (size_t i = 0; i < m * 2 * (n_samp); i += 2) { 
    z[i] = 0;    
    z[i + 1] = 0;
  }                  
  for (int k = 0; k < m; k++) {
    for (size_t i = 0; i < 2 * (n_samp - 1); i += 2) {                  
      c[k*(2*len) + i] = 0;    
      c[k*(2*len) + i + 1] = 0;
    }
    memcpy(c + k*(2*len) + 2 * (n_samp - 1), x, 2 * n_samp * sizeof(cedr_re_flt_type));
    c[k*(2*len) + 2 * len - 2] = 0;
    c[k*(2*len) + 2 * len - 1] = 0;
    memcpy(d + k*(2*len), y, 2 * n_samp * sizeof(cedr_re_flt_type));
    memcpy(d + k*(2*len) + 2 * n_samp, z + k*len, 2 * (n_samp) * sizeof(cedr_re_flt_type)); 
  }  
  cedr_re_flt_type *X1 = (cedr_re_flt_type*) malloc(m * 2 * len * sizeof(cedr_re_flt_type));
  cedr_re_flt_type *X2 = (cedr_re_flt_type*) malloc(m * 2 * len * sizeof(cedr_re_flt_type));
  cedr_re_flt_type *corr_freq = (cedr_re_flt_type*) malloc(m * 2 * len * sizeof(cedr_re_flt_type));
  cedr_re_flt_type *corr = (cedr_re_flt_type*) malloc(m * 2 * (2 * n_samples) * sizeof(cedr_re_flt_type)); // array for the output of matched filter
  cedr_re_flt_type *mf = (cedr_re_flt_type*) malloc((2 * n_samples) * m * 2 * sizeof(cedr_re_flt_type)); // build a 2D array for the output of the // matched filter

  int num_ffts = 2*n_samples;


  struct timespec start, end;
  std::vector<long long> elapsed_times;
  long long elapsed_time;
  int num_measurements = 100;
  for(int iii=0;iii<100;iii++){
  clock_gettime(CLOCK_MONOTONIC_RAW, &start);

  int first_p1 = 0;
  int last_p1 = m;
  //auto phase1 = taskflow.for_each_index(std::ref(first_p1), std::ref(last_p1), 1, [c, &X1, len, argc, argv](int k){
  #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)    
  for(int k=first_p1; k<last_p1; k++){
    cedr_cmplx_flt_type *fft_inp_x1 = (cedr_cmplx_flt_type*) &c[k*2*len];
    cedr_cmplx_flt_type *fft_out_x1 = (cedr_cmplx_flt_type*) &X1[k*2*len];
    bool is_fwd = true;
    CEDR_FFT_flt(fft_inp_x1, fft_out_x1, len, is_fwd, .argc=argc, .argv=argv);
  //}).name("phase1");
   }

  int first_p2 = 0;
  int last_p2 = m;
  //auto phase2 = taskflow.for_each_index(std::ref(first_p2), std::ref(last_p2), 1, [d, &X2, len, argc, argv](int k){
   #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)    
   for(int k=first_p2; k<last_p2; k++){
    cedr_cmplx_flt_type *fft_inp_x2 = (cedr_cmplx_flt_type*) &d[k*2*len];
    cedr_cmplx_flt_type *fft_out_x2 = (cedr_cmplx_flt_type*) &X2[k*2*len];
    bool is_fwd = true;
    CEDR_FFT_flt(fft_inp_x2, fft_out_x2, len, is_fwd, .argc=argc, .argv=argv);
//  }).name("phase2");
	  }

  int first_ct = 0;
  int last_ct = m;
  //auto corr_task = taskflow.for_each_index(std::ref(first_ct), std::ref(last_ct), 1, [X1, X2, &corr_freq, len](int k){
   //#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)    
   for(int k=first_ct; k<last_ct; k++){
    for (size_t i = 0; i < 2 * len; i += 2) {
      corr_freq[k*2*len + i] = (X1[k*2*len + i] * X2[k*2*len + i]) + (X1[k*2*len + i + 1] * X2[k*2*len + i + 1]);
      corr_freq[k*2*len + i + 1] = (X1[k*2*len + i + 1] * X2[k*2*len + i]) - (X1[k*2*len + i] * X2[k*2*len + i + 1]);
    }
 // }).name("zips");
	  }

  int first_p3 = 0;
  int last_p3 = m;
 // auto phase3 = taskflow.for_each_index(std::ref(first_p3), std::ref(last_p3), 1, [corr_freq, &corr, len, argc, argv](int k){
   #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)    
   for(int k=first_p3; k<last_p3; k++){
    cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) &corr_freq[k*2*len];
    cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) &corr[k*2*len];
    bool is_fwd = false;
    CEDR_FFT_flt(fft_inp, fft_out, len, is_fwd, .argc=argc, .argv=argv);
  //}).name("phase3");
	  }

  int first_a2D = 0;
  int last_a2D = m;
  //auto array2D = taskflow.for_each_index(std::ref(first_a2D), std::ref(last_a2D), 1, [corr, &mf, n_samples](int k){
   //#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)    
   for(int k=first_a2D; k<last_a2D; k++){
    // put the output into a new 2D array
    for(int n = 0; n < 2 * (2 * n_samples); n += 2) {
      mf[n / 2 + (2 * k) * (2 * n_samples)] = corr[k*2*2*n_samples + n];
      mf[n / 2 + (2 * k + 1) * (2 * n_samples)] = corr[k*2*2*n_samples + n + 1];
    }
  //}).name("array2D");
		}

  int first_p4 = 0;
  int last_p4 = num_ffts;
//  auto phase4 = taskflow.for_each_index(std::ref(first_p4), std::ref(last_p4), 1, [mf, &q, m, num_ffts, argc, argv](int k){
   #pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS)    
   for(int k=first_p4; k<last_p4; k++){
    bool is_fwd = true;
    cedr_re_flt_type *l = (cedr_re_flt_type*) malloc(2 * m * sizeof(cedr_re_flt_type));
    for (size_t o = 0; o < 2 * m; o++) {
      l[o] = mf[k + o * num_ffts];
    }
    cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) l;
    cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) q[k];
    CEDR_FFT_flt(fft_inp, fft_out, m, is_fwd, .argc=argc, .argv=argv);
    free(l);
 // }).name("phase4");
    }

    cedr_re_flt_type max = 0, a, b;
    double PRI = 1.27e-4;
    cedr_re_flt_type **r = (cedr_re_flt_type**) malloc(2*n_samples* sizeof(cedr_re_flt_type*));
    cedr_re_flt_type *f  = (cedr_re_flt_type*)  malloc(m * (2 * n_samples) * sizeof(cedr_re_flt_type));
    for(int x = 0; x < num_ffts; x++){
      r[x] = (cedr_re_flt_type*) malloc(m *sizeof(cedr_re_flt_type));
      for (int y = 0; y < 2 * m; y += 2) {
        r[x][y / 2] = sqrt(
            q[x][y] * q[x][y] +
            q[x][y + 1] * q[x][y + 1]); // calculate the absolute value of the output
      }
      fftshift(r[x], m);
      for (int z = 0; z < m; z++) {
        f[x + z * num_ffts] = r[x][z]; // put the elements of output into corresponding location of the 2D array
        if (r[x][z] > max) {
          max = r[x][z];
          a = z + 1;
          b = x + 1;
        }
      }
    }

    cedr_re_flt_type rg, dp;
    rg = (b - n_samples) / (n_samples - 1) * PRI;
    dp = (a - (m + 1) / 2) / (m - 1) / PRI;

    /*FILE *fp;
    fp = fopen("./output/pulse_doppler_output.txt", "w");
    if (fp != NULL) {
      fprintf(fp, "Doppler shift = %lf, time delay = %lf\n", dp, rg);
      fclose(fp);
    }
    printf("Doppler shift = %lf, time delay = %lf\n", dp, rg);
    printf("[Pulse Doppler] Execution is complete...\n");*/
    free(f);
    for(int x = 0; x < num_ffts; x++) free(r[x]);
    free(r);

    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
    elapsed_time = (end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec);
    elapsed_times.push_back(elapsed_time);
  }
  auto min_it = std::min_element(elapsed_times.begin(), elapsed_times.end());
  auto max_it = std::max_element(elapsed_times.begin(), elapsed_times.end());
  long long min = *min_it;
  long long max = *max_it;
  size_t min_index = std::distance(elapsed_times.begin(), min_it);
  size_t max_index = std::distance(elapsed_times.begin(), max_it);
  double average = std::accumulate(elapsed_times.begin(), elapsed_times.end(), 0.0) / num_measurements;
  std::sort(elapsed_times.begin(), elapsed_times.end());
  double median = num_measurements % 2 == 0
                  ? (elapsed_times[num_measurements / 2 - 1] + elapsed_times[num_measurements / 2]) / 2.0
                  : elapsed_times[num_measurements / 2];
  printf("Average: %.0f ns\n", average);
  printf("Median: %.0f ns\n", median);
  printf("Min(%ld): %lld ns\n", min_index, min);
  printf("Max(%ld): %lld ns\n", max_index, max);
  printf("%.0f\n%.0f\n%lld\n%lld\n", average, median, min, max);
  
  free(mf);
  free(corr);

  free(c);
  free(d);
  free(z);
  free(X1);
  free(X2);
  free(corr_freq);
}

void swap(cedr_re_flt_type *v1, cedr_re_flt_type *v2) {
  cedr_re_flt_type tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

void fftshift(cedr_re_flt_type *data, cedr_re_flt_type count) {
  int k = 0;
  int c = (cedr_re_flt_type)floor((float)count / 2);
  // For odd and for even numbers of element use different algorithm
  if ((int)count % 2 == 0) {
    for (k = 0; k < c; k++)
      swap(&data[k], &data[k + c]);
  } else {
    cedr_re_flt_type tmp = data[0];
    for (k = 0; k < c; k++) {
      data[k] = data[c + k + 1];
      data[c + k + 1] = data[k + 1];
    }
    data[c] = tmp;
  }
}

int main(int argc, char *argv[]) {
  ARGC = argc;
  ARGV = argv;

  size_t m = 128;    // number of pulses
  size_t n_samples = 64; // length of single pulse
  double PRI = 1.27e-4;
  int i, j, k, n, x, y, z, o;

  cedr_re_flt_type *p = (cedr_re_flt_type*) malloc(m * 2 * n_samples * sizeof(cedr_re_flt_type)); // array for pulse with noise
  cedr_re_flt_type *pulse = (cedr_re_flt_type*) malloc(2 * n_samples * sizeof(cedr_re_flt_type)); // array for the original pulse

  // Read the original pulse
  FILE *fp;
  fp = fopen(PDPULSE, "r");
  for (i = 0; i < 2 * n_samples; i++) {
    if(fscanf(fp, "%f", &pulse[i]) < 0){
      pulse[i] = 0.0;
    }
  }
  fclose(fp);

  // Run the input samples through the matched filter
  fp = fopen(PDPS, "r"); // read the multiple pulses with noise and delay
  for (k = 0; k < m * 2 * n_samples; k++) {
      if(fscanf(fp, "%f", &p[k]) < 0){
        p[k] = 0.0;
      }
  }
  fclose(fp);

  cedr_re_flt_type **q = (cedr_re_flt_type**) malloc(2*n_samples* sizeof(cedr_re_flt_type*));
  
  for(size_t i = 0; i < 2*n_samples; i++){
      q[i] = (cedr_re_flt_type*) malloc(2 * m *sizeof(cedr_re_flt_type));
  }

  /* matched filter */
  //xcorr(p, pulse, n_samples, corr);
  general_func(p, pulse, n_samples, q, m, n_samples);

  free(pulse);
  free(p);
  free(q);
}
