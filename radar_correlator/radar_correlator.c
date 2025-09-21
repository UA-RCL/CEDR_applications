#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "libcedr.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <numeric>

int main(int argc, char** argv)
{
    size_t len, x_count, y_count, i, j, k;
    size_t n_samples; 
    size_t dft_size; 

    double lag;
    double T;
    double B;
    double sampling_rate;
    double max_corr;
    double index;

    double *c, *d;
    double *X1, *X2, *corr_freq;
    double *time, *received, *dftMatrix, *indftMatrix, *corr, *gen_wave;

    FILE *fp;

    int row, column, row2, column2;

    n_samples = 256;
    dft_size = 2 * n_samples;
    len = 2 * n_samples;

    T = 0.000512;
    B = 500000;
    sampling_rate = 1000;
    max_corr = 0;
    index = 0;

    time = (double*) malloc(n_samples*sizeof(double));
    received = (double*) malloc(2 * n_samples*sizeof(double));
    dftMatrix = (double*) malloc(2 * dft_size * dft_size * sizeof(double));
    indftMatrix = (double*) malloc(2 * dft_size * dft_size * sizeof(double));
    corr = (double*) malloc(2 * len * sizeof(double));
    gen_wave = (double*) malloc(2 * n_samples * sizeof(double));

    fp = fopen("./input/time_input.txt","r");
    if (fp == NULL) { printf("Unable to open time_input.txt!\n"); return 1; }
    int scan_check;
    for(i = 0; i < n_samples; i++) {
        scan_check = fscanf(fp, "%lf", &time[i]);
    }
    fclose(fp);
    fp = NULL;

    for (i = 0; i < len; i += 2) {
        gen_wave[i] = sin(M_PI * B / T * pow(time[i / 2], 2));
        gen_wave[i + 1] = cos(M_PI * B / T * pow(time[i / 2], 2));
    }

    fp = fopen("./input/received_input.txt","r");
    if (fp == NULL) { printf("Unable to open received_input.txt!\n"); return 1; }
    for(i = 0; i < len; i++) {
        scan_check = fscanf(fp,"%lf", &received[i]);
    }
    fclose(fp);
    fp = NULL;

    c = (double*) malloc(2 * len * sizeof(double));
    d = (double*) malloc(2 * len * sizeof(double));

    x_count = 0;
    y_count = 0;
    for (i = 0; i < 2 * len; i += 2) {
        if (i/2 > n_samples - 1) {
            c[i] = gen_wave[x_count];
            c[i + 1] = gen_wave[x_count + 1];
            x_count += 2;
        } else {
            c[i] = 0;
            c[i + 1] = 0;
        }
        if (i > n_samples) {
            d[i] = 0;
            d[i + 1] = 0;
        } else {
            d[i] = received[y_count];
            d[i + 1] = received[y_count + 1];
            y_count += 2;
        }
    }

    X1 = (double*) malloc(2 * len * sizeof(double));
    X2 = (double*) malloc(2 * len * sizeof(double));
    corr_freq = (double*) malloc(2 * len * sizeof(double));

    struct timespec start, end;
    std::vector<long long> elapsed_times;
    long long elapsed_time;
    int num_measurements = 1000;
    for(int iii=0;iii<1000;iii++){
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    #pragma omp parallel sections num_threads(2)
		{
      #pragma omp section
      {
        cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
        cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));

        for (size_t i = 0; i < len; i++) {
          fft_inp[i].re = (cedr_re_flt_type) c[2*i];
          fft_inp[i].im = (cedr_re_flt_type) c[2*i+1];
        }
        bool direction=true;
        CEDR_FFT_flt(fft_inp, fft_out, len, direction, .argc=argc, .argv=argv);
        //CEDR_FFT_flt(fft_inp, fft_out, len, direction, .argc=argc, .argv=argv, .config=confp1, .task_id=0);

        for (size_t i = 0; i < len; i++) {
          X1[2*i] = (double) fft_out[i].re;
          X1[2*i+1] = (double) fft_out[i].im;
        }

        free(fft_inp);
        free(fft_out);
      }
      
      #pragma omp section
      {
        cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
        cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));

        for (size_t i = 0; i < len; i++) {
          fft_inp[i].re = (cedr_re_flt_type) d[2*i];
          fft_inp[i].im = (cedr_re_flt_type) d[2*i+1];
        }
        bool direction=true;
        CEDR_FFT_flt(fft_inp, fft_out, len, direction, .argc=argc, .argv=argv);
        //CEDR_FFT_flt(fft_inp, fft_out, len, direction, .argc=argc, .argv=argv, .config=confp2, .task_id=0);

        for (size_t i = 0; i < len; i++) {
          X2[2*i] = (double) fft_out[i].re;
          X2[2*i+1] = (double) fft_out[i].im;
        }

        free(fft_inp);
        free(fft_out);
      }
		}

    for (size_t i = 0; i < 2 * len; i += 2) {
        corr_freq[i] = (X1[i] * X2[i]) + (X1[i + 1] * X2[i + 1]);
        corr_freq[i + 1] = (X1[i + 1] * X2[i]) - (X1[i] * X2[i + 1]);
    }

    // lol scope
    {
      cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
      cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));

      for (size_t i = 0; i < len; i++) {
        fft_inp[i].re = (cedr_re_flt_type) corr_freq[2*i];
        fft_inp[i].im = (cedr_re_flt_type) corr_freq[2*i+1];
      }
      bool direction=false;

      CEDR_FFT_flt(fft_inp, fft_out, len, direction, .argc=argc, .argv=argv);
      //CEDR_FFT_flt(fft_inp, fft_out, len, direction, .argc=argc, .argv=argv, .config=confp3, .task_id=0);

      for (size_t i = 0; i < len; i++) {
        corr[2*i] = (double) fft_out[i].re;
        corr[2*i+1] = (double) fft_out[i].im;
      }

      free(fft_inp);
      free(fft_out);
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

    for (i = 0; i < 2 * len; i += 2)
    {
        // Only finding maximum of real part of correlation
        // (corr[i] > max_corr) && (max_corr = corr[i], index = i / 2);
        if (corr[i] > max_corr) {
            max_corr = corr[i];
            index = i / 2;
        }
    }

    lag = (n_samples - index) / sampling_rate;
    //printf("Radar correlator complete; lag Value is: %lf\n", lag);
       
    free(time);    
    free(received);
    free(dftMatrix);
    free(indftMatrix);
    free(corr);
    free(gen_wave);
    free(c);
    free(d);
    free(X1);
    free(X2);
    free(corr_freq);

    return 0;
}
