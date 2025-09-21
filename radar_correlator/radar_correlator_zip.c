#include <math.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>

#include "libcedr.h"

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
    cedr_cmplx_flt_type *X1, *X2, *corr_freq;
    double *time, *received, *dftMatrix, *indftMatrix, *corr, *gen_wave;

    clock_t begin, end;

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

    X1 = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
    X2 = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
    corr_freq = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
   
    // lol scope
    {
      cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));

      for (size_t i = 0; i < len; i++) {
        fft_inp[i].re = (cedr_re_flt_type) c[2*i];
        fft_inp[i].im = (cedr_re_flt_type) c[2*i+1];
      }

      CEDR_FFT_flt(fft_inp, X1, len, true,.argc=argc, .argv=argv);
      free(fft_inp);
    }
    // lol scope
    {
      cedr_cmplx_flt_type *fft_inp = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));

      for (size_t i = 0; i < len; i++) {
        fft_inp[i].re = (cedr_re_flt_type) d[2*i];
        fft_inp[i].im = (cedr_re_flt_type) d[2*i+1];
      }

      CEDR_FFT_flt(fft_inp, X2, len, true,.argc=argc, .argv=argv);
      free(fft_inp);
    }
    for (size_t i = 0; i < len; i++) {
      X2[i].im = -X2[i].im;
    }
    
    CEDR_ZIP_flt(X1, X2, corr_freq, len, ZIP_MULT,.argc=argc, .argv=argv);

    // lol scope
    {
      cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));

      CEDR_FFT_flt(corr_freq, fft_out, len, false,.argc=argc, .argv=argv);

      for (size_t i = 0; i < len; i++) {
        corr[2*i] = (double) fft_out[i].re;
        corr[2*i+1] = (double) fft_out[i].im;
      }

      free(fft_out);
    }

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
    printf("Radar correlator complete; lag Value is: %lf\n", lag);
       
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