#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>                       
#include "libcedr.h"
#include <omp.h>  
                                          
#define PROGPATH "./input/"               
#define PDPULSE PROGPATH "input_pd_pulse.txt"
#define PDPS PROGPATH "input_pd_ps.txt"   
#define OUTPUT PROGPATH "output_pd_f.txt" 
                                          
#define KERN_ENTER(str)
#define KERN_EXIT(str) 
 
int ARGC;
char** ARGV;      
/* Function Declarations */               
void xcorr(double *, double *, size_t, double *);
void swap(double *, double *);  
void fftshift(double *, double);
                       
/* Function Definitions */                
void xcorr(double *x, double *y, size_t n_samp, double *corr) {
  size_t len;    
  len = 2 * n_samp;  
  double *c = malloc(2 * len * sizeof(double));
  double *d = malloc(2 * len * sizeof(double));
                                          
  size_t x_count = 0;
  size_t y_count = 0;

  double *z = malloc(2 * (n_samp) * sizeof(double));    
  for (size_t i = 0; i < 2 * (n_samp); i += 2) { 
    z[i] = 0;    
    z[i + 1] = 0;
  }                  
  for (size_t i = 0; i < 2 * (n_samp - 1); i += 2) {                  
    c[i] = 0;    
    c[i + 1] = 0;
  }                                       
  memcpy(c + 2 * (n_samp - 1), x, 2 * n_samp * sizeof(double));   
  c[2 * len - 2] = 0;
  c[2 * len - 1] = 0;
  memcpy(d, y, 2 * n_samp * sizeof(double));   
  memcpy(d + 2 * n_samp, z, 2 * (n_samp) * sizeof(double)); 
  cedr_cmplx_flt_type *X1 = malloc(len * sizeof(cedr_cmplx_flt_type));        
  cedr_cmplx_flt_type *X2 = malloc(len * sizeof(cedr_cmplx_flt_type));        
  cedr_cmplx_flt_type *corr_freq = malloc(len * sizeof(cedr_cmplx_flt_type));
 
  // gsl_fft(c, X1, len) and gsl_fft(d, X2, len) (independent);   
  // lol scope
  {

    cedr_cmplx_flt_type *fft_inp_x1 = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
    cedr_cmplx_flt_type *fft_inp_x2 = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
    bool is_fwd = true;

    for (size_t i = 0; i < len; i++) {
      fft_inp_x1[i].re = (cedr_re_flt_type) c[2*i];
      fft_inp_x1[i].im = (cedr_re_flt_type) c[2*i+1];
      fft_inp_x2[i].re = (cedr_re_flt_type) d[2*i];
      fft_inp_x2[i].im = (cedr_re_flt_type) d[2*i+1];
    }



		#pragma omp parallel sections num_threads(2)
		{
				#pragma omp section
				{
				    CEDR_FFT_flt(fft_inp_x1, X1, len, is_fwd,.argc=ARGC, .argv=ARGV);
				}
				
				#pragma omp section
				{
				    CEDR_FFT_flt(fft_inp_x2, X2, len, is_fwd,.argc=ARGC, .argv=ARGV);
				}
		}


    free(fft_inp_x1);
    free(fft_inp_x2);
  }
  for (size_t i = 0; i < len; i++) {
    X2[i].im = -X2[i].im;
  }
  
  free(c);
  free(d);
  free(z);
  KERN_ENTER(make_label("ZIP[multiply][complex][float64][%d]", len));
  //CEDR_ZIP_flt(X1, X2, corr_freq, len, ZIP_MULT);
  CEDR_ZIP_flt(X1, X2, corr_freq, len, ZIP_MULT,.argc=ARGC, .argv=ARGV);  
  KERN_EXIT(make_label("ZIP[multiply][complex][float64][%d]", len));
  free(X1);
  free(X2);

  // gsl_ifft(corr_freq, corr, len);
  // lol scope
  {
    cedr_cmplx_flt_type *fft_out = (cedr_cmplx_flt_type*) malloc(len * sizeof(cedr_cmplx_flt_type));
    bool is_fwd = false;

    CEDR_FFT_flt(corr_freq, fft_out, len, is_fwd,.argc=ARGC,.argv=ARGV);

    for (size_t i = 0; i < len; i++) {
      corr[2*i] = (double) fft_out[i].re;
      corr[2*i+1] = (double) fft_out[i].im;
    }

    free(fft_out);
  }

  free(corr_freq);
}

void swap(double *v1, double *v2) {
  double tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

void fftshift(double *data, double count) {
  int k = 0;
  int c = (double)floor((float)count / 2);
  // For odd and for even numbers of element use different algorithm
  if ((int)count % 2 == 0) {
    for (k = 0; k < c; k++)
      swap(&data[k], &data[k + c]);
  } else {
    double tmp = data[0];
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

  double *mf = malloc((2 * n_samples) * m * 2 *
                      sizeof(double)); // build a 2D array for the output of the
                                       // matched filter
  double *p =
      malloc(2 * n_samples * sizeof(double)); // array for pulse with noise
  double *pulse =
      malloc(2 * n_samples * sizeof(double)); // array for the original pulse
  double *corr =
      malloc(2 * (2 * n_samples) *
             sizeof(double)); // array for the output of matched filter

  // Read the original pulse
  FILE *fp;
  fp = fopen(PDPULSE, "r");
  for (i = 0; i < 2 * n_samples; i++) {
    fscanf(fp, "%lf", &pulse[i]);
  }
  fclose(fp);

  // Run the input samples through the matched filter
  fp = fopen(PDPS, "r"); // read the multiple pulses with noise and delay
  for (k = 0; k < m; k++) {
    for (j = 0; j < 2 * n_samples; j++) {
      fscanf(fp, "%lf", &p[j]);
    }

    /* matched filter */
    xcorr(p, pulse, n_samples, corr);

    /* put the output into a new 2D array */
    for (n = 0; n < 2 * (2 * n_samples); n += 2) {
      mf[n / 2 + (2 * k) * (2 * n_samples)] = corr[n];
      mf[n / 2 + (2 * k + 1) * (2 * n_samples)] = corr[n + 1];
    }
  }
  fclose(fp);
  free(p);
  free(pulse);
  free(corr);

  double **q = malloc(2*n_samples* sizeof(double));
  double **r = malloc(2*n_samples* sizeof(double));
  double *l = malloc(2 * m * sizeof(double));
  double *f = malloc(m * (2 * n_samples) * sizeof(double));
  double max = 0, a, b;
  

  /* FFT */
  // gsl_fft(l,q,m);
  // lol scope
  {
  int num_ffts = 2*n_samples;
  bool is_fwd = true;
  cedr_cmplx_flt_type **fft_inp = (cedr_cmplx_flt_type**) malloc(num_ffts * sizeof(cedr_cmplx_flt_type));
  cedr_cmplx_flt_type **fft_out = (cedr_cmplx_flt_type**) malloc(num_ffts * sizeof(cedr_cmplx_flt_type));
  
  for (x = 0; x < num_ffts; x++) {
    
    for (o = 0; o < 2 * m; o++) {
      l[o] = mf[x + o * num_ffts];
    }
    
      fft_inp[x] = (cedr_cmplx_flt_type*) malloc(m * sizeof(cedr_cmplx_flt_type));
      fft_out[x] = (cedr_cmplx_flt_type*) malloc(m * sizeof(cedr_cmplx_flt_type));

      for (size_t i = 0; i < m; i++) {
        fft_inp[x][i].re = (cedr_re_flt_type) l[2*i];
        fft_inp[x][i].im = (cedr_re_flt_type) l[2*i+1];
      }
  }
  // exit for loop to do non-blocking api calls over array of fft_inp
  // CEDR_FFT_flt(fft_inp, fft_out, m, true /* is_forward_transform? */);
	#pragma omp parallel for schedule(SCHEDULER) num_threads(THREADS) 
  for (size_t i = 0; i < num_ffts; i++) {
    CEDR_FFT_flt(fft_inp[i], fft_out[i], m, is_fwd,.argc=ARGC,.argv=ARGV);
    //CEDR_FFT_flt_nb(&fft_inp[i], &fft_out[i], &m, &is_fwd, &barrier);
  }



  for(x = 0; x < num_ffts; x++){
      q[x] = malloc(2 * m *sizeof(double));
      for (size_t i = 0; i < m; i++) {
        q[x][2*i] = (double) fft_out[x][i].re;
        q[x][2*i+1] = (double) fft_out[x][i].im;
      }
  }

  free(fft_inp);
  free(fft_out);
    
  for(x = 0; x < num_ffts; x++){
    r[x] = malloc(m *sizeof(double));
    for (y = 0; y < 2 * m; y += 2) {
      r[x][y / 2] = sqrt(
          q[x][y] * q[x][y] +
          q[x][y + 1] * q[x][y + 1]); // calculate the absolute value of the output
    }
    fftshift(r[x], m);

    for (z = 0; z < m; z++) {
      f[x + z * num_ffts] = r[x][z]; // put the elements of output into corresponding location of the 2D array
      if (r[x][z] > max) {
        max = r[x][z];
        a = z + 1;
        b = x + 1;
      }
    }
  }
  }
  
  free(mf);
  free(q);
  free(r);
  free(l);

  double rg, dp;
  rg = (b - n_samples) / (n_samples - 1) * PRI;
  dp = (a - (m + 1) / 2) / (m - 1) / PRI;

  fp = fopen("./output/pulse_doppler_output.txt", "w");
  if (fp != NULL) {
    fprintf(fp, "Doppler shift = %lf, time delay = %lf\n", dp, rg);
    fclose(fp);
  }
  printf("Doppler shift = %lf, time delay = %lf\n", dp, rg);
  printf("[Pulse Doppler] Execution is complete...\n");
}
