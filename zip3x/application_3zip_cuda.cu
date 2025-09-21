#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "libcedr.h"
#include <time.h>
#include <mutex>
#include <inttypes.h>
#include <cstring>
#include <cuda_runtime.h>

using namespace std;
#define SEC2NANOSEC 1000000000


// Function for launching ZIP on GPU
__global__ void vector_mult(const cedr_cmplx_flt_type* x, const cedr_cmplx_flt_type* y, cedr_cmplx_flt_type* z, int len) {
  int id = threadIdx.x + blockIdx.x * blockDim.x;
  cedr_re_flt_type r1=x[id].re, r2=y[id].re;
  cedr_re_flt_type i1=x[id].im, i2=y[id].im;
  if (id < len) {
    z[id].re = r1*r2 - i1*i2;
    z[id].im = r1*i2 + r2*i1;
  }
}

int main(int argc, char** argv) {
  #if defined(TIME_CAPTURE)
  const char* filename_time = "results_3zip_cuda.csv";  // Replace with your desired filename
  // Open the file in append mode
  FILE* file = fopen(filename_time, "w");
  if (file == NULL) {
    printf("Failed to open the file.\n");
    return 1;
  }
  fprintf(file, "Size\tAllocation Time\tZIP Execution Time\tDeallocation Time\n");
  struct timespec start_time{};
  struct timespec end_time{};
  long long start, end;
  long long total_time;
  #endif
  int M; // Repeated experiments
  if (argc == 2){
    M = atoi(argv[1]);
  } else{
    M = 100;
  }
  cout << "[APP] Launched app main function with M=" << M << "!" << endl;

  // Allocate input and output memories and generate input
  cedr_cmplx_flt_type *ref_A0, *ref_A1;
  cedr_cmplx_flt_type *ref_B0, *ref_B1;
  cedr_cmplx_flt_type *ref_C0, *ref_C1;
  cedr_cmplx_flt_type *ref_D;
  cedr_cmplx_flt_type *dev_A0, *dev_A1;
  cedr_cmplx_flt_type *dev_B0, *dev_B1;
  cedr_cmplx_flt_type *dev_C0, *dev_C1;
  cedr_cmplx_flt_type *dev_D;


  // CEDR_ZIP API call argument setup
  size_t ZIP_SIZE;  // = N;

  for (size_t N = 64; N <= 65536*2; N=N*2){
    cout << "[APP] Running experiments for size:" << N << endl;
    #if defined(TIME_CAPTURE)
    fprintf(file, "%ld\t", N);
    #endif
    ZIP_SIZE = N;

    // GPU stuff
    const int threadsPerBlock = 512;
    const int blocksPerGrid = (ZIP_SIZE + threadsPerBlock - 1) / threadsPerBlock;

    #if defined(TIME_CAPTURE)
    clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
    #endif

    ref_A0 = new cedr_cmplx_flt_type[N];
    ref_A1 = new cedr_cmplx_flt_type[N];
    ref_B0 = new cedr_cmplx_flt_type[N];
    ref_B1 = new cedr_cmplx_flt_type[N];
    ref_C0 = new cedr_cmplx_flt_type[N];
    ref_C1 = new cedr_cmplx_flt_type[N];
    ref_D  = new cedr_cmplx_flt_type[N];

    cudaMalloc((void**)&dev_A0, N*sizeof(cedr_cmplx_flt_type));
    cudaMalloc((void**)&dev_A1, N*sizeof(cedr_cmplx_flt_type));
    cudaMalloc((void**)&dev_B0, N*sizeof(cedr_cmplx_flt_type));
    cudaMalloc((void**)&dev_B1, N*sizeof(cedr_cmplx_flt_type));
    cudaMalloc((void**)&dev_C0, N*sizeof(cedr_cmplx_flt_type));
    cudaMalloc((void**)&dev_C1, N*sizeof(cedr_cmplx_flt_type));
    cudaMalloc((void**)&dev_D, N*sizeof(cedr_cmplx_flt_type));

    #if defined(TIME_CAPTURE)
    clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
    end = end_time.tv_nsec + end_time.tv_sec*SEC2NANOSEC;
    start = start_time.tv_nsec + start_time.tv_sec*SEC2NANOSEC;
    fprintf(file, "%lld\t", end-start);
    #endif

    for (int i = 0; i < N; i++) {
      ref_A0[i].re = i;
      ref_A0[i].im = i*2;
      ref_A1[i].re = i*3;
      ref_A1[i].im = i*4;
      ref_B0[i].re = i;
      ref_B0[i].im = i*2;
      ref_B1[i].re = i*3;
      ref_B1[i].im = i*4;
      ref_D[i].re  = 0;
      ref_D[i].im  = 0;
    }

    #if defined(TIME_CAPTURE)
      total_time = 0;
    #endif
    for (int ii=0; ii<M; ii++){
      #if defined(TIME_CAPTURE)
        clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
      #endif

      //#pragma omp parallel sections num_threads(2)
      //{
      //  #pragma omp section
      //  {
          cudaMemcpy(dev_A0, ref_A0, N*sizeof(cedr_cmplx_flt_type), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_B0, ref_B0, N*sizeof(cedr_cmplx_flt_type), cudaMemcpyHostToDevice);
          vector_mult<<<blocksPerGrid, threadsPerBlock>>>(dev_A0, dev_B0, dev_C0, ZIP_SIZE);
      //  }
      //  #pragma omp section
      //  {
          cudaMemcpy(dev_A1, ref_A1, N*sizeof(cedr_cmplx_flt_type), cudaMemcpyHostToDevice);
          cudaMemcpy(dev_B1, ref_B1, N*sizeof(cedr_cmplx_flt_type), cudaMemcpyHostToDevice);
          vector_mult<<<blocksPerGrid, threadsPerBlock>>>(dev_A1, dev_B1, dev_C1, ZIP_SIZE);
      //  }
      //}
      // Launch 
      vector_mult<<<blocksPerGrid, threadsPerBlock>>>(dev_C0, dev_C1, dev_D, ZIP_SIZE);
      cudaMemcpy(ref_D, dev_D, N*sizeof(cedr_cmplx_flt_type), cudaMemcpyDeviceToHost);
      //cudaDeviceSynchronize();
      #if defined(TIME_CAPTURE)
        clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
        start = start_time.tv_nsec + start_time.tv_sec * SEC2NANOSEC;
        end = end_time.tv_nsec + end_time.tv_sec * SEC2NANOSEC;
        if (ii >= 0.1 * M){
          total_time += end - start;
        }
      #endif
    }
    #if defined(TIME_CAPTURE)
      fprintf(file, "%lf\t", total_time/(float)(0.9 * M));
    #endif

    
    // Free the allocated memory
    #if defined(TIME_CAPTURE)
      clock_gettime(CLOCK_MONOTONIC_RAW, &start_time);
    #endif
    delete[] ref_A0;
    delete[] ref_A1;
    delete[] ref_B0;
    delete[] ref_B1;
    delete[] ref_C0;
    delete[] ref_C1;
    delete[] ref_D;
    cudaFree(dev_A0);
    cudaFree(dev_A1);
    cudaFree(dev_B0);
    cudaFree(dev_B1);
    cudaFree(dev_C0);
    cudaFree(dev_C1);
    cudaFree(dev_D);
    #if defined(TIME_CAPTURE)
      clock_gettime(CLOCK_MONOTONIC_RAW, &end_time);
      start = start_time.tv_nsec + start_time.tv_sec * SEC2NANOSEC;
      end = end_time.tv_nsec + end_time.tv_sec * SEC2NANOSEC;
      total_time = end - start;
      fprintf(file, "%lld\n", total_time);
    #endif
  }

  cout << "[APP] Exiting app main function!" << endl;

  return 0;
}
