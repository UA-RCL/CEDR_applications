#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>
#include <complex>
#include "libcedr.h"
#include <time.h>
#include <mutex>
#include <inttypes.h>
#include <stdlib.h>
#include <cstring>

#include <ctime>
#include <numeric>

using namespace std;
#define SEC2NANOSEC 1000000000

int main(int argc, char** argv) {
  int M; // Repeated experiments
  if (argc == 2){
    M = atoi(argv[1]);
  } else{
    M = 10000;
  }
  #ifdef DEBUG
    cout << "[APP] Launched app main function with M=" << M << "!" << endl;
  #endif

  // Allocate input and output memories and generate input
  cedr_cmplx_flt_type *ref_A0, *ref_A1;
  cedr_cmplx_flt_type *ref_B0, *ref_B1;
  cedr_cmplx_flt_type *ref_C0, *ref_C1;
  cedr_cmplx_flt_type *ref_D;

  // CEDR_ZIP API call argument setup
  size_t ZIP_SIZE;
  int N = 1024*4;
  #ifdef DEBUG
    cout << "[APP] Running experiments for size:" << N << endl;
  #endif
  ZIP_SIZE = N;
  ref_A0 = new cedr_cmplx_flt_type[N];
  ref_A1 = new cedr_cmplx_flt_type[N];
  ref_B0 = new cedr_cmplx_flt_type[N];
  ref_B1 = new cedr_cmplx_flt_type[N];
  ref_C0 = new cedr_cmplx_flt_type[N];
  ref_C1 = new cedr_cmplx_flt_type[N];
  ref_D  = new cedr_cmplx_flt_type[N];

  #ifdef DEBUG
  cout << "[APP] Input values:" << endl;
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
    #ifdef DEBUG
    cout << ref_input1[i].re << " + " << ref_input1[i].im << "i; \t";
    cout << ref_input2[i].re << " + " << ref_input2[i].im << "i; \t";
    #endif
  }

  #if defined(TIME_CAPTURE)
    struct timespec start, end;
    std::vector<long long> elapsed_times;
    long long elapsed_time;
    int num_measurements = 1;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  #endif

  for(int iii=0;iii<M;iii++){
    #if defined(TIME_CAPTURE)
      clock_gettime(CLOCK_MONOTONIC_RAW, &start);
    #endif
    /* Create DAG Nodes*/
    #pragma omp parallel sections num_threads(2)
    {
      #pragma omp section
      {
        CEDR_ZIP_flt(ref_A0, ref_B0, ref_C0, ZIP_SIZE, ZIP_MULT, .argc=argc, .argv=argv);
      }
      #pragma omp section
      {
        CEDR_ZIP_flt(ref_A1, ref_B1, ref_C1, ZIP_SIZE, ZIP_MULT, .argc=argc, .argv=argv);
      }
    }
    CEDR_ZIP_flt(ref_C0, ref_C1, ref_D,  ZIP_SIZE, ZIP_MULT, .argc=argc, .argv=argv);

    #if defined(TIME_CAPTURE)
      clock_gettime(CLOCK_MONOTONIC_RAW, &end);
      elapsed_time = (end.tv_sec - start.tv_sec) * 1e9 + (end.tv_nsec - start.tv_nsec);
      elapsed_times.push_back(elapsed_time);
    #endif
  }
  #if defined(TIME_CAPTURE)
    printf("Start:%0.f\nEnd:%0.f\n", start.tv_sec * 1e9 + start.tv_nsec, end.tv_sec * 1e9 + end.tv_nsec);
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
  #endif

  #ifdef DEBUG
    cout << "[APP] Reference FFT values:" << endl;
    for (int i = 0; i < N; i++) {
      cout << ref_output[i].re << " + " << ref_output[i].im << "i; \t";
    }
  #endif
  
  // Free the allocated memory
  delete[] ref_A0;
  delete[] ref_A1;
  delete[] ref_B0;
  delete[] ref_B1;
  delete[] ref_C0;
  delete[] ref_C1;
  delete[] ref_D;

  #ifdef DEBUG
  cout << "[APP] Exiting app main function!" << endl;
  #endif

  return 0;
}
