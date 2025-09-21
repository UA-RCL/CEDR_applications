#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <string.h>
#include <unistd.h>
#include "common.h"
#include "IFFT_FFT.h"
#include "scrambler_descrambler.h"
#include "CyclicPrefix.h"
#include "Preamble_ST_LG.h"
#include "baseband_lib.h"
#include "interleaver_deintleaver.h"
#include "qpsk_Mod_Demod.h"
#include "datatypeconv.h"
#include "baseband_lib.h"
#include "pilot.h"
#include "fft_hs.h"
#include "txData.h"
#include "rf_interface.h"
#include "rfInf.h"

#include <sys/stat.h>
#include <sys/mman.h>
#include <fcntl.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <errno.h>
#include <limits.h>

#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <numeric>

#include <complex>
#include "libcedr.h"
#include "libtaskflow.h"
#include <omp.h>          
#include <taskflow/taskflow.hpp>
#include <taskflow/algorithm/for_each.hpp>
#include "viterbi.h"
static char         attr[1024];
int                 fft_id;

//typedef double complex cplx;
typedef std::complex<double> cplx;
const std::complex<double> I (0.0, 1.0);

void random_wait_time(int random_wait_time_in_us);

#define SEC2NANOSEC 1000000000

int _txOption;
int hw_fft_busy;

int fft_a53_count = 0;
int fft_acc_count = 0;

unsigned int *base_addr;

void readConfig();

int main(int argc,char** argv) {
    
   char rate = PUNC_RATE_1_2;
   int encoderId;
   int frameLength;
   unsigned char* inbit;
   unsigned char* scram;
   signed char* enc_out;
   signed char* enc_dep_out;
   unsigned char* intl_out;
   double* sig_real;
   double* sig_img;
   float *in_ifft;
   comp_t** pilot_out;

   // cedr_FFT input-output
   cedr_cmplx_flt_type **ifft_input, **ifft_output;
   

   //comp_t cyclic_out[SYM_NUM*TOTAL_LEN];
   comp_t* cyclic_out;
   //comp_t pre_out[SYM_NUM*TOTAL_LEN + PREAMBLE_LEN + 2048];
   comp_t* pre_out;
   int i, j;
   //comp_t txdata[TX_DATA_LEN];
   comp_t* txdata;

   int user_data_len;

   FILE *cfp, *fp;

   struct timespec qpr_start, qpr_end;

   for (i = 0; i < 1; i++) {}

   clock_gettime(CLOCK_MONOTONIC, &qpr_start);

   inbit = (unsigned char*) calloc(1024, sizeof(unsigned char));
   scram = (unsigned char*) calloc(SYM_NUM*1024, sizeof(unsigned char));
   enc_out = (signed char*) calloc(SYM_NUM*OUTPUT_LEN, sizeof(signed char));
   enc_dep_out = (signed char*) calloc(SYM_NUM*OUTPUT_LEN, sizeof(signed char));
   intl_out = (unsigned char*) calloc(SYM_NUM*OUTPUT_LEN, sizeof(unsigned char));
   sig_real = (double*) calloc(SYM_NUM*OUTPUT_LEN, sizeof(double));
   sig_img = (double*) calloc(SYM_NUM*OUTPUT_LEN, sizeof(double));
   in_ifft = (float*) calloc(SYM_NUM*FFT_N*2, sizeof(float));
   
   pilot_out = (comp_t**)malloc(sizeof(comp_t*) * SYM_NUM);
   for (i = 0; i < SYM_NUM; i++) {
    pilot_out[i] = (comp_t*) calloc(FFT_N, sizeof(comp_t));
   }
   
   ifft_input = (cedr_cmplx_flt_type**) malloc(sizeof(cedr_cmplx_flt_type*) * SYM_NUM);
   ifft_output = (cedr_cmplx_flt_type**) malloc(sizeof(cedr_cmplx_flt_type*) * SYM_NUM);

   cyclic_out = (comp_t*) calloc(SYM_NUM*TOTAL_LEN, sizeof(comp_t));
   pre_out = (comp_t*) calloc(SYM_NUM*TOTAL_LEN + PREAMBLE_LEN + 2048, sizeof(comp_t));
   txdata = (comp_t*) calloc(TX_DATA_LEN, sizeof(comp_t));
   
   // Object Instatiation
   init_viterbiEncoder();
   encoderId = get_viterbiEncoder();
   set_viterbiEncoder(encoderId);
   #ifdef TARGET
      create_rfInf(TXMODE, RFCARD, INT_2BYTE, NO_SMPERR, 0);
   #else
      create_rfInf(TXMODE, DATFILE, INT_2BYTE, NO_SMPERR, 0);
   #endif


   readConfig();
   int txOption = _txOption;
   struct timespec start1, end1;
   float exec_time;

   int frame_count = 0;

   tf::Taskflow taskflow;
   std::map<std::string, cedr_task_config_t*> task_configs;
  
   taskflow.name("WiFi-TX");

    //printf("[wifi-tx] Generating input data\n");
    // input data generation
    auto t_DataGen = taskflow.emplace([txOption, &inbit, &user_data_len]() {
      user_data_len = txDataGen(txOption, inbit, SYM_NUM);
    }).name("txDataGen-non_kernel");

    // transmitter chain
    int first_scrambler = 0;
    int last_scrambler = SYM_NUM;
    auto t_scrambler = taskflow.for_each_index(std::ref(first_scrambler), std::ref(last_scrambler), 1, [inbit, &scram](int i){
      scrambler(USR_DAT_LEN, &inbit[i*USR_DAT_LEN], &scram[i*1024]);
    }).name("scrambler-non_kernel");

    int first_viterbi = 0;
    int last_viterbi = SYM_NUM;
    auto t_viterbi = taskflow.for_each_index(std::ref(first_viterbi), std::ref(last_viterbi), 1, [encoderId, scram, &enc_out, rate, &enc_dep_out](int i){
      #ifndef SCRAMBLER_CE_HW
        //###############################################################
        //## SW Convolution Encoder
        //###############################################################
        viterbi_encoding(encoderId, &scram[i*1024], &enc_out[i*OUTPUT_LEN]);

        //###############################################################
        //## SW Viterbi Puncturing
        //###############################################################
        viterbi_puncturing(rate, &enc_out[i*OUTPUT_LEN], &enc_dep_out[i*OUTPUT_LEN]);
      #endif
    }).name("viterbi-non_kernel");

    //###############################################################
    //## Interleaver
    //###############################################################
    int first_interleaver = 0;
    int last_interleaver = SYM_NUM;
    #ifdef SCRAMBLER_CE_HW
      auto t_interleaver = taskflow.for_each_index(std::ref(first_interleaver), std::ref(last_interleaver), 1, [scram, &intl_out](int i){
        interleaver(&scram[i*1024], OUTPUT_LEN, &intl_out[i*OUTPUT_LEN]);
    #else
      auto t_interleaver = taskflow.for_each_index(std::ref(first_interleaver), std::ref(last_interleaver), 1, [enc_dep_out, &intl_out](int i){
        interleaver(&enc_dep_out[i*OUTPUT_LEN], OUTPUT_LEN, &intl_out[i*OUTPUT_LEN]);
    #endif
    }).name("interleaver-non_kernel");

    int first_qpsk = 0;
    int last_qpsk = SYM_NUM;
    auto t_qpsk = taskflow.for_each_index(std::ref(first_qpsk), std::ref(last_qpsk), 1, [intl_out, &sig_real, &sig_img, &in_ifft](int i){
      MOD_QPSK(OUTPUT_LEN, &intl_out[i*OUTPUT_LEN], &sig_real[i*OUTPUT_LEN], &sig_img[i*OUTPUT_LEN], &in_ifft[i*FFT_N*2]);
    }).name("qpsk-non_kernel");

    //###############################################################
    //## Pilot Insertion
    //###############################################################
      
    int first_pilotInsertion = 0;
    int last_pilotInsertion = SYM_NUM;
    auto t_pilotInsertion = taskflow.for_each_index(std::ref(first_pilotInsertion), std::ref(last_pilotInsertion), 1, [in_ifft, &pilot_out](int i){
      comp_t *ifft_in = (comp_t *) &in_ifft[i*FFT_N*2];
      //printf("[wifi-tx] Inserting pilot\n");
      pilotInsertion(ifft_in, pilot_out[i]);
    }).name("pilotInsertion-non_kernel");
      
    //###############################################################
    //## Inverse FFT: Input data initialization
    //###############################################################

    int first_dataConv = 0;
    int last_dataConv = SYM_NUM;
    auto t_dataConv = taskflow.for_each_index(std::ref(first_dataConv), std::ref(last_dataConv), 1, [pilot_out, &ifft_input, &ifft_output](int i){
      int random_wait_time_in_us; // = 15;
      float random_wait;
      
      fft_id = 1;

      hw_fft_busy = 1;
      
      ifft_input[i] = (cedr_cmplx_flt_type*) malloc(sizeof(cedr_cmplx_flt_type) * FFT_N);
      ifft_output[i] = (cedr_cmplx_flt_type*) malloc(sizeof(cedr_cmplx_flt_type) * FFT_N);
      for (int k = 0; k < FFT_N; k++){
        ifft_input[i][k].re = (cedr_re_flt_type) pilot_out[i][k].real;
        ifft_input[i][k].im = (cedr_re_flt_type) pilot_out[i][k].imag;
      }
    }).name("dataConv-non_kernel");
    
    //###############################################################
    //## Inverse FFT: IFFT computation
    //###############################################################

    size_t size = FFT_N;
    bool forwardTrans = false;

    // Parallel cedr_FFT portion ----------------------------------------------------------------

    // Non-blocking API call to cedr_FFT
    int first_ifft = 0;
    int last_ifft = SYM_NUM;
    task_configs["ifft-CEDR_FFT_flt-" + std::to_string(SYM_NUM)] = NULL;
    cedr_task_config_t *&confp = task_configs["ifft-CEDR_FFT_flt-" + std::to_string(SYM_NUM)];
    auto t_ifft = taskflow.for_each_index(std::ref(first_ifft), std::ref(last_ifft), 1, [ifft_input, &ifft_output, size, forwardTrans, argc, argv, &confp](int i){
      CEDR_FFT_flt(ifft_input[i], ifft_output[i], size, forwardTrans, .argc=argc, .argv=argv);
      //CEDR_FFT_flt(ifft_input[i], ifft_output[i], size, forwardTrans, .argc=argc, .argv=argv, .config=confp, .task_id=i);
    }).name("ifft-CEDR_FFT_flt-" + std::to_string(SYM_NUM));

    //--------------------------------------------------------------------------------------------*/

    int first_fftshift = 0;
    int last_fftshift = SYM_NUM;
    auto t_fftshift = taskflow.for_each_index(std::ref(first_fftshift), std::ref(last_fftshift), 1, [&ifft_input, &ifft_output, &pilot_out](int i){
      free(ifft_input[i]);
    
      // FFT-shift ##########################################################
      cplx buf[FFT_N], tmp;
      for (int k = 0; k < FFT_N; k++){
        buf[k] = (double) ifft_output[i][k].re + (double) ifft_output[i][k].im * I;
      }
      free(ifft_output[i]);

      int n2 = FFT_N/2;
      buf[0] = buf[0]/(double)FFT_N;
      buf[n2] = buf[n2]/(double)FFT_N;
      for(int l=1; l<n2; l++) {
        tmp = buf[l]/(double)FFT_N;
        buf[l] = buf[FFT_N-l]/(double)FFT_N;
        buf[FFT_N-l] = tmp;
      }

      for (int k=0; k<FFT_N; k++){
        pilot_out[i][k].real = (float)buf[k].real();
        pilot_out[i][k].imag = (float)buf[k].imag();
      }
      //###################################################
    }).name("fftshift-non_kernel");

    int first_crc = 0;
    int last_crc = SYM_NUM;
    auto t_crc = taskflow.for_each_index(std::ref(first_crc), std::ref(last_crc), 1, [pilot_out, &cyclic_out](int i){
      //###############################################################
      //## CRC
      //###############################################################
      cyclicPrefix(pilot_out[i], &cyclic_out[i*(TOTAL_LEN)], FFT_N, CYC_LEN);

    }).name("crc-non_kernel");
  
    auto t_preamble = taskflow.emplace([cyclic_out, &pre_out]() {
      // zero padding
      for(int i=0; i<256; i++) { // 512 zero pad
        pre_out[i].real = pre_out[i].imag = 0;
      }
      // data payload
      preamble(cyclic_out, &pre_out[256], SYM_NUM*(TOTAL_LEN)); // 322 preamble + SYM_NUM*80
      // total = 256 + 322 + SYM_NUM*(64+16)
    }).name("preamble-non_kernel");

    auto t_frameDup = taskflow.emplace([pre_out, &txdata]() {
      // frame duplications
      int frameLength = 256 + PREAMBLE_LEN + SYM_NUM*(TOTAL_LEN);//in complex number
      int i;
      for(i=0; i<TX_DATA_LEN - frameLength; i+=frameLength) {
        for(int j=0; j<frameLength; j++) {
            txdata[i+j].real = pre_out[j].real;
            txdata[i+j].imag = pre_out[j].imag;
        }
      }

      for( ; i<TX_DATA_LEN; i++) {
        txdata[i].real = 0;
        txdata[i].imag = 0;
      }
    }).name("frameDup-non_kernel");

    auto t_txDataDump = taskflow.emplace([txOption, txdata]() {
      // send data
      //txDataDump(pre_out, SYM_NUM*(64+16)+322);

      #ifdef TARGET
      if(txOption == 0) {
        printf("press enter to transmit: ");
        getchar();
      }
      for(int i=0; i<3; i++) {
        rfInfWrite(txdata, TX_DATA_LEN);
      }
      #else
      //printf("- RF file dump!!\n");
      //rfInfWrite(txdata, TX_DATA_LEN);
      FILE *txdata_file = fopen("./output/wifi_tx_output.txt", "w");
        if (txdata_file != NULL) {
        for(int i = 0; i < TX_DATA_LEN; i++){
          fprintf(txdata_file,"%f %f\n", txdata[i].real, txdata[i].imag);
        }
        fclose(txdata_file);
      }
    
      //printf("- completion!!\n\n");
      #endif
    }).name("txDataDump-non_kernel");

  t_scrambler.succeed(t_DataGen);
  t_viterbi.succeed(t_scrambler);
  t_interleaver.succeed(t_viterbi);
  t_qpsk.succeed(t_interleaver);
  t_pilotInsertion.succeed(t_qpsk);
  t_dataConv.succeed(t_pilotInsertion);
  t_ifft.succeed(t_dataConv);
  t_fftshift.succeed(t_ifft);
  t_crc.succeed(t_fftshift);
  t_preamble.succeed(t_crc);
  t_frameDup.succeed(t_preamble);
  t_txDataDump.succeed(t_frameDup);
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
  printf("Average: %.0f ns\n", average);
  printf("Median: %.0f ns\n", median);
  printf("Min: %lld ns\n", min);
  printf("Max: %lld ns\n", max);
  printf("%.0f\n%.0f\n%lld\n%lld\n", average, median, min, max);
  
  
   //fclose(txdata_file);

   // Free out memories
   free(inbit);
   free(scram);
   free(enc_out);
   free(enc_dep_out);
   free(intl_out);
   free(sig_real);
   free(sig_img);
   free(in_ifft);
   for (i = 0; i < SYM_NUM; i++) {
    free(pilot_out[i]);
   }
   free(pilot_out);
   
   free(ifft_input);
   free(ifft_output);

   free(cyclic_out);
   free(pre_out);
   free(txdata);

  // printf("[nk] Non-kernel thread execution is complete...\n");
  // printf("[WiFi-TX] Execution is complete...\n");
  return 0;
}

void random_wait_time(int random_wait_time_in_us) {
    for(int k = 0; k < random_wait_time_in_us; k++) {
        for(int i = 0; i < 170; i++);
    }
}

void readConfig() {

    FILE *cfp;
    char buf[1024];

    cfp = fopen("./input/tx.cfg", "r");
    if(cfp == NULL) {
        char currWorkingDir[PATH_MAX];
        printf("fail to open config file: %s\n", strerror(errno));
        getcwd(currWorkingDir, PATH_MAX);
        printf("current working dir: %s\n", currWorkingDir);
        exit(1);
    }

    fgets(buf, 1024, cfp);
    sscanf(buf, "%d", &_txOption);
    // printf("- %s\n", (txOption == 0) ? "Tx fixed string" : "Tx variable string");
}
