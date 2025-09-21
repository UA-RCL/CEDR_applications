ARCH ?= x86

INCLUDES =  -I../../../../libcedr/
FLAGS = -shared -fPIC -O3
STANDALONE_FLAGS = -DCPU_ONLY -O3

ifeq ($(ARCH),x86)
CXX=g++
CC=gcc
STANDALONE_LIB = -lm -lgsl -lgslcblas -lpthread -lstdc++ -Wno-unused-result
ASM_TIMER_FLAG = 
endif

ifeq ($(ARCH),aarch64)
CXX=aarch64-linux-gnu-g++
CC=aarch64-linux-gnu-gcc
FLAGS += -DUSE_ASM_TIMERS
STANDALONE_LIB = ../../../../libcedr/aarch64_lib/libgsl.a ../../../../libcedr/aarch64_lib/libgslcblas.a -lm -lpthread -lstdc++ -Wno-unused-result -DUSE_ASM_TIMERS
ASM_TIMER_FLAG = -DUSE_ASM_TIMERS
endif

CPU_ONLY_SOURCES = ../../../../libcedr/cpu/*.cpp

APP_NAME = wifi-tx
C_SOURCES = ../src/crc.c ../src/txData.c ../src/baseband_lib.c ../src/viterbi.c \
	../src/interleaver_deintleaver.c ../src/scrambler_descrambler.c ../src/qpsk_Mod_Demod.c \
	../src/CyclicPrefix.c ../src/Preamble_ST_LG.c ../src/diag.c ../src/rfInf.c ../src/pilot.c
CPP_SOURCES = ../src/TX.cpp ../src/fft_hs.cpp
CPP_SOURCES_NB = ../src/TX_nb.cpp ../src/fft_hs.cpp
CPP_SOURCES_TASKFLOW = ../src/TX_taskflow_nb.cpp ../src/fft_hs.cpp

SOURCES = $(C_SOURCES) $(CPP_SOURCES)
SOURCES_NB = $(C_SOURCES) $(CPP_SOURCES_NB)
SOURCES_TASKFLOW = $(C_SOURCES) $(CPP_SOURCES_TASKFLOW)

omp_sch ?= dynamic
omp_n ?= 4
OMP_FLAGS= -fopenmp -DSCHEDULER=$(omp_sch) -DTHREADS=$(omp_n)
#MAIN_FUNC_TIME?=-DTIME_CAPTURE
default: blocking nonblocking

blocking:
	$(CXX) $(INCLUDES) $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(FLAGS) $(SOURCES) -o $(APP_NAME)-$(ARCH).so
	$(CXX) $(INCLUDES) $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) $(SOURCES) $(CPU_ONLY_SOURCES) -o $(APP_NAME)-$(ARCH).out ${STANDALONE_LIB}
	
nonblocking:	
	$(CXX) $(INCLUDES) $(OMP_FLAGS) $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(FLAGS) $(SOURCES_NB) -o $(APP_NAME)_nb-$(ARCH).so
	$(CXX) $(INCLUDES) $(OMP_FLAGS) $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) $(SOURCES_NB) $(CPU_ONLY_SOURCES) -o $(APP_NAME)_nb-$(ARCH).out ${STANDALONE_LIB}
	
SHELL := /bin/bash
taskflow: INCLUDES += -I../../../../taskflow
taskflow: CPU_ONLY_SOURCES += ../../../../libcedr/taskflow/*.cpp
taskflow: FLAGS += -std=c++17
taskflow: STANDALONE_FLAGS += -std=c++17
taskflow: $(SOURCES_TASKFLOW)
	$(CXX) $(INCLUDES) $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(FLAGS) $(SOURCES_TASKFLOW) -o $(APP_NAME)_taskflow_nb-$(ARCH).so
	$(CXX) $(INCLUDES) -DDAG_PARSE $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) $(SOURCES_TASKFLOW) $(CPU_ONLY_SOURCES) -o $(APP_NAME)_taskflow_nb-$(ARCH).out ${STANDALONE_LIB}

clean:
	rm -f *.so
	rm -f *.o
	rm -f *.out	
