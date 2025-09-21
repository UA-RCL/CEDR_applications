ARCH ?= x86

INCLUDES =  -I../../../libcedr/
FLAGS = -DMEASURE_OVERHEAD -shared -fPIC -O3 #-DDEBUG #Use -DMEASURE_OVERHEAD to measure App-based (blocking app only) API executions in cycles and print them into a csv file
STANDALONE_FLAGS = -DMEASURE_OVERHEAD -DCPU_ONLY -O3 #-DDEBUG

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
STANDALONE_LIB = ../../../libcedr/aarch64_lib/libgsl.a ../../../libcedr/aarch64_lib/libgslcblas.a -lm -lpthread -lstdc++ -Wno-unused-result -DUSE_ASM_TIMERS
ASM_TIMER_FLAG = -DUSE_ASM_TIMERS
endif

CPU_ONLY_SOURCES = ../../../libcedr/cpu/*.cpp
SOURCES_CEDR = application_3zip_cedr.cpp
SOURCES_ZIP_TASLFLOW = application_3zip_taskflow.cpp	
SOURCE_CUDA = application_3zip_cuda.cu

#OpenMP Configurations:
omp_sch ?= static
omp_n ?= 8
OMP_FLAGS= -fopenmp -DSCHEDULER=$(omp_sch) -DTHREADS=$(omp_n)
MAIN_FUNC_TIME?=-DTIME_CAPTURE

cedr: $(SOURCES_CEDR)
	$(CXX) $(INCLUDES) $(ASM_TIMER_FLAG) $(FLAGS) $(OMP_FLAGS) $(SOURCES_CEDR) -o zip3x_cedr_nb-$(ARCH).so
	$(CXX) $(INCLUDES) $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) $(OMP_FLAGS) $(SOURCES_CEDR) $(CPU_ONLY_SOURCES) -o zip3x_cedr_nb-$(ARCH).out ${STANDALONE_LIB}

cuda: $(SOURCE_CUDA)
	nvcc $(INCLUDES) -Xcompiler $(OMP_FLAGS) -DTIME_CAPTURE application_3zip_cuda.cu -o 3zip_cuda.out

cuda_opt: $(SOURCE_CUDA)
	nvcc $(INCLUDES) -Xcompiler $(OMP_FLAGS) -O3 -DTIME_CAPTURE application_3zip_cuda.cu -o 3zip_cuda-opt.out

SHELL := /bin/bash
taskflow: INCLUDES += -I../../../taskflow
taskflow: CPU_ONLY_SOURCES += ../../../libcedr/taskflow/*.cpp
taskflow: FLAGS += -std=c++17
taskflow: STANDALONE_FLAGS += -std=c++17
taskflow: $(SOURCES_ZIP_TASLFLOW)
	$(CXX) $(INCLUDES) $(ASM_TIMER_FLAG) $(FLAGS) $(SOURCES_ZIP_TASLFLOW) -o zip3x_taskflow_nb-$(ARCH).so
	$(CXX) $(INCLUDES) -DDAG_PARSE $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) $(SOURCES_ZIP_TASLFLOW) $(CPU_ONLY_SOURCES) -o zip3x_taskflow_nb-$(ARCH).out ${STANDALONE_LIB}

clean:
	rm -f *.so
	rm -f *.o
	rm -f *.out

