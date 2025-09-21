ARCH ?= x86

INCLUDES =  -I../../../libcedr/
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
	STANDALONE_LIB = ../../../libcedr/aarch64_lib/libgsl.a ../../../libcedr/aarch64_lib/libgslcblas.a -lm -lpthread -lstdc++ -Wno-unused-result -DUSE_ASM_TIMERS
	ASM_TIMER_FLAG = -DUSE_ASM_TIMERS
endif

CPU_ONLY_SOURCES = ../../../libcedr/cpu/*.cpp

DEFINES = -DPRINT_IMAGES

omp_sch ?= dynamic
omp_n ?= 4
OMP_FLAGS= -fopenmp -DSCHEDULER=$(omp_sch) -DTHREADS=$(omp_n)
APP_NAME = lane_detection

2Dconv: 
	$(CXX) $(INCLUDES) -DCONV -DPRINT_IMAGES $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(FLAGS)  track_2Dconv.cpp -o $(APP_NAME)-2Dconv-$(ARCH).so
	$(CXX) $(INCLUDES) -DCONV -DPRINT_IMAGES $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) track_2Dconv.cpp $(CPU_ONLY_SOURCES) -o $(APP_NAME)-2Dconv-$(ARCH).out ${STANDALONE_LIB}

uint8: 
	$(CXX) $(INCLUDES) -DCONV -DPRINT_IMAGES $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(FLAGS)  track_2Dconv_uint8.cpp -o $(APP_NAME)-2Dconv-uint8-$(ARCH).so
	$(CXX) $(INCLUDES) -DCONV -DPRINT_IMAGES $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) track_2Dconv_uint8.cpp $(CPU_ONLY_SOURCES) -o $(APP_NAME)-2Dconv-uint8-$(ARCH).out ${STANDALONE_LIB}

SHELL := /bin/bash
taskflow: INCLUDES += -I../../../taskflow
taskflow: CPU_ONLY_SOURCES += ../../../libcedr/taskflow/*.cpp
taskflow: FLAGS += -std=c++17
taskflow: STANDALONE_FLAGS += -std=c++17
taskflow: $(SOURCES_TASKFLOW)
	$(CXX) $(INCLUDES) -DCONV -DPRINT_IMAGES $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(FLAGS)  track_2Dconv_taskflow.cpp -o $(APP_NAME)-2Dconv-taskflow-$(ARCH).so
	$(CXX) $(INCLUDES) -DCONV -DPRINT_IMAGES -DDAG_PARSE $(MAIN_FUNC_TIME) $(ASM_TIMER_FLAG) $(STANDALONE_FLAGS) track_2Dconv_taskflow.cpp $(CPU_ONLY_SOURCES) -o $(APP_NAME)-2Dconv-taskflow-$(ARCH).out ${STANDALONE_LIB}

clean:
	-rm -f *.so
	-rm -f *.out
	-find . -maxdepth 1 -type f -name '*.png' ! -name 'image.png' -exec rm -f {} +
	-rm -f o*.jpeg

