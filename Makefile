.PHONY: all clean

THREADS	    		?= 2

INPUT_CONFIG		:= ./data/2006/2006_000000000000.cfg
OUTPUT_CONFIG		:= ./data/2006/output_2006
OUTPUT				:= ./data/2006/output_2006_000000016000_Temperature.stt  #md5sum: 0c071cd864046d3c6aaf30997290ad6c   /   704a4a65d1890589e952b155d53b110d
STEPS				:= 1000
REDUCE_INTERVL		:= 1000
THICKNESS_THRESHOLD := 1.0  #resulting in 16000 steps



CXX					:= g++
NVCC				:= nvcc


STSRCS					:= $(shell find src/standard -name '*.cpp')
STHDRS					:= $(shell find src/standard -name '*.h')
CUSRCS					:= $(shell find src/cuda -name '*.cu')
CUHDRS					:= $(shell find src/cuda -name '*.h')
CUSRCS_STRAIGHTFORWARD	:= $(shell find src/cuda-straightforward -name '*.cu')
CUHDRS_STRAIGHTFORWARD	:= $(shell find src/cuda-straightforward -name '*.h')
CUSRCS_TILED_HALO		:= $(shell find src/cuda-tiled-halo -name '*.cu')
CUHDRS_TILED_HALO		:= $(shell find src/cuda-tiled-halo -name '*.h')
CUSRCS_TILED_NO_HALO	:= $(shell find src/cuda-tiled-no-halo -name '*.cu')
CUHDRS_TILED_NO_HALO	:= $(shell find src/cuda-tiled-no-halo -name '*.h')

SERIAL					:= sciara-fv2-serial
PARALLEL				:= sciara-fv2-parallel
CUDA_STRAIGHTFORWARD	:= sciara-fv2-cuda-straightforward
CUDA_TILED_HALO			:= sciara-fv2-cuda-tiled-halo
CUDA_TILED_NO_HALO		:= sciara-fv2-cuda-tiled-no-halo

CXXFLAGS	:=-O3
CUFLAGS		:=-O3 -I src/cuda


all: $(SERIAL) $(PARALLEL) $(CUDA)

serial: $(SERIAL)
parallel: $(PARALLEL)
cuda-straightforward: $(CUDA_STRAIGHTFORWARD)
cuda-tiled-halo: $(CUDA_TILED_HALO)
cuda-tiled-no-halo: $(CUDA_TILED_NO_HALO)


$(SERIAL): $(STSRCS) $(STHDRS)
	$(CXX) $(CXXFLAGS) -o $(SERIAL) $(STSRCS) $(LIBS)

$(PARALLEL): $(STSRCS) $(STHDRS)
	$(CXX) $(CXXFLAGS) -fopenmp -o $(PARALLEL) $(STSRCS) $(LIBS)

$(CUDA_STRAIGHTFORWARD): $(CUSRCS) $(CUHDRS) $(CUSRCS_STRAIGHTFORWARD) $(CUHDRS_STRAIGHTFORWARD)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_STRAIGHTFORWARD) -o $@ $(LIBS)

$(CUDA_TILED_HALO): $(CUSRCS) $(CUHDRS) $(CUSRCS_TILED_HALO) $(CUHDRS_TILED_HALO)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_TILED_HALO) -o $@ $(LIBS)

$(CUDA_TILED_NO_HALO): $(CUSRCS) $(CUHDRS) $(CUSRCS_TILED_NO_HALO) $(CUHDRS_TILED_NO_HALO)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_TILED_NO_HALO) -o $@ $(LIBS)


run-serial: $(SERIAL)
	./$(SERIAL) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "0c071cd864046d3c6aaf30997290ad6c"

run-parallel: $(PARALLEL)
	OMP_NUM_THREADS=$(THREADS) ./$(PARALLEL) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "0c071cd864046d3c6aaf30997290ad6c"

run-cuda-straightforward: $(CUDA_STRAIGHTFORWARD)
	./$(CUDA_STRAIGHTFORWARD) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "0c071cd864046d3c6aaf30997290ad6c"

run-cuda-tiled-halo: $(CUDA_TILED_HALO)
	./$(CUDA_TILED_HALO) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "0c071cd864046d3c6aaf30997290ad6c"

run-cuda-tiled-no-halo: $(CUDA_TILED_NO_HALO)
	./$(CUDA_TILED_NO_HALO) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "0c071cd864046d3c6aaf30997290ad6c"

clean:
	$(RM) $(SERIAL) $(PARALLEL) $(CUDA_STRAIGHTFORWARD) $(CUDA_TILED_HALO) $(CUDA_TILED_NO_HALO) $(OUTPUT_CONFIG)*