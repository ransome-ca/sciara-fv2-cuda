.PHONY: all clean

THREADS	    		?= 8

INPUT_CONFIG		:= ./data/2006/2006_000000000000.cfg
OUTPUT_CONFIG		:= ./data/2006/output_2006
OUTPUT				:= ./data/2006/output_2006_000000016000_Temperature.stt  #md5sum: 0c071cd864046d3c6aaf30997290ad6c   /   704a4a65d1890589e952b155d53b110d
STEPS				:= 1000
REDUCE_INTERVL		:= 1000
THICKNESS_THRESHOLD := 1.0  #resulting in 16000 steps



CXX					:= g++
NVCC				:= nvcc


EXTRA_CUFLAGS       ?=
EXTRA_CXXFLAGS      ?=


STSRCS						:= $(shell find src/standard -name '*.cpp')
STHDRS						:= $(shell find src/standard -name '*.h')
CUSRCS						:= $(shell find src/cuda -name '*.cu')
CUHDRS						:= $(shell find src/cuda -name '*.h')
CUSRCS_STRAIGHTFORWARD		:= $(shell find src/cuda-straightforward -name '*.cu')
CUHDRS_STRAIGHTFORWARD		:= $(shell find src/cuda-straightforward -name '*.h')
CUSRCS_TILED_HALO			:= $(shell find src/cuda-tiled-halo -name '*.cu')
CUHDRS_TILED_HALO			:= $(shell find src/cuda-tiled-halo -name '*.h')
CUSRCS_MULTI_GPU			:= $(shell find src/cuda-multi-gpu -name '*.cu')
CUHDRS_MULTI_GPU			:= $(shell find src/cuda-multi-gpu -name '*.h')
CUSRCS_TILED_NO_HALO		:= $(shell find src/cuda-tiled-no-halo -name '*.cu')
CUHDRS_TILED_NO_HALO		:= $(shell find src/cuda-tiled-no-halo -name '*.h')

SERIAL						:= sciara-fv2-serial
PARALLEL					:= sciara-fv2-parallel
CUDA_STRAIGHTFORWARD		:= sciara-fv2-cuda-straightforward
CUDA_TILED_HALO				:= sciara-fv2-cuda-tiled-halo
CUDA_MULTI_GPU				:= sciara-fv2-cuda-multi-gpu
CUDA_TILED_NO_HALO			:= sciara-fv2-cuda-tiled-no-halo


DEFINES		:=
CXXFLAGS	:=-O3 $(addprefix -D,$(DEFINES))
CUFLAGS		:=-O3 -I src/cuda -gencode=arch=compute_52,code="compute_52" -fmad=true --maxrregcount=64 $(addprefix -D,$(DEFINES)) $(EXTRA_CUFLAGS)


all: $(SERIAL) $(PARALLEL) $(CUDA_STRAIGHTFORWARD) $(CUDA_TILED_HALO) $(CUDA_MULTI_GPU) $(CUDA_TILED_NO_HALO)

serial: $(SERIAL)
parallel: $(PARALLEL)
cuda-straightforward: $(CUDA_STRAIGHTFORWARD)
cuda-tiled-halo: $(CUDA_TILED_HALO)
cuda-multi-gpu: $(CUDA_MULTI_GPU)
cuda-tiled-no-halo: $(CUDA_TILED_NO_HALO)


$(SERIAL): $(STSRCS) $(STHDRS)
	$(CXX) $(CXXFLAGS) -o $(SERIAL) $(STSRCS) $(LIBS)

$(PARALLEL): $(STSRCS) $(STHDRS)
	$(CXX) $(CXXFLAGS) -fopenmp -o $(PARALLEL) $(STSRCS) $(LIBS)

$(CUDA_STRAIGHTFORWARD): $(CUSRCS) $(CUHDRS) $(CUSRCS_STRAIGHTFORWARD) $(CUHDRS_STRAIGHTFORWARD)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_STRAIGHTFORWARD) -o $@ $(LIBS)

$(CUDA_TILED_HALO): $(CUSRCS) $(CUHDRS) $(CUSRCS_TILED_HALO) $(CUHDRS_TILED_HALO)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_TILED_HALO) -o $@ $(LIBS)

$(CUDA_MULTI_GPU): $(CUSRCS) $(CUHDRS) $(CUSRCS_MULTI_GPU) $(CUHDRS_MULTI_GPU)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_MULTI_GPU) -o $@ $(LIBS) -Xcompiler -fopenmp

$(CUDA_TILED_NO_HALO): $(CUSRCS) $(CUHDRS) $(CUSRCS_TILED_NO_HALO) $(CUHDRS_TILED_NO_HALO)
	$(NVCC) $(CUFLAGS) $(CUSRCS) $(CUSRCS_TILED_NO_HALO) -o $@ $(LIBS)


run-serial: $(SERIAL)
	./$(SERIAL) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "704a4a65d1890589e952b155d53b110d"

run-parallel: $(PARALLEL)
	OMP_NUM_THREADS=$(THREADS) ./$(PARALLEL) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "704a4a65d1890589e952b155d53b110d"

run-cuda-straightforward: $(CUDA_STRAIGHTFORWARD)
	$(PROF) ./$(CUDA_STRAIGHTFORWARD) $(THREADS) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "704a4a65d1890589e952b155d53b110d"

run-cuda-tiled-halo: $(CUDA_TILED_HALO)
	$(PROF) ./$(CUDA_TILED_HALO) $(THREADS) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "704a4a65d1890589e952b155d53b110d"

run-cuda-multi-gpu: $(CUDA_MULTI_GPU)
	$(PROF) ./$(CUDA_MULTI_GPU) $(THREADS) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "704a4a65d1890589e952b155d53b110d"

run-cuda-tiled-no-halo: $(CUDA_TILED_NO_HALO)
	$(PROF) ./$(CUDA_TILED_NO_HALO) $(THREADS) $(INPUT_CONFIG) $(OUTPUT_CONFIG) $(STEPS) $(REDUCE_INTERVL) $(THICKNESS_THRESHOLD) && md5sum $(OUTPUT) && echo "704a4a65d1890589e952b155d53b110d"



debug-cuda-straightforward: $(CUDA_STRAIGHTFORWARD)
	cuda-gdb ./$(CUDA_STRAIGHTFORWARD)

debug-cuda-tiled-halo: $(CUDA_TILED_HALO)
	cuda-gdb ./$(CUDA_TILED_HALO)

debug-cuda-multi-gpu: $(CUDA_MULTI_GPU)
	cuda-gdb ./$(CUDA_MULTI_GPU)

clean:
	$(RM) $(SERIAL) $(PARALLEL) $(CUDA_STRAIGHTFORWARD) $(CUDA_TILED_HALO) $(CUDA_MULTI_GPU) $(CUDA_TILED_NO_HALO) $(OUTPUT_CONFIG)*



METRICS 	?= flop_count_dp

profile-cuda-straightforward: $(CUDA_STRAIGHTFORWARD)
	nvprof --log-file roofline.txt --print-summary $(addprefix --kernels ,$(KERNELS)) $(addprefix -m ,$(METRICS)) ./$(CUDA_STRAIGHTFORWARD)

profile-cuda-tiled-halo: $(CUDA_TILED_HALO)
	nvprof --log-file roofline.txt --print-summary $(addprefix --kernels ,$(KERNELS)) $(addprefix -m ,$(METRICS)) ./$(CUDA_TILED_HALO)

profile-cuda-tiled-no-halo: $(CUDA_TILED_NO_HALO)
	nvprof --log-file roofline.txt --print-summary $(addprefix --kernels ,$(KERNELS)) $(addprefix -m ,$(METRICS)) ./$(CUDA_TILED_NO_HALO)

profile-cuda-multi-gpu: $(CUDA_MULTI_GPU)
	nvprof --log-file roofline.txt --print-summary $(addprefix --kernels ,$(KERNELS)) $(addprefix -m ,$(METRICS)) ./$(CUDA_MULTI_GPU)