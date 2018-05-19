CUDA_INSTALL_PATH = /usr/local/cuda-7.5
MPI_INSTALL_PATH = /usr

NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
GCC = $(MPI_INSTALL_PATH)/bin/gcc

LDFLAGS = -L$(CUDA_INSTALL_PATH)/lib64
LIB = -lcudart -lcurand
LGTK = `pkg-config --cflags --libs gtk+-2.0`

CFILES = 4_model_gtk_gpu.c
CUFILES = 4_model_GPU_kernels.cu
OBJECTS = 4_model_gtk_gpu.o 4_model_GPU_kernels.o 
EXECNAME = a

all:
	$(GCC) -c $(CFILES)  $(LGTK)
	$(NVCC) -c $(CUFILES)
	$(GCC) -o $(EXECNAME) $(LDFLAGS) $(LIB) $(LGTK) $(OBJECTS)

clean:
	rm -f *.o $(EXECNAME)

