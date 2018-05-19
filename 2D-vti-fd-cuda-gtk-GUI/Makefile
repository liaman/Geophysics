CUDA_INSTALL_PATH = /usr/local/cuda-7.5
MPI_INSTALL_PATH = /usr

NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
GCC = $(MPI_INSTALL_PATH)/bin/gcc

LDFLAGS = -L$(CUDA_INSTALL_PATH)/lib64
LIB = -lcudart -lcurand
LGTK = `pkg-config --cflags --libs gtk+-2.0`

CFILES = GTK_VTI.c
CUFILES = GPU_VTImodel_kernels.cu
OBJECTS = GTK_VTI.o GPU_VTImodel_kernels.o 
EXECNAME = a

all:
	$(GCC) -c $(CFILES)  $(LGTK)
	$(NVCC) -c $(CUFILES)
	$(GCC) -o $(EXECNAME) $(LDFLAGS) $(LIB) $(LGTK) $(OBJECTS)

clean:
	rm -f *.o $(EXECNAME)

