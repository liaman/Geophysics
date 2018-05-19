# gpu_23D_vti_fd
Finite Difference Forward Simulation Software for 2D 3D VTI Media Based on GPU Acceleration Computation on Linux

## dependence & envrioment
* Linux
* gcc
* cuda7.5+
* gtk+-2.0 or gtk+-3.0

## Compiled and Run
* Makefile
```shell
$ make
$./binaryname
```
## Makefile
```shell
CUDA_INSTALL_PATH = /usr/local/cuda-7.5
GCC_INSTALL_PATH = /usr

NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
GCC = $(GCC_INSTALL_PATH)/bin/gcc
MPICC = mpicc

LCUDA = -L$(CUDA_INSTALL_PATH)/lib64
LIB = -lcudart -lcurand -w -lm
LGTK = `pkg-config --cflags --libs gtk+-2.0`

CFILES = CPU_vti2dfd_kernels.c
GTKFILES = GTK_main.c 
GPU2DFILES = GPU_vti2dfd_kernels.cu
GPU3DFILES = GPU_vti3dfd_kernels.cu
OBJECTS = GTK_main.o GPU_vti2dfd_kernels.o CPU_vti2dfd_kernels.o GPU_vti3dfd_kernels.o
EXECNAME = Tvtigpufd

all:
	$(GCC) -c $(GTKFILES)  $(LGTK)
	$(NVCC) -c $(GPU2DFILES)
	$(NVCC) -c $(GPU3DFILES)
	$(GCC) -c $(CFILES) 
	$(MPICC) -o $(EXECNAME) $(LCUDA) $(LIB) $(LGTK) $(OBJECTS)

	rm -f *.o *~
```

