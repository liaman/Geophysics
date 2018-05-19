CUDA_INSTALL_PATH = /usr/local/cuda-7.5
MPI_INSTALL_PATH = /home/leonvel/software/mpi/mpich

NVCC = $(CUDA_INSTALL_PATH)/bin/nvcc
MPICC = $(MPI_INSTALL_PATH)/bin/mpicc

LDFLAGS = -L$(CUDA_INSTALL_PATH)/lib64
LIB = -lcudart -lcurand

CFILES = Toa_3dfd_mpi.c
CUFILES = Toa_3dfd_cuda.cu
OBJECTS = Toa_3dfd_mpi.o Toa_3dfd_cuda.o 
EXECNAME = fd

#CFILES = test.c
#CUFILES = test_cuda.cu
#OBJECTS = test.o test_cuda.o 
#EXECNAME = test

all:
	$(MPICC) -c $(CFILES)
	$(NVCC) -c $(CUFILES)
	$(MPICC) -o $(EXECNAME) $(LDFLAGS) $(LIB) $(OBJECTS)

clean:
	rm -f *.o $(EXECNAME)

