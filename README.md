# Geophysical-data-processing-methods
Geophysical data processing methods, Program code implementation of geophysical data processing methods, including parallel algorithms (cuda, mpich, openmp, etc.) with finite difference forward modeling, inverse time migration, and full waveform inversion

## Key Words 关键词

### Geophysical Medium Key Words-地球物理介质关键词
* 2D (two dimension, 二维)
* 3D (three dimension, 三维)
* VTI (水平横向各向同性介质, a kind of anisotropic）
* TTI (倾斜横向各向同性介质, s kind of anisotropic)
* ISO (各向同性介质, isotropic)

### Geophysical Key Words-地球物理方法关键词
* FD (Finite Difference, 有限差分)
* RTM (Reverse Time Migration, 逆时偏移)
* FWI (Full Waveform Inversion, 全波形反演)
* TOMO (Tomography, 层析成像)
* RayTracing (射线追踪)

### Geophysical Gather Key Words-地球物理道集关键词
* ADCIGs (Angel Domain Common Imaging Gather, 角度域共成像点道集)

### Geophysical file format Key Words-地球物理文件格式关键词
* binary/data (裸数据)
* SU (Seismic Unix格式文件.su)
* segy (seg-Y格式数据.segy或者.sgy)

### Calculation Tool Key Words-计算机关键词
* HPC (High-performance computing, 高性能运算)
* omp (openmp, 多cpu线程)
* mpi (mpich or openmpi, 多核/多计算节点)
* cuda (NVIDIA CUDA-Toolkit, version>=7.5, 多gpu线程)

### GUI Tool Key Words-图形用户界面开发关键词
* GUI (Graphical User Interface, 图形用户界面)
* gtk (GIMP Toolkit, version>=2.0, such as "gtk+-2.0" or "gtk+-3.0")
* Qt (Qt-Creater Toolkit, C++开发工具)
* software (软件)

### Other Key Words-其他关键字
* utils （工具)
* mix (混合编程, mix-mpi+cuda表示mpi和cuda混合编程)

## 详情请见每个文件夹中的README.md文件

### Dependence && Envrionment
* OS:  linux
* Compiler:  gcc, nvcc, mpicc, javac
* Software:  gcc, cuda, mpich/openmpi, openmp, Qt-Creater, gtk+-2.0/3.0, JDK, make/cmake
### some usage example
* to use gtk
```c
#include<gtk/gtk.h>
```
* to use cuda
```c
#include<cuda_runtime.h>
```
* to use mpich/openmpi
```c
#include<mpi.h>
```
* to use openmp
```c
#include<omp.h>
#pragma omp parallel for//such as
```

