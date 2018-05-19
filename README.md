# Geophysical-data-processing-methods
Geophysical data processing methods, Program code implementation of geophysical data processing methods, including parallel algorithms (cuda, mpich, openmp, etc.) with finite difference forward modeling, inverse time migration, and full waveform inversion
# Content
##FD
##RTM
##FWI
## Envrionment
* OS: 
linux
* Compiler:
gcc, nvcc, mpicc, javac
* Software: 
gcc, cuda, mpich/openmpi, openmp, Qt-Creater, gtk+-2.0/3.0, JDK, make/cmake
### ALL dependence
* Basic: 
gcc
* HPC: 
mpich, openmp, cuda, 
* GUI: 
gtk, Qt, java Swing/AWT
#### such as 
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
## FD-finite-difference
* VTI（Vertical transverse isotropy）
* TTI（Inclined transversely isotropic）
* Isotropic

## RTM-reverse-time-migration
* VTI（Vertical transverse isotropy）
* TTI（Inclined transversely isotropic）
* Isotropic

## FWI-full-waveform-inversion
* Isotropic
