# Geophysical-data-processing-methods
Geophysical data processing methods, Program code implementation of geophysical data processing methods, including parallel algorithms (cuda, mpich, openmp, etc.) with finite difference forward modeling, inverse time migration, and full waveform inversion

# Key Words
* 2D (two dimension)
* 3D (three dimension)
* FD (Finite Difference)
* RTM (Reverse Time Migration)
* FWI (Full Waveform Inversion)
* VTI (水平横向各向同性介质, a kind of anisotropic）
* TTI (倾斜横向各向同性介质, s kind of anisotropic)
* ISO (各向同性介质, isotropic)
* mpi (mpich or openmpi)
* cuda (NVIDIA CUDA-Toolkit, version>=7.5, 多gpu线程)
* gtk (GIMP Toolkit, version>=2.0, such as "gtk+-2.0" or "gtk+-3.0")
* omp (openmp, 多cpu线程)
* software (软件)
* GUI (Graphical User Interface, 图形用户界面)

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
