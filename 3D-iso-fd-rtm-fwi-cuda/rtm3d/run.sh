#!/bin/sh


ns=10

nt=1201
kt=600
n1=200
npd=50 
mm=4

# Complie
nvcc -o rtm3d -lcublas -lcufft -lcusparse gpurtm3d.cu -L/home/leonvel/software/madagascar/madagascar-1.7/lib -L/usr/lib64/atlas -lrsf++ -llapack -lcblas -lumfpack -lcholmod -lamd -lcamd -lcolamd -lccolamd -lrsfpwd -lrsf -lm -lf77blas -lcblas -latlas -llapack -lcblas -lgomp -lfftw3f -lfftw3f_threads -I /home/leonvel/software/madagascar/madagascar-1.7/include
# Run 
./rtm3d < vel.rsf   snap=snap.rsf shot=shot.rsf gradient=grad.rsf verb=n nt=$nt kt=$kt dt=0.001 fm=30 ns=$ns szbeg=1 sxbeg=10 sybeg=10 jsz=0 jsx=20 jsy=20 

#show
#sfbyte < Fw.rsf bar=gbar.rsf gainpanel=all allpos=y |sfgrey3 bar=gbar.rsf scalebar=y color=j movie=1 frame1=100 frame2=100 frame3=100 flat=n |sfpen&
ximage < grad.rsf@ n1=200 perc=99&
ximage < snap.rsf@ n1=200 perc=99&
ximage < shot.rsf@ n1=$nt perc=99&
