#!/bin/sh

gcc -o half half.c
./half

nx=100
ny=200
nz=200


dx=5
dy=5
dz=5


size1=20
size2=10
size3=15

axiscolor='red'
gridcolor='red'
cchar=' '
#'cmap=rgb1 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0' 
labelsize=90
labelfont='Helvetica'
perc=99
#
pscube  <  a_snap_new.dat n1=$nz d1=$dz n2=$nx d2=$dx n3=$ny d3=$dy \
           label1='depth[m]-Z' label2='width[m]-X' label3='width[m]-Y'  \
           $cchar  size1=$size1 size2=$size2 size3=$size3 labelfont=$labelfont \
           gridcolor=$gridcolor axescolor=$axiscolor labelsize=$labelsize perc=$perc > a_eps_snap.eps





