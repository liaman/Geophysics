#!/bin/sh

     FN='../fd/layer_shot_v2000e0.5d0.5ref.dat'
   file='../fd/layer_shot_v2000e0.5d0.5refmaxValue.dat'
outfile='../fd/layer_shot_v2000e0.5d0.5refmaxValue.eps'

echo "gcc"
gcc -o a maxValueLine.c
echo "run"

./a $FN $file

echo "done"

n1=11001
n2=2000

d1=0.0005
d2=0.005

label2='width[km]'
label1='time[s]'

title='shot'

#red blue green
cchar='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

#white red blue
#cchar='wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0'

#migration
psimage < $file n1=$n1 d1=$d1 n2=$n2 d2=$d2 n2tic=20 \
       label1=$label1 label2=$label2 labelsize=24 title=$title \
       width=8 height=8 $cchar  perc=99 bclip=1 wclip=0 brgb=1.6,0.0,0.7 grgb=0.0,0.5,1.0 wrgb=1.0,1.0,1.0 \
      > $outfile

echo "gimp"
gimp $outfile &
