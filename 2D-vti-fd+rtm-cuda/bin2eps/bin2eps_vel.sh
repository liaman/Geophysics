#!/bin/sh

n1=1000
n2=2000

d1=0.005
d2=0.005

   file='../trial/waxian_v_2000_1000.dat'
outfile='../trial/waxian_v_2000_1000.eps'

label2='width[km]'
label1='depth[km]'

title='media'

#red blue green
cchar='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

#white red blue
#cchar='wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0'

#migration
psimage < $file n1=$n1 d1=$d1 n2=$n2 d2=$d2 n2tic=20 \
       label1=$label1 label2=$label2 title=$title \
       width=8 height=4 $cchar > $outfile

gimp $outfile &
