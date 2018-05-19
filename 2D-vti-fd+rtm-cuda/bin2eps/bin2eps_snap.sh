#!/bin/sh

n1=500
n2=500

d1=0.005
d2=0.005

   file='snap_e0.5d0.0-0.1-0.2-0.3-0.4-0.5.dat'
outfile='snap_e0.5d0.0-0.1-0.2-0.3-0.4-0.5.eps'


label2='width[km]'
label1='depth[km]'

title='snap'

#red blue green
cchar='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

#white red blue
#cchar='wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0'

#migration
psimage < $file n1=$n1 d1=$d1 n2=$n2 d2=$d2 n2tic=20 \
       label1=$label1 label2=$label2 labelsize=24 title=$title \
       width=8 height=8 $cchar  perc=99 > $outfile

gimp $outfile &
