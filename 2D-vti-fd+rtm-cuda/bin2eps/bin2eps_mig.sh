#!/bin/sh

n1=1000
n2=2000

d1=0.005
d2=0.005

   file='../trial/waxian_3_stkadcigs.dat'
outfile='../trial/waxian_3_stkadcigs.eps'


label2='width[km]'
label1='depth[km]'

title='migration'

#red blue green
cchar='wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

#white red blue
#cchar='wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0'


#migration
psimage < $file n1=$n1 d1=$d1 n2=$n2 d2=$d2   \
       label1=$label1 label2=$label2 title=$title \
       width=8 height=4 $cchar perc=98 > $outfile

gimp $outfile &
