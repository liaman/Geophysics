#!/bin/sh

n1=1000
n2=2000

d1=0.005
d2=0.005

   file='../trial/waxian_real_adcigs_chouxi.dat'
outfile='../trial/waxian_real_adcigs_chouxi_300_600_1000_1400_1700.eps'


label2='angle[degree]'
label1='depth[km]'

title='ADCIGs'

#red blue green
cchar='cmap=rgb1 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

#white red blue
#cchar='wrgb=1.0,1.0,1.0 grgb=1.0,0.0,0.0 brgb=0,0,1.0'


#adcigs single
#psimage < $file n1=$n1 d1=$d1 n2=$n2 f2=0 d2=$d2 x2beg=0 x2end=70 \
#      label1=$label1 label2=$label2 title=$title \
#       width=3 height=8 $cchar perc=99 > $outfile

#adcigs multi
psimage < $file n1=$n1 d1=$d1  \
       title=$title label1=$label1 \
       width=4 height=4 $cchar perc=98 >$outfile

gimp $outfile &
