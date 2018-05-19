#!/bin/sh

n1=6001
n2=1000

d1=0.0005
d2=0.005


file='../trial/sshengli_real_shot_obs.dat'
      outfile='sshengli_real_shot_obs.eps'

label2='distance[km]'
label1='time[s]'

title='shot'

#red blue green
cchar='cmap=rgb1 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0'

#migration
psimage < $file n1=$n1 d1=$d1 n2=$n2 d2=$d2 \
       titile=$title label1=$label1 label2=$label2 title=$title \
       width=8 height=8 $cchar perc=99 > $outfile

gimp $outfile &
