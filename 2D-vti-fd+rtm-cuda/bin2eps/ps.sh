#!/bin/sh

nx=1001
nxa=770
nz=301


dx=5
dz=5


width=20
depth=6

cchar=''
#'cmap=rgb1 wrgb=1.0,0,0 grgb=1.0,1.0,1.0 brgb=0,0,1.0' 
labelsize=30

perc=99
#
psimage  <  waxian_vel_1001_301.dat n1=$nz d1=$dz n2=$nx d2=$dx label1='depth[m]' label2='width[m]'  \
             $cchar  width=$width height=$depth  labelsize=$labelsize perc=$perc > eps_migration.eps



#psimage  <  thrust_vel_711_300.bin n1=$nz d1=$dz n2=$nx d2=$dx label1='depth[m]' label2='width[m]'  \
#             $cchar  width=$width height=$depth legend=1 lstyle=vertright  lheight=5 labelsize=$labelsize perc=$perc > eps_thrust_vel_711_300.eps

#psimage  <  thrust_epsilon_711_300.bin n1=$nz d1=$dz n2=$nx d2=$dx label1='depth[m]' label2='width[m]'  \
#             $cchar  width=$width height=$depth legend=1 lstyle=vertright lheight=5 labelsize=$labelsize perc=$perc > eps_thrust_epsilon_711_300.eps

#psimage  <  thrust_delta_711_300.bin n1=$nz d1=$dz n2=$nx d2=$dx label1='depth[m]' label2='width[m]'  \
#             $cchar  width=$width height=$depth legend=1 lstyle=vertright lheight=5 labelsize=$labelsize perc=$perc > eps_thrust_delta_711_300.eps



