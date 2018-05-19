#!/bin/sh

n1=2500
n2=2500

d1=5
d2=5

na=60

binary='layers_v0.80e0.50d0.50_adcigs.dat'

suaddhead <   $binary    ns=$n1   |  
sushw   key=cdp  a=0   b=1  c=0       j=$na  |
sushw   key=dt       a=$d1  >   tmp.su 

segyhdrs <tmp.su |segywrite tape=tmp.segy endian=0

rm tmp.su
