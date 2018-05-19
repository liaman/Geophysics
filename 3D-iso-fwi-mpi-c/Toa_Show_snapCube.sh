#!/bin/sh


n11=100
n22=100
n33=100

bclip11=6.34e+09 
wclip11=-8.27e+09

pscube < a2.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a2" >a2.eps

pscube < a3.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a3" >a3.eps

pscube < a4.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a4" >a4.eps

pscube < a5.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a5" >a5.eps

pscube < a6.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a6" >a6.eps

pscube < a7.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a7" >a7.eps

pscube < a8.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a8" >a8.eps

pscube < a9.dat n1=$n11 n2=$n22 n3=$n33 bclip=$bclip11 wclip=$wclip11 perc=99  title="a9" >a9.eps

