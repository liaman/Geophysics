#! /bin/sh

dt=0.005              # time sample interval in ray tracing             
nt=401        # number of time samples in ray tracing           
                                                                        
fz=0                   #first depth sample in velocity                  
nz=200        # number of depth samples in velocity             
dz=10         #depth interval in velocity      
                
fx=0          #         first lateral sample in velocity                
nx=200        # number of lateral samples in velocity           
dx=10         #lateral interval in velocity  
      
nxs=10
dxs=100


#smooth2 n1=$nz n2=$nx r1=5 r2=5 < vel_600_300.dat > vel_600_300smooth.bin

xmovie < vel200200.dat n1=$nz n2=$nx title="Vel Model" &

# use rayt2d to generate traveltime tables from the smooth Marmousi model
rayt2d dt=$dt nt=$nt  fz=$fz nz=$nz dz=$dz fx=$fx nx=$nx dx=$dx \
nxs=$nxs dxs=$dxs vfile=vel200200.dat tfile=veltime.bin 


# view traveltime tables with xmovie
suaddhead < veltime.bin ns=$nz | 
suximage  n1=$nz   title="Smooth frame"  &

# use rayt2d to generate traveltime tables from the hard Marmousi model
#rayt2d dt=$dt nt=$nt  fz=$fz nz=$nz dz=$dz fx=$fx nx=$nx dx=$dx \
#nxs=$nxs dxs=$dxs vfile=marmhard.bin tfile=marmhardtime.bin 

# view traveltime tables with xmovie
#suaddhead < marmhardtime.bin ns=122 | sushw key=dt a=1 | sugain trap=1 |
#suximage n1=122 n2=384 loop=1 title="hard frame=%g" width=768 height=244 &

exit 0
