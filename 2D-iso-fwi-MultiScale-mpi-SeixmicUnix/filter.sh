suaddhead< fault_shot_obs.dat ns=3000 ntracl=560 |
sushw key=dt a=600 |
sufilter f=1,3,20,30 amps=1,1,0,0 |
sustrip > fault_shot_obs_filter.dat 


ximage < fault_shot_obs_filter.dat n1=3000 perc=99&


