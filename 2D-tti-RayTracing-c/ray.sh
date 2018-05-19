gcc -o ray rayt2d_no_surf_rslocation.c -lm 

./ray

suaddhead < time.dat ns=100 |
suximage perc=99&

#suaddhead < tvfile ns=200 |     #走时表变化的输出文件
#suximage perc=99&

#suaddhead < csfile ns=200 |
#suximage perc=99&
