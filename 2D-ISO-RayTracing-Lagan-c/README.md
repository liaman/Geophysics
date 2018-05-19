# lagan-raytracing-iso
Isotropic medium lagan law ray tracing

## ste_self_new_2.c
* 该文件是源代码文件

## shot_location.txt
* 该文件为炮点位置文件，其中第一列是炮点横坐标，第二列是炮点纵坐标，第三列是第一条射线方向

## tomo-2.par
* 参数文件，格式如下：
```
  iflag_surface=  1
    FN_velmodel=  vel_600_300.dat
   FN_elevation=  step1_elevation1.txt
             nx=  600
       layer_nz=  300
             dx=  2
             dz=  2
         s_step=  0.1
          nray0=  180
          nshot=  1
         dangle=  -1.000000
    start_angle=  180.0
```

## Compiled and Run
```shell
$ gcc -o raytracing step_self_new_2.c -lm -w
$ ./raytracing
```
