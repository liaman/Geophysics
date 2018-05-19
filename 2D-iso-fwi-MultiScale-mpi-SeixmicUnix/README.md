# 基于MPICH的三维多尺度波形反演
## Scale_Filter_Smooth.c
对炮记录进行平滑处理
## Scale_PRP_model_rw_mu.c
对应尺度的正演模拟，需要读取子波
## Scale_changegrid.c
网格延拓或限制
## Scale_changev.c

## Scale_filter.sh
用Seismic Unix进行滤波
## Scale_fwi_sl_zm_rw_mu.c
多尺度波形反演源代码，需要读取对应尺度子波
## Scale_getVel_from_velRecord.c
从一些列输出数据中提取切片
## Scale_read_nloop_vel.c
同上
## Scale_smooth.sh
平滑shell
## Scale_vel_model.c Scale_vel_model2.c
生成模型
## filter.sh
滤波
## smooth.f95
平滑的fortran程序
