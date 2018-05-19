#include<stdio.h>
#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
void main()
{
   int i,j,k;
   int nx,nz;
   float **vp;
   

   char FN1[250]={"fault_shot_obs.dat"};
   char FN2[250]={"fault_shot_obs_280_bigscale.dat"};

   nx=22400;
   nz=3000;

   vp=alloc2float(nz,nx);

   FILE *fp;
   fp=fopen(FN1,"rb");
   for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
         {
        fread(&vp[i][j],4L,1,fp);
           vp[i][j]=vp[i][j];
         }
   fclose(fp);
 
   FILE *fpc;
   fpc=fopen(FN2,"wb");
   for(i=0;i<11200;i++)
      for(j=0;j<nz;j++)
         fwrite(&vp[i*2][j],4L,1,fpc);
   fclose(fpc);

}

