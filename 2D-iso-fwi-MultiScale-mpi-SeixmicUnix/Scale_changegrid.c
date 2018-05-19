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
   float **vp,**vp2;
   

   char FN1[250]={"fault_vel_280_180_ite100.dat"};
   char FN2[250]={"fault_vel_560_360_ite100.dat"};

   nx=280;
   nz=180;

   vp=alloc2float(nz,nx);

   FILE *fp;
   fp=fopen(FN1,"rb");
   for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
         fread(&vp[i][j],4L,1,fp);
   fclose(fp);
 //×××××××××××××××××××××××××××××××××××××××××××××××××××××//速度限制因子
 /*  FILE *fpc;              
   fpc=fopen(FN2,"wb");
   for(i=0;i<nx/2;i++)
      for(j=0;j<nz/2;j++)
         fwrite(&vp[i*2][j*2],4L,1,fpc);
   fclose(fpc);   */


//**************************************************//波场限制因子
  /*  FILE *fpc;              
   fpc=fopen(FN2,"wb");
   for(i=0;i<nx/2;i++)
      for(j=0;j<nz;j++)
         fwrite(&vp[i*2][j],4L,1,fpc);
   fclose(fpc);   */

 //×××××××××××××××××××××××××××××××××××××××××××××××××××××//速度延拓因子

   vp2=alloc2float(nz*2,nx*2);  
   for(i=0;i<nx;i++)
      for(j=0;j<nz;j++)
         vp2[i*2][j*2]=vp[i][j];
   for(i=0;i<nx*2-1;i++)
      for(j=0;j<nz*2-1;j++)
      {
         if((i%2!=0)&&(j%2!=0))
         {
            vp2[i][j]=(vp2[i-1][j-1]+vp2[i+1][j-1]+vp2[i-1][j+1]+vp2[i+1][j+1])/4.0;
         }
         else if((i%2!=0)&&(j%2==0))
         {
            vp2[i][j]=(vp2[i-1][j]+vp2[i+1][j])/2.0;
         }
         else if((i%2==0)&&(j%2!=0))
         {
            vp2[i][j]=(vp2[i][j-1]+vp2[i][j+1])/2.0;
         }
      }
   for(i=0;i<nx*2-1;i++)
      vp2[i][nz*2-1]=vp2[i][nz*2-2];
   for(j=0;j<nz*2-1;j++)
      vp2[nx*2-1][j]=vp2[nx*2-2][j];
   vp2[nx*2-1][nz*2-1]=vp2[nx*2-2][nz*2-2];

   FILE *fpc;              
   fpc=fopen(FN2,"wb");
   for(i=0;i<nx*2;i++)
      for(j=0;j<nz*2;j++)
         fwrite(&vp2[i][j],4L,1,fpc);
   fclose(fpc);  

   











}

