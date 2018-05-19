#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
void main()
{
  float ***vp,*layer,val;
  int i,j,k,ij,l,lay;
  int nx,ny,nz;
 
  char FN1[250]={"model_del_201201201.dat"};
  FILE *fpvel;
  fpvel=fopen(FN1,"wb");

  nx=201;
  ny=201;
  nz=201; 

       vp=alloc3float(nz,ny,nx);       zero3float(vp,nz,ny,nx);


       layer=alloc1float(nx);
       zero1float(layer,nx);
/************************ waxian   ******************************/
             for(k=0;k<nx;k++) 
                 {
                if(k<25) layer[k]=100;
                if(k>=25&&k<75) layer[k]=0.5*(k-25)+100;
                if(k>=75&&k<125) layer[k]=125;
                if(k>=125&&k<175) layer[k]=-0.5*(k-125)+125;
                if(k>=175) layer[k]=100;

                 }
  for(j=0;j<ny;j++)
  {
    for(i=0;i<nx;i++)
    {
       for(k=0;k<nz;k++)
       {
             vp[i][j][k]=0.0;//3000.0;//3000;
          /*   if(k>30)vp[i][j][k]=0.05;//3200;
             if(k>60)vp[i][j][k]=0.1;//3400;
             if(k>layer[i])vp[i][j][k]=0.15;//3600;
             if(k>160)vp[i][j][k]=0.25;//3800; */
                if(k>60)vp[i][j][k]=0.15;//3400.0;
                if(k>120)vp[i][j][k]=0.3;//3800.0;

       }
    }
  }

/************************ abnormal cube   ******************************/
 /* l=8;

  for(j=0;j<ny;j++)
  {
    for(i=0;i<nx;i++)
    {
       for(k=0;k<nz;k++)
       {
             vp[i][j][k]=0.0;//3000;
             if(k>80)vp[i][j][k]=0.15;//0.2;//3200;
             if((i>nx/2-l&&i<nx/2+l)&&(j>ny/2-l&&j<ny/2+l)&&(k>nz*3/5-l&&k<nz*3/5+l))vp[i][j][k]=0.2;//0.3;//4000;

       }
    }
  }*/
/************************ thrust   ******************************/



/*

  for(j=0;j<ny;j++)
  {
    for(i=0;i<nx;i++)
    {
       for(k=0;k<nz;k++)
       {
             vp[i][j][k]=0.0;//3000;
             if(k>60)vp[i][j][k]=0.0;//3140;
             if(k>180)vp[i][j][k]=0.0;//3500;
       }
    }
  }

   val=0.081;//3325.0;

for(ij=0;ij<nx+ny;ij++)
  for(j=0;j<ny;j++)
  {
    for(i=0;i<nx;i++)
    {
       if(i+j==ij)
         //if(((i-j)*(i-j)<10000))
       for(k=0;k<nz;k++)
       {
          if(ij<=100&&k>140&&k<180)vp[i][j][k]=val;

          if(ij>=100&&ij<=150&&k>=(140-0.2*(ij-100))&&k<=180)vp[i][j][k]=val;

          if(ij>=150&&ij<=200&&k>=(140-0.2*(ij-100))&&k<=(180-0.2*(ij-150)))vp[i][j][k]=val;

          if(ij>=200&&ij<=250&&k>=(120-0.4*(ij-200))&&k<=(180-0.2*(ij-150)))vp[i][j][k]=val;

          if(ij>=250&&ij<=300&&k>=(100-0.4*(ij-250))&&k<=(160-0.4*(ij-250)))vp[i][j][k]=val;

          if(ij>=300&&ij<=320&&k>=(80-0.8*(ij-300))&&k<=(160-0.4*(ij-250)))vp[i][j][k]=val;

          if(ij>=320&&ij<=350&&k>=(80-0.8*(ij-300))&&k<=(132-0.4*(ij-320)))vp[i][j][k]=val;

          if(ij>=350&&ij<=430&&k>=40&&k<=(120-0.8*(ij-350)))vp[i][j][k]=val;

          if(ij>=430&&ij<=500&&k>=40&&k<=(56-2*(ij-430)))vp[i][j][k]=val;








       }
    }
  }
*/
/************************end*****************************/









      for(j=0;j<ny;j++)
      {
        for(i=0;i<nx;i++)
	{   //if(i==j)
	  for(k=0;k<nz;k++)
	  {
        //  vp[i][j][k]=0.3;
	    fwrite(&vp[i][j][k],4L,1,fpvel);
	  }
	}
      }
  fclose(fpvel);
  free3float(vp);printf("done\n");
  exit(0);
}








