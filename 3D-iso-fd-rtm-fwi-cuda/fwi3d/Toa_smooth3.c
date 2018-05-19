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
  float ***vel,vmin=9999999.9,vmax=-999999.9;
  int i,j,k;
  int nx,ny,nz,nn,mmz,n;
 
  char FN1[250]={"vel201202203.dat"};
  char FN2[250]={"vel201202203initial.dat"};

  FILE *fpvel;
  fpvel=fopen(FN1,"rb");
  FILE *fpsmooth;
  fpsmooth=fopen(FN2,"wb");

  nx=201;
  ny=202;
  nz=203;

  nn=50;
  mmz=10; 

  vel=alloc3float(nz,ny,nx);
  zero3float(vel,nz,ny,nx);

    for(j=0;j<ny;j++)
      {
        for(i=0;i<nx;i++)
	{ 
	  for(k=0;k<nz;k++)
	  {
	    fread(&vel[i][j][k],4L,1,fpvel);
	  }
	}
      }


/********************************************** smooth3 **********************************/
for(n=0;n<nn;n++)
{
  printf("smooth ite=%d's circle!\n",n);
   //************* inner
   for(i=1;i<nx-1;i++)
   {
      for(j=1;j<ny-1;j++)
      {
         for(k=mmz;k<nz-1;k++)
         {
            vel[i][j][k]=((vel[i-1][j-1][k-1]+vel[i-1][j][k-1]+vel[i-1][j+1][k-1]+vel[i][j-1][k-1]+vel[i][j][k-1]+vel[i][j+1][k-1]+vel[i+1][j-1][k-1]+vel[i+1][j][k-1]+vel[i+1][j+1][k-1])+(vel[i-1][j-1][k]+vel[i-1][j][k]+vel[i-1][j+1][k]+vel[i][j-1][k]+vel[i][j][k]+vel[i][j+1][k]+vel[i+1][j-1][k]+vel[i+1][j][k]+vel[i+1][j+1][k])+(vel[i-1][j-1][k+1]+vel[i-1][j][k+1]+vel[i-1][j+1][k+1]+vel[i][j-1][k+1]+vel[i][j][k+1]+vel[i][j+1][k+1]+vel[i+1][j-1][k+1]+vel[i+1][j][k+1]+vel[i+1][j+1][k+1]))/27;
         }
      }
   } 
   //************* up and down
   for(j=1;j<ny-1;j++)
   {
      for(i=1;i<nx-1;i++)
      {
         vel[i][j][mmz]=((vel[i-1][j-1][mmz]+vel[i-1][j][mmz]+vel[i-1][j+1][mmz]+vel[i][j-1][mmz]+vel[i][j][mmz]+vel[i][j+1][mmz]+vel[i+1][j-1][mmz]+vel[i+1][j][mmz]+vel[i+1][j+1][mmz])+(vel[i-1][j-1][mmz+1]+vel[i-1][j][mmz+1]+vel[i-1][j+1][mmz+1]+vel[i][j-1][mmz+1]+vel[i][j][mmz+1]+vel[i][j+1][mmz+1]+vel[i+1][j-1][mmz+1]+vel[i+1][j][mmz+1]+vel[i+1][j+1][mmz+1]))/18;
         vel[i][j][nz-1]=((vel[i-1][j-1][nz-1]+vel[i-1][j][nz-1]+vel[i-1][j+1][nz-1]+vel[i][j-1][nz-1]+vel[i][j][nz-1]+vel[i][j+1][nz-1]+vel[i+1][j-1][nz-1]+vel[i+1][j][nz-1]+vel[i+1][j+1][nz-1])+(vel[i-1][j-1][nz-2]+vel[i-1][j][nz-2]+vel[i-1][j+1][nz-2]+vel[i][j-1][nz-2]+vel[i][j][nz-2]+vel[i][j+1][nz-2]+vel[i+1][j-1][nz-2]+vel[i+1][j][nz-2]+vel[i+1][j+1][nz-2]))/18;
      }
   }
   //************* left and right
   for(j=1;j<ny-1;j++)
   {
      for(k=mmz;k<nz-1;k++)
      {
         vel[0][j][k]=((vel[0][j-1][k-1]+vel[0][j-1][k]+vel[0][j-1][k+1]+vel[0][j][k-1]+vel[0][j][k]+vel[0][j][k+1]+vel[0][j+1][k-1]+vel[0][j+1][k]+vel[0][j+1][k+1])+(vel[1][j-1][k-1]+vel[1][j-1][k]+vel[1][j-1][k+1]+vel[1][j][k-1]+vel[1][j][k]+vel[1][j][k+1]+vel[1][j+1][k-1]+vel[1][j+1][k]+vel[1][j+1][k+1]))/18;
         vel[nx-1][j][k]=((vel[nx-1][j-1][k-1]+vel[nx-1][j-1][k]+vel[nx-1][j-1][k+1]+vel[nx-1][j][k-1]+vel[nx-1][j][k]+vel[nx-1][j][k+1]+vel[nx-1][j+1][k-1]+vel[nx-1][j+1][k]+vel[nx-1][j+1][k+1])+(vel[nx-2][j-1][k-1]+vel[nx-2][j-1][k]+vel[nx-2][j-1][k+1]+vel[nx-2][j][k-1]+vel[nx-2][j][k]+vel[nx-2][j][k+1]+vel[nx-2][j+1][k-1]+vel[nx-2][j+1][k]+vel[nx-2][j+1][k+1]))/18;
      } 
   }
   //************** front and back
   for(i=1;i<nx-1;i++)
   {
      for(k=mmz;k<nz;k++)
      {
         vel[i][0][k]=((vel[i-1][0][k-1]+vel[i-1][0][k]+vel[i-1][0][k+1]+vel[i][0][k-1]+vel[i][0][k]+vel[i][0][k+1]+vel[i+1][0][k-1]+vel[i+1][0][k]+vel[i+1][0][k+1])+(vel[i-1][1][k-1]+vel[i-1][1][k]+vel[i-1][1][k+1]+vel[i][1][k-1]+vel[i][1][k]+vel[i][1][k+1]+vel[i+1][1][k-1]+vel[i+1][1][k]+vel[i+1][1][k+1]))/18;
         vel[i][ny-1][k]=((vel[i-1][ny-1][k-1]+vel[i-1][ny-1][k]+vel[i-1][ny-1][k+1]+vel[i][ny-1][k-1]+vel[i][ny-1][k]+vel[i][ny-1][k+1]+vel[i+1][ny-1][k-1]+vel[i+1][ny-1][k]+vel[i+1][ny-1][k+1])+(vel[i-1][ny-2][k-1]+vel[i-1][ny-2][k]+vel[i-1][ny-2][k+1]+vel[i][ny-2][k-1]+vel[i][ny-2][k]+vel[i][ny-2][k+1]+vel[i+1][ny-2][k-1]+vel[i+1][ny-2][k]+vel[i+1][ny-2][k+1]))/18;
      }
   }
   
   //********* the eight point
   vel[0][0][mmz]=(vel[1][0][mmz]+vel[0][1][mmz]+vel[1][1][mmz]+vel[0][0][mmz]+vel[1][0][mmz+1]+vel[0][1][mmz+1]+vel[1][1][mmz+1]+vel[0][0][mmz+1])/8;

   vel[0][0][nz-1]=(vel[1][0][nz-1]+vel[0][1][nz-1]+vel[1][1][nz-1]+vel[0][0][nz-1]+vel[1][0][nz-2]+vel[0][1][nz-2]+vel[1][1][nz-2]+vel[0][0][nz-2])/8;

   vel[nx-1][0][mmz]=(vel[nx-2][0][mmz]+vel[nx-1][1][mmz]+vel[nx-2][1][mmz]+vel[nx-1][0][mmz]+vel[nx-2][0][mmz+1]+vel[nx-1][1][mmz+1]+vel[nx-2][1][mmz+1]+vel[nx-1][0][mmz+1])/8;

   vel[nx-1][0][nz-1]=(vel[nx-2][0][nz-1]+vel[nx-1][1][nz-1]+vel[nx-2][1][nz-1]+vel[nx-1][0][nz-1]+vel[nx-2][0][nz-2]+vel[nx-1][1][nz-2]+vel[nx-2][1][nz-2]+vel[nx-1][0][nz-2])/8;

   vel[0][ny-1][mmz]=(vel[1][ny-1][mmz]+vel[0][ny-2][mmz]+vel[1][ny-2][mmz]+vel[0][ny-1][mmz]+vel[1][ny-1][mmz+1]+vel[0][ny-2][mmz+1]+vel[1][ny-2][mmz+1]+vel[0][ny-1][mmz+1])/8;

   vel[0][ny-1][nz-1]=(vel[1][ny-1][nz-1]+vel[0][ny-2][nz-1]+vel[1][ny-2][nz-1]+vel[0][ny-1][nz-1]+vel[1][ny-1][nz-2]+vel[0][ny-2][nz-2]+vel[0][ny-1][nz-2]+vel[1][ny-2][nz-2])/8;

   vel[nx-1][ny-1][mmz]=(vel[nx-2][ny-1][mmz]+vel[nx-1][ny-2][mmz]+vel[nx-2][ny-2][mmz]+vel[nx-1][ny-1][mmz]+vel[nx-2][ny-1][mmz+1]+vel[nx-1][ny-2][mmz+1]+vel[nx-2][ny-2][mmz+1]+vel[nx-1][ny-1][mmz+1])/8;

   vel[nx-1][ny-1][nz-1]=(vel[nx-2][ny-1][nz-1]+vel[nx-1][ny-2][nz-1]+vel[nx-2][ny-2][nz-1]+vel[nx-1][ny-1][nz-1]+vel[nx-2][ny-1][nz-2]+vel[nx-1][ny-2][nz-2]+vel[nx-2][ny-2][nz-2]+vel[nx-1][ny-1][nz-2])/8;

   
   for(k=mmz;k<nz-1;k++)
   {
      //*************  x=0 & y=0
      vel[0][0][k]=((vel[1][0][k]+vel[0][1][k]+vel[1][1][k]+vel[0][0][k])+(vel[1][0][k-1]+vel[0][1][k-1]+vel[1][1][k-1]+vel[0][0][k-1])+(vel[1][0][k+1]+vel[0][1][k+1]+vel[1][1][k+1]+vel[0][0][k+1]))/12;
      //*************  x=0 & y=ny-1
      vel[0][ny-1][k]=((vel[1][ny-1][k]+vel[0][ny-2][k]+vel[1][ny-2][k]+vel[0][ny-1][k])+(vel[1][ny-1][k-1]+vel[0][ny-2][k-1]+vel[1][ny-2][k-1]+vel[0][ny-1][k-1])+(vel[1][ny-1][k+1]+vel[0][ny-2][k+1]+vel[1][ny-2][k+1]+vel[0][ny-1][k+1]))/12;
      //*************  x=nx-1 & y=0
      vel[nx-1][0][k]=((vel[nx-2][0][k]+vel[nx-1][1][k]+vel[nx-2][1][k]+vel[nx-1][0][k])+(vel[nx-2][0][k-1]+vel[nx-1][1][k-1]+vel[nx-2][1][k-1]+vel[nx-1][0][k-1])+(vel[nx-2][0][k+1]+vel[nx-1][1][k+1]+vel[nx-2][1][k+1]+vel[nx-1][0][k+1]))/12;
      //*************  x=nx-1 & y=ny-1
      vel[nx-1][ny-1][k]=((vel[nx-2][ny-1][k]+vel[nx-1][ny-2][k]+vel[nx-2][ny-2][k]+vel[nx-1][ny-1][k])+(vel[nx-2][ny-1][k-1]+vel[nx-1][ny-2][k-1]+vel[nx-2][ny-2][k-1]+vel[nx-1][ny-1][k-1])+(vel[nx-2][ny-1][k+1]+vel[nx-1][ny-2][k+1]+vel[nx-2][ny-2][k+1]+vel[nx-1][ny-1][k+1]))/12;
   } 
   for(i=1;i<nx-1;i++)
   {
      //************** y=0 & z=mmz
      vel[i][0][mmz]=((vel[i][1][mmz]+vel[i][0][mmz+1]+vel[i][1][mmz+1]+vel[i][0][mmz])+(vel[i-1][1][mmz]+vel[i-1][0][mmz+1]+vel[i-1][1][mmz+1]+vel[i-1][0][mmz])+(vel[i+1][1][mmz]+vel[i+1][0][mmz+1]+vel[i+1][1][mmz+1]+vel[i+1][0][mmz]))/12;
      //************** y=0 & z=nz-1
      vel[i][0][nz-1]=((vel[i][1][nz-1]+vel[i][0][nz-2]+vel[i][1][nz-2]+vel[i][0][nz-1])+(vel[i-1][1][nz-1]+vel[i-1][0][nz-2]+vel[i-1][1][nz-2]+vel[i-1][0][nz-1])+(vel[i+1][1][nz-1]+vel[i+1][0][nz-2]+vel[i+1][1][nz-2]+vel[i+1][0][nz-1]))/12;
      //************** y=ny-1 & z=mmz
      vel[i][ny-1][mmz]=((vel[i][ny-2][mmz]+vel[i][ny-1][mmz+1]+vel[i][ny-2][mmz+1]+vel[i][ny-1][mmz])+(vel[i-1][ny-2][mmz]+vel[i-1][ny-1][mmz+1]+vel[i-1][ny-2][mmz+1]+vel[i-1][ny-1][mmz])+(vel[i+1][ny-2][mmz]+vel[i+1][ny-1][mmz+1]+vel[i+1][ny-2][mmz+1]+vel[i+1][ny-1][mmz]))/12;
      //************** y=ny-1 & z=nz-1
      vel[i][ny-1][nz-1]=((vel[i][ny-2][nz-1]+vel[i][ny-1][nz-2]+vel[i][ny-2][nz-2]+vel[i][ny-1][nz-1])+(vel[i-1][ny-2][nz-1]+vel[i-1][ny-1][nz-2]+vel[i-1][ny-2][nz-2]+vel[i-1][ny-1][nz-1])+(vel[i+1][ny-2][nz-1]+vel[i+1][ny-1][nz-2]+vel[i+1][ny-2][nz-2]+vel[i+1][ny-1][nz-1]))/12;
   }
   for(j=1;j<ny-1;j++)
   {
      //************** x=0 & z=mmz
      vel[0][j][mmz]=((vel[1][j][mmz]+vel[0][j][mmz+1]+vel[1][j][mmz+1]+vel[0][j][mmz])+(vel[1][j-1][mmz]+vel[0][j-1][mmz+1]+vel[1][j-1][mmz+1]+vel[0][j-1][mmz])+(vel[1][j+1][mmz]+vel[0][j+1][mmz+1]+vel[1][j+1][mmz+1]+vel[0][j+1][mmz]))/12;
      //************** x=0 & z=nz-1
      vel[0][j][nz-1]=((vel[1][j][nz-1]+vel[0][j][nz-2]+vel[1][j][nz-2]+vel[0][j][nz-1])+(vel[1][j-1][nz-1]+vel[0][j-1][nz-2]+vel[1][j-1][nz-2]+vel[0][j-1][nz-1])+(vel[1][j+1][nz-1]+vel[0][j+1][nz-2]+vel[1][j+1][nz-2]+vel[0][j+1][nz-1]))/12;
      //************** x=nx-1 & z=mmz
      vel[nx-1][j][mmz]=((vel[nx-2][j][mmz]+vel[nx-1][j][mmz+1]+vel[nx-2][j][mmz+1]+vel[nx-1][j][mmz])+(vel[nx-2][j-1][mmz]+vel[nx-1][j-1][mmz+1]+vel[nx-2][j-1][mmz+1]+vel[nx-1][j-1][mmz])+(vel[nx-2][j+1][mmz]+vel[nx-1][j+1][mmz+1]+vel[nx-2][j+1][mmz+1]+vel[nx-1][j+1][mmz]))/12;
      //************** x=nx-1 &z=nz-1
      vel[nx-1][j][nz-1]=((vel[nx-2][j][nz-1]+vel[nx-1][j][nz-2]+vel[nx-2][j][nz-2]+vel[nx-1][j][nz-1])+(vel[nx-2][j-1][nz-1]+vel[nx-1][j-1][nz-2]+vel[nx-2][j-1][nz-2]+vel[nx-1][j-1][nz-1])+(vel[nx-2][j+1][nz-1]+vel[nx-1][j+1][nz-2]+vel[nx-2][j+1][nz-2]+vel[nx-1][j+1][nz-1]))/12;
   }



}








/*************************************************************************************************/ 
      for(j=0;j<ny;j++)
      {
        for(i=0;i<nx;i++)
	{ 
	  for(k=0;k<nz;k++)
	  {
            if(vel[i][j][k]<vmin)vmin=vel[i][j][k];
            if(vel[i][j][k]>vmax)vmax=vel[i][j][k];
            vel[i][j][k]=vel[i][j][k];
	    fwrite(&vel[i][j][k],4L,1,fpsmooth);
	  }
	}
      }
printf("vmin=%.2f, vmax=%.2f \n",vmin,vmax);
  fclose(fpvel);
  free3float(vel);
}








