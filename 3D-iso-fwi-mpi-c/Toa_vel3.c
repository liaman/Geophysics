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
  float ***vel;
  int i,j,k;
  int nx,ny,nz;
 
  char FN1[250]={"vel200200200.dat"};
  FILE *fpvel;
  fpvel=fopen(FN1,"wb");

  nx=200;
  ny=200;
  nz=200; 

  vel=alloc3float(nz,ny,nx);
  zero3float(vel,nz,ny,nx);

  for(j=0;j<ny;j++)
  {
    for(i=0;i<nx;i++)
    {
       for(k=0;k<nz;k++)
       {
            vel[i][j][k]=3000;
             if(k>40)
               vel[i][j][k]=3200;
             if(k>70)
               vel[i][j][k]=3700;
             if(k>100)
               vel[i][j][k]=4000;
             if(k>140)
               vel[i][j][k]=4200;
             if(k>180)
               vel[i][j][k]=4500;
          /*       if(k>80)
               vel[i][j][k]=3200;
               if(k>120)
               vel[i][j][k]=3300;
               if(k>160)
               vel[i][j][k]=3400;
            if((k-(0.15)*(i+j))>55)
               vel[i][j][k]=3300;
            if((k-(0.20)*(-1*i+(200-0.5*j)))>125)
               vel[i][j][k]=3400;  
            if((k-(0.12)*(1.2*i+0.7*j))>150)
               vel[i][j][k]=3500; */

           //if((i>23&&i<27)&&(j>23&&j<27)&&(k>23&&k<27)) vel[i][j][k]=4000;
       }
    }
  }

 /*  for(i=0;i<nx;i++)
  {
    for(j=0;j<ny;j++)
    {
       for(k=40;k<nz;k++)
       {
            vel[i][j][k]=3000;
       }
    }
  }
 */
 // printf("============%f\n",vel[nx-1][ny-1][nz-1]);
  /* for(i=0;i<nx;i++)
     for(j=0;j<ny;j++)
       for(k=0;k<nz;k++)
       {
         if((pow((i-nx/4),2)+pow((j-ny/4),2)+pow((k-nz/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx*3/4),2)+pow((j-ny/4),2)+pow((k-nz/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx/4),2)+pow((j-ny*3/4),2)+pow((k-nz/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx/4),2)+pow((j-ny/4),2)+pow((k-nz*3/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx*3/4),2)+pow((j-ny*3/4),2)+pow((k-nz/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx/4),2)+pow((j-ny*3/4),2)+pow((k-nz*3/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx*3/4),2)+pow((j-ny/4),2)+pow((k-nz*3/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx*3/4),2)+pow((j-ny*3/4),2)+pow((k-nz*3/4),2))<=400)
           { vel[i][j][k]=3500;}
         if((pow((i-nx*2/4),2)+pow((j-ny*2/4),2)+pow((k-nz*2/4),2))<=400)
           { vel[i][j][k]=3500;}
       }    */

 /*     FILE *fp555;
      fp555=fopen("vel.txt","a");
      
      for(i=0;i<nx;i++)
        for(j=0;j<ny;j++)
          for(k=0;k<nz;k++)
            { fprintf(fp555,"vel[%2d][%2d][%2d] = %f \n",i,j,k,vel[i][j][k]);
              if(vel[i][j][k]==(i+j+k))
                 fprintf(fp555,"***\n");
            }
      fclose(fp555);
 */
      for(j=0;j<ny;j++)
      {
        for(i=0;i<nx;i++)
	{ 
	  for(k=0;k<nz;k++)
	  {
	    fwrite(&vel[i][j][k],4L,1,fpvel);
	  }
	}
      }
  fclose(fpvel);
  free3float(vel);printf("done\n");
}








