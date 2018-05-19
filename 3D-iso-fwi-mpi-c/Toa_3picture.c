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
  int nx,ny,nz,px,py,pz;
 
  char FN1[250]={"multi_grad_zm.dat"};

  char FN2[250]={"a2.dat"};
  char FN3[250]={"a3.dat"};
  char FN4[250]={"a4.dat"};
  char FN5[250]={"a5.dat"};
  char FN6[250]={"a6.dat"};
  char FN7[250]={"a7.dat"};
  char FN8[250]={"a8.dat"};
  char FN9[250]={"a9.dat"};

  FILE *fpvel;
  fpvel=fopen(FN1,"rb");

  nx=200;   px=100;
  ny=200;   py=100;
  nz=200;   pz=100;

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
  fclose(fpvel);


  fpvel=fopen(FN2,"wb");
    for(j=0;j<py;j++)
     for(i=0;i<px;i++)
       for(k=0;k<pz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN3,"wb");
    for(j=0;j<py;j++)
     for(i=px;i<nx;i++)
       for(k=0;k<pz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN4,"wb");
    for(j=py;j<ny;j++)
     for(i=px;i<nx;i++)
       for(k=0;k<pz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN5,"wb");
    for(j=py;j<ny;j++)
     for(i=0;i<px;i++)
       for(k=0;k<pz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN6,"wb");
    for(j=0;j<py;j++)
     for(i=0;i<px;i++)
       for(k=pz;k<nz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN7,"wb");
    for(j=0;j<py;j++)
     for(i=px;i<nx;i++)
       for(k=pz;k<nz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN8,"wb");
    for(j=py;j<ny;j++)
     for(i=px;i<nx;i++)
       for(k=pz;k<nz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);

  fpvel=fopen(FN9,"wb");
    for(j=py;j<ny;j++)
     for(i=0;i<px;i++)
       for(k=pz;k<nz;k++)
	    fwrite(&vel[i][j][k],4L,1,fpvel);
  fclose(fpvel);


     
  free3float(vel);printf("done\n");
}








