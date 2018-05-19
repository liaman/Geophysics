#include<stdio.h>
#include<stdlib.h>

void main()
{

  int ix,iy,iz,nx,ny,nz;
  float val1,val2;

  FILE *in1,*in2,*out;

  in1=fopen("snap_front013.dat","rb");
  in2=fopen("snap_back013.dat","rb");
 out=fopen("poynting_y013.dat","wb");


  nx=nz=ny=200;

  for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
      for(iz=0;iz<nz;iz++)
       {

         fread(&val1,4L,1,in1);
         fread(&val2,4L,1,in2);
          val1=val1+val2;
         fwrite(&val1,4L,1,out);







       }


  fclose(in1);
  fclose(in2);
  fclose(out);





















}
