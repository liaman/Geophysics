#include<stdio.h>
#include<stdlib.h>

void main()
{

  int ix,iy,iz,nx,ny,nz;
  float val;

  FILE *in,*out;

  in=fopen("poynting_y013.dat","rb");
 out=fopen("a_snap_new.dat","wb");


  nx=nz=ny=200;

  for(iy=0;iy<ny;iy++)
    for(ix=0;ix<nx;ix++)
      for(iz=0;iz<nz;iz++)
       {

         fread(&val,4L,1,in);

         if(ix<100)fwrite(&val,4L,1,out);







       }


  fclose(in);
  fclose(out);





















}
