#include<stdio.h>
#include<stdlib.h>


void main()
{
  int i,is,it,ix,ns,nt,nx;
  float **sobs,**scal,**scom;


   FILE *fpsobs,*fpscal,*fpscom;

   fpsobs=fopen("vel201202203.dat","rb");
   fpscal=fopen("vel201202203initial.dat","rb");
   fpscom=fopen("com.dat","wb");


  ns=  202;
  nx=  201;
  nt=  203;


  sobs=(float **)malloc(nt*sizeof(float *));
  for(i=0;i<nt;i++)sobs[i]=(float *)malloc(nx*sizeof(float));

  scal=(float **)malloc(nt*sizeof(float *));
  for(i=0;i<nt;i++)scal[i]=(float *)malloc(nx*sizeof(float));

  scom=(float **)malloc(nt*sizeof(float *));
  for(i=0;i<nt;i++)scom[i]=(float *)malloc(nx*sizeof(float));



  for(is=0;is<ns;is++)
  {
    for(ix=0;ix<nx;ix++)
     {
      for(it=0;it<nt;it++)
       {
           fread(&sobs[ix][it],sizeof(float),1,fpsobs);
           fread(&scal[ix][it],sizeof(float),1,fpscal);

           scom[ix][it]=sobs[ix][it]-scal[ix][it];

           fwrite(&scom[ix][it],sizeof(float),1,fpscom);
       }
     }

  }

  fclose(fpsobs);
  fclose(fpscal);
  fclose(fpscom);










}






