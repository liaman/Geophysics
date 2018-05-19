#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
main()
{

   void smooth(int nx,int nz,float **in,float **out,int nsm);


        int i,j,k;
        int nx,nz,nsm;
        float **value,**value_out;
 



       char FN1[250]={"fault_shot_obs_280_bigscale.dat"};
       char FN2[250]={"fault_shot_obs_280_bigscale_smooth.dat"};


      nx=11200;
     nz=3000;

     nsm=2500;



     value=alloc2float(nz,nx);   zero2float(value,nz,nx);
     value_out=alloc2float(nz,nx);   zero2float(value_out,nz,nx);



      FILE *fp1;
      fp1=fopen(FN1,"rb");
        for(i=0;i<nx;i++)
          for(j=0;j<nz;j++)
              fread(&value[i][j],FSIZE,1,fp1);
       fclose(fp1);

      
       smooth(nx,nz,value,value_out,nsm);


      fp1=fopen(FN2,"wb");
        for(i=0;i<nx;i++)
          for(j=0;j<nz;j++)
              fwrite(&value_out[i][j],FSIZE,1,fp1);
       fclose(fp1);

     free2float(value);
      free2float(value_out);

}
 void smooth(int nx,int nz,float **in,float **out,int nsm)
{

    int i,j,k;

     float **temp;

      temp=alloc2float(nz,nx);  zero2float(temp,nz,nx);

   for(i=0;i<nx;i++)
    for(j=0;j<nz;j++)
      temp[i][j]=in[i][j];

    for(i=0;i<nx;i++)
     {
      if(i%100==0)
         printf(" i==%d/%d \n",i,nx);
     for(k=0;k<nsm;k++)
      {
      for(j=0;j<nz;j++)
        { 
       // if(j==0)
          // temp[i][j]=( temp[i][j]*2 + temp[i][j] ) / 3.0;
        if(j>0&&j<(nz-1))
           temp[i][j]=( temp[i][j-1] + temp[i][j]*2 + temp[i][j+1] ) / 4.0;
        if(j==(nz-1))
           temp[i][j]=( temp[i][j-1] + temp[i][j]*2 ) / 3.0;

 
        }

      }
     }

   for(i=0;i<nx;i++)
    for(j=0;j<nz;j++)
      out[i][j]=temp[i][j];

     free2float(temp);
}

















