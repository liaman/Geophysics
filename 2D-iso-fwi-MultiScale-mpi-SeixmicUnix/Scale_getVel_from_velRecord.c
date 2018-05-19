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



        int i,j,k;
        int nx,nz,ite;
        float **value;
 



       char FN1[250]={"Velocity_IteRecord.dat"};
       char FN2[250]={"fault_vel_560_360_ite35.dat"};


      nx=560;
     nz=360;

     ite=35;



     value=alloc2float(nz,nx);   zero2float(value,nz,nx);



      FILE *fp1;
      fp1=fopen(FN1,"rb");
       fseek(fp1,(ite-1)*nx*nz*FSIZE,0);
        for(i=0;i<nx;i++)
          for(j=0;j<nz;j++)
              fread(&value[i][j],FSIZE,1,fp1);
       fclose(fp1);

      

      fp1=fopen(FN2,"wb");
        for(i=0;i<nx;i++)
          for(j=0;j<nz;j++)
              fwrite(&value[i][j],FSIZE,1,fp1);
       fclose(fp1);

     free2float(value);

}

















