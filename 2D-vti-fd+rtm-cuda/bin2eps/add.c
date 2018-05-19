#include<stdio.h>
#include<stdlib.h>

void main(){


          char FN1[250]={"snap_e0.00d0.00_nsv_max.dat"};      
          char FN2[250]={"snap_e0.10d0.00_nsv_max.dat"};
          char FN3[250]={"snap_e0.20d0.00_nsv_max.dat"};
          char FN4[250]={"snap_e0.30d0.00_nsv_max.dat"};
          char FN5[250]={"snap_e0.40d0.00_nsv_max.dat"};
          char FN6[250]={"snap_e0.50d0.00_nsv_max.dat"};


          char FNn[250]={"snap_e0.5_max.dat"};


        int nx, nz, i;

        nx = 500;
        nz = 500;


        float *a = malloc(sizeof(float)*nx*nz);
        float *b = malloc(sizeof(float)*nx*nz);

        FILE *fp;

        fp = fopen(FN1,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);
        for(i=0;i<nx*nz;i++)
            b[i] += a[i];

        fp = fopen(FN2,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);
        for(i=0;i<nx*nz;i++)
            b[i] += a[i];

        fp = fopen(FN3,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);
        for(i=0;i<nx*nz;i++)
            b[i] += a[i];

        fp = fopen(FN4,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);
        for(i=0;i<nx*nz;i++)
            b[i] += a[i];

        fp = fopen(FN5,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);
        for(i=0;i<nx*nz;i++)
            b[i] += a[i];

        fp = fopen(FN6,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);
        for(i=0;i<nx*nz;i++)
            b[i] += a[i];

        for(i=0;i<nx*nz;i++)
            if(b[i] >=1)b[i] = 1;


        fp = fopen(FNn,"wb");
        fwrite(b,4,nx*nz,fp);
        fclose(fp);






}
