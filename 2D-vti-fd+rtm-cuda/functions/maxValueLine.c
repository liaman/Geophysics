#include<stdio.h>
#include<stdlib.h>

void main(int argc, char *argv[]){


          char *FN1=argv[1];//{"../fd/layer_shot_v2000e0.0d0.0ref.dat"};   


          char *FNn=argv[2];//{"../fd/layer_shot_v2000e0.0d0.0refmaxValue.dat"};


        int nx, nz, i;

        nx = 2000;
        nz = 11001;


        float *a = malloc(sizeof(float)*nx*nz);
        float *b = malloc(sizeof(float)*nx*nz);
        float max = 9999.9;

        FILE *fp;

        fp = fopen(FN1,"rb");
        fread(a,4,nx*nz,fp);
        fclose(fp);

        for(i=0;i<nx*nz;i++)
           if(max > a[i])
                max = a[i];

        for(i=0;i<nx*nz;i++)
           if(a[i] < max*0.16){
                b[i] = 1;
                printf("%d\n",i);
            }


        fp = fopen(FNn,"wb");
        fwrite(b,4,nx*nz,fp);
        fclose(fp);






}
