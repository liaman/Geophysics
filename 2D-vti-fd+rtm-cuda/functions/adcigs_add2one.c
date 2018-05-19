#include<stdio.h>
#include<stdlib.h>

void main(){


    int nx, nz, na, ix, iz, ia, i;

    char FN1[100] = {"../layers_update_layer0_v0.80e1.00d1.00_adcigs.dat"};
    char FN2[100] = {"../layers_update_layer0_v0.80e1.00d1.00_adcigs_add.dat"};


    nx = 2500;
    nz = 2500;
    na = 60;

    float *in = (float*)malloc(sizeof(float)*nz*na);

    float *tmp = (float*)malloc(sizeof(float)*nz*na);

    FILE *fp = fopen(FN1,"rb");

    for(ix=0;ix<nx;ix++){

        printf("ix = %d\n",ix);

        fread(in,sizeof(float),na*nz,fp);
        for(i=0;i<nz*na;i++){

            tmp[i] += (in[i]);
        }

    }
    fclose(fp);
    fp = fopen(FN2,"wb");
    fwrite(tmp,sizeof(float),na*nz,fp);
    fclose(fp);


}
