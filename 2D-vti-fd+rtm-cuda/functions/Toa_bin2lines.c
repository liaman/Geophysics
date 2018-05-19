#include<stdio.h>
#include<stdlib.h>

void main(){

    int nx, nz, ix, iz, id;

    char FN1[100] = {"v4_2500_2500.dat"};
    char FN2[100] = {"linev4.txt"};

    nx = 2500;
    nz = 2500;

    float *v = malloc(sizeof(float)*nx*nz);

    FILE *fp = fopen(FN1,"rb");

    fread(v, sizeof(float), nx*nz, fp);

    fclose(fp);

    fp = fopen(FN2, "w");

    for(iz = 0; iz < nz; iz ++) {

        fprintf(fp, "%f\n",v[iz + nz * 1250]);
    }
    fclose(fp);

}
