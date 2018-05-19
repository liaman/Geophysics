#include<stdio.h>
#include<stdlib.h>

void main(){


        int ix, iz, id;
        int nx, nz;

        char FNv[100] = {"vti_v_10000_2500.dat"};
        char FNe[100] = {"vti_e_10000_2500.dat"};
        char FNd[100] = {"vti_d_10000_2500.dat"};


        nx = 10000;
        nz = 2500;

        float *v, *e, *d;

        v = malloc(sizeof(float)*nx*nz);
        e = malloc(sizeof(float)*nx*nz);
        d = malloc(sizeof(float)*nx*nz);

        for(ix = 0; ix < nx; ix ++){

                for(iz = 0; iz < nz; iz ++){

                        id = ix*nz + iz;

                        v[id] = 2000.0;
                        e[id] = 0.0;
                        d[id] = 0.0;

                        if(iz > 500){

                                v[id] = 2600;
                                e[id] = 0.07;
                                d[id] = 0.02;
                        }
                        if(iz > 1000){

                                v[id] = 3300;
                                e[id] = 0.12;
                                d[id] = 0.05;
                        }

                        if(iz > 1500){

                                v[id] = 3800;
                                e[id] = 0.19;
                                d[id] = 0.08;
                        }
                        if(iz > 2000){

                                v[id] = 4500;
                                e[id] = 0.25;
                                d[id] = 0.1;
                        }

                }
        }

        FILE *fpv, *fpe, *fpd;
        fpv = fopen(FNv,"wb");
        fpe = fopen(FNe,"wb");
        fpd = fopen(FNd,"wb");

        fwrite(v, nx*nz, sizeof(float), fpv);
        fwrite(e, nx*nz, sizeof(float), fpe);
        fwrite(d, nx*nz, sizeof(float), fpd);


        fclose(fpv);
        fclose(fpe);
        fclose(fpd);


        free(v);
        free(e);
        free(d);


}
