#include<stdio.h>
#include<stdlib.h>

void main(){


        int ix, iz, id;
        int nx, nz;

        char FNv[100] = {"../fd/layer_v_2000_2.dat"};
        char FNe[100] = {"../fd/layer_e_0.5.dat"};
        char FNd[100] = {"../fd/layer_d_0.5.dat"};


        nx = 2000;
        nz = 210;

        float *v, *e, *d;
        int *wlayer;

        v = malloc(sizeof(float)*nx*nz);
        e = malloc(sizeof(float)*nx*nz);
        d = malloc(sizeof(float)*nx*nz);

        wlayer = malloc(sizeof(int)*nx);

        for(ix = 0; ix < nx; ix ++){

            if(ix < 400)
                wlayer[ix] = 800;
            else if(ix >=400 && ix < 800)
                wlayer[ix] = 800+(ix-400);
            else if(ix >=800 && ix < 1200)
                wlayer[ix] = 1200;
            else if(ix >=1200 && ix < 1600)
                wlayer[ix] = 1200 -(ix - 1200);
            else 
                wlayer[ix] = 800;
        }

        for(ix = 0; ix < nx; ix ++){

                for(iz = 0; iz < nz; iz ++){

                        id = ix*nz + iz;

                        v[id] = 2000.0;
                        e[id] = 0.5;
                        d[id] = 0.5;

                        if(iz > 200){

                                v[id] = 2500;
                                //e[id] = 0.5;
                                //d[id] = 0.5;
                        }
                   /*     if(iz > wlayer[ix]/2){

                                v[id] = 3600;
                                e[id] = 0.2;
                                d[id] = 0.15;
                        }

                        if(iz > 800){

                                v[id] = 4000;
                                e[id] = 0.3;
                                d[id] = 0.2;
                        } */
                     /*   if(iz > 2000){

                                v[id] = 4500;
                                e[id] = 0.25;
                                d[id] = 0.1;
                        }*/

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
