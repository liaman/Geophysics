

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>



int main(int argc, char* argv[])
{
    int nz,nx,dcdp,na,ix,iz,ia,id,ido,nna;
    float *in,*out;

    char FN1[250]={"../trial/waxian_0_adcigs.dat"};
    char FN2[250]={"../trial/waxian_0_adcigs_chouxi.dat"};

    FILE *fp1, *fp2;
    fp1  = fopen (FN1,"rb");
    fp2 = fopen (FN2,"wb");


    nx=2000;
    nz=1000;
    na=70;
    dcdp=500;

    nna=70;


     in = (float*)malloc(nx*nz*na*sizeof(float));printf("----\n");
    out = (float*)malloc(nx/dcdp*nz*nna*sizeof(float));printf("----\n");

    fread(in,sizeof(float),nx*nz*na,fp1);printf("----2\n");
    for (ix=0; ix<nx; ix++) 
    {
           for (ia=0; ia<na; ia++) 
            {
               for (iz=0; iz<nz; iz++) 
                 {
                      id=ix*na*nz+ia*nz+iz;
                      if(ix%dcdp==0&&ia<nna)
                        {
                            ido=ix/dcdp*nna*nz+ia*nz+iz;
                              out[ido]=in[id];
                        }
                  }
              }
    }
    fwrite(out,sizeof(float),nx/dcdp*nz*nna,fp2);
    exit(0);



}
