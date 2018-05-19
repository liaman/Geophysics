

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>



int main(int argc, char* argv[])
{
    int nz,nx,dcdp,na,ix,iz,ia,id,ido,nna,i,j;
    float *in;

    char FN1[250]={"../trial/waxian_adcigs.dat"};
    char FN2[250]={"../trial/waxian_adcigs_chouxi.dat"};




    nx=2000;
    nz=1000;
    na=70;



    in = (float*)malloc(nz*na*sizeof(float));


    FILE *fpin, *fpout;
    fpin  = fopen (FN1,"rb");
    fpout = fopen (FN2,"wb");

    
    fseek(fpin,sizeof(float)*nz*na*300,0);
    fread(in,sizeof(float),nz*na,fpin);
    fwrite(in,sizeof(float),na*nz,fpout);

    fseek(fpin,sizeof(float)*nz*na*600,0);
    fread(in,sizeof(float),nz*na,fpin);
    fwrite(in,sizeof(float),na*nz,fpout);

    fseek(fpin,sizeof(float)*nz*na*1000,0);
    fread(in,sizeof(float),nz*na,fpin);
    fwrite(in,sizeof(float),na*nz,fpout);

    fseek(fpin,sizeof(float)*nz*na*1400,0);
    fread(in,sizeof(float),nz*na,fpin);
    fwrite(in,sizeof(float),na*nz,fpout);

    fseek(fpin,sizeof(float)*nz*na*1700,0);
    fread(in,sizeof(float),nz*na,fpin);
    fwrite(in,sizeof(float),na*nz,fpout);



    fclose(fpin);
    fclose(fpout);






}
