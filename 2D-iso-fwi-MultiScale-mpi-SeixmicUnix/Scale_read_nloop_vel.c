#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "/home/rongtao/gpfs03/hc/cjbsegy.h"
#include "/home/rongtao/gpfs03/hc/fft.c"
#include "/home/rongtao/gpfs03/hc/alloc.c"
#include "/home/rongtao/gpfs03/hc/complex.c"
main()
{
    float **vp;
    int nx,nz,nloop,i,j,k;
    
    char FN1[250]={"record_velocity.dat"};
    char FN2[250]={"vel_nloop.dat"};

    nx=200;
    nz=180;
    nloop=140;
    scanf("%d",&nloop);
    vp=alloc2float(nz,nx);
    zero2float(vp,nz,nx);

    FILE *fp;
    fp=fopen(FN1,"rb");
    fseek(fp,(nloop-1)*nx*nz*4L,0);
    for(i=0;i<nx;i++)
       for(j=0;j<nz;j++)
          fread(&vp[i][j],4L,1,fp);
    fclose(fp);


    fp=fopen(FN2,"wb");
    for(i=0;i<nx;i++)
       for(j=0;j<nz;j++)
          fwrite(&vp[i][j],4L,1,fp);
    fclose(fp);

}
