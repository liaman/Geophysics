#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
int main()
{    



	int i,j,k,nx,nz;
	float **vp,**vp1;
/*******************************************/
/*****change the name and the parameter*****/
        char FN1[250]={"shot_all.dat"};
        char FN2[250]={"shot_dir.dat"};

        char FN3[250]={"shot_ref.dat"};
        nx=883*3;
        nz=5001; 

	 vp=alloc2float(nz,nx);	 vp1=alloc2float(nz,nx);
        zero2float(vp,nz,nx);


	FILE *fp1,*fp2;
	fp1=fopen(FN1,"rb");
	fp2=fopen(FN2,"rb");
	for(i=0;i<nx;i++)
	   for(j=0;j<nz;j++)
	   {
            fread(&vp[i][j],4,1,fp1);
            fread(&vp1[i][j],4,1,fp2);vp[i][j]-=vp1[i][j];
           }
	fclose(fp1); 
	fclose(fp2); 
            
     
/*********the velocity model ENDing**********/
/********************************************/

	fp1=fopen(FN3,"wb");
	for(i=0;i<nx;i++)
	   for(j=0;j<nz;j++)
	      fwrite(&vp[i][j],4,1,fp1);
	fclose(fp1); 

    free2float(vp);

        return 0;
}












