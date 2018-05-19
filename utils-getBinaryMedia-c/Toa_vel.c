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
	float **vp,*layer;
/*******************************************/
/*****change the name and the parameter*****/
        char FN1[250]={"3layer_theta_601_301.dat"};
        nx=601;
        nz=301; 

	 vp=alloc2float(nz,nx);
        zero2float(vp,nz,nx);
        layer=alloc1float(nx);
        zero1float(layer,nx);
/********************************************/
/*************the velocity model*************/
             for(k=0;k<nx;k++) 
                 {
                if(k<150) layer[k]=150;
                if(k>=150&&k<250) layer[k]=0.5*(k-150)+150;
                if(k>=250&&k<350) layer[k]=200;
                if(k>=350&&k<450) layer[k]=-0.5*(k-350)+200;
                if(k>=450) layer[k]=150;

                 }         /* nx=601 */   


      for(i=0;i<nx;i++)
        {
	 for(j=0;j<nz;j++)
	  {
             vp[i][j]=30;
             if(j>=100)vp[i][j]=45;
             if(j>=200)vp[i][j]=60;
	  }
	}    

                                    
     
/*********the velocity model ENDing**********/
/********************************************/
	FILE *fp1;
	fp1=fopen(FN1,"wb");
	for(i=0;i<nx;i++)
	   for(j=0;j<nz;j++)
	      fwrite(&vp[i][j],4,1,fp1);
	fclose(fp1); 

    free2float(vp);
    free1float(layer);
        return 0;
}












