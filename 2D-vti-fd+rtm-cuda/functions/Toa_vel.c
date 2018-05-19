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
        char FN1[250]={"sshengli_vel_1000_550.dat"};
        char FN2[250]={"sshengli_epsilon_1000_550.dat"};

        nx=1000;
        nz=550; 

	 vp=alloc2float(nz,nx);
        zero2float(vp,nz,nx);
        layer=alloc1float(nx);
        zero1float(layer,nx);
/********************************************/
/*************the velocity model*************/
             for(k=0;k<nx;k++) 
                 {
                if(k<250) layer[k]=150;
                if(k>=250&&k<350) layer[k]=0.5*(k-250)+150;
                if(k>=350&&k<450) layer[k]=200;
                if(k>=450&&k<550) layer[k]=-0.5*(k-450)+200;
                if(k>=550) layer[k]=150;

                 }         /* nx=601 */   
             /*  for(k=0;k<nx;k++) 
                 {
                if(k<150) layer[k]=150;
                if(k>=150&&k<250) layer[k]=0.5*(k-150)+150;
                if(k>=250&&k<350) layer[k]=200;
                if(k>=350&&k<450) layer[k]=-0.5*(k-350)+200;
                if(k>=450&&k<550) layer[k]=150;

                if(k>=550&&k<650) layer[k]=0.5*(k-550)+150;
                if(k>=650&&k<750) layer[k]=200;
                if(k>=750&&k<850) layer[k]=-0.5*(k-750)+200;
                if(k>=850) layer[k]=150;

                 }*/  /* nx=1001 */ 
          /*   for(k=0;k<nx;k++) 
                 {
                if(k<150) layer[k]=150;
                if(k>=150&&k<450) layer[k]=200;
                if(k>=450) layer[k]=150;

                 }*/



	FILE *fp1;
	fp1=fopen(FN1,"rb");
	for(i=0;i<nx;i++)
	   for(j=0;j<nz;j++)
	      fread(&vp[i][j],4,1,fp1);
	fclose(fp1);                
     


	for(i=0;i<nx;i++)
	   for(j=0;j<nz;j++)
          {

            if(vp[i][j]==2000)vp[i][j]=0.0;
            if(vp[i][j]==2200)vp[i][j]=0.15;
            if(vp[i][j]==2400)vp[i][j]=0.2;
            if(vp[i][j]==2650)vp[i][j]=0.25;
            if(vp[i][j]==2900)vp[i][j]=0.3;
            if(vp[i][j]==3200)vp[i][j]=0.35;
            if(vp[i][j]==3600)vp[i][j]=0.4;


         }
/*********the velocity model ENDing**********/
/********************************************/

	fp1=fopen(FN2,"wb");
	for(i=0;i<nx;i+=1)
	   for(j=0;j<nz;j+=1)
	      fwrite(&vp[i][j],4,1,fp1);
	fclose(fp1); 

    free2float(vp);
    free1float(layer);
        return 0;
}
/************************a************************/












