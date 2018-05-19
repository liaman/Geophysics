#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>

#define pi 3.141592653


/*************func**************/    
void mute_directwave(int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float vp, float epsilon, float *shot,int tt)
{
    int mu_t,mu_nt,id;
    float mu_x,mu_z,mu_t0;



   for(id = 0; id < nx*nt; id ++)
   {
        int ix=id/nt;
        int it=id%nt;

        mu_x=dx*abs(ix-fs-(is-1)*ds);
        mu_z=dz*zs;
        mu_t0=sqrt(pow(mu_x,2)+pow(mu_z,2))/(vp*sqrtf(1+2*epsilon));
        mu_t=(int)(2.0/(dt*favg));
        mu_nt=(int)(mu_t0/dt)+mu_t+tt;

           if((it>(int)(mu_t0/dt)-tt)&&(it<mu_nt))
              shot[id]=0.0;
   }
}


void main(){

	int is, it, nx, nz, nt;
	int ns, ds, fs, zs;
	float dx, dz, dt, favg;

          nx=2666;              
	   nz=1200;         favg=15;     

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=16001;    
          dt=0.0005;
     
          ns=260;       
          fs=nx/ns/2;      
          ds=nx/ns;
          zs=1;  

	   char FN1[250]={"trial/shot_obs.dat"};//shot obs
	   char FN2[250]={"trial/shot_obs_muteDir.dat"};//shot obs

        FILE *fp1 = fopen(FN1,"rb");
        FILE *fp2 = fopen(FN2,"wb");

        float *shot;

        shot = malloc(sizeof(float)*nx*nt);

        for(is = 1; is <= ns; is ++){

                printf("is = %d\n",is);

                fread(shot,sizeof(float),nx*nt,fp1);

                mute_directwave(nx,nt,dt,favg,dx,dz,fs,ds,zs,is,2000, 0.0, shot,190);

                fwrite(shot,sizeof(float),nx*nt,fp2);

        }

        fclose(fp1);
        fclose(fp2);












}
