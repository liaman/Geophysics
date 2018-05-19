//#########################################################
//##         2D Acoustic VTI Medium Forward   
//##  Ps : P + sv wave and get rid of sv        
//##                                     Rong Tao 
//########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "mpi.h"
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"

#define pi 3.141592653

main(int argc,char *argv[])
{
void model_vti(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
             float vdx,float vdz,float favg,float tmax,float dt,float dtout,
             float pfac,char FN1[],char FN2[],char FN3[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
             int is,float **p_cal,float Circle_Radius,
             int mm,int wtype,int hsx,int myid,float *mu_v);       
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs_sxd,int ds_sxd,int zs_sxd,int is,
                     float mu_v,float **p_cal,int tt);        
	int i,j,k,is,nx,nz,nt,vnx,vnz,i_start,i_end,mm,wtype,hsx;
	int ns_sxd,ds_sxd,fs_sxd,zs_sxd,npd;
	float dx,dz,vdx,vdz,tmax,dt,dtout,pfac,favg;
	int myid,numprocs,Circle_Radius,flag_mu;
      float mu_v;
      float **p_cal;


/******** ranks,wavelet,receivers,mute direct **********/
        mm=4;wtype=1;hsx=1;flag_mu=1;
/******************* dat document **********************/
        char FN1[250]={"vel601301.dat"};
        char FN2[250]={"epsilu601301.dat"};
        char FN3[250]={"deta601301.dat"};
	  char FN4[250]={"shot.dat"};

/***************** parameters **************************/

        nx=601;         npd=50;      tmax=1.5;
	  nz=301;         favg=40;     pfac=1000.0;

 	  dx=5.0;   
        dz=5.0;   
     
	  nt=1501;    
        dt=1.0;
     
        ns_sxd=1;       
        fs_sxd=300;      
        ds_sxd=100;
        zs_sxd=1;     


        Circle_Radius=15;//for get rid of SV
/*************************v*****************************/
      vdz=dz;vdx=dx;vnx=nx;vnz=nz;dtout=dt;   
/************************Loop start**************************/

      FILE *fp3;
      fp3=fopen(FN4,"wb");

	p_cal=alloc2float(nt,nx);


/*******************MPI************************/
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
/*******************MPI***********************/	
      MPI_Barrier(MPI_COMM_WORLD);
   if(myid==0)
    {
        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("---   \n");
    }                                                      
/********************IS Loop start**************/
   for(is=1+myid;is<=ns_sxd;is=is+numprocs)	
    {     
      if(myid==0)
        {
         printf("---   IS========%d  \n",is);
         printf("---   The forward is start  !  \n");
        }
      zero2float(p_cal,nt,nx);
	
      model_vti(nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,FN1,FN2,FN3,
            ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_cal,Circle_Radius,
            mm,wtype,hsx,myid,&mu_v);    
/**************** mute direct wave **************/ 
      mute_directwave(flag_mu,nx,nt,dt,favg,dx,dz,fs_sxd,ds_sxd,zs_sxd,is,mu_v,p_cal,85);


      fseek(fp3,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	  for(j=0;j<nt;j++)
	    fwrite(&p_cal[i][j],4L,1,fp3);
	
      MPI_Barrier(MPI_COMM_WORLD);
    } 
/*****************IS Loop end******************/ 		     
   if(myid==0)printf("---   The forward is over    \n"); 
   if(myid==0)printf("---   Complete!!!!!!!!! \n");  

   fclose(fp3);
   free2float(p_cal);
/******************MPI************************/   
   MPI_Finalize();
}
/***********************************func********************************************/
void model_vti(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
           float vdx,float vdz,float favg,float tmax,float dt,float dtout,
           float pfac,char FN1[],char FN2[],char FN3[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
           int is,float **p_cal,float Circle_Radius,
           int mm,int wtype,int hsx,int myid,float *mu_v)
/*******************************************************
Function for VTI medium modeling,2016.9.24

 Ps:  the function of modeling following:
      
          du/dt=1/rho*dp/dx , 
          dw/dt=1/rho*dq/dz ,  
          dp/dt=rho*vpx^2*du/dx+rho*vp0*vpn*dw/dz ,
          dq/dt=rho*vp0*vpn*du/dx+rho*vp0^2*dw/dz ,
                     vpx^2=vp0^2*(1+2*epsilu);
                     vpn^2=vp0^2*(1+2*deta);
********************************************************/
{
void cal_c(int mm,float c[]);
void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
                  float favg,float **s,int wtype,int npd,int is,int ds_sxd);
void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                  float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
                  float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
                  float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,
                  float *coffz2);
void update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
                  float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
                  float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
                  float **s,float **vp,float c[],int npd,
                  float **tpx1,float **tpx0,float **tpz1,float **tpz0,
                  float **tqx1,float **tqx0,float **tqz1,float **tqz0,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                  float **deta,float **epsilu,int fs_sxd,int ds_sxd,int zs_sxd,int is,
                  float Circle_Radius);
void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,
                  int npd,float tmax,float favg,float dtout,float dt,float **vp0,float ndtt);
void pad_vv(int nx,int nz,int npd,float **ee);
void read_vrho(char FN1[],char FN2[],char FN3[],int nx,int nz,float **vv,float **epsilu,float **deta,
               float **rho0,int npd);
void current_shot(float **vp0,float **rho0,float **vp,float **rho,
                  int nx,int nz,int npd,int vnx,int vnz,int ds_sxd,int is);
void initial_coffe(float dt,float d0,int nx,int nz,
                  float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);	

	  int i,j;
	  int ntout,it;
	  float t,ndtt,d0;

	  
	  FILE *fp1;

	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;
     
	  float **vp0;
	  float **rho0;
	  float **vp;
	  float **rho;

	  float **up0;    float **up1;
	  float **wp0;    float **wp1;
        float **uq0;    float **uq1;
	  float **wq0;    float **wq1;

	  float **tpx0;   float **tqx0;
	  float **tpx1;   float **tqx1;
	  float **tpz0;   float **tqz0;
	  float **tpz1;   float **tqz1;

	  float **txx0;   float **txx1;
        float **tzz0;   float **tzz1;

	  float **s,**epsilu,**deta;
     
	  float c[mm];
	  

          cal_c(mm,c);   

      
	  vp0=alloc2float(nz+2*npd,nx+2*npd);                                               
	  rho0=alloc2float(nz+2*npd,nx+2*npd);                                            

        zero2float(vp0,nz+2*npd,nx+2*npd); 
        zero2float(rho0,nz+2*npd,nx+2*npd); 

      epsilu=alloc2float(vnz+2*npd,vnx+2*npd);
      deta=alloc2float(vnz+2*npd,vnx+2*npd);
      zero2float(epsilu,vnz+2*npd,vnx+2*npd);
      zero2float(deta,vnz+2*npd,vnx+2*npd);
/********************************************************/
       read_vrho(FN1,FN2,FN3,vnx,vnz,vp0,epsilu,deta,rho0,npd); 
              
/********************************************************/
      
             pad_vv(nx,nz,npd,epsilu);
             pad_vv(nx,nz,npd,deta);
/********************************************************/
       *mu_v=vp0[1+npd][1+npd]*sqrtf((1+2*epsilu[1+npd][1+npd]));/////////////  
/********************************************************/
      
	 vp=alloc2float(nz+2*npd,nx+2*npd);                                              
	 rho=alloc2float(nz+2*npd,nx+2*npd);   

	 up0=alloc2float(nz+2*npd,nx+2*npd);
	 up1=alloc2float(nz+2*npd,nx+2*npd);
	 wp0=alloc2float(nz+2*npd,nx+2*npd);
	 wp1=alloc2float(nz+2*npd,nx+2*npd);  
	 uq0=alloc2float(nz+2*npd,nx+2*npd);
	 uq1=alloc2float(nz+2*npd,nx+2*npd);
	 wq0=alloc2float(nz+2*npd,nx+2*npd);
	 wq1=alloc2float(nz+2*npd,nx+2*npd); 

	 txx0=alloc2float(nz+2*npd,nx+2*npd);
	 txx1=alloc2float(nz+2*npd,nx+2*npd);
        tzz0=alloc2float(nz+2*npd,nx+2*npd);
	 tzz1=alloc2float(nz+2*npd,nx+2*npd);

	 tpx0=alloc2float(nz+2*npd,nx+2*npd);
	 tpx1=alloc2float(nz+2*npd,nx+2*npd);
	 tpz0=alloc2float(nz+2*npd,nx+2*npd);
	 tpz1=alloc2float(nz+2*npd,nx+2*npd);
	 tqx0=alloc2float(nz+2*npd,nx+2*npd);
	 tqx1=alloc2float(nz+2*npd,nx+2*npd);
	 tqz0=alloc2float(nz+2*npd,nx+2*npd);
	 tqz1=alloc2float(nz+2*npd,nx+2*npd);
      
	 s=alloc2float(nz+2*npd,nx+2*npd);   


       d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);
       dt=dt/1000;
/***********************************************************/

        float *coffx1;float *coffx2;float *coffz1;float *coffz2;
        float *acoffx1;float *acoffx2;float *acoffz1;float *acoffz2;
        coffx1=alloc1float(nx+2*npd);
        coffx2=alloc1float(nx+2*npd);
	  coffz1=alloc1float(nz+2*npd); 
        coffz2=alloc1float(nz+2*npd);
	  acoffx1=alloc1float(nx+2*npd);
	  acoffx2=alloc1float(nx+2*npd);
	  acoffz1=alloc1float(nz+2*npd);
	  acoffz2=alloc1float(nz+2*npd);
        zero1float(coffx1,nx+2*npd);
        zero1float(coffx2,nx+2*npd);
        zero1float(acoffx1,nx+2*npd);
        zero1float(acoffx2,nx+2*npd);
        zero1float(coffz1,nz+2*npd);
        zero1float(coffz2,nz+2*npd);
        zero1float(acoffz1,nz+2*npd);
        zero1float(acoffz2,nz+2*npd);

        initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,
                         acoffx1,acoffx2,acoffz1,acoffz2,npd);

/***********************************************************/
	  ndtt=(int)ndtt;
/*******************zero************************/  
           zero2float(p_cal,nt,nx); 
           zero2float(up0,nz+2*npd,nx+2*npd); 
           zero2float(up1,nz+2*npd,nx+2*npd);   
           zero2float(wp0,nz+2*npd,nx+2*npd);  
           zero2float(wp1,nz+2*npd,nx+2*npd); 
           zero2float(uq0,nz+2*npd,nx+2*npd); 
           zero2float(uq1,nz+2*npd,nx+2*npd);   
           zero2float(wq0,nz+2*npd,nx+2*npd);  
           zero2float(wq1,nz+2*npd,nx+2*npd);  

           zero2float(txx0,nz+2*npd,nx+2*npd); 
           zero2float(txx1,nz+2*npd,nx+2*npd);   
           zero2float(tzz0,nz+2*npd,nx+2*npd);  
           zero2float(tzz1,nz+2*npd,nx+2*npd); 

           zero2float(tpx0,nz+2*npd,nx+2*npd); 
           zero2float(tpx1,nz+2*npd,nx+2*npd);   
           zero2float(tpz0,nz+2*npd,nx+2*npd);  
           zero2float(tpz1,nz+2*npd,nx+2*npd); 

           zero2float(tqx0,nz+2*npd,nx+2*npd); 
           zero2float(tqx1,nz+2*npd,nx+2*npd);   
           zero2float(tqz0,nz+2*npd,nx+2*npd);  
           zero2float(tqz1,nz+2*npd,nx+2*npd); 
			 
           current_shot(vp0,rho0,vp,rho,nx,nz,npd,vnx,vnz,ds_sxd,is);//
        
           pad_vv(nx,nz,npd,vp); 

           pad_vv(nx,nz,npd,rho);
       
          
	  for(i=0;i<=nx+2*npd-1;i++)
	   {
		for(j=0;j<=nz+2*npd-1;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
		   rho[i][j]=1.0/rho[i][j];
		}
	   }
        
	      FILE *fpsnap,*fpsnap1;
              fpsnap=fopen("snap-sv.dat","wb");
              fpsnap1=fopen("snap1-sv.dat","wb");	
	  

    for(it=0,t=0.0;it<nt;it++,t+=dt)
     { 
       
      if(it%100==0&&myid==0)printf("---   is===%d   it===%d\n",is,it);

	ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,npd,is,ds_sxd);
      update_vel(nx,nz,npd,mm,dt,dx,dz,up0,wp0,uq0,wq0,txx0,tzz0,
                up1,wp1,uq1,wq1,txx1,tzz1,rho,c,coffx1,coffx2,coffz1,coffz2);
      update_stress(nx,nz,dt,dx,dz,mm,up0,wp0,uq0,wq0,txx0,tzz0,
                up1,wp1,uq1,wq1,txx1,tzz1,s,vp,c,npd,
                tpx1,tpx0,tpz1,tpz0,tqx1,tqx0,tqz1,tqz0,
                acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,
                fs_sxd,ds_sxd,zs_sxd,is,Circle_Radius);
       
	  for(i=npd;i<npd+nx;i++)  
	  {   
		p_cal[i-npd][it]=tzz1[i][npd+hsx-1]+txx1[i][npd+hsx-1];
                //p_cal[i-npd][it]=txx1[i][npd+hsx-1];
	  }


	  for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			up0[i][j]=up1[i][j];
			wp0[i][j]=wp1[i][j];
                  uq0[i][j]=uq1[i][j];
			wq0[i][j]=wq1[i][j];

			tpx0[i][j]=tpx1[i][j];
			tpz0[i][j]=tpz1[i][j];	
                  tqx0[i][j]=tqx1[i][j];
			tqz0[i][j]=tqz1[i][j];
		
			txx0[i][j]=txx1[i][j];
                  tzz0[i][j]=tzz1[i][j];
		}
	   }

           if((is==1)&&(it%5==0)&&(myid==0))
           {
              fseek(fpsnap,(int)(it/5)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&txx1[i][j],4L,1,fpsnap);
           
              
              fseek(fpsnap1,(int)(it/5)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&tzz1[i][j],4L,1,fpsnap1);
           }
     }//it loop end
/**********************close************************/ 
          fclose(fpsnap);fclose(fpsnap1);
/**********************free*************************/        
          free1float(coffx1);free1float(coffx2);
          free1float(coffz1);free1float(coffz2);
          free1float(acoffx1);free1float(acoffx2);
          free1float(acoffz1);free1float(acoffz2);

          free2float(up0);   free2float(up1);   free2float(uq0);   free2float(uq1);
          free2float(wp0);   free2float(wp1);   free2float(wq0);   free2float(wq1);

          free2float(txx0);  free2float(txx1);  free2float(tzz0);  free2float(tzz1);

          free2float(tpx0);  free2float(tpx1);  free2float(tpz0);  free2float(tpz1);
          free2float(tqx0);  free2float(tqx1);  free2float(tqz0);  free2float(tqz1);

          free2float(s);
	    
}                                                  
/************************************func***************************************/
void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
              float favg,float **s,int wtype,int npd,int is,int ds_sxd)
{
float get_wavelet(float ts,float favg,int wtype);

	    int i,j,ixs,izs,x,z;
	    float tdelay,ts,source,fs;
      
       zero2float(s,nz+2*npd,nx+2*npd);     
	 tdelay=1.0/favg;
       ts=t-tdelay;
       fs=xsn+(is-1)*ds_sxd;
       if(t<=2*tdelay)
         {
          source=get_wavelet(ts,favg,wtype);            
	    ixs = (int)(fs+0.5)+npd-1;
          izs = (int)(zsn+0.5)+npd-1;
         for(j=izs-3;j<=izs+3;j++)
	    { 
		 for(i=ixs-3;i<=ixs+3;i++)
		  {  
		    x=i-ixs;z=j-izs;
                s[i][j]=pfac*source*exp(-z*z-x*x);
		  }
	    }
         
	}
}
/**********************************func**************************************/
float get_wavelet(float ts,float favg,int wtype)
 {

	float x,xx,source;

      source=0.0;
	if(wtype==1)//ricker wavelet
	  {
		  x=favg*ts;
		  xx=x*x;
	        source=(1-2*pi*pi*(xx))*exp(-(pi*pi*xx));
	  }
	else if(wtype==2)//derivative of gaussian
	  {
		  x=(-4)*favg*favg*pi*pi/log(0.1);
		  source=(-2)*pi*pi*ts*exp(-x*ts*ts);
          }
      else if(wtype==3)//derivative of gaussian
          {
              x=(-1)*favg*favg*pi*pi/log(0.1);
              source=exp(-x*ts*ts);
          }
      return (source);
}
/**************************************func******************************************/
void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
           float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
           float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
           float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
		 int ii,i,j,im;
		 float dtxx,dtxz,dtx,dtz,xx,zz;

		 dtx=dt/dx;
		 dtz=dt/dz;
		 for(j=mm;j<=(2*npd+nz-mm-1);j++)
		 { 
			 for(i=mm;i<=(2*npd+nx-mm-1);i++)
			 {
                     xx=0.0;
                     zz=0.0;
			   for(im=0;im<mm;im++)
                            {
                        xx+=c[im]*(txx0[i+im+1][j]-txx0[i-im][j]);
                        zz+=c[im]*(tzz0[i][j+im+1]-tzz0[i][j-im]);
                            }
                     up1[i][j]=coffx2[i]*up0[i][j]-coffx1[i]*dtx*rho[i][j]*xx;
                     wq1[i][j]=coffz2[j]*wq0[i][j]-coffz1[j]*dtz*rho[i][j]*zz;
			 }
		 }
}
/*************************************func**********************************************/
void update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
            float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
            float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
            float **s,float **vp,float c[],int npd,
            float **tpx1,float **tpx0,float **tpz1,float **tpz0,
            float **tqx1,float **tqx0,float **tqz1,float **tqz0,
            float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
            float **deta,float **epsilu,int xsn,int ds_sxd,int zsn,int is,float Circle_Radius)
{
		 int i,j,ii,im,ix,iz;
		 float dtx,dtz,dux,dwz,xx,zz;
             int fs,ixs,izs,CR;
             float **deta1,**epsilu1;

            fs=xsn+(is-1)*ds_sxd;
            ixs=(int)(fs+0.5)+npd-1;
            izs=(int)(zsn+0.5)+npd-1;

            CR=Circle_Radius;///////////////////////

            epsilu1=alloc2float(nz+2*npd,nx+2*npd);
            deta1=alloc2float(nz+2*npd,nx+2*npd);
            zero2float(epsilu1,nz+2*npd,nx+2*npd);
            zero2float(deta1,nz+2*npd,nx+2*npd);

		 dtx=dt/dx;
		 dtz=dt/dz; 
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
		   for(j=mm;j<=(2*npd+nz-mm-1);j++)
		    {
              /**** get the smooth circle to get rid of SV wave ****/
                  ix=i-ixs;
                  iz=j-izs;
                  if((ix*ix+iz*iz)<=CR*CR)
                        {
                    if((ix*ix+iz*iz)<=(CR*CR/16))
                           { 
                       epsilu1[i][j]=0.0;
                       deta1[i][j]=0.0;
                    }else{
                       epsilu1[i][j]=0.5*(1-cos(pi*((pow((ix*ix+iz*iz),0.5)-CR/4.0)*4.0/(CR*3.0-1))))*epsilu[i][j];
                       deta1[i][j]=0.5*(1-cos(pi*((pow((ix*ix+iz*iz),0.5)-CR/4.0)*4.0/(CR*3.0-1))))*deta[i][j];   
                           }
                  }else{
                       epsilu1[i][j]=epsilu[i][j];
                       deta1[i][j]=deta[i][j]; 
                        }          
                  xx=0.0;
                  zz=0.0;
                  for(im=0;im<mm;im++)
                        {
                    xx+=c[im]*(up1[i+im][j]-up1[i-im-1][j]);
                    zz+=c[im]*(wq1[i][j+im]-wq1[i][j-im-1]);
                        }
                 tpx1[i][j]=acoffx2[i]*tpx0[i][j]-acoffx1[i]*vp[i][j]*(1+2*epsilu1[i][j])*dtx*xx;  
                 tpz1[i][j]=acoffz2[j]*tpz0[i][j]-acoffz1[j]*vp[i][j]*(pow((1+2*deta1[i][j]),0.5))*dtz*zz;
                 tqx1[i][j]=acoffx2[i]*tqx0[i][j]-acoffx1[i]*vp[i][j]*(pow((1+2*deta1[i][j]),0.5))*dtx*xx;
                 tqz1[i][j]=acoffz2[j]*tqz0[i][j]-acoffz1[j]*vp[i][j]*dtz*zz;
                 txx1[i][j]=tpx1[i][j]+tpz1[i][j]+s[i][j];
                 tzz1[i][j]=tqx1[i][j]+tqz1[i][j]+s[i][j];
		    }
		 }
}                      
/***************************************func********************************************/
float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,
                   float tmax,float favg,float dtout,float dt,float **vp0,float ndtt)
{
		 int i,j;
		 float vpmax,vpmin,H_min;
		 float dt_max,dx_max,dz_max,d0;

		 vpmax=vp0[npd][npd];
		 vpmin=vp0[npd][npd];
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 if(vpmax<vp0[i][j]) vpmax=vp0[i][j];
				 if(vpmin>vp0[i][j]) vpmin=vp0[i][j];
			 }
		 }
		 d0=3.0*vpmax*log(100000.0)/(2.0*npd*dx);
		 if(dx<dz) H_min=dx;
		 else H_min=dz;
/********determine time sampling interval to ensure stability*****/
		 dt_max=0.5*1000*H_min/vpmax;
             dx_max=vpmin/favg*0.2;
             dz_max=dx_max;

            if(dx_max<dx)
                { 
               printf("dx_max===%f, vpmin===%f, favg===%f \n",dx_max,vpmin,favg);
		   printf("YOU NEED HAVE TO REDEFINE DX ! \n");
               exit(0);
		 }
             if(dz_max<dz)
		 {
		   printf("YOU NEED HAVE TO REDEFINE DZ ! \n");
               exit(0);
		 }
	       if(dt_max<dt)
		 {
               printf("dt_max===%f, H_min===%f, vpmax===%f \n",dt_max,H_min,vpmax);
		   printf("YOU NEED HAVE TO REDEFINE dt ! \n");
               exit(0);
		 }
         
             return d0;
}
/**************************func*************************************/
void pad_vv(int nx,int nz,int npd,float **ee)
{
		 int i,j;


/*****pad left side                    */
            for(j=npd;j<=(nz+npd-1);j++)
	    {
              for(i=0;i<=npd-1;i++)
              { 
               ee[i][j]=ee[npd][j];
              }
	    }
       
/*****pad right side                    */
            for(j=npd;j<=(nz+npd-1);j++)
		{
              for(i=nx+npd;i<=(nx+2*npd-1);i++)
              {
                ee[i][j]=ee[nx+npd-1][j];
              }
		}
/*****pad upper side                    */
            for(j=0;j<=(npd-1);j++)
		{
              for(i=0;i<=(nx+2*npd-1);i++)
              {
                ee[i][j]=ee[i][npd];
              }
		}
/*****lower side                        */
            for(j=nz+npd;j<=(nz+2*npd-1);j++)
		{
              for(i=0;i<=(nx+2*npd-1);i++)
              {
                ee[i][j]=ee[i][nz+npd-1];
              }
		}	
}
/**************************func*************************************/
void read_vrho(char FN1[],char FN2[],char FN3[],int nx,int nz,float **vv,float **epsilu,float **deta,
               float **rho0,int npd)
{

		 int i,j;
		
		 FILE *fp1,*fp2,*fp3;
		 fp1=fopen(FN1,"rb");
		 fp2=fopen(FN2,"rb");
		 fp3=fopen(FN3,"rb");
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 fread(&vv[i][j],4,1,fp1);
				 fread(&epsilu[i][j],4,1,fp2);
				 fread(&deta[i][j],4,1,fp3);

			 }
		 }
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 rho0[i][j]=1.0;
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
}
/**************************func*************************************/
void current_shot(float **vp0,float **rho0,float **vp,float **rho,
                  int nx,int nz,int npd,int vnx,int vnz,int ds_sxd,int is)
{
      
            
                   int ivstart,ivend;
			 int i,ix,iz;
                   is=1;
                  // ivstart=1+(is-1)*ds_sxd;
			// ivend=nx+(is-1)*ds_sxd;
                   ivstart=1;
			 ivend=nx;

			 if(ivstart<=0)
			 {
				 printf("ivstart less than zero \n");
				 exit(0);
			 }
			 if(ivend>vnx)
			 {
				 printf("ivend great than Vnx \n");
				 exit(0);
			 }

			 for(ix=npd+ivstart-1;ix<ivend+npd;ix++)
			 {
				 for(iz=npd;iz<nz+npd;iz++)
				 {
				  vp[ix-ivstart+1][iz]=vp0[ix][iz];
                          rho[ix-ivstart+1][iz]=rho0[ix][iz];
				 }
			 }
}
/**************************func*************************************/
void initial_coffe(float dt,float d0,int nx,int nz,
                   float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                   float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd)
{		
		 int i,j;
         

		 for(i=0;i<npd;i++)
		 {   
			 coffx1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz2[i]=coffz1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);

	
		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 coffx1[i]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
		 }
		 for(i=npd+nz;i<nz+2*npd;i++)
		 {
			 coffz1[i]=1/(1+(dt*d0*pow((0.5+i-nz-npd)/npd,2))/2);
			 coffz2[i]=coffz1[i]*(1-(dt*d0*pow((0.5+i-nz-npd)/npd,2))/2);
		 }

		 for(i=npd;i<npd+nx;i++)
		 {
			 coffx1[i]=1.0;
			 coffx2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 coffz1[i]=1.0;
			 coffz2[i]=1.0;
		 }

		 
		 for(i=0;i<npd;i++)    
		 {    
			 acoffx1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffx2[i]=coffx1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz2[i]=coffz1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);

		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 acoffx1[i]=1/(1+(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
			 acoffx2[i]=coffx1[i]*(1-(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
		 }
		 for(i=npd+nz;i<nz+2*npd;i++)
		 {
			 acoffz1[i]=1/(1+(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
			 acoffz2[i]=coffz1[i]*(1-(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
		 }

		 for(i=npd;i<npd+nx;i++)
		 {
			 acoffx1[i]=1.0;
			 acoffx2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 acoffz1[i]=1.0;
			 acoffz2[i]=1.0;
		 }	       
}
/**************************func****************************/                                                  
void cal_c(int mm,float c[])                                             
{                                                      
	if(mm==2)
	{
        c[0]=1.125;
        c[1]=-0.04166667;
	}
	if(mm==3)
	{
	  c[0]=1.1718750;
        c[1]=-0.065104167;
        c[2]=0.0046875;
	}
	if(mm==4)
	{
	  c[0]=1.196289;
        c[1]=-0.0797526;
        c[2]=0.009570313;
        c[3]=-0.0006975447;
	}
	if(mm==5)
	{
	  c[0]=1.211243;
        c[1]=-0.08972168;
        c[2]=0.01384277;
        c[3]=-0.00176566;
        c[4]=0.0001186795;
	}
	if(mm==6)
	{
        c[0]=1.2213364;
        c[1]=-0.096931458;
        c[2]=0.017447662;
        c[3]=-0.0029672895;
        c[4]=0.0003590054;
        c[5]=-0.000021847812;
   	}
 	if(mm==7)
  	{
        c[0]=1.2286062;
        c[1]=-0.10238385;
        c[2]=0.020476770;
        c[3]=-0.0041789327;
        c[4]=0.00068945355;
        c[5]=-0.000076922503;
        c[6]=0.0000042365148;
        }
      if(mm==8)
        {
        c[0]=1.2340911;
        c[1]=-0.10664985;
        c[2]=0.023036367;
        c[3]=-0.0053423856;
        c[4]=0.0010772712;
        c[5]=-0.00016641888;
        c[6]=0.000017021711;
        c[7]=-0.00000085234642;
   	}                                                                         
}  
/**************************func****************************/    
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs_sxd,int ds_sxd,int zs_sxd,int is,
                     float mu_v,float **p_cal,int tt)
{
  int i,j,mu_t,mu_nt;
  float mu_x,mu_z,mu_t0;

    if(flag_mu)   
     for(i=0;i<nx;i++)
       {
        mu_x=dx*abs(i-fs_sxd-(is-1)*ds_sxd);
        mu_z=dz*zs_sxd;
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_z,2))/mu_v;
        mu_t=(int)(2.0/(dt/1000*favg));
        mu_nt=(int)(mu_t0/dt*1000)+mu_t+tt;
        for(j=0;j<nt;j++)if(j<mu_nt)
           p_cal[i][j]=0.0;
       }else{}
}


