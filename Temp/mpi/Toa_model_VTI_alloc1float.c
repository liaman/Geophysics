//#########################################################
//##         2D Acoustic VTI Medium Forward   
//##  Ps : P + sv wave and get rid of sv        
//##                                     Rong Tao 
/****************************
Function for VTI medium modeling,2017.2.5

 Ps:  the function of modeling following:
      
          du/dt=1/rho*dp/dx , 
          dw/dt=1/rho*dq/dz ,  
          dp/dt=rho*vpx^2*du/dx+rho*vp*vpn*dw/dz ,
          dq/dt=rho*vp*vpn*du/dx+rho*vp^2*dw/dz ,
                     vpx^2=vp^2*(1+2*epsilu);
                     vpn^2=vp^2*(1+2*deta);
****************************/
//########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "/home/rongtao/gpfs03/hc/cjbsegy.h"
#include "/home/rongtao/gpfs03/hc/fft.c"
#include "/home/rongtao/gpfs03/hc/alloc.c"
#include "/home/rongtao/gpfs03/hc/complex.c"

#define pi 3.141592653

/******************func********************/
void ptsource(float pfac,float xsn,float zsn,int nx,int nz,int nnx,int nnz,float dt,float t,
              float favg,float *s,int wtype,int npd,int is,int ds);
float get_wavelet(float ts,float favg,int wtype);
void update_vel(int nx,int nz,int nnx,int nnz,int npd,int mm,float dt,float dx,float dz,
           float *u0,float *w0,float *u1,float *w1,float *P,float *Q,
           float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2);
void update_stress(int nx,int nz,int nnx,int nnz,float dt,float dx,float dz,int mm,
            float *u1,float *w1,float *P,float *Q,float *s,float *vp,float c[],int npd,
            float *px1,float *px0,float *pz1,float *pz0,float *qx1,float *qx0,float *qz1,float *qz0,
            float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
            float *deta,float *epsilu,int xsn,int ds,int zsn,int is,float Circle_iso,int flag_SV);        
float get_constant(float dx,float dz,int nx,int nz,int nnx,int nnz,int nt,int npd,float favg,float dt,float *vp);
void pad_vv(int nx,int nz,int nnx,int nnz,int npd,float *ee);
void read_file(char FN1[],char FN2[],char FN3[],int nx,int nz,int nnx,int nnz,float *vv,float *epsilu,float *deta,int npd);
void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                   float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);                                              
void cal_c(int mm,float c[]);         
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,float dx,float dz,int fs,int ds,int zs,int is,
                     float mu_v,float *p_cal,int tt);
//a########################################################################
int main(int argc,char *argv[])
{
	int i, j, k, is, it, nx, nz, nnx, nnz, nt, mm, wtype, hsx;
	int ns, ds, fs, zs, npd;
	float dx, dz, dt, t, d0, pfac, favg;
	int Circle_iso, flag_mu, flag_SV;
        float mu_v;
        float *p_cal;

/**** ranks,wavelet,receivers,mute direct *****/
          mm=4;wtype=1;hsx=1;npd=30;flag_mu=1;flag_SV=0;
/********** dat document ***********/
          char FN1[250]={"vel5050.dat"};
          char FN2[250]={"epsilu5050.dat"};
          char FN3[250]={"deta5050.dat"};
	  char FN4[250]={"shot.dat"};

/********* parameters *************/

          nx=50;              
	  nz=50;         favg=20;     pfac=10.0;

 	  dx=5.0;   
          dz=5.0;   
     
	  nt=701;    
          dt=0.0005;
     
          ns=1;       
          fs=25;      
          ds=10;
          zs=25;     


          Circle_iso=15;
/*************v***************/ 
          nnx=nx+2*npd;
          nnz=nz+2*npd;
/************Loop start*************/

          FILE *fp3;
          fp3=fopen(FN4,"wb");
	  FILE *fp1;
	  p_cal=alloc1float(nt*nx);


	  float *vp, *epsilu, *deta;
	  float *u0, *u1, *px0, *qx0, *px1, *qx1;
          float *w0, *w1, *pz0, *qz0, *pz1, *qz1;
	  float *P, *Q, *s;
	  float c[mm];

          cal_c(mm,c);

   	vp=alloc1float(nnx*nnz); 
        epsilu=alloc1float(nnx*nnz);
        deta=alloc1float(nnx*nnz);
        read_file(FN1,FN2,FN3,nx,nz,nnx,nnz,vp,epsilu,deta,npd); 
              
/****************************/
        pad_vv(nx,nz,nnx,nnz,npd,epsilu);
        pad_vv(nx,nz,nnx,nnz,npd,deta);
        pad_vv(nx,nz,nnx,nnz,npd,vp); 
/****************************/
        mu_v=vp[npd]*sqrtf((1+2*epsilu[npd]));printf("surface vel >> %.2f\n",mu_v);

/****************************/
	 u0=alloc1float(nnx*nnz);	 u1=alloc1float(nnx*nnz);
	 w0=alloc1float(nnx*nnz);	 w1=alloc1float(nnx*nnz); 
	 P=alloc1float(nnx*nnz);         Q=alloc1float(nnx*nnz);
	 px0=alloc1float(nnx*nnz);	 px1=alloc1float(nnx*nnz);
	 pz0=alloc1float(nnx*nnz);	 pz1=alloc1float(nnx*nnz);
	 qx0=alloc1float(nnx*nnz);	 qx1=alloc1float(nnx*nnz);
	 qz0=alloc1float(nnx*nnz);	 qz1=alloc1float(nnx*nnz);
	 s=alloc1float(nnx*nnz);   

        d0=get_constant(dx,dz,nx,nz,nnx,nnz,nt,npd,favg,dt,vp);
/******************************/

        float *coffx1,*coffx2,*coffz1,*coffz2,*acoffx1,*acoffx2,*acoffz1,*acoffz2;
        coffx1=alloc1float(nnx);        coffx2=alloc1float(nnx);
	coffz1=alloc1float(nnz);        coffz2=alloc1float(nnz); 
	acoffx1=alloc1float(nnx);	acoffx2=alloc1float(nnx);
	acoffz1=alloc1float(nnz);	acoffz2=alloc1float(nnz);

        initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/**********zero************/  
	FILE *fpsnap,*fpsnap1;
        fpsnap=fopen("snap-sv.dat","wb");
        fpsnap1=fopen("snap1-sv.dat","wb");
        printf("--------------------------------------------------------\n");
        printf("---   \n");                                                     
/**********IS Loop start*******/
   for(is=1;is<=ns;is++)	
    {     
         printf("---   IS========%d  \n",is);
           zero1float(p_cal,nt*nx);

	   zero1float(u0,nnx*nnz);        zero1float(u1,nnx*nnz); 
           zero1float(w0,nnx*nnz);        zero1float(w1,nnx*nnz); 
           zero1float(P,nnx*nnz);         zero1float(Q,nnx*nnz); 
           zero1float(px0,nnx*nnz);       zero1float(px1,nnx*nnz); 
           zero1float(pz0,nnx*nnz);       zero1float(pz1,nnx*nnz); 
           zero1float(qx0,nnx*nnz);       zero1float(qx1,nnx*nnz); 
           zero1float(qz0,nnx*nnz);       zero1float(qz1,nnx*nnz); 

     for(it=0,t=dt;it<nt;it++,t+=dt)
     { 
       
      //if(it%100==0&&is==1)printf("---   is===%d   it===%d\n",is,it);

	ptsource(pfac,fs,zs,nx,nz,nnx,nnz,dt,t,favg,s,wtype,npd,is,ds);
        update_vel(nx,nz,nnx,nnz,npd,mm,dt,dx,dz,u0,w0,u1,w1,P,Q,c,coffx1,coffx2,coffz1,coffz2);
        update_stress(nx,nz,nnx,nnz,dt,dx,dz,mm,u1,w1,P,Q,s,vp,c,npd,px1,px0,pz1,pz0,qx1,qx0,qz1,qz0,
                      acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,fs,ds,zs,is,Circle_iso,flag_SV);
       
	/*  for(i=npd;i<npd+nx;i++)  
	  {   
		p_cal[i-npd it]=tzz1[i npd+hsx-1]+txx1[i npd+hsx-1];
                //p_cal[i-npd it]=txx1[i npd+hsx-1];
	  }  */


		for(i=0;i<nnx*nnz;i++)
		{
			u0[i]=u1[i];     w0[i]=w1[i];
			px0[i]=px1[i];   pz0[i]=pz1[i];
                        qx0[i]=qx1[i];   qz0[i]=qz1[i];
		}

           if((is==1)&&(it%50==0))
           {
              fseek(fpsnap,(int)(it/50)*(nnx)*(nnz)*4L,0);
              fwrite(P,4L,nnx*nnz,fpsnap);
           
              
              fseek(fpsnap1,(int)(it/50)*(nnx)*(nnz)*4L,0);
              fwrite(Q,4L,nnx*nnz,fpsnap1);
           }
     }//it loop end
    } 
/*********IS Loop end*********/ 		     
   printf("---   The forward is over    \n"); 
   printf("---   Complete!!!!!!!!! \n");  

   fclose(fp3);
   free1float(p_cal);
/***********close************/ 
          fclose(fpsnap);fclose(fpsnap1);
/***********free*************/        
          free1float(coffx1);free1float(coffx2);
          free1float(coffz1);free1float(coffz2);
          free1float(acoffx1);free1float(acoffx2);
          free1float(acoffz1);free1float(acoffz2);

          free1float(u0);   free1float(u1);  
          free1float(w0);   free1float(w1);

          free1float(P);  free1float(Q);

          free1float(px0);  free1float(px1);  free1float(pz0);  free1float(pz1);
          free1float(qx0);  free1float(qx1);  free1float(qz0);  free1float(qz1);

          free1float(s);
}
/******************func********************/
void ptsource(float pfac,float xsn,float zsn,int nx,int nz,int nnx,int nnz,float dt,float t,
              float favg,float *s,int wtype,int npd,int is,int ds)
{
float get_wavelet(float ts,float favg,int wtype);

       int i,j,ixs,izs,x,z;
       float tdelay,ts,source,fs;
      
       zero1float(s,nnx*nnz);     
       tdelay=1.0/favg;
       ts=t-tdelay;
       fs=xsn+(is-1)*ds;
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
                    s[i*nnz+j]=pfac*source*exp(-z*z-x*x);
		  }
	    }
       }
}
/*****************func*******************/
float get_wavelet(float ts,float favg,int wtype)
 {
	float x,xx,source;

        source=0.0;
	if(wtype==1)//ricker wavelet
	{
          x=favg*ts;
          xx=x*x;
          source=(1-2*pi*pi*(xx))*exp(-(pi*pi*xx));
	}else if(wtype==2){//derivative of gaussian
          x=(-4)*favg*favg*pi*pi/log(0.1);
          source=(-2)*pi*pi*ts*exp(-x*ts*ts);
        }else if(wtype==3){//derivative of gaussian
          x=(-1)*favg*favg*pi*pi/log(0.1);
          source=exp(-x*ts*ts);
        }
        return (source);
}
/*******************func*********************/
void update_vel(int nx,int nz,int nnx,int nnz,int npd,int mm,float dt,float dx,float dz,
           float *u0,float *w0,float *u1,float *w1,float *P,float *Q,
           float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
		 int ix,iz,im,id;
		 float dtx,dtz,xx,zz;

		 dtx=dt/dx;
		 dtz=dt/dz;
                 for(id=mm;id<nnx*nnz-mm;id++)
                 {
                   ix=id/nnz;
                   iz=id%nnz;
                   if(ix>=mm&&ix<(nnx-mm)&&iz>=mm&&iz<(nnz-mm))
                   {
                     xx=0.0;

                     zz=0.0;
	             for(im=0;im<mm;im++)
                     {
                        xx+=c[im]*(P[id+(im+1)*nnz]-P[id-im*nnz]);
                        zz+=c[im]*(Q[id+im+1]      -Q[id-im]);
                     }
                     u1[id]=coffx2[ix]*u0[id]-coffx1[ix]*dtx*xx;
                     w1[id]=coffz2[iz]*w0[id]-coffz1[iz]*dtz*zz;
                   }
                 }
}
/*******************func***********************/
void update_stress(int nx,int nz,int nnx,int nnz,float dt,float dx,float dz,int mm,
            float *u1,float *w1,float *P,float *Q,float *s,float *vp,float c[],int npd,
            float *px1,float *px0,float *pz1,float *pz0,float *qx1,float *qx0,float *qz1,float *qz0,
            float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
            float *deta,float *epsilu,int xsn,int ds,int zsn,int is,float Circle_iso,int flag_SV)
{
		 int i,j,ii,im,ix,iz,rx,rz,id;
		 float dtx,dtz, xx,zz;
                 int fs,ixs,izs,CR;

            float *deta1,*epsilu1;

            fs=xsn+(is-1)*ds;
            ixs=(int)(fs+0.5)+npd-1;
            izs=(int)(zsn+0.5)+npd-1;

            CR=Circle_iso;///////////////////////

            epsilu1=alloc1float(nnx*nnz);
            deta1=alloc1float(nnx*nnz);

                 dtx=dt/dx;
		 dtz=dt/dz;
                 for(id=mm;id<nnx*nnz-mm;id++)
                 {
                   ix=id/nnz;
                   iz=id%nnz;

              /** get the smooth circle to get rid of SV wave **/
                  rx=ix-ixs;
                  rz=iz-izs;
               if(flag_SV){
                  if((rx*rx+rz*rz)<=CR*CR){
                       if((rx*rx+rz*rz)<=(CR*CR/16)){ 
                              epsilu1[id]=0.0;
                              deta1[id]=0.0;
                       }else{
                              epsilu1[id]=0.5*(1-cos(pi*((pow((rx*rx+rz*rz),0.5)-CR/4.0)*4.0/(CR*3.0-1))))*epsilu[id];
                              deta1[id]  =0.5*(1-cos(pi*((pow((rx*rx+rz*rz),0.5)-CR/4.0)*4.0/(CR*3.0-1))))*deta[id];   
                       }
                  }else{
                       epsilu1[id]=epsilu[id];
                       deta1[id]  =deta[id]; 
                  }  
               }else{
                  epsilu1[id]=epsilu[id];
                  deta1[id]  =deta[id]; 
               }  
              /** get the smooth circle to get rid of SV wave **/


                   if(ix>=mm&&ix<(nnx-mm)&&iz>=mm&&iz<(nnz-mm))
                   {
                     xx=0.0;
                     zz=0.0;
	             for(im=0;im<mm;im++)
                     {
                        xx+=c[im]*(u1[id+im*nnz]-u1[id-(im+1)*nnz]);
                        zz+=c[im]*(w1[id+im]    -w1[id-im-1]);
                     }
                     px1[id]=acoffx2[ix]*px0[id]-acoffx1[ix]*vp[id]*vp[id]*(1+2*epsilu1[id])*dtx*xx;
                     pz1[id]=acoffz2[iz]*pz0[id]-acoffz1[iz]*vp[id]*vp[id]*sqrtf(1+2*deta1[id])*dtz*zz;
                     qx1[id]=acoffx2[ix]*qx0[id]-acoffx1[ix]*vp[id]*vp[id]*sqrtf(1+2*deta1[id])*dtx*xx;
                     qz1[id]=acoffz2[iz]*qz0[id]-acoffz1[iz]*vp[id]*vp[id]*dtz*zz;

                     P[id]=px1[id]+pz1[id]+s[id];
                     Q[id]=qx1[id]+qz1[id]+s[id];
                   }
                 }
}                      
/********************func**********************/
float get_constant(float dx,float dz,int nx,int nz,int nnx,int nnz,int nt,int npd,float favg,float dt,float *vp)
{
		 int i,j,id;
		 float vpmax,vpmin,H_min;
		 float dt_max,dx_max,dz_max,d0;

		 vpmax=vp[npd];
		 vpmin=vp[npd];
		 for(id=npd;id<nnx*nnz;id++)
                 {
			if(vpmax<vp[id]) vpmax=vp[id];
			if(vpmin>vp[id]) vpmin=vp[id];
                 }
		 d0=3.0*vpmax*log(100000.0)/(2.0*npd*dx);
		 if(dx<dz) H_min=dx;
		 else H_min=dz;
/****determine time sampling interval to ensure stability***/
		 dt_max=0.5*H_min/vpmax;
                 dx_max=vpmin/favg*0.2;
                 dz_max=dx_max;

                if(dx_max<dx)
                { 
                   printf("dx_max=%f, vpmin=%f, favg=%f \n",dx_max,vpmin,favg);
		   printf("Redefine <dx> !\n");

                   exit(0);
		}
                if(dz_max<dz)
		{
		   printf("Redefine <dz> !\n");
                   exit(0);
		}
	        if(dt_max<dt)
		{
                   printf("dt_max=%f, H_min=%f, vpmax=%f \n",dt_max,H_min,vpmax);
		   printf("Redefine <dt> !\n");
                   exit(0);
		}
             return d0;
}
/*************func*******************/
void pad_vv(int nx,int nz,int nnx,int nnz,int npd,float *ee)
{
     int ix,iz,id;
 
     for(id=0;id<nnx*nnz;id++)
     {
       ix=id/nnz;
       iz=id%nnz;
       if(ix<npd){
           ee[id]=ee[npd*nnz+iz];  //left
       }else if(ix>=nnx-npd){
           ee[id]=ee[(nnx-npd-1)*nnz+iz];//right
       }
     }
     for(id=0;id<nnx*nnz;id++)
     {
       ix=id/nnz;
       iz=id%nnz;
       if(iz<npd){
           ee[id]=ee[ix*nnz+npd];//up
       }else if(iz>=nnz-npd){
           ee[id]=ee[ix*nnz+nnz-npd-1];//down
       }
       if(ee[id]==0){printf("ee[%d][%d]==0\n",ix,iz);exit(0);}
     }
}
/*************func*******************/
void read_file(char FN1[],char FN2[],char FN3[],int nx,int nz,int nnx,int nnz,float *vv,float *epsilu,float *deta,int npd)
{
		 int i,j,id;
		
		 FILE *fp1,*fp2,*fp3;
		 if((fp1=fopen(FN1,"rb"))==NULL)printf("error open <%s>!\n",FN1);
		 if((fp2=fopen(FN2,"rb"))==NULL)printf("error open <%s>!\n",FN2);
		 if((fp3=fopen(FN3,"rb"))==NULL)printf("error open <%s>!\n",FN3);
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
                            id=i*nnz+j;
				 fread(&vv[id],4L,1,fp1);
				 fread(&epsilu[id],4L,1,fp2);
				 fread(&deta[id],4L,1,fp3);
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
}
/*************func*******************/
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
/*************func**************/                                                  
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
/*************func**************/    
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float mu_v,float *p_cal,int tt)
{
  int i,j,mu_t,mu_nt;
  float mu_x,mu_z,mu_t0;

    if(flag_mu)   
     for(i=0;i<nx;i++)
       {
        mu_x=dx*abs(i-fs-(is-1)*ds);
        mu_z=dz*zs;
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_z,2))/mu_v;
        mu_t=(int)(2.0/(dt/1000*favg));
        mu_nt=(int)(mu_t0/dt*1000)+mu_t+tt;
        for(j=0;j<nt;j++)if(j<mu_nt)
           p_cal[i]=0.0;
       }else{}
}
