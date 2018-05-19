//#########################################################
//##         2D Acoustic VTI Medium Forward   
//##  Ps : P + sv wave and get rid of sv        
//##                                     Rong Tao 
//########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include "mpi.h"
#include "/home/rongtao/gpfs03/hc/cjbsegy.h"
#include "/home/rongtao/gpfs03/hc/fft.c"
#include "/home/rongtao/gpfs03/hc/alloc.c"
#include "/home/rongtao/gpfs03/hc/complex.c"

#define pi 3.141592653
float stencil[5]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};
main(int argc,char *argv[])
{
void model_vti(int nx,int nz,int nt,int npd,float dx,float dz,float favg,float tmax,float dt,
             float pfac,char FN1[],char FN2[],char FN3[],int ns,int ds,int fs,int zs,
             int is,float **p_cal,float **P_up,float **P_do,float **P_le,float **P_ri,float **Q_up,float **Q_do,float **Q_le,float **Q_ri,
             float Circle_Radius,int mm,int wtype,int hsx,int myid,float *mu_v,float **P,float **Q);       
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float mu_v,float **p_cal,int tt);        
	int i,j,k,is,nx,nz,nt,i_start,i_end,mm,wtype,hsx;
	int ns,ds,fs,zs,npd;
	float dx,dz,tmax,dt,pfac,favg;
	int myid,numprocs,Circle_Radius,flag_mu;
        float mu_v;
        float **p_cal,**P_up,**P_do,**P_le,**P_ri,**Q_up,**Q_do,**Q_le,**Q_ri;
        float **P,**Q;


/******** ranks,wavelet,receivers,mute direct **********/
        mm=4;wtype=1;hsx=1;flag_mu=1;
/******************* dat document **********************/
        char FN1[250]={"vel5151.dat"};
        char FN2[250]={"epsilu5151.dat"};
        char FN3[250]={"deta5151.dat"};
	char FN4[250]={"shot.dat"};

/***************** parameters **************************/

        nx=51;         npd=50;      tmax=1.5;
	nz=51;         favg=40;     pfac=1000.0;

 	dx=5.0;   
        dz=5.0;   
     
	nt=301;    
        dt=1.0;
     
        ns=1;       
        fs=nx/ns/2;      
        ds=nx/ns;
        zs=1;     


        Circle_Radius=15;//for get rid of SV

/************************Loop start**************************/

      FILE *fp3;
      fp3=fopen(FN4,"wb");

	P=alloc2float(nz,nx);	        Q=alloc2float(nz,nx);

	p_cal=alloc2float(nt,nx);
	P_up=alloc2float(nt,nx);	Q_up=alloc2float(nt,nx);
	P_do=alloc2float(nt,nx);	Q_do=alloc2float(nt,nx);
	P_le=alloc2float(nt,nz);	Q_le=alloc2float(nt,nz);
	P_ri=alloc2float(nt,nz);	Q_ri=alloc2float(nt,nz);


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
   for(is=1+myid;is<=ns;is=is+numprocs)	
    {     
      if(myid==0)
        {
         printf("---   IS========%d  \n",is);
         printf("---   The forward is start  !  \n");
        }
      zero2float(p_cal,nt,nx);
      zero2float(P_up,nt,nx);      zero2float(Q_up,nt,nx);
      zero2float(P_do,nt,nx);      zero2float(Q_do,nt,nx);
      zero2float(P_le,nt,nz);      zero2float(Q_le,nt,nz);
      zero2float(P_ri,nt,nz);      zero2float(Q_ri,nt,nz);
	
      model_vti(nx,nz,nt,npd,dx,dz,favg,tmax,dt,pfac,FN1,FN2,FN3,
            ns,ds,fs,zs,is,p_cal,P_up,P_do,P_le,P_ri,Q_up,Q_do,Q_le,Q_ri,
            Circle_Radius,mm,wtype,hsx,myid,&mu_v,P,Q);    
/**************** mute direct wave **************/ 
      mute_directwave(flag_mu,nx,nt,dt,favg,dx,dz,fs,ds,zs,is,mu_v,p_cal,25);


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

   free2float(P);      free2float(Q);
   free2float(p_cal);
   free2float(P_up);   free2float(Q_up);
   free2float(P_do);   free2float(Q_do);
   free2float(P_le);   free2float(Q_le);
   free2float(P_ri);   free2float(Q_ri);
/******************MPI************************/   
   MPI_Finalize();exit(0);
}
/***********************************func********************************************/
void model_vti(int nx,int nz,int nt,int npd,float dx,float dz,float favg,float tmax,float dt,
           float pfac,char FN1[],char FN2[],char FN3[],int ns,int ds,int fs,int zs,
           int is,float **p_cal,float **P_up,float **P_do,float **P_le,float **P_ri,float **Q_up,float **Q_do,float **Q_le,float **Q_ri,
           float Circle_Radius,int mm,int wtype,int hsx,int myid,float *mu_v,float **P,float **Q)
/*******************************************************
Function for VTI medium modeling,2017.1.31
    
                               Rong Tao
********************************************************/
{
void add_source(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
                  float favg,float **s,int wtype,int npd,int is,int ds);
void fd_vti_2d_PQ(int flag,int nx,int nz,float dt,float dx,float dz,int mm,float **P0,float **Q0,float **P1,float **Q1,
                  float **vp,int npd,float **deta,float **epsilu,int fs,int ds,int zs,int is,float Circle_Radius, float **s);
void pad_vv(int nx,int nz,int npd,float **ee);
void read_file(char FN1[],char FN2[],char FN3[],int nx,int nz,float **vv,float **epsilu,float **deta,int npd);
void absorb_bndr(float **P0,float **P1,float **Q0,float **Q1,int nx,int nz,int npd,float qp) ;

	  int i,j;
	  int it;
	  float t,temp;

     
	  float **vp;

	  float **P0;   float **P1;
          float **Q0;   float **Q1;

	  float **s,**epsilu,**deta;
	  
      
	  vp=alloc2float(nz+2*npd,nx+2*npd);                                              

        zero2float(vp,nz+2*npd,nx+2*npd); 

      epsilu=alloc2float(nz+2*npd,nx+2*npd);
      deta=alloc2float(nz+2*npd,nx+2*npd);
      zero2float(epsilu,nz+2*npd,nx+2*npd);
      zero2float(deta,nz+2*npd,nx+2*npd);
/********************************************************/
       read_file(FN1,FN2,FN3,nx,nz,vp,epsilu,deta,npd); 
              
/********************************************************/
      
             pad_vv(nx,nz,npd,epsilu);
             pad_vv(nx,nz,npd,deta);
/********************************************************/
       *mu_v=vp[1+npd][1+npd]*sqrtf((1+2*epsilu[1+npd][1+npd]));/////////////  
/********************************************************/
	 P0=alloc2float(nz+2*npd,nx+2*npd);
	 P1=alloc2float(nz+2*npd,nx+2*npd);
         Q0=alloc2float(nz+2*npd,nx+2*npd);
	 Q1=alloc2float(nz+2*npd,nx+2*npd);

      
	 s=alloc2float(nz+2*npd,nx+2*npd);   

       dt=dt/1000;
/*******************zero************************/  
           zero2float(p_cal,nt,nx); 
 

           zero2float(P0,nz+2*npd,nx+2*npd); 
           zero2float(P1,nz+2*npd,nx+2*npd);   
           zero2float(Q0,nz+2*npd,nx+2*npd);  
           zero2float(Q1,nz+2*npd,nx+2*npd); 

           pad_vv(nx,nz,npd,vp); 
       
          
	  for(i=0;i<=nx+2*npd-1;i++)
		for(j=0;j<=nz+2*npd-1;j++)
		   vp[i][j]=vp[i][j]*vp[i][j];
        
	      FILE *fpsnap,*fpsnap1;
              fpsnap=fopen("snap-sv.dat","wb");
              fpsnap1=fopen("snap1-sv.dat","wb");	
	  




    for(it=0,t=0.0;it<nt;it++,t+=dt)
     { 
       
      if(it%50==0&&myid==0)printf("---   is=%d, it=%4d, \n",is,it);

	add_source(pfac,fs,zs,nx,nz,dt,t,favg,s,wtype,npd,is,ds);
        fd_vti_2d_PQ(1,nx,nz,dt,dx,dz,mm,P0,Q0,P1,Q1,vp,npd,deta,epsilu,
                  fs,ds,zs,is,Circle_Radius,s);
        absorb_bndr(P0,P1,Q0,Q1,nx,nz,npd,-0.1);
	  for(i=npd;i<npd+nx;i++)  
	  {   
		p_cal[i-npd][it]=Q0[i][npd+hsx-1]+P0[i][npd+hsx-1];

		P_up[i-npd][it]=P0[i][npd];		Q_up[i-npd][it]=Q0[i][npd];
		P_do[i-npd][it]=P0[i][npd+nz-1];	Q_do[i-npd][it]=Q0[i][npd+nz-1];
	  }for(j=npd;j<npd+nz;j++){   
		P_le[j-npd][it]=P0[npd][j];		Q_le[j-npd][it]=Q0[npd][j];
		P_ri[j-npd][it]=P0[npd+nx-1][j];	Q_ri[j-npd][it]=Q0[npd+nx-1][j];
	  }


	  for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			temp=P1[i][j];P1[i][j]=P0[i][j];P0[i][j]=temp;
			temp=Q1[i][j];Q1[i][j]=Q0[i][j];Q0[i][j]=temp;
		}
	   }

           if((is==1)&&(it%20==0)&&(myid==0))
           {
              fseek(fpsnap,(int)(it/20)*(nx+2*npd)*(nz+2*npd)*4L,0);
              for(i=0;i<nx+2*npd;i++)
                 for(j=0;j<nz+2*npd;j++)
                    fwrite(&P1[i][j],4L,1,fpsnap);
           
              
              fseek(fpsnap1,(int)(it/20)*(nx+2*npd)*(nz+2*npd)*4L,0);
              for(i=0;i<nx+2*npd;i++)
                 for(j=0;j<nz+2*npd;j++)
                    fwrite(&Q1[i][j],4L,1,fpsnap1);
           }
     }//it loop end
	  for(j=0;j<nz;j++)
	    for(i=0;i<nx;i++)
               {P[i][j]=P0[i+npd][j+npd];Q[i][j]=Q0[i+npd][j+npd];}
/**********************close************************/ 
          fclose(fpsnap);fclose(fpsnap1);
/**********************free*************************/        
          free2float(P0);  free2float(P1);  free2float(Q0);  free2float(Q1);

          free2float(s);
	    
}  
/************************************func***************************************/      
void absorb_bndr(float **P0,float **P1,float **Q0,float **Q1,int nx,int nz,int npd,float qp) 
{
   int ix,iz;
   for(ix=0;ix<nx+2*npd;ix++){
      for(iz=0;iz<nz+2*npd;iz++){
         if(ix<npd){
               P0[ix][iz]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*P0[ix][iz];
               P1[ix][iz]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*P1[ix][iz];
               Q0[ix][iz]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*Q0[ix][iz];
               Q1[ix][iz]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*Q1[ix][iz];
         }else if(ix>=nx+npd){
               P0[ix][iz]=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 )*P0[ix][iz];
               P1[ix][iz]=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 )*P1[ix][iz];
               Q0[ix][iz]=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 )*Q0[ix][iz];
               Q1[ix][iz]=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 )*Q1[ix][iz];
         }if(iz<npd){
               P0[ix][iz]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*P0[ix][iz];
               P1[ix][iz]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*P1[ix][iz];
               Q0[ix][iz]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*Q0[ix][iz];
               Q1[ix][iz]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*Q1[ix][iz];
         }else if(iz>=nz+npd){
               P0[ix][iz]=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 )*P0[ix][iz];
               P1[ix][iz]=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 )*P1[ix][iz];
               Q0[ix][iz]=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 )*Q0[ix][iz];
               Q1[ix][iz]=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 )*Q1[ix][iz];
         }
      }
   }
}                       
/************************************func***************************************/
void add_source(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
              float favg,float **s,int wtype,int npd,int is,int ds)
{
float get_wavelet(float ts,float favg,int wtype);

	    int i,j,ixs,izs,x,z;
	    float tdelay,ts,source,fs;
      
       zero2float(s,nz+2*npd,nx+2*npd);     
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

/*************************************func**********************************************/
void fd_vti_2d_PQ(int flag,int nx,int nz,float dt,float dx,float dz,int mm,float **P0,float **Q0,float **P1,float **Q1,
            float **vp,int npd,float **deta,float **epsilu,int xsn,int ds,int zsn,int is,float Circle_Radius,float **s)
{
          
		 int i,j,ii,im,ix,iz;
		 float dt2x2,dt2z2,xx,zz,ss;
       /*      int fs,ixs,izs,CR;
             float **deta1,**epsilu1;

            fs=xsn+(is-1)*ds;
            ixs=(int)(fs+0.5)+npd-1;
            izs=(int)(zsn+0.5)+npd-1;

            CR=Circle_Radius;///////////////////////

            epsilu1=alloc2float(nz+2*npd,nx+2*npd);
            deta1=alloc2float(nz+2*npd,nx+2*npd);
            zero2float(epsilu1,nz+2*npd,nx+2*npd);
            zero2float(deta1,nz+2*npd,nx+2*npd);*/

		 dt2x2=dt*dt/(dx*dx);
		 dt2z2=dt*dt/(dz*dz); 
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
		   for(j=mm;j<=(2*npd+nz-mm-1);j++)
		    {
              /**** get the smooth circle to get rid of SV wave ****/
              /*    ix=i-ixs;
                  iz=j-izs;
                  if((ix*ix+iz*iz)<=CR*CR){
                       if((ix*ix+iz*iz)<=(CR*CR/16)){ 
                          epsilu1[i][j]=0.0;
                          deta1[i][j]=0.0;
                       }else{
                          epsilu1[i][j]=0.5*(1-cos(pi*((pow((ix*ix+iz*iz),0.5)-CR/4.0)*4.0/(CR*3.0-1))))*epsilu[i][j];
                          deta1[i][j]=0.5*(1-cos(pi*((pow((ix*ix+iz*iz),0.5)-CR/4.0)*4.0/(CR*3.0-1))))*deta[i][j];   
                       }
                  }else{
                          epsilu1[i][j]=epsilu[i][j];
                          deta1[i][j]=deta[i][j]; 
                       } */      
 
                      xx=stencil[0]*P1[i][j];
                      zz=stencil[0]*Q1[i][j];
                      for(im=1;im<=mm;im++)
                        {
                        xx+=stencil[im]*(P1[i+im][j]+P1[i-im][j]);
                        zz+=stencil[im]*(Q1[i][j+im]+Q1[i][j-im]);
                        }
                       xx*=dt2x2;
                       zz*=dt2z2;

                       if(flag){ss=s[i][j];}else{ss=0.0;}

                       P0[i][j] = 2.0*P1[i][j] - P0[i][j] 
                                  + xx*vp[i][j]*(1+2*epsilu[i][j]) + zz*vp[i][j] + ss;

                       Q0[i][j] = 2.0*Q1[i][j] - Q0[i][j] 
                                  + xx*vp[i][j]*(1+2* deta [i][j]) + zz*vp[i][j] + ss;
		    }
		 }
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
void read_file(char FN1[],char FN2[],char FN3[],int nx,int nz,float **vv,float **epsilu,float **deta,int npd)
{

		 int i,j;
		
		 FILE *fp1,*fp2,*fp3;
		 if((fp1=fopen(FN1,"rb"))==NULL)printf("Open %s eror!\n",FN1);
		 if((fp2=fopen(FN2,"rb"))==NULL)printf("Open %s eror!\n",FN2);
		 if((fp3=fopen(FN3,"rb"))==NULL)printf("Open %s eror!\n",FN3);
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 fread(&vv[i][j],4,1,fp1);
				 fread(&epsilu[i][j],4,1,fp2);
				 fread(&deta[i][j],4,1,fp3);

			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
}

/**************************func****************************/    
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float mu_v,float **p_cal,int tt)
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
           p_cal[i][j]=0.0;
       }else{}
}


