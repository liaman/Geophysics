//#########################################################
//##         2D Acoustic VTI Medium RTM
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
void backward_corr(int nx,int nz,int nt,int npd,float dx,float dz,float favg,float tmax,float dt,
           float pfac,char FN1[],char FN2[],char FN3[],int ns,int ds,int fs,int zs,
           int is,float **p_cal,float **P_up,float **P_do,float **P_le,float **P_ri,float **Q_up,float **Q_do,float **Q_le,float **Q_ri,
           int mm,int wtype,int hsx,int myid,float **P,float **Q,float **mig_is,float **mig_ns0,
           float ***adcigs_is,float ***adcigs_ns0,int na);
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float mu_v,float **p_cal,int tt);
void laplace_filter(int adj, int nz, int nx, float **in, float **out);
	int i,j,k,is,nx,nz,nt,i_start,i_end,mm,wtype,hsx,na,ia;
	int ns,ds,fs,zs,npd;
	float dx,dz,tmax,dt,pfac,favg;
	int myid,numprocs,Circle_Radius,flag_mu;
        float mu_v;
        float **p_cal,**P_up,**P_do,**P_le,**P_ri,**Q_up,**Q_do,**Q_le,**Q_ri;
        float **P,**Q,**mig_is,**mig_ns0,**mig_ns,***adcigs_is,***adcigs_ns0,***adcigs_ns;


/******** ranks,wavelet,receivers,mute direct **********/
        mm=4;wtype=1;hsx=1;flag_mu=1;na=90;
/******************* dat document **********************/
        char FN1[250]={"vel100100.dat"};
        char FN2[250]={"epsilu100100.dat"};
        char FN3[250]={"deta100100.dat"};
	char FN4[250]={"shot.dat"};
	char FN5[250]={"migration.dat"};
	char FN6[250]={"adcigs.dat"};

/***************** parameters **************************/

        nx=100;         npd=50;      tmax=1.5;
	nz=100;         favg=40;     pfac=1000.0;

 	dx=5.0;
        dz=5.0;

	nt=501;
        dt=1.0;

        ns=10;
        fs=5;
        ds=10;
        zs=1;


        Circle_Radius=15;//for get rid of SV

/************************Loop start**************************/

      FILE *fp4,*fp5,*fp6;
      fp4=fopen(FN4,"wb");
      fp5=fopen(FN5,"wb");
      fp6=fopen(FN6,"wb");

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
   mig_is =alloc2float(nz,nx);
   mig_ns =alloc2float(nz,nx); zero2float(mig_ns, nz,nx);
   mig_ns0=alloc2float(nz,nx); zero2float(mig_ns0,nz,nx);
   adcigs_is=alloc3float(nz,nx,na);   zero3float(adcigs_is,nz,nx,na);
   adcigs_ns0=alloc3float(nz,nx,na);  zero3float(adcigs_ns0,nz,nx,na);
   adcigs_ns=alloc3float(nz,nx,na);   zero3float(adcigs_ns,nz,nx,na);
/********************IS Loop start**************/
   for(is=1+myid;is<=ns;is=is+numprocs)
    {
         zero2float(mig_is, nz,nx);
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


      fseek(fp4,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	  for(j=0;j<nt;j++)
	    fwrite(&p_cal[i][j],4L,1,fp4);


      zero2float(mig_is,nz,nx);
      zero3float(adcigs_is,nz,nx,na);
      backward_corr(nx,nz,nt,npd,dx,dz,favg,tmax,dt,pfac,FN1,FN2,FN3,
            ns,ds,fs,zs,is,p_cal,P_up,P_do,P_le,P_ri,Q_up,Q_do,Q_le,Q_ri,
            mm,wtype,hsx,myid,P,Q,mig_is,mig_ns0,adcigs_is,adcigs_ns0,na);


      MPI_Barrier(MPI_COMM_WORLD);
    } //IS loop end
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Reduce(mig_ns0[0], mig_ns[0], nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(adcigs_ns0[0][0], adcigs_ns[0][0], na*nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
   laplace_filter(1, nz, nx, mig_ns, mig_ns0);
	for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++)
	    fwrite(&mig_ns0[i][j],4L,1,fp5);

      //for(ia=0;ia<na;ia++) laplace_filter(1, nz, nx, adcigs_ns[ia], adcigs_ns0[ia]);
      for(i=0;i<nx;i++)
       for(ia=0;ia<na;ia++)
	  for(j=0;j<nz;j++)
	    fwrite(&adcigs_ns[ia][i][j],4L,1,fp6);

/**********************************************/
   if(myid==0)printf("---   The forward is over    \n");
   if(myid==0)printf("---   Complete!!!!!!!!! \n");

   fclose(fp4);
   fclose(fp5);
   fclose(fp6);

   free2float(P);      free2float(Q);
   free2float(p_cal);
   free2float(P_up);   free2float(Q_up);
   free2float(P_do);   free2float(Q_do);
   free2float(P_le);   free2float(Q_le);
   free2float(P_ri);   free2float(Q_ri);
   free2float(mig_is);   free2float(mig_ns);   free2float(mig_ns0);
   free3float(adcigs_is);   free3float(adcigs_ns);   free3float(adcigs_ns0);
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

	    /*  FILE *fp1,*fp2;
              fp1=fopen("snap-sv.dat","wb");
              fp2=fopen("snap1-sv.dat","wb");	*/





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

       /*    if((is==1)&&(it%20==0)&&(myid==0))
           {
              fseek(fp1,(int)(it/20)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&P1[i][j],4L,1,fp1);


              fseek(fp2,(int)(it/20)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&Q1[i][j],4L,1,fp2);
           }*/
     }//it loop end
	  for(j=0;j<nz;j++)
	    for(i=0;i<nx;i++)
               {P[i][j]=P0[i+npd][j+npd];Q[i][j]=Q0[i+npd][j+npd];}
/**********************close************************/
         // fclose(fp1);fclose(fp2);
/**********************free*************************/
          free2float(P0);  free2float(P1);  free2float(Q0);  free2float(Q1);

          free2float(s);

}
/***********************************func********************************************/
void backward_corr(int nx,int nz,int nt,int npd,float dx,float dz,float favg,float tmax,float dt,
           float pfac,char FN1[],char FN2[],char FN3[],int ns,int ds,int fs,int zs,
           int is,float **p_obs,float **P_up,float **P_do,float **P_le,float **P_ri,float **Q_up,float **Q_do,float **Q_le,float **Q_ri,
           int mm,int wtype,int hsx,int myid,float **P,float **Q,float **mig_is,float **mig_ns0,float ***adcigs_is,float ***adcigs_ns0,int na)
/*******************************************************
Function for VTI backward propagation,2017.1.31

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

	  int i,j,ia,im;
	  int it;
	  float t,temp;


	  float **vp,**illum;

	  float **s_P0;   float **s_P1;
          float **s_Q0;   float **s_Q1;

	  float **g_P0;   float **g_P1;
          float **g_Q0;   float **g_Q1;

	  float **epsilu,**deta;

        float Ssx, Ssz, Sgx, Sgz, b1, b2, a, s_u=0.0, s_w=0.0, g_u=0.0, g_w=0.0;

	illum=alloc2float(nz,nx); zero2float(illum,nz,nx);

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
/********************************************************/
	 s_P0=alloc2float(nz+2*npd,nx+2*npd);
	 s_P1=alloc2float(nz+2*npd,nx+2*npd);
         s_Q0=alloc2float(nz+2*npd,nx+2*npd);
	 s_Q1=alloc2float(nz+2*npd,nx+2*npd);
	 g_P0=alloc2float(nz+2*npd,nx+2*npd);
	 g_P1=alloc2float(nz+2*npd,nx+2*npd);
         g_Q0=alloc2float(nz+2*npd,nx+2*npd);
	 g_Q1=alloc2float(nz+2*npd,nx+2*npd);

       dt=dt/1000;
/*******************zero************************/


           zero2float(s_P0,nz+2*npd,nx+2*npd);
           zero2float(s_P1,nz+2*npd,nx+2*npd);
           zero2float(s_Q0,nz+2*npd,nx+2*npd);
           zero2float(s_Q1,nz+2*npd,nx+2*npd);
           zero2float(g_P0,nz+2*npd,nx+2*npd);
           zero2float(g_P1,nz+2*npd,nx+2*npd);
           zero2float(g_Q0,nz+2*npd,nx+2*npd);
           zero2float(g_Q1,nz+2*npd,nx+2*npd);

           pad_vv(nx,nz,npd,vp);


	  for(i=0;i<=nx+2*npd-1;i++)
		for(j=0;j<=nz+2*npd-1;j++)
		   vp[i][j]=vp[i][j]*vp[i][j];

	    /*  FILE *fp1,*fp2;
              fp1=fopen("chonggouP.dat","wb");
              fp2=fopen("chonggouQ.dat","wb");	*/


	 for(j=0;j<nz;j++)
	    for(i=0;i<nx;i++)
               {s_P1[i+npd][j+npd]=P[i][j];s_Q1[i+npd][j+npd]=Q[i][j];}
    for(it=nt-1; it>=0; it--)
     {//it loop start

      if(it%50==0&&myid==0)printf("---   is=%d, it=%4d, \n",is,it);
          for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			temp=s_P0[i][j];s_P0[i][j]=s_P1[i][j];s_P1[i][j]=temp;
			temp=s_Q0[i][j];s_Q0[i][j]=s_Q1[i][j];s_Q1[i][j]=temp;
		}
	   }
      for(i=0;i<nx;i++)
        {
         s_P1[npd+i][npd]      =  P_up[i][it];
         s_P1[npd+i][npd+nz-1] =  P_do[i][it];
         s_Q1[npd+i][npd]      =  Q_up[i][it];
         s_Q1[npd+i][npd+nz-1] =  Q_do[i][it];
        }
      for(j=0;j<nz;j++)
        {
         s_P1[npd][npd+j]      =  P_le[j][it];
         s_P1[npd+nx-1][npd+j] =  P_ri[j][it];
         s_Q1[npd][npd+j]      =  Q_le[j][it];
         s_Q1[npd+nx-1][npd+j] =  Q_ri[j][it];
        }
        fd_vti_2d_PQ(0,nx,nz,dt,dx,dz,mm,s_P0,s_Q0,s_P1,s_Q1,vp,npd,deta,epsilu,
                  fs,ds,zs,is,1,NULL);
        absorb_bndr(s_P0,s_P1,s_Q0,s_Q1,nx,nz,npd,-0.1);

      /*  if((is==1)&&(it%20==0)&&(myid==0))
           {
              fseek(fp1,(int)(it/20)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&s_P1[i][j],4L,1,fp1);


              fseek(fp2,(int)(it/20)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&s_Q1[i][j],4L,1,fp2);   }*/
      for(i=0;i<nx;i++)
        {
         g_P1[npd+i][npd]      =  p_obs[i][it];
         g_Q1[npd+i][npd]      =  p_obs[i][it];
        }
        fd_vti_2d_PQ(0,nx,nz,dt,dx,dz,mm,g_P0,g_Q0,g_P1,g_Q1,vp,npd,deta,epsilu,
                  fs,ds,zs,is,1,NULL);
        absorb_bndr(g_P0,g_P1,g_Q0,g_Q1,nx,nz,npd,-0.1);
          for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			temp=g_P0[i][j];g_P0[i][j]=g_P1[i][j];g_P1[i][j]=temp;
			temp=g_Q0[i][j];g_Q0[i][j]=g_Q1[i][j];g_Q1[i][j]=temp;
		}
	   }
     /*   if((is==1)&&(it%20==0)&&(myid==0))
           {
              fseek(fp1,(int)(it/20)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&g_P1[i][j],4L,1,fp1);


              fseek(fp2,(int)(it/20)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&g_Q1[i][j],4L,1,fp2);
          }   */

        for(i=npd;i<nx+npd;i++)
           for(j=npd;j<nz+npd;j++)
           {
              mig_is[i-npd][j-npd]+=s_P1[i][j]*g_P1[i][j];
              illum[i-npd][j-npd]+=pow(s_P1[i][j],2);

                      s_u=stencil[0]*s_P1[i][j];
                      s_w=stencil[0]*s_Q1[i][j];
                      g_u=stencil[0]*g_P1[i][j];
                      g_w=stencil[0]*g_Q1[i][j];
                      for(im=1;im<=mm;im++)
                        {
                        s_u+=stencil[im]*(s_P1[i+im][j]+s_P1[i-im][j]);
                        s_w+=stencil[im]*(s_Q1[i][j+im]+s_Q1[i][j-im]);
                        g_u+=stencil[im]*(g_P1[i+im][j]+g_P1[i-im][j]);
                        g_w+=stencil[im]*(g_Q1[i][j+im]+g_Q1[i][j-im]);
                        }
                       s_u/=dx;
                       s_w/=dz;
                       g_u/=dx;
                       g_w/=dz;

              Ssx=-s_P1[i][j]*s_u;
              Ssz=-s_Q1[i][j]*s_w;
              Sgx= g_P1[i][j]*g_u;
              Sgz= g_Q1[i][j]*g_w;

              b1=Ssz*Ssz+Ssx*Ssx;
              b2=Sgz*Sgz+Sgx*Sgx;

              a= 0.5*acosf( (Ssx*Sgx+Ssz*Sgz) / (sqrtf(b1*b2)*(1+0.0)) );

              ia=(int)(  a*180/pi  );

              if ((ia>=0)&&(ia<na))
                 adcigs_is[ia][i-npd][j-npd]+=(g_P1[i][j]*s_P1[i][j])*pow(cos(ia*pi/180.0),3);
           }

     }//it loop end

        for(i=0;i<nx;i++)
           for(j=0;j<nz;j++)
           {
              mig_ns0[i][j]+=mig_is[i][j]/illum[i][j];
             for(ia=0;ia<na;ia++)
              adcigs_ns0[ia][i][j]+=adcigs_is[ia][i][j];
           }

/**********************close************************/
        //  fclose(fp1);fclose(fp2);
/**********************free*************************/
          free2float(s_P0);  free2float(s_P1);  free2float(s_Q0);  free2float(s_Q1);
          free2float(g_P0);  free2float(g_P1);  free2float(g_Q0);  free2float(g_Q1);

}
/************************************func***************************************/
void laplace_filter(int adj, int nz, int nx, float **in, float **out)
/*< linear operator, come from Madagascar Mlaplac2>*/
{
    int iz,ix;
    for (ix=0; ix < nx; ix++)for (iz=0; iz < nz; iz++) out[ix][iz]=0.0;
    for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
	    if (iz > 0) {
		if (adj) {
		    out[ix][iz-1] -= in[ix][iz];//out[j-1] -= in[j];
		    out[ix][iz]   += in[ix][iz];//out[j]   += in[j];
		} else {
		    out[ix][iz] += in[ix][iz] - in[ix][iz-1];//out[j] += in[j] - in[j-1];
		}
	    }
	    if (iz < nz-1) {
		if (adj) {
		    out[ix][iz+1] -= in[ix][iz];//out[j+1] -= in[j];
		    out[ix][iz]   += in[ix][iz];//out[j]   += in[j];
		} else {
		    out[ix][iz] += in[ix][iz] - in[ix][iz+1];//out[j] += in[j] - in[j+1];
		}
	    }

	    if (ix > 0) {
		if (adj) {
		    out[ix-1][iz] -= in[ix][iz];//out[j-nz] -= in[j];
		    out[ix][iz]    += in[ix][iz];//out[j]    += in[j];
		} else {
		    out[ix][iz] += in[ix][iz] - in[ix-1][iz];//out[j] += in[j] - in[j-nz];
		}
	    }
	    if (ix < nx-1) {
		if (adj) {
		    out[ix+1][iz] -= in[ix][iz];//out[j+nz] -= in[j];
		    out[ix][iz]    += in[ix][iz];//out[j]    += in[j];
		} else {
		    out[ix][iz] += in[ix][iz] - in[ix+1][iz];//out[j] += in[j] - in[j+nz];
		}
	    }
	}
    }
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


