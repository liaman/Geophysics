//#############################################################a
//##s
//##s       2D Acoustic VTI Medium RTM  & Pick Up ADCIG
//##s
//##s----------------------------------------------------------
//##s    PML , P&SV , VTI , Acoustic , RTM , -SV , Adcig,
//##s    Migration , Both sides receive , SU or dat(v,e,d) , 
//##s    Adcig( poynting , s-r , offset-domain ),
//##s 
//##s----------------------------------------------------------
//##s                      Rong Tao
//##s                      2016
//############################################################a
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "mpi.h"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/alloc.h"
#include "/home/Toa/hc/complex.c"
#include "/home/Toa/hc/complex.h"
#include "/home/Toa/hc/bhdr.h"
#include "/home/Toa/hc/hdr.h"
#include "/home/Toa/hc/fft.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/mute_direct.h"
#include "/home/Toa/hc/cjbsegy.h"
/********* SEG-Y header *********/
typedef cjbsegy segy;
/********** SU & SEG-Y **********/
#define SU_NKEYS        80      /* Number of key header words           */
#define HDRBYTES        240     /* Bytes in the trace header            */
#define EBCBYTES        3200    /* Bytes in the card image EBCDIC block */
#define BNYBYTES        400     /* Bytes in the binary coded block      */
#define SEGY_HDRBYTES   240     /* Bytes in the tape trace header       */
#define SEGY_NKEYS      71      /* Number of mandated header fields     */
#define BHED_NKEYS      27      /* Number of mandated binary fields     */

#define pi 3.1415926535898
/********* SEG-Y header *********/
       Y_3200   y3200;
       bhed     head_400;
       cjbsegy  tr, vtr;
/********* SEG-Y header *********/
/********** SU function *********/
void swap_short_2(short *tni2);
void swap_u_short_2(unsigned short *tni2);
void swap_int_4(int *tni4);
void swap_u_int_4(unsigned int *tni4);
void swap_long_4(long *tni4);
void swap_u_long_4(unsigned long *tni4);
void swap_float_4(float *tnf4);
void swap_double_8(double *tndd8);
void swaphval(segy *tr, int index);
/********* SU function *********/

/* MAIN */
main(int argc,char *argv[])
{
/*************************************** function ********************************************/
void model_vti_get_boundry(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
             float vdx,float vdz,float favg,float tmax,float dt,float dtout,float pfac,
             float **vp,float **epsilu,float **deta,float **rho,
             int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,int is,
             float **p_cal_x,float **p_cal_z,
             float **p_top_x,float **p_bottom_x,float **p_left_x,float **p_right_x,
             float **p_top_z,float **p_bottom_z,float **p_left_z,float **p_right_z,int _Circle_,
             int mm,int wtype,int hsx,int myid,float *mu_v,int flag_snap,int seismic);
void mute_directwave(int flag_mu,int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs_sxd,int ds_sxd,int zs_sxd,int is,
                     float mu_v,float **p_cal,int tt);
void RTM_corr_adcig(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
           float vdx,float vdz,float favg,float tmax,float dt,float dtout,
           float pfac,float **vp,float **epsilu,float **deta,float **rho,char FN5[],
           int ns_sxd,int ds_sxd,int fs_sxd,int ds_initial,int fs_initial,int zs_sxd,int is,int flag_cdp,
           float **p_top_x,float **p_bottom_x,float **p_left_x,float **p_right_x,////////////////////////////////
           float **p_top_z,float **p_bottom_z,float **p_left_z,float **p_right_z,float **p_obs_x,float **p_obs_z,
           int mm,int wtype,int hsx,int myid,float **mig_is,float **mig_ns0,
           float ***adcig_is,float ***adcig_ns0,float ***Ixhz_is,float ***Ixhz_ns0,int nh,
           int flag_snap,int seismic,int flag_adcig);
void smooth1float(float *v,int r,int n);
void smooth2float(int nx,int rx,int nz,int rz,float **v);
void pad_vv(int nx,int nz,int npd,float **ee);
void read_v_e_d_r(char FN1[],char FN2[],char FN3[],int nx,int nz,float **vv,float **epsilu,float **deta,
               float **rho0,int npd,int seismic);

    /*************************** parameter statement ***********************/

	int i,j,k,l,m,ih,is,nx,nz,nt,vnx,vnz,i_start,i_end,mm,wtype,hsx,ia,nxs,flag_cdp,flag_adcig,nh;
	int ns_sxd,ds_sxd,fs_sxd,zs_sxd,fs,ds,npd;
	float dx,dz,vdx,vdz,tmax,dt,dtout,pfac,favg;
	int myid,numprocs,_Circle_,flag_mu,flag_snap,seismic,nangle,dangle,fangle;
      float mu_v;

      /***************** wave float **************/

      float **p_cal_x,**p_cal_z,**p_top_x,**p_bottom_x,**p_left_x,**p_right_x;
      float **p_top_z,**p_bottom_z,**p_left_z,**p_right_z;

      /***************** image float *************/

      float **mig_is,**mig_ns,**mig_ns0,***adcig_is,***adcig_ns,***adcig_ns0,**adcig2d;
      float ***Ixhz_ns0,***Ixhz_ns,***Ixhz_is;

	float **vp,**rho,**deta,**epsilu;
	float **vps,**rhos,**detas,**epsilus;


//a###########################################################################
//####                 input the parameter and confirm                    ####
//a###########################################################################
/************************ dat document **********************/

/* Input velocity        */char FN1[250]={"vel_1801_862.dat"};
/* Input epsilu          */char FN2[250]={"eps_1801_862.dat"};
/* Input deta            */char FN3[250]={"del_1801_862.dat"};
/* Cal data              */char FN4[250]={"shot_cal.dat"};
/* Obs data              */char FN5[250]={"shot_obs.dat"};
/* Migration             */char FN6[250]={"mig_ns.dat"};
/* IS tempor migration   */char FN7[250]={"mig_is_tempor.dat"};
/* Adcig initial         */char FN8[250]={"adcig.dat"};
/* Adcig smooth          */char FN9[250]={"adcig_smooth.dat"};
/* Migration adcig stack */char FN10[250]={"mig_stack_adcig.dat"};
/* Half offset I(x,h,z)  */char FN11[250]={"Image_x_h_z.dat"};

/*************************** type ***************************/

/* Wavelet type (1-3)       */ wtype=1;
/* 1-mute , 0->don't mute   */ flag_mu=1;
/* 1-snap , 0->nosnap       */ flag_snap=0;
/* 1.dat,  2.su  (v,e,d)    */ seismic=1;

/*************************** cdp ****************************/

/* Activate "nxs"(=1)       */ flag_cdp=1;  /* 1-activate ;0-invaliable */
/* CDP of each shot         */ nxs=501;     /* Must be odd number */

/*************************** cdp ****************************/

/* choice adcig type        */ flag_adcig=1;  /* 1-poynting */
                                              /* 2-source-receivers domain */
                                              /* 3-half-offset-domain */
/*************************** wave ***************************/

              hsx=1;mm=4;npd=20;_Circle_=15;

/******************** observation system ********************/

              favg=20;    pfac=1000.0;

              nx=1801;     dx=10.0;
              nz=862;      dz=5.0;

              tmax=6.5;
              dt=0.5;//ms
              nt=11001;

              ns_sxd=250;
              fs_sxd=255;
              ds_sxd=5;
              zs_sxd=1;

              nangle=90;
              fangle=0;
              dangle=1;

/*************************v****************************/

      vdz=dz;vdx=dx;vnx=nx;vnz=nz;dtout=dt;

/************************FILE**************************/

      FILE *fp4,*fp6,*fp7,*fp8,*fp9,*fp10,*fp11;
      fp4=fopen(FN4,"wb");    /* Cal data              */
      fp6=fopen(FN6,"wb");    /* Migration             */
      fp7=fopen(FN7,"wb");    /* IS tempor migration   */
      fp8=fopen(FN8,"wb");    /* Adcig initial         */
      fp9=fopen(FN9,"wb");    /* Adcig smooth          */
      fp10=fopen(FN10,"wb");  /* Migration adcig stack */
   if(flag_adcig==3)
      fp11=fopen(FN11,"wb");  /* Image half-offset     */

/************************** read_file ******************************/

       vp = alloc2float(nz+2*npd,nx+2*npd);   zero2float(vp,nz+2*npd,nx+2*npd);
      rho = alloc2float(nz+2*npd,nx+2*npd);   zero2float(rho,nz+2*npd,nx+2*npd);
   epsilu = alloc2float(nz+2*npd,nx+2*npd);   zero2float(epsilu,nz+2*npd,nx+2*npd);
     deta = alloc2float(nz+2*npd,nx+2*npd);   zero2float(deta,nz+2*npd,nx+2*npd);

       read_v_e_d_r(FN1,FN2,FN3,vnx,vnz,vp,epsilu,deta,rho,npd,seismic);

        pad_vv(nx,nz,npd,vp);
        pad_vv(nx,nz,npd,rho);
        pad_vv(nx,nz,npd,epsilu);
        pad_vv(nx,nz,npd,deta);

   /** determine NXS and nh **/
     if(flag_cdp==1){
        nxs=nxs;
        nh=nxs/2+1;
     }else{
        nxs=nx;
        nh=0;
       }

       vps = alloc2float(nz+2*npd,nxs+2*npd);   zero2float(vps,nz+2*npd,nxs+2*npd);
      rhos = alloc2float(nz+2*npd,nxs+2*npd);   zero2float(rhos,nz+2*npd,nxs+2*npd);
   epsilus = alloc2float(nz+2*npd,nxs+2*npd);   zero2float(epsilus,nz+2*npd,nxs+2*npd);
     detas = alloc2float(nz+2*npd,nxs+2*npd);   zero2float(detas,nz+2*npd,nxs+2*npd);

/********************* alloc ***********************/
	p_cal_x=alloc2float(nt,nxs);
	p_cal_z=alloc2float(nt,nxs);

	p_top_x=alloc2float(nt,nxs);  p_bottom_x=alloc2float(nt,nxs);
	p_top_z=alloc2float(nt,nxs);	p_bottom_z=alloc2float(nt,nxs);
	
	p_left_x=alloc2float(nt,nz);	p_right_x=alloc2float(nt,nz);
	p_left_z=alloc2float(nt,nz);  p_right_z=alloc2float(nt,nz);

/********************* image alloc *************************/
      mig_is=alloc2float(nz,nxs);        zero2float(mig_is,nz,nxs);
      adcig_is=alloc3float(nz,nxs,90);   zero3float(adcig_is,nz,nxs,90);
      /* half-offset domain image alloc */
   if(flag_adcig==3)
    {  Ixhz_is=alloc3float(nz,nh,nxs);       zero3float(Ixhz_is,nz,nh,nxs); }

      mig_ns=alloc2float(nz,nx);         zero2float(mig_ns,nz,nx);
      mig_ns0=alloc2float(nz,nx);        zero2float(mig_ns0,nz,nx);
      adcig_ns=alloc3float(nz,nx,90);    zero3float(adcig_ns,nz,nx,90);
      adcig_ns0=alloc3float(nz,nx,90);   zero3float(adcig_ns0,nz,nx,90);
      
      /* half-offset domain image alloc */
   if(flag_adcig==3)
   {  Ixhz_ns=alloc3float(nz,nh,nx);        zero3float(Ixhz_ns,nz,nh,nx);
      Ixhz_ns0=alloc3float(nz,nh,nx);       zero3float(Ixhz_ns0,nz,nh,nx);}

      adcig2d=alloc2float(nz,90);        zero2float(adcig2d,nz,90);

/*******************MPI************************/
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

/*******************MPI***********************/
      MPI_Barrier(MPI_COMM_WORLD);
   if(myid==0)
    {
        printf("--------------------------------------------------------\n");
        printf("    ####################################################\n");
        printf("    ###         Please confirm the parameter !       ###\n");
        printf("    ####################################################\n");
        printf("    ###                                              ###\n");
        printf("    ###    mm=%d, wtype=%d, hsx=%d, flag_mu=%d\n",mm,wtype,hsx,flag_mu);
        printf("    ###\n");
        printf("    ###    nx=vnx=%3d, dx=vdx=%.1f\n",nx,dx);
        printf("    ###    nz=vnz=%3d, dz=vdz=%.1f\n",nz,dz);
        printf("    ###\n");
        printf("    ###    npd=%3d,      favg=%.1f,   pfac=%.1f\n",npd,favg,pfac);
        printf("    ###    tmax=%.2f(s), dt=%.2f(ms), nt=%d\n",tmax,dt,nt);
        printf("    ###\n");
        printf("    ###    ns=%3d\n",ns_sxd);
        printf("    ###    fs=%3d\n",fs_sxd);
        printf("    ###    ds=%3d\n",ds_sxd);
        printf("    ###    zs=%3d\n",zs_sxd);
        printf("    ###\n");
        printf("    ###    _Circle_=%3d\n",_Circle_);
        printf("    ###\n"); 
        printf("    ###    numprocs= %2d\n",numprocs);
        printf("    ###                                              ###\n");
        printf("    ####################################################\n");
        //system("pause");
    }if( ( flag_cdp==1 ) && ( nx<nxs ) ){
        if(myid==0)
           {
          printf("\n\n\n");
          printf("    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
          printf("    $$$   nx must more than nxs!   $$$\n");
          printf("    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");printf("\n\n\n");
           }
        exit(0);
    }if( ( flag_cdp==1 ) && ( ( fs_sxd<(nxs/2+1) ) || ( (fs_sxd+(ns_sxd-1)*ds_sxd+nxs/2)>nx ) ) ){
        if(myid==0)
           {
          printf("\n\n\n");
          printf("    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
          printf("    $$$   Receivers location out of model boundary !   $$$\n");
          printf("    $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");printf("\n\n\n");
           }
        exit(0);
     }
/************************IS Loop start**************************/
/************************IS Loop start**************************/
   for(is=1+myid;is<=ns_sxd;is+=numprocs)
    {
      if(myid==0)
        {
         printf("--------------------------------------------------------\n");
         printf("--------------------------------------------------------\n");
         printf("---   IS========%d  \n",is);
         printf("---   The forward is start  !  \n");
        }
      zero2float(p_cal_x,nt,nxs);  zero2float(p_cal_z,nt,nxs);
      zero2float(p_top_x,nt,nxs);  zero2float(p_bottom_x,nt,nxs);
      zero2float(p_top_z,nt,nxs);  zero2float(p_bottom_z,nt,nxs);
      zero2float(p_left_x,nt,nz);  zero2float(p_right_x,nt,nz);
      zero2float(p_left_z,nt,nz);  zero2float(p_right_z,nt,nz);
 /* determine IS vp,rho,deta,epsilu */
   if(flag_cdp==1){
     k=fs_sxd+(is-1)*ds_sxd;
     for(i=k-nxs/2+npd;i<=k+nxs/2+npd;i++)
       for(j=npd;j<nz+npd;j++)
         {
         vps[i-k+nxs/2][j]=vp[i][j];
         rhos[i-k+nxs/2][j]=rho[i][j];
         detas[i-k+nxs/2][j]=deta[i][j];
         epsilus[i-k+nxs/2][j]=epsilu[i][j];
         }
     ds=0;
     fs=nxs/2+1;
    }else{
     for(i=npd;i<nx+npd;i++)
       for(j=npd;j<nz+npd;j++)
         {
         vps[i][j]=vp[i][j];
         rhos[i][j]=rho[i][j];
         detas[i][j]=deta[i][j];
         epsilus[i][j]=epsilu[i][j];
         }
      fs=fs_sxd;
      ds=ds_sxd;
     }
        pad_vv(nxs,nz,npd,vps);
        pad_vv(nxs,nz,npd,rhos);
        pad_vv(nxs,nz,npd,epsilus);
        pad_vv(nxs,nz,npd,detas);
/**************** model vti forwarding**************/
      model_vti_get_boundry(nxs,nz,nxs,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,
                           vps,epsilus,detas,rhos,ns_sxd,ds,fs,zs_sxd,is,
                           p_cal_x,p_cal_z,
                           p_top_x,p_bottom_x,p_left_x,p_right_x,
                           p_top_z,p_bottom_z,p_left_z,p_right_z,
                           _Circle_,mm,wtype,hsx,myid,&mu_v,flag_snap,seismic);
      if(myid==0)printf("---   The forward is over  !  \n");
        /* make the vp back to real */
	  for(i=0;i<=nxs+2*npd-1;i++)
	   {
		for(j=0;j<=nz+2*npd-1;j++)
		{
               rhos[i][j]=1.0/rhos[i][j];
		   vps[i][j]=sqrtf(vps[i][j]/rhos[i][j]);
		   
		}
	   }
/**************** mute direct wave **************/
      mute_directwave(flag_mu,nxs,nt,dt,favg,dx,dz,fs,ds,zs_sxd,is,mu_v,p_cal_x,285);
      mute_directwave(flag_mu,nxs,nt,dt,favg,dx,dz,fs,ds,zs_sxd,is,mu_v,p_cal_z,285);

/**************** output p_cal **************/
      fseek(fp4,(is-1)*nxs*nt*4L,0);
	for(i=0;i<nxs;i++)
	  for(j=0;j<nt;j++)
	    fwrite(&p_cal_x[i][j],4L,1,fp4);

/**************** RTM correlation **************/
      if(myid==0)
        {
         printf("--------------------------------------------\n");
         printf("---\n");
         printf("---   The RTM correlation is start  !  \n");
        }
      MPI_Barrier(MPI_COMM_WORLD);
      zero2float(mig_is,nz,nxs);
      zero3float(adcig_is,nz,nxs,90);
   if(flag_adcig==3)
    {  zero3float(Ixhz_is,nz,nh,nxs); }

      RTM_corr_adcig(nxs,nz,nxs,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,
                     vps,epsilus,detas,rhos,FN5,ns_sxd,ds,fs,ds_sxd,fs_sxd,zs_sxd,is,flag_cdp,                                     
                     p_top_x,p_bottom_x,p_left_x,p_right_x,
                     p_top_z,p_bottom_z,p_left_z,p_right_z,p_cal_x,p_cal_z, /* p_cal->p_obs */
                     mm,wtype,hsx,myid,mig_is,mig_ns0,adcig_is,adcig_ns0,Ixhz_is,Ixhz_ns0,nh,
                     flag_snap,seismic,flag_adcig);
/**************** output tempor migration **************/
      fseek(fp7,(is-1)*nxs*nz*4L,0);
	for(i=0;i<nxs;i++)
	  for(j=0;j<nz;j++)
	    fwrite(&mig_is[i][j],4L,1,fp7);

      MPI_Barrier(MPI_COMM_WORLD);
    }//is loop ending
/*****************IS Loop end******************/

   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Reduce(mig_ns0[0], mig_ns[0], nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
   MPI_Reduce(adcig_ns0[0][0], adcig_ns[0][0], 90*nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
 if(flag_adcig==3)
   MPI_Reduce(Ixhz_ns0[0][0], Ixhz_ns[0][0], nh*nx*nz, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);


   if(myid==0)
    {
      for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++)
	    fwrite(&mig_ns[i][j],4L,1,fp6);

      for(i=0;i<nx;i++)
       for(ia=0;ia<90;ia++)
	  for(j=0;j<nz;j++)
	    fwrite(&adcig_ns[ia][i][j],4L,1,fp8);

   if(flag_adcig==3)
      for(i=0;i<nx;i++)
       for(ih=0;ih<nh;ih++)
	  for(j=0;j<nz;j++)
	    fwrite(&Ixhz_ns[i][ih][j],4L,1,fp11);

/************** smooth the adcig & stack ***************/
    zero2float(mig_ns0,nz,nx);
     for(i=0;i<nx;i++)
      {
       for(ia=0;ia<90;ia++)
	   for(j=0;j<nz;j++)
              adcig2d[ia][j]=adcig_ns[ia][i][j];
       smooth2float(90,10,nz,0,adcig2d);
       for(ia=0;ia<90;ia++)
	   for(j=0;j<nz;j++)
            {
              adcig_ns[ia][i][j]=adcig2d[ia][j];
              mig_ns0[i][j]+=adcig2d[ia][j];
            }
      }
      for(i=0;i<nx;i++)
       for(ia=0;ia<90;ia++)
	  for(j=0;j<nz;j++)
	    fwrite(&adcig_ns[ia][i][j],4L,1,fp9);
      for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++)
	    fwrite(&mig_ns0[i][j],4L,1,fp10);
    }
   MPI_Barrier(MPI_COMM_WORLD);
   if(myid==0)printf("---   The RTM correlation is over !   \n");
   if(myid==0)printf("---   Complete!!!!!!!!! \n");


/************** fclose ***************/
   fclose(fp4);
   fclose(fp6);
   fclose(fp7);
   fclose(fp8);
   fclose(fp9);
   fclose(fp10);
   if(flag_adcig==3)
     fclose(fp11);
/************** free ***************/
   free2float(rho);   free2float(vp);  free2float(epsilu);  free2float(deta);
   free2float(rhos);  free2float(vps); free2float(epsilus); free2float(detas);

   free2float(p_cal_x);      free2float(p_cal_z);
   free2float(p_top_x);      free2float(p_top_z);
   free2float(p_bottom_x);   free2float(p_bottom_z);
   free2float(p_left_x);     free2float(p_left_z);
   free2float(p_right_x);    free2float(p_right_z);

   free2float(mig_is);       
   free2float(mig_ns);
   free2float(mig_ns0);

   free3float(adcig_is);
   free3float(adcig_ns);
   free3float(adcig_ns0);
   free2float(adcig2d);

   if(flag_adcig==3)
  {   free3float(Ixhz_is);
      free3float(Ixhz_ns);
      free3float(Ixhz_ns0); }
/******************MPI************************/
   MPI_Finalize();
}
/***********************************func********************************************/
void model_vti_get_boundry(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
           float vdx,float vdz,float favg,float tmax,float dt,float dtout,float pfac,
           float **vp,float **epsilu,float **deta,float **rho,
           int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,int is,
           float **p_cal_x,float **p_cal_z,
           float **p_top_x,float **p_bottom_x,float **p_left_x,float **p_right_x,
           float **p_top_z,float **p_bottom_z,float **p_left_z,float **p_right_z,
           int _Circle_,int mm,int wtype,int hsx,int myid,float *mu_v,int flag_snap,int seismic)
/*****************************************************a**
Function for VTI medium modeling,2016.9.24

 Ps:  the function of modeling following:

          du/dt=1/rho*dp/dx ,
          dw/dt=1/rho*dq/dz ,
          dp/dt=rho*vpx^2*du/dx+rho*vp0*vpn*dw/dz ,
          dq/dt=rho*vp0*vpn*du/dx+rho*vp0^2*dw/dz ,
                     vpx^2=vp0^2*(1+2*epsilu);
                     vpn^2=vp0^2*(1+2*deta);

                                               Rong Tao
******************************************************a**/
{
void cal_c(int mm,float c[]);
void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
                  float favg,float **s,int wtype,int npd,int is,int ds_sxd);
void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                  float **u0,float **w0,float **txx0,float **tzz0,
                  float **u1,float **w1,float **txx1,float **tzz1,
                  float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2);
void update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
                  float **u0,float **w0,float **txx0,float **tzz0,
                  float **u1,float **w1,float **txx1,float **tzz1,
                  float **s,float **vp,float c[],int npd,
                  float **tpx1,float **tpx0,float **tpz1,float **tpz0,
                  float **tqx1,float **tqx0,float **tqz1,float **tqz0,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                  float **deta,float **epsilu,int fs_sxd,int ds_sxd,int zs_sxd,int is,int _Circle_);
float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,
                  int npd,float tmax,float favg,float dtout,float dt,float **vp,float ndtt,int myid);
void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);

	  int i,j;
	  int ntout,it;
	  float t,ndtt,d0;



	  float **u0;    float **u1;
	  float **w0;    float **w1;

	  float **tpx0;   float **tqx0;
	  float **tpx1;   float **tqx1;
	  float **tpz0;   float **tqz0;
	  float **tpz1;   float **tqz1;

	  float **txx0;   float **txx1;
        float **tzz0;   float **tzz1;

	  float **s;

	  float c[mm];

	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;


        cal_c(mm,c);


/********************************************************/
       *mu_v=vp[1+npd][1+npd]*sqrtf((1+2*epsilu[1+npd][1+npd]));/////////////
/********************************************************/

	 u0=alloc2float(nz+2*npd,nx+2*npd);
	 u1=alloc2float(nz+2*npd,nx+2*npd);
	 w0=alloc2float(nz+2*npd,nx+2*npd);
	 w1=alloc2float(nz+2*npd,nx+2*npd);

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


       d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp,ndtt,myid);
       dt=dt/1000;
/**************************** boundry *******************************/
        float *coffx1;  float *coffx2;  float *coffz1;  float *coffz2;
        float *acoffx1; float *acoffx2; float *acoffz1; float *acoffz2;
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

        initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/***********************************************************/
	  ndtt=(int)ndtt;
/*******************zero************************/
           zero2float(p_cal_x,nt,nx);
           zero2float(p_cal_z,nt,nx);

           zero2float(u0,nz+2*npd,nx+2*npd);
           zero2float(u1,nz+2*npd,nx+2*npd);
           zero2float(w0,nz+2*npd,nx+2*npd);
           zero2float(w1,nz+2*npd,nx+2*npd);

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

           zero2float(s,nz+2*npd,nx+2*npd);






	  for(i=0;i<=nx+2*npd-1;i++)
	   {
		for(j=0;j<=nz+2*npd-1;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
		   rho[i][j]=1.0/rho[i][j];
		}
	   }

	      FILE *fpsnap,*fpsnap1;
          if((is==1)&&(flag_snap))
            {
              fpsnap=fopen("snap_forward_txx.dat","wb");
              fpsnap1=fopen("snap_forward_tzz.dat","wb");
            }


    for(it=0,t=0.0;it<nt;it++,t+=dt)
     {

      if(it%100==0&&myid==0)printf("---   FOR   is===%d   it===%d\n",is,it);

	ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,npd,is,ds_sxd);
      update_vel(nx,nz,npd,mm,dt,dx,dz,u0,w0,txx0,tzz0,
                u1,w1,txx1,tzz1,rho,c,coffx1,coffx2,coffz1,coffz2);
      update_stress(nx,nz,dt,dx,dz,mm,u0,w0,txx0,tzz0,
                u1,w1,txx1,tzz1,s,vp,c,npd,
                tpx1,tpx0,tpz1,tpz0,tqx1,tqx0,tqz1,tqz0,
                acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,
                fs_sxd,ds_sxd,zs_sxd,is,_Circle_);

	  for(i=npd;i<npd+nx;i++)
	  {
		p_cal_x[i-npd][it]     =   txx1[i][npd+hsx-1];//////
		p_cal_z[i-npd][it]     =   tzz1[i][npd+hsx-1];//////
            p_top_x[i-npd][it]     =   txx1[i][npd];
            p_bottom_x[i-npd][it]  =   txx1[i][npd+nz-1];
            p_top_z[i-npd][it]     =   tzz1[i][npd];
            p_bottom_z[i-npd][it]  =   tzz1[i][npd+nz-1];
	  }
	  for(j=npd;j<npd+nz;j++)
	  {
            p_left_x[j-npd][it]    =   txx1[npd][j];
            p_right_x[j-npd][it]   =   txx1[npd+nx-1][j];
            p_left_z[j-npd][it]    =   tzz1[npd][j];
            p_right_z[j-npd][it]   =   tzz1[npd+nx-1][j];
	  }


	  for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			u0[i][j]=u1[i][j];
			w0[i][j]=w1[i][j];

			tpx0[i][j]=tpx1[i][j];
			tpz0[i][j]=tpz1[i][j];
                  tqx0[i][j]=tqx1[i][j];
			tqz0[i][j]=tqz1[i][j];

			txx0[i][j]=txx1[i][j];
                  tzz0[i][j]=tzz1[i][j];
		}
	   }

           if((is==1)&&(it%100==0)&&(flag_snap))
           {
              fseek(fpsnap,(int)(it/100)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&txx1[i][j],4L,1,fpsnap);


              fseek(fpsnap1,(int)(it/100)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&tzz1[i][j],4L,1,fpsnap1);
           }
     }//it loop end
/**********************close************************/
          if((is==1)&&(flag_snap))
            {
          fclose(fpsnap);fclose(fpsnap1);
            }
/**********************free*************************/
          free1float(coffx1);free1float(coffx2);
          free1float(coffz1);free1float(coffz2);
          free1float(acoffx1);free1float(acoffx2);
          free1float(acoffz1);free1float(acoffz2);

          free2float(u0);   free2float(u1);
          free2float(w0);   free2float(w1);

          free2float(txx0);  free2float(txx1);  free2float(tzz0);  free2float(tzz1);

          free2float(tpx0);  free2float(tpx1);  free2float(tpz0);  free2float(tpz1);
          free2float(tqx0);  free2float(tqx1);  free2float(tqz0);  free2float(tqz1);

          free2float(s);     
}
/************************************func***************************************/
void RTM_corr_adcig(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
           float vdx,float vdz,float favg,float tmax,float dt,float dtout,float pfac,
           float **vp,float **epsilu,float **deta,float **rho,char FN5[],
           int ns_sxd,int ds_sxd,int fs_sxd,int ds_initial,int fs_initial,int zs_sxd,int is,int flag_cdp,
           float **p_top_x,float **p_bottom_x,float **p_left_x,float **p_right_x,
           float **p_top_z,float **p_bottom_z,float **p_left_z,float **p_right_z,float **p_obs_x,float **p_obs_z,
           int mm,int wtype,int hsx,int myid,float **mig_is,float **mig_ns0,
           float ***adcig_is,float ***adcig_ns0,float ***Ixhz_is,float ***Ixhz_ns0,int nh,
           int flag_snap,int seismic,int flag_adcig)
/************************************************************a*
  function for RTM
    PS:correlation image condition
        pick up the adcig and stack
                                          Rong Tao
*************************************************************b*/
{
void cal_c(int mm,float c[]);
void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                  float **s_u0,float **s_w0,float **s_txx0,float **s_tzz0,
                  float **s_u1,float **s_w1,float **s_txx1,float **s_tzz1,
                  float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2);
void inv_update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
                  float **s_u0,float **s_w0,float **s_txx0,float **s_tzz0,
                  float **s_u1,float **s_w1,float **s_txx1,float **s_tzz1,
                  float **vp,float c[],int npd,
                  float **s_tpx1,float **s_tpx0,float **s_tpz1,float **s_tpz0,
                  float **s_tqx1,float **s_tqx0,float **s_tqz1,float **s_tqz0,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                  float **deta,float **epsilu,int fs_sxd,int ds_sxd,int zs_sxd,int is);
float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,
                  int npd,float tmax,float favg,float dtout,float dt,float **vp,float ndtt,int myid);
void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);

	  int i,j,k,ih;
	  int ntout,it,ia;
	  float t,ndtt,d0,atan_s,atan_g,angle,sx,sz,gx,gz,b1,b2,a,a1,a2;

/****************** source wavefiled **************/
	  float **s_u0;    float **s_u1;
	  float **s_w0;    float **s_w1;

	  float **s_tpx0;   float **s_tqx0;
	  float **s_tpx1;   float **s_tqx1;
	  float **s_tpz0;   float **s_tqz0;
	  float **s_tpz1;   float **s_tqz1;

	  float **s_txx0;   float **s_txx1;
      float **s_tzz0;   float **s_tzz1;
/***************** reciever wavefiled *************/
	  float **g_u0;    float **g_u1;
	  float **g_w0;    float **g_w1;

	  float **g_tpx0;   float **g_tqx0;
	  float **g_tpx1;   float **g_tqx1;
	  float **g_tpz0;   float **g_tqz0;
	  float **g_tpz1;   float **g_tqz1;

	  float **g_txx0;   float **g_txx1;
        float **g_tzz0;   float **g_tzz1;
/***************** Lighting source ***************/
        float **source_x,**source_z,***sourceI;

/***************** float end  *************/

	  float c[mm];


	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;


        cal_c(mm,c);

/********************************************************/
/********************************************************/

/****************** source wavefiled *********************** reciever wavefiled *************/
	 s_u0=alloc2float(nz+2*npd,nx+2*npd);	 g_u0=alloc2float(nz+2*npd,nx+2*npd);
	 s_u1=alloc2float(nz+2*npd,nx+2*npd);	 g_u1=alloc2float(nz+2*npd,nx+2*npd);
	 s_w0=alloc2float(nz+2*npd,nx+2*npd);	 g_w0=alloc2float(nz+2*npd,nx+2*npd);
	 s_w1=alloc2float(nz+2*npd,nx+2*npd); 	 g_w1=alloc2float(nz+2*npd,nx+2*npd);

	 s_txx0=alloc2float(nz+2*npd,nx+2*npd);	 g_txx0=alloc2float(nz+2*npd,nx+2*npd);
	 s_txx1=alloc2float(nz+2*npd,nx+2*npd);	 g_txx1=alloc2float(nz+2*npd,nx+2*npd);
       s_tzz0=alloc2float(nz+2*npd,nx+2*npd);    g_tzz0=alloc2float(nz+2*npd,nx+2*npd);
	 s_tzz1=alloc2float(nz+2*npd,nx+2*npd);	 g_tzz1=alloc2float(nz+2*npd,nx+2*npd);

	 s_tpx0=alloc2float(nz+2*npd,nx+2*npd);	 g_tpx0=alloc2float(nz+2*npd,nx+2*npd);
	 s_tpx1=alloc2float(nz+2*npd,nx+2*npd);	 g_tpx1=alloc2float(nz+2*npd,nx+2*npd);
	 s_tpz0=alloc2float(nz+2*npd,nx+2*npd);	 g_tpz0=alloc2float(nz+2*npd,nx+2*npd);
	 s_tpz1=alloc2float(nz+2*npd,nx+2*npd);	 g_tpz1=alloc2float(nz+2*npd,nx+2*npd);

	 s_tqx0=alloc2float(nz+2*npd,nx+2*npd);	 g_tqx0=alloc2float(nz+2*npd,nx+2*npd);
	 s_tqx1=alloc2float(nz+2*npd,nx+2*npd);	 g_tqx1=alloc2float(nz+2*npd,nx+2*npd);
	 s_tqz0=alloc2float(nz+2*npd,nx+2*npd);	 g_tqz0=alloc2float(nz+2*npd,nx+2*npd);
	 s_tqz1=alloc2float(nz+2*npd,nx+2*npd);	 g_tqz1=alloc2float(nz+2*npd,nx+2*npd);
/*************zero source wavefiled ****************** zero reciever wavefiled *************/
        zero2float(s_u0,nz+2*npd,nx+2*npd);      zero2float(g_u0,nz+2*npd,nx+2*npd);
        zero2float(s_u1,nz+2*npd,nx+2*npd);      zero2float(g_u1,nz+2*npd,nx+2*npd);
        zero2float(s_w0,nz+2*npd,nx+2*npd);      zero2float(g_w0,nz+2*npd,nx+2*npd);
        zero2float(s_w1,nz+2*npd,nx+2*npd);      zero2float(g_w1,nz+2*npd,nx+2*npd);

        zero2float(s_txx0,nz+2*npd,nx+2*npd);    zero2float(g_txx0,nz+2*npd,nx+2*npd);
        zero2float(s_txx1,nz+2*npd,nx+2*npd);    zero2float(g_txx1,nz+2*npd,nx+2*npd);
        zero2float(s_tzz0,nz+2*npd,nx+2*npd);    zero2float(g_tzz0,nz+2*npd,nx+2*npd);
        zero2float(s_tzz1,nz+2*npd,nx+2*npd);    zero2float(g_tzz1,nz+2*npd,nx+2*npd);

        zero2float(s_tpx0,nz+2*npd,nx+2*npd);    zero2float(g_tpx0,nz+2*npd,nx+2*npd);
        zero2float(s_tpx1,nz+2*npd,nx+2*npd);    zero2float(g_tpx1,nz+2*npd,nx+2*npd);
        zero2float(s_tpz0,nz+2*npd,nx+2*npd);    zero2float(g_tpz0,nz+2*npd,nx+2*npd);
        zero2float(s_tpz1,nz+2*npd,nx+2*npd);    zero2float(g_tpz1,nz+2*npd,nx+2*npd);

        zero2float(s_tqx0,nz+2*npd,nx+2*npd);    zero2float(g_tqx0,nz+2*npd,nx+2*npd);
        zero2float(s_tqx1,nz+2*npd,nx+2*npd);    zero2float(g_tqx1,nz+2*npd,nx+2*npd);
        zero2float(s_tqz0,nz+2*npd,nx+2*npd);    zero2float(g_tqz0,nz+2*npd,nx+2*npd);
        zero2float(s_tqz1,nz+2*npd,nx+2*npd);    zero2float(g_tqz1,nz+2*npd,nx+2*npd);
/************************* Lighting source alloc & zero ******************************/
         source_x=alloc2float(nz,nx);       zero2float(source_x,nz,nx);
         source_z=alloc2float(nz,nx);       zero2float(source_z,nz,nx);

         sourceI=alloc3float(nz,2*nh+1,nx); zero3float(sourceI,nz,2*nh+1,nx);


       d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp,ndtt,myid);
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

        initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/***********************************************************/
	  ndtt=(int)ndtt;

	  for(i=0;i<=nx+2*npd-1;i++)
	   {
		for(j=0;j<=nz+2*npd-1;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
		   rho[i][j]=1.0/rho[i][j];
		}
	   }

	      FILE *fpsnap,*fpsnap1,*fpsnap2,*fpsnap3;
          if((is==1)&&(flag_snap))
            {
              fpsnap=fopen("snap_construction_txx.dat","wb");
              fpsnap1=fopen("snap_construction_tzz.dat","wb");
              fpsnap2=fopen("snap_backpropagation_txx.dat","wb");
              fpsnap3=fopen("snap_backpropagation_tzz.dat","wb");
	    }

    for(it=nt-1;it>=0;it--)
     {

      if(it%100==0&&myid==0)printf("---   RTM   is===%d   it===%d\n",is,it);
/************************ construction waveform start *************************/
      for(i=0;i<nx;i++)
        {
         s_txx0[npd+i][npd]=p_top_x[i][it];
         s_txx0[npd+i][npd+nz-1]=p_bottom_x[i][it];
         s_tzz0[npd+i][npd]=p_top_z[i][it];
         s_tzz0[npd+i][npd+nz-1]=p_bottom_z[i][it];
        }
      for(j=0;j<nz;j++)
        {
         s_txx0[npd][npd+j]=p_left_x[j][it];
         s_txx0[npd+nx-1][npd+j]=p_right_x[j][it];
         s_tzz0[npd][npd+j]=p_left_z[j][it];
         s_tzz0[npd+nx-1][npd+j]=p_right_z[j][it];
        }

      update_vel(nx,nz,npd,mm,dt,dx,dz,s_u0,s_w0,s_txx0,s_tzz0,
                s_u1,s_w1,s_txx1,s_tzz1,rho,c,coffx1,coffx2,coffz1,coffz2);
      inv_update_stress(nx,nz,dt,dx,dz,mm,s_u0,s_w0,s_txx0,s_tzz0,
                s_u1,s_w1,s_txx1,s_tzz1,vp,c,npd,
                s_tpx1,s_tpx0,s_tpz1,s_tpz0,s_tqx1,s_tqx0,s_tqz1,s_tqz0,
                acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,fs_sxd,ds_sxd,zs_sxd,is);
	  for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			s_u0[i][j]=s_u1[i][j];
			s_w0[i][j]=s_w1[i][j];

			s_tpx0[i][j]=s_tpx1[i][j];
			s_tpz0[i][j]=s_tpz1[i][j];
                  s_tqx0[i][j]=s_tqx1[i][j];
			s_tqz0[i][j]=s_tqz1[i][j];

			s_txx0[i][j]=s_txx1[i][j];
                  s_tzz0[i][j]=s_tzz1[i][j];
		}
	   }
/************************ construction waveform end ************************/
/************************** back propagation start *************************/
      for(i=0;i<nx;i++)
        {
         g_txx0[npd+i][npd]=p_obs_x[i][it];/////////////////
         g_tzz0[npd+i][npd]=p_obs_z[i][it];/////////////////
        }
      update_vel(nx,nz,npd,mm,dt,dx,dz,g_u0,g_w0,g_txx0,g_tzz0,
                g_u1,g_w1,g_txx1,g_tzz1,rho,c,coffx1,coffx2,coffz1,coffz2);
      inv_update_stress(nx,nz,dt,dx,dz,mm,g_u0,g_w0,g_txx0,g_tzz0,
                g_u1,g_w1,g_txx1,g_tzz1,vp,c,npd,
                g_tpx1,g_tpx0,g_tpz1,g_tpz0,g_tqx1,g_tqx0,g_tqz1,g_tqz0,
                acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,fs_sxd,ds_sxd,zs_sxd,is);
	  for(j=0;j<nz+2*npd;j++)
	  {
		for(i=0;i<nx+2*npd;i++)
		{
			g_u0[i][j]=g_u1[i][j];
			g_w0[i][j]=g_w1[i][j];

			g_tpx0[i][j]=g_tpx1[i][j];
			g_tpz0[i][j]=g_tpz1[i][j];
                  g_tqx0[i][j]=g_tqx1[i][j];
			g_tqz0[i][j]=g_tqz1[i][j];

			g_txx0[i][j]=g_txx1[i][j];
                  g_tzz0[i][j]=g_tzz1[i][j];
		}
	  }
/************************** back propagation end *************************/
/************************** snap *************************/
        if((is==1)&&(it%100==0)&&(flag_snap))
           {
              fseek(fpsnap,(int)(it/100)*(nx)*(nz)*4L,0);
              fseek(fpsnap1,(int)(it/100)*(nx)*(nz)*4L,0);
              fseek(fpsnap2,(int)(it/100)*(nx)*(nz)*4L,0);
              fseek(fpsnap3,(int)(it/100)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                       {
                    fwrite(&s_txx1[i][j],4L,1,fpsnap);
                    fwrite(&s_tzz1[i][j],4L,1,fpsnap1);
                    fwrite(&g_txx1[i][j],4L,1,fpsnap2);
                    fwrite(&g_tzz1[i][j],4L,1,fpsnap3);
                       }
           }

/************************** image condition start *************************/
/************************** image condition start *************************/
        for(i=npd;i<nx+npd;i++)
            for(j=npd;j<nz+npd;j++)
                {
                ia=0;
          /*############ poynting adcig ############*///poynting method is good.
              if(flag_adcig==1){
                 sx=-s_txx0[i][j]*s_u0[i][j];
                 sz=-s_tzz0[i][j]*s_w0[i][j];
                 gx= g_txx0[i][j]*g_u0[i][j];
                 gz= g_tzz0[i][j]*g_w0[i][j];
    /*   sx=-s_txx0[i][j]*(1.125*(s_txx0[i+1][j]-s_txx0[i][j])-0.0416666666667*(s_txx0[i+2][j]-s_txx0[i-1][j]));
       sz=-s_tzz0[i][j]*((s_tzz0[i][j+1]-s_tzz0[i][j])-0.0416666666667*(s_tzz0[i][j+2]-s_tzz0[i][j-1]));
       gx= g_txx0[i][j]*(1.125*(g_txx0[i+1][j]-g_txx0[i][j])-0.0416666666667*(g_txx0[i+2][j]-g_txx0[i-1][j]));
       gz= g_tzz0[i][j]*((g_tzz0[i][j+1]-g_tzz0[i][j])-0.0416666666667*(g_tzz0[i][j+2]-g_tzz0[i][j-1]));*/
                 b1=sx*sx+sz*sz;
                 b2=gx*gx+gz*gz;
                 a =sx*gx+sz*gz;
                 a1=a/(sqrtf(b1*b2)+0.0000000000001);
                 if(a1>=-1&&a1<=1)
                       {
                     a2=0.5*acosf(a1)*180/pi;
                     ia=(int)(a2/1.0);
                     if(ia==90)ia-=1;
                       }
          /*############ source-receivers adcig ############*///imageing is very bad, pass this method.
              }else if(flag_adcig==2){
                 sx=-s_txx0[i][j]*s_u0[i][j];
                 sz=-s_tzz0[i][j]*s_w0[i][j];
                 gx= g_txx0[i][j]*g_u0[i][j];
                 gz= g_tzz0[i][j]*g_w0[i][j];
                 b1=sqrtf(sx*sx+sz*sz);
                 b2=sqrtf(gx*gx+gz*gz);
                 sz=sx/b1;
                 gz=gx/b2;
                 ia=0;
                 if(sz>=-1&&sz<=1&&gz>=-1&&gz<=1)
                       {
                     a2=0.5*( asinf(sz) + asinf(gz) )*180/pi;
                     ia=(int)(a2/1.0);
                     if(ia<0)ia=-1*ia;
                     if(ia==90)ia-=1;
                       }
                       }
          /*############ half-offset domain adcig ############*///this method is very compelicated, pass this method.
              else if(flag_adcig==3){
                 if(flag_cdp==1){/* both sides */
                    
                for(ih=-nh/2;ih<nh/2+1;ih++)
                     {
                  if( (i+ih-npd)<nx && (i-ih-npd)>0 && (i+ih-npd)>0 && (i-ih-npd)<nx  )
                        {
                    Ixhz_is[i-npd][ih+nh/2][j-npd]+=(g_txx0[i+ih][j]*s_txx0[i-ih][j]);//+g_tzz0[i+ih][j]*s_tzz0[i-ih][j]);

                    //sourceI[i-npd][ih+nh][j-npd]+=s_txx0[i+ih][j]*s_txx0[i-ih][j];
                    
                        }
                  //if(sourceI[i-npd][ih+nh][j-npd]==0)sourceI[i-npd][ih+nh][j-npd]=1.0;
                     }

                 }else if(flag_cdp==0){
                     printf("Don't support this image condition!\n");
                     printf("When flag_cdp==0,flag_adcig!=3  .\n");
                      }
                  
                     }
          /*#######################end############################*/
                 
                adcig_is[ia][i-npd][j-npd]+=(g_txx0[i][j]*s_txx0[i][j])//+g_tzz0[i][j]*s_tzz0[i][j])
                                                  *pow(cos(ia*pi/180.0),3);
                      

                source_x[i-npd][j-npd]+=pow(s_txx0[i][j],2);
                source_z[i-npd][j-npd]+=pow(s_tzz0[i][j],2);
                if(source_x[i-npd][j-npd]==0)source_x[i-npd][j-npd]=1.0;
                if(source_z[i-npd][j-npd]==0)source_z[i-npd][j-npd]=1.0;

                mig_is[i-npd][j-npd]+=(g_txx0[i][j]*s_txx0[i][j]);//+g_tzz0[i][j]*s_tzz0[i][j]);
                  

                }
/************************** image condition over *************************/
/************************** image condition over *************************/
       //pickup_adcig();



     }//nt->0 loop end

/************************* is stack to ns **************************/
        for(i=npd;i<nx+npd;i++)
            for(j=npd;j<nz+npd;j++)
                {
               for(ia=0;ia<90;ia++)
                    {
                   adcig_is[ia][i-npd][j-npd]/=(source_x[i-npd][j-npd]);//+source_z[i-npd][j-npd]);
                   if(flag_cdp==1){
                      adcig_ns0[ia][i+fs_initial+(is-1)*ds_initial-nx/2-1-npd][j-npd]+=adcig_is[ia][i-npd][j-npd];
                   }else{
                      adcig_ns0[ia][i-npd][j-npd]+=adcig_is[ia][i-npd][j-npd];
                         }
                    }
               mig_is[i-npd][j-npd]/=(source_x[i-npd][j-npd]);//+source_z[i-npd][j-npd]);
               if(flag_cdp==1){
                  mig_ns0[i+fs_initial+(is-1)*ds_initial-nx/2-1-npd][j-npd]+=mig_is[i-npd][j-npd];

                  if(flag_adcig==3)
                    for(ih=0;ih<nh;ih++)
                      Ixhz_ns0[i+fs_initial+(is-1)*ds_initial-nx/2-1-npd][ih][j-npd]+=Ixhz_is[i-npd][ih][j-npd];
               }else{
                  mig_ns0[i-npd][j-npd]+=mig_is[i-npd][j-npd];
                   }
                }

/**********************close************************/
          if((is==1)&&(flag_snap))
            {
          fclose(fpsnap);   fclose(fpsnap1);
          fclose(fpsnap2);  fclose(fpsnap3);
            }
/**********************free*************************/
          free1float(coffx1);  free1float(coffx2);
          free1float(coffz1);  free1float(coffz2);
          free1float(acoffx1); free1float(acoffx2);
          free1float(acoffz1); free1float(acoffz2);
/********************** free source alloc ************************/
          free2float(s_u0);    free2float(s_u1);
          free2float(s_w0);    free2float(s_w1);
          free2float(s_txx0);  free2float(s_txx1);  free2float(s_tzz0);  free2float(s_tzz1);
          free2float(s_tpx0);  free2float(s_tpx1);  free2float(s_tpz0);  free2float(s_tpz1);
          free2float(s_tqx0);  free2float(s_tqx1);  free2float(s_tqz0);  free2float(s_tqz1);
/********************** free receiver alloc ************************/
          free2float(g_u0);    free2float(g_u1);
          free2float(g_w0);    free2float(g_w1);
          free2float(g_txx0);  free2float(g_txx1);  free2float(g_tzz0);  free2float(g_tzz1);
          free2float(g_tpx0);  free2float(g_tpx1);  free2float(g_tpz0);  free2float(g_tpz1);
          free2float(g_tqx0);  free2float(g_tqx1);  free2float(g_tqz0);  free2float(g_tqz1);
/********************** free -------- alloc ************************/
          free2float(source_x);free2float(source_z);free3float(sourceI);
          
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
       if(fs>nx){
          printf("###########################\n");
          printf("#Shot location(%f) >> nx(%d)\n",fs,nx);
          printf("###########################\n");
          exit(0);
         }
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
           float **u0,float **w0,float **txx0,float **tzz0,
           float **u1,float **w1,float **txx1,float **tzz1,
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
                     u1[i][j]=coffx2[i]*u0[i][j]-coffx1[i]*dtx*rho[i][j]*xx;
                     w1[i][j]=coffz2[j]*w0[i][j]-coffz1[j]*dtz*rho[i][j]*zz;
			 }
		 }
}
/*************************************func**********************************************/
void update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
            float **u0,float **w0,float **txx0,float **tzz0,
            float **u1,float **w1,float **txx1,float **tzz1,
            float **s,float **vp,float c[],int npd,
            float **tpx1,float **tpx0,float **tpz1,float **tpz0,
            float **tqx1,float **tqx0,float **tqz1,float **tqz0,
            float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
            float **deta,float **epsilu,int xsn,int ds_sxd,int zsn,int is,int _Circle_)
{
		 int i,j,ii,im,ix,iz;
		 float dtx,dtz,dux,dwz,xx,zz;
             int fs,ixs,izs,CR;
             float **deta1,**epsilu1;

            fs=xsn+(is-1)*ds_sxd;
            ixs=(int)(fs+0.5)+npd-1;
            izs=(int)(zsn+0.5)+npd-1;

            CR=_Circle_;///////////////////////

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
                    xx+=c[im]*(u1[i+im][j]-u1[i-im-1][j]);
                    zz+=c[im]*(w1[i][j+im]-w1[i][j-im-1]);
                        }
                 tpx1[i][j]=acoffx2[i]*tpx0[i][j]-acoffx1[i]*vp[i][j]*(1+2*epsilu1[i][j])*dtx*xx;
                 tpz1[i][j]=acoffz2[j]*tpz0[i][j]-acoffz1[j]*vp[i][j]*(pow((1+2*deta1[i][j]),0.5))*dtz*zz;
                 tqx1[i][j]=acoffx2[i]*tqx0[i][j]-acoffx1[i]*vp[i][j]*(pow((1+2*deta1[i][j]),0.5))*dtx*xx;
                 tqz1[i][j]=acoffz2[j]*tqz0[i][j]-acoffz1[j]*vp[i][j]*dtz*zz;
                 txx1[i][j]=tpx1[i][j]+tpz1[i][j]+s[i][j];
                 tzz1[i][j]=tqx1[i][j]+tqz1[i][j]+s[i][j];
		    }
		 }
             free2float(deta1);free2float(epsilu1);
}
/*************************************func**********************************************/
void inv_update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
            float **u0,float **w0,float **txx0,float **tzz0,
            float **u1,float **w1,float **txx1,float **tzz1,
            float **vp,float c[],int npd,
            float **tpx1,float **tpx0,float **tpz1,float **tpz0,
            float **tqx1,float **tqx0,float **tqz1,float **tqz0,
            float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
            float **deta,float **epsilu,int xsn,int ds_sxd,int zsn,int is)
{
		 int i,j,ii,im,ix,iz;
		 float dtx,dtz,dux,dwz,xx,zz;
             float **deta1,**epsilu1;

		 dtx=dt/dx;
		 dtz=dt/dz;
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
		   for(j=mm;j<=(2*npd+nz-mm-1);j++)
		    {
                  xx=0.0;
                  zz=0.0;
                  for(im=0;im<mm;im++)
                        {
                    xx+=c[im]*(u1[i+im][j]-u1[i-im-1][j]);
                    zz+=c[im]*(w1[i][j+im]-w1[i][j-im-1]);
                        }
                 tpx1[i][j]=acoffx2[i]*tpx0[i][j]-acoffx1[i]*vp[i][j]*(1+2*epsilu[i][j])*dtx*xx;
                 tpz1[i][j]=acoffz2[j]*tpz0[i][j]-acoffz1[j]*vp[i][j]*(pow((1+2*deta[i][j]),0.5))*dtz*zz;
                 tqx1[i][j]=acoffx2[i]*tqx0[i][j]-acoffx1[i]*vp[i][j]*(pow((1+2*deta[i][j]),0.5))*dtx*xx;
                 tqz1[i][j]=acoffz2[j]*tqz0[i][j]-acoffz1[j]*vp[i][j]*dtz*zz;
                 txx1[i][j]=tpx1[i][j]+tpz1[i][j];
                 tzz1[i][j]=tqx1[i][j]+tqz1[i][j];
		    }
		 }
}
/***************************************func********************************************/
float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,
                   float tmax,float favg,float dtout,float dt,float **vp,float ndtt,int myid)
{
		 int i,j;
		 float vpmax,vpmin,H_min;
		 float dt_max,dx_max,dz_max,d0;

		 vpmax=vp[npd][npd];
		 vpmin=vp[npd][npd];
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 if(vpmax<vp[i][j]) vpmax=vp[i][j];
				 if(vpmin>vp[i][j]) vpmin=vp[i][j];
			 }
		 }
		 d0=3.0*vpmax*log(100000.0)/(2.0*npd*dx);
		 if(dx<dz) H_min=dx;
		 else H_min=dz;
/********determine time sampling interval to ensure stability*****/
		 dt_max=0.5*1000*H_min/vpmax;
             dx_max=vpmin/favg*0.2;
             dz_max=dx_max;

            if((dx_max<dx)&&(myid==0))
                {
               printf("dx_max===%f, vpmin===%f, favg===%f \n",dx_max,vpmin,favg);
		   printf("YOU NEED HAVE TO REDEFINE DX ! \n");
               //exit(0);
		 }
             if((dz_max<dz)&&(myid==0))
		 {
		   printf("YOU NEED HAVE TO REDEFINE DZ ! \n");
               //exit(0);
		 }
	       if((dt_max<dt)&&(myid==0))
		 {
               printf("dt_max===%f, H_min===%f, vpmax===%f \n",dt_max,H_min,vpmax);
		   printf("YOU NEED HAVE TO REDEFINE dt ! \n");
               //exit(0);
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
void read_v_e_d_r(char FN1[],char FN2[],char FN3[],int nx,int nz,float **vv,float **epsilu,float **deta,
               float **rho0,int npd,int seismic)
{

	int i,j,k;

		FILE *fp1,*fp2,*fp3;
		if((fp1=fopen(FN1,"rb"))==NULL)printf("Open <%s> error!\n",FN1);
		if((fp2=fopen(FN2,"rb"))==NULL)printf("Open <%s> error!\n",FN2);
		if((fp3=fopen(FN3,"rb"))==NULL)printf("Open <%s> error!\n",FN3);
            rewind(fp1);
            rewind(fp2);
            rewind(fp3);
		 for(i=npd;i<nx+npd;i++)
		 {
                 if(seismic==2)  fseek(fp1,((i-npd+1)*HDRBYTES+(i-npd)*nz*FSIZE),0);
                 if(seismic==2)  fseek(fp2,((i-npd+1)*HDRBYTES+(i-npd)*nz*FSIZE),0);
                 if(seismic==2)  fseek(fp3,((i-npd+1)*HDRBYTES+(i-npd)*nz*FSIZE),0);
		     for(j=npd;j<nz+npd;j++)
		      {
		       fread(&vv[i][j],FSIZE,1,fp1);//vv[i][j] *= 0.75;
                   if(seismic==2) swap_float_4(&vv[i][j]);
		       fread(&epsilu[i][j],FSIZE,1,fp2);
                   if(seismic==2) swap_float_4(&epsilu[i][j]);
		       fread(&deta[i][j],FSIZE,1,fp3);//if(epsilu[i][j]<deta[i][j])epsilu[i][j]=deta[i][j];
                   if(seismic==2) swap_float_4(&deta[i][j]);
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
			 acoffx2[i]=acoffx1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz2[i]=acoffz1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);

		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 acoffx1[i]=1/(1+(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
			 acoffx2[i]=acoffx1[i]*(1-(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
		 }
		 for(i=npd+nz;i<nz+2*npd;i++)
		 {
			 acoffz1[i]=1/(1+(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
			 acoffz2[i]=acoffz1[i]*(1-(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
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
        for(j=0;j<mu_nt;j++)
           p_cal[i][j]=0.0;
       }else{}
}
/***************smooth2float****************/
void smooth2float(int nx,int rx,int nz,int rz,float **v)
{
  void smooth1float(float *v,int r,int n);

      float *f=NULL;
      int nmax;	/* max of nz and nx */
	int ix, iz;	/* counters */

	nmax = (nz<nx)?nx:nz;

	f = alloc1float(nmax);

        for(iz=0; iz<nz; ++iz)
           {
	 	for(ix=0; ix<nx; ++ix)
                {
			f[ix] = v[ix][iz];
		}
            smooth1float(f,rx,nx);
	 	for(ix=0; ix<nx; ++ix)
			v[ix][iz] = f[ix];
	    }
         for(ix=0; ix<nx; ++ix)
            {
	 	for(iz=0; iz<nz; ++iz)
                {
			f[iz] = v[ix][iz];
		}
              smooth1float(f,rz,nz);
	 	for(iz=0; iz<nz; ++iz)
			v[ix][iz] = f[iz];
	   }

}
/***************smooth1float****************/
void smooth1float(float *v,int r,int n)
{
  int i,j,ir;
  float *v0;
   v0=alloc1float(n);
   zero1float(v0,n);

  for(i=0;i<n;i++)
    v0[i]=v[i];
  for(ir=0;ir<r;ir++)
   {
    for(i=1;i<n-1;i++)
        v[i]=(v0[i-1]+v0[i]+v0[i+1])/3.0;
    v[0]=(v0[0]+v0[1])/2.0;
    v[n-1]=(v0[n-1]+v0[n-2])/2.0;
    for(i=0;i<n;i++)
       v0[i]=v[i];
   }
  for(i=0;i<n;i++)
     v[i]=v0[i];
  free1float(v0);

}
/*************** SU header conversion ****************/
/*************** SU header conversion ****************/
void swap_short_2(short *tni2)
/**************************************************************************
swap_short_2		swap a short integer
***************************************************************************/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));
}
void swap_u_short_2(unsigned short *tni2)
/**************************************************************************
swap_u_short_2		swap an unsigned short integer
***************************************************************************/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));
}
void swap_int_4(int *tni4)
/**************************************************************************
swap_int_4		swap a 4 byte integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
void swap_u_int_4(unsigned int *tni4)
/**************************************************************************
swap_u_int_4		swap an unsigned integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
void swap_long_4(long *tni4)
/**************************************************************************
swap_long_4		swap a long integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
void swap_u_long_4(unsigned long *tni4)
/**************************************************************************
swap_u_long_4		swap an unsigned long integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
void swap_float_4(float *tnf4)
/**************************************************************************
swap_float_4		swap a float
***************************************************************************/
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
void swap_double_8(double *tndd8)
/**************************************************************************
swap_double_8		swap a double
***************************************************************************/
{
  char *tnd8=(char *)tndd8;
  char tnc;

  tnc= *tnd8;
  *tnd8= *(tnd8+7);
  *(tnd8+7)=tnc;

  tnc= *(tnd8+1);
  *(tnd8+1)= *(tnd8+6);
  *(tnd8+6)=tnc;

  tnc= *(tnd8+2);
  *(tnd8+2)= *(tnd8+5);
  *(tnd8+5)=tnc;

  tnc= *(tnd8+3);
  *(tnd8+3)= *(tnd8+4);
  *(tnd8+4)=tnc;
}
/*************** SU header conversion ****************/
void swaphval(segy *tr, int index)
{
        register char *tp= (char *) tr;
        switch(*(hdr[index].type)) {
        case 'h': swap_short_2((short*)(tp + hdr[index].offs));
        break;
        case 'u': swap_u_short_2((unsigned short*)(tp + hdr[index].offs));
        break;
        case 'i': swap_int_4((int*)(tp + hdr[index].offs));
        break;
        case 'p': swap_u_int_4((unsigned int*)(tp + hdr[index].offs));
        break;
        case 'l': swap_long_4((long*)(tp + hdr[index].offs));
        break;
        case 'v': swap_u_long_4((unsigned long*)(tp + hdr[index].offs));
        break;
        case 'f': swap_float_4((float*)(tp + hdr[index].offs));
        break;
        case 'd': swap_double_8((double*)(tp + hdr[index].offs));
        break;
        default: err("%s: %s: unsupported data type", __FILE__, __LINE__);
        break;
        }
}

