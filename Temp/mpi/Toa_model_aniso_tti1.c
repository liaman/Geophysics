//#########################################################
//##
//##         2D Acoustic TTI Medium Forward 
//##              
//##  Ps : P + sv wave        
//##                                     Rong Tao 
//##                 
//##
//##  Ps:  du/dt=1/rho*dp/dx , dw/dt=1/rho*dq/dz
//##      
//##       dp/dt=rho*vpx^2*du/dx+rho*vp0*vpn*dw/dz
//##
//##       dq/dt=rho*vp0*vpn*du/dx+rho*vp0^2*dw/dz
//##
//##                    vpx^2=vp0^2*(1+2*epsilu);
//##                    vpn^2=vp0^2*(1+2*deta);
//##
//##                    d/dx'=cos@*d/dx+sin@*d/dz
//##                    d/dz'=-sin@*d/dx+cos@*d/dz
//##
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
main(int argc,char *argv[])
{


        
        void model(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
                   float vdx,float vdz,float favg,float tmax,float dt,float dtout,
                   float pfac,char FN1[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
                   int is,float **p_cal,float deta,float epsilu,float ceta);
	void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);        
        void freea(float **a,int m);
        void zero(float **a,int nx,int nz);
              
//	char FN1[250]={"vmar_350X230w.dat"};
//	char FN3[250]={"shot_mar_350X4001.dat"};
        char FN1[250]={"vel_300_200.dat"};
	char FN3[250]={"shot_vel.dat"};

        FILE *fp3;
        fp3=fopen(FN3,"wb");
        
	int i,j,k,is,nx,nz,nt,vnx,vnz,i_start,i_end;
	int ns_sxd,ds_sxd,fs_sxd,zs_sxd,npd;
	
	float dx,dz,vdx,vdz,tmax,dt,dtout,pfac,favg,deta,epsilu,ceta;
	
	int myid,numprocs,count;

/**********************************************************/
        nx=300;         npd=50;      tmax=0.5;
	nz=200;         favg=50;     pfac=1000.0;
 	                dx=5.0;   
	vnx=300;        dz=5.0;   
	vnz=200;        vdx=5.0;
	nt=1001;         vdz=5.0;
        ns_sxd=1;       dt=0.5;
        ds_sxd=10;      dtout=0.5;
        fs_sxd=150;
        zs_sxd=100;     

        deta=0.10;
        epsilu=0.24;
        ceta=90;  // degree

        ceta=ceta*3.141593/180.0;
          
 
        float **p_cal;

	p_cal=alloc2float(nt,nx);

/*******************************************/
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
/******************************************/	
      MPI_Barrier(MPI_COMM_WORLD);
/******************************************************************************/
    if(myid==0)
    {
        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("-----sin(ceta)===%f-----\n",sin(ceta)); 
        printf("\n");
        printf("\n");
        printf("\n");   
    }
/*************************************/
    MPI_Barrier(MPI_COMM_WORLD);                                                       
/**************begin the loop******************************/
    for(is=1+myid;is<=ns_sxd;is=is+numprocs)	
    {     
        if(myid==0)
        {
             printf("begin the loop************ IS========%d  \n",is);
        }
 
/***************************************/
        zero2float(p_cal,nt,nx);
/***************************************/
   
        if(myid==0)
	{
	  printf("the model is start    \n");
	}
	
        model(nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,FN1,
               ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_cal,deta,epsilu,ceta);        
            
        fseek(fp3,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	{ 
	  for(j=0;j<nt;j++)
	  {
	    fwrite(&p_cal[i][j],4L,1,fp3);
	  }
	}
	
		     
	if(myid==0)
	{
	  printf("the model is over    \n");
	}
 	
        MPI_Barrier(MPI_COMM_WORLD);

    } 
/************** the loop is done***************************/ 
      if(myid==0)
      {
          printf("the loop is done\n");
      }

      MPI_Barrier(MPI_COMM_WORLD);
  
      fclose(fp3);
/******************************************/ 
      free2float(p_cal);
      if(myid==0)
      {
          printf("complete!!!!!!!!! \n");
          printf("**********************************************\n");
      }            
      MPI_Finalize();
}
/************************************************************************/
void zero(float **a,int nx,int nz)
{
	int i,j;
  	for(i=0;i<nx;i++)
	{
		for(j=0;j<nz;j++)
		{
			a[i][j]=0.0;
		}
	}
}
/*******************************************************************************/

void model(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,
           float vdx,float vdz,float favg,float tmax,float dt,float dtout,
           float pfac,char FN1[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
           int is,float **p_cal,float deta,float epsilu,float ceta)

{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
                  float favg,float **s,int wtype,float pi,int npd,int is,int ds_sxd);
        void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                  float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
                  float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
                  float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,
                  float *coffz2,float ceta);
        void update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
                  float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
                  float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
                  float **s,float **vp,float c[],int npd,
                  float **tpx1,float **tpx0,float **tpz1,float **tpz0,
                  float **tqx1,float **tqx0,float **tqz1,float **tqz0,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                  float deta,float epsilu,float ceta);
        void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
        float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,
                  int npd,float tmax,float favg,float dtout,float dt,float **vp0,float ndtt);
        void pad_vv(int nx,int nz,int npd,float **ee);
	void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);
	void current_shot(float **vp0,float **rho0,float **vp,float **rho,
                  int nx,int nz,int npd,int vnx,int vnz,int ds_sxd,int is);
	void initial_coffe(float dt,float d0,int nx,int nz,
                  float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);	

      


	  int mm=4; 
          int hsx=1;
	  int i,j;

	  int ntout,wtype,it,ifx,ilx,jfz,jlz;
	  float pi,t,ndtt,d0;

	  
	  FILE *fp1;

          wtype=1;
         
          
	  pi=3.141593;
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

	  float **s;
     
	  float c[4];
	  

          cal_c(mm,c);   

      
	  vp0=alloc2float(nz+2*npd,nx+2*npd);                                               
	  rho0=alloc2float(nz+2*npd,nx+2*npd);                                            

          zero2float(vp0,nz+2*npd,nx+2*npd); 
          zero2float(rho0,nz+2*npd,nx+2*npd); 
                  
/********************************************************/
          read_file(FN1,vnx,vnz,vp0,rho0,npd);          
          
      
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
        
	  for(i=0;i<nx+2*npd;i++)
	  {
		 coffx1[i]=0.0;
		 coffx2[i]=0.0;
		 acoffx1[i]=0.0;
		 acoffx2[i]=0.0;
	  }

	  for(i=0;i<nz+2*npd;i++)
	  {
		 coffz1[i]=0.0;
		 coffz2[i]=0.0;
		 acoffz1[i]=0.0;
		 acoffz2[i]=0.0;
	  }

          initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,
                        acoffx1,acoffx2,acoffz1,acoffz2,npd);
       
/***********************************************************/
	   ndtt=(int)ndtt;



       //  printf("IS============%d \n",is);
/***********************************************/  
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
              fpsnap=fopen("snap.dat","wb");
              fpsnap1=fopen("snap1.dat","wb");

    for(it=0,t=0.0;it<nt;it++,t+=dt)
    { 

         if((int)(it%100)==0)
          printf("is===%d   it===%d\n",is,it);

	  ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd);
          update_vel(nx,nz,npd,mm,dt,dx,dz,up0,wp0,uq0,wq0,txx0,tzz0,
                     up1,wp1,uq1,wq1,txx1,tzz1,rho,c,coffx1,coffx2,coffz1,coffz2,ceta);
          update_stress(nx,nz,dt,dx,dz,mm,up0,wp0,uq0,wq0,txx0,tzz0,
                     up1,wp1,uq1,wq1,txx1,tzz1,s,vp,c,npd,tpx1,tpx0,tpz1,tpz0,
                     tqx1,tqx0,tqz1,tqz0,acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,ceta);
       


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

           if((is==1)&&(it%5==0))
           {
              
    /*          fseek(fpsnap,(int)(it/10)*nx*nz*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&txx1[i][j],4L,1,fpsnap);
           
              
              fseek(fpsnap1,(int)(it/10)*nx*nz*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&tzz1[i][j],4L,1,fpsnap1);     */

              fseek(fpsnap,(int)(it/5)*(nx)*(nz)*4L,0);
              for(i=npd;i<nx+npd;i++)
                 for(j=npd;j<nz+npd;j++)
                    fwrite(&txx1[i][j],4L,1,fpsnap);
           
              
              fseek(fpsnap1,(int)(it/5)*(nx+2*npd)*(nz+2*npd)*4L,0);
              for(i=0;i<nx+2*npd;i++)
                 for(j=0;j<nz+2*npd;j++)
                    fwrite(&tzz1[i][j],4L,1,fpsnap1);
              
           }


     }
          fclose(fpsnap);fclose(fpsnap1);
/***********************************************/        
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

/*******************************************************/
                                                       //
                                                       //
void cal_c(int mm,float c[])                           //                   
{                                                      //
                                                       //
                                                       //
                                                       //
                                                       //
	  c[0]=1.196289;                               //
          c[1]=-0.0797526;                             //
          c[2]=0.009570313;                            //
          c[3]=-0.000697544;                           //
	                                                 //
                                                       //
}                                                      //




void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,
              float favg,float **s,int wtype,float pi,int npd,int is,int ds_sxd)
{

	    float get_source(float ts,float favg,int wtype);

	    int i,j,ixs,izs,x,z;
	    float tdelay,ts,source,fs;
      
	    for(j=0;j<=nz+2*npd-1;j++)
	    {
		   for(i=0;i<=nx+2*npd-1;i++)
		   {
			s[i][j]=0.0;
		   }
	    }
          
	    tdelay=1.0/favg;
            ts=t-tdelay;
            fs=xsn+(is-1)*ds_sxd;
       if(t<=2*tdelay)
       {
            source=get_source(ts,favg,wtype);            
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
/*float get_source(float ts,float favg,int wtype)
{
        int i,j,nt,ntfft;
        int nf1,nf2,nf3,nf4;
	    float x,pi;
	    float source=0.0;
	    float  f1,f2,dw,fw;
        float tmpp=0.0;
        float dt=0.5/1000;
        float tdelay;   
        float    *trace;
        complex  *ctrace;
        tdelay=1.0/favg;
        nt=4001;

        f1=5.0;
        f2=6.0; 
        pi=3.14;     
        dw=2.0*pi/(ntfft*dt);
        fw=2.0*pi*f1;
        nf1=fw/dw+0.5;
        fw=2.0*pi*f2;
        nf2=fw/dw+0.5;
        
        ntfft=npfa(nt);
        nf3=ntfft-1;
        nf4=ntfft;
        trace=alloc1float(nt);
        ctrace=alloc1complex(ntfft);
        zero1float(trace,nt); 
        zero1complex(ctrace,ntfft);   
            
          for(i=0;i<nt;i++)
          {
          if(((i+1)*dt)<=2*tdelay)
          {
           x=(favg*((i+1)*dt-tdelay))*(favg*((i+1)*dt-tdelay));            
           trace[i]=(1-2*pi*pi*(x))*exp(-(pi*pi*x));
           }
           else
           trace[i]=0.0;
          } 
          
         for(i=0;i<nt;i++)
         {   
           ctrace[i]=cmplx(trace[i],0.0);
         }
         for(i=nt;i<ntfft;i++)
         { 
           ctrace[i]=cmplx(0.0, 0.0);	
         }	    
         pfacc(1, ntfft, ctrace);                    
               
	    for(i=0;i<ntfft;i++)
	    {
          if(i>=nf1&&i<=nf2)
          {
              tmpp=0.54+0.46*cos(PI*(i-nf1)/(nf2-nf1)-PI);
              ctrace[i]=crmul(ctrace[i],tmpp);
          }
          else if(i>=nf3&&i<=nf4)
          {
              tmpp=0.54+0.46*cos(PI*(nf3-i)/(nf4-nf3));
              ctrace[i]=crmul(ctrace[i],tmpp);
          }
          else if(i<nf1)
              ctrace[i]=cmplx(0.0, 0.0);
	    }  
	    
	    pfacc (-1, ntfft, ctrace); 
	    for(i=0;i<ntfft;i++)
	    {
	        ctrace[i]=crmul(ctrace[i],1.0/ntfft);
	    }      
		    j=((ts+tdelay)/dt)-1;
            source=ctrace[j].r;
       
	    return (source);


}*/
/*
float get_source(float ts,float favg,int wtype)
 {
		  float x,xx,pi;
		  float source=0.0,pi2;
		  pi=4*atan(1.0);
		  pi2=pi*pi;

		  if(wtype==1)//ricker wavelet
		  {
		      x=favg*ts;
		      xx=x*x;
		      source=(1-2*pi2*xx)*pow(exp,-pi2*xx);
		  }
		  if(wtype==2)//derivative of gaussian
		  {
		      x=(-4)*favg*favg*pi2/log(0.1);
		      source=(-2)*pi2*(ts)*pow(exp,-x*ts*ts);
          }
          if(wtype==3)//derivative of gaussian
          {
              x=(-1)*favg*favg*pi2/log(0.1);
              source=pow(exp,-x*ts*ts);
          }

		  return (source);

}
*/
float get_source(float ts,float favg,int wtype)
{
	    float x,pi;
	    float source=0.0;
       
	    pi=3.14;
            x=(favg*(ts))*(favg*(ts));
		   
            source=(1-2*pi*pi*(x))*exp(-(pi*pi*x));
	    return (source);


}



void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
                float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
                float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,
                float *coffz2,float ceta)
{
		 int ii,i,j;
		 float dtxx,dtxz,dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;
         
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
			 for(j=mm;j<=(2*npd+nz-mm-1);j++)
			 {
                   //  du/dt=dp/dx-dp/dz
                           up1[i][j]=coffx2[i]*coffz2[j]*up0[i][j] - coffx1[i]*coffz1[j]*rho[i][j]*( cos(ceta)*dtx*(c[0]*(txx0[i+1][j]-txx0[i][j])+c[1]*(txx0[i+2][j]-txx0[i-1][j])+c[2]*(txx0[i+3][j]-txx0[i-2][j])+c[3]*(txx0[i+4][j]-txx0[i-3][j])) - sin(ceta)*dtz*(c[0]*(txx0[i][j+1]-txx0[i][j])+c[1]*(txx0[i][j+2]-txx0[i][j-1])+c[2]*(txx0[i][j+3]-txx0[i][j-2])+c[3]*(txx0[i][j+4]-txx0[i][j-3])) );

                           wq1[i][j]=coffz2[j]*coffx2[i]*wq0[i][j] - coffz1[j]*coffx1[i]*rho[i][j]*( cos(ceta)*dtz*(c[0]*(tzz0[i][j+1]-tzz0[i][j])+c[1]*(tzz0[i][j+2]-tzz0[i][j-1])+c[2]*(tzz0[i][j+3]-tzz0[i][j-2])+c[3]*(tzz0[i][j+4]-tzz0[i][j-3])) + sin(ceta)*dtx*(c[0]*(tzz0[i+1][j]-tzz0[i][j])+c[1]*(tzz0[i+2][j]-tzz0[i-1][j])+c[2]*(tzz0[i+3][j]-tzz0[i-2][j])+c[3]*(tzz0[i+4][j]-tzz0[i-3][j])) );    

                   //  du/dt=dp/dx-dq/dz   ****** so wrong 
                      /*     up1[i][j]=coffx2[i]*up0[i][j] - coffx1[i]*rho[i][j]*cos(ceta)*dtx*(c[0]*(txx0[i+1][j]-txx0[i][j])+c[1]*(txx0[i+2][j]-txx0[i-1][j])+c[2]*(txx0[i+3][j]-txx0[i-2][j])+c[3]*(txx0[i+4][j]-txx0[i-3][j])) + coffz1[j]*rho[i][j]*sin(ceta)*dtz*(c[0]*(tzz0[i][j+1]-tzz0[i][j])+c[1]*(tzz0[i][j+2]-tzz0[i][j-1])+c[2]*(tzz0[i][j+3]-tzz0[i][j-2])+c[3]*(tzz0[i][j+4]-tzz0[i][j-3]));

                           wq1[i][j]=coffz2[j]*wq0[i][j] - coffz1[j]*rho[i][j]*cos(ceta)*dtz*(c[0]*(tzz0[i][j+1]-tzz0[i][j])+c[1]*(tzz0[i][j+2]-tzz0[i][j-1])+c[2]*(tzz0[i][j+3]-tzz0[i][j-2])+c[3]*(tzz0[i][j+4]-tzz0[i][j-3])) - coffx1[i]*rho[i][j]*sin(ceta)*dtx*(c[0]*(txx0[i+1][j]-txx0[i][j])+c[1]*(txx0[i+2][j]-txx0[i-1][j])+c[2]*(txx0[i+3][j]-txx0[i-2][j])+c[3]*(txx0[i+4][j]-txx0[i-3][j]));   */

			 }
		 }


         //合并到上面了，这部分没有删除的东西
	/*	 for(j=mm;j<=(2*npd+nz-mm-1);j++)
		 { 
			 for(i=mm;i<=(2*npd+nx-mm-1);i++)
			 {
			   up1[i][j]=coffx2[i]*up0[i][j]-coffx1[i]*rho[i][j]*( cos(ceta)*dtx*(c[0]*(txx0[i+1][j]-txx0[i][j])+c[1]*(txx0[i+2][j]-txx0[i-1][j])+c[2]*(txx0[i+3][j]-txx0[i-2][j])+c[3]*(txx0[i+4][j]-txx0[i-3][j]))) + coffz1[j]*rho[i][j]*(sin(ceta)*dtz*(c[0]*(txx0[i][j+1]-txx0[i][j])+c[1]*(txx0[i][j+2]-txx0[i][j-1])+c[2]*(txx0[i][j+3]-txx0[i][j-2])+c[3]*(txx0[i][j+4]-txx0[i][j-3])));
			 }
		 }       */
		

}
     

void update_stress(int nx,int nz,float dt,float dx,float dz,int mm,
              float **up0,float **wp0,float **uq0,float **wq0,float **txx0,float **tzz0,
              float **up1,float **wp1,float **uq1,float **wq1,float **txx1,float **tzz1,
              float **s,float **vp,float c[],int npd,
              float **tpx1,float **tpx0,float **tpz1,float **tpz0,
              float **tqx1,float **tqx0,float **tqz1,float **tqz0,
              float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
              float deta,float epsilu,float ceta)
{

		 float dux,dwz;
		 int i,j,ii;
		 float dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;
                     
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
			 for(j=mm;j<=(2*npd+nz-mm-1);j++)
			 {
                     //  dp/dt=du/dx-du/dz
                           tpx1[i][j]=acoffx2[i]*acoffz2[j]*tpx0[i][j]-acoffx1[i]*acoffz1[j]*vp[i][j]*(1+2*epsilu)*( cos(ceta)*dtx*(c[0]*(up1[i][j]-up1[i-1][j])+c[1]*(up1[i+1][j]-up1[i-2][j])+c[2]*(up1[i+2][j]-up1[i-3][j])+c[3]*(up1[i+3][j]-up1[i-4][j])) - sin(ceta)*dtz*(c[0]*(up1[i][j]-up1[i][j-1])+c[1]*(up1[i][j+1]-up1[i][j-2])+c[2]*(up1[i][j+2]-up1[i][j-3])+c[3]*(up1[i][j+3]-up1[i][j-4])) );	
 
                           tpz1[i][j]=acoffz2[j]*acoffx2[i]*tpz0[i][j]-acoffz1[j]*acoffx1[i]*vp[i][j]*(pow((1+2*deta),0.5))*( cos(ceta)*dtz*(c[0]*(wq1[i][j]-wq1[i][j-1])+c[1]*(wq1[i][j+1]-wq1[i][j-2])+c[2]*(wq1[i][j+2]-wq1[i][j-3])+c[3]*(wq1[i][j+3]-wq1[i][j-4])) + sin(ceta)*dtx*(c[0]*(wq1[i][j]-wq1[i-1][j])+c[1]*(wq1[i+1][j]-wq1[i-2][j])+c[2]*(wq1[i+2][j]-wq1[i-3][j])+c[3]*(wq1[i+3][j]-wq1[i-4][j])) );	         

                           tqx1[i][j]=acoffx2[i]*acoffz2[j]*tqx0[i][j]-acoffx1[i]*acoffz1[j]*vp[i][j]*(pow((1+2*deta),0.5))*( cos(ceta)*dtx*(c[0]*(up1[i][j]-up1[i-1][j])+c[1]*(up1[i+1][j]-up1[i-2][j])+c[2]*(up1[i+2][j]-up1[i-3][j])+c[3]*(up1[i+3][j]-up1[i-4][j])) - sin(ceta)*dtz*(c[0]*(up1[i][j]-up1[i][j-1])+c[1]*(up1[i][j+1]-up1[i][j-2])+c[2]*(up1[i][j+2]-up1[i][j-3])+c[3]*(up1[i][j+3]-up1[i][j-4])) );      

                           tqz1[i][j]=acoffx2[i]*acoffz2[j]*tqz0[i][j]-acoffx1[i]*acoffz1[j]*vp[i][j]*( cos(ceta)*dtz*(c[0]*(wq1[i][j]-wq1[i][j-1])+c[1]*(wq1[i][j+1]-wq1[i][j-2])+c[2]*(wq1[i][j+2]-wq1[i][j-3])+c[3]*(wq1[i][j+3]-wq1[i][j-4])) + sin(ceta)*dtx*(c[0]*(wq1[i][j]-wq1[i-1][j])+c[1]*(wq1[i+1][j]-wq1[i-2][j])+c[2]*(wq1[i+2][j]-wq1[i-3][j])+c[3]*(wq1[i+3][j]-wq1[i-4][j])) );    

                     //  dp/dt=du/dx-dw/dz ********* so wrong 
                       /*    tpx1[i][j]=acoffx2[i]*tpx0[i][j]-acoffx1[i]*vp[i][j]*(1+2*epsilu)*cos(ceta)*dtx*(c[0]*(up1[i][j]-up1[i-1][j])+c[1]*(up1[i+1][j]-up1[i-2][j])+c[2]*(up1[i+2][j]-up1[i-3][j])+c[3]*(up1[i+3][j]-up1[i-4][j])) + acoffz1[j]*vp[i][j]*(1+2*epsilu)*sin(ceta)*dtz*(c[0]*(wq1[i][j]-wq1[i][j-1])+c[1]*(wq1[i][j+1]-wq1[i][j-2])+c[2]*(wq1[i][j+2]-wq1[i][j-3])+c[3]*(wq1[i][j+3]-wq1[i][j-4]));	
 
                           tpz1[i][j]=acoffz2[j]*tpz0[i][j]-acoffz1[j]*vp[i][j]*(pow((1+2*deta),0.5))*cos(ceta)*dtz*(c[0]*(wq1[i][j]-wq1[i][j-1])+c[1]*(wq1[i][j+1]-wq1[i][j-2])+c[2]*(wq1[i][j+2]-wq1[i][j-3])+c[3]*(wq1[i][j+3]-wq1[i][j-4])) -acoffx1[i]*vp[i][j]*sin(ceta)*dtx*(c[0]*(up1[i][j]-up1[i-1][j])+c[1]*(up1[i+1][j]-up1[i-2][j])+c[2]*(up1[i+2][j]-up1[i-3][j])+c[3]*(up1[i+3][j]-up1[i-4][j]));	         

                           tqx1[i][j]=acoffx2[i]*tqx0[i][j]-acoffx1[i]*vp[i][j]*(pow((1+2*deta),0.5))*cos(ceta)*dtx*(c[0]*(up1[i][j]-up1[i-1][j])+c[1]*(up1[i+1][j]-up1[i-2][j])+c[2]*(up1[i+2][j]-up1[i-3][j])+c[3]*(up1[i+3][j]-up1[i-4][j])) + acoffz1[j]*vp[i][j]*(pow((1+2*deta),0.5))*sin(ceta)*dtz*(c[0]*(wq1[i][j]-wq1[i][j-1])+c[1]*(wq1[i][j+1]-wq1[i][j-2])+c[2]*(wq1[i][j+2]-wq1[i][j-3])+c[3]*(wq1[i][j+3]-wq1[i][j-4]));      

                           tqz1[i][j]=acoffz2[j]*tqz0[i][j]-acoffz1[j]*vp[i][j]*cos(ceta)*dtz*(c[0]*(wq1[i][j]-wq1[i][j-1])+c[1]*(wq1[i][j+1]-wq1[i][j-2])+c[2]*(wq1[i][j+2]-wq1[i][j-3])+c[3]*(wq1[i][j+3]-wq1[i][j-4])) - acoffx1[i]*vp[i][j]*sin(ceta)*dtx*(c[0]*(up1[i][j]-up1[i-1][j])+c[1]*(up1[i+1][j]-up1[i-2][j])+c[2]*(up1[i+2][j]-up1[i-3][j])+c[3]*(up1[i+3][j]-up1[i-4][j]));	   */
     
                           txx1[i][j]=tpx1[i][j]+tpz1[i][j]+s[i][j];
                           tzz1[i][j]=tqx1[i][j]+tqz1[i][j]+s[i][j];

			 }
		 }
      

}           


           

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


/*====== determine time sampling interval to ensure stability====*/

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
		 {printf("dt_max===%f, H_min===%f, vpmax===%f \n",dt_max,H_min,vpmax);
		   printf("YOU NEED HAVE TO REDEFINE dt ! \n");
                   exit(0);
		 }
         
                 return d0;



}


void zero2d(int nx,int nz,float **vv0,int npd)
{
		 int i,j;
		 
		 for(i=0;i<nx+2*npd;i++)
		 {
			 for(j=0;j<nz+2*npd;j++)
			 {
				 vv0[i][j]=0.0;
			 }
		 }

}


void zero2d2(int nx,int nz,int npd,float **vv)
{
		 int i,j;


		 for(i=0;i<=nz+2*npd-1;i++)
		 {
			 for(j=0;j<=nx+2*npd-1;j++)
			 {
				 vv[i][j]=0.0;
			 }
		 }

}
    


void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,
              float **u1,float **w1,float **txx1,float **vp,float **rho)
{
		 int i,j,k;
 

		   for(i=0;i<nz+2*npd;i++)
		   {
			  for(j=0;j<nx+2*npd;j++)
			  {
				 vp[i][j]=0.0;
				 rho[i][j]=0.0;
				
				 u0[i][j]=0.0;
                                 w0[i][j]=0.0;
                                 txx0[i][j]=0.0;
				 u1[i][j]=0.0;
                                 w1[i][j]=0.0;
                                 txx1[i][j]=0.0;
			  }
		   }

}

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

void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd)
{

		 int i,j;
		
		 FILE *fp1;
		 fp1=fopen(FN1,"rb");
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 fread(&vv[i][j],4,1,fp1);

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
}


void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,
                  int npd,int vnx,int vnz,int ds_sxd,int is)
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



void freea(float **a,int m)
{
	int i;

	for(i=0;i<m;i++)
	{
		free(a[i]);
	}
	free(a);
}


