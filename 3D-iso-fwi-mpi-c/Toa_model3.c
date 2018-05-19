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

        void model_step(int nx,int ny,int nz,int vnx,int vny,int vnz,int nt,int npd,float dx,float dy,float dz,float vdx,float vdy,float vdz,float favg,float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int xs_sxd,int ys_sxd,int zs_sxd,int is,float ***p_cal);
	void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd);        

       
              
//	char FN1[250]={"vmar_350X230w.dat"};
//	char FN3[250]={"shot_mar_350X4001.dat"};
        char FN1[250]={"vel_100_100_100.dat"};
	char FN3[250]={"shot_100_100_601_obs.dat"};

        FILE *fp3;
        fp3=fopen(FN3,"wb");
        
	int i,j,k,is,nx,ny,nz,nt,vnx,vny,vnz,i_start,i_end;
	int ns_sxd,zs_sxd,npd;
	int xs_sxd[10],ys_sxd[10];
	float dx,dy,dz,vdx,vdy,vdz,tmax,dt,dtout,pfac,favg;
	
	int myid,numprocs,count;


        nx=100;          npd=5;   tmax=0.3;
        ny=100;          favg=30;  pfac=1000.0;
	nz=100;         
 	
        vnx=100;         dx=5.0;   
	vny=100;         dy=5.0;   
	vnz=100;         dz=5.0; 
                        vdx=5.0;
                        vdy=5.0;
	nt=601;         vdz=5.0;
        dt=0.5;
        dtout=0.5;
         
        ns_sxd=1;
  /*      xs_sxd[0]=10;xs_sxd[1]=20;xs_sxd[2]=30;xs_sxd[3]=40;
        xs_sxd[4]=50;xs_sxd[5]=60;xs_sxd[6]=70;
        ys_sxd[0]=10;ys_sxd[1]=20;ys_sxd[2]=30;ys_sxd[3]=40;
        ys_sxd[4]=50;ys_sxd[5]=60;ys_sxd[6]=70;   */
        xs_sxd[0]=50;
        ys_sxd[0]=50;
        zs_sxd=50;
             
        
        float ***p_cal;
      
	p_cal=alloc3float(nt,ny,nx);

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
        zero3float(p_cal,nt,ny,nx);
/***************************************/
   
        if(myid==0)
	{
	  printf("the model is start    \n");
	}
	
        model_step(nx,ny,nz,vnx,vny,vnz,nt,npd,dx,dy,dz,vdx,vdy,vdz,favg,tmax,dt,dtout,pfac,FN1,ns_sxd,xs_sxd[is-1],ys_sxd[is-1],zs_sxd,is,p_cal);        
            
        fseek(fp3,(is-1)*nx*ny*nt*4L,0);
        for(j=0;j<ny;j++)
	  for(i=0;i<nx;i++)
	  { 
	    for(k=0;k<nt;k++)
	    {
	      fwrite(&p_cal[i][j][k],4L,1,fp3);
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
      free3float(p_cal);
      if(myid==0)
      {
          printf("complete!!!!!!!!! \n");
          printf("**********************************************\n");
      }            
      MPI_Finalize();
}

/*******************************************************************************/


void model_step(int nx,int ny,int nz,int vnx,int vny,int vnz,int nt,int npd,float dx,float dy,float dz,float vdx,float vdy,float vdz,float favg,float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int xs_sxd,int ys_sxd,int zs_sxd,int is,float ***p_cal)
{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float ysn,float zsn,int nx,int ny,int nz,float dt,float t,float favg,float ***s,int wtype,float pi,int npd,int is);
        void update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,float ***u0,float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***rho,float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,float *coffz2);
        void update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***s,float ***vp,float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,float ***tz1,float ***tz0,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2);
        void abs_bc(float ***u1,float ***w1,float ***txx1,int nx,int ny,int nz,int npd,float absbord[]);
        float get_constant(float dx,float dy,float dz,int nx,int ny,int nz,int nt,int ntout,int npd,float tmax,float favg,float dtout,float dt,float ***vp0,float ndtt);
        void pad_vv(int nx,int ny,int nz,int npd,float ***ee);
	void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd);
	void current_shot(float ***vp0,float ***rho0,float ***vp,float ***rho,int nx,int ny,int nz,int npd,int vnx,int vny,int vnz,int is);	
        void initial_coffe(float dt,float d0,int nx,int ny,int nz,float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,float *coffz2,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2,int npd);
       // void initial_coffe2(float dt,float d0,int nx,int ny,int nz,float **coffxy1,float **coffxy2,float **coffxz1,float **coffxz2,float **coffyz1,float **coffyz2,float **acoffxy1,float **acoffxy2,float **acoffxz1,float **acoffxz2,float **acoffyz1,float **acoffyz2,int npd);
         void freea(float **a,int m);
      

          
	  int mm=4; 
          int hsx=1;
	  int i,j,k;

	  int ntout,wtype,it,ifx,ilx,ify,ily,jfz,jlz;
	  float pi,t,ndtt,d0;
	  
	  
	  FILE *fp1;

          wtype=1;
          

	  pi=3.141593;
	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;
     
	  float ***vp0;
	  float ***rho0;
	  float ***vp;
	  float ***rho;

	  float ***u0;
          float ***v0;
	  float ***w0;

          float ***u1;
          float ***v1;
	  float ***w1;

	  float ***txx0; 
          float ***txx1;
 
          float ***txx_snap;
	  
	  
	  float ***tx0;
	  float ***tx1;
          float ***ty0;
	  float ***ty1;
	  float ***tz0;
	  float ***tz1;
	  float ***s;
     
	  float c[4];
	  

          cal_c(mm,c);   

          vp0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          rho0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	  zero3float(vp0,vnz+2*npd,vny+2*npd,vnx+2*npd);  
          zero3float(rho0,vnz+2*npd,vny+2*npd,vnx+2*npd);                                              
     

    
                   
          read_file(FN1,vnx,vny,vnz,vp0,rho0,npd);          
      
          vp=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	  rho=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);  
                                            
          u0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          u1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          v0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          v1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          w0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          w1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 

          txx0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          txx1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 
         
          txx_snap=alloc3float(vnz,vny,vnx);

          tx0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          tx1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          ty0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          ty1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	  tz0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 
          tz1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 
          s=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	 
	 
	  


          d0=get_constant(dx,dy,dz,nx,ny,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);
          dt=dt/1000;
/**********************************************************************************************************************/
/***************************************************  wrong coffe  *******************************************************************/
/**********************************************************************************************************************/


   /*       float **coffxy1;
          float **coffxy2;
          float **coffxz1;
          float **coffxz2;
          float **coffyz1;
          float **coffyz2;
          float **acoffxy1;
          float **acoffxy2;
          float **acoffxz1;
          float **acoffxz2;
          float **acoffyz1;
          float **acoffyz2;

          coffxy1=alloc2float(vny+2*npd,vnx+2*npd);
          coffxy2=alloc2float(vny+2*npd,vnx+2*npd);
          coffxz1=alloc2float(vnz+2*npd,vnx+2*npd);
          coffxz2=alloc2float(vnz+2*npd,vnx+2*npd);
          coffyz1=alloc2float(vnz+2*npd,vny+2*npd);
          coffyz2=alloc2float(vnz+2*npd,vny+2*npd);

          acoffxy1=alloc2float(vny+2*npd,vnx+2*npd);
          acoffxy2=alloc2float(vny+2*npd,vnx+2*npd);
          acoffxz1=alloc2float(vnz+2*npd,vnx+2*npd);
          acoffxz2=alloc2float(vnz+2*npd,vnx+2*npd);
          acoffyz1=alloc2float(vnz+2*npd,vny+2*npd);
          acoffyz2=alloc2float(vnz+2*npd,vny+2*npd);

          zero2float(coffxy1,vny+2*npd,vnx+2*npd);
          zero2float(coffxy2,vny+2*npd,vnx+2*npd);
          zero2float(coffxz1,vnz+2*npd,vnx+2*npd);
          zero2float(coffxz2,vnz+2*npd,vnx+2*npd);
          zero2float(coffyz1,vnz+2*npd,vny+2*npd);
          zero2float(coffyz2,vnz+2*npd,vny+2*npd);

          zero2float(acoffxy1,vny+2*npd,vnx+2*npd);
          zero2float(acoffxy2,vny+2*npd,vnx+2*npd);
          zero2float(acoffxz1,vnz+2*npd,vnx+2*npd);
          zero2float(acoffxz2,vnz+2*npd,vnx+2*npd);
          zero2float(acoffyz1,vnz+2*npd,vny+2*npd);
          zero2float(acoffyz2,vnz+2*npd,vny+2*npd);
    */

/**********************************************************************************************************************/
/***************************************************  write coffe  *******************************************************************/
/**********************************************************************************************************************/
          float *coffx1;float *coffx2;float *coffy1;float *coffy2;float *coffz1;float *coffz2;
          float *acoffx1;float *acoffx2;float *acoffy1;float *acoffy2;float *acoffz1;float *acoffz2;
          coffx1=alloc1float(vnx+2*npd);
          coffx2=alloc1float(vnx+2*npd);
          coffy1=alloc1float(vny+2*npd);
          coffy2=alloc1float(vny+2*npd);
	  coffz1=alloc1float(vnz+2*npd);
          coffz2=alloc1float(vnz+2*npd);

	  acoffx1=alloc1float(vnx+2*npd);
	  acoffx2=alloc1float(vnx+2*npd);
          acoffy1=alloc1float(vny+2*npd);
	  acoffy2=alloc1float(vny+2*npd);
	  acoffz1=alloc1float(vnz+2*npd);
	  acoffz2=alloc1float(vnz+2*npd);

          zero1float(coffx1,vnx+2*npd);
          zero1float(coffx2,vnx+2*npd);
          zero1float(coffz1,vnz+2*npd);
          zero1float(coffz2,vnz+2*npd);
          zero1float(coffy1,vny+2*npd);
          zero1float(coffy2,vny+2*npd);

          zero1float(acoffx1,vnx+2*npd);
          zero1float(acoffx2,vnx+2*npd);
          zero1float(acoffz1,vnz+2*npd);
          zero1float(acoffz2,vnz+2*npd);
          zero1float(acoffy1,vny+2*npd);
          zero1float(acoffy2,vny+2*npd);


	 
          initial_coffe(dt,d0,nx,ny,nz,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2,acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2,npd);
        //  initial_coffe2(dt,d0,nx,ny,nz,coffxy1,coffxy2,coffxz1,coffxz2,coffyz1,coffyz2,acoffxy1,acoffxy2,acoffxz1,acoffxz2,acoffyz1,acoffyz2,npd);

/***********************************************************/
	   ndtt=(int)ndtt;



       //  printf("IS============%d \n",is);
/***********************************************/ 
           zero3float(p_cal,nt,ny,nx);
           zero3float(tx0,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(tx1,vnz+2*npd,vny+2*npd,vnx+2*npd);  
           zero3float(ty0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(ty1,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(tz0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(tz1,vnz+2*npd,vny+2*npd,vnx+2*npd); 

           zero3float(u0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(u1,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(v0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(v1,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(w0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(w1,vnz+2*npd,vny+2*npd,vnx+2*npd);

           zero3float(txx0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(txx1,vnz+2*npd,vny+2*npd,vnx+2*npd);
       
           zero3float(txx_snap,vnz,vny,vnx);

           zero3float(vp,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(rho,vnz+2*npd,vny+2*npd,vnx+2*npd);
 

			 
           current_shot(vp0,rho0,vp,rho,nx,ny,nz,npd,vnx,vny,vnz,is);//*********there is no transform the vp and  rho

           pad_vv(nx,ny,nz,npd,vp); 

           pad_vv(nx,ny,nz,npd,rho);
       
      
	   for(i=0;i<=nx+2*npd-1;i++)
	   {
	     for(j=0;j<=ny+2*npd-1;j++)
	     {
                for(k=0;k<=nz+2*npd-1;k++)
                {  
		   vp[i][j][k]=rho[i][j][k]*(vp[i][j][k]*vp[i][j][k]);
                   rho[i][j][k]=1.0/rho[i][j][k];
		}		
	     }
	   }
		
	   t=0.0;

    for(it=0;it<nt;it++)
    { 
          if(it%50==0)
          printf("is==%3d/%d  ,  it===%4d/%d\n",is,ns_sxd,it,nt);
	  t=t+dt;
	  ptsource(pfac,xs_sxd,ys_sxd,zs_sxd,nx,ny,nz,dt,t,favg,s,wtype,pi,npd,is);
          update_vel(nx,ny,nz,npd,mm,dt,dx,dy,dz,u0,v0,w0,txx0,u1,v1,w1,txx1,rho,c,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
          update_txx(nx,ny,nz,dt,dx,dy,dz,mm,u0,v0,w0,txx0,u1,v1,w1,txx1,s,vp,c,npd,tx1,tx0,ty1,ty0,tz1,tz0,acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2);
       

          
	  for(i=npd;i<npd+nx;i++)  
	  {   
              for(j=npd;j<npd+ny;j++)
              {
		p_cal[i-npd][j-npd][it]=txx1[i][j][npd+hsx-1];
              }
	  }
          

	  for(i=0;i<nx+2*npd;i++)
	  {
		for(j=0;j<ny+2*npd;j++)
		{
                     for(k=0;k<nz+2*npd;k++)
                     {
			u0[i][j][k]=u1[i][j][k];
                        v0[i][j][k]=v1[i][j][k];
			w0[i][j][k]=w1[i][j][k];
			tx0[i][j][k]=tx1[i][j][k];
                        ty0[i][j][k]=ty1[i][j][k];
			tz0[i][j][k]=tz1[i][j][k];			
			txx0[i][j][k]=txx1[i][j][k];
                     }
		}
	   }

           if(it==150&&is==1)
           {
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it150.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }

           if(it==200&&is==1)
           {
              zero3float(txx_snap,vnz,vny,vnx);
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it200.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }
           
           if(it==250&&is==1)
           {
              zero3float(txx_snap,vnz,vny,vnx);
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it250.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }

           if(it==300&&is==1)
           {
              zero3float(txx_snap,vnz,vny,vnx);
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it300.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }

           if(it==350&&is==1)
           {
              zero3float(txx_snap,vnz,vny,vnx);
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it350.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }

           if(it==400&&is==1)
           {
              zero3float(txx_snap,vnz,vny,vnx);
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it400.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }

           
   
          

         
	

          

     }
        
/***********************************************/        
      	  free3float(vp0);
          free3float(rho0);
      	  free3float(vp);
          free3float(rho);

          free3float(u0);
          free3float(v0);
          free3float(w0);
          free3float(txx0);
          free3float(u1);
          free3float(v1);
          free3float(w1);
          free3float(txx1);

          free3float(txx_snap);

          free3float(tx0);
          free3float(tx1);
          free3float(ty0);
          free3float(ty1);
          free3float(tz0);
          free3float(tz1);
		  
	  free3float(s);
       
          free1float(coffx1);free1float(coffx2);
          free1float(coffz1);free1float(coffz2);
          free1float(coffy1);free1float(coffy2);
          free1float(acoffx1);free1float(acoffx2);
          free1float(acoffz1);free1float(acoffz2);
          free1float(acoffy1);free1float(acoffy2);

  /*        free2float(coffxy1);free2float(coffxy2);
          free2float(coffxz1);free2float(coffxz2);
          free2float(coffyz1);free2float(coffyz2);
          free2float(acoffxy1);free2float(acoffxy2);
          free2float(acoffxz1);free2float(acoffxz2);
          free2float(acoffyz1);free2float(acoffyz2);

   */   
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




void ptsource(float pfac,float xsn,float ysn,float zsn,int nx,int ny,int nz,float dt,float t,float favg,float ***s,int wtype,float pi,int npd,int is)
{

	    float get_source(float ts,float favg,int wtype);

	    int i,j,k,ixs,iys,izs,x,y,z;
	    float tdelay,ts,source,fs;
           
            

	    for(i=0;i<=nx+2*npd-1;i++)
	    {
		   for(j=0;j<=ny+2*npd-1;j++)
		   {
                       for(k=0;k<nz+2*npd;k++)
			  s[i][j][k]=0.0*dt*favg;
		   }
	    }

	    tdelay=1.0/favg;
            ts=t-tdelay;
            
           // fs=xsn;
       if(t<=2*tdelay)
       {
            source=get_source(ts,favg,wtype);            
	    ixs = (int)(xsn+0.5)+npd-1;
            iys = (int)(ysn+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;

            for(i=ixs-3;i<=ixs+3;i++)
	    { 
		  for(j=iys-3;j<=iys+3;j++)
		  {  
                      for(k=izs-3;k<=izs+3;k++)
                      {
			  x=i-ixs;
                          y=j-iys;
                          z=k-izs;
                          s[i][j][k]=pfac*source*exp(-z*z-y*y-x*x);
                      }
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



void update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,float ***u0,float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***rho,float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,float *coffz2)
{
		 int ii,i,j,k;
		 float dtxx,dtxz,dtx,dty,dtz;


                 

		 dtx=dt/dx;
                 dty=dt/dy;
		 dtz=dt/dz;
         
		 for(k=mm;k<=(2*npd+nz-mm-1);k++)
		 {
		   for(j=mm;j<=(2*npd+ny-mm-1);j++)
	           {
                       for(i=mm;i<=(2*npd+nx-mm-1);i++)
                       {
			   u1[i][j][k]=coffx2[i]*u0[i][j][k]-coffx1[i]*dtx*rho[i][j][k]*(c[0]*(txx0[i+1][j][k]-txx0[i][j][k])+c[1]*(txx0[i+2][j][k]-txx0[i-1][j][k])+c[2]*(txx0[i+3][j][k]-txx0[i-2][j][k])+c[3]*(txx0[i+4][j][k]-txx0[i-3][j][k]));
                       }
		   }
		 }

                 for(k=mm;k<=(2*npd+nz-mm-1);k++)
		 {
		   for(i=mm;i<=(2*npd+nx-mm-1);i++)
	           {
                       for(j=mm;j<=(2*npd+ny-mm-1);j++)
                       {
			   v1[i][j][k]=coffy2[j]*v0[i][j][k]-coffy1[j]*dty*rho[i][j][k]*(c[0]*(txx0[i][j+1][k]-txx0[i][j][k])+c[1]*(txx0[i][j+2][k]-txx0[i][j-1][k])+c[2]*(txx0[i][j+3][k]-txx0[i][j-2][k])+c[3]*(txx0[i][j+4][k]-txx0[i][j-3][k]));
                       }
		   }
		 }

         
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 { 
		   for(j=mm;j<=(2*npd+ny-mm-1);j++)
		   {
                       for(k=mm;k<=(2*npd+nz-mm-1);k++)
                       {
			   w1[i][j][k]=coffz2[k]*w0[i][j][k]-coffz1[k]*dtz*rho[i][j][k]*(c[0]*(txx0[i][j][k+1]-txx0[i][j][k])+c[1]*(txx0[i][j][k+2]-txx0[i][j][k-1])+c[2]*(txx0[i][j][k+3]-txx0[i][j][k-2])+c[3]*(txx0[i][j][k+4]-txx0[i][j][k-3]));
                       }
		   }
		 }
		

}


void update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***s,float ***vp,float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,float ***tz1,float ***tz0,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2)
{

		 float dux,dwz;
		 int i,j,k,ii;
		 float dtx,dty,dtz;


                 

		 dtx=dt/dx;
                 dty=dt/dy;
		 dtz=dt/dz;

		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
		    for(j=mm;j<=(2*npd+ny-mm-1);j++)
		    {
		       for(k=mm;k<=(2*npd+nz-mm-1);k++)
                       {
                            
		    	   tx1[i][j][k]=acoffx2[i]*tx0[i][j][k]-acoffx1[i]*vp[i][j][k]*dtx*(c[0]*(u1[i][j][k]-u1[i-1][j][k])+c[1]*(u1[i+1][j][k]-u1[i-2][j][k])+c[2]*(u1[i+2][j][k]-u1[i-3][j][k])+c[3]*(u1[i+3][j][k]-u1[i-4][j][k]));

                           ty1[i][j][k]=acoffy2[j]*ty0[i][j][k]-acoffy1[j]*vp[i][j][k]*dty*(c[0]*(v1[i][j][k]-v1[i][j-1][k])+c[1]*(v1[i][j+1][k]-v1[i][j-2][k])+c[2]*(v1[i][j+2][k]-v1[i][j-3][k])+c[3]*(v1[i][j+3][k]-v1[i][j-4][k]));

                           tz1[i][j][k]=acoffz2[k]*tz0[i][j][k]-acoffz1[k]*vp[i][j][k]*dtz*(c[0]*(w1[i][j][k]-w1[i][j][k-1])+c[1]*(w1[i][j][k+1]-w1[i][j][k-2])+c[2]*(w1[i][j][k+2]-w1[i][j][k-3])+c[3]*(w1[i][j][k+3]-w1[i][j][k-4]));
                    
			   txx1[i][j][k]=tx1[i][j][k]+ty1[i][j][k]+tz1[i][j][k]+s[i][j][k];
                       }
	            }
		 }

}           


           

float get_constant(float dx,float dy,float dz,int nx,int ny,int nz,int nt,int ntout,int npd,float tmax,float favg,float dtout,float dt,float ***vp0,float ndtt)
{
		 int i,j,k;
		 float vpmax,vpmin,H_min;
		 float dt_max,dx_max,dz_max,dy_max,d0;


                 

		 vpmax=vp0[npd][npd][npd];
		 vpmin=vp0[npd][npd][npd];
                 for(j=npd;j<ny+npd;j++)
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(k=npd;k<nz+npd;k++)
			 {
				 if(vpmax<vp0[i][j][k]) vpmax=vp0[i][j][k];
				 if(vpmin>vp0[i][j][k]) vpmin=vp0[i][j][k];
			 }
		 }

             
		 d0=3.0*vpmax*log(100000.0)/(2.0*npd*dx);

                 H_min=dx;
		 if(dy<H_min) H_min=dy;
                 if(dz<H_min) H_min=dz;
		 


/*====== determine time sampling interval to ensure stability====*/

		 dt_max=0.5*1000*H_min/vpmax;
                 dx_max=vpmin/favg*0.2;
                 dy_max=dx_max;
                 dz_max=dx_max;

             
                if(dx_max<dx)
                { 
                   printf("dx_max===%f, vpmin===%f, favg===%f \n",dx_max,vpmin,favg);
		   printf("YOU NEED HAVE TO REDEFINE DX ! \n");
                   exit(0);
		 }
                 if(dy_max<dy)
                { 
                   printf("dy_max===%f, vpmin===%f, favg===%f \n",dy_max,vpmin,favg);
		   printf("YOU NEED HAVE TO REDEFINE DY ! \n");
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


void pad_vv(int nx,int ny,int nz,int npd,float ***ee)
{
		 int i,j,k;

                   

/*****pad left side                    */
          for(j=npd;j<npd+ny;j++)
            for(k=npd;k<nz+npd;k++)
	    {
              for(i=0;i<npd;i++)
              { 
                 ee[i][j][k]=ee[npd][j][k];
              }
	    }
       
/*****pad right side                    */
          for(j=npd;j<npd+ny;j++)
             for(k=npd;k<nz+npd;k++)
	     {
                for(i=nx+npd;i<nx+2*npd;i++)
                 {
                  ee[i][j][k]=ee[nx+npd-1][j][k];
                 }
	     }
/*****pad front side                    */
          for(i=0;i<npd*2+nx;i++)
             for(k=npd;k<nz+npd;k++)
	     {
                for(j=0;j<npd;j++)
                 {
                  ee[i][j][k]=ee[i][npd][k];
                 }
	     }
/*****pad back side                    */
          for(i=0;i<npd*2+nx;i++)
             for(k=npd;k<nz+npd;k++)
	     {
                for(j=ny+npd;j<ny+2*npd;j++)
                 {
                  ee[i][j][k]=ee[i][ny+npd-1][k];
                 }
	     }

/*****pad upper side                    */
            for(i=0;i<nx+2*npd;i++)
	    {
              for(j=0;j<ny+2*npd;j++)
              {
                for(k=0;k<npd;k++)
                   ee[i][j][k]=ee[i][j][npd];
              }
	    }
/*****lower side                        */
            for(i=0;i<nx+2*npd;i++)
	    {
              for(j=0;j<ny+2*npd;j++)
              {
                for(k=npd+nz;k<npd*2+nz;k++)
                   ee[i][j][k]=ee[i][j][npd+nz-1];
              }
	    }
		
		
}

void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd)
{

		 int i,j,k;
		

                 

		 FILE *fp1;
		 fp1=fopen(FN1,"rb");
                 for(j=npd;j<ny+npd;j++)
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(k=npd;k<nz+npd;k++)
			 {
				 fread(&vv[i][j][k],4,1,fp1);

			 }
		 }
                 for(j=npd;j<ny+npd;j++)
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(k=npd;k<nz+npd;k++)
			 {
				 rho0[i][j][k]=1.0;
			 }
		 }
		 fclose(fp1);
}


void current_shot(float ***vp0,float ***rho0,float ***vp,float ***rho,int nx,int ny,int nz,int npd,int vnx,int vny,int vnz,int is)
{
      
            
                         int ivstart,ivend;
			 int i,ix,iy,iz;
                         is=1;

                        

                     //    ivstart=1+(is-1)*ds_sxd;
		    //	 ivend=nx+(is-1)*ds_sxd;
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
                         for(iy=npd;iy<ny+npd;iy++)
			 for(ix=npd;ix<nx+npd;ix++)
			 {
				 for(iz=npd;iz<nz+npd;iz++)
				 {
				  vp[ix][iy][iz]=vp0[ix][iy][iz];
                                  rho[ix][iy][iz]=rho0[ix][iy][iz];
				  
				 }
			 }

}
/**************************************************************************************************************/
void initial_coffe(float dt,float d0,int nx,int ny,int nz,float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,float *coffz2,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2,int npd)
{		
		 int i,j,k;
         

		 for(i=0;i<npd;i++)
		 {   
			 coffx1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);
                         coffy1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffy2[i]=coffy1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz2[i]=coffz1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);

	
		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 coffx1[i]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
		 }
        
                 for(i=npd+ny;i<ny+2*npd;i++)
		 {
			 coffy1[i]=1/(1+(dt*d0*pow((0.5+i-ny-npd)/npd,2))/2);
			 coffy2[i]=coffy1[i]*(1-(dt*d0*pow((0.5+i-ny-npd)/npd,2))/2);
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
                 for(i=npd;i<npd+ny;i++)
		 {
			 coffy1[i]=1.0;
			 coffy2[i]=1.0;
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
                         acoffy1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffy2[i]=acoffy1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz2[i]=acoffz1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);

		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 acoffx1[i]=1/(1+(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
			 acoffx2[i]=acoffx1[i]*(1-(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
		 }
                 for(i=npd+ny;i<ny+2*npd;i++)
		 {
			 acoffy1[i]=1/(1+(dt*d0*pow(((1+i-ny-npd)*1.0)/npd,2))/2);
			 acoffy2[i]=acoffy1[i]*(1-(dt*d0*pow(((1+i-ny-npd)*1.0)/npd,2))/2);
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
                 for(i=npd;i<npd+ny;i++)
		 {
			 acoffy1[i]=1.0;
			 acoffy2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 acoffz1[i]=1.0;
			 acoffz2[i]=1.0;
		 }
			       
}

/**************************************************************************************************************/
/********************************************   wrong coffe  ******************************************************************/
/**************************************************************************************************************/
/*
void initial_coffe2(float dt,float d0,int nx,int ny,int nz,float **coffxy1,float **coffxy2,float **coffxz1,float **coffxz2,float **coffyz1,float **coffyz2,float **acoffxy1,float **acoffxy2,float **acoffxz1,float **acoffxz2,float **acoffyz1,float **acoffyz2,int npd)
{		
		 int i,j,k;
         
        */          

/******************************************************************************************************************/	
/****************************************************  coff start   ***********************************************/
/******************************************************************************************************************/
/*		 for(i=0;i<npd;i++)//************all coff
		 {   
                     for(j=0;j<npd;j++)
                     {
			 coffxy1[i][j]=1/(1+(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);

                         coffxz1[i][j]=1/(1+(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
			 coffxz2[i][j]=coffxz1[i][j]*(1-(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);

			 coffyz1[i][j]=1/(1+(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
			 coffyz2[i][j]=coffyz1[i][j]*(1-(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
                     }	
		 }

*/
/*******************************************************  coffxy   *******************************************/
/*		 for(i=npd+nx;i<nx+2*npd;i++) 
		 {
                     for(j=npd+ny;j<ny+2*npd;j++)
                     {
			 coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i+j-nx-ny-2*npd)/(npd),2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+i+j-nx-ny-2*npd)/(npd),2))/2);
                     }
                     for(j=0;j<npd;j++)
                     {
                         coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i-j-nx)/npd,2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+i-j-nx)/npd,2))/2);
                     }
                     for(j=npd;j<ny+npd;j++)
                     {
                         coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
                     }
		 }
                 for(j=ny+npd;j<ny+2*npd;j++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+j-i-ny)/npd,2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+j-i-ny)/npd,2))/2);
                     }
                     for(i=npd;i<nx+npd;i++)
                     {
                         coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(j=0;j<npd;j++)
                     {
                         coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+j-npd)/npd,2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+j-npd)/npd,2))/2);
                     }
                 }
                 for(j=npd;j<npd+ny;j++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         coffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i-npd)/npd,2))/2);
			 coffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow((0.5+i-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(j=npd;j<npd+ny;j++)
                     {
                         coffxy1[i][j]=1.0;
			 coffxy2[i][j]=1.0;
                     }
                 }
*/
/*******************************************************  coffxy over  *******************************************/

/*******************************************************  coffxz   *******************************************/
/*		 for(i=npd+nx;i<nx+2*npd;i++) 
		 {
                     for(k=npd+nz;k<nz+2*npd;k++)
                     {
			 coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i+k-nx-nz-2*npd)/(npd),2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+i+k-nx-nz-2*npd)/(npd),2))/2);
                     }
                     for(k=0;k<npd;k++)
                     {
                         coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i-k-nx)/npd,2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+i-k-nx)/npd,2))/2);
                     }
                     for(k=npd;k<nz+npd;k++)
                     {
                         coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
                     }
		 }
                 for(k=nz+npd;k<nz+2*npd;k++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+k-i-nz)/npd,2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+k-i-nz)/npd,2))/2);
                     }
                     for(i=npd;i<nx+npd;i++)
                     {
                         coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(k=0;k<npd;k++)
                     {
                         coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+k-npd)/npd,2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+k-npd)/npd,2))/2);
                     }
                 }
                 for(k=npd;k<npd+nz;k++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         coffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i-npd)/npd,2))/2);
			 coffxz2[i][k]=coffxz1[i][k]*(1-(dt*d0*pow((0.5+i-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(k=npd;k<npd+nz;k++)
                     {
                         coffxz1[i][k]=1.0;
			 coffxz2[i][k]=1.0;
                     }
                 }
*/
/*******************************************************  coffxz over  *******************************************/

/*******************************************************  coffyz   *******************************************/
/*		 for(j=npd+ny;j<ny+2*npd;j++) 
		 {
                     for(k=npd+nz;k<nz+2*npd;k++)
                     {
			 coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j+k-ny-nz-2*npd)/(npd),2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+j+k-ny-nz-2*npd)/(npd),2))/2);
                     }
                     for(k=0;k<npd;k++)
                     {
                         coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j-k-ny)/npd,2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+j-k-ny)/npd,2))/2);
                     }
                     for(k=npd;k<nz+npd;k++)
                     {
                         coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
                     }
		 }
                 for(k=nz+npd;k<nz+2*npd;k++)
                 {
                     for(j=0;j<npd;j++)
                     {
                         coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+k-j-nz)/npd,2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+k-j-nz)/npd,2))/2);
                     }
                     for(j=npd;j<ny+npd;j++)
                     {
                         coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
                     }
                 }
                 for(j=npd;j<npd+ny;j++)
                 {
                     for(k=0;k<npd;k++)
                     {
                         coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+k-npd)/npd,2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+k-npd)/npd,2))/2);
                     }
                 }
                 for(k=npd;k<npd+nz;k++)
                 {
                     for(j=0;j<npd;j++)
                     {
                         coffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j-npd)/npd,2))/2);
			 coffyz2[j][k]=coffyz1[j][k]*(1-(dt*d0*pow((0.5+j-npd)/npd,2))/2);
                     }
                 }
                 for(j=npd;j<npd+ny;j++)
                 {
                     for(k=npd;k<npd+nz;k++)
                     {
                         coffyz1[j][k]=1.0;
			 coffyz2[j][k]=1.0;
                     }
                 }
*/
/*******************************************************  coffyz over  *******************************************/
		


/******************************************************************************************************************/	
/****************************************************  acoff start   ***********************************************/
/******************************************************************************************************************/
/*		 for(i=0;i<npd;i++)//************all acoff
		 {   
                     for(j=0;j<npd;j++)
                     {
			 acoffxy1[i][j]=1/(1+(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);

                         acoffxz1[i][j]=1/(1+(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
			 acoffxz2[i][j]=acoffxz1[i][j]*(1-(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);

			 acoffyz1[i][j]=1/(1+(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
			 acoffyz2[i][j]=acoffyz1[i][j]*(1-(dt*d0*pow((2*npd-0.5-i-j)/(npd),2))/2);
                     }	
		 }

*/
/*******************************************************  acoffxy   *******************************************/
/*		 for(i=npd+nx;i<nx+2*npd;i++) 
		 {
                     for(j=npd+ny;j<ny+2*npd;j++)
                     {
			 acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i+j-nx-ny-2*npd)/(npd),2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+i+j-nx-ny-2*npd)/(npd),2))/2);
                     }
                     for(j=0;j<npd;j++)
                     {
                         acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i-j-nx)/npd,2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+i-j-nx)/npd,2))/2);
                     }
                     for(j=npd;j<ny+npd;j++)
                     {
                         acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
                     }
		 }
                 for(j=ny+npd;j<ny+2*npd;j++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+j-i-ny)/npd,2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+j-i-ny)/npd,2))/2);
                     }
                     for(i=npd;i<nx+npd;i++)
                     {
                         acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(j=0;j<npd;j++)
                     {
                         acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+j-npd)/npd,2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+j-npd)/npd,2))/2);
                     }
                 }
                 for(j=npd;j<npd+ny;j++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         acoffxy1[i][j]=1/(1+(dt*d0*pow((0.5+i-npd)/npd,2))/2);
			 acoffxy2[i][j]=acoffxy1[i][j]*(1-(dt*d0*pow((0.5+i-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(j=npd;j<npd+ny;j++)
                     {
                         acoffxy1[i][j]=1.0;
			 acoffxy2[i][j]=1.0;
                     }
                 }
*/
/*******************************************************  acoffxy over  *******************************************/

/*******************************************************  acoffxz   *******************************************/
/*		 for(i=npd+nx;i<nx+2*npd;i++) 
		 {
                     for(k=npd+nz;k<nz+2*npd;k++)
                     {
			 acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i+k-nx-nz-2*npd)/(npd),2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+i+k-nx-nz-2*npd)/(npd),2))/2);
                     }
                     for(k=0;k<npd;k++)
                     {
                         acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i-k-nx)/npd,2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+i-k-nx)/npd,2))/2);
                     }
                     for(k=npd;k<nz+npd;k++)
                     {
                         acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
                     }
		 }
                 for(k=nz+npd;k<nz+2*npd;k++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+k-i-nz)/npd,2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+k-i-nz)/npd,2))/2);
                     }
                     for(i=npd;i<nx+npd;i++)
                     {
                         acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(k=0;k<npd;k++)
                     {
                         acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+k-npd)/npd,2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+k-npd)/npd,2))/2);
                     }
                 }
                 for(k=npd;k<npd+nz;k++)
                 {
                     for(i=0;i<npd;i++)
                     {
                         acoffxz1[i][k]=1/(1+(dt*d0*pow((0.5+i-npd)/npd,2))/2);
			 acoffxz2[i][k]=acoffxz1[i][k]*(1-(dt*d0*pow((0.5+i-npd)/npd,2))/2);
                     }
                 }
                 for(i=npd;i<npd+nx;i++)
                 {
                     for(k=npd;k<npd+nz;k++)
                     {
                         acoffxz1[i][k]=1.0;
			 acoffxz2[i][k]=1.0;
                     }
                 }
*/
/*******************************************************  acoffxz over  *******************************************/

/*******************************************************  acoffyz   *******************************************/
/*		 for(j=npd+ny;j<ny+2*npd;j++) 
		 {
                     for(k=npd+nz;k<nz+2*npd;k++)
                     {
			 acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j+k-ny-nz-2*npd)/(npd),2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+j+k-ny-nz-2*npd)/(npd),2))/2);
                     }
                     for(k=0;k<npd;k++)
                     {
                         acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j-k-ny)/npd,2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+j-k-ny)/npd,2))/2);
                     }
                     for(k=npd;k<nz+npd;k++)
                     {
                         acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+j-ny-npd)/npd,2))/2);
                     }
		 }
                 for(k=nz+npd;k<nz+2*npd;k++)
                 {
                     for(j=0;j<npd;j++)
                     {
                         acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+k-j-nz)/npd,2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+k-j-nz)/npd,2))/2);
                     }
                     for(j=npd;j<ny+npd;j++)
                     {
                         acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+k-nz-npd)/npd,2))/2);
                     }
                 }
                 for(j=npd;j<npd+ny;j++)
                 {
                     for(k=0;k<npd;k++)
                     {
                         acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+k-npd)/npd,2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+k-npd)/npd,2))/2);
                     }
                 }
                 for(k=npd;k<npd+nz;k++)
                 {
                     for(j=0;j<npd;j++)
                     {
                         acoffyz1[j][k]=1/(1+(dt*d0*pow((0.5+j-npd)/npd,2))/2);
			 acoffyz2[j][k]=acoffyz1[j][k]*(1-(dt*d0*pow((0.5+j-npd)/npd,2))/2);
                     }
                 }
                 for(j=npd;j<npd+ny;j++)
                 {
                     for(k=npd;k<npd+nz;k++)
                     {
                         acoffyz1[j][k]=1.0;
			 acoffyz2[j][k]=1.0;
                     }
                 }
*/
/*******************************************************  acoffyz over  *******************************************/
	 
/*		 for(i=0;i<npd;i++)    //******all acoff
		 {    
                     for(j=0;j<npd;i++)
                     {
			 acoffxy1[i][j]=1/(1+(dt*d0*pow(((2*npd-i-j)*1.0)/npd,2))/2);
			 acoffxy2[i][j]=coffxy1[i][j]*(1-(dt*d0*pow(((2*npd-i-j)*1.0)/npd,2))/2);
                         acoffxz1[i][j]=1/(1+(dt*d0*pow(((2*npd-i-j)*1.0)/npd,2))/2);
			 acoffxz2[i][j]=coffxz1[i][j]*(1-(dt*d0*pow(((2*npd-i-j)*1.0)/npd,2))/2);
			 acoffyz1[i][j]=1/(1+(dt*d0*pow(((2*npd-i-j)*1.0)/npd,2))/2);
			 acoffyz2[i][j]=coffyz1[i][j]*(1-(dt*d0*pow(((2*npd-i-j)*1.0)/npd,2))/2);
                     }
		 }
*/
/****************************************************  acoff  ***********************************************/
/*		 for(i=npd+nx;i<nx+2*npd;i++)
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
*/
/*		       
}
*/

