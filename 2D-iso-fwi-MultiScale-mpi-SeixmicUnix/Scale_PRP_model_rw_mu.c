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


        void readfile(int nx,int nt,int is,char FN1[],float **p_real);
        void model_step(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,float favg,
                        float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
                        int is,float **p_cal,int wavelet);
	  void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);        
        void freea(float **a,int m);
        void zero(float **a,int nx,int nz);
              

        char FN1[250]={"fault_vel_280_180.dat"};
  
	  char FN3[250]={"fault_shot_obs_280_3000_250_rw.dat"};

        char FN4[250]={"Temp_vel_for_mute.dat"};
        

        FILE *fp3;
        fp3=fopen(FN3,"wb");
        
	int i,j,k,is,nx,nz,nt,vnx,vnz,i_start,i_end,mz,vmz,wavelet;
	int ns_sxd,ds_sxd,fs_sxd,zs_sxd,npd;
	
	float dx,dz,vdx,vdz,tmax,dt,dtout,pfac,favg;
	
	int myid,numprocs,count;

/*************parameter****************/  

        nx=280;         npd=50;   tmax=1.8;
	  nz=180;         favg=15;  pfac=250.0;
 	                  dx=10.0;   
	  vnx=280;        dz=10.0;   
	  vnz=180;        vdx=10.0;
	  nt=3000;        vdz=10.0;
        ns_sxd=1;       dt=0.6;
        ds_sxd=5;       dtout=0.6;
        fs_sxd=40;
        zs_sxd=1;        
 
        mz=30;//++++for cute  

           wavelet=1;//wavelet==0:use initial wavelet
                     //wavelet==1:read wavelet


        vmz=mz;

        float **vtempor,**vmute;

          vtempor=(float **)calloc((vnx),sizeof(float *));        
	  {
		  for(i=0;i<vnx;i++)                           
		  {
		  vtempor[i]=(float *)calloc((vnz),sizeof(float)); 
		  }
	  }   
        vmute=(float **)calloc((vnx),sizeof(float *));        
	  {
		  for(i=0;i<vnx;i++)                           
		  {
		  vmute[i]=(float *)calloc((vmz),sizeof(float)); 
		  }
	  }      
                 FILE *fptempor;
		 fptempor=fopen(FN1,"rb");
		 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<nz;j++)
			 {
				 fread(&vtempor[i][j],4,1,fptempor);

			 }
		 }
                 fclose(fptempor);
                 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<mz;j++)
			 {
				 vmute[i][j]=vtempor[i][j];

			 }
		 }
                 fptempor=fopen(FN4,"wb");
		 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<mz;j++)
			 {
				 fwrite(&vmute[i][j],4,1,fptempor);

			 }
		 }
                 fclose(fptempor);
        float **p_cal;

	p_cal=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		p_cal[i]=(float *)calloc(nt,sizeof(float));
	}  
        float **p_m;

	p_m=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		p_m[i]=(float *)calloc(nt,sizeof(float));
	} 
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
             printf(" begin the loop  !  \n");
        }
 
/***************************************/
        zero(p_cal,nx,nt);
        zero(p_m,nx,nt);
/***************************************/
   
        if(myid==0)
	{
	  printf("1. the model is start  !  \n");
	}
	
        model_step(nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,FN1,
                     ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_cal,wavelet);   
         if(myid==0)
	{
	  printf("2. start mute the direct wave  !\n");
	}
        model_step(nx,mz,vnx,vmz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,FN4,
                     ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_m,wavelet);//++++++  
           
        for(i=0;i<nx;i++)//++++++++++++
	{ 
	  for(j=0;j<nt;j++)
	  {
	    p_cal[i][j]=p_cal[i][j]-p_m[i][j];
	  }
	}
        
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
	  printf("3. the model is over    \n");
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
      freea(p_cal,nx);
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

void model_step(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,
                float favg,float tmax,float dt,float dtout,float pfac,char FN1[],
                 int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,int is,float **p_cal,int wavelet)

{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                    int wtype,float pi,int npd,int is,int ds_sxd);
      void ptsource_wavelet(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                    int wtype,float pi,int npd,int is,int ds_sxd,int nt,int it);//++++++++++++++
      void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                     float **u0,float **w0,float **txx0,float **u1,float **w1,float **txx1,
                     float **rho,float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2);
      void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                     float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,
                     float **tx1,float **tx0,float **tz1,float **tz0,float *acoffx1,float *acoffx2,float *acoffz1,
                     float *acoffz2); 
      void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,
                     float favg,float dtout,float dt,float **vp0,float ndtt);
      void zero(float **a,int nx,int nz);
      void zero2d(int nx,int nz,float **vv0,int npd);
	void zero2d2(int nx,int nz,int npd,float **vv);
	void zero2d3(int nx,int nz,int npd,float **vv);
	void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,
                  float **u1,float **w1,float **txx1,float **vp,float **rho);
      void pad_vv(int nx,int nz,int npd,float **ee);
	void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);
	void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,
                          int npd,int vnx,int vnz,int ds_sxd,int is);
	
	void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                                  float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);	
      void freea(float **a,int m);
      


	  int mm=4; 
          int hsx=1;
	  int i,j;

	  int ntout,wtype,it,ifx,ilx,jfz,jlz;
	  float pi,t,ndtt,d0;
	  float iabsorb[4];
	  
	  FILE *fp1;

          wtype=1;
          iabsorb[0]=1.0;
          iabsorb[1]=1.0;
          iabsorb[2]=1.0;
          iabsorb[3]=1.0;

	  pi=3.1415927;
	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;
     
	  float **vp0;
	  float **rho0;
	  float **vp;
	  float **rho;
	  float **u0;
	  float **w0;
	  float **txx0; 
	  float **u1;
	  float **w1;
	  float **txx1;
	  float **tx0;
	  float **tx1;
	  float **tz0;
	  float **tz1;
	  float **s;
     
	  float c[4];
	  

          cal_c(mm,c);   

      
	  vp0=(float **)calloc((vnx+2*npd),sizeof(float *));        
	  {
		  for(i=0;i<(vnx+2*npd);i++)                           
		  {
		  vp0[i]=(float *)calloc((vnz+2*npd),sizeof(float)); 
		  }
	  }                                                 
     

    
	  rho0=(float **)calloc((vnx+2*npd),sizeof(float *));       
	  {
		  for(i=0;i<(vnx+2*npd);i++)                           
		  {
		  rho0[i]=(float *)calloc((vnz+2*npd),sizeof(float)); 
		  }
	  }                                               

     

          zero2d(vnx,vnz,vp0,npd);                     
          zero2d(vnx,vnz,rho0,npd);                    
          read_file(FN1,vnx,vnz,vp0,rho0,npd);          
      
      
	  vp=(float **)calloc((nz+2*npd),sizeof(float *));       
	  {
		  for(i=0;i<nz+2*npd;i++)                            
		  {
		  vp[i]=(float *)calloc((nx+2*npd),sizeof(float));  
		  }
	  }                                                

      
	  rho=(float **)calloc((nz+2*npd),sizeof(float *));       
	  {
		  for(i=0;i<nz+2*npd;i++)                          
		  {
		  rho[i]=(float *)calloc((nx+2*npd),sizeof(float)); 
		  }
	  }                                              


	 u0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }
	  

	 u1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

      
	 w0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }
	  

	 w1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }                                                                         


	 txx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 txx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }
	 
	 tx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }
	 
	 tx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }
	 
	 tz0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }
	 
	 tz1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }	 	 	 	 	 	 	 
	 
	 

      
	  s=(float **)calloc((nz+2*npd),sizeof(float *));       
	  { 
		  for(i=0;i<nz+2*npd;i++)                          
		  {
		  s[i]=(float *)calloc((nx+2*npd),sizeof(float)); 
		  }
	  }   


          d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);
          dt=dt/1000;
/***********************************************************/

          float *coffx1;float *coffx2;float *coffz1;float *coffz2;float *acoffx1;float *acoffx2;float *acoffz1;float *acoffz2;
          coffx1=(float *)calloc((nx+2*npd),sizeof(float));
          coffx2=(float *)calloc((nx+2*npd),sizeof(float));
	  coffz1=(float *)calloc((nz+2*npd),sizeof(float));  
          coffz2=(float *)calloc((nz+2*npd),sizeof(float));
	  acoffx1=(float *)calloc((nx+2*npd),sizeof(float));
	  acoffx2=(float *)calloc((nx+2*npd),sizeof(float));
	  acoffz1=(float *)calloc((nz+2*npd),sizeof(float));
	  acoffz2=(float *)calloc((nz+2*npd),sizeof(float));

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

          initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/***********************************************************/
	   ndtt=(int)ndtt;



       //  printf("IS============%d \n",is);
/***********************************************/  
           zero(p_cal,nx,nt); 
           zero2d2(nx,nz,npd,tx0);
	   zero2d2(nx,nz,npd,tx1);
	   zero2d2(nx,nz,npd,tz0);
	   zero2d2(nx,nz,npd,tz1);              
	   zero3d(nx,nz,nt,npd,u0,w0,txx0,u1,w1,txx1,vp,rho);
			 
           current_shot(vp0,rho0,vp,rho,nx,nz,npd,vnx,vnz,ds_sxd,is);//

           pad_vv(nx,nz,npd,vp); 

           pad_vv(nx,nz,npd,rho);
       
      
	   for(i=0;i<=nz+2*npd-1;i++)
	   {
		for(j=0;j<=nx+2*npd-1;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
				
		}
	   }
		
	   t=0.0;

    for(it=0;it<nt;it++)
    { 

         if(it%200==0&&is==1)
                printf("it===%d\n",it);

	  t=t+dt;
         if(wavelet==0)
	    ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd);
         if(wavelet==1)
          ptsource_wavelet(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd,nt,it);//+++++++++++++++++++++++++++
          update_vel(nx,nz,npd,mm,dt,dx,dz,u0,w0,txx0,u1,w1,txx1,rho,c,coffx1,coffx2,coffz1,coffz2);
          update_txx(nx,nz,dt,dx,dz,mm,u0,w0,txx0,u1,w1,txx1,s,vp,c,npd,tx1,tx0,tz1,tz0,acoffx1,acoffx2,acoffz1,acoffz2);
       


	  for(j=npd;j<npd+nx;j++)  
	  {   
		p_cal[j-npd][it]=txx1[npd+hsx-1][j];
     
	  }


	  for(i=0;i<nz+2*npd;i++)
	  {
		for(j=0;j<nx+2*npd;j++)
		{
			u0[i][j]=u1[i][j];
			w0[i][j]=w1[i][j];
			tx0[i][j]=tx1[i][j];
			tz0[i][j]=tz1[i][j];			
			txx0[i][j]=txx1[i][j];
		}
	   }

          
/*		
          if(nz>90)
          {
           FILE *fp3;
           fp3=fopen("snap.dat","wb");
            if(is==1)
	   {

		if(it%20==0)  
                 {

                    fseek(fp3,(int)(it/20)*(nx+2*npd)*(nz+2*npd)*4L,0);
                    for(i=0;i<(nx+2*npd);i++)
		    {
                        for(j=0;j<(nz+2*npd);j++)
                        {
                            fwrite(&txx1[i][j],4L,1,fp3);
                        }
                    }

                }

            }
            fclose(fp3);
            }

*/
     }
        
/***********************************************/        
      	  
	  freea(vp0,(vnx+2*npd));
	  freea(rho0,(vnx+2*npd));
	  freea(vp,(nz+2*npd));
	  freea(rho,(nz+2*npd));
	  freea(u0,(nz+2*npd));
	  freea(w0,(nz+2*npd));
	  freea(txx0,(nz+2*npd));
	  freea(u1,(nz+2*npd));
	  freea(w1,(nz+2*npd));
	  freea(txx1,(nz+2*npd));
	  freea(tx0,(nz+2*npd));
	  freea(tx1,(nz+2*npd));
	  freea(tz0,(nz+2*npd));
	  freea(tz1,(nz+2*npd));	  
	  freea(s,(nz+2*npd));
	  free(coffx1);free(coffx2);free(coffz1);free(coffz2);
	  free(acoffx1);free(acoffx2);free(acoffz1);free(acoffz2);
	    
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




void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,
                float **s,int wtype,float pi,int npd,int is,int ds_sxd)
{
            if((t==dt)&&(is==1))
               printf("we use the initial source wavelet !\n");

	    float get_source(float ts,float favg,int wtype);

	    int i,j,ixs,izs,x,z;
	    float tdelay,ts,source,fs;
      
	    for(i=0;i<=nz+2*npd-1;i++)
	    {
		   for(j=0;j<=nx+2*npd-1;j++)
		   {
			s[i][j]=0.0*dt*favg;
		   }
	    }

	    tdelay=1.0/favg;
            ts=t-2.5*tdelay;//++++++++++change
            fs=xsn+(is-1)*ds_sxd;
           if(t<=10*tdelay)//++++++++++change
           {
            source=get_source(ts,favg,wtype);  
          
	    ixs = (int)(fs+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;

       /*     for(i=izs-3;i<=izs+3;i++)
	    { 
		  for(j=ixs-3;j<=ixs+3;j++)
		  {                           *///+++++++++++++++++++change
            for(i=izs;i<=izs;i++)
	    { 
		  for(j=ixs;j<=ixs;j++)//+++++++++++++++++change
		  {  
			  x=j-ixs;z=i-izs;
                          s[i][j]=pfac*source*exp(-z*z-x*x);

                      /*    FILE *fpwave;
                          fpwave=fopen("wavelet.txt","a");
                          fprintf(fpwave,"%f \n",s[i][j]);
                          fclose(fpwave);  */
		  }
	    }
            
	}

}

void ptsource_wavelet(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,
                    float **s,int wtype,float pi,int npd,int is,int ds_sxd,int nt,int it)

{

            if((t==dt)&&(is==1))
               printf("we use the filter or smooth wavelet !\n");

	    float get_source(float ts,float favg,int wtype);
            float read_wavelet(char FN[],int nt,int it);//+++++++++++++++++++++++
	    int i,j,ixs,izs,x,z;
	    float tdelay,ts,source,fs;
            float wavelet=0.0;//+++++++++++++++++++++++++++
      
	    for(i=0;i<=nz+2*npd-1;i++)
	    {
		   for(j=0;j<=nx+2*npd-1;j++)
		   {
			s[i][j]=0.0;
		   }
	    }

	    tdelay=1.0/favg;
            ts=t-2.5*tdelay;//++++++++++change
            fs=xsn+(is-1)*ds_sxd;
           if(t<=10*tdelay)//++++++++++change
           {
            //source=get_source(ts,favg,wtype);   
            wavelet=read_wavelet("wavelet_filter.dat",nt,it);//++++++++++++++++++++
	    ixs = (int)(fs+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;

            /*     for(i=izs-3;i<=izs+3;i++)
	    { 
		  for(j=ixs-3;j<=ixs+3;j++)
		  {                           *///++++++++++++++++++++change
            for(i=izs;i<=izs;i++)
	    { 
		  for(j=ixs;j<=ixs;j++)//+++++++++++++++++change
		  {  
			  x=j-ixs;z=i-izs;
                          s[i][j]=pfac*wavelet*exp(-z*z-x*x);//++++++++++
                          //s[i][j]=pfac*wavelet;//++++++++++

                       /*   FILE *fpwave;
                          fpwave=fopen("wavelet.txt","a");
                          fprintf(fpwave,"%f \n",s[i][j]);
                          fclose(fpwave);   */
		  }
	    }
	   }

}
float read_wavelet(char FN[],int nt,int it)//+++++++++++++++++++++
{
    int i,j;
    float wavelet=0.0;
    FILE *fp1;
    fp1=fopen(FN,"rb");
    fseek(fp1,((it-1)*4L),0);//++++++you can change 
    fread(&wavelet,4,1,fp1);
    fclose(fp1);
    return (wavelet);

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



void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,float **u0,float **w0,float **txx0,
               float **u1,float **w1,float **txx1,float **rho,float c[],
               float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
		 int ii,i,j;
		 float dtxx,dtxz,dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;
         
		 for(i=mm;i<=(2*npd+nz-mm-1);i++)
		 {
			 for(j=mm;j<=(2*npd+nx-mm-1);j++)
			 {

			   u1[i][j]=coffx2[j]*u0[i][j]-coffx1[j]*dtx*rho[i][j]*(c[0]*(txx0[i][j+1]-txx0[i][j])+c[1]*(txx0[i][j+2]-txx0[i][j-1])+c[2]*(txx0[i][j+3]-txx0[i][j-2])+c[3]*(txx0[i][j+4]-txx0[i][j-3]));

			 }
		 }


         
		 for(j=mm;j<=(2*npd+nx-mm-1);j++)
		 { 
			 for(i=mm;i<=(2*npd+nz-mm-1);i++)
			 {

			   w1[i][j]=coffz2[i]*w0[i][j]-coffz1[i]*dtz*rho[i][j]*(c[0]*(txx0[i+1][j]-txx0[i][j])+c[1]*(txx0[i+2][j]-txx0[i-1][j])+c[2]*(txx0[i+3][j]-txx0[i-2][j])+c[3]*(txx0[i+4][j]-txx0[i-3][j]));


			 }
		 }
		

}


void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,float **tx1,
               float **tx0,float **tz1,float **tz0,float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2)
{

		 float dux,dwz;
		 int i,j,ii;
		 float dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;

		 for(i=mm;i<=(2*npd+nz-mm-1);i++)
		 {
			 for(j=mm;j<=(2*npd+nx-mm-1);j++)
			 {
		         
		    	   tx1[i][j]=acoffx2[j]*tx0[i][j]-acoffx1[j]*vp[i][j]*dtx*(c[0]*(u1[i][j]-u1[i][j-1])+c[1]*(u1[i][j+1]-u1[i][j-2])+c[2]*(u1[i][j+2]-u1[i][j-3])+c[3]*(u1[i][j+3]-u1[i][j-4]));

                     tz1[i][j]=acoffz2[i]*tz0[i][j]-acoffz1[i]*vp[i][j]*dtz*(c[0]*(w1[i][j]-w1[i-1][j])+c[1]*(w1[i+1][j]-w1[i-2][j])+c[2]*(w1[i+2][j]-w1[i-3][j])+c[3]*(w1[i+3][j]-w1[i-4][j]));
                    
			   txx1[i][j]=tx1[i][j]+tz1[i][j]+s[i][j];

			 }
		 }

}           


           

float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,float favg,float dtout,float dt,float **vp0,float ndtt)
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
    


void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,float **u1,float **w1,float **txx1,
            float **vp,float **rho)
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
            for(i=npd;i<=(nz+npd-1);i++)
		{
              for(j=0;j<=npd-1;j++)
              { 
               ee[i][j]=ee[i][npd];
              }
		}
       
/*****pad right side                    */
            for(i=npd;i<=(nz+npd-1);i++)
		{
              for(j=nx+npd;j<=(nx+2*npd-1);j++)
              {
                ee[i][j]=ee[i][nx+npd-1];
              }
		}
/*****pad upper side                    */
            for(i=0;i<=(npd-1);i++)
		{
              for(j=0;j<=(nx+2*npd-1);j++)
              {
                ee[i][j]=ee[npd][j];
              }
		}
/*****lower side                        */
            for(i=nz+npd;i<=(nz+2*npd-1);i++)
		{
              for(j=0;j<=(nx+2*npd-1);j++)
              {
                ee[i][j]=ee[nz+npd-1][j];
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


void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,int npd,int vnx,int vnz,int ds_sxd,int is)
{
      
            
                         int ivstart,ivend;
			 int i,ix,iz;
                         is=1;
                         ivstart=1+(is-1)*ds_sxd;
			 ivend=nx+(is-1)*ds_sxd;

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
				  vp[iz][ix-ivstart+1]=vp0[ix][iz];
                                  rho[iz][ix-ivstart+1]=rho0[ix][iz];
				  
				 }
			 }

}


void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
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



void freea(float **a,int m)
{
	int i;

	for(i=0;i<m;i++)
	{
		free(a[i]);
	}
	free(a);
}


