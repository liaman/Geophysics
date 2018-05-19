//a#########################################################
//a##         2D Elastic VTI Medium Forward   
//a##  Ps : P0 + sv wave and get rid of sv        
//a##       GPU(CUDA)  
//a##
//a##/*a***************************
//a##Function for VTI medium modeling,2017.2.13
//a##
//a## Ps:  the function of modeling following:
//a##      
//a##          du/dt=1/rho*dp/dx , 
//a##          dw/dt=1/rho*dq/dz ,  
//a##          dp/dt=rho*vpx^2*du/dx+rho*vp*vpn*dw/dz ,
//a##          dq/dt=rho*vp*vpn*du/dx+rho*vp^2*dw/dz ,
//a##                     vpx^2=vp^2*(1+2*epsilon);
//a##                     vpn^2=vp^2*(1+2*delta);
//a##*********a*******************/
//a##                        programming by Rong Tao
//a##                                     Rong Tao 
//a##                                  2017.2.14
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

#define pi 3.141592653

#define mm 4

__constant__ float c[mm]={1.196289,-0.0797526,0.009570313,-0.0006975447};
__constant__ float stencil[5]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};

//a################################################################################
void check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	printf("Cuda error: %s: %s\n", msg, cudaGetErrorString (err)); 
	exit(0);   
    }
}
//a################################################################################
__global__ void add_source(float pfac,float xsn,float zsn,int nx,int nz,int nnx,int nnz,float dt,float t,
                        float favg,int wtype,int npd,int is,int ds,float *P,float *Q)
/*< generate ricker wavelet with time deley >*/
{
       int ixs,izs;
       float x_,xx_,tdelay,ts,source=0.0,fs; 
  
       tdelay=1.0/favg;
       ts=t-tdelay;
       fs=xsn+(is-1)*ds;

	if(wtype==1)//ricker wavelet
	{
          x_=favg*ts;
          xx_=x_*x_;
          source=(1-2*pi*pi*(xx_))*exp(-(pi*pi*xx_));
	}else if(wtype==2){//derivative of gaussian
          x_=(-4)*favg*favg*pi*pi/log(0.1);
          source=(-2)*pi*pi*ts*exp(-x_*ts*ts);
        }else if(wtype==3){//derivative of gaussian
          x_=(-1)*favg*favg*pi*pi/log(0.1);
          source=exp(-x_*ts*ts);
        }

       if(t<=2*tdelay)
       {         
	     ixs = (int)(fs+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;
            P[ixs*nnz+izs]+=pfac*source;
            Q[ixs*nnz+izs]+=pfac*source;
       }
}
/*******************func*********************/
/*******************func***********************/
__global__ void VTI_FD(int nx,int nz,int nnx,int nnz,float dt,float dx,float dz,
                           float *P0,float *Q0,float *P1,float *Q1,float *vp,float *vs,int npd,
                           float *delta,float *epsilon)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

	int im,ix,iz;
	float dtx,dtz, xx,zz,xz,zx,ee,dd,c11,c13,c33,c44;

        ix=id/nnz;
        iz=id%nnz;

               dtx=dt*dt/(dx*dx);
		 dtz=dt*dt/(dz*dz);
               if(id>=0&&id<nnx*nnz)
                 {
                      ee=epsilon[id];
                      dd=delta[id];

                   if(ix>=mm&&ix<(nnx-mm-1)&&iz>=mm&&iz<(nnz-mm-1))
                     {


                      xx=stencil[0]*P1[id];
                      xz=stencil[0]*P1[id];
                      zz=stencil[0]*Q1[id];
                      zx=stencil[0]*Q1[id];

	             for(im=1;im<=mm;im++)
                       {
                        xx+=stencil[im]*(P1[id+im*nnz]+P1[id-im*nnz]);
                        xz+=stencil[im]*(P1[id+im]    +P1[id-im]);
                        zz+=stencil[im]*(Q1[id+im]    +Q1[id-im]);
                        zx+=stencil[im]*(Q1[id+im*nnz]+Q1[id-im*nnz]);

                       }
                       xx*=dtx;
                       xz*=dtz;
                       zz*=dtz;
                       zx*=dtx;

                     c11=vp[id]*vp[id]*(1+2*ee);
                     c13=vp[id]*vp[id]*(1+2*dd);
                     //c13=sqrtf((2*dd*vp[id]*vp[id]+(vp[id]*vp[id]-vs[id]*vs[id]))*(vp[id]*vp[id]-vs[id]*vs[id]))-vs[id]*vs[id];
                     c33=vp[id]*vp[id];
                     c44=vs[id]*vs[id];

                       P0[id] = 2.0*P1[id] - P0[id]  + xx*c11 + zz*c33 + c44*(xz - zz);

                       Q0[id] = 2.0*Q1[id] - Q0[id]  + xx*c13 + zz*c33 - c44*(xx - zx);
                   }
                 }
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
       //if(ee[id]==0){printf("ee[%d][%d]==0.0\n",ix,iz);exit(0);}
     }
}
/*************func*******************/
void read_file(char FN1[],char FN2[],char FN3[],int nx,int nz,int nnx,int nnz,float *vv,float *epsilon,float *delta,float *vs,int npd)
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
				 fread(&vv[id],4L,1,fp1);//vv[id]=2000;
                                                            
				 fread(&epsilon[id],4L,1,fp2);//epsilon[id]=0.24;
				 fread(&delta[id],4L,1,fp3);//delta[id]=0.1;

                                                             // vs[id]=sqrt(vv[id]*vv[id]*(epsilon[id]-delta[id])/0.75);
                                                             // vs[id]=vv[id]/2.0;
                                                              vs[id]=0.0;
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
}

/*************func*******************/
__global__ void shot_record(int nnx, int nnz, int nx, int nz, int npd, int it, int nt, float *P0, float *shot)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<nx)
            {
               shot[it+nt*id]=P0[npd+nnz*(id+npd)];
            }       
}
/*************func**************/    
__global__ void mute_directwave(int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float *vp,float *epsilon,float *shot,int tt)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

    int mu_t,mu_nt;
    float mu_x,mu_z,mu_t0;

    int ix=id/nt;
    int it=id%nt;

   if(id<nx*nt)
   {
        mu_x=dx*abs(ix-fs-(is-1)*ds);
        mu_z=dz*zs;
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_z,2))/(vp[1]*sqrtf(1+2*epsilon[1]));
        mu_t=(int)(2.0/(dt*favg));
        mu_nt=(int)(mu_t0/dt)+mu_t+tt;

           if(it<mu_nt)
              shot[id]=0.0;
   }
}/************************************func***************************************/      
__global__ void absorb_bndr(float *P0,float *P1,float *Q0,float *Q1,int nx,int nz,int nnz,int npd,float qp) 
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix,iz;

        ix=id/nnz;
        iz=id%nnz;

   if(id<(nx+2*npd)*nnz-1)
   {
         if(ix<npd){
               P0[id]*=( qp*pow((npd-ix)/(1.0*npd),2) + 1 );
               P1[id]*=( qp*pow((npd-ix)/(1.0*npd),2) + 1 );
               Q0[id]*=( qp*pow((npd-ix)/(1.0*npd),2) + 1 );
               Q1[id]*=( qp*pow((npd-ix)/(1.0*npd),2) + 1 );
         }else if(ix>=nx+npd){
               P0[id]*=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 );
               P1[id]*=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 );
               Q0[id]*=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 );
               Q1[id]*=( qp*pow((ix-npd-nx)/(1.0*npd),2) + 1 );
         }if(iz<npd){
               P0[id]*=( qp*pow((npd-iz)/(1.0*npd),2) + 1 );
               P1[id]*=( qp*pow((npd-iz)/(1.0*npd),2) + 1 );
               Q0[id]*=( qp*pow((npd-iz)/(1.0*npd),2) + 1 );
               Q1[id]*=( qp*pow((npd-iz)/(1.0*npd),2) + 1 );
         }else if(iz>=nz+npd){
               P0[id]*=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 );
               P1[id]*=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 );
               Q0[id]*=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 );
               Q1[id]*=( qp*pow((iz-npd-nz)/(1.0*npd),2) + 1 );
         }
    }
}   
//a########################################################################
int main(int argc,char *argv[])
{
	int is, it, nx, nz, nnx, nnz, nt, wtype;
	int ns, ds, fs, zs, npd;
	float dx, dz, dt, t, pfac, favg;

	float *v, *e, *d, *s;
	float *vp, *epsilon, *delta, *vs;
	float *P0, *Q0, *P1, *Q1, *shot_Dev, *shot_Hos, *buffer;

       clock_t start, end;
/*************wavelet\boundary**************/
          wtype=1;npd=50;
/********** dat document ***********/
          char FN1[250]={"BP_vel.dat"};
          char FN2[250]={"BP_epsilon.dat"};
          char FN3[250]={"BP_delta.dat"};
	   char FN4[250]={"BP_shot_test.dat"};
	   char FN5[250]={"BP_snap_test.dat"};

/********aaa************/  
	 FILE *fpsnap, *fpshot;
        fpshot=fopen(FN4,"wb");
        fpsnap=fopen(FN5,"wb");

 
/********* parameters *************/

          nx=801;              
	   nz=601;         favg=30;     pfac=10.0;

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=1501;    
          dt=0.0005;
     
          ns=1;       
          fs=401;      
          ds=0;
          zs=201;     
/*************v***************/ 
          nnx=nx+2*npd;
          nnz=nz+2*npd;
/************a*************/


    	 v=(float*)malloc(nnz*nnx*sizeof(float));
    	 e=(float*)malloc(nnz*nnx*sizeof(float));
    	 d=(float*)malloc(nnz*nnx*sizeof(float));
    	 s=(float*)malloc(nnz*nnx*sizeof(float));
    	 shot_Hos=(float*)malloc(nt*nx*sizeof(float));
        read_file(FN1,FN2,FN3,nx,nz,nnx,nnz,v,e,d,s,npd);
/****************************/
        pad_vv(nx,nz,nnx,nnz,npd,e);
        pad_vv(nx,nz,nnx,nnz,npd,d);
        pad_vv(nx,nz,nnx,nnz,npd,v); 
        pad_vv(nx,nz,nnx,nnz,npd,s); 

        cudaSetDevice(0);// initialize device, default device=0;
	 check_gpu_error("Failed to initialize device!");

/****************************/
        cudaMalloc(&vp, nnz*nnx*sizeof(float));
        cudaMalloc(&vs, nnz*nnx*sizeof(float));
        cudaMalloc(&epsilon, nnz*nnx*sizeof(float));
        cudaMalloc(&delta, nnz*nnx*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(vs, s, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilon, e, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(delta, d, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
/****************************/

        cudaMalloc(&P0, nnz*nnx*sizeof(float));     cudaMalloc(&Q0, nnz*nnx*sizeof(float));
        cudaMalloc(&P1, nnz*nnx*sizeof(float));     cudaMalloc(&Q1, nnz*nnx*sizeof(float));

        cudaMalloc(&shot_Dev, nx*nt*sizeof(float));
/******************************/
	 check_gpu_error("Failed to allocate memory for variables!");
        printf("--------------------------------------------------------\n");
        printf("---   \n");   
        start = clock();                                  
/**********IS Loop start*******/
   for(is=1;is<=ns;is++)	
    {     
         printf("---   IS=%3d  \n",is);

     cudaMemset(P0, 0, nnz*nnx*sizeof(float));      cudaMemset(Q0, 0, nnz*nnx*sizeof(float));
     cudaMemset(P1, 0, nnz*nnx*sizeof(float));      cudaMemset(Q1, 0, nnz*nnx*sizeof(float));

     cudaMemset(shot_Dev, 0, nt*nx*sizeof(float));

     for(it=0,t=dt;it<nt;it++,t+=dt)
     { 
      if(it%100==0&&is==1)printf("---   is===%d   it===%d\n",is,it);

	 add_source<<<1,1>>>(pfac,fs,zs,nx,nz,nnx,nnz,dt,t,favg,wtype,npd,is,ds,P0,Q0);
        VTI_FD<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,P0,Q0,P1,Q1,vp,vs,npd,delta,epsilon);
        buffer=P0;P0=P1;P1=buffer; buffer=Q0;Q0=Q1;Q1=buffer;
        absorb_bndr<<<(nnx*nnz+511)/512, 512>>>(P0,P1,Q0,Q1,nx,nz,nnz,npd,-0.25);
        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npd, it, nt, P0, shot_Dev);


           if((is==1)&&(it%100==0)&&it!=0)
            {
	       cudaMemcpy(e, P0, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
             // fseek(fpsnap,(int)(it/1200)*(nnx)*(nnz)*4L,0);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }
     }//it loop end
    //  mute_directwave<<<(nx*nt+511)/512, 512>>>(nx,nt,dt,favg,dx,dz,fs,ds,zs,is,vp,epsilon,shot_Dev,30);
      cudaMemcpy(shot_Hos, shot_Dev, nt*nx*sizeof(float), cudaMemcpyDeviceToHost);
      fseek(fpshot,(is-1)*nt*nx*sizeof(float),0);
      fwrite(shot_Hos,sizeof(float),nt*nx,fpshot);

    }//is loop end
    end = clock();
/*********IS Loop end*********/ 		     
   printf("---   The forward is over    \n"); 
   printf("---   Complete!!!!!!!!! \n");  
   printf("total %d shots: %f (s)\n", ns, ((float)(end-start))/CLOCKS_PER_SEC);



/***********close************/ 
          fclose(fpsnap);   fclose(fpshot);
/***********free*************/ 

       cudaFree(P0);            cudaFree(Q0);
       cudaFree(P1);            cudaFree(Q1);

       cudaFree(shot_Dev);

       cudaFree(vp);
       cudaFree(vs);
       cudaFree(epsilon);
       cudaFree(delta);
/***************host free*****************/
	free(v);	free(e);	free(d);    free(s);
       free(shot_Hos);
}

