//a#########################################################
//a##         2D Elastic TTI Medium Forward   
//a##  Ps : P0 + sv wave and get rid of sv        
//a##       GPU(CUDA)  
//a##
//a##                        programming by Rong Tao
//a##                                     Rong Tao 
//a##                                   2017.4.8
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

#define pi 3.141592653

#define mm 4

__constant__ float c1[2]={0.0,0.5};
__constant__ float c4[5]={0.0,0.8,-0.2,0.038095,-0.0035714};
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
       float x_,x2_,tdelay,ts,source=0.0,fs; 
  
       tdelay=1.0/favg;
       ts=t-tdelay;
       fs=xsn+(is-1)*ds;

	if(wtype==1)//ricker wavelet
	{
          x_=favg*ts;
          x2_=x_*x_;
          source=(1-2*pi*pi*(x2_))*exp(-(pi*pi*x2_));
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
                           float *delta,float *epsilon,float *theta,int fs,int ds,int zs,int is,bool SV)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

	int im,ix,iz,imx,imz,rx,rz,R=15,r=4;
	float dttxx,dttzz,dttxz,Pxx,Qzz,Pzz,Qxx,Pxz,Qxz,ee,dd,c11,c13,c33,c44,cc,ss,s2;

        ix=id/nnz;
        iz=id%nnz;

               dttxx=dt*dt/(dx*dx);
		 dttzz=dt*dt/(dz*dz);
		 dttxz=dt*dt/(dx*dz);

               if(id>=0&&id<nnx*nnz)
                 {
/************************i****************************************/
/************************iso circle start*************************/
                   rx=ix-(fs+(is-1)*ds+npd);
                   rz=iz-(zs+npd);
                   if(SV){
                       if((rx*rx+rz*rz)<=R*R){
                           if((rx*rx+rz*rz)<=r*r){
                               ee = 0.0;
                               dd = 0.0;
                           }else{
                               ee = 0.5*(1-cos(pi*((sqrtf(rx*rx+rz*rz)-r)*4.0/(R*3.0-1))))*epsilon[id];
                               dd = 0.5*(1-cos(pi*((sqrtf(rx*rx+rz*rz)-r)*4.0/(R*3.0-1))))*delta[id]; 
                              }
                       }else{
                          ee=epsilon[id];
                          dd=delta[id];
                          }
                   }else{
                      ee=epsilon[id];
                      dd=delta[id];
                     }
/************************ iso circle end *************************/
/************************i****************************************/

                   if(ix>=mm&&ix<(nnx-mm-1)&&iz>=mm&&iz<(nnz-mm-1))
                     {

/************************ Pxx Pzz Qzz Qzx ********************************/
                      Pxx=stencil[0]*P1[id];
                      Pzz=stencil[0]*P1[id];
                      Qzz=stencil[0]*Q1[id];
                      Qxx=stencil[0]*Q1[id];

	             for(im=1;im<=mm;im++)
                       {
                        Pxx+=stencil[im]*(P1[id+im*nnz]+P1[id-im*nnz]);
                        Pzz+=stencil[im]*(P1[id+im]    +P1[id-im]);
                        Qzz+=stencil[im]*(Q1[id+im]    +Q1[id-im]);
                        Qxx+=stencil[im]*(Q1[id+im*nnz]+Q1[id-im*nnz]);
                       }
/************************ Pxz Qxz ********************************/
                      Pxz=0.0;
                      Qxz=0.0; 
	             for(imz=0;imz<=mm;imz++)
                      {
	                 for(imx=0;imx<=mm;imx++)
                           {
                            Pxz+=c4[imx]*c4[imz]*(P1[id+imx*nnz+imz]+P1[id-imx*nnz-imz]-P1[id-imx*nnz+imz]-P1[id+imx*nnz-imz]);
                            Qxz+=c4[imx]*c4[imz]*(Q1[id+imx*nnz+imz]+Q1[id-imx*nnz-imz]-Q1[id-imx*nnz+imz]-Q1[id+imx*nnz-imz]);
                           }

                      }
/****a*****************************************************************/

                       Pxx*=dttxx;
                       Pzz*=dttzz;
                       Pxz*=dttxz;
                       Qzz*=dttzz;
                       Qxx*=dttxx;
                       Qxz*=dttxz;

/****a*****************************************************************/
                     cc=cos(theta[id])*cos(theta[id]);
                     ss=sin(theta[id])*sin(theta[id]);
                     s2=sin(theta[id]*2.0);
/****a*****************************************************************/
                     c11=vp[id]*vp[id]*(1+2*ee);
                     c13=vp[id]*vp[id]*(1+2*dd);
                     //c13=sqrtf((2*dd*vp[id]*vp[id]+(vp[id]*vp[id]-vs[id]*vs[id]))*(vp[id]*vp[id]-vs[id]*vs[id]))-vs[id]*vs[id];
                     c33=vp[id]*vp[id];
                     c44=vs[id]*vs[id];
/****a*****************************************************************/

                       P0[id] = 2.0*P1[id] - P0[id]  + c11*(cc*Pxx + ss*Pzz - s2*Pxz)
                                                     + c33*(ss*Qxx + cc*Qzz + s2*Qxz)
                                                     + c44*(ss*Pxx + cc*Pzz + s2*Pxz - ss*Qxx - cc*Qzz - s2*Qxz);

                       Q0[id] = 2.0*Q1[id] - Q0[id]  + c13*(cc*Pxx + ss*Pzz - s2*Pxz)
                                                     + c33*(ss*Qxx + cc*Qzz + s2*Qxz)
                                                     - c44*(cc*Pxx + ss*Pzz - s2*Pxz - cc*Qxx - ss*Qzz + s2*Qxz);
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
void read_file(char FN1[],char FN2[],char FN3[],char FN4[],int nx,int nz,int nnx,int nnz,float *vv,float *epsilon,float *delta,float *theta,
               float *vs,int npd)
{
		 int i,j,id;
		
		 FILE *fp1,*fp2,*fp3,*fp4;
		 if((fp1=fopen(FN1,"rb"))==NULL)printf("error open <%s>!\n",FN1);
		 if((fp2=fopen(FN2,"rb"))==NULL)printf("error open <%s>!\n",FN2);
		 if((fp3=fopen(FN3,"rb"))==NULL)printf("error open <%s>!\n",FN3);
		 if((fp4=fopen(FN4,"rb"))==NULL)printf("error open <%s>!\n",FN4);
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
                            id=i*nnz+j;
				 fread(&vv[id],4L,1,fp1);//vv[id]=2000;
                                                            
				 fread(&epsilon[id],4L,1,fp2);//epsilon[id]=0.24;
				 fread(&delta[id],4L,1,fp3);//delta[id]=0.1;
				 fread(&theta[id],4L,1,fp4);theta[id]*=pi/180.0;//theta[id]=45.0;

                                                              vs[id]=sqrt(vv[id]*vv[id]*(epsilon[id]-delta[id])/999999);
                                                              //vs[id]=vv[id]/2.0;
                                                             // vs[id]=0.0;
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
		 fclose(fp4);
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

	float *v, *e, *d, *s, *th;
	float *vp, *epsilon, *delta, *theta, *vs;
	float *P0, *Q0, *P1, *Q1, *shot_Dev, *shot_Hos, *buffer;

       clock_t start, end;
/*************wavelet\boundary**************/
          wtype=1;npd=50;
/********** dat document ***********/
          char FN1[250]={"2layer_vel_1001_601.dat"};
          char FN2[250]={"2layer_epsilon_1001_601.dat"};
          char FN3[250]={"2layer_delta_1001_601.dat"};
          char FN4[250]={"2layer_theta_1001_601.dat"};
	   char FN5[250]={"thrust_shot_stable.dat"};
	   char FN6[250]={"thrust_snap_stable.dat"};
/********aaa************/  
	 FILE *fpsnap, *fpshot;
        fpshot=fopen(FN5,"wb");
        fpsnap=fopen(FN6,"wb");

 
/********* parameters *************/

          nx=1001;              
	   nz=601;         favg=40;     pfac=10.0;

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=2501;    
          dt=0.0005;
     
          ns=1;       
          fs=nx/2;      
          ds=0;
          zs=1;     
/*************v***************/ 
          nnx=nx+2*npd;
          nnz=nz+2*npd;
/************a*************/


    	 v=(float*)malloc(nnz*nnx*sizeof(float));
    	 e=(float*)malloc(nnz*nnx*sizeof(float));
    	 d=(float*)malloc(nnz*nnx*sizeof(float));
    	 s=(float*)malloc(nnz*nnx*sizeof(float));
    	 th=(float*)malloc(nnz*nnx*sizeof(float));
    	 shot_Hos=(float*)malloc(nt*nx*sizeof(float));
        read_file(FN1,FN2,FN3,FN4,nx,nz,nnx,nnz,v,e,d,th,s,npd);
/****************************/
        pad_vv(nx,nz,nnx,nnz,npd,e);
        pad_vv(nx,nz,nnx,nnz,npd,d);
        pad_vv(nx,nz,nnx,nnz,npd,v);
        pad_vv(nx,nz,nnx,nnz,npd,th); 
        pad_vv(nx,nz,nnx,nnz,npd,s); 

        cudaSetDevice(0);// initialize device, default device=0;
	 check_gpu_error("Failed to initialize device!");

/****************************/
        cudaMalloc(&vp, nnz*nnx*sizeof(float));
        cudaMalloc(&theta, nnz*nnx*sizeof(float));
        cudaMalloc(&vs, nnz*nnx*sizeof(float));
        cudaMalloc(&epsilon, nnz*nnx*sizeof(float));
        cudaMalloc(&delta, nnz*nnx*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(vs, s, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilon, e, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(delta, d, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(theta, th, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
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
        VTI_FD<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,P0,Q0,P1,Q1,vp,vs,npd,delta,epsilon,theta,fs,ds,zs,is,false);
        buffer=P0;P0=P1;P1=buffer; buffer=Q0;Q0=Q1;Q1=buffer;
        absorb_bndr<<<(nnx*nnz+511)/512, 512>>>(P0,P1,Q0,Q1,nx,nz,nnz,npd,-0.25);
        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npd, it, nt, P0, shot_Dev);


           if((is==1)&&(it%300==0)&&it!=0)
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
       cudaFree(theta);
       cudaFree(vs);
       cudaFree(epsilon);
       cudaFree(delta);
/***************host free*****************/
	free(v);	free(e);	free(d);    free(s);   free(th);
       free(shot_Hos);
}

