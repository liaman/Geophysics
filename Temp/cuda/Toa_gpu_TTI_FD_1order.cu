//a#########################################################
//a##         2D Acoustic TTI Medium Forward   
//a##  Ps : P + sv wave and get rid of sv        
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
//a##                                   2017.2.21
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

#define pi 3.141592653

#define mm 4

__device__ float d0;

__constant__ float c[mm]={1.196289,-0.0797526,0.009570313,-0.0006975447};

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
                        float favg,int wtype,int npml,int is,int ds,float *P,float *Q)
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
	     ixs = (int)(fs+0.5)+npml-1;
            izs = (int)(zsn+0.5)+npml-1;
            P[ixs*nnz+izs]+=pfac*source;
            Q[ixs*nnz+izs]+=pfac*source;
       }
}
/*******************func*********************/
__global__ void update_vel(int nx,int nz,int nnx,int nnz,int npml,float dt,float dx,float dz,
                           float *u0,float *w0,float *u1,float *w1,float *P,float *Q,
                           float *coffx1,float *coffx2,float *coffz1,float *coffz2,float *theta)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

	int ix,iz,im;
	float dtx,dtz,Px,Pz,Qz,Qx;

        ix=id/nnz;
        iz=id%nnz;

		 dtx=dt/dx;
		 dtz=dt/dz;
               if(id>=mm&&id<nnx*nnz-mm)
                 {
                   if(ix>=mm&&ix<(nnx-mm)&&iz>=mm&&iz<(nnz-mm))
                    {
                     Px=0.0;
                     Pz=0.0;
                     Qz=0.0;
                     Qx=0.0;
	             for(im=0;im<mm;im++)
                      {
                        Px+=c[im]*(P[id+(im+1)*nnz]-P[id-im*nnz]);
                        Pz+=c[im]*(P[id+im+1]      -P[id-im]);
                        Qz+=c[im]*(Q[id+im+1]      -Q[id-im]);
                        Qx+=c[im]*(Q[id+(im+1)*nnz]-Q[id-im*nnz]);
                      }
                     u1[id]=coffx2[ix]*coffz2[iz]*u0[id]-coffx1[ix]*coffz1[iz]*(cos(theta[id])*dtx*Px-sin(theta[id])*dtz*Pz);
                     w1[id]=coffx2[ix]*coffz2[iz]*w0[id]-coffx1[ix]*coffz1[iz]*(sin(theta[id])*dtx*Qx+cos(theta[id])*dtz*Qz);
                   }
                 }
}
/*******************func***********************/
__global__ void update_stress(int nx,int nz,int nnx,int nnz,float dt,float dx,float dz,
                           float *u1,float *w1,float *P,float *Q,float *vp,int npml,
                           float *px1,float *px0,float *pz1,float *pz0,float *qx1,float *qx0,float *qz1,float *qz0,
                           float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                           float *delta,float *epsilon,float *theta,int fs,int ds,int zs,int is,bool SV)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

	int im,ix,iz,rx,rz,R=15,r=4;
	float dtx,dtz, ux,uz,wz,wx,ee,dd;

        ix=id/nnz;
        iz=id%nnz;

               dtx=dt/dx;
		 dtz=dt/dz;
               if(id>=mm&&id<nnx*nnz-mm)
                 {
/************************i****************************************/
/************************iso circle start*************************/
                   rx=ix-(fs+(is-1)*ds+npml);
                   rz=iz-(zs+npml);
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
                   if(ix>=mm&&ix<(nnx-mm)&&iz>=mm&&iz<(nnz-mm))
                     {
                     ux=0.0;
                     uz=0.0;
                     wz=0.0;
                     wx=0.0;
	             for(im=0;im<mm;im++)
                       {
                        ux+=c[im]*(u1[id+im*nnz]-u1[id-(im+1)*nnz]);
                        uz+=c[im]*(u1[id+im]    -u1[id-im-1]);
                        wz+=c[im]*(w1[id+im]    -w1[id-im-1]);
                        wx+=c[im]*(w1[id+im*nnz]-w1[id-(im+1)*nnz]);
                       }
                     px1[id]=acoffx2[ix]*acoffz2[iz]*px0[id]
                            -acoffx1[ix]*acoffz1[iz]*vp[id]*vp[id]*(1+2*ee)*     (cos(theta[id])*dtx*ux-sin(theta[id])*dtz*uz);
                     pz1[id]=acoffx2[ix]*acoffz2[iz]*pz0[id]
                            -acoffx1[ix]*acoffz1[iz]*vp[id]*vp[id]*sqrtf(1+2*dd)*(sin(theta[id])*dtx*wx+cos(theta[id])*dtz*wz);
                     qx1[id]=acoffx2[ix]*acoffz2[iz]*qx0[id]
                            -acoffx1[ix]*acoffz1[iz]*vp[id]*vp[id]*sqrtf(1+2*dd)*(cos(theta[id])*dtx*ux-sin(theta[id])*dtz*uz);
                     qz1[id]=acoffx2[ix]*acoffz2[iz]*qz0[id]
                            -acoffx1[ix]*acoffz1[iz]*vp[id]*vp[id]*              (sin(theta[id])*dtx*wx+cos(theta[id])*dtz*wz);

                     P[id]=px1[id]+pz1[id];
                     Q[id]=qx1[id]+qz1[id];
                   }
                 }
}                      
/********************func**********************/
__global__ void get_d0(float dx,float dz,int nnx,int nnz,int npml,float *vp)
{
   d0=10.0*vp[nnx*nnz/2]*log(100000.0)/(2.0*npml*((dx+dz)/2.0));
}
/*************func*******************/
void pad_vv(int nx,int nz,int nnx,int nnz,int npml,float *ee)
{
     int ix,iz,id;
 
    for(id=0;id<nnx*nnz;id++)
     {
       ix=id/nnz;
       iz=id%nnz;
       if(ix<npml){
           ee[id]=ee[npml*nnz+iz];  //left
       }else if(ix>=nnx-npml){
           ee[id]=ee[(nnx-npml-1)*nnz+iz];//right
       }
     }
    for(id=0;id<nnx*nnz;id++)
     {
       ix=id/nnz;
       iz=id%nnz;
       if(iz<npml){
           ee[id]=ee[ix*nnz+npml];//up
       }else if(iz>=nnz-npml){
           ee[id]=ee[ix*nnz+nnz-npml-1];//down
       }
      // if(ee[id]==0){printf("ee[%d][%d]==0.0\n",ix,iz);exit(0);}
     }
}
/*************func*******************/
void read_file(char FN1[],char FN2[],char FN3[],char FN4[],int nx,int nz,int nnx,int nnz,
               float *vv,float *epsilon,float *delta,float *theta,int npml)
{
		 int i,j,id;
		
		 FILE *fp1,*fp2,*fp3,*fp4;
		 if((fp1=fopen(FN1,"rb"))==NULL)printf("error open <%s>!\n",FN1);
		 if((fp2=fopen(FN2,"rb"))==NULL)printf("error open <%s>!\n",FN2);
		 if((fp3=fopen(FN3,"rb"))==NULL)printf("error open <%s>!\n",FN3);
		 if((fp4=fopen(FN4,"rb"))==NULL)printf("error open <%s>!\n",FN4);
		 for(i=npml;i<nx+npml;i++)
		 {
			 for(j=npml;j<nz+npml;j++)
			 {
                            id=i*nnz+j;
				 fread(&vv[id],4L,1,fp1);//vv[id]=vv[npml*nnz+npml];
				 fread(&epsilon[id],4L,1,fp2);//epsilon[id]=0;
				 fread(&delta[id],4L,1,fp3);//delta[id]=0;
				 fread(&theta[id],4L,1,fp4);theta[id]*=pi/180.0;//theta[id]=0;
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
		 fclose(fp4);
}
/*************func*******************/
__global__ void initial_coffe(float dt,int nn,float *coff1,float *coff2,float *acoff1,float *acoff2,int npml)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<nn+2*npml)
            {
		 if(id<npml)
		 {   
			 coff1[id]=1.0/(1.0+(dt*d0*pow((npml-0.5-id)/npml,2.0))/2.0);
			 coff2[id]=coff1[id]*(1.0-(dt*d0*pow((npml-0.5-id)/npml,2.0))/2.0);

			 acoff1[id]=1.0/(1.0+(dt*d0*pow(((npml-id)*1.0)/npml,2.0))/2.0);
			 acoff2[id]=acoff1[id]*(1.0-(dt*d0*pow(((npml-id)*1.0)/npml,2.0))/2.0);

		 }else if(id>=npml&&id<npml+nn){

			 coff1[id]=1.0;
			 coff2[id]=1.0;

			 acoff1[id]=1.0;
			 acoff2[id]=1.0;

		 }else{

			 coff1[id]=1.0/(1.0+(dt*d0*pow((0.5+id-nn-npml)/npml,2.0))/2.0);
			 coff2[id]=coff1[id]*(1.0-(dt*d0*pow((0.5+id-nn-npml)/npml,2.0))/2.0);

			 acoff1[id]=1.0/(1.0+(dt*d0*pow(((id-nn-npml)*1.0)/npml,2.0))/2.0);
			 acoff2[id]=acoff1[id]*(1.0-(dt*d0*pow(((id-nn-npml)*1.0)/npml,2.0))/2.0);
		 }	
            }       
}
/*************func*******************/
__global__ void shot_record(int nnx, int nnz, int nx, int nz, int npml, int it, int nt, float *P, float *shot)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<nx)
            {
               shot[it+nt*id]=P[npml+nnz*(id+npml)];
            }       
}
/*************func**************/    
/*************func**************/    
__global__ void mute_directwave(int nx,int nt,float dt,float favg,
                     float dx,float dz,int fs,int ds,int zs,int is,
                     float *vp,float *epsilu,float *theta,float *shot,int tt)
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
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_z,2))/(vp[1]*(  1+(sqrtf(1+2*epsilu[1])-1) * cos(theta[1])  )   );
        mu_t=(int)(2.0/(dt*favg));
        mu_nt=(int)(mu_t0/dt)+mu_t+tt;

           if((it>(int)(mu_t0/dt)-tt)&&(it<mu_nt))
              shot[id]=0.0;
   }
}
//a########################################################################
int main(int argc,char *argv[])
{
	int is, it, nx, nz, nnx, nnz, nt, wtype;
	int ns, ds, fs, zs, npml;
	float dx, dz, dt, t, pfac, favg;

	float *v, *e, *d, *b;
	float *vp, *epsilon, *delta, *theta;
	float *u0, *u1, *px0, *qx0, *px1, *qx1;
       float *w0, *w1, *pz0, *qz0, *pz1, *qz1;
	float *P, *Q, *shot_Dev, *shot_Hos;

       float *coffx1,*coffx2,*coffz1,*coffz2,*acoffx1,*acoffx2,*acoffz1,*acoffz2;

       clock_t start, end;
/*************wavelet\boundary**************/
          wtype=1;npml=50;
/********** dat document ***********/
          char FN1[250]={"thrust_vel_711_300.dat"};
          char FN2[250]={"thrust_epsilon_711_300.dat"};
          char FN3[250]={"thrust_delta_711_300.dat"};
          char FN4[250]={"thrust_theta_711_300.dat"};
	   char FN5[250]={"thrust_shot_unstable.dat"};
	   char FN6[250]={"thrust_snap_unstable.dat"};

/********aaa************/  
	 FILE *fpsnap, *fpshot;
        fpshot=fopen(FN5,"wb");
        fpsnap=fopen(FN6,"wb");

 
/********* parameters *************/


          nx=600;              
	   nz=300;         favg=40;     pfac=10.0;

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=801;    
          dt=0.0005;
     
          ns=1;       
          fs=200;      
          ds=0;
          zs=140;  
  
/*************v***************/ 
          nnx=nx+2*npml;
          nnz=nz+2*npml;
/************a*************/


    	 v=(float*)malloc(nnz*nnx*sizeof(float));
    	 e=(float*)malloc(nnz*nnx*sizeof(float));
    	 d=(float*)malloc(nnz*nnx*sizeof(float));
    	 b=(float*)malloc(nnz*nnx*sizeof(float));
    	 shot_Hos=(float*)malloc(nt*nx*sizeof(float));
        read_file(FN1,FN2,FN3,FN4,nx,nz,nnx,nnz,v,e,d,b,npml);
/****************************/
        pad_vv(nx,nz,nnx,nnz,npml,e);
        pad_vv(nx,nz,nnx,nnz,npml,d);
        pad_vv(nx,nz,nnx,nnz,npml,v); 
        pad_vv(nx,nz,nnx,nnz,npml,b); 

        cudaSetDevice(0);// initialize device, default device=0;
	 check_gpu_error("Failed to initialize device!");

/****************************/
        cudaMalloc(&vp, nnz*nnx*sizeof(float));
        cudaMalloc(&epsilon, nnz*nnx*sizeof(float));
        cudaMalloc(&delta, nnz*nnx*sizeof(float));
        cudaMalloc(&theta, nnz*nnx*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilon, e, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(delta, d, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(theta, b, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
/****************************/
        cudaMalloc(&u0, nnz*nnx*sizeof(float));    cudaMalloc(&u1, nnz*nnx*sizeof(float));
        cudaMalloc(&w0, nnz*nnx*sizeof(float));    cudaMalloc(&w1, nnz*nnx*sizeof(float));

        cudaMalloc(&P, nnz*nnx*sizeof(float));     cudaMalloc(&Q, nnz*nnx*sizeof(float));

        cudaMalloc(&px0, nnz*nnx*sizeof(float));   cudaMalloc(&px1, nnz*nnx*sizeof(float));
        cudaMalloc(&pz0, nnz*nnx*sizeof(float));   cudaMalloc(&pz1, nnz*nnx*sizeof(float));
        cudaMalloc(&qx0, nnz*nnx*sizeof(float));   cudaMalloc(&qx1, nnz*nnx*sizeof(float));
        cudaMalloc(&qz0, nnz*nnx*sizeof(float));   cudaMalloc(&qz1, nnz*nnx*sizeof(float));

        cudaMalloc(&coffx1, nnx*sizeof(float));     cudaMalloc(&coffx2, nnx*sizeof(float));
        cudaMalloc(&coffz1, nnz*sizeof(float));     cudaMalloc(&coffz2, nnz*sizeof(float));
        cudaMalloc(&acoffx1, nnx*sizeof(float));    cudaMalloc(&acoffx2, nnx*sizeof(float));
        cudaMalloc(&acoffz1, nnz*sizeof(float));    cudaMalloc(&acoffz2, nnz*sizeof(float));

        cudaMalloc(&shot_Dev, nx*nt*sizeof(float));
/******************************/
	 check_gpu_error("Failed to allocate memory for variables!");

        get_d0<<<1, 1>>>(dx, dz, nnx, nnz, npml, vp);
        initial_coffe<<<(nnx+511)/512, 512>>>(dt,nx,coffx1,coffx2,acoffx1,acoffx2,npml);
        initial_coffe<<<(nnz+511)/512, 512>>>(dt,nz,coffz1,coffz2,acoffz1,acoffz2,npml);



        printf("--------------------------------------------------------\n");
        printf("---   \n");   
        start = clock();                                  
/**********IS Loop start*******/
   for(is=1;is<=ns;is++)	
    {     
         printf("---   IS=%3d  \n",is);

     cudaMemset(u0, 0, nnz*nnx*sizeof(float));     cudaMemset(u1, 0, nnz*nnx*sizeof(float));
     cudaMemset(w0, 0, nnz*nnx*sizeof(float));     cudaMemset(w1, 0, nnz*nnx*sizeof(float));

     cudaMemset(P, 0, nnz*nnx*sizeof(float));      cudaMemset(Q, 0, nnz*nnx*sizeof(float));

     cudaMemset(px0, 0, nnz*nnx*sizeof(float));    cudaMemset(px1, 0, nnz*nnx*sizeof(float));
     cudaMemset(pz0, 0, nnz*nnx*sizeof(float));    cudaMemset(pz1, 0, nnz*nnx*sizeof(float));
     cudaMemset(qx0, 0, nnz*nnx*sizeof(float));    cudaMemset(qx1, 0, nnz*nnx*sizeof(float));
     cudaMemset(qz0, 0, nnz*nnx*sizeof(float));    cudaMemset(qz1, 0, nnz*nnx*sizeof(float));

     cudaMemset(shot_Dev, 0, nt*nx*sizeof(float));

     for(it=0,t=dt;it<nt;it++,t+=dt)
     { 
      if(it%1000==0&&is==1)printf("---   is===%d   it===%d\n",is,it);

	 add_source<<<1,1>>>(pfac,fs,zs,nx,nz,nnx,nnz,dt,t,favg,wtype,npml,is,ds,P,Q);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,u0,w0,u1,w1,P,Q,coffx1,coffx2,coffz1,coffz2,theta);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,u1,w1,P,Q,vp,npml,px1,px0,pz1,pz0,qx1,qx0,qz1,qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,delta,epsilon,theta,fs,ds,zs,is,false);
        u0=u1; w0=w1; px0=px1; pz0=pz1; qx0=qx1; qz0=qz1; 

        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, P, shot_Dev);


           if((is==1)&&(it%100==0)&&it!=0)
            {
	       cudaMemcpy(e, P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              //fseek(fpsnap,(int)(it/100)*(nnx)*(nnz)*4L,0);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }
     }//it loop end
     // mute_directwave<<<(nx*nt+511)/512, 512>>>(nx,nt,dt,favg,dx,dz,fs,ds,zs,is,vp,epsilu,theta,shot_Dev,30);
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
       cudaFree(coffx1);       cudaFree(coffx2);
       cudaFree(coffz1);       cudaFree(coffz2);
       cudaFree(acoffx1);      cudaFree(acoffx2);
       cudaFree(acoffz1);      cudaFree(acoffz2);

       cudaFree(u0);           cudaFree(u1);
       cudaFree(w0);           cudaFree(w1);

       cudaFree(P);            cudaFree(Q);

       cudaFree(px0);          cudaFree(px1);
       cudaFree(pz0);          cudaFree(pz1);
       cudaFree(qx0);          cudaFree(qx1);
       cudaFree(qz0);          cudaFree(qz1);

       cudaFree(shot_Dev);

       cudaFree(vp);
       cudaFree(epsilon);
       cudaFree(delta);
       cudaFree(theta);
/***************host free*****************/
	free(v);	free(e);	free(d);free(b);
       free(shot_Hos);
}

