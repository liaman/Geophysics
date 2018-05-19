#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>
#include<cuda_runtime.h>

#define pi 3.141592653

__device__ float d0;

__constant__ float c[4]={1.196289,-0.0797526,0.009570313,-0.0006975447};

void check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	printf("Cuda error: %s: %s\n", msg, cudaGetErrorString (err)); 
	exit(0);   
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
       if(ee[id]==0){printf("ee[%d][%d]==0.0\n",ix,iz);exit(0);}
     }
}
/*************func*******************/
void read_file(const char FN1[],const char FN2[],const char FN3[],
               int nx,int nz,int nnx,int nnz,float *vv,float *epsilu,float *deta,int npd)
{
		 int i,j,id,vmax=0.0;
		
		 FILE *fp1,*fp2,*fp3;
		 if((fp1=fopen(FN1,"rb"))==NULL)printf("error open <%s>!\n",FN1);
		 if((fp2=fopen(FN2,"rb"))==NULL)printf("error open <%s>!\n",FN2);
		 if((fp3=fopen(FN3,"rb"))==NULL)printf("error open <%s>!\n",FN3);
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
                            id=i*nnz+j;
				 fread(&vv[id],4L,1,fp1);if(vmax<vv[id])vmax=vv[id];
				 fread(&epsilu[id],4L,1,fp2);
				 fread(&deta[id],4L,1,fp3);
			 }
		 }
		 fclose(fp1);printf("vmax=%d\n",vmax);
		 fclose(fp2);
		 fclose(fp3);
}
/********************func**********************/
__global__ void get_d0(float dx_,float dz_,int nnx,int nnz,int npd,float *vp)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
       if(id<1)d0=10.0*vp[nnx*nnz/2]*log(100000.0)/(2.0*npd*((dx_+dz_)/2.0));
}
/*************func*******************/
__global__ void initial_coffe(float dt_,int nn,float *coff1,float *coff2,float *acoff1,float *acoff2,int npd)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<nn+2*npd)
            {
		 if(id<npd)
		 {   
			 coff1[id]=1.0/(1.0+(dt_*d0*pow((npd-0.5-id)/npd,2.0))/2.0);
			 coff2[id]=coff1[id]*(1.0-(dt_*d0*pow((npd-0.5-id)/npd,2.0))/2.0);

			 acoff1[id]=1.0/(1.0+(dt_*d0*pow(((npd-id)*1.0)/npd,2.0))/2.0);
			 acoff2[id]=acoff1[id]*(1.0-(dt_*d0*pow(((npd-id)*1.0)/npd,2.0))/2.0);

		 }else if(id>=npd&&id<npd+nn){

			 coff1[id]=1.0;
			 coff2[id]=1.0;

			 acoff1[id]=1.0;
			 acoff2[id]=1.0;

		 }else{

			 coff1[id]=1.0/(1.0+(dt_*d0*pow((0.5+id-nn-npd)/npd,2.0))/2.0);
			 coff2[id]=coff1[id]*(1.0-(dt_*d0*pow((0.5+id-nn-npd)/npd,2.0))/2.0);

			 acoff1[id]=1.0/(1.0+(dt_*d0*pow(((id-nn-npd)*1.0)/npd,2.0))/2.0);
			 acoff2[id]=acoff1[id]*(1.0-(dt_*d0*pow(((id-nn-npd)*1.0)/npd,2.0))/2.0);
		 }	
            }       
}
/*************func*******************/
__global__ void shot_record(int nnx, int nnz, int nx, int nz, int npd, int it, int nt_, float *P, float *shot)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<nx)
            {
               shot[it+nt_*id]=P[npd+nnz*(id+npd)];
            }       
}
/*************func**************/    
__global__ void mute_directwave(int nx,int nt_,float dt_,float favg_,
                     float dx_,float dz_,int fs,int ds,int zs,int is,
                     float *vp,float *epsilu,float *shot,int tt)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

    int mu_t,mu_nt_;
    float mu_x,mu_z,mu_t0;

    int ix=id/nt_;
    int it=id%nt_;

   if(id<nx*nt_)
   {
        mu_x=dx_*abs(ix-fs-(is-1)*ds);
        mu_z=dz_*zs;
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_z,2))/(vp[1]*sqrtf(1+2*epsilu[1]));
        mu_t=(int)(2.0/(dt_*favg_));
        mu_nt_=(int)(mu_t0/dt_)+mu_t+tt;

           if(it<mu_nt_)
              shot[id]=0.0;
   }
}
//a################################################################################
__global__ void add_source(float pfac,float xsn,float zsn,int nx,int nz,int nnx,int nnz,float dt_,float t,
                        float favg_,int wtype,int npd,int is,int ds,float *P,float *Q)
/*< generate ricker wavelet with time deley >*/
{
       int ixs,izs;
       float x_,xx_,tdelay,ts,source=0.0,fs; 
  
       tdelay=1.0/favg_;
       ts=t-tdelay;
       fs=xsn+(is-1)*ds;

	if(wtype==1)//ricker wavelet
	{
          x_=favg_*ts;
          xx_=x_*x_;
          source=(1-2*pi*pi*(xx_))*exp(-(pi*pi*xx_));
	}else if(wtype==2){//derivative of gaussian
          x_=(-4)*favg_*favg_*pi*pi/log(0.1);
          source=(-2)*pi*pi*ts*exp(-x_*ts*ts);
        }else if(wtype==3){//derivative of gaussian
          x_=(-1)*favg_*favg_*pi*pi/log(0.1);
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
__global__ void update_vel(int nx,int nz,int nnx,int nnz,int npd,int mm,float dt_,float dx_,float dz_,
                           float *u0,float *w0,float *u1,float *w1,float *P,float *Q,
                           float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

	int ix,iz,im;
	float dtx,dtz,xx,zz;

        ix=id/nnz;
        iz=id%nnz;

		 dtx=dt_/dx_;
		 dtz=dt_/dz_;
               if(id>=mm&&id<nnx*nnz-mm)
                 {
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
__global__ void update_stress(int nx,int nz,int nnx,int nnz,float dt_,float dx_,float dz_,
                           float *u1,float *w1,float *P,float *Q,float *vp,int npd,int mm,
                           float *px1,float *px0,float *pz1,float *pz0,float *qx1,float *qx0,float *qz1,float *qz0,
                           float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                           float *deta,float *epsilu,int fs,int ds,int zs,int is,bool SV)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

	int im,ix,iz,rx,rz,R=15,r=4;
	float dtx,dtz, xx,zz,ee,dd;

        ix=id/nnz;
        iz=id%nnz;

               dtx=dt_/dx_;
		 dtz=dt_/dz_;
               if(id>=mm&&id<nnx*nnz-mm)
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
                               ee = 0.5*(1-cos(pi*((sqrtf(rx*rx+rz*rz)-r)*4.0/(R*3.0-1))))*epsilu[id];
                               dd = 0.5*(1-cos(pi*((sqrtf(rx*rx+rz*rz)-r)*4.0/(R*3.0-1))))*deta[id]; 
                              }
                       }else{
                          ee=epsilu[id];
                          dd=deta[id];
                          }
                   }else{
                      ee=epsilu[id];
                      dd=deta[id];
                     }
/************************ iso circle end *************************/
/************************i****************************************/
                   if(ix>=mm&&ix<(nnx-mm)&&iz>=mm&&iz<(nnz-mm))
                     {
                     xx=0.0;
                     zz=0.0;
	             for(im=0;im<mm;im++)
                       {
                        xx+=c[im]*(u1[id+im*nnz]-u1[id-(im+1)*nnz]);
                        zz+=c[im]*(w1[id+im]    -w1[id-im-1]);
                       }
                     px1[id]=acoffx2[ix]*px0[id]-acoffx1[ix]*vp[id]*vp[id]*(1+2*ee)*dtx*xx;
                     pz1[id]=acoffz2[iz]*pz0[id]-acoffz1[iz]*vp[id]*vp[id]*sqrtf(1+2*dd)*dtz*zz;
                     qx1[id]=acoffx2[ix]*qx0[id]-acoffx1[ix]*vp[id]*vp[id]*sqrtf(1+2*dd)*dtx*xx;
                     qz1[id]=acoffz2[iz]*qz0[id]-acoffz1[iz]*vp[id]*vp[id]*dtz*zz;

                     P[id]=px1[id]+pz1[id];
                     Q[id]=qx1[id]+qz1[id];
                   }
                 }
} 
//###################################model#######################################
extern "C"  void model(int nx, int nz,int dx,int dz,int npd,int mm,
           const char FNv[],const char FNe[],const char FNd[],
           int favg,int ns,int fs,int ds,int zs,
           const char FNshot[],const char FNsnap[],int nt, int dt,int run_count)
{
  float dx_,dz_,favg_,dt_,pfac;

  dx_=(float)dx;
  dz_=(float)dz;
  favg_=(float)favg;

  printf("##### model start #####\n");
  printf("#  nx=%2d, dx=%.2f, npd=%d\n",nx,dx_,npd);
  printf("#  nz=%2d, dz=%.2f, mm=%d\n",nz,dz_,mm);
  printf("#     vel=<%s>\n",FNv);
  printf("#  epsilu=<%s>\n",FNe);
  printf("#    deta=<%s>\n",FNd);
  printf("#  favg=%.2f\n",favg_);
  printf("#  ns=%3d\n",ns);
  printf("#  fs=%3d\n",fs);
  printf("#  ds=%3d\n",ds);
  printf("#  zs=%3d\n",zs);
  printf("#    shot=<%s>\n",FNshot);
  printf("#    snap=<%s>\n",FNsnap);


	 FILE *fpsnap, *fpshot;
        fpshot=fopen(FNshot,"wb");
        fpsnap=fopen(FNsnap,"wb");


	int is, it, nnx, nnz, nt_, wtype;
	float *v, *e, *d, t;
	float *vp, *epsilu, *deta;
	float *u0, *u1, *px0, *qx0, *px1, *qx1;
       float *w0, *w1, *pz0, *qz0, *pz1, *qz1;
	float *P, *Q, *shot_Dev, *shot_Hos;

       float *coffx1,*coffx2,*coffz1,*coffz2,*acoffx1,*acoffx2,*acoffz1,*acoffz2;

       clock_t start, end;

          wtype=1;
	   nt_=nt;    
          dt_=(float)(dt*1.0/1000000);
          pfac=10.0;

          nnx=nx+2*npd;
          nnz=nz+2*npd;

    	 v=(float*)malloc(nnz*nnx*sizeof(float));
    	 e=(float*)malloc(nnz*nnx*sizeof(float));
    	 d=(float*)malloc(nnz*nnx*sizeof(float));
    	 shot_Hos=(float*)malloc(nt_*nx*sizeof(float));
        read_file(FNv,FNe,FNd,nx,nz,nnx,nnz,v,e,d,npd);
        pad_vv(nx,nz,nnx,nnz,npd,e);
        pad_vv(nx,nz,nnx,nnz,npd,d);
        pad_vv(nx,nz,nnx,nnz,npd,v); 

        cudaSetDevice(0);// initialize device, default device=0;
	 if(run_count==0)check_gpu_error("Failed to initialize device!");

/****************************/
        cudaMalloc(&vp, nnz*nnx*sizeof(float));
        cudaMalloc(&epsilu, nnz*nnx*sizeof(float));
        cudaMalloc(&deta, nnz*nnx*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilu, e, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(deta, d, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
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

        cudaMalloc(&shot_Dev, nx*nt_*sizeof(float));

	 if(run_count==0)check_gpu_error("Failed to allocate memory for variables!");


        get_d0<<<1, 1>>>(dx_, dz_, nnx, nnz, npd, vp);

        initial_coffe<<<(nnx+511)/512, 512>>>(dt_,nx,coffx1,coffx2,acoffx1,acoffx2,npd);
        initial_coffe<<<(nnz+511)/512, 512>>>(dt_,nz,coffz1,coffz2,acoffz1,acoffz2,npd);

        printf("--------------------------------------------------------\n");
        printf("---   \n");   
        start = clock(); 

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

     cudaMemset(shot_Dev, 0, nt_*nx*sizeof(float));

     for(it=0,t=dt_;it<nt_;it++,t+=dt_)
     { 
     // if(it%100==0&&is==1)printf("---   is===%d   it===%d\n",is,it);

	 add_source<<<1,1>>>(pfac,fs,zs,nx,nz,nnx,nnz,dt_,t,favg_,wtype,npd,is,ds,P,Q);     
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npd,mm,dt_,dx_,dz_,u0,w0,u1,w1,P,Q,coffx1,coffx2,coffz1,coffz2);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt_,dx_,dz_,u1,w1,P,Q,vp,npd,mm,px1,px0,pz1,pz0,qx1,qx0,qz1,qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,fs,ds,zs,is,true);
        u0=u1; w0=w1; px0=px1; pz0=pz1; qx0=qx1; qz0=qz1; 

        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npd, it, nt_, P, shot_Dev);


           if((is==1)&&(it%50==0))
            {
	       cudaMemcpy(e, P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fseek(fpsnap,(int)(it/50)*(nnx)*(nnz)*4L,0);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }
     }//it loop end
      mute_directwave<<<(nx*nt_+511)/512, 512>>>(nx,nt_,dt_,favg_,dx_,dz_,fs,ds,zs,is,vp,epsilu,shot_Dev,100);
      cudaMemcpy(shot_Hos, shot_Dev, nt_*nx*sizeof(float), cudaMemcpyDeviceToHost);
      fseek(fpshot,(is-1)*nt_*nx*sizeof(float),0);
      fwrite(shot_Hos,sizeof(float),nt_*nx,fpshot);

    }

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
       cudaFree(epsilu);
       cudaFree(deta);

       
/***************host free*****************/
	free(v);	free(e);	free(d);
       free(shot_Hos);


     //  exit(0);

}
