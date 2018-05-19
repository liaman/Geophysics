//a#########################################################
//a##         3D Acoustic VTI Medium RTM 
//a##    
//a##  Ps :GPU(CUDA)  ,-SV    ,3D laplace filter
//a##
//a##/*a***************************
//a##Function for VTI medium modeling,
//a##
//a## Ps:  the function of modeling following:
//a##      
//a##          du/dt=1/rho*dp/dx , 
//a##          dv/dt=1/rho*dp/dy , 
//a##          dw/dt=1/rho*dq/dz ,  
//a##          dp/dt=rho*vpx^2*(du/dx+dv/dy)+rho*vp*vpn*dw/dz ,
//a##          dq/dt=rho*vp*vpn*(du/dx+dv/dy)+rho*vp^2*dw/dz ,
//a##                     vpx^2=vp^2*(1+2*epsilu);
//a##                     vpn^2=vp^2*(1+2*deta);
//a##  
//a##*********a*******************/
//a##
//a##                                   code by Rong Tao 
//a##                            
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

#define pi 3.141592653

#define BlockSize1 16// tile size in 1st-axis
#define BlockSize2 16// tile size in 2nd-axis

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
__global__ void add_source(float pfac,int fsx,int fsy,int sz,int nx,int ny,int nz,int nnx,int nny,int nnz,float dt,float t,
                        float favg,int wtype,int npml,int is,int dsx,int dsy,float *P,float *Q,int nsx)
/*< generate ricker wavelet with time deley >*/
{
       int ixs,iys,izs;
       float x_,xx_,tdelay,ts,source=0.0,sx,sy; 
  
       tdelay=1.0/favg;
       ts=t-tdelay;

      // sx=fsx+is%nsx*dsx;
      // sy=fsy+is/nsx*dsy;

       sx=fsx+is*dsx;
       sy=fsy+is*dsy;

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
	     ixs = sx+npml-1;
	     iys = sy+npml-1;
            izs = sz+npml-1;
            P[izs+ixs*nnz+iys*nnz*nnx]+=pfac*source;
            Q[izs+ixs*nnz+iys*nnz*nnx]+=pfac*source;
       }
}
/*******************func*********************/
__global__ void update_vel(int nx,int ny,int nz,int nnx,int nny,int nnz,int npml,float dt,float dx,float dy,float dz,
                           float *u0,float *v0,float *w0,float *u1,float *v1,float *w1,float *P,float *Q,
                           float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,float *coffz2)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy,im;
	float dtx,dty,dtz,xx,yy,zz;

		 dtx=dt/dx;
		 dty=dt/dy;
		 dtz=dt/dz;

       for(iy=0;iy<nny;iy++)
        {
               id=iz+ix*nnz+iy*nnz*nnx;
               if(id>=mm&&id<nnx*nny*nnz-mm)
                 {
                   if(ix>=mm&&ix<(nnx-mm)&&iy>=mm&&iy<(nny-mm)&&iz>=mm&&iz<(nnz-mm))
                    {
                     xx=0.0;
                     yy=0.0;
                     zz=0.0;
	             for(im=0;im<mm;im++)
                      {
                        yy+=c[im]*(P[id+(im+1)*nnz*nnx] - P[id-im*nnz*nnx]);
                        xx+=c[im]*(P[id+(im+1)*nnz]     - P[id-im*nnz]);
                        zz+=c[im]*(Q[id+im+1]           - Q[id-im]);
                      }
                     u1[id]=coffx2[ix]*u0[id]-coffx1[ix]*dtx*xx;
                     v1[id]=coffy2[iy]*v0[id]-coffy1[iy]*dty*yy;
                     w1[id]=coffz2[iz]*w0[id]-coffz1[iz]*dtz*zz;
                   }
                 }
        }  



}
/*******************func***********************/
__global__ void update_stress(int nx,int ny,int nz,int nnx,int nny,int nnz,float dt,float dx,float dy,float dz,
                           float *u1,float *v1,float *w1,float *P,float *Q,float *vp,int npml,
                           float *px1,float *px0,float *py1,float *py0,float *pz1,float *pz0,
                           float *qx1,float *qx0,float *qy1,float *qy0,float *qz1,float *qz0,
                           float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2,
                           float *deta,float *epsilu,int fsx,int dsx,int fsy,int dsy,int zs,int is,int nsx,bool SV)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy,im,rx,ry,rz,R=15,r=4;
	float dtx,dty,dtz,xx,yy,zz,ee,dd;

		 dtx=dt/dx;
		 dty=dt/dy;
		 dtz=dt/dz;

       for(iy=0;iy<nny;iy++)
        {
               id=iz+ix*nnz+iy*nnz*nnx;
               if(id>=mm&&id<nnx*nnz*nny-mm)
                 {
/************************i****************************************/
/************************iso circle start*************************/
                  // rx=ix-(fsx+is%nsx*dsx+npml-1);
                  // ry=iy-(fsy+is/nsx*dsy+npml-1);
                   rx=ix-(fsx+is*dsx+npml-1);
                   ry=iy-(fsy+is*dsy+npml-1);
                   rz=iz-(zs+npml-1);
                   if(SV){
                       if((rx*rx+ry*ry+rz*rz)<=R*R){
                           if((rx*rx+ry*ry+rz*rz)<=r*r){
                               ee = 0.0;
                               dd = 0.0;
                           }else{
                               ee = 0.5*(1-cos(pi*((sqrtf(rx*rx+ry*ry+rz*rz)-r)*4.0/(R*3.0-1))))*epsilu[id];
                               dd = 0.5*(1-cos(pi*((sqrtf(rx*rx+ry*ry+rz*rz)-r)*4.0/(R*3.0-1))))*deta[id]; 
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
                   if(ix>=mm&&ix<(nnx-mm)&&iy>=mm&&iy<(nny-mm)&&iz>=mm&&iz<(nnz-mm))
                     {
                     xx=0.0;
                     yy=0.0;
                     zz=0.0;
	             for(im=0;im<mm;im++)
                       {
                        yy+=c[im]*(v1[id+im*nnz*nnx] - v1[id-(im+1)*nnz*nnx]);
                        xx+=c[im]*(u1[id+im*nnz]     - u1[id-(im+1)*nnz]);
                        zz+=c[im]*(w1[id+im]         - w1[id-im-1]);
                       }
                     px1[id]=acoffx2[ix]*px0[id] - acoffx1[ix]*vp[id]*vp[id]*(1+2*ee)*dtx*xx;
                     py1[id]=acoffy2[iy]*py0[id] - acoffy1[iy]*vp[id]*vp[id]*(1+2*ee)*dty*yy;
                     pz1[id]=acoffz2[iz]*pz0[id] - acoffz1[iz]*vp[id]*vp[id]*sqrtf(1+2*dd)*dtz*zz;

                     qx1[id]=acoffx2[ix]*qx0[id] - acoffx1[ix]*vp[id]*vp[id]*sqrtf(1+2*dd)*dtx*xx;
                     qy1[id]=acoffy2[iy]*qy0[id] - acoffy1[iy]*vp[id]*vp[id]*sqrtf(1+2*dd)*dty*yy;
                     qz1[id]=acoffz2[iz]*qz0[id] - acoffz1[iz]*vp[id]*vp[id]*dtz*zz;

                     P[id]=px1[id]+py1[id]+pz1[id];
                     Q[id]=qx1[id]+qy1[id]+qz1[id];
                   }
                 }
         }
}                      
/********************func**********************/
__global__ void get_d0(float dx,float dy,float dz,int nnx,int nny,int nnz,int npml,float *vp)
{
   d0=10.0*vp[nny*nnx*nnz/2]*log(100000.0)/(2.0*npml*((dx+dy+dz)/3.0));
}
/*************func*******************/
void pad_vv(int nx,int ny,int nz,int nnx,int nny,int nnz,int npml,float *ee)
{
     int ix,iy,iz,id;
 
	    for(iy=0;iy<nny;iy++)
		 for(ix=0;ix<nnx;ix++)
		 {
			 for(iz=0;iz<nnz;iz++)
			 {
                             id=iz+ix*nnz+iy*nnz*nnx;

                             if(ix<npml){
                                ee[id]=ee[iz+npml*nnz+iy*nnz*nnx];  //left
                             }else if(ix>=nnx-npml){
                                ee[id]=ee[iz+(nnx-npml-1)*nnz+iy*nnz*nnx];//right
                                 }
			 }
		 }
	    for(iy=0;iy<nny;iy++)
		 for(ix=0;ix<nnx;ix++)
		 {
			 for(iz=0;iz<nnz;iz++)
			 {
                             id=iz+ix*nnz+iy*nnz*nnx;

                             if(iy<npml){
                                ee[id]=ee[iz+ix*nnz+npml*nnz*nnx];  //front
                             }else if(iy>=nny-npml){
                                ee[id]=ee[iz+ix*nnz+(nny-npml-1)*nnz*nnx];//back
                                 }
			 }
		 }
	    for(iy=0;iy<nny;iy++)
		 for(ix=0;ix<nnx;ix++)
		 {
			 for(iz=0;iz<nnz;iz++)
			 {
                             id=iz+ix*nnz+iy*nnz*nnx;

                             if(iz<npml){
                                ee[id]=ee[npml+ix*nnz+iy*nnz*nnx];  //up
                             }else if(iz>=nnz-npml){
                                ee[id]=ee[nnz-npml-1+ix*nnz+iy*nnz*nnx];//down
                                 }
			 }
		 }

}
/*************func*******************/
void read_file(char FN1[],char FN2[],char FN3[],int nx,int ny,int nz,int nnx,int nny,int nnz,float *vv,float *epsilu,float *deta,int npml)
{
		 int ix,iy,iz,id;
		
		 FILE *fp1,*fp2,*fp3;
		 if((fp1=fopen(FN1,"rb"))==NULL){printf("error open <%s>!\n",FN1);exit(0);}
		 if((fp2=fopen(FN2,"rb"))==NULL){printf("error open <%s>!\n",FN2);exit(0);}
		 if((fp3=fopen(FN3,"rb"))==NULL){printf("error open <%s>!\n",FN3);exit(0);}

	    for(iy=npml;iy<ny+npml;iy++)
		 for(ix=npml;ix<nx+npml;ix++)
		 {
			 for(iz=npml;iz<nz+npml;iz++)
			 {
                             id=iz+ix*nnz+iy*nnz*nnx;
				 fread(&vv[id],4L,1,fp1);
				 fread(&epsilu[id],4L,1,fp2);
				 fread(&deta[id],4L,1,fp3);
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
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
__global__ void shot_record(int nnx,int nny, int nnz, int nx,int ny, int nz, int npml, int it, int nt,
                            float *P, float *Q, float *shot, bool flag)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

        int ix=id%nx;
        int iy=id/nx;

           if(id<nx*ny)
            {
              if(flag)
                {
                 shot[it+nt*ix+nt*nx*iy]=P[npml+nnz*(ix+npml)+nnz*nnx*(iy+npml)]
                                        +Q[npml+nnz*(ix+npml)+nnz*nnx*(iy+npml)];
              }else{
                 P[npml-1+nnz*(ix+npml)+nnz*nnx*(iy+npml)]=shot[it+nt*ix+nt*nx*iy];
                 Q[npml-1+nnz*(ix+npml)+nnz*nnx*(iy+npml)]=shot[it+nt*ix+nt*nx*iy];
                }
            }       
}
/*************func**************/ 
void window3d(float *a, float *b, int nz, int nx, int ny, int nnz, int nnx, int npml)
/*< window a 3d subvolume >*/
{
	int iz, ix, iy;
	
	for(iy=0; iy<ny; iy++)
	for(ix=0; ix<nx; ix++)
	for(iz=0; iz<nz; iz++)
	{
		a[iz+nz*ix+nz*nx*iy]=b[(iz+npml)+nnz*(ix+npml)+nnz*nnx*(iy+npml)];
	}
}
/*************func**************/    
__global__ void mute_directwave(int nx,int ny,int nt,float dt,float favg, float dx,float dy,float dz,int fsx,int fsy,int dsx,int dsy,
                                int zs,int is, float *vp,float *epsilu,float *shot,int tt,int nsx)
{

    const int ix = blockIdx.x * blockDim.x + threadIdx.x;
    const int iy = blockIdx.y * blockDim.y + threadIdx.y;

    int id,it;
    int mu_t,mu_nt;
    float mu_x,mu_y,mu_z,mu_t0;

       for(it=0;it<nt;it++)
        {
          id=it+ix*nt+iy*nx*nt;
          if(ix<nx&&iy<ny&&it<nt)
            {
            //  mu_x=dx*abs(ix-fsx-(is%nsx)*dsx);
            //  mu_y=dy*abs(iy-fsy-(is/nsx)*dsy);
              mu_x=dx*abs(ix-fsx-is*dsx);
              mu_y=dy*abs(iy-fsy-is*dsy);
              mu_z=dz*zs;
              mu_t0=sqrtf(pow(mu_x,2)+pow(mu_y,2)+pow(mu_z,2))/(vp[1]*sqrtf(1+2*epsilu[1]));
              mu_t=(int)(2.0/(dt*favg));
              mu_nt=(int)(mu_t0/dt)+mu_t+tt;

           if((it>(int)(mu_t0/dt)-tt)&&(it<mu_nt))
                    shot[id]=0.0;
            }
        }
/*    int id=threadIdx.x+blockDim.x*blockIdx.x;

    int mu_t,mu_nt;
    float mu_x,mu_y,mu_z,mu_t0;

    int ix=(id/nt)%nx;
    int iy=(id/nt)/nx;
    int it=id%nt;

   if(id<nx*ny*nt)
   {
        mu_x=dx*abs(ix-fsx-(is%nsx)*dsx);
        mu_y=dy*abs(iy-fsy-(is/nsx)*dsy);
        mu_z=dz*zs;
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_y,2)+pow(mu_z,2))/(vp[1]*sqrtf(1+2*epsilu[1]));
        mu_t=(int)(2.0/(dt*favg));
        mu_nt=(int)(mu_t0/dt)+mu_t+tt;

           if(it<mu_nt)
              shot[id]=0.0;
   }  */
}
/*************func*******************/
__global__ void wavefield_bndr(int nnx, int nny, int nnz, int nx, int ny, int nz, int npml, int it, int nt, 
                               float *P, float *Q, float *P_bndr, float *Q_bndr, bool flag)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;
        int ix,iy,iz;

           if(id<2*nx*ny+2*nz*ny+2*nx*nz)
            {
            if(flag)/////////////////////////////////save boundary
             {
              if(id<nx*ny){//up

               ix=id%nx;
               iy=id/nx;
               P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=P[npml-1+nnz*(ix+npml)+nnz*nnx*(iy+npml)];
               Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=Q[npml-1+nnz*(ix+npml)+nnz*nnx*(iy+npml)];

              }else if(id>=nx*ny&&id<(2*nx*ny)){//down

               ix=(id-nx*ny)%nx;
               iy=(id-nx*ny)/nx;
               P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=P[npml+nz+nnz*(ix+npml)+nnz*nnx*(iy+npml)];
               Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=Q[npml+nz+nnz*(ix+npml)+nnz*nnx*(iy+npml)];

              }else if(id>=(2*nx*ny)&&id<(2*nx*ny+nz*ny)){//left

               iz=(id-2*nx*ny)%nz;
               iy=(id-2*nx*ny)/nz;
               P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=P[npml+iz+nnz*(npml-1)+nnz*nnx*(iy+npml)];
               Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=Q[npml+iz+nnz*(npml-1)+nnz*nnx*(iy+npml)];

              }else if(id>=(2*nx*ny+nz*ny)&&id<(2*nx*ny+2*nz*ny)){//right

               iz=(id-2*nx*ny-nz*ny)%nz;
               iy=(id-2*nx*ny-nz*ny)/nz;
               P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=P[npml+iz+nnz*(nx+npml)+nnz*nnx*(iy+npml)];
               Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=Q[npml+iz+nnz*(nx+npml)+nnz*nnx*(iy+npml)];

              }else if(id>=(2*nx*ny+2*nz*ny)&&id<(2*nx*ny+2*nz*ny+nx*nz)){//front

               iz=(id-2*nx*ny-2*nz*ny)%nz;
               ix=(id-2*nx*ny-2*nz*ny)/nz;
               P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=P[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml-1)];
               Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=Q[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml-1)];

                }else if(id>=(2*nx*ny+2*nz*ny+nx*nz)&&id<(2*nx*ny+2*nz*ny+2*nx*nz)){//back

               iz=(id-2*nx*ny-2*nz*ny-nx*nz)%nz;
               ix=(id-2*nx*ny-2*nz*ny-nx*nz)/nz;
               P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=P[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml+ny)];
               Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id]=Q[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml+ny)];

                }

             }else{

              if(id<nx*ny){//up

               ix=id%nx;
               iy=id/nx;
               P[npml-1+nnz*(ix+npml)+nnz*nnx*(iy+npml)]=P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];
               Q[npml-1+nnz*(ix+npml)+nnz*nnx*(iy+npml)]=Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];

              }else if(id>=nx*ny&&id<(2*nx*ny)){//down

               ix=(id-nx*ny)%nx;
               iy=(id-nx*ny)/nx;
               P[npml+nz+nnz*(ix+npml)+nnz*nnx*(iy+npml)]=P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];
               Q[npml+nz+nnz*(ix+npml)+nnz*nnx*(iy+npml)]=Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];

              }else if(id>=(2*nx*ny)&&id<(2*nx*ny+nz*ny)){//left

               iz=(id-2*nx*ny)%nz;
               iy=(id-2*nx*ny)/nz;
               P[npml+iz+nnz*(npml-1)+nnz*nnx*(iy+npml)]=P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];
               Q[npml+iz+nnz*(npml-1)+nnz*nnx*(iy+npml)]=Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];

              }else if(id>=(2*nx*ny+nz*ny)&&id<(2*nx*ny+2*nz*ny)){//right

               iz=(id-2*nx*ny-nz*ny)%nz;
               iy=(id-2*nx*ny-nz*ny)/nz;
               P[npml+iz+nnz*(nx+npml)+nnz*nnx*(iy+npml)]=P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];
               Q[npml+iz+nnz*(nx+npml)+nnz*nnx*(iy+npml)]=Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];

              }else if(id>=(2*nx*ny+2*nz*ny)&&id<(2*nx*ny+2*nz*ny+nx*nz)){//front

               iz=(id-2*nx*ny-2*nz*ny)%nz;
               ix=(id-2*nx*ny-2*nz*ny)/nz;
               P[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml-1)]=P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];
               Q[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml-1)]=Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];

                }else if(id>=(2*nx*ny+2*nz*ny+nx*nz)&&id<(2*nx*ny+2*nz*ny+2*nx*nz)){//back

               iz=(id-2*nx*ny-2*nz*ny-nx*nz)%nz;
               ix=(id-2*nx*ny-2*nz*ny-nx*nz)/nz;
               P[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml+ny)]=P_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];
               Q[npml+iz+nnz*(ix+npml)+nnz*nnx*(npml+ny)]=Q_bndr[it*(2*nx*ny+2*nz*ny+2*nx*nz)+id];

                }
             }
            }       
}
/*************func**************/    
__global__ void cal_migration(int nnx, int nny, int nnz,int nx, int ny, int nz, int npml, float *migration, float *s, float *g)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

       for(iy=0;iy<ny;iy++)
        {
           id=iz+ix*nz+iy*nz*nx;
           if(ix<nx&&iy<ny&&iz<nz)
              migration[id]+=s[iz+npml+nnz*(ix+npml)+nnx*nnz*(iy+npml)]*g[iz+npml+nnz*(ix+npml)+nnx*nnz*(iy+npml)];
        }
}
/*************func**************/    
__global__ void cal_illumination(int nnx, int nny, int nnz,int nx, int ny, int nz, int npml, float *illumination, float *P, float *Q)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

       for(iy=0;iy<ny;iy++)
        {
           id=iz+ix*nz+iy*nz*nx;
           if(ix<nx&&iy<ny&&iz<nz)
              illumination[id]+=P[iz+npml+nnz*(ix+npml)+nnx*nnz*(iy+npml)]*P[iz+npml+nnz*(ix+npml)+nnx*nnz*(iy+npml)]
                               +Q[iz+npml+nnz*(ix+npml)+nnx*nnz*(iy+npml)]*Q[iz+npml+nnz*(ix+npml)+nnx*nnz*(iy+npml)];
        }
}
/*************func**************/    
__global__ void migration_illum(int nnx, int nny, int nnz,int nx, int ny, int nz, int npml, float *illumination, float *migration)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

       for(iy=0;iy<ny;iy++)
        {
           id=iz+ix*nz+iy*nz*nx;
           if(ix<nx&&iy<ny&&iz<nz)
               if(illumination[id]!=0)
                     migration[id]/=illumination[id];
        }
}
/*************func**************/
void laplace_3d_filter(int adj, int nz, int nx,int ny, float *in, float *out)
/*< linear operator, come from Madagascar Mlaplac2>*/
{
    int iz,ix,iy,j;
    for (j=0;j<nx*nz*ny;j++) out[j]=0.0;

  for(iy=0;iy<ny;iy++)
    for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {

	    j = iz+ix*nz+iy*nx*nz;
	    if (iz > 0) {
		if (adj) {
		    out[j-1] -= in[j];
		    out[j]   += in[j];
		} else {
		    out[j] += in[j] - in[j-1];
		}
	    }
	    if (iz < nz-1) {
		if (adj) {
		    out[j+1] -= in[j];
		    out[j]   += in[j];
		} else {
		    out[j] += in[j] - in[j+1];
		}
	    }
	    if (ix > 0) {
		if (adj) {
		    out[j-nz] -= in[j];
		    out[j]    += in[j];
		} else {
		    out[j] += in[j] - in[j-nz];
		}
	    }
	    if (ix < nx-1) {
		if (adj) {
		    out[j+nz] -= in[j];
		    out[j]    += in[j];
		} else {
		    out[j] += in[j] - in[j+nz];
		}
	    }
	    if (iy > 0) {
		if (adj) {
		    out[j-nz*nx] -= in[j];
		    out[j]    += in[j];
		} else {
		    out[j] += in[j] - in[j-nz*nx];
		}
	    }
	    if (iy < ny-1) {
		if (adj) {
		    out[j+nz*nx] -= in[j];
		    out[j]    += in[j];
		} else {
		    out[j] += in[j] - in[j+nz*nx];
		}
	    }
	}
    }
}
/*************func**************/    
__global__ void Poynting_Adcigs(int nnx, int nny, int nnz, int nx, int ny, int nz, int npml, int na, int da,float *adcigs,int dcdp, 
                           float *s_P, float *s_Q, float *s_u, float *s_v, float *s_w, 
                           float *g_P, float *g_Q, float *g_u, float *g_v, float *g_w)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy,ia;

    float Ssx, Ssy, Ssz, Sgx, Sgy, Sgz, b1, b2, a;

       for(iy=0;iy<(int)(ny/dcdp);iy++)
        {

           if(ix<(int)(nx/dcdp)&&iz<nz)
            {
               ia=0;
               Ssx=-s_P[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]*s_u[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)];
               Ssy=-s_P[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]*s_v[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)];
               Ssz=-s_Q[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]*s_w[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)];
               Sgx= g_P[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]*g_u[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)];
               Sgy= g_P[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]*g_v[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)];
               Sgz= g_Q[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]*g_w[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)];

               b1= Ssx*Ssx + Ssy*Ssy + Ssz*Ssz;
               b2= Sgx*Sgx + Sgy*Sgy + Sgz*Sgz;
                a=(Ssx*Sgx + Ssy*Sgy + Ssz*Sgz)/(sqrtf(b1*b2)*(1 - 0.1));

               if(a>=-1&&a<=1)
                 {
                   a=0.5*acosf(a)*180.0/pi;
                   ia=(int)(a/(da*1.0));
                   if(ia<na)
                     {
                        id=iz+ia*nz+ix*na*nz+iy*nz*na*((int)(nx/dcdp));
                        adcigs[id] += s_P[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]
                                     *g_P[iz+npml+nnz*(ix*dcdp+npml)+nnz*nnx*(iy*dcdp+npml)]
                                     *cosf(ia*pi/180.0)*cosf(ia*pi/180.0)*cosf(ia*pi/180.0);
                     }
                 }
            }

        }
}
/*************func**************/    
__global__ void adcigs_illum(int nnx, int nny, int nnz, int nx, int ny, int nz, int npml, int na, float *adcigs,int dcdp, float *illum)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy,ia;

       for(iy=0;iy<(int)(ny/dcdp);iy++)
        {
          for(ia=0;ia<na;ia++)
           {

             if(ix<(int)(nx/dcdp)&&iz<nz)
               {
               
                        id=iz+ia*nz+ix*na*nz+iy*nz*na*((int)(nx/dcdp));
                        if(illum[iz+ix*dcdp*nz+iy*dcdp*nx*nz]!=0)
                             adcigs[id] /=illum[iz+ix*dcdp*nz+iy*dcdp*nx*nz];
               }
           }
       }
}
//a########################################################################
int main(int argc,char *argv[])
{
	int is, it, nx, ny, nz, nnx, nny, nnz, nt, wtype, na, da, dcdp, nxa, nya;
	int ns, nsx, dsx, fsx, dsy, fsy, zs, npml;
	float dx, dy, dz, dt, t, pfac, favg;

	float *v, *e, *d;
	float *vp, *epsilu, *deta;
	float *s_u0, *s_u1, *s_px0, *s_qx0, *s_px1, *s_qx1;
	float *s_v0, *s_v1, *s_py0, *s_qy0, *s_py1, *s_qy1;
       float *s_w0, *s_w1, *s_pz0, *s_qz0, *s_pz1, *s_qz1;
	float *s_P, *s_Q;
	float *g_u0, *g_u1, *g_px0, *g_qx0, *g_px1, *g_qx1;
	float *g_v0, *g_v1, *g_py0, *g_qy0, *g_py1, *g_qy1;
       float *g_w0, *g_w1, *g_pz0, *g_qz0, *g_pz1, *g_qz1;
	float *g_P, *g_Q;
       float *s_P_bndr, *s_Q_bndr;
       float *shot_Dev, *shot_Hos;
       float *migration, *illumination, *adcigs;
       float *Atemp;


       float *coffx1,*coffx2,*coffy1,*coffy2,*coffz1,*coffz2;
       float *acoffx1,*acoffx2,*acoffy1,*acoffy2,*acoffz1,*acoffz2;

       clock_t start, end, is_t0, is_t1;
/*************wavelet\boundary**************/
          wtype=1;npml=20;
/********** dat document ***********/
          char FN1[250]={"waxian_vel_201201201.dat"};
          char FN2[250]={"waxian_eps_201201201.dat"};
          char FN3[250]={"waxian_del_201201201.dat"};
	   char FN4[250]={"waxian_shot.dat"};
	   char FN5[250]={"waxian_snap.dat"};
	   char FN6[250]={"waxian_migration.dat"};
	   char FN7[250]={"waxian_migration_laplace.dat"};
	   char FN8[250]={"waxian_illumination.dat"};
	   char FN9[250]={"waxian_adcigs.dat"};

/********aaa************/  
	 FILE *fpsnap, *fpshot, *fpmig,*fpmigla, *fpillum, *fpadcigs;
        fpshot=fopen(FN4,"wb");
        fpsnap=fopen(FN5,"wb");
        fpmig=fopen(FN6,"wb");
        fpmigla=fopen(FN7,"wb");
        fpillum=fopen(FN8,"wb");
        fpadcigs=fopen(FN9,"wb");

 
/********* parameters *************/

          nx=201; 
          ny=201;              
	   nz=201;         favg=60;     pfac=10.0;

 	   dx=5.0;  
 	   dy=5.0;   
          dz=5.0;   
     
	   nt=1501;    
          dt=0.0005;
     
          ns=625;          nsx=25;  
         // fsx=nx/nsx/2;    dsx=nx/nsx;         
         // fsy=ny/(ns/nsx)/2;    dsy=ny/(ns/nsx);
          fsx=4;//nx/ns/2;    
          dsx=8;//nx/ns;         
          fsy=4;//ny/ns/2;//200;//100;//ny/ns/2;    
          dsy=8;//ny/ns;//0;//ny/ns;
          zs=1;    

          na=65; 
          da=1;
          dcdp=1;
/*************v***************/ 
          nnx=nx+2*npml;
          nny=ny+2*npml;
          nnz=nz+2*npml;
          nxa=(int)(nx/dcdp);
          nya=(int)(ny/dcdp);
/************a*************/
    	 Atemp=(float*)malloc(nz*nxa*nya*na*sizeof(float));

    	 v=(float*)malloc(nnz*nnx*nny*sizeof(float));
    	 e=(float*)malloc(nnz*nnx*nny*sizeof(float));
    	 d=(float*)malloc(nnz*nnx*nny*sizeof(float));
    	 shot_Hos=(float*)malloc(nt*nx*ny*sizeof(float));
        read_file(FN1,FN2,FN3,nx,ny,nz,nnx,nny,nnz,v,e,d,npml);
/****************************/
        pad_vv(nx,ny,nz,nnx,nny,nnz,npml,e);
        pad_vv(nx,ny,nz,nnx,nny,nnz,npml,d);
        pad_vv(nx,ny,nz,nnx,nny,nnz,npml,v);

        cudaSetDevice(0);// initialize device, default device=0;
	 check_gpu_error("Failed to initialize device!");

	dim3 dimg, Xdimg, Adimg, dimb;
	dimg.x=(nnz+BlockSize1-1)/BlockSize1;
	dimg.y=(nnx+BlockSize2-1)/BlockSize2;
	Xdimg.x=(nnx+BlockSize1-1)/BlockSize1;
	Xdimg.y=(nny+BlockSize2-1)/BlockSize2;
	Adimg.x=(nz+BlockSize1-1)/BlockSize1;
	Adimg.y=(nxa+BlockSize2-1)/BlockSize2;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;
/****************************/
        cudaMalloc(&vp, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&epsilu, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&deta, nnz*nnx*nny*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilu, e, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(deta, d, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
/****************************/
        cudaMalloc(&s_u0, nnz*nnx*nny*sizeof(float));    cudaMalloc(&s_u1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_v0, nnz*nnx*nny*sizeof(float));    cudaMalloc(&s_v1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_w0, nnz*nnx*nny*sizeof(float));    cudaMalloc(&s_w1, nnz*nnx*nny*sizeof(float));

        cudaMalloc(&s_P, nnz*nnx*nny*sizeof(float));     cudaMalloc(&s_Q, nnz*nnx*nny*sizeof(float));

        cudaMalloc(&s_px0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&s_px1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_py0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&s_py1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_pz0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&s_pz1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_qx0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&s_qx1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_qy0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&s_qy1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&s_qz0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&s_qz1, nnz*nnx*nny*sizeof(float));

        cudaMalloc(&g_u0, nnz*nnx*nny*sizeof(float));    cudaMalloc(&g_u1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_v0, nnz*nnx*nny*sizeof(float));    cudaMalloc(&g_v1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_w0, nnz*nnx*nny*sizeof(float));    cudaMalloc(&g_w1, nnz*nnx*nny*sizeof(float));

        cudaMalloc(&g_P, nnz*nnx*nny*sizeof(float));     cudaMalloc(&g_Q, nnz*nnx*nny*sizeof(float));

        cudaMalloc(&g_px0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&g_px1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_py0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&g_py1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_pz0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&g_pz1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_qx0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&g_qx1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_qy0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&g_qy1, nnz*nnx*nny*sizeof(float));
        cudaMalloc(&g_qz0, nnz*nnx*nny*sizeof(float));   cudaMalloc(&g_qz1, nnz*nnx*nny*sizeof(float));

        cudaMalloc(&coffx1, nnx*sizeof(float));     cudaMalloc(&coffx2, nnx*sizeof(float));
        cudaMalloc(&coffy1, nnx*sizeof(float));     cudaMalloc(&coffy2, nnx*sizeof(float));
        cudaMalloc(&coffz1, nnz*sizeof(float));     cudaMalloc(&coffz2, nnz*sizeof(float));
        cudaMalloc(&acoffx1, nnx*sizeof(float));    cudaMalloc(&acoffx2, nnx*sizeof(float));
        cudaMalloc(&acoffy1, nnx*sizeof(float));    cudaMalloc(&acoffy2, nnx*sizeof(float));
        cudaMalloc(&acoffz1, nnz*sizeof(float));    cudaMalloc(&acoffz2, nnz*sizeof(float));

        cudaMalloc(&s_P_bndr, nt*(2*nx*ny+2*nz*ny+2*nx*nz)*sizeof(float));
        cudaMalloc(&s_Q_bndr, nt*(2*nx*ny+2*nz*ny+2*nx*nz)*sizeof(float));

        cudaMalloc(&migration, nz*nx*ny*sizeof(float)); 
        cudaMalloc(&illumination, nz*nx*ny*sizeof(float)); 
        cudaMalloc(&adcigs, nz*nxa*nya*na*sizeof(float)); 

        cudaMalloc(&shot_Dev, nx*ny*nt*sizeof(float));
/******************************/
	 check_gpu_error("Failed to allocate memory for variables!");

        get_d0<<<1, 1>>>(dx,dy,dz,nnx,nny,nnz,npml,vp);
        initial_coffe<<<(nnx+511)/512, 512>>>(dt,nx,coffx1,coffx2,acoffx1,acoffx2,npml);
        initial_coffe<<<(nny+511)/512, 512>>>(dt,ny,coffy1,coffy2,acoffy1,acoffy2,npml);
        initial_coffe<<<(nnz+511)/512, 512>>>(dt,nz,coffz1,coffz2,acoffz1,acoffz2,npml);

        cudaMemset(migration, 0, nz*nx*ny*sizeof(float)); 
        cudaMemset(illumination, 0, nz*nx*ny*sizeof(float)); 
        cudaMemset(adcigs, 0, nz*na*nxa*nya*sizeof(float)); 

        printf("--------------------------------------------------------\n");
        printf("---   \n");   
        start = clock();                                  
/**********IS Loop start*******/
   for(is=0;is<ns;is++)	
    {     
       //  printf("---   IS=%3d  \n",is);
     is_t0 = clock(); 

     cudaMemset(s_u0, 0, nnz*nnx*nny*sizeof(float));     cudaMemset(s_u1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_v0, 0, nnz*nnx*nny*sizeof(float));     cudaMemset(s_v1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_w0, 0, nnz*nnx*nny*sizeof(float));     cudaMemset(s_w1, 0, nnz*nnx*nny*sizeof(float));

     cudaMemset(s_P, 0, nnz*nnx*nny*sizeof(float));      cudaMemset(s_Q, 0, nnz*nnx*nny*sizeof(float));

     cudaMemset(s_px0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(s_px1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_py0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(s_py1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_pz0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(s_pz1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_qx0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(s_qx1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_qy0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(s_qy1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(s_qz0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(s_qz1, 0, nnz*nnx*nny*sizeof(float));

     cudaMemset(g_u0, 0, nnz*nnx*nny*sizeof(float));     cudaMemset(g_u1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_v0, 0, nnz*nnx*nny*sizeof(float));     cudaMemset(g_v1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_w0, 0, nnz*nnx*nny*sizeof(float));     cudaMemset(g_w1, 0, nnz*nnx*nny*sizeof(float));

     cudaMemset(g_P, 0, nnz*nnx*nny*sizeof(float));      cudaMemset(g_Q, 0, nnz*nnx*nny*sizeof(float));

     cudaMemset(g_px0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(g_px1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_py0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(g_py1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_pz0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(g_pz1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_qx0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(g_qx1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_qy0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(g_qy1, 0, nnz*nnx*nny*sizeof(float));
     cudaMemset(g_qz0, 0, nnz*nnx*nny*sizeof(float));    cudaMemset(g_qz1, 0, nnz*nnx*nny*sizeof(float));

     cudaMemset(s_P_bndr, 0, nt*(2*nx*ny+2*nz*ny+2*nx*nz)*sizeof(float));
     cudaMemset(s_Q_bndr, 0, nt*(2*nx*ny+2*nz*ny+2*nx*nz)*sizeof(float));

     cudaMemset(shot_Dev, 0, nt*nx*ny*sizeof(float));

     for(it=0,t=dt;it<nt;it++,t+=dt)
     { 
      if(it%100==0)printf("---   IS===%d   it===%d\n",is,it);
        add_source<<<1,1>>>(pfac,fsx,fsy,zs,nx,ny,nz,nnx,nny,nnz,dt,t,favg,wtype,npml,is,dsx,dsy,s_P,s_Q,nsx);
        update_vel<<<dimg,dimb>>>(nx,ny,nz,nnx,nny,nnz,npml,dt,dx,dy,dz,
                                 s_u0,s_v0,s_w0,s_u1,s_v1,s_w1,s_P,s_Q,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
        update_stress<<<dimg,dimb>>>(nx,ny,nz,nnx,nny,nnz,dt,dx,dy,dz,s_u1,s_v1,s_w1,s_P,s_Q,vp,npml,
                                     s_px1,s_px0,s_py1,s_py0,s_pz1,s_pz0,s_qx1,s_qx0,s_qy1,s_qy0,s_qz1,s_qz0,
                                     acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2,deta,epsilu, 
                                     fsx, dsx, fsy, dsy,zs, is, nsx, true);
        s_u0=s_u1; s_v0=s_v1; s_w0=s_w1; s_px0=s_px1; s_py0=s_py1; s_pz0=s_pz1; s_qx0=s_qx1; s_qy0=s_qy1; s_qz0=s_qz1; 

        wavefield_bndr<<<(2*nx*ny+2*nz*ny+2*nx*nz+511)/512, 512>>>(nnx,nny,nnz,nx,ny,nz,npml,it,nt,s_P,s_Q,s_P_bndr,s_Q_bndr,true);
        shot_record<<<(nx*ny+511)/512, 512>>>(nnx,nny, nnz, nx, ny, nz, npml, it, nt, s_P, s_Q, shot_Dev, true);

        cal_illumination<<<dimg,dimb>>>(nnx,nny,nnz,nx,ny,nz,npml,illumination,s_P,s_Q);

    /*       if((is==0)&&(it!=0&&it%100==0))
            {
	       cudaMemcpy(e, s_P, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
            //  fseek(fpsnap,(int)(it/100-1)*(nx*ny*nz)*4L,0);
              window3d(v, e, nz, nx, ny, nnz, nnx, npml);
              fwrite(v,4L,nx*nz*ny,fpsnap);
            }   */
     }//it loop end

      mute_directwave<<<Xdimg,dimb>>>(nx,ny,nt,dt,favg,dx,dy,dz,fsx,fsy,dsx,dsy,zs,is,vp,epsilu,shot_Dev,70,nsx);
      
   //   if(is==0){   
          cudaMemcpy(shot_Hos, shot_Dev, nt*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
         // fseek(fpshot,is*nt*nx*ny*sizeof(float),0);
          fwrite(shot_Hos,sizeof(float),nt*nx*ny,fpshot);
    //   }

    for(it=nt-1;it>=0;it--)
     { 
      if(it%100==0)printf("---   IS===%d   it===%d\n",is,it);
        wavefield_bndr<<<(2*nx*ny+2*nz*ny+2*nx*nz+511)/512, 512>>>(nnx,nny,nnz,nx,ny,nz,npml,it,nt,s_P,s_Q,s_P_bndr,s_Q_bndr,false);
        update_vel<<<dimg,dimb>>>(nx,ny,nz,nnx,nny,nnz,npml,dt,dx,dy,dz,
                                 s_u0,s_v0,s_w0,s_u1,s_v1,s_w1,s_P,s_Q,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
        update_stress<<<dimg,dimb>>>(nx,ny,nz,nnx,nny,nnz,dt,dx,dy,dz,s_u1,s_v1,s_w1,s_P,s_Q,vp,npml,
                                     s_px1,s_px0,s_py1,s_py0,s_pz1,s_pz0,s_qx1,s_qx0,s_qy1,s_qy0,s_qz1,s_qz0,
                                     acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2,deta,epsilu, 
                                     fsx, dsx, fsy, dsy,zs, is, nsx, true);
        s_u0=s_u1; s_v0=s_v1; s_w0=s_w1; s_px0=s_px1; s_py0=s_py1; s_pz0=s_pz1; s_qx0=s_qx1; s_qy0=s_qy1; s_qz0=s_qz1; 

      /*     if((is==0)&&(it!=0&&it%100==0))
            {
	       cudaMemcpy(e, s_P, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
             // fseek(fpsnap,(int)(it/100-1)*(nx*ny*nz)*4L,0);
              window3d(v, e, nz, nx, ny, nnz, nnx, npml);
              fwrite(v,4L,nx*nz*ny,fpsnap);
            }    */

        shot_record<<<(nx*ny+511)/512, 512>>>(nnx,nny, nnz, nx, ny, nz, npml, it, nt, g_P, g_Q, shot_Dev, false);
        update_vel<<<dimg,dimb>>>(nx,ny,nz,nnx,nny,nnz,npml,dt,dx,dy,dz,
                                 g_u0,g_v0,g_w0,g_u1,g_v1,g_w1,g_P,g_Q,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
        update_stress<<<dimg,dimb>>>(nx,ny,nz,nnx,nny,nnz,dt,dx,dy,dz,g_u1,g_v1,g_w1,g_P,g_Q,vp,npml,
                                     g_px1,g_px0,g_py1,g_py0,g_pz1,g_pz0,g_qx1,g_qx0,g_qy1,g_qy0,g_qz1,g_qz0,
                                     acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2,deta,epsilu, 
                                     fsx, dsx, fsy, dsy,zs, is, nsx, true);
        g_u0=g_u1; g_v0=g_v1; g_w0=g_w1; g_px0=g_px1; g_py0=g_py1; g_pz0=g_pz1; g_qx0=g_qx1; g_qy0=g_qy1; g_qz0=g_qz1; 

    /*       if((is==0)&&(it!=0&&it%100==0))
            {
	       cudaMemcpy(e, g_P, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
            //  fseek(fpsnap,(int)(it/100-1)*(nx*ny*nz)*4L,0);
              window3d(v, e, nz, nx, ny, nnz, nnx, npml);
              fwrite(v,4L,nx*nz*ny,fpsnap);
            }    */
        cal_illumination<<<dimg,dimb>>>(nnx,nny,nnz,nx,ny,nz,npml,illumination,g_P,g_Q);
        cal_migration<<<dimg,dimb>>>(nnx,nny,nnz,nx,ny,nz,npml,migration,s_P,g_P);

        Poynting_Adcigs<<<Adimg,dimb>>>(nnx,nny,nnz,nx,ny,nz,npml,na,da,adcigs,dcdp,
                                                s_P,s_Q,s_u0,s_v0,s_w0,g_P,g_Q,g_u0,g_v0,g_w0);

     }//it loop end
    is_t1 = clock(); 
    printf("IS=%3d: %f (min)\n", is, ((float)(is_t1-is_t0))/60.0/CLOCKS_PER_SEC);
    }//is loop end

   migration_illum<<<dimg,dimb>>>(nnx,nny,nnz,nx,ny,nz,npml,illumination,migration);
   adcigs_illum<<<Adimg,dimb>>>(nnx,nny,nnz,nx,ny,nz,npml,na,adcigs,dcdp,illumination);

   cudaMemcpy(v, illumination, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(v,4L,nx*nz*ny,fpillum);

   cudaMemcpy(v, migration, nz*nx*ny*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(v,4L,nx*nz*ny,fpmig);
   laplace_3d_filter(1, nz, nx, ny, v, e);
   fwrite(e,4L,nx*nz*ny,fpmigla);

   cudaMemcpy(Atemp, adcigs, nz*nxa*nya*na*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(Atemp,sizeof(float),nz*nxa*nya*na,fpadcigs);


    end = clock();
/*********IS Loop end*********/ 		     
   printf("---   The forward is over    \n"); 
   printf("---   Complete!!!!!!!!! \n");  
   printf("total %d shots: %f (min)\n", ns, ((float)(end-start))/60.0/CLOCKS_PER_SEC);



/***********close************/ 
          fclose(fpsnap);   fclose(fpshot);   fclose(fpmig);fclose(fpmigla); fclose(fpillum); fclose(fpadcigs);
/***********free*************/ 
       cudaFree(coffx1);       cudaFree(coffx2);
       cudaFree(coffz1);       cudaFree(coffz2);
       cudaFree(acoffx1);      cudaFree(acoffx2);
       cudaFree(acoffz1);      cudaFree(acoffz2);

       cudaFree(s_u0);           cudaFree(s_u1);
       cudaFree(s_v0);           cudaFree(s_v1);
       cudaFree(s_w0);           cudaFree(s_w1);

       cudaFree(s_P);            cudaFree(s_Q);

       cudaFree(s_px0);          cudaFree(s_px1);
       cudaFree(s_py0);          cudaFree(s_py1);
       cudaFree(s_pz0);          cudaFree(s_pz1);
       cudaFree(s_qx0);          cudaFree(s_qx1);
       cudaFree(s_qy0);          cudaFree(s_qy1);
       cudaFree(s_qz0);          cudaFree(s_qz1);

       cudaFree(g_u0);           cudaFree(g_u1);
       cudaFree(g_v0);           cudaFree(g_v1);
       cudaFree(g_w0);           cudaFree(g_w1);

       cudaFree(g_P);            cudaFree(g_Q);

       cudaFree(g_px0);          cudaFree(g_px1);
       cudaFree(g_py0);          cudaFree(g_py1);
       cudaFree(g_pz0);          cudaFree(g_pz1);
       cudaFree(g_qx0);          cudaFree(g_qx1);
       cudaFree(g_qy0);          cudaFree(g_qy1);
       cudaFree(g_qz0);          cudaFree(g_qz1);

       cudaFree(s_P_bndr);
       cudaFree(s_Q_bndr);

       cudaFree(shot_Dev);

       cudaFree(migration);
       cudaFree(illumination);
       cudaFree(adcigs);
/***************host free*****************/
	free(v);	free(e);	free(d);
       free(shot_Hos);  free(Atemp);
}

