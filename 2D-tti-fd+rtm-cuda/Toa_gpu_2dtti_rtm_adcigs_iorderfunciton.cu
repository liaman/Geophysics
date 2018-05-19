//a#########################################################
//a##         2D Acoustic TTI Medium RTM   
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
//a##                     vpx^2=vp^2*(1+2*epsilu);
//a##                     vpn^2=vp^2*(1+2*deta);
//a##*********a*******************/
//a##
//a##                                  code by Rong Tao 
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
/*************func**************/
void laplace_filter(int adj, int nz, int nx, float *in, float *out)
/*< linear operator, come from Madagascar Mlaplac2>*/
{
    int iz,ix,j;
    for (j=0;j<nx*nz;j++) out[j]=0.0;
    for (ix=0; ix < nx; ix++) {
	for (iz=0; iz < nz; iz++) {
	    j = iz+ix*nz;
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
	}
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
                           float *deta,float *epsilu,float *theta,int fs,int ds,int zs,int is,bool SV)
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
   d0=3.0*vp[nnx*nnz/2]*log(100000.0)/(2.0*npml*((dx+dz)/2.0));
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
               float *vv,float *epsilu,float *deta,float *theta,int npml)
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
				 fread(&vv[id],4L,1,fp1);
				 fread(&epsilu[id],4L,1,fp2);
				 fread(&deta[id],4L,1,fp3);
				 fread(&theta[id],4L,1,fp4);theta[id]*=pi/180.0;
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
			 coff1[id]=1.0/(1.0+(dt*d0*pow((npml-id)*1.0/npml,2.0))/2.0);
			 coff2[id]=coff1[id]*(1.0-(dt*d0*pow((npml-id)*1.0/npml,2.0))/2.0);

			 acoff1[id]=1.0/(1.0+(dt*d0*pow(((npml-id)*1.0)/npml,2.0))/2.0);
			 acoff2[id]=acoff1[id]*(1.0-(dt*d0*pow(((npml-id)*1.0)/npml,2.0))/2.0);

		 }else if(id>=npml&&id<npml+nn){

			 coff1[id]=1.0;
			 coff2[id]=1.0;

			 acoff1[id]=1.0;
			 acoff2[id]=1.0;

		 }else{

			 coff1[id]=1.0/(1.0+(dt*d0*pow((id-nn-npml)*1.0/npml,2.0))/2.0);
			 coff2[id]=coff1[id]*(1.0-(dt*d0*pow((id-nn-npml)*1.0/npml,2.0))/2.0);

			 acoff1[id]=1.0/(1.0+(dt*d0*pow(((id-nn-npml)*1.0)/npml,2.0))/2.0);
			 acoff2[id]=acoff1[id]*(1.0-(dt*d0*pow(((id-nn-npml)*1.0)/npml,2.0))/2.0);
		 }	
            }      
}
/*************func*******************/
__global__ void shot_record(int nnx, int nnz, int nx, int nz, int npml, int it, int nt, float *P, float *shot, bool flag)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<nx)
            {
             if(flag){
               shot[it+nt*id]=P[npml+nnz*(id+npml)];
             }else{
               P[npml+nnz*(id+npml)]=shot[it+nt*id];
              }
            }       
}
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
}/*************func*******************/
__global__ void wavefield_bndr(int nnx, int nnz, int nx, int nz, int npml, int it, int nt, 
                               float *P, float *Q, float *P_bndr, float *Q_bndr, bool flag)
{		
	 int id=threadIdx.x+blockDim.x*blockIdx.x;

           if(id<2*nx+2*nz)
            {
            if(flag)/////////////////////////////////save boundary
             {
              if(id<nx){//up

               P_bndr[it*(2*nx+2*nz)+id]=P[npml-1+nnz*(id+npml)];
               Q_bndr[it*(2*nx+2*nz)+id]=Q[npml-1+nnz*(id+npml)];

              }else if(id>=nx&&id<(2*nx)){//down
   
               P_bndr[it*(2*nx+2*nz)+id]=P[npml+nz+1+nnz*(id-nx+npml)];
               Q_bndr[it*(2*nx+2*nz)+id]=Q[npml+nz+1+nnz*(id-nx+npml)];


              }else if(id>=(2*nx)&&id<(2*nx+nz)){//left

               P_bndr[it*(2*nx+2*nz)+id]=P[id-2*nx+npml+nnz*(npml-1)];
               Q_bndr[it*(2*nx+2*nz)+id]=Q[id-2*nx+npml+nnz*(npml-1)];

              }else if(id>=(2*nx+nz)){//right

               P_bndr[it*(2*nx+2*nz)+id]=P[id-2*nx-nz+npml+nnz*(npml+nx+1)];
               Q_bndr[it*(2*nx+2*nz)+id]=Q[id-2*nx-nz+npml+nnz*(npml+nx+1)];

                }
            }else{/////////////////////////////add boundary
              if(id<nx){//up

               P[npml-1+nnz*(id+npml)]=P_bndr[it*(2*nx+2*nz)+id];
               Q[npml-1+nnz*(id+npml)]=Q_bndr[it*(2*nx+2*nz)+id];

              }else if(id>=nx&&id<(2*nx)){//down
   
               P[npml+nz+1+nnz*(id-nx+npml)]=P_bndr[it*(2*nx+2*nz)+id];
               Q[npml+nz+1+nnz*(id-nx+npml)]=Q_bndr[it*(2*nx+2*nz)+id];


              }else if(id>=(2*nx)&&id<(2*nx+nz)){//left

               P[id-2*nx+npml+nnz*(npml-1)]=P_bndr[it*(2*nx+2*nz)+id];
               Q[id-2*nx+npml+nnz*(npml-1)]=Q_bndr[it*(2*nx+2*nz)+id];

              }else if(id>=(2*nx+nz)){//right

               P[id-2*nx-nz+npml+nnz*(npml+nx+1)]=P_bndr[it*(2*nx+2*nz)+id];
               Q[id-2*nx-nz+npml+nnz*(npml+nx+1)]=Q_bndr[it*(2*nx+2*nz)+id];

                }
             }
            }       
}
/*************func**************/    
__global__ void cal_migration(int nnx, int nnz, int nz, int npml, float *migration, float *s, float *g)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix=id/nz;
    int iz=id%nz;

   if(id<nnx*nnz)
   {
      migration[id]+=s[iz+npml+nnz*(ix+npml)]*g[iz+npml+nnz*(ix+npml)];
   }
}
/*************func**************/    
__global__ void cal_illumination(int nnx, int nnz, int nz, int npml, float *illumination, float *P, float *Q)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix=id/nz;
    int iz=id%nz;

   if(id<nnx*nnz)
   {
      illumination[id]+=P[iz+npml+nnz*(ix+npml)]*P[iz+npml+nnz*(ix+npml)]
                         +Q[iz+npml+nnz*(ix+npml)]*Q[iz+npml+nnz*(ix+npml)];
      if(illumination[id]==0)illumination[id]=1.0;
   }
}
/*************func**************/    
__global__ void migration_illum(int nx, int nz, int npml, float *migration, float *illumination)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

   if(id<nx*nz)
   {
      migration[id]/=illumination[id];
   }
}
/*************func**************/    
__global__ void Poynting_Adcigs(int nnz, int nx, int nz, int npml, int na, int da, int dcdp, float *adcigs, 
                           float *s_P, float *s_Q, float *s_u, float *s_w, 
                           float *g_P, float *g_Q, float *g_u, float *g_w)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix=id/nz*dcdp;
    int iz=id%nz;

    int ia=0;

    float Ssx=-s_P[iz+npml+nnz*(ix+npml)]*s_u[iz+npml+nnz*(ix+npml)];
    float Ssz=-s_Q[iz+npml+nnz*(ix+npml)]*s_w[iz+npml+nnz*(ix+npml)];
    float Sgx= g_P[iz+npml+nnz*(ix+npml)]*g_u[iz+npml+nnz*(ix+npml)];
    float Sgz= g_Q[iz+npml+nnz*(ix+npml)]*g_w[iz+npml+nnz*(ix+npml)];

    float b1= Ssx*Ssx + Ssz*Ssz;
    float b2= Sgx*Sgx + Sgz*Sgz;
    float  a=(Ssx*Sgx + Ssz*Sgz)/(sqrtf(b1*b2)*(1 - 0.1));

   if(id<nx/dcdp*nz)
   {
     if(a>=-1&&a<=1)
      {
         a=0.5*acosf(a)*180.0/pi;
         ia=(int)(a/(da*1.0));
         if(ia<na)
          {
             adcigs[iz+nz*ia+nz*na*(id/nz)] += s_P[iz+npml+nnz*(ix+npml)]*g_P[iz+npml+nnz*(ix+npml)]
                                                *cosf(ia*pi/180.0)*cosf(ia*pi/180.0)*cosf(ia*pi/180.0);
          }
      }
   }
}
/*************func**************/    
__global__ void adcigs_illum(int nx, int nz, int na, int da, int dcdp, float *adcigs, float *illumination)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix=id/(nz*na)*dcdp;
    int iz=id%nz;

   if(id<nx*nz/dcdp*na)
   {
      adcigs[id]/=illumination[iz+nz*ix];
   }
}
//a########################################################################
int main(int argc,char *argv[])
{

	int is, it, nx, nz, nnx, nnz, nt, wtype, na, da, dcdp, nxa;
	int ns, ds, fs, zs, npml;
	float dx, dz, dt, t, pfac, favg;
       float *coffx1,*coffx2,*coffz1,*coffz2,*acoffx1,*acoffx2,*acoffz1,*acoffz2;
	float *v, *e, *d, *b;
	float *vp, *epsilu, *deta, *theta;
	float *s_u0, *s_u1, *s_px0, *s_qx0, *s_px1, *s_qx1;
       float *s_w0, *s_w1, *s_pz0, *s_qz0, *s_pz1, *s_qz1;
	float *s_P, *s_Q;
	float *g_u0, *g_u1, *g_px0, *g_qx0, *g_px1, *g_qx1;
       float *g_w0, *g_w1, *g_pz0, *g_qz0, *g_pz1, *g_qz1;
	float *g_P, *g_Q;
       float *shot_Dev, *shot_Hos, *s_P_bndr, *s_Q_bndr;

       float *migration, *illumination, *adcigs;
       float *Atemp;

       clock_t start, end, is_t0, is_t1;
/*************wavelet\boundary**************/
          wtype=1;npml=50;
/********** dat document ***********/
      /*    char FN1[250]={"waxian_vel_601_301.dat"};
          char FN2[250]={"waxian_epsilon_601_301.dat"};
          char FN3[250]={"waxian_delta_601_301.dat"};
          char FN4[250]={"waxian_theta_601_301.dat"};
	   char FN5[250]={"waxian_45_shot_cal.dat"};
	   char FN6[250]={"waxian_45_snap.dat"};
	   char FN7[250]={"waxian_45_migration.dat"};
	   char FN8[250]={"waxian_45_illumination.dat"};
	   char FN9[250]={"waxian_45_adcigs.dat"};*/

     /*     char FN1[250]={"layers_vel_601_301.dat"};
          char FN2[250]={"layers_epsilon_601_301.dat"};
          char FN3[250]={"layers_delta_601_301.dat"};
          char FN4[250]={"layers_theta_601_301.dat"};
	   char FN5[250]={"layers_shot_cal.dat"};
	   char FN6[250]={"layers_snap.dat"};
	   char FN7[250]={"layers_migration.dat"};
	   char FN8[250]={"layers_illumination.dat"};
	   char FN9[250]={"layers_adcigs65.dat"}; */

     /*     char FN1[250]={"thrust_vel_711_300.dat"};
          char FN2[250]={"thrust_epsilon_711_300.dat"};
          char FN3[250]={"thrust_delta_711_300.dat"};
          char FN4[250]={"thrust_theta_711_300_theta10-20-30.dat"};
	   char FN5[250]={"thrust_shot_cal.dat"};
	   char FN6[250]={"thrust_snap.dat"};
	   char FN7[250]={"thrust_migration.dat"};
	   char FN8[250]={"thrust_illumination.dat"};
	   char FN9[250]={"thrust_adcigs.dat"};*/

    /*      char FN1[250]={"hesh_vel_883_301.dat"};
          char FN2[250]={"hesh_epsilon_883_301.dat"};
          char FN3[250]={"hesh_delta_883_301.dat"};
          char FN4[250]={"hesh_theta_883_301_0.75.dat"};
	   char FN5[250]={"hesh_shot_cal.dat"};
	   char FN6[250]={"hesh_snap.dat"};
	   char FN7[250]={"hesh_migration.dat"};
	   char FN8[250]={"hesh_illumination.dat"};
	   char FN9[250]={"hesh_adcigs.dat"}; */

   /*       char FN1[250]={"yongan_vel_650_301.dat"};
          char FN2[250]={"epsilon_const0.3.dat"};
          char FN3[250]={"delta_const0.2.dat"};
          char FN4[250]={"theta_const30.dat"};
	   char FN5[250]={"yongan_shot_cal.dat"};
	   char FN6[250]={"yongan_snap.dat"};
	   char FN7[250]={"yongan_migration.dat"};
	   char FN8[250]={"yongan_illumination.dat"};
	   char FN9[250]={"yongan_adcigs.dat"};*/

          char FN1[250]={"BP_vel.dat"};
          char FN2[250]={"BP_epsilon.dat"};
          char FN3[250]={"BP_delta.dat"};
          char FN4[250]={"BP_theta.dat"};
	   char FN5[250]={"BP_shot_cal.dat"};
	   char FN6[250]={"BP_snap.dat"};
	   char FN7[250]={"BP_migration.dat"};
	   char FN8[250]={"BP_illumination.dat"};
	   char FN9[250]={"BP_adcigs.dat"};
/********aaa************/  
	 FILE *fpsnap, *fpshot, *fpmig, *fpillum, *fpadcigs;
        fpshot=fopen(FN5,"wb");
        fpsnap=fopen(FN6,"wb");
        fpmig=fopen(FN7,"wb");
        fpillum=fopen(FN8,"wb");
        fpadcigs=fopen(FN9,"wb");

 
/********* parameters *************/

          nx=1700;              
	   nz=601;         favg=20;     pfac=10.0;

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=9001;    
          dt=0.0005;
     
          ns=340;       
          fs=nx/ns/2;      
          ds=nx/ns;
          zs=1; 
  
          na=65; 
          da=1;  
          dcdp=1;
/*************v***************/ 
          nxa=(int)(nx/dcdp);
          nnx=nx+2*npml;
          nnz=nz+2*npml;
/************a*************/
    	 Atemp=(float*)malloc(nz*nx/dcdp*na*sizeof(float));

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
        cudaMalloc(&epsilu, nnz*nnx*sizeof(float));
        cudaMalloc(&deta, nnz*nnx*sizeof(float));
        cudaMalloc(&theta, nnz*nnx*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilu, e, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(deta, d, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(theta, b, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
/****************************/
        cudaMalloc(&s_u0, nnz*nnx*sizeof(float));    cudaMalloc(&s_u1, nnz*nnx*sizeof(float));
        cudaMalloc(&s_w0, nnz*nnx*sizeof(float));    cudaMalloc(&s_w1, nnz*nnx*sizeof(float));

        cudaMalloc(&s_P, nnz*nnx*sizeof(float));     cudaMalloc(&s_Q, nnz*nnx*sizeof(float));

        cudaMalloc(&s_px0, nnz*nnx*sizeof(float));   cudaMalloc(&s_px1, nnz*nnx*sizeof(float));
        cudaMalloc(&s_pz0, nnz*nnx*sizeof(float));   cudaMalloc(&s_pz1, nnz*nnx*sizeof(float));
        cudaMalloc(&s_qx0, nnz*nnx*sizeof(float));   cudaMalloc(&s_qx1, nnz*nnx*sizeof(float));
        cudaMalloc(&s_qz0, nnz*nnx*sizeof(float));   cudaMalloc(&s_qz1, nnz*nnx*sizeof(float));

        cudaMalloc(&g_u0, nnz*nnx*sizeof(float));    cudaMalloc(&g_u1, nnz*nnx*sizeof(float));
        cudaMalloc(&g_w0, nnz*nnx*sizeof(float));    cudaMalloc(&g_w1, nnz*nnx*sizeof(float));

        cudaMalloc(&g_P, nnz*nnx*sizeof(float));     cudaMalloc(&g_Q, nnz*nnx*sizeof(float));

        cudaMalloc(&g_px0, nnz*nnx*sizeof(float));   cudaMalloc(&g_px1, nnz*nnx*sizeof(float));
        cudaMalloc(&g_pz0, nnz*nnx*sizeof(float));   cudaMalloc(&g_pz1, nnz*nnx*sizeof(float));
        cudaMalloc(&g_qx0, nnz*nnx*sizeof(float));   cudaMalloc(&g_qx1, nnz*nnx*sizeof(float));
        cudaMalloc(&g_qz0, nnz*nnx*sizeof(float));   cudaMalloc(&g_qz1, nnz*nnx*sizeof(float));

        cudaMalloc(&coffx1, nnx*sizeof(float));     cudaMalloc(&coffx2, nnx*sizeof(float));
        cudaMalloc(&coffz1, nnz*sizeof(float));     cudaMalloc(&coffz2, nnz*sizeof(float));
        cudaMalloc(&acoffx1, nnx*sizeof(float));    cudaMalloc(&acoffx2, nnx*sizeof(float));
        cudaMalloc(&acoffz1, nnz*sizeof(float));    cudaMalloc(&acoffz2, nnz*sizeof(float));

        cudaMalloc(&shot_Dev, nx*nt*sizeof(float));
        cudaMalloc(&s_P_bndr, nt*(2*nx+2*nz)*sizeof(float));
        cudaMalloc(&s_Q_bndr, nt*(2*nx+2*nz)*sizeof(float));

        cudaMalloc(&migration, nz*nx*sizeof(float)); 
        cudaMalloc(&illumination, nz*nx*sizeof(float));
        cudaMalloc(&adcigs, nz*na*nxa*sizeof(float));
/******************************/
	 check_gpu_error("Failed to allocate memory for variables!");

        get_d0<<<1, 1>>>(dx, dz, nnx, nnz, npml, vp);
        initial_coffe<<<(nnx+511)/512, 512>>>(dt,nx,coffx1,coffx2,acoffx1,acoffx2,npml);
        initial_coffe<<<(nnz+511)/512, 512>>>(dt,nz,coffz1,coffz2,acoffz1,acoffz2,npml);

        cudaMemset(migration, 0, nz*nx*sizeof(float)); 
        cudaMemset(illumination, 0, nz*nx*sizeof(float));
        cudaMemset(adcigs, 0, nz*na*nxa*sizeof(float));

        printf("--------------------------------------------------------\n");
        printf("---   \n");   
        start = clock();                                  
/**********IS Loop start*******/
   for(is=1;is<=ns;is++)	
    {     
         printf("---   IS=%3d,  ",is);
        is_t0 = clock(); 
     cudaMemset(s_u0, 0, nnz*nnx*sizeof(float));     cudaMemset(s_u1, 0, nnz*nnx*sizeof(float));
     cudaMemset(s_w0, 0, nnz*nnx*sizeof(float));     cudaMemset(s_w1, 0, nnz*nnx*sizeof(float));

     cudaMemset(s_P, 0, nnz*nnx*sizeof(float));      cudaMemset(s_Q, 0, nnz*nnx*sizeof(float));

     cudaMemset(s_px0, 0, nnz*nnx*sizeof(float));    cudaMemset(s_px1, 0, nnz*nnx*sizeof(float));
     cudaMemset(s_pz0, 0, nnz*nnx*sizeof(float));    cudaMemset(s_pz1, 0, nnz*nnx*sizeof(float));
     cudaMemset(s_qx0, 0, nnz*nnx*sizeof(float));    cudaMemset(s_qx1, 0, nnz*nnx*sizeof(float));
     cudaMemset(s_qz0, 0, nnz*nnx*sizeof(float));    cudaMemset(s_qz1, 0, nnz*nnx*sizeof(float));

     cudaMemset(g_u0, 0, nnz*nnx*sizeof(float));     cudaMemset(g_u1, 0, nnz*nnx*sizeof(float));
     cudaMemset(g_w0, 0, nnz*nnx*sizeof(float));     cudaMemset(g_w1, 0, nnz*nnx*sizeof(float));

     cudaMemset(g_P, 0, nnz*nnx*sizeof(float));      cudaMemset(g_Q, 0, nnz*nnx*sizeof(float));

     cudaMemset(g_px0, 0, nnz*nnx*sizeof(float));    cudaMemset(g_px1, 0, nnz*nnx*sizeof(float));
     cudaMemset(g_pz0, 0, nnz*nnx*sizeof(float));    cudaMemset(g_pz1, 0, nnz*nnx*sizeof(float));
     cudaMemset(g_qx0, 0, nnz*nnx*sizeof(float));    cudaMemset(g_qx1, 0, nnz*nnx*sizeof(float));
     cudaMemset(g_qz0, 0, nnz*nnx*sizeof(float));    cudaMemset(g_qz1, 0, nnz*nnx*sizeof(float));

     cudaMemset(shot_Dev, 0, nt*nx*sizeof(float));
     cudaMemset(s_P_bndr, 0, nt*(2*nx+2*nz)*sizeof(float));
     cudaMemset(s_Q_bndr, 0, nt*(2*nx+2*nz)*sizeof(float));

     for(it=0,t=dt;it<nt;it++,t+=dt)
     { 
     // if(it%100==0&&is==1)printf("---   is===%d   it===%d\n",is,it);
       /*a#####################a*/
       /*a##     Forward     ##a*/
       /*a#####################a*/
	 add_source<<<1,1>>>(pfac,fs,zs,nx,nz,nnx,nnz,dt,t,favg,wtype,npml,is,ds,s_P,s_Q);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,s_u0,s_w0,s_u1,s_w1,s_P,s_Q,coffx1,coffx2,coffz1,coffz2,theta);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,s_u1,s_w1,s_P,s_Q,vp,npml,
                                                  s_px1,s_px0,s_pz1,s_pz0,s_qx1,s_qx0,s_qz1,s_qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,theta,fs,ds,zs,is,true);
        s_u0=s_u1; s_w0=s_w1; s_px0=s_px1; s_pz0=s_pz1; s_qx0=s_qx1; s_qz0=s_qz1; 

        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, s_P, shot_Dev, true);
        wavefield_bndr<<<((2*nx+2*nz)+511)/512,512>>>(nnx, nnz, nx, nz, npml, it, nt, s_P, s_Q, s_P_bndr, s_Q_bndr, true);
        cal_illumination<<<(nx*nz+511)/512, 512>>>(nnx, nnz, nz, npml, illumination, s_P, s_Q);
        /*   if((is==1)&&(it%50==0))
            {
	       cudaMemcpy(e, s_P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fseek(fpsnap,(int)(it/50)*(nnx)*(nnz)*4L,0);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }   */
     }//it loop end

      mute_directwave<<<(nx*nt+511)/512, 512>>>(nx,nt,dt,favg,dx,dz,fs,ds,zs,is,vp,epsilu,theta,shot_Dev,30);
      cudaMemcpy(shot_Hos, shot_Dev, nt*nx*sizeof(float), cudaMemcpyDeviceToHost);
      fseek(fpshot,(is-1)*nt*nx*sizeof(float),0);
      fwrite(shot_Hos,sizeof(float),nt*nx,fpshot);

     for(it=nt-1;it>=0;it--)
     { 
     // if(it%100==0&&is==1)printf("---   is===%d   it===%d\n",is,it);
       /*a#####################a*/
       /*a##  Reconstruction ##a*/
       /*a#####################a*/
        wavefield_bndr<<<((2*nx+2*nz)+511)/512,512>>>(nnx, nnz, nx, nz, npml, it, nt, s_P, s_Q, s_P_bndr, s_Q_bndr, false);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,s_u0,s_w0,s_u1,s_w1,s_P,s_Q,coffx1,coffx2,coffz1,coffz2,theta);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,s_u1,s_w1,s_P,s_Q,vp,npml,
                                                  s_px1,s_px0,s_pz1,s_pz0,s_qx1,s_qx0,s_qz1,s_qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,theta,fs,ds,zs,is,false);
        s_u0=s_u1; s_w0=s_w1; s_px0=s_px1; s_pz0=s_pz1; s_qx0=s_qx1; s_qz0=s_qz1; 

      /*     if((is==1)&&(it%50==0))
            {
	       cudaMemcpy(e, s_P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fseek(fpsnap,(int)(it/50)*(nnx)*(nnz)*4L,0);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }  */

       /*a#####################a*/
       /*a##     Backward    ##a*/
       /*a#####################a*/
        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, g_P, shot_Dev, false);
        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, g_Q, shot_Dev, false);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,g_u0,g_w0,g_u1,g_w1,g_P,g_Q,coffx1,coffx2,coffz1,coffz2,theta);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,g_u1,g_w1,g_P,g_Q,vp,npml,
                                                  g_px1,g_px0,g_pz1,g_pz0,g_qx1,g_qx0,g_qz1,g_qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,deta,epsilu,theta,fs,ds,zs,is,false);
        g_u0=g_u1; g_w0=g_w1; g_px0=g_px1; g_pz0=g_pz1; g_qx0=g_qx1; g_qz0=g_qz1; 

    /*       if((is==1)&&(it%50==0))
            {
	       cudaMemcpy(e, g_P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fseek(fpsnap,(int)(it/50)*(nnx)*(nnz)*4L,0);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }   */

        cal_migration<<<(nx*nz+511)/512, 512>>>(nnx, nnz, nz, npml, migration, s_P, g_P);

        Poynting_Adcigs<<<(nxa*nz+511)/512, 512>>>(nnz, nx, nz, npml, na, da, dcdp, adcigs, 
                                                       s_P, s_Q, s_u0, s_w0, g_P, g_Q, g_u0, g_w0);


       }//it loop end



      is_t1 = clock();    printf("%.2f (min)\n", ns, ((float)(is_t1-is_t0))/60.0/CLOCKS_PER_SEC);
    }//is loop end

   migration_illum<<<(nx*nz+511)/512, 512>>>(nx, nz, npml, migration, illumination);
   adcigs_illum<<<(nxa*nz*na+511)/512, 512>>>(nx, nz, na, da, dcdp, adcigs, illumination);


   cudaMemcpy(e, migration, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
   laplace_filter(1,nz,nx,e,d);
   fwrite(d,sizeof(float),nx*nz,fpmig);

   cudaMemcpy(e, illumination, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(e,sizeof(float),nx*nz,fpillum);

   cudaMemcpy(Atemp, adcigs, nz*nxa*na*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(Atemp,sizeof(float),nz*nxa*na,fpadcigs);

    end = clock();
/*********IS Loop end*********/ 		     
   printf("---   The forward is over    \n"); 
   printf("---   Complete!!!!!!!!! \n");  
   printf("total %d shots: %.2f (min)\n", ns, ((float)(end-start))/60.0/CLOCKS_PER_SEC);



/***********close************/ 
          fclose(fpsnap);   fclose(fpshot);  fclose(fpmig);  fclose(fpillum);  fclose(fpadcigs);
/***********free*************/ 
       cudaFree(coffx1);       cudaFree(coffx2);
       cudaFree(coffz1);       cudaFree(coffz2);
       cudaFree(acoffx1);      cudaFree(acoffx2);
       cudaFree(acoffz1);      cudaFree(acoffz2);

       cudaFree(s_u0);           cudaFree(s_u1);
       cudaFree(s_w0);           cudaFree(s_w1);

       cudaFree(s_P);            cudaFree(s_Q);

       cudaFree(s_px0);          cudaFree(s_px1);
       cudaFree(s_pz0);          cudaFree(s_pz1);
       cudaFree(s_qx0);          cudaFree(s_qx1);
       cudaFree(s_qz0);          cudaFree(s_qz1);

       cudaFree(g_u0);           cudaFree(g_u1);
       cudaFree(g_w0);           cudaFree(g_w1);

       cudaFree(g_P);            cudaFree(g_Q);

       cudaFree(g_px0);          cudaFree(g_px1);
       cudaFree(g_pz0);          cudaFree(g_pz1);
       cudaFree(g_qx0);          cudaFree(g_qx1);
       cudaFree(g_qz0);          cudaFree(g_qz1);

       cudaFree(shot_Dev);
       cudaFree(s_P_bndr);        cudaFree(s_Q_bndr);

       cudaFree(migration); 
       cudaFree(illumination);
       cudaFree(adcigs);

       cudaFree(vp);
       cudaFree(epsilu);
       cudaFree(deta);
       cudaFree(theta);
/***************host free*****************/
	free(v);	free(e);	free(d);free(b);
       free(shot_Hos);   free(Atemp);
}

