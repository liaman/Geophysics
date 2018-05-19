//a#########################################################
//a##         2D Acoustic VTI Medium RTM   
//a##  Ps : P + sv wave and get rid of sv        
//a##       GPU(CUDA) ,poynting adcigs, read shot
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
//a##                        first code: 2017.2.15
//a##                    adcigs process: 2017.4.13
//a##
//a##                                  code by Rong Tao 
//a##                                     
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include <cuda_runtime.h>

#define pi 3.141592653

#define mm 4

//__constant__ float c[mm]={1.125,-0.04166667};/*mm==2*/
//__constant__ float c[mm]={1.1718750,-0.065104167,0.0046875};/*mm==3*/
__constant__ float c[mm]={1.196289,-0.0797526,0.009570313,-0.0006975447};/*mm==4*/
//__constant__ float c[mm]={1.211243,-0.08972168,0.01384277,-0.00176566,0.0001186795};/*mm==5*/


__device__ float d0;
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
/*************func**************/
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
                           float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

	int ix,iz,im;
	float dtx,dtz,xx,zz;

        ix=id/nnz;
        iz=id%nnz;

		 dtx=dt/dx;
		 dtz=dt/dz;
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
__global__ void update_stress(int nx,int nz,int nnx,int nnz,float dt,float dx,float dz,
                           float *u1,float *w1,float *P,float *Q,float *vp,int npml,
                           float *px1,float *px0,float *pz1,float *pz0,float *qx1,float *qx0,float *qz1,float *qz0,
                           float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,
                           float *delta,float *epsilon,int fs,int ds,int zs,int is,bool SV)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

	int im,ix,iz,rx,rz,R=15,r=5;
	float dtx,dtz, xx,zz,ee,dd;

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
     }
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
/*************func*******************/
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

           if((it>(int)(mu_t0/dt)-tt)&&(it<mu_nt))
              shot[id]=0.0;
   }
}/*************func**************/    
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
__global__ void migration_illum(int nx, int nz, int npml, float *migration, float *illumination)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;

   if(id<nx*nz)
   {
      migration[id]/=illumination[id];//*illumination[id];
   }
}
/*************func**************/    
__global__ void Poynting_Adcigs(int nnz, int nx, int nz, int npml, int na, int da, float *adcigs, 
                           float *s_P, float *s_Q, float *s_u, float *s_w, 
                           float *g_P, float *g_Q, float *g_u, float *g_w)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix=id/nz;
    int iz=id%nz;

    int ia=0;

    float Ssx=-s_P[iz+npml+nnz*(ix+npml)]*s_u[iz+npml+nnz*(ix+npml)];
    float Ssz=-s_Q[iz+npml+nnz*(ix+npml)]*s_w[iz+npml+nnz*(ix+npml)];
    float Sgx= g_P[iz+npml+nnz*(ix+npml)]*g_u[iz+npml+nnz*(ix+npml)];
    float Sgz= g_Q[iz+npml+nnz*(ix+npml)]*g_w[iz+npml+nnz*(ix+npml)];

    float b1= Ssx*Ssx + Ssz*Ssz;
    float b2= Sgx*Sgx + Sgz*Sgz;
    float  a=(Ssx*Sgx + Ssz*Sgz)/(sqrtf(b1*b2)*(1 - 0.1));

   if(id<nx*nz)
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
__global__ void adcigs_illum(int nx, int nz, int na, int da, float *adcigs, float *illumination)
{
    int id=threadIdx.x+blockDim.x*blockIdx.x;
    int ix=id/(nz*na);
    int iz=id%nz;

   if(id<nx*nz*na)
   {
      adcigs[id]/=illumination[iz+nz*ix];//*illumination[iz+nz*ix];
   }
}
/*************func**************/ 
void stk_adcigs(int nx,int nz,int na,float *adcigs,float *migration)
{
   int ix,iz,ia,id,ido;
   float stk;
   float *temp;

   temp=(float*)malloc(nz*nx*sizeof(float));

   for (ix=0; ix<nx; ix++)  {
       for (iz=0; iz<nz; iz++)  {
           stk=0.0;
           for (ia=0; ia<na; ia++)  {
                id=ix*na*nz+ia*nz+iz;
                stk+=adcigs[id];
             }
           ido=ix*nz+iz;
           temp[ido]=stk;
        }
   }
   laplace_filter(1,nz,nx,temp,migration);
}
/*************func**************/ 
void adcigs_smiled(int nx,int nz,int na,int dcdp,float *adcigs)   
{
   int ix,iz,ia,id,ido;
   float *temp;

   temp=(float*)malloc(nz*nx/dcdp*na*sizeof(float));
   for (ix=0; ix<nx; ix++)  {
       for (ia=0; ia<na; ia++)  {
           for (iz=0; iz<nz; iz++)  {
                id=ix*na*nz+ia*nz+iz;
                if(ix%dcdp==0) {
                      ido=ix/dcdp*na*nz+ia*nz+iz;
                      temp[ido]=adcigs[id];
                      adcigs[ido]=temp[ido];
                  }
             }
        }
   }

}
/*************func*******************/
void read_file(char FN1[],char FN2[],char FN3[],int nx,int nz,int nnx,int nnz,float dx,float dz,float favg,float dt,
               float *v,float *e,float *d,int npml)
{
		 int i,j,id;
               float vmax,  vmin,emax, emin, dmax, dmin,    H_min, dt_max, dxz_max, C, tmp;

		
		 FILE *fp1,*fp2,*fp3;
		 if((fp1=fopen(FN1,"rb"))==NULL){printf("error open <%s>!\n",FN1);exit(0);}
		 if((fp2=fopen(FN2,"rb"))==NULL){printf("error open <%s>!\n",FN2);exit(0);}
		 if((fp3=fopen(FN3,"rb"))==NULL){printf("error open <%s>!\n",FN3);exit(0);}

               vmin= 999999.9;
               vmax=-999999.9;
		 for(i=npml;i<nx+npml;i++)
		 {
			 for(j=npml;j<nz+npml;j++)
			 {
                            id=i*nnz+j;
				 fread(&v[id],4L,1,fp1);/* test the paras *///v[id]*=1.00;
				 fread(&e[id],4L,1,fp2);/* test the paras *///e[id]*=1.00;
				 fread(&d[id],4L,1,fp3);/* test the paras *///d[id]*=0.00;

                            if(vmax<v[id]) vmax = v[id];
			       if(vmin>v[id]) vmin = v[id];
                            if(emax<e[id]) emax = e[id];
			       if(emin>e[id]) emin = e[id];
                            if(dmax<d[id]) dmax = d[id];
			       if(dmin>d[id]) dmin = d[id];
			 }
		 }
		 fclose(fp1);
		 fclose(fp2);
		 fclose(fp3);
        printf("------------------------------------\n---\n");
        printf("---   Vmax=%.2f, Vmin=%.2f\n",vmax,vmin);
        printf("---   Emax=%.4f, Emin=%.4f\n",emax,emin);
        printf("---   Dmax=%.4f, Dmin=%.4f\n---\n",dmax,dmin);
        /*********boundary*********/
        pad_vv(nx,nz,nnx,nnz,npml,e);
        pad_vv(nx,nz,nnx,nnz,npml,d);
        pad_vv(nx,nz,nnx,nnz,npml,v); 
      
       H_min=dx<dz?dx:dz;
       dt_max = 0.5*H_min/vmin;
       dxz_max = vmax/favg*0.2;

       if(dxz_max<dz||dxz_max<dx){printf("---   You need have to redefine DX and DZ ! \n");exit(0);}
	if(dt_max<dt){printf("---   You need have to redefine DT ! \n");exit(0);}
       if (   favg >= vmin/( 5.0*(dx>dz?dx:dz) )   ||   favg >= vmin/( 5.0*(dx>dz?dx:dz) )  )
	             {printf("---   Non-dispersion relation not satisfied! \n");exit(0);}

	else if ( mm == 2 )     C = 0.857;
	else if ( mm == 3 )     C = 0.8;
	else if ( mm == 4 )     C = 0.777;
	else if ( mm == 5 )     C = 0.759;

       tmp = dt*vmax*sqrtf( 1.0/(dx*dx)+1.0/(dz*dz) );
       if ( tmp >= C){ printf("---   Stability condition not satisfied! tmp = %f, C = %f\n",tmp,C);exit(0);}
}
//a########################################################################
//a##                          Main Function                             ##
//a########################################################################
int main(int argc,char *argv[])
{
	int is, it, nx, nz, nnx, nnz, nt, wtype, na, da, dcdp;
	int ns, ds, fs, zs, npml;
	float dx, dz, dt, t, pfac, favg;
       float *coffx1,*coffx2,*coffz1,*coffz2,*acoffx1,*acoffx2,*acoffz1,*acoffz2;
	float *v, *e, *d;
	float *vp, *epsilon, *delta;
	float *s_u0, *s_u1, *s_px0, *s_qx0, *s_px1, *s_qx1;
       float *s_w0, *s_w1, *s_pz0, *s_qz0, *s_pz1, *s_qz1;
	float *g_u0, *g_u1, *g_px0, *g_qx0, *g_px1, *g_qx1;
       float *g_w0, *g_w1, *g_pz0, *g_qz0, *g_pz1, *g_qz1;
	float *s_P, *s_Q, *g_P, *g_Q, *shot_Dev, *shot_Hos, *P_bndr, *Q_bndr;
       float *migration, *illumination, *adcigs;
       float *Atemp;
       bool read;


       clock_t start, end;
/*************wavelet\boundary**************/
          wtype=1;npml=20;
/********** dat document ***********/

          char FN1[250]={"waxian_3layer/waxian_vel_1001_301.dat"};
          char FN2[250]={"waxian_3layer/waxian_epsilon_1001_301.dat"};
          char FN3[250]={"waxian_3layer/waxian_delta_1001_301.dat"};
	   char FN4[250]={"waxian_3layer/waxian_shot_obs.dat"};//shot obs
	   char FN5[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_shot_cal.dat"};//shot cal
	   char FN6[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_snap.dat"};//snap
	   char FN7[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_migration.dat"};
	   char FN8[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_illumination.dat"};
	   char FN9[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_adcigs.dat"}; 
	  char FN10[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_stkadcigs.dat"};
	  char FN11[250]={"waxian_3layer_2/waxian_v1.00e1.00d0.00_smiled_adcigs.dat"};

/********* parameters *************/
          read=true;/* true: read shot; flase: use right shot record */
/********* parameters *************/
          nx=1001;              
	   nz=301;         favg=60;     pfac=1000.0;

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=2501;    
          dt=0.0005;
     
          ns=500;       
          fs=nx/ns/2;      
          ds=nx/ns;
          zs=1;   

          na=70; 
          da=1;  
          dcdp=25;

     /*     char FN1[250]={"sshengli_duankuai/sshengli_vel_1000_550.dat"};      
          char FN2[250]={"sshengli_duankuai/sshengli_epsilon_1000_550.dat"};
          char FN3[250]={"sshengli_duankuai/sshengli_delta_1000_550.dat"}; 
	   char FN4[250]={"sshengli_duankuai/sshengli_shot_obs.dat"};//shot obs
	   char FN5[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_shot_cal.dat"};//shot cal
	   char FN6[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_snap.dat"};//snap
	   char FN7[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_migration.dat"};
	   char FN8[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_illumination.dat"};
	   char FN9[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_adcigs.dat"}; 
	  char FN10[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_stkadcigs.dat"};
	  char FN11[250]={"sshengli_duankuai/sshengli_v1.00e0.60d0.60_smiled_adcigs.dat"};*/

/********* parameters *************/
        //  read=true;/* true: read shot; flase: use right shot record */
/********* parameters *************/
     /*     nx=1000;              
	   nz=550;         favg=40;     pfac=1000.0;

 	   dx=5.0;   
          dz=5.0;   
     
	   nt=4501;    
          dt=0.0005;
     
          ns=500;       
          fs=nx/ns/2;      
          ds=nx/ns;
          zs=1;   

          na=70; 
          da=1;  
          dcdp=25;  */

/********aaa************/  
	 FILE *fpsnap, *fpobs, *fpcal, *fpmig, *fpillum, *fpadcigs, *fpadcigs2, *fpstk;
        if((fpobs=fopen(FN4,"rb"))==NULL){printf("error open <%s>!\n",FN4);exit(0);}
        fpcal=fopen(FN5,"wb");
        fpsnap=fopen(FN6,"wb");
        fpmig=fopen(FN7,"wb");
        fpillum=fopen(FN8,"wb");
        fpadcigs=fopen(FN9,"wb");
        fpstk=fopen(FN10,"wb");
        fpadcigs2=fopen(FN11,"wb");
/*************v***************/ 
          nnx=nx+2*npml;
          nnz=nz+2*npml;
/************a*************/
    	 Atemp=(float*)malloc(nz*nx*na*sizeof(float));

    	 v=(float*)malloc(nnz*nnx*sizeof(float));
    	 e=(float*)malloc(nnz*nnx*sizeof(float));
    	 d=(float*)malloc(nnz*nnx*sizeof(float));
    	 shot_Hos=(float*)malloc(nt*nx*sizeof(float));
        read_file(FN1,FN2,FN3,nx,nz,nnx,nnz,dx,dz,favg,dt,v,e,d,npml);
/****************************/

        cudaSetDevice(0);// initialize device, default device=0;
	 check_gpu_error("Failed to initialize device!");

/****************************/
        cudaMalloc(&vp, nnz*nnx*sizeof(float));
        cudaMalloc(&epsilon, nnz*nnx*sizeof(float));
        cudaMalloc(&delta, nnz*nnx*sizeof(float));
	 cudaMemcpy(vp, v, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(epsilon, e, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	 cudaMemcpy(delta, d, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
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
        cudaMalloc(&P_bndr, nt*(2*nx+2*nz)*sizeof(float));
        cudaMalloc(&Q_bndr, nt*(2*nx+2*nz)*sizeof(float));

        cudaMalloc(&migration, nz*nx*sizeof(float)); 
        cudaMalloc(&illumination, nz*nx*sizeof(float));
        cudaMalloc(&adcigs, nz*na*nx*sizeof(float));
/******************************/
	 check_gpu_error("Failed to allocate memory for variables!");

        get_d0<<<1, 1>>>(dx, dz, nnx, nnz, npml, vp);
        initial_coffe<<<(nnx+511)/512, 512>>>(dt,nx,coffx1,coffx2,acoffx1,acoffx2,npml);
        initial_coffe<<<(nnz+511)/512, 512>>>(dt,nz,coffz1,coffz2,acoffz1,acoffz2,npml);

        cudaMemset(migration, 0, nz*nx*sizeof(float)); 
        cudaMemset(illumination, 0, nz*nx*sizeof(float));
        cudaMemset(adcigs, 0, nz*na*nx*sizeof(float));
        printf("--------------------------------------------------------\n");
        printf("---");   
        start = clock();                                  
/**********IS Loop start*******/
   for(is=1;is<=ns;is++)	
    {     
         printf("\n---   IS=%3d  ",is);

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
     cudaMemset(P_bndr, 0, nt*(2*nx+2*nz)*sizeof(float));
     cudaMemset(Q_bndr, 0, nt*(2*nx+2*nz)*sizeof(float));

/*a***********************************Forward*******************************************/
     for(it=0,t=dt;it<nt;it++,t+=dt)
     { 
      //if(it==0)printf(" > F >",is,it);
       /*a#####################a*/
       /*a##     Forward     ##a*/
       /*a#####################a*/
	 add_source<<<1,1>>>(pfac,fs,zs,nx,nz,nnx,nnz,dt,t,favg,wtype,npml,is,ds,s_P,s_Q);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,
                                               s_u0,s_w0,s_u1,s_w1,s_P,s_Q,coffx1,coffx2,coffz1,coffz2);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,s_u1,s_w1,s_P,s_Q,vp,npml,
                                                  s_px1,s_px0,s_pz1,s_pz0,s_qx1,s_qx0,s_qz1,s_qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,delta,epsilon,fs,ds,zs,is,true);
        s_u0=s_u1; s_w0=s_w1; s_px0=s_px1; s_pz0=s_pz1; s_qx0=s_qx1; s_qz0=s_qz1; 

        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, s_P, shot_Dev, true);
        wavefield_bndr<<<((2*nx+2*nz)+511)/512,512>>>(nnx, nnz, nx, nz, npml, it, nt, s_P, s_Q, P_bndr, Q_bndr, true);
        cal_illumination<<<(nx*nz+511)/512, 512>>>(nnx, nnz, nz, npml, illumination, s_P, s_Q);

     /*      if((is==1)&&(it%300==0))
            {
	       cudaMemcpy(e, s_P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }*/
     }//it loop end
      mute_directwave<<<(nx*nt+511)/512, 512>>>(nx,nt,dt,favg,dx,dz,fs,ds,zs,is,vp,epsilon,shot_Dev,20);
      cudaMemcpy(shot_Hos, shot_Dev, nt*nx*sizeof(float), cudaMemcpyDeviceToHost);
      //fseek(fpcal,(is-1)*nt*nx*sizeof(float),0);
      //fwrite(shot_Hos,sizeof(float),nt*nx,fpcal);

    if(read){
         fseek(fpobs,(is-1)*nt*nx*sizeof(float),0);
         fread(shot_Hos,sizeof(float),nt*nx,fpobs);
         cudaMemcpy(shot_Dev, shot_Hos, nt*nx*sizeof(float), cudaMemcpyHostToDevice);
    }
/*a***********************************Backward*******************************************/
     for(it=nt-1;it>=0;it--)
     { 
     // if(it==0)printf("  B ",is,it);
       /*a#####################a*/
       /*a##  Reconstruction ##a*/
       /*a#####################a*/
        wavefield_bndr<<<((2*nx+2*nz)+511)/512,512>>>(nnx, nnz, nx, nz, npml, it, nt, s_P, s_Q, P_bndr, Q_bndr, false);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,
                                               s_u0,s_w0,s_u1,s_w1,s_P,s_Q,coffx1,coffx2,coffz1,coffz2);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,s_u1,s_w1,s_P,s_Q,vp,npml,
                                                  s_px1,s_px0,s_pz1,s_pz0,s_qx1,s_qx0,s_qz1,s_qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,delta,epsilon,fs,ds,zs,is,false);
        s_u0=s_u1; s_w0=s_w1; s_px0=s_px1; s_pz0=s_pz1; s_qx0=s_qx1; s_qz0=s_qz1; 

         /*  if((is==1)&&(it%300==0))
            {
	       cudaMemcpy(e, s_P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }*/
       /*a#####################a*/
       /*a##     Backward    ##a*/
       /*a#####################a*/
        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, g_P, shot_Dev, false);
        shot_record<<<(nx+511)/512, 512>>>(nnx, nnz, nx, nz, npml, it, nt, g_Q, shot_Dev, false);
        update_vel<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,npml,dt,dx,dz,
                                               g_u0,g_w0,g_u1,g_w1,g_P,g_Q,coffx1,coffx2,coffz1,coffz2);
        update_stress<<<(nnx*nnz+511)/512, 512>>>(nx,nz,nnx,nnz,dt,dx,dz,g_u1,g_w1,g_P,g_Q,vp,npml,
                                                  g_px1,g_px0,g_pz1,g_pz0,g_qx1,g_qx0,g_qz1,g_qz0,
                                                  acoffx1,acoffx2,acoffz1,acoffz2,delta,epsilon,fs,ds,zs,is,false);
        g_u0=g_u1; g_w0=g_w1; g_px0=g_px1; g_pz0=g_pz1; g_qx0=g_qx1; g_qz0=g_qz1; 

        /*   if((is==1)&&(it%300==0))
            {
	       cudaMemcpy(e, g_P, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
              fwrite(e,4L,nnx*nnz,fpsnap);
            }*/
        cal_migration<<<(nx*nz+511)/512, 512>>>(nnx, nnz, nz, npml, migration, s_P, g_P);

        Poynting_Adcigs<<<(nx*nz+511)/512, 512>>>(nnz, nx, nz, npml, na, da, adcigs, 
                                                       s_P, s_Q, s_u0, s_w0, g_P, g_Q, g_u0, g_w0);
     }//it loop end

   }//is loop end

   migration_illum<<<(nx*nz+511)/512, 512>>>(nx, nz, npml, migration, illumination);
   adcigs_illum<<<(nx*nz*na+511)/512, 512>>>(nx, nz, na, da, adcigs, illumination);
   /* output multi-shot migration */
   cudaMemcpy(e, migration, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
   laplace_filter(1,nz,nx,e,d);
   fwrite(d,sizeof(float),nx*nz,fpmig);
   /* output multi-shot illumination */
   cudaMemcpy(e, illumination, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(e,sizeof(float),nx*nz,fpillum);
   /* output multi-shot adcigs */
   cudaMemcpy(Atemp, adcigs, nz*nx*na*sizeof(float), cudaMemcpyDeviceToHost);
   fwrite(Atemp,sizeof(float),nz*nx*na,fpadcigs);
   /* output adcigs stk migration */
   stk_adcigs(nx,nz,na,Atemp,d);
   fwrite(d,sizeof(float),nx*nz,fpstk);
   /* output smiled adcigs */
   adcigs_smiled(nx,nz,na,dcdp,Atemp);
   fwrite(Atemp,sizeof(float),nz*nx/dcdp*na,fpadcigs2);

   end = clock();
/*********IS Loop end*********/ 		     
   printf("\n---   Complete!!!!!!!!! \n");  
   printf("total %d shots: %f (min)\n", ns, ((float)(end-start))/60.0/CLOCKS_PER_SEC);



/***********close************/ 
     fclose(fpsnap);   fclose(fpobs);  fclose(fpmig);
     fclose(fpillum);  fclose(fpadcigs);fclose(fpstk);
     fclose(fpadcigs2);
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

       cudaFree(P_bndr);        cudaFree(Q_bndr);

       cudaFree(migration); 
       cudaFree(illumination);
       cudaFree(adcigs);
/***************host free*****************/
	free(v);	free(e);	free(d);
       free(shot_Hos);    free(Atemp);
}

