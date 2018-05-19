//a##################################################################
//a##
//a##                    CUDA-VTI-RTM-ADCIGs
//a##
//a##----------------------------------------------------------------
//a## Features:
//a##    Read initial models & shots derive Migrations and ADCIGs, 
//a## use poynting vector method to calculate angle of reflection.
//a## This is a CUDA code, initial code comes from "/Madagascar/user
//a## /pyang/Mgpu3dfd.cu". That code don't include boundary condition
//a## and it's a 3-D forward.
//a##----------------------------------------------------------------
//a##
//a##
//a##              | npml |mm|    nx    |mm| npml |
//a##           -- 0------------------------------> nx+2*mm+2*npml
//a##         npml |             npml             |
//a##           -- |   ------------------------   |
//a##           mm |   |          mm          |   |
//a##           -- |   |   ----------------   |   |
//a##              |   |   |              |   |   |
//a##              |   |   |              |   |   |
//a##           nz |   |   |     ved      |   |   |
//a##              |   |   |              |   |   |
//a##              |   |   |              |   |   |
//a##           -- |   |   ----------------   |   |
//a##           mm |   |                      |   |
//a##           -- |   ------------------------   |
//a##         npml |                              |
//a##           -- |-------------------------------
//a##       nz+2*mm+2*npml         
//a##
//a##---------------------------------------------------------------
//a## Ps: some of function you can search in Madagascar/user/pyang
//a##
//a##
//a##
//a##---------------------------------------------------------------
//a##                                            Rong Tao
//a##################################################################
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

#define PI 	3.141592653

#define BlockSize1 32// tile size in 1st-axis
#define BlockSize2 32// tile size in 2nd-axis

#define mm      4    // half of the order in space
#define npml     50   // absorbing boundry condition wield

//a################################################################################
__global__ void cuda_cal_c(float *c)
/*< get source >*/
{
  int id=threadIdx.x+blockIdx.x*blockDim.x;
    if(id<1){
       if(mm==2){
        c[0]=1.125;
        c[1]=-0.04166667;
	}else if(mm==3){
	 c[0]=1.1718750;
        c[1]=-0.065104167;
        c[2]=0.0046875;
	}else if(mm==4){
	 c[0]=1.196289;
        c[1]=-0.0797526;
        c[2]=0.009570313;
        c[3]=-0.0006975447;
	}else if(mm==5){
	 c[0]=1.211243;
        c[1]=-0.08972168;
        c[2]=0.01384277;
        c[3]=-0.00176566;
        c[4]=0.0001186795;
	}else if(mm==6){
        c[0]=1.2213364;
        c[1]=-0.096931458;
        c[2]=0.017447662;
        c[3]=-0.0029672895;
        c[4]=0.0003590054;
        c[5]=-0.000021847812;
   	}else if(mm==7){
        c[0]=1.2286062;
        c[1]=-0.10238385;
        c[2]=0.020476770;
        c[3]=-0.0041789327;
        c[4]=0.00068945355;
        c[5]=-0.000076922503;
        c[6]=0.0000042365148;
       }else if(mm==8){
        c[0]=1.2340911;
        c[1]=-0.10664985;
        c[2]=0.023036367;
        c[3]=-0.0053423856;
        c[4]=0.0010772712;
        c[5]=-0.00016641888;
        c[6]=0.000017021711;
        c[7]=-0.00000085234642;
   	}  
     }           
}
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
__constant__ float stencil[mm+1]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};
//a################################################################################
__global__ void cuda_ricker_wavelet(float *wlt, float fm, float dt, int nt)
/*< generate ricker wavelet with time deley >*/
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
    	if (it<nt){
	  float tmp = PI*fm*fabsf(it*dt-1.0/fm);//delay the wavelet to exhibit all waveform
	  tmp *=tmp;
	  wlt[it]= (1.0-2.0*tmp)*expf(-tmp);// ricker wavelet at time: t=nt*dt
	}
}
//a################################################################################
__global__ void cuda_set_s(int *szx, int fsz, int fsx, int dsz, int dsx, int ns, int nz, int nx)
/*< set the positions of sources  in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npml;
    	if (id<ns) szx[id]=(fsz+id*dsz+mm+npml)+nnz*(fsx+id*dsx+mm+npml);
}
//a################################################################################
__global__ void cuda_set_g(int *gzx, int ng, int nz, int nx)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npml;
       int ix=id%nx;
    	if (id<ng) gzx[id]=(mm+npml)+nnz*(ix*1+mm+npml);
}
//a################################################################################
__global__ void cuda_trans_x2tx(float *x, float *tx, int it, int nt, int ng, bool flag)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng){
          if(flag){
             tx[it+id*nt]=x[id];
           }else{
             x[id]=tx[it+id*nt];
            }
        }
}
//a################################################################################
__global__ void cuda_absorb_bndr(float *p1,float *p2,float *p3,float *p4,int nz,int nx,float qp)
/*< absorb boundry condition >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id;
	int nnz=nz+2*mm+2*npml;

          id=iz+ix*nnz;
            /*< left & right (0<x<nx) >*/
             if ( ix < npml ){
               p1[id]=( qp*pow((npml-ix)/(1.0*npml),2) + 1 )*p1[id];
               p2[id]=( qp*pow((npml-ix)/(1.0*npml),2) + 1 )*p2[id];
               p3[id]=( qp*pow((npml-ix)/(1.0*npml),2) + 1 )*p3[id];
               p4[id]=( qp*pow((npml-ix)/(1.0*npml),2) + 1 )*p4[id];
             }else if ( ix >= 2*mm + npml + nx ){
               p1[id]=( qp*pow((ix-2*mm-npml-nx)/(1.0*npml),2) + 1 )*p1[id];
               p2[id]=( qp*pow((ix-2*mm-npml-nx)/(1.0*npml),2) + 1 )*p2[id];
               p3[id]=( qp*pow((ix-2*mm-npml-nx)/(1.0*npml),2) + 1 )*p3[id];
               p4[id]=( qp*pow((ix-2*mm-npml-nx)/(1.0*npml),2) + 1 )*p4[id];
              }
            /*< up & down (0<z<nz) >*/
             if ( iz < npml ){
               p1[id]=( qp*pow((npml-iz)/(1.0*npml),2) + 1 )*p1[id];
               p2[id]=( qp*pow((npml-iz)/(1.0*npml),2) + 1 )*p2[id];
               p3[id]=( qp*pow((npml-iz)/(1.0*npml),2) + 1 )*p3[id];
               p4[id]=( qp*pow((npml-iz)/(1.0*npml),2) + 1 )*p4[id];
             }else if ( iz >= 2*mm + npml + nz ){
               p1[id]=( qp*pow((iz-2*mm-npml-nz)/(1.0*npml),2) + 1 )*p1[id]; 
               p2[id]=( qp*pow((iz-2*mm-npml-nz)/(1.0*npml),2) + 1 )*p2[id]; 
               p3[id]=( qp*pow((iz-2*mm-npml-nz)/(1.0*npml),2) + 1 )*p3[id]; 
               p4[id]=( qp*pow((iz-2*mm-npml-nz)/(1.0*npml),2) + 1 )*p4[id];   
               }
}
//a################################################################################
//__global__ void cuda_initial_PML(float *coff1, float *coff2, int nx, int nz, float dx, float dz, float dt, float vmax)
/*< PML boundry condition >*/
/*{
    int id=threadIdx.x+blockIdx.x*blockDim.x;

    float d0=3.0*vmax*log(100000.0)/(2.0*npml*(dx+dz)/2);
    int iz=id%(nz+2*mm+2*npml);
    int ix=id/(nz+2*mm+2*npml);

   if(id<(nx+2*npml*2*mm)*(nz+2*npml*2*mm))
   {
      if(iz<=npml){
            coff1[id]=1/(1+(dt*d0*pow((npml-0.5-iz)/npml,2))/2);
            coff2[id]=coff1[id]*(1-(dt*d0*pow((npml-0.5-iz)/npml,2))/2);
      }else if(iz>=nz+2*mm+npml){
            coff1[id]=1/(1+(dt*d0*pow((0.5+iz-nz-2*mm-npml)/npml,2))/2);
	     coff2[id]=coff1[id]*(1-(dt*d0*pow((0.5+iz-nz-2*mm-npml)/npml,2))/2);
      }if(ix<=npml&&(ix<=iz)&&ix<=(nz+2*mm+2*npml-iz)){
            coff1[id]=1/(1+(dt*d0*pow((npml-0.5-ix)/npml,2))/2);
            coff2[id]=coff1[id]*(1-(dt*d0*pow((npml-0.5-ix)/npml,2))/2);
      }else if((ix>=nx+2*mm+npml)&&(ix>=iz+nx-nz)&&(iz>=nx+2*mm+2*npml-ix)){
            coff1[id]=1/(1+(dt*d0*pow((0.5+ix-nx-2*mm-npml)/npml,2))/2);
	     coff2[id]=coff1[id]*(1-(dt*d0*pow((0.5+ix-nx-2*mm-npml)/npml,2))/2);
      }if(ix>=npml&&ix<=(npml+nx+2*mm)&&iz>=npml&&iz<=(npml+nz+2*mm)){
            coff1[id]=1.0;
	     coff2[id]=1.0;
      }
   }        
}*/
//a################################################################################
__global__ void cuda_record(float *p, float *seis, int *gx, int ng, bool flag)//++++++++++++
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) {
           if(flag){
                seis[id]=p[gx[id]];
            }else{
                p[gx[id]]=seis[id];
             }
        }
}
//a################################################################################
__global__ void cuda_add_source(bool add, float *p, float *source, int *szx, int ns)
/*< add/subtract sources: length of source[]=ns, index stored in szxy[] >*/
{
  int id=threadIdx.x+blockIdx.x*blockDim.x;

  if(id<ns){
    if(add){
      p[szx[id]]+=source[id];
    }else{
      p[szx[id]]-=source[id];
    }
  }
}
//a################################################################################
__global__ void cuda_step_fd2d(float *p0, float *p1, float *q0, float *q1, float *vv, float *vx, float *vn, 
                               float _dz2, float _dx2,int nz, int nx, 
                               bool forward, int *szx, int r, int R)
/*< step forward: 3-D FD, order=8 >*/
{
    bool validr = true;
    bool validw = true;
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix
    const int ltid1 = threadIdx.x;//ithreadz
    const int ltid2 = threadIdx.y;//ithreadx
    const int work1 = blockDim.x;//nblockz
    const int work2 = blockDim.y;//nblockx
    __shared__ float tile[BlockSize2 + 2 * mm][BlockSize1 + 2 * mm];
    __shared__ float tile2[BlockSize2 + 2 * mm][BlockSize1 + 2 * mm];

    float vvv, vvx, vvn;

    const int stride2 = nz + 2 * mm + 2 * npml;
    int inIndex = 0;
    int outIndex = 0;

    // Advance inputIndex to start of inner volume
    inIndex += (mm ) * stride2 + mm ;
    // Advance inputIndex to target element
    inIndex += ix * stride2 + iz;

    float current, current2;
    const int t1 = ltid1 + mm;
    const int t2 = ltid2 + mm;
    // Check in bounds
    if ((iz >= nz + mm + 2*npml) ||(ix >= nx + mm + 2*npml)) validr = false;
    if ((iz >= nz + 2*npml) ||(ix >= nx + 2*npml)) validw = false;

    if (validr) {current = p1[inIndex];current2 = q1[inIndex];}

    outIndex = inIndex;
    __syncthreads();

    if (ltid2 < mm){

       tile[ltid2][t1]                  = p1[outIndex - mm * stride2];
       tile[ltid2 + work2 + mm][t1]     = p1[outIndex + work2 * stride2];
       tile2[ltid2][t1]                  = q1[outIndex - mm * stride2];
       tile2[ltid2 + work2 + mm][t1]     = q1[outIndex + work2 * stride2];

    }if (ltid1 < mm){// Halo left & right

       tile[t2][ltid1]                  = p1[outIndex - mm];
       tile[t2][ltid1 + work1 + mm]     = p1[outIndex + work1];
       tile2[t2][ltid1]                  = q1[outIndex - mm];
       tile2[t2][ltid1 + work1 + mm]     = q1[outIndex + work1];
     }

    tile[t2][t1] = current;
    tile2[t2][t1] = current2;
    __syncthreads();

   // Compute the output value
    float c2, c3;
    c2=stencil[0]*current;
    c3=stencil[0]*current2;        

    for (int i=1; i <= mm ; i++){
	c2 +=stencil[i]*(tile[t2-i][t1]+ tile[t2+i][t1]);//x
       c3 +=stencil[i]*(tile2[t2][t1-i]+ tile2[t2][t1+i]);//z
     }
    c2*=_dx2;
    c3*=_dz2;	
    if (validw){
     if(!forward){
          vvv=vv[outIndex];
          vvx=vx[outIndex];
          vvn=vn[outIndex];
     }else{

       int iix=outIndex/stride2;
       int iiz=outIndex%stride2;

       int sx=*szx/stride2;
       int sz=*szx%stride2;

       int d=(int)sqrtf( ((iix-sx)*(iix-sx)) + ((iiz-sz)*(iiz-sz)) );

           if(d<=r){
               vvv=vv[outIndex];
               vvx=vv[outIndex];
               vvn=vv[outIndex];
           }else if(d>r&&d<=R){
               vvv=vv[outIndex];
               vvx=vv[outIndex]+(vx[outIndex]-vv[outIndex])/2.0*(1+cos( PI/(R-r)*(d-r) +PI ));
               vvn=vv[outIndex]+(vn[outIndex]-vv[outIndex])/2.0*(1+cos( PI/(R-r)*(d-r) +PI ));
           }else{
               vvv=vv[outIndex];
               vvx=vx[outIndex];
               vvn=vn[outIndex];
            }
      }

       p0[outIndex]=2.0*p1[outIndex]-p0[outIndex]
                   +vvv*(c3)+vvx*(c2);

       q0[outIndex]=2.0*q1[outIndex]-q0[outIndex]
                   +vvv*(c3)+vvn*(c2);

    }
}
//a################################################################################
//__global__ void cuda_iso_source_ring(float *vv, float *vx, float *vn, int nx, int nz, int *szx, int r, int R)
///*< Isotropic source ring >*/
//{
//	int id=threadIdx.x+blockDim.x*blockIdx.x;
//
//       int nnx=nx+2*mm+2*npml;
//       int nnz=nz+2*mm+2*npml;
//
//       int ix=id/nnz;
//       int iz=id%nnz;
//
//       int sx=*szx/nnz;
//       int sz=*szx%nnz;
//
//       int d=(int)sqrtf( ((ix-sx)*(ix-sx)) + ((iz-sz)*(iz-sz)) );
//
//    	if (id<nnx*nnz){
//           if(d<=r){
//               vx[id]=vv[id];
//               vn[id]=vv[id];
//            }else if(d>r&&d<=R){
//               vx[id]=vv[id]+(vx[id]-vv[id])/2.0*(1+cos( PI/(R-r)*(d-r) +PI ));
//               vn[id]=vv[id]+(vn[id]-vv[id])/2.0*(1+cos( PI/(R-r)*(d-r) +PI ));
//             }else{}
//        }
//
//}
//a################################################################################
void velocity_transform(float *v0, float*vv, float*vx, float*vn, float dt, int nz, int nx, float *vmax, float *vmin, float *vmute)
 /*< velocit2 transform: vv=v0*dt; vv<--vv^2 >*/
{
  int i1, i2, nnz, nnx;
  float tmp;

  nnz=nz+2*mm+2*npml;
  nnx=nx+2*mm+2*npml;
  *vmax=v0[0];
  *vmin=v0[0];
  *vmute=v0[0]*sqrtf(1+2*0.3);
  // inner zone
    for(i2=0; i2<nx; i2++){//x
      for(i1=0; i1<nz; i1++){//z
       if(*vmax<v0[i1+nz*i2])*vmax=v0[i1+nz*i2];
       if(*vmin>v0[i1+nz*i2])*vmin=v0[i1+nz*i2];
	tmp=v0[i1+nz*i2]*dt;
	vv[(i1+mm+npml)+nnz*(i2+mm+npml)]=tmp*tmp;
	vx[(i1+mm+npml)+nnz*(i2+mm+npml)]=tmp*tmp*(1+2*0.3);
	vn[(i1+mm+npml)+nnz*(i2+mm+npml)]=tmp*tmp*(1+2*0.2);
      }
    }
    //top & down 
	for(i2=0; i2<nnx; i2++){//x
	    for (i1=0; i1<mm+npml; i1++){//z
		vv[i1+nnz*i2]=vv[mm+npml+nnz*i2];
		vv[(nnz-i1-1)+nnz*i2]=vv[(nnz-mm-npml-1)+nnz*i2];
		vx[i1+nnz*i2]=vx[mm+npml+nnz*i2];
		vx[(nnz-i1-1)+nnz*i2]=vx[(nnz-mm-npml-1)+nnz*i2];
		vn[i1+nnz*i2]=vn[mm+npml+nnz*i2];
		vn[(nnz-i1-1)+nnz*i2]=vn[(nnz-mm-npml-1)+nnz*i2];
	    }
	}
    //left & right
	for(i2=0; i2<mm+npml; i2++){//x
	    for (i1=0; i1<nnz; i1++){//z
		vv[i1+nnz*i2]=vv[i1+nnz*(mm+npml)];
		vv[i1+nnz*(nnx-i2-1)]=vv[i1+nnz*(nnx-mm-npml-1)];
		vx[i1+nnz*i2]=vx[i1+nnz*(mm+npml)];
		vx[i1+nnz*(nnx-i2-1)]=vx[i1+nnz*(nnx-mm-npml-1)];
		vn[i1+nnz*i2]=vn[i1+nnz*(mm+npml)];
		vn[i1+nnz*(nnx-i2-1)]=vn[i1+nnz*(nnx-mm-npml-1)];
	    }
	}
}
//a################################################################################
void window3d(float *a, float *b, int nz, int nx)
/*< window a 3d subvolume >*/
{
	int i1, i2, nnz;
	nnz=nz+2*mm+ 2*npml;//z
	
	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	{
          a[i1+nz*i2]=b[(i1+mm+npml)+nnz*(i2+mm+npml)];
	}
}
//a################################################################################
__global__ void cuda_set_cooLR(int *left, int *right, int nx, int nz)
/*< set the positions of  left & right in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npml;
    	if (id<nz){
            left[id]=(nnz+1)*(mm+npml)+id;
           right[id]=nnz*(mm+npml+nx)+mm+npml+id;
        }
}
//a################################################################################
__global__ void cuda_set_cooUD(int *up, int *down, int nx, int nz)
/*< set the positions of  up & down in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npml;
    	if (id<nx){
             up[id]=(nnz+1)*(mm+npml)+id*nnz;
           down[id]=(nnz+1)*(mm+npml)+id*nnz+nz;
        }
}
//a#############################################################################################
__global__ void cuda_s_bndr(float *s_p_bndr, float *s_q_bndr, float *p, float *q, 
                            int *left, int *right, int *up, int *down, int nz, int nx, bool write)
/*< write boundaries out or read them into wavefield variables p>*/
{
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	if(write){
		if(id<nz){ /* left  boundary */
                    s_p_bndr[id]=p[left[id]]; 
                    s_q_bndr[id]=q[left[id]];         
		}else if((id>=nz)&&(id<(2*nz))){ /* right boundary */
                    s_p_bndr[id]=p[right[id-nz]]; 
                    s_q_bndr[id]=q[right[id-nz]];  
		}else if(id>=(2*nz)&&(id<(2*nz+nx))){  /* up    boundary */
                    s_p_bndr[id]=p[up[id-2*nz]];  
                    s_q_bndr[id]=q[up[id-2*nz]];  
		}else if(id>=(2*nz+nx)&&id<(2*nz+2*nx)){ /* down boundary */
                    s_p_bndr[id]=p[down[id-2*nz-nx]];
                    s_q_bndr[id]=q[down[id-2*nz-nx]];
                }
	}else{
		if(id<nz){ /* left  boundary */
                    p[left[id]]=s_p_bndr[id]; 
                    q[left[id]]=s_q_bndr[id];         
		}else if((id>=nz)&&(id<(2*nz))){ /* right boundary */
                    p[right[id-nz]]=s_p_bndr[id]; 
                    q[right[id-nz]]=s_q_bndr[id];  
		}else if(id>=(2*nz)&&(id<(2*nz+nx))){  /* up    boundary */
                    p[up[id-2*nz]]=s_p_bndr[id];  
                    q[up[id-2*nz]]=s_q_bndr[id];  
		}else if(id>=(2*nz+nx)&&id<(2*nz+2*nx)){ /* down boundary */
                    p[down[id-2*nz-nx]]=s_p_bndr[id];
                    q[down[id-2*nz-nx]]=s_q_bndr[id];
                }
        }
}
//a################################################################################
__global__ void cuda_cal_corr(float *mig, float *s, float *g, int nx, int nz)
/*< correlation imaging condition >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<(nx+2*mm+2*npml)*(nz+2*mm+2*npml)){
              mig[id]+=s[id]*g[id];
        }
}
//a################################################################################
__global__ void cuda_cal_illum_matrix(float *illum, float *wave, int nx, int nz)
/*< illumination matrix >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<(nx+2*mm+2*npml)*(nz+2*mm+2*npml)){
              illum[id]+=wave[id]*wave[id];
        }
}
//a################################################################################
__global__ void cuda_illumination(float *mig_ns, float *mig_is, float *illum_ns, float *illum_is, int nx, int nz)
/*< illumination matrix >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<(nx+2*mm+2*npml)*(nz+2*mm+2*npml)){
              mig_is[id]/=illum_is[id];
              mig_ns[id]+=mig_is[id];
              illum_ns[id]+=illum_is[id];
        }
}
//a################################################################################
__global__ void cuda_mute_direct(float *p, int nx, int nz, float dx, float dz, int nt, float dt, float fm, 
                                 float vmute, int *coo_zx, int tt)
/*< illumination matrix >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

       int ix=id/nt;
       int it=id%nt;

       int sx=*coo_zx/(nz+2*mm+2*npml)-mm-npml;
       int sz=*coo_zx%(nz+2*mm+2*npml)-mm-npml;

       int t0=(int)sqrtf((dx*(ix-sx)*dx*(ix-sx))+(dz*sz*dz*sz));
       int t1=(int)(t0/vmute/dt);
       int t2=(int)(2.0/(dt*fm));

    	if (id<nx*nt){
          if(it <= t1+t2+tt){
              p[id]=0.0;     
           }else{}
        }

}
//a################################################################################
__global__ void cuda_mut_v(float *s_vv, float *s_vx, float *s_vn, float *mut_vv, float *mut_vx, float *mut_vn, 
                           int nx, int nz)
/*< copy velocity to mute velocity >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

       int nnx=nx+2*mm+2*npml;
       int nnz=nz+2*mm+2*npml;

    	if (id<nnx*nnz){
           mut_vv[id]=s_vv[0];
           mut_vx[id]=s_vx[0];
           mut_vn[id]=s_vn[0];
        }
}
//a################################################################################
__global__ void cuda_difference(float *p1, float *p2, int nx, int nt)
/*< get wavefiled(x-t) difference >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

    	if (id<nx*nt) p1[id]-=p2[id];
}
//a################################################################################
__global__ void cuda_poynting_adcigs(float *adcigs_Dev, float *s_p0, float *s_q0, float *g_p0, float *g_q0, 
                                     int nx, int nz, int na, float _dx2, float _dz2, float *illum_is, int *angle_count_Dev, 
                                     float *c)
/*< poynting vector get ADCIGs >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

       int nnz=nz+2*mm+2*npml;
       int nnx=nx+2*mm+2*npml;

       int iz=id%nnz-mm-npml;
       int ix=id/nnz-mm-npml;

       float Ssx, Ssz, Sgx, Sgz, tmp, s_u=0.0, s_w=0.0, g_u=0.0, g_w=0.0;

       s_u=stencil[0]*s_p0[id]*sqrtf(_dx2);
       s_w=stencil[0]*s_q0[id]*sqrtf(_dz2);
       g_u=stencil[0]*g_p0[id]*sqrtf(_dx2);
       g_w=stencil[0]*g_q0[id]*sqrtf(_dz2);

       for(int im=1; im<=mm; im++){
            s_u+=stencil[im]*(s_p0[id+im*nnz]+s_p0[id-im*nnz])*sqrtf(_dx2);
            s_w+=stencil[im]*(s_p0[id+im]    +s_p0[id-im])      *sqrtf(_dz2);
            g_u+=stencil[im]*(g_p0[id+im*nnz]+g_p0[id-im*nnz])*sqrtf(_dx2);
            g_w+=stencil[im]*(g_p0[id+im]    +g_p0[id-im])      *sqrtf(_dz2);
        }

       Ssx=-s_p0[id]*s_u;
       Ssz=-s_q0[id]*s_w;
       Sgx= g_p0[id]*g_u;
       Sgz= g_q0[id]*g_w;

       float b1=Ssz*Ssz+Ssx*Ssx;
       float b2=Sgz*Sgz+Sgx*Sgx;

       float a= 0.5*acosf( (Ssx*Sgx+Ssz*Sgz) / (sqrtf(b1*b2)*(1+0.0)) );

       int ia=(int)(  a*180/PI  );

    	if (((ia>=0)&&(ia<na))&&iz>-1&&ix>-1&&id<nnz*nnx&&(iz+nz*ia+nz*na*ix)<nz*na*nx){
              angle_count_Dev[ia]++;
              tmp=s_p0[id]*g_p0[id];//+s_q0[id]*g_q0[id];
              //tmp=tmp*expf(-(a-ia)*(a-ia)/ ( 2.0/9.0 ) );
              adcigs_Dev[iz+nz*ia+nz*na*ix]+=tmp/illum_is[id];
        }
}
//a################################################################################
__global__ void cuda_smooth_adcigs(float *adcigs, int nx, int nz, int na, int nsmooth)
/*< poynting vector get ADCIGs >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

       int iz=id%nz;
       int ix=id/nz;

    if (id<nz*nx)
      for(int in=0; in<nsmooth; in++){
         for(int ia=0; ia<na; ia++){
            if((iz+nz*ia+nz*na*ix)<nz*na*nx)
              if(ia==0){
                adcigs[iz+nz*ia+nz*na*ix]=(adcigs[iz+nz*na*ix]+adcigs[iz+nz+nz*na*ix])/2.0;
              }else if(ia==na-1){
                adcigs[iz+nz*ia+nz*na*ix]=(adcigs[iz+nz*(na-1)+nz*na*ix]+adcigs[iz+nz*(na-2)+nz*na*ix])/2.0;
              }else{
                adcigs[iz+nz*ia+nz*na*ix]=(adcigs[iz+nz*(ia-1)+nz*na*ix]+adcigs[iz+nz*(ia+1)+nz*na*ix]+adcigs[iz+nz*ia+nz*na*ix])/3.0;
               }
          }
       }  
}
//a################################################################################
//__global__ void cuda_uw(float *p, float *q, float _dz2, float _dx2,int nz, int nx, float dt, float *u, float *w)
///*< get zhidian velocity u,w >*/
//{
//	int id=threadIdx.x+blockDim.x*blockIdx.x;
//
//       int nnz=nz+2*mm+2*npml;
//       int nnx=nx+2*mm+2*npml;
//       float diffx, diffz;
//
//       if(id<nnx*nnz){
//           diffx=1.125*(p[id+nnz]-p[id])-0.0416666666667*(p[id+2*nnz]-p[id+nnz]);
//           diffz=1.125*(q[id+1]-q[id])-0.0416666666667*(q[id+2]-q[id+1]);
//           u[id]+=1000*dt*sqrtf(_dx2)*diffx;
//           w[id]+=1000*dt*sqrtf(_dz2)*diffz;
//        }
//}
//a################################################################################
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
//a###                                                                         ####
//a###                             Main Function                               ####
//a###                                                                         ####
//a################################################################################
int main(int argc, char* argv[])
{
	int nz, nx, nnz, nnx, ns, nt, kt, it, is, fsz, fsx,  dsz, dsx, ng, i, na, r, R;
	int *szx, *gzx, *coo_left, *coo_right, *coo_up, *coo_down, *angle_count_Dev, *angle_count_Hos;
	float dz, dx,  fm, dt, _dz2, _dx2, vmax=0.0, vmin=0.0, vmute=0.0;
	float *temp, *v0, *vv, *vx, *vn, *s_wlt, *s_vv, *s_vx, *s_vn, *s_p0, *s_p1, *s_q0, *s_q1, *g_p0, *g_p1, *g_q0, *g_q1, *ptr;
       float *cal_it_Dev, *cal_nt_Dev, *cal_nt_Hos, *obs_it_Dev, *obs_nt_Dev, *obs_nt_Hos, *s_p_bndr, *s_q_bndr;
       float *c;
       float *mig_is, *mig_ns, *illum_is, *illum_ns;
       float *mut_vv, *mut_vx, *mut_vn, *mut_p0, *mut_p1, *mut_q0, *mut_q1, *mut_nt_Dev;
       float *adcigs_Dev, *adcigs_Hos;
       bool flag_snap, flag_laplace, flag_adcigs_smooth;
//a##################################################
	char     FNvel[250]={"vel_600_300.dat"};
       char    FNsnap[250]={"snap.dat"};
       char FNshotcal[250]={"shot_cal.dat"};
       char FNshotobs[250]={"shot_cal.dat"};
       char   FNmigis[250]={"mig_is.dat"};
       char   FNmigns[250]={"mig_ns.dat"};
       char FNillumis[250]={"illum_is.dat"};
       char FNillumns[250]={"illum_ns.dat"};
       char  FNadcigs[250]={"adcigs.dat"};
//a##################################################
             flag_snap=true;
          flag_laplace=true;
    flag_adcigs_smooth=true;
//a##################################################
       r=10; R=20;
//a##################################################
       fm=20;     

    	nx=300;   dx=10;
    	nz=300;   dz=10;
    	
   	nt=2501;   kt=50;    dt=0.001;

   	ns=10;
       fsx=nx/ns/2;    dsx=nx/ns;
       fsz=1;    dsz=0;

       na=70;

//a##################################################
       FILE *fpvel, *fpsnap, *fpshotcal, *fpshotobs, *fpmigis, *fpmigns, *fpillumis, *fpillumns, *fpadcigs;
       if((fpvel=fopen(FNvel,"rb"))==NULL) printf("ERROR:open %s error!\n",FNvel);
       if(flag_snap) fpsnap=fopen(FNsnap,"wb");
       fpshotcal=fopen(FNshotcal,"wb");
       fpshotobs=fopen(FNshotobs,"rb");
         fpmigis=fopen(FNmigis,"wb");
         fpmigns=fopen(FNmigns,"wb");
       fpillumis=fopen(FNillumis,"wb");
       fpillumns=fopen(FNillumns,"wb");
        fpadcigs=fopen(FNadcigs,"wb");

//a##################################################
	_dz2=1.0/(dz*dz);
	_dx2=1.0/(dx*dx);
	nnz=nz+2*mm+2*npml;
	nnx=nx+2*mm+2*npml;
       ng=nx;
//a##################################################
    	v0=(float*)malloc(nz*nx*sizeof(float));
    	temp=(float*)malloc(nz*nx*sizeof(float));
    	vv=(float*)malloc(nnz*nnx*sizeof(float));
    	vx=(float*)malloc(nnz*nnx*sizeof(float));
    	vn=(float*)malloc(nnz*nnx*sizeof(float));
    	cal_nt_Hos=(float*)malloc(ng*nt*sizeof(float));
    	obs_nt_Hos=(float*)malloc(ng*nt*sizeof(float));
    	adcigs_Hos=(float*)malloc(nz*na*nx*sizeof(float));
    	angle_count_Hos=(int*)malloc(na*sizeof(int));
	fread(v0, sizeof(float), nz*nx, fpvel);
	velocity_transform(v0, vv, vx, vn, dt, nz, nx, &vmax, &vmin, &vmute);
//a##################################################
       printf("###############################\n");
       printf("##  vmin=%.2f, vmax=%.2f\n",vmin,vmax);
       printf("###############################\n");
//a##################################################
    	cudaSetDevice(0);// initialize device, default device=0;
	check_gpu_error("Failed to initialize device!");
//a##################################################
	dim3 dimg, dimb;
	dimg.x=(nz+2*npml+2*mm+BlockSize1-1)/BlockSize1;
	dimg.y=(nx+2*npml+2*mm+BlockSize2-1)/BlockSize2;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;
//a##################################################
	cudaMalloc(&s_wlt, nt*sizeof(float));
	cudaMalloc(&c, mm*sizeof(float));
	cudaMalloc(&s_vv, nnz*nnx*sizeof(float));
	cudaMalloc(&s_vx, nnz*nnx*sizeof(float));
	cudaMalloc(&s_vn, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_vv, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_vx, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_vn, nnz*nnx*sizeof(float));
	cudaMalloc(&s_p0, nnz*nnx*sizeof(float));
	cudaMalloc(&s_p1, nnz*nnx*sizeof(float));
	cudaMalloc(&s_q0, nnz*nnx*sizeof(float));
	cudaMalloc(&s_q1, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_p0, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_p1, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_q0, nnz*nnx*sizeof(float));
	cudaMalloc(&mut_q1, nnz*nnx*sizeof(float));
	cudaMalloc(&g_p0, nnz*nnx*sizeof(float));
	cudaMalloc(&g_p1, nnz*nnx*sizeof(float));
	cudaMalloc(&g_q0, nnz*nnx*sizeof(float));
	cudaMalloc(&g_q1, nnz*nnx*sizeof(float));
//a##################################################
	cudaMalloc(&mig_is, nnz*nnx*sizeof(float));
	cudaMalloc(&mig_ns, nnz*nnx*sizeof(float));
	cudaMalloc(&illum_is, nnz*nnx*sizeof(float));
	cudaMalloc(&illum_ns, nnz*nnx*sizeof(float));
//a##################################################
	cudaMalloc(&adcigs_Dev, nz*na*nx*sizeof(float));
	cudaMalloc(&angle_count_Dev, na*sizeof(int));
//a##################################################
	cudaMalloc(&szx, ns*sizeof(int));
	cudaMalloc(&gzx, ng*sizeof(int));
	cudaMalloc(&coo_left , nz*sizeof(int));
	cudaMalloc(&coo_right, nz*sizeof(int));
	cudaMalloc(&coo_up   , nx*sizeof(int));
	cudaMalloc(&coo_down , nx*sizeof(int));
//a##################################################
	cudaMalloc(&cal_it_Dev, ng*sizeof(float));	
	cudaMalloc(&cal_nt_Dev, ng*nt*sizeof(float));
	cudaMalloc(&obs_it_Dev, ng*sizeof(float));	
	cudaMalloc(&obs_nt_Dev, ng*nt*sizeof(float));
	cudaMalloc(&mut_nt_Dev, ng*nt*sizeof(float));
	cudaMalloc(&s_p_bndr, (2*nx+2*nz)*nt*sizeof(float)); 
	cudaMalloc(&s_q_bndr, (2*nx+2*nz)*nt*sizeof(float)); 
//a##################################################
	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(s_wlt, fm, dt, nt);
	check_gpu_error("Failed to allocate memory for variables!");
//a##################################################
	cudaMemcpy(s_vv, vv, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(s_vx, vx, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(s_vn, vn, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
//a##################################################
       cuda_mut_v<<<(nnx*nnz+511)/512, 512>>>(s_vv, s_vx, s_vn, mut_vv, mut_vx, mut_vn, nx, nz);
//a##################################################
       cuda_cal_c<<<1, 1>>>(c);
	cuda_set_s<<<1, ns>>>(szx, fsz, fsx, dsz, dsx, ns, nz, nx);
	cuda_set_g<<<(ng+511)/512,512>>>(gzx, ng, nz, nx);
       cuda_set_cooLR<<<(nz+511)/512,512>>>(coo_left, coo_right, nx, nz);
       cuda_set_cooUD<<<(nx+511)/512,512>>>(coo_up,   coo_down,  nx, nz);
//a##################################################
	clock_t t0, t1;
	t0 = clock();
//a##################################################
	cudaMemset(adcigs_Dev, 0, nz*na*nx*sizeof(float));
	cudaMemset(mig_ns, 0, nnz*nnx*sizeof(float));
	cudaMemset(illum_ns, 0, nnz*nnx*sizeof(float));
	cudaMemset(angle_count_Dev, 0, na*sizeof(int));
     for(is=0; is<ns; is++){
	  cudaMemset(s_p0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_p1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_q0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_q1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(mut_p0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(mut_p1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(mut_q0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(mut_q1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(g_p0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(g_p1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(g_q0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(g_q1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(mig_is, 0, nnz*nnx*sizeof(float));
	  cudaMemset(illum_is, 0, nnz*nnx*sizeof(float));
	  cudaMemset(cal_it_Dev, 0, ng*sizeof(float));
	  cudaMemset(cal_nt_Dev, 0, ng*nt*sizeof(float));
	  cudaMemset(obs_it_Dev, 0, ng*sizeof(float));
	  cudaMemset(obs_nt_Dev, 0, ng*nt*sizeof(float));
	  cudaMemset(obs_it_Dev, 0, ng*sizeof(float));
	  cudaMemset(mut_nt_Dev, 0, ng*nt*sizeof(float));
//a##################################################
         //cuda_iso_source_ring<<<(nnz*nnx+511)/512,512>>>(s_vv, s_vx, s_vn, nx, nz, &szx[is], 5, 15);
//a####################################################################################################
         printf("##  >>  is =%3d\n",is);
	  for(it=0; it<nt; it++){
	    cuda_add_source<<<1,1>>>(true, s_p1, &s_wlt[it], &szx[is], 1);
	    cuda_add_source<<<1,1>>>(true, s_q1, &s_wlt[it], &szx[is], 1);
	    cuda_step_fd2d<<<dimg,dimb>>>(s_p0, s_p1, s_q0, s_q1, s_vv, s_vx, s_vn, _dz2, _dx2, nz, nx,
                                        true, &szx[is], r, R);
           cuda_absorb_bndr<<<dimg,dimb>>>(s_p0, s_p1, s_q0, s_q1, nz, nx, -0.25);
           cuda_s_bndr<<<((2*nx+2*nz)+511)/512,512>>>(&s_p_bndr[it*(2*nx+2*nz)], &s_q_bndr[it*(2*nx+2*nz)], 
                                                      s_p0, s_q0, coo_left, coo_right, coo_up, coo_down, nz, nx, true);
	    ptr=s_p0; s_p0=s_p1; s_p1=ptr;
	    ptr=s_q0; s_q0=s_q1; s_q1=ptr;
//a##################################################
           cuda_cal_illum_matrix<<<(nnx*nnz+511)/512, 512>>>(illum_is, s_p0, nx, nz);
//a##################################################
	    cuda_record<<<(ng+511)/512, 512>>>(s_p0, cal_it_Dev, gzx, ng, true);
           cuda_trans_x2tx<<<(ng+511)/512, 512>>>(cal_it_Dev, cal_nt_Dev, it, nt, ng, true);
//a####################################################################################################
	    cuda_add_source<<<1,1>>>(true, mut_p1, &s_wlt[it], &szx[is], 1);
	    cuda_add_source<<<1,1>>>(true, mut_q1, &s_wlt[it], &szx[is], 1);
	    cuda_step_fd2d<<<dimg,dimb>>>(mut_p0, mut_p1, mut_q0, mut_q1, mut_vv, mut_vx, mut_vn, _dz2, _dx2, nz, nx,
                                        true, &szx[is], r, R);
           cuda_absorb_bndr<<<dimg,dimb>>>(mut_p0, mut_p1, mut_q0, mut_q1, nz, nx, -0.25);

	    ptr=mut_p0; mut_p0=mut_p1; mut_p1=ptr;
	    ptr=mut_q0; mut_q0=mut_q1; mut_q1=ptr;
//a##################################################
	    cuda_record<<<(ng+511)/512, 512>>>(mut_p0, cal_it_Dev, gzx, ng, true);
           cuda_trans_x2tx<<<(ng+511)/512, 512>>>(cal_it_Dev, mut_nt_Dev, it, nt, ng, true);
//a##################################################
	    if(is==0&&it!=0&&it%kt==0&&flag_snap){
	      cudaMemcpy(vv, s_p0, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
	      window3d(temp, vv, nz, nx);
	      fwrite(temp, sizeof(float),nz*nx, fpsnap);	  
             }
            }
//a##################################################
         //cuda_mute_direct<<<(ng*nt+511)/512,512>>>(cal_nt_Dev, nx, nz, dx, dz, nt, dt, fm, vmute, &szx[is], 30);
//a##################################################
           cuda_difference<<<(ng*nt+511)/512, 512>>>(cal_nt_Dev, mut_nt_Dev, ng, nt);
//a##################################################
	  cudaMemcpy(cal_nt_Hos, cal_nt_Dev, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
         fseek(fpshotcal,is*ng*nt*sizeof(float),0);
	  fwrite(cal_nt_Hos, sizeof(float), ng*nt, fpshotcal);
//a##################################################
         for(i=0;i<ng*nt;i++)
               obs_nt_Hos[i]=cal_nt_Hos[i];
	  cudaMemcpy(obs_nt_Dev, obs_nt_Hos, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
//a####################################################################################################
	  cudaMemset(s_p0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_p1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_q0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_q1, 0, nnz*nnx*sizeof(float));
//a##################################################
	  for(it=nt-1; it>=0; it--){
//a##################################################
	    ptr=s_p0; s_p0=s_p1; s_p1=ptr;
	    ptr=s_q0; s_q0=s_q1; s_q1=ptr;
           cuda_s_bndr<<<((2*nx+2*nz)+511)/512,512>>>(&s_p_bndr[it*(2*nx+2*nz)], &s_q_bndr[it*(2*nx+2*nz)], 
                                                      s_p1, s_q1, coo_left, coo_right, coo_up, coo_down, nz, nx, false);
           cuda_step_fd2d<<<dimg,dimb>>>(s_p0, s_p1, s_q0, s_q1, s_vv, s_vx, s_vn, _dz2, _dx2, nz, nx,
                                        false, NULL, NULL, NULL);
           cuda_absorb_bndr<<<dimg,dimb>>>(s_p0, s_p1, s_q0, s_q1, nz, nx, -0.25);
//a##################################################
	  //  if(is==0&&it!=0&&it%kt==0&&flag_snap){
	  //    cudaMemcpy(vv, s_p0, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
	  //    window3d(temp, vv, nz, nx);
	  //    fwrite(temp, sizeof(float),nz*nx, fpsnap);	  
          //   }
//a##################################################
           cuda_trans_x2tx<<<(ng+511)/512, 512>>>(obs_it_Dev, obs_nt_Dev, it, nt, ng, false);
	    cuda_record<<<(ng+511)/512, 512>>>(g_p1, obs_it_Dev, gzx, ng, false);
	    cuda_record<<<(ng+511)/512, 512>>>(g_q1, obs_it_Dev, gzx, ng, false);
           cuda_step_fd2d<<<dimg,dimb>>>(g_p0, g_p1, g_q0, g_q1, s_vv, s_vx, s_vn, _dz2, _dx2, nz, nx,
                                        false, NULL, NULL, NULL);
           cuda_absorb_bndr<<<dimg,dimb>>>(g_p0, g_p1, g_q0, g_q1, nz, nx, -0.25);
	    ptr=g_p0; g_p0=g_p1; g_p1=ptr;
	    ptr=g_q0; g_q0=g_q1; g_q1=ptr;
           cuda_cal_illum_matrix<<<(nnx*nnz+511)/512, 512>>>(illum_is, g_p0, nx, nz);
//a##################################################
	  //  if(is==0&&it!=0&&it%kt==0&&flag_snap){
	  //    cudaMemcpy(vv, g_p0, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
	  //    window3d(temp, vv, nz, nx);
	  //    fwrite(temp, sizeof(float),nz*nx, fpsnap);	  
          //   }
//a##################################################
           cuda_cal_corr<<<(nnx*nnz+511)/512, 512>>>(mig_is, s_p1, g_p1, nx, nz);
//a##################################################
           cuda_poynting_adcigs<<<(nnx*nnz+511)/512, 512>>>(adcigs_Dev, s_p0, s_q0, g_p0, g_q0, nx, nz, na, _dx2, _dz2, 
                                                            illum_is, angle_count_Dev, c);
          }
//a##################################################
          cudaMemcpy(vv, illum_is, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
	   window3d(temp, vv, nz, nx);
          fseek(fpillumis,is*nx*nz*sizeof(float),0);
	   fwrite(temp, sizeof(float),nz*nx, fpillumis);
//a##################################################
          cuda_illumination<<<(nnx*nnz+511)/512, 512>>>(mig_ns, mig_is, illum_ns, illum_is, nx, nz);
//a##################################################
          cudaMemcpy(vv, mig_is, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
	   window3d(temp, vv, nz, nx);
          fseek(fpmigis,is*nx*nz*sizeof(float),0);
	   fwrite(temp, sizeof(float),nz*nx, fpmigis);
      }//end of IS loop
//a##################################################
     cudaMemcpy(vv, mig_ns, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
     window3d(temp, vv, nz, nx);
     if(flag_laplace){
          laplace_filter(1, nz, nx, temp, v0);
          fwrite(v0, sizeof(float),nz*nx, fpmigns);
     }else{ fwrite(temp, sizeof(float),nz*nx, fpmigns); }
//a##################################################
     cudaMemcpy(vv, illum_ns, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
     window3d(temp, vv, nz, nx);
     fwrite(temp, sizeof(float),nz*nx, fpillumns);
//a##################################################
     if(flag_adcigs_smooth)cuda_smooth_adcigs<<<(nx*nz+511)/512, 512>>>(adcigs_Dev, nx, nz, na, 10);
     cudaMemcpy(adcigs_Hos, adcigs_Dev, nz*na*nx*sizeof(float), cudaMemcpyDeviceToHost);
     fwrite(adcigs_Hos, sizeof(float),nz*na*nx, fpadcigs);
//a##################################################
     cudaMemcpy(angle_count_Hos, angle_count_Dev, na*sizeof(int), cudaMemcpyDeviceToHost);
     for(i=0;i<na;i++) printf("The number of %2d degree is %d.\n",i+1,angle_count_Hos[i]);
//a##################################################


     t1 = clock();
     printf("total %d shots: %f (s)\n", ns, ((float)(t1-t0))/CLOCKS_PER_SEC);

	/* free memory on device */
	cudaFree(c);
	cudaFree(s_wlt);
	cudaFree(s_vv);
	cudaFree(s_vx);
	cudaFree(s_vn);
	cudaFree(s_p0);
	cudaFree(s_p1);
	cudaFree(s_q0);
	cudaFree(s_q1);
	cudaFree(mut_vv);
	cudaFree(mut_vx);
	cudaFree(mut_vn);
	cudaFree(mut_p0);
	cudaFree(mut_p1);
	cudaFree(mut_q0);
	cudaFree(mut_q1);
	cudaFree(g_p0);
	cudaFree(g_p1);
	cudaFree(g_q0);
	cudaFree(g_q1);
	cudaFree(mig_is);
	cudaFree(mig_ns);
	cudaFree(illum_is);
	cudaFree(illum_ns);
	cudaFree(szx);
	cudaFree(gzx);
	cudaFree(cal_it_Dev);
	cudaFree(cal_nt_Dev);
	cudaFree(obs_it_Dev);
	cudaFree(obs_nt_Dev);
	cudaFree(mut_nt_Dev);
	cudaFree(adcigs_Dev);
	cudaFree(angle_count_Dev);

	free(v0);
	free(temp);
	free(vv);
	free(vx);
	free(vn);
	free(cal_nt_Hos);
	free(obs_nt_Hos);
	free(adcigs_Hos);
	free(angle_count_Hos);

       fclose(fpvel);
       if(flag_snap) fclose(fpsnap);
       fclose(fpshotcal);
       fclose(fpshotobs);
       fclose(fpmigis);
       fclose(fpmigns);
       fclose(fpillumis);
       fclose(fpillumns);
       fclose(fpadcigs);

    	exit (0);
}

