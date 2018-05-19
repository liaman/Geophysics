//a##################################################################
//a##
//a##                    CUDA-VTI-FD
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

#define BlockSize1 16// tile size in 1st-axis
#define BlockSize2 16// tile size in 2nd-axis

#define mm      4    // half of the order in space
#define npd     50   // absorbing boundry condition wield
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
	int nnz=nz+2*mm+2*npd;
    	if (id<ns) szx[id]=(fsz+id*dsz+mm+npd)+nnz*(fsx+id*dsx+mm+npd);
}
//a################################################################################
__global__ void cuda_set_g(int *gzx, int ng, int nz, int nx)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
       int ix=id%nx;
    	if (id<ng) gzx[id]=(mm+npd)+nnz*(ix*1+mm+npd);
}
//a################################################################################
__global__ void cuda_trans_x2tx(float *x, float *tx, int it, int nt, int ng)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) tx[it+id*nt]+=x[id];
       
}
//a################################################################################
__global__ void cuda_absorb_bndr(float *p,int nz,int nx,float qp)
/*< absorb boundry condition >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id;
	int nnz=nz+2*mm+2*npd;

          id=iz+ix*nnz;
            /*< left & right (0<x<nx) >*/
             if ( ix < npd )
               p[id]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*p[id];
             else if ( ix >= 2*mm + npd + nx )
               p[id]=( qp*pow((ix-2*mm-npd-nx)/(1.0*npd),2) + 1 )*p[id];
            /*< up & down (0<z<nz) >*/
             if ( iz < npd )
               p[id]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*p[id];
             else if ( iz >= 2*mm + npd + nz )
               p[id]=( qp*pow((iz-2*mm-npd-nz)/(1.0*npd),2) + 1 )*p[id];   
}
//a################################################################################
__global__ void cuda_initial_PML(float *coff1, float *coff2, int nx, int nz, float dx, float dz, float dt, float vmax)
/*< PML boundry condition >*/
{
    int id=threadIdx.x+blockIdx.x*blockDim.x;

    float d0=6.0*vmax*log(100000.0)/(2.0*npd*(dx+dz)/2);
    int iz=id%(nz+2*mm+2*npd);
    int ix=id/(nz+2*mm+2*npd);

   if(id<(nx+2*npd*2*mm)*(nz+2*npd*2*mm))
   {
      if(iz<npd){
            coff1[id]=1/(1+(dt*d0*pow((npd-0.5-iz)/npd,2))/2);
            coff2[id]=coff1[id]*(1-(dt*d0*pow((npd-0.5-iz)/npd,2))/2);
      }else if(iz>=nz+2*mm+npd){
            coff1[id]=1/(1+(dt*d0*pow((0.5+iz-nz-2*mm-npd)/npd,2))/2);
	     coff2[id]=coff1[id]*(1-(dt*d0*pow((0.5+iz-nz-2*mm-npd)/npd,2))/2);
      }if(ix<npd){
            coff1[id]=1/(1+(dt*d0*pow((npd-0.5-ix)/npd,2))/2);
            coff2[id]=coff1[id]*(1-(dt*d0*pow((npd-0.5-ix)/npd,2))/2);
      }else if(ix>=nx+2*mm+npd){
            coff1[id]=1/(1+(dt*d0*pow((0.5+ix-nx-2*mm-npd)/npd,2))/2);
	     coff2[id]=coff1[id]*(1-(dt*d0*pow((0.5+ix-nx-2*mm-npd)/npd,2))/2);
      }if(ix>=npd&&ix<(npd+nx+2*mm)&&iz>=npd&&iz<(npd+nz+2*mm)){
            coff1[id]=1.0;
	     coff2[id]=1.0;
      }
   }        
}
//a################################################################################
__global__ void cuda_record(float *p, float *seis, int *gx, int ng)//++++++++++++
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) seis[id]=p[gx[id]];
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
                               float _dz2, float _dx2,int nz, int nx, float *coff1, float *coff2)
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

    const int stride2 = nz + 2 * mm + 2 * npd;
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
    if ((iz >= nz + mm + 2*npd) ||(ix >= nx + mm + 2*npd)) validr = false;
    if ((iz >= nz + 2*npd) ||(ix >= nx + 2*npd)) validw = false;

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
       p0[outIndex]=coff1[outIndex]*2.0*p1[outIndex]-coff2[outIndex]*p0[outIndex]+vv[outIndex]*(c3)+vx[outIndex]*(c2);
       q0[outIndex]=coff1[outIndex]*2.0*q1[outIndex]-coff2[outIndex]*q0[outIndex]+vv[outIndex]*(c3)+vn[outIndex]*(c2);
     }
}
//a################################################################################
void velocity_transform(float *v0, float*vv, float*vx, float*vn, float dt, int nz, int nx, float *vmax)
 /*< velocit2 transform: vv=v0*dt; vv<--vv^2 >*/
{
  int i1, i2, nnz, nnx;
  float tmp;

  nnz=nz+2*mm+2*npd;
  nnx=nx+2*mm+2*npd;
  *vmax=v0[0];;
  // inner zone
    for(i2=0; i2<nx; i2++){//x
      for(i1=0; i1<nz; i1++){//z
       if(*vmax<v0[i1+nz*i2])*vmax=v0[i1+nz*i2];
	tmp=v0[i1+nz*i2]*dt;
	vv[(i1+mm+npd)+nnz*(i2+mm+npd)]=tmp*tmp;
	vx[(i1+mm+npd)+nnz*(i2+mm+npd)]=tmp*tmp*(1+2*0.3);
	vn[(i1+mm+npd)+nnz*(i2+mm+npd)]=tmp*tmp*(1+2*0.2);
      }
    }
    //top & down 
	for(i2=0; i2<nnx; i2++){//x
	    for (i1=0; i1<mm+npd; i1++){//z
		vv[i1+nnz*i2]=vv[mm+npd+nnz*i2];
		vv[(nnz-i1-1)+nnz*i2]=vv[(nnz-mm-npd-1)+nnz*i2];
		vx[i1+nnz*i2]=vx[mm+npd+nnz*i2];
		vx[(nnz-i1-1)+nnz*i2]=vx[(nnz-mm-npd-1)+nnz*i2];
		vn[i1+nnz*i2]=vn[mm+npd+nnz*i2];
		vn[(nnz-i1-1)+nnz*i2]=vn[(nnz-mm-npd-1)+nnz*i2];
	    }
	}
    //left & right
	for(i2=0; i2<mm+npd; i2++){//x
	    for (i1=0; i1<nnz; i1++){//z
		vv[i1+nnz*i2]=vv[i1+nnz*(mm+npd)];
		vv[i1+nnz*(nnx-i2-1)]=vv[i1+nnz*(nnx-mm-npd-1)];
		vx[i1+nnz*i2]=vx[i1+nnz*(mm+npd)];
		vx[i1+nnz*(nnx-i2-1)]=vx[i1+nnz*(nnx-mm-npd-1)];
		vn[i1+nnz*i2]=vn[i1+nnz*(mm+npd)];
		vn[i1+nnz*(nnx-i2-1)]=vn[i1+nnz*(nnx-mm-npd-1)];
	    }
	}
}
//a################################################################################
void window3d(float *a, float *b, int nz, int nx)
/*< window a 3d subvolume >*/
{
	int i1, i2, nnz;
	nnz=nz+2*mm+ 2*npd;//z
	
	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	{
          a[i1+nz*i2]=b[(i1+mm+npd)+nnz*(i2+mm+npd)];
	}
}
//a################################################################################
int main(int argc, char* argv[])
{

	int nz, nx, nnz, nnx, ns, nt, kt, it, is, fsz, fsx,  dsz, dsx, ng;
	int *s_szx,*s_gzx;
	float dz, dx,  fm, dt, _dz2, _dx2, vmax=0;
	float *v0, *vv, *vx, *vn, *s_wlt, *s_vv, *s_vx, *s_vn, *s_p0, *s_p1, *s_q0, *s_q1, *ptr;
       float *s_dcal, *s_dcal0, *s_dcal1;
       float *coff1, *coff2;
//a##################################################
	char FNvel[250]={"vel_600_300_const.dat"};
       char FNsnap[250]={"snap.dat"};
       char FNshot[250]={"shot.dat"};
//a##################################################
       FILE *fpvel,*fpsnap,*fpshot;
       fpvel=fopen(FNvel,"rb");
       fpsnap=fopen(FNsnap,"wb");
       fpshot=fopen(FNshot,"wb");
//a##################################################
       fm=40;     

    	nx=600;   dx=5;
    	nz=300;   dz=5;
    	
   	nt=1301;   kt=150;    dt=0.0005;

   	ns=5;
       fsx=100;dsx=100;
       fsz=150;dsz=0;
//a##################################################

	_dz2=1.0/(dz*dz);
	_dx2=1.0/(dx*dx);
	nnz=nz+2*mm+2*npd;
	nnx=nx+2*mm+2*npd;
       ng=nx;

    	v0=(float*)malloc(nz*nx*sizeof(float));
    	vv=(float*)malloc(nnz*nnx*sizeof(float));
    	vx=(float*)malloc(nnz*nnx*sizeof(float));
    	vn=(float*)malloc(nnz*nnx*sizeof(float));
    	s_dcal1=(float*)malloc(ng*nt*sizeof(float));//++++++++++

	fread(v0, sizeof(float), nz*nx, fpvel);// read velocit2 model v0
	velocity_transform(v0, vv, vx, vn, dt, nz, nx, &vmax);

    	cudaSetDevice(0);// initialize device, default device=0;
	check_gpu_error("Failed to initialize device!");

	dim3 dimg, dimb;
	dimg.x=(nz+2*npd+2*mm+BlockSize1-1)/BlockSize1;//BlockSize1=16;
	dimg.y=(nx+2*npd+2*mm+BlockSize2-1)/BlockSize2;//BlockSize2=16;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;

	/* allocate memory on device */
	cudaMalloc(&s_wlt, nt*sizeof(float));
	cudaMalloc(&s_vv, nnz*nnx*sizeof(float));
	cudaMalloc(&s_vx, nnz*nnx*sizeof(float));
	cudaMalloc(&s_vn, nnz*nnx*sizeof(float));
	cudaMalloc(&s_p0, nnz*nnx*sizeof(float));
	cudaMalloc(&s_p1, nnz*nnx*sizeof(float));
	cudaMalloc(&s_q0, nnz*nnx*sizeof(float));
	cudaMalloc(&s_q1, nnz*nnx*sizeof(float));
	cudaMalloc(&s_szx, ns*sizeof(int));
	cudaMalloc(&s_gzx, ng*sizeof(int));
	cudaMalloc(&s_dcal, ng*sizeof(float));	
	cudaMalloc(&s_dcal0, ng*nt*sizeof(float)); 

	cudaMalloc(&coff1, nnz*nnx*sizeof(float));
	cudaMalloc(&coff2, nnz*nnx*sizeof(float));
	cudaMemset(coff1, 0, nnz*nnx*sizeof(float));
	cudaMemset(coff2, 0, nnz*nnx*sizeof(float));

       cuda_initial_PML<<<(nnx*nnz+511)/512, 512>>>(coff1, coff2, nx, nz, dx, dz, dt, vmax);
	check_gpu_error("Failed to allocate memory for variables!");

	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(s_wlt, fm, dt, nt);
	cudaMemcpy(s_vv, vv, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(s_vx, vx, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(s_vn, vn, nnz*nnx*sizeof(float), cudaMemcpyHostToDevice);
	cuda_set_s<<<1, ns>>>(s_szx, fsz, fsx, dsz, dsx, ns, nz, nx);
	cuda_set_g<<<(ng+511)/512,512>>>(s_gzx, ng, nz, nx);

	clock_t t0, t1;
	t0 = clock();
     for(is=0; is<ns; is++){printf("\nIS=%2d << it=",is);
	  cudaMemset(s_p0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_p1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_q0, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_q1, 0, nnz*nnx*sizeof(float));
	  cudaMemset(s_dcal, 0, ng*sizeof(float));
	  cudaMemset(s_dcal0, 0, ng*nt*sizeof(float));
	  for(it=0; it<nt; it++){
           if(it%kt==0)  printf("%d-",it);
	    cuda_add_source<<<1,1>>>(true, s_p1, &s_wlt[it], &s_szx[is], 1);
	    cuda_add_source<<<1,1>>>(true, s_q1, &s_wlt[it], &s_szx[is], 1);
	    cuda_step_fd2d<<<dimg,dimb>>>(s_p0, s_p1, s_q0, s_q1, s_vv, s_vx, s_vn, _dz2, _dx2, nz, nx, coff1, coff2);
	    ptr=s_p0; s_p0=s_p1; s_p1=ptr;
	    ptr=s_q0; s_q0=s_q1; s_q1=ptr;

	    cuda_record<<<(ng+511)/512, 512>>>(s_p0, s_dcal, s_gzx, ng);
           cuda_trans_x2tx<<<(ng+511)/512, 512>>>(s_dcal, s_dcal0, it, nt, ng);

	    if(is==0&&it!=0&&it%kt==0){
	      cudaMemcpy(vv, s_p0, nnz*nnx*sizeof(float), cudaMemcpyDeviceToHost);
	      window3d(v0, vv, nz, nx);
	      fwrite(v0, sizeof(float),nz*nx, fpsnap);	  
             }
	  }printf(" >>");

	  cudaMemcpy(s_dcal1, s_dcal0, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
         fseek(fpshot,is*ng*nt*sizeof(float),0);
	  fwrite(s_dcal1, sizeof(float), ng*nt, fpshot);
      }
     t1 = clock();
     printf("total %d shots: %f (s)\n", ns, ((float)(t1-t0))/CLOCKS_PER_SEC);

	/* free memory on device */
	cudaFree(s_wlt);
	cudaFree(s_vv);
	cudaFree(s_vx);
	cudaFree(s_vn);
	cudaFree(s_p0);
	cudaFree(s_p1);
	cudaFree(s_q0);
	cudaFree(s_q1);
	cudaFree(s_szx);
	cudaFree(s_gzx);
	cudaFree(s_dcal);
	cudaFree(s_dcal0);

	free(v0);
	free(vv);
	free(vx);
	free(vn);
	free(s_dcal1);

    	exit (0);
}
