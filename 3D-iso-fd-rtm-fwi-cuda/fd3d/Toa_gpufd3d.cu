/* GPU-based finite difference on 3-D grid
*/
/*
  Copyright (C) 2014  Xi'an Jiaotong University (Pengliang Yang)

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  Reference:
    Micikevicius, Paulius. "3D finite difference computation on GPUs
    using CUDA." Proceedings of 2nd Workshop on General Purpose 
    Processing on Graphics Processing Units. ACM, 2009.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include <rsf.h>
}

#ifndef PI
#define PI 	SF_PI
#endif
#define BlockSize1 16// tile size in 1st-axis
#define BlockSize2 16// tile size in 2nd-axis
#define mm      4    // half of the order in space
#define npd     50   // absorbing boundry condition wield

void sf_check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); 
	exit(0);   
    }
}

__constant__ float stencil[mm+1]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};

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



__global__ void cuda_set_s(int *szxy, int szbeg, int sxbeg, int sybeg, int jsz, int jsx, int jsy, int ns, int nz, int nx, int ny)
/*< set the positions of sources  in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
    	if (id<ns) szxy[id]=(szbeg+id*jsz+mm+npd)+nnz*(sxbeg+id*jsx+mm+npd)+nnz*nnx*(sybeg+id*jsy+mm+npd);
}

__global__ void cuda_set_g(int *gzxy, int ng, int nz, int nx, int ny)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
       int iy=id/nx;
       int ix=id%nx;
    	if (id<ng) gzxy[id]=(mm+npd)+nnz*(ix*1+mm+npd)+nnz*nnx*(iy*1+mm+npd);
       
}
__global__ void cuda_trans_xy2txy(float *xy, float *txy, int it, int nt, int ng)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) txy[it+id*nt]+=xy[id];
       
}

__global__ void cuda_absorb_bndr(float *d_p,int nz,int nx,int ny,float qp)
/*< absorb boundry condition >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
	int nny=ny+2*mm+2*npd;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
            /*< front & back (0<y<ny) >*/
             if ( iy < npd )
               d_p[id]=( qp*pow((npd-iy)/(1.0*npd),2) + 1 )*d_p[id];
             else if ( iy >= 2*mm + npd + ny )
               d_p[id]=( qp*pow((iy-2*mm-npd-ny)/(1.0*npd),2) + 1 )*d_p[id];
            /*< left & right (0<x<nx) >*/
             if ( ix < npd )
               d_p[id]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*d_p[id];
             else if ( ix >= 2*mm + npd + nx )
               d_p[id]=( qp*pow((ix-2*mm-npd-nx)/(1.0*npd),2) + 1 )*d_p[id];
            /*< up & down (0<z<nz) >*/
             if ( iz < npd )
               d_p[id]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*d_p[id];
             else if ( iz >= 2*mm + npd + nz )
               d_p[id]=( qp*pow((iz-2*mm-npd-nz)/(1.0*npd),2) + 1 )*d_p[id];
        }
       
}
__global__ void cuda_record(float *p, float *seis, int *gxz, int ng)//++++++++++++
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) seis[id]=p[gxz[id]];
}


__global__ void cuda_add_source(bool add, float *p, float *source, int *szxy, int ns)
/*< add/subtract sources: length of source[]=ns, index stored in szxy[] >*/
{
  int id=threadIdx.x+blockIdx.x*blockDim.x;

  if(id<ns){
    if(add){
      p[szxy[id]]+=source[id];
    }else{
      p[szxy[id]]-=source[id];
    }
  }
}

__global__ void cuda_step_fd3d(float *p0, float *p1, float *vv, float _dz2, float _dx2, float _dy2, int n1, int n2, int n3)
/*< step forward: 3-D FD, order=8 >*/
{
    bool validr = true;
    bool validw = true;
    const int gtid1 = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int gtid2 = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix
    const int ltid1 = threadIdx.x;//ithreadz
    const int ltid2 = threadIdx.y;//ithreadx
    const int work1 = blockDim.x;//nblockz
    const int work2 = blockDim.y;//nblockx
    __shared__ float tile[BlockSize2 + 2 * mm][BlockSize1 + 2 * mm];//tile[16+2*mm][16+2*mm]

    const int stride2 = n1 + 2 * mm + 2 * npd;//n1=nz
    const int stride3 = stride2 * (n2 + 2 * mm + 2 * npd);//n2=nx   stride3=(nz+2*mm)*(nx+2*mm)

    int inIndex = 0;
    int outIndex = 0;

    // Advance inputIndex to start of inner volume
    inIndex += (mm ) * stride2 + mm ;// inIndex=mm*(nz+2*mm+2*npd)+mm;

    // Advance inputIndex to target element
    inIndex += gtid2 * stride2 + gtid1; // inIndex=mm*(nz+2*mm)+mm+ix*(nz+2*mm+2*npd)+iz;:igrid

    float infront[mm];
    float behind[mm];
    float current;

    const int t1 = ltid1 + mm;
    const int t2 = ltid2 + mm;

    // Check in bounds
    if ((gtid1 >= n1 + mm + 2*npd) ||(gtid2 >= n2 + mm + 2*npd)) validr = false;
    if ((gtid1 >= n1 + 2*npd) ||(gtid2 >= n2 + 2*npd)) validw = false;

    // Preload the "infront" and "behind" data
    for (int i = mm -2 ; i >= 0 ; i--)//change 'mm-2' to 'mm-1'+++++++++++++++++++
    {
        if (validr) behind[i] = p1[inIndex];
        inIndex += stride3;//stride3=(nz+2*mm)*(nx+2*mm)
    }

    if (validr)	current = p1[inIndex];

    outIndex = inIndex;
    inIndex += stride3;//stride3=(nz+2*mm)*(nx+2*mm)

    for (int i = 0 ; i < mm ; i++)
    {
	if (validr) infront[i] = p1[inIndex];
        inIndex += stride3;//stride3=(nz+2*mm)*(nx+2*mm)
    }

    // Step through the zx-planes

    for (int i3 = mm ; i3 < n3 + 2*npd + mm ; i3++)
    {
        // Advance the slice (move the thread-front)
        for (int i = mm - 1 ; i > 0 ; i--) behind[i] = behind[i - 1];

        behind[0] = current;
        current = infront[0];

        for (int i = 0 ; i < mm - 1 ; i++) infront[i] = infront[i + 1];

        if (validr) infront[mm - 1] = p1[inIndex];

        inIndex += stride3;
        outIndex += stride3;
        __syncthreads();

        // Update the data slice in the local tile
        // Halo above & below
        if (ltid2 < mm)
        {
          /*   tile[ithread][ithread+mm]=p1[igrid - mm*(nz+2*mm)]  */
            tile[ltid2][t1]                  = p1[outIndex - mm * stride2];//t1 = ltid1 + mm;
            tile[ltid2 + work2 + mm][t1] = p1[outIndex + work2 * stride2];
        }

        // Halo left & right
        if (ltid1 < mm)
        {
            tile[t2][ltid1]                  = p1[outIndex - mm];
            tile[t2][ltid1 + work1 + mm] = p1[outIndex + work1];
        }

        tile[t2][t1] = current;
        __syncthreads();

        // Compute the output value
		float c1, c2, c3;
			c1=c2=c3=stencil[0]*current;        

        for (int i=1; i <= mm ; i++)
        {
			c1 +=stencil[i]*(tile[t2][t1-i]+ tile[t2][t1+i]);//z
			c2 +=stencil[i]*(tile[t2-i][t1]+ tile[t2+i][t1]);//x
			c3 +=stencil[i]*(infront[i-1]  + behind[i-1]  ); //y
        }
			c1*=_dz2;	
			c2*=_dx2;
			c3*=_dy2;
	
	
        if (validw) p0[outIndex]=2.0*p1[outIndex]-p0[outIndex]+vv[outIndex]*(c1+c2+c3);
    }

}

void velocity_transform(float *v0, float*vv, float dt, int n1, int n2, int n3)
 /*< velocit2 transform: vv=v0*dt; vv<--vv^2 >*/
{
  int i1, i2, i3, nn1, nn2, nn3;
  float tmp;

  nn1=n1+2*mm+2*npd;
  nn2=n2+2*mm+2*npd;
  nn3=n3+2*mm+2*npd;

  // inner zone
  for(i3=0; i3<n3; i3++){//y
    for(i2=0; i2<n2; i2++){//x
      for(i1=0; i1<n1; i1++){//z
	tmp=v0[i1+n1*i2+n1*n2*i3]*dt;
	vv[(i1+mm+npd)+nn1*(i2+mm+npd)+nn1*nn2*(i3+mm+npd)]=tmp*tmp;
      }
    }
  }  
    //top & down 
    for(i3=0; i3<nn3; i3++){//y
	for(i2=0; i2<nn2; i2++){//x
	    for (i1=0; i1<mm+npd; i1++){//z
		vv[i1+nn1*i2+nn1*nn2*i3]=vv[mm+npd+nn1*i2+nn1*nn2*i3];
		vv[(nn1-i1-1)+nn1*i2+nn1*nn2*i3]=vv[(nn1-mm-npd-1)+nn1*i2+nn1*nn2*i3];
	    }
	}
    }

    //left & right
    for(i3=0; i3<nn3; i3++){//y
	for(i2=0; i2<mm+npd; i2++){//x
	    for (i1=0; i1<nn1; i1++){//z
		vv[i1+nn1*i2+nn1*nn2*i3]=vv[i1+nn1*(mm+npd)+nn1*nn2*i3];
		vv[i1+nn1*(nn2-i2-1)+nn1*nn2*i3]=vv[i1+nn1*(nn2-mm-npd-1)+nn1*nn2*i3];
	    }
	}
    }
    //front & back
    for(i3=0; i3<mm+npd; i3++){//y
	for(i2=0; i2<nn2; i2++){//x
	    for(i1=0; i1<nn1; i1++){//z
		vv[i1+nn1*i2+nn1*nn2*i3]=vv[i1+nn1*i2+nn1*nn2*(mm+npd)];
		vv[i1+nn1*i2+nn1*nn2*(nn3-1-i3)]=vv[i1+nn1*i2+nn1*nn2*(nn3-mm-npd-1)];
	    }
	}
    }
}


void window3d(float *a, float *b, int n1, int n2, int n3)
/*< window a 3d subvolume >*/
{
	int i1, i2, i3, nn1, nn2;
	nn1=n1+2*mm+ 2*npd;//z
	nn2=n2+2*mm+ 2*npd;//x
	
	for(i3=0; i3<n3; i3++)
	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	{
		a[i1+n1*i2+n1*n2*i3]=b[(i1+mm+npd)+nn1*(i2+mm+npd)+nn1*nn2*(i3+mm+npd)];
	}
}


int main(int argc, char* argv[])
{
	bool verb;
	int nz, nx, ny, nnz, nnx, nny, ns, nt, kt, it, is, szbeg, sxbeg, sybeg, jsz, jsx, jsy, ng;
	int *d_szxy,*d_gzxy;
	float dz, dx, dy, fm, dt, _dz2, _dx2, _dy2;
	float *v0, *vv, *d_wlt, *d_vv, *d_p0, *d_p1, *ptr;
       float *d_dcal,*d_dcal0,*d_dcal1;//+++++++++++++

	sf_file Fv, Fw, Fshot;

    	sf_init(argc,argv);
	Fv=sf_input("in");
	Fw=sf_output("snap");
	Fshot=sf_output("shot");

    	if (!sf_getbool("verb",&verb)) verb=false; /* verbosit2 */
    	if (!sf_histint(Fv,"n1",&nz)) sf_error("No n1= in input");
    	if (!sf_histint(Fv,"n2",&nx)) sf_error("No n2= in input");
    	if (!sf_histint(Fv,"n3",&ny)) sf_error("No n3= in input");
    	if (!sf_histfloat(Fv,"d1",&dz)) sf_error("No d1= in input");
    	if (!sf_histfloat(Fv,"d2",&dx)) sf_error("No d2= in input");
    	if (!sf_histfloat(Fv,"d3",&dy)) sf_error("No d3= in input");
   	if (!sf_getint("nt",&nt))  sf_error("nt required");
	/* total number of time steps */
    	if (!sf_getint("kt",&kt)) sf_error("kt required");
	/* record wavefield at time kt */
   	if (!sf_getfloat("dt",&dt))  sf_error("dt required");
	/* time sampling interval */
   	if (!sf_getfloat("fm",&fm))  fm=20;
	/* dominant frequency of Ricker wavelet */
   	if (!sf_getint("ns",&ns))  ns=1;
	/* number of sources */
	if (!sf_getint("szbeg",&szbeg)) sf_error("No szbeg");
	/* source beginning of z-axis */
	if (!sf_getint("sxbeg",&sxbeg)) sf_error("No sxbeg");
	/* source beginning of x-axis */
	if (!sf_getint("sybeg",&sybeg)) sf_error("No sybeg");
	/* source beginning of y-axis */
	if (!sf_getint("jsz",&jsz)) sf_error("No jsz");
	/* source jump interval in z-axis */
	if (!sf_getint("jsx",&jsx)) sf_error("No jsx");
	/* source jump interval in x-axis */
	if (!sf_getint("jsy",&jsy)) sf_error("No jsy");
	/* source jump interval in y-axis */

	sf_putint(Fw,"n1",nz);
	sf_putint(Fw,"n2",nx);
	sf_putint(Fw,"n3",ny);



	_dz2=1.0/(dz*dz);
	_dx2=1.0/(dx*dx);
	_dy2=1.0/(dy*dy);

	nnz=nz+2*mm+2*npd;
	nnx=nx+2*mm+2*npd;
	nny=ny+2*mm+2*npd;
       //ng=nx*ny;//+++++++++++++
       ng=nx*ny;

	sf_putint(Fshot,"n1",nt);
	sf_putint(Fshot,"n2",nx);
	sf_putint(Fshot,"n3",ny);

    	v0=(float*)malloc(nz*nx*ny*sizeof(float));
    	vv=(float*)malloc(nnz*nnx*nny*sizeof(float));
    	d_dcal1=(float*)malloc(ng*nt*sizeof(float));//++++++++++

	sf_floatread(v0, nz*nx*ny, Fv);// read velocit2 model v0
	velocity_transform(v0, vv, dt, nz, nx, ny);// init

    	cudaSetDevice(0);// initialize device, default device=0;
	sf_check_gpu_error("Failed to initialize device!");

	dim3 dimg, dimb;
	dimg.x=(nz+2*npd+2*mm+BlockSize1-1)/BlockSize1;//BlockSize1=16;
	dimg.y=(nx+2*npd+2*mm+BlockSize2-1)/BlockSize2;//BlockSize2=16;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;

	/* allocate memory on device */
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_vv, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_p0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_p1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_szxy, ns*sizeof(int));
	cudaMalloc(&d_gzxy, ng*sizeof(int));//++++++++++
	cudaMalloc(&d_dcal, ng*sizeof(float));	/* calculated/synthetic seismic data *///++++++++++++++
	cudaMalloc(&d_dcal0, ng*nt*sizeof(float)); //++++++++++++++


	sf_check_gpu_error("Failed to allocate memory for variables!");

	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(d_wlt, fm, dt, nt);
	cudaMemcpy(d_vv, vv, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
	cuda_set_s<<<1, ns>>>(d_szxy, szbeg, sxbeg, sybeg, jsz, jsx, jsy, ns, nz, nx, ny);
	cuda_set_g<<<(ng+511)/512,512>>>(d_gzxy, ng, nz, nx, ny);//++++++++++++geophone

	float mstimer;
	clock_t t0, t1;
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);
	for(is=0; is<ns; is++){
	  cudaEventRecord(start);

	  cudaMemset(d_p0, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(d_p1, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(d_dcal, 0, ng*sizeof(float));//++++++++++
	  cudaMemset(d_dcal0, 0, ng*nt*sizeof(float));//++++++++++
	    //  sf_floatwrite(vv, nnz*nnx*nny, Fw);


	  for(it=0; it<nt; it++){
	    cuda_add_source<<<1,1>>>(true, d_p1, &d_wlt[it], &d_szxy[is], 1);
	    cuda_step_fd3d<<<dimg,dimb>>>(d_p0, d_p1, d_vv, _dz2, _dx2, _dy2, nz, nx, ny);//<<<[(nz+16-1)/16,(nx+16-1)/16],[16,16]>>>
	    ptr=d_p0; d_p0=d_p1; d_p1=ptr;//toggle buffers

           cuda_absorb_bndr<<<dimg,dimb>>>(d_p0, nz, nx, ny, -0.25);
           cuda_absorb_bndr<<<dimg,dimb>>>(d_p1, nz, nx, ny, -0.25);

	    cuda_record<<<(ng+511)/512, 512>>>(d_p0, d_dcal, d_gzxy, ng);//++++++++++
           cuda_trans_xy2txy<<<(ng+511)/512, 512>>>(d_dcal, d_dcal0, it, nt, ng);



	    if(it==kt){
	      t0 = clock();

	      cudaMemcpy(vv, d_p0, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	      window3d(v0, vv, nz, nx, ny);
	      sf_floatwrite(v0, nz*nx*ny, Fw);	  

 	      t1 = clock();
 	      sf_warning("save the volume: %f (s)", ((float)(t1-t0))/CLOCKS_PER_SEC); 	
	    }

          if(it%100==0)
	     sf_warning("is=%2d, it=%d",is,it);
	  }

	    cudaMemcpy(d_dcal1, d_dcal0, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
	    sf_floatwrite(d_dcal1, ng*nt, Fshot);

	  cudaEventRecord(stop);
          cudaEventSynchronize(stop);
  	  cudaEventElapsedTime(&mstimer, start, stop);
    	  sf_warning("%d shot finished: %g (s)",is+1, mstimer*1.e-3);
	}
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	/* free memory on device */
	cudaFree(d_wlt);
	cudaFree(d_vv);
	cudaFree(d_p0);
	cudaFree(d_p1);
	cudaFree(d_szxy);
	cudaFree(d_gzxy);
	cudaFree(d_dcal);//+++++++++++
	cudaFree(d_dcal0);//+++++++++++
	free(v0);
	free(vv);
	free(d_dcal1);

    	exit (0);
}
