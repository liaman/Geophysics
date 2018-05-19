//############################################
//#
//#
//#
//#
//##############################################
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

#define PI 3.141592653

#define EPS 0.000000001

#define BlockSize1 16// tile size in 1st-axis
#define BlockSize2 16// tile size in 2nd-axis
#define BlockSize 512

#define mm      4    // half of the order in space
#define npd     20   // absorbing boundry condition wield


__device__ float s, t, r;
//a#############################################################################################
__constant__ float stencil[mm+1]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};
//a#############################################################################################
__global__ void cuda_step_fd3d(float *p0, float *p1, float *VV, float _dz2, float _dx2, float _dy2, int n1, int n2, int n3, 
                               float dt, float *pdt2, bool pdt)
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
	 if (validw&&pdt) pdt2[outIndex]= (c1+c2+c3);
	
        if (validw) p0[outIndex]=2.0*p1[outIndex]-p0[outIndex]+VV[outIndex]*VV[outIndex]*dt*dt*(c1+c2+c3);
    }

}

//a#############################################################################################
void check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	printf ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); 
	exit(0);   
    }
}
//a#############################################################################################
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
//a#############################################################################################
void velocity_transform(float *v, float*vv, float dt, int n1, int n2, int n3)
 /*< velocit2 transform: vv=v*dt; vv<--vv^2 >*/
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
	tmp=v[i1+n1*i2+n1*n2*i3];
	vv[(i1+mm+npd)+nn1*(i2+mm+npd)+nn1*nn2*(i3+mm+npd)]=tmp;
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
//a#############################################################################################
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
//a#############################################################################################
__global__ void cuda_record(float *P, float *seis, int *gxz, int ng, int it, int nt, bool record)//++++++++++++
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng)
        {
        if(record) seis[it+id*nt]=P[gxz[id]];
        else  P[gxz[id]]=seis[it+id*nt];
        }
}
//a#############################################################################################
__global__ void cuda_cal_illum(float *s, float *p, int nz, int nx, int ny)
/*< calculate the source lighting matrix >*/
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
          if(id<nnz*nnx*nny) s[id]+=p[id]*p[id];
        }    
}
//a#############################################################################################
__global__ void cuda_illum(float *g1, float *illum, int nz, int nx, int ny)
/*< source lighting >*/
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
          if(id<nnz*nnx*nny&&illum[id]!=0) g1[id]/=illum[id];
        }    
}
//a#############################################################################################
__global__ void cuda_sum(float *ns, float *is, int nz, int nx, int ny)
/*< source lighting >*/
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
          if(id<nnz*nnx*nny) ns[id]+=is[id];
        }    
}
//a#############################################################################################
__global__ void cuda_cal_g1(float *g1, float *s, float *g, int nz, int nx, int ny)
/*< calculate is g1 >*/
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
          if(id<nnz*nnx*nny) g1[id]+=s[id]*g[id];
        }    
}
//a#############################################################################################
__global__ void cuda_absorb_bndr(float *P,float *Q,int nz,int nx,int ny,float qp)
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
             if ( iy < npd ){
               P[id]=( qp*pow((npd-iy)/(1.0*npd),2) + 1 )*P[id];
               Q[id]=( qp*pow((npd-iy)/(1.0*npd),2) + 1 )*Q[id];
             }else if ( iy >= 2*mm + npd + ny ){
               P[id]=( qp*pow((iy-2*mm-npd-ny)/(1.0*npd),2) + 1 )*P[id];
               Q[id]=( qp*pow((iy-2*mm-npd-ny)/(1.0*npd),2) + 1 )*Q[id];
               }
            /*< left & right (0<x<nx) >*/
             if ( ix < npd ){
               P[id]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*P[id];
               Q[id]=( qp*pow((npd-ix)/(1.0*npd),2) + 1 )*Q[id];
             }else if ( ix >= 2*mm + npd + nx ){
               P[id]=( qp*pow((ix-2*mm-npd-nx)/(1.0*npd),2) + 1 )*P[id];
               Q[id]=( qp*pow((ix-2*mm-npd-nx)/(1.0*npd),2) + 1 )*Q[id];
               }
            /*< up & down (0<z<nz) >*/
             if ( iz < npd ){
               P[id]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*P[id];
               Q[id]=( qp*pow((npd-iz)/(1.0*npd),2) + 1 )*Q[id];
             }else if ( iz >= 2*mm + npd + nz ){
               P[id]=( qp*pow((iz-2*mm-npd-nz)/(1.0*npd),2) + 1 )*P[id];
               Q[id]=( qp*pow((iz-2*mm-npd-nz)/(1.0*npd),2) + 1 )*Q[id];
               }
        }    
}
//a#############################################################################################
__global__ void cuda_cal_residuals(float *obj, float *cal, float *obs, float *com, int nn, int nx, int ny, int nt)
{
    const int it = blockIdx.x * blockDim.x + threadIdx.x;//0--nt's thread:it
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix
	int id, iy;
       if(it<nt){
        for(iy=0;iy<ny;iy++) 
         {
            id=it+ix*nt+iy*nt*nx;
            if (id<nn) 
              com[id]=cal[id] - obs[id];
              *obj+=com[id]*com[id];
         }
       } 
    	

}
//a#############################################################################################
__global__ void cuda_ricker_wavelet(float *wlt, float favg, float dt, int nt, float pfac)
/*< generate ricker wavelet with time deley >*/
{
	int it=threadIdx.x+blockDim.x*blockIdx.x;
    	if (it<nt){
	  float tmp = PI*favg*fabsf(it*dt-1.0/favg);//delay the wavelet to exhibit all waveform
	  tmp *=tmp;
	  wlt[it]= (1.0-2.0*tmp)*expf(-tmp);// ricker wavelet at time: t=nt*dt
	}
}
//a#############################################################################################
__global__ void cuda_set_s(int *szxy, int fsz, int fsx, int fsy, int dsz, int dsx, int dsy, int ns, int nsx, int nz, int nx, int ny)
/*< set the positions of sources  in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;

       int ixs=id%nsx;
       int iys=id/nsx;

    	if (id<ns) szxy[id]=(fsz+mm+npd)+nnz*(fsx+ixs*dsx+mm+npd)+nnz*nnx*(fsy+iys*dsy+mm+npd);
}
//a#############################################################################################
__global__ void cuda_set_up_do(int *gzxy, int *up, int *down, int ng, int nz, int nx, int ny)
/*< set the positions of  geophones & down in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
       int iy=id/nx;
       int ix=id%nx;
    	if (id<ng){
           gzxy[id]=(mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(iy+mm+npd);
           up[id]=(mm+npd-1)+nnz*(ix+mm+npd)+nnz*nnx*(iy+mm+npd);
           down[id]=(nz+mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(iy+mm+npd);
        }
       
}
//a#############################################################################################
__global__ void cuda_set_fr_ba(int *front, int *back, int ng, int nz, int nx, int ny)
/*< set the positions of  front & back in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
       int ix=id/nz;
       int iz=id%nz;
    	if (id<ng){
           front[id]=(iz+mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(mm+npd-1);
           back[id]=(iz+mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(ny+mm+npd);
        }
}
//a#############################################################################################
__global__ void cuda_set_le_ri(int *left, int *right,int ng, int nz,int nx, int ny)
/*< set the positions of  left & right in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
       int iy=id/nz;
       int iz=id%nz;
    	if (id<ng){
           left[id]=(iz+mm+npd)+nnz*(mm+npd-1)+nnz*nnx*(iy+mm+npd);
           right[id]=(iz+mm+npd)+nnz*(nx+mm+npd)+nnz*nnx*(iy+mm+npd);
        }
}
//a#############################################################################################
__global__ void cuda_save_bndr(float *bndr, float *p0, int *front, int *back, int *left, int *right, int *up, int *down, 
                               int nz, int nx, int ny, bool write)//(2*nz*nx+2*nz*ny+nx*ny)
/*< write boundaries out or read them into wavefield variables p>*/
{
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	if(write){
		if(id<nz*nx) 
                    bndr[id]=p0[front[id]];                /* front boundary */
		else if((id>=nz*nx)&&(id<2*nz*nx)) 
                    bndr[id]=p0[back[id-nz*nx]];           /* back  boundary */
		else if((id>=2*nz*nx)&&(id<(2*nz*nx+nz*ny))) 
                    bndr[id]=p0[left[id-2*nz*nx]];         /* left  boundary */
		else if((id>=(2*nz*nx+nz*ny))&&(id<(2*nz*nx+2*nz*ny))) 
                    bndr[id]=p0[right[id-2*nz*nx-nz*ny]];  /* right boundary */
		else if((id>=(2*nz*nx+2*nz*ny))&&(id<(2*nz*nx+2*nz*ny+nx*ny))) 
                    bndr[id]=p0[up[id-2*nz*nx-2*nz*ny]];   /* up    boundary */
		else if((id>=(2*nz*nx+2*nz*ny+nx*ny))&&(id<(2*nz*nx+2*nz*ny+2*nx*ny))) 
                    bndr[id]=p0[down[id-2*nz*nx-2*nz*ny-nx*ny]];/* down boundary */
	}else{
		if(id<nz*nx) 
                    p0[front[id]]=bndr[id];                /* front boundary */
		else if((id>=nz*nx)&&(id<2*nz*nx)) 
                    p0[back[id-nz*nx]]=bndr[id];           /* back  boundary */
		else if((id>=2*nz*nx)&&(id<(2*nz*nx+nz*ny))) 
                    p0[left[id-2*nz*nx]]=bndr[id];         /* left  boundary */
		else if((id>=(2*nz*nx+nz*ny))&&(id<(2*nz*nx+2*nz*ny))) 
                    p0[right[id-2*nz*nx-nz*ny]]=bndr[id];  /* right boundary */
		else if((id>=(2*nz*nx+2*nz*ny))&&(id<(2*nz*nx+2*nz*ny+nx*ny))) 
                    p0[up[id-2*nz*nx-2*nz*ny]]=bndr[id];   /* up    boundary */
		else if((id>=(2*nz*nx+2*nz*ny+nx*ny))&&(id<(2*nz*nx+2*nz*ny+2*nx*ny))) 
                    p0[down[id-2*nz*nx-2*nz*ny-nx*ny]]=bndr[id];   /* down boundary */
	}
}
//a#############################################################################################
__global__ void cuda_scale_gradient(float *g1, float *VV, float *illum, int nnx, int nny, int nnz, bool precon)
/*< scale g1 >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz)
           {
		float a=VV[id];
		if (precon) a*=sqrtf(illum[id]+EPS);/*precondition with residual wavefield illum*/
		g1[id]*=2.0/a;
           }

        }  

}
//a#############################################################################################
__global__ void cuda_bell_smoothz(float *g1, int rbell, int nnx, int nny, int nnz)
/*< smoothing with gaussian function >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int i,id,iy;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz)
           {
		float s=0.0;
              for(i=-rbell; i<=rbell; i++) if(iz+i>=0 && iz+i<nnz) s+=expf(-(2.0*i*i)/rbell)*g1[id+i];
              g1[id]=s;
           }

        }  
}
//a#############################################################################################
__global__ void cuda_bell_smoothx(float *g1, int rbell, int nnx, int nny, int nnz)
/*< smoothing with gaussian function >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int i,id,iy;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz)
           {
		float s=0.0;
              for(i=-rbell; i<=rbell; i++) if(ix+i>=0 && ix+i<nnx) s+=expf(-(2.0*i*i)/rbell)*g1[id+i*nnz];
              g1[id]=s;
           }

        }  
}
//a#############################################################################################
__global__ void cuda_bell_smoothy(float *g1, int rbell, int nnx, int nny, int nnz)
/*< smoothing with gaussian function >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int i,id,iy;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz)
           {
		float s=0.0;
              for(i=-rbell; i<=rbell; i++) if(iy+i>=0 && iy+i<nny) s+=expf(-(2.0*i*i)/rbell)*g1[id+i*nnz*nnx];
              g1[id]=s;
           }

        }  
}
//a#############################################################################################
__global__ void cuda_cal_beta_step1(float *g0, float *g1, float *cg, int nnx, int nny, int nnz)
/*< calculate beta for nonlinear conjugate gradient algorithm 
configuration requirement: <<<1,BlockSize>>> >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

       s=0.0,t=0.0,r=0.0;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz) 
           {
              float a=g0[id];
		float b=g1[id];
		float c=cg[id];

		/* HS: Hestenses-Stiefel NLCG algorithm */
		s += b*(b-a);	// numerator of HS
		t += c*(b-a);	// denominator of HS,DY
		r += b*b;	// numerator of DY

           }
        }  
}
//a#############################################################################################
__global__ void cuda_cal_beta_step2(float *beta, int nnx, int nny, int nnz)
/*< set the positions of  geophones & down in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;

       if(id<1)
       {
	   float beta_HS=0.0;
	   float beta_DY=0.0;
	   if(t!=0) 
	   {
	 	beta_HS=s/t; 
		beta_DY=r/t;
	   } 
	   *beta=max(0.0, min(beta_HS, beta_DY));/* Hybrid HS-DY method combined with iteration restart */
        }
}
//a#############################################################################################
__global__ void cuda_cal_conjgrad(float *g1, float *cg, float beta, int nnx, int nny, int nnz)
/*< calculate nonlinear conjugate gradient >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;
       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz) 
           {
            cg[id] = -g1[id]+beta*cg[id];
           }
        }  
}
//a#############################################################################################
__global__ void cuda_cal_epsilon(float *VV, float *cg, float *epsil, int N)
/*< calculate estimated stepsize (epsil) according to Taratola's method
configuration requirement: <<<1, Block_Size>>> >*/ 
{
    	__shared__ float sdata[BlockSize];/* find max(|vv(:)|) */
	__shared__ float tdata[BlockSize];/* find max(|cg(:)|) */
    	int tid = threadIdx.x;
    	sdata[tid] = 0.0f;
    	tdata[tid] = 0.0f;
	for(int s=0; s<(N+BlockSize-1)/BlockSize; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		float a=(id<N)?fabsf(VV[id]):0.0f;
		float b=(id<N)?fabsf(cg[id]):0.0f;
		sdata[tid]= max(sdata[tid], a);
		tdata[tid]= max(tdata[tid], b);
	} 
    	__syncthreads();

    	/* do reduction in shared mem */
    	for(int s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s)	{sdata[tid]=max(sdata[tid], sdata[tid+s]);tdata[tid]=max(tdata[tid], tdata[tid+s]);} 
		__syncthreads();
    	}  
   	if (tid < 32)
   	{
		if (blockDim.x >=  64) { sdata[tid] =max(sdata[tid],sdata[tid + 32]);tdata[tid]=max(tdata[tid], tdata[tid+32]);}
		if (blockDim.x >=  32) { sdata[tid] =max(sdata[tid],sdata[tid + 16]);tdata[tid]=max(tdata[tid], tdata[tid+16]);}
		if (blockDim.x >=  16) { sdata[tid] =max(sdata[tid],sdata[tid + 8]);tdata[tid]=max(tdata[tid], tdata[tid+8]);}
		if (blockDim.x >=   8) { sdata[tid] =max(sdata[tid],sdata[tid + 4]);tdata[tid]=max(tdata[tid], tdata[tid+4]);}
		if (blockDim.x >=   4) { sdata[tid] =max(sdata[tid],sdata[tid + 2]);tdata[tid]=max(tdata[tid], tdata[tid+2]);}
		if (blockDim.x >=   2) { sdata[tid] =max(sdata[tid],sdata[tid + 1]);tdata[tid]=max(tdata[tid], tdata[tid+1]);}
    	}

    	if (tid == 0) { if(tdata[0]>EPS) *epsil=0.01*sdata[0]/tdata[0]; else *epsil=0.0;}
}
//a#############################################################################################
__global__ void cuda_com2derr(float *com, float *derr, int nx, int ny, int nt)
{
    const int it = blockIdx.x * blockDim.x + threadIdx.x;//0--nt's thread:it
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix
	int id, iy;
       if(it<nt){
        for(iy=0;iy<ny;iy++) 
         {
            id=it+ix*nt+iy*nt*nx;
            if (id<nx*ny*nt) derr[id]=com[id];
         }
       } 
}
//a#############################################################################################
__global__ void cuda_cal_vtmp(float *VVtmp, float *VV, float *cg, float epsil, int nnx, int nny, int nnz)
/*< calculate temporary velocity >*/ 
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;
       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz) 
           {
             VVtmp[id] =VV[id] + epsil*cg[id];
           }
        }  
}
//a#############################################################################################
__global__ void cuda_sum_alpha12(float *alpha1, float *alpha2, float *cal, float *obs, float *derr, 
                                 int nx, int ny, int nz, int nt)
{
    const int it = blockIdx.x * blockDim.x + threadIdx.x;//0--nt's thread:it
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix
	int id, iy;
       if(it<nt)
       {
        for(iy=0;iy<ny;iy++) 
         {
            id=it+ix*nt+iy*nt*nx;
            if (id<nx*ny*nt) 
              {
                float c=derr[id];
                float a=obs[id]+c;/* since f(mk)-dobs[id]=derr[id], thus f(mk)=b+c; */
                float b=cal[id]-a;/* f(mk+epsil*cg)-f(mk) */
                alpha1[ix+nx*iy]-=b*c;
                alpha2[ix+nx*iy]+=b*b;
 
              }
         }
       } 
}
//a#############################################################################################
__global__ void cuda_cal_alpha(float *alpha, float *alpha1, float *alpha2, float epsil, int ng)
/*< calculate searched stepsize (alpha) according to Taratola's method
configuration requirement: <<<1, Block_Size>>> >*/ 
{
  	__shared__ float sdata[BlockSize];
	__shared__ float tdata[BlockSize];
    	int tid=threadIdx.x;
    	sdata[tid]=0.0f;
	tdata[tid]=0.0f;
	for(int s=0; s<(ng+BlockSize-1)/BlockSize; s++)
	{
		int id=s*blockDim.x+threadIdx.x;
		float a=(id<ng)?alpha1[id]:0.0f;
		float b=(id<ng)?alpha2[id]:0.0f;
		sdata[tid] +=a;	
		tdata[tid] +=b;	
	} 
    	__syncthreads();

    	/* do reduction in shared mem */
    	for(int s=blockDim.x/2; s>32; s>>=1) 
    	{
		if (threadIdx.x < s) { sdata[tid] += sdata[tid + s];tdata[tid] += tdata[tid + s]; } __syncthreads();
    	}
   	if (tid < 32)
   	{
		if (blockDim.x >=  64) { sdata[tid] += sdata[tid + 32]; tdata[tid] += tdata[tid + 32];}
		if (blockDim.x >=  32) { sdata[tid] += sdata[tid + 16]; tdata[tid] += tdata[tid + 16];}
		if (blockDim.x >=  16) { sdata[tid] += sdata[tid +  8]; tdata[tid] += tdata[tid +  8];}
		if (blockDim.x >=   8) { sdata[tid] += sdata[tid +  4]; tdata[tid] += tdata[tid +  4];}
		if (blockDim.x >=   4) { sdata[tid] += sdata[tid +  2]; tdata[tid] += tdata[tid +  2];}
		if (blockDim.x >=   2) { sdata[tid] += sdata[tid +  1]; tdata[tid] += tdata[tid +  1];}
    	}
     
    	if (tid == 0) 
        {
             if(tdata[0]>EPS) *alpha=epsil*sdata[0]/(tdata[0]+EPS); 
             else *alpha=0.0;
        }
}
//a#############################################################################################
__global__ void cuda_update_vel(float *VV, float *cg, float alpha, int nnx, int nny, int nnz)
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

       for(iy=0;iy<nny;iy++)
        {
          id=iz+ix*nnz+iy*nnz*nnx;
          if(id<nnx*nny*nnz) VV[id]=VV[id]+alpha*cg[id];
        }  
}
//a#############################################################################################
//a###                                                                                       ###
//a###                                 Main Function                                         ###
//a###                                                                                       ###
//a#############################################################################################
int main(int argc, char* argv[])
{

	int nz, nx, ny, nnz, nnx, nny, ns, nsx, nt, it, is, fsz, fsx, fsy, dsz, dsx, dsy, ng, iter, niter;
	int *coo_source, *coo_receivers, *coo_up, *coo_down, *coo_front, *coo_back, *coo_left, *coo_right;
	float dz, dx, dy, favg, dt, _dz2, _dx2, _dy2, pfac;
	float *v, *vv, *wavelet, *VV, *VVtmp, *s_P0, *s_P1, *ptr, *g_P0, *g_P1, *s_Ptt;
       float *p_cal, *p_IO, *p_obs, *p_com, *s_bndr, *p_derr;
       float *g0, *g1, *cg, *illum, *pars;
       float obj1, obj, beta, epsil, alpha, *alpha1, *alpha2, *objval;
 //a######################################
	char FNvel[250]={"vel201202203initial.dat"};
	char FNsobs[250]={"shot_obs.dat"};
	char FNscal[250]={"shot_cal.dat"};
       char FNscom[250]={"shot_com.dat"};
	char FNgrad[250]={"gradient.dat"};
       char FNillum[250]={"illum.dat"};
	char FNupdatevel[250]={"velupdate.dat"};
	char FNlastvel[250]={"vellastIter.dat"};
	char FNobjs[250]={"objections.txt"};
//a######################################
    	nx=201;    dx=10;
    	ny=1;    dy=10;
    	nz=203;    dz=10;

   	nt=1501;     favg=20;   pfac=100;
       dt=0.001;  

	ns=5;    nsx=5; 

	fsx=10;   dsx=40;
	fsy=1;   dsy=0;
	fsz=1;     dsz=0;

       niter=50;
//a######################################
  FILE *fpvel, *fpscal, *fpsobs, *fpgrad, *fpscom, *fpillum, *fpupdatevel, *fplastvel, *fpobjs;
  if((fpvel=fopen(FNvel,"rb"))==NULL){printf("###   < %s > read error!\n",FNvel);exit(0);}
if((fpsobs=fopen(FNsobs,"rb"))==NULL){printf("###   < %s > read error!\n",FNsobs);exit(0);}
    fpscal=fopen(FNscal,"wb");
    fpscom=fopen(FNscom,"wb");
    fpgrad=fopen(FNgrad,"wb");
   fpillum=fopen(FNillum,"wb");
   fpupdatevel=fopen(FNupdatevel,"wb");
   fplastvel=fopen(FNlastvel,"wb");
  fpobjs=fopen(FNobjs,"w");
//a######################################
	_dz2=1.0/(dz*dz);
	_dx2=1.0/(dx*dx);
	_dy2=1.0/(dy*dy);

	nnz=nz+2*mm+2*npd;
	nnx=nx+2*mm+2*npd;
	nny=ny+2*mm+2*npd;

       ng=nx*ny;
//a######################################
    	v=(float*)malloc(nz*nx*ny*sizeof(float));
    	vv=(float*)malloc(nnz*nnx*nny*sizeof(float));
    	p_IO=(float*)malloc(ng*nt*sizeof(float));
    	objval=(float*)malloc(niter*sizeof(float));

	memset(p_IO, 0, ng*nt*sizeof(float));
	memset(objval, 0, niter*sizeof(float));

	fread(v, sizeof(float), nz*nx*ny, fpvel);
	velocity_transform(v, vv, dt, nz, nx, ny);

       /*< initialize device, default device=0 >*/
    	cudaSetDevice(0);
	check_gpu_error("Failed to initialize device!");

	dim3 dimg, dimb, dimt;
	dimg.x=(nz+2*npd+2*mm+BlockSize1-1)/BlockSize1;
	dimg.y=(nx+2*npd+2*mm+BlockSize2-1)/BlockSize2;
	dimt.x=(nt+BlockSize1-1)/BlockSize1;
	dimt.y=(nx+BlockSize2-1)/BlockSize2;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;

	/* allocate memory on device */
       /*< wavelet & velocity >*/
	cudaMalloc(&wavelet, nt*sizeof(float));
	cudaMalloc(&VV, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&VVtmp, nnz*nnx*nny*sizeof(float));
       /*< forward & backward & receivers wavefield >*/
	cudaMalloc(&s_P0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&s_P1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&g_P0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&g_P1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&s_Ptt, nnz*nnx*nny*sizeof(float));
       /*< shot & receivers location >*/
	cudaMalloc(&coo_source, ns*sizeof(int));
	cudaMalloc(&coo_receivers, ng*sizeof(int));
       /*< boundary location >*/
	cudaMalloc(&coo_up , nx*ny*sizeof(int));
	cudaMalloc(&coo_down , nx*ny*sizeof(int));
	cudaMalloc(&coo_front, nx*nz*sizeof(int));
	cudaMalloc(&coo_back , nx*nz*sizeof(int));
	cudaMalloc(&coo_left , ny*nz*sizeof(int));
	cudaMalloc(&coo_right, ny*nz*sizeof(int));
       /*< calculated/synthetic seismic data (it & nt & 6's boundary) >*/
	cudaMalloc(&p_cal, ng*nt*sizeof(float)); 
	cudaMalloc(&p_obs, ng*nt*sizeof(float)); 
	cudaMalloc(&p_com, ng*nt*sizeof(float)); 
	cudaMalloc(&p_derr, ng*nt*sizeof(float)); 
	cudaMalloc(&alpha1, ng*sizeof(float)); 
	cudaMalloc(&alpha2, ng*sizeof(float)); 
	cudaMalloc(&s_bndr, nt*(2*nz*nx+2*nz*ny+2*nx*ny)*sizeof(float));
       /*< The is & ns gradient ,lighting matrix >*/
	cudaMalloc(&g0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&g1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&cg, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&illum, nnz*nnx*nny*sizeof(float));
       cudaMemset(g1, 0, nnz*nnx*nny*sizeof(float));
       cudaMemset(cg, 0, nnz*nnx*nny*sizeof(float));
       /* d_pars[0]: obj; d_pars[1]: beta; d_pars[2]: epsilon; d_pars[3]: alpha; */
       cudaMalloc(&pars, 4*sizeof(float));	
	cudaMemset(pars, 0, 4*sizeof(float));

	check_gpu_error("Failed to allocate memory for variables!");

	cuda_ricker_wavelet<<<(nt+BlockSize-1)/BlockSize, BlockSize>>>(wavelet, favg, dt, nt, pfac);
	cudaMemcpy(VV, vv, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
       /*< shot location >*/
	cuda_set_s<<<1, ns>>>(coo_source, fsz, fsx, fsy, dsz, dsx, dsy, ns, nsx, nz, nx, ny);
       /*< receivers(up),down,front,back,left,right location >*/
	cuda_set_up_do<<<(nx*ny+BlockSize-1)/BlockSize,BlockSize>>>(coo_receivers, coo_up,coo_down, nx*ny, nz, nx, ny); 
	cuda_set_fr_ba<<<(nz*nx+BlockSize-1)/BlockSize,BlockSize>>>(coo_front, coo_back,  nz*nx, nz, nx, ny);
	cuda_set_le_ri<<<(nz*ny+BlockSize-1)/BlockSize,BlockSize>>>(coo_left,  coo_right, nz*ny, nz, nx, ny);

	clock_t iter_t0, iter_t1, is_t0, is_t1, ns_t0, ns_t1;

       printf("##########################################\n");
       printf("###\n");
       for(iter=0; iter<niter; iter++)
        {
              iter_t0=clock();
              printf("########## Iter =%3d  ##########  \n###\n",iter+1);
		cudaMemcpy(g0, g1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToDevice);
	       cudaMemset(g1, 0, nnz*nnx*nny*sizeof(float));
	       cudaMemset(illum, 0, nnz*nnx*nny*sizeof(float));
	       cudaMemset(p_derr, 0, ng*nt*sizeof(float));
              cudaMemset(alpha1, 0, ng*sizeof(float));
              cudaMemset(alpha2, 0, ng*sizeof(float));
              cudaMemset(pars, 0, 4*sizeof(float));

              rewind(fpscal);
              rewind(fpsobs);
              rewind(fpscom);
              rewind(fpillum);
              rewind(fpgrad);

              ns_t0=clock();
	       for(is=0; is<ns; is++)
                {   
                   is_t0=clock();
	            cudaMemset(s_P0, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(s_P1, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(g_P0, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(g_P1, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(s_Ptt, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(p_cal, 0, ng*nt*sizeof(float));
	            cudaMemset(p_obs, 0, ng*nt*sizeof(float));
	            cudaMemset(p_com, 0, ng*nt*sizeof(float));
	            cudaMemset(s_bndr, 0, nt*(2*nz*nx+2*nz*ny+2*nx*ny)*sizeof(float));

	            for(it=0; it<nt; it++)
                      {
                          //if(it%400==0) printf("For: is=%2d, it=%d\n",is,it);
	                   cuda_add_source<<<1,1>>>(true, s_P1, &wavelet[it], &coo_source[is], 1);
	                   cuda_step_fd3d<<<dimg,dimb>>>(s_P0, s_P1, VV, _dz2, _dx2, _dy2, nz, nx, ny, dt, NULL, false);
	                   ptr=s_P0; s_P0=s_P1; s_P1=ptr;
                          cuda_absorb_bndr<<<dimg,dimb>>>(s_P0, s_P1, nz, nx, ny, -0.25);

                          cuda_save_bndr<<<((2*nz*nx+2*nz*ny+2*nx*ny)+BlockSize-1)/BlockSize,BlockSize>>>(
                                                           &s_bndr[it*(2*nz*nx+2*nz*ny+2*nx*ny)], 
                                                           s_P0, coo_front, coo_back, coo_left, coo_right, coo_up, coo_down,
                                                           nz, nx, ny, true);
                          cuda_cal_illum<<<dimg,dimb>>>(illum, s_P0, nz, nx, ny);

	                   cuda_record<<<(ng+BlockSize-1)/BlockSize, BlockSize>>>(s_P0, p_cal, coo_receivers, ng, it, nt, true);

	             }//it loop end

	              cudaMemcpy(p_IO, p_cal, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
	                        fwrite(p_IO, sizeof(float), ng*nt, fpscal);

                    fseek(fpsobs,is*ng*nt*sizeof(float),0);
	             fread(p_IO, sizeof(float), ng*nt, fpsobs);
	             cudaMemcpy(p_obs, p_IO, ng*nt*sizeof(float), cudaMemcpyHostToDevice);

                    cuda_cal_residuals<<<dimt, dimb>>>(&pars[0], p_cal, p_obs, p_com, ng*nt, nx, ny, nt);
                    if(is==0)cuda_com2derr<<<dimt, dimb>>>(p_com, p_derr, nx, ny, nt);
		      cudaMemcpy(&obj, &pars[0], sizeof(float), cudaMemcpyDeviceToHost);

	            if(is==(ns/2+1)){   cudaMemcpy(p_IO, p_com, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
                                fseek(fpscom,is*ng*nt*sizeof(float),0);
	                         fwrite(p_IO, sizeof(float), ng*nt, fpscom); }

	             for(it=nt-1; it>-1; it--)
                       {
                          //if(it%400==0) printf("Back: is=%2d, it=%d\n",is,it);
	                   ptr=s_P0; s_P0=s_P1; s_P1=ptr;
	                   cuda_save_bndr<<<((2*nz*nx+2*nz*ny+2*nx*ny)+BlockSize-1)/BlockSize,BlockSize>>>(
                                                           &s_bndr[it*(2*nz*nx+2*nz*ny+2*nx*ny)], 
                                                           s_P1, coo_front, coo_back, coo_left, coo_right, coo_up, coo_down, 
                                                           nz, nx, ny, false);
	                   cuda_step_fd3d<<<dimg,dimb>>>(s_P0, s_P1, VV, _dz2, _dx2, _dy2, nz, nx, ny, dt, s_Ptt, true);
                          cuda_absorb_bndr<<<dimg,dimb>>>(s_P0, s_P1, nz, nx, ny, -0.25);


                          cuda_record<<<(ng+BlockSize-1)/BlockSize, BlockSize>>>(g_P1, p_com, coo_receivers, ng, it, nt, false);
	                   cuda_step_fd3d<<<dimg,dimb>>>(g_P0, g_P1, VV, _dz2, _dx2, _dy2, nz, nx, ny, dt, NULL, false);
	                   ptr=g_P0; g_P0=g_P1; g_P1=ptr;
                          cuda_absorb_bndr<<<dimg,dimb>>>(g_P0, g_P1, nz, nx, ny, -0.25);

                          cuda_cal_illum<<<dimg,dimb>>>(illum, g_P1, nz, nx, ny);
                          cuda_cal_g1<<<dimg,dimb>>>(g1, s_Ptt, g_P1, nz, nx, ny);
	              }// it loop end
           
                     is_t1=clock();
                     printf("###   IS:(%2d) %.2f(min);\n",is,((float)(is_t1-is_t0))/60000000.0);
	       }//IS loop end
              ns_t1=clock();
              printf("###   Cal gradient: %.2f (min)\n",((float)(ns_t1-ns_t0))/60000000.0);

              cudaMemcpy(vv, illum, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
              window3d(v, vv, nz, nx, ny);
              fwrite(v, sizeof(float), nz*nx*ny, fpillum);
		/* compute the gradient of FWI by scaling, precondition incorporated here */
              cuda_scale_gradient<<<dimg,dimb>>>(g1, VV, illum, nnx, nny, nnz, true);
		/* Gaussian smoothing for the sharp gradient */
              cuda_bell_smoothz<<<dimg,dimb>>>(g1, 2, nnx, nny, nnz);
              cuda_bell_smoothx<<<dimg,dimb>>>(g1, 2, nnx, nny, nnz);
              cuda_bell_smoothy<<<dimg,dimb>>>(g1, 2, nnx, nny, nnz);

		/* calculate the factor beta in conjugate gradient method */
              if (iter>0) 
                {
                    cuda_cal_beta_step1<<<dimg,dimb>>>(g0, g1, cg, nnx, nny, nnz);
                    cuda_cal_beta_step2<<<1,1>>>(&pars[1], nnx, nny, nnz);
                }
		cudaMemcpy(&beta, &pars[1], sizeof(float), cudaMemcpyDeviceToHost);
		/* compute the conjugate gradient */
              cuda_cal_conjgrad<<<dimg,dimb>>>(g1, cg, beta, nnx, nny, nnz);

              cudaMemcpy(vv, cg, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
              window3d(v, vv, nz, nx, ny);
              fwrite(v, sizeof(float), nz*nx*ny, fpgrad);

		/* estimate epsilon according to equation 11 */
		cuda_cal_epsilon<<<1, BlockSize>>>(VV, cg, &pars[2], nnx*nnz*nny);
		cudaMemcpy(&epsil, &pars[2], sizeof(float), cudaMemcpyDeviceToHost);

		/* obtain a tentative velocity model to estimate a good stepsize alpha */
              cuda_cal_vtmp<<<dimg,dimb>>>(VVtmp, VV, cg, epsil, nnx, nny, nnz);

             ns_t0=clock();
             printf("###   Cal alpha:");
	      for(is=0; is<1; is++)
               {   
	            cudaMemset(s_P0, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(s_P1, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(p_cal, 0, ng*nt*sizeof(float));

                   fseek(fpsobs,is*ng*nt*sizeof(float),0);
	            fread(p_IO, sizeof(float), ng*nt, fpsobs);
	            cudaMemcpy(p_obs, p_IO, ng*nt*sizeof(float), cudaMemcpyHostToDevice);

	            for(it=0; it<nt; it++)
                      {
	                   cuda_add_source<<<1,1>>>(true, s_P1, &wavelet[it], &coo_source[is], 1);
	                   cuda_step_fd3d<<<dimg,dimb>>>(s_P0, s_P1, VVtmp, _dz2, _dx2, _dy2, nz, nx, ny, dt, NULL, false);
	                   ptr=s_P0; s_P0=s_P1; s_P1=ptr;
                          cuda_absorb_bndr<<<dimg,dimb>>>(s_P0, s_P1, nz, nx, ny, -0.25);
	                   cuda_record<<<(ng+BlockSize-1)/BlockSize, BlockSize>>>(s_P0, p_cal, coo_receivers, ng, it, nt, true);
	             }//it loop end
                    cuda_sum_alpha12<<<dimt, dimb>>>(alpha1, alpha2, p_cal, p_obs, p_derr, nx, ny, nz, nt);
              }//is loop end
     
		cuda_cal_alpha<<<1,BlockSize>>>(&pars[3], alpha1, alpha2, epsil, ng);
		cudaMemcpy(&alpha, &pars[3], sizeof(float), cudaMemcpyDeviceToHost);

              ns_t1=clock();printf(" %.2f (min)\n",((float)(ns_t1-ns_t0))/60000000.0);
		/* update the velocity model according to previous velocity, conjugate gradient and estimated stepsize */
              cuda_update_vel<<<dimg,dimb>>>(VV, cg, alpha, nnx, nny, nnz);
	       cudaMemcpy(vv, VV, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	       window3d(v, vv, nz, nx, ny);
	       fwrite(v, sizeof(float),nz*nx*ny, fpupdatevel);

		/* compute the normalized objective function */
		if(iter==0) 	{obj1=obj; objval[iter]=1.0;}
		else		objval[iter]=obj/obj1;

              iter_t1=clock();
              printf("###   objval=%f, beta=%f, epsil=%.2f, alpha=%.2f : %.2f(min)\n",
                                            objval[iter],beta,epsil,alpha,((float)(iter_t1-iter_t0))/60000000.0);
              fprintf(fpobjs,"iter=%3d, obj=%f;\n",iter+1,objval[iter]);
 
	       cudaMemcpy(vv, VV, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	       window3d(v, vv, nz, nx, ny);
              rewind(fplastvel);
	       fwrite(v, sizeof(float),nz*nx*ny, fplastvel);

        }//iter loop end
        printf("##################################\n");

	/* free memory on device */
	cudaFree(wavelet);
	cudaFree(VV);
	cudaFree(VVtmp);
      /*< wavefield(x-y-z) >*/
	cudaFree(s_P0);
	cudaFree(s_P1);
	cudaFree(g_P0);
	cudaFree(g_P1);
	cudaFree(s_Ptt);
      /*< location >*/
	cudaFree(coo_source);
	cudaFree(coo_receivers);
	cudaFree(coo_front);
	cudaFree(coo_back);
	cudaFree(coo_left);
	cudaFree(coo_right);
	cudaFree(coo_down);
	cudaFree(coo_up);
      /*< gradient >*/
	cudaFree(g0);
	cudaFree(g1);
	cudaFree(cg);
	cudaFree(illum);
	cudaFree(pars);
      /*< wavefield(t-x-y-z) >*/
	cudaFree(p_cal);
	cudaFree(p_obs);
	cudaFree(p_com);
	cudaFree(p_derr);
	cudaFree(alpha1);
	cudaFree(alpha2);
	cudaFree(s_bndr);
      /*< free alloc >*/
	free(v);
	free(vv);
	free(p_IO);
	free(objval);

       fclose(fpvel);
       fclose(fpscal);
       fclose(fpsobs);
       fclose(fpscom);
       fclose(fpgrad);
       fclose(fpillum);
       fclose(fpupdatevel);
       fclose(fplastvel);
       fclose(fpobjs);

    	exit (0);
}


