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
#define npd     50   // absorbing boundry condition wield


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
__global__ void cuda_illum(float *migration, float *illum, int nz, int nx, int ny)
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
          if(id<nnz*nnx*nny&&illum[id]!=0) migration[id]/=illum[id];
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
__global__ void cuda_cal_migration(float *migration, float *s, float *g, int nz, int nx, int ny)
/*< calculate is migration >*/
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
          if(id<nnz*nnx*nny) migration[id]+=s[id]*g[id];
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
__global__ void cuda_scale_gradient(float *migration, float *VV, float *illum, int nnx, int nny, int nnz, bool precon)
/*< scale migration >*/
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
		migration[id]*=2.0/a;
           }

        }  

}
//a#############################################################################################
__global__ void cuda_bell_smoothz(float *migration, int rbell, int nnx, int nny, int nnz)
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
              for(i=-rbell; i<=rbell; i++) if(iz+i>=0 && iz+i<nnz) s+=expf(-(2.0*i*i)/rbell)*migration[id+i];
              migration[id]=s;
           }

        }  
}
//a#############################################################################################
__global__ void cuda_bell_smoothx(float *migration, int rbell, int nnx, int nny, int nnz)
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
              for(i=-rbell; i<=rbell; i++) if(ix+i>=0 && ix+i<nnx) s+=expf(-(2.0*i*i)/rbell)*migration[id+i*nnz];
              migration[id]=s;
           }

        }  
}
//a#############################################################################################
__global__ void cuda_bell_smoothy(float *migration, int rbell, int nnx, int nny, int nnz)
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
              for(i=-rbell; i<=rbell; i++) if(iy+i>=0 && iy+i<nny) s+=expf(-(2.0*i*i)/rbell)*migration[id+i*nnz*nnx];
              migration[id]=s;
           }

        }  
}



/*************func**************/    
__global__ void mute_directwave(int nx,int ny,int nt,float dt,float favg, float dx,float dy,float dz,int fsx,int fsy,int dsx,int dsy,
                                int zs,int is, float *vp,float *shot,int tt,int nsx)
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
              mu_x=dx*abs(ix-fsx-(is%nsx)*dsx);
              mu_y=dy*abs(iy-fsy-(is/nsx)*dsy);
              mu_z=dz*zs;
              mu_t0=sqrtf(pow(mu_x,2)+pow(mu_y,2)+pow(mu_z,2))/(vp[1]);
              mu_t=(int)(2.0/(dt*favg));
              mu_nt=(int)(mu_t0/dt)+mu_t+tt;

                 if(it<mu_nt)
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

//a#############################################################################################
//a###                                                                                       ###
//a###                                 Main Function                                         ###
//a###                                                                                       ###
//a#############################################################################################
int main(int argc, char* argv[])
{

	int nz, nx, ny, nnz, nnx, nny, ns, nsx, nt, it, is, fsz, fsx, fsy, dsz, dsx, dsy, ng;
	int *coo_source, *coo_receivers, *coo_up, *coo_down, *coo_front, *coo_back, *coo_left, *coo_right;
	float dz, dx, dy, favg, dt, _dz2, _dx2, _dy2, pfac;
	float *v, *vv, *wavelet, *VV,  *s_P0, *s_P1, *ptr, *g_P0, *g_P1;
       float *p_cal, *p_IO,  *s_bndr;
       float *migration, *illum;

 //a######################################
	char FNvel[250]={"vel201202203.dat"};

	char FNscal[250]={"shot_cal.dat"};
	char FNgrad[250]={"migration.dat"};
       char FNillum[250]={"illumination.dat"};


//a######################################
    	nx=201;    dx=10;
    	ny=201;    dy=10;
    	nz=203;    dz=10;

   	nt=1501;     favg=20;   pfac=100;
       dt=0.001;  

	ns=1;    nsx=1; 

	fsx=100;   dsx=100;
	fsy=50;   dsy=100;
	fsz=1;     dsz=0;

//a######################################
  FILE *fpvel, *fpscal,  *fpmig,  *fpillum;
  if((fpvel=fopen(FNvel,"rb"))==NULL){printf("###   < %s > read error!\n",FNvel);exit(0);}
    fpscal=fopen(FNscal,"wb");
    fpmig=fopen(FNgrad,"wb");
   fpillum=fopen(FNillum,"wb");


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


	memset(p_IO, 0, ng*nt*sizeof(float));


	fread(v, sizeof(float), nz*nx*ny, fpvel);
	velocity_transform(v, vv, dt, nz, nx, ny);

       /*< initialize device, default device=0 >*/
    	cudaSetDevice(0);
	check_gpu_error("Failed to initialize device!");

	dim3 dimg, dimb,Xdimg, dimt;
	Xdimg.x=(nnx+BlockSize1-1)/BlockSize1;
	Xdimg.y=(nny+BlockSize2-1)/BlockSize2;

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
       /*< forward & backward & receivers wavefield >*/
	cudaMalloc(&s_P0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&s_P1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&g_P0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&g_P1, nnz*nnx*nny*sizeof(float));
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

	cudaMalloc(&s_bndr, nt*(2*nz*nx+2*nz*ny+2*nx*ny)*sizeof(float));
       /*< The is & ns gradient ,lighting matrix >*/

	cudaMalloc(&migration, nnz*nnx*nny*sizeof(float));

	cudaMalloc(&illum, nnz*nnx*nny*sizeof(float));
       cudaMemset(migration, 0, nnz*nnx*nny*sizeof(float));




	check_gpu_error("Failed to allocate memory for variables!");

	cuda_ricker_wavelet<<<(nt+BlockSize-1)/BlockSize, BlockSize>>>(wavelet, favg, dt, nt, pfac);
	cudaMemcpy(VV, vv, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
       /*< shot location >*/
	cuda_set_s<<<1, ns>>>(coo_source, fsz, fsx, fsy, dsz, dsx, dsy, ns, nsx, nz, nx, ny);
       /*< receivers(up),down,front,back,left,right location >*/
	cuda_set_up_do<<<(nx*ny+BlockSize-1)/BlockSize,BlockSize>>>(coo_receivers, coo_up,coo_down, nx*ny, nz, nx, ny); 
	cuda_set_fr_ba<<<(nz*nx+BlockSize-1)/BlockSize,BlockSize>>>(coo_front, coo_back,  nz*nx, nz, nx, ny);
	cuda_set_le_ri<<<(nz*ny+BlockSize-1)/BlockSize,BlockSize>>>(coo_left,  coo_right, nz*ny, nz, nx, ny);

	clock_t  is_t0, is_t1, ns_t0, ns_t1;

       printf("##########################################\n");
       printf("###\n");



	       cudaMemset(migration, 0, nnz*nnx*nny*sizeof(float));
	       cudaMemset(illum, 0, nnz*nnx*nny*sizeof(float));



              rewind(fpscal);


              rewind(fpillum);
              rewind(fpmig);

              ns_t0=clock();
	       for(is=0; is<ns; is++)
                {   printf("###  is= %d\n",is);
                   is_t0=clock();
	            cudaMemset(s_P0, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(s_P1, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(g_P0, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(g_P1, 0, nnz*nnx*nny*sizeof(float));
	            cudaMemset(p_cal, 0, ng*nt*sizeof(float));

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


                     mute_directwave<<<Xdimg,dimb>>>(nx,ny,nt,dt,favg,dx,dy,dz,fsx,fsy,dsx,dsy,fsz,is,VV,p_cal,70,nsx);

	              cudaMemcpy(p_IO, p_cal, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
	                        fwrite(p_IO, sizeof(float), ng*nt, fpscal);



	             for(it=nt-1; it>-1; it--)
                       {
                          //if(it%400==0) printf("Back: is=%2d, it=%d\n",is,it);
	                   ptr=s_P0; s_P0=s_P1; s_P1=ptr;
	                   cuda_save_bndr<<<((2*nz*nx+2*nz*ny+2*nx*ny)+BlockSize-1)/BlockSize,BlockSize>>>(
                                                           &s_bndr[it*(2*nz*nx+2*nz*ny+2*nx*ny)], 
                                                           s_P1, coo_front, coo_back, coo_left, coo_right, coo_up, coo_down, 
                                                           nz, nx, ny, false);
	                   cuda_step_fd3d<<<dimg,dimb>>>(s_P0, s_P1, VV, _dz2, _dx2, _dy2, nz, nx, ny, dt, NULL, true);
                          cuda_absorb_bndr<<<dimg,dimb>>>(s_P0, s_P1, nz, nx, ny, -0.25);


                          cuda_record<<<(ng+BlockSize-1)/BlockSize, BlockSize>>>(g_P1, p_cal, coo_receivers, ng, it, nt, false);
	                   cuda_step_fd3d<<<dimg,dimb>>>(g_P0, g_P1, VV, _dz2, _dx2, _dy2, nz, nx, ny, dt, NULL, false);
	                   ptr=g_P0; g_P0=g_P1; g_P1=ptr;
                          cuda_absorb_bndr<<<dimg,dimb>>>(g_P0, g_P1, nz, nx, ny, -0.25);

                          cuda_cal_illum<<<dimg,dimb>>>(illum, g_P1, nz, nx, ny);
                          cuda_cal_migration<<<dimg,dimb>>>(migration, s_P1, g_P1, nz, nx, ny);
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
              cuda_scale_gradient<<<dimg,dimb>>>(migration, VV, illum, nnx, nny, nnz, true);

              cudaMemcpy(vv, migration, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
              window3d(v, vv, nz, nx, ny);
             // laplace_3d_filter(1, nz, nx, ny, v, vv);
              fwrite(v, sizeof(float), nz*nx*ny, fpmig);

        printf("##################################\n");

	/* free memory on device */
	cudaFree(wavelet);
	cudaFree(VV);
      /*< wavefield(x-y-z) >*/
	cudaFree(s_P0);
	cudaFree(s_P1);
	cudaFree(g_P0);
	cudaFree(g_P1);
      /*< location >*/
	cudaFree(coo_source);
	cudaFree(coo_receivers);
	cudaFree(coo_front);
	cudaFree(coo_back);
	cudaFree(coo_left);
	cudaFree(coo_right);
	cudaFree(coo_down);
	cudaFree(coo_up);

	cudaFree(migration);

	cudaFree(illum);

      /*< wavefield(t-x-y-z) >*/
	cudaFree(p_cal);


	cudaFree(s_bndr);
      /*< free alloc >*/
	free(v);
	free(vv);
	free(p_IO);


       fclose(fpvel);
       fclose(fpscal);


       fclose(fpmig);
       fclose(fpillum);


    	exit (0);
}


