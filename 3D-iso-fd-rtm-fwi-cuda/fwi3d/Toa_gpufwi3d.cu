//
//
/*                   3D acousic FWI           */
//
//
/*  Initial modeling Code (fd3d) in 2014, 
   Xi'an Jiaotong University (by Pengliang Yang) */
//
//
/*  Ps :
             absorb boundary(mm-npd-nxyz-npd-mm), 




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
//a#############################################################################################
__constant__ float stencil[mm+1]={-205.0/72.0,8.0/5.0,-1.0/5.0,8.0/315.0,-1.0/560.0};
//a#############################################################################################
void sf_check_gpu_error (const char *msg);
void velocity_transform(float *v0, float*vv, float dt, int n1, int n2, int n3);
void window3d(float *a, float *b, int n1, int n2, int n3);
//a#############################################################################################
__global__ void cuda_ricker_wavelet(float *wlt, float fm, float dt, int nt);
__global__ void cuda_set_s(int *szxy, int szbeg, int sxbeg, int sybeg, int jsz, int jsx, int jsy, int ns, int nz, int nx, int ny);
__global__ void cuda_set_up_do(int *gzxy, int *down, int ng, int nz, int nx, int ny);
__global__ void cuda_set_fr_ba(int *front, int *back, int ng, int nz, int nx, int ny);
__global__ void cuda_set_le_ri(int *left, int *right,int ng, int nz,int nx, int ny);
__global__ void cuda_save_bndr(float *bndr,float *p0, int *front, int *back, int *left, int *right, int *up, int *down, 
                               int nz, int nx, int ny, bool write);
__global__ void cuda_trans_xy2txy(float *xy, float *txy, int it, int nt, int ng, bool output);
__global__ void cuda_absorb_bndr(float *d_p,int nz,int nx,int ny,float qp);
__global__ void cuda_free_surface(float *p0,float *p1,int nz,int nx,int ny);
__global__ void cuda_record(float *p, float *seis, int *gxz, int ng, bool record);
/*__global__ void cuda_load_receiver(bool mute,float *bndr, float *p0, int *up, int nz, int nx, int ny, 
                                   int it, float fm, float dt, float dz,float dx, float dy, int is, int *szxy, float v_mu);*/
__global__ void cuda_add_source(bool add, float *p, float *source, int *szxy, int ns);
__global__ void cuda_step_fd3d(float *p0, float *p1, float *vv, float _dz2, float _dx2, float _dy2, int n1, int n2, int n3, 
                               float dt, float *d_pdt2, bool pdt);
__global__ void cuda_cal_gradient(float *g_is, float *s, float *g, int nz, int nx, int ny);
__global__ void cuda_light_matrix(float *s, float *p, int nz, int nx, int ny);
__global__ void cuda_IS_lighting(float *g_is, float *source, int nz, int nx, int ny);
__global__ void cuda_sum(float *ns, float *is, int nz, int nx, int ny);
__global__ void cuda_cal_residuals(float *cal, float *obs, float *com, int nn, int nx, int ny, int nt, float obj);

//a#############################################################################################
int main(int argc, char* argv[])
{
	bool verb;
	int nz, nx, ny, nnz, nnx, nny, ns, nt, kt, it, is, szbeg, sxbeg, sybeg, jsz, jsx, jsy, ng;
	int *d_szxy, *d_gzxy, *bndr_down, *bndr_front, *bndr_back, *bndr_left, *bndr_right;
	float dz, dx, dy, fm, dt, _dz2, _dx2, _dy2, obj;
	float *v0, *vv, *d_wlt, *d_vv, *d_p0, *d_p1, *ptr, *d_g0, *d_g1, *d_pdt2; //, v_mu;
       float *p_cal_it, *p_cal_nt, *p_IO, *p_obs_it, *p_obs_nt, *p_com, *d_bndr;
       float *g_is, *g_ns, *source;


	sf_file Fv, Fsnap, Fscal, Fg, Fsobs, Fscom, Fobj;

    	sf_init(argc,argv);
	Fv=sf_input("vel");
	Fsobs=sf_input("shot_obs");
	Fsnap=sf_output("snap");
	Fscal=sf_output("shot_cal");
	Fscom=sf_output("shot_com");
	Fg=sf_output("gradient");
	Fobj=sf_output("objection");

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

	sf_putint(Fsnap,"n1",nz);
	sf_putint(Fsnap,"n2",nx);
	sf_putint(Fsnap,"n3",ny);
	sf_putint(Fg,"n1",nz);
	sf_putint(Fg,"n2",nx);
	sf_putint(Fg,"n3",ny);



	_dz2=1.0/(dz*dz);
	_dx2=1.0/(dx*dx);
	_dy2=1.0/(dy*dy);

	nnz=nz+2*mm+2*npd;
	nnx=nx+2*mm+2*npd;
	nny=ny+2*mm+2*npd;

       ng=nx*ny;


	sf_putint(Fscal,"n1",nt);
	sf_putint(Fscal,"n2",nx);
	sf_putint(Fscal,"n3",ny);
	sf_putint(Fscom,"n1",nt);
	sf_putint(Fscom,"n2",nx);
	sf_putint(Fscom,"n3",ny);

    	v0=(float*)malloc(nz*nx*ny*sizeof(float));
    	vv=(float*)malloc(nnz*nnx*nny*sizeof(float));
    	p_IO=(float*)malloc(ng*nt*sizeof(float));

	sf_floatread(v0, nz*nx*ny, Fv);


      // v_mu=v0[1];
	velocity_transform(v0, vv, dt, nz, nx, ny);

       /*< initialize device, default device=0 >*/
    	cudaSetDevice(0);
	sf_check_gpu_error("Failed to initialize device!");

	dim3 dimg, dimb, dimt;
	dimg.x=(nz+2*npd+2*mm+BlockSize1-1)/BlockSize1;
	dimg.y=(nx+2*npd+2*mm+BlockSize2-1)/BlockSize2;
	dimt.x=(nt+BlockSize1-1)/BlockSize1;
	dimt.y=(nx+BlockSize2-1)/BlockSize2;
	dimb.x=BlockSize1;
	dimb.y=BlockSize2;

	/* allocate memory on device */
       /*< wavelet & velocity >*/
	cudaMalloc(&d_wlt, nt*sizeof(float));
	cudaMalloc(&d_vv, nnz*nnx*nny*sizeof(float));
       /*< forward & backward & receivers wavefield >*/
	cudaMalloc(&d_p0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_p1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_g0, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_g1, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&d_pdt2, nnz*nnx*nny*sizeof(float));
       /*< shot & receivers location >*/
	cudaMalloc(&d_szxy, ns*sizeof(int));
	cudaMalloc(&d_gzxy, ng*sizeof(int));
       /*< boundary location >*/
	cudaMalloc(&bndr_down , nx*ny*sizeof(int));
	cudaMalloc(&bndr_front, nx*nz*sizeof(int));
	cudaMalloc(&bndr_back , nx*nz*sizeof(int));
	cudaMalloc(&bndr_left , ny*nz*sizeof(int));
	cudaMalloc(&bndr_right, ny*nz*sizeof(int));
       /*< calculated/synthetic seismic data (it & nt & 6's boundary) >*/
	cudaMalloc(&p_cal_it, ng*sizeof(float));
	cudaMalloc(&p_obs_it, ng*sizeof(float));
	cudaMalloc(&p_cal_nt, ng*nt*sizeof(float)); 
	cudaMalloc(&p_obs_nt, ng*nt*sizeof(float)); 
	cudaMalloc(&p_com, ng*nt*sizeof(float)); 
	cudaMalloc(&d_bndr, nt*(2*nz*nx+2*nz*ny+2*nx*ny)*sizeof(float));
       /*< The is & ns gradient ,lighting matrix >*/
	cudaMalloc(&g_is, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&g_ns, nnz*nnx*nny*sizeof(float));
	cudaMalloc(&source, nnz*nnx*nny*sizeof(float));


	sf_check_gpu_error("Failed to allocate memory for variables!");

	cuda_ricker_wavelet<<<(nt+511)/512, 512>>>(d_wlt, fm, dt, nt);
	cudaMemcpy(d_vv, vv, nnz*nnx*nny*sizeof(float), cudaMemcpyHostToDevice);
       /*< shot location >*/
	cuda_set_s<<<1, ns>>>(d_szxy, szbeg, sxbeg, sybeg, jsz, jsx, jsy, ns, nz, nx, ny);
       /*< receivers(up),down,front,back,left,right location >*/
	cuda_set_up_do<<<(nx*ny+511)/512,512>>>(d_gzxy,     bndr_down,  nx*ny, nz, nx, ny); 
	cuda_set_fr_ba<<<(nz*nx+511)/512,512>>>(bndr_front, bndr_back,  nz*nx, nz, nx, ny);
	cuda_set_le_ri<<<(nz*ny+511)/512,512>>>(bndr_left,  bndr_right, nz*ny, nz, nx, ny);

	float mstimer;
	clock_t t0, t1;
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);

       obj=0.0;

       /*< Multi gradient >*/
	cudaMemset(g_ns, 0, nnz*nnx*nny*sizeof(float));

	sf_seek(Fscal, 0L, SEEK_SET);
	sf_seek(Fscom, 0L, SEEK_SET);
       /*< IS shot loop start ! >*/
	for(is=0; is<ns; is++){
	  cudaEventRecord(start);
         /*< matrix set '0' >*/
	  cudaMemset(d_p0, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(d_p1, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(p_cal_it, 0, ng*sizeof(float));
	  cudaMemset(p_cal_nt, 0, ng*nt*sizeof(float));
	  cudaMemset(p_obs_it, 0, ng*sizeof(float));
	  cudaMemset(p_obs_nt, 0, ng*nt*sizeof(float));
	  cudaMemset(p_com,    0, ng*nt*sizeof(float));
	  cudaMemset(d_bndr, 0, nt*(2*nz*nx+2*nz*ny+2*nx*ny)*sizeof(float));

	  cudaMemset(g_is, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(source, 0, nnz*nnx*nny*sizeof(float));


         /*<##################################################a>*/
         /*<###               Forward  Start               ###>*/
         /*<##################################################a>*/
	  for(it=0; it<nt; it++){
           if(it%100==0) sf_warning("forw: is=%2d, it=%d",is,it);

            /*< FD3D >*/
	    cuda_add_source<<<1,1>>>(true, d_p1, &d_wlt[it], &d_szxy[is], 1);
	    cuda_step_fd3d<<<dimg,dimb>>>(d_p0, d_p1, d_vv, _dz2, _dx2, _dy2, nz, nx, ny, dt, d_pdt2, false);
	    ptr=d_p0; d_p0=d_p1; d_p1=ptr;
           cuda_absorb_bndr<<<dimg,dimb>>>(d_p0, nz, nx, ny, -0.25);
           cuda_absorb_bndr<<<dimg,dimb>>>(d_p1, nz, nx, ny, -0.25);
           cuda_save_bndr<<<((2*nz*nx+2*nz*ny+2*nx*ny)+511)/512,512>>>(&d_bndr[it*(2*nz*nx+2*nz*ny+2*nx*ny)], d_p0, 
                                                bndr_front, bndr_back, bndr_left, bndr_right, d_gzxy, bndr_down, nz, nx, ny, true);

           cuda_light_matrix<<<dimg,dimb>>>(source, d_p0, nz, nx, ny);

           //cuda_free_surface<<<dimg,dimb>>>(d_p0, d_p1, nz, nx, ny);


            /*< record & (xy)2(txy) >*/
	    cuda_record<<<(ng+511)/512, 512>>>(d_p0, p_cal_it, d_gzxy, ng, true);
           cuda_trans_xy2txy<<<(ng+511)/512, 512>>>(p_cal_it, p_cal_nt, it, nt, ng, true);


	/*    if(is==0&&it%300==0){
	      t0 = clock();

	      cudaMemcpy(vv, d_p0, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	      window3d(v0, vv, nz, nx, ny);
	      sf_floatwrite(v0, nz*nx*ny, Fsnap);	  

 	      t1 = clock();
 	      sf_warning("save the volume: %f (s)", ((float)(t1-t0))/CLOCKS_PER_SEC); 	
	     } */

	   }//it loop end

         /*< write common shot gathers >*/
	  cudaMemcpy(p_IO, p_cal_nt, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
	  sf_floatwrite(p_IO, ng*nt, Fscal);


	  cudaMemset(d_g0, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(d_g1, 0, nnz*nnx*nny*sizeof(float));
	  cudaMemset(d_pdt2, 0, nnz*nnx*nny*sizeof(float));

	  sf_floatread(p_IO, ng*nt, Fsobs);
	  cudaMemcpy(p_obs_nt, p_IO, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
         cuda_cal_residuals<<<dimt, dimb>>>(p_cal_nt, p_obs_nt, p_com, ng*nt, nx, ny, nt, obj);


	  cudaMemcpy(p_IO, p_com, ng*nt*sizeof(float), cudaMemcpyDeviceToHost);
	  sf_floatwrite(p_IO, ng*nt, Fscom);
         /*<##################################################a>*/
         /*<###              Backward  Start               ###>*/
         /*<##################################################a>*/
	  for(it=nt-1; it>-1; it--){
           if(it%100==0) sf_warning("back: is=%2d, it=%d",is,it);
            /*< FD3D >*/
	    ptr=d_p0; d_p0=d_p1; d_p1=ptr;
	    cuda_save_bndr<<<((2*nz*nx+2*nz*ny+2*nx*ny)+511)/512,512>>>(&d_bndr[it*(2*nz*nx+2*nz*ny+2*nx*ny)], d_p1, 
                                                bndr_front, bndr_back, bndr_left, bndr_right, d_gzxy, bndr_down, nz, nx, ny, false);
	    cuda_step_fd3d<<<dimg,dimb>>>(d_p0, d_p1, d_vv, _dz2, _dx2, _dy2, nz, nx, ny, dt, d_pdt2, true);
           cuda_absorb_bndr<<<dimg,dimb>>>(d_p0, nz, nx, ny, -0.25);
           cuda_absorb_bndr<<<dimg,dimb>>>(d_p1, nz, nx, ny, -0.25);

           //cuda_free_surface<<<dimg,dimb>>>(d_p0, d_p1, nz, nx, ny);

            /*< FD3D >*/
	    //cuda_load_receiver<<<(nx*ny+511)/512,512>>>(true, &d_bndr[it*(2*nz*nx+2*nz*ny+2*nx*ny)], d_g1, d_gzxy, nz, nx, ny, 
           //                                             it, fm, dt, dz, dx, dy, is, d_szxy, v_mu);
           
           cuda_trans_xy2txy<<<(ng+511)/512, 512>>>(p_obs_it, p_com, it, nt, ng, false);
           cuda_record<<<(ng+511)/512, 512>>>(d_g1, p_obs_it, d_gzxy, ng, false);

	    cuda_step_fd3d<<<dimg,dimb>>>(d_g0, d_g1, d_vv, _dz2, _dx2, _dy2, nz, nx, ny, dt, d_pdt2, false);
	    ptr=d_g0; d_g0=d_g1; d_g1=ptr;
           cuda_absorb_bndr<<<dimg,dimb>>>(d_g0, nz, nx, ny, -0.25);
           cuda_absorb_bndr<<<dimg,dimb>>>(d_g1, nz, nx, ny, -0.25);


           cuda_light_matrix<<<dimg,dimb>>>(source, d_g1, nz, nx, ny);

         /*  if(is==0&&it%300==0){
	      t0 = clock();
	      cudaMemcpy(vv, d_g1, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	      window3d(v0, vv, nz, nx, ny);
	      sf_floatwrite(v0, nz*nx*ny, Fsnap);	  

 	      t1 = clock();
 	      sf_warning("save the volume: %f (s)", ((float)(t1-t0))/CLOCKS_PER_SEC); 	
	     } */
            /*< calculate the is gradient >*/
           cuda_cal_gradient<<<dimg,dimb>>>(g_is, d_pdt2, d_g1, nz, nx, ny);

	  }//back propagation it loop over


         /*< source lighting >*/
         cuda_IS_lighting<<<dimg,dimb>>>(g_is, source, nz, nx, ny);
         cuda_sum<<<dimg,dimb>>>(g_ns, g_is, nz, nx, ny);

	  cudaMemcpy(vv, g_is, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	  window3d(v0, vv, nz, nx, ny);
	  sf_floatwrite(v0, nz*nx*ny, Fg);

	  cudaEventRecord(stop);
         cudaEventSynchronize(stop);
  	  cudaEventElapsedTime(&mstimer, start, stop);
    	  sf_warning("%d shot finished: %g (s)",is+1, mstimer*1.e-3);

	}/*< IS shot loop end ! >*/

       sf_floatwrite(&obj,1,Fobj);

	  cudaMemcpy(vv, g_ns, nnz*nnx*nny*sizeof(float), cudaMemcpyDeviceToHost);
	  window3d(v0, vv, nz, nx, ny);
	  sf_floatwrite(v0, nz*nx*ny, Fg);

	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	/* free memory on device */
	cudaFree(d_wlt);
	cudaFree(d_vv);
      /*< wavefield(x-y-z) >*/
	cudaFree(d_p0);
	cudaFree(d_p1);
	cudaFree(d_g0);
	cudaFree(d_g1);
	cudaFree(d_pdt2);
      /*< location >*/
	cudaFree(d_szxy);
	cudaFree(d_gzxy);
	cudaFree(bndr_front);
	cudaFree(bndr_back);
	cudaFree(bndr_left);
	cudaFree(bndr_right);
	cudaFree(bndr_down);
      /*< gradient >*/
	cudaFree(g_is);
	cudaFree(g_ns);
	cudaFree(source);
      /*< wavefield(t-x-y-z) >*/
	cudaFree(p_cal_it);
	cudaFree(p_cal_nt);
	cudaFree(p_obs_it);
	cudaFree(p_obs_nt);
	cudaFree(p_com);
	cudaFree(d_bndr);
      /*< free alloc >*/
	free(v0);
	free(vv);
	free(p_IO);

    	exit (0);
}

//a#############################################################################################
__global__ void cuda_step_fd3d(float *p0, float *p1, float *vv, float _dz2, float _dx2, float _dy2, int n1, int n2, int n3, 
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
	 if (validw&&pdt) pdt2[outIndex]=( dt/sqrtf(vv[outIndex]) ) *(c1+c2+c3);
	
        if (validw) p0[outIndex]=2.0*p1[outIndex]-p0[outIndex]+vv[outIndex]*(c1+c2+c3);
    }

}

//a#############################################################################################
void sf_check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); 
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
//a#############################################################################################

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
__global__ void cuda_record(float *p, float *seis, int *gxz, int ng, bool record)//++++++++++++
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng)
        {
        if(record) seis[id]=p[gxz[id]];
        else  p[gxz[id]]=seis[id];
        }
}
//a#############################################################################################
__global__ void cuda_light_matrix(float *s, float *p, int nz, int nx, int ny)
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
__global__ void cuda_IS_lighting(float *g_is, float *source, int nz, int nx, int ny)
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
          if(id<nnz*nnx*nny&&source[id]!=0) g_is[id]/=source[id];
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
__global__ void cuda_cal_gradient(float *g_is, float *s, float *g, int nz, int nx, int ny)
/*< calculate is gradient >*/
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
          if(id<nnz*nnx*nny) g_is[id]+=s[id]*g[id];
        }    
}
//a#############################################################################################
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
//a#############################################################################################
__global__ void cuda_free_surface(float *p0,float *p1,int nz,int nx,int ny)
/*< absorb boundry condition >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
	int nny=ny+2*mm+2*npd;
       
       if(iz<mm+npd){
        for(iy=0;iy<nny;iy++) {
            id=iz+ix*nnz+iy*nnz*nnx;
            p0[id]=-p0[id+2*(mm+npd)-iz];
            p1[id]=-p1[id+2*(mm+npd)-iz];
        }}else if(iz==mm+npd){ p0[id]=0.0;p1[id]=0.0;}  
}
//a#############################################################################################
__global__ void cuda_trans_xy2txy(float *xy, float *txy, int it, int nt, int ng, bool output)
/*< set the positions of  geophones in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng)
        {
          if(output) txy[it+id*nt]=xy[id];
          else xy[id]=txy[it+id*nt];
        }
       
}
//a#############################################################################################
__global__ void cuda_cal_residuals(float *cal, float *obs, float *com, int nn, int nx, int ny, int nt, float obj)
{
    const int it = blockIdx.x * blockDim.x + threadIdx.x;//0--nt's thread:it
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix
	int id, iy;
       if(it<nt){
        for(iy=0;iy<ny;iy++) {
            id=it+ix*nt+iy*nt*nx;
            if (id<nn) 
              com[id]=cal[id] - obs[id];
              obj+=com[id]*com[id];
            }} 
    	

}
//a#############################################################################################
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
//a#############################################################################################
__global__ void cuda_set_s(int *szxy, int szbeg, int sxbeg, int sybeg, int jsz, int jsx, int jsy, int ns, int nz, int nx, int ny)
/*< set the positions of sources  in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
    	if (id<ns) szxy[id]=(szbeg+id*jsz+mm+npd)+nnz*(sxbeg+id*jsx+mm+npd)+nnz*nnx*(sybeg+id*jsy+mm+npd);
}
//a#############################################################################################
__global__ void cuda_set_up_do(int *gzxy, int *down, int ng, int nz, int nx, int ny)
/*< set the positions of  geophones & down in whole domain >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;
       int iy=id/nx;
       int ix=id%nx;
    	if (id<ng){
           gzxy[id]=(mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(iy+mm+npd);
           down[id]=(nz+mm+npd-1)+nnz*(ix+mm+npd)+nnz*nnx*(iy+mm+npd);
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
           front[id]=(iz+mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(mm+npd);
           back[id]=(iz+mm+npd)+nnz*(ix+mm+npd)+nnz*nnx*(ny+mm+npd-1);
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
           left[id]=(iz+mm+npd)+nnz*(mm+npd)+nnz*nnx*(iy+mm+npd);
           right[id]=(iz+mm+npd)+nnz*(nx+mm+npd-1)+nnz*nnx*(iy+mm+npd);
        }
}
//a#############################################################################################
/*__global__ void cuda_load_receiver(bool mute,float *bndr, float *p0, int *up, int nz, int nx, int ny, int it, float fm, 
                                   float dt, float dz, float dx, float dy,int is, int *szxy, float v_mu)*/
/*< load receivers wavefield >*/
/*{
	int id=threadIdx.x+blockIdx.x*blockDim.x;
	int nnz=nz+2*mm+2*npd;
	int nnx=nx+2*mm+2*npd;

       int sy= szxy[is]/(nnz*nnx)-mm-npd;
       int sx=(szxy[is]%(nnz*nnx))/nnz-mm-npd;
       int sz=(szxy[is]%(nnz*nnx))%nnz-mm-npd;

       int iy= up[id]/(nnz*nnx)-mm-npd;
       int ix=(up[id]%(nnz*nnx))/nnz-mm-npd;
       int iz=(up[id]%(nnz*nnx))%nnz-mm-npd;

		if(id<nx*ny) 
                {
                  p0[up[id]]=bndr[id+2*nz*nx+2*nz*ny]; */  /* up    boundary */
               /*   if(mute)
                    if(it<sqrtf(pow((ix-sx)*dx,2) + pow((iy-sy)*dy,2) + pow((iz-sz)*dz,2))/v_mu/dt+2.0/(dt*fm)+50)
                          p0[up[id]]=0.0;
                }

}*/
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
