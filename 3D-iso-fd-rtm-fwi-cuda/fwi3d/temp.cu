//a#############################################################################################
__global__ void cuda_cal_beta(float *beta, float *g0, float *g1, float *cg, int N)
/*< calculate beta for nonlinear conjugate gradient algorithm 
configuration requirement: <<<1,BlockSize>>> >*/
{
    const int iz = blockIdx.x * blockDim.x + threadIdx.x;//0--nz's thread:iz
    const int ix = blockIdx.y * blockDim.y + threadIdx.y;//0--nx's thread:ix

       int id,iy;

      float s=0.0f,t=0.0f,r=0.0f;

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
       __syncthreads();

       if(iz==0)
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
