/* 2-D finite-difference Laplacian operation. */

#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<string.h>

void laplac2_lop(int adj, int nz, int nx, float *in, float *out)
/*< linear operator >*/
{
    int iz,ix,j;

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

int main(int argc, char* argv[])
{
    int nz,nx,n3,i3;
    float *in,*out;
    int adj;
    char FN1[250]={"migration_stkADCIGs0-65.dat"};
    char FN2[250]={"migration_stkADCIGs0-65_laplace.dat"};

    FILE *fp1, *fp2;
    fp1  = fopen (FN1,"rb");
    fp2 = fopen (FN2,"wb");

    nz=301;
    nx=601;printf("----\n");
    n3 = 1;

    adj=1; /* adjoint flag */printf("----\n");

     in = (float*)malloc(nx*nz*sizeof(float));printf("----\n");
    out = (float*)malloc(nx*nz*sizeof(float));printf("----\n");

    fread(in,sizeof(float),nx*nz,fp1);printf("----2\n");
    for (i3=0; i3<n3; i3++) {
       laplac2_lop ( adj, nz, nx, in,  out );
    }
    fwrite(out,sizeof(float),nx*nz,fp2);
    exit(0);



}
