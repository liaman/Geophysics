

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

/*************func**************/ 
void adcigs_laplace_filter(int nx,int nz,int na,float *adcigs,int flag)
{
   int ix,iz,ia,id,ido;
   float *in, *out,*temp;

    	 in=(float*)malloc(nz*nx*sizeof(float));
    	out=(float*)malloc(nz*nx*sizeof(float));
    	temp=(float*)malloc(nz*nx*na*sizeof(float));

   if(flag==1||flag==3)
   for (ia=0; ia<na; ia++){
       for (ix=0; ix<nx; ix++){
           for (iz=0; iz<nz; iz++){
                   id=ix*na*nz+ia*nz+iz;
                   ido=ix*nz+iz;
                   in[ido]=adcigs[id];
             }
        }
       laplac2_lop( 1, nz, nx, in,  out );
       for (ix=0; ix<nx; ix++)  {
          for (iz=0; iz<nz; iz++) {
                   id=ix*na*nz+ia*nz+iz;
                   ido=ix*nz+iz;
                   temp[id]+=out[ido];
                   if(flag==1)adcigs[id]=temp[id];
            }
        }
   }
   if(flag!=1)
   for (ia=na-1; ia>=0; ia--) {
       for (ix=0; ix<nx; ix++) {
          for (iz=0; iz<nz; iz++) {
                   id=ix*na*nz+ia*nz+iz;
                   ido=ix*nz+iz;
                   in[ido]=adcigs[id];
           }
       }
      laplac2_lop( 1, nz, nx, in,  out );
      for (ix=0; ix<nx; ix++) {
          for (iz=0; iz<nz; iz++) {
                   id=ix*na*nz+ia*nz+iz;
                   ido=ix*nz+iz;
                   temp[id]+=out[ido];
                   if(flag==2||flag==3) adcigs[id]=temp[id];
            }
       }
   }
}

int main(int argc, char* argv[])
{
    int nz,nx,dcdp,na,ix,iz,ia,id,ido,nna;
    float *in,*out,*mig0,*mig1;

    char FN1[250]={"adcigs.dat"};
    char FN2[250]={"adcigs_laplace.dat"};

    FILE *fp1, *fp2;
    fp1  = fopen (FN1,"rb");
    fp2 = fopen (FN2,"wb");


    nx=801;
    nz=301;
    na=65;



     in = (float*)malloc(nx*nz*na*sizeof(float));printf("----\n");
    out = (float*)malloc(nx*nz*na*sizeof(float));printf("----\n");
   mig0 = (float*)malloc(nx*nz*sizeof(float));printf("----\n");
   mig1 = (float*)malloc(nx*nz*sizeof(float));printf("----\n");


    fread(in,sizeof(float),nx*nz*na,fp1);printf("----2\n");
    adcigs_laplace_filter(nx,nz,na,in,2);/*1:(0-na);2:(na-0);3:(0-na-0)*/
/*for (ia=0; ia<na; ia++) 
{
    for (ix=0; ix<nx; ix++) 
    {

               for (iz=0; iz<nz; iz++) 
                 {
                      id=ix*na*nz+ia*nz+iz;
                      ido=ix*nz+iz;
                      mig0[ido]=in[id];

                  }
    }
       laplac2_lop ( 1, nz, nx, mig0,  mig1 );printf("----3\n");

    for (ix=0; ix<nx; ix++) 
    {

               for (iz=0; iz<nz; iz++) 
                 {
                      id=ix*na*nz+ia*nz+iz;
                      ido=ix*nz+iz;
                      out[id]=mig1[ido];

                  }

    }
}*/printf("----4\n");
    fwrite(in,sizeof(float),nx*nz*na,fp2);
    exit(0);



}
