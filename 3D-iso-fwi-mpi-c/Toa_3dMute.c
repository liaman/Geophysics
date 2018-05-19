#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
void main()
{
  float ***val,dx,dy,dz;
  int i,j,k,is;
  int nx,ny,nt;

  int mu_t,mu_nt;
  float mu_x,mu_y,mu_z,mu_t0,favg=60,dt=0.4;
 
  char FN1[250]={"shot.dat"};
  char FN2[250]={"shot_mute.dat"};


  int xs_sxd[100],ys_sxd[100];
	int ns_sxd,zs_sxd;

  nx=70;   dx=dy=dz=5;
  ny=70;   
  nt=725;   

                ns_sxd=28;
 
                for(i=0;i<14;i++)
                {  xs_sxd[i]=5+i*5;
                   ys_sxd[i]=1;}
                for(i=14;i<28;i++)
                {  xs_sxd[i]=69;
                   ys_sxd[i]=0+(i-14)*5;}
                zs_sxd=1;  


  val=alloc3float(nt,ny,nx);
  zero3float(val,nt,ny,nx);

  FILE *fp1,*fp2;
  fp1=fopen(FN1,"rb");
  fp2=fopen(FN2,"wb");
for(is=1;is<=ns_sxd;is++)
{
  printf("is==%d\n",is);
     fseek(fp1,(is-1)*nx*ny*nt*4L,0);
      for(j=0;j<ny;j++)
       for(i=0;i<nx;i++)
	  for(k=0;k<nt;k++)
	    fread(&val[i][j][k],4L,1,fp1);

      for(j=0;j<ny;j++)
       for(i=0;i<nx;i++)
        {
        mu_x=dx*abs(i-xs_sxd[is-1]);
        mu_y=dy*abs(j-ys_sxd[is-1]);
        mu_z=dz*zs_sxd;
        mu_t0=sqrtf(pow(mu_x,2)+pow(mu_y,2)+pow(mu_z,2))/3000;
        mu_t=(int)(2.0/(dt/1000*favg));
        mu_nt=(int)(mu_t0/dt*1000)+mu_t+20;
        for(k=0;k<nt;k++)if(k<mu_nt)
           val[i][j][k]=0.0;
        }
     fseek(fp2,(is-1)*nx*ny*nt*4L,0);
      for(j=0;j<ny;j++)
       for(i=0;i<nx;i++)
	  for(k=0;k<nt;k++)
	    fwrite(&val[i][j][k],4L,1,fp2);


}
  fclose(fp1);
  fclose(fp2);





     
  free3float(val);printf("done\n");
}







