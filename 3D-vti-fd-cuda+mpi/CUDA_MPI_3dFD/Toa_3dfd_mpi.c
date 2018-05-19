//a#########################################################
//a##    3D iso acoustic fd :MPI + CUDA 
//a##                       code by Rong Tao
//a#########################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "mpi.h"

void cuda_3dfd(FILE *fpvel, FILE *fpsnap, FILE *fpshot, int is, int ns, int myid, int nx, int ny, int nz, 
              float dx, float dy, float dz, int sxbeg, int sybeg, int szbeg, int jsx, int jsy, int jsz, 
              int nt, int kt, float dt, float fm);

main(int argc,char *argv[])
{

	int myid, numprocs;
	int is, ns, nx, ny, nz, szbeg, sxbeg, sybeg, jsz, jsx, jsy, nt, kt;
       float dx, dy, dz, dt, fm;

       clock_t t0, t1;


	char FNvel[250]={"vel201202203.dat"};
       char FNsnap[250]={"snap.dat"};
       char FNshot[250]={"shot.dat"};

       FILE *fpvel,*fpsnap,*fpshot;
       fpsnap=fopen(FNsnap,"wb");
       fpshot=fopen(FNshot,"wb+");


    	nx=201;   dx=10;   sxbeg=10;   jsx=9;
    	ny=202;   dy=10;   sybeg=10;   jsy=9;
    	nz=203;   dz=10;   szbeg=1;    jsz=0;
    	
       nt=1301;kt=300;dt=0.001;

   	fm=15;	

   	ns=4;

/*******************************************/
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
      MPI_Barrier(MPI_COMM_WORLD);


      if(myid==0)t0 = clock();
   
      for(is=myid;is<ns;is=is+numprocs)
      {
       fpvel=fopen(FNvel,"rb");

       if(myid==0){
           printf("\n################################\n");
           printf("MPIforward: is==%d\n",is);}

       cuda_3dfd(fpvel, fpsnap, fpshot, is, ns, myid, nx, ny, nz, dx, dy, dz, 
                         sxbeg, sybeg, szbeg, jsx, jsy, jsz, nt, kt, dt, fm);
       fclose(fpvel);
      }

     if(myid==0)t1 = clock();
     if(myid==0)printf("####### MPI_Totally: %f (s)\n", ((float)(t1-t0)/1000000.0)); 

     MPI_Barrier(MPI_COMM_WORLD);


     fclose(fpsnap);
     fclose(fpshot);


     MPI_Finalize();
}

