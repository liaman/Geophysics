#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "mpi.h"
#include "time.h"
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"

main(int argc,char *argv[])
{


        void readfile(int nx,int nt,int is,char FN1[],float **p_real);
        void compare(int nx,int nt,int is,char FN1[],float **pp_real,float **p_initial,float **p_com);
	  void model(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,float favg,
                 float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
                   int is,int myid,float **p_initial,float **p_b,float **p_l,float **p_r,float **ws0,int wavelet);
        void model_mute_directwave(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,
                   float favg,float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int ds_sxd,int fs_sxd,
                    int zs_sxd,int is,int myid,float **p_initial,int wavelet);//++++++++++++

	  void back_grad(int myid,int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,
                   float favg,float tmax,float dt,float dtout,char FN1[],int ns_sxd,int ds_sxd,int is,int fs_sxd,
                   int zs_sxd,float **p_com,float **p_t,float **p_b,float **p_l,float **p_r,float **g);
	  float step(int nx,int nz,int vnx,int vnz,int nt,int is,int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,int myid,
                   int numprocs,float factor[],char FN1[],char FN2[],char FN4[],char FN5[],int npd,float dx,float dz,
                   float vdx,float vdz,float favg,float tmax,float dt,float dtout,float pfac,int wavelet,int mz,int vmz);
        void read_velfile(int vnx,int vnz,char FN1[],float **v);
        void zero(float **a,int nx,int nz);
        void window(int vnx,int vnz,int m,int n,float **g,float **g1);
        float square_sum(int nx,int nz,float **a);
        void freea(float **a,int m);

	char FN1[250]={"fault_shot_obs.dat"};
	char FN2[250]={"fault_vel_560_360_ite100.dat"};
	char FN3[250]={"fault_shot_cal.dat"};
	char FN4[250]={"multi_gradient.dat"};
	char FN5[250]={"cg2_gradient.dat"};


	int i,j,k,is,nx,nz,mz,vmz,nt,vnx,vnz,i_start,i_end,wavelet;
	int ns_sxd,ds_sxd,fs_sxd,zs_sxd,npd;

	float dx,dz,vdx,vdz,tmax,dt,dtout,pfac,favg;
	float wsmax;
	float a,alph,beta;float factor[1];
	int myid,numprocs,count;

       clock_t time_start,time_end;
       double time_total;

         float E_com1=0.0,E_com11=0.0,E_com=0.0;





   /*     nx=100;         npd=50;   tmax=0.7;
	nz=100;         favg=20;  pfac=1000.0;
 	                dx=5.0;
	vnx=100;        dz=5.0;
	vnz=100;        vdx=5.0;
	nt=1401;        vdz=5.0;
        ns_sxd=11;       dt=0.5;
        ds_sxd=10;       dtout=0.5;
        fs_sxd=1;
        zs_sxd=1;

        mz=vmz=9;//for mute
   */
         printf("-------------------------\n");


         nx=560;         npd=50;   tmax=1.8;
	  nz=360;         favg=15;  pfac=1000.0;
 	                  dx=5.0;   
	  vnx=560;        dz=5.0;   
	  vnz=360;        vdx=5.0;
	  nt=3000;        vdz=5.0;
        ns_sxd=40;       dt=0.6;
        ds_sxd=10;       dtout=0.6;
        fs_sxd=80;
        zs_sxd=2;        
 
        mz=vmz=20;//++++for cute  

           wavelet=0;//wavelet==0:use initial wavelet
                     //wavelet==1:read wavelet



	float **g;
	float **g0;
	float **ggg0;
	float **prp;
	float **new_g;
	float **v_initial;
        float **p_initial;
        float **p_b;
        float **p_l;
        float **p_r;
        float **p_com;
	float **pp_real;
	float **ws;
	float **ws0;
        float **p_m;//+++++++++++++++++++++

/******************add the tempor vel_model**********///+++++++++++++++++++++++(
        float **vtempor,**vmute;

          vtempor=(float **)calloc((vnx),sizeof(float *));
	  {
		  for(i=0;i<vnx;i++)
		  {
		  vtempor[i]=(float *)calloc((vnz),sizeof(float));
		  }
	  }
        vmute=(float **)calloc((vnx),sizeof(float *));
	  {
		  for(i=0;i<vnx;i++)
		  {
		  vmute[i]=(float *)calloc((vmz),sizeof(float));
		  }
	  }

           FILE *fptempor;
		 fptempor=fopen(FN2,"rb");
		 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<nz;j++)
			 {
				 fread(&vtempor[i][j],4,1,fptempor);

			 }
		 }

                 fclose(fptempor);
                 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<mz;j++)
			 {
				 vmute[i][j]=vtempor[i][j];

			 }
		 }

                 fptempor=fopen("Temp_vel_for_mute_direction.dat","wb");
		 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<mz;j++)
			 {
				 fwrite(&vmute[i][j],4,1,fptempor);

			 }
		 }
                 fclose(fptempor);


	p_m=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		p_m[i]=(float *)calloc(nt,sizeof(float));
	}
/*********************************************************///+++++++++++++++++++++++++++++)
	p_initial=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		p_initial[i]=(float *)calloc(nt,sizeof(float));
	}

	p_b=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		p_b[i]=(float *)calloc(nt,sizeof(float));
	}

	p_l=(float **)calloc(nz,sizeof(float *));
	for(i=0;i<nz;i++)
	{
		p_l[i]=(float *)calloc(nt,sizeof(float));
	}

	p_r=(float **)calloc(nz,sizeof(float *));
	for(i=0;i<nz;i++)
	{
		p_r[i]=(float *)calloc(nt,sizeof(float));
	}

	p_com=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		p_com[i]=(float *)calloc(nt,sizeof(float));
	}

	pp_real=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		pp_real[i]=(float *)calloc(nt,sizeof(float));
	}


	v_initial=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		v_initial[i]=(float *)calloc(vnz,sizeof(float));
	}

	g=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		g[i]=(float *)calloc(vnz,sizeof(float));
	}

	new_g=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		new_g[i]=(float *)calloc(vnz,sizeof(float));
	}

	ggg0=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		ggg0[i]=(float *)calloc(vnz,sizeof(float));
	}

	prp=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		prp[i]=(float *)calloc(vnz,sizeof(float));
	}


	g0=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		g0[i]=(float *)calloc(vnz,sizeof(float));
	}

	ws=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		ws[i]=(float *)calloc(vnz,sizeof(float));
	}
	ws0=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		ws0[i]=(float *)calloc(vnz,sizeof(float));
	}
/*******************************************/
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);



	zero(v_initial,vnx,vnz);
	zero(g,vnx,vnz);
	zero(g0,vnx,vnz);
	zero(p_initial,nx,nt);
	zero(new_g,vnx,vnz);
	zero(prp,vnx,vnz);
	zero(ws,vnx,vnz);
	zero(ws0,vnx,vnz);

	for(i=0;i<vnx;i++)
	{
	    for(j=0;j<vnz;j++)
	    {
	        ggg0[i][j]=1.0;
	    }
	}

/******************************************/
      MPI_Barrier(MPI_COMM_WORLD);
      zero(new_g,vnx,vnz);
/******************************************************************************/
      count=39;


      FILE *fpjj;
      fpjj=fopen("Velocity_IteRecord.dat","wb");

do
{

    E_com1=0.0;
    E_com11=0.0;
     E_com=0.0;

    if(myid==0)
    {
      time_start=clock();


      fseek(fpjj,(count-1)*nx*nz*4L,0);
                 for(i=0;i<nx;i++)
		 {
			 for(j=0;j<nz;j++)
			 {
				 fwrite(&v_initial[i][j],4L,1,fpjj);
			 }
		 }


        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("---\n");
        printf("---\n");
        printf("---\n");
        printf("---   The count=%d iterative Begin !!!!!!!! \n",count);
    }
/*************************************/




    FILE *fp3;
    fp3=fopen(FN3,"wb");
    FILE *fp2;
    fp2=fopen("shot_compare_m.dat","wb");
    FILE *fp4;
    fp4=fopen("Temp_grad_m.dat","wb");
    FILE *fp444;
    fp444=fopen("Temp_ws2_m.dat","wb");
/*************************************/
    a=0.0;
    factor[0]=0.0;
    zero(v_initial,vnx,vnz);
    zero(g,vnx,vnz);
    zero(g0,vnx,vnz);
    zero(ws,vnx,vnz);
    zero(ws0,vnx,vnz);
    MPI_Barrier(MPI_COMM_WORLD);
/**************begin the loop******************************/
    for(is=1+myid;is<=ns_sxd;is=is+numprocs)
    {
        if(myid==0)
        {
             printf("---   begin the loop ! \n");
        }
        MPI_Barrier(MPI_COMM_WORLD);

         E_com1=0.0;
/***************************************/
        zero(p_initial,nx,nt);
        zero(p_m,nx,nt);//++++++++++++++++++++++
        zero(pp_real,nx,nt);
        zero(p_b,nx,nt);
        zero(p_l,nz,nt);
        zero(p_r,nz,nt);
	zero(p_com,nx,nt);
	zero(g0,vnx,vnz);
	zero(ws0,vnx,vnz);
/***************************************/

        if(myid==0)
	{
	  printf("---   the  model  is start  !   \n");
	}
        model(nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,FN2,
                ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,myid,p_initial,p_b,p_l,p_r,ws0,wavelet);
        if(myid==0)
	{
	  printf("---   mute the direct wave is start !    \n");
	}
        model_mute_directwave(nx,mz,vnx,vmz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,"Temp_vel_for_mute_direction.dat",
                ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,myid,p_m,wavelet);//+++++++

         if(myid==0)
        {
             printf("---   mute the direct wave is over !  \n");
        }

        for(i=0;i<nx;i++)//++++++++++++mute the direct wave
	{
	  for(j=0;j<nt;j++)
	  {
	    p_m[i][j]=p_initial[i][j]-p_m[i][j];
	  }
	}

        fseek(fp3,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	{
	  for(j=0;j<nt;j++)
	  {
	    fwrite(&p_m[i][j],4,1,fp3);//++++++++++++change 'p_initial' to 'p_m'
	  }
	}

	fseek(fp444,(is-1)*nx*nz*4L,0);
	for(i=0;i<nx;i++)
	{
	  for(j=0;j<nz;j++)
	  {
	    fwrite(&ws0[i][j],4,1,fp444);
	  }
	}

	if(myid==0)
	{
	  printf("---   the  model  is over  !   \n");
          printf("---   the compare is start  !   \n");
	}

        MPI_Barrier(MPI_COMM_WORLD);
        compare(nx,nt,is,FN1,pp_real,p_m,p_com);//++++++++++++change 'p_initial' to 'p_m'

        fseek(fp2,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			fwrite(&p_com[i][j],4,1,fp2);
		}
	}
          E_com1=square_sum(nx,nt,p_com);
           E_com11=E_com11+E_com1;


	if(myid==0)
	{
          printf("---   the compare is over  !  \n");
          printf("---   the back propagation is start  !    \n");
	}
/**************************************************************************/
        back_grad(myid,nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,FN2,
                   ns_sxd,ds_sxd,is,fs_sxd,zs_sxd,p_com,p_initial,p_b,p_l,p_r,g0);



	if(myid==0)
	{
          printf("---   the back propagation is over  !   \n");
	}



	fseek(fp4,(is-1)*vnx*vnz*4L,0);
	for(i=0;i<vnx;i++)
	{
	  for(j=0;j<vnz;j++)
	  {
	     fwrite(&g0[i][j],4,1,fp4);
	  }
	}

    }
/**********************************************************************************************/
/****************  if you do not use the complicated step**************/
/***************** then will out put the shot_compare enerage**********/
/**********************************************************************************************/

      MPI_Barrier(MPI_COMM_WORLD);

      MPI_Reduce(&E_com11,&E_com,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);


/************** the loop is done***************************/
      if(myid==0)
      {

           FILE *fpE;
           fpE=fopen("E_shot_compare.txt","a");
           fprintf(fpE,"%f\n",E_com);
           fclose(fpE);



          printf("---   the loop is done\n");
      }

      MPI_Barrier(MPI_COMM_WORLD);

      fclose(fp2);
      fclose(fp3);
      fclose(fp4);
      fclose(fp444);
/******************************************/

   if(myid==0)
   {
      FILE *fp5;
      fp5=fopen("Temp_grad_m.dat","rb");
      FILE *fp555;
      fp555=fopen("Temp_ws2_m.dat","rb");
      FILE *fp6;
      fp6=fopen(FN4,"wb");
      FILE *fp666;
      fp666=fopen("ws2_m.dat","wb");

      for(k=0;k<ns_sxd;k++)
      {
          zero(g0,vnx,vnz);
          zero(ws0,vnx,vnz);
          fseek(fp5,k*vnx*vnz*4L,0);
          for(i=0;i<vnx;i++)
          {
               for(j=0;j<vnz;j++)
               {
                    fread(&g0[i][j],4,1,fp5);
               }
          }
          fseek(fp555,k*vnx*vnz*4L,0);
          for(i=0;i<vnx;i++)
          {
               for(j=0;j<vnz;j++)
               {
                    fread(&ws0[i][j],4,1,fp555);
               }
          }

          for(i=0;i<vnx;i++)
          {
               for(j=0;j<vnz;j++)
               {
                   g[i][j]=g[i][j]+g0[i][j];
               }
          }
          for(i=0;i<vnx;i++)
          {
               for(j=0;j<vnz;j++)
               {
                   ws[i][j]=ws[i][j]+ws0[i][j];
               }
          }
      }
      wsmax=0.0;
      for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{
			if(wsmax<fabs(ws[i][j]))
			{
			   wsmax=fabs(ws[i][j]);
		      }
		}
	}

      for(i=0;i<vnx;i++)
      {
          for(j=0;j<vnz;j++)
          {
               //ws[i][j]=sqrt(ws[i][j]);
               ws[i][j]=ws[i][j]/wsmax;
               g[i][j]=g[i][j]/ws[i][j];
          }
      }
       for(i=0;i<vnx;i++)
      {
          for(j=0;j<vnz;j++)
          {
               fwrite(&g[i][j],4,1,fp6);
               fwrite(&ws[i][j],4,1,fp666);
          }
      }
      fclose(fp5);
      fclose(fp555);
      fclose(fp6);
      fclose(fp666);

      printf("---   complete the global gradient  !  \n");


      //////////////conjugate process////////////////
      alph=0.0;
      zero(prp,vnx,vnz);


      for(i=0;i<vnx;i++)
      {
            for(j=0;j<vnz;j++)
            {
                prp[i][j]=(g[i][j]-ggg0[i][j])*g[i][j];
            }
      }

      for(i=0;i<vnx;i++)
      {
            for(j=0;j<vnz;j++)
            {
               alph=alph+prp[i][j];
            }
      }


      beta=alph/(square_sum(vnx,vnz,ggg0));

      for(i=0;i<vnx;i++)
      {
            for(j=0;j<vnz;j++)
            {
               new_g[i][j]=g[i][j]+beta*new_g[i][j];
            }
      }


      for(i=0;i<vnx;i++)
      {
            for(j=0;j<vnz;j++)
            {
                ggg0[i][j]=g[i][j];
            }
      }
   }

      if(myid==0)
      {
         FILE *fp123;
         fp123=fopen(FN5,"wb");
         for(i=0;i<vnx;i++)
         {
                for(j=0;j<vnz;j++)
                {
                   fwrite(&new_g[i][j],4,1,fp123);
                }
         }

         fclose(fp123);
      }


      MPI_Barrier(MPI_COMM_WORLD);
/******************************step " a " ************************************/

    // a=step(nx,nz,vnx,vnz,nt,is,ns_sxd,ds_sxd,fs_sxd,zs_sxd,myid,numprocs,factor,FN2,FN5,FN1,FN3,
     //           npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,wavelet,mz,vmz);
      //试探步长改为2/1000和4/1000

      if(myid==0)
      {
         read_velfile(vnx,vnz,FN2,v_initial);
         zero(g0,vnx,vnz);
         window(vnx,vnz,20,vnz,new_g,g0);
   /**************************************************************///+++++++++++++++++get the easy "a"
         float vamin,gamax;
         vamin=9999999.0;
         gamax=g0[vnx-1][vnz-1];
       for(i=0;i<vnx;i++)
         {
                for(j=32;j<vnz;j++)
                {
                   if(vamin>v_initial[i][j])
                        vamin=v_initial[i][j];
                  if(gamax<g0[i][j])
                        gamax=g0[i][j];
                }
         }
        a=1.0*vamin/((280+1.5*count)*gamax);//+++++++++++++++++++++
       FILE *fpaaa;
      fpaaa=fopen("new_A_m.txt","a");
      fprintf(fpaaa,"vamin=%f,gamax=%f,a=%f \n",vamin,gamax,a);
      fclose(fpaaa);
   /**************************************************************///+++++++++++++++++get the easy "a"
         for(i=0;i<vnx;i++)
         {
             for(j=0;j<vnz;j++)
             {
                 v_initial[i][j]=v_initial[i][j]-a*g0[i][j];
             }
         }


         FILE *fp7;
         fp7=fopen(FN2,"wb");
         for(i=0;i<vnx;i++)
         {
             for(j=0;j<vnz;j++)
             {
                 fwrite(&v_initial[i][j],4,1,fp7);
             }
         }

         fclose(fp7);
      }

      count=count+1;



  if(myid==0)
     {
         FILE *fptime;
       fptime=fopen("Time_Record.txt","a");
       time_end=clock();
        time_total=(double)(time_end-time_start) / CLOCKS_PER_SEC;
        fprintf(fptime,"count=%4d  ,   ite_time=%f \n",count,time_total/60.0);
          fclose(fptime);
      }


}
while(1>0.5);//+++++++++++++++++++change
//while(fabs(factor[0])>0.5);

      fclose(fpjj);
       // fclose(fptime);

      freea(g,vnx);
      freea(g0,vnx);
      freea(ws,vnx);
      freea(ws0,vnx);
      freea(ggg0,vnx);
      freea(prp,vnx);
      freea(new_g,vnx);
      freea(v_initial,vnx);
      freea(p_initial,nx);
      freea(p_b,nx);
      freea(p_l,nz);
      freea(p_r,nz);
      freea(p_com,nx);
      freea(pp_real,nx);

      count=count-1;
      if(myid==0)
      {
          printf("---   complete the count=%d ite!!!!!!!!! \n",count);
          printf("---\n");
          printf("---\n");
          printf("---\n");
      }
	MPI_Finalize();


}


/**************************************************************************************/
void compare(int nx,int nt,int is,char FN1[],float **pp_real,float **p_initial,float **p_com)
{
	void freea(float **a,int m);
	int i,j;

	FILE *fp1;

	for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
	        pp_real[i][j]=0.0;
		}
	}

	fp1=fopen(FN1,"rb");
        fseek(fp1,(is-1)*nx*nt*4L,0);

	for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			fread(&pp_real[i][j],4,1,fp1);

		}
	}


	for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
		   p_com[i][j]=p_initial[i][j]-pp_real[i][j];

		}
	}

	fclose(fp1);


}

/************************************************************************/
void zero(float **a,int nx,int nz)
{
	int i,j;
  	for(i=0;i<nx;i++)
	{
		for(j=0;j<nz;j++)
		{
			a[i][j]=0.0;
		}
	}
}
/**********************************************************************/
void read_velfile(int vnx,int vnz,char FN1[],float **v)
{
	int i,j;
	FILE *fp1;

	fp1=fopen(FN1,"rb");

	for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{
			fread(&v[i][j],4,1,fp1);
		}
	}

	fclose(fp1);

}

/*******************************************************************************************/
void model(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,float favg,
           float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
            int is,int myid,float **p_initial,float **p_b,float **p_l,float **p_r,float **ws0,int wavelet)

{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                    int wtype,float pi,int npd,int is,int ds_sxd);
	void ptsource_wavelet(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                    int wtype,float pi,int npd,int is,int ds_sxd,int nt,int it);//++++++++++++++
      void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,float **u0,float **w0,float **txx0,
                    float **u1,float **w1,float **txx1,float **rho,float c[],
                    float *coffx1,float *coffx2,float *coffz1,float *coffz2);
      void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                    float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,float **tx1,
                    float **tx0,float **tz1,float **tz0,
                    float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2);
      void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,float favg,
                      float dtout,float dt,float **vp0,float ndtt);
      void zero(float **a,int nx,int nz);
      void zero2d(int nx,int nz,float **vv0,int npd);
	void zero2d2(int nx,int nz,int npd,float **vv);
	void zero2d3(int nx,int nz,int npd,float **vv);
	void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,
                  float **u1,float **w1,float **txx1,float **vp,float **rho);
      void pad_vv(int nx,int nz,int npd,float **ee);
	void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);
	void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,int npd,
                      int vnx,int vnz,int ds_sxd,int is);
	void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                           float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);

        void freea(float **a,int m);




/***********************************main****************/

	int i,j;
        int mm=4;
	int ntout,wtype,it,ifx,ilx,jfz,jlz,hsx;
	float pi,hsz,vsx,t,ndtt,d0;
	float iabsorb[4];

        hsx=1;

        wtype=1;
        iabsorb[0]=1.0;
        iabsorb[1]=1.0;
        iabsorb[2]=1.0;
        iabsorb[3]=1.0;

	pi=3.14;
	ndtt=dtout/dt;
	ntout=(int)(tmax*1000/dtout+0.5)+1;



	  float **vp0;
	  float **rho0;
	  float **vp;
	  float **rho;
	  float **u0;
	  float **w0;
	  float **txx0;
	  float **u1;
	  float **w1;
	  float **txx1;
	  float **tx0;
	  float **tx1;
	  float **tz0;
	  float **tz1;
	  float **s;

	  float c[4];


          cal_c(mm,c);


	  vp0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  vp0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



	  rho0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  rho0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



        zero2d(vnx,vnz,vp0,npd);

        zero2d(vnx,vnz,rho0,npd);

        read_file(FN1,vnx,vnz,vp0,rho0,npd);



	  vp=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  vp[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


	  rho=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  rho[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


	 u0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 u1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 w0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	  s=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  s[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }

	 d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);

/***************************************************************/

         float *coffx1;float *coffx2;float *coffz1;float *coffz2;float *acoffx1;float *acoffx2;float *acoffz1;float *acoffz2;
         coffx1=(float *)calloc((nx+2*npd),sizeof(float));
         coffx2=(float *)calloc((nx+2*npd),sizeof(float));
	 coffz1=(float *)calloc((nz+2*npd),sizeof(float));
         coffz2=(float *)calloc((nz+2*npd),sizeof(float));
	 acoffx1=(float *)calloc((nx+2*npd),sizeof(float));
	 acoffx2=(float *)calloc((nx+2*npd),sizeof(float));
	 acoffz1=(float *)calloc((nz+2*npd),sizeof(float));
	 acoffz2=(float *)calloc((nz+2*npd),sizeof(float));

	 for(i=0;i<nx+2*npd;i++)
	 {
		 coffx1[i]=0.0;
		 coffx2[i]=0.0;
		 acoffx1[i]=0.0;
		 acoffx2[i]=0.0;
	 }

	 for(i=0;i<nz+2*npd;i++)
	 {
		 coffz1[i]=0.0;
		 coffz2[i]=0.0;
		 acoffz1[i]=0.0;
		 acoffz2[i]=0.0;
	 }


       dt=dt/1000;

       initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/************************************************************/


       ndtt=(int)ndtt;

/******************begin wave propagation******************/
      zero2d2(nx,nz,npd,tx0);
      zero2d2(nx,nz,npd,tx1);
      zero2d2(nx,nz,npd,tz0);
      zero2d2(nx,nz,npd,tz1);
      zero3d(nx,nz,nt,npd,u0,w0,txx0,u1,w1,txx1,vp,rho);

      current_shot(vp0,rho0,vp,rho,nx,nz,npd,vnx,vnz,ds_sxd,is);

      pad_vv(nx,nz,npd,vp);

      pad_vv(nx,nz,npd,rho);


       for(i=0;i<nz+2*npd;i++)
       {
		for(j=0;j<nx+2*npd;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
	           rho[i][j]=1.0/rho[i][j];
		}
	}



      t=0.0;



      for(it=0;it<nt;it++)
      {

            t=t+dt;

/********************************************************/
          if(wavelet==0)
            ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd);
          if(wavelet==1)
            ptsource_wavelet(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd,nt,it);//+++++++++++++++++++++++++++
            update_vel(nx,nz,npd,mm,dt,dx,dz,u0,w0,txx0,u1,w1,txx1,rho,c,coffx1,coffx2,coffz1,coffz2);
            update_txx(nx,nz,dt,dx,dz,mm,u0,w0,txx0,u1,w1,txx1,s,vp,c,npd,tx1,tx0,tz1,tz0,acoffx1,acoffx2,acoffz1,acoffz2);


//////////////////////////////////////////////////////////////
		for(j=npd;j<npd+nx;j++)
		{
			p_initial[j-npd][it]=txx1[npd+1-1][j];

		}


		for(j=npd;j<npd+nx;j++)
		{
			p_b[j-npd][it]=txx1[npd+nz-1][j];

		}

		for(i=npd;i<npd+nz;i++)
		{
			p_l[i-npd][it]=txx1[i][npd+1-1];

		}

		for(i=npd;i<npd+nz;i++)
		{
			p_r[i-npd][it]=txx1[i][npd+nx-1];

		}
/////////////////////////////////////////////////////////////////

		for(i=0;i<nz+2*npd;i++)
		{
			for(j=0;j<nx+2*npd;j++)
			{
				u0[i][j]=u1[i][j];
				w0[i][j]=w1[i][j];
			    	tx0[i][j]=tx1[i][j];
			    	tz0[i][j]=tz1[i][j];
				txx0[i][j]=txx1[i][j];
			}
		}

       for(i=npd;i<npd+nz;i++)
		{
	    	      for(j=npd;j<nx+npd;j++)
			{
			      ws0[j-npd][i-npd]=ws0[j-npd][i-npd]+(txx1[i][j]*txx1[i][j]*txx1[i][j]*txx1[i][j]);
			}
		}
      }

/**********************release*******************/

	  free(coffx1);free(coffx2);free(coffz1);free(coffz2);
	  free(acoffx1);free(acoffx2);free(acoffz1);free(acoffz2);
	  freea(vp0,(vnx+2*npd));
	  freea(rho0,(vnx+2*npd));
	  freea(vp,(nz+2*npd));
	  freea(rho,(nz+2*npd));
	  freea(u0,(nz+2*npd));
	  freea(w0,(nz+2*npd));
	  freea(txx0,(nz+2*npd));
	  freea(u1,(nz+2*npd));
	  freea(w1,(nz+2*npd));
	  freea(txx1,(nz+2*npd));
	  freea(tx0,(nz+2*npd));
	  freea(tx1,(nz+2*npd));
	  freea(tz0,(nz+2*npd));
	  freea(tz1,(nz+2*npd));
	  freea(s,(nz+2*npd));

}

/****************************************************************************************/
void model_mute_directwave(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,
                float favg,float tmax,float dt,float dtout,float pfac,char FN1[],int ns_sxd,int ds_sxd,
                int fs_sxd,int zs_sxd,int is,int myid,float **p_initial,int wavelet)

{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                     int wtype,float pi,int npd,int is,int ds_sxd);
	void ptsource_wavelet(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                     int wtype,float pi,int npd,int is,int ds_sxd,int nt,int it);//++++++++++++++
      void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,float **u0,float **w0,float **txx0,
                     float **u1,float **w1,float **txx1,float **rho,float c[],float *coffx1,float *coffx2,
                     float *coffz1,float *coffz2);
      void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                     float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,
                     float **tx1,float **tx0,float **tz1,float **tz0,
                      float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2);
      void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,float favg,
                          float dtout,float dt,float **vp0,float ndtt);
      void zero(float **a,int nx,int nz);
      void zero2d(int nx,int nz,float **vv0,int npd);
	void zero2d2(int nx,int nz,int npd,float **vv);
	void zero2d3(int nx,int nz,int npd,float **vv);
	void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,
                       float **u1,float **w1,float **txx1,float **vp,float **rho);
      void pad_vv(int nx,int nz,int npd,float **ee);
	void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);
	void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,int npd,
                            int vnx,int vnz,int ds_sxd,int is);
	void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                           float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);

      void freea(float **a,int m);





/***********************************main****************/

	int i,j;
        int mm=4;
	int ntout,wtype,it,ifx,ilx,jfz,jlz,hsx;
	float pi,hsz,vsx,t,ndtt,d0;
	float iabsorb[4];

        hsx=1;

        wtype=1;
        iabsorb[0]=1.0;
        iabsorb[1]=1.0;
        iabsorb[2]=1.0;
        iabsorb[3]=1.0;

	pi=3.14;
	ndtt=dtout/dt;
	ntout=(int)(tmax*1000/dtout+0.5)+1;



	  float **vp0;
	  float **rho0;
	  float **vp;
	  float **rho;
	  float **u0;
	  float **w0;
	  float **txx0;
	  float **u1;
	  float **w1;
	  float **txx1;
	  float **tx0;
	  float **tx1;
	  float **tz0;
	  float **tz1;
	  float **s;

	  float c[4];


          cal_c(mm,c);


	  vp0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  vp0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



	  rho0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  rho0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



        zero2d(vnx,vnz,vp0,npd);

        zero2d(vnx,vnz,rho0,npd);

        read_file(FN1,vnx,vnz,vp0,rho0,npd);



	  vp=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  vp[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


	  rho=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  rho[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


	 u0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 u1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 w0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	  s=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  s[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }

	 d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);

/***************************************************************/

         float *coffx1;float *coffx2;float *coffz1;float *coffz2;float *acoffx1;float *acoffx2;float *acoffz1;float *acoffz2;
         coffx1=(float *)calloc((nx+2*npd),sizeof(float));
         coffx2=(float *)calloc((nx+2*npd),sizeof(float));
	 coffz1=(float *)calloc((nz+2*npd),sizeof(float));
         coffz2=(float *)calloc((nz+2*npd),sizeof(float));
	 acoffx1=(float *)calloc((nx+2*npd),sizeof(float));
	 acoffx2=(float *)calloc((nx+2*npd),sizeof(float));
	 acoffz1=(float *)calloc((nz+2*npd),sizeof(float));
	 acoffz2=(float *)calloc((nz+2*npd),sizeof(float));

	 for(i=0;i<nx+2*npd;i++)
	 {
		 coffx1[i]=0.0;
		 coffx2[i]=0.0;
		 acoffx1[i]=0.0;
		 acoffx2[i]=0.0;
	 }

	 for(i=0;i<nz+2*npd;i++)
	 {
		 coffz1[i]=0.0;
		 coffz2[i]=0.0;
		 acoffz1[i]=0.0;
		 acoffz2[i]=0.0;
	 }


       dt=dt/1000;

       initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/************************************************************/


       ndtt=(int)ndtt;

/******************begin wave propagation******************/
      zero2d2(nx,nz,npd,tx0);
      zero2d2(nx,nz,npd,tx1);
      zero2d2(nx,nz,npd,tz0);
      zero2d2(nx,nz,npd,tz1);
      zero3d(nx,nz,nt,npd,u0,w0,txx0,u1,w1,txx1,vp,rho);

      current_shot(vp0,rho0,vp,rho,nx,nz,npd,vnx,vnz,ds_sxd,is);

      pad_vv(nx,nz,npd,vp);

      pad_vv(nx,nz,npd,rho);


       for(i=0;i<nz+2*npd;i++)
       {
		for(j=0;j<nx+2*npd;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
	           rho[i][j]=1.0/rho[i][j];
		}
	}



      t=0.0;



      for(it=0;it<nt;it++)
      {

            t=t+dt;

/********************************************************/
          if(wavelet==0)
            ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd);
          if(wavelet==1)
            ptsource_wavelet(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd,nt,it);//+++++++++++++++++++++++++++
            update_vel(nx,nz,npd,mm,dt,dx,dz,u0,w0,txx0,u1,w1,txx1,rho,c,coffx1,coffx2,coffz1,coffz2);
            update_txx(nx,nz,dt,dx,dz,mm,u0,w0,txx0,u1,w1,txx1,s,vp,c,npd,tx1,tx0,tz1,tz0,acoffx1,acoffx2,acoffz1,acoffz2);


//////////////////////////////////////////////////////////////
		for(j=npd;j<npd+nx;j++)
		{
			p_initial[j-npd][it]=txx1[npd+1-1][j];

		}
/////////////////////////////////////////////////////////////////

		for(i=0;i<nz+2*npd;i++)
		{
			for(j=0;j<nx+2*npd;j++)
			{
				u0[i][j]=u1[i][j];
				w0[i][j]=w1[i][j];
			    	tx0[i][j]=tx1[i][j];
			    	tz0[i][j]=tz1[i][j];
				txx0[i][j]=txx1[i][j];
			}
		}
      }

/**********************release*******************/

	  free(coffx1);free(coffx2);free(coffz1);free(coffz2);
	  free(acoffx1);free(acoffx2);free(acoffz1);free(acoffz2);
	  freea(vp0,(vnx+2*npd));
	  freea(rho0,(vnx+2*npd));
	  freea(vp,(nz+2*npd));
	  freea(rho,(nz+2*npd));
	  freea(u0,(nz+2*npd));
	  freea(w0,(nz+2*npd));
	  freea(txx0,(nz+2*npd));
	  freea(u1,(nz+2*npd));
	  freea(w1,(nz+2*npd));
	  freea(txx1,(nz+2*npd));
	  freea(tx0,(nz+2*npd));
	  freea(tx1,(nz+2*npd));
	  freea(tz0,(nz+2*npd));
	  freea(tz1,(nz+2*npd));
	  freea(s,(nz+2*npd));

}

/*******************************************************************************************/
void back_grad(int myid,int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,
              float favg,float tmax,float dt,float dtout,char FN1[],int ns_sxd,int ds_sxd,int is,
              int fs_sxd,int zs_sxd,float **p_com,float **p_t,float **p_b,float **p_l,float **p_r,float **g)
{

      void cal_c(int mm,float c[]);

      void for_update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                       float **u0,float **w0,float **txx0,float **u1,float **w1,float **txx1,float **rho,
                       float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2);

      void for_update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                      float **u1,float **w1,float **txx1,float **vp,float c[],int npd,
                      float **tx1,float **tx0,float **tz1,float **tz0,
                      float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2);

      void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,
                      float **u0,float **w0,float **txx0,float **u1,float **w1,float **txx1,float **rho,
                      float c[],float *coffx1,float *coffx2,float *coffz1,float *coffz2);

      void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                      float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,
                       float **tx1,float **tx0,float **tz1,float **tz0,
                      float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2);

      void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,float favg,
                       float dtout,float dt,float **vp0,float ndtt);
      void zero(float **a,int nx,int nz);
      void zero2d(int nx,int nz,float **vv0,int npd);
      void zero2d2(int nx,int nz,int npd,float **vv);
      void zero2d3(int nx,int nz,int npd,float **vv);
      void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,
                         float **u1,float **w1,float **txx1,float **vp,float **rho);
      void pad_vv(int nx,int nz,int npd,float **ee);
      void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);
      void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,int npd,
                          int vnx,int vnz,int ds_sxd,int is);
      void read_velfile(int vnx,int vnz,char FN1[],float **v);
	void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                           float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);
      void freea(float **a,int m);


	int mm=4;
	int hsx=1;
	int i,j;

	int wtype,it,ifx,ilx,jfz,jlz,ntout;
	float pi,d0,ndtt;
	float iabsorb[4];

      wtype=1;
      iabsorb[0]=1.0;
      iabsorb[1]=1.0;
      iabsorb[2]=1.0;
      iabsorb[3]=1.0;

	pi=3.14;
	ndtt=dtout/dt;
        ntout=(int)(tmax*1000/dtout+0.5)+1;



	  float **vp0;
	  float **rho0;
	  float **vp;
	  float **rho;
	  float **u0_f;    float **u0_b;
	  float **w0_f;    float **w0_b;
	  float **txx0_f;  float **txx0_b;
	  float **u1_f;    float **u1_b;
	  float **w1_f;    float **w1_b;
	  float **txx1_f;  float **txx1_b;
	  float **tx0_f;   float **tx0_b;
	  float **tx1_f;   float **tx1_b;
	  float **tz0_f;   float **tz0_b;
	  float **tz1_f;   float **tz1_b;
          float **s;

	  float c[4];


          cal_c(mm,c);


	  vp0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  vp0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



	  rho0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  rho0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



        zero2d(vnx,vnz,vp0,npd);

        zero2d(vnx,vnz,rho0,npd);

        read_file(FN1,vnx,vnz,vp0,rho0,npd);


	  vp=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  vp[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }




	  rho=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  rho[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }

	  s=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		   s[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


/**************************************************************/
	 u0_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u0_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 u1_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u1_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w0_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w0_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w1_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w1_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx0_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx0_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 txx1_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx1_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx0_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx0_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx1_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx1_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz0_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz0_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz1_f=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz1_f[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }



/*********************************************************************/


	 u0_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u0_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 u1_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u1_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w0_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w0_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w1_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w1_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx0_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx0_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 txx1_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx1_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx0_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx0_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx1_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx1_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz0_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz0_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz1_b=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz1_b[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);
/*********************************************************************/

/***************the declaration sterm of gradient**************/

	float **f0;
	float **b0;
	float **f1;
	float **b1;
	float **f2;
	float **b2;
	float **t;
	float **a;
	float **vvv;

	f0=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		f0[i]=(float *)calloc(nz,sizeof(float));
	}

	b0=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		b0[i]=(float *)calloc(nz,sizeof(float));
	}

	f1=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		f1[i]=(float *)calloc(nz,sizeof(float));
	}

	b1=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		b1[i]=(float *)calloc(nz,sizeof(float));
	}

	f2=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		f2[i]=(float *)calloc(nz,sizeof(float));
	}

	b2=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		b2[i]=(float *)calloc(nz,sizeof(float));
	}


	t=(float **)calloc(nx,sizeof(float *));
	for(i=0;i<nx;i++)
	{
		t[i]=(float *)calloc(nz,sizeof(float));
	}

	a=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		a[i]=(float *)calloc(vnz,sizeof(float));
	}

	vvv=(float **)calloc(vnx,sizeof(float *));
	for(i=0;i<vnx;i++)
	{
		vvv[i]=(float *)calloc(vnz,sizeof(float));
	}

        dt=dt/1000;

	zero(f0,nx,nz);
	zero(b0,nx,nz);
	zero(f1,nx,nz);
	zero(b1,nx,nz);
	zero(f2,nx,nz);
	zero(b2,nx,nz);
	zero(t,nx,nz);
	zero(a,vnx,vnz);
	zero(vvv,vnx,vnz);
	zero(s,(nz+2*npd),(nx+2*npd));
	read_velfile(vnx,vnz,FN1,vvv);
/*************************************************************/

        float *coffx1;float *coffx2;float *coffz1;float *coffz2;float *acoffx1;float *acoffx2;float *acoffz1;float *acoffz2;
        coffx1=(float *)calloc((nx+2*npd),sizeof(float));
        coffx2=(float *)calloc((nx+2*npd),sizeof(float));
	coffz1=(float *)calloc((nz+2*npd),sizeof(float));
        coffz2=(float *)calloc((nz+2*npd),sizeof(float));
	acoffx1=(float *)calloc((nx+2*npd),sizeof(float));
	acoffx2=(float *)calloc((nx+2*npd),sizeof(float));
	acoffz1=(float *)calloc((nz+2*npd),sizeof(float));
	acoffz2=(float *)calloc((nz+2*npd),sizeof(float));

	for(i=0;i<nx+2*npd;i++)
	{
		coffx1[i]=0.0;
		coffx2[i]=0.0;
		acoffx1[i]=0.0;
		acoffx2[i]=0.0;
	}

	for(i=0;i<nz+2*npd;i++)
	{
		coffz1[i]=0.0;
		coffz2[i]=0.0;
		acoffz1[i]=0.0;
		acoffz2[i]=0.0;
	}

        initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/*************************************************************/


          ndtt=(int)ndtt;

//      for(is=1;is<=ns_sxd;is++)
        //printf("back propagation IS===%d \n",is);

          zero2d2(nx,nz,npd,tx0_f);
	  zero2d2(nx,nz,npd,tx1_f);
	  zero2d2(nx,nz,npd,tz0_f);
	  zero2d2(nx,nz,npd,tz1_f);
          zero2d2(nx,nz,npd,tx0_b);
	  zero2d2(nx,nz,npd,tx1_b);
	  zero2d2(nx,nz,npd,tz0_b);
	  zero2d2(nx,nz,npd,tz1_b);

	  zero3d(nx,nz,nt,npd,u0_f,w0_f,txx0_f,u1_f,w1_f,txx1_f,vp,rho);
	  zero3d(nx,nz,nt,npd,u0_b,w0_b,txx0_b,u1_b,w1_b,txx1_b,vp,rho);

          current_shot(vp0,rho0,vp,rho,nx,nz,npd,vnx,vnz,ds_sxd,is);

          pad_vv(nx,nz,npd,vp);
          pad_vv(nx,nz,npd,rho);


	  for(i=0;i<=nz+2*npd-1;i++)
	  {
		for(j=0;j<=nx+2*npd-1;j++)
		{
		  vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
		  rho[i][j]=1.0/rho[i][j];
		}
	  }



/*************************************************/
        for(it=nt-1;it>=0;it--)
	{


/*******************input the forward P********************/
/********************and forward to back*******************/
	    	for(j=npd;j<npd+nx;j++)
	        {
		    	txx1_f[npd][j]=p_t[j-npd][it];
		}
	    	for(j=npd;j<npd+nx;j++)
		{
		    	txx1_f[npd+nz-1][j]=p_b[j-npd][it];
		}
	    	for(i=npd;i<npd+nz;i++)
		{
		    	txx1_f[i][npd]=p_l[i-npd][it];
		}
	    	for(i=npd;i<npd+nz;i++)
		{
		    	txx1_f[i][npd+nx-1]=p_r[i-npd][it];
		}


              for_update_vel(nx,nz,npd,mm,dt,dx,dz,u0_f,w0_f,txx0_f,u1_f,w1_f,txx1_f,rho,c,coffx1,coffx2,coffz1,coffz2);
              for_update_txx(nx,nz,dt,dx,dz,mm,u0_f,w0_f,txx0_f,u1_f,w1_f,txx1_f,vp,c,npd,
                             tx1_f,tx0_f,tz1_f,tz0_f,acoffx1,acoffx2,acoffz1,acoffz2);
            //abs_bc(u0_f,w0_f,txx0_f,nx,nz,npd,iabsorb);

		for(i=npd;i<npd+nz;i++)
		{
	    	      for(j=npd;j<nx+npd;j++)
			{
			      f0[j-npd][i-npd]=txx0_f[i][j];
			}
		}

//if(myid==0)  printf("the forward model!!!!!!!!!!!!!! \n");


	        for(i=0;i<=nz+2*npd-1;i++)
		{
			for(j=0;j<=nx+2*npd-1;j++)
			{
				u1_f[i][j]=u0_f[i][j];
				w1_f[i][j]=w0_f[i][j];
			    	tx1_f[i][j]=tx0_f[i][j];
			    	tz1_f[i][j]=tz0_f[i][j];
				txx1_f[i][j]=txx0_f[i][j];
			}
		}



/*******************input the residual P***************************/
/********************and back to back******************************/

	   	for(j=npd;j<npd+nx;j++)
	        {
		    s[npd+hsx-1][j]=-p_com[j-npd][it];
		}


            update_vel(nx,nz,npd,mm,dt,dx,dz,u0_b,w0_b,txx0_b,u1_b,w1_b,txx1_b,rho,c,coffx1,coffx2,coffz1,coffz2);
            update_txx(nx,nz,dt,dx,dz,mm,u0_b,w0_b,txx0_b,u1_b,w1_b,txx1_b,s,vp,c,npd,
                              tx1_b,tx0_b,tz1_b,tz0_b,acoffx1,acoffx2,acoffz1,acoffz2);
            //abs_bc(u1_b,w1_b,txx1_b,nx,nz,npd,iabsorb);



		for(i=npd;i<npd+nz;i++)
		{
	    	      for(j=npd;j<nx+npd;j++)
			{
			     b0[j-npd][i-npd]=txx1_b[i][j];
			}
		}

	     //if(myid==0)  printf("the back forward model!!!!!!!!!!!!!! \n");

	        for(i=0;i<=nz+2*npd-1;i++)
		{
			for(j=0;j<=nx+2*npd-1;j++)
			{
				u0_b[i][j]=u1_b[i][j];
				w0_b[i][j]=w1_b[i][j];
			    	tx0_b[i][j]=tx1_b[i][j];
			    	tz0_b[i][j]=tz1_b[i][j];
				txx0_b[i][j]=txx1_b[i][j];
			}
		}
/*****************begin the gradient*********************/


            if(it<=nt-2)
            {
                  if(it<=nt-3)
                  {
                        for(i=1;i<nx-1;i++)
                        {
                           for(j=1;j<nz-1;j++)
                           {
                              t[i][j]=t[i][j]+(f2[i][j]-2*f1[i][j]+f0[i][j])*(b0[i][j])/(dt*dt);
                            }
                        }
                  }

                  for(i=0;i<nx;i++)
                  {
                        for(j=0;j<nz;j++)
                        {
                             f2[i][j]=f1[i][j];
                        }
                  }
            }

            for(i=0;i<nx;i++)
            {
                  for(j=0;j<nz;j++)
                  {
                        f1[i][j]=f0[i][j];
                  }
            }

	}

/**************************the time loop  over***********************/

/***************process the gradient****************/
/*
	for(i=0;i<nx;i++)
	{
		for(j=0;j<nz;j++)
		{
			a[i+(is-1)*ds_sxd][j]=t[i][j];
		}
	}

*/
	for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{
			g[i][j]=2*t[i][j]/(pow(vvv[i][j],3));
		}
	}
/**************************************************/


	  freea(s,(nz+2*npd));
	  freea(vp0,(vnx+2*npd));
	  freea(rho0,(vnx+2*npd));
	  freea(vp,(nz+2*npd));
	  freea(rho,(nz+2*npd));
	  freea(u0_f,(nz+2*npd));
	  freea(w0_f,(nz+2*npd));
	  freea(txx0_f,(nz+2*npd));
	  freea(u1_f,(nz+2*npd));
	  freea(w1_f,(nz+2*npd));
	  freea(txx1_f,(nz+2*npd));
	  freea(tx0_f,(nz+2*npd));
	  freea(tx1_f,(nz+2*npd));
	  freea(tz0_f,(nz+2*npd));
	  freea(tz1_f,(nz+2*npd));
	  freea(u0_b,(nz+2*npd));
	  freea(w0_b,(nz+2*npd));
	  freea(txx0_b,(nz+2*npd));
	  freea(u1_b,(nz+2*npd));
	  freea(w1_b,(nz+2*npd));
	  freea(txx1_b,(nz+2*npd));
	  freea(tx0_b,(nz+2*npd));
	  freea(tx1_b,(nz+2*npd));
	  freea(tz0_b,(nz+2*npd));
	  freea(tz1_b,(nz+2*npd));
	  freea(f0,nx);
	  freea(b0,nx);
	  freea(f1,nx);
	  freea(b1,nx);
	  freea(f2,nx);
	  freea(b2,nx);
	  freea(t,nx);
	  freea(a,vnx);
	  freea(vvv,vnx);
	  free(coffx1);free(coffx2);free(coffz1);free(coffz2);
	  free(acoffx1);free(acoffx2);free(acoffz1);free(acoffz2);




}
/***************main end********************/



/*******************************************************************************/

float step(int nx,int nz,int vnx,int vnz,int nt,int is,int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,int myid,
              int numprocs,float factor[],char FN1[],char FN2[],char FN4[],char FN5[],int npd,float dx,float dz,
              float vdx,float vdz,float favg,float tmax,float dt,float dtout,float pfac,int wavelet,int mz,int vmz)
{

        void model_step(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,
                      float favg,float tmax,float dt,float dtout,float pfac,int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
                       int is,float **p_cal,int wavelet,int vflag);
        void step_zero(int nx,int nz,float **vv0);
        float square_sum(int nx,int nt,float **a);
        void freea(float **a,int m);
        void window(int vnx,int vnz,int m,int n,float **g,float **g1);
	FILE *fp1;
	FILE *fp2;
	FILE *fp3;
	FILE *fp4;
	FILE *fp5;
      FILE *fptemp;


	float a,b,E0,E00,E000,E1,g0,vmin,gmax;
	int i,j,k,m,n;
	float a0[2],E11,E111;
        float L[2],ratio;
	a=0.0;
	ratio=0.0;
	a0[0]=0.0;
	a0[1]=0.0;
	L[0]=0.0;
	L[1]=0.0;
	g0=0.0;
        E0=0.0;
	E00=0.0;
	E000=0.0;
        E1=0.0;
        E11=0.0;
        E111=0.0;
        b=0.0;
        n=0;
/*******************************/
        for(i=0;i<2;i++)
        {
          a0[i]=0.0;
          L[i]=0.0;
        }

/******************************/
	float **v0;
	float **g;
	float **g1;
	float **v1;
	float **p_cal;
	float **p_real;
	float **p_initial;
      float **p_direct;


	fp1=fopen(FN1,"rb");
	fp2=fopen(FN2,"rb");

	fp4=fopen(FN4,"rb");
	fp5=fopen(FN5,"rb");

	v0=(float **)calloc(vnx,sizeof(float *));
	{
		for(i=0;i<vnx;i++)
		{
			v0[i]=(float *)calloc(vnz,sizeof(float));
		}
	}

	v1=(float **)calloc(vnx,sizeof(float *));
	{
		for(i=0;i<vnx;i++)
		{
			v1[i]=(float *)calloc(vnz,sizeof(float));
		}
	}

	g=(float **)calloc(vnx,sizeof(float *));
	{
		for(i=0;i<vnx;i++)
		{
			g[i]=(float *)calloc(vnz,sizeof(float));
		}
	}

	g1=(float **)calloc(vnx,sizeof(float *));
	{
		for(i=0;i<vnx;i++)
		{
			g1[i]=(float *)calloc(vnz,sizeof(float));
		}
	}


	p_cal=(float **)calloc(nx,sizeof(float *));
	{
		for(i=0;i<nx;i++)
		{
			p_cal[i]=(float *)calloc(nt,sizeof(float));
		}
	}


	p_real=(float **)calloc(nx,sizeof(float *));
	{
		for(i=0;i<nx;i++)
		{
			p_real[i]=(float *)calloc(nt,sizeof(float));
		}
	}


	p_initial=(float **)calloc(nx,sizeof(float *));
	{
		for(i=0;i<nx;i++)
		{
			p_initial[i]=(float *)calloc(nt,sizeof(float));
		}
	}
	p_direct=(float **)calloc(nx,sizeof(float *));
	{
		for(i=0;i<nx;i++)
		{
			p_direct[i]=(float *)calloc(nt,sizeof(float));
		}
	}


	step_zero(vnx,vnz,v0);
        step_zero(vnx,vnz,v1);
	step_zero(vnx,vnz,g);
	step_zero(vnx,vnz,g1);
        step_zero(nx,nt,p_cal);
	step_zero(nx,nt,p_real);
	step_zero(nx,nt,p_initial);
        step_zero(nx,nt,p_direct);



        for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{
			fread(&v0[i][j],4,1,fp1);
			fread(&g[i][j],4,1,fp2);

		}
	}


        MPI_Barrier(MPI_COMM_WORLD);

        fclose(fp1);
	  fclose(fp2);

        vmin=0.0;
        gmax=0.0;
        vmin=v0[0][0];
	  window(vnx,vnz,20,vnz,g,g1);

	for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{
			if(vmin>v0[i][j])
			{
			   vmin=v0[i][j];
			}
		}
	}

	for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{
			if(gmax<fabs(g1[i][j]))
			{
			   gmax=fabs(g1[i][j]);
		      }
		}
	}

      step_zero(vnx,vnz,v1);


      MPI_Barrier(MPI_COMM_WORLD);



/*******************calculate the step ratio****************************/

   for(m=1;m<=2;m++)
   {
      if(myid==0)
      {
         printf("---   the m===%d \n",m);
      }

      E00=0.0;
      E000=0.0;
      E11=0.0;
      E111=0.0;
      a0[m-1]=0.0;


	a0[m-1]=2*m*vmin/(1000*gmax);


	for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{

		  v1[i][j]=v0[i][j]-a0[m-1]*g1[i][j];

		}

	}

	if(myid==0)
	{
	     fp3=fopen("vel_for_step.dat","wb");
        fptemp=fopen("Temp_vel_for_step.dat","wb");
	     for(i=0;i<vnx;i++)
	     {
		    for(j=0;j<vnz;j++)
		    {
			fwrite(&v1[i][j],4,1,fp3);
		    }
	     }
           for(i=0;i<vnx;i++)
	     {
		   for(j=0;j<vmz;j++)
		    {
			fwrite(&v1[i][j],4,1,fptemp);
		    }
	     }

	     fclose(fp3);
           fclose(fptemp);
	}


      MPI_Barrier(MPI_COMM_WORLD);



    //////////////////////////////////////////////////
    for(is=1+myid;is<=ns_sxd;is=is+numprocs)
    {
        MPI_Barrier(MPI_COMM_WORLD);

        E0=0.0;
	E1=0.0;
	step_zero(nx,nt,p_initial);
	step_zero(nx,nt,p_real);
	step_zero(nx,nt,p_cal);

        fseek(fp4,(is-1)*nx*nt*4L,0);
        fseek(fp5,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			fread(&p_real[i][j],4,1,fp4);
			fread(&p_initial[i][j],4,1,fp5);
		}
	}


        for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			p_initial[i][j]=p_initial[i][j]-p_real[i][j];
		}
	}
	E0=square_sum(nx,nt,p_initial);
	E00=E00+E0;




        model_step(nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,
                       ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_cal,wavelet,1);
        model_step(nx,mz,vnx,vmz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,
                       ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_direct,wavelet,2);



        for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			p_cal[i][j]=p_cal[i][j]-p_real[i][j]-p_direct[i][j];

		}
	}



	E1=square_sum(nx,nt,p_cal);

	E11=E11+E1;



    }
    //////////////////////////////////////////////////

      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Reduce(&E00,&E000,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);
      MPI_Reduce(&E11,&E111,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);

      if(myid==0)
      {
        E000=E000/2;
        E111=E111/2;
        L[m-1]=(E111-E000)/(a0[m-1]*a0[m-1]);

        printf("---   E000========%f  of m==%d  \n",E000,m);
        printf("---   E111========%f  of m==%d  \n",E111,m);

      }

   }

    if(myid==0)
    {
        ratio=(L[1]-L[0])*a0[0]*a0[1]/(a0[1]-a0[0]);

        g0=square_sum(vnx,vnz,g1);

        ratio=ratio/g0;

        printf("---   ratio====%f \n",ratio);

        a=ratio*g0*a0[1]*a0[1]/(2*(E111-E000+a0[1]*ratio*g0));

        FILE *fpaaa;
      fpaaa=fopen("new_A.txt","a");
      fprintf(fpaaa,"%f  \n",a);
      fclose(fpaaa);
    }

/******************* the step ratio over****************************/


    MPI_Bcast(&a,1,MPI_FLOAT,0,MPI_COMM_WORLD);

///////////////////////////////////////////////////////

    E11=0.0;
    step_zero(vnx,vnz,v1);
    step_zero(vnx,vnz,v0);

    if(myid==0)
    {
         FILE *fp11;
         fp11=fopen(FN1,"rb");
         for(i=0;i<vnx;i++)
         {
		     for(j=0;j<vnz;j++)
		     {
			     fread(&v0[i][j],4,1,fp11);
		     }
         }
         fclose(fp11);


         for(i=0;i<vnx;i++)
         {
		     for(j=0;j<vnz;j++)
		     {

		          v1[i][j]=v0[i][j]-a*g1[i][j];;

		     }

         }

         FILE *fp33;
         FILE *fptemp3;
         fp33=fopen("vel_for_step.dat","wb");
      fptemp3=fopen("Temp_vel_for_step.dat","wb");
         for(i=0;i<vnx;i++)
         {
		     for(j=0;j<vnz;j++)
		     {
		          fwrite(&v1[i][j],4,1,fp33);
		     }
         }
         for(i=0;i<vnx;i++)
         {
		     for(j=0;j<vmz;j++)
		     {
		          fwrite(&v1[i][j],4,1,fptemp3);
		     }
         }
         fclose(fp33);
         fclose(fptemp3);
    }


    MPI_Barrier(MPI_COMM_WORLD);


    for(is=1+myid;is<=ns_sxd;is=is+numprocs)
    {

	E1=0.0;
	step_zero(nx,nt,p_real);
	step_zero(nx,nt,p_cal);

        fseek(fp4,(is-1)*nx*nt*4L,0);
	for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			fread(&p_real[i][j],4,1,fp4);
		}
	}

        model_step(nx,nz,vnx,vnz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,
                    ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_cal,wavelet,1);
        model_step(nx,mz,vnx,vmz,nt,npd,dx,dz,vdx,vdz,favg,tmax,dt,dtout,pfac,
                    ns_sxd,ds_sxd,fs_sxd,zs_sxd,is,p_direct,wavelet,2);


        for(i=0;i<nx;i++)
	{
		for(j=0;j<nt;j++)
		{
			p_cal[i][j]=p_cal[i][j]-p_real[i][j]-p_direct[i][j];

		}
	}

	E1=square_sum(nx,nt,p_cal);

	E11=E11+E1;

    }

    MPI_Reduce(&E11,&E111,1,MPI_FLOAT,MPI_SUM,0,MPI_COMM_WORLD);

///////////////////////////////////////////////////////////

    if(myid==0)
    {
      E111=E111/2;
      factor[0]=100*(E000-E111)/E000;
      printf("---   new E111=======%f a=========%f  \n",E111,a);

      FILE *fp555;
      fp555=fopen("new_E_m.txt","a");
      fprintf(fp555,"%f \n",E111);
      fclose(fp555);


    }

    MPI_Bcast(factor,1,MPI_FLOAT,0,MPI_COMM_WORLD);



	fclose(fp4);
        fclose(fp5);

	freea(v0,vnx);
	freea(g,vnx);
	freea(g1,vnx);
	freea(v1,vnx);
	freea(p_cal,nx);
	freea(p_real,nx);
	freea(p_initial,nx);

	MPI_Barrier(MPI_COMM_WORLD);
	return a;

}


/*******************************/
float square_sum(int nx,int nz,float **a)
{
	int i,j;
	float s=0.0;
	for(i=0;i<nx;i++)
	{
		for(j=0;j<nz;j++)
		{
			s=s+a[i][j]*a[i][j];
		}
	}
	return s;
}

/*******************************/
void step_zero(int nx,int nz,float **vv0)
{
	int i,j;

	for(i=0;i<nx;i++)
	{
	    for(j=0;j<nz;j++)
	    {
	      vv0[i][j]=0.0;
	    }
	}
}
/**********************************/
void window(int vnx,int vnz,int m,int n,float **g,float **g1)
{
      int i,j;

	for(i=0;i<vnx;i++)
	{
		for(j=0;j<vnz;j++)
		{

		  g1[i][j]=0.0;

		}

	}

	for(i=0;i<vnx;i++)
	{
		for(j=m;j<n;j++)
		{

		  g1[i][j]=g[i][j];

		}

	}
}

/************************************/







/*******************************************************************************/

void model_step(int nx,int nz,int vnx,int vnz,int nt,int npd,float dx,float dz,float vdx,float vdz,float favg,
             float tmax,float dt,float dtout,float pfac,int ns_sxd,int ds_sxd,int fs_sxd,int zs_sxd,
              int is,float **p_cal,int wavelet,int vflag)

{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,int wtype,
                     float pi,int npd,int is,int ds_sxd);
      void ptsource_wavelet(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                       int wtype,float pi,int npd,int is,int ds_sxd,int nt,int it);//++++++++++++++
      void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,float **u0,float **w0,float **txx0,
                  float **u1,float **w1,float **txx1,float **rho,float c[],
                 float *coffx1,float *coffx2,float *coffz1,float *coffz2);
      void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                   float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,
                   float **tx1,float **tx0,float **tz1,float **tz0,
                   float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2);
      void abs_bc(float **u1,float **w1,float **txx1,int nx,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,float favg,
                            float dtout,float dt,float **vp0,float ndtt);
      void zero(float **a,int nx,int nz);
      void zero2d(int nx,int nz,float **vv0,int npd);
	void zero2d2(int nx,int nz,int npd,float **vv);
	void zero2d3(int nx,int nz,int npd,float **vv);
	void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,
                     float **u1,float **w1,float **txx1,float **vp,float **rho);
      void pad_vv(int nx,int nz,int npd,float **ee);
	void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd);
	void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,int npd,
                         int vnx,int vnz,int ds_sxd,int is);

      void initial_coffe(float dt,float d0,int nx,int nz,
                       float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                       float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd);
      void freea(float **a,int m);



	  int mm=4;
          int hsx=1;
	  int i,j;

	  int ntout,wtype,it,ifx,ilx,jfz,jlz;
	  float pi,t,ndtt,d0;
	  float iabsorb[4];

	  FILE *fp1;FILE *fp2;
	  char FN1[250];

          if(vflag==1)
              strcpy(FN1,"vel_for_step.dat");
          if(vflag==2)
              strcpy(FN1,"Temp_vel_for_step.dat");

          wtype=1;
          iabsorb[0]=1.0;
          iabsorb[1]=1.0;
          iabsorb[2]=1.0;
          iabsorb[3]=1.0;

	  pi=3.14;
	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;

	  float **vp0;
	  float **rho0;
	  float **vp;
	  float **rho;
	  float **u0;
	  float **w0;
	  float **txx0;
	  float **u1;
	  float **w1;
	  float **txx1;
	  float **tx0;
	  float **tx1;
	  float **tz0;
	  float **tz1;
	  float **s;

	  float c[4];


          cal_c(mm,c);


	  vp0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  vp0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



	  rho0=(float **)calloc((vnx+2*npd),sizeof(float *));
	  {
		  for(i=0;i<(vnx+2*npd);i++)
		  {
		  rho0[i]=(float *)calloc((vnz+2*npd),sizeof(float));
		  }
	  }



          zero2d(vnx,vnz,vp0,npd);
          zero2d(vnx,vnz,rho0,npd);
          read_file(FN1,vnx,vnz,vp0,rho0,npd);


	  vp=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  vp[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


	  rho=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  rho[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


	 u0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 u1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     u1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 w1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     w1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }


	 txx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 txx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     txx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tx1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tx1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz0=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz0[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }

	 tz1=(float **)calloc((nz+2*npd),sizeof(float *));
	 {
	   for(i=0;i<nz+2*npd;i++)
	   {
	     tz1[i]=(float *)calloc((nx+2*npd),sizeof(float));
	   }
	 }




	  s=(float **)calloc((nz+2*npd),sizeof(float *));
	  {
		  for(i=0;i<nz+2*npd;i++)
		  {
		  s[i]=(float *)calloc((nx+2*npd),sizeof(float));
		  }
	  }


          d0=get_constant(dx,dz,nx,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);
          dt=dt/1000;
/***********************************************************/

          float *coffx1;float *coffx2;float *coffz1;float *coffz2;float *acoffx1;float *acoffx2;float *acoffz1;float *acoffz2;
          coffx1=(float *)calloc((nx+2*npd),sizeof(float));
          coffx2=(float *)calloc((nx+2*npd),sizeof(float));
	  coffz1=(float *)calloc((nz+2*npd),sizeof(float));
          coffz2=(float *)calloc((nz+2*npd),sizeof(float));
	  acoffx1=(float *)calloc((nx+2*npd),sizeof(float));
	  acoffx2=(float *)calloc((nx+2*npd),sizeof(float));
	  acoffz1=(float *)calloc((nz+2*npd),sizeof(float));
	  acoffz2=(float *)calloc((nz+2*npd),sizeof(float));

	  for(i=0;i<nx+2*npd;i++)
	  {
		 coffx1[i]=0.0;
		 coffx2[i]=0.0;
		 acoffx1[i]=0.0;
		 acoffx2[i]=0.0;
	  }

	  for(i=0;i<nz+2*npd;i++)
	  {
		 coffz1[i]=0.0;
		 coffz2[i]=0.0;
		 acoffz1[i]=0.0;
		 acoffz2[i]=0.0;
	  }

          initial_coffe(dt,d0,nx,nz,coffx1,coffx2,coffz1,coffz2,acoffx1,acoffx2,acoffz1,acoffz2,npd);

/***********************************************************/
	   ndtt=(int)ndtt;



       //  printf("IS============%d \n",is);
/***********************************************/
           zero(p_cal,nx,nt);
           zero2d2(nx,nz,npd,tx0);
	   zero2d2(nx,nz,npd,tx1);
	   zero2d2(nx,nz,npd,tz0);
	   zero2d2(nx,nz,npd,tz1);
	   zero3d(nx,nz,nt,npd,u0,w0,txx0,u1,w1,txx1,vp,rho);

           current_shot(vp0,rho0,vp,rho,nx,nz,npd,vnx,vnz,ds_sxd,is);

           pad_vv(nx,nz,npd,vp);

           pad_vv(nx,nz,npd,rho);


	   for(i=0;i<=nz+2*npd-1;i++)
	   {
		for(j=0;j<=nx+2*npd-1;j++)
		{
		   vp[i][j]=rho[i][j]*(vp[i][j]*vp[i][j]);
		   rho[i][j]=1.0/rho[i][j];
		}
	   }

	   t=0.0;

    for(it=0;it<nt;it++)
    {

	  t=t+dt;
         if(wavelet==0)
	    ptsource(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd);
         if(wavelet==1)
          ptsource_wavelet(pfac,fs_sxd,zs_sxd,nx,nz,dt,t,favg,s,wtype,pi,npd,is,ds_sxd,nt,it);//+++++++++++++++++++++++++++
          update_vel(nx,nz,npd,mm,dt,dx,dz,u0,w0,txx0,u1,w1,txx1,rho,c,coffx1,coffx2,coffz1,coffz2);
          update_txx(nx,nz,dt,dx,dz,mm,u0,w0,txx0,u1,w1,txx1,s,vp,c,npd,tx1,tx0,tz1,tz0,acoffx1,acoffx2,acoffz1,acoffz2);



	  for(j=npd;j<npd+nx;j++)
	  {
		p_cal[j-npd][it]=txx1[npd+hsx-1][j];

	  }


	  for(i=0;i<nz+2*npd;i++)
	  {
		for(j=0;j<nx+2*npd;j++)
		{
			u0[i][j]=u1[i][j];
			w0[i][j]=w1[i][j];
			tx0[i][j]=tx1[i][j];
			tz0[i][j]=tz1[i][j];
			txx0[i][j]=txx1[i][j];
		}
	   }


     }

/***********************************************/

	  freea(vp0,(vnx+2*npd));
	  freea(rho0,(vnx+2*npd));
	  freea(vp,(nz+2*npd));
	  freea(rho,(nz+2*npd));
	  freea(u0,(nz+2*npd));
	  freea(w0,(nz+2*npd));
	  freea(txx0,(nz+2*npd));
	  freea(u1,(nz+2*npd));
	  freea(w1,(nz+2*npd));
	  freea(txx1,(nz+2*npd));
	  freea(tx0,(nz+2*npd));
	  freea(tx1,(nz+2*npd));
	  freea(tz0,(nz+2*npd));
	  freea(tz1,(nz+2*npd));
	  freea(s,(nz+2*npd));
	  free(coffx1);free(coffx2);free(coffz1);free(coffz2);
	  free(acoffx1);free(acoffx2);free(acoffz1);free(acoffz2);

}



////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/************************************************************************************/

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void for_update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,float **u0,float **w0,float **txx0,
                  float **u1,float **w1,float **txx1,float **rho,float c[],
                  float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
		 int ii,i,j;
		 float dtxx,dtxz,dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;

		 for(i=mm-1;i<=(2*npd+nz-mm-1);i++)
		 {
			 for(j=mm-1;j<=(2*npd+nx-mm-1);j++)
			 {

			   u0[i][j]=coffx2[j]*u1[i][j]+coffx1[j]*dtx*rho[i][j]*(c[0]*(txx1[i][j+1]-txx1[i][j])+c[1]*(txx1[i][j+2]-txx1[i][j-1])+c[2]*(txx1[i][j+3]-txx1[i][j-2])+c[3]*(txx1[i][j+4]-txx1[i][j-3]));

			 }
		 }



		 for(j=mm;j<=(2*npd+nx-mm-1);j++)
		 {
			 for(i=mm;i<=(2*npd+nz-mm-1);i++)
			 {

			   w0[i][j]=coffz2[i]*w1[i][j]+coffz1[i]*dtz*rho[i][j]*(c[0]*(txx1[i+1][j]-txx1[i][j])+c[1]*(txx1[i+2][j]-txx1[i-1][j])+c[2]*(txx1[i+3][j]-txx1[i-2][j])+c[3]*(txx1[i+4][j]-txx1[i-3][j]));

			 }
		 }


}


void for_update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
                 float **u1,float **w1,float **txx1,float **vp,float c[],int npd,float **tx1,float **tx0,
                 float **tz1,float **tz0,float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2)
{

		 float dux,dwz;
		 int i,j,ii;
		 float dtx,dtz;



		 dtx=dt/dx;
		 dtz=dt/dz;

		 for(i=mm;i<=(2*npd+nz-mm-1);i++)
		 {
			 for(j=mm;j<=(2*npd+nx-mm-1);j++)
			 {

		    	   tx0[i][j]=acoffx2[j]*tx1[i][j]+acoffx1[j]*vp[i][j]*dtx*(c[0]*(u0[i][j]-u0[i][j-1])+c[1]*(u0[i][j+1]-u0[i][j-2])+c[2]*(u0[i][j+2]-u0[i][j-3])+c[3]*(u0[i][j+3]-u0[i][j-4]));

                           tz0[i][j]=acoffz2[i]*tz1[i][j]+acoffz1[i]*vp[i][j]*dtz*(c[0]*(w0[i][j]-w0[i-1][j])+c[1]*(w0[i+1][j]-w0[i-2][j])+c[2]*(w0[i+2][j]-w0[i-3][j])+c[3]*(w0[i+3][j]-w0[i-4][j]));

			   txx0[i][j]=tx0[i][j]+tz0[i][j];


			 }
		 }

}


/*******************************************************/
                                                       //
                                                       //
void cal_c(int mm,float c[])                           //
{                                                      //
                                                       //
                                                       //
                                                       //
                                                       //
	  c[0]=1.196289;                               //
          c[1]=-0.0797526;                             //
          c[2]=0.009570313;                            //
          c[3]=-0.000697544;                           //
	                                                 //
                                                       //
}                                                      //




void ptsource(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,int wtype,
                float pi,int npd,int is,int ds_sxd)
{

            if((t==dt)&&(is==1))
	{

          printf("---   ****  get_source  ****\n");

	}

	    float get_source(float ts,float favg,int wtype);

	    int i,j,ixs,izs,x,z;
	    float tdelay,ts,source,fs;

	    for(i=0;i<=nz+2*npd-1;i++)
	    {
		   for(j=0;j<=nx+2*npd-1;j++)
		   {
			s[i][j]=0.0*dt*favg;
		   }
	    }

	    tdelay=1.0/favg;
            ts=t-2.5*tdelay;
            fs=xsn+(is-1)*ds_sxd;
           if(t<=10*tdelay)
       {
            source=get_source(ts,favg,wtype);

	    ixs = (int)(fs+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;

      /*      for(i=izs-3;i<=izs+3;i++)
	    {
		  for(j=ixs-3;j<=ixs+3;j++)
		  {        */
             for(i=izs;i<=izs;i++)//+++++++++++++++++++change
	    {
		  for(j=ixs;j<=ixs;j++)
		  {
			  x=j-ixs;z=i-izs;
                          s[i][j]=pfac*source*exp(-z*z-x*x);
		  }
	    }
	}

}

/****************************************************************************///+++++++++++++++++(
void ptsource_wavelet(float pfac,float xsn,float zsn,int nx,int nz,float dt,float t,float favg,float **s,
                int wtype,float pi,int npd,int is,int ds_sxd,int nt,int it)

{
            if((t==dt)&&(is==1))
	{

          printf("---   **** read wavelet ****\n");

	}


	      float get_source(float ts,float favg,int wtype);
            float read_wavelet(char FN[],int nt,int it);//+++++++++++++++++++++++
	      int i,j,ixs,izs,x,z;
	      float tdelay,ts,source,fs;
            float wavelet=0.0;//+++++++++++++++++++++++++++

	    for(i=0;i<=nz+2*npd-1;i++)
	    {
		   for(j=0;j<=nx+2*npd-1;j++)
		   {
			s[i][j]=0.0;
		   }
	    }

	    tdelay=1.0/favg;
            ts=t-2.5*tdelay;//++++++++++change
            fs=xsn+(is-1)*ds_sxd;
           if(t<=10*tdelay)//++++++++++change
           {
            //source=get_source(ts,favg,wtype);
            wavelet=read_wavelet("wavelet_filter.dat",nt,it);//++++++++++++++++++++you should change the name to you used
	    ixs = (int)(fs+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;


            /*      for(i=izs-3;i<=izs+3;i++)
	    {
		  for(j=ixs-3;j<=ixs+3;j++)
		  {        */
             for(i=izs;i<=izs;i++)//+++++++++++++++++++change
	    {
		  for(j=ixs;j<=ixs;j++)
		  {
			  x=j-ixs;z=i-izs;
                          s[i][j]=pfac*wavelet*exp(-z*z-x*x);//++++++++++
                          //s[i][j]=pfac*wavelet;//++++++++++

                       /*   FILE *fpwave;
                          fpwave=fopen("wavelet.txt","a");
                          fprintf(fpwave,"%f \n",s[i][j]);
                          fclose(fpwave);   */
		  }
	    }
	   }

}
float read_wavelet(char FN[],int nt,int it)//+++++++++++++++++++++
{
    int i,j;
    float wavelet=0.0;
    FILE *fp1;
    fp1=fopen(FN,"rb");
    fseek(fp1,((it-1)*4L),0);//++++++you can change
    fread(&wavelet,4,1,fp1);
    fclose(fp1);
    return (wavelet);

}
/****************************************************************************///+++++++++++++++++)
/*float get_source(float ts,float favg,int wtype)
{
        int i,j,nt,ntfft;
        int nf1,nf2,nf3,nf4;
	    float x,pi;
	    float source=0.0;
	    float  f1,f2,dw,fw;
        float tmpp=0.0;
        float dt=0.6/1000;
        float tdelay;
        float    *trace;
        complex  *ctrace;
        tdelay=1.0/favg;
        nt=2501;

        f1=5.0;
        f2=6.0;
        pi=3.14;
        dw=2.0*pi/(ntfft*dt);
        fw=2.0*pi*f1;
        nf1=fw/dw+0.5;
        fw=2.0*pi*f2;
        nf2=fw/dw+0.5;

        ntfft=npfa(nt);
        nf3=ntfft-1;
        nf4=ntfft;
        trace=alloc1float(nt);
        ctrace=alloc1complex(ntfft);
        zero1float(trace,nt);
        zero1complex(ctrace,ntfft);

          for(i=0;i<nt;i++)
          {
          if(((i+1)*dt)<=2*tdelay)
          {
           x=(favg*((i+1)*dt-tdelay))*(favg*((i+1)*dt-tdelay));
           trace[i]=(1-2*pi*pi*(x))*exp(-(pi*pi*x));
           }
           else
           trace[i]=0.0;
          }

         for(i=0;i<nt;i++)
         {
           ctrace[i]=cmplx(trace[i],0.0);
         }
         for(i=nt;i<ntfft;i++)
         {
           ctrace[i]=cmplx(0.0, 0.0);
         }
         pfacc(1, ntfft, ctrace);

	    for(i=0;i<ntfft;i++)
	    {
          if(i>=nf1&&i<=nf2)
          {
              tmpp=0.54+0.46*cos(PI*(i-nf1)/(nf2-nf1)-PI);
              ctrace[i]=crmul(ctrace[i],tmpp);
          }
          else if(i>=nf3&&i<=nf4)
          {
              tmpp=0.54+0.46*cos(PI*(nf3-i)/(nf4-nf3));
              ctrace[i]=crmul(ctrace[i],tmpp);
          }
          else if(i<nf1)
              ctrace[i]=cmplx(0.0, 0.0);
	    }

	    pfacc (-1, ntfft, ctrace);
	    for(i=0;i<ntfft;i++)
	    {
	        ctrace[i]=crmul(ctrace[i],1.0/ntfft);
	    }
		    j=((ts+tdelay)/dt)-1;
            source=ctrace[j].r;

	    return (source);


}*/

float get_source(float ts,float favg,int wtype)
{
	    float x,pi;
	    float source=0.0;

	    pi=3.14;
            x=(favg*(ts))*(favg*(ts));

            source=(1-2*pi*pi*(x))*exp(-(pi*pi*x));
	    return (source);


}



void update_vel(int nx,int nz,int npd,int mm,float dt,float dx,float dz,float **u0,float **w0,float **txx0,
                float **u1,float **w1,float **txx1,float **rho,float c[],
                float *coffx1,float *coffx2,float *coffz1,float *coffz2)
{
		 int ii,i,j;
		 float dtxx,dtxz,dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;

		 for(i=mm-1;i<=(2*npd+nz-mm-1);i++)
		 {
			 for(j=mm-1;j<=(2*npd+nx-mm-1);j++)
			 {

			   u1[i][j]=coffx2[j]*u0[i][j]-coffx1[j]*dtx*rho[i][j]*(c[0]*(txx0[i][j+1]-txx0[i][j])+c[1]*(txx0[i][j+2]-txx0[i][j-1])+c[2]*(txx0[i][j+3]-txx0[i][j-2])+c[3]*(txx0[i][j+4]-txx0[i][j-3]));

			 }
		 }



		 for(j=mm;j<=(2*npd+nx-mm-1);j++)
		 {
			 for(i=mm;i<=(2*npd+nz-mm-1);i++)
			 {

			   w1[i][j]=coffz2[i]*w0[i][j]-coffz1[i]*dtz*rho[i][j]*(c[0]*(txx0[i+1][j]-txx0[i][j])+c[1]*(txx0[i+2][j]-txx0[i-1][j])+c[2]*(txx0[i+3][j]-txx0[i-2][j])+c[3]*(txx0[i+4][j]-txx0[i-3][j]));


			 }
		 }


}


void update_txx(int nx,int nz,float dt,float dx,float dz,int mm,float **u0,float **w0,float **txx0,
               float **u1,float **w1,float **txx1,float **s,float **vp,float c[],int npd,
               float **tx1,float **tx0,float **tz1,float **tz0,float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2)
{

		 float dux,dwz;
		 int i,j,ii;
		 float dtx,dtz;

		 dtx=dt/dx;
		 dtz=dt/dz;

		 for(i=mm;i<=(2*npd+nz-mm-1);i++)
		 {
			 for(j=mm;j<=(2*npd+nx-mm-1);j++)
			 {

		    	   tx1[i][j]=acoffx2[j]*tx0[i][j]-acoffx1[j]*vp[i][j]*dtx*(c[0]*(u1[i][j]-u1[i][j-1])+c[1]*(u1[i][j+1]-u1[i][j-2])+c[2]*(u1[i][j+2]-u1[i][j-3])+c[3]*(u1[i][j+3]-u1[i][j-4]));

                     tz1[i][j]=acoffz2[i]*tz0[i][j]-acoffz1[i]*vp[i][j]*dtz*(c[0]*(w1[i][j]-w1[i-1][j])+c[1]*(w1[i+1][j]-w1[i-2][j])+c[2]*(w1[i+2][j]-w1[i-3][j])+c[3]*(w1[i+3][j]-w1[i-4][j]));

			   txx1[i][j]=tx1[i][j]+tz1[i][j]+s[i][j];

			 }
		 }

}




float get_constant(float dx,float dz,int nx,int nz,int nt,int ntout,int npd,float tmax,float favg,
                   float dtout,float dt,float **vp0,float ndtt)
{
		 int i,j;
		 float vpmax,vpmin,H_min;
		 float dt_max,dx_max,dz_max,d0;

		 vpmax=vp0[npd][npd];
		 vpmin=vp0[npd][npd];

		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 if(vpmax<vp0[i][j]) vpmax=vp0[i][j];
				 if(vpmin>vp0[i][j]) vpmin=vp0[i][j];
			 }
		 }


		 d0=3.0*vpmax*log(100000.0)/(2.0*npd*dx);


		 if(dx<dz) H_min=dx;
		 else H_min=dz;


/*====== determine time sampling interval to ensure stability====*/

		 dt_max=0.5*1000*H_min/vpmax;
                 dx_max=vpmin/favg*0.2;
                 dz_max=dx_max;


                if(dx_max<dx)
                {
                   printf("---   dx_max===%f, vpmin===%f, favg===%f \n",dx_max,vpmin,favg);
		   printf("---   YOU NEED HAVE TO REDEFINE DX ! \n");
                   exit(0);
		 }
                 if(dz_max<dz)
		 {
		   printf("---   YOU NEED HAVE TO REDEFINE DZ ! \n");
                   exit(0);
		 }
	         if(dt_max<dt)
		 {
		   printf("---   YOU NEED HAVE TO REDEFINE dt ! \n");
                   exit(0);
		 }

                 return d0;



}


void zero2d(int nx,int nz,float **vv0,int npd)
{
		 int i,j;

		 for(i=0;i<nx+2*npd;i++)
		 {
			 for(j=0;j<nz+2*npd;j++)
			 {
				 vv0[i][j]=0.0;
			 }
		 }

}


void zero2d2(int nx,int nz,int npd,float **vv)
{
		 int i,j;


		 for(i=0;i<=nz+2*npd-1;i++)
		 {
			 for(j=0;j<=nx+2*npd-1;j++)
			 {
				 vv[i][j]=0.0;
			 }
		 }

}



void zero3d(int nx,int nz,int nt,int npd,float **u0,float **w0,float **txx0,float **u1,float **w1,float **txx1,
              float **vp,float **rho)
{
		 int i,j,k;


		   for(i=0;i<nz+2*npd;i++)
		   {
			  for(j=0;j<nx+2*npd;j++)
			  {
				 vp[i][j]=0.0;
				 rho[i][j]=0.0;

				 u0[i][j]=0.0;
                                 w0[i][j]=0.0;
                                 txx0[i][j]=0.0;
				 u1[i][j]=0.0;
                                 w1[i][j]=0.0;
                                 txx1[i][j]=0.0;
			  }
		   }

}

void pad_vv(int nx,int nz,int npd,float **ee)
{
		 int i,j;


/*****pad left side                    */
            for(i=npd;i<=(nz+npd-1);i++)
		{
              for(j=0;j<=npd-1;j++)
              {
               ee[i][j]=ee[i][npd];
              }
		}

/*****pad right side                    */
            for(i=npd;i<=(nz+npd-1);i++)
		{
              for(j=nx+npd;j<=(nx+2*npd-1);j++)
              {
                ee[i][j]=ee[i][nx+npd-1];
              }
		}
/*****pad upper side                    */
            for(i=0;i<=(npd-1);i++)
		{
              for(j=0;j<=(nx+2*npd-1);j++)
              {
                ee[i][j]=ee[npd][j];
              }
		}
/*****lower side                        */
            for(i=nz+npd;i<=(nz+2*npd-1);i++)
		{
              for(j=0;j<=(nx+2*npd-1);j++)
              {
                ee[i][j]=ee[nz+npd-1][j];
              }
		}


}

void read_file(char FN1[],int nx,int nz,float **vv,float **rho0,int npd)
{

		 int i,j;

		 FILE *fp1;
		 fp1=fopen(FN1,"rb");
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 fread(&vv[i][j],4,1,fp1);

			 }
		 }
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(j=npd;j<nz+npd;j++)
			 {
				 rho0[i][j]=1.0;
			 }
		 }
		 fclose(fp1);
}


void current_shot(float **vp0,float **rho0,float **vp,float **rho,int nx,int nz,int npd,int vnx,int vnz,int ds_sxd,int is)
{


                         int ivstart,ivend;
			 int i,ix,iz;
                         is=1;
                         ivstart=1+(is-1)*ds_sxd;
			 ivend=nx+(is-1)*ds_sxd;

			 if(ivstart<=0)
			 {
				 printf("---   ivstart less than zero \n");
				 exit(0);
			 }
			 if(ivend>vnx)
			 {
				 printf("---   ivend great than Vnx \n");
				 exit(0);
			 }

			 for(ix=npd+ivstart-1;ix<ivend+npd;ix++)
			 {
				 for(iz=npd;iz<nz+npd;iz++)
				 {
				  vp[iz][ix-ivstart+1]=vp0[ix][iz];
                                  rho[iz][ix-ivstart+1]=rho0[ix][iz];

				 }
			 }

}


void initial_coffe(float dt,float d0,int nx,int nz,float *coffx1,float *coffx2,float *coffz1,float *coffz2,
                   float *acoffx1,float *acoffx2,float *acoffz1,float *acoffz2,int npd)
{
		 int i,j;


		 for(i=0;i<npd;i++)
		 {
			 coffx1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz2[i]=coffz1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);


		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 coffx1[i]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
		 }
		 for(i=npd+nz;i<nz+2*npd;i++)
		 {
			 coffz1[i]=1/(1+(dt*d0*pow((0.5+i-nz-npd)/npd,2))/2);
			 coffz2[i]=coffz1[i]*(1-(dt*d0*pow((0.5+i-nz-npd)/npd,2))/2);
		 }

		 for(i=npd;i<npd+nx;i++)
		 {
			 coffx1[i]=1.0;
			 coffx2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 coffz1[i]=1.0;
			 coffz2[i]=1.0;
		 }


		 for(i=0;i<npd;i++)
		 {
			 acoffx1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffx2[i]=coffx1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz2[i]=coffz1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);

		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 acoffx1[i]=1/(1+(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
			 acoffx2[i]=coffx1[i]*(1-(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
		 }
		 for(i=npd+nz;i<nz+2*npd;i++)
		 {
			 acoffz1[i]=1/(1+(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
			 acoffz2[i]=coffz1[i]*(1-(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
		 }

		 for(i=npd;i<npd+nx;i++)
		 {
			 acoffx1[i]=1.0;
			 acoffx2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 acoffz1[i]=1.0;
			 acoffz2[i]=1.0;
		 }

}



void freea(float **a,int m)
{
	int i;

	for(i=0;i<m;i++)
	{
		free(a[i]);
	}
	free(a);
}


