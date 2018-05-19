#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "mpi.h"
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
main(int argc,char *argv[])
{

        void model(int nx,int ny,int nz,int vnx,int vny,int vnz,int nt,int npd,
                  float dx,float dy,float dz,float vdx,float vdy,float vdz,float favg,float tmax,float dt,
                  float dtout,float pfac,char FN1[],int ns_sxd,int xs_sxd,int ys_sxd,int zs_sxd,int is,
                  float ***p_initial,float ***p_down,float ***p_left,float ***p_right,float ***p_front,
                  float ***p_back,float ***ws0,int myid);
        void shot_compare(int nx,int ny,int nt,int is,char FN1[],float ***p_obs,float ***p_initial,float ***p_com);
        void back_grad(int myid,int nx,int ny,int nz,int vnx,int vny,int vnz,int nt,int npd,
                      float dx,float dy,float dz,float vdx,float vdy,float vdz,float favg,float tmax,float dt,
                      float dtout,char FN1[],int ns_sxd,int xs_sxd,int ys_sxd,int zs_sxd,int is,float ***p_com,
                      float ***p_up,float ***p_down,float ***p_left,float ***p_right,float ***p_front,
                      float ***p_back,float ***g,float ***g0_mig,float ***g0_tom);
	  void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd);  
        float square_sum(int nx,int ny,int nz,float ***a); 
        void read_velfile(int vnx,int vny,int vnz,char FN1[],float ***v);  
        void window(int vnx,int vny,int vnz,int m,int n,float ***g,float ***g1);   

       
              
/*******************************************************/
        char FN1[250]={"shot707070725.dat"};//real shot
	  char FN2[250]={"vel707070initial.dat"};//initial OR uodate velocity
        char FN3[250]={"shot_initial.dat"};//shot of initial velocity
        char FN4[250]={"multi_gradient_zm.dat"};//multi_gradient
        char FN8[250]={"multi_gradient_mig_zm.dat"};//multi_gradient-migration
        char FN9[250]={"multi_gradient_tom_zm.dat"};//multi_gradient-tomography
        char FN5[250]={"cg_gradient.dat"};//conjuage gradient

        char FN6[250]={"Vel_Recording1.dat"};//recording iteration velocity

        char FN7[250]={"vel707070.dat"};//real velocity for test the gradient "back_grad"

        char FN10[250]={"multi_grad_no_zm.dat"};//multi_gradient
        char FN11[250]={"multi_grad_no_zm_mig.dat"};//multi_gradient-migration
        char FN12[250]={"multi_grad_no_zm_tom.dat"};//multi_gradient-tomography         



/*******************************************************/     
        
	int i,j,k,is,nx,ny,nz,nt,vnx,vny,vnz,i_start,i_end,l,ipoynting;
	int ns_sxd,zs_sxd,npd;
	int xs_sxd[100],ys_sxd[100];
	float dx,dy,dz,vdx,vdy,vdz,tmax,dt,dtout,pfac,favg;
	float wsmax;
        float a,alph,beta;float factor[1];
	
	int myid,numprocs,count;

/*******************************************************/
      ipoynting=1; //1->poynting ; other->usual
                         


        nx=70;          npd=30;   tmax=3.0;
        ny=70;          favg=60;  pfac=10000.0;
	  nz=70;         
 	
        vnx=70;         dx=5.0;   
	  vny=70;         dy=5.0;   
	  vnz=70;         dz=5.0; 
                         vdx=5.0;
                         vdy=5.0;
	  nt=725;         vdz=5.0;
        dt=0.4;
        dtout=0.4;
         
        ns_sxd=28;
    /*    xs_sxd[0]=10;xs_sxd[1]=20;xs_sxd[2]=30;
        xs_sxd[3]=40;xs_sxd[4]=50;xs_sxd[5]=60;xs_sxd[6]=70;
        ys_sxd[0]=10;ys_sxd[1]=20;ys_sxd[2]=30;
        ys_sxd[3]=40;ys_sxd[4]=50;ys_sxd[5]=60;ys_sxd[6]=70;  */
   /*      for(i=0;i<8;i++)
           for(j=0;j<8;j++)
               {
             xs_sxd[i*8+j]=(i+1)*10;
               ys_sxd[i*8+j]=(j+1)*10+i;}
          zs_sxd=1; */

         /*for(i=0;i<8;i++)
           for(j=0;j<8;j++)
               {
             xs_sxd[i*8+j]=(i+1)*5;
               ys_sxd[i*8+j]=(j+1)*5+i;}*/
                for(i=0;i<14;i++)
                {  xs_sxd[i]=5+i*5;
                   ys_sxd[i]=1;}
                for(i=14;i<28;i++)
                {  xs_sxd[i]=69;
                   ys_sxd[i]=0+(i-14)*5;}
                zs_sxd=1;  
             
/*******************************************************/        

        float ***g;
        float ***g_mig;
        float ***g_tom;
	  float ***g0;
	  float ***g0_mig;
	  float ***g0_tom;

     	float ***ggg0;
	float ***prp;
	float ***new_g;
	float ***v_initial;
        float ***p_initial;
        float ***p_down;
        float ***p_left;
        float ***p_right;
        float ***p_front;
        float ***p_back;
        float ***p_com;
	float ***p_obs;
	float ***ws;
	float ***ws0;
/*******************************************************/      
	
        p_initial=alloc3float(nt,ny,nx);
        p_down=alloc3float(nt,ny,nx);

        p_left=alloc3float(nt,nz,ny);
        p_right=alloc3float(nt,nz,ny);

        p_front=alloc3float(nt,nz,nx);
        p_back=alloc3float(nt,nz,nx);

        p_com=alloc3float(nt,ny,nx);
        p_obs=alloc3float(nt,ny,nx);

        g=alloc3float(nz,ny,nx);
        g_mig=alloc3float(nz,ny,nx);
        g_tom=alloc3float(nz,ny,nx);
        g0=alloc3float(nz,ny,nx);
        g0_tom=alloc3float(nz,ny,nx);
        g0_mig=alloc3float(nz,ny,nx);
        g0_tom=alloc3float(nz,ny,nx);
        ggg0=alloc3float(nz,ny,nx);
        prp=alloc3float(nz,ny,nx);
        new_g=alloc3float(nz,ny,nx);

        ws=alloc3float(nz,ny,nx);
        ws0=alloc3float(nz,ny,nx);
        v_initial=alloc3float(nz,ny,nx);


/*******************************************/
      MPI_Init(&argc,&argv);
      MPI_Comm_rank(MPI_COMM_WORLD,&myid);
      MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
/******************************************/	


        zero3float(v_initial,nz,ny,nx);
        zero3float(g,nz,ny,nx);
        zero3float(g_mig,nz,ny,nx);
        zero3float(g_tom,nz,ny,nx);
        zero3float(g0,nz,ny,nx);
        zero3float(g0_mig,nz,ny,nx);
        zero3float(g0_tom,nz,ny,nx);
        zero3float(new_g,nz,ny,nx);
        zero3float(prp,nz,ny,nx);
        zero3float(ws,nz,ny,nx);
        zero3float(ws0,nz,ny,nx);
        zero3float(p_initial,nt,ny,nx);

        for(i=0;i<vnx;i++)
          for(j=0;j<vny;j++)
	    for(k=0;k<vnz;k++)
	        ggg0[i][j][k]=1.0;


        MPI_Barrier(MPI_COMM_WORLD);
        zero3float(new_g,nz,ny,nx);

/******************************************************************************/
        count=79;
  
       FILE *fpvelrecording,*fpdebug;
      fpvelrecording=fopen(FN6,"wb");
     fpdebug=fopen("Debug.txt","w");

do
{
  
    if(myid==0)
    {

      
        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("--------------------------------------------------------\n");
        printf("\n");
        printf("\n");
        printf("\n");   
        printf("The count=%d iterative Begin !!!!!!!! \n",count);
        fprintf(fpdebug,"------------------------------------------------\n");
        fprintf(fpdebug,"------------------------------------------------\n");
        fprintf(fpdebug,"---   Count = %d is start .\n",count);
    }
      //goto loop1;
/*************************************/
    FILE *fp3;
    fp3=fopen(FN3,"wb");  
    FILE *fpcomp;
    fpcomp=fopen("shot_compare.dat","wb");   
    FILE *fptgrad,*fptgrad_mig,*fptgrad_tom;
    fptgrad=fopen("tempor_grad.dat","wb");
    fptgrad_mig=fopen("tempor_grad_mig.dat","wb");
    fptgrad_tom=fopen("tempor_grad_tom.dat","wb");
    FILE *fptws;
    fptws=fopen("tempor_ws.dat","wb"); 
/*************************************/ 
    a=0.0;  
   // factor[0]=1.0;


    zero3float(v_initial,nz,ny,nx);
    zero3float(g,nz,ny,nx);
    zero3float(g_mig,nz,ny,nx);
    zero3float(g_tom,nz,ny,nx);
    zero3float(g0,nz,ny,nx);
    zero3float(g0_mig,nz,ny,nx);
    zero3float(g0_tom,nz,ny,nx);
    zero3float(ws,nz,ny,nx);
    zero3float(ws0,nz,ny,nx);
/*************************************/
    MPI_Barrier(MPI_COMM_WORLD);                                                       
/**************begin the loop******************************/
    for(is=1+myid;is<=ns_sxd;is=is+numprocs)	
    {     
        if(myid==0)
        {
             printf("begin ************ IS========%d  \n",is);
        fprintf(fpdebug,"---   begin ************ IS========%d  \n",is);
        }
 
/***************************************/

        zero3float(p_obs,nt,ny,nx);

        zero3float(p_initial,nt,ny,nx);
        zero3float(p_down,nt,ny,nx);
        zero3float(p_left,nt,nz,ny);
        zero3float(p_right,nt,nz,ny);
        zero3float(p_front,nt,nz,nx);
        zero3float(p_back,nt,nz,nx);

        zero3float(p_com,nt,ny,nx);
        zero3float(g0,nz,ny,nx);
        zero3float(g0_mig,nz,ny,nx);
        zero3float(g0_tom,nz,ny,nx);
        zero3float(ws0,nz,ny,nx);

/***************************************/
   
        if(myid==0)
	{
	  printf("the model is start  !    \n");
        fprintf(fpdebug,"---   the model is start  !    \n");
	}
	
        model(nx,ny,nz,vnx,vny,vnz,nt,npd,dx,dy,dz,vdx,vdy,vdz,favg,tmax,dt,dtout,pfac,FN2,ns_sxd,
            xs_sxd[is-1],ys_sxd[is-1],zs_sxd,is,p_initial,p_down,p_left,p_right,p_front,p_back,ws0,myid);        
            
        fseek(fp3,(is-1)*nx*ny*nt*4L,0);
        for(j=0;j<ny;j++)
	  for(i=0;i<nx;i++)
	  { 
	    for(k=0;k<nt;k++)
	    {
	      fwrite(&p_initial[i][j][k],4L,1,fp3);
	    }
	  }
         




        fseek(fptws,(is-1)*nx*ny*nz*4L,0);
        for(j=0;j<ny;j++)
	  for(i=0;i<nx;i++)
	  { 
	    for(k=0;k<nz;k++)
	    {
	      fwrite(&ws0[i][j][k],4L,1,fptws);
	    }
	  }
	
     
	if(myid==0)
	{
	  printf("the model is over  !    \n");
          printf("the compare is start  !    \n");
        fprintf(fpdebug,"---   the model is over  !\n");
        fprintf(fpdebug,"---   the compare is start  !\n");

	}
 	
        MPI_Barrier(MPI_COMM_WORLD);
        shot_compare(nx,ny,nt,is,FN1,p_obs,p_initial,p_com);

        fseek(fpcomp,(is-1)*nx*ny*nt*4L,0);
        for(j=0;j<ny;j++)
	  for(i=0;i<nx;i++)
	  { 
	    for(k=0;k<nt;k++)
	    {
	      fwrite(&p_com[i][j][k],4L,1,fpcomp);
	    }
	  }

        if(myid==0)
	{
          printf("the compare is over  !    \n");
          printf("the back propagation is start  !    \n");
        fprintf(fpdebug,"---   the compare is over  !\n");
        fprintf(fpdebug,"---   the back propagation is start  !\n");
	}

/******************************************************************************************/  
/***********************************  back propagation   **********************************/  
/******************************************************************************************/   
        back_grad(myid,nx,ny,nz,vnx,vny,vnz,nt,npd,dx,dy,dz,vdx,vdy,vdz,favg,tmax,dt,dtout,FN2,ns_sxd,
                 xs_sxd[is-1],ys_sxd[is-1],zs_sxd,is,p_com,p_initial,p_down,p_left,p_right,p_front,p_back,
                 g0,g0_mig,g0_tom);


        if(myid==0)
	{
          printf("the back propagation is over    \n");
        fprintf(fpdebug,"---   the back propagation is over  !\n");
	} 

        fseek(fptgrad,(is-1)*vnx*vny*vnz*4L,0);
        fseek(fptgrad_mig,(is-1)*vnx*vny*vnz*4L,0);
        fseek(fptgrad_tom,(is-1)*vnx*vny*vnz*4L,0);
        for(j=0;j<vny;j++)
	  for(i=0;i<vnx;i++)
	  { 
	    for(k=0;k<vnz;k++)
	    {
	      fwrite(&g0[i][j][k],4L,1,fptgrad);
	      fwrite(&g0_mig[i][j][k],4L,1,fptgrad_mig);
	      fwrite(&g0_tom[i][j][k],4L,1,fptgrad_tom);
	    }
	  }
        
      




    } //***the IS loop ending 

      if(myid==0)
      {
          printf("the loop is done\n");
        fprintf(fpdebug,"---   the IS loop is over  !\n");
      }

      MPI_Barrier(MPI_COMM_WORLD);
  
      fclose(fp3);
      fclose(fpcomp);
      fclose(fptgrad);
      fclose(fptgrad_mig);
      fclose(fptgrad_tom);
      fclose(fptws);

      if(myid==0)
      {         
      
         FILE *fptgrad2,*fptgrad2_mig,*fptgrad2_tom;
         fptgrad2=fopen("tempor_grad.dat","rb");
         fptgrad2_mig=fopen("tempor_grad_mig.dat","rb");
         fptgrad2_tom=fopen("tempor_grad_tom.dat","rb");
         FILE *fptws2;
         fptws2=fopen("tempor_ws.dat","rb");
         FILE *fpmgrad,*fpmgrad_mig,*fpmgrad_tom,*fpmg,*fpmg_mig,*fpmg_tom;
         fpmgrad=fopen(FN4,"wb");
         fpmgrad_mig=fopen(FN8,"wb");
         fpmgrad_tom=fopen(FN9,"wb");
         fpmg=fopen(FN10,"wb");
         fpmg_mig=fopen(FN11,"wb");
         fpmg_tom=fopen(FN12,"wb");
         FILE *fpmws;
         fpmws=fopen("multi_ws.dat","wb");
    
         for(l=0;l<ns_sxd;l++)
         {
            zero3float(g0,nz,ny,nx);
            zero3float(g0_mig,nz,ny,nx);
            zero3float(g0_tom,nz,ny,nx);
            zero3float(ws0,nz,ny,nx);
            fseek(fptgrad2,l*vnx*vny*vnz*4L,0);
            fseek(fptgrad2_mig,l*vnx*vny*vnz*4L,0);
            fseek(fptgrad2_tom,l*vnx*vny*vnz*4L,0);
            for(j=0;j<vny;j++)
            {
	       for(i=0;i<vnx;i++)
	       { 
	         for(k=0;k<vnz;k++)
	         {
	           fread(&g0[i][j][k],4L,1,fptgrad2);
	           fread(&g0_mig[i][j][k],4L,1,fptgrad2_mig);
	           fread(&g0_tom[i][j][k],4L,1,fptgrad2_tom);
	         }
	       }
            }
             
            fseek(fptws2,l*vnx*vny*vnz*4L,0);
            for(j=0;j<vny;j++)
            {
	       for(i=0;i<vnx;i++)
	       { 
	         for(k=0;k<vnz;k++)
	         {
	           fread(&ws0[i][j][k],4L,1,fptws2);
	         }
	       }
            }       
    
            for(j=0;j<vny;j++)
            {
	       for(i=0;i<vnx;i++)
	       { 
	         for(k=0;k<vnz;k++)
	         {
	            g[i][j][k]=g[i][j][k]+g0[i][j][k];
	            g_mig[i][j][k]=g_mig[i][j][k]+g0_mig[i][j][k];
	            g_tom[i][j][k]=g_tom[i][j][k]+g0_tom[i][j][k];
                  ws[i][j][k]=ws[i][j][k]+ws0[i][j][k];
	         }
	       }
            }    
         } //  end of "l" (the ns_sxd loop)   
                   
         wsmax=0.0;
         for(j=0;j<vny;j++)
         {
            for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         if(wsmax<fabs(ws[i][j][k]))
                      wsmax=fabs(ws[i][j][k]);
	      }
	    }
         }
/********************* no ZM's gradient******************/
        for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         fwrite(&g[i][j][k],4L,1,fpmg);
	         fwrite(&g_mig[i][j][k],4L,1,fpmg_mig);
	         fwrite(&g_tom[i][j][k],4L,1,fpmg_tom);
               fwrite(&ws[i][j][k],4L,1,fpmws);
	      }
	    }
         }
/************************* ZM *************************/
         for(j=0;j<vny;j++)
         {
            for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
                  ws[i][j][k]=ws[i][j][k]/wsmax;
                  g[i][j][k]=g[i][j][k]/ws[i][j][k];
                  g_mig[i][j][k]=g_mig[i][j][k]/ws[i][j][k];
                  g_tom[i][j][k]=g_tom[i][j][k]/ws[i][j][k];
	      }
	    }
         }
/********************* ZM's gradient******************/
         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         fwrite(&g[i][j][k],4L,1,fpmgrad);
	         fwrite(&g_mig[i][j][k],4L,1,fpmgrad_mig);
	         fwrite(&g_tom[i][j][k],4L,1,fpmgrad_tom);
               fwrite(&ws[i][j][k],4L,1,fpmws);
	      }
	    }
         }
/********************* ZM's ending ******************/       
         fclose(fptgrad2);
         fclose(fptgrad2_mig);
         fclose(fptgrad2_tom);
         fclose(fptws2);
         fclose(fpmgrad);
         fclose(fpmgrad_mig);
         fclose(fpmgrad_tom);
         fclose(fpmg);
         fclose(fpmg_mig);
         fclose(fpmg_tom);
         fclose(fpmws);

         printf("complete the global gradient  !\n");
        fprintf(fpdebug,"---   complete the global gradient  !\n");

/******************************************************************/
/************************  conjugate process***********************/
/******************************************************************/
         alph=0.0;
         zero3float(prp,nz,ny,nx);

         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         prp[i][j][k]=(g[i][j][k]-ggg0[i][j][k])*g[i][j][k];
                 alph=alph+prp[i][j][k];
	      }
	    }
         }
       
         beta=alph/(square_sum(vnx,vny,vnz,ggg0));
       
         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         new_g[i][j][k]=g[i][j][k]+beta*new_g[i][j][k];
                 ggg0[i][j][k]=g[i][j][k];
	      }
	    }
         }
         printf("complete the cg gradient  !\n");
        fprintf(fpdebug,"---   complete the cg gradient  !\n");

      }//  end of the "if(myid==0)"

      if(myid==0)
      {
         FILE *fpcg;
         fpcg=fopen(FN5,"wb");
         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         fwrite(&new_g[i][j][k],4L,1,fpcg);
	      }
	    }
         }
 
         fclose(fpcg);
      }
      MPI_Barrier(MPI_COMM_WORLD);

/******************************************************************************************/  
/***********************************  get the  step   *************************************/  
/******************************************************************************************/  \
 //loop1:
      if(myid==0)
      {
       
         read_velfile(vnx,vny,vnz,FN2,v_initial);
         zero3float(g0,nz,ny,nx);
         window(vnx,vny,vnz,10,vnz,new_g,g0);
         float vamin,gamax;
         vamin=v_initial[10][10][10];
         gamax=g0[vnx-10][vny-10][vnz-10];
       
         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	          if(vamin>v_initial[i][j][k])
                        vamin=v_initial[i][j][k];
                  if(gamax<g0[i][j][k])
                        gamax=g0[i][j][k];
	      }
	    }
         }
        
         a=1.0*vamin/((400+count)*gamax);//+++++++++++++++++++++ 
  
         FILE *fpaaa;
         fpaaa=fopen("NEW_a.txt","a");
         fprintf(fpaaa,"vamin=%f,gamax=%f,a=%f \n",vamin,gamax,a);
        fprintf(fpdebug,"---   this loop >> vamin = %f,gamax = %f,a = %f \n",vamin,gamax,a);
         fclose(fpaaa); 
         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	          v_initial[i][j][k]=v_initial[i][j][k]-a*g0[i][j][k];
                  //v_initial[i][j][k]=v_initial[i][j][k];
	      }
	    }
         }
        
         FILE *fpv;
         fpv=fopen(FN2,"wb");
         for(j=0;j<vny;j++)
         {
	    for(i=0;i<vnx;i++)
	    { 
	      for(k=0;k<vnz;k++)
	      {
	         fwrite(&v_initial[i][j][k],4L,1,fpv);
            
	      }
	    }
         }
     
         fclose(fpv);
         printf("complete the update velocity  !\n");
        fprintf(fpdebug,"---   complete the update velocity  !\n");

     //fseek(fpvelrecording,vnx*vny*vnz*4L,1);
      for(j=0;j<ny;j++)
	  for(i=0;i<nx;i++)
	    for(k=0;k<nz;k++)
               fwrite(&v_initial[i][j][k],4L,1,fpvelrecording);
      }
      count=count+1;
 
}
//while(fabs(factor[0])>0.5);
while(count<200);
/******************************************/ 
      free3float(g);
      free3float(g_mig);
      free3float(g_tom);
      free3float(g0);
      free3float(g0_mig);
      free3float(g0_tom);
      free3float(ws);
      free3float(ws0);
      free3float(ggg0);
      free3float(prp);
      free3float(new_g);

  
      free3float(v_initial);
      free3float(p_initial);
      free3float(p_down);
      free3float(p_left);
      free3float(p_right);
      free3float(p_front);
      free3float(p_back);

      free3float(p_com);
      free3float(p_obs);


      count=count-1;

      if(myid==0)
      {
          printf("complete the count==%d ite  !!!!!!!!! \n",count);
          printf("**********************************************\n");
        fprintf(fpdebug,"---   complete the count ===== %d ite  !\n",count);
        fprintf(fpdebug,"--------------------------------------------------\n");
      }       
      fclose(fpvelrecording);
      fclose(fpdebug);

      MPI_Finalize();
}

/**************************************************************************************/
void shot_compare(int nx,int ny,int nt,int is,char FN1[],float ***p_obs,float ***p_initial,float ***p_com)
{
	
	int i,j,k;
	
	FILE *fp1;      

	zero3float(p_obs,nt,ny,nx);			
      
	if((fp1=fopen(FN1,"rb"))==NULL){exit(0);printf("The %s open error!",FN1);};	         
        fseek(fp1,(is-1)*nx*ny*nt*4L,0);

        for(j=0;j<ny;j++)
        {
	    for(i=0;i<nx;i++)
	    {
		for(k=0;k<nt;k++)
		{
			fread(&p_obs[i][j][k],4L,1,fp1);
		}
	    }
        }
	

	for(i=0;i<nx;i++)
	{
	   for(j=0;j<ny;j++)
	   {
              for(k=0;k<nt;k++)
              {
		   p_com[i][j][k]=p_initial[i][j][k]-p_obs[i][j][k];
	      }
	   }
	}

	fclose(fp1);

	
}


/*******************************************************************************/


void model(int nx,int ny,int nz,int vnx,int vny,int vnz,int nt,int npd,float dx,float dy,float dz,
          float vdx,float vdy,float vdz,float favg,float tmax,float dt,float dtout,float pfac,char FN1[],
          int ns_sxd,int xs_sxd,int ys_sxd,int zs_sxd,int is,float ***p_initial,float ***p_down,
          float ***p_left,float ***p_right,float ***p_front,float ***p_back,float ***ws0,int myid)
{


	void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float ysn,float zsn,int nx,int ny,int nz,float dt,float t,
                    float favg,float ***s,int wtype,float pi,int npd,int is);
      void update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,
                    float ***u0,float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,
                    float ***txx1,float ***rho,float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,
                    float *coffz1,float *coffz2);
      void update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,float ***v0,
                    float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***s,
                    float ***vp,float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,
                    float ***tz1,float ***tz0,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,
                    float *acoffz1,float *acoffz2);
      void abs_bc(float ***u1,float ***w1,float ***txx1,int nx,int ny,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dy,float dz,int nx,int ny,int nz,int nt,int ntout,int npd,float tmax,
                         float favg,float dtout,float dt,float ***vp0,float ndtt);
      void pad_vv(int nx,int ny,int nz,int npd,float ***ee);
	void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd);
	void current_shot(float ***vp0,float ***rho0,float ***vp,float ***rho,int nx,int ny,int nz,int npd,
                        int vnx,int vny,int vnz,int is);	
      void initial_coffe(float dt,float d0,int nx,int ny,int nz,float *coffx1,float *coffx2,float *coffy1,
                        float *coffy2,float *coffz1,float *coffz2,float *acoffx1,float *acoffx2,float *acoffy1,
                        float *acoffy2,float *acoffz1,float *acoffz2,int npd);
      void freea(float **a,int m);
      

          
	  int mm=4; 
          int hsx=1;
	  int i,j,k;

	  int ntout,wtype,it,ifx,ilx,ify,ily,jfz,jlz;
	  float pi,t,ndtt,d0;
	  
	  
	  FILE *fp1;

          wtype=1;
          

	  pi=3.141593;
	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;
     
	  float ***vp0;
	  float ***rho0;
	  float ***vp;
	  float ***rho;

	  float ***u0;
          float ***v0;
	  float ***w0;

          float ***u1;
          float ***v1;
	  float ***w1;

	  float ***txx0; 
          float ***txx1;
 
          float ***txx_snap;
	  
	  
	  float ***tx0;
	  float ***tx1;
          float ***ty0;
	  float ***ty1;
	  float ***tz0;
	  float ***tz1;
	  float ***s;
     
	  float c[4];
	  

          cal_c(mm,c);   

          vp0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          rho0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	  zero3float(vp0,vnz+2*npd,vny+2*npd,vnx+2*npd);  
          zero3float(rho0,vnz+2*npd,vny+2*npd,vnx+2*npd);                                              
     

    
                   
          read_file(FN1,vnx,vny,vnz,vp0,rho0,npd);          
      
          vp=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	  rho=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);  
                                            
          u0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          u1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          v0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          v1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          w0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          w1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 

          txx0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          txx1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 
         
          txx_snap=alloc3float(vnz,vny,vnx);

          tx0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          tx1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          ty0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          ty1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	  tz0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 
          tz1=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd); 
          s=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	 

          d0=get_constant(dx,dy,dz,nx,ny,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);
          dt=dt/1000;

          float *coffx1;float *coffx2;float *coffy1;float *coffy2;float *coffz1;float *coffz2;
          float *acoffx1;float *acoffx2;float *acoffy1;float *acoffy2;float *acoffz1;float *acoffz2;
          coffx1=alloc1float(vnx+2*npd);
          coffx2=alloc1float(vnx+2*npd);
          coffy1=alloc1float(vny+2*npd);
          coffy2=alloc1float(vny+2*npd);
	  coffz1=alloc1float(vnz+2*npd);
          coffz2=alloc1float(vnz+2*npd);

	  acoffx1=alloc1float(vnx+2*npd);
	  acoffx2=alloc1float(vnx+2*npd);
          acoffy1=alloc1float(vny+2*npd);
	  acoffy2=alloc1float(vny+2*npd);
	  acoffz1=alloc1float(vnz+2*npd);
	  acoffz2=alloc1float(vnz+2*npd);

          zero1float(coffx1,vnx+2*npd);
          zero1float(coffx2,vnx+2*npd);
          zero1float(coffz1,vnz+2*npd);
          zero1float(coffz2,vnz+2*npd);
          zero1float(coffy1,vny+2*npd);
          zero1float(coffy2,vny+2*npd);

          zero1float(acoffx1,vnx+2*npd);
          zero1float(acoffx2,vnx+2*npd);
          zero1float(acoffz1,vnz+2*npd);
          zero1float(acoffz2,vnz+2*npd);
          zero1float(acoffy1,vny+2*npd);
          zero1float(acoffy2,vny+2*npd);


	 
          initial_coffe(dt,d0,nx,ny,nz,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2,
                       acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2,npd);

/***********************************************************/
	   ndtt=(int)ndtt;


       //  printf("IS============%d \n",is);
/***********************************************/ 
           zero3float(p_initial,nt,ny,nx);
           zero3float(p_down,nt,ny,nx);
           zero3float(p_left,nt,nz,ny);
           zero3float(p_right,nt,nz,ny);
           zero3float(p_front,nt,nz,nx);
           zero3float(p_back,nt,nz,nx);

           zero3float(tx0,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(tx1,vnz+2*npd,vny+2*npd,vnx+2*npd);  
           zero3float(ty0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(ty1,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(tz0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(tz1,vnz+2*npd,vny+2*npd,vnx+2*npd); 

           zero3float(u0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(u1,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(v0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(v1,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(w0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(w1,vnz+2*npd,vny+2*npd,vnx+2*npd);

           zero3float(txx0,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(txx1,vnz+2*npd,vny+2*npd,vnx+2*npd);
       
           zero3float(txx_snap,vnz,vny,vnx);

           zero3float(vp,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(rho,vnz+2*npd,vny+2*npd,vnx+2*npd);
 

			 
           current_shot(vp0,rho0,vp,rho,nx,ny,nz,npd,vnx,vny,vnz,is);//*********there is no transform the vp and  rho

           pad_vv(nx,ny,nz,npd,vp); 

           pad_vv(nx,ny,nz,npd,rho);
       
      
	   for(i=0;i<=nx+2*npd-1;i++)
	   {
	     for(j=0;j<=ny+2*npd-1;j++)
	     {
                for(k=0;k<=nz+2*npd-1;k++)
                {  
		       vp[i][j][k]=rho[i][j][k]*(vp[i][j][k]*vp[i][j][k]);
                   rho[i][j][k]=1.0/rho[i][j][k];
		}		
	     }
	   }
		
	   t=0.0;

    for(it=0;it<nt;it++)
    { 
          if(it%50==0&&myid==0)
          printf("the model is====%d  ,  it===%d\n",is,it);
	    t=t+dt;
	    ptsource(pfac,xs_sxd,ys_sxd,zs_sxd,nx,ny,nz,dt,t,favg,s,wtype,pi,npd,is);
          update_vel(nx,ny,nz,npd,mm,dt,dx,dy,dz,u0,v0,w0,txx0,u1,v1,w1,txx1,rho,
                     c,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
          update_txx(nx,ny,nz,dt,dx,dy,dz,mm,u0,v0,w0,txx0,u1,v1,w1,txx1,s,vp,c,npd,
                  tx1,tx0,ty1,ty0,tz1,tz0,acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2);
       

          
	  for(i=npd;i<npd+nx;i++)  
	  {   
              for(j=npd;j<npd+ny;j++)
              {
		    p_initial[i-npd][j-npd][it]=txx1[i][j][npd+hsx-1];
                p_down[i-npd][j-npd][it]=txx1[i][j][npd+nz-1];
              }
	  }

          for(j=npd;j<npd+ny;j++)  
	  {   
              for(k=npd;k<npd+nz;k++)
              {
		    p_left[j-npd][k-npd][it]=txx1[npd+hsx-1][j][k];
                p_right[j-npd][k-npd][it]=txx1[npd+nx-1][j][k];
              }
	  }

          for(i=npd;i<npd+nx;i++)  
	  {   
              for(k=npd;k<npd+nz;k++)
              {
		    p_front[i-npd][k-npd][it]=txx1[i][npd+hsx-1][k];
                p_back[i-npd][k-npd][it]=txx1[i][npd+ny-1][k];
              }
	  }
          

	  for(i=0;i<nx+2*npd;i++)
	  {
		for(j=0;j<ny+2*npd;j++)
		{
                     for(k=0;k<nz+2*npd;k++)
                     {
			      u0[i][j][k]=u1[i][j][k];
                        v0[i][j][k]=v1[i][j][k];
			      w0[i][j][k]=w1[i][j][k];
			      tx0[i][j][k]=tx1[i][j][k];
                        ty0[i][j][k]=ty1[i][j][k];
			      tz0[i][j][k]=tz1[i][j][k];			
			      txx0[i][j][k]=txx1[i][j][k];
                     }
		}
	   }

           for(i=npd;i<npd+nx;i++)
	   {
	    	 for(j=npd;j<ny+npd;j++)
		 {
                      for(k=npd;k<nz+npd;k++)
                      {
			      ws0[i-npd][j-npd][k-npd]=ws0[i-npd][j-npd][k-npd]
                                                +(txx1[i][j][k]*txx1[i][j][k]*txx1[i][j][k]*txx1[i][j][k]);
                      }
		 }
	   }   
/******************************************** the snap dat output   ******************/
 /*          if(it==150&&is==1)
           {
              for(i=0;i<nx;i++)
	      {
		 for(j=0;j<ny;j++)
		 {
                     for(k=0;k<nz;k++)
                     {
			txx_snap[i][j][k]=txx1[i+npd][j+npd][k+npd];
                     }
		 }
	      }
              FILE *fpsnap;
              fpsnap=fopen("snap_it150.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&txx_snap[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }

  */
           

     }
        
/***********************************************/        
      	  free3float(vp0);
          free3float(rho0);
      	  free3float(vp);
          free3float(rho);

          free3float(u0);
          free3float(v0);
          free3float(w0);
          free3float(txx0);
          free3float(u1);
          free3float(v1);
          free3float(w1);
          free3float(txx1);

          free3float(txx_snap);

          free3float(tx0);
          free3float(tx1);
          free3float(ty0);
          free3float(ty1);
          free3float(tz0);
          free3float(tz1);
		  
	  free3float(s);
       
          free1float(coffx1);free1float(coffx2);
          free1float(coffz1);free1float(coffz2);
          free1float(coffy1);free1float(coffy2);
          free1float(acoffx1);free1float(acoffx2);
          free1float(acoffz1);free1float(acoffz2);
          free1float(acoffy1);free1float(acoffy2);

 
}
/**************************************************************************************************************/
void back_grad(int myid,int nx,int ny,int nz,int vnx,int vny,int vnz,int nt,int npd,float dx,float dy,float dz,
               float vdx,float vdy,float vdz,float favg,float tmax,float dt,float dtout,char FN1[],int ns_sxd,
               int xs_sxd,int ys_sxd,int zs_sxd,int is,float ***p_com,float ***p_up,float ***p_down,
               float ***p_left,float ***p_right,float ***p_front,float ***p_back,
               float ***g,float ***g_mig,float ***g_tom)
{


      void cal_c(int mm,float c[]);
	void ptsource(float pfac,float xsn,float ysn,float zsn,int nx,int ny,int nz,float dt,float t,float favg,
                    float ***s,int wtype,float pi,int npd,int is);
      void update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,float ***u0,
                    float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,
                    float ***rho,float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,
                    float *coffz2);
      void for_update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,float ***u0,
                    float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,
                    float ***rho,float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,
                    float *coffz2);
      void update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,float ***v0,
                    float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***s,
                    float ***vp,float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,
                    float ***tz1,float ***tz0,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,
                    float *acoffz1,float *acoffz2);
      void for_update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,
                    float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,
                    float ***vp,float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,float ***tz1,
                    float ***tz0,float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,
                    float *acoffz2);
      void abs_bc(float ***u1,float ***w1,float ***txx1,int nx,int ny,int nz,int npd,float absbord[]);
      float get_constant(float dx,float dy,float dz,int nx,int ny,int nz,int nt,int ntout,int npd,float tmax,
                         float favg,float dtout,float dt,float ***vp0,float ndtt);
      void pad_vv(int nx,int ny,int nz,int npd,float ***ee);
	void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd);
	void current_shot(float ***vp0,float ***rho0,float ***vp,float ***rho,int nx,int ny,int nz,int npd,
                        int vnx,int vny,int vnz,int is);	
      void initial_coffe(float dt,float d0,int nx,int ny,int nz,float *coffx1,float *coffx2,float *coffy1,
                         float *coffy2,float *coffz1,float *coffz2,float *acoffx1,float *acoffx2,float *acoffy1,
                         float *acoffy2,float *acoffz1,float *acoffz2,int npd);
      void freea(float **a,int m);
      

          
	  int mm=4; 
          int hsx=1;
	  int i,j,k;

	  int ntout,wtype,it,ifx,ilx,ify,ily,jfz,jlz;
	  float pi,ndtt,d0;
	  
	  
	  FILE *fp1;

          wtype=1;
          
          

	  pi=3.141593;
	  ndtt=dtout/dt;
	  ntout=(int)(1000*tmax/dtout+0.5)+1;
     
	  float ***vp0;
	  float ***rho0;
	  float ***vp;
	  float ***rho;

	    float ***u0_f;      float ***u0_b;
          float ***v0_f;      float ***v0_b;
	    float ***w0_f;      float ***w0_b;

          float ***u1_f;      float ***u1_b;
          float ***v1_f;      float ***v1_b;
	    float ***w1_f;      float ***w1_b;

	    float ***txx0_f;    float ***txx0_b;
          float ***txx1_f;    float ***txx1_b;
	  
	    float ***tx0_f;     float ***tx0_b;
	    float ***tx1_f;     float ***tx1_b;
          float ***ty0_f;     float ***ty0_b;
	    float ***ty1_f;     float ***ty1_b;
	    float ***tz0_f;     float ***tz0_b;
	    float ***tz1_f;     float ***tz1_b;
	    float ***s;

          float ***txx_snap;
     
	  float c[4];
	  
         

          cal_c(mm,c);   

          vp0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          rho0=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	    zero3float(vp0,vnz+2*npd,vny+2*npd,vnx+2*npd);  
          zero3float(rho0,vnz+2*npd,vny+2*npd,vnx+2*npd);  
        
/***************************************************************/                   
          read_file(FN1,vnx,vny,vnz,vp0,rho0,npd);          
        
          vp=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	    rho=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);  
                                            
          u0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);    u0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          u1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);    u1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          v0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);    v0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          v1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);    v1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          w0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);    w0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          w1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);    w1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);

          txx0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);  txx0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          txx1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);  txx1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
         
          tx0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);   tx0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          tx1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);   tx1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          ty0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);   ty0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          ty1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);   ty1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
	    tz0_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);   tz0_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);
          tz1_f=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);   tz1_b=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);

          s=alloc3float(vnz+2*npd,vny+2*npd,vnx+2*npd);

          txx_snap=alloc3float(vnz,vny,vnx);
	
/*****************************************************************************************/
          d0=get_constant(dx,dy,dz,nx,ny,nz,nt,ntout,npd,tmax,favg,dtout,dt,vp0,ndtt);

	  float ***f0;
	  float ***b0;	
	  float ***f1;
	  float ***b1;
	  float ***f2;
	  float ***b2;
	  float ***t;
	  float ***a;
	  float ***vvv;
/******************** poynting wave propagation alloc ****************/
        float ***f_z,***f_up0,***f_up1,***f_up2;
        float ***f_down0,***f_down1,***f_down2;
        float ***b_z,***b_up0,***b_down0;
          
	  float ***t_tom;
	  float ***t_mig;



          f_z      =  alloc3float(nz,ny,nx);
          f_up0    =  alloc3float(nz,ny,nx);
          f_up1    =  alloc3float(nz,ny,nx);
          f_up2    =  alloc3float(nz,ny,nx);
          f_down0  =  alloc3float(nz,ny,nx);
          f_down1  =  alloc3float(nz,ny,nx);
          f_down2  =  alloc3float(nz,ny,nx);
          b_z      =  alloc3float(nz,ny,nx);
          b_up0    =  alloc3float(nz,ny,nx);
          b_down0  =  alloc3float(nz,ny,nx);

          t_tom      =  alloc3float(nz,ny,nx);
          t_mig      =  alloc3float(nz,ny,nx);


          zero3float(f_z,nz,ny,nx);
          zero3float(f_up0,nz,ny,nx);
          zero3float(f_up1,nz,ny,nx);
          zero3float(f_up2,nz,ny,nx);
          zero3float(f_down0,nz,ny,nx);
          zero3float(f_down1,nz,ny,nx);
          zero3float(f_down2,nz,ny,nx);
          zero3float(b_z,nz,ny,nx);
          zero3float(b_up0,nz,ny,nx);
          zero3float(b_down0,nz,ny,nx);

          zero3float(t_tom,nz,ny,nx);
          zero3float(t_mig,nz,ny,nx);

          zero3float(g_tom,nz,ny,nx);
          zero3float(g_mig,nz,ny,nx);
/******************** poynting wave propagation alloc ****************/

          f0=alloc3float(nz,ny,nx);
          f1=alloc3float(nz,ny,nx);
          f2=alloc3float(nz,ny,nx);
          b0=alloc3float(nz,ny,nx);
          b1=alloc3float(nz,ny,nx);
          b2=alloc3float(nz,ny,nx);
          t=alloc3float(nz,ny,nx);
          a=alloc3float(vnz,vny,vnx);
          vvv=alloc3float(vnz,vny,vnx);
          
          zero3float(f0,nz,ny,nx);
          zero3float(f1,nz,ny,nx);
          zero3float(f2,nz,ny,nx);
          zero3float(b0,nz,ny,nx);
          zero3float(b1,nz,ny,nx);
          zero3float(b2,nz,ny,nx);
          zero3float(t,vnz,vny,vnx);
          zero3float(a,vnz,vny,vnx);
          zero3float(vvv,vnz,vny,vnx);
          
          for(i=npd;i<vnx+npd;i++)
          {   
             for(j=npd;j<vny+npd;j++)
             {
                for(k=npd;k<vnz+npd;k++)
                {
                   vvv[i-npd][j-npd][k-npd]=vp0[i][j][k];
                }
             }
          }
                   

          dt=dt/1000;
        
          float *coffx1;float *coffx2;float *coffy1;float *coffy2;float *coffz1;float *coffz2;
          float *acoffx1;float *acoffx2;float *acoffy1;float *acoffy2;float *acoffz1;float *acoffz2;
          coffx1=alloc1float(vnx+2*npd);
          coffx2=alloc1float(vnx+2*npd);
          coffy1=alloc1float(vny+2*npd);
          coffy2=alloc1float(vny+2*npd);
	    coffz1=alloc1float(vnz+2*npd);
          coffz2=alloc1float(vnz+2*npd);
          
	    acoffx1=alloc1float(vnx+2*npd);
	    acoffx2=alloc1float(vnx+2*npd);
          acoffy1=alloc1float(vny+2*npd);
	    acoffy2=alloc1float(vny+2*npd);
	    acoffz1=alloc1float(vnz+2*npd);
	    acoffz2=alloc1float(vnz+2*npd);
          
          zero1float(coffx1,vnx+2*npd);
          zero1float(coffx2,vnx+2*npd);
          zero1float(coffz1,vnz+2*npd);
          zero1float(coffz2,vnz+2*npd);
          zero1float(coffy1,vny+2*npd);
          zero1float(coffy2,vny+2*npd);
          
          zero1float(acoffx1,vnx+2*npd);
          zero1float(acoffx2,vnx+2*npd);
          zero1float(acoffz1,vnz+2*npd);
          zero1float(acoffz2,vnz+2*npd);
          zero1float(acoffy1,vny+2*npd);
          zero1float(acoffy2,vny+2*npd);
          
/********************************************************************************************************/
          initial_coffe(dt,d0,nx,ny,nz,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2,
                        acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2,npd);

	   ndtt=(int)ndtt;

          
       //  printf("IS============%d \n",is);
/***********************************************/ 

           zero3float(tx0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);    zero3float(tx0_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(tx1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);    zero3float(tx1_b,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(ty0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);    zero3float(ty0_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(ty1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);    zero3float(ty1_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(tz0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);    zero3float(tz0_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(tz1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);    zero3float(tz1_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
      
           zero3float(u0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);     zero3float(u0_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(u1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);     zero3float(u1_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(v0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);     zero3float(v0_b,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(v1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);     zero3float(v1_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(w0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);     zero3float(w0_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(w1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);     zero3float(w1_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
      
           zero3float(txx0_f,vnz+2*npd,vny+2*npd,vnx+2*npd);   zero3float(txx0_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
           zero3float(txx1_f,vnz+2*npd,vny+2*npd,vnx+2*npd);   zero3float(txx1_b,vnz+2*npd,vny+2*npd,vnx+2*npd);
       
           zero3float(txx_snap,vnz,vny,vnx);

           zero3float(vp,vnz+2*npd,vny+2*npd,vnx+2*npd); 
           zero3float(rho,vnz+2*npd,vny+2*npd,vnx+2*npd);
 
           
			 
           current_shot(vp0,rho0,vp,rho,nx,ny,nz,npd,vnx,vny,vnz,is);//*********there is no transform the vp and  rho
       
           pad_vv(nx,ny,nz,npd,vp); 

           pad_vv(nx,ny,nz,npd,rho);
        
      
	   for(i=0;i<=nx+2*npd-1;i++)
	   {
	     for(j=0;j<=ny+2*npd-1;j++)
	     {
                for(k=0;k<=nz+2*npd-1;k++)
                {  
		       vp[i][j][k]=rho[i][j][k]*(vp[i][j][k]*vp[i][j][k]);
                   rho[i][j][k]=1.0/rho[i][j][k];
		}		
	     }
	   }
		
/**********************************************************************************************************************/
/********************************************** the reserve it is start ***********************************************/
/**********************************************************************************************************************/
    for(it=nt-1;it>=0;it--)
    { 
         if(it%50==0&&myid==0)
          printf("the back_propagation is====%d  ,  it===%d\n",is,it);
     
/********************************** the six bondary propagation ******************************************/
          for(i=npd;i<nx+npd;i++)
          {
             for(j=npd;j<ny+npd;j++)
             {
                 txx1_f[i][j][npd]=p_up[i-npd][j-npd][it];
                 txx1_f[i][j][npd+nz-1]=p_down[i-npd][j-npd][it];
             }
          }
          
          for(j=npd;j<ny+npd;j++)
          {
             for(k=npd;k<nz+npd;k++)
             {
                 txx0_f[npd][j][k]=p_left[j-npd][k-npd][it];
                 txx0_f[npd+nx-1][j][k]=p_right[j-npd][k-npd][it];
             }
          }
          for(i=npd;i<nx+npd;i++)
          {
             for(k=npd;k<nz+npd;k++)
             {
                 txx0_f[i][npd][k]=p_front[i-npd][k-npd][it];
                 txx0_f[i][npd+ny-1][k]=p_back[i-npd][k-npd][it];
             }
          }
	             
          for_update_vel(nx,ny,nz,npd,mm,dt,dx,dy,dz,u0_f,v0_f,w0_f,txx0_f,u1_f,v1_f,w1_f,txx1_f,
                        rho,c,coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
          for_update_txx(nx,ny,nz,dt,dx,dy,dz,mm,u0_f,v0_f,w0_f,txx0_f,u1_f,v1_f,w1_f,txx1_f,
                        vp,c,npd,tx1_f,tx0_f,ty1_f,ty0_f,tz1_f,tz0_f,acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2);
                     
          zero3float(f_z,nz,ny,nx);
          zero3float(f_up0,nz,ny,nx);
          zero3float(f_down0,nz,ny,nx);

          for(i=npd;i<nx+npd;i++)
          { 
             for(j=npd;j<ny+npd;j++)
             { 
                for(k=npd;k<nz+npd;k++)
                {
                    f0[i-npd][j-npd][k-npd]=txx1_f[i][j][k];
                    f_z[i-npd][j-npd][k-npd]=(w1_f[i][j][k]-w0_f[i][j][k])*(txx1_f[i][j][k]-txx0_f[i][j][k])/(dt*dt);
                    if(f_z[i-npd][j-npd][k-npd]>0.0){
                          f_up0[i-npd][j-npd][k-npd]=f0[i-npd][j-npd][k-npd];
                          f_down0[i-npd][j-npd][k-npd]=0.0;
                    }else{
                          f_up0[i-npd][j-npd][k-npd]=0.0;
                          f_down0[i-npd][j-npd][k-npd]=f0[i-npd][j-npd][k-npd];
                          }
                                  
                }
             }
          }
         

          for(i=0;i<nx+2*npd;i++)
	  {
		for(j=0;j<ny+2*npd;j++)
		{
                     for(k=0;k<nz+2*npd;k++)
                     {
			      u0_f[i][j][k]=u1_f[i][j][k];
                        v0_f[i][j][k]=v1_f[i][j][k];
			      w0_f[i][j][k]=w1_f[i][j][k];
			      tx0_f[i][j][k]=tx1_f[i][j][k];
                        ty0_f[i][j][k]=ty1_f[i][j][k];
			      tz0_f[i][j][k]=tz1_f[i][j][k];			
			      txx0_f[i][j][k]=txx1_f[i][j][k];
                     }
		}
	   }

/********************************** get back propagation wavefield ***************************************/
          
          for(i=npd;i<nx+npd;i++)
          {
             for(j=npd;j<ny+npd;j++)
             {
                 s[i][j][npd+hsx-1]=-p_com[i-npd][j-npd][it];
             }
          }
                
          update_vel(nx,ny,nz,npd,mm,dt,dx,dy,dz,u0_b,v0_b,w0_b,txx0_b,u1_b,v1_b,w1_b,txx1_b,rho,c,
                      coffx1,coffx2,coffy1,coffy2,coffz1,coffz2);
          update_txx(nx,ny,nz,dt,dx,dy,dz,mm,u0_b,v0_b,w0_b,txx0_b,u1_b,v1_b,w1_b,txx1_b,
                     s,vp,c,npd,tx1_b,tx0_b,ty1_b,ty0_b,tz1_b,tz0_b,acoffx1,acoffx2,acoffy1,acoffy2,acoffz1,acoffz2);
                  
          zero3float(b_z,nz,ny,nx);
          zero3float(b_up0,nz,ny,nx);
          zero3float(b_down0,nz,ny,nx);

          for(i=npd;i<nx+npd;i++)
          { 
             for(j=npd;j<ny+npd;j++)
             { 
                for(k=npd;k<nz+npd;k++)
                {
                    b0[i-npd][j-npd][k-npd]=txx1_b[i][j][k];
                    b_z[i-npd][j-npd][k-npd]=(w1_b[i][j][k]-w0_b[i][j][k])*(txx1_b[i][j][k]-txx0_b[i][j][k])/(dt*dt);
                    if(b_z[i-npd][j-npd][k-npd]>0.0){
                          b_up0[i-npd][j-npd][k-npd]=b0[i-npd][j-npd][k-npd];
                          b_down0[i-npd][j-npd][k-npd]=0.0;
                    }else{
                          b_up0[i-npd][j-npd][k-npd]=0.0;
                          b_down0[i-npd][j-npd][k-npd]=b0[i-npd][j-npd][k-npd];
                          }
                }
             }
          }
               
	  for(i=0;i<nx+2*npd;i++)
	  {
		for(j=0;j<ny+2*npd;j++)
		{
                     for(k=0;k<nz+2*npd;k++)
                     {
			      u0_b[i][j][k]=u1_b[i][j][k];
                        v0_b[i][j][k]=v1_b[i][j][k];
			      w0_b[i][j][k]=w1_b[i][j][k];
			      tx0_b[i][j][k]=tx1_b[i][j][k];
                        ty0_b[i][j][k]=ty1_b[i][j][k];
			      tz0_b[i][j][k]=tz1_b[i][j][k];			
			      txx0_b[i][j][k]=txx1_b[i][j][k];
                     }
		}
	   }

        if(it<=nt-2)
            {
            if(it<=nt-3)
                  {
                  for(i=0;i<nx-0;i++)
                        {
                    for(j=0;j<ny-0;j++)
                           {   
                      for(k=0;k<nz-0;k++)
                              {                          
                          t[i][j][k]=t[i][j][k]+(f2[i][j][k]-2*f1[i][j][k]+f0[i][j][k])*(b0[i][j][k])/(dt*dt);
                          t_mig[i][j][k]=t_mig[i][j][k]
                                        +(f_up2[i][j][k]-2*f_up1[i][j][k]+f_up0[i][j][k])*b_down0[i][j][k]/(dt*dt)
                                        +(f_down2[i][j][k]-2*f_down1[i][j][k]+f_down0[i][j][k])*b_up0[i][j][k]/(dt*dt);
                          t_tom[i][j][k]=t_tom[i][j][k]
                                        +(f_up2[i][j][k]-2*f_up1[i][j][k]+f_up0[i][j][k])*b_up0[i][j][k]/(dt*dt)
                                        +(f_down2[i][j][k]-2*f_down1[i][j][k]+f_down0[i][j][k])*b_down0[i][j][k]/(dt*dt);
                              }
                           }
                        }
                  }
                  
             for(i=0;i<nx;i++)
                  {
                for(j=0;j<ny;j++)
                     {
                  for(k=0;k<nz;k++)
                        {
                     f2[i][j][k]=f1[i][j][k];
                     f_up2[i][j][k]=f_up1[i][j][k];
                     f_down2[i][j][k]=f_down1[i][j][k];
                        }
                     }
                  }
            }
            
         for(i=0;i<nx;i++)
            {
           for(j=0;j<ny;j++)
               {
             for(j=0;j<ny;j++)
                  {
                f1[i][j][k]=f0[i][j][k];
                f_up1[i][j][k]=f_up0[i][j][k];
                f_down1[i][j][k]=f_down0[i][j][k];
                  }
               }
            }   
             
            
/******************************************** the snap dat output   ******************/
       
           if(it==150&&is==1)
           {
              FILE *fpsnap;
              fpsnap=fopen("snapb_it150.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&b0[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }
           if(it==250&&is==1)
           {
              FILE *fpsnap;
              fpsnap=fopen("snapb_it250.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&b0[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }
           if(it==300&&is==1)
           {
              FILE *fpsnap;
              fpsnap=fopen("snapb_it300.dat","wb");
              for(j=0;j<ny;j++)
	      {
		 for(i=0;i<nx;i++)
		 {
                     for(k=0;k<nz;k++)
                     {
			fwrite(&b0[i][j][k],4L,1,fpsnap);
                     }
		 }
	      }
              fclose(fpsnap);
           }
  
           

     }//The nt -> 0 loop ending
          
          for(i=0;i<vnx;i++)
	  {
	     for(j=0;j<vny;j++)
	     {
                for(k=0;k<vnz;k++)
                {
	            g[i][j][k]=2*t[i][j][k]/(pow(vvv[i][j][k],3));
	            g_tom[i][j][k]=2*t_tom[i][j][k]/(pow(vvv[i][j][k],3));
	            g_mig[i][j][k]=2*t_mig[i][j][k]/(pow(vvv[i][j][k],3));
                  //g_poynting[i][j][k]=g_mig[i][j][k]+g_tom[i][j][k];
                }
             }
	  } 
      
/***********************************************/        
          free3float(vp0);
          free3float(rho0);
          free3float(vp);
          free3float(rho);

          free3float(u0_f);      free3float(u0_b);
          free3float(v0_f);      free3float(v0_b);
          free3float(w0_f);      free3float(w0_b);
          free3float(txx0_f);    free3float(txx0_b);
          free3float(u1_f);      free3float(u1_b);
          free3float(v1_f);      free3float(v1_b);
          free3float(w1_f);      free3float(w1_b);
          free3float(txx1_f);    free3float(txx1_b);

          free3float(tx0_f);     free3float(tx0_b);
          free3float(tx1_f);     free3float(tx1_b);
          free3float(ty0_f);     free3float(ty0_b);
          free3float(ty1_f);     free3float(ty1_b);
          free3float(tz0_f);     free3float(tz0_b);
          free3float(tz1_f);     free3float(tz1_b);

          free3float(f_z);          
          free3float(f_up0);       free3float(f_up1);         free3float(f_up2);
          free3float(f_down0);     free3float(f_down1);       free3float(f_down2); 
          free3float(b_z);          
          free3float(b_up0);       free3float(b_down0);  
          free3float(t_tom);       free3float(t_mig);  
        //  free3float(g_tom);       free3float(g_mig);      

          free3float(txx_snap);
		  
	    free3float(s);
       
          free1float(coffx1);free1float(coffx2);
          free1float(coffz1);free1float(coffz2);
          free1float(coffy1);free1float(coffy2);
          free1float(acoffx1);free1float(acoffx2);
          free1float(acoffz1);free1float(acoffz2);
          free1float(acoffy1);free1float(acoffy2);

 
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



/**************************************************************************************************************/
void ptsource(float pfac,float xsn,float ysn,float zsn,int nx,int ny,int nz,float dt,float t,float favg,
               float ***s,int wtype,float pi,int npd,int is)
{

	    float get_source(float ts,float favg,int wtype);

	    int i,j,k,ixs,iys,izs,x,y,z;
	    float tdelay,ts,source,fs;
           
            

	    for(i=0;i<=nx+2*npd-1;i++)
	    {
		   for(j=0;j<=ny+2*npd-1;j++)
		   {
                       for(k=0;k<nz+2*npd;k++)
			  s[i][j][k]=0.0*dt*favg;
		   }
	    }

	    tdelay=1.0/favg;
            ts=t-tdelay;
            
           // fs=xsn;
       if(t<=2*tdelay)
       {
            source=get_source(ts,favg,wtype);            
	    ixs = (int)(xsn+0.5)+npd-1;
            iys = (int)(ysn+0.5)+npd-1;
            izs = (int)(zsn+0.5)+npd-1;

            for(i=ixs-3;i<=ixs+3;i++)
	    { 
		  for(j=iys-3;j<=iys+3;j++)
		  {  
                      for(k=izs-3;k<=izs+3;k++)
                      {
			  x=i-ixs;
                          y=j-iys;
                          z=k-izs;
                          s[i][j][k]=pfac*source*exp(-z*z-y*y-x*x);
                      }
		  }
	    }
	}

}
/*float get_source(float ts,float favg,int wtype)
{
        int i,j,nt,ntfft;
        int nf1,nf2,nf3,nf4;
	    float x,pi;
	    float source=0.0;
	    float  f1,f2,dw,fw;
        float tmpp=0.0;
        float dt=0.5/1000;
        float tdelay;   
        float    *trace;
        complex  *ctrace;
        tdelay=1.0/favg;
        nt=4001;

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
/*
float get_source(float ts,float favg,int wtype)
 {
		  float x,xx,pi;
		  float source=0.0,pi2;
		  pi=4*atan(1.0);
		  pi2=pi*pi;

		  if(wtype==1)//ricker wavelet
		  {
		      x=favg*ts;
		      xx=x*x;
		      source=(1-2*pi2*xx)*pow(exp,-pi2*xx);
		  }
		  if(wtype==2)//derivative of gaussian
		  {
		      x=(-4)*favg*favg*pi2/log(0.1);
		      source=(-2)*pi2*(ts)*pow(exp,-x*ts*ts);
          }
          if(wtype==3)//derivative of gaussian
          {
              x=(-1)*favg*favg*pi2/log(0.1);
              source=pow(exp,-x*ts*ts);
          }

		  return (source);

}
*/


/**************************************************************************************************************/
float get_source(float ts,float favg,int wtype)
{
	    float x,pi;
	    float source=0.0;
       
            

	    pi=3.14;
            x=(favg*(ts))*(favg*(ts));
		   
            source=(1-2*pi*pi*(x))*exp(-(pi*pi*x));
	    return (source);


}


/**************************************************************************************************************/
void update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,float ***u0,float ***v0,
                float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***rho,
                float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,float *coffz1,float *coffz2)
{
		 int ii,i,j,k;
		 float dtxx,dtxz,dtx,dty,dtz;


                 

		 dtx=dt/dx;
                 dty=dt/dy;
		 dtz=dt/dz;
         
		 for(k=mm;k<=(2*npd+nz-mm-1);k++)
		 {
		   for(j=mm;j<=(2*npd+ny-mm-1);j++)
	           {
                       for(i=mm;i<=(2*npd+nx-mm-1);i++)
                       {
			   u1[i][j][k]=coffx2[i]*u0[i][j][k]-coffx1[i]*dtx*rho[i][j][k]*(c[0]*(txx0[i+1][j][k]-txx0[i][j][k])+c[1]*(txx0[i+2][j][k]-txx0[i-1][j][k])+c[2]*(txx0[i+3][j][k]-txx0[i-2][j][k])+c[3]*(txx0[i+4][j][k]-txx0[i-3][j][k]));
                       }
		   }
		 }

                 for(k=mm;k<=(2*npd+nz-mm-1);k++)
		 {
		   for(i=mm;i<=(2*npd+nx-mm-1);i++)
	           {
                       for(j=mm;j<=(2*npd+ny-mm-1);j++)
                       {
			   v1[i][j][k]=coffy2[j]*v0[i][j][k]-coffy1[j]*dty*rho[i][j][k]*(c[0]*(txx0[i][j+1][k]-txx0[i][j][k])+c[1]*(txx0[i][j+2][k]-txx0[i][j-1][k])+c[2]*(txx0[i][j+3][k]-txx0[i][j-2][k])+c[3]*(txx0[i][j+4][k]-txx0[i][j-3][k]));
                       }
		   }
		 }

         
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 { 
		   for(j=mm;j<=(2*npd+ny-mm-1);j++)
		   {
                       for(k=mm;k<=(2*npd+nz-mm-1);k++)
                       {
			   w1[i][j][k]=coffz2[k]*w0[i][j][k]-coffz1[k]*dtz*rho[i][j][k]*(c[0]*(txx0[i][j][k+1]-txx0[i][j][k])+c[1]*(txx0[i][j][k+2]-txx0[i][j][k-1])+c[2]*(txx0[i][j][k+3]-txx0[i][j][k-2])+c[3]*(txx0[i][j][k+4]-txx0[i][j][k-3]));
                       }
		   }
		 }
		

}

/**************************************************************************************************************/
void for_update_vel(int nx,int ny,int nz,int npd,int mm,float dt,float dx,float dy,float dz,float ***u0,
                    float ***v0,float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,
                    float ***rho,float c[],float *coffx1,float *coffx2,float *coffy1,float *coffy2,
                    float *coffz1,float *coffz2)
{
		 int ii,i,j,k;
		 float dtxx,dtxz,dtx,dty,dtz;


                 

		 dtx=dt/dx;
                 dty=dt/dy;
		 dtz=dt/dz;
         
		 for(k=mm;k<=(2*npd+nz-mm-1);k++)
		 {
		   for(j=mm;j<=(2*npd+ny-mm-1);j++)
	           {
                       for(i=mm;i<=(2*npd+nx-mm-1);i++)
                       {
			   u1[i][j][k]=coffx2[i]*u0[i][j][k]-coffx1[i]*dtx*rho[i][j][k]*(c[0]*(txx0[i+1][j][k]-txx0[i][j][k])+c[1]*(txx0[i+2][j][k]-txx0[i-1][j][k])+c[2]*(txx0[i+3][j][k]-txx0[i-2][j][k])+c[3]*(txx0[i+4][j][k]-txx0[i-3][j][k]));
                       }
		   }
		 }

                 for(k=mm;k<=(2*npd+nz-mm-1);k++)
		 {
		   for(i=mm;i<=(2*npd+nx-mm-1);i++)
	           {
                       for(j=mm;j<=(2*npd+ny-mm-1);j++)
                       {
			   v1[i][j][k]=coffy2[j]*v0[i][j][k]-coffy1[j]*dty*rho[i][j][k]*(c[0]*(txx0[i][j+1][k]-txx0[i][j][k])+c[1]*(txx0[i][j+2][k]-txx0[i][j-1][k])+c[2]*(txx0[i][j+3][k]-txx0[i][j-2][k])+c[3]*(txx0[i][j+4][k]-txx0[i][j-3][k]));
                       }
		   }
		 }

         
		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 { 
		   for(j=mm;j<=(2*npd+ny-mm-1);j++)
		   {
                       for(k=mm;k<=(2*npd+nz-mm-1);k++)
                       {
			   w1[i][j][k]=coffz2[k]*w0[i][j][k]-coffz1[k]*dtz*rho[i][j][k]*(c[0]*(txx0[i][j][k+1]-txx0[i][j][k])+c[1]*(txx0[i][j][k+2]-txx0[i][j][k-1])+c[2]*(txx0[i][j][k+3]-txx0[i][j][k-2])+c[3]*(txx0[i][j][k+4]-txx0[i][j][k-3]));
                       }
		   }
		 }
		

}

/**************************************************************************************************************/
void update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,float ***v0,
               float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***s,float ***vp,
               float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,float ***tz1,float ***tz0,
               float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2)
{

		 float dux,dwz;
		 int i,j,k,ii;
		 float dtx,dty,dtz;


                 

		 dtx=dt/dx;
                 dty=dt/dy;
		 dtz=dt/dz;

		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
		    for(j=mm;j<=(2*npd+ny-mm-1);j++)
		    {
		       for(k=mm;k<=(2*npd+nz-mm-1);k++)
                       {
                            
		    	   tx1[i][j][k]=acoffx2[i]*tx0[i][j][k]-acoffx1[i]*vp[i][j][k]*dtx*(c[0]*(u1[i][j][k]-u1[i-1][j][k])+c[1]*(u1[i+1][j][k]-u1[i-2][j][k])+c[2]*(u1[i+2][j][k]-u1[i-3][j][k])+c[3]*(u1[i+3][j][k]-u1[i-4][j][k]));

                           ty1[i][j][k]=acoffy2[j]*ty0[i][j][k]-acoffy1[j]*vp[i][j][k]*dty*(c[0]*(v1[i][j][k]-v1[i][j-1][k])+c[1]*(v1[i][j+1][k]-v1[i][j-2][k])+c[2]*(v1[i][j+2][k]-v1[i][j-3][k])+c[3]*(v1[i][j+3][k]-v1[i][j-4][k]));

                           tz1[i][j][k]=acoffz2[k]*tz0[i][j][k]-acoffz1[k]*vp[i][j][k]*dtz*(c[0]*(w1[i][j][k]-w1[i][j][k-1])+c[1]*(w1[i][j][k+1]-w1[i][j][k-2])+c[2]*(w1[i][j][k+2]-w1[i][j][k-3])+c[3]*(w1[i][j][k+3]-w1[i][j][k-4]));
                    
			   txx1[i][j][k]=tx1[i][j][k]+ty1[i][j][k]+tz1[i][j][k]+s[i][j][k];
                       }
	            }
		 }

}           

/**************************************************************************************************************/
void for_update_txx(int nx,int ny,int nz,float dt,float dx,float dy,float dz,int mm,float ***u0,float ***v0,
               float ***w0,float ***txx0,float ***u1,float ***v1,float ***w1,float ***txx1,float ***vp,
               float c[],int npd,float ***tx1,float ***tx0,float ***ty1,float ***ty0,float ***tz1,float ***tz0,
               float *acoffx1,float *acoffx2,float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2)
{

		 float dux,dwz;
		 int i,j,k,ii;
		 float dtx,dty,dtz;


                 

		 dtx=dt/dx;
                 dty=dt/dy;
		 dtz=dt/dz;

		 for(i=mm;i<=(2*npd+nx-mm-1);i++)
		 {
		    for(j=mm;j<=(2*npd+ny-mm-1);j++)
		    {
		       for(k=mm;k<=(2*npd+nz-mm-1);k++)

                       {
                            
		    	   tx1[i][j][k]=acoffx2[i]*tx0[i][j][k]-acoffx1[i]*vp[i][j][k]*dtx*(c[0]*(u1[i][j][k]-u1[i-1][j][k])+c[1]*(u1[i+1][j][k]-u1[i-2][j][k])+c[2]*(u1[i+2][j][k]-u1[i-3][j][k])+c[3]*(u1[i+3][j][k]-u1[i-4][j][k]));

                           ty1[i][j][k]=acoffy2[j]*ty0[i][j][k]-acoffy1[j]*vp[i][j][k]*dty*(c[0]*(v1[i][j][k]-v1[i][j-1][k])+c[1]*(v1[i][j+1][k]-v1[i][j-2][k])+c[2]*(v1[i][j+2][k]-v1[i][j-3][k])+c[3]*(v1[i][j+3][k]-v1[i][j-4][k]));

                           tz1[i][j][k]=acoffz2[k]*tz0[i][j][k]-acoffz1[k]*vp[i][j][k]*dtz*(c[0]*(w1[i][j][k]-w1[i][j][k-1])+c[1]*(w1[i][j][k+1]-w1[i][j][k-2])+c[2]*(w1[i][j][k+2]-w1[i][j][k-3])+c[3]*(w1[i][j][k+3]-w1[i][j][k-4]));
                    
			   txx1[i][j][k]=tx1[i][j][k]+ty1[i][j][k]+tz1[i][j][k];
                       }
	            }
		 }

}           


           
/**************************************************************************************************************/
float get_constant(float dx,float dy,float dz,int nx,int ny,int nz,int nt,int ntout,int npd,float tmax,
                   float favg,float dtout,float dt,float ***vp0,float ndtt)
{
		 int i,j,k;
		 float vpmax,vpmin,H_min;
		 float dt_max,dx_max,dz_max,dy_max,d0;


                 

		 vpmax=vp0[npd][npd][npd];
		 vpmin=vp0[npd][npd][npd];
                 for(j=npd;j<ny+npd;j++)
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(k=npd;k<nz+npd;k++)
			 {
				 if(vpmax<vp0[i][j][k]) vpmax=vp0[i][j][k];
				 if(vpmin>vp0[i][j][k]) vpmin=vp0[i][j][k];
			 }
		 }

             
		 d0=3.0*vpmax*log(100000.0)/(2.0*npd*dx);

                 H_min=dx;
		 if(dy<H_min) H_min=dy;
                 if(dz<H_min) H_min=dz;
		 


/*====== determine time sampling interval to ensure stability====*/

		 dt_max=0.5*1000*H_min/vpmax;
                 dx_max=vpmin/favg*0.2;
                 dy_max=dx_max;
                 dz_max=dx_max;

             
                if(dx_max<dx)
                { 
                   printf("dx_max===%f, vpmin===%f, favg===%f \n",dx_max,vpmin,favg);
		   printf("YOU NEED HAVE TO REDEFINE DX ! \n");
                   exit(0);
		 }
                 if(dy_max<dy)
                { 
                   printf("dy_max===%f, vpmin===%f, favg===%f \n",dy_max,vpmin,favg);
		   printf("YOU NEED HAVE TO REDEFINE DY ! \n");
                   exit(0);
		 }
                 if(dz_max<dz)
		 {
		   printf("YOU NEED HAVE TO REDEFINE DZ ! \n");
                   exit(0);
		 }
	         if(dt_max<dt)
		 {printf("dt_max===%f, H_min===%f, vpmax===%f \n",dt_max,H_min,vpmax);
		   printf("YOU NEED HAVE TO REDEFINE dt ! \n");
                   exit(0);
		 }
         
                 return d0;



}

/**************************************************************************************************************/
void pad_vv(int nx,int ny,int nz,int npd,float ***ee)
{
		 int i,j,k;

                   

/*****pad left side                    */
          for(j=npd;j<npd+ny;j++)
            for(k=npd;k<nz+npd;k++)
	    {
              for(i=0;i<npd;i++)
              { 
                 ee[i][j][k]=ee[npd][j][k];
              }
	    }
       
/*****pad right side                    */
          for(j=npd;j<npd+ny;j++)
             for(k=npd;k<nz+npd;k++)
	     {
                for(i=nx+npd;i<nx+2*npd;i++)
                 {
                  ee[i][j][k]=ee[nx+npd-1][j][k];
                 }
	     }
/*****pad front side                    */
          for(i=0;i<npd*2+nx;i++)
             for(k=npd;k<nz+npd;k++)
	     {
                for(j=0;j<npd;j++)
                 {
                  ee[i][j][k]=ee[i][npd][k];
                 }
	     }
/*****pad back side                    */
          for(i=0;i<npd*2+nx;i++)
             for(k=npd;k<nz+npd;k++)
	     {
                for(j=ny+npd;j<ny+2*npd;j++)
                 {
                  ee[i][j][k]=ee[i][ny+npd-1][k];
                 }
	     }

/*****pad upper side                    */
            for(i=0;i<nx+2*npd;i++)
	    {
              for(j=0;j<ny+2*npd;j++)
              {
                for(k=0;k<npd;k++)
                   ee[i][j][k]=ee[i][j][npd];
              }
	    }
/*****lower side                        */
            for(i=0;i<nx+2*npd;i++)
	    {
              for(j=0;j<ny+2*npd;j++)
              {
                for(k=npd+nz;k<npd*2+nz;k++)
                   ee[i][j][k]=ee[i][j][npd+nz-1];
              }
	    }
		
		
}
/**************************************************************************************************************/
void read_file(char FN1[],int nx,int ny,int nz,float ***vv,float ***rho0,int npd)
{

		 int i,j,k;
		

                 

		 FILE *fp1;
		if(( fp1=fopen(FN1,"rb"))==NULL){exit(0);printf("The %s open error!",FN1);};
                 for(j=npd;j<ny+npd;j++)
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(k=npd;k<nz+npd;k++)
			 {
				 fread(&vv[i][j][k],4,1,fp1);

			 }
		 }
                 for(j=npd;j<ny+npd;j++)
		 for(i=npd;i<nx+npd;i++)
		 {
			 for(k=npd;k<nz+npd;k++)
			 {
				 rho0[i][j][k]=1.0;
			 }
		 }
		 fclose(fp1);
}

/**************************************************************************************************************/
void current_shot(float ***vp0,float ***rho0,float ***vp,float ***rho,int nx,int ny,int nz,int npd,
                  int vnx,int vny,int vnz,int is)
{
      
            
                         int ivstart,ivend;
			 int i,ix,iy,iz;
                         is=1;

                        

                     //    ivstart=1+(is-1)*ds_sxd;
		    //	 ivend=nx+(is-1)*ds_sxd;
                         ivstart=1;
		    	 ivend=nx;
			 if(ivstart<=0)
			 {
				 printf("ivstart less than zero \n");
				 exit(0);
			 }
			 if(ivend>vnx)
			 {
				 printf("ivend great than Vnx \n");
				 exit(0);
			 }
                         for(iy=npd;iy<ny+npd;iy++)
			 for(ix=npd;ix<nx+npd;ix++)
			 {
				 for(iz=npd;iz<nz+npd;iz++)
				 {
				  vp[ix][iy][iz]=vp0[ix][iy][iz];
                                  rho[ix][iy][iz]=rho0[ix][iy][iz];
				  
				 }
			 }

}
/**************************************************************************************************************/
void initial_coffe(float dt,float d0,int nx,int ny,int nz,float *coffx1,float *coffx2,float *coffy1,
                   float *coffy2,float *coffz1,float *coffz2,float *acoffx1,float *acoffx2,
                   float *acoffy1,float *acoffy2,float *acoffz1,float *acoffz2,int npd)
{		
		 int i,j,k;
         

		 for(i=0;i<npd;i++)
		 {   
			 coffx1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);
                         coffy1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffy2[i]=coffy1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz1[i]=1/(1+(dt*d0*pow((npd-0.5-i)/npd,2))/2);
			 coffz2[i]=coffz1[i]*(1-(dt*d0*pow((npd-0.5-i)/npd,2))/2);

	
		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 coffx1[i]=1/(1+(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
			 coffx2[i]=coffx1[i]*(1-(dt*d0*pow((0.5+i-nx-npd)/npd,2))/2);
		 }
        
                 for(i=npd+ny;i<ny+2*npd;i++)
		 {
			 coffy1[i]=1/(1+(dt*d0*pow((0.5+i-ny-npd)/npd,2))/2);
			 coffy2[i]=coffy1[i]*(1-(dt*d0*pow((0.5+i-ny-npd)/npd,2))/2);
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
                 for(i=npd;i<npd+ny;i++)
		 {
			 coffy1[i]=1.0;
			 coffy2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 coffz1[i]=1.0;
			 coffz2[i]=1.0;
		 }

		 
		 for(i=0;i<npd;i++)    
		 {    
			 acoffx1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffx2[i]=acoffx1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
                         acoffy1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffy2[i]=acoffy1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz1[i]=1/(1+(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);
			 acoffz2[i]=acoffz1[i]*(1-(dt*d0*pow(((npd-i)*1.0)/npd,2))/2);

		 }

		 for(i=npd+nx;i<nx+2*npd;i++)
		 {
			 acoffx1[i]=1/(1+(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
			 acoffx2[i]=acoffx1[i]*(1-(dt*d0*pow(((1+i-nx-npd)*1.0)/npd,2))/2);
		 }
                 for(i=npd+ny;i<ny+2*npd;i++)
		 {
			 acoffy1[i]=1/(1+(dt*d0*pow(((1+i-ny-npd)*1.0)/npd,2))/2);
			 acoffy2[i]=acoffy1[i]*(1-(dt*d0*pow(((1+i-ny-npd)*1.0)/npd,2))/2);
		 }
		 for(i=npd+nz;i<nz+2*npd;i++)
		 {
			 acoffz1[i]=1/(1+(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
			 acoffz2[i]=acoffz1[i]*(1-(dt*d0*pow(((1+i-nz-npd)*1.0)/npd,2))/2);
		 }

		 for(i=npd;i<npd+nx;i++)
		 {
			 acoffx1[i]=1.0;
			 acoffx2[i]=1.0;
		 }
                 for(i=npd;i<npd+ny;i++)
		 {
			 acoffy1[i]=1.0;
			 acoffy2[i]=1.0;
		 }
		 for(i=npd;i<npd+nz;i++)
		 {
			 acoffz1[i]=1.0;
			 acoffz2[i]=1.0;
		 }
			       
}
/**************************************************************/
float square_sum(int nx,int ny,int nz,float ***a)
{
	int i,j,k;
	float s=0.0;
	for(i=0;i<nx;i++)
	{
	   for(j=0;j<ny;j++)
	   {
              for(k=0;k<nz;k++)
              {
			s=s+a[i][j][k]*a[i][j][k];
              }
	   }
	}
	return s;
}
/**************************************************************/
void read_velfile(int vnx,int vny,int vnz,char FN1[],float ***v)
{
	int i,j,k;
	FILE *fp1;

	if((fp1=fopen(FN1,"rb"))==NULL){exit(0);printf("The %s open error!",FN1);};
        for(j=0;j<vny;j++)
	for(i=0;i<vnx;i++)
	{
		for(k=0;k<vnz;k++)
		{
			fread(&v[i][j][k],4,1,fp1);
		}
	}

	fclose(fp1);
	
}
/**********************************/ 
void window(int vnx,int vny,int vnz,int m,int n,float ***g,float ***g1)
{
      int i,j,k;
      
	for(j=0;j<vny;j++)
	for(i=0;i<vnx;i++)
	{
		for(k=0;k<vnz;k++)
		{
			g1[i][j][k]=0;
		}
	}     
  
        for(j=0;j<vny;j++)
	for(i=0;i<vnx;i++)
	{
		for(k=m;k<n;k++)
		{
			g1[i][j][k]=g[i][j][k];
		}
	}
	
}










