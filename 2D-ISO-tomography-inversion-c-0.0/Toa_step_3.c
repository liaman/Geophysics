//##################################################################
//##
//##          STEP3:       LSQR    
//##
//##     Ps:  Calculate the x matrix of the AX=B ,by lsqr method
//##          *Initial code comes from Doc.Qin Ning
//##                
//##                                   2016.4.19   RongTao
//##                       
//##################################################################
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
  void step3_keyin(int *nray,int *nx,int *nz,float *derr,float *drr,int *itmax,
                  float *dampx,float *dampz);
  void step3_readpar(int *nray,int *nx,int *nz,float *derr,float *drr,int *itmax,
                  float *dampx,float *dampz);
  void step3_lsqr_method(int nray,int ngrid,int nx,int nz,int n_sum,float derr,
                  float drr,float dampx,float dampz,int itmax);

    long Lmax=9999999;
    float *Buff,*a,*b,*c,*d,*e;
    int i,j,L1,L2,L3,L4,L5,L6;

    int nray,nx,nz;
    int n_sum;
    int ngrid,itmax,flag;

    float derr,drr,dampx,dampz;


    Buff=alloc1float(Lmax);   zero1float(Buff,Lmax);
    a=alloc1float(Lmax);      zero1float(a,Lmax);
    b=alloc1float(Lmax);      zero1float(b,Lmax);
    c=alloc1float(Lmax);      zero1float(c,Lmax);
    d=alloc1float(Lmax);      zero1float(d,Lmax);
    e=alloc1float(Lmax);      zero1float(e,Lmax);

     printf("Ready to run this program ?\n");
     printf("Make sure <tomo-3.par> is done !\n");
     system("pause");

//===================================
//  Read parameter  *.par file
//===================================
//   parameter instruction
//===================================
//      nray:ray number of all shot
//      nx:number of vel-x
//      nz:number of grid
//===================================
     printf("The program is started ! \n");
     printf("IF YOU INPUT THE OPERATING PARAMETERS FROM\n");
     printf(" Key in : 1\n");
     printf("Read par: 2\n");

     //scanf("%d",&flag);
     flag=2;

     if(flag==1)
     {
       step3_keyin(&nray,&nx,&nz,&derr,&drr,&itmax,
                  &dampx,&dampz);
     }else if(flag==2){
       step3_readpar(&nray,&nx,&nz,&derr,&drr,&itmax,
                  &dampx,&dampz);
     }
     else
     { printf("Please input the right flag number!\n");exit(1);}

     ngrid=(nx-1)*(nz-1);

     FILE *fp;
     if((fp=fopen("step2_kk_ngrid_all_ray1.txt","r"))==NULL)
      {printf("Open file error ! <%s>\n","step2_kk_ngrid_all_ray1.txt");exit(1);}

     fscanf(fp,"%d\n",&n_sum);
     fclose(fp);

      printf("*****************************\n");
      printf("**  derr= %f\n",derr);
      printf("**  drr= %f\n",drr);
      printf("**  nray=%d\n",nray);
      printf("**  nx=%d\n",nx);
      printf("**  nz=%d\n",nz);
      printf("**  itmax=%d\n",itmax);
      printf("**  dampx=%f\n",dampx);
      printf("**  dampz=%f\n",dampz);
      printf("**  ngrid=%d\n",ngrid);
      printf("**  n_sum=%d\n",n_sum);
      printf("*****************************\n");
      printf("Ready to continue :\n\n\n");
      system("pause");

      L1=1;
      L2=L1+nx*nz;
      L3=L2+(nray+2*ngrid);
      L4=L3+(nray+2*ngrid);
      L5=L4+(n_sum+4*ngrid);
      L6=L5+(n_sum+4*ngrid);

      printf("*****************************\n");
      printf("**   The LSQR is begin !   **\n");
      printf("*****************************\n\n");

      step3_lsqr_method(nray,ngrid,nx,nz,n_sum,derr,drr,dampx,dampz,itmax);

      printf("*****************************\n");
      printf("**  The program is over !  **\n");
      printf("*****************************\n\n");

}

//*****************************************************************************//
//============================================================================
//                LSQR method
//============================================================================
//   Parameter instructions:
//============================================================================
//                              vel:velocity
//                              delta_t:residual of time
//                              nray:number of rays
//                              row,col:parameter of data  compression
//                             length:length of rays in each grid
//============================================================================
void step3_lsqr_method(int nray,int ngrid,int nx,int nz,int n_sum,float derr,
                  float drr,float dampx,float dampz,int itmax)
{
  void step3_regularization(float *row,int *col,float *length,int nray,
                  int ngrid,int n_sum,int nx,float dampx,float dampz); 
  void step3_lsqr_tomo(float *delta_t,float *row,int *col,float *length,
                  float *delta_s,int nray,int ngrid,int n_sum,
                  float derr,float drr,int itmax);
  void step3_determin_igrid(int ix,int iz,int nx,int *igrid);


      float **vel,**s_update,*delta_t,*row;
      float *length,*delta_s;
      int *col;

      int i,j,ix,iz,igrid;

      char fn1[250]={"step1_velmodel1.dat"}; 
      char fn2[250]={"step1_dt1.txt"}; 
      char fn3[250]={"step2_row_ngrid_per_ray_add1.txt"}; 
      char fn4[250]={"step2_igrid1.txt"}; 
      char fn5[250]={"step2_length1.txt"}; 
      char fn6[250]={"step3_vel_update1.dat"}; 
      char fn7[250]={"step3_delta_s1.dat"}; 

      FILE *fp1,*fp2,*fp3,*fp4,*fp5,*fp6,*fp7;

   vel=alloc2float(nz+1,nx+1);           zero2float(vel,nz+1,nx+1); 
   s_update=alloc2float(nz+1,nx+1);      zero2float(s_update,nz+1,nx+1);
   delta_t=alloc1float(nray+2*ngrid+1);  zero1float(delta_t,nray+2*ngrid+1);
   row=alloc1float(nray+2*ngrid+1);      zero1float(row,nray+2*ngrid+1);    
   col=alloc1int(n_sum+4*ngrid+1);       zero1int(col,n_sum+4*ngrid+1);
   length=alloc1float(n_sum+4*ngrid+1);  zero1float(length,n_sum+4*ngrid+1);
   delta_s=alloc1float(ngrid+1);         zero1float(delta_s,ngrid+1);

//========== read the 'velmodel.dat' ========== 1
      if((fp1=fopen(fn1,"rb"))==NULL)
      {printf("Open file error < %s > !\n",fn1);exit(1);}
      for(i=1;i<=nx;i++)
         for(j=1;j<=nz;j++)
            fread(&vel[i][j],4L,1,fp1);
      fclose(fp1);

//========== read the 'dt.txt' ========== 2
      if((fp2=fopen(fn2,"r"))==NULL)
      {printf("Open file error < %s > !\n",fn2);exit(1);}
      for(i=1;i<=nray;i++)
         fscanf(fp2,"%f\n",&delta_t[i]);
      fclose(fp2);

//========== read the 'row.txt' ========== 3
      if((fp3=fopen(fn3,"r"))==NULL)
      {printf("Open file error < %s > !\n",fn3);exit(1);}     
      for(i=1;i<=nray;i++)
         fscanf(fp3,"%f\n",&row[i]);
      fclose(fp3);

//========== read the 'col.txt' ========== 4
      if((fp4=fopen(fn4,"r"))==NULL)
      {printf("Open file error < %s > !\n",fn4);exit(1);}     
      for(i=1;i<=n_sum;i++)
         fscanf(fp4,"%d\n",&col[i]);
      fclose(fp4);
      
//========== read the 'length.txt' ========== 5
      if((fp5=fopen(fn5,"r"))==NULL)
      {printf("Open file error < %s > !\n",fn5);exit(1);}     
      for(i=1;i<=n_sum;i++)
         fscanf(fp5,"%f\n",&length[i]);
      fclose(fp5);

      //printf("vel=%f\n,delta_t=%e\n,row=%f\n,col=%f\n,length=%e\n",
      //        vel[1][1],delta_t[nray],row[nray],col[n_sum],length[n_sum]);

      step3_regularization(row,col,length,nray,ngrid,n_sum,nx,dampx,dampz);

      step3_lsqr_tomo(delta_t,row,col,length,delta_s,nray,ngrid,n_sum,
                      derr,drr,itmax);

      for(ix=1;ix<=nx;ix++)
      {
         for(iz=1;iz<=nz;iz++)
         {
            step3_determin_igrid(ix,iz,nx,&igrid);
            s_update[ix][iz]=delta_s[igrid];
            vel[ix][iz]=1000.0/(1000.0/vel[ix][iz]+delta_s[igrid]);
         }
      }

      fp6=fopen(fn6,"wb");
      for(ix=1;ix<=nx;ix++)
         for(iz=1;iz<=nz;iz++)
            fwrite(&vel[ix][iz],4L,1,fp6);
      fclose(fp6);

      fp7=fopen(fn7,"wb");
      for(ix=1;ix<=nx;ix++)
         for(iz=1;iz<=nz;iz++)
            fwrite(&s_update[ix][iz],4L,1,fp7);
      fclose(fp7);
 
      printf("nx=%d\nnz=%d\n",nx,nz);

}

//*****************************************************************************//
void step3_determin_igrid(int ix,int iz,int nx,int *igrid)
{
     if(iz==1&&ix==1)
        *igrid=1;
     else if(iz==1&&ix!=1)
        *igrid=ix-1;
     else if(ix==1&&iz!=1)
        *igrid=(iz-2)*(nx-1)+ix;
     else
        *igrid=(iz-2)*(nx-1)+(ix-1);
}
//*****************************************************************************//
void step3_lsqr_tomo(float *delta_t,float *row,int *col,float *length,
                  float *delta_s,int nray,int ngrid,int n_sum,
                  float derr,float drr,int itmax)
{
   void step3_normal(int num,float *xx,float *xx_norm);
   void step3_get_vv(float *vv,float *row,int  *col,float *length,
                     float *delta_t,int ngrid,int nray,int n_sum);
   void step3_get_deltat(float *delta_t,float *row,int *col,float *length,
                     float *vv,int ngrid,int nray,int n_sum);

     int iray,igrid,iterate_num,i,j;
     float beta,alpha,gama;
     float rhobar,phibar;
     float t1,t2;
     float beta_start;
     float rr_ab0,rr_re0;
     float rr_relative;
     float *vv,*ww;

     vv=alloc1float(ngrid+1);    zero1float(vv,ngrid+1);
     ww=alloc1float(ngrid+1);    zero1float(ww,ngrid+1);

//========================
//      STEP_1
//========================
     zero1float(delta_s,ngrid+1);

     printf("ngrid = %d \n",ngrid);
     printf("nray = %d \n",nray);
     

     step3_normal(nray+2*ngrid,delta_t,&beta);

     step3_get_vv(vv,row,col,length,delta_t,ngrid,nray,n_sum);

     step3_normal(ngrid,vv,&alpha);

     for(igrid=1;igrid<=ngrid;igrid++)
        ww[igrid]=vv[igrid];

     phibar=beta;
     rhobar=alpha;
     beta_start=beta;
     rr_relative=phibar/beta_start;

     printf("*****************************\n");
     printf("** The initial error is:\n");
     printf("** rr_absolute=%e\n",phibar);
     printf("** rr_relative=%e\n",rr_relative);
     printf("*****************************\n");
     printf("** beta = %e \n\n\n",beta);
     //system("pause");

     if(phibar<derr||rr_relative<derr)
     {
        printf("Enough precise!\n");
        goto loop_200;
     }
//========================
//      STEP_2
//========================
     for(iterate_num=1;iterate_num<=itmax;iterate_num++)
     {
        printf("ite_num=%5d/%d\n",iterate_num,itmax);
        rr_ab0=phibar;
        rr_re0=rr_relative;
        
        for(iray=1;iray<=nray+2*ngrid;iray++)
            delta_t[iray]=-alpha*delta_t[iray];

        step3_get_deltat(delta_t,row,col,length,vv,ngrid,nray,n_sum);

        step3_normal(nray+2*ngrid,delta_t,&beta);

        for(igrid=1;igrid<=ngrid;igrid++)
           vv[igrid]=-beta*vv[igrid];

        step3_get_vv(vv,row,col,length,delta_t,ngrid,nray,n_sum);

        step3_normal(ngrid,vv,&alpha);

        gama=sqrt(rhobar*rhobar+beta*beta);
        t1=rhobar*phibar/(gama*gama);
    	t2=-alpha*beta/(gama*gama);
        rhobar=-rhobar*alpha/gama;
    	phibar=phibar*beta/gama;

        for(igrid=1;igrid<=ngrid;igrid++)
        {
           delta_s[igrid]=delta_s[igrid]+ww[igrid]*t1;
    	   ww[igrid]=vv[igrid]+ww[igrid]*t2;
        }
        rr_relative=phibar/beta_start;
 
        if(phibar<derr||rr_relative<derr)
        {
           printf("******** loop end ********\n");
           goto loop_200;
        }
        else if((fabs(rr_ab0-phibar)<drr)
              &&(fabs(rr_re0-rr_relative)<drr))
        {
           printf("******** loop end ********\n");
           goto loop_200;
        }
     }
loop_200:
      printf("*****************************\n");
      printf("** iterate_num=%d \n",iterate_num);
      printf("** rr_absolute=%e \n",phibar);
      printf("** rr_relative=%e \n",rr_relative);
      printf("*****************************\n\n\n");
}
//*****************************************************************************//
void step3_get_deltat(float *delta_t,float *row,int *col,float *length,
                     float *vv,int ngrid,int nray,int n_sum)
{
   void step3_get_av(float *av,float *row,int *col,float *length,
                     float *vv,int ngrid,int nray,int n_sum); 
     float *av;
     int iray;

     av=alloc1float(nray+2*ngrid+1);   zero1float(av,nray+2*ngrid+1);

     step3_get_av(av,row,col,length,vv,ngrid,nray,n_sum);

     for(iray=1;iray<=nray+2*ngrid;iray++)
        delta_t[iray]+=av[iray];
}
//*****************************************************************************//
void step3_get_av(float *av,float *row,int *col,float *length,
                     float *vv,int ngrid,int nray,int n_sum)
{
     int i,iray;

     row[0]=0.0;
   
     for(iray=1;iray<=nray+2*ngrid;iray++)
     {
        av[iray]=0.0;
        for(i=row[iray-1]+1;i<=row[iray];i++)
        {
           av[iray]+=length[i]*vv[col[i]];
        }
     }
}
//*****************************************************************************//
void step3_get_vv(float *vv,float *row,int  *col,float *length,
                     float *delta_t,int ngrid,int nray,int n_sum)
{
   void step3_get_atb(float *atb,float *row,int *col,float *length,
                      float *delta_t,int ngrid,int nray,int n_sum);
     float *atb;
     int igrid;

     atb=alloc1float(ngrid+1);    zero1float(atb,ngrid+1);

     step3_get_atb(atb,row,col,length,delta_t,ngrid,nray,n_sum);

     for(igrid=1;igrid<=ngrid;igrid++)
        vv[igrid]+=atb[igrid];

     free1float(atb);

}
//*****************************************************************************//
void step3_get_atb(float *atb,float *row,int *col,float *length,
                      float *delta_t,int ngrid,int nray,int n_sum)
{
     int igrid,iray,i;

     row[0]=0.0;
     for(iray=1;iray<=nray+2*ngrid;iray++)
        for(i=row[iray-1]+1;i<=row[iray];i++)
            atb[col[i]]+=length[i]*delta_t[iray];
         
}
//*****************************************************************************//
//=====================================
//      CACULATE NORMAL
//=====================================
void step3_normal(int num,float *xx,float *xx_norm)
{
    int inum;

     *xx_norm=0.0;
     for(inum=1;inum<=num;inum++)
        *xx_norm+=xx[inum]*xx[inum];
     *xx_norm=sqrt(*xx_norm);
     for(inum=1;inum<=num;inum++)
        xx[inum]=xx[inum]/(*xx_norm);
}
//*****************************************************************************//
//======================================================================
//==                 Metrix regularization
//======================================================================
void step3_regularization(float *row,int *col,float *length,int nray,
                  int ngrid,int n_sum,int nx,float dampx,float dampz)
{
   void step3_add_x_regular1(float *row,int *col,float *length,
                  int igrid,int nray,int ngrid,int n_sum,float dampx);
   void step3_add_x_regular2(float *row,int *col,float *length,
                  int igrid,int nray,int ngrid,int n_sum,float dampx);
   void step3_add_z_regular1(float *row,int *col,float *length,
              int igrid,int nray,int ngrid,int n_sum,int nx,float dampz);
   void step3_add_z_regular2(float *row,int *col,float *length,
              int igrid,int nray,int ngrid,int n_sum,int nx,float dampz);
   

     int igrid;

//====================================================
//==           lateral one_order derivative
//==             type regularation matrix
//====================================================
     for(igrid=1;igrid<=ngrid;igrid++)
     {
         if((igrid%(nx-1))!=0)
             step3_add_x_regular1(row,col,length,igrid,
                                  nray,ngrid,n_sum,dampx);
         else 
             step3_add_x_regular2(row,col,length,igrid,
                                  nray,ngrid,n_sum,dampx);
     } 
//====================================================
//==           vertical one_order derivative
//==             type regularation matrix
//====================================================
     for(igrid=1;igrid<=ngrid;igrid++)
     {
         if((igrid+nx-1)<=ngrid)
             step3_add_z_regular1(row,col,length,igrid,
                                  nray,ngrid,n_sum,nx,dampz);
         else 
             step3_add_z_regular2(row,col,length,igrid,
                                  nray,ngrid,n_sum,nx,dampz);
     } 
}
//*****************************************************************************//
//========================================
//    Lateral one order derivation
//            regulation
//       mod(igrid,nx-1)!0
//========================================
void step3_add_x_regular1(float *row,int *col,float *length,
                  int igrid,int nray,int ngrid,int n_sum,float dampx)
{
     row[nray+igrid]=row[nray+igrid-1]+2.0;
     col[n_sum+2*igrid-1]=igrid;
     col[n_sum+2*igrid]=igrid+1;
     length[n_sum+2*igrid-1]=dampx;
     length[n_sum+2*igrid]=-dampx;
}
//*****************************************************************************//
//========================================
//    Lateral one order derivation
//            regulation
//       mod(igrid,nx-1)==0
//========================================
void step3_add_x_regular2(float *row,int *col,float *length,
                  int igrid,int nray,int ngrid,int n_sum,float dampx)
{
        row[nray+igrid]=row[nray+igrid-1]+2.0;
	col[n_sum+2*igrid-1]=igrid-1;
	col[n_sum+2*igrid]=igrid;
	length[n_sum+2*igrid-1]=dampx;
	length[n_sum+2*igrid]=-dampx;
}
//*****************************************************************************//
//========================================
//    Vateral one order derivation
//            regulation
//       (igrid+nx-1)<=ngrid
//========================================
void step3_add_z_regular1(float *row,int *col,float *length,
                int igrid,int nray,int ngrid,int n_sum,int nx,float dampz)
{
        row[nray+ngrid+igrid]=row[nray+ngrid+igrid-1]+2.0;
        col[n_sum+2*ngrid+2*igrid-1]=igrid;
        col[n_sum+2*ngrid+2*igrid]=igrid+nx-1;
        length[n_sum+2*ngrid+2*igrid-1]=dampz;
        length[n_sum+2*ngrid+2*igrid]=-dampz;
}
//*****************************************************************************//
//========================================
//    Vateral one order derivation
//            regulation
//       (igrid+nx-1)>ngrid
//========================================
void step3_add_z_regular2(float *row,int *col,float *length,
                int igrid,int nray,int ngrid,int n_sum,int nx,float dampz)
{
        row[nray+ngrid+igrid]=row[nray+ngrid+igrid-1]+2.0;
        col[n_sum+2*ngrid+2*igrid-1]=igrid-nx+1;
        col[n_sum+2*ngrid+2*igrid]=igrid;
        length[n_sum+2*ngrid+2*igrid-1]=dampz;
        length[n_sum+2*ngrid+2*igrid]=-dampz;
}
//*****************************************************************************//
void step3_readpar(int *nray,int *nx,int *nz,float *derr,float *drr,int *itmax,
                  float *dampx,float *dampz)
{
    FILE *fp;
    char cc[250];
     if((fp=fopen("tomo-3.par","r"))==NULL)
      {printf("Open file error !\n");exit(1);}
      
      fscanf(fp,"%s%f\n",cc,derr);
      fscanf(fp,"%s%f\n",cc,drr);
      fscanf(fp,"%s%d\n",cc,nray);
      fscanf(fp,"%s%d\n",cc,nx);
      fscanf(fp,"%s%d\n",cc,nz);
      fscanf(fp,"%s%d\n",cc,itmax);
      fscanf(fp,"%s%f\n",cc,dampx);
      fscanf(fp,"%s%f\n",cc,dampz);

      fclose(fp);
}
//*****************************************************************************//
void step3_keyin(int *nray,int *nx,int *nz,float *derr,float *drr,int *itmax,
                  float *dampx,float *dampz)
{
      FILE *fp;
      fp=fopen("LSQR_regular.par","w");

      printf("***************************************************\n");
      printf("*             Parameters Inputing                 *\n");
      printf("***************************************************\n");


        printf("input error of the absolute or relative rr: derr\n");
        scanf("%f",derr);
        fprintf(fp,"%f\n",*derr);

        printf("input error of the neighboring rr: drr\n");
        scanf("%f",drr);
        fprintf(fp,"%f\n",*drr);

        printf("input ray number: nray\n");
        scanf("%d",nray);
        fprintf(fp,"%d\n",*nray);

        printf("input the the number of lateral sample_point: nx\n");
        scanf("%d",nx);
        fprintf(fp,"%d\n",*nx);

        printf("input the the number of vertical sample_point: nz\n");
        scanf("%d",nz);
        fprintf(fp,"%d\n",*nz);

        printf("input the max iterate number: itmax\n");
        scanf("%d",itmax);
        fprintf(fp,"%d\n",*itmax);

        printf("input lateral regular factor: dampx\n");
        scanf("%f",dampx);
        fprintf(fp,"%f\n",*dampx);

        printf("input vertical regular factor: dampz\n");
        scanf("%f",dampz);
        fprintf(fp,"%f\n",*dampz);
     
     fclose(fp);

}
