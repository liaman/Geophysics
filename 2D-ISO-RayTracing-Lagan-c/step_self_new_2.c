//##################################################################
//##
//##          STEP2:       Raytrcing_2D  
//##
//##     Ps:  adapt to plane and fluctuate
//##          *Initial code comes from Doc.Qin Ning
//##                
//##                         plane :     2016.4.13   RongTao
//##                       fluctuate:    2016.4.14   RongTao
//##################################################################
#include<stdio.h>
#include<malloc.h>
#include<math.h>
#include<stdlib.h>
#include "hc/cjbsegy.h"
#include "hc/fft.c"
#include "hc/alloc.c"
#include "hc/complex.c"

void main()
{
   void step2_keyin(int *iflag_surface,char fn1[],char fn2[],int *nx,int *nz,
                   float *dx,float *dz,float *s,int *nray0,int *nshot,
                   float *dangle,float *start_angle);
   void step2_readpar(int *iflag_surface,char fn1[],char fn2[],int *nx,int *nz,
                   float *dx,float *dz,float *s,int *nray0,int *nshot,
                   float *dangle,float *start_angle);
   void step2_Raytracing_2D_plane(float *a,float *b,float *c,char fn1[],
                   int nx,int nz,float dx,float dz,float s,float pvel,float svel,
                   int nray,int nshot,int nray0,int ngrid,float dangle,
                   float start_angle); 
   void step2_read_elev(char fn2[],int nx,int *elev);
   void step2_get_maxelev(int *elev,int nx,float *max_elev);
  
   void step2_Raytracing_2D_fluctuate(float *a,float *b,float *c,int *elev,
                   char fn1[],
                   int nx,int nz,float dx,float dz,float s,float pvel,float svel,
                   int nray,int nshot,int nray0,int ngrid,float dangle,
                   float start_angle,float max_elev);


    int svel=1700,pvel=2000;
    long Lmax=9999999;
    float *Buff,*a,*b,*c;
    float max_elev,iflag;
    int iflag_surface,flag;
    int ngrid;
    int nx,nz,nray0,nshot,nray;
    float dx,dz,s,dangle,start_angle;

    int i,j,L1,L2,L3,L4,L5;
   
    char fn1[250],fn2[250];

    Buff=alloc1float(Lmax);   zero1float(Buff,Lmax);
    a=alloc1float(Lmax);      zero1float(a,Lmax);
    b=alloc1float(Lmax);      zero1float(b,Lmax);
    c=alloc1float(Lmax);      zero1float(c,Lmax);
 
    printf("The program is started ! \n");
//===================================
//  Read parameter  *.par file
//===================================
//   parameter instruction
//===================================
//      nray0:ray number of each shot
//      nshot:number of shots
//      ngrid:number of grid
//      nray:total number of ray
//===================================

      printf("Ready to run this program ?\n");
      printf("Make sure <tomo-2.par> is done !\n");
      system("pause");    
 

     printf("IF YOU INPUT THE OPERATING PARAMETERS FROM\n");
     printf(" Key in : 1\n");
     printf("Read par: 2\n");

     //scanf("%d",&flag);
     flag=2;
     if(flag==1)
     {
       step2_keyin(&iflag_surface,fn1,fn2,&nx,&nz,&dx,&dz,&s,
                  &nray0,&nshot,&dangle,&start_angle);
     }else if(flag==2){
       step2_readpar(&iflag_surface,fn1,fn2,&nx,&nz,&dx,&dz,&s,
                  &nray0,&nshot,&dangle,&start_angle);
     }
     else
     { printf("Please input the right flag number!\n");exit(1);}

     nray=nray0*nshot;
     ngrid=(nx-1)*(nz-1);

      printf("*****************************\n");
      printf("**  iflag_surface=%d\n",iflag_surface);
      printf("**  fn1= %s\n",fn1);
      printf("**  fn2= %s\n",fn2);
      printf("**  nx=%d\n",nx);
      printf("**  nz=%d\n",nz);
      printf("**  dx=%f\n",dx);
      printf("**  dz=%f\n",dz);
      printf("**  s=%f\n",s);
      printf("**  nray0=%d\n",nray0);
      printf("**  nshot=%d\n",nshot);
      printf("**  dangle=%f\n",dangle);
      printf("**  start_angle=%f\n",start_angle);
      printf("**  nray=%d\n",nray);
      printf("**  ngrid=%d\n",ngrid);
      printf("*****************************\n\n\n");
      printf("Ready to continue :\n");
      system("pause");

     FILE *fptc;
     fptc=fopen("step2_touched_cell1.txt","w");
     fprintf(fptc,"%8d\n",ngrid);
      for(i=1;i<=ngrid;i++)
         fprintf(fptc,"%5d\n",i);
      fclose(fptc); 

//====================some parameter no use in here 
     L1=1;
     L2=L1+nx*nz;
     L3=L2+ngrid;
     L4=L3+ngrid;
     L5=L4+nx;

//==================================================
//   After get parameter define the arrays
//==================================================
//     elev: "elevation.txt"
//     max_elev: the max elev
//              others is same as plane surface  
//==================================================
    
        int *elev;
        elev=alloc1int(nx+1);
        zero1int(elev,nx+1);

        max_elev=0.0;
//==================================================
//   Change to adapt fluctuate surface
//==================================================
//     if surface is plane: key in iflag_surface=1
//               fluctuate: key in iflag_surface=2
//==================================================

     if(iflag_surface==1)
     {
        step2_Raytracing_2D_plane(a,b,c,fn1,nx,nz,dx,dz,s,pvel,svel,nray,
                                 nshot,nray0,ngrid,dangle,start_angle);
     }else if(iflag_surface==2) 
     {
        step2_read_elev(fn2,nx,elev);

        step2_get_maxelev(elev,nx,&max_elev);

        step2_Raytracing_2D_fluctuate(a,b,c,elev,fn1,nx,nz,dx,dz,s,pvel,svel,nray,
                                 nshot,nray0,ngrid,dangle,start_angle,max_elev);



     }
     
     printf("The program has finished!\n");

     free1float(Buff);
     free1float(a);
     free1float(b);
     free1float(c);

     free1int(elev);


}
//*******************************************************************************//
void step2_Raytracing_2D_fluctuate(float *a,float *b,float *c,int *elev,char fn1[],
                   int nx,int nz,float dx,float dz,float s,float pvel,float svel,
                   int nray,int nshot,int nray0,int ngrid,float dangle,
                   float start_angle,float max_elev)
//===============================================================================//
//===============================================================================//
//==            Subroutine of raytrcing_2D fluctuate surface                   ==//
//===============================================================================//
//==                          parameter instruction                            ==//
//===============================================================================//
//==                                        fn11:step2_raypath_x_z1.txt         ==//
//==                                        fn22:step2_igrid1.txt               ==//
//==                                        fn33:step2_length1.txt              ==//
//==                                        fn44:step2_row_ngrid_per_ray_add1.txt==//
//==                                        fn55:step2_kk_ngrid_all_ray1.txt     ==//
//==                                        fn66:step2_nj_raynum_per_grid1.txt   ==//
//==                                        vel_1:velocity model               ==//
//===============================================================================//
//===============================================================================//
{
  void step2_determin_shotpoint(float *p_x0,float *p_z0,int nshot,float dx,
                                float dz,float *dip);
  void step2_cal_gridpoint(int *ip_lux,int *ip_ldx,int *ip_rux,int *ip_rdx,
                           int *ip_luz,int *ip_ldz,int *ip_ruz,int *ip_rdz,
                           float dx,float dz,float p_x,float p_z);
  void step2_cal_gridvel(float *v0,float *l_x,float *l_z,float **vel_1,
                         int nx,int nz,float dx,float dz,
                         int ip_lux,int ip_ldx,int ip_rux,int ip_rdx,
                         int ip_luz,int ip_ldz,int ip_ruz,int ip_rdz);
  void step2_cal_path(float p_x,float p_z,float s,
                      float n_x,float n_z,float l_x,float l_z,
                      float *p_xend,float *p_zend,float *n_xnew,float *n_znew,
                      float *time,float v0);
  void step2_determin_gridnum(float p_x,float p_z,int nx,float dx,float dz,
                              int *igrid);


//=========================================same an plane
     float pai=3.141592653;
     int nout=1;

     float p_x,p_z,n_x,n_z,l_x,l_z,n_x0,n_z0;
     float p_xend,p_zend;
     float n_xnew,n_znew;

     float *p_x0,*p_z0;
     float **vel_1;
     float *length_g;
     float *nj;
     float *dip;

     int ix,iz,i,j,k;
     int ip_lux,ip_ldx,ip_rux,ip_rdx;
     int ip_luz,ip_ldz,ip_ruz,ip_rdz;
     int istep,igrid,iray,irayon,ishot;
     int *row;

     float angle;
     float beta,sita;

     int kk_all,mnmn;

     float anglefircor,time,v0;
//=========================================add parameter
     int ixl,ixr;

     float zzl,zzr,slope,zb;
     float rmin;
//===============================================


     char fn11[250]={"step2_raypath_x_z1.txt"};
     char fn22[250]={"step2_igrid1.txt"};
     char fn33[250]={"step2_length1.txt"};
     char fn44[250]={"step2_row_ngrid_per_ray_add1.txt"};
     char fn55[250]={"step2_kk_ngrid_all_ray1.txt"};
     char fn66[250]={"step2_nj_raynum_per_grid1.txt"};

     
     p_x0=alloc1float(nshot+1);        zero1float(p_x0,nshot+1);
     p_z0=alloc1float(nshot+1);        zero1float(p_z0,nshot+1);
     vel_1=alloc2float(nz+1,nx+1);     zero2float(vel_1,nz+1,nx+1);
     length_g=alloc1float(ngrid+1);    zero1float(length_g,ngrid+1);
     nj=alloc1float(ngrid+1);          zero1float(nj,ngrid+1);
     dip=alloc1float(nshot+1);         zero1float(dip,nshot+1);

     row=alloc1int(nray+1);            zero1int(row,nray+1);

     printf("********************************************\n");
     printf("** The Raytracing_2D program has started! **\n");
     printf("********************************************\n");
     printf("** nx=%d\n",nx);
     printf("** nz=%d\n",nz);
     printf("** ngrid3=%d\n",ngrid);
     printf("********************************************\n\n\n");
     
     FILE *fpvel,*fpray,*fpcol,*fplen,*fprow,*fpkk,*fpnj;

     if((fpvel=fopen(fn1,"rb"))==NULL)
     { printf("Open file error ! <%s>\n",fn1);exit(1);}

     fpray=fopen(fn11,"w");
     fpcol=fopen(fn22,"w");
     fplen=fopen(fn33,"w");
     fprow=fopen(fn44,"w");
      fpkk=fopen(fn55,"w");
      fpnj=fopen(fn66,"w");

     FILE *fp_i_len;
     fp_i_len=fopen("step2_igrid_length1.txt","w");
     
     printf("nz=%d\n",nz);

     for(ix=1;ix<=nx;ix++)
     {
        for(iz=1;iz<=nz;iz++)
        {
           fread(&vel_1[ix][iz],4L,1,fpvel);
        }
        //printf("ix=%d,iz=%d,vel_1=%f\n",ix,iz,vel_1[ix][iz]);
     }

     for(ix=1;ix<=nx;ix++)
        for(iz=1;iz<=nz;iz++)
            vel_1[ix][iz]=vel_1[ix][iz];

     for(igrid=1;igrid<=ngrid;igrid++)
     {
        nj[igrid]=0.0;
        //printf("igrid = %d/%d\n",igrid,ngrid);
     }
     printf("===================================\n");

     step2_determin_shotpoint(p_x0,p_z0,nshot,dx,dz,dip);

     kk_all=0;
     row[0]=0;
//=======================================
//      Cycle bu number of ray
//=======================================
//       parameter instruction
//=======================================
//       length_g:length of ray
//       ishot:shot number
//       dip:slope of CDP points
//       sita:exit angle
//       beta:angel of incidence
//       p_x:abscissa of CDP points
//       p_z:ordinate of CDP points
//       angle:sita add slope angle of
//             stratum
//       n_x:azimuth of x direction
//       n_z:azimuth of z direction
//======================================= 
    for(iray=1;iray<=nray;iray++)
    {
       if((iray%100)==0)
          printf("iray = %d / %d \n",iray,nray);

       for(igrid=1;igrid<=ngrid;igrid++)
           length_g[igrid]=0.0;

       irayon=iray%nray0;
  
       if(irayon==0)
       {
          ishot=(iray-irayon)/nray0;
       /*   beta=start_angle+(nray0-1)*dangle:
          sita=asin(sin(beta*pai/180.0)*svel/pvel)*180.0/pai;
          angle=sita+atan(dip[ishot])*180.0/pai;    */    //** hu code
          anglefircor=start_angle+atan(dip[ishot])*180.0/pai;
          angle=anglefircor+(nray0-1)*dangle;            //** qin code

       }else
       {
          ishot=(iray-irayon)/nray0+1;
       /*   beta=start_angle+(irayon-1)*dangle:
          sita=asin(sin(beta*pai/180.0)*svel/pvel)*180.0/pai;
          angle=sita+atan(dip[ishot])*180.0/pai;    */    //** hu code  
          anglefircor=start_angle+atan(dip[ishot])*180.0/pai;
          angle=anglefircor+(irayon-1)*dangle;            //** qin code 
       }
       
       p_x=p_x0[ishot];
       p_z=p_z0[ishot];

       istep=0;

       angle=angle*pai/180.0;
 
       n_x0=cos(angle);
       n_z0=sin(angle);
       n_x=n_x0;
       n_z=n_z0;

       time=0.0;
//for goto:
loop_100:

       do
       {
         istep+=1;
         if((istep%nout)==0)
           fprintf(fpray,"%f  %f\n",p_x,p_z);

         step2_cal_gridpoint(&ip_lux,&ip_ldx,&ip_rux,&ip_rdx,
                             &ip_luz,&ip_ldz,&ip_ruz,&ip_rdz,
                              dx,dz,p_x,p_z);

         step2_cal_gridvel(&v0,&l_x,&l_z,vel_1,nx,nz,dx,dz,
                           ip_lux,ip_ldx,ip_rux,ip_rdx,
                           ip_luz,ip_ldz,ip_ruz,ip_rdz);

         step2_cal_path(p_x,p_z,s,n_x,n_z,l_x,l_z,
                       &p_xend,&p_zend,&n_xnew,&n_znew,&time,v0);
    
         step2_determin_gridnum(p_x,p_z,nx,dx,dz,&igrid);
         
         length_g[igrid]+=s;
         if(igrid>ngrid) 
         {
            printf("Error : igrid > ngrid \n");
            printf("You should break the program ! \n");
            system("pause");
         }
              p_x=p_xend;
	      p_z=p_zend;
	      n_x=n_xnew;
	      n_z=n_znew;
       
         if(p_z<=((max_elev-1)*dz*1.0))
             break;

         //===========================================
         //=       < Judge if out of boundary >      =
         //=  Caculate node(p_xend,p_zend) is out of = 
         //=        boundary condition               =  
         //===========================================
       }while((p_xend>0.0)&&(p_xend<((nx-1)*dx))&&
               (p_zend>0.0)&&(p_zend<(nz-1)*dz));
       

//======================================= add for fluctuate begin
       if(p_z<=((max_elev-1)*dz*1.0))
       {
          if((p_xend>0.0)&&(p_xend<((nx-1)*dx))&&
               (p_zend>0.0)&&(p_zend<(nz-1)*dz))
          {
              ixl=(int)(p_x/dx)+1;
	      ixr=ixl+1;
	      zzl=(elev[ixl]-1)*dz;
	      zzr=(elev[ixr]-1)*dz;
	      slope=(zzr-zzl)/dx;
	      zb=zzl-(zzr-zzl)*(ixl-1);
              
              if(p_z>(slope*p_x+zb))
              {
                 rmin=fabs(p_z-slope*p_x-zb)/sqrt(slope*slope+1);
                 
                 if(rmin<0.5*s)
                   printf("Reach the surface !\n");
                 else goto loop_100;
              }
              else  printf("Stop the surface !\n");
          }
       }
//======================================= add for fluctuate  end 

       k=0;
       for(igrid=1;igrid<=ngrid;igrid++)
       {
          if(length_g[igrid]!=0.0)
          {
             nj[igrid]=nj[igrid]+1.0;
             k+=1;
             kk_all+=1;
             fprintf(fpcol,"%d\n",igrid);                      // "igrid.txt"
             fprintf(fplen,"%e\n",length_g[igrid]);            // "length.txt"
             fprintf(fp_i_len,"%5d  %lf\n",igrid,length_g[igrid]); // "igrid_length1.txt"
          }
       }
       mnmn=-1;
       fprintf(fp_i_len,"%5d\n",mnmn);
       row[iray]=row[iray-1]+k;
       fprintf(fprow,"%5d\n",row[iray]);
       //***********************//
       //** end of loop iray  **//   
       //***********************//
     }//loop iray end

     fclose(fpray);
     fclose(fpcol);
     fclose(fplen);
     fclose(fprow);
     fclose(fp_i_len);

     printf("kk_all = %d\n",kk_all);
     fprintf(fpkk,"%5d\n",kk_all);
     
     fclose(fpkk);

     for(igrid=1;igrid<=ngrid;igrid++)
        fprintf(fpnj,"%e\n",nj[igrid]);

     fclose(fpnj);
     fclose(fpvel);

     printf("The Raytracing_2D program has finished!\n");

     free1float(p_x0);
     free1float(p_z0);
     free2float(vel_1);
     free1float(length_g);
     free1float(nj);
     free1float(dip);

     free1int(row);

}
//*******************************************************************************//
void step2_get_maxelev(int *elev,int nx,float *max_elev)
{
     int ix;
     *max_elev=-999999;
     for(ix=1;ix<=nx;ix++)
     {
        if(*max_elev<elev[ix])
           *max_elev=elev[ix];
     }
}
//*******************************************************************************//
void step2_read_elev(char fn2[],int nx,int *elev)
{
     int ix;
 
     FILE *fpelev;
     if((fpelev=fopen(fn2,"r"))==NULL)
     { printf("Open file error ! <%s>\n",fn2);exit(1);}
     for(ix=1;ix<=nx;ix++)
         fscanf(fpelev,"%d\n",&elev[ix]);

     fclose(fpelev);
}

void step2_Raytracing_2D_plane(float *a,float *b,float *c,char fn1[],int nx,
                   int nz,float dx,float dz,float s,float pvel,float svel,
                   int nray,int nshot,int nray0,int ngrid,float dangle,
                   float start_angle)
//*******************************************************************************//
//===============================================================================//
//===============================================================================//
//==              Subroutine of raytrcing_2D plane surface                     ==//
//===============================================================================//
//==                          parameter instruction                            ==//
//===============================================================================//
//==                                        fn11:step2_raypath_x_z1.txt         ==//
//==                                        fn22:step2_igrid1.txt               ==//
//==                                        fn33:step2_length1.txt              ==//
//==                                        fn44:step2_row_ngrid_per_ray_add1.txt==//
//==                                        fn55:step2_kk_ngrid_all_ray1.txt     ==//
//==                                        fn66:step2_nj_raynum_per_grid1.txt  ==//
//==                                        vel_1:velocity model               ==//
//===============================================================================//
//===============================================================================//
{
  void step2_determin_shotpoint(float *p_x0,float *p_z0,int nshot,float dx,
                                float dz,float *dip);
  void step2_cal_gridpoint(int *ip_lux,int *ip_ldx,int *ip_rux,int *ip_rdx,
                           int *ip_luz,int *ip_ldz,int *ip_ruz,int *ip_rdz,
                           float dx,float dz,float p_x,float p_z);
  void step2_cal_gridvel(float *v0,float *l_x,float *l_z,float **vel_1,
                         int nx,int nz,float dx,float dz,
                         int ip_lux,int ip_ldx,int ip_rux,int ip_rdx,
                         int ip_luz,int ip_ldz,int ip_ruz,int ip_rdz);
  void step2_cal_path(float p_x,float p_z,float s,
                      float n_x,float n_z,float l_x,float l_z,
                      float *p_xend,float *p_zend,float *n_xnew,float *n_znew,
                      float *time,float v0);
  void step2_determin_gridnum(float p_x,float p_z,int nx,float dx,float dz,
                              int *igrid);


     float pai=3.141592653;
     int nout=1;

     float p_x,p_z,n_x,n_z,l_x,l_z,n_x0,n_z0;
     float p_xend,p_zend;
     float n_xnew,n_znew;

     float *p_x0,*p_z0;
     float **vel_1;
     float *length_g;
     float *nj;
     float *dip;

     int ix,iz,i,j,k;
     int ip_lux,ip_ldx,ip_rux,ip_rdx;
     int ip_luz,ip_ldz,ip_ruz,ip_rdz;
     int istep,igrid,iray,irayon,ishot;
     int *row;

     float angle;
     float beta,sita;

     int kk_all,mnmn;

     float anglefircor,time,time_ray,v0;


     char fn11[250]={"step2_raypath_x_z1.txt"};
     char fn22[250]={"step2_igrid1.txt"};
     char fn33[250]={"step2_length1.txt"};
     char fn44[250]={"step2_row_ngrid_per_ray_add1.txt"};
     char fn55[250]={"step2_kk_ngrid_all_ray1.txt"};
     char fn66[250]={"step2_nj_raynum_per_grid1.txt"};
     char fn77[250]={"step2_raypath_x_z_t1.txt"};

     p_x0=alloc1float(nshot+1);        zero1float(p_x0,nshot+1);
     p_z0=alloc1float(nshot+1);        zero1float(p_z0,nshot+1);
     vel_1=alloc2float(nz+1,nx+1);     zero2float(vel_1,nz+1,nx+1);
     length_g=alloc1float(ngrid+1);    zero1float(length_g,ngrid+1);
     nj=alloc1float(ngrid+1);          zero1float(nj,ngrid+1);
     dip=alloc1float(nshot+1);         zero1float(dip,nshot+1);

     row=alloc1int(nray+1);            zero1int(row,nray+1);

     printf("********************************************\n");
     printf("** The Raytracing_2D program has started! **\n");
     printf("********************************************\n");
     printf("** nx=%d\n",nx);
     printf("** nz=%d\n",nz);
     printf("** ngrid3=%d\n",ngrid);
     printf("********************************************\n\n\n");
     
     FILE *fpvel,*fpray,*fpcol,*fplen,*fprow,*fpkk,*fpnj,*fptime;

     if((fpvel=fopen(fn1,"rb"))==NULL)
     { printf("Open file error ! <%s>\n",fn1);exit(1);}

     fpray=fopen(fn11,"w");
     fpcol=fopen(fn22,"w");
     fplen=fopen(fn33,"w");
     fprow=fopen(fn44,"w");
      fpkk=fopen(fn55,"w");
      fpnj=fopen(fn66,"w");
    fptime=fopen(fn77,"w");

     FILE *fp_i_len;
     fp_i_len=fopen("step2_igrid_length1.txt","w");
     
     printf("nz=%d\n",nz);

     for(ix=1;ix<=nx;ix++)
     {
        for(iz=1;iz<=nz;iz++)
        {
           fread(&vel_1[ix][iz],4L,1,fpvel);
        }
        printf("ix=%d,iz=%d,vel_1=%f\n",ix,iz,vel_1[ix][iz]);
     }

     for(ix=1;ix<=nx;ix++)
        for(iz=1;iz<=nz;iz++)
            vel_1[ix][iz]=vel_1[ix][iz];

     for(igrid=1;igrid<=ngrid;igrid++)
     {
        nj[igrid]=0.0;
        printf("igrid = %d/%d\n",igrid,ngrid);
     }
     printf("===================================\n");

     step2_determin_shotpoint(p_x0,p_z0,nshot,dx,dz,dip);

     kk_all=0;
     row[0]=0;
//=======================================
//      Cycle bu number of ray
//=======================================
//       parameter instruction
//=======================================
//       length_g:length of ray
//       ishot:shot number
//       dip:slope of CDP points
//       sita:exit angle
//       beta:angel of incidence
//       p_x:abscissa of CDP points
//       p_z:ordinate of CDP points
//       angle:sita add slope angle of
//             stratum
//       n_x:azimuth of x direction
//       n_z:azimuth of z direction
//======================================= 
    for(iray=1;iray<=nray;iray++)
    {
       if((iray%100)==0)
          printf("iray = %d / %d \n",iray,nray);

       for(igrid=1;igrid<=ngrid;igrid++)
           length_g[igrid]=0.0;

       irayon=iray%nray0;
  
       if(irayon==0)
       {
          ishot=(iray-irayon)/nray0;
       /*   beta=start_angle+(nray0-1)*dangle:
          sita=asin(sin(beta*pai/180.0)*svel/pvel)*180.0/pai;
          angle=sita+atan(dip[ishot])*180.0/pai;    */    //** hu code
          anglefircor=start_angle+atan(dip[ishot])*180.0/pai;
          angle=anglefircor+(nray0-1)*dangle;            //** qin code

       }else
       {
          ishot=(iray-irayon)/nray0+1;
       /*   beta=start_angle+(irayon-1)*dangle:
          sita=asin(sin(beta*pai/180.0)*svel/pvel)*180.0/pai;
          angle=sita+atan(dip[ishot])*180.0/pai;    */    //** hu code  
          anglefircor=start_angle+atan(dip[ishot])*180.0/pai;
          angle=anglefircor+(irayon-1)*dangle;            //** qin code 
       }
       
       p_x=p_x0[ishot];
       p_z=p_z0[ishot];

       istep=0;

       angle=angle*pai/180.0;
 
       n_x0=cos(angle);
       n_z0=sin(angle);
       n_x=n_x0;
       n_z=n_z0;

       time=0.0;
       time_ray=0.0;

       do
       {
         istep+=1;
         if((istep%nout)==0)
         {
           fprintf(fpray,"%f  %f\n",p_x,p_z);
           fprintf(fptime,"%f  %f  %f\n",p_x,p_z,time_ray);
         }

         step2_cal_gridpoint(&ip_lux,&ip_ldx,&ip_rux,&ip_rdx,
                             &ip_luz,&ip_ldz,&ip_ruz,&ip_rdz,
                              dx,dz,p_x,p_z);

         step2_cal_gridvel(&v0,&l_x,&l_z,vel_1,nx,nz,dx,dz,
                           ip_lux,ip_ldx,ip_rux,ip_rdx,
                           ip_luz,ip_ldz,ip_ruz,ip_rdz);

         step2_cal_path(p_x,p_z,s,n_x,n_z,l_x,l_z,
                       &p_xend,&p_zend,&n_xnew,&n_znew,&time,v0);
    
         step2_determin_gridnum(p_x,p_z,nx,dx,dz,&igrid);
         
         time_ray+=time;
         length_g[igrid]+=s;
         if(igrid>ngrid) 
         {
            printf("Error : igrid > ngrid \n");
            exit(1);
         }
              p_x=p_xend;
	      p_z=p_zend;
	      n_x=n_xnew;
	      n_z=n_znew;
         //===========================================
         //=       < Judge if out of boundary >      =
         //=  Caculate node(p_xend,p_zend) is out of = 
         //=        boundary condition               =  
         //===========================================
       }while((p_xend>0.0)&&(p_xend<((nx-1)*dx))&&
               (p_zend>0.0)&&(p_zend<(nz-1)*dz));
       
       k=0;
       for(igrid=1;igrid<=ngrid;igrid++)
       {
          if(length_g[igrid]!=0.0)
          {
             nj[igrid]=nj[igrid]+1.0;
             k+=1;
             kk_all+=1;
             fprintf(fpcol,"%d\n",igrid);                      // "igrid.txt"
             fprintf(fplen,"%e\n",length_g[igrid]);            // "length.txt"
             fprintf(fp_i_len,"%5d  %lf\n",igrid,length_g[igrid]); // "igrid_length1.txt"
          }
       }
       mnmn=-1;
       fprintf(fp_i_len,"%5d\n",mnmn);
       row[iray]=row[iray-1]+k;
       fprintf(fprow,"%5d\n",row[iray]);
       //***********************//
       //** end of loop iray  **//   
       //***********************//
     }//loop iray end

     fclose(fpray);
     fclose(fpcol);
     fclose(fplen);
     fclose(fprow);
     fclose(fp_i_len);
     fclose(fptime);

     printf("kk_all = %d\n",kk_all);
     fprintf(fpkk,"%5d\n",kk_all);
     
     fclose(fpkk);

     for(igrid=1;igrid<=ngrid;igrid++)
        fprintf(fpnj,"%e\n",nj[igrid]);

     fclose(fpnj);
     fclose(fpvel);

     printf("The Raytracing_2D program has finished!\n");

}

//********************************************************************************//
void step2_determin_gridnum(float p_x,float p_z,int nx,float dx,float dz,
                              int *igrid)
{
  void determin_igrid(int ix,int iz,int nx,int *igrid);

    int igridxnum,igridznum;

    igridxnum=(int)(p_x/dx)+1;
    igridznum=(int)(p_z/dz)+1;
 
    determin_igrid(igridxnum,igridznum,nx,igrid);

}
//********************************************************************************//
void determin_igrid(int ix,int iz,int nx,int *igrid)
{
    if((iz==1)&&(ix==1))
       *igrid=1;
    else if((iz==1)&&(ix!=1))
       *igrid=ix-1;
    else if((iz!=1)&&(ix==1))
       *igrid=(iz-2)*(nx-1)+ix;
    else
       *igrid=(iz-2)*(nx-1)+(ix-1);
}
//********************************************************************************//
void step2_cal_path(float p_x,float p_z,float s,float n_x,float n_z,float l_x,float l_z,
                    float *p_xend,float *p_zend,float *n_xnew,float *n_znew,
                    float *time,float v0)
{
     float dotmult_ln,dotmult_ll;

     dotmult_ln=l_x*n_x+l_z*n_z;
     dotmult_ll=l_x*l_x+l_z*l_z;
	   
      *p_xend=p_x+n_x*s*(1+dotmult_ln*0.5*s/v0)-0.5*l_x*s*s/v0
              -n_x*s*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(6.0*v0*v0);

      *p_zend=p_z+n_z*s*(1+dotmult_ln*0.5*s/v0)-0.5*l_z*s*s/v0
              -n_z*s*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(6.0*v0*v0);

      *n_xnew=n_x*(1+dotmult_ln*s/v0)-l_x*s/v0
              -n_x*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(2.0*v0*v0);

      *n_znew=n_z*(1+dotmult_ln*s/v0)-l_z*s/v0
              -n_z*s*s*(dotmult_ll-dotmult_ln*dotmult_ln)/(2.0*v0*v0);

      *time=s/v0*(1+s*s*(dotmult_ll+dotmult_ln*dotmult_ln)
            /(6*v0*v0)-dotmult_ln*s/(2*v0));

}
//********************************************************************************//
void step2_cal_gridvel(float *v0,float *l_x,float *l_z,float **vel,
                         int nx,int nz,float dx,float dz,
                         int ip_lux,int ip_ldx,int ip_rux,int ip_rdx,
                         int ip_luz,int ip_ldz,int ip_ruz,int ip_rdz)
{
     float l_x1,l_x2,l_z1,l_z2;

     l_x1=(vel[ip_rux][ip_ruz]-vel[ip_lux][ip_luz])/dx;
     l_x2=(vel[ip_rdx][ip_rdz]-vel[ip_ldx][ip_ldz])/dx;
     *l_x=(l_x1+l_x2)/2.0;

     l_z1=(vel[ip_ldx][ip_ldz]-vel[ip_lux][ip_luz])/dz;
     l_z2=(vel[ip_rdx][ip_rdz]-vel[ip_rux][ip_ruz])/dz;
     *l_z=(l_z1+l_z2)/2.0;

     *v0=vel[ip_rux][ip_ruz]+vel[ip_lux][ip_luz]+vel[ip_rdx][ip_rdz]
           +vel[ip_ldx][ip_ldz];
     *v0=*v0/4.0; 
}
//********************************************************************************//
void step2_cal_gridpoint(int *ip_lux,int *ip_ldx,int *ip_rux,int *ip_rdx,
                           int *ip_luz,int *ip_ldz,int *ip_ruz,int *ip_rdz,
                           float dx,float dz,float p_x,float p_z)
//#################### use the pointer 
{
          *ip_lux=(int)(p_x/dx)+1;
	  *ip_luz=(int)(p_z/dz)+1;

	  *ip_ldx=*ip_lux;
	  *ip_ldz=*ip_luz+1;

	  *ip_rux=*ip_lux+1;
	  *ip_ruz=*ip_luz;

	  *ip_rdx=*ip_lux+1;
	  *ip_rdz=*ip_luz+1;
}
//********************************************************************************//
void step2_determin_shotpoint(float *p_x0,float *p_z0,int nshot,float dx,
                                float dz,float *dip)
{
      int *shotx,*shotz;
      int ishot;

      shotx=alloc1int(nshot+1);   zero1int(shotx,nshot+1);
      shotz=alloc1int(nshot+1);   zero1int(shotz,nshot+1);
      
      FILE *fp;
      if((fp=fopen("step1_shot_location1.txt","r"))==NULL)
      { printf("Open file error ! <step1_shot_location1.txt>\n");exit(1);}
     
      for(ishot=1;ishot<=nshot;ishot++)
         fscanf(fp,"%d%d%e\n",&shotx[ishot],&shotz[ishot],&dip[ishot]);

      

      for(ishot=1;ishot<=nshot;ishot++)
      {
         p_x0[ishot]=(shotx[ishot]-1)*dx*1.0;
         p_z0[ishot]=(shotz[ishot]-1)*dz*1.0;
      }

      free1int(shotx);
      free1int(shotz);
      //printf("------%d,%d,%e,%f,%f \n",shotx[1],shotz[1],dip[1],p_x0[1],p_z0[1]);
}
//********************************************************************************//
void step2_readpar(int *iflag_surface,char fn1[],char fn2[],int *nx,int *nz,
                   float *dx,float *dz,float *s,int *nray0,int *nshot,
                   float *dangle,float *start_angle)
{
      FILE *fp;
     char cc[250];

     if((fp=fopen("tomo-2.par","r"))==NULL)
      {printf("Open file error !\n");exit(1);}
      
      fscanf(fp,"%s%d\n",cc,iflag_surface);
      fscanf(fp,"%s%s\n",cc,fn1);
      fscanf(fp,"%s%s\n",cc,fn2);
      fscanf(fp,"%s%d\n",cc,nx);
      fscanf(fp,"%s%d\n",cc,nz);
      fscanf(fp,"%s%f\n",cc,dx);
      fscanf(fp,"%s%f\n",cc,dz);
      fscanf(fp,"%s%f\n",cc,s);
      fscanf(fp,"%s%d\n",cc,nray0);
      fscanf(fp,"%s%d\n",cc,nshot);
      fscanf(fp,"%s%f\n",cc,dangle);
      fscanf(fp,"%s%f\n",cc,start_angle);

      fclose(fp);
}
//*******************************************************************************//
void step2_keyin(int *iflag_surface,char fn1[],char fn2[],int *nx,int *nz,
                 float *dx,float *dz,float *s,int *nray0,int *nshot,
                 float *dangle,float *start_angle)
{ 
      FILE *fp;
      fp=fopen("Raytracing_undulation.par","w");

      printf("***************************************************\n");
      printf("*             Parameters Inputing                 *\n");
      printf("***************************************************\n");

        printf("Input the type of the surface : iflag_surface\n");
        printf("If plane, key in 1---If fluctuate, key in 2\n");
        scanf("%d",iflag_surface);
        fprintf(fp,"%d\n",*iflag_surface);

        printf("input the velocity field filename: fn1\n");
        scanf("%s",fn1);
        fprintf(fp,"%s\n",fn1);

        printf("input the elevation filename: fn2\n");
        scanf("%s",fn2);
        fprintf(fp,"%s\n",fn2);

        printf("input the the number of lateral sample_point: nx\n");
        scanf("%d",nx);
        fprintf(fp,"%d\n",*nx);

        printf("input the the number of vertical sample_point: nz\n");
        scanf("%d",nz);
        fprintf(fp,"%d\n",*nz);

        printf("input the the interval of lateral sample_point: dx\n");
        scanf("%f",dx);
        fprintf(fp,"%f\n",*dx);

        printf("input the the interval of vertical sample_point: dz\n");
        scanf("%f",dz);
        fprintf(fp,"%f\n",*dz);

        printf("input the Langan step length: s\n");
        scanf("%f",s);
        fprintf(fp,"%f\n",*s);

        printf("input the ray number of each shot: nray0\n");
        scanf("%d",nray0);
        fprintf(fp,"%d\n",*nray0);

        printf("input the number of shots: nshot\n");
        scanf("%d",nshot);
        fprintf(fp,"%d\n",*nshot);

        printf("input the angle interval: dangle\n");
        scanf("%f",dangle);
        fprintf(fp,"%f\n",*dangle);

        printf("input the initial ray direction: start_angle\n");
        scanf("%f",start_angle);
        fprintf(fp,"%f\n",*start_angle);
     
     fclose(fp);
}
