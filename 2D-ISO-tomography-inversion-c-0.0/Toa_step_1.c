//##################################################################
//##
//##          STEP1:       GET dt & Velmodel
//##
//##     Ps:  get the dt and velmodel from adcig and migration
//##          *Initial code comes from Doc.Qin Ning
//##
//##                             2016.6.25   RongTao
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
   void step1_keyin(char fnvel[],char fn1[],char fn2[],int *iflag,float *factor,int *sign,
              int *nx,int *nz,float *dx,float *dz,float *svel,float *pvel,
              int *ncdp,int *nangle,float *dangle,float *anglestart,int *xstart,
              int *xend,float *dzmax,int *nmove,int * ntra_per_adcig ,int *nzu,int *nzd);
   void step1_readpar(char fnvel[],char fn1[],char fn2[],int *iflag,float *factor,int *sign,
              int *nx,int *nz,float *dx,float *dz,float *svel,float *pvel,
              int *ncdp,int *nangle,float *dangle,float *anglestart,int *xstart,
              int *xend,float *dzmax,int *nmove,int * ntra_per_adcig ,int *nzu,int *nzd);
   void step1_cal_dt_velmodel(char fnvel[],char fn1[],char fn2[],int iflag,float factor,int sign,
              int nx,int nz,float dx,float dz,float svel,float pvel,
              int ncdp,int nangle,float dangle,float anglestart,int xstart,
              int xend,float dzmax,int nmove,int  ntra_per_adcig ,int nzu,int nzd,
              int ntrace_adcig_use ,int ntrcid_mig_use ,int dshot);


      int num=1, ncdpstep_cig =1, dshot= 10,l1=3000,l2=1000;
    float **valcig,**adcig_data,**valmig,**mig_data,**v;
    float *deltt,*depth,*z,autval[10][200],aut[21];
    float *deltz,angl_select[10];
      int *layer,ind[10];
    float *layerz,gama[10],*layer3,k;
    float *x0,*xx,*dy,*ddy,*h,*s,*ds,*dds;
      int *topz;
      int flag,ntrace_adcig_use ,ntrcid_mig_use ;
     char fnvel[250],fn1[250],fn2[250];
      int i, ntra_per_adcig ,nzu,nzd,xstart,nangle, mid_per_adcig ;
      int iflag,nx,nz,ncdp,xend,sign;
      int nmove;
    float factor,dx,dz,svel,pvel,dangle,anglestart,dzmax;

     printf("Ready to run this program ?\n");
     printf("Make sure <tomo-1.par> is done !\n");
     system("pause");

    printf("IF YOU INPUT THE OPERATING PARAMETERS FROM\n");
    printf(" Key in : 1\n");
    printf("Read par: 2\n");
    

    //scanf("%d",&flag);
    flag=2;

    if(flag==1)
    {
       step1_keyin(fnvel,fn1,fn2,&iflag,&factor,&sign,&nx,&nz,&dx,&dz,
            &svel,&pvel,&ncdp,&nangle,&dangle,&anglestart,
             &xstart,&xend,&dzmax,&nmove,& ntra_per_adcig ,&nzu,&nzd);
    }else if(flag==2){
       step1_readpar(fnvel,fn1,fn2,&iflag,&factor,&sign,&nx,&nz,&dx,&dz,
            &svel,&pvel,&ncdp,&nangle,&dangle,&anglestart,
             &xstart,&xend,&dzmax,&nmove,& ntra_per_adcig ,&nzu,&nzd);
    }
    else
    { printf("Please input the right flag number!\n");exit(1);}


      ntrace_adcig_use =(xend-xstart+1)/num;
      ntrcid_mig_use =ntrace_adcig_use /nangle;

      printf("*****************************\n");
      printf("*  fnvel= %s\n",fnvel);
      printf("*  fn1= %s\n",fn1);
      printf("*  fn2= %s\n",fn2);
      printf("*  iflag=%d\n",iflag);
      printf("*  factor=%f\n",factor);
      printf("*  sign=%d\n",sign);
      printf("*  nx=%d\n",nx);
      printf("*  nz=%d\n",nz);
      printf("*  dx=%f\n",dx);
      printf("*  dz=%f\n",dz);
      printf("*  pvel=%f\n",pvel);
      printf("*  ncdp=%d\n",ncdp);
      printf("*  nangle=%d\n",nangle);
      printf("*  dangle=%f\n",dangle);
      printf("*  anglestart=%f\n",anglestart);
      printf("*  xstart,xend=%d %d\n",xstart,xend);
      printf("*  dzmax=%f\n",dzmax);
      printf("*  nmove=%d\n",nmove);
      printf("*   ntra_per_adcig =%d\n", ntra_per_adcig );
      printf("*  ntrace_adcig_use =%d\n",ntrace_adcig_use );
      printf("*****************************\n\n\n");
      printf("Ready to continue :\n");
      system("pause");

      step1_cal_dt_velmodel(fnvel,fn1,fn2,iflag,factor,sign,nx,nz,dx,dz,
                     svel,pvel,ncdp,nangle,dangle,anglestart,
                     xstart,xend,dzmax,nmove, ntra_per_adcig ,nzu,nzd,
                     ntrace_adcig_use ,ntrcid_mig_use ,dshot);


      printf("The program is finished !\n");




}


//**********************************************************************************//
void step1_cal_dt_velmodel(char fnvel[],char fn1[],char fn2[],int iflag,float factor,int sign,
              int nx,int nz,float dx,float dz,float svel,float pvel,
              int ncdp,int nangle,float dangle,float anglestart,int xstart,
              int xend,float dzmax,int nmove,int  ntra_per_adcig ,int nzu,int nzd,
              int ntrace_adcig_use ,int ntrcid_mig_use ,int dshot)
{
  void step1_espl1(float *x0,float *layerz,int ncdp,float dy1,float dyn,float *xx,
              int i2,float *dy,float *ddy,float *s,float *ds,float *dds,float t,
              float *h);
  void step1_corautopick(int i,int  ntra_per_adcig ,int nzu,int nzd,float *z,float **valcig,
              int *ind,int xstart,int nangle,int  mid_per_adcig ,float dz,float z0,int nx,int nz,
                   int ncdpstep_cig);   

      float pi=3.1415927;
      int num=1, ncdpstep_cig =1;
      float veltop=10,velbottom=10;

      char fn6[250]={"step1_dt1.txt"};
      char fn7[250]={"step1_dz1.txt"};
      char fn8[250]={"step1_z1.txt"};
      
      float **valcig,**adcig_data;
      float *deltt,valmax;
      float **mig_data,**valmig;
      float *depth,valmigmax,*z,**v,**v_initial;
      float *deltz,*angl_select,*layer3;
      int *layer,*ind;
      float *layerz,*gama,x,y,z0,*x0,xx[3],*dy,*ddy,*h,s[3],ds[3],dds[3],dy1,dyn,t;
      int zmin,zmax,nmax;
      int *topz;

      int i,j,k,ix,iz;
      int  mid_per_adcig ,icdp;
 
      float faj,gamave,midd,angle;
      float vsum;

/****************************   Data Initialization ***********************/
          zero1float(xx,3);zero1float(s,3);
          zero1float(ds,3);zero1float(dds,3);

           valcig=alloc2float(nz+1,nx+1);      zero2float(valcig,nz+1,nx+1);
      adcig_data=alloc2float(nz+1,nx+1);      zero2float(adcig_data,nz+1,nx+1);

           deltt=alloc1float(ntrace_adcig_use +1);     zero1float(deltt,ntrace_adcig_use +1);
           depth=alloc1float(ntrace_adcig_use +1);     zero1float(depth,ntrace_adcig_use +1);
           deltz=alloc1float(ntrace_adcig_use +1);     zero1float(deltz,ntrace_adcig_use +1);

        mig_data=alloc2float(nz+1,ncdp+1);    zero2float(mig_data,nz+1,ncdp+1);
          valmig=alloc2float(nz+1,ncdp+1);    zero2float(valmig,nz+1,ncdp+1);
               v=alloc2float(nz+1,ncdp+1);    zero2float(v,nz+1,ncdp+1);
       v_initial=alloc2float(nz+1,ncdp+1);    zero2float(v_initial,nz+1,ncdp+1);

               z=alloc1float( ntra_per_adcig +1);       zero1float(z, ntra_per_adcig +1);       
            angl_select=alloc1float( ntra_per_adcig +1);       zero1float(angl_select, ntra_per_adcig +1);
            gama=alloc1float( ntra_per_adcig +1);       zero1float(gama, ntra_per_adcig +1);
             ind=alloc1int( ntra_per_adcig +1);         zero1int(ind, ntra_per_adcig +1);

          layer3=alloc1float(ncdp+1);       zero1float(layer3,ncdp+1);
          layerz=alloc1float(ncdp+1);       zero1float(layerz,ncdp+1);
              x0=alloc1float(ncdp+1);       zero1float(x0,ncdp+1);
              dy=alloc1float(ncdp+1);       zero1float(dy,ncdp+1);
             ddy=alloc1float(ncdp+1);       zero1float(ddy,ncdp+1);
               h=alloc1float(ncdp+1);       zero1float(h,ncdp+1);
 
           layer=alloc1int(ncdp+1);         zero1int(layer,ncdp+1);
            topz=alloc1int(ncdp+1);         zero1int(topz,ncdp+1);
/***************************************************************************/
//=============================================
//  Read data from simply migrated section
//   and compare each other to get layers
//=============================================
//         parameter instruction
//=============================================
//     ncdp:number of CDP points
//     nz:number of interval samples
//     mig_data:data of migrated image
//     nzu:start of window near layer
//     nzd:end of window near layer
//     layer(i):depth of choosen layer
//     zmax:the most depth location of the layer
//==============================================       
      dy1=0.0;
      dyn=0.0;
      for(i=1;i<=ncdp;i++)
         x0[i]=(i-1)*dx;
      
      xx[1]=2*dx;
      xx[2]=3*dx;

      FILE *fp1;
      if((fp1=fopen(fn1,"rb"))==NULL)
       { printf("Open file error ! <%s>\n",fn1);exit(1);}
      for(ix=1;ix<=ncdp;ix++)
         for(iz=1;iz<=nz;iz++)
            fread(&mig_data[ix][iz],4L,1,fp1);
      fclose(fp1);
      

      for(ix=1;ix<=ncdp;ix++)
         for(iz=1;iz<=nz;iz++)
            valmig[ix][iz]=fabs(mig_data[ix][iz]);
      
      FILE *fpposition;
      if((fpposition=fopen("position.txt","r"))==NULL)
       { printf("Open file error ! <position.txt>\n");exit(1);}

      for(i=1;i<=ncdp;i++)
      {
         fscanf(fpposition,"%d %f\n",&k,&layer3[i]);
         layer[i]=(int)(layer3[i]/dz+0.5);
      }
      fclose(fpposition);      

      for(i=1;i<=ncdp;i++)
      {
         nzu=layer[i]-1;
         nzd=layer[i]+1;
         valmigmax=valmig[i][nzu];
        // layer[i]=nzu;
         for(iz=nzu;iz<=nzd;iz++)
         {
            if(valmig[i][iz]>valmigmax)
            {
               valmigmax=valmig[i][iz];
               layer[i]=iz;
            }
         }
      }

//*****************************************************//
//**                   smooth layer                  **//
//*****************************************************//
 /*    for(i=1;i<=50;i++)
         layer[i]=layer[51];

      for(i=1;i<=50;i++)
         layer[ncdp+i-50]=layer[ncdp-51];     */  

     for(i=1;i<=1+xstart/nangle* ncdpstep_cig;i++)
         layer[i]=layer[1+xstart/nangle* ncdpstep_cig+1+1];

      for(i=1+xend/nangle* ncdpstep_cig;i<=ncdp;i++)
         layer[i]=layer[1+xend/nangle* ncdpstep_cig];      

      layer[1]=(int)((layer[3]+layer[1]+layer[2])/3+0.5);
      layer[2]=layer[1];
      layer[ncdp]=(int)((layer[ncdp-1]+layer[ncdp-2]+layer[ncdp])/3+0.5);
      layer[ncdp-1]=layer[ncdp];

      for(i=3;i<=ncdp-2;i++)
         layer[i]=(int)((layer[i-1]+layer[i]+layer[i+1])/3+0.5);
//******************************************************//

      FILE *fplayer2;
      fplayer2=fopen("layer2.txt","w");

      for(i=1;i<=ncdp;i++)
         fprintf(fplayer2,"%5d %5d\n",i,layer[i]);

      fclose(fplayer2); 

      zmax=layer[1];

      for(i=1;i<=ncdp;i++)
         if(layer[i]>zmax)
            zmax=layer[i];
      
//===================================================
//          elevation data
//===================================================
//       parameter instruction
//===================================================
//      topz():elevation data
//      zmin:the most minimum elevation in
//           the data
//===================================================           
      FILE *fplayer1;
      if((fplayer1=fopen("layer1.txt","r"))==NULL)
       { printf("Open file error ! <layer1.txt>\n");exit(1);}
      for(i=1;i<=ncdp;i++)
      {
         fscanf(fplayer1,"%d %d\n",&j,&topz[i]);
      }
      fclose(fplayer1);

      zmin=topz[i];
      for(i=1;i<=ncdp;i++)
         if(topz[i]<zmin)
            zmin=topz[i];

      if(iflag==2)
      {
        FILE *fpelev;
        fpelev=fopen("step1_elevation.txt","w");
        for(i=1;i<=ncdp;i++)
           fprintf(fpelev,"%5d\n",topz[i]-zmin);
        fclose(fpelev);
      }
//====================================================
//            get velocity field
//====================================================
//          parameter instruction
//====================================================
//        v(i,j):velocity
//        veltop=10
//        velbottom=10
//        pvel:p-wave velocity
//====================================================      
       nmax=zmax-zmin;
      printf("******************************************\n");
      printf("* Please remenber the zmin and zmax!!  **\n");
      printf("******************************************\n");
      printf("* zmax=%d,zmin=%d,\n",zmax,zmin);
      printf("* dimension of velmodel is %d*%d,\n",ncdp,nmax); 
      printf("******************************************\n\n\n");
      system("pause");

      FILE *fpvel;
      
      if((fpvel=fopen(fnvel,"rb"))==NULL)
       { printf("Open file error ! <%s>\n",fnvel);exit(1);}
      for(i=1;i<=ncdp;i++)
      {
         for(j=1;j<=nz;j++)
         {
                fread(&v_initial[i][j],FSIZE,1,fpvel);
         }
      }

       k=0;
       vsum=0.0;
       pvel=0.0;
 
      for(i=1;i<=ncdp;i++)
      {
         for(j=1;j<=nz;j++)
         {
            if(j<=topz[i])
               v[i][j]=veltop;
            else if((j>topz[i])&&(j<=layer[i]))
                {
               v[i][j]=v_initial[i][j];
               k++;
               vsum+=v_initial[i][j];
                }
            else 
               v[i][j]=velbottom;
         }
      }

       pvel=vsum*1.0/(k*1.0);
       printf("pvel=%f\n",pvel);
       system("pause");

      FILE *fpvelmodel;
      fpvelmodel=fopen("step1_velmodel1.dat","wb");
      for(i=1;i<=ncdp;i++)
         for(j=zmin+1;j<=zmax;j++)
            fwrite(&v[i][j],4L,1,fpvelmodel);
      fclose(fpvelmodel);
      printf("Get the < step1_velmodel.dat > !\n");
      system("pause");
     

//=====================================================
//    slop of reflection points and smooth
//=====================================================
//            parameter instruction
//=====================================================
//        layerz(i):depth of reflection points
//        dy(i):slop of reflection points on layer
//=====================================================
       for(i=1;i<=ncdp;i++)
          layerz[i]=(layer[i]-1)*dz;
      
       step1_espl1(x0,layerz,ncdp,dy1,dyn,xx,2,dy,ddy,s,ds,dds,t,h);

       for(i=1;i<=ncdp;i++)
          if(fabs(dy[i])>=1.0)
             dy[i]=dy[i]/fabs(dy[i])*1.0;

//**********************   smooth dy  *************************//
       dy[1]=(dy[4]+dy[1]+dy[2]+dy[3])/4.0;
       dy[2]=dy[1];
       dy[3]=dy[1];
       dy[ncdp]=(dy[ncdp-1]+dy[ncdp-2]+dy[ncdp-3]+dy[ncdp])/4.0;
       dy[ncdp-1]=dy[ncdp];
       dy[ncdp-2]=dy[ncdp];
       for(i=4;i<=ncdp-3;i++)
          dy[i]=(dy[i-2]+dy[i-1]+dy[i]+dy[i+1]+dy[i+2])/5.0;

       FILE *fpsl;
       fpsl=fopen("step1_shot_location1.txt","w");
       for(i=1+xstart/nangle* ncdpstep_cig ;i<=xend/nangle* ncdpstep_cig ;i+=num* ncdpstep_cig )
           fprintf(fpsl,"%5d  %5d  %e\n",i,layer[i]-zmin,dy[i]);
        fclose(fpsl);
//=======================================================
//   Read data from adcig and select traces in adcig
//=======================================================
//              parameter instruction
//=======================================================
//        adcig_data(i,iz):data from ADCIG
//        valcig(i,iz):data without header information
//         ntra_per_adcig :traces selected from each ADCIG(odd number)
//        angl_select():angle of selected trace in ADCIGS
//        ind():trace number of selected traces
//        anglestart:start angle of ADCIGS
//        nmove:traces or angles moved towards 0
//              because of bad quality in depth
//=======================================================
        FILE *fpadcig;
        if((fpadcig=fopen(fn2,"rb"))==NULL)
        printf("Open file error ! <%s>\n",fn2);
        for(ix=1;ix<=nx;ix++)
           for(iz=1;iz<=nz;iz++)
           {
              fread(&adcig_data[ix][iz],4L,1,fpadcig);
              valcig[ix][iz]=adcig_data[ix][iz];
           }
        fclose(fpadcig);

        mid_per_adcig = ntra_per_adcig /2+1;
        angl_select[ mid_per_adcig ]=0.0;
        ind[ mid_per_adcig ]=nangle/2+1;
        for(i=1;i<= mid_per_adcig -1;i++)
        {
           angl_select[i]=anglestart+i-1+1+nmove;
           ind[i]=i+nmove;
        }
        for(i= mid_per_adcig +1;i<= ntra_per_adcig ;i++)
        {
           angl_select[i]=-1.0*angl_select[2* mid_per_adcig -i];
           ind[i]=2*ind[ mid_per_adcig ]-ind[2* mid_per_adcig -i];
        }
//===================================================
//     Get depth compared to the middle trace
//       in adcigs through cross-correlation
//===================================================
//            parameter instructions
//===================================================
//           valmax:maxium of middle angel in ADCIGS
//           z0:depth of layer in the middle of angel
//           corautopick:a subroutine which is related
//                       to cross-correlation and get
//                       depth compared to middle
//                       angel
//====================================================
        for(i=1;i<=ntrcid_mig_use ;i++)
        {
           if(i%10==0)
              printf(" i/ntrcid_mig_use =%5d/%d \n",i,ntrcid_mig_use );
           nzu=layer[(xstart-1)* ncdpstep_cig /nangle+(i-1)* ncdpstep_cig ]-1;
           nzd=layer[(xstart-1)* ncdpstep_cig /nangle+(i-1)* ncdpstep_cig ]+1;
    
           valmax=fabs(valcig[ind[ mid_per_adcig ]+xstart-(xstart/ncdpstep_cig)%nangle
                           +(i-1)*num*nangle][nzu]);
           z0=(nzu-1)*dz;

           for(k=nzu+1;k<=nzd;k++)
           {
             if(fabs(valcig[ind[ mid_per_adcig ]+xstart-(xstart/ncdpstep_cig)%nangle
                                      +(i-1)*num*nangle][k])>valmax)
             {
               valmax=fabs(valcig[ind[ mid_per_adcig ]+xstart-(xstart/ncdpstep_cig)%nangle
                                              +(i-1)*num*nangle][k]);
               z0=(k-1)*dz;
             }
           }

           step1_corautopick(i, ntra_per_adcig ,nzu,nzd,z,valcig,ind,xstart,
                          nangle, mid_per_adcig ,dz,z0,nx,nz,ncdpstep_cig);
           
           for(j=1;j<= ntra_per_adcig ;j++)
              printf("z =  %f *******************\n",z[j]);
           
           if(sign==1)
           { 
              for(j=1;j<= ntra_per_adcig ;j++)
              {
                 if(z[j]<=z[ mid_per_adcig ])
                   { z[j]=z[ mid_per_adcig ];}
                 else if(z[j]>(z[ mid_per_adcig ]+dzmax))
                   { z[j]=z[ mid_per_adcig ]+dzmax;}
              } 
           }else if(sign==-1)
           { 
              for(j=1;j<= ntra_per_adcig ;j++)
              {
                 if(z[j]>z[ mid_per_adcig ])
                   { z[j]=z[ mid_per_adcig ];}
                 else if(z[j]<=(z[ mid_per_adcig ]-dzmax))
                   { z[j]=z[ mid_per_adcig ]-dzmax;}
              } 
           }
//=====================================================
//           Calculate gama
//=====================================================
//         parameter instructions
//=====================================================
//       faj:radian deriven from angel
//       gama:ratio between migrated depth
//            and real depth
//       gamave:residual vurvature
//=====================================================
            gama[ mid_per_adcig ]=0.0;
            for(j=1;j<= mid_per_adcig -1;j++)
            {
               x=z[j]*z[j]/(z[ mid_per_adcig ]*z[ mid_per_adcig ])-1;
               faj=angl_select[j]*3.141593/180.0;
               y=tan(faj)*tan(faj);
               gama[j]=1.0-x/y;
               gama[j]=sqrt(1.0/gama[j]);
               //printf("gama=%f\n",gama[j]);
            }
           // system("pause");
            for(j= mid_per_adcig +1;j<= ntra_per_adcig ;j++)
            {
               x=z[j]*z[j]/(z[ mid_per_adcig ]*z[ mid_per_adcig ])-1;
               faj=angl_select[j]*3.141593/180.0;
               y=tan(faj)*tan(faj);
               gama[j]=1.0-x/y;
               gama[j]=sqrt(1.0/gama[j]);
               //printf("gama=%f\n",gama[j]);
            }
            gamave=0;
            for(j=1;j<= ntra_per_adcig ;j++)
               gamave+=gama[j];
            gamave=gamave/( ntra_per_adcig -1);

            if((sign==1)&&(gamave<1.0))
                printf("model adcig probelm!\n");
            if((sign!=1)&&(gamave>1.0))
                printf("model adcig probelm!\n");
            
//================================================
//              Get deltaz
//================================================
//         parameter instructions
//================================================
//        midd:the middle angle of ADCIG
//        deltz:residual depth
//================================================
            x=0.0;
            midd=nangle/2+1;
            for(j=1;j<=nangle;j++)
            {
               y=1.0*(j-midd)*3.141593/180.0;
               x=1.0+tan(y)*tan(y)*(1.0-1.0/(gamave*gamave));
               depth[(i-1)*nangle+j]=(int)(z[ mid_per_adcig ]*sqrt(x)+0.5);
            }
            for(j=1;j<=nangle;j++)
               deltz[(i-1)*nangle+j]=depth[(i-1)*nangle+j]-z[ mid_per_adcig ];
            for(j=1;j<= ntra_per_adcig ;j++)
              printf("z = %f gama= %f++++++++++++++\n",z[j],gama[j]);
            // system("pause");
        }//loop i=1:ntrcid_mig_use   ending
        
        FILE *fp7,*fp8;
        
        fp7=fopen(fn7,"w");
        fp8=fopen(fn8,"w");
        for(i=1;i<=ntrcid_mig_use *nangle;i+=1)
        {
           fprintf(fp7,"%f\n",deltz[i]);
           fprintf(fp8,"%f\n",depth[i]);
        }
        fclose(fp7);
        fclose(fp8);
        
//==================================================
//               Get dt
//==================================================
        for(i=1;i<=ntrcid_mig_use *nangle;i+=1)
        {
           j=i%nangle;
           icdp=(int)(1+xstart/nangle+num*(i-1)/nangle);

           if(j==0)
             angle=anglestart+(nangle-1)*dangle;
           else
             angle=anglestart+(j-1)*dangle;

           angle=angle/180.0*3.14159265;
           deltt[i]=factor*deltz[i]*cos(angle)*cos(atan(dy[icdp]))/pvel;
        }

        FILE *fp6;
        
        fp6=fopen(fn6,"w");
        for(i=1;i<=ntrcid_mig_use *nangle;i+=1)
        {
           fprintf(fp6,"%e\n",deltt[i]);
        }
        fclose(fp6);

//==================================================
//==       free all alloc
//==================================================
        free2float(valcig);
        free2float(adcig_data);

        free1float(deltt);
        free1float(depth);
        free1float(deltz);

        free1float(z);
        free1float(angl_select);
        free1float(gama);
        free1float(layer3);
        free1float(layerz);
        free1float(x0);
        free1float(dy);
        free1float(ddy);
        free1float(h);

        free1int(ind);
        free1int(layer);
        free1int(topz);



}
//**********************************************************************************//
void step1_corautopick(int i,int  ntra_per_adcig ,int nzu,int nzd,float *z,float **valcig,
              int *ind,int xstart,int nangle,int  mid_per_adcig ,float dz,float z0,int nx,int nz,
                   int ncdpstep_cig)
{

      int j,k,l;
      float *aut,**autval,max;

      aut=alloc1float(21+1);   zero1float(aut,21+1);
      autval=alloc2float(nzd-nzu+1+20+1, ntra_per_adcig +1);
      zero2float(autval,nzd-nzu+1+20+1, ntra_per_adcig +1);

      for(k=1;k<= ntra_per_adcig ;k++)
      {
         for(j=1;j<=10;j++)
            autval[k][j]=0.0;
         for(j=11;j<=nzd-nzu+1+10;j++)
         {
            //autval[k][j]=valcig[ind[k]+xstart-1+(i-1)*10*nangle][j+nzu-10];
            autval[k][j]=valcig[ind[k]+xstart-1+(i-1)*nangle][j+nzu-10];
            printf("autval = %e \n",autval[k][j]);
                
         }
         for(j=nzd-nzu+1+10+1;j<=nzd-nzu+1+20;j++)
            autval[k][j]=0.0;
         
      }
      for(k=1;k<= ntra_per_adcig ;k++)
      {
         //for(j=-10;j<=10;j++)
         for(j=0;j<=20;j++)
         {
            aut[j]=0.0;
            for(l=1+10;l<=nzd-nzu+1+10;l++)
            {
               aut[j]=aut[j]+autval[ mid_per_adcig ][l]*autval[k][l+j-10];
               printf("aut = %e\n",aut[j]);
            }
         }
         max=aut[0];
         z[k]=z0-10*dz;
         //for(j=-10;j<=10;j++)
         for(j=0;j<=20;j++)
         {
            if(max<aut[j])
            {
               max=aut[j];
               z[k]=z0+(j-10)*dz;
               
            }
         }
         printf("z = %f-----------------\n",z[k]);
      }
 //***end of autopickprogram
      free1float(aut);
      free2float(autval);
}


//**********************************************************************************//
void step1_espl1(float *x,float *y,int n,float dy1,float dyn,float *xx,
              int m,float *dy,float *ddy,float *s,float *ds,float *dds,float t,
              float *h)
{
     float alpha,beta,h1,h0;
     int i,j;

     dy[1]=0.0;
     h[1]=dy1;
     h0=x[2]-x[1];
     
     for(j=2;j<=n-1;j++)
     {
        h1=x[j+1]-x[j];
        alpha=h0/(h0+h1);
        beta=(1.0-alpha)*(y[j]-y[j-1])/h0;
        beta=3.0*(beta+alpha*(y[j+1]-y[j])/h1);
        dy[j]=-1.0*alpha/(2.0+(1.0-alpha)*dy[j-1]);
        h[j]=(beta-(1.0-alpha)*h[j-1]);
        h[j]=h[j]/(2.0+(1.0-alpha)*dy[j-1]);
        h0=h1;
     }
     dy[n]=dyn;
     for(j=n-1;j>=1;j--)
        dy[j]=dy[j]*dy[j+1]+h[j];
     for(j=1;j<=n-1;j++)
        h[j]=x[j+1]-x[j];
     for(j=1;j<=n-1;j++)
     {
        h1=h[j]*h[j];
        ddy[j]=6.0*(y[j+1]-y[j])/h1-2.0*(2.0*dy[j]+dy[j+1])/h[j];
     }
     h1=h[n-1]*h[n-1];
     ddy[n]=6.0*(y[n-1]-y[n])/h1+2.0*(2.0*dy[n]+dy[n-1])/h[n-1];
     
     for(i=1,t=0.0;i<=n-1;i++,t+=h1)
     {
        h1=0.5*h[i]*(y[i]+y[i+1]);
        h1=h1-h[i]*h[i]*h[i]*(ddy[i]+ddy[i+1])/24.0;
     }

     for(j=1;j<=m;j++)
     {
        if(xx[j]>=x[n])
           i=n-1;
        else
        {
           i=1;
           if(xx[j]>x[i+1])
           {
             do{i+=1;}
             while(xx[j]>x[i+1]);
           }
        }
       // printf("i=%d\n",i);
         h1=(x[i+1]-xx[j])/h[i];
         s[j]=(3.0*h1*h1-2.0*h1*h1*h1)*y[i];
         s[j]=s[j]+h[i]*(h1*h1-h1*h1*h1)*dy[i];
         ds[j]=6.0*(h1*h1-h1)*y[i]/h[i];
         ds[j]=ds[j]+(3.0*h1*h1-2.0*h1)*dy[i];
         dds[j]=(6.0-12.0*h1)*y[i]/(h[i]*h[i]);
         dds[j]=dds[j]+(2.0-6.0*h1)*dy[i]/h[i];
         h1=(xx[j]-x[i])/h[i];
         s[j]=s[j]+(3.0*h1*h1-2.0*h1*h1*h1)*y[i+1];
         s[j]=s[j]-h[i]*(h1*h1-h1*h1*h1)*dy[i+1];
         ds[j]=ds[j]-6.0*(h1*h1-h1)*y[i+1]/h[i];
         ds[j]=ds[j]+(3.0*h1*h1-2.0*h1)*dy[i+1];
         dds[j]=dds[j]+(6.0-12.0*h1)*y[i+1]/(h[i]*h[i]);
         dds[j]=dds[j]-(2.0-6.0*h1)*dy[i+1]/h[i];
     }
}
//**********************************************************************************//
void step1_readpar(char fnvel[],char fn1[],char fn2[],int *iflag,float *factor,int *sign,
              int *nx,int *nz,
              float *dx,float *dz,float *svel,float *pvel,int *ncdp,int *nangle,
              float *dangle,float *anglestart,int *xstart,int *xend,float *dzmax,
              int *nmove,int * ntra_per_adcig ,int *nzu,int *nzd)
{
      FILE *fp;
        char cc[250];

     if((fp=fopen("tomo-1.par","r"))==NULL)
      {printf("Open file error !\n");exit(1);}
      
      fscanf(fp,"%s%s\n",cc,fnvel);
      fscanf(fp,"%s%s\n",cc,fn1);
      fscanf(fp,"%s%s\n",cc,fn2);
      fscanf(fp,"%s%d\n",cc,iflag);
      fscanf(fp,"%s%f\n",cc,factor);
      fscanf(fp,"%s%d\n",cc,sign);
      fscanf(fp,"%s%d\n",cc,nx);
      fscanf(fp,"%s%d\n",cc,nz);
      fscanf(fp,"%s%f\n",cc,dx);
      fscanf(fp,"%s%f\n",cc,dz);
     // fscanf(fp,"%s%f\n",cc,pvel);
      fscanf(fp,"%s%d\n",cc,ncdp);
      fscanf(fp,"%s%d\n",cc,nangle);
      fscanf(fp,"%s%f\n",cc,dangle);
      fscanf(fp,"%s%f\n",cc,anglestart);
      fscanf(fp,"%s%d%d\n",cc,xstart,xend);
      fscanf(fp,"%s%f\n",cc,dzmax);
      fscanf(fp,"%s%d\n",cc,nmove);
      fscanf(fp,"%s%d\n",cc, ntra_per_adcig );

      fclose(fp);
}
//*********************************************************************************//
void step1_keyin(char fnvel[],char fn1[],char fn2[],int *iflag,float *factor,int *sign,int *nx,int *nz,
           float *dx,float *dz,float *svel,float *pvel,int *ncdp,int *nangle,
           float *dangle,float *anglestart,int *xstart,int *xend,float *dzmax,
           int *nmove,int * ntra_per_adcig ,int *nzu,int *nzd)
{ 
      FILE *fp;
      fp=fopen("step1_get_dt_and_velmodel-1.par","w");

      printf("***************************************************\n");
      printf("*             Parameters Inputing                 *\n");
      printf("***************************************************\n");
        printf("input the velocity porfile filename : fnvel\n");
        scanf("%s",fnvel);
        fprintf(fp,"%s\n",fnvel);
        printf("input the migration porfile filename : fn1\n");
        scanf("%s",fn1);
        fprintf(fp,"%s\n",fn1);
        printf("input the adcigs filename : fn2\n");
        scanf("%s",fn2);
        fprintf(fp,"%s\n",fn2);
        printf("input iflag=1:the first layer,else iflag=2  :iflag\n");
        scanf("%d",iflag);
        fprintf(fp,"%d\n",*iflag);
        printf("input factor\n");
        scanf("%f",factor);
        fprintf(fp,"%f\n",*factor);
        printf("Adcigs sign=1 :sign\n");
        scanf("%d",sign);
        fprintf(fp,"%d\n",*sign);
        printf("input the total traces of Adcigs : nx\n");
        scanf("%d",nx);
        fprintf(fp,"%d\n",*nx);
        printf("input the number of vertical sample_point in Adcigs: nz\n");
        scanf("%d",nz);
        fprintf(fp,"%d\n",*nz);
        printf("input the interval of lateral sample_point :dx\n");
        scanf("%f",dx);
        fprintf(fp,"%f\n",*dx);
        printf("input the interval of vertical sample_point :dz\n");
        scanf("%f",dz);
        fprintf(fp,"%f\n",*dz);
        printf("input the p&s_velcity :svel,pvel\n");
        scanf("%f %f",svel,pvel);
        fprintf(fp,"%f %f\n",*svel,*pvel);
        printf("input the number of cdps in the migration: ncdp\n");
        scanf("%d",ncdp);
        fprintf(fp,"%d\n",*ncdp);
        printf("input total traces in each ADCIG:nangle\n");
        scanf("%d",nangle);
        fprintf(fp,"%d\n",*nangle);
        printf("input:dangle\n");
        scanf("%f",dangle);
        fprintf(fp,"%f\n",*dangle);
        printf("input start angle:anglestart\n");
        scanf("%f",anglestart);
        fprintf(fp,"%f\n",*anglestart);
        printf("input start and end traces:xstart,xend\n");
        scanf("%d %d",xstart,xend);
        fprintf(fp,"%d %d\n",*xstart,*xend);
        printf("input dzmax\n");
        scanf("%f",dzmax);
        fprintf(fp,"%f\n",*dzmax);
        printf("input the number of traces moved:nmove\n");
        scanf("%d",nmove);
        fprintf(fp,"%d\n",*nmove);
        printf("input traces fitting in each ADCIG: ntra_per_adcig \n");
        scanf("%d", ntra_per_adcig );
        fprintf(fp,"%d\n",* ntra_per_adcig );
        printf("input the window size:nzu,nzd'\n");
        scanf("%d %d",nzu,nzd);
        fprintf(fp,"%d %d\n",*nzu,*nzd);
     
     fclose(fp);
}

