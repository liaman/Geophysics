//##################################################################
//##
//##          STEP4:       Update velmodel  
//##
//##     Ps:  
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
  void step4_keyin(char fn0[],char fn3[],char fn77[],int *allcdp,int *nz,
                  int *zmin,int *zmax,int *pro_flag,int *iflag);
  void step4_readpar(char fn0[],char fn3[],char fn77[],int *allcdp,int *nz,
                  int *zmin,int *zmax,int *pro_flag,int *iflag);
  void step4_vel_update(char fn0[],char fn3[],char fn77[],int allcdp,
                  int nz,int zmin,int zmax,int iflag,int pro_flag);


    int allcdp,nz,zmin,zmax,pro_flag,iflag;

    int flag;
    char fn0[250],fn3[250],fn77[250];

     printf("Ready to run this program ?\n");
     printf("Make sure <tomo-4.par> is done !\n");
     system("pause");

    printf("IF YOU INPUT THE OPERATING PARAMETERS FROM\n");
    printf(" Key in : 1\n");
    printf("Read par: 2\n");

    //scanf("%d",&flag);
    flag=2;

    if(flag==1)
    {
       step4_keyin(fn0,fn3,fn77,&allcdp,&nz,&zmin,&zmax,&pro_flag,&iflag);
    }else if(flag==2){
       step4_readpar(fn0,fn3,fn77,&allcdp,&nz,&zmin,&zmax,&pro_flag,&iflag);
    }
    else
    { printf("Please input the right flag number!\n");exit(1);}
  
      printf("*****************************\n");
      printf("*  fn0= %s\n",fn0);
      printf("*  fn3= %s\n",fn3);
      printf("*  fn77= %s\n",fn77);
      printf("*  iflag=%d\n",iflag);
      printf("*  pro_flag=%d\n*\n",pro_flag);
      printf("*  allcdp=%d\n",allcdp);
      printf("*  nz=%d\n",nz);
      printf("*  zmin=%d\n",zmin);
      printf("*  zmax=%d\n",zmax);
      printf("*****************************\n\n\n");
      printf("Ready to continue :\n");
      system("pause");


     step4_vel_update(fn0,fn3,fn77,allcdp,nz,zmin,zmax,iflag,pro_flag);

     printf("End of the program !\n");



}
//*****************************************************************************//
void step4_vel_update(char fn0[],char fn3[],char fn77[],int allcdp,
                  int nz,int zmin,int zmax,int iflag,int pro_flag)
{
    float **vel,**vel_date;
    int *vellayer;
    float **vel1,**ds;
    int *layer1,*layer2;

    float vv;
    int i,j,ix,iz;

    char fn[250]={"step3_delta_s1.dat"};
    char fn1[250]={"step1_velmodel1.dat"};
    char fn4[250]={"layer1.txt"};
    char fn5[250]={"layer2.txt"};

  vel=alloc2float(nz+1,allcdp+1);        zero2float(vel,nz+1,allcdp+1); 
  vel_date=alloc2float(nz+1,allcdp+1);   zero2float(vel_date,nz+1,allcdp+1);  
  vellayer=alloc1int(allcdp+1);          zero1int(vellayer,allcdp+1);
  ds=alloc2float(zmax-zmin+1,allcdp+1);  zero2float(ds,zmax-zmin+1,allcdp+1);
  vel1=alloc2float(zmax-zmin+1,allcdp+1);zero2float(vel1,zmax-zmin+1,allcdp+1);
  layer1=alloc1int(allcdp+1);            zero1int(layer1,allcdp+1);
  layer2=alloc1int(allcdp+1);            zero1int(layer2,allcdp+1);
  

//========================================
//==         Read interface data
//========================================
    FILE *fp4;//****layer1.txt
    if((fp4=fopen(fn4,"r"))==NULL)
      {printf("Open file error < %s > !\n",fn4);exit(1);}
    for(i=1;i<=allcdp;i++)
    {
       fscanf(fp4,"%d%d\n",&j,&layer1[i]);
       if(i%100==0)
          printf("j=%5d ,   layer1=%4d\n",j,layer1[i]);
    }
    fclose(fp4);

    FILE *fp5;//****layer2.txt
    if((fp5=fopen(fn5,"r"))==NULL)
      {printf("Open file error < %s > !\n",fn5);exit(1);}
    for(i=1;i<=allcdp;i++)
    {
       fscanf(fp5,"%d%d\n",&j,&layer2[i]);
       if(i%100==0)
          printf("j=%5d ,   layer2=%4d\n",j,layer2[i]);
    }
    fclose(fp5);
    
    FILE *fp1;//****velmodel.dat
    if((fp1=fopen(fn1,"rb"))==NULL)
      {printf("Open file error < %s > !\n",fn1);exit(1);}
      for(ix=1;ix<=allcdp;ix++)
         for(iz=1;iz<=zmax-zmin;iz++)
            fread(&vel1[ix][iz],4L,1,fp1);
    fclose(fp1);

    FILE *fp;//****delta_s.dat
    if((fp=fopen(fn,"rb"))==NULL)
      {printf("Open file error < %s > !\n",fn);exit(1);}
      for(ix=1;ix<=allcdp;ix++)
         for(iz=1;iz<=zmax-zmin;iz++)
            fread(&ds[ix][iz],4L,1,fp);
    fclose(fp);

    
    for(ix=1;ix<=allcdp;ix++)
       for(iz=1;iz<=zmax-zmin;iz++)
          vel1[ix][iz]=1.0/(1.0/vel1[ix][iz]+ds[ix][iz]);


    printf("allcdp = %d ,\nzmax-zmin = %d \n",allcdp,zmax-zmin);

//========================================
//==         Read initial veldata
//========================================
    FILE *fp0;//****vel_initial.dat
    if((fp0=fopen(fn0,"rb"))==NULL)
      {printf("Open file error < %s > !\n",fn0);exit(1);}
      for(ix=1;ix<=allcdp;ix++)
         for(iz=1;iz<=nz;iz++)
          { fread(&vel_date[ix][iz],4L,1,fp0);
            vel[ix][iz]=vel_date[ix][iz]; }
    fclose(fp0);

//========================================
//== Get average velocity of 'the layer'
//========================================
    i=0;
    vv=0.0;
    for(ix=1;ix<=allcdp;ix++)
     for(iz=1;iz<=zmax-zmin;iz++)
      if((iz<=(layer2[ix]-zmin))&&(iz>(layer1[ix]-zmin)))
      {
         vv+=vel1[ix][iz];
         i+=1;
      }     
     vv=vv/i;

     printf("average vel: %f\nnum = %d\n",vv,i);       

    for(ix=1;ix<=allcdp;ix++)
     for(iz=zmin+1;iz<=zmax;iz++)
      if((iz<=layer2[ix])&&(iz>layer1[ix]))
        vel[ix][iz]=vv;

//========================================
//==       Select vel interface
//========================================
    if(pro_flag==1)
    {
      for(ix=1;ix<=allcdp;ix++)
      {
        for(iz=layer2[ix];iz<=nz;iz++)
          if((vel[ix][iz]!=vv)&&(vel[ix][iz]!=vel[ix][iz+1]))
          {
             vellayer[ix]=iz;
             goto loop_10;
          }
      }

loop_10:
      for(ix=1;ix<=allcdp;ix++)
        for(iz=layer2[ix];iz<=vellayer[ix];iz++)
           vel[ix][iz]=vel[ix][vellayer[ix]+1];
    }//***end if

//========================================
//==      Velocity update
//========================================
    FILE *fp3;//****vel update
    fp3=fopen(fn3,"wb");
      for(ix=1;ix<=allcdp;ix++)
         for(iz=1;iz<=nz;iz++)
            fwrite(&vel[ix][iz],4L,1,fp3);
    fclose(fp3);

}
//*****************************************************************************//
void step4_readpar(char fn0[],char fn3[],char fn77[],int *allcdp,int *nz,
                  int *zmin,int *zmax,int *pro_flag,int *iflag)
{
    FILE *fp;
    char cc[250];
     if((fp=fopen("tomo-4.par","r"))==NULL)
      {printf("Open file error !\n");exit(1);}
      
      fscanf(fp,"%s%s\n",cc,fn0);
      fscanf(fp,"%s%s\n",cc,fn3);
      fscanf(fp,"%s%s\n",cc,fn77);
      fscanf(fp,"%s%d\n",cc,iflag);
      fscanf(fp,"%s%d\n",cc,pro_flag);
      fscanf(fp,"%s%d\n",cc,allcdp);
      fscanf(fp,"%s%d\n",cc,nz);
      fscanf(fp,"%s%d\n",cc,zmin);
      fscanf(fp,"%s%d\n",cc,zmax);
    fclose(fp); 
}
//*****************************************************************************//
void step4_keyin(char fn0[],char fn3[],char fn77[],int *allcdp,int *nz,
                  int *zmin,int *zmax,int *pro_flag,int *iflag)
{
      FILE *fp;
      fp=fopen("update_vel.par","w");

      printf("***************************************************\n");
      printf("*             Parameters Inputing                 *\n");
      printf("***************************************************\n");
        printf("input the filename of initial velocity field: fn0\n");
        scanf("%s",fn0);
        fprintf(fp,"%s\n",fn0);
        printf("input the filename of updated velocity field: fn3\n");
        scanf("%s",fn3);
        fprintf(fp,"%s\n",fn3);
        printf("input the filename of updated velocity field: fn77\n");
        scanf("%s",fn3);
        fprintf(fp,"%s\n",fn3);
        printf("the 1st layer--iflag=1,the deeper layer--iflag=2:iflag\n");
        scanf("%d",iflag);
        fprintf(fp,"%d\n",*iflag);
        printf("input the flag of interlayer process:  pro_flag\n");
        scanf("%d",pro_flag);
        fprintf(fp,"%d\n",*pro_flag);
        printf("input number of cdp:allcdp\n");
        scanf("%d",allcdp);
        fprintf(fp,"%d\n",*allcdp);
        printf("input the interval of vertical sample_point: nz\n");
        scanf("%d",nz);
        fprintf(fp,"%d\n",*nz);
        printf("input the minimum vertical sample_point :zmin\n");
        scanf("%d",zmin);
        fprintf(fp,"%d\n",*zmin);
        printf("input the maximum vertical sample_point :zmax\n");
        scanf("%d",zmax);
        fprintf(fp,"%d\n",*zmax);
     fclose(fp);

}
