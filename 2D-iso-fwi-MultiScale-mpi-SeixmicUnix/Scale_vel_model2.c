#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
int main()
{    
	int i,j,k,nx,nz,flag;
	float **vp;
/*******************************************/
/*****change the name and the parameter*****/
        char FN1[250]={"fault_vel_560_360.dat"};
        nx=560;
        nz=360; 
        flag=1;//如果不等于1,则输出梯度场层模型
/******************************************/
	vp=alloc2float(nz,nx);
       
 
           for(i=0;i<nx;i++)
              for(j=0;j<nz;j++)
                { vp[i][j]=0.0;}
 
/********************************************/
/*************the vel and rho*************/
     

/***************the fault model***************/
      if(flag==1)
      {
        for(i=0;i<nx;i++)
	  for(j=0;j<nz;j++)
             {   
                 if(j<=80)
                    {vp[i][j]=2000;}
                 else 
                    {vp[i][j]=2300;}
                 if(i-50<=0.9*j+50)
                 {
                     if(j>=140)
                         vp[i][j]=2800;
                     if(j>=190)
                                 vp[i][j]=3100;
                      if(j>=200)
                         vp[i][j]=3500;
                   

                 }
                 if(i-50>=0.9*j+50)
                 {
                     if(j>=185)
                         vp[i][j]=2800;
                      if(j>=235)
                          vp[i][j]=3100;
                      if(j>=245)
                         vp[i][j]=3500;
                     
                      
                 }
           
                 if(pow((i-280),2)+pow((j-460),2)<=36000)
                    {
                        vp[i][j]=3800;
                        if(j<300)
                         {vp[i][j]=2800;}
                        else if(j<335&&j>=290)
                         {vp[i][j]=3300;}
                      }
             }
        }
/*********************得到梯度速度场***********************/
        if(flag!=1)
        { 
          for(i=0;i<nx;i++)
	     for(j=0;j<nz;j++)
             {   
                 if(j<=30)
                    {vp[i][j]=2000;}
                 else if(j>30&&j<=50)
                    {vp[i][j]=2500;} 
                 else if(j>50&&j<=70)
                    {vp[i][j]=2600;}  
                 else if(j>70&&j<=90)
                    {vp[i][j]=2700;}
                 else if(j>90&&j<=110)
                    {vp[i][j]=2800;}
                 else if(j>110&&j<=130)
                    {vp[i][j]=3000;}
                 else if(j>130&&j<=150)
                    {vp[i][j]=3250;}
                 else if(j>150&&j<=180)
                    {vp[i][j]=3500;}
             }
          
        }



/**********************************×××××××××××××**********/ 
      


/********************************************/
	FILE *fp1;
	fp1=fopen(FN1,"wb");
 
	for(i=0;i<nx;i++)
	   for(j=0;j<nz;j++)
	     { fwrite(&vp[i][j],4,1,fp1);}
  
        fclose(fp1);
       free2float(vp);
       printf("The vel  model is done!\n");
        return 0;
}

