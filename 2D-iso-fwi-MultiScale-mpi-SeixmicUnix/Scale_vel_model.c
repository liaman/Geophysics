#include<stdio.h>
#include<math.h>
#include<stdlib.h>
int main()
{    
	int i,j,k,vnx,vnz;
	float **vp;
/*******************************************/
/*****change the name and the parameter*****/
        char FN1[250]={"mar_vel_737_750.dat"};
        char FN2[250]={"mar_vel_600_700.dat"};
        vnx=737;
        vnz=750; 

	vp=(float **)calloc(vnz,sizeof(float *));
           for(i=0;i<vnz;i++)
		vp[i]=(float *)calloc(vnx,sizeof(float));

        for(i=0;i<vnx;i++)
             for(j=0;j<vnz;j++)
		 vp[i][j]=0.0;

/*********the velocity model ENDing**********/
/********************************************/
	FILE *fp1;
	fp1=fopen(FN1,"rb");
	for(i=0;i<vnx;i++)
	   for(j=0;j<vnz;j++)
	      fread(&vp[i][j],4,1,fp1);
	fclose(fp1); 

        for(i=0;i<vnx;i++)
             for(j=0;j<40;j++)
		 vp[i][j]=1600;        


        FILE *fp2;
	fp2=fopen(FN2,"wb");
	for(i=87;i<vnx-50;i++)
	   for(j=0;j<700;j++)
	      fwrite(&vp[i][j],4,1,fp2);
	fclose(fp2); 



        return 0;
}

