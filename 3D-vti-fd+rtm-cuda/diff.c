#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include "/home/Toa/hc/cjbsegy.h"
#include "/home/Toa/hc/fft.c"
#include "/home/Toa/hc/alloc.c"
#include "/home/Toa/hc/complex.c"
int main()
{    



	int i,nx,ny,nt;
	float val1,val2,val3;
/*******************************************/
/*****change the name and the parameter*****/
        char FN1[250]={"waxian_shot.dat"};
        char FN2[250]={"waxian_shot_direct.dat"};
        char FN3[250]={"waxian_shot_mutedirect.dat"};


        nx=301;
        ny=301; 
        nt=1101;



	FILE *fp1,*fp2,*fp3;
	fp1=fopen(FN1,"rb");
	fp2=fopen(FN2,"rb");
	fp3=fopen(FN3,"wb");


	   for(i=0;i<nx*ny*nt;i++)
	   {
            fread(&val1,4,1,fp1);
            fread(&val2,4,1,fp2);
               val3=val1-val2;

            fwrite(&val3,4,1,fp3);
           }

     fclose(fp1);fclose(fp3);fclose(fp2);


        return 0;
}












