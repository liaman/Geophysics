#include<stdio.h>
#include<stdlib.h>

void main(){


          char FN1[250]={"tianjin_196/eps_1801_431.dat"};      
          char FN2[250]={"tianjin_196/eps_1801_862.dat"};


        float *a = malloc(sizeof(float)*1801*431);
        float *b = malloc(sizeof(float)*1801*862);

        FILE *fp = fopen(FN1,"rb");

        fread(a,4,1801*431,fp);

        fclose(fp);

        fp = fopen(FN2,"wb");


        int i, j, id;

        for(i=0;i<1801;i++){

            for(j=0;j<431;j++){

                id = i*431 + j;

                b[id*2] = a[id];

            }
        }

        for(i=0;i<1801;i++){

            for(j=0;j<862;j++){

                id = i*862 + j;

                if((id%2)!=0)
                    b[id] = b[id-1];

            }
        }
        fwrite(b,4,1801*862,fp);

       fclose(fp);






}
