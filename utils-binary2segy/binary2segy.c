#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "hc/bhdr.h"
#include "hc/hdr.h"
#include "hc/cjbsegy.h"
#include "math.h"

typedef cjbsegy segy;

#define SU_NKEYS        80      /* Number of key header words */
#define HDRBYTES        240     /* Bytes in the trace header */
#define EBCBYTES        3200    /* Bytes in the card image EBCDIC block */
#define BNYBYTES        400     /* Bytes in the binary coded block      */
#define SEGY_HDRBYTES   240     /* Bytes in the tape trace header       */
#define SEGY_NKEYS      71      /* Number of mandated header fields     */
#define BHED_NKEYS      27      /* Number of mandated binary fields     */
/***** SEG-Y header *********/
       Y_3200   y3200;
       bhed     head_400;
       cjbsegy  tr, vtr;
/***** SU function *********/
void swap_short_2(short *tni2);
void swap_u_short_2(unsigned short *tni2);
void swap_int_4(int *tni4);
void swap_u_int_4(unsigned int *tni4);
void swap_long_4(long *tni4);
void swap_u_long_4(unsigned long *tni4);
void swap_float_4(float *tnf4);
void swap_double_8(double *tndd8);
void swaphval(segy *tr, int index);


/**
 *  Complier
 *
 * -bash-4.1$ gcc -o a dat2su.c -lm 
 * -bash-4.1$ ./a
 */

void main()
{
    int i,j,k,ntrace,ns,cc;
    char FN1[250]={"v1.su"};//你所需要的su道头的输入文件
    char FN2[250]={"v2.dat"};//你输入的二进制裸数据
    char FN3[250]={"v2.su"};//输出的文件：v1.su的道头+v2.dat的数据体

    ntrace = 1700;//最大的道数：三维情况等于nx*ny*炮数, 二维情况等于nx*炮数 or ng*炮数
    
    ns = 601;//每一道的纵向采样点数
    
    FILE *fp1 = fopen(FN1,"rb");
    FILE *fp2 = fopen(FN2,"rb");
    FILE *fp3 = fopen(FN3,"wb");
    
    float *trace = (float*)malloc(sizeof(float)*ns);
    
    int itrace;
    for(itrace = 0; itrace < ntrace; itrace ++)
    {
        if(itrace%200 == 0) 
            printf("itrace = %5d/%d;\n",itrace, ntrace);
        
        fread(&tr, sizeof(cjbsegy), 1, fp1);
        
        /**
         *  don't do this swaphval function
         */
        //for (k = 0; k < SU_NKEYS; ++k) 
        //    swaphval(&tr,k);
        
        fread(trace, sizeof(float), ns, fp1);
        fread(trace, sizeof(float), ns, fp2);
        
        for(j=0;j<ns;j++)
            swap_float_4(&trace[j]);
        
        fwrite(&tr, sizeof(cjbsegy), 1, fp3);
        fwrite(trace, sizeof(float), ns, fp3);
    }

    free(trace);

}//end of main


void swap_short_2(short *tni2)
/**************************************************************************
swap_short_2		swap a short integer
***************************************************************************/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));
}

void swap_u_short_2(unsigned short *tni2)
/**************************************************************************
swap_u_short_2		swap an unsigned short integer
***************************************************************************/
{
 *tni2=(((*tni2>>8)&0xff) | ((*tni2&0xff)<<8));
}

void swap_int_4(int *tni4)
/**************************************************************************
swap_int_4		swap a 4 byte integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}

void swap_u_int_4(unsigned int *tni4)
/**************************************************************************
swap_u_int_4		swap an unsigned integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}

void swap_long_4(long *tni4)
/**************************************************************************
swap_long_4		swap a long integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}

void swap_u_long_4(unsigned long *tni4)
/**************************************************************************
swap_u_long_4		swap an unsigned long integer
***************************************************************************/
{
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}

void swap_float_4(float *tnf4)
/**************************************************************************
swap_float_4		swap a float
***************************************************************************/
{
 int *tni4=(int *)tnf4;
 *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
	    ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}

void swap_double_8(double *tndd8)
/**************************************************************************
swap_double_8		swap a double
***************************************************************************/
{
  char *tnd8=(char *)tndd8;
  char tnc;

  tnc= *tnd8;
  *tnd8= *(tnd8+7);
  *(tnd8+7)=tnc;

  tnc= *(tnd8+1);
  *(tnd8+1)= *(tnd8+6);
  *(tnd8+6)=tnc;

  tnc= *(tnd8+2);
  *(tnd8+2)= *(tnd8+5);
  *(tnd8+5)=tnc;

  tnc= *(tnd8+3);
  *(tnd8+3)= *(tnd8+4);
  *(tnd8+4)=tnc;
}
void swaphval(segy *tr, int index)
{
        register char *tp= (char *) tr;

        switch(*(hdr[index].type)) {
        case 'h': swap_short_2((short*)(tp + hdr[index].offs));
        break;
        case 'u': swap_u_short_2((unsigned short*)(tp + hdr[index].offs));
        break;
        case 'i': swap_int_4((int*)(tp + hdr[index].offs));
        break;
        case 'p': swap_u_int_4((unsigned int*)(tp + hdr[index].offs));
        break;
        case 'l': swap_long_4((long*)(tp + hdr[index].offs));
        break;
        case 'v': swap_u_long_4((unsigned long*)(tp + hdr[index].offs));
        break;
        case 'f': swap_float_4((float*)(tp + hdr[index].offs));
        break;
        case 'd': swap_double_8((double*)(tp + hdr[index].offs));
        break;
        default: err("%s: %s: unsupported data type", __FILE__, __LINE__);
        break;
        }

}

