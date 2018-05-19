#ifndef ALLOC_H
#define ALLOC_H

#include "cwp.h"
#include "complex.h"

/* allocate and free multi-dimensional arrays */
void *alloc1 (size_t n1, size_t size);
void *realloc1 (void *v, size_t n1, size_t size);
void **alloc2 (size_t n1, size_t n2, size_t size);
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6, 
                   size_t size);

void free1 (void *p);
void free2 (void **p);
void free3 (void ***p);
void free4 (void ****p);
void free5 (void *****p);
void free6 (void ******p);

int *alloc1int (size_t n1);
int *realloc1int (int *v, size_t n1);
int **alloc2int (size_t n1, size_t n2);
int ***alloc3int (size_t n1, size_t n2, size_t n3);
float *alloc1float (size_t n1);
float *realloc1float (float *v, size_t n1);
float **alloc2float (size_t n1, size_t n2);
float ***alloc3float (size_t n1, size_t n2, size_t n3);

float ****alloc4float (size_t n1, size_t n2, size_t n3, size_t n4);
void free4float (float ****p);
float *****alloc5float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5float (float *****p);
float ******alloc6float (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6);
void free6float (float ******p);
int ****alloc4int (size_t n1, size_t n2, size_t n3, size_t n4);
void free4int (int ****p);
int *****alloc5int (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
void free5int (int *****p);

unsigned char ******alloc6uchar(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5, size_t n6);
unsigned char *****alloc5uchar(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5);
void free6uchar(unsigned char ******p);
void free5uchar(unsigned char *****p);

unsigned short ******alloc6ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5, size_t n6);
unsigned short *****alloc5ushort(size_t n1,size_t n2,size_t n3,size_t n4,
        size_t n5);
unsigned short ***alloc3ushort(size_t n1,size_t n2,size_t n3);
unsigned short **alloc2ushort(size_t n1,size_t n2);

void free2ushort(unsigned short **p);
void free3ushort(unsigned short ***p);
void free5ushort(unsigned short *****p);
void free6ushort(unsigned short ******p);

double *alloc1double (size_t n1);
double *realloc1double (double *v, size_t n1);
double **alloc2double (size_t n1, size_t n2);
double ***alloc3double (size_t n1, size_t n2, size_t n3);
complex *alloc1complex (size_t n1);
complex *realloc1complex (complex *v, size_t n1);
complex **alloc2complex (size_t n1, size_t n2);
complex ***alloc3complex (size_t n1, size_t n2, size_t n3);
complex ****alloc4complex (size_t n1, size_t n2, size_t n3, size_t n4);
complex *****alloc5complex (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);

void free1int (int *p);
void free2int (int **p);
void free3int (int ***p);

void free1float (float *p);
void free2float (float **p);
void free3float (float ***p);

void free1double (double *p);
void free2double (double **p);
void free3double (double ***p);

void free1complex (complex *p);
void free2complex (complex **p);
void free3complex (complex ***p);
void free4complex (complex ****p);
void free5complex (complex *****p);

void zero1int(int *p, size_t n1);
void zero2int(int **p, size_t n1, size_t n2);
void zero3int(int ***p, size_t n1, size_t n2, size_t n3);

void zero2ushort(unsigned short **p, size_t n1, size_t n2);

void zero3ushort(unsigned short ***p, size_t n1, size_t n2, size_t n3);

void zero1float(float *p, size_t n1);
void zero2float(float **p, size_t n1, size_t n2);
void zero3float(float ***p, size_t n1, size_t n2, size_t n3);
void zero4float(float ****p, size_t n1, size_t n2, size_t n3, size_t n4);

void zero1double (double *p, size_t n1);
void zero2double (double **p, size_t n1, size_t n2);
void zero3double (double ***p, size_t n1, size_t n2, size_t n3);

void zero1complex (complex *p, size_t n1);
void zero2complex (complex **p, size_t n1, size_t n2);
void zero3complex (complex ***p, size_t n1, size_t n2, size_t n3);
void zero4complex (complex ****p, size_t n1, size_t n2, size_t n3, size_t n4);
void zero5complex (complex *****p, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);

#endif
