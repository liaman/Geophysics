/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */


#include "cwp.h"

/*********************** self documentation **********************/
/*****************************************************************************
EALLOC - Allocate and free multi-dimensional arrays with error reports.

ealloc1			allocate a 1d array
erealloc1		reallocate a 1d array
ealloc2			allocate a 2d array
ealloc3			allocate a 3d array
ealloc4			allocate a 4d array
ealloc5                 allocate a 5d array
ealloc6                 allocate a 6d array
ealloc1int		allocate a 1d array of ints
erealloc1int		reallocate a 1d array of ints
ealloc2int		allocate a 2d array of ints
ealloc3int		allocate a 3d array of ints
ealloc4int              allocate a 4d array of ints
ealloc5int              allocate a 5d array of ints
ealloc1float		allocate a 1d array of floats
erealloc1float		reallocate a 1d array of floats
ealloc2float		allocate a 2d array of floats
ealloc3float		allocate a 3d array of floats
ealloc4float            allocate a 4d array of floats
ealloc5float            allocate a 5d array of floats
ealloc6float            allocate a 6d array of floats
ealloc5ushort           allocate a 5d array of unsigned shorts
ealloc5uchar            allocate a 5d array of unsigned chars
ealloc6ushort           allocate a 6d array of unsigned shorts
ealloc6uchar            allocate a 6d array of unsigned chars
ealloc1double		allocate a 1d array of doubles
erealloc1double		reallocate a 1d array of doubles
ealloc2double		allocate a 2d array of doubles
ealloc3double		allocate a 3d array of doubles
ealloc1complex		allocate a 1d array of complex values	
erealloc1complex	reallocate a 1d array of complex values
ealloc2complex		allocate a 2d array of complex values
ealloc3complex		allocate a 3d array of complex values

*****************************************************************************
Function Prototypes:
void *ealloc1 (size_t n1, size_t size);
void *erealloc1 (void *v, size_t n1, size_t size);
void **ealloc2 (size_t n1, size_t n2, size_t size);
void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size);
void ****ealloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****ealloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******ealloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, 
                   size_t n6, size_t size);

int *ealloc1int(size_t n1);
int *erealloc1int(int *v, size_t n1);
int **ealloc2int(size_t n1, size_t n2);
int ***ealloc3int(size_t n1, size_t n2, size_t n3);
int ****ealloc4int(size_t n1, size_t n2, size_t n3, size_t n4);
int *****ealloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);

float *ealloc1float(size_t n1);
float *erealloc1float(float *v, size_t n1);
float **ealloc2float(size_t n1, size_t n2);
float ***ealloc3float(size_t n1, size_t n2, size_t n3);
float ****ealloc4float(size_t n1, size_t n2, size_t n3, size_t n4);
float *****ealloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
float ******ealloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, 
                        size_t n6);

unsigned short *****ealloc5ushort(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
unsigned char *****ealloc5uchar(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5);
unsigned short ******ealloc6ushort(size_t n1, size_t n2, size_t n3, size_t n4, 
      size_t n5, size_t n6);
unsigned char ******ealloc6uchar(size_t n1, size_t n2, size_t n3, size_t n4, 
      size_t n5, size_t n6);

double *ealloc1double(size_t n1);
double *erealloc1double(double *v, size_t n1);
double **ealloc2double(size_t n1, size_t n2);
double ***ealloc3double(size_t n1, size_t n2, size_t n3);
complex *ealloc1complex(size_t n1);
complex *erealloc1complex(complex *v, size_t n1);
complex **ealloc2complex(size_t n1, size_t n2);
complex ***ealloc3complex(size_t n1, size_t n2, size_t n3);

*****************************************************************************
Notes:
These routines simply call those in ../../cwp/lib/alloc.c and issue
an error message via the syserr routine if the underlying malloc
came up empty.  See alloc.c for notes on the routines.

The routines in ../../cwp/lib/alloc.c were written by Dave Hale
(Zhaobo Meng added 4d (except alloc4), 5d and 6d functions).

*****************************************************************************
Author: Jack Cohen, Center for Wave Phenomena
Zhaobo Meng added 4d (except ealloc4), 5d and 6d functions
*****************************************************************************/
/**************** end self doc ********************************/

#include "par.h"
#define ERROR	NULL

/* allocate a 1-d array */
void *ealloc1 (size_t n1, size_t size)
{
	void *p;

	if (ERROR == (p=alloc1(n1, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* re-allocate a 1-d array */
void *erealloc1 (void *v, size_t n1, size_t size)
{
	void *p;

	if (ERROR == (p=realloc1(v, n1, size)))
		syserr("%s: realloc failed", __FILE__);
	return p;
}

/* allocate a 2-d array */
void **ealloc2 (size_t n1, size_t n2, size_t size)
{
	void **p;

	if (ERROR == (p=alloc2(n1, n2, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 3-d array */
void ***ealloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
	void ***p;

	if (ERROR == (p=alloc3(n1, n2, n3, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 4-d array */
void ****ealloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
	void ****p;

	if (ERROR == (p=alloc4(n1, n2, n3, n4, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 5-d array */
void *****ealloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
	void *****p;

	if (ERROR == (p=alloc5(n1, n2, n3, n4, n5, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 6-d array */
void ******ealloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, 
                    size_t n6, size_t size)
{
	void ******p;

	if (ERROR == (p=alloc6(n1, n2, n3, n4, n5, n6, size)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 1-d array of ints */
int *ealloc1int(size_t n1)
{
	int *p;

	if (ERROR == (p=alloc1int(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* re-allocate a 1-d array of ints */
int *erealloc1int(int *v, size_t n1)
{
	int *p;

	if (ERROR == (p=realloc1int(v,n1)))
		syserr("%s: realloc failed", __FILE__);
	return p;
}


/* allocate a 2-d array of ints */
int **ealloc2int(size_t n1, size_t n2)
{
	int **p;

	if (ERROR == (p=alloc2int(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 3-d array of ints */
int ***ealloc3int(size_t n1, size_t n2, size_t n3)
{
	int ***p;

	if (ERROR == (p=alloc3int(n1, n2, n3)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 4-d array of ints */
int ****ealloc4int(size_t n1, size_t n2, size_t n3, size_t n4)
{
	int ****p;

	if (ERROR == (p=alloc4int(n1, n2, n3,n4)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 5-d array of ints */
int *****ealloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
	int *****p;

	if (ERROR == (p=alloc5int(n1, n2, n3, n4, n5)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 1-d array of floats */
float *ealloc1float(size_t n1)
{
	float *p;

	if (ERROR == (p=alloc1float(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* re-allocate a 1-d array of floats */
float *erealloc1float(float *v, size_t n1)
{
	float *p;

	if (ERROR == (p=realloc1float(v, n1)))
		syserr("%s: realloc failed", __FILE__);
	return p;
}


/* allocate a 2-d array of floats */
float **ealloc2float(size_t n1, size_t n2)
{
	float **p;

	if (ERROR == (p=alloc2float(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 3-d array of floats */
float ***ealloc3float(size_t n1, size_t n2, size_t n3)
{
	float ***p;

	if (ERROR == (p=alloc3float(n1, n2, n3)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 4-d array of floats */
float ****ealloc4float(size_t n1, size_t n2, size_t n3, size_t n4)
{
	float ****p;

	if (ERROR == (p=alloc4float(n1, n2, n3, n4)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 5-d array of floats */
float *****ealloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
	float *****p;

	if (ERROR == (p=alloc5float(n1, n2, n3, n4, n5)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 6-d array of floats */
float ******ealloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5,
                         size_t n6)
{
	float ******p;

	if (ERROR == (p=alloc6float(n1, n2, n3, n4, n5, n6)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}

/* allocate a 5-d array of unsigned shorts */
unsigned short *****ealloc5ushort(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
        unsigned short *****p;

        if (ERROR == (p=alloc5ushort(n1, n2, n3, n4, n5)))
                syserr("%s: malloc failed", __FILE__);
        return p;
}

/* allocate a 6-d array of unsigned shorts */
unsigned short ******ealloc6ushort(size_t n1, size_t n2, size_t n3, size_t n4, 
      size_t n5, size_t n6)
{
        unsigned short ******p;

        if (ERROR == (p=alloc6ushort(n1, n2, n3, n4, n5, n6)))
                syserr("%s: malloc failed", __FILE__);
        return p;
}

/* allocate a 5-d array of unsigned chars */
unsigned char *****ealloc5uchar(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
        unsigned char *****p;

        if (ERROR == (p=alloc5uchar(n1, n2, n3, n4, n5)))
                syserr("%s: malloc failed", __FILE__);
        return p;
}

/* allocate a 1-d array of doubles */
double *ealloc1double(size_t n1)
{
	double *p;

	if (ERROR == (p=alloc1double(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* re-allocate a 1-d array of doubles */
double *erealloc1double(double *v, size_t n1)
{
	double *p;

	if (ERROR == (p=realloc1double(v, n1)))
		syserr("%s: realloc failed", __FILE__);
	return p;
}


/* allocate a 2-d array of doubles */
double **ealloc2double(size_t n1, size_t n2)
{
	double **p;

	if (ERROR == (p=alloc2double(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 3-d array of doubles */
double ***ealloc3double(size_t n1, size_t n2, size_t n3)
{
	double ***p;

	if (ERROR == (p=alloc3double(n1, n2, n3)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 1-d array of complexs */
complex *ealloc1complex(size_t n1)
{
	complex *p;

	if (ERROR == (p=alloc1complex(n1)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* re-allocate a 1-d array of complexs */
complex *erealloc1complex(complex *v, size_t n1)
{
	complex *p;

	if (ERROR == (p=realloc1complex(v, n1)))
		syserr("%s: realloc failed", __FILE__);
	return p;
}


/* allocate a 2-d array of complexs */
complex **ealloc2complex(size_t n1, size_t n2)
{
	complex **p;

	if (ERROR == (p=alloc2complex(n1, n2)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}


/* allocate a 3-d array of complexs */
complex ***ealloc3complex(size_t n1, size_t n2, size_t n3)
{
	complex ***p;

	if (ERROR == (p=alloc3complex(n1, n2, n3)))
		syserr("%s: malloc failed", __FILE__);
	return p;
}



#ifdef TEST
main()
{
	short   *hv, **hm;
	int     *iv, **im;
	float   *fv, **fm;
	double  *dv, **dm;
	size_t i1, i2;
	size_t n1, n2;

	scanf("%d %*[^\n]", &n1);
	scanf("%d %*[^\n]", &n2);

	/* Exercise 1-d routines */
	hv = (short *) ealloc1(n1, sizeof(short));
	iv = ealloc1int(n1);
	fv = ealloc1float(n1);
	dv = ealloc1double(n1);

	for (i1 = 0; i1 < n1; ++i1) {
		hv[i1] = i1;
		iv[i1] = i1;
		fv[i1]  = (float) i1;
		dv[i1] = (double) i1;
	}

	printf("short vector:\n");
	for (i1 = 0; i1 < n1; ++i1) {
		printf("hv[%d] = %hd\n", i1, hv[i1]);
	}
	putchar('\n');

	printf("int vector:\n");
	for (i1 = 0; i1 < n1; ++i1) {
		printf("iv[%d] = %d\n", i1, iv[i1]);
	}
	putchar('\n');

	printf("float vector:\n");
	for (i1 = 0; i1 < n1; ++i1) {
		printf("fv[%d] = %f\n", i1, fv[i1]);
	}
	putchar('\n');

	printf("double vector:\n");
	for (i1 = 0; i1 < n1; ++i1) {
		printf("dv[%d] = %lf\n", i1, dv[i1]);
	}
	putchar('\n');


	free1(hv);
	free1int(iv);
	free1float(fv);
	free1double(dv);


	/* Exercise 2-d routines */
	hm = (short **) ealloc2(n1, n2, sizeof(short));
	im = ealloc2int(n1, n2);
	fm = ealloc2float(n1, n2);
	dm = ealloc2double(n1, n2);

	for (i2 = 0; i2 < n2; ++i2) {
		for (i1 = 0; i1 < n1; ++i1) {
			hm[i2][i1] = i1 + 2*i2;
			im[i2][i1] = i1 + 2*i2;
			fm[i2][i1] = (float) (i1 + 2*i2);
			dm[i2][i1] = (double) (i1 + 2*i2);
		}
	}

	printf("short matrix:\n");
	for (i2 = 0; i2 < n2; ++i2) {
		for (i1 = 0; i1 < n1; ++i1) {
			printf("hm[%d, %d] = %hd ", i2, i1, hm[i2][i1]);
		}
		putchar('\n');
	}
	putchar('\n');

	printf("int matrix:\n");
	for (i2 = 0; i2 < n2; ++i2) {
		for (i1 = 0; i1 < n1; ++i1) {
			printf("im[%d, %d] = %d ", i2, i1, im[i2][i1]);
		}
		putchar('\n');
	}
	putchar('\n');

	printf("float matrix:\n");
	for (i2 = 0; i2 < n2; ++i2) {
		for (i1 = 0; i1 < n1; ++i1) {
			printf("fm[%d, %d] = %f ", i2, i1, fm[i2][i1]);
		}
		putchar('\n');
	}
	putchar('\n');

	printf("double matrix:\n");
	for (i2 = 0; i2 < n2; ++i2) {
		for (i1 = 0; i1 < n1; ++i1) {
			printf("dm[%d, %d] = %lf ", i2, i1, dm[i2][i1]);
		}
		putchar('\n');
	}
	putchar('\n');

	free2((void*)hm);
	free2int(im);
	free2float(fm);
	free2double(dm);

	return EXIT_SUCCESS;
}
#endif
