/*********************** self documentation **********************/
/******************************************************************
ALLOC - Allocate and free multi-dimensional arrays

alloc1		allocate a 1-d array
realloc1	re-allocate a 1-d array
free1		free a 1-d array
alloc2		allocate a 2-d array
free2		free a 2-d array
alloc3		allocate a 3-d array
free3		free a 3-d array
alloc4		allocate a 4-d array
free4		free a 4-d array
alloc5		allocate a 5-d array
free5		free a 5-d array
alloc6		allocate a 6-d array
free6		free a 6-d arrayalloc1int	
allocate a 1-d array of ints
realloc1int	re-allocate a 1-d array of ints
free1int	free a 1-d array of ints
alloc2int	allocate a 2-d array of ints
free2int	free a 2-d array of ints
alloc3int	allocate a 3-d array of ints
free3int	free a 3-d array of ints
alloc1float	allocate a 1-d array of floats
realloc1float	re-allocate a 1-d array of floats
free1float	free a 1-d array of floats
alloc2float	allocate a 2-d array of floats
free2float	free a 2-d array of floats
alloc3float	allocate a 3-d array of floats
free3float	free a 3-d array of floats
alloc4float	allocate a 4-d array of floats 
free4float      free a 4-d array of floats 
alloc5float     allocate a 5-d array of floats 
free5float      free a 5-d array of floats 
alloc6float     allocate a 6-d array of floats 
free6float      free a 6-d array of floats 
alloc4int       allocate a 4-d array of ints 
free4int        free a 4-d array of ints 
alloc5int       allocate a 5-d array of ints 
free5int        free a 5-d array of ints 
alloc5uchar	allocate a 5-d array of unsigned chars 
free5uchar	free a 5-d array of unsiged chars 
alloc2ushort    allocate a 2-d array of unsigned shorts 
free2ushort     free a 2-d array of unsiged shorts
alloc3ushort    allocate a 3-d array of unsigned shorts 
free3ushort     free a 3-d array of unsiged shorts
alloc5ushort    allocate a 5-d array of unsigned shorts 
free5ushort     free a 5-d array of unsiged shorts
alloc6ushort    allocate a 6-d array of unsigned shorts 
free6ushort     free a 6-d array of unsiged shorts
alloc1double	allocate a 1-d array of doubles
realloc1double	re-allocate a 1-d array of doubles
free1double	free a 1-d array of doubles
alloc2double	allocate a 2-d array of doubles
free2double	free a 2-d array of doubles
alloc3double	allocate a 3-d array of doubles
free3double	free a 3-d array of doubles
alloc1complex	allocate a 1-d array of complexs
realloc1complex	re-allocate a 1-d array of complexs
free1complex	free a 1-d array of complexs
alloc2complex	allocate a 2-d array of complexs
free2complex	free a 2-d array of complexs
alloc3complex	allocate a 3-d array of complexs
free3complex	free a 3-d array of complexs
alloc4complex	allocate a 4-d array of complexs
free4complex	free a 4-d array of complexs
alloc5complex	allocate a 5-d array of complexs
free5complex	free a 5-d array of complexs

zero1int        initialize the 1-d int array with zero
zero2int        initialize the 2-d int array with zero
zero3int        initialize the 3-d int array with zero

zero1float      initialize the 1-d float array with zero
zero2float      initialize the 2-d float array with zero
zero3float      initialize the 3-d float array with zero
zero4float      initialize the 4-d float array with zero

zero1double     initialize the 1-d double array with zero
zero2double     initialize the 2-d double array with zero
zero3double     initialize the 3-d double array with zero

zero1complex    initialize the 1-d complex array with zero
zero2complex    initialize the 2-d complex array with zero
zero3complex    initialize the 3-d complex array with zero
zero4complex    initialize the 4-d complex array with zero
zero5complex    initialize the 5-d complex array with zero

******************************************************************************
Notes:
The functions defined below are intended to simplify manipulation
of multi-dimensional arrays in scientific programming in C.  These
functions are useful only because true multi-dimensional arrays
in C cannot have variable dimensions (as in FORTRAN).  For example,
the following function IS NOT valid in C:
	void badFunc(a,n1,n2)
	float a[n2][n1];
	{
		a[n2-1][n1-1] = 1.0;
	}
However, the following function IS valid in C:
	void goodFunc(a,n1,n2)
	float **a;
	{
		a[n2-1][n1-1] = 1.0;
	}
Therefore, the functions defined below do not allocate true
multi-dimensional arrays, as described in the C specification.
Instead, they allocate and initialize pointers (and pointers to 
pointers) so that, for example, a[i2][i1] behaves like a 2-D array.

The array dimensions are numbered, which makes it easy to add 
functions for arrays of higher dimensions.  In particular,
the 1st dimension of length n1 is always the fastest dimension,
the 2nd dimension of length n2 is the next fastest dimension,
and so on.  Note that the 1st (fastest) dimension n1 is the 
first argument to the allocation functions defined below, but 
that the 1st dimension is the last subscript in a[i2][i1].
(This is another important difference between C and Fortran.)

The allocation of pointers to pointers implies that more storage
is required than is necessary to hold a true multi-dimensional array.
The fraction of the total storage allocated that is used to hold 
pointers is approximately 1/(n1+1).  This extra storage is unlikely
to represent a significant waste for large n1.

The functions defined below are significantly different from similar 
functions described by Press et al, 1988, Numerical Recipes in C.
In particular, the functions defined below:
	(1) Allocate arrays of arbitrary size elements.
	(2) Allocate contiguous storage for arrays.
	(3) Return NULL if allocation fails (just like malloc).
	(4) Do not provide arbitrary lower and upper bounds for arrays.

Contiguous storage enables an allocated multi-dimensional array to
be passed to a C function that expects a one-dimensional array.
For example, to allocate and zero an n1 by n2 two-dimensional array
of floats, one could use
	a = alloc2(n1,n2,sizeof(float));
	zeroFloatArray(n1*n2,a[0]);
where zeroFloatArray is a function defined as
	void zeroFloatArray(int n, float *a)
	{
		int i;
		for (i=0; i<n; i++)
			a[i] = 0.0;
	}

Internal error handling and arbitrary array bounds, if desired,
should be implemented in functions that call the functions defined 
below, with the understanding that these enhancements may limit 
portability.
**************************************************************************/

#include "alloc.h"

/* allocate a 1-d array */
void *alloc1 (size_t n1, size_t size)
{
	void *p;

	if ((p=malloc(n1*size))==NULL)
		return NULL;
	return p;
}

/* re-allocate a 1-d array */
void *realloc1(void *v, size_t n1, size_t size)
{
	void *p;

	if ((p=realloc(v,n1*size))==NULL)
		return NULL;
	return p;
}

/* free a 1-d array */
void free1 (void *p)
{
	free(p);
}

/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
	size_t i2;
	void **p;

	if ((p=(void**)malloc(n2*sizeof(void*)))==NULL) 
		return NULL;
	if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
		free(p);
		return NULL;
	}
	for (i2=0; i2<n2; i2++)
		p[i2] = (char*)p[0]+size*n1*i2;
	return p;
}

/* free a 2-d array */
void free2 (void **p)
{
	free(p[0]);
	free(p);
}

/* allocate a 3-d array */
void ***alloc3 (size_t n1, size_t n2, size_t n3, size_t size)
{
	size_t i3,i2;
	void ***p;

	if ((p=(void***)malloc(n3*sizeof(void**)))==NULL)
		return NULL;
	if ((p[0]=(void**)malloc(n3*n2*sizeof(void*)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void*)malloc(n3*n2*n1*size))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i3=0; i3<n3; i3++) {
		p[i3] = p[0]+n2*i3;
		for (i2=0; i2<n2; i2++)
			p[i3][i2] = (char*)p[0][0]+size*n1*(i2+n2*i3);
	}
	return p;
}

/* free a 3-d array */
void free3 (void ***p)
{
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* allocate a 4-d array */
void ****alloc4 (size_t n1, size_t n2, size_t n3, size_t n4, size_t size)
{
	size_t i4,i3,i2;
	void ****p;

	if ((p=(void****)malloc(n4*sizeof(void***)))==NULL)
		return NULL;
	if ((p[0]=(void***)malloc(n4*n3*sizeof(void**)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void**)malloc(n4*n3*n2*sizeof(void*)))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void*)malloc(n4*n3*n2*n1*size))==NULL) {
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i4=0; i4<n4; i4++) {
		p[i4] = p[0]+i4*n3;
		for (i3=0; i3<n3; i3++) {
			p[i4][i3] = p[0][0]+n2*(i3+n3*i4);
			for (i2=0; i2<n2; i2++)
				p[i4][i3][i2] = (char*)p[0][0][0]+
						size*n1*(i2+n2*(i3+n3*i4));
		}
	}
	return p;
}

/* free a 4-d array */
void free4 (void ****p)
{
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 5-d array */
void *****alloc5 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size)
{
	size_t i5,i4,i3,i2;
	void *****p;

	if ((p=(void*****)malloc(n5*sizeof(void****)))==NULL)
		return NULL;
	if ((p[0]=(void****)malloc(n5*n4*sizeof(void***)))==NULL) {
		free(p);
		return NULL;
	}
	if ((p[0][0]=(void***)malloc(n5*n4*n3*sizeof(void**)))==NULL) {
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void**)malloc(n5*n4*n3*n2*sizeof(void*)))==NULL) {
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0][0]=(void*)malloc(n5*n4*n3*n2*n1*size))==NULL) {
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	for (i5=0; i5<n5; i5++) {
		p[i5] = p[0]+n4*i5;
		for (i4=0; i4<n4; i4++) {
			p[i5][i4] = p[0][0]+n3*(i4+n4*i5);
			for (i3=0; i3<n3; i3++) {
				p[i5][i4][i3] = p[0][0][0]+n2*(i3+n3*(i4+n4*i5));
				for (i2=0; i2<n2; i2++)
					p[i5][i4][i3][i2] = (char*)p[0][0][0][0]+
						size*n1*(i2+n2*(i3+n3*(i4+n4*i5)));
			}
		}
	}
	return p;
}

/* free a 5-d array */
void free5 (void *****p)
{
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* The following two functions were added by Zhaobo Meng, Jan. 1997*/
/* allocate a 6-d array */
void ******alloc6 (size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6,
                  size_t size)
{
	size_t i6,i5,i4,i3,i2;
	void ******p;

	if ((p=(void******)malloc(n6*sizeof(void*****)))==NULL)
		return NULL;

	if ((p[0]=(void*****)malloc(n6*n5*sizeof(void****)))==NULL) {
                free(p);
		return NULL;
        }

	if ((p[0][0]=(void****)malloc(n6*n5*n4*sizeof(void***)))==NULL) {
		free(p[0]);
                free(p);
		return NULL;
	}
	if ((p[0][0][0]=(void***)malloc(n6*n5*n4*n3*sizeof(void**)))==NULL) {
		free(p[0][0]);
                free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0][0]=(void**)malloc(n6*n5*n4*n3*n2*sizeof(void*)))==NULL) {
	        free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}
	if ((p[0][0][0][0][0]=(void*)malloc(n6*n5*n4*n3*n2*n1*size))==NULL) {
	        free(p[0][0][0][0]);
		free(p[0][0][0]);
		free(p[0][0]);
		free(p[0]);
		free(p);
		return NULL;
	}

        for (i6=0; i6<n6; i6++) {
                p[i6] = p[0]+n5*i6;
	        for (i5=0; i5<n5; i5++) {
	                p[i6][i5] = p[0][0]+n4*(i5+n5*i6);
			for (i4=0; i4<n4; i4++) {
			        p[i6][i5][i4] = p[0][0][0]+n3*(i4+n4*(i5+n5*i6));
				for (i3=0; i3<n3; i3++) {
				        p[i6][i5][i4][i3] = p[0][0][0][0]
					      +n2*(i3+n3*(i4+n4*(i5+n5*i6)));
					for (i2=0; i2<n2; i2++)
					        p[i6][i5][i4][i3][i2] = 
						      (char*)p[0][0][0][0][0]+
				       size*n1*(i2+n2*(i3+n3*(i4+n4*(i5+n5*i6))));
			        }
		        }
	        }
        }
	return p;
}

/* free a 6-d array */
void free6 (void ******p)
{
        free(p[0][0][0][0][0]);
	free(p[0][0][0][0]);
	free(p[0][0][0]);
	free(p[0][0]);
	free(p[0]);
	free(p);
}

/* allocate a 1-d array of ints */
int *alloc1int(size_t n1)
{
	return (int*)alloc1(n1,sizeof(int));
}

/* re-allocate a 1-d array of ints */
int *realloc1int(int *v, size_t n1)
{
	return (int*)realloc1(v,n1,sizeof(int));
}

/* free a 1-d array of ints */
void free1int(int *p)
{
	free1(p);
}

/* allocate a 2-d array of ints */
/*  n1: fast dimension; n2: slow dimension */
int **alloc2int(size_t n1, size_t n2)
{
	return (int**)alloc2(n1,n2,sizeof(int));
}

/* free a 2-d array of ints */
void free2int(int **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of ints */
int ***alloc3int(size_t n1, size_t n2, size_t n3)
{
	return (int***)alloc3(n1,n2,n3,sizeof(int));
}

/* free a 3-d array of ints */
void free3int(int ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of floats */
float *alloc1float(size_t n1)
{
	return (float*)alloc1(n1,sizeof(float));
}

/* re-allocate a 1-d array of floats */
float *realloc1float(float *v, size_t n1)
{
	return (float*)realloc1(v,n1,sizeof(float));
}

/* free a 1-d array of floats */
void free1float(float *p)
{
	free1(p);
}

/* allocate a 2-d array of floats */
/*  n1: fast dimension; n2: slow dimension */
float **alloc2float(size_t n1, size_t n2)
{
	return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of floats */
float ***alloc3float(size_t n1, size_t n2, size_t n3)
{
	return (float***)alloc3(n1,n2,n3,sizeof(float));
}

/* free a 3-d array of floats */
void free3float(float ***p)
{
	free3((void***)p);
}

/* allocate a 4-d array of floats, added by Zhaobo Meng, 1997 */
float ****alloc4float(size_t n1, size_t n2, size_t n3, size_t n4)
{
        return (float****)alloc4(n1,n2,n3,n4,sizeof(float));
}

/* free a 4-d array of floats, added by Zhaobo Meng, 1997 */
void free4float(float ****p)
{
        free4((void****)p);
}

/* allocate a 5-d array of floats, added by Zhaobo Meng, 1997 */
float *****alloc5float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
        return (float*****)alloc5(n1,n2,n3,n4,n5,sizeof(float));
}

/* free a 5-d array of floats, added by Zhaobo Meng, 1997 */
void free5float(float *****p)
{
        free5((void*****)p);
}

/* allocate a 6-d array of floats, added by Zhaobo Meng, 1997 */
float ******alloc6float(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6)
{
        return (float******)alloc6(n1,n2,n3,n4,n5,n6,sizeof(float));
}

/* free a 6-d array of floats, added by Zhaobo Meng, 1997 */
void free6float(float ******p)
{
        free6((void******)p);
}

/* allocate a 4-d array of ints, added by Zhaobo Meng, 1997 */
int ****alloc4int(size_t n1, size_t n2, size_t n3, size_t n4)
{
        return (int****)alloc4(n1,n2,n3,n4,sizeof(int));
}

/* free a 4-d array of ints, added by Zhaobo Meng, 1997 */
void free4int(int ****p)
{
        free4((void****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
int *****alloc5int(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
        return (int*****)alloc5(n1,n2,n3,n4,n5,sizeof(int));
}

/* free a 5-d array of ints, added by Zhaobo Meng, 1997 */
void free5int(int *****p)
{
        free5((void*****)p);
}

/* allocate a 5-d array of chars, added by Zhaobo Meng, 1997 */
unsigned char *****alloc5uchar(size_t n1, size_t n2, size_t n3, 
	size_t n4, size_t n5)
{
        return (unsigned char*****)alloc5(n1,n2,n3,n4,n5,sizeof(unsigned char));
}

/* free a 5-d array of chars, added by Zhaobo Meng, 1997 */
void free5uchar(unsigned char *****p)
{
        free5((void*****)p);
}

/* allocate a 5-d array of ints, added by Zhaobo Meng, 1997 */
unsigned short *****alloc5ushort(size_t n1, size_t n2, size_t n3,
        size_t n4, size_t n5)
{
        return (unsigned short*****)alloc5(n1,n2,n3,n4,n5,sizeof(unsigned short));
}
/* allocate a 3-d array of ints, added by Meng, 1997 */
unsigned short ***alloc3ushort(size_t n1, size_t n2, size_t n3)
{
        return (unsigned short***)alloc3(n1,n2,n3,sizeof(unsigned short));
}

/* allocate a 2-d array of ints, added by Meng, 1997 */
unsigned short **alloc2ushort(size_t n1, size_t n2)
{
        return (unsigned short**)alloc2(n1,n2,sizeof(unsigned short));
}

/* free a 5-d array of shorts, added by Zhaobo Meng, 1997 */
void free5ushort(unsigned short *****p)
{
        free5((void*****)p);
}
/* free a 3-d array of shorts, added by Zhaobo Meng, 1997 */
void free3ushort(unsigned short ***p)
{
        free3((void***)p);
}

/* free a 2-d array of shorts, added by Zhaobo Meng, 1997 */
void free2ushort(unsigned short **p)
{
        free2((void**)p);
}
/* allocate a 6-d array of ints, added by Zhaobo Meng, 1997 */
unsigned short ******alloc6ushort(size_t n1, size_t n2, size_t n3,
        size_t n4, size_t n5, size_t n6)
{
        return (unsigned short******)alloc6(n1,n2,n3,n4,n5,n6,sizeof(unsigned short));
}

/* free a 6-d array of shorts, added by Zhaobo Meng, 1997 */
void free6ushort(unsigned short ******p)
{
        free6((void******)p);
}

/* allocate a 1-d array of doubles */
double *alloc1double(size_t n1)
{
	return (double*)alloc1(n1,sizeof(double));
}

/* re-allocate a 1-d array of doubles */
double *realloc1double(double *v, size_t n1)
{
	return (double*)realloc1(v,n1,sizeof(double));
}


/* free a 1-d array of doubles */
void free1double(double *p)
{
	free1(p);
}

/* allocate a 2-d array of doubles */
/*  n1: fast dimension; n2: slow dimension */
double **alloc2double(size_t n1, size_t n2)
{
	return (double**)alloc2(n1,n2,sizeof(double));
}

/* free a 2-d array of doubles */
void free2double(double **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of doubles */
double ***alloc3double(size_t n1, size_t n2, size_t n3)
{
	return (double***)alloc3(n1,n2,n3,sizeof(double));
}

/* free a 3-d array of doubles */
void free3double(double ***p)
{
	free3((void***)p);
}

/* allocate a 1-d array of complexs */
complex *alloc1complex(size_t n1)
{
	return (complex*)alloc1(n1,sizeof(complex));
}

/* re-allocate a 1-d array of complexs */
complex *realloc1complex(complex *v, size_t n1)
{
	return (complex*)realloc1(v,n1,sizeof(complex));
}

/* free a 1-d array of complexs */
void free1complex(complex *p)
{
	free1(p);
}

/* allocate a 2-d array of complexs */
/*  n1: fast dimension; n2: slow dimension */
complex **alloc2complex(size_t n1, size_t n2)
{
	return (complex**)alloc2(n1,n2,sizeof(complex));
}

/* free a 2-d array of complexs */
void free2complex(complex **p)
{
	free2((void**)p);
}

/* allocate a 3-d array of complexs */
complex ***alloc3complex(size_t n1, size_t n2, size_t n3)
{
	return (complex***)alloc3(n1,n2,n3,sizeof(complex));
}

/* free a 3-d array of complexs */
void free3complex(complex ***p)
{
	free3((void***)p);
}

/* allocate a 4-d array of complexs */
complex ****alloc4complex(size_t n1, size_t n2, size_t n3, size_t n4)
{
	return (complex****)alloc4(n1,n2,n3,n4,sizeof(complex));
}

/* allocate a 5-d array of complexs */
complex *****alloc5complex(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
	return (complex*****)alloc5(n1,n2,n3,n4,n5,sizeof(complex));
}

/* free a 4-d array of complexs */
void free4complex(complex ****p)
{
	free4((void****)p);
}
/* free a 5-d array of complexs */
void free5complex(complex *****p)
{
	free5((void*****)p);
}
void zero1int(int *p, size_t n1)
{
     int i;
     for(i=0;i<n1;i++) p[i]=0;
}
void zero2int(int **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0;
}
void zero3int(int ***p, size_t n1, size_t n2, size_t n3)
{
     int i, j, k;
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0;
}
void zero3ushort(unsigned short ***p, size_t n1, size_t n2, size_t n3)
{
     int i, j, k;
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0;
}
void zero2ushort(unsigned short **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
            p[i][j]=0;
}
void zero1float(float *p, size_t n1)
{
     int i;
     for(i=0;i<n1;i++) p[i]=0.0;
}
void zero2float(float **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0.0;
}
void zero3float(float ***p, size_t n1, size_t n2, size_t n3)
{
     int i, j, k;
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0.0;
}
void zero4float(float ****p, size_t n1, size_t n2, size_t n3, size_t n4)
{
     int m, i, j, k;
     for(m=0;m<n4;m++) 
       for(i=0;i<n3;i++) 
          for(j=0;j<n2;j++)
            for(k=0;k<n1;k++)
               p[m][i][j][k]=0.0;
}
void zero1double(double *p, size_t n1)
{
     int i;
     for(i=0;i<n1;i++) p[i]=0.0;
}
void zero2double(double **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=0.0;
}
void zero3double(double ***p, size_t n1, size_t n2, size_t n3)
{
     int i, j, k;
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=0.0;
}
void zero1complex(complex *p, size_t n1)
{
     int i;
     for(i=0;i<n1;i++) p[i]=cmplx(0.0, 0.0);
}
void zero2complex(complex **p, size_t n1, size_t n2)
{
     int i, j;
     for(i=0;i<n2;i++) 
       for(j=0;j<n1;j++)
         p[i][j]=cmplx(0.0, 0.0);
}
void zero3complex(complex ***p, size_t n1, size_t n2, size_t n3)
{
     int i, j, k;
     for(i=0;i<n3;i++) 
       for(j=0;j<n2;j++)
          for(k=0;k<n1;k++)
            p[i][j][k]=cmplx(0.0, 0.0);
}
void zero4complex(complex ****p, size_t n1, size_t n2, size_t n3, size_t n4)
{
     int m, i, j, k;

     for(m=0;m<n4;m++) 
       for(i=0;i<n3;i++) 
         for(j=0;j<n2;j++)
           for(k=0;k<n1;k++)
              p[m][i][j][k]=cmplx(0.0, 0.0);
}
void zero5complex(complex *****p, size_t n1, size_t n2, size_t n3, size_t n4, size_t n5)
{
     int l, m, i, j, k;

     for(l=0;l<n5;l++) 
       for(m=0;m<n4;m++) 
         for(i=0;i<n3;i++) 
           for(j=0;j<n2;j++)
             for(k=0;k<n1;k++)
                p[l][m][i][j][k]=cmplx(0.0, 0.0);
}
/**************************************************************/
#ifdef TEST1
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
	hv = (short *) alloc1(n1, sizeof(short));
	iv = alloc1int(n1);
	fv = alloc1float(n1);
	dv = alloc1double(n1);

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
	hm = (short *) alloc2(n1, n2, sizeof(short));
	im = alloc2int(n1, n2);
	fm = alloc2float(n1, n2);
	dm = alloc2double(n1, n2);

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

	free2(hm);
	free2int(im);
	free2float(fm);
	free2double(dm);

	exit(0);
}
#endif
