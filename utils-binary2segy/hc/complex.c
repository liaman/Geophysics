/* Copyright (c) Colorado School of Mines, 2005.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
COMPLEX - Functions to manipulate complex numbers

cadd	add two complex numbers
csub	subtract two complex numbers
cmul	multiply two complex numbers
cdiv	divide two complex numbers
cmplx	make a complex number from two real numbers
conjg	complex conjugate of a complex number 
cneg	negate a complex number
cinv	invert a complex number
csqrt	complex square root of a complex number
cexp	complex exponential of a complex number
crmul	multiply a complex number by a real number 
rcabs	real magnitude of a complex number

******************************************************************************
Structure:
typedef struct _complexStruct {  complex number
	float r,i;
} complex;

******************************************************************************
Function Prototypes:
complex cadd (complex a, complex b);
complex csub (complex a, complex b);
complex cmul (complex a, complex b);
complex cdiv (complex a, complex b);
float rcabs (complex z);
complex cmplx (float re, float im);
complex conjg (complex z);
complex cneg (complex z);
complex cinv (complex z);
complex csqrt (complex z);
complex cexp (complex z);
complex crmul (complex a, float x);

******************************************************************************
Notes:
The function "rcabs" was originally called "fcabs". This produced
a collision on some systems so a new name was chosen.

******************************************************************************
Reference:
Adapted from Press et al, 1988, Numerical Recipes in C (Appendix E).

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
Modified:  Dave Hale, Colorado School of Mines, 04/26/90
	Added function cinv().
*****************************************************************************/
/**************** end self doc ********************************/

#include "complex.h"

complex cadd(complex a, complex b)
{
	complex c;
	c.r = a.r+b.r;
	c.i = a.i+b.i;
	return c;
}

complex csub(complex a, complex b)
{
	complex c;
	c.r = a.r-b.r;
	c.i = a.i-b.i;
	return c;
}

complex cmul(complex a, complex b)
{
	complex c;
	c.r = a.r*b.r-a.i*b.i;
	c.i = a.i*b.r+a.r*b.i;
	return c;
}

complex cdiv(complex a, complex b)
{
	complex c;
	float r,den;
	if (fabs(b.r)>=fabs(b.i)) {
		r = b.i/b.r;
		den = b.r+r*b.i;
		c.r = (a.r+r*a.i)/den;
		c.i = (a.i-r*a.r)/den;
	} else {
		r = b.r/b.i;
		den = b.i+r*b.r;
		c.r = (a.r*r+a.i)/den;
		c.i = (a.i*r-a.r)/den;
	}
	return c;
}

complex cmplx(float re, float im)
{
	complex c;
	c.r = re;
	c.i = im;
	return c;
}

complex conjg(complex z)
{
	complex c;
	c.r = z.r;
	c.i = -z.i;
	return c;
}

complex cneg(complex z)
{
	complex c;
	c.r = -z.r;
	c.i = -z.i;
	return c;
}

complex cinv(complex z)
{
	complex c;
	float s;
	s = 1.0/(z.r*z.r+z.i*z.i);
	c.r = z.r*s;
	c.i = -z.i*s;
	return c;
}

complex ccsqrt(complex z)//+++++++=change the "c" to "cc"
{
	complex c;
	float x,y,w,r;
	if (z.r==0.0 && z.i==0.0) {
		c.r = c.i = 0.0;
		return c;
	} else {
		x = fabs(z.r);
		y = fabs(z.i);
		if (x>=y) {
			r = y/x;
			w = sqrt(x)*sqrt(0.5*(1.0+sqrt(1.0+r*r)));
		} else {
			r = x/y;
			w = sqrt(y)*sqrt(0.5*(r+sqrt(1.0+r*r)));
		}
		if (z.r>=0.0) {
			c.r = w;
			c.i = z.i/(2.0*w);
		} else {
			c.i = (z.i>=0.0) ? w : -w;
			c.r = z.i/(2.0*c.i);
		}
		return c;
	}
}

complex ccexp(complex z)//+++++++=change the "c" to "cc"
{
	float a;
	complex c;
	a = exp(z.r);
	c.r = a*cos(z.i);
	c.i = a*sin(z.i);
	return c;
}

complex crmul(complex a, float x)
{
	complex c;
	c.r = x*a.r;
	c.i = x*a.i;
	return c;
}

float rcabs(complex z)
{
	float x,y,ans,temp;
	x = fabs(z.r);
	y = fabs(z.i);
	if (x==0.0)
		ans = y;
	else if (y==0.0)
		ans = x;
	else if (x>y) {
		temp = y/x;
		ans = x*sqrt(1.0+temp*temp);
	} else {
		temp =x/y;
		ans = y*sqrt(1.0+temp*temp);
	}
	return ans;
}

