#ifndef COMPLEX_H
#define COMPLEX_H

/* TYPEDEFS */
#ifdef CRAY
typedef struct _complexStruct { /* complex number */
	float r,i;
}  cwp_complex;
typedef struct _dcomplexStruct { /* double-precision complex number */
	double r,i;
}  cwp_dcomplex;
#define complex cwp_complex
#define dcomplex cwp_dcomplex
#define cadd cwp_cadd
#define csub cwp_csub
#define cmul cwp_cmul
#define cdiv cwp_cdiv
#define rcabs cwp_rcabs
#define cmplx cwp_cmplx
#define conjg cwp_conjg
#define cneg cwp_cneg
#define cinv cwp_cinv
#define csqrt cwp_csqrt
#define cexp cwp_cexp
#define crmul cwp_crmul
#define cipow cwp_cipow
#define crpow cwp_crpow
#define rcpow cwp_rcpow
#define ccpow cwp_ccpow
#define ccos cwp_ccos
#define csin cwp_csin
#define ccosh cwp_ccosh
#define csinh cwp_csinh
#define cexp1 cwp_cexp1
#define clog cwp_clog

#else
 
  #ifndef __cplusplus /* if not C++, define the C struct complex */
     #ifndef complex
     typedef struct _complexStruct 
     { /* complex number */
   	     float r,i;
     } complex;
     #endif/* complex */
  #else /* if C++, define the C++ class complex */
     #include "Complex.h"
  #endif /* C++ */

#endif

/* FUNCTION PROTOTYPES */

#ifndef __cplusplus /* if not C++, declare C complex functions */
/* complex number manipulation */
complex cadd (complex a, complex b);
complex csub (complex a, complex b);
complex cmul (complex a, complex b);
complex cdiv (complex a, complex b);
float rcabs (complex z);
complex cmplx (float re, float im);
complex conjg (complex z);
complex cneg (complex z);
complex cinv (complex z);
complex ccsqrt (complex z); //+++++++=change the "c" to "cc"
complex ccexp (complex z);  //+++++++=change the "c" to "cc"
complex crmul (complex a, float x);

/* complex functions */
complex cipow(complex a, int p);
complex crpow(complex a, float p);
complex rcpow(float a, complex p);
complex ccpow (complex a, complex p);
complex cccos(complex a);//+++++++=change the "c" to "cc"
complex ccsin(complex a);//+++++++=change the "c" to "cc"
complex cccosh(complex a);//+++++++=change the "c" to "cc"
complex ccsinh(complex a);//+++++++=change the "c" to "cc"
complex cexp1(complex a);
complex cclog(complex a);//+++++++=change the "c" to "cc"
#endif

#endif /* Complex.h */
