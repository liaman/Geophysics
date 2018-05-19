/* fft.h  -- FFT transform head file */
/**********************************************************************
PFAFFT - Functions to perform Prime Factor (PFA) FFT's, in place

npfa		return valid n for complex-to-complex PFA
npfar		return valid n for real-to-complex/complex-to-real PFA
npfao		return optimal n for complex-to-complex PFA
npfaro		return optimal n for real-to-complex/complex-to-real PFA
pfacc		1D PFA complex to complex
pfacr		1D PFA complex to real
pfarc		1D PFA real to complex
pfamcc		multiple PFA complex to real
pfa2cc		2D PFA complex to complex
pfa2cr		2D PFA complex to real
pfa2rc		2D PFA real to complex
**********************************************************************/
#ifndef FFT_H
#define FFT_H

#include "cwp.h"
#include "complex.h"

/* Prime Factor FFTs */
int npfa (int nmin);
int npfao (int nmin, int nmax);
int npfar (int nmin);
int npfaro (int nmin, int nmax);
void pfacc (int isign, int n, complex z[]);
void pfarc (int isign, int n, float rz[], complex cz[]);
void pfacr (int isign, int n, complex cz[], float rz[]);
void pfa2cc (int isign, int idim, int n1, int n2, complex z[]);
void pfa2rc (int isign, int idim, int n1, int n2, float rz[], complex cz[]);
void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], float rz[]);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);

#endif
