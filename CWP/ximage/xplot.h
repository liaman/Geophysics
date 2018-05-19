/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* include file for X graphics */

#ifndef XPLOT_H
#define XPLOT_H


/* INCLUDES */
#include "par.h"      /* required for portability, do not remove */
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>

/* DEFINES */

/* axes drawing */
#define NONE 0
#define DOT 1
#define DASH 2
#define SOLID 3
#define NORMAL 0
#define SEISMIC 1

/* 256 pixel values for truecolor model*/
extern unsigned long truecolor_pixel[256];


/* FUNCTION PROTOTYPES */

/* windows */
Window xNewWindow (Display *dpy, int x, int y, int width, int height,
	int border, int background, char *name);

/* axes drawing */
void xDrawAxesBox (Display *dpy, Window win,
	int x, int y, int width, int height,
	float x1beg, float x1end, float p1beg, float p1end,
	float d1num, float f1num, int n1tic, int grid1, char *label1,
	float x2beg, float x2end, float p2beg, float p2end,
	float d2num, float f2num, int n2tic, int grid2, char *label2,
	char *labelfont, char *title, char *titlefont, 
	char *axescolor, char *titlecolor, char *gridcolor,
	int style);
void xSizeAxesBox (Display *dpy, Window win, 
	char *labelfont, char *titlefont, int style,
	int *x, int *y, int *width, int *height);

/* images */
XImage *xNewImage (Display *dpy, unsigned long pmin, unsigned long pmax,
	int width, int height, float blank, unsigned char *bytes);

/* rubberbanding box */
void xRubberBox (Display *dpy, Window win, XEvent event,
	int *x, int *y, int *width, int *height);

/* colormaps */
Status xCreateRGBDefaultMap (Display *dpy, XStandardColormap *scmap);
unsigned long xGetFirstPixel (Display *dpy);
unsigned long xGetLastPixel (Display *dpy);
Colormap xCreateRGBColormap (Display *dpy, Window win,
	 char *str_cmap, int verbose);
Colormap xCreateHSVColormap (Display *dpy, Window win,
	 char *str_cmap, int verbose);
Colormap xCreateGrayColormap (Display *dpy, Window win);
Colormap xCreateHueColormap (Display *dpy, Window win);
void xDrawLegendBox(Display *dpy, Window win,
        int x, int y, int width, int height,
        float bclip, float wclip, char *units, char *legendfont,
        char *labelfont, char *title, char *titlefont,
        char *axescolor, char *titlecolor, char *gridcolor,
        int style);

/* for xcontour */
void xContour(Display *dpy, Window win,GC gcc, GC gcl, 
	       float *cp,int nx, float x[], int ny, float y[], float z[], 
	       char lcflag,char *lcf,char *lcc, float *w, int nplaces);

/* curve drawing */
void xDrawCurve(Display *dpy, Window win,
		int x, int y, int width, int height,
		float x1beg, float x1end, float p1beg, float p1end,
		float x2beg, float x2end, float p2beg, float p2end,
		float *x1curve, float *x2curve, int ncurve,
		char *curvecolor, int style);

#endif /* XPLOT_H */
