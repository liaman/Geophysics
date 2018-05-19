/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/******************************************************************************
Xtcwp.h:  header file for X cwp library
*******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 08/28/90
******************************************************************************/

#ifndef XTCWP_H
#define XTCWP_H

/* INCLUDES */
#include "par.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/Xatom.h>
#include <X11/Intrinsic.h>
#include <X11/IntrinsicP.h>
#include <X11/StringDefs.h>
#include <X11/CoreP.h>


/* DATA STRUCTURES */
typedef struct {
	GC gc;
	float xshift,yshift;		/* x and y shifts */
	float xscale,yscale;		/* x and y scale factors */
	int clip;			/* clip if non-zero */
	float xmin,ymin,xmax,ymax;	/* clip rectangle */
} *FGC;
typedef struct {
	float fx,fy;
} FXPoint;
typedef struct {
	float fx,fy;
	float fwidth,fheight;
} FXRectangle;

/* MACROS */
#define MapFX(fgc,fx) ((int)(((fgc)->xshift)+(fx)*((fgc)->xscale)))
#define MapFY(fgc,fy) ((int)(((fgc)->yshift)+(fy)*((fgc)->yscale)))
#define MapFWidth(fgc,fwidth) ((int)((fwidth)*((fgc)->xscale)))
#define MapFHeight(fgc,fheight) ((int)((fheight)*((fgc)->yscale)))
#define MapFAngle(fgc,fangle) ((int)((fangle)*64.0))
#define MapX(fgc,x) (((x)-((fgc)->xshift))/((fgc)->xscale))
#define MapY(fgc,y) (((y)-((fgc)->yshift))/((fgc)->yscale))
#define MapWidth(fgc,width) ((width)/((fgc)->xscale))
#define MapHeight(fgc,height) ((height)/((fgc)->yscale))
#define MapAngle(fgc,angle) ((angle)/64.0)

/* RESOURCES */
#define XtcwpNaxesGrid "axesGrid"
#define XtcwpCAxesGrid "AxesGrid"
#define XtcwpRAxesGrid "AxesGrid"
#define XtcwpNONE 0
#define XtcwpDOT 1
#define XtcwpDASH 2
#define XtcwpSOLID 3
#define XtcwpNaxesStyle "axesStyle"
#define XtcwpCAxesStyle "AxesStyle"
#define XtcwpRAxesStyle "AxesStyle"
#define XtcwpNORMAL 0
#define XtcwpSEISMIC 1
#define XtcwpRFloat "Float"

/* CALLBACK REASONS */
#define XtcwpCR_RESIZE 1
#define XtcwpCR_EXPOSE 2
#define XtcwpCR_INPUT 3

/* FUNCTION PROTOTYPES */
int FMapFX (FGC fgc, float fx);
int FMapFY (FGC fgc, float fy);
int FMapFWidth (FGC fgc, float fwidth);
int FMapFHeight (FGC fgc, float fheight);
int FMapFAngle (FGC fgc, float fangle);
void FMapFPoint (FGC fgc, float fx, float fy, int *x_return, int *y_return);
void FMapFPoints (FGC fgc, FXPoint fpoints[], int npoints, 
	XPoint points_return[]);
float FMapX (FGC fgc, int x);
float FMapY (FGC fgc, int y);
float FMapWidth (FGC fgc, int width);
float FMapHeight (FGC fgc, int height);
float FMapAngle (FGC fgc, int angle);
void FMapPoint (FGC fgc, int x, int y, float *fx_return, float *fy_return);
void FMapPoints (FGC fgc, XPoint points[], int npoints, 
	FXPoint fpoints_return[]);
void FSetGC (FGC fgc, GC gc);
void FSetMap (FGC fgc, int x, int y, int width, int height,
	float fx, float fy, float fwidth, float fheight);
void FSetClipRectangle(FGC fgc, float fxa, float fya, float fxb, float fyb);
void FClipOn (FGC fgc);
void FClipOff (FGC fgc);
int FClipPoint (FGC fgc, float fx, float fy);
int FClipLine (FGC fgc, float fx1, float fy1, float fx2, float fy2,
	float *fx1c, float *fy1c, float *fx2c, float *fy2c);
int FClipRectangle (FGC fgc, float fx, float fy, float fwidth, float fheight,
	float *fxc, float *fyc, float *fwidthc, float *fheightc);
FGC FXCreateFGC (GC gc, int x, int y, int width, int height,
	float fx, float fy, float fwidth, float fheight);
void FXFreeFGC (FGC fgc);
void FXDrawPoint (Display *display, Drawable d, FGC fgc, float fx, float fy);
void FXDrawPoints (Display *display, Drawable d, FGC fgc, 
	FXPoint fpoints[], int npoints, int mode);
void FXDrawLine (Display *display, Drawable d, FGC fgc,
	float fx1, float fy1, float fx2, float fy2);
void FXDrawLines (Display *display, Drawable d, FGC fgc,
	FXPoint fpoints[], int npoints, int mode);
void FXDrawRectangle (Display *display, Drawable d, FGC fgc, 
	float fx, float fy, float fwidth, float fheight);
void FXDrawArc (Display *display, Drawable d, FGC fgc,
	float fx, float fy, float fwidth, float fheight, 
	float fangle1, float fangle2);
void FXDrawString (Display *display, Drawable d, FGC fgc, 
	float fx, float fy, char *string, int length);
void FXFillRectangle (Display *display, Drawable d, FGC fgc, 
	float fx, float fy, float fwidth, float fheight);
void XtcwpStringToFloat (XrmValue *args, int *nargs, 
	XrmValue *fromVal, XrmValue *toVal);
Status XtcwpCreateRGBDefaultMap (Display *dpy, XStandardColormap *scmap);
unsigned long XtcwpGetFirstPixel (Display *dpy);
unsigned long XtcwpGetLastPixel (Display *dpy);
Colormap XtcwpCreateRGBColormap (Display *dpy, Window win);
Colormap XtcwpCreateGrayColormap (Display *dpy, Window win);
Colormap XtcwpCreateHueColormap (Display *dpy, Window win,
	float fhue, float lhue, float sat, float bright);
Colormap XtcwpCreateSatColormap (Display *dpy, Window win,
	float fhue, float lhue, float wfrac, float bright);
void XtcwpRubberBox (Display *dpy, Window win, XEvent event,
	int *x, int *y, int *width, int *height);
void XtcwpDrawString90 (Display *dpy, Drawable d, GC gc,
	int x, int y, char *string, int count);

#endif /* XTCWP_H */
