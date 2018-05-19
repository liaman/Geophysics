/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/******************************************************************************
AxesP.h:  Private header file for Axes Widget
*******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 08/28/90
******************************************************************************/

#ifndef AXESP_H
#define AXESP_H

#include "Xtcwp/Axes.h"

typedef struct _XtcwpAxesClassPart {
	int ignore;
} XtcwpAxesClassPart;

typedef struct _XtcwpAxesClassRec {
	CoreClassPart core_class;
	XtcwpAxesClassPart axes_class;
} XtcwpAxesClassRec;

extern XtcwpAxesClassRec XtcwpaxesClassRec;

typedef struct _XtcwpAxesPart {
	int style;			/* normal or seismic */
	Position x,y;			/* axes box upper-left corner */
	Dimension width,height;		/* axes box dimensions */
	float x1beg,x1end;		/* axes values for dimension 1 */
	float x2beg,x2end;		/* axes values for dimension 2 */
	float p1beg,p1end;		/* axes pads for dimension 1 */
	float p2beg,p2end;		/* axes pads for dimension 2 */
	int grid1,grid2;		/* none, dot, dash, or solid */
	int n1tic,n2tic;		/* tics per numbered tic */
	char *label1,*label2;		/* axes labels */
	char *title;			/* axes title */
	Pixel axescolor;		/* for box, tics, and labels */
	Pixel gridcolor;		/* for grid lines */
	Pixel titlecolor;		/* for grid lines */
	Font labelfont;			/* font for axes labels */
	Font titlefont;			/* font for title */
	XtCallbackList resize;		/* callback list */
	XtCallbackList expose;		/* callback list */
	XtCallbackList input;		/* callback list */
} XtcwpAxesPart;

typedef struct _XtcwpAxesRec {
   CorePart core;
   XtcwpAxesPart axes;
} XtcwpAxesRec;

#endif /* AXESP_H */
