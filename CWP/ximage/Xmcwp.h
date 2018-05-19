/* INCLUDES */

#if (XT_REVISION >= 5)
#include <X11/Xlocale.h>
#endif
#include <X11/Intrinsic.h>
#include <X11/cursorfont.h>
#include <Xm/Xm.h>

#if (XT_REVISION >= 5)
#include <Xm/XmStrDefs.h>
#endif

#include <Xm/MainW.h>
#include <Xm/Frame.h>
#include <Xm/Form.h>
#include <Xm/LabelG.h>
#include <Xm/DrawingA.h>
#include <Xm/PushBG.h>
#include <Xm/ToggleBG.h>
#include "Xtcwp/Xtcwp.h"


/* TYPEDEFS */

typedef enum {
	DRAW,	/* sample nearest mouse is set to mouse position */
	NEGATE, /* sample nearest mouse is negated */
	ZERO,	/* samples nearest mouse is zeroed */
	NONE	/* editing disabled */
} EditMode;

typedef struct SamplesStruct {
	Widget frame; /* XmFrame contains label and drawing area */
	Widget da; /* XmDrawingArea in which samples are drawn */
	GC gcDraw,gcErase; /* graphics contexts for drawing and erasing */
	int n; /* number of samples */
	float *data; /* sample values */
	int origin; /* index of samples origin (0 to n-1) */
	float plotValue; /* sample value plotted at half-maximum */
	int editMode; /* determines editing of samples */
	void (*editDone)(struct SamplesStruct *s); /* called when edit done */
	float width,base,scale,radius; /* parameters for drawing samples */
} Samples;


/* FUNCTION PROTOTYPES */

Samples *samplesCreate (Widget parent, char *title,
	void (*editDone)(Samples *s));	
void samplesDraw (Samples *s);
void samplesSetN (Samples *s, int n);
void samplesSetData (Samples *s, float *d);
void samplesSetPlotValue (Samples *s, float pv);
void samplesSetEditMode (Samples *s, EditMode m);
void samplesSetOrigin (Samples *s, int i);

/* MISCELLANEOUS DEFINES */
#ifndef CHARSET
#define CHARSET ((XmStringCharSet)XmSTRING_DEFAULT_CHARSET)
#endif

/* FUNCTION PROTOTYPES */
Widget XtcwpCreateStringRadioButtons (Widget parent, char *label,
	int nstrings, char **strings, int first,
	void (*callback)(int selected, void *clientdata), void *clientdata);

