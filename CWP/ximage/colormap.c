/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* COLORMAP: $Revision: 1.19 $ ; $Date: 2011/11/21 17:02:44 $	*/

/*********************** self documentation **********************/
/*****************************************************************************
COLORMAP - Functions to manipulate X colormaps:

xCreateRGBDefaultMap	create XA_RGB_DEFAULT_MAP property of root window if
			it does not already exist
xGetFirstPixel		return first pixel in range of contiguous pixels in
			XA_RGB_DEFAULT_MAP
xGetLastPixel		return last pixel in range of contiguous pixels in
			XA_RGB_DEFAULT_MAP
xCreateHSVColormap	create a 2 ramp colormap (HSV - Model)
xCreateRGBColormap	create a 2 ramp colormap (RGB - Model)

******************************************************************************
Function Prototypes:
Status xCreateRGBDefaultMap (Display *dpy, XStandardColormap *scmap);
unsigned long xGetFirstPixel (Display *dpy);
unsigned long xGetLastPixel (Display *dpy);
Colormap xCreateRGBColormap (Display *dpy, Window win,
	char *str_cmap, int verbose)
Colormap xCreateHSVColormap (Display *dpy, Window win,
	char *str_cmap, int verbose)

******************************************************************************
xCreateRGBDefaultMap:
Input:
dpy		display

Output:
scmap		the standard colormap structure

******************************************************************************
xGetFirstPixel, xGetLastPixel:
Input:
dpy		display

******************************************************************************
Notes:
PROBLEM
-------

Most mid-range display devices today support what X calls
the "PseudoColor visual".  Typically, only 256 colors (or gray
levels) may be displayed simultaneously.  Although these 256 colors
may be chosen from a much larger (4096 or more) set of available
colors, only 256 colors can appear on a display at one time.

These 256 colors are indexed by pixel values in a table called
the colormap.  Each window can have its own colormap, but only
one colormap can be installed in the display hardware at a time.
(Again, only 256 colors may be displayed at one time.)  The window 
manager is responsible for installing a window's colormap when that 
window becomes the key window.

Many of the applications we are likely to write require a large,
contiguous range of pixels (entries in the colormap).  In this
range, we must be able to:
(1) given a color (or gray), determine the corresponding pixel.
(2) given a pixel, determine the corresponding color (or gray).
An example would be an imaging application that uses a gray scale
to display images in shades of gray between black and white.
Such applications are also likely to require a few additional colors
for drawing axes, text, etc.

The problem is to coordinate the use of the limited number of
256 simultaneous colors so that windows for different applications 
appear reasonable, even when their particular colormaps are not
installed in the display hardware.  For example, we might expect 
an analog xclock's hands to be visible even when xclock's window
is not the key window, when its colormap is not installed.

We should ensure that the range of contiguous pixels used by one
application (perhaps for imaging) does not conflict with the pixels
used by other applications to draw text, clock hands, etc.


SOLUTION
--------

Applications that do not require special colormaps should simply
use the default colormap inherited from the root window when new
top-level windows are created.

Applications that do require a special colormap MUST create their
own colormap.  They must not assume that space will be available
in the default colormap for a contiguous range of read/write pixels,
because the server or window manager may have already allocated
these pixels as read-only.  Even if sufficient pixels are available
in the default colormap, they should not be allocated by a single
application.  The default colormap should be used only for windows
requiring a limited number of typical colors, such as red, yellow, etc.

Applications that require a contiguous range of read/write pixels
should allocate these pixels in their window's private colormaps.
They should determine which contiguous pixels to allocate from 
parameters in the standard colormap XA_RGB_DEFAULT_MAP.  In particular,
the first pixel in the range of contiguous pixels should be 
	base_pixel
and the last pixel in the range should be 
	base_pixel+red_max*red_mult+green_max*green_mult+blue_max*blue_mult,
where base_pixel, red_max, etc. are members in the XStandardColormap
structure.  On an 8-bit display, this range will typically provide 216
contiguous pixels, which may be set to a gray scale, color scale, or
whatever.  This leaves 40 colors for drawing text, axes, etc.

If the XA_RGB_DEFAULT_MAP does not exist, it should be created to 
consist of various colors composed of an equal number of reds, 
greens, and blues.  For example, if 216 colors are to be allocated,
then red_max=green_max=blue_max=5, red_mult=36, green_mult=6, and
blue_mult=1.  Because of the difficulty in forcing a particular 
pixel to correspond to a particular color in read-only color cells,
these 216 colors will likely be read/write color cells unless
created by the X server.  In any case, these 216 colors should not
be modified by any application.  In creating custom colormaps, the
only use of XA_RGB_DEFAULT_MAP should be in determining which 216
pixels to allocate for contiguous pixels.

In creating a custom colormap for a window, the application should
initialize this colormap to the colors already contained in the
window's colormap, which was inherited initially from its parent.
This will ensure that typical colors already allocated by other
applications will be consistent with pixels used by the application
requiring the custom colormap.  Ideally, windows might have
different colormaps, but the only differences would be in the
range of contiguous colors used for imaging, rendering, etc.
Ideally, the pixels corresponding to colors used to draw text, 
axes, etc. would be consistent for all windows.

Unfortunately, it is impractical to maintain complete consistency 
among various private colormaps.  For example, suppose a custom
colormap is created for a window before other applications have
had the opportunity to allocate their colors from the default
colormap.  Then, when the window with the custom colormap becomes
the key window, the windows of the other applications may be 
displayed with false colors, since the colormap of the key window
may not contain the true colors.  The colors used by the other 
applications did not exist when the custom colormap was created.
One solution to this problem might be to initially allocate a
set of "common" colors in the default colormap before launching
any applications.  This will increase the likelihood that typical
colors will be consistent among various colormaps.

Functions are provided below to
(1) create the standard colormap XA_RGB_DEFAULT_MAP, if it does not exist,
(2) determine the first and last pixels in the contiguous range of pixels,
(3) create some common private colormaps 

xCreateRGBDefaultMap:
This function returns 0 if the XA_RGB_DEFAULT_MAP property does not exist
and cannot be created.  At least 8 contiguous color cells must be free
in the default colormap to create the XA_RGB_DEFAULT_MAP.  If created, the
red_max, green_max, and blue_max values returned in scmap will be equal.

xGetFirstPixel, xGetLastPixel:
If it does not already exist, XA_RGB_DEFAULT_MAP will be created.
If XA_RGB_DEFAULT_MAP does not exist and cannot be created, then
this function returns 0.

xCreateRGBColormap, xCreateHSVColormap:
The returned colormap is only created; the window's colormap attribute
is not changed, and the colormap is not installed by this function.
The returned colormap is a copy of the window's current colormap, but 
with an RGB color scale allocated in the range of contiguous cells
determined by XA_RGB_DEFAULT_MAP.  If it does not already exist,
XA_RGB_DEFAULT_MAP will be created.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/30/90
*****************************************************************************/
/**************** end self doc ********************************/

#include "xplot.h"
#include "cwpcmaps.h"

unsigned long truecolor_pixel[256];

/* functions defined and used internally */
static float rgbvalue (float n1, float n2, float hue);
static void hsv2rgb (float h, float s, float v, float *r, float *g, float *b);

Status xCreateRGBDefaultMap (Display *dpy, XStandardColormap *scmap)
/*****************************************************************************
create XA_RGB_DEFAULT_MAP property of root window if it does not already exist
******************************************************************************
Input:
dpy		display

Output:
scmap		the standard colormap structure
******************************************************************************
Notes:
This function returns 0 if the XA_RGB_DEFAULT_MAP property does not exist
and cannot be created.  At least 8 contiguous color cells must be free
in the default colormap to create the XA_RGB_DEFAULT_MAP.  If created, the
red_max, green_max, and blue_max values returned in scmap will be equal.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/29/90
*****************************************************************************/
{
	Screen *scr=XDefaultScreenOfDisplay(dpy);
	Window root=XRootWindowOfScreen(scr);
	Colormap cmap;
	XColor color;
	int i,ncells;
	unsigned long npixels;
	unsigned long bpixel,epixel,pixel1,pixel2,imax,rmult,gmult,bmult;
	unsigned long pixel[4096];

#ifndef  __osf__

	/* grab the server */
	XGrabServer(dpy); 

#endif

	/* if XA_RGB_DEFAULT_MAP does not exist, then */
	if (!XGetStandardColormap(dpy,root,scmap,XA_RGB_DEFAULT_MAP)) {
		
		/* use default colormap */
		cmap = DefaultColormapOfScreen(scr);

		/* determine largest number of contiguous free color cells */
		ncells = CellsOfScreen(scr);
		while(ncells && 
			!XAllocColorCells(dpy,cmap,True,NULL,0,pixel,ncells))
			ncells--;
		
		/* determine beginning and ending pixel of contiguous cells */
		for (i=1,bpixel=epixel=pixel1=pixel2=pixel[0]; i<ncells; i++) {
			if (pixel[i]==pixel[i-1]+1)
				pixel2 = pixel[i];
			else
				pixel1 = pixel2 = pixel[i];
			if (pixel2-pixel1>=epixel-bpixel) {
				bpixel = pixel1;
				epixel = pixel2;
			}
		}
		
		/* number of pixels must be at least 8 */
		npixels = epixel-bpixel+1;
		if (npixels<8) {
#ifndef __osf__
			XUngrabServer(dpy);
#endif
			return 0;
		}
		
		/* force number of contiguous cells to be an integer cubed */
		for (i=2,imax=0; i*i*i<=npixels; i++,imax++);
		npixels = (imax+1)*(imax+1)*(imax+1);
		bpixel = epixel-npixels+1;
		
		/* free cells not in contiguous range */
		for (i=0; i<ncells; i++)
			if (pixel[i]<bpixel || pixel[i]>epixel)
				XFreeColors(dpy,cmap,&pixel[i],1,0);

		/* store colors in contiguous range of allocated cells */
		rmult = (imax+1)*(imax+1);
		gmult = imax+1;
		bmult = 1;
		for (i=0; i<npixels; i++) {
			color.pixel = bpixel+i;
			color.red = (unsigned short) (i/rmult);
			color.green = (unsigned short) ((i-color.red*rmult)/gmult);
			color.blue = (unsigned short) (i-color.red*rmult-color.green*gmult);
			color.red *= 65535/imax;
			color.green *= 65535/imax;
			color.blue *= 65535/imax;
			color.flags = DoRed|DoGreen|DoBlue;
			XStoreColor(dpy,cmap,&color);
		}
		
		/* set standard colormap */
		scmap->colormap = cmap;
		scmap->red_max = imax;
		scmap->green_max = imax;
		scmap->blue_max = imax;
		scmap->red_mult = rmult;
		scmap->green_mult = gmult;
		scmap->blue_mult = bmult;
		scmap->base_pixel = bpixel;
		XSetStandardColormap(dpy,root,scmap,XA_RGB_DEFAULT_MAP);
	}

	/* ungrab the server before returning */
#ifndef __osf__
	XUngrabServer(dpy);
#endif
	return 1;
}

unsigned long xGetFirstPixel (Display *dpy)
/*****************************************************************************
return first pixel in range of contiguous pixels in XA_RGB_DEFAULT_MAP
******************************************************************************
Input:
dpy		display
******************************************************************************
Notes:
If it does not already exist, XA_RGB_DEFAULT_MAP will be created.
If XA_RGB_DEFAULT_MAP does not exist and cannot be created, then
this function returns 0.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/29/90
*****************************************************************************/
{
	Screen *scr=XDefaultScreenOfDisplay(dpy);
	Window root=XRootWindowOfScreen(scr);
	XStandardColormap scmap;
	
	/* if XA_RGB_DEFAULT_MAP does not exist, create it */
	if (!XGetStandardColormap(dpy,root,&scmap,XA_RGB_DEFAULT_MAP))
		if (!xCreateRGBDefaultMap(dpy,&scmap))
			return 0;
	
	/* return first pixel in range of contiguous pixels */
	return scmap.base_pixel; 
}

unsigned long xGetLastPixel (Display *dpy)
/*****************************************************************************
return last pixel in range of contiguous pixels in XA_RGB_DEFAULT_MAP
******************************************************************************
Input:
dpy		display
******************************************************************************
Notes:
If it does not already exist, XA_RGB_DEFAULT_MAP will be created.
If XA_RGB_DEFAULT_MAP does not exist and cannot be created, then
this function returns 0.
******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 09/29/90
*****************************************************************************/
{
	Screen *scr=XDefaultScreenOfDisplay(dpy);
	Window root=XRootWindowOfScreen(scr);
	XStandardColormap scmap;
	
	/* if XA_RGB_DEFAULT_MAP does not exist, create it */
	if (!XGetStandardColormap(dpy,root,&scmap,XA_RGB_DEFAULT_MAP))
		if (!xCreateRGBDefaultMap(dpy,&scmap))
			return 0;
	
	/* return last pixel in range of contiguous pixels */
	return scmap.base_pixel+
		scmap.red_max*scmap.red_mult+
		scmap.green_max*scmap.green_mult+
		scmap.blue_max*scmap.blue_mult;
}


/*
 * new  colormaps, the source code is adapted from the routines above
 *
 * Toralf Foerster
 * Baltic Sea Research Institute
 * Rostock Warnemuende
 * Germany
 * 1996
 *
 * email : toralf.foerster@io-warnemuende.de
 *
 * There are many color models, we use 2:
 * the RGB - Model and the HSV - Model (often called HLS).
 * The RGB has the advantage of easy defining the values of the 3
 * basic colors : Red, Green and Blue, the HSV is more closed to 
 * human understanding of colors.
 *
 * The created colormap goes from color[0] over color[1] to color[2]
 * So we have 2 ramps between 3 color points.
 */

Colormap xCreateRGBColormap (Display *dpy, Window win,
			char *str_cmap, int verbose)
{
	Screen *scr=XDefaultScreenOfDisplay(dpy);
	/* Window root=XRootWindowOfScreen(scr);  --unused? */
	Colormap cmap,wcmap;
	XColor color;
	XWindowAttributes wa;
	int i,ncells;
	unsigned long npixels;
	unsigned long bpixel,epixel,pixel[4096];
	unsigned int depth,sr;

/*

# define RGB_BLACK	{0x00, 0x00, 0x00}
# define RGB_WHITE	{0xff, 0xff, 0xff}
# define RGB_GRAY	{0x80, 0x80, 0x80}

# define RGB_ORANGE	{0xff, 0x80, 0x00}

# define RGB_RED	{0xe0, 0x00, 0x50}
# define RGB_BLUE	{0x00, 0x40, 0xc0}
# define RGB_GREEN	{0x06, 0x5b, 0x3f}
# define RGB_BROWN	{0x72, 0x5b, 0x3f}
# define RGB_REDBROWN	{0xa0, 0x40, 0x00}

# define RGB_GRAY2	{0xb0, 0xb0, 0xb0}

# define RGB_LGRAY	{0xf0, 0xf0, 0xf0}
# define RGB_LBLUE	{0x55, 0x9c, 0xe0}
# define RGB_YELLOW	{0xd0, 0xb0, 0x20}

	float c_rgb [][3][3]  =	{
	{ RGB_BLACK,    RGB_GRAY,   RGB_WHITE  },

	{ RGB_RED,      RGB_LGRAY,  RGB_BLUE  },
	{ RGB_RED,      RGB_LGRAY,  RGB_GREEN },
	{ RGB_BROWN,    RGB_LGRAY,  RGB_BLUE  },
	{ RGB_BROWN,    RGB_LGRAY,  RGB_GREEN },
	{ RGB_REDBROWN, RGB_LGRAY,  RGB_BLUE  },
	{ RGB_REDBROWN, RGB_LGRAY,  RGB_GREEN },
	{ RGB_ORANGE,   RGB_LGRAY,  RGB_BLUE  },
	{ RGB_ORANGE,   RGB_LGRAY,  RGB_GREEN },
	{ RGB_BROWN,    RGB_GRAY2,  RGB_GREEN },
	{ RGB_BROWN,    RGB_GRAY2,  RGB_BLUE  },
	{ RGB_BROWN,    RGB_YELLOW, RGB_BLUE  }
	};
*/
	
	static int	c_nr = -1;
	unsigned long	max_cmap, half, ih;

	sr=DefaultScreen(dpy);
	depth=(unsigned int)DefaultDepth(dpy,sr);

	if( verbose == 1 ) {
		warn("\n in colormap, depth=%d\n",depth);
	}

	/* determine window's current colormap */
	XGetWindowAttributes(dpy,win,&wa);
	wcmap = wa.colormap;


	max_cmap = sizeof (c_rgb) / sizeof (float[3][3]);
	/* We got the specific number of the cmap from the string	*/
	if (STREQ (str_cmap, "rgb_up"))
		c_nr++;
	else if (STREQ (str_cmap, "rgb_down"))
		c_nr--;
	else	{
		if (strlen (str_cmap) > 3)	{
			str_cmap[0] = str_cmap[1] = str_cmap[2] = ' ';
			c_nr = atoi (str_cmap);
			if (c_nr < 0 || c_nr >= max_cmap)	{
				warn ("\"cmap=rgb%i\" not installed !", c_nr);
				c_nr = 0;
				warn (" using : \"cmap=rgb%i\"", c_nr);
			}
		}
	}

	/* cycle through the cmaps					*/
	while (c_nr < 0)
		c_nr += max_cmap;
		
	while (c_nr >= max_cmap)
		c_nr -= max_cmap;
		
	if (verbose == 1)
		warn (" using : \"cmap=rgb%i\"", c_nr);


if(depth<=8){
	
	/* determine beginning and ending pixels in contiguous range	*/
	bpixel = xGetFirstPixel(dpy);
	epixel = xGetLastPixel(dpy);
	if (epixel<=bpixel) return None;
	
	
	/* create new colormap and allocate all cells read/write */
	cmap = XCreateColormap(dpy,win,DefaultVisualOfScreen(scr),AllocNone);
	ncells = CellsOfScreen(scr);
	XAllocColorCells(dpy,cmap,True,NULL,0,pixel,ncells);
	
	/* copy color cells from window's colormap to new colormap */
	for (i=0; i<ncells; ++i) {
		if (i<bpixel || i>epixel) {
			color.pixel = i;
			XQueryColor(dpy,wcmap,&color);
			XFreeColors(dpy,cmap,&pixel[i],1,0);
			XAllocColor(dpy,cmap,&color);
		}
	}
	
	/* build scale in contiguous cells in new colormap */
	npixels = epixel-bpixel+1;
	half = npixels / 2;
	

	/* Build the 1st ramp						*/
	for (ih = 0; ih < half; ++ih) {
		color.pixel = bpixel + ih;
		color.red   = c_rgb[c_nr][0][0] +
			(c_rgb[c_nr][1][0] - c_rgb[c_nr][0][0]) * ((float) ih)/((float) half);
		color.green = c_rgb[c_nr][0][1] +
			(c_rgb[c_nr][1][1] - c_rgb[c_nr][0][1]) * ((float) ih)/((float) half);
		color.blue  = c_rgb[c_nr][0][2] +
			(c_rgb[c_nr][1][2] - c_rgb[c_nr][0][2]) * ((float) ih)/((float) half);
		
		color.red   *= 257.0;
		color.green *= 257.0;
		color.blue  *= 257.0;
		
		color.flags = DoRed|DoGreen|DoBlue;
		XStoreColor(dpy,cmap,&color);
	}

	/* Build the 2nd ramp						*/
	for (ih=half; ih<npixels; ++ih) {
		color.pixel = bpixel+ih;
		color.red   = c_rgb[c_nr][1][0] +
			(c_rgb[c_nr][2][0] - c_rgb[c_nr][1][0]) * ((float) (ih-half))/((float) half);
		color.green = c_rgb[c_nr][1][1] +
			(c_rgb[c_nr][2][1] - c_rgb[c_nr][1][1]) * ((float) (ih-half))/((float) half);
		color.blue  = c_rgb[c_nr][1][2] +
			(c_rgb[c_nr][2][2] - c_rgb[c_nr][1][2]) * ((float) (ih-half))/((float) half);
		
		color.red   *= 257.0;
		color.green *= 257.0;
		color.blue  *= 257.0;
		
		color.flags = DoRed|DoGreen|DoBlue;
		XStoreColor(dpy,cmap,&color);
	}

	/* return colormap */
	return cmap;
}
else{


        /* Build the 1st ramp                                           */
        for (ih = 0; ih < 128; ++ih) {
                color.red   = c_rgb[c_nr][0][0] +
                        (c_rgb[c_nr][1][0] - c_rgb[c_nr][0][0]) * ((float) ih)/((float) 128);
                color.green = c_rgb[c_nr][0][1] +
                        (c_rgb[c_nr][1][1] - c_rgb[c_nr][0][1]) * ((float) ih)/((float) 128);
                color.blue  = c_rgb[c_nr][0][2] +
                        (c_rgb[c_nr][1][2] - c_rgb[c_nr][0][2]) * ((float) ih)/((float) 128);
                
                color.red   *= 257.0;
                color.green *= 257.0;
                color.blue  *= 257.0;
        
                XAllocColor(dpy,wcmap,&color);
		truecolor_pixel[ih]=(unsigned long)color.pixel;
        }


        /* Build the 2nd ramp                                           */
        for (ih=128; ih<256; ++ih) {
                color.red   = c_rgb[c_nr][1][0] +
                        (c_rgb[c_nr][2][0] - c_rgb[c_nr][1][0]) * ((float) (ih-128))/((float) 128);
                color.green = c_rgb[c_nr][1][1] +
                        (c_rgb[c_nr][2][1] - c_rgb[c_nr][1][1]) * ((float) (ih-128))/((float) 128);
                color.blue  = c_rgb[c_nr][1][2] +
                        (c_rgb[c_nr][2][2] - c_rgb[c_nr][1][2]) * ((float) (ih-128))/((float) 128);
        
                color.red   *= 257.0;
                color.green *= 257.0;
                color.blue  *= 257.0;
 
                XAllocColor(dpy,wcmap,&color);
		truecolor_pixel[ih]=(unsigned long)color.pixel;
        }

	return wcmap;

}


}


#include "cwpcmaps.h"

Colormap xCreateHSVColormap (Display *dpy, Window win,
			char *str_cmap, int verbose)
{
	Screen *scr=XDefaultScreenOfDisplay(dpy);
	/* Window root=XRootWindowOfScreen(scr); --unused? */
	Colormap cmap,wcmap;
	XColor color;
	XWindowAttributes wa;
	int i,ncells;
	unsigned long npixels;
	unsigned long bpixel,epixel,pixel[4096];

/*
# define HSV_BLACK	{  0.0, 0.00, 0.00}
# define HSV_GRAY	{  0.0, 0.00, 0.50}
# define HSV_WHITE	{  0.0, 0.00, 1.00}

# define HSV_HUE1	{240.0, 1.00, 0.50}
# define HSV_HUE2	{120.0, 1.00, 0.50}
# define HSV_HUE3	{  0.0, 1.00, 0.50}

# define HSV_DRED	{  0.0, 1.00, 0.50}
# define HSV_BROWN	{ 30.0, 1.00, 0.30}
# define HSV_GREEN	{140.0, 1.00, 0.50}
# define HSV_BLUE	{240.0, 1.00, 0.70}

# define HSV_YELLOW	{ 70.0, 1.00, 0.50}

	float c_hsv [][3][3]  =	{
	{ HSV_WHITE,  HSV_GRAY,   HSV_BLACK },
	{ HSV_HUE1,   HSV_HUE2,   HSV_HUE3  },
	{ HSV_HUE3,   HSV_HUE2,   HSV_HUE1  },
	{ HSV_BROWN,  HSV_GREEN,  HSV_BLUE  },
	{ HSV_DRED,  HSV_WHITE,  HSV_BLUE  },
	{ HSV_BLUE,  HSV_WHITE,  HSV_DRED  },
	{ HSV_WHITE,  HSV_DRED,  HSV_BLUE  },
	{ HSV_WHITE,  HSV_GREEN,  HSV_BLUE  },
	{ HSV_BLUE,  HSV_DRED,  HSV_WHITE  },
	{ HSV_BLUE,  HSV_GREEN,  HSV_WHITE  },
	{ HSV_BLUE,  HSV_WHITE,  HSV_GREEN  },
	{ HSV_YELLOW,  HSV_DRED,  HSV_BROWN  },
	{ HSV_BROWN,  HSV_DRED,  HSV_YELLOW  },
	{ HSV_DRED,  HSV_YELLOW,  HSV_BROWN  }
	};
*/

	static int	c_nr = -1;
	unsigned long	max_cmap, half, ih;
	float r,g,b, h,s,v;	/* Red,Green,Blue, Hue,Sat,Value	*/
	unsigned int depth, sr;

        sr=DefaultScreen(dpy);
        depth=(unsigned int)DefaultDepth(dpy,sr);
	

	/* determine window's current colormap */
	XGetWindowAttributes(dpy,win,&wa);
	wcmap = wa.colormap;

	max_cmap = sizeof (c_hsv) / sizeof (float[3][3]);
	
	/* We got the specific number of the cmap from the string	*/
	if (STREQ (str_cmap, "hsv_up"))
		c_nr++;
	else if (STREQ (str_cmap, "hsv_down"))
		c_nr--;
	else	{
		if (strlen (str_cmap) > 3)	{
			str_cmap[0] = str_cmap[1] = str_cmap[2] = ' ';
			c_nr = atoi (str_cmap);
			if (c_nr < 0 || c_nr >= max_cmap)	{
				warn ("\"cmap=hsv%i\" not installed !", c_nr);
				c_nr = 0;
				warn (" using : \"cmap=hsv%i\"", c_nr);
			}
		}
	}

	/* cycle through the cmaps					*/
	while (c_nr < 0)
		c_nr += max_cmap;
		
	while (c_nr >= max_cmap)
		c_nr -= max_cmap;
		
	if (verbose == 1)
		warn (" using : \"cmap=hsv%i\"", c_nr);


if(depth<=8){
	
	/* determine beginning and ending pixels in contiguous range	*/
	bpixel = xGetFirstPixel(dpy);
	epixel = xGetLastPixel(dpy);
	if (epixel<=bpixel) return None;
	
	
	/* create new colormap and allocate all cells read/write */
	cmap = XCreateColormap(dpy,win,DefaultVisualOfScreen(scr),AllocNone);
	ncells = CellsOfScreen(scr);
	XAllocColorCells(dpy,cmap,True,NULL,0,pixel,ncells);
	
	/* copy color cells from window's colormap to new colormap */
	for (i=0; i<ncells; ++i) {
		if (i<bpixel || i>epixel) {
			color.pixel = i;
			XQueryColor(dpy,wcmap,&color);
			XFreeColors(dpy,cmap,&pixel[i],1,0);
			XAllocColor(dpy,cmap,&color);
		}
	}
	
	/* build scale in contiguous cells in new colormap */
	npixels = epixel-bpixel+1;
	half = npixels / 2;
	
	
	/* Build the 1st ramp						*/
	for (ih = 0; ih < half; ++ih) {
		color.pixel = bpixel + ih;
		
		h = c_hsv[c_nr][0][0] +
			(c_hsv[c_nr][1][0] - c_hsv[c_nr][0][0]) * ((float) ih) /((float) half);
		s = c_hsv[c_nr][0][1] +
			(c_hsv[c_nr][1][1] - c_hsv[c_nr][0][1]) * ((float) ih) / ((float) half);
		v = c_hsv[c_nr][0][2] +
			(c_hsv[c_nr][1][2] - c_hsv[c_nr][0][2]) * ((float) ih) / ((float) half);
		
		hsv2rgb (h, s, v, &r, &g, &b);
		color.red   = 65535.0 * r;
		color.green = 65535.0 * g;
		color.blue  = 65535.0 * b;
		
		color.flags = DoRed|DoGreen|DoBlue;
		XStoreColor(dpy,cmap,&color);
	}

	/* Build the 2nd ramp						*/
	for (ih = half; ih < npixels; ++ih) {
		color.pixel = bpixel + ih;
		
		h = c_hsv[c_nr][1][0] +
			(c_hsv[c_nr][2][0] - c_hsv[c_nr][1][0]) * ((float) (ih-half))/((float) half);
		s = c_hsv[c_nr][1][1] +
			(c_hsv[c_nr][2][1] - c_hsv[c_nr][1][1]) * ((float) (ih-half))/((float) half);
		v = c_hsv[c_nr][1][2] +
			(c_hsv[c_nr][2][2] - c_hsv[c_nr][1][2]) * ((float) (ih-half))/((float) half);
		
		hsv2rgb (h, s, v, &r, &g, &b);
		color.red   = 65535.0 * r;
		color.green = 65535.0 * g;
		color.blue  = 65535.0 * b;
		
		color.flags = DoRed|DoGreen|DoBlue;
		XStoreColor(dpy,cmap,&color);
	}
		
	/* return colormap */
	return cmap;

}else{

        /* Build the 1st ramp                                           */
        for (ih = 0; ih < 128; ++ih) {
                h = c_hsv[c_nr][0][0] + 
                        (c_hsv[c_nr][1][0] - c_hsv[c_nr][0][0]) * ((float) ih) /((float) 128);
                s = c_hsv[c_nr][0][1] +
                        (c_hsv[c_nr][1][1] - c_hsv[c_nr][0][1]) * ((float) ih) / ((float) 128);
                v = c_hsv[c_nr][0][2] +
                        (c_hsv[c_nr][1][2] - c_hsv[c_nr][0][2]) * ((float) ih) / ((float) 128);
        
                hsv2rgb (h, s, v, &r, &g, &b);
                color.red   = 65535.0 * r;
                color.green = 65535.0 * g;
                color.blue  = 65535.0 * b;
                XAllocColor(dpy,wcmap,&color);
                truecolor_pixel[ih]=(unsigned long)color.pixel;
        }

        /* Build the 2nd ramp                                           */
        for (ih = 128; ih < 256; ++ih) {
                h = c_hsv[c_nr][1][0] +
                        (c_hsv[c_nr][2][0] - c_hsv[c_nr][1][0]) * ((float) (ih-128))/((float) 128);
                s = c_hsv[c_nr][1][1] +
                        (c_hsv[c_nr][2][1] - c_hsv[c_nr][1][1]) * ((float) (ih-128))/((float) 128);
                v = c_hsv[c_nr][1][2] +
                        (c_hsv[c_nr][2][2] - c_hsv[c_nr][1][2]) * ((float) (ih-128))/((float) 128);
                
                hsv2rgb (h, s, v, &r, &g, &b);
                color.red   = 65535.0 * r;
                color.green = 65535.0 * g;
                color.blue  = 65535.0 * b;

                XAllocColor(dpy,wcmap,&color);
                truecolor_pixel[ih]=(unsigned long)color.pixel;
        }

        /* return colormap */
        return wcmap;


}


}


/* internal functions to convert HSV to RGB */
static float rgbvalue (float n1, float n2, float hue)
{
	while (hue > 360.0)
		hue -= 360.0;
	while (hue < 0.0)
		hue += 360.0;
	
	if (hue < 60.0)
		return n1 + (n2 - n1) * hue / 60.0;
	else if (hue<180.0)
		return n2;
	else if (hue < 240.0)
		return n1 + (n2 - n1) * (240.0 - hue) / 60.0;
	else
		return n1;
}

/*
 * variable	range
 * --------	------------
 *	h	0.0 .. 360.0
 *	s	0.0 .. 1.0
 *	v	0.0 .. 1.0
 */
static void hsv2rgb (float h, float s, float v, float *r, float *g, float *b)
{
	float m1,m2;
	/* float rgbvalue (float,float,float);*/
	
	if (v <= 0.5)
		m2 = v * (1.0 + s);
	else
		m2 = v + s - v * s;
	m1 = 2 * v - m2;
	if (s == 0.0) {
                *r = *g = *b = v;
	} else {
		*r = rgbvalue(m1, m2, h + 120.0);
		*g = rgbvalue(m1, m2, h);
		*b = rgbvalue(m1, m2, h - 120.0);

		if (*r > 1.0)
			*r = 1.0;
		if (*g > 1.0)
			*g = 1.0;
		if (*b > 1.0)
			*b = 1.0;
	}
}
	
/* test program - compile with "cc -DTEST colormap.c $INCS $LIBS ..." */
#ifdef TEST
#include <stdio.h>
#define X 100
#define Y 100
#define WIDTH 256
#define HEIGHT 64
main()
{
	Display *dpy;
	Window root,win;
	Colormap cmap;
	XStandardColormap scmap;
	XColor color,junk;
	XImage *image;
	XEvent event;
	GC gc;
	int scr,i;
	unsigned long black,white,pmin,pmax;
	char *data;
	
	/* connect to X server */
	dpy = XOpenDisplay(NULL);
	if ((dpy=XOpenDisplay(NULL))==NULL) {
		fprintf(stderr,"Cannot open display!\n");
		exit(-1);
	}
	scr = DefaultScreen(dpy);
	root = RootWindow(dpy,scr);
	black = BlackPixel(dpy,scr);
	white = WhitePixel(dpy,scr);
	
	/* create and map window */
	win = XCreateSimpleWindow(dpy,root,X,Y,WIDTH,HEIGHT,4,black,white);
	cmap = xCreateRGBColormap(dpy,win, "rgb0", 1);
	XSetWindowColormap(dpy,win,cmap);
	XMapWindow(dpy,win);
	
	/* determine range of contiguous pixels from standard colormap */
	if (!xCreateRGBDefaultMap(dpy,&scmap)) {
		fprintf(stderr,"Cannot create standard colormap!\n");
		exit(-1);
	}
	pmin = xGetFirstPixel(dpy);
	pmax = xGetLastPixel(dpy);
	
	/* create image */
	data = (char*)malloc(WIDTH*HEIGHT);
	for (i=0; i<WIDTH*HEIGHT; ++i)
		data[i] = pmin+(pmax-pmin)*(i%WIDTH)/WIDTH;
	image = XCreateImage(dpy,DefaultVisual(dpy,scr),
		DefaultDepth(dpy,scr),ZPixmap,
		0,data,WIDTH,HEIGHT,BitmapPad(dpy),WIDTH);
	gc = XCreateGC(dpy,win,0,NULL);
	XAllocNamedColor(dpy,cmap,"red",&color,&junk);
	XSetForeground(dpy,gc,color.pixel);
	
	/* set event mask */
	XSelectInput(dpy,win,ExposureMask);

	/* loop forever */
	XPutImage(dpy,win,gc,image,0,0,0,0,WIDTH,HEIGHT);
	while(True) {
		XNextEvent(dpy,&event);
		while (XCheckTypedEvent(dpy,Expose,&event));
		XPutImage(dpy,win,gc,image,0,0,0,0,WIDTH,HEIGHT);
		XDrawLine(dpy,win,gc,0,0,WIDTH,HEIGHT);
	}

	/* close display */
	XCloseDisplay(dpy);
}
#endif /* TEST */
