/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
WINDOW - Function to create a window in X-windows graphics

xNewWindow	Create a new window and return the window ID

******************************************************************************
Function Prototype:
Window xNewWindow (Display *dpy, int x, int y, int width, int height,
	int border, int background, char *name);

******************************************************************************
Input:
dpy		display pointer
x		x in pixels of upper left corner
y		y in pixels of upper left corner
width		width in pixels
height		height in pixels
border		border pixel
background	background pixel
name		name of window (also used for icon)

******************************************************************************
Notes:
The parent window is the root window.
The border_width is 4 pixels.

******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/06/90
*****************************************************************************/
/**************** end self doc ********************************/

#include "xplot.h"


Window
xNewWindow (Display *dpy, int x, int y, int width, int height,
	int border, int background, char *name)
/*****************************************************************************
Create a new window and return the window ID
******************************************************************************
Input:
dpy		display pointer
x		x in pixels of upper left corner
y		y in pixels of upper left corner
width		width in pixels
height		height in pixels
border		border pixel
background	background pixel
name		name of window (also used for icon)
******************************************************************************
Notes:
The parent window is the root window.
The border_width is 4 pixels.
******************************************************************************
Author:		Dave Hale, Colorado School of Mines, 01/06/90
*****************************************************************************/
{
	Window root,win;
	XSizeHints size_hints;
	XWMHints wm_hints;
	int scr,border_w=4;

	/* get screen and root window */
	scr = DefaultScreen(dpy);
	root = RootWindow(dpy,scr);

	/* create window */
	win = XCreateSimpleWindow(dpy,root,x,y,width,height,
		border_w,border,background);

	/* set window properties for window manager */
	size_hints.flags = USPosition|USSize;
	size_hints.x = x;
	size_hints.y = y;
	size_hints.width = width;
	size_hints.height = height;
	XSetStandardProperties(dpy,win,name,name,None,0,0,&size_hints);
	wm_hints.flags = InputHint;
	wm_hints.input = True;
	XSetWMHints(dpy,win,&wm_hints);

	/* return window ID */
	return win;
}
