/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* xwindow.h for xpicker */

#ifndef X_WINDOW_H
#define X_WINDOW_H
#ifndef EGSTERN
#define EGSTERN extern
#endif/*EGSTERN*/


#define CELLS 2 /** this is how many color cells we want to grab **/
/***************************************************/
/** These definitions are for the garnish buttons **/
/***************************************************/
#define SHADOW 30000
#define RADIO_SHADOW 15000
#define SHRINK 0
#define EXPAND 1
#define RIGHT 3
#define LEFT 2
#define UP 1 
#define DOWN 0
#define RELATIVE 1
#define FIXED 2
#define TOP 3
#define BOTTOM 14
/***************************************************/
#define COLOR_FLAGS DoRed | DoGreen | DoBlue
#define XEVENT_MASK ButtonPressMask | KeyPressMask | ExposureMask | PointerMotionMask | ButtonReleaseMask | KeyReleaseMask | EnterWindowMask | LeaveWindowMask | FocusChangeMask

EGSTERN Display     		*display;
EGSTERN Visual			*visual;
EGSTERN Window    	        window,win2;
EGSTERN Window                  inwin;
EGSTERN GC 	      		gc;
EGSTERN Colormap    		colormap;
EGSTERN XColor      		generic, foreground_color, background_color;
EGSTERN XEvent      		event;
EGSTERN Pixmap                  draw;
EGSTERN Pixmap                  highlight;
EGSTERN Pixmap                  escherknot;
EGSTERN XImage			*image;
EGSTERN KeySym      		key;
EGSTERN XSizeHints  		hint;
EGSTERN XWMHints                wmhint;
EGSTERN Cursor                  standard_cursor, busy_cursor;
EGSTERN Cursor                  grid_cursor, mutate_cursor;
EGSTERN Status      		result;
EGSTERN Font 			font;
EGSTERN XFontStruct             *font_struct;
EGSTERN Time                    time_stamp;
EGSTERN unsigned long		width,height;
EGSTERN int			xhot,yhot,status,x;
EGSTERN unsigned long 		plane_mask,depth;
EGSTERN unsigned long 		pixel;
EGSTERN unsigned long 		foreground, background;
EGSTERN int			bitmap_pad;
EGSTERN int         		screen;
EGSTERN int                     cells;
EGSTERN int                     current_pixel;
EGSTERN char                    text[11]; /** the text read from keyboard **/

#endif /* end of X_WINDOW_H */

