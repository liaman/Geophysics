/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* GARNISH: $Revision: 1.8 $ ; $Date: 2011/11/21 17:02:44 $	*/

#include "par.h"
#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include "xwindow.h"
#include "garnish.h"
#include "picking.h"

/*****************************************************************************/
/** This will draw shaded buttons of either "fixed" length or the length of **/
/** the string to be written in them.  The choice is made through style     **/
/*****************************************************************************/
void
DrawButton(window,style,state,x,y,textheight, textwidth, fixed, text, color)
     Window			window;
     int			style,state,x,y,textheight,textwidth,fixed;
     char			text[];
     XColor			color;
{
  int 				w;

  XFlush(display);
  XSetFillStyle(display, gc, FillSolid);
  XSetLineAttributes(display, gc, 3,LineSolid,CapRound,JoinMiter);
  XClearArea(display, window, x,y,textwidth*((int) strlen(text))+6,textheight + 6,False);
  generic.flags = DoRed | DoBlue | DoGreen;
  generic.red = color.red; generic.green = color.green; generic.blue = color.blue;
  if (generic.red > 0) generic.red -= SHADOW;
  if (generic.blue > 0) generic.blue -= SHADOW;
  if (generic.green > 0) generic.green -= SHADOW;
  (void) XAllocColor(display, colormap, &generic);
  if (style == RELATIVE)
    w = textwidth*((int) strlen(text));
  else
    w = fixed;
  if (state == UP) {
    XSetForeground(display, gc, color.pixel);
    XFillRectangle(display, window, gc, x+3, y+3,(w), textheight+12);
    XSetForeground(display, gc, generic.pixel);
    XDrawRectangle(display, window, gc, x,y,(w+6), textheight+14);
    XSetForeground(display, gc, BlackPixel(display, DefaultScreen(display)));
    XSetBackground(display, gc, color.pixel);
    XDrawImageString(display, window, gc, x+3,y+textheight+3,text,(int) strlen(text));
  }
  else {
    XSetForeground(display, gc, generic.pixel);
    XFillRectangle(display, window, gc, x+3, y+3,(w), textheight+12);
    XSetForeground(display, gc, color.pixel);
    XDrawRectangle(display, window, gc, x,y,(w+6), textheight+14);
    XSetForeground(display, gc, BlackPixel(display, DefaultScreen(display)));
    XSetBackground(display, gc, generic.pixel);
    XDrawImageString(display, window, gc, x+4,y+textheight+3,text,(int) strlen(text));
  }
  XSync(display,0);
  XFlush(display);
}


void DrawRadio(int state, int x, int y, XColor color, int size) 
{
  generic.flags = DoRed | DoBlue | DoGreen;
  generic.red =color.red;generic.green=color.green;generic.blue = color.blue;
  if (generic.red > 0) generic.red -= RADIO_SHADOW;
  if (generic.blue > 0) generic.blue -= RADIO_SHADOW;
  if (generic.green > 0) generic.green -= RADIO_SHADOW;
  XAllocColor(display, colormap, &generic);
  if (state == UP) {
    XSetForeground(display, gc, color.pixel);
    XFillRectangle(display, window, gc, x, y, size,size);
    XSetForeground(display, gc, generic.pixel);
    XFillRectangle(display, window, gc, x+2, y+2,size-4,size-4);
  }
  else {
    XSetForeground(display, gc, generic.pixel);
    XFillRectangle(display, window, gc, x, y, size,size);
    XSetForeground(display, gc, color.pixel);
    XFillRectangle(display, window, gc, x+2, y+2,size-4,size-4);
  }
}

/*********************************************************************/
/** This procedure draws a shadowed border box of 'thick' thickness **/
/** Placed at the coordinates indicated by 'xcorner', 'ycorner'.The **/
/** shadowed box will be 'length' long and 'height' high.           **/
/*********************************************************************/
int DrawShadowedBox(Window win,XColor color, int thick, int xcorner,
		     int ycorner, int length, int height)
{
  int count, increment,max=0;
  static int first = True;
  static XColor colors[MAX_BOX_THICK];
  
  if (thick > MAX_BOX_THICK) {
    fprintf(stderr, "ERROR : shadow box thickness specified too large.\n");
    return -1;
  }
  if (first) {
    first = False;
    count = 0;
    increment = 30000/(MAX_BOX_THICK/2);
    while (count < MAX_BOX_THICK) {
      if (count < MAX_BOX_THICK/2) {
	if (color.red > 1) color.red += increment;
	if (color.green > 1) color.green += increment;
	if (color.blue > 1) color.blue += increment;
      }
      else {
	if (color.red > 1) color.red -= increment; 
	if (color.green > 1) color.green -= increment;
	if (color.blue > 1) color.blue -= increment;
      }
      colors[count] = color;
      XAllocColor(display, colormap, &colors[count]);
      count++;
    }
  }
  increment = 0;
  count = 0; 
  while (count < thick) {
    XSetForeground(display, gc, colors[increment].pixel);
    
    /** Draw left border **/
    XDrawLine(display, win, gc, xcorner+count, ycorner, 
	      xcorner+count, ycorner+height+thick*2);
    /** Draw top border **/
    XDrawLine(display, win, gc, xcorner, ycorner+count, 
	      xcorner+length+thick*2, ycorner+count);
    /** Draw right border **/
    XDrawLine(display, win, gc, xcorner+length+thick+count, ycorner,
	      xcorner+thick+length+count, ycorner+height+thick*2);
    /** Draw bottom border **/
    XDrawLine(display, win, gc, xcorner, ycorner+height+count+thick,
	      xcorner+length+thick*2, ycorner+height+count+thick);
    if (count <= thick/2) {
      increment ++;
      max = increment;
    }
    else {
      if (increment == max)
	increment = MAX_BOX_THICK/2+(MAX_BOX_THICK/2-increment);
      else
	increment ++;
    }
    count++;
  }
  return EXIT_SUCCESS;
}
 
/*************************************************************************/
/** This procedure sets up a radio button set.  The programmer sends in **/
/** the number of buttons desired, the space between buttons, and the   **/
/** size of the buttons.  Note that currently you can only get          **/
/** vertically distributed radio buttons.                               **/
/*************************************************************************/
void DisplayRadioButtonSet(int x, int y, int number, int space, int size,
		      int state[], XColor color)
{
  int loop;

  XSetForeground(display, gc, color.pixel);
  XSetLineAttributes(display, gc, 1,LineSolid,CapRound,JoinMiter);
  for(loop = 0; loop < number; loop++) {
    XSetForeground(display, gc, background);
    XFillArc(display, window, gc, x, y, size,size,0,64*360);
    XSetForeground(display, gc, color.pixel);
    if (state[loop] == UP) 
      XDrawArc(display, window, gc, x, y, size,size,0,64*360);
    else
      XFillArc(display, window, gc, x, y, size,size,0,64*360);
    y += (size+space);
  }
}
    


/********************************************************************/
/** This procedure checks the radio button set specified with the  **/
/** coordinates given.  It returns 0 on failure, or 1-n depending  **/
/** on which radio button the coordinates corrospond to.           **/
/********************************************************************/
int CheckRadioButtonSet(x,y,number,space,size,xpos,ypos)
     int x,y,number,space,size,xpos,ypos;
{
  int ret,diff,extra;

  if (xpos < x || xpos > x+size) return 0; /** Not within the x-range **/
  if (ypos < y) return 0;  /** Too high **/
  if (ypos > y+(number*size)+(space*(number-1))) return 0; /** Too low **/

  diff = ypos-y;
  ret = diff/(size+space)+1;
  extra = diff%size;
  if (extra > size) return 0; /** coordinate was in the space between **/
  return ret;
}  

void NewButton(Window win, int x, int y, int width, int height, int state,
	  XColor color, XColor textcolor, char text[], int c_width,
	  int c_height, int mode)
{
  static XColor dark;
  static XColor medium;
  static XColor light;
  int tx,ty;
  int bw,bh;
  int loop;
  int index = 0;
  static int count = 0;
  static XColor colors[5][4];

  XSetFillStyle(display, gc, FillSolid);
  for (loop = 0; loop < count; loop++) {
    if (colors[loop][0].pixel == color.pixel) {
      index = loop;
      loop = count+1;
    }
  }
  if (loop == count) {
    count++;
    colors[loop][0] = color;
    dark.red = color.red; dark.green = color.green; dark.blue = color.blue; 
    if (dark.red > 1) dark.red -= 20000;
    if (dark.green > 1) dark.green -= 20000;
    if (dark.blue > 1) dark.blue -= 20000;
    dark.flags = COLOR_FLAGS;
    colors[loop][1] = dark;
    XAllocColor(display, colormap, &colors[loop][1]);
    dark = colors[loop][1];
    
    medium = dark;
    if (medium.red > 1) medium.red += 16000;
    if (medium.green > 1) medium.green += 16000;
    if (medium.blue > 1) medium.blue += 16000;
    colors[loop][2] = medium;
    XAllocColor(display, colormap, &colors[loop][2]);
    medium = colors[loop][2];
    

    light.red = color.red; light.green = color.green; light.blue = color.blue;
    if (light.red > 1) light.red += 20000;
    if (light.green > 1) light.green += 20000;
    if (light.blue > 1) light.blue += 20000;
    light.flags = COLOR_FLAGS;
    colors[loop][3] = light;
    XAllocColor(display, colormap, &colors[loop][3]);
    light = colors[loop][3];
  }
  else {
    dark = colors[index][1];
    medium = colors[index][2];
    light = colors[index][3];
  }
  bw = 3;
  bh = 3;
  
  XSetLineAttributes(display, gc, 1,LineSolid,CapRound,JoinMiter);
  /** Draw the main text rectangle **/
  if (state == UP) 
    XSetForeground(display, gc, color.pixel);
  else
    XSetForeground(display, gc, medium.pixel);
  XFillRectangle(display, win, gc, x,y,width,height);
  
  /** Place the text in the rectangle **/
  if (((int) strlen(text)) > 0) {
    if (mode == FIXED)
      tx = 0;
    else
      tx = width/2-(((int) strlen(text))*c_width)/2;
    ty = height - (height - c_height)/2;
    XSetForeground(display, gc, textcolor.pixel);
    XDrawString(display, win, gc, x+tx,y+ty,text,(int) strlen(text));
  }

  /** Draw the upper highlight border **/
  if (state == UP)
    XSetForeground(display, gc, light.pixel);
  else
    XSetForeground(display, gc, dark.pixel);
  XFillRectangle(display, win, gc, x-bw, y-bh+1,bw,height+bh);
  for(loop = 0; loop < bh; loop++) 
    XDrawLine(display, win, gc, x, y-loop,x+width+loop,y-loop);
  if (state == UP)
    XSetForeground(display, gc, dark.pixel);
  else
    XSetForeground(display, gc, light.pixel);
  XFillRectangle(display, win, gc, x-bw, y+height,width+bw*2,bh);
  for(loop = 0; loop < bw; loop++)
    XDrawLine(display, win, gc, x+width+loop, y-loop,x+width+loop,
	      y+height);
}

void
DrawScrollBar(Window win,int x,int y,int width,int height,int pos,XColor color)
{
  XColor accentuate;
  XColor foreground;

  if ((((unsigned int) (color.red)) + ((unsigned int) color.blue) +
	((unsigned int) color.green)) < ((unsigned int) 90000))
    foreground.red = foreground.blue = foreground.green = 0;
  else
    foreground.red = foreground.blue = foreground.green = 65000;
  XAllocColor(display, colormap, &foreground);
  accentuate.red = color.red;
  accentuate.green = color.green;
  accentuate.blue = color.blue;
  accentuate.flags = COLOR_FLAGS;
  if (accentuate.red > 1) accentuate.red -= 15000;
  if (accentuate.green > 1) accentuate.green -= 15000;
  if (accentuate.blue > 1) accentuate.blue -= 15000;
  XAllocColor(display, colormap, &accentuate);
  XSetForeground(display, gc, color.pixel);
  XFillRectangle(display, win, gc, x,y,width,height);
  /** Draw the Arrows Boxes **/
  XSetForeground(display, gc, accentuate.pixel);
  XFillRectangle(display, win, gc, x,y,20,width);
  XFillRectangle(display, win, gc, x,y+height-20,20,width);
  XSetForeground(display, gc, foreground.pixel);
  XSetLineAttributes(display, gc, 1,LineSolid,CapRound,JoinMiter);
  /** Draw up arrow **/
  XDrawLine(display, win, gc, x+width/2,y+18,x+width/2,y+2);
  XDrawLine(display, win, gc, x+width/2-4,y+6,x+width/2,y+2);
  XDrawLine(display, win, gc, x+width/2+4,y+6,x+width/2,y+2);
  /** Draw down arrow **/
  XDrawLine(display, win, gc, x+width/2,y+height-2,x+width/2,y+height-18);
  XDrawLine(display, win, gc, x+width/2-4,y+height-6,x+width/2,y+height-2);
  XDrawLine(display, win, gc, x+width/2+4,y+height-6,x+width/2,y+height-2);
  /** Draw the position marker **/
  XSetForeground(display, gc, foreground.pixel);
  XFillRectangle(display, win, gc, x+2,y+20+pos*(height-40)/40,
		 width-4,(height-40)/40);
}

void RefreshTextLine(TextSet *set, int which)
{
  int x, y,l;
  char stemp[80];
  int state;
  XColor fore,back;

  if (which == set->currently_pointed) 
    state = DOWN;
  else
    state = UP;
  fore.pixel = set->foreground;
  back.pixel = set->background;
  fore.flags = back.flags = COLOR_FLAGS;
  XQueryColor(display, colormap, &fore);
  XQueryColor(display, colormap, &back);
  x = set->startx; y = set->starty;
  y += which*(set->space + set->char_height);
  XSetForeground(display,  gc, set->background);
/*  XFillRectangle(display, set->win, gc, set->startx-3, y-set->char_height,
		 set->max*set->char_width+6,set->char_height+4); */
  XSetFont(display, gc, set->font);
  XSetForeground(display, gc, set->foreground);
  XSetBackground(display, gc, set->background);
/*
  XDrawLine(display, set->win, gc, x, y+1, x+set->max*set->char_width, y+1);
*/
  l = (int) strlen(&set->text[which][set->view_start[which]]);
  memset(stemp, 0, sizeof(stemp));
  if (l < set->max)
    strncpy(stemp, &set->text[which][set->view_start[which]],
	    l);
  else
    strncpy(stemp, &set->text[which][set->view_start[which]],
	    set->max);
  NewButton(set->win, x, y-set->char_height/*-set->space/3*/, 
	    set->max*set->char_width,
	    set->char_height, state, back, fore, 
	    stemp, set->char_width, set->char_height,FIXED);
  XSetForeground(display, gc, set->foreground);
  if (which == set->currently_pointed) {
    int loop;

    x += (set->char_pos-set->view_start[set->currently_pointed])*
      set->char_width;
    y -= 2;
    for (loop = 0; loop < 4; loop++) {
      XDrawLine(display, set->win, gc, x-loop, y+loop, x,y);
      XDrawLine(display, set->win, gc, x-loop, y+loop, x+loop,y+loop);
      XDrawLine(display, set->win, gc, x+loop, y+loop, x, y);
    }
  }
}

/*****************************************************************************/
/** This procedure refreshes the entire test set.                           **/
/*****************************************************************************/
void RefreshTextSet(TextSet *set)
{
  int loop;

  for (loop = 0; loop < set->number; loop++) 
    RefreshTextLine(set, loop);
}

void CheckTextSet(TextSet *set, int x, int y)
{
  int dist; /* distance of pointer from start of line */

  dist = x-set->startx;
  if (x < set->startx || x > set->startx+set->max*set->char_width)
    return;
  if (y < set->starty - set->char_height ||
      y > set->starty + (set->char_height*set->number-1)+
      (set->space*(set->number-1)))
    return;
  if (y > set->starty-set->char_height && y < set->starty)  {
    int i;

    i = set->currently_pointed;
    set->currently_pointed = 0;
    
    set->char_pos = set->view_start[set->currently_pointed] +
      (dist/set->char_width);
    if (set->char_pos > (((int) strlen(set->text[set->currently_pointed]))-1))
      set->char_pos = (int) strlen(set->text[set->currently_pointed]);
    if (i >= 0)
      RefreshTextLine(set, i);
    if (set->number != 1)
      RefreshTextLine(set, 0);
  }
  else {
    int i;
    
    if (set->number < 2) 
      return;
    y = y - set->starty;
    y = y / (set->char_height+set->space);
    i = set->currently_pointed;
    if (y >= set->number-1)  return;
    set->currently_pointed = y+1;
    set->char_pos = set->view_start[set->currently_pointed]
      + (dist/set->char_width);
    if (set->char_pos > (((int) strlen(set->text[set->currently_pointed]))-1))
      set->char_pos = (int) strlen(set->text[set->currently_pointed]);
    if (i >= 0 && i < set->number)
      RefreshTextLine(set, i);
    if (y+1 < set->number)
    RefreshTextLine(set, y+1);
  }
  /*
   * As we may have several text sets operational at once, we need to
   * remember which one was last indicated.
   */
  CurrentTextSet = set;
}

int SetTextSetLine(TextSet *set, int which) 
{
  if (which > -1 && which < set->number)
    set->currently_pointed = which;
  else
    return 0;
  return 1;
}

int GetTextSetLine(set) 
     TextSet *set;
{
  if (!set) return -1;
  return set->currently_pointed;
}

void GetTextLineString(TextSet *set, int which, char *str)
{
  memset(str, 0,  strlen(str)); strcpy(str, set->text[which]);
}

void AddTextSetString(TextSet *set, char *str)
{
  if (set->currently_pointed == -1) {
    fprintf(stderr, "ERROR : Currently pointed is -1!\n");
    exit(0);
  }
  strcpy(set->text[set->currently_pointed], str);
  RefreshTextLine(set, set->currently_pointed);
}

TextSet *QueryCurrentTextSet()
{
  return CurrentTextSet;
}

void SetCurrentTextSet(TextSet *set, int dir)
{
  if (CurrentTextSet != NULL)
    RefreshTextSet(CurrentTextSet);
  if (CurrentTextSet == set)
    return;
  if (CurrentTextSet != NULL)
    CurrentTextSet->currently_pointed = -1;
  if (CurrentTextSet != NULL)
    RefreshTextSet(CurrentTextSet);
  if (dir == UP)
    set->currently_pointed = set->number-1;
  else
    set->currently_pointed = 0;
  if (CurrentTextSet != NULL)
    set->char_pos = CurrentTextSet->char_pos;
  else
    set->char_pos = 0;
  CurrentTextSet = set;
  if (set->text[set->currently_pointed][0] == '\0')
    set->char_pos = 0;
  else if (set->char_pos > (((int) strlen(set->text[set->currently_pointed]))-1))
    set->char_pos = (int) strlen(set->text[set->currently_pointed]);
  RefreshTextSet(set);
}


/*****************************************************************************/
/** This procedure adds a character to the currently pointed text line in   **/
/** the indicated set.  If the set is NULL then the CurrentTextSet is used. **/
/*****************************************************************************/
void AddTextSetChar(TextSet *set, char c)
{
  char stemp[180];

  if (!set) 
    set = CurrentTextSet;
  if (set->currently_pointed == -1)
    return;

  memset(stemp, 0, sizeof(stemp));
  strncpy(stemp, set->text[set->currently_pointed], set->char_pos);
  stemp[set->char_pos] = c;
  strcat(stemp, &set->text[set->currently_pointed][set->char_pos]);
  strcpy(set->text[set->currently_pointed], stemp);
  set->char_pos++;
  if (set->char_pos > set->view_start[set->currently_pointed]+set->max)
    set->view_start[set->currently_pointed]++;
  RefreshTextLine(set, set->currently_pointed);
}

void TextLineBack(TextSet *set)
{
  if (!set) 
    set = CurrentTextSet;
  if (set->char_pos == 0) return;
  if (set->char_pos == set->view_start[set->currently_pointed]) {
    set->char_pos --;
    set->view_start[set->currently_pointed] --;
  }
  else
    set->char_pos--;
  RefreshTextLine(set, set->currently_pointed);
}

void TextLineForward(TextSet *set)
{
  if (!set) 
    set = CurrentTextSet;
  if (set->char_pos > (((int) strlen(set->text[set->currently_pointed]))-1))
    return;
  set->char_pos++;
  if ((set->char_pos - set->view_start[set->currently_pointed]) > set->max) 
    set->view_start[set->currently_pointed]++;
  RefreshTextLine(set, set->currently_pointed);
}

int TextLineUp(TextSet *set)
{
  if (!set) 
    set = CurrentTextSet;
  if (set->currently_pointed == 0)
    return -1;
  set->currently_pointed--;
  if (set->char_pos >  (((int) strlen(set->text[set->currently_pointed]))-1))
    set->char_pos = (int) strlen(set->text[set->currently_pointed]);
  RefreshTextLine(set, set->currently_pointed+1);
  RefreshTextLine(set, set->currently_pointed);
  return EXIT_SUCCESS;
}

int TextLineDown(TextSet *set)
{
  if (!set) 
    set = CurrentTextSet;
  if (set->currently_pointed < set->number-1) {
    set->currently_pointed ++;
    if (set->char_pos >  (((int) strlen(set->text[set->currently_pointed]))-1))
      set->char_pos = (int) strlen(set->text[set->currently_pointed]);
    RefreshTextLine(set, set->currently_pointed-1);
    RefreshTextLine(set, set->currently_pointed);
    return 0;
  }
  return -1;
}

void DeleteTextSetChar(TextSet *set)
{
  char stemp[180];
  
  
  if (!set) 
    set = CurrentTextSet;
  if (!set->char_pos) return;
  set->char_pos--;
  memset(stemp, 0, sizeof(stemp));
  if (set->char_pos - set->view_start[set->currently_pointed] < 0)
    set->view_start[set->currently_pointed]--;
  strncpy(stemp, set->text[set->currently_pointed], set->char_pos);
  strcat(stemp, &set->text[set->currently_pointed][set->char_pos+1]);
  strcpy(set->text[set->currently_pointed], stemp);
  RefreshTextLine(set, set->currently_pointed);
}

TextSet *CreateTextSet(Window wind, int startx, int starty, int number,
		       int space, int max, Font font, int char_width,
		       int char_height, unsigned long foreground,
		       unsigned long background)
{
  int loop;
  TextSet *ret;  

  ret = (TextSet *) malloc(sizeof(TextSet));
  ret->win         = wind;
  ret->startx      = startx;
  ret->starty      = starty;
  ret->number      = number;
  ret->space       = space;
  ret->max         = max;
  ret->font        = font;
  ret->char_width  = char_width;
  ret->char_height = char_height;
  ret->currently_pointed = 1;
  ret->char_pos    = 0;
  ret->foreground  = foreground;
  ret->background  = background;
  for (loop = 0; loop < number; loop ++) {
    ret->text[loop][0] = '\0';
    ret->view_start[loop]  = 0;
  }
  CurrentTextSet=ret;
  return ret;
}

