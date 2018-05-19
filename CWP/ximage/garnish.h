/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

#ifndef _GARNISH
#define _GARNISH
#ifndef EGSTERN
#define EGSTERN extern
#endif/*EGSTERN*/

#define MAX_BOX_THICK 30

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

EGSTERN XColor generic;

typedef struct _textset {
  Window win;                /* window text set resides in                   */
  Font font;                 /* font for use in the text set                 */
  int startx;                /* startign x position of the text set          */
  int starty;                /* starting y position of the text set          */
  int number;                /* number of text lines in the set              */
  int space;                 /* space between the text lines                 */
  int max;                   /* max characters displayable at a time         */
  int char_width;            /* character width                              */
  int char_height;           /* character height                             */
  int currently_pointed;     /* current text line being pointed to           */
  int char_pos;              /* current "pointer" position on the line       */
  unsigned long foreground;  /* forground pixel                              */
  unsigned long background;  /* background pixel                             */
  char text[40][200];        /* Text for the text lines                      */
  int view_start[40];        /* position in string of first visible char     */
} TextSet;

EGSTERN TextSet *CurrentTextSet;     /* Currently indicated Text Set                 */

/* Prototypes */
void TextLineBack(TextSet *set);
void TextLineForward(TextSet *set);
int TextLineUp(TextSet *set);
int TextLineDown(TextSet *set);
void DeleteTextSetChar(TextSet *set);
void AddTextSetChar(TextSet *set, char c);
TextSet *CreateTextSet(Window wind, int startx, int starty, int number,
		       int space, int max, Font font, int char_width,
		       int char_height, unsigned long foreground,
		       unsigned long background);
int SetTextSetLine(TextSet *set, int which);
void SetCurrentTextSet(TextSet *set, int dir);
void AddTextSetString(TextSet *set, char *str);
void RefreshTextSet(TextSet *set);
void NewButton(Window win, int x, int y, int width, int height, int state,
	  XColor color, XColor textcolor, char text[], int c_width,
	  int c_height, int mode);
void GetTextLineString(TextSet *set, int which, char *str);
void DrawRadio(int state, int x, int y, XColor color, int size);
int DrawShadowedBox(Window win,XColor color, int thick, int xcorner,
		     int ycorner, int length, int height);
void DisplayRadioButtonSet(int x, int y, int number, int space, int size,
		      int state[], XColor color);
void DrawScrollBar(Window win,int x,int y,int width,int height,int pos,
		   XColor color);
void RefreshTextLine(TextSet *set, int which);
void CheckTextSet(TextSet *set, int x, int y);

#endif
