/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */

#ifndef _PICKING_H
#define _PICKING_H
#ifndef EGSTERN
#define EGSTERN extern
#endif

#define COLOR_FLAGS DoRed | DoGreen | DoBlue
#define COMMAND_WIDTH 100
#define BUTTON_HEIGHT 30
#define BUTTON_WIDTH  90
#define BUTTON_BRIGHTNESS 60000
#define FONT_NAME "9x15"
#define PICK_MODE 1
#define REGULAR_MODE 0
#define ADD_MODE 1
#define DELETE_MODE 0
#define DRAW_FLAG 1
#define ERASE_FLAG 0

EGSTERN int char_width, char_height;
EGSTERN XColor grey_color, black_color,red_color,blue_color;
EGSTERN unsigned long grey_pixel,black_pixel,red_pixel,blue_pixel;
EGSTERN GC blue_r_gc;
EGSTERN GC red_r_gc;

typedef struct pick_tag {
  float x2;
  float time;
  int picked;
} pick_t;

#endif /* end PICKING */
