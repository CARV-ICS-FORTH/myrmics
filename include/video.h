// ============================================================================
// Copyright (c) 2013, FORTH-ICS / CARV 
//                     (Foundation for Research & Technology -- Hellas,
//                      Institute of Computer Science,
//                      Computer Architecture & VLSI Systems Laboratory)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// 
//     http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// 
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : XUP board video in/out functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: video.h,v $
// CVS revision  : $Revision: 1.23 $
// Last modified : $Date: 2013/03/22 12:25:40 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _VIDEO_H
#define _VIDEO_H

// XUP board DRAM
#define VID_DRAM_BASE                 ((unsigned int) 0x00000000)
#define VID_DRAM_SIZE                 (256 * 1024 * 1024) // 256 MB
                                                   // 0x0FFFFFFC last DRAM adr

#define VID_BUF_SIZE                  (4 * 1024 * 1024) // 4 MB per frame
#define VID_NUM_BUFFERS               8

// Video input registers
#define VID_IN_CONTROL                ((unsigned int) 0x40100000)
#define VID_IN_FB_ADR                 ((unsigned int) 0x40100004)
#define VID_IN_FB_INDICES             ((unsigned int) 0x40100008)
#define VID_IN_STATUS                 ((unsigned int) 0x4010000C)
#define VID_IN_HOR_PIPE               ((unsigned int) 0x40100010)
#define VID_IN_HOR_SYNC               ((unsigned int) 0x40100014)
#define VID_IN_HOR_FRONT_PORCH        ((unsigned int) 0x40100018)
#define VID_IN_HOR_BACK_PORCH         ((unsigned int) 0x4010001C)
#define VID_IN_HOR_PIXEL_TIME         ((unsigned int) 0x40100020)
#define VID_IN_VER_SYNC               ((unsigned int) 0x40100024)
#define VID_IN_VER_FRONT_PORCH        ((unsigned int) 0x40100028)
#define VID_IN_VER_BACK_PORCH         ((unsigned int) 0x4010002C)
#define VID_IN_VER_FRAME              ((unsigned int) 0x40100030)
                                       
// Video output registers              
#define VID_OUT_CONTROL               ((unsigned int) 0x40110000)
#define VID_OUT_FB_ADR                ((unsigned int) 0x40110004)
#define VID_OUT_FB_INDICES            ((unsigned int) 0x40110008)
#define VID_OUT_STATUS                ((unsigned int) 0x4011000C)
#define VID_OUT_HOR_PIPE              ((unsigned int) 0x40110010)
#define VID_OUT_HOR_SYNC              ((unsigned int) 0x40110014)
#define VID_OUT_HOR_FRONT_PORCH       ((unsigned int) 0x40110018)
#define VID_OUT_HOR_BACK_PORCH        ((unsigned int) 0x4011001C)
#define VID_OUT_HOR_PIXEL_TIME        ((unsigned int) 0x40110020)
#define VID_OUT_VER_SYNC              ((unsigned int) 0x40110024)
#define VID_OUT_VER_FRONT_PORCH       ((unsigned int) 0x40110028)
#define VID_OUT_VER_BACK_PORCH        ((unsigned int) 0x4011002C)
#define VID_OUT_VER_FRAME             ((unsigned int) 0x40110030)

// I2C
#define VID_I2C_SL                    ((unsigned int) 0x40120000)
#define VID_I2C_SH                    ((unsigned int) 0x40120004)
#define VID_I2C_CT                    ((unsigned int) 0x40120008)
#define VID_I2C_TX_RX                 ((unsigned int) 0x4012000C)
#define VID_I2C_CR_SR                 ((unsigned int) 0x40120010)

#define VID_I2C_VID_IN_SLV_ADR        0x4C
#define VID_I2C_VID_OUT_SLV_ADR       0x76


// Chosen setup for 800x600 in/out frames
#define VID_FULL_WIDTH         800     // Camera 16:9 width = Full frame width
#define VID_STRIP_HEIGHT       448     // Camera 16:9 height
#define VID_FULL_HEIGHT        600     // Full frame height


// Exported core functions
extern void vid_init(int my_bid, int my_cid);

extern void vid_out_reset(int my_bid, int my_cid);
extern void vid_out_set_frame_buffer(int my_bid, int my_cid, void *framebuf);
extern void vid_out_set_frame_mask(int my_bid, int my_cid, unsigned int mask);

extern void vid_in_reset(int my_bid, int my_cid);
extern void vid_in_set_frame_buffer(int my_bid, int my_cid, void *framebuf);
extern void vid_in_set_frame_mask(int my_bid, int my_cid, unsigned int mask);
extern int  vid_in_hold_last_frame(int my_bid, int my_cid, 
                                   unsigned int full_vid_in_mask);

extern void vid_font_putstr(unsigned int *buf, int buf_width, int buf_x, 
                            int buf_y, int font_large, int condense, 
                            unsigned int fg_color, unsigned int bg_color, 
                            int transparency, char *str);

extern void vid_blit_logo(unsigned int *buf, int buf_width, int buf_x, 
                          int buf_y, int which_logo, unsigned int *colors);

extern void vid_draw_line(unsigned int *buf, int buf_width, int buf_height,
                          int x1, int y1, int x2, int y2, unsigned int color,
                          int transparency);
extern void vid_draw_circle(unsigned int *buf, int buf_width, int buf_height, 
                            int x0, int y0, int radius, unsigned int color);
extern void vid_draw_floodfill(unsigned int *buf, int buf_width, int buf_height,
                               int x, int y, unsigned int color, 
                               unsigned int boundary);


// Video demo, MPI version
extern void demo_mpi();

// Video demo, Myrmics version
#define DEMO_MYRMICS_REDKEY_RED_MIN     190  // Define these to detect a 
#define DEMO_MYRMICS_REDKEY_GREEN_MAX   5    // red ball in the camera
#define DEMO_MYRMICS_REDKEY_BLUE_MAX    5 

#define DEMO_MYRMICS_GREENKEY_RED_MAX   10   // Define these to detect a 
#define DEMO_MYRMICS_GREENKEY_GREEN_MIN 120  // green ball in the camera
#define DEMO_MYRMICS_GREENKEY_BLUE_MAX  40

#define DEMO_MYRMICS_BALL_SKETCH_RADIUS 40   // In-out sketching ball radius
#define DEMO_MYRMICS_RED_SKETCH_COLOR   0x0000FF00 // Red ball color (BBGGRRxx)
#define DEMO_MYRMICS_RED_SKETCH_FILL    0x0000DF00 // Red ball fill color
#define DEMO_MYRMICS_GREEN_SKETCH_COLOR 0x00AA0000 // Green ball color (BBGGRRxx)
#define DEMO_MYRMICS_GREEN_SKETCH_FILL  0x008A0000 // Green ball fill color

#define DEMO_MYRMICS_BOUNDARY_COLOR     0xFFFFFF00 // Tile boundary clr (BBGGRRxx)
#define DEMO_MYRMICS_BOUNDARY_TRANSP    48         // Tile boundary transparency

#define DEMO_MYRMICS_MIN_TASKS          1
#define DEMO_MYRMICS_MAX_TASKS          1024

#define DEMO_MYRMICS_OUT_UPDATE_MSEC    100  // Min msec to update output frame
#define DEMO_MYRMICS_START_TIMEOUT      10   // Seconds to wait before starting
#define DEMO_MYRMICS_END_TIMEOUT        3    // Seconds to wait before restarting

#define DEMO_MYRMICS_IN_CNT_XUP         1    // Video input core counter to use
                                             // for XUP transactions

#define DEMO_MYRMICS_OUT_CNT_XUP        20   // Video output core counter to
                                             // use for XUP transactions

typedef struct {

  int   red_present;
  int   red_x;
  int   red_y;

  int   green_present;
  int   green_x;
  int   green_y;

  int   pad[10];        // Fill this so that the struct is 64-B aligned

} DemoMyrmicsInScene;

extern void demo_myrmics_input_loop();
extern void demo_myrmics_output_loop();

extern void demo_myrmics_passthrough();

extern void demo_myrmics_pix_buf_coords(int tile_id, int num_tiles,
                                        int *ret_origin_x, int *ret_origin_y,
                                        int *ret_width, int *ret_height);
extern void demo_myrmics_init(unsigned int *worker_base, int *red_present, 
                              int *red_x, int *red_y, int *green_present, 
                              int *green_x, int *green_y, int *scene, 
                              int *num_tasks);
//extern void demo_myrmics_get_new_rod(DemoMyrmicsState *state);
extern void demo_myrmics_draw_output(unsigned int *pix_buf, int width, 
                                     int height, int origin_x, int origin_y,
                                     unsigned int vid_out_buf_base, int tile_id,
                                     int num_tiles);

extern void (**demo_myrmics_task_table)();

#endif
