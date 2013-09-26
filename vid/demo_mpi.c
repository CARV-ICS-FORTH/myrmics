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
// Abstract      : XUP-based video demo using MPI
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: demo_mpi.c,v $
// CVS revision  : $Revision: 1.17 $
// Last modified : $Date: 2013/04/03 07:59:53 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <video.h>
#include <fmpi.h>


#define DEMO_MPI_REDKEY_RED_MIN     170    // Define these to detect a 
#define DEMO_MPI_REDKEY_GREEN_MAX   10     // red ball in the camera
#define DEMO_MPI_REDKEY_BLUE_MAX    10

#define DEMO_MPI_GREENKEY_RED_MAX   10     // Define these to detect a 
#define DEMO_MPI_GREENKEY_GREEN_MIN 120    // green ball in the camera
#define DEMO_MPI_GREENKEY_BLUE_MAX  40

#define DEMO_MPI_AREA_WIDTH     (VID_FULL_WIDTH - 2 * DEMO_MPI_BORDER - \
                                 2 * DEMO_MPI_BALL_RADIUS)
#define DEMO_MPI_AREA_HEIGHT    (VID_STRIP_HEIGHT - 2 * DEMO_MPI_BORDER - \
                                 2 * DEMO_MPI_BALL_RADIUS)

#define DEMO_MPI_BORDER         4          // Border to draw in pixels
#define DEMO_MPI_BORDER_COLOR   0xA0A0A000 // Border color in BBGGRRxx format
#define DEMO_MPI_CORE_SEP_COLOR 0x60606000 // Core separator color in BBGGRRxx

#define DEMO_MPI_ROD_COLOR      0x00C00000 // Rod color in BBGGRRxx format
#define DEMO_MPI_ROD_FILL_COLOR 0x00300000 // Rod fill color in BBGGRRxx format
#define DEMO_MPI_ROD_RADIUS     80         // Drawn rod thickness in pixels

#define DEMO_MPI_BALL_RADIUS    3          // Ball radius in pixels

#define DEMO_MPI_MAX_ABS_VEL    8          // Maximum absolute velocity
#define DEMO_MPI_COL_EXTRA_VEL  1.1f       // Added velocity after collision

#define DEMO_MPI_MIN_BALLS      4          // Minimum possible balls
#define DEMO_MPI_MAX_BALLS      512        // Maximum possible balls
#define DEMO_MPI_SCENE_ARRAY    342        // DEMO_MPI_MAX_BALLS * 2/3 round up
#define DEMO_MPI_MAX_XFER_BALLS 16         // Maximum balls allowed to be
                                           // transferred among neighbors in a
                                           // single time step
#define DEMO_MPI_MAX_WORKERS    512        // Maximum possible workers

#define DEMO_MPI_CMD_BUFS       4          // Number of command sets for video
                                           // input core
#define DEMO_MPI_OUT_BUFS       16         // Number of buffer sets for video
                                           // output core

#define DEMO_MPI_IN_CNT_XUP     1          // Video input core counter to use
                                           // for XUP transactions

#define DEMO_MPI_OUT_CNT_BASE   0          // Video output core counter base
                                           // to receive incoming buffers.
                                           // DEMO_MPI_OUT_BUFS will be used
                                           // in total, starting from
                                           // DEMO_MPI_OUT_CNT_BASE.
#define DEMO_MPI_OUT_CNT_XUP    20         // Video output core counter to
                                           // use for XUP transactions


// ===========================================================================
// This needs the CFG_512_MB setup (512 MB). The 4 ARM arm cores (ARM0 board)
// must be present, but they communicate outside the MPI setup.
// ===========================================================================

typedef struct {

  float pos_x;
  float pos_y;

  float vel_x;
  float vel_y;

  int   valid;  // used only during transfers

} DemoMPIBall;

typedef struct {

  int   num_slaves;
  float blk_width;
  float blk_height;
  int   tiles_x;
  int   tiles_y;
  int   rod_present;
  int   rod_x0;
  int   rod_y0;
  int   init_balls;
  int   init_num_balls;

  int   pad[6];   // Make sure struct is a multiple of 64 B

} DemoMPICommand;


typedef struct {

  unsigned int  balls           // Array of balls, 2/3 integer per ball:
        [DEMO_MPI_SCENE_ARRAY]; // pos 0:    x1 | y0 | x0
                                // pos 1:    y2 | x2 | y1
                                // pos 2:    x3 | y2 | x2
                                // pos 3:    y4 | x4 | y3
                                // ...
                                // (where each xN/yN is 10 bits)
  
  int           active_balls;   // Number of balls to draw from the array
  
  int           pad[9];         // Make sure struct is a multiple of 64 B

} DemoMPIOutScene;

// ===========================================================================
// ===========================================================================
static inline int cnt_id(int frame) {

  switch (frame % DEMO_MPI_CMD_BUFS) {
    case 0: return NOC_COUNTER_WAKEUP0;
    case 1: return NOC_COUNTER_WAKEUP1;
    case 2: return NOC_COUNTER_WAKEUP2;
    case 3: return NOC_COUNTER_WAKEUP3;
  }
  // Must be below 4, we only have 4 counters not used by the MPI library
  ar_abort();
}


// ===========================================================================
// ===========================================================================
static inline unsigned int scene_buf(unsigned int scene_base, int set,
                                     int rank) {

  return scene_base + (set * DEMO_MPI_MAX_WORKERS + rank) * 
                                                sizeof(DemoMPIOutScene);
}

// ===========================================================================
// ===========================================================================
int demo_mpi_user_menu(int *num_balls, int *config) {

  int   change;
  int   c;
  char  *menu = "\r\n"
                "Current settings: %d balls, %d worker cores\r\n"
                "\r\n"
                "Select an option:\r\n"
                "[1] Increase number of balls\r\n"
                "[2] Decrease number of balls\r\n"
                "[3] Increase worker cores\r\n"
                "[4] Decrease worker cores\r\n";
  char  *err_b = "\r\n** Error: Number of balls must be between %d and %d **\r";
  char  *err_c = "\r\n** Error: Worker cores must be between 1 and 512 **\r";

  change = 0;

  if (*config == -1) {
    *config = 0;
    *num_balls = DEMO_MPI_MIN_BALLS;
    kt_printf(menu, *num_balls, kt_int_pow(2, *config));
    change = 1;
  }

  c = ar_uart_get_char();
  if ((c >= '1') && (c <= '4')) {
    change = 1;
    switch (c) {

      case '1':
        if (*num_balls < DEMO_MPI_MAX_BALLS) {
          *num_balls *= 2;
        }
        else {
          change = 0;
          kt_printf(err_b, DEMO_MPI_MIN_BALLS, DEMO_MPI_MAX_BALLS);
        }
        break;

      case '2':
        if (*num_balls > DEMO_MPI_MIN_BALLS) {
          *num_balls /= 2;
        }
        else {
          change = 0;
          kt_printf(err_b, DEMO_MPI_MIN_BALLS, DEMO_MPI_MAX_BALLS);
        }
        break;

      case '3':
        if (*config < 9) {
          (*config)++;
        }
        else {
          change = 0;
          kt_printf(err_c);
        }
        break;

      case '4':
        if (*config > 0) {
          (*config)--;
        }
        else {
          change = 0;
          kt_printf(err_c);
        }
        break;
    }
    
    kt_printf(menu, *num_balls, kt_int_pow(2, *config));
  }

  return change;
}


// ===========================================================================
// ===========================================================================
void demo_mpi_input_loop(int my_bid, int my_cid) {

  unsigned int            *buf[2];
  volatile DemoMPICommand *cmd[DEMO_MPI_CMD_BUFS];
  unsigned int            buf_base[DEMO_MPI_MAX_WORKERS + 1];
  int                     cur;
  int                     in;
  int                     frame;
  int                     cur_cmd;
  int                     i;
  int                     x;
  int                     y;
  int                     x0;
  int                     y0;
  int                     min_x;
  int                     min_y;
  int                     max_x;
  int                     max_y;
  unsigned int            red;
  unsigned int            green;
  unsigned int            blue;
  unsigned int            *p;
  unsigned int            word;
  int                     peer_bid;
  int                     peer_cid;
  int                     num_balls;
  int                     config;
  int                     change;


  // Initialize XUP video
  vid_init(my_bid, my_cid);
  vid_in_set_frame_mask(my_bid, my_cid,  0xF0);
  vid_out_set_frame_mask(my_bid, my_cid, 0x01);

  // Initialize two video buffers
  for (i = 0; i < 2; i++) {
    buf[i] = kt_malloc(VID_FULL_WIDTH * VID_FULL_HEIGHT * 4);
  }

  // Initialize command buffers and command counters
  for (i = 0; i < DEMO_MPI_CMD_BUFS; i++) {
    cmd[i] = kt_malloc(sizeof(DemoMPICommand));
    ar_cnt_set(my_cid, cnt_id(i), 0);
  }


  // Discover buffer bases for all MPI cores plus the video out core, 
  // so we can push commands to them
  for (i = 0; i < DEMO_MPI_MAX_WORKERS + 1; i++) {
    word = ar_mbox_get(my_cid);
    peer_bid = (word >> 8) & 0xFF;
    peer_cid = (word >> 0) & 0xFF;
    if (peer_bid < 64) {
      buf_base[peer_bid * 8 + peer_cid] = ar_mbox_get(my_cid);
    }
    else if (peer_bid == AR_ARM0_BID) {
      ar_assert(peer_cid == 1);
      buf_base[DEMO_MPI_MAX_WORKERS] = ar_mbox_get(my_cid);
    }
    else {
      ar_abort();
    }
  }

  // Loop
  cur = 0;
  frame = 0;
  config = -1;
  num_balls = -1;
  while (1) {

    // Find last stable input frame and hold it.
    in = vid_in_hold_last_frame(my_bid, my_cid, 0xF0);
    ar_assert(in <= 7);
    if (in < 4) {
      // No frame yet
      continue;
    }

    // Copy it from XUP to our working buffer
    ar_receive_xup_frame(my_bid, my_cid, (unsigned int) buf[cur], 
                         VID_DRAM_BASE + (in * VID_BUF_SIZE), 
                         DEMO_MPI_IN_CNT_XUP);
    while (ar_cnt_get(my_cid, DEMO_MPI_IN_CNT_XUP)) {
      ;
    }

    // Return frame to video input, we've now copied it
    vid_in_set_frame_mask(my_bid, my_cid, 0xF0);

    // Do a keying pass to locate a red ball. Do the analysis on pixel 
    // by pixel and line by line, so we don't devote too much time.
    min_x = 9999;
    min_y = 9999;
    max_x = -1;
    max_y = -1;
    p = buf[cur] + VID_FULL_WIDTH * (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 +
                   VID_FULL_WIDTH / 4;
    for (y = (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2;
         y < (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + VID_STRIP_HEIGHT;
         y += 1, p += 0 * VID_FULL_WIDTH + VID_FULL_WIDTH / 2) {
      for (x = VID_FULL_WIDTH / 4; x < (VID_FULL_WIDTH * 3) / 4;
           x += 1, p += 1) {
        red   = (*p >> 8 ) & 0xFF;
        green = (*p >> 16) & 0xFF;
        blue  = (*p >> 24) & 0xFF;
        if ((red   <= DEMO_MPI_GREENKEY_RED_MAX) &&
            (green >= DEMO_MPI_GREENKEY_GREEN_MIN) &&
            (blue  <= DEMO_MPI_GREENKEY_BLUE_MAX)) {
          if (x < min_x) {
            min_x = x;
          }
          if (y < min_y) {
            min_y = y;
          }
          if (x > max_x) {
            max_x = x;
          }
          if (y > max_y) {
            max_y = y;
          }
        }
      }
    }

    // Can we reuse the command buffer?
    if (frame >= DEMO_MPI_CMD_BUFS) {
      while (ar_cnt_get(my_cid, cnt_id(frame - DEMO_MPI_CMD_BUFS)) != 
             DEMO_MPI_MAX_WORKERS + 1) {
        ;
      }
      ar_cnt_set(my_cid, cnt_id(frame - DEMO_MPI_CMD_BUFS), 0);
    }

    // Prepare command. The camera mirrors the image, so to
    // match the rod as the user would expect it, we also mirror it. 
    // Additionally, we translate the coordinates from the full frame system
    // to the in-border drawing area.
    cur_cmd = frame % DEMO_MPI_CMD_BUFS;
    if (min_x < 9999) {

      x = ((max_x + min_x) / 2 - VID_FULL_WIDTH / 4) * 2;
      y = (max_y + min_y) / 2;

      // Mirror
      x0 = VID_FULL_WIDTH - x;
      y0 = y;

      // Translate and clamp
      x0 -= DEMO_MPI_BORDER;
      y0 -= (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + DEMO_MPI_BORDER;
      if (x0 < DEMO_MPI_ROD_RADIUS) {
        x0 = DEMO_MPI_ROD_RADIUS;
      }
      if (y0 < DEMO_MPI_ROD_RADIUS) {
        y0 = DEMO_MPI_ROD_RADIUS;
      }
      if (x0 >= DEMO_MPI_AREA_WIDTH - DEMO_MPI_ROD_RADIUS) {
        x0 = DEMO_MPI_AREA_WIDTH - 1 - DEMO_MPI_ROD_RADIUS;
      }
      if (y0 >= DEMO_MPI_AREA_HEIGHT - DEMO_MPI_ROD_RADIUS) {
        y0 = DEMO_MPI_AREA_HEIGHT - 1 - DEMO_MPI_ROD_RADIUS;
      }

      // Fill buffer
      cmd[cur_cmd]->rod_present = 1;
      cmd[cur_cmd]->rod_x0 = x0;
      cmd[cur_cmd]->rod_y0 = y0;
    }
    else {
      cmd[cur_cmd]->rod_present = 0;
    }

    // Get worker configuration and number of balls from a user menu
    change = demo_mpi_user_menu(&num_balls, &config);

    // Add slave configuration to the command
    switch (config) {

      case 0:
        cmd[cur_cmd]->num_slaves = 1;
        cmd[cur_cmd]->tiles_x    = 1;
        cmd[cur_cmd]->tiles_y    = 1;
        break;
        
      case 1:
        cmd[cur_cmd]->num_slaves = 2;
        cmd[cur_cmd]->tiles_x    = 2;
        cmd[cur_cmd]->tiles_y    = 1;
        break;

      case 2:
        cmd[cur_cmd]->num_slaves = 4;
        cmd[cur_cmd]->tiles_x    = 2;
        cmd[cur_cmd]->tiles_y    = 2;
        break;
        
      case 3:
        cmd[cur_cmd]->num_slaves = 8;
        cmd[cur_cmd]->tiles_x    = 4;
        cmd[cur_cmd]->tiles_y    = 2;
        break;

      case 4:
        cmd[cur_cmd]->num_slaves = 16;
        cmd[cur_cmd]->tiles_x    = 4;
        cmd[cur_cmd]->tiles_y    = 4;
        break;

      case 5:
        cmd[cur_cmd]->num_slaves = 32;
        cmd[cur_cmd]->tiles_x    = 8;
        cmd[cur_cmd]->tiles_y    = 4;
        break;

      case 6:
        cmd[cur_cmd]->num_slaves = 64;
        cmd[cur_cmd]->tiles_x    = 8;
        cmd[cur_cmd]->tiles_y    = 8;
        break;

      case 7:
        cmd[cur_cmd]->num_slaves = 128;
        cmd[cur_cmd]->tiles_x    = 16;
        cmd[cur_cmd]->tiles_y    = 8;
        break;

      case 8:
        cmd[cur_cmd]->num_slaves = 256;
        cmd[cur_cmd]->tiles_x    = 16;
        cmd[cur_cmd]->tiles_y    = 16;
        break;

      case 9:
        cmd[cur_cmd]->num_slaves = 512;
        cmd[cur_cmd]->tiles_x    = 32;
        cmd[cur_cmd]->tiles_y    = 16;
        break;

      default:
        ar_abort();
    }

    cmd[cur_cmd]->blk_width  = (float) DEMO_MPI_AREA_WIDTH / 
                               (float) cmd[cur_cmd]->tiles_x;
    cmd[cur_cmd]->blk_height = (float) DEMO_MPI_AREA_HEIGHT /
                               (float) cmd[cur_cmd]->tiles_y;
    
    // Was there any change? Tell the workers to initialize their balls.
    if (change) {
      cmd[cur_cmd]->init_balls = 1;
      cmd[cur_cmd]->init_num_balls = num_balls / cmd[cur_cmd]->num_slaves;
      if (!cmd[cur_cmd]->init_num_balls) {
        cmd[cur_cmd]->init_num_balls = -num_balls; // mark 1 ball until ranks 
                                                   // reach num_balls
      }
    }
    else {
      cmd[cur_cmd]->init_balls = 0;
    }


    // Send command to all MPI cores, as well as the video output core
    for (i = 0; i < DEMO_MPI_MAX_WORKERS; i++) {
      while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
        ;
      }
      ar_dma_with_ack(my_cid,
                      my_bid, my_cid, (unsigned int) cmd[cur_cmd],
                      i / 8, i % 8, 
                      buf_base[i] + cur_cmd * sizeof(DemoMPICommand),
                      i / 8, i % 8, cnt_id(frame),
                      sizeof(DemoMPICommand), 0, 0, 0);
    }
    while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
      ;
    }
    ar_dma_with_ack(my_cid,
                    my_bid, my_cid, (unsigned int) cmd[cur_cmd],
                    AR_ARM0_BID, 1, 
                    buf_base[DEMO_MPI_MAX_WORKERS] + cur_cmd * 
                                                        sizeof(DemoMPICommand),
                    AR_ARM0_BID, 1, cnt_id(frame),
                    sizeof(DemoMPICommand), 0, 0, 0);

    // Change buffer
    cur = (cur + 1) % 2;
    frame++;
  }
}


// ===========================================================================
// ===========================================================================
void demo_mpi_output_loop(int my_bid, int my_cid) {

  unsigned int                  *buf[3];
  int                           cur_buf;
  int                           cur_out;
  int                           cur_cmd;
  int                           cur_scene;
  int                           x;
  int                           y;
  int                           bx;
  int                           by;
  char                          msg[128];
  int                           frame;
  float                         fps;
  unsigned int                  elapsed;
  unsigned int                  *p;
  unsigned int                  *q;
  int                           i;
  int                           j;
  int                           idx;
  unsigned int                  *ball_buf[DEMO_MPI_MAX_WORKERS];
  int                           total_balls;
  unsigned int                  *user_ball_buf;
  volatile DemoMPICommand       *cmd;
  unsigned int                  cmd_base;
  unsigned int                  scene_base;
  volatile DemoMPIOutScene      *scene;
  unsigned int                  word;
  int                           peer_bid;
  int                           peer_cid;
  int                           r;
  int                           g;
  int                           b;
  unsigned int                  col_outer;
  unsigned int                  col_inner;
  unsigned int                  seed;
  unsigned int                  forth_colors[] = {0x20A02000, 0x18901800,
                                                  0x10601000, 0x00000000};


  // Allocate clean video buffers
  for (i = 0; i < 3; i++) {
    buf[i] = kt_zalloc(VID_FULL_WIDTH * VID_FULL_HEIGHT * 4);
  }

  // Allocate incoming command buffers
  cmd_base = (unsigned int) kt_malloc(DEMO_MPI_CMD_BUFS * sizeof(DemoMPICommand));

  // Allocate incoming scene buffers
  scene_base = (unsigned int) kt_malloc(DEMO_MPI_OUT_BUFS * 
                                        DEMO_MPI_MAX_WORKERS * 
                                        sizeof(DemoMPIOutScene));

  // Initialize our command counters
  for (i = 0; i < DEMO_MPI_CMD_BUFS; i++) {
    ar_cnt_set(my_cid, cnt_id(i), 0);
  }

  // Initialize our scene counters
  for (i = 0; i < DEMO_MPI_OUT_BUFS; i++) {
    ar_cnt_set(my_cid, DEMO_MPI_OUT_CNT_BASE + i, 0);
  }

  // Send our command buffer base to the video input core
  ar_mbox_send2(my_cid, AR_ARM0_BID, 0, (my_bid << 8) | my_cid, cmd_base);

  // Put static stamps and borders on the video buffers
  for (i = 0; i < 3; i++) {
    vid_blit_logo(buf[i], VID_FULL_WIDTH, 0, 524, 1, forth_colors);
    
    vid_font_putstr(buf[i], VID_FULL_WIDTH, 82, 535, 0, 1, 0x20A02000, 0, 0, 
                    "512-core Prototype MPI Demonstration");
    vid_font_putstr(buf[i], VID_FULL_WIDTH, 82, 547, 0, 1, 0x20A02000, 0, 0, 
                    "FORTH-ICS/CARV (c) 2012-2013");

    // Upper border
    p = buf[i] + VID_FULL_WIDTH * (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2;
    for (y = 0; y < DEMO_MPI_BORDER; y++) {
      for (x = 0; x < VID_FULL_WIDTH; x++) {
        *p++ = DEMO_MPI_BORDER_COLOR;
      }
    }

    // Left/right borders
    for (y = 0; y < VID_STRIP_HEIGHT - 2 * DEMO_MPI_BORDER; y++) {
      for (x = 0; x < DEMO_MPI_BORDER; x++) {
        *p++ = DEMO_MPI_BORDER_COLOR;
      }
      p += VID_FULL_WIDTH - 2 * DEMO_MPI_BORDER;
      for (x = 0; x < DEMO_MPI_BORDER; x++) {
        *p++ = DEMO_MPI_BORDER_COLOR;
      }
    }

    // Lower border
    for (y = 0; y < DEMO_MPI_BORDER; y++) {
      for (x = 0; x < VID_FULL_WIDTH; x++) {
        *p++ = DEMO_MPI_BORDER_COLOR;
      }
    }
  }

  // Prepare balls in predrawn buffers. Each ball gets a random color.
  seed = 1;
  for (i = 0; i < DEMO_MPI_MAX_WORKERS; i++) {
    r = ((seed = kt_rand(seed)) % 128) + 128;
    g = ((seed = kt_rand(seed)) % 128) + 128;
    b = ((seed = kt_rand(seed)) % 128) + 128;
    col_outer = (b << 24) | (g << 16) | (r << 8);
    r -= 32;
    g -= 32;
    b -= 32;
    col_inner = (b << 24) | (g << 16) | (r << 8);

    ball_buf[i] = kt_zalloc((2 * DEMO_MPI_BALL_RADIUS + 1) * 
                            (2 * DEMO_MPI_BALL_RADIUS + 1) * sizeof(int));

    vid_draw_circle(ball_buf[i], 2 * DEMO_MPI_BALL_RADIUS + 1, 
                    2 * DEMO_MPI_BALL_RADIUS + 1, 
                    DEMO_MPI_BALL_RADIUS, DEMO_MPI_BALL_RADIUS, 
                    DEMO_MPI_BALL_RADIUS, col_outer);
    vid_draw_floodfill(ball_buf[i], 2 * DEMO_MPI_BALL_RADIUS + 1, 
                       2 * DEMO_MPI_BALL_RADIUS + 1, 
                       DEMO_MPI_BALL_RADIUS, DEMO_MPI_BALL_RADIUS, 
                       col_inner, col_outer);
  }

  // Prepare the big, user-controlled, green ball 
  user_ball_buf = kt_zalloc((2 * DEMO_MPI_ROD_RADIUS + 1) * 
                            (2 * DEMO_MPI_ROD_RADIUS + 1) * sizeof(int));

  vid_draw_circle(user_ball_buf, 2 * DEMO_MPI_ROD_RADIUS + 1, 
                  2 * DEMO_MPI_ROD_RADIUS + 1, 
                  DEMO_MPI_ROD_RADIUS, DEMO_MPI_ROD_RADIUS, 
                  DEMO_MPI_ROD_RADIUS, col_outer);
  vid_draw_floodfill(user_ball_buf, 2 * DEMO_MPI_ROD_RADIUS + 1, 
                     2 * DEMO_MPI_ROD_RADIUS + 1, 
                     DEMO_MPI_ROD_RADIUS, DEMO_MPI_ROD_RADIUS, 
                     col_inner, col_outer);

  // Distribute our scene buffer base to all workers
  for (i = 0; i < DEMO_MPI_MAX_WORKERS; i++) {
    word = ar_mbox_get(my_cid);
    peer_bid = (word >> 8) & 0xFF;
    peer_cid = (word >> 0) & 0xFF;
    ar_assert(peer_bid < 64);
    ar_mslot_send(my_cid, peer_bid, peer_cid, scene_base);
  }

  // Loop
  frame = 0;
  cur_scene = 0;
  cur_cmd = 0;
  cur_buf = 0;
  cur_out = 0;
  elapsed = 1;
  while (1) {
    
    // Wait for a command from the input, so we know the number of active slaves
    while (ar_cnt_get(my_cid, cnt_id(cur_cmd)) != sizeof(DemoMPICommand)) {
      ;
    }

    cmd = (DemoMPICommand *) (cmd_base + cur_cmd * sizeof(DemoMPICommand));

    // Wait for a scene batch to arrive
    while (ar_cnt_get(my_cid, DEMO_MPI_OUT_CNT_BASE + cur_scene) != 
           cmd->num_slaves * sizeof(DemoMPIOutScene)) {
      ;
    }

    ar_cnt_set(my_cid, DEMO_MPI_OUT_CNT_BASE + cur_scene, 0);
    
    

    // Draw core boundaries
    for (x = cmd->blk_width; 
         x <= (cmd->tiles_x - 1) * cmd->blk_width; 
         x += cmd->blk_width) {
      vid_draw_line(buf[cur_buf], VID_FULL_WIDTH, VID_FULL_HEIGHT, 
                    x + DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS, 
                    (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + DEMO_MPI_BORDER,
                    x + DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS, 
                    (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + 
                                        VID_STRIP_HEIGHT - DEMO_MPI_BORDER - 1,
                    DEMO_MPI_CORE_SEP_COLOR, 0);
    }
    for (y = cmd->blk_height; 
         y <= (cmd->tiles_y - 1) * cmd->blk_height; 
         y += cmd->blk_height) {
      vid_draw_line(buf[cur_buf], VID_FULL_WIDTH, VID_FULL_HEIGHT, 
                    DEMO_MPI_BORDER,
                    y + (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + 
                        DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS,
                    VID_FULL_WIDTH - DEMO_MPI_BORDER - 1,
                    y + (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 +
                        DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS,
                    DEMO_MPI_CORE_SEP_COLOR, 0);
    }

    // Draw rod, if it exists
    if (cmd->rod_present) {

      // Translate rod coordinates
      cmd->rod_x0 += DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS;
      cmd->rod_y0 += (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + 
                          DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS;

      // Blit the user ball from the predrawn buffer
      p = buf[cur_buf] + (cmd->rod_y0 - DEMO_MPI_ROD_RADIUS) * VID_FULL_WIDTH +
                         (cmd->rod_x0 - DEMO_MPI_ROD_RADIUS);
      q = user_ball_buf;
      for (y = 0; 
           y < 2 * DEMO_MPI_ROD_RADIUS + 1; 
           y++, p += VID_FULL_WIDTH - 2 * DEMO_MPI_ROD_RADIUS - 1) {
        for (x = 0; 
             x < 2 * DEMO_MPI_ROD_RADIUS + 1; 
             x++, p++, q++) {
          if (*q) {
            *p = *q;
          }
        }
      }
    }

    // Draw each worker's balls 
    total_balls = 0;
    for (j = 0; j < cmd->num_slaves; j++) {

      scene = (DemoMPIOutScene *) scene_buf(scene_base, cur_scene, j);
      for (i = 0, idx = 0; i < scene->active_balls; i++, total_balls++) {

        // Extract ball info and translate coordinates
        switch (i % 3) {
          case 0: 
            bx = scene->balls[idx] & 0x3FF; 
            by = (scene->balls[idx] >> 10) & 0x3FF; 
            break;
          case 1: 
            bx = (scene->balls[idx++] >> 20) & 0x3FF; 
            by = scene->balls[idx] & 0x3FF; 
            break;
          case 2: 
            bx = (scene->balls[idx] >> 10) & 0x3FF; 
            by = (scene->balls[idx++] >> 20) & 0x3FF; 
            break;
          default:
            ar_abort();
        }
        ar_assert(bx < DEMO_MPI_AREA_WIDTH);
        ar_assert(by < DEMO_MPI_AREA_HEIGHT);
        bx += DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS;
        by += (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 +
              DEMO_MPI_BORDER + DEMO_MPI_BALL_RADIUS;

        // Blit the ball from the predrawn buffer for this worker
        p = buf[cur_buf] + (by - DEMO_MPI_BALL_RADIUS) * VID_FULL_WIDTH +
                           (bx - DEMO_MPI_BALL_RADIUS);
        q = ball_buf[j];
        for (y = 0; 
             y < 2 * DEMO_MPI_BALL_RADIUS + 1; 
             y++, p += VID_FULL_WIDTH - 2 * DEMO_MPI_BALL_RADIUS - 1) {
          for (x = 0; 
               x < 2 * DEMO_MPI_BALL_RADIUS + 1; 
               x++, p++, q++) {
            if (*q) {
              *p = *q;
            }
          }
        }
      }
    }

    // Stamp buffer
    kt_sprintf(msg, "%d balls, %d workers     ", 
               total_balls, cmd->num_slaves);
    vid_font_putstr(buf[cur_buf], VID_FULL_WIDTH, 10, 10, 1, 0, 0x00FFFF00, 
                    0, 0, msg);
//kt_printf("%s\r\n", msg);

    fps = 1000.0 / (float) elapsed;
    kt_sprintf(msg, "Frame time: %d msec | FPS: %d.%02d   ", 
               elapsed,
               (int) fps, (int) ((fps - (int) fps) * 100.0));
    vid_font_putstr(buf[cur_buf], VID_FULL_WIDTH, 10, 42, 1, 0, 0x00FFFF00, 
                    0, 0, msg);
//kt_printf("%s\r\n", msg);

    // Signal to the input core we're ready for a new command
    ar_cnt_set(my_cid, cnt_id(cur_cmd), 0);
    ar_cnt_incr(my_cid, AR_ARM0_BID, 0, cnt_id(cur_cmd), 1);

    // Send current buffer to XUP
    ar_send_xup_frame(my_bid, my_cid, (unsigned int) buf[cur_buf], 
                      VID_DRAM_BASE + (cur_out * VID_BUF_SIZE), 
                      DEMO_MPI_OUT_CNT_XUP);

    // Clean next buffer while XUP transmission is under way
    p = buf[(cur_buf + 1) % 3] + VID_FULL_WIDTH *
                  ((VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 +
                   DEMO_MPI_BORDER) + DEMO_MPI_BORDER;
    for (y = 0; y < DEMO_MPI_AREA_HEIGHT + 2 * DEMO_MPI_BALL_RADIUS; y++) {
      for (x = 0; x < DEMO_MPI_AREA_WIDTH + 2 * DEMO_MPI_BALL_RADIUS; x++) {
        *p++ = 0;
      }
      p += 2 * DEMO_MPI_BORDER;
    }

    // Wait for XUP frame to arrive
    while (ar_cnt_get(my_cid, DEMO_MPI_OUT_CNT_XUP)) {
      ;
    }

    // Display frame
    vid_out_set_frame_mask(my_bid, my_cid, 1 << cur_out);

    // Next
    cur_out = (cur_out + 1) % 4;

    cur_buf = (cur_buf + 1) % 3;

    cur_cmd++;
    if (cur_cmd >= DEMO_MPI_CMD_BUFS) {
      cur_cmd = 0;
    }

    cur_scene++; 
    if (cur_scene >= DEMO_MPI_OUT_BUFS) {
      cur_scene = 0;
    }

    // Update stats
    frame++;
    elapsed = ar_timer_get_msec();
    ar_timer_reset();
  }
}


// ===========================================================================
// ===========================================================================
void demo_mpi_new_pos(float old_pos_x, float old_pos_y, float old_vel_x,
                      float old_vel_y, float *ret_pos_x, float *ret_pos_y, 
                      float *ret_vel_x, float *ret_vel_y) {

  // Compute ideal next position and velocity, if nothing else interferes
  *ret_pos_x = old_pos_x + old_vel_x / ((float) DEMO_MPI_MAX_ABS_VEL);
  *ret_pos_y = old_pos_y + old_vel_y / ((float) DEMO_MPI_MAX_ABS_VEL);

  *ret_vel_x = old_vel_x;
  *ret_vel_y = old_vel_y;

  // Does the new position hit onto any wall? Undo the move and inverse
  // the velocity.
  if ((*ret_pos_x < 0.0f) || (*ret_pos_x >= (float) DEMO_MPI_AREA_WIDTH)) {
    *ret_vel_x = -*ret_vel_x;
    *ret_pos_x = old_pos_x;
  }
  if ((*ret_pos_y < 0.0f) || (*ret_pos_y >= (float) DEMO_MPI_AREA_HEIGHT)) {
    *ret_vel_y = -*ret_vel_y;
    *ret_pos_y = old_pos_y;
  }
}


// ===========================================================================
// ===========================================================================
void demo_mpi_collision(float *my_pos_x, float *my_pos_y, float *my_vel_x,
                        float *my_vel_y, float my_radius, float peer_pos_x, 
                        float peer_pos_y, float peer_vel_x, 
                        float peer_vel_y, float peer_radius) {
  float dist_x;
  float dist_y;
  float dist;
  float margin;
  float my_proj_initial;
  float peer_proj_initial;
  float my_proj_final;


  // Are we colliding?
  dist_x = *my_pos_x - peer_pos_x;
  dist_y = *my_pos_y - peer_pos_y;
  dist = ar_float_sqrt(dist_x * dist_x + dist_y * dist_y);
  if (dist == 0.0f) {
    // Correct corner case where the two balls merge into one
    dist_x = 1.0f;
    dist_y = 0.0f;
    dist = 1.0f;
  }
  margin = dist - my_radius - peer_radius;
  if (margin > 0.0f) {
    return;
  }

  // Undo move
  //my_pos_x = state->balls[state->idx][idx]->pos_x;
  //my_pos_y = state->balls[state->idx][idx]->pos_y;

  // Compute the projection of the two velocities onto the axis defined
  // by the two centres of the balls
  dist_x /= dist;
  dist_y /= dist;
  my_proj_initial   = *my_vel_x  * dist_x + *my_vel_y  * dist_y;
  peer_proj_initial = peer_vel_x * dist_x + peer_vel_y * dist_y;

  // After the collision, projected velocities are inversed. We also add 
  // a little margin, to make sure we avoid balls that merge together if 
  // their velocities are rounded badly.
  my_proj_final   = peer_proj_initial * DEMO_MPI_COL_EXTRA_VEL;

  // Correct velocity vector by the amount of the collision
  *my_vel_x += (my_proj_final - my_proj_initial) * dist_x;
  *my_vel_y += (my_proj_final - my_proj_initial) * dist_y;

  // Correct our position to be outside of the other ball. We do this by
  // half if our peer is another ball, because it'll do the same correction.
  // We correct fully if our peer is the rod, because it won't move by itself.
  if (peer_radius == DEMO_MPI_BALL_RADIUS) {
    margin /= 2.0f;
  }
  *my_pos_x -= margin * dist_x;
  *my_pos_y -= margin * dist_y;

  // Clamp to min/max velocities
  if (*my_vel_x < -DEMO_MPI_MAX_ABS_VEL) {
    *my_vel_x = -DEMO_MPI_MAX_ABS_VEL;
  }
  if (*my_vel_y < -DEMO_MPI_MAX_ABS_VEL) {
    *my_vel_y = -DEMO_MPI_MAX_ABS_VEL;
  }
  if (*my_vel_x > DEMO_MPI_MAX_ABS_VEL) {
    *my_vel_x = DEMO_MPI_MAX_ABS_VEL;
  }
  if (*my_vel_y > DEMO_MPI_MAX_ABS_VEL) {
    *my_vel_y = DEMO_MPI_MAX_ABS_VEL;
  }
  
  // Clamp to min/max allowed ball positions
  if (*my_pos_x < 0.0f) {
    *my_pos_x = 0.0f;
  }
  if (*my_pos_x >= (float) DEMO_MPI_AREA_WIDTH) {
    *my_pos_x = (float) DEMO_MPI_AREA_WIDTH - 1.0f;
  }
  if (*my_pos_y < 0.0f) {
    *my_pos_y = 0.0f;
  }
  if (*my_pos_y >= (float) DEMO_MPI_AREA_HEIGHT) {
    *my_pos_y = (float) DEMO_MPI_AREA_HEIGHT - 1.0f;
  }
}

// ===========================================================================
// ===========================================================================
void demo_mpi_do_step(int rank, int north_rank, int east_rank, int south_rank, 
                      int west_rank, float my_min_x, float my_min_y, 
                      float my_max_x, float my_max_y, DemoMPIBall *my_balls, 
                      int *my_num_balls, DemoMPIBall **comm_balls,
                      int rod_present, int rod_x0, int rod_y0) {

  float         my_new_pos_x;
  float         my_new_pos_y;
  float         my_new_vel_x;
  float         my_new_vel_y;
  float         peer_new_pos_x;
  float         peer_new_pos_y;
  float         peer_new_vel_x;
  float         peer_new_vel_y;
  MPI_Status    status;
  MPI_Request   send_reqs[4];
  MPI_Request   recv_reqs[4];
  int           xfer[4];
  int           wait_recvs;
  int           wait_sends;
  DemoMPIBall   tmp_balls[DEMO_MPI_MAX_BALLS];
  int           i;
  int           j;
  int           k;


  // Do enough intermediate time steps without drawing, so that we advance at
  // most one pixel per repetition. This simplifies collision detection,
  // without the need to check if fast balls "jump" others.
  for (j = 0; j < DEMO_MPI_MAX_ABS_VEL; j++) {

    // Post receives for each of our neighbors who want to send us balls
    wait_recvs = 0;
    if (north_rank > -1) {
      MPI_Irecv(comm_balls[4], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall), 
                MPI_CHAR, north_rank, j * 4 + 0, MPI_COMM_WORLD, 
                &recv_reqs[wait_recvs++]);
    }
    if (east_rank > -1) {
      MPI_Irecv(comm_balls[5], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall), 
                MPI_CHAR, east_rank, j * 4 + 1, MPI_COMM_WORLD, 
                &recv_reqs[wait_recvs++]);
    }
    if (south_rank > -1) {
      MPI_Irecv(comm_balls[6], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall), 
                MPI_CHAR, south_rank, j * 4 + 2, MPI_COMM_WORLD, 
                &recv_reqs[wait_recvs++]);
    }
    if (west_rank > -1) {
      MPI_Irecv(comm_balls[7], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall), 
                MPI_CHAR, west_rank, j * 4 + 3, MPI_COMM_WORLD, 
                &recv_reqs[wait_recvs++]);
    }

    // Advance positions of all my balls, checking for collisions among them
    for (i = 0; i < *my_num_balls; i++) {

      // Compute our next position and velocity, taking into account if
      // we hit any wall. 
      demo_mpi_new_pos(my_balls[i].pos_x, my_balls[i].pos_y, 
                       my_balls[i].vel_x, my_balls[i].vel_y, 
                       &my_new_pos_x, &my_new_pos_y,
                       &my_new_vel_x, &my_new_vel_y);

      // Compute my new velocities if we're colliding with the rod
      if (rod_present) {
        demo_mpi_collision(&my_new_pos_x, &my_new_pos_y, &my_new_vel_x,
                           &my_new_vel_y, DEMO_MPI_BALL_RADIUS, 
                           rod_x0, rod_y0, -my_new_vel_x, 
                           -my_new_vel_y, DEMO_MPI_ROD_RADIUS);
      }

      // For all other balls than the i one
      for (k = 0; k < *my_num_balls; k++) {
        if (k == i) {
          continue;
        }

        // Compute the other ball's next position
        demo_mpi_new_pos(my_balls[k].pos_x, my_balls[k].pos_y, 
                         my_balls[k].vel_x, my_balls[k].vel_y, 
                         &peer_new_pos_x, &peer_new_pos_y,
                         &peer_new_vel_x, &peer_new_vel_y);

        // Compute my new velocities if we're colliding with the other ball
        demo_mpi_collision(&my_new_pos_x, &my_new_pos_y, &my_new_vel_x,
                           &my_new_vel_y, DEMO_MPI_BALL_RADIUS, 
                           peer_new_pos_x, peer_new_pos_y, peer_new_vel_x, 
                           peer_new_vel_y, DEMO_MPI_BALL_RADIUS);
      }

      // Assign new values to tmp array. We do that so all balls interact
      // with each other without moving some of them in between (this could
      // affect collision detection)
      tmp_balls[i].pos_x = my_new_pos_x;
      tmp_balls[i].pos_y = my_new_pos_y;
      tmp_balls[i].vel_x = my_new_vel_x;
      tmp_balls[i].vel_y = my_new_vel_y;
    }

    // Transfer values to the official array
    for (i = 0; i < *my_num_balls; i++) {
      my_balls[i].pos_x = tmp_balls[i].pos_x;
      my_balls[i].pos_y = tmp_balls[i].pos_y;
      my_balls[i].vel_x = tmp_balls[i].vel_x;
      my_balls[i].vel_y = tmp_balls[i].vel_y;
    }

    // Find which of our balls need to change rank
    xfer[0] = 0; // north
    xfer[1] = 0; // east
    xfer[2] = 0; // south
    xfer[3] = 0; // west
    for (i = *my_num_balls - 1; i >= 0; i--) {
      if ((my_balls[i].pos_y < my_min_y) && (north_rank > -1)) {
        comm_balls[0][xfer[0]] = my_balls[i];
        comm_balls[0][xfer[0]++].valid = 1;
        ar_assert(xfer[0] < DEMO_MPI_MAX_XFER_BALLS);
        if (i != *my_num_balls - 1) {
          my_balls[i] = my_balls[*my_num_balls - 1];
        }
        (*my_num_balls)--;
        continue;
      }
      if ((my_balls[i].pos_x > my_max_x) && (east_rank > -1)) {
        comm_balls[1][xfer[1]] = my_balls[i];
        comm_balls[1][xfer[1]++].valid = 1;
        ar_assert(xfer[1] < DEMO_MPI_MAX_XFER_BALLS);
        if (i != *my_num_balls - 1) {
          my_balls[i] = my_balls[*my_num_balls - 1];
        }
        (*my_num_balls)--;
        continue;
      }
      if ((my_balls[i].pos_y > my_max_y) && (south_rank > -1)) {
        comm_balls[2][xfer[2]] = my_balls[i];
        comm_balls[2][xfer[2]++].valid = 1;
        ar_assert(xfer[2] < DEMO_MPI_MAX_XFER_BALLS);
        if (i != *my_num_balls - 1) {
          my_balls[i] = my_balls[*my_num_balls - 1];
        }
        (*my_num_balls)--;
        continue;
      }
      if ((my_balls[i].pos_x < my_min_x) && (west_rank > -1)) {
        comm_balls[3][xfer[3]] = my_balls[i];
        comm_balls[3][xfer[3]++].valid = 1;
        ar_assert(xfer[3] < DEMO_MPI_MAX_XFER_BALLS);
        if (i != *my_num_balls - 1) {
          my_balls[i] = my_balls[*my_num_balls - 1];
        }
        (*my_num_balls)--;
        continue;
      }
    }

    // Send them, marking the next ball after the last of each set "invalid"
    wait_sends = 0;
    if (north_rank > -1) {
      comm_balls[0][xfer[0]].valid = 0;
      MPI_Isend(comm_balls[0], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall),
                MPI_CHAR, north_rank, j * 4 + 2, MPI_COMM_WORLD, 
                &send_reqs[wait_sends++]);
    }
    if (east_rank > -1) {
      comm_balls[1][xfer[1]].valid = 0;
      MPI_Isend(comm_balls[1], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall),
                MPI_CHAR, east_rank, j * 4 + 3, MPI_COMM_WORLD, 
                &send_reqs[wait_sends++]);
    }
    if (south_rank > -1) {
      comm_balls[2][xfer[2]].valid = 0;
      MPI_Isend(comm_balls[2], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall),
                MPI_CHAR, south_rank, j * 4 + 0, MPI_COMM_WORLD, 
                &send_reqs[wait_sends++]);
    }
    if (west_rank > -1) {
      comm_balls[3][xfer[3]].valid = 0;
      MPI_Isend(comm_balls[3], DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall),
                MPI_CHAR, west_rank, j * 4 + 1, MPI_COMM_WORLD, 
                &send_reqs[wait_sends++]);
    }

    // See what our neighbors sent us
    MPI_Waitall(wait_recvs, recv_reqs, &status);

    if (north_rank > -1) {
      for (i = 0; i < DEMO_MPI_MAX_XFER_BALLS; i++) {
        if (!comm_balls[4][i].valid) {
          break;
        }
        ar_assert(*my_num_balls < DEMO_MPI_MAX_BALLS - 1);
        my_balls[(*my_num_balls)++] = comm_balls[4][i];
      }
    }
    if (east_rank > -1) {
      for (i = 0; i < DEMO_MPI_MAX_XFER_BALLS; i++) {
        if (!comm_balls[5][i].valid) {
          break;
        }
        ar_assert(*my_num_balls < DEMO_MPI_MAX_BALLS - 1);
        my_balls[(*my_num_balls)++] = comm_balls[5][i];
      }
    }
    if (south_rank > -1) {
      for (i = 0; i < DEMO_MPI_MAX_XFER_BALLS; i++) {
        if (!comm_balls[6][i].valid) {
          break;
        }
        ar_assert(*my_num_balls < DEMO_MPI_MAX_BALLS - 1);
        my_balls[(*my_num_balls)++] = comm_balls[6][i];
      }
    }
    if (west_rank > -1) {
      for (i = 0; i < DEMO_MPI_MAX_XFER_BALLS; i++) {
        if (!comm_balls[7][i].valid) {
          break;
        }
        ar_assert(*my_num_balls < DEMO_MPI_MAX_BALLS - 1);
        my_balls[(*my_num_balls)++] = comm_balls[7][i];
      }
    }

    // Wait for our posted sends, before we reuse the arrays for the next rep
    MPI_Waitall(wait_sends, send_reqs, &status);
  }
}


// ===========================================================================
// ===========================================================================
void demo_mpi_worker_loop(int rank, int my_bid, int my_cid) {

  volatile DemoMPIOutScene      *scene[DEMO_MPI_OUT_BUFS];
  int                           cur_cmd;
  int                           cur_scene;
  int                           frame;
  int                           i;
  volatile DemoMPICommand       *cmd;
  unsigned int                  cmd_base;
  unsigned int                  scene_base;
  DemoMPIBall                   my_balls[DEMO_MPI_MAX_BALLS];
  DemoMPIBall                   *comm_balls[8];
  int                           my_num_balls;
  int                           my_tile_x;
  int                           my_tile_y;
  int                           north_rank;
  int                           east_rank;
  int                           west_rank;
  int                           south_rank;
  float                         my_min_x;
  float                         my_max_x;
  float                         my_min_y;
  float                         my_max_y;
  unsigned int                  seed;
  int                           idx;


  
  // Allocate incoming command buffers
  cmd_base = (unsigned int) kt_malloc(DEMO_MPI_CMD_BUFS * sizeof(DemoMPICommand));

  // Allocate outgoing scene buffers
  for (i = 0; i < DEMO_MPI_OUT_BUFS; i++) {
    scene[i] = kt_malloc(sizeof(DemoMPIOutScene));
  }

  // Allocate four plus four communication buffers
  for (i = 0; i < 8; i++) {
    comm_balls[i] = kt_malloc(DEMO_MPI_MAX_XFER_BALLS * sizeof(DemoMPIBall));
  }

  // Initialize our command counters
  for (i = 0; i < DEMO_MPI_CMD_BUFS; i++) {
    ar_cnt_set(my_cid, cnt_id(i), 0);
  }

  // Send our command buffer base to the video input core
  ar_mbox_send2(my_cid, AR_ARM0_BID, 0, (my_bid << 8) | my_cid, cmd_base);

  // Learn the scene buffer base of the video output core
  ar_mbox_send(my_cid, AR_ARM0_BID, 1, (my_bid << 8) | my_cid);
  scene_base = ar_mslot_get(my_cid);


  // Loop
  cur_cmd = 0;
  cur_scene = 0;
  frame = 0;
  my_num_balls = 0;
  my_tile_x = -1;
  my_tile_y = -1;
  my_min_x = -1.0f;
  my_max_x = -1.0f;
  my_min_y = -1.0f;
  my_max_y = -1.0f;
  north_rank = -1;
  east_rank = -1;
  south_rank = -1;
  west_rank = -1;
  seed = rank;
  while (1) {

    // Wait for a command from the input
    while (ar_cnt_get(my_cid, cnt_id(cur_cmd)) != sizeof(DemoMPICommand)) {
      ;
    }
    cmd = (DemoMPICommand *) (cmd_base + cur_cmd * sizeof(DemoMPICommand));


    // Are we active for the current frame?
    if (rank < cmd->num_slaves) {

      // Init balls?
      if (cmd->init_balls) {

        if (cmd->init_num_balls < 0) {
          my_num_balls = (rank < -cmd->init_num_balls) ? 1 : 0;
        }
        else {
          my_num_balls = cmd->init_num_balls;
        }

        my_tile_x = rank % cmd->tiles_x;
        my_tile_y = rank / cmd->tiles_x;
        ar_assert(my_tile_y < cmd->tiles_y);

        west_rank  = (my_tile_x > 0)                ? rank - 1 : -1;
        east_rank  = (my_tile_x < cmd->tiles_x - 1) ? rank + 1 : -1;
        north_rank = (my_tile_y > 0)                ? rank - cmd->tiles_x : -1;
        south_rank = (my_tile_y < cmd->tiles_y - 1) ? rank + cmd->tiles_x : -1;

        my_min_x = my_tile_x * cmd->blk_width;
        my_max_x = (float) ((my_tile_x + 1) * cmd->blk_width) - 0.0001f;
        my_min_y = my_tile_y * cmd->blk_height;
        my_max_y = (float) ((my_tile_y + 1) * cmd->blk_height) - 0.0001f;

        for (i = 0; i < my_num_balls; i++) {
          my_balls[i].pos_x = ((seed = kt_rand(seed)) % 
                               (int) (my_max_x - my_min_x)) + my_min_x;
          my_balls[i].pos_y = ((seed = kt_rand(seed)) % 
                               (int) (my_max_y - my_min_y)) + my_min_y;
          my_balls[i].vel_x = ((float) ((seed = kt_rand(seed)) %
                                        (200 * DEMO_MPI_MAX_ABS_VEL)) -
                               (100.0f * DEMO_MPI_MAX_ABS_VEL)) / 100.0f;
          my_balls[i].vel_y = ((float) ((seed = kt_rand(seed)) %
                                        (200 * DEMO_MPI_MAX_ABS_VEL)) -
                               (100.0f * DEMO_MPI_MAX_ABS_VEL)) / 100.0f;
        }
      }

      // Process
      demo_mpi_do_step(rank, north_rank, east_rank, south_rank, west_rank,
                       my_min_x, my_min_y, my_max_x, my_max_y, 
                       my_balls, &my_num_balls, comm_balls, cmd->rod_present,
                       cmd->rod_x0, cmd->rod_y0);

      // Prepare scene to output
      scene[cur_scene]->active_balls = my_num_balls;
      for (i = 0, idx = 0; i < my_num_balls; i++) {
        scene[cur_scene]->balls[i] = (int) my_balls[i].pos_x | 
                                     (((int) my_balls[i].pos_y) << 10);
        switch (i % 3) {
          case 0: 
            scene[cur_scene]->balls[idx] = (int) my_balls[i].pos_x |
                                           (((int) my_balls[i].pos_y) << 10);
            break;
          case 1: 
            scene[cur_scene]->balls[idx++] |= (((int) my_balls[i].pos_x) << 20);
            scene[cur_scene]->balls[idx] = (int) my_balls[i].pos_y;
            break;
          case 2: 
            scene[cur_scene]->balls[idx++] |= (((int) my_balls[i].pos_x) << 10) |
                                              (((int) my_balls[i].pos_y) << 20);
            break;
          default:
            ar_abort();
        }
      }

      // Send to video output core
      while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
        ;
      }
      ar_dma_with_ack(my_cid,
                      my_bid, my_cid, (unsigned int) scene[cur_scene],
                      AR_ARM0_BID, 1, scene_buf(scene_base, cur_scene, rank),
                      AR_ARM0_BID, 1, DEMO_MPI_OUT_CNT_BASE + cur_scene,
                      sizeof(DemoMPIOutScene), 0, 0, 0);
    }

    // Signal to the input core we're ready for a new command
    ar_cnt_set(my_cid, cnt_id(cur_cmd), 0);
    while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
      ;
    }
    ar_cnt_incr(my_cid, AR_ARM0_BID, 0, cnt_id(cur_cmd), 1);

    // Next
    cur_cmd++;
    if (cur_cmd >= DEMO_MPI_CMD_BUFS) {
      cur_cmd = 0;
    }
    cur_scene++;
    if (cur_scene >= DEMO_MPI_OUT_BUFS) {
      cur_scene = 0;
    }
    frame++;
  }
}



// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
void demo_mpi() {

  int                   rank;
  int                   my_bid;
  int                   my_cid;


  // Sanity checks
  ar_assert(sizeof(DemoMPICommand) % 64 == 0);
  ar_assert(sizeof(DemoMPIOutScene) % 64 == 0);

  // Who are we?
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();

  // ARM cores
  if (my_bid == AR_ARM0_BID) {
    if (my_cid == 0) {
      demo_mpi_input_loop(my_bid, my_cid);
    }
    else if (my_cid == 1) {
      demo_mpi_output_loop(my_bid, my_cid);
    }
    else {
      while (1) {
        ;
      }
    }
  }

  // Initialize MPI
  //MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // MPI Workers
  demo_mpi_worker_loop(rank, my_bid, my_cid);
}


