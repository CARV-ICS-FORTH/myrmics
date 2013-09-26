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
// Abstract      : XUP-based video demo using Myrmics. This file contains the
//                 code run by the two ARM cores for video input and output,
//                 as well as any low-level functionality that needs
//                 arch- and/or context-level access, which will confuse 
//                 SCOOP and is separated here.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: demo_myrmics_vid.c,v $
// CVS revision  : $Revision: 1.16 $
// Last modified : $Date: 2013/04/03 07:59:53 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <memory_management.h>
#include <video.h>
#include <noc.h>


// ===========================================================================
// Simple video passthrough (gets a XUP input frame and pushes it to the 
// XUP output). This is a stand-alone routine, not used for demo, but used
// as a test for XUP functionality and Myrmics integration.
// ===========================================================================
void demo_myrmics_passthrough() {

  Context                       *context;
  int                           my_bid;
  int                           my_cid;
  unsigned int                  *buf;
  int                           out_cur;
  int                           in;


  // Who am I?
  context = mm_get_context(ar_get_core_id());
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();

  // Initialize XUP video
  vid_init(my_bid, my_cid);
  vid_in_set_frame_mask(my_bid, my_cid,  0xF0);
  vid_out_set_frame_mask(my_bid, my_cid, 0x01);

  // Initialize video buffer
  buf = kt_malloc(VID_FULL_WIDTH * VID_FULL_HEIGHT * 4);

  // Loop
  out_cur = 0;
  while (1) {

    // Find last stable input frame and hold it.
    in = vid_in_hold_last_frame(my_bid, my_cid, 0xF0);
    ar_assert(in <= 7);
    if (in < 4) {
      // No frame yet
      continue;
    }

    // Copy it from XUP to our working buffer
    ar_receive_xup_frame(my_bid, my_cid, (unsigned int) buf, 
                         VID_DRAM_BASE + (in * VID_BUF_SIZE), 
                         NOC_COUNTER_WAKEUP3);
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP3)) {
      ;
    }

    // Return frame to video input, we've now copied it
    vid_in_set_frame_mask(my_bid, my_cid, 0xF0);

    // Send frame to XUP
    ar_send_xup_frame(my_bid, my_cid, (unsigned int) buf, 
                      VID_DRAM_BASE + (out_cur * VID_BUF_SIZE), 
                      NOC_COUNTER_WAKEUP3);
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP3)) {
      ;
    }

    // Display it
    vid_out_set_frame_mask(my_bid, my_cid, 1 << out_cur);

    // Change vid_out target for next time
    out_cur = (out_cur + 1) % 4;
  }
}


// ===========================================================================
// ===========================================================================
int demo_myrmics_user_menu(int force_print, int *config, int *num_tasks, 
                           int *scene) {

  int   key_hit;
  int   c;
  char  *menu = "\r\n"
                "Current settings: %d tasks, %d workers, %s scheduling, scene %d\r\n"
                "\r\n"
                "Select an option:\r\n"
                "[1] Increase tasks\r\n"
                "[2] Decrease tasks\r\n"
                "[3] Increase worker cores\r\n"
                "[4] Decrease worker cores\r\n"
                "[5] Toggle flat/hierarchical scheduling\r\n"
                "[6] Increase scene complexity\r\n"
                "[7] Decrease scene complexity\r\n";
  char  *err_t = "\r\n** Error: Number of tasks must be between %d and %d **\r";
  char  *err_w = "\r\n** Error: Worker cores must be between 1 and 512 **\r";
  char  *err_s = "\r\n** Error: Scene complexity must be between 1 and 3 **\r";

  if (*config == -1) {
    *config = 4;     // 16 workers, flat
    *num_tasks = 16; // 16 tasks
    *scene = 0;      // simple scene
  }
  if (force_print) {
    kt_printf(menu, *num_tasks, kt_int_pow(2, *config & 0xF),
              (*config & 0x10) ? "hier" : "flat", *scene + 1);
  }

  c = ar_uart_get_char();
  if (c > -1) {
    key_hit = 1;
  }
  else {
    key_hit = 0;
  }


  if ((c >= '1') && (c <= '7')) {
    switch (c) {

      case '1':
        if (*num_tasks < DEMO_MYRMICS_MAX_TASKS) {
          *num_tasks *= 2;
        }
        else {
          kt_printf(err_t, DEMO_MYRMICS_MIN_TASKS, DEMO_MYRMICS_MAX_TASKS);
        }
        break;

      case '2':
        if (*num_tasks > DEMO_MYRMICS_MIN_TASKS) {
          *num_tasks /= 2;
        }
        else {
          kt_printf(err_t, DEMO_MYRMICS_MIN_TASKS, DEMO_MYRMICS_MAX_TASKS);
        }
        break;

      case '3':
        if ((*config & 0xF) < 9) {
          *config = (*config & 0x10) | ((*config & 0xF) + 1);
        }
        else {
          kt_printf(err_w);
        }
        break;

      case '4':
        if ((*config & 0xF) > 0) {
          *config = (*config & 0x10) | ((*config & 0xF) - 1);
        }
        else {
          kt_printf(err_w);
        }
        break;

      case '5':
        *config = (0x10 - (*config & 0x10)) | (*config & 0xF);
        break;

      case '6':
        if (*scene < 2) {
          (*scene)++;
        }
        else {
          kt_printf(err_s);
        }
        break;

      case '7':
        if (*scene > 0) {
          (*scene)--;
        }
        else {
          kt_printf(err_s);
        }
        break;

    }
    
    kt_printf(menu, *num_tasks, kt_int_pow(2, *config & 0xF),
              (*config & 0x10) ? "hier" : "flat", *scene + 1);
  }

  return key_hit;
}



// ===========================================================================
// ===========================================================================
void demo_myrmics_input_loop() {

  Context                       *context;
  int                           my_bid;
  int                           my_cid;
  volatile DemoMyrmicsInScene   *comm_buf;
  unsigned int                  *buf;
  int                           in;
  int                           x;
  int                           y;
  unsigned int                  *p;
  unsigned int                  red;
  unsigned int                  green;
  unsigned int                  blue;
  unsigned int                  peer_buf;
  int                           red_min_x;
  int                           red_min_y;
  int                           red_max_x;
  int                           red_max_y;
  int                           green_min_x;
  int                           green_min_y;
  int                           green_max_x;
  int                           green_max_y;


  // Who am I?
  context = mm_get_context(ar_get_core_id());
  my_bid = context->vid_demo_in_bid;
  my_cid = context->vid_demo_in_cid;

  // Initialize XUP video
  vid_init(my_bid, my_cid);
  vid_in_set_frame_mask(my_bid, my_cid,  0xF0);
  vid_out_set_frame_mask(my_bid, my_cid, 0x01);

  // Initialize video buffer
  buf = kt_malloc(VID_FULL_WIDTH * VID_FULL_HEIGHT * 4);

  // Allocate a communication buffer
  comm_buf = kt_malloc(sizeof(DemoMyrmicsInScene));


  // Loop
  while (1) {

    // Find last stable input frame and hold it.
    in = vid_in_hold_last_frame(my_bid, my_cid, 0xF0);
    ar_assert(in <= 7);
    if (in < 4) {
      // No frame yet
      continue;
    }

    // Copy it from XUP to our working buffer
    ar_receive_xup_frame(my_bid, my_cid, (unsigned int) buf, 
                         VID_DRAM_BASE + (in * VID_BUF_SIZE), 
                         DEMO_MYRMICS_IN_CNT_XUP);
    while (ar_cnt_get(my_cid, DEMO_MYRMICS_IN_CNT_XUP)) {
      ;
    }

    // Return frame to video input, we've now copied it
    vid_in_set_frame_mask(my_bid, my_cid, 0xF0);


    // Do a keying pass to locate a red and a green ball. Do the analysis on
    // pixel by pixel and line by line, so we don't devote too much time.
    red_min_x = 9999;
    red_min_y = 9999;
    red_max_x = -1;
    red_max_y = -1;

    green_min_x = 9999;
    green_min_y = 9999;
    green_max_x = -1;
    green_max_y = -1;

    p = buf + VID_FULL_WIDTH * (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 +
              VID_FULL_WIDTH / 4;

    for (y = (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2;
         y < (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 + VID_STRIP_HEIGHT;
         y += 1, p += 0 * VID_FULL_WIDTH + VID_FULL_WIDTH / 2) {
      for (x = VID_FULL_WIDTH / 4; x < (VID_FULL_WIDTH * 3) / 4; 
           x += 1, p += 1) {

        red   = (*p >> 8 ) & 0xFF;
        green = (*p >> 16) & 0xFF;
        blue  = (*p >> 24) & 0xFF;

        if ((red   >= DEMO_MYRMICS_REDKEY_RED_MIN) && 
            (green <= DEMO_MYRMICS_REDKEY_GREEN_MAX) && 
            (blue  <= DEMO_MYRMICS_REDKEY_BLUE_MAX)) {
          if (x < red_min_x) {
            red_min_x = x;
          }
          if (y < red_min_y) {
            red_min_y = y;
          }
          if (x > red_max_x) {
            red_max_x = x;
          }
          if (y > red_max_y) {
            red_max_y = y;
          }
        }

        if ((red   <= DEMO_MYRMICS_GREENKEY_RED_MAX) && 
            (green >= DEMO_MYRMICS_GREENKEY_GREEN_MIN) && 
            (blue  <= DEMO_MYRMICS_GREENKEY_BLUE_MAX)) {
          if (x < green_min_x) {
            green_min_x = x;
          }
          if (y < green_min_y) {
            green_min_y = y;
          }
          if (x > green_max_x) {
            green_max_x = x;
          }
          if (y > green_max_y) {
            green_max_y = y;
          }
        }
      }
    }

    // Anything in mailbox?
    if (!(ar_mbox_status_get(my_cid) & 0xFFFF)) {
      continue;
    }


    // Dequeue video output core buffer
    peer_buf = ar_mbox_get(my_cid);

    // Fill the communication buffer. The camera mirrors the image, so to
    // match the rod as the user would expect it, we also mirror it. 
    // Additionally, we translate the coordinates from the full frame system
    // to the in-border drawing area.
    if (red_min_x < 9999) {

      x = ((red_max_x + red_min_x) / 2 - VID_FULL_WIDTH / 4) * 2;
      y = (red_max_y + red_min_y) / 2;
      
      // Mirror
      x = VID_FULL_WIDTH - x;

      // Translate and clamp
      y -= (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2;
      if (x < DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        x = DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }
      if (y < DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        y = DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }
      if (x >= VID_FULL_WIDTH - DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        x = VID_FULL_WIDTH - 1 - DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }
      if (y >= VID_STRIP_HEIGHT - DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        y = VID_STRIP_HEIGHT - 1 - DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }

      // Fill buffer
      comm_buf->red_present = 1;
      comm_buf->red_x = x;
      comm_buf->red_y = y;
    }
    else {
      comm_buf->red_present = 0;
    }

    if (green_min_x < 9999) {

      x = ((green_max_x + green_min_x) / 2 - VID_FULL_WIDTH / 4) * 2;
      y = (green_max_y + green_min_y) / 2;
      
      // Mirror
      x = VID_FULL_WIDTH - x;

      // Translate and clamp
      y -= (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2;
      if (x < DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        x = DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }
      if (y < DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        y = DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }
      if (x >= VID_FULL_WIDTH - DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        x = VID_FULL_WIDTH - 1 - DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }
      if (y >= VID_STRIP_HEIGHT - DEMO_MYRMICS_BALL_SKETCH_RADIUS) {
        y = VID_STRIP_HEIGHT - 1 - DEMO_MYRMICS_BALL_SKETCH_RADIUS;
      }

      // Fill buffer
      comm_buf->green_present = 1;
      comm_buf->green_x = x;
      comm_buf->green_y = y;
    }
    else {
      comm_buf->green_present = 0;
    }

    // Transfer it
    ar_dma_no_ack(my_cid,
                  my_bid, my_cid, (unsigned int) comm_buf,
                  context->vid_demo_out_bid, context->vid_demo_out_cid, peer_buf,
                  sizeof(DemoMyrmicsInScene), 0, 0, 0);

    // Send video output a mail to notify that the transfer is complete. This
    // is safe because it's a DMA going after the previous DMA, and we're 
    // talking about ARM cores (ARM cache is disabled).
    ar_mslot_send(my_cid, context->vid_demo_out_bid, context->vid_demo_out_cid,
                  0xDEADBEEF);
  }
}


// ===========================================================================
// ===========================================================================
void demo_myrmics_output_loop() {

  Context                       *context;
  int                           my_bid;
  int                           my_cid;
  unsigned int                  *buf;
  unsigned int                  worker_base;
  unsigned int                  word;
  unsigned int                  peer_bid;
  unsigned int                  peer_cid;
  int                           tile_width;
  int                           tile_height;
  int                           tile_id;
  int                           num_tiles;
  unsigned int                  *src;
  unsigned int                  *dst;
  int                           origin_x;
  int                           origin_y;
  int                           tile_size;
  int                           tiles_done;
  unsigned int                  start;
  unsigned int                  elapsed;
  unsigned int                  prev_xup;
  char                          msg[128];
  volatile DemoMyrmicsInScene   *in_scene[2];
  int                           x;
  int                           y;
  int                           z;
  int                           cur_in_scene;
  unsigned int                  *p;
  unsigned int                  *q;
  int                           i;
  int                           core_config;
  int                           num_tasks;
  int                           which_scene;
  int                           w;
  int                           s;
  unsigned int                  *ball_buf[2];
  unsigned int                  forth_colors[] = {0x20A02000, 0x10801000,
                                                  0x00400000, 0x00000000};
  unsigned int                  myrmics_colors[] = {0x10801000, 0x08500800,
                                                    0x00300000, 0x00000000};



  // Who am I?
  context = mm_get_context(ar_get_core_id());
  my_bid = context->vid_demo_out_bid;
  my_cid = context->vid_demo_out_cid;

  // Allocate clean video buffer
  buf = kt_zalloc(VID_FULL_WIDTH * VID_FULL_HEIGHT * 4);

  // Allocate a place for the workers to send their rendered tiles. We assume
  // 1024 max tiles, each rendering an 25*14 = 350 pixels tile, moved to 
  // 384 pixels to accommodate 64B-line padding.
  worker_base = (unsigned int) kt_malloc(384 * 4 * 1024);

  // Allocate two clean buffers to communicate with the video input core
  for (i = 0; i < 2; i++) {
    in_scene[i] = kt_zalloc(sizeof(DemoMyrmicsInScene));
  }

  // Put static stamps on the video buffer
  vid_blit_logo(buf, VID_FULL_WIDTH, 0, 524, 1, forth_colors);

  vid_font_putstr(buf, VID_FULL_WIDTH, 82, 535, 0, 0, 0x20A02000, 0, 0, 
                  "Myrmics Runtime System Demonstration");
  vid_font_putstr(buf, VID_FULL_WIDTH, 82, 547, 0, 0, 0x20A02000, 0, 0, 
                  "FORTH-ICS/CARV (c) 2012-2013");
    
  vid_blit_logo(buf, VID_FULL_WIDTH, 692, 524, 0, myrmics_colors);

  // Always display the same XUP frame
  vid_out_set_frame_mask(my_bid, my_cid, 1 << 0);

  // Show the empty XUP frame
  ar_send_xup_frame(my_bid, my_cid, (unsigned int) buf, 
                    VID_DRAM_BASE + (0 * VID_BUF_SIZE), 
                    DEMO_MYRMICS_OUT_CNT_XUP);
  while (ar_cnt_get(my_cid, DEMO_MYRMICS_OUT_CNT_XUP)) {
    ;
  }

  // Prepare the red and green ball sprites for the preparatory loop
  for (i = 0; i < 2; i++) {
    ball_buf[i] = kt_zalloc((2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1) * 
                            (2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1) * 
                            sizeof(int));

    vid_draw_circle(ball_buf[i], 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1, 
                    2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1, 
                    DEMO_MYRMICS_BALL_SKETCH_RADIUS, 
                    DEMO_MYRMICS_BALL_SKETCH_RADIUS, 
                    DEMO_MYRMICS_BALL_SKETCH_RADIUS, 
                    (i == 0) ? DEMO_MYRMICS_RED_SKETCH_COLOR :
                               DEMO_MYRMICS_GREEN_SKETCH_COLOR);

    vid_draw_floodfill(ball_buf[i], 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1, 
                       2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1, 
                       DEMO_MYRMICS_BALL_SKETCH_RADIUS, 
                       DEMO_MYRMICS_BALL_SKETCH_RADIUS, 
                       (i == 0) ? DEMO_MYRMICS_RED_SKETCH_FILL :
                                  DEMO_MYRMICS_GREEN_SKETCH_FILL, 
                       (i == 0) ? DEMO_MYRMICS_RED_SKETCH_COLOR :
                                  DEMO_MYRMICS_GREEN_SKETCH_COLOR);
  }

  // Get default values
  core_config = -1;
  demo_myrmics_user_menu(0, &core_config, &num_tasks, &which_scene);

start:

  // Print the menu
  demo_myrmics_user_menu(1, &core_config, &num_tasks, &which_scene);

  // Preparatory loop
  cur_in_scene = 0;
  ar_timer_reset();
  start = ar_free_timer_get_ticks();
  while (1) {

    // User input?
    if (demo_myrmics_user_menu(0, &core_config, &num_tasks, &which_scene)) {
      // Restart the timeout when user types anything
      start = ar_free_timer_get_ticks();
    }



    // Is it time to draw a new one?
    if (ar_timer_get_msec() < DEMO_MYRMICS_OUT_UPDATE_MSEC) {
      continue;
    }
    ar_timer_reset();

    // Get a new input scene
    ar_mbox_send(my_cid, context->vid_demo_in_bid, context->vid_demo_in_cid,
                 (unsigned int) in_scene[cur_in_scene]);

    // Wait until it arrives
    ar_assert(ar_mslot_get(my_cid) == 0xDEADBEEF);

    // Clear previous ball positions
    for (i = 0; i < 2; i++) {
      if (i == 0) {
        if (!in_scene[1 - cur_in_scene]->red_present) {
          continue;
        }
        x = in_scene[1 - cur_in_scene]->red_x;
        y = in_scene[1 - cur_in_scene]->red_y;
      }
      else {
        if (!in_scene[1 - cur_in_scene]->green_present) {
          continue;
        }
        x = in_scene[1 - cur_in_scene]->green_x;
        y = in_scene[1 - cur_in_scene]->green_y;
      }

      p = buf + (y - DEMO_MYRMICS_BALL_SKETCH_RADIUS) * VID_FULL_WIDTH +
                (x - DEMO_MYRMICS_BALL_SKETCH_RADIUS);
      for (y = 0; 
           y < 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1; 
           y++, p += VID_FULL_WIDTH - 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS - 1) {
        for (x = 0; 
             x < 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1; 
             x++, p++) {
          *p = 0;
        }
      }
    }

    // Fire Myrmics here, before scene coordinates translation
    elapsed = (ar_free_timer_get_ticks() - start) / 10000000;
    if (DEMO_MYRMICS_START_TIMEOUT - elapsed <= 0) {
      break;
    }


    // Translate new positions and blit red and green balls
    for (i = 0; i < 2; i++) {
      if (i == 0) {
        if (!in_scene[cur_in_scene]->red_present) {
          continue;
        }
        in_scene[cur_in_scene]->red_y += (VID_FULL_HEIGHT-VID_STRIP_HEIGHT) / 2;
        x = in_scene[cur_in_scene]->red_x;
        y = in_scene[cur_in_scene]->red_y;
      }
      else {
        if (!in_scene[cur_in_scene]->green_present) {
          continue;
        }
        in_scene[cur_in_scene]->green_y += (VID_FULL_HEIGHT-VID_STRIP_HEIGHT) / 2;
        x = in_scene[cur_in_scene]->green_x;
        y = in_scene[cur_in_scene]->green_y;
      }

      p = buf + (y - DEMO_MYRMICS_BALL_SKETCH_RADIUS) * VID_FULL_WIDTH +
                (x - DEMO_MYRMICS_BALL_SKETCH_RADIUS);
      q = ball_buf[i];
      for (y = 0; 
           y < 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1; 
           y++, p += VID_FULL_WIDTH - 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS - 1) {
        for (x = 0; 
             x < 2 * DEMO_MYRMICS_BALL_SKETCH_RADIUS + 1; 
             x++, p++, q++) {
          if (*q) {
            *p = *q;
          }
        }
      }
    }

    // Print messages
    kt_sprintf(msg, "Select ball positions, change parameters");
    vid_font_putstr(buf, VID_FULL_WIDTH, 10, 10, 1, 0, 0x00FFFF00, 0, 0, msg);

    kt_sprintf(msg, "Starting render in %d...                   ", 
               DEMO_MYRMICS_START_TIMEOUT - elapsed);
//kt_printf("%s\r\n", msg);
    vid_font_putstr(buf, VID_FULL_WIDTH, 10, 42, 1, 0, 0x00FFFF00, 0, 0, msg);

    // Update XUP frame
    ar_send_xup_frame(my_bid, my_cid, (unsigned int) buf, 
                      VID_DRAM_BASE + (0 * VID_BUF_SIZE), 
                      DEMO_MYRMICS_OUT_CNT_XUP);
    while (ar_cnt_get(my_cid, DEMO_MYRMICS_OUT_CNT_XUP)) {
      ;
    }

    // New scene
    cur_in_scene = 1 - cur_in_scene;
  }


  // -------------------------------------------------------------------------

  // Inform all cores (except video in/out) for selected configuration
  for (x = AR_FORMIC_MIN_X; x <= AR_FORMIC_MAX_X; x++) {
    for (y = AR_FORMIC_MIN_Y; y <= AR_FORMIC_MAX_Y; y++) {
      for (z = AR_FORMIC_MIN_Z; z <= AR_FORMIC_MAX_Z; z++) {
        for (i = 0; i < AR_FORMIC_CORES_PER_BOARD; i++) {
          ar_mslot_send(my_cid, (x << 4) | (y << 2) | z, i, core_config);
        }
      }
    }
  }
  if (AR_ARM0_BID != -1) {
    for (i = 0; i < AR_ARM0_CORES_PER_BOARD; i++) {
      if (((AR_ARM0_BID != context->vid_demo_in_bid) || 
           (i != context->vid_demo_in_cid)) &&
          ((AR_ARM0_BID != my_bid) || (i != my_cid))) {
        ar_mslot_send(my_cid, AR_ARM0_BID, i, core_config);
      }
    }
  }
  if (AR_ARM1_BID != -1) {
    for (i = 0; i < AR_ARM1_CORES_PER_BOARD; i++) {
      if (((AR_ARM1_BID != context->vid_demo_in_bid) || 
           (i != context->vid_demo_in_cid)) &&
          ((AR_ARM1_BID != my_bid) || (i != my_cid))) {
        ar_mslot_send(my_cid, AR_ARM1_BID, i, core_config);
      }
    }
  }

  // -------------------------------------------------------------------------

  // Stamp top marquee
  w = kt_int_pow(2, core_config & 0xF);
  if (core_config & 0x10) {
    s = ((core_config & 0xF) <= 4) ? 1 :
        ((core_config & 0xF) == 5) ? 2 : 4;
  }
  else {
    s = 1;
  }
  kt_sprintf(msg, "%d worker core%s, %d scheduler core%s          ", 
             w, (w > 1) ? "s" : "",
             s, (s > 1) ? "s" : "");
  vid_font_putstr(buf, VID_FULL_WIDTH, 10, 10, 1, 0, 0x00FFFF00, 0, 0, msg);


  // Working loop
  start = ar_free_timer_get_ticks();
  tiles_done = 0;
  num_tiles = -1;
  prev_xup = 0;
  while (1) {

    // Any mail?
    if ((ar_mbox_status_get(my_cid) & 0xFFFF) > 0) {

      word = ar_mbox_get(my_cid);
//kt_printf("vo: got mail 0x%08X\r\n", word);
      tile_id   = (word >> 12) & 0xFFF;
      num_tiles = (word >> 24) & 0xFF;

      // Request to send our buffer base, ball positions or selected scene?
      if (tile_id == 0xFFF) {
        peer_bid = (word >> 0) & 0xFF;
        peer_cid = (word >> 8) & 0xF;
        if (num_tiles == 0xFF) {
          ar_mslot_send(my_cid, peer_bid, peer_cid, worker_base);
          continue;
        }
        else if (num_tiles == 0xFE) {
          ar_mslot_send(my_cid, peer_bid, peer_cid, 
                        in_scene[cur_in_scene]->red_present |
                        (in_scene[cur_in_scene]->red_x << 4) | 
                        (in_scene[cur_in_scene]->red_y << 16));
          continue;
        }
        else if (num_tiles == 0xFD) {
          ar_mslot_send(my_cid, peer_bid, peer_cid, 
                        in_scene[cur_in_scene]->green_present |
                        (in_scene[cur_in_scene]->green_x << 4) | 
                        (in_scene[cur_in_scene]->green_y << 16));
          continue;
        }
        else if (num_tiles == 0xFC) {
          ar_mslot_send(my_cid, peer_bid, peer_cid, 
                        which_scene | (num_tasks << 8));
          continue;
        }
        else {
          ar_abort();
        }
      }

      // Got a new tile. Find out its origin and size.
      num_tiles = kt_int_pow(2, num_tiles);
      demo_myrmics_pix_buf_coords(tile_id, num_tiles,
                                  &origin_x, &origin_y, 
                                  &tile_width, &tile_height);

      // Copy it from the worker buffer base to the video buffer
      tile_size = ((tile_width * tile_height * 4 + 63) / 64) * 64;
      src = (unsigned int *) (worker_base + tile_id * tile_size);
      dst = buf + VID_FULL_WIDTH * ((VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2 +
                                    origin_y) + origin_x;
  //kt_printf("vo: tile %d of %d, base 0x%08X, size %dx%d, org %d,%d\r\n", tile_id, num_tiles, src, tile_width, tile_height, origin_x, origin_y);

      for (y = 0; y < tile_height; y++, dst += VID_FULL_WIDTH - tile_width) {
        for (x = 0; x < tile_width; x++) {
          *dst++ = *src++;
        }
      }
      tiles_done++;
    }
    
    // Compute time and proceed for XUP update only when above a threshold
    elapsed = (ar_free_timer_get_ticks() - start) / 10000;
    if (elapsed - prev_xup > DEMO_MYRMICS_OUT_UPDATE_MSEC) {
      prev_xup = elapsed;
    }
    else {
      continue;
    }

    // Stamp buffer
    kt_sprintf(msg, "Time: %d.%d sec, tiles rendered: %d      ", 
               elapsed / 1000, (elapsed % 1000) / 100, tiles_done);
    vid_font_putstr(buf, VID_FULL_WIDTH, 10, 42, 1, 0, 0x00FFFF00, 0, 0, msg);
//kt_printf("%s\r\n", msg);


    // Update XUP frame
    ar_send_xup_frame(my_bid, my_cid, (unsigned int) buf, 
                      VID_DRAM_BASE + (0 * VID_BUF_SIZE), 
                      DEMO_MYRMICS_OUT_CNT_XUP);
    while (ar_cnt_get(my_cid, DEMO_MYRMICS_OUT_CNT_XUP)) {
      ;
    }

    // Finished?
    if (tiles_done == num_tiles) {
      break;
    }
  }

  // Wait a few seconds
  ar_timer_busy_wait_msec(1000 * DEMO_MYRMICS_END_TIMEOUT);

  // Clear the buffer
  p = buf + VID_FULL_WIDTH * (VID_FULL_HEIGHT - VID_STRIP_HEIGHT) / 2;
  for (y = 0; y < VID_STRIP_HEIGHT; y++) {
    for (x = 0; x < VID_FULL_WIDTH; x++) {
      *p++ = 0;
    }
  }

  // Restart
  goto start;
}

// ===========================================================================
// ===========================================================================
void demo_myrmics_init(unsigned int *worker_base, int *red_present, int *red_x,
                       int *red_y, int *green_present, int *green_x, 
                       int *green_y, int *scene, int *num_tasks) {
  
  Context               *context;
  int                   my_bid;
  int                   my_cid;
  unsigned int          word;


  // Get context
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Send a mail to video output core to request the base
  ar_mbox_send(my_cid, context->vid_demo_out_bid, context->vid_demo_out_cid,
               (my_bid >> 0) | (my_cid << 8) | (0xFFF << 12) | (0xFF << 24));

  // Wait until we get an answer
  *worker_base = ar_mslot_get(my_cid);

  // Using the same way, learn about the red ball position...
  ar_mbox_send(my_cid, context->vid_demo_out_bid, context->vid_demo_out_cid,
               (my_bid >> 0) | (my_cid << 8) | (0xFFF << 12) | (0xFE << 24));
  word = ar_mslot_get(my_cid);
  *red_present = word & 0x1;
  if (*red_present) {
    *red_x = (word >> 4) & 0x3FF;
    *red_y = (word >> 16) & 0x3FF;
  }

  // ... the green ball position...
  ar_mbox_send(my_cid, context->vid_demo_out_bid, context->vid_demo_out_cid,
               (my_bid >> 0) | (my_cid << 8) | (0xFFF << 12) | (0xFD << 24));
  word = ar_mslot_get(my_cid);
  *green_present = word & 0x1;
  if (*green_present) {
    *green_x = (word >> 4) & 0x3FF;
    *green_y = (word >> 16) & 0x3FF;
  }

  // ... and finally the selected scene and number of tasks.
  ar_mbox_send(my_cid, context->vid_demo_out_bid, context->vid_demo_out_cid,
               (my_bid >> 0) | (my_cid << 8) | (0xFFF << 12) | (0xFC << 24));
  word = ar_mslot_get(my_cid);
  *scene = word & 0xFF;
  *num_tasks = (word >> 8) & 0xFFFF;
}


// ===========================================================================
// ===========================================================================
void demo_myrmics_pix_buf_coords(int tile_id, int num_tiles,
                                 int *ret_origin_x, int *ret_origin_y, 
                                 int *ret_width, int *ret_height) {

  int tiles_x;
  int tiles_y;
  
  switch (num_tiles) {
    case 1:    tiles_x = 1;  tiles_y = 1;  break;
    case 2:    tiles_x = 2;  tiles_y = 1;  break;
    case 4:    tiles_x = 2;  tiles_y = 2;  break;
    case 8:    tiles_x = 4;  tiles_y = 2;  break;
    case 16:   tiles_x = 4;  tiles_y = 4;  break;
    case 32:   tiles_x = 8;  tiles_y = 4;  break;
    case 64:   tiles_x = 8;  tiles_y = 8;  break;
    case 128:  tiles_x = 16; tiles_y = 8;  break;
    case 256:  tiles_x = 16; tiles_y = 16; break;
    case 512:  tiles_x = 32; tiles_y = 16; break;
    case 1024: tiles_x = 32; tiles_y = 32; break;
    default: ar_abort();
  }

  *ret_width = VID_FULL_WIDTH / tiles_x;
  *ret_height = VID_STRIP_HEIGHT / tiles_y;

  if (ret_origin_x) {
    *ret_origin_x = (tile_id % tiles_x) * *ret_width;
  }
  if (ret_origin_y) {
    *ret_origin_y = (tile_id / tiles_x) * *ret_height;
  }
}


// ===========================================================================
// ===========================================================================
#if 0
void demo_myrmics_get_new_rod(DemoMyrmicsState *state) {
  Context               *context;
  int                   my_bid;
  int                   my_cid;
  ListNode              *cnt_node;
  int                   cnt;


  // Get context
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Get a hardware counter temporarily
  cnt_node = kt_list_head(context->noc_cnt_free);
  ar_assert(cnt_node);
  cnt = (int) cnt_node->data;

  // Initialize it to wait for a DemoMyrmicsInScene transfer
  ar_cnt_set(my_cid, cnt, -sizeof(DemoMyrmicsInScene));

  // Wait until DMA engine can support 1 more DMA
  while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
    ;
  }

  // Send to video input core a 2-word message, indicating who we are,
  // what's our hardware counter and what's our comm buffer address
  ar_mbox_send2(my_cid, context->vid_demo_in_bid, context->vid_demo_in_cid,
                (my_bid << 24) | (my_cid << 16) | (cnt << 0),
                (unsigned int) state->in_com_buf);

  // Wait until counter says the transfer is done
  while (ar_cnt_get(my_cid, cnt)) {
    ;
  }

  // Transfer useful fields
  state->rod_present = state->in_com_buf->rod_present;
  state->rod_x0 = state->in_com_buf->rod_x0;
  state->rod_y0 = state->in_com_buf->rod_y0;
}
#endif


// ===========================================================================
// ===========================================================================
void demo_myrmics_draw_output(unsigned int *pix_buf, int width, int height,
                              int origin_x, int origin_y, 
                              unsigned int vid_out_buf_base, int tile_id,
                              int num_tiles) {

  Context               *context;
  int                   my_bid;
  int                   my_cid;
  ListNode              *cnt_node;
  int                   my_cnt;
  int                   tile_size;
  unsigned int          dst;


  // Get context
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);


  // Compute video output core destination address
  tile_size = ((width * height * 4 + 63) / 64) * 64;
  dst = vid_out_buf_base + tile_id * tile_size;

  // Get a hardware counter temporarily
  cnt_node = kt_list_head(context->noc_cnt_free);
  ar_assert(cnt_node);
  my_cnt = (int) cnt_node->data;

  // Wait until DMA engine can support 2 more DMAs
  while ((ar_ni_status_get(my_cid) & 0xFF) < 2) {
    ;
  }

  // Transfer our buffer to the video output core
  ar_cnt_set(my_cid, my_cnt, -tile_size);
  ar_dma_with_ack(my_cid,
                  my_bid, my_cid, (unsigned int) pix_buf,
                  context->vid_demo_out_bid, context->vid_demo_in_cid, dst,
                  my_bid, my_cid, my_cnt,
                  tile_size, 0, 0, 0);
  
  // Wait for transfer to finish
  while (ar_cnt_get(my_cid, my_cnt)) {
    ;
  }

//kt_printf("%d: transfer to 0x%08X done\r\n", context->pr_core_id, dst);

  // Send mail to notify
  ar_mbox_send(my_cid, context->vid_demo_out_bid, context->vid_demo_out_cid,
               (tile_id << 12) | (kt_int_log2(num_tiles) << 24));

}
