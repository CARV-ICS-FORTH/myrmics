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
// Abstract      : Multilevel kernel benchmark: does RGB->YUV conversion on
//                 data split into 0, 1-, 2- or 3-level deep regions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: multilevel_myrmics.c,v $
// CVS revision  : $Revision: 1.3 $
// Last modified : $Date: 2013/02/20 08:45:17 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>

// ===========================================================================
// ===========================================================================
// void multilevel_myrmics_wait_region(rid_t r) {
//   // Do nothing
// }


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_do_level0(unsigned char *buf, int buf_size, int delay) {

  float r;
  float g;
  float b;
  float y;
  float u;
  float v;
  int   i;


  //for (i = 0; i < buf_size; i += 3) {

  //  // Read pixel buffer
  //  r = (float) buf[i];
  //  g = (float) buf[i+1];
  //  b = (float) buf[i+2];

  //  // Apply contrast filter
  //  r = (r - 128.0F) / 0.63F + 64.0F;
  //  g = (g - 128.0F) / 0.63F + 64.0F;
  //  b = (b - 128.0F) / 0.63F + 64.0F;

  //  // Convert to YUV
  //  y =  0.256788F * r + 0.504129F * g + 0.097906F * b + 16.0F;
  //  u = -0.148223F * r - 0.290993F * g + 0.439216F * b + 128.0F;
  //  v =  0.439216F * r - 0.367788F * g - 0.071427F * b + 128.0F;

  //  // Clamp values
  //  if (y < 0.0F) y = 0.0F;
  //  if (y > 255.0F) y = 255.0F;
  //  if (u < 0.0F) y = 0.0F;
  //  if (u > 255.0F) y = 255.0F;
  //  if (v < 0.0F) y = 0.0F;
  //  if (v > 255.0F) y = 255.0F;

  //  // Write pixel buffer
  //  buf[i]   = (unsigned char) y;
  //  buf[i+1] = (unsigned char) u;
  //  buf[i+2] = (unsigned char) v;
  //}

  sys_timer_busy_wait_cycles(delay);
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_do_level1(rid_t r, unsigned char **buf, int l0_size, 
                                  int buf_size, int delay) {
  unsigned char *scoop_buf;
  int           i;


  for (i = 0; i < l0_size; i++) {

    scoop_buf = buf[i];

    #pragma myrmics task inout(scoop_buf) in(buf_size, delay)
    multilevel_myrmics_do_level0(scoop_buf, buf_size, delay);
  }
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_do_level2(rid_t r, rid_t *regions, 
                                  unsigned char ***buf, int l1_size, 
                                  int l0_size, int buf_size, int delay) {

  rid_t         scoop_region;
  unsigned char **scoop_buf;
  int           i;


  // For all level-1 regions
  for (i = 0; i < l1_size; i++) {

    scoop_region = regions[i];
    scoop_buf    = buf[i];

    #pragma myrmics task region inout(scoop_region) \
                         in(scoop_buf) safe(scoop_buf) \
                         in(l0_size, buf_size, delay)
    multilevel_myrmics_do_level1(scoop_region, scoop_buf, l0_size, buf_size,
                                 delay);
  }
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_do_level3(rid_t r, rid_t *l2_regions, 
                                  rid_t **l1_regions, unsigned char ****buf, 
                                  int l2_size, int l1_size, int l0_size, 
                                  int buf_size, int delay) {

  rid_t         scoop_l2_region;
  rid_t         *scoop_l1_region;
  unsigned char ***scoop_buf;
  int           i;


  // For all level-2 regions
  for (i = 0; i < l2_size; i++) {

    scoop_l2_region = l2_regions[i];
    scoop_l1_region = l1_regions[i];
    scoop_buf       = buf[i];

    #pragma myrmics task region inout(scoop_l2_region) \
                         in(scoop_l1_region, scoop_buf) \
                         safe(scoop_l1_region, scoop_buf) \
                         in(l1_size, l0_size, buf_size, delay)
    multilevel_myrmics_do_level2(scoop_l2_region, scoop_l1_region, scoop_buf, 
                                 l1_size, l0_size, buf_size, delay);
  }
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_do_level4(rid_t r, rid_t *l3_regions,
                                  rid_t **l2_regions, rid_t ***l1_regions, 
                                  unsigned char *****buf, 
                                  int l3_size, int l2_size, int l1_size, 
                                  int l0_size, int buf_size, int delay) {

  rid_t         scoop_l3_region;
  rid_t         *scoop_l2_region;
  rid_t         **scoop_l1_region;
  unsigned char ****scoop_buf;
  int           i;


  // For all level-3 regions
  for (i = 0; i < l3_size; i++) {

    scoop_l3_region = l3_regions[i];
    scoop_l2_region = l2_regions[i];
    scoop_l1_region = l1_regions[i];
    scoop_buf       = buf[i];

    #pragma myrmics task region inout(scoop_l3_region) \
                         in(scoop_l2_region, scoop_l1_region, scoop_buf) \
                         safe(scoop_l2_region, scoop_l1_region, scoop_buf) \
                         in(l2_size, l1_size, l0_size, buf_size, delay)
    multilevel_myrmics_do_level3(scoop_l3_region, scoop_l2_region, 
                                 scoop_l1_region, scoop_buf, 
                                 l2_size, l1_size, l0_size, buf_size, delay);
  }
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_time(rid_t unused, unsigned int time_start,
                             int l3_size, int l2_size, int l1_size, 
                             int l0_size) {

  unsigned int          time_stop;
  unsigned int          time;

  time_stop = sys_free_timer_get_ticks();
  if (time_stop > time_start) {
    time = time_stop - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_stop);
  }

  printf(
      "%3d wrk, %3d sch: %10u cycles (%6u msec) [%d x %d x %d x %d tasks]\r\n", 
      sys_get_num_workers(), sys_get_num_schedulers(),
      time, time / 10000,
      l3_size, l2_size, l1_size, l0_size);
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_l0_main(int l0_size, int buf_size, int num_reps, 
                                int delay) {

  rid_t                 r;
  unsigned char         **buf;
  unsigned int          time_start;
  int                   i;
  int                   j;


  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Create arrays
  buf = sys_alloc(l0_size * sizeof(unsigned char *), r);

  if (l0_size > 1) {
    sys_balloc(buf_size * sizeof(unsigned char), r, l0_size, buf);
  }
  else {
    buf[0] = sys_alloc(buf_size * sizeof(unsigned char), r);
  }


  // Print we're starting
  printf("Tasks of %d cycles starting in %d L0 tile(s)\r\n", 
         delay, l0_size);


  // Start time
  time_start = sys_free_timer_get_ticks();

  // For all reps
  for (i = 0; i < num_reps; i++) {

    // Use one-level-up function to compute this wave
    #pragma myrmics task region inout(r) in(buf) safe(buf) \
                         in(l0_size, buf_size, delay)
    multilevel_myrmics_do_level1(r, buf, l0_size, buf_size, delay);

    // // Emulate transfer to framebuffer
    // #pragma myrmics task region in(r)
    // multilevel_myrmics_wait_region(r);
  }

  // Stop time
  int l3_size = 0;
  int l2_size = 0;
  int l1_size = 0;
  #pragma myrmics task region in(r) \
                       in(time_start, l3_size, l2_size, l1_size, l0_size)
  multilevel_myrmics_time(r, time_start, l3_size, l2_size, l1_size, l0_size);
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_l1_main(int l1_size, int l0_size, int buf_size,
                                int num_reps, int delay) {

  rid_t                 r;
  rid_t                 *regions;
  unsigned char         ***buf;
  unsigned int          time_start;
  int                   i;
  int                   j;
  int                   k;


  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Create regions
  regions = sys_alloc(l1_size * sizeof(rid_t), r);
  for (i = 0; i < l1_size; i++) {
    regions[i] = sys_ralloc(r, 0); // lowest level
  }

  // Create arrays
  buf = sys_alloc(l1_size * sizeof(unsigned char **), r);
  for (i = 0; i < l1_size; i++) {
    buf[i] = sys_alloc(l0_size * sizeof(unsigned char *), regions[i]);

    if (l0_size > 1) {
      sys_balloc(buf_size * sizeof(float), regions[i], l0_size, buf[i]);
    }
    else {
      buf[i][0] = sys_alloc(buf_size * sizeof(float), regions[i]);
    }
  }


  // Print we're starting
  printf("Tasks of %d cycles starting in %d L1 x %d L0 tile(s)\r\n", 
         delay, l1_size, l0_size);

  // Start time
  time_start = sys_free_timer_get_ticks();

  // For all reps
  for (i = 0; i < num_reps; i++) {

    // Use one-level-up function to compute this wave
    #pragma myrmics task region inout(r) in(regions, buf) safe(regions, buf) \
                         in(l1_size, l0_size, buf_size, delay)
    multilevel_myrmics_do_level2(r, regions, buf, l1_size, l0_size, buf_size,
                                 delay);

    // // Emulate transfer to framebuffer
    // #pragma myrmics task region in(r)
    // multilevel_myrmics_wait_region(r);
  }


  // Stop time
  int l3_size = 0;
  int l2_size = 0;
  #pragma myrmics task region in(r) \
                       in(time_start, l3_size, l2_size, l1_size, l0_size)
  multilevel_myrmics_time(r, time_start, l3_size, l2_size, l1_size, l0_size);

}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_l2_main(int l2_size, int l1_size, int l0_size, 
                                int buf_size, int num_reps, int delay) {

  rid_t                 r;
  rid_t                 *l2_regions;
  rid_t                 **l1_regions;
  unsigned char         ****buf;
  unsigned int          time_start;
  int                   i;
  int                   j;
  int                   k;


  // Create all-holding region
  r = sys_ralloc(0, 99); // top level

  // Create regions
  l2_regions = sys_alloc(l2_size * sizeof(rid_t), r);
  l1_regions = sys_alloc(l2_size * sizeof(rid_t *), r);
  for (i = 0; i < l2_size; i++) {
    l2_regions[i] = sys_ralloc(r, 0); // mid level
    l1_regions[i] = sys_alloc(l1_size * sizeof(rid_t), l2_regions[i]);
    for (j = 0; j < l1_size; j++) {
      l1_regions[i][j] = sys_ralloc(l2_regions[i], 0); // bottom level
    }
  }

  // Create arrays
  buf = sys_alloc(l2_size * sizeof(unsigned char ***), r);
  for (i = 0; i < l2_size; i++) {
    buf[i] = sys_alloc(l1_size * sizeof(unsigned char **), l2_regions[i]);

    for (j = 0; j < l1_size; j++) {
      buf[i][j] = sys_alloc(l0_size * sizeof(unsigned char *), 
                            l1_regions[i][j]);

      if (l0_size > 1) {
        sys_balloc(buf_size * sizeof(float), l1_regions[i][j], l0_size, 
                   buf[i][j]);
      }
      else {
        buf[i][j][0] = sys_alloc(buf_size * sizeof(float), l1_regions[i][j]);
      }
    }
  }


  // Print we're starting
  printf("Tasks of %d cycles starting in %d L2 x %d L1 x %d L0 tile(s)\r\n", 
         delay, l2_size, l1_size, l0_size);

  // Start time
  time_start = sys_free_timer_get_ticks();

  // For all reps
  for (i = 0; i < num_reps; i++) {

    // Use one-level-up function to compute this wave
    #pragma myrmics task region inout(r) in(l2_regions, l1_regions, buf) \
                         safe(l2_regions, l1_regions, buf) \
                         in(l2_size, l1_size, l0_size, buf_size, delay)
    multilevel_myrmics_do_level3(r, l2_regions, l1_regions, buf, 
                                 l2_size, l1_size, l0_size, buf_size, delay);

    // // Emulate transfer to framebuffer
    // #pragma myrmics task region in(r)
    // multilevel_myrmics_wait_region(r);
  }


  // Stop time
  int l3_size = 0;
  #pragma myrmics task region in(r) \
                       in(time_start, l3_size, l2_size, l1_size, l0_size)
  multilevel_myrmics_time(r, time_start, l3_size, l2_size, l1_size, l0_size);
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics_l3_main(int l3_size, int l2_size, int l1_size, 
                                int l0_size, int buf_size, int num_reps, 
                                int delay) {

  rid_t                 r;
  rid_t                 *l3_regions;
  rid_t                 **l2_regions;
  rid_t                 ***l1_regions;
  unsigned char         *****buf;
  unsigned int          time_start;
  int                   i;
  int                   j;
  int                   k;


  // Create all-holding region
  r = sys_ralloc(0, 99); // top level

  // Create regions
  l3_regions = sys_alloc(l3_size * sizeof(rid_t), r);
  l2_regions = sys_alloc(l2_size * sizeof(rid_t *), r);
  l1_regions = sys_alloc(l1_size * sizeof(rid_t **), r);
  for (i = 0; i < l3_size; i++) {
    l3_regions[i] = sys_ralloc(r, 0); // high level
    l2_regions[i] = sys_alloc(l2_size * sizeof(rid_t), l3_regions[i]);
    l1_regions[i] = sys_alloc(l2_size * sizeof(rid_t *), l3_regions[i]);
    for (j = 0; j < l2_size; j++) {
      l2_regions[i][j] = sys_ralloc(l3_regions[i], 0); // mid level
      l1_regions[i][j] = sys_alloc(l1_size * sizeof(rid_t), l2_regions[i][j]);
      for (k = 0; k < l1_size; k++) {
        l1_regions[i][j][k] = sys_ralloc(l2_regions[i][j], 0); // low level
      }
    }
  }

  // Create arrays
  buf = sys_alloc(l3_size * sizeof(unsigned char ****), r);
  for (i = 0; i < l3_size; i++) {
    buf[i] = sys_alloc(l2_size * sizeof(unsigned char ***), l3_regions[i]);

    for (j = 0; j < l2_size; j++) {
      buf[i][j] = sys_alloc(l1_size * sizeof(unsigned char **), 
                            l2_regions[i][j]);

      for (k = 0; k < l1_size; k++) {
        buf[i][j][k] = sys_alloc(l0_size * sizeof(unsigned char *),
                                 l1_regions[i][j][k]);

        if (l0_size > 1) {
          sys_balloc(buf_size * sizeof(float), l1_regions[i][j][k], l0_size, 
                     buf[i][j][k]);
        }
        else {
          buf[i][j][k][0] = sys_alloc(buf_size * sizeof(float), 
                                      l1_regions[i][j][k]);
        }
      }
    }
  }


  // Print we're starting
  printf(
    "Tasks of %d cycles starting in %d L3 x %d L2 x %d L1 x %d L0 tile(s)\r\n", 
    delay, l3_size, l2_size, l1_size, l0_size);

  // Start time
  time_start = sys_free_timer_get_ticks();

  // For all reps
  for (i = 0; i < num_reps; i++) {

    // Use one-level-up function to compute this wave
    #pragma myrmics task region inout(r) \
                         in(l3_regions, l2_regions, l1_regions, buf) \
                         safe(l3_regions, l2_regions, l1_regions, buf) \
                         in(l3_size, l2_size, l1_size, l0_size, buf_size, delay)
    multilevel_myrmics_do_level4(r, l3_regions, l2_regions, l1_regions, buf, 
                                 l3_size, l2_size, l1_size, l0_size, buf_size, 
                                 delay);

    // // Emulate transfer to framebuffer
    // #pragma myrmics task region in(r)
    // multilevel_myrmics_wait_region(r);
  }


  // Stop time
  #pragma myrmics task region in(r) \
                       in(time_start, l3_size, l2_size, l1_size, l0_size)
  multilevel_myrmics_time(r, time_start, l3_size, l2_size, l1_size, l0_size);
}


// ===========================================================================
// ===========================================================================
void multilevel_myrmics(int l3_size, int l2_size, int l1_size, int l0_size, 
                        int blk_size, int num_reps, int delay) {

  if (l3_size > 0) {
    multilevel_myrmics_l3_main(l3_size, l2_size, l1_size, l0_size, blk_size, 
                               num_reps, delay);
  }
  else if (l2_size > 0) {
    multilevel_myrmics_l2_main(l2_size, l1_size, l0_size, blk_size, num_reps, 
                               delay);
  }
  else if (l1_size > 0) {
    multilevel_myrmics_l1_main(l1_size, l0_size, blk_size, num_reps, delay);
  }
  else if (l0_size > 0) {
    multilevel_myrmics_l0_main(l0_size, blk_size, num_reps, delay);
  }
  else {
    sys_abort();
  }
}

