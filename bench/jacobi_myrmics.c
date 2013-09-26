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
// Author        : Iakovos Mavroidis / Spyros Lyberis
// Abstract      : 2D Jacobi kernel, row-wise distribution, Myrmics version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: jacobi_myrmics.c,v $
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2013/01/23 15:22:18 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>

//#define DELAY_SUBREGION         500
//#define DELAY_SUBGRID           100
//#define DELAY_CHECKSUM          500
//#define EXTRA_SANITY_CHECKS

// ===========================================================================
// ===========================================================================
void jacobi_myrmics_do_checksum(float *cur_top, float *cur_core, 
                                float *cur_bot, int tile_rows, int num_cols, 
                                float *checksum) {
  int   i;
  int   j;
  float element;
  float old_checksum;

  old_checksum = *checksum;

  for (i = 0; i < tile_rows; i++) {
    for (j = 0; j < num_cols; j++) {
      element = (i == 0)             ? cur_top[j] :
                (i == tile_rows - 1) ? cur_bot[j] :
                                       cur_core[(i - 1) * num_cols + j];
      sys_assert((element >= 0.0) && (element <= 1.0));
      *checksum += element;
    }
  }
  printf("Partial checksum: %lf [%+lf]\r\n", 
         *checksum, *checksum - old_checksum);
}


// ===========================================================================
// ===========================================================================
void jacobi_myrmics_print_checksum(int num_rows, int num_cols, 
                                   unsigned int time_start, float *checksum) {
  float         average;
  unsigned int  time_stop;
  unsigned int  time;


#ifdef DELAY_CHECKSUM
  sys_timer_busy_wait_msec(DELAY_CHECKSUM);
#endif

  // Compute elapsed time
  time_stop = sys_free_timer_get_ticks();
  if (time_stop > time_start) {
    time = time_stop - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_stop);
  }

  printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);


  // Compute average
  average = *checksum / (float) (num_rows * num_cols);

  printf("Checksum: %lf\r\nAverage: %lf\r\n", *checksum, average);
}


// ===========================================================================
// ===========================================================================
void jacobi_myrmics_print_time(rid_t unused, unsigned int time_start) {
  unsigned int  time_stop;
  unsigned int  time;


  // Compute elapsed time
  time_stop = sys_free_timer_get_ticks();
  if (time_stop > time_start) {
    time = time_stop - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_stop);
  }

  printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);
}


// ===========================================================================
// ===========================================================================
void jacobi_myrmics_subgrid_compute(float *prev_bot, // old buf, prev tile bot
                                    float *nxt_top,  // old buf, next tile top
                                    float *cur_top,  // old buf, curr tile top
                                    float *cur_core, // old buf, curr tile core 
                                    float *cur_bot,  // old buf, curr tile bot 
                                    float *new_top,  // new buf, curr tile top
                                    float *new_core, // new buf, curr tile core
                                    float *new_bot,  // new buf, curr tile bot
                                    int start_row, int end_row, 
                                    int tile_rows, int num_cols) {

  int   i;
  int   j;
  float *north;
  float *center;
  float *south;
  float *new;


  //printf("%d: %X %X %X %X %X\r\n", sys_get_worker_id(), prev_bot, nxt_top, cur_top, cur_core, cur_bot);

#ifdef DELAY_SUBGRID
  sys_timer_busy_wait_msec(DELAY_SUBGRID);
#endif


  // Do extensive sanity checks
#ifdef EXTRA_SANITY_CHECKS
  for (i = 0; i < num_cols; i++) {
    if (prev_bot) {
      if (!((prev_bot[i] >= 0.0) && (prev_bot[i] <= 1.0))) {
        printf("%d: Error %X[%d] = %f\r\n", sys_get_worker_id(), prev_bot, i, prev_bot[i]);
        sys_abort();
      }
    }
    if (nxt_top) {
      if (!((nxt_top[i]  >= 0.0) && (nxt_top[i]  <= 1.0))) {
        printf("%d: Error %X[%d] = %f\r\n", sys_get_worker_id(), nxt_top, i, nxt_top[i]);
        sys_abort();
      }
    }
    if (!((cur_top[i]  >= 0.0) && (cur_top[i]  <= 1.0))) {
      printf("%d: Error %X[%d] = %f\r\n", sys_get_worker_id(), cur_top, i, cur_top[i]);
      sys_abort();
    }
    if (!((cur_bot[i]  >= 0.0) && (cur_bot[i]  <= 1.0))) {
      printf("%d: Error %X[%d] = %f\r\n", sys_get_worker_id(), cur_bot, i, cur_bot[i]);
      sys_abort();
    }
    for (j = 0; j < tile_rows - 2; j++) {
      if (!((cur_core[j * num_cols + i] >= 0.0) && 
            (cur_core[j * num_cols + i] <= 1.0))) {
        printf("%d: Error %X[%d] = %f\r\n", sys_get_worker_id(), cur_core, j * num_cols + i, cur_core[j * num_cols + i]);
        sys_abort();
      }
    }
  }
#endif

  for (i = start_row; i <= end_row; i++) {

    if (i == 0) {
      sys_assert(prev_bot);

      north  = prev_bot;
      center = cur_top;
      south  = cur_core + (i - 0) * num_cols;
      
      new    = new_top;
    }
    else if (i == 1) {
      north  = cur_top;
      center = cur_core + (i - 1) * num_cols;
      south  = cur_core + (i - 0) * num_cols;
      
      new    = new_core + (i - 1) * num_cols;
    }
    else if (i == tile_rows - 2) {
      north  = cur_core + (i - 2) * num_cols;
      center = cur_core + (i - 1) * num_cols;
      south  = cur_bot;
      
      new    = new_core + (i - 1) * num_cols;
    }
    else if (i == tile_rows - 1) {
      sys_assert(nxt_top);

      north  = cur_core + (i - 2) * num_cols;
      center = cur_bot;
      south  = nxt_top;
      
      new    = new_bot;
    }
    else {
      north  = cur_core + (i - 2) * num_cols;
      center = cur_core + (i - 1) * num_cols;
      south  = cur_core + (i - 0) * num_cols;
      
      new    = new_core + (i - 1) * num_cols;
    }

    for (j = 1; j < num_cols - 1; j++) {
      new[j] = (north[j] + center[j - 1] + center[j + 1] + south[j]) * 0.25;
    }
  }
}


// ===========================================================================
// ===========================================================================
void jacobi_myrmics_subregion_compute(rid_t cur_region, float **cur_top, 
                                      float **cur_core, float **cur_bot,
                                      float *ext_prev_bot, float *ext_nxt_top,
                                      rid_t new_region, float **new_top, 
                                      float **new_core, float **new_bot,
                                      int region, int num_regions, 
                                      int region_tiles, int tile_rows, 
                                      int num_cols) {

  int           start_row;
  int           end_row;
  float         *prev_bot_row;
  float         *nxt_top_row;
  float         *cur_core_row;
  float         *cur_top_row;
  float         *cur_bot_row;
  float         *new_core_row;
  float         *new_top_row;
  float         *new_bot_row;
  int           tile;

#ifdef DELAY_SUBREGION
  sys_timer_busy_wait_msec(DELAY_SUBREGION);
#endif

  for (tile = 0; tile < region_tiles; tile++) {

    if ((tile > 0) && (region > 0)) {
      start_row = 0;
    }
    else {
      start_row = 1;
    }

    if ((tile < region_tiles - 1) && (region < num_regions - 1)) {
      end_row = tile_rows - 1;
    }
    else {
      end_row = tile_rows - 2;
    }

    prev_bot_row = (tile > 0)                ? cur_bot[tile-1] : ext_prev_bot;
    nxt_top_row  = (tile < region_tiles - 1) ? cur_top[tile+1] : ext_nxt_top;
    cur_top_row  =                             cur_top[tile];
    cur_core_row =                             cur_core[tile];
    cur_bot_row  =                             cur_bot[tile];
    new_top_row  =                             new_top[tile];
    new_core_row =                             new_core[tile];
    new_bot_row  =                             new_bot[tile];

//printf("prv 0x%X nxt 0x%X cur 0x%X 0x%X 0x%X new 0x%X 0x%X 0x%X\r\n", prev_bot_row, nxt_top_row, cur_top_row, cur_core_row, cur_bot_row, new_top_row, new_core_row, new_bot_row);

    #pragma myrmics task in(prev_bot_row, nxt_top_row, \
                            cur_top_row, cur_core_row, cur_bot_row) \
                         inout(new_top_row, new_core_row, new_bot_row) \
                         in(start_row, end_row, tile_rows, num_cols)
    jacobi_myrmics_subgrid_compute(prev_bot_row, nxt_top_row, 
                                   cur_top_row, cur_core_row, cur_bot_row, 
                                   new_top_row, new_core_row, new_bot_row,
                                   start_row, end_row, 
                                   tile_rows, num_cols);
  }
}


// ===========================================================================
// ===========================================================================
void jacobi_myrmics_subregion_copy_boundaries(float *prv_bot, float *nxt_top, 
                                              float *cur_prv_bot, 
                                              float *cur_nxt_top, 
                                              int num_cols) {
  int i;

  if (prv_bot) {
    for (i = 0; i < num_cols; i++) {
      cur_prv_bot[i] = prv_bot[i];
    }
  }

  if (nxt_top) {
    for (i = 0; i < num_cols; i++) {
      cur_nxt_top[i] = nxt_top[i];
    }
  }
}


// ===========================================================================
// ===========================================================================
void jacobi_myrmics(int num_regions, int num_tiles, int num_rows, 
                    int num_cols, int num_iter, int verify) {

  rid_t         r;
  rid_t         *grid0_regions;
  rid_t         *grid1_regions;
  float         ***grid0_core;
  float         ***grid0_top;
  float         ***grid0_bot;
  float         **grid0_prv_bot;
  float         **grid0_nxt_top;
  float         ***grid1_core;
  float         ***grid1_top;
  float         ***grid1_bot;
  float         **grid1_prv_bot;
  float         **grid1_nxt_top;
  float         *checksum;
  int           tile_rows;
  int           region_tiles;
  int           i;
  int           j;
  int           start_row;
  int           end_row;
  int           tile;
  int           region;
  int           rep;
  unsigned int  time_start;
  float         *prv_bot;
  float         *nxt_top;
  rid_t         cur_region;
  float         **cur_top;
  float         **cur_core;
  float         **cur_bot;
  float         *cur_prv_bot;
  float         *cur_nxt_top;
  rid_t         new_region;
  float         **new_top;
  float         **new_core;
  float         **new_bot;
  float         *cur_top_row;
  float         *cur_core_row;
  float         *cur_bot_row;


  // Sanity checks
  if (num_rows % num_tiles) {
    printf("%d rows not divisible by %d tiles\r\n", num_rows, num_tiles);
    return;
  }
  tile_rows = num_rows / num_tiles;
  if (tile_rows < 4) {
    printf("rows per tile must be >= 4 (they are %d now)\r\n", tile_rows);
    return;
  }
  if (num_tiles % num_regions) {
    printf("%d tiles not divisible by %d regions\r\n", num_tiles, 
              num_regions);
    return;
  }
  region_tiles = num_tiles / num_regions;


  // Infomercial
  printf(
      "Jacobi 2D-row of %d x %d starting in %d tile(s) and %d region(s)\r\n",
      num_rows, num_cols, num_tiles, num_regions);

  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Create regions
  grid0_regions = sys_alloc(num_regions * sizeof(rid_t), r);
  grid1_regions = sys_alloc(num_regions * sizeof(rid_t), r);
  for (i = 0; i < num_regions; i++) {
    grid0_regions[i] = sys_ralloc(r, 0); // lowest level
  }
  for (i = 0; i < num_regions; i++) {
    grid1_regions[i] = sys_ralloc(r, 0); // lowest level
  }

  // Allocate buffers
  grid0_core    = sys_alloc(num_regions * sizeof(float **), r);
  grid0_top     = sys_alloc(num_regions * sizeof(float **), r);
  grid0_bot     = sys_alloc(num_regions * sizeof(float **), r);
  grid0_prv_bot = sys_alloc(num_regions * sizeof(float *), r);
  grid0_nxt_top = sys_alloc(num_regions * sizeof(float *), r);

  grid1_core    = sys_alloc(num_regions * sizeof(float **), r);
  grid1_top     = sys_alloc(num_regions * sizeof(float **), r);
  grid1_bot     = sys_alloc(num_regions * sizeof(float **), r);
  grid1_prv_bot = sys_alloc(num_regions * sizeof(float *), r);
  grid1_nxt_top = sys_alloc(num_regions * sizeof(float *), r);

  for (i = 0; i < num_regions; i++) {
    grid0_core[i] = sys_alloc(region_tiles * sizeof(float *), grid0_regions[i]);
    grid0_top[i]  = sys_alloc(region_tiles * sizeof(float *), grid0_regions[i]);
    grid0_bot[i]  = sys_alloc(region_tiles * sizeof(float *), grid0_regions[i]);

    grid1_core[i] = sys_alloc(region_tiles * sizeof(float *), grid1_regions[i]);
    grid1_top[i]  = sys_alloc(region_tiles * sizeof(float *), grid1_regions[i]);
    grid1_bot[i]  = sys_alloc(region_tiles * sizeof(float *), grid1_regions[i]);

    if (region_tiles > 1) {
      sys_balloc((tile_rows - 2) * num_cols * sizeof(float), grid0_regions[i], 
                 region_tiles, grid0_core[i]);
      sys_balloc(num_cols * sizeof(float), grid0_regions[i], 
                 region_tiles, grid0_top[i]);
      sys_balloc(num_cols * sizeof(float), grid0_regions[i], 
                 region_tiles, grid0_bot[i]);

      sys_balloc((tile_rows - 2) * num_cols * sizeof(float), grid1_regions[i], 
                 region_tiles, grid1_core[i]);
      sys_balloc(num_cols * sizeof(float), grid1_regions[i], 
                 region_tiles, grid1_top[i]);
      sys_balloc(num_cols * sizeof(float), grid1_regions[i], 
                 region_tiles, grid1_bot[i]);
    }
    else {
      grid0_core[i][0] = sys_alloc((tile_rows - 2) * num_cols * sizeof(float), 
                                   grid0_regions[i]);
      grid0_top[i][0] = sys_alloc(num_cols * sizeof(float), grid0_regions[i]);
      grid0_bot[i][0] = sys_alloc(num_cols * sizeof(float), grid0_regions[i]);

      grid1_core[i][0] = sys_alloc((tile_rows - 2) * num_cols * sizeof(float), 
                                   grid1_regions[i]);
      grid1_top[i][0] = sys_alloc(num_cols * sizeof(float), grid1_regions[i]); 
      grid1_bot[i][0] = sys_alloc(num_cols * sizeof(float), grid1_regions[i]); 
    }

    grid0_prv_bot[i] = sys_alloc(num_cols * sizeof(float), grid0_regions[i]);
    grid0_nxt_top[i] = sys_alloc(num_cols * sizeof(float), grid0_regions[i]);

    grid1_prv_bot[i] = sys_alloc(num_cols * sizeof(float), grid1_regions[i]);
    grid1_nxt_top[i] = sys_alloc(num_cols * sizeof(float), grid1_regions[i]);
  }

  checksum = sys_alloc(sizeof(float), 0);


  // Initialize grid and extra row buffers to 0.0, grid borders to 1.0
  for (region = 0; region < num_regions; region++) {
    for (tile = 0; tile < region_tiles; tile++) {
      for (j = 0; j < num_cols; j++) {

        // Initialize top line buffers
        if ((tile == 0) || (j == 0) || (j == num_cols - 1)) {
            grid0_top[region][tile][j] = 1.0;
            grid1_top[region][tile][j] = 1.0;
        } else {
            grid0_top[region][tile][j] = 0.0;
            grid1_top[region][tile][j] = 0.0;
        }

        // Initialize bottom line buffers
        if ((tile == num_tiles - 1) || (j == 0) || (j == num_cols - 1)) {
            grid0_bot[region][tile][j] = 1.0;
            grid1_bot[region][tile][j] = 1.0;
        } else {
            grid0_bot[region][tile][j] = 0.0;
            grid1_bot[region][tile][j] = 0.0;
        }

        // Initialize core buffers
        for (i = 0; i < tile_rows - 2; i++) {
          if ((j == 0) || (j == num_cols - 1)) {
            grid0_core[region][tile][i * num_cols + j] = 1.0;
            grid1_core[region][tile][i * num_cols + j] = 1.0;
          }
          else {
            grid0_core[region][tile][i * num_cols + j] = 0.0;
            grid1_core[region][tile][i * num_cols + j] = 0.0;
          }
        }
      }
    }
  }


  // Start time
  time_start = sys_free_timer_get_ticks();


  // Compute
  for (rep = 0; rep < num_iter; rep++) {

    for (region = 0; region < num_regions; region++) {
    
      if (rep % 2) {
        prv_bot     = (region > 0) ? grid1_bot[region-1][region_tiles-1] : NULL;
        nxt_top     = (region < num_regions-1) ? grid1_top[region+1][0] : NULL;

        cur_region  = grid1_regions[region];
        cur_top     = grid1_top[region];
        cur_core    = grid1_core[region];
        cur_bot     = grid1_bot[region];
        cur_prv_bot = grid1_prv_bot[region];
        cur_nxt_top = grid1_nxt_top[region];

        new_region  = grid0_regions[region];
        new_top     = grid0_top[region];
        new_core    = grid0_core[region];
        new_bot     = grid0_bot[region];
      }
      else {
        prv_bot     = (region > 0) ? grid0_bot[region-1][region_tiles-1] : NULL;
        nxt_top     = (region < num_regions-1) ? grid0_top[region+1][0] : NULL;

        cur_region  = grid0_regions[region];
        cur_top     = grid0_top[region];
        cur_core    = grid0_core[region];
        cur_bot     = grid0_bot[region];
        cur_prv_bot = grid0_prv_bot[region];
        cur_nxt_top = grid0_nxt_top[region];

        new_region  = grid1_regions[region];
        new_top     = grid1_top[region];
        new_core    = grid1_core[region];
        new_bot     = grid1_bot[region];
      }


      #pragma myrmics task in(prv_bot, nxt_top) \
                           inout(cur_prv_bot, cur_nxt_top) \
                           in(num_cols)
      jacobi_myrmics_subregion_copy_boundaries(prv_bot, nxt_top, 
                                               cur_prv_bot, cur_nxt_top,
                                               num_cols);


      #pragma myrmics task region in(cur_region) \
                           in(cur_top, cur_core, cur_bot, \
                              cur_prv_bot, cur_nxt_top) \
                           safe(cur_top, cur_core, cur_bot, \
                                cur_prv_bot, cur_nxt_top) \
                           region inout(new_region) \
                           in(new_top, new_core, new_bot) \
                           safe(new_top, new_core, new_bot) \
                           in(region, num_regions, region_tiles, \
                              tile_rows, num_cols)
      jacobi_myrmics_subregion_compute(cur_region, cur_top, cur_core, cur_bot,
                                       cur_prv_bot, cur_nxt_top,
                                       new_region, new_top, new_core, new_bot,
                                       region, num_regions, region_tiles, 
                                       tile_rows, num_cols);

    }
  }




  // Compute final average and checksum 
  *checksum = 0.0;
  if (verify == 1) {
    for (region = 0; region < num_regions; region++) {
      for (tile = 0; tile < region_tiles; tile++) {
        if (num_iter % 2) {
          cur_top_row  = grid1_top[region][tile];
          cur_core_row = grid1_core[region][tile];
          cur_bot_row  = grid1_bot[region][tile];
        }
        else {
          cur_top_row  = grid0_top[region][tile];
          cur_core_row = grid0_core[region][tile];
          cur_bot_row  = grid0_bot[region][tile];
        }

        #pragma myrmics task in(cur_top_row, cur_core_row, cur_bot_row) \
                             in(tile_rows, num_cols) inout(checksum)
        jacobi_myrmics_do_checksum(cur_top_row, cur_core_row, cur_bot_row, 
                                   tile_rows, num_cols, checksum);
      }
    }

    // Print the time and checksum
    #pragma myrmics task in(num_rows, num_cols, time_start, checksum)
    jacobi_myrmics_print_checksum(num_rows, num_cols, time_start, checksum);
  }
  else if (verify == -1) {
    // Just print the time
    #pragma myrmics task region in(r) in(time_start)
    jacobi_myrmics_print_time(r, time_start);
  }

}

