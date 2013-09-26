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
// Abstract      : Matrix multiplication benchmark
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: matrix_mult_myrmics.c,v $
// CVS revision  : $Revision: 1.11 $
// Last modified : $Date: 2013/01/16 17:01:02 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>

// ===========================================================================
// ===========================================================================
void matrix_mult_myrmics_do_submatrix(float *a, float *b, float *c, 
                                      int blk_size) {
  int i;
  int j;
  int k;

  //printf("%d: do_submatrix: %08X x %08X -> %08X\r\n", sys_get_worker_id(), a, b, c);

  for (i = 0; i < blk_size; i++) {
    for (j = 0; j < blk_size; j++) {
      for (k = 0; k < blk_size; k++) {
        c[i * blk_size + j] += a[i * blk_size + k] * b[k * blk_size + j];
      }
    }
  }
}


// ===========================================================================
// ===========================================================================
void matrix_mult_myrmics_do_region(rid_t a_rid, rid_t b_rid, rid_t c_rid, 
                                   float **a, float **b, float **c, 
                                   int tile_size, int blk_size) {
  int   phase;
  int   tile_row;
  int   tile_col;
  float *a_tile;
  float *b_tile;
  float *c_tile;


  // For all the phases in the algorithm
  for (phase = 0; phase < tile_size; phase++) {

    // For all tiles
    for (tile_row = 0; tile_row < tile_size; tile_row++) {
      for (tile_col = 0; tile_col < tile_size; tile_col++) { 

        a_tile = a[tile_row * tile_size + phase];
        b_tile = b[phase * tile_size + tile_col];
        c_tile = c[tile_row * tile_size + tile_col];

        #pragma myrmics task in(a_tile, b_tile) inout(c_tile) \
                             in(blk_size)
        matrix_mult_myrmics_do_submatrix(a_tile, b_tile, c_tile, blk_size);
      }
    }
  }
}


// ===========================================================================
// ===========================================================================
void matrix_mult_myrmics_copy_region(rid_t src_region, rid_t dst_region,
                                     float **src, float **dst, 
                                     int tile_size, int blk_size) {
  int i;
  int j;

  for (i = 0; i < tile_size * tile_size; i++) {
    for (j = 0; j < blk_size * blk_size; j++) {
      dst[i][j] = src[i][j];
    }
  }
}


// ===========================================================================
// ===========================================================================
static inline int reg_idx(int row, int col, int region_size, int whole_size) {

  int reg_row = (row * region_size) / whole_size;
  int reg_col = (col * region_size) / whole_size;

  return (reg_row * region_size + reg_col);
}

// ===========================================================================
// ===========================================================================
static inline int tile_idx(int row, int col, int tile_size, int tile_blk_size) {
  
  row %= tile_blk_size;
  col %= tile_blk_size;

  int tile_row = (row * tile_size) / tile_blk_size;
  int tile_col = (col * tile_size) / tile_blk_size;

  return (tile_row * tile_size + tile_col);
}

// ===========================================================================
// ===========================================================================
static inline int blk_idx(int row, int col, int blk_size) {
  
  int blk_row = row % blk_size;
  int blk_col = col % blk_size;

  return (blk_row * blk_size + blk_col);
}


// ===========================================================================
// ===========================================================================
void matrix_mult_myrmics_verify(rid_t r, float ***a, float ***b, float ***c,
                                int region_size, int tile_size, int blk_size) {

  int   whole_size;
  int   tile_blk_size;
  float verify_val;
  int   i;
  int   j;
  int   k;

  printf("Multiplication finished, verifying...\r\n");

  whole_size = region_size * tile_size * blk_size;
  tile_blk_size = tile_size * blk_size;

  for (i = 0; i < whole_size; i++) {
    for (j = 0; j < whole_size; j++) {
      verify_val = 0;
      for (k = 0; k < whole_size; k++) {
        
        // We want: verify_val += A[i][k] * B[k][j];

        verify_val += a[reg_idx (i, k, region_size, whole_size)]
                       [tile_idx(i, k, tile_size, tile_blk_size)]
                       [blk_idx (i, k, blk_size)] *

                      b[reg_idx (k, j, region_size, whole_size)]
                       [tile_idx(k, j, tile_size, tile_blk_size)]
                       [blk_idx (k, j, blk_size)];
      }
      // Check C[i][j]
      if (c[reg_idx (i, j, region_size, whole_size)]
           [tile_idx(i, j, tile_size, tile_blk_size)]
           [blk_idx (i, j, blk_size)] != verify_val) {
        printf("Verification [1;31mFAILED[0m at C[%d, %d]\r\n", i, j);
        sys_abort();
      }
    }
  }

  printf("Verification [1;32mPASSED[0m\r\n");
}


// ===========================================================================
// ===========================================================================
void matrix_mult_myrmics_checksum(rid_t r, float ***c, int num_regions, 
                                  int tile_size_sq, int blk_size_sq) {

  float checksum;
  int   i;
  int   j;
  int   k;

  printf("Multiplication finished, checksumming...\r\n");

  checksum = 0.0;
  for (i = 0; i < num_regions; i++) {
    for (j = 0; j < tile_size_sq; j++) {
      for (k = 0; k < blk_size_sq; k++) {
        
        checksum += c[i][j][k];
      }
    }
  }

  printf("Checksum: %f\r\n", checksum);
}


// ===========================================================================
// ===========================================================================
void matrix_mult_myrmics_time(rid_t unused, unsigned int time_start) {

  unsigned int          time_stop;
  unsigned int          time;

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
// matrix_mult_myrmics()        Myrmics version of dense matrix multiplication.
//                              Tiles are laid out in a 2D space
//                              tile_size * tile_size layout. Each tile
//                              has a portion of matrices A and B (each such
//                              portion blk_size * blk_size numbers) and
//                              computes a portion of matrix C = A * B.
// ===========================================================================
// * INPUTS
//   int region_size            Number of region rows and columns. If set to
//                              1, run is effectively non-hierarchical (the
//                              matrix_mult_myrmics_do_region() function is
//                              called only once).
//   int tile_size              Number of tile rows and columns per region
//   int blk_size               Number of rows and columns per tile
//                              The total size of A, B and C is
//                              region_size ^ 2 * tile_size ^ 2 * blk_size ^ 2 
//                              elements
//   int verify                 -1: don't verify the multiplication
//                               1: a final task does the verification
//                               2: a final task just takes a checksum
//   int dummy_regions          -1: don't create any dummy regions
//                              >0: create this many dummy regions per
//                                  c/wrk_a/wrk_b loop. Use it to force
//                                  the handling of c/wrk_a/wrk_b on the same
//                                  leaf scheduler (see code for details)
// ===========================================================================
void matrix_mult_myrmics(int region_size, int tile_size, int blk_size, 
                         int verify, int dummy_regions) {

  rid_t                 r;
  rid_t                 *a_regions;
  rid_t                 *b_regions;
  rid_t                 *c_regions;
  rid_t                 **wrk_a_regions;
  rid_t                 **wrk_b_regions;
  int                   num_regions;
  int                   tile_size_sq;
  int                   blk_size_sq;
  int                   whole_size;
  int                   phase;
  float                 ***a;
  float                 ***b;
  float                 ***c;
  float                 ****wrk_a;
  float                 ****wrk_b;
  unsigned int          time_start;
  float                 **cur_a;
  float                 **cur_b;
  float                 **cur_c;
  float                 **cur_wrk_a;
  float                 **cur_wrk_b;
  rid_t                 a_region;
  rid_t                 b_region;
  rid_t                 c_region;
  rid_t                 wrk_a_region;
  rid_t                 wrk_b_region;
  int                   idx_a;
  int                   idx_b;
  int                   idx_c;
  int                   region_row;
  int                   region_col;
  int                   i;
  int                   j;
  int                   k;


  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level

  // Basic sizes
  num_regions = region_size * region_size;
  tile_size_sq = tile_size * tile_size;
  blk_size_sq = blk_size * blk_size;

  // Create regions
  a_regions = sys_alloc(num_regions * sizeof(rid_t), r);
  b_regions = sys_alloc(num_regions * sizeof(rid_t), r);
  c_regions = sys_alloc(num_regions * sizeof(rid_t), r);

  wrk_a_regions = sys_alloc(2 * sizeof(rid_t *), r);
  wrk_a_regions[0] = sys_alloc(num_regions * sizeof(rid_t), r);
  wrk_a_regions[1] = sys_alloc(num_regions * sizeof(rid_t), r);

  wrk_b_regions = sys_alloc(2 * sizeof(rid_t *), r);
  wrk_b_regions[0] = sys_alloc(num_regions * sizeof(rid_t), r);
  wrk_b_regions[1] = sys_alloc(num_regions * sizeof(rid_t), r);

  for (i = 0; i < num_regions; i++) {
    a_regions[i] = sys_ralloc(r, 99); // highest level
    b_regions[i] = sys_ralloc(r, 99); // highest level
  }

  // We use separate loops for all low-level regions, to ensure that
  // c_regions[y] is on the same scheduler with wrk_a_regions[*][y] and
  // wrk_b_regions[*][y] for all y values. To ensure this, num regions (i.e.
  // region_size squared) must be set to equal to the number of L0 schedulers
  // (for two-level hierarchies). If this is not possible, use dummy_regions
  // to create extra dummies.
  for (i = 0; i < num_regions; i++) {
    c_regions[i] = sys_ralloc(r, 0); // lowest level
  }
  for (i = 0; i < dummy_regions; i++) {
    sys_ralloc(r, 0); // lowest level, discard the outcome
  }

  for (i = 0; i < 2; i++) {
    for (j = 0; j < num_regions; j++) {
      wrk_a_regions[i][j] = sys_ralloc(r, 0); // lowest level
    }
    for (j = 0; j < dummy_regions; j++) {
      sys_ralloc(r, 0); // lowest level, discard the outcome
    }

    for (j = 0; j < num_regions; j++) {
      wrk_b_regions[i][j] = sys_ralloc(r, 0); // lowest level
    }
    for (j = 0; j < dummy_regions; j++) {
      sys_ralloc(r, 0); // lowest level, discard the outcome
    }
  }
  
  // Create arrays
  a     = sys_alloc(num_regions * sizeof(float **), r);
  b     = sys_alloc(num_regions * sizeof(float **), r);
  c     = sys_alloc(num_regions * sizeof(float **), r);

  wrk_a = sys_alloc(2 * sizeof(float ***), r);
  wrk_a[0] = sys_alloc(num_regions * sizeof(float **), r);
  wrk_a[1] = sys_alloc(num_regions * sizeof(float **), r);

  wrk_b = sys_alloc(2 * sizeof(float ***), r);
  wrk_b[0] = sys_alloc(num_regions * sizeof(float **), r);
  wrk_b[1] = sys_alloc(num_regions * sizeof(float **), r);

  for (i = 0; i < num_regions; i++) {
    a[i]        = sys_alloc(tile_size_sq*sizeof(float *), a_regions[i]);
    b[i]        = sys_alloc(tile_size_sq*sizeof(float *), b_regions[i]);
    c[i]        = sys_alloc(tile_size_sq*sizeof(float *), c_regions[i]);
    wrk_a[0][i] = sys_alloc(tile_size_sq*sizeof(float *), wrk_a_regions[0][i]);
    wrk_a[1][i] = sys_alloc(tile_size_sq*sizeof(float *), wrk_a_regions[1][i]);
    wrk_b[0][i] = sys_alloc(tile_size_sq*sizeof(float *), wrk_b_regions[0][i]);
    wrk_b[1][i] = sys_alloc(tile_size_sq*sizeof(float *), wrk_b_regions[1][i]);

    if (tile_size_sq > 1) {
      sys_balloc(blk_size_sq * sizeof(float), a_regions[i],
                 tile_size_sq, a[i]);
      sys_balloc(blk_size_sq * sizeof(float), b_regions[i],
                 tile_size_sq, b[i]);
      sys_balloc(blk_size_sq * sizeof(float), c_regions[i], 
                 tile_size_sq, c[i]);
      sys_balloc(blk_size_sq * sizeof(float), wrk_a_regions[0][i], 
                 tile_size_sq, wrk_a[0][i]);
      sys_balloc(blk_size_sq * sizeof(float), wrk_a_regions[1][i], 
                 tile_size_sq, wrk_a[1][i]);
      sys_balloc(blk_size_sq * sizeof(float), wrk_b_regions[0][i], 
                 tile_size_sq, wrk_b[0][i]);
      sys_balloc(blk_size_sq * sizeof(float), wrk_b_regions[1][i], 
                 tile_size_sq, wrk_b[1][i]);
    }
    else {
      a[i][0]        = sys_alloc(blk_size_sq*sizeof(float), 
                                 a_regions[i]);
      b[i][0]        = sys_alloc(blk_size_sq*sizeof(float), 
                                 b_regions[i]);
      c[i][0]        = sys_alloc(blk_size_sq*sizeof(float), 
                                 c_regions[i]);
      wrk_a[0][i][0] = sys_alloc(blk_size_sq*sizeof(float), 
                                 wrk_a_regions[0][i]);
      wrk_a[1][i][0] = sys_alloc(blk_size_sq*sizeof(float), 
                                 wrk_a_regions[1][i]);
      wrk_b[0][i][0] = sys_alloc(blk_size_sq*sizeof(float), 
                                 wrk_b_regions[0][i]);
      wrk_b[1][i][0] = sys_alloc(blk_size_sq*sizeof(float), 
                                 wrk_b_regions[1][i]);
    }
  }


  // Initialize A and B with some values, zero out C
  for (i = 0; i < num_regions; i++) {
    for (j = 0; j < tile_size_sq; j++) {
      for (k = 0; k < blk_size_sq; k++) {
        a[i][j][k] = i + j + k;
        b[i][j][k] = i - j - k;
        c[i][j][k] = 0.0F;
      }
    }
  }

  // Print we're starting
  whole_size = region_size * tile_size * blk_size;
  printf("Matrix multiplication of %d x %d starting split into %d tile(s)\r\n", 
         whole_size, whole_size, num_regions * tile_size_sq);


  // Start time
  time_start = sys_free_timer_get_ticks();


  // For all the phases in the algorithm
  for (phase = 0; phase < region_size; phase++) {

    // For all regions
    for (region_row = 0; region_row < region_size; region_row++) {
      for (region_col = 0; region_col < region_size; region_col++) { 

        idx_a = region_row * region_size + phase;
        idx_b = phase * region_size + region_col;
        idx_c = region_row * region_size + region_col;

        a_region     = a_regions[idx_a];
        b_region     = b_regions[idx_b];
        c_region     = c_regions[idx_c];
        wrk_a_region = wrk_a_regions[phase & 1][idx_c];
        wrk_b_region = wrk_b_regions[phase & 1][idx_c];

        cur_a     = a[idx_a];
        cur_b     = b[idx_b];
        cur_c     = c[idx_c];
        cur_wrk_a = wrk_a[phase & 1][idx_c];
        cur_wrk_b = wrk_b[phase & 1][idx_c];

        #pragma myrmics task region in(a_region) \
                             region inout(wrk_a_region) \
                             in(cur_a, cur_wrk_a) safe(cur_a, cur_wrk_a) \
                             in(tile_size, blk_size)
        matrix_mult_myrmics_copy_region(a_region, wrk_a_region,
                                        cur_a, cur_wrk_a,
                                        tile_size, blk_size);

        #pragma myrmics task region in(b_region) \
                             region inout(wrk_b_region) \
                             in(cur_b, cur_wrk_b) safe(cur_b, cur_wrk_b) \
                             in(tile_size, blk_size)
        matrix_mult_myrmics_copy_region(b_region, wrk_b_region,
                                        cur_b, cur_wrk_b,
                                        tile_size, blk_size);

        #pragma myrmics task region in(wrk_a_region, wrk_b_region) \
                             region inout(c_region) \
                             in(cur_wrk_a, cur_wrk_b, cur_c) \
                             safe(cur_wrk_a, cur_wrk_b, cur_c) \
                             in(tile_size, blk_size)
        matrix_mult_myrmics_do_region(wrk_a_region, wrk_b_region, c_region, 
                                      cur_wrk_a, cur_wrk_b, cur_c,
                                      tile_size, blk_size);
      }
    }
  }


  // Stop time
  #pragma myrmics task region in(r) in(time_start)
  matrix_mult_myrmics_time(r, time_start);


  // Gather A, B and C's and verify or do a checksum
  if (verify == 1) {
     
    #pragma myrmics task region in(r) \
                         in(a, b, c) safe(a, b, c) \
                         in(region_size, tile_size, blk_size)
    matrix_mult_myrmics_verify(r, a, b, c,
                               region_size, tile_size, blk_size);
  }
  else if (verify == 2) {
     
    #pragma myrmics task region in(r) in(c) safe(c) \
                         in(num_regions, tile_size_sq, blk_size_sq)
    matrix_mult_myrmics_checksum(r, c, num_regions, tile_size_sq, blk_size_sq);
  }

}
