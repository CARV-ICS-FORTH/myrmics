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
// Abstract      : Parallel bitonic sort, Myrmics version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: bitonic_myrmics.c,v $
// CVS revision  : $Revision: 1.7 $
// Last modified : $Date: 2013/01/24 10:20:00 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>


// ===========================================================================
// ===========================================================================
void bitonic_myrmics_merge_even(int *src1, int *src2, int *dst, 
                                int tile_size, int verify) {

  int i;
  int j = 0;
  int k = 0;

  //printf("%d: merge_even src1 = 0x%08X, src2 = 0x%08X, dst = 0x%08X\r\n",
  //        sys_get_worker_id(), src1, src2, dst);

  if (verify == 1) {
    for (i = 1; i < tile_size; i++) {
      if (src1[i] < src1[i-1]) {
        printf("%d: merge_even [31;1mFAILED[0m for src1[%d]\r\n", 
               sys_get_worker_id(), i);
        sys_abort();
      }
      if (src2[i] < src2[i-1]) {
        printf("%d: merge_even [31;1mFAILED[0m for src2[%d]\r\n", 
               sys_get_worker_id(), i);
        sys_abort();
      }
    }
  }

  for (i = 0; i < tile_size; i++) {
    if (src1[j] <= src2[k]) {
      dst[i] = src1[j++];
    } 
    else {
      dst[i] = src2[k++];
    }
  }
}

// ===========================================================================
// ===========================================================================
void bitonic_myrmics_merge_odd(int *src1, int *src2, int *dst, 
                               int tile_size, int verify) {

  int i;
  int j = tile_size - 1;
  int k = tile_size - 1;

  //printf("%d: merge_odd src1 = 0x%08X, src2 = 0x%08X, dst = 0x%08X\r\n",
  //        sys_get_worker_id(), src1, src2, dst);

  if (verify == 1) {
    for (i = tile_size - 2; i >= 0; i--) {
      if (src1[i] > src1[i+1]) {
        printf("%d: merge_odd [31;1mFAILED[0m for src1[%d]\r\n", 
               sys_get_worker_id(), i);
        sys_abort();
      }
      if (src2[i] > src2[i+1]) {
        printf("%d: merge_odd [31;1mFAILED[0m for src2[%d]\r\n", 
               sys_get_worker_id(), i);
        sys_abort();
      }
    }
  }

  for (i = tile_size - 1; i >= 0; i--) {
    if (src1[j] >= src2[k]) {
      dst[i] = src1[j--];
    }
    else {
      dst[i] = src2[k--];
    }
  }
}


// ===========================================================================
// ===========================================================================
void bitonic_myrmics_region_merge(rid_t src1_id, int **src1_ptr, 
                                  rid_t src2_id, int **src2_ptr,
                                  rid_t dst_id,  int **dst_ptr, 
                                  int src1_i, int src2_i, int mask, int mask2,
                                  int region_size, int tile_size, int verify) {

  int *src1;
  int *src2;
  int *dst;
  int src1_j;
  int src2_j;
  int src1_ij;
  int src2_ij;


  for (src1_j = 0; src1_j < region_size; src1_j++) {

    // Select tiles
    if (mask2 < region_size) {
      src2_j = src1_j ^ mask2;
    }
    else {
      src2_j = src1_j;
    }

    src1 = src1_ptr[src1_j];
    src2 = src2_ptr[src2_j];
    dst  = dst_ptr[src1_j];

//printf("merge 0x%X 0x%X 0x%X [%d %d, %d %d]\r\n", src1, src2, dst, src1_i, src1_j, src2_i, src2_j);

    // Merge
    src1_ij = src1_i * region_size + src1_j;
    src2_ij = src2_i * region_size + src2_j;
    if ((((src1_ij & mask) == 0) && (src1_ij < src2_ij)) ||
        (((src1_ij & mask) != 0) && (src1_ij > src2_ij))) {

      #pragma myrmics task in(src1, src2) inout(dst) in(tile_size, verify)
      bitonic_myrmics_merge_even(src1, src2, dst, tile_size, verify);
    }
    else {

      #pragma myrmics task in(src1, src2) inout(dst) in(tile_size, verify)
      bitonic_myrmics_merge_odd(src1, src2, dst, tile_size, verify);
    }
  }
}


// ===========================================================================
// ===========================================================================
void bitonic_myrmics_tile_sort(int *tile, int tile_size) {

  int i;
  int j;
  int k;
  int l;
  int tmp;

  //printf("sorting tile 0x%X, %d elements\r\n", tile, tile_size);

  for (k = 2; k <= tile_size; k *= 2) {
    for (j = k / 2; j > 0; j /= 2) {
      for (i = 0; i < tile_size; i++) {
        l = i ^ j;
        if (l < i) {
          i = i + j - 1;
        }
        else if ((tile[i] > tile[l]) == !(i & k)) {
          sys_assert(i >= 0);
          sys_assert(i < tile_size);
          sys_assert(l >= 0);
          sys_assert(l < tile_size);
          tmp = tile[i];
          tile[i] = tile[l];
          tile[l] = tmp;
        }
      }
    }
  }
}


// ===========================================================================
// ===========================================================================
void bitonic_myrmics_region_copy(rid_t src_id, int **src_ptr, 
                                 rid_t dst_id, int **dst_ptr, 
                                 int region_size, int tile_size) {
  int i;
  int j;

  for (i = 0; i < region_size; i++) {
    for (j = 0; j < tile_size; j++) {
      dst_ptr[i][j] = src_ptr[i][j];
    }
  }
}


// ===========================================================================
// ===========================================================================
void bitonic_myrmics_region_sort(rid_t region_id, int **region_ptr, 
                                 int region_size, int tile_size) {
  int i;
  int *tile;

  for (i = 0; i < region_size; i++) {
    tile = region_ptr[i];

//printf("tile_sort tile ptr 0x%X\r\n", tile);

    #pragma myrmics task inout(tile) in(tile_size)
    bitonic_myrmics_tile_sort(tile, tile_size);
  }
}


// ===========================================================================
// ===========================================================================
void bitonic_myrmics_time(rid_t unused, unsigned int time_start) {

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
// ===========================================================================
void bitonic_myrmics_verify(rid_t r, int ***array, int num_regions,
                            int region_size, int num_tiles, int tile_size) {

  int val;
  int val_prev = -1;
  int i;
  int j;
  int k;

  printf("Sorting finished, verifying...\r\n");

  for (i = 0; i < num_regions; i++) {
    for (j = 0; j < region_size; j++) {
      for (k = 0; k < tile_size; k++) {

        if ((i != 0 || j != 0 || k != 0) && (array[i][j][k] < val_prev)) {
          printf("Verification [31;1mFAILED[0m: "
                 "array[%d][%d][%d] = %d, array[%d][%d][%d] = %d\r\n", 
                 i, 
                 j, 
                 k, 
                 array[i][k], 
                 (k > 0) ? i : (j > 0) ? j : i - 1,
                 (k > 0) ? j : (j > 0) ? j - 1 : j,
                 (k > 0) ? k - 1 : tile_size - 1,
                 val_prev);
          //sys_abort();
          return;
        }
        val_prev = array[i][j][k];
      }
    }
  }

  printf("Verification [32;1mPASSED[0m\r\n");
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
void bitonic_myrmics(int num_regions, int num_tiles, int num_elements, 
                     int verify, int dummy_regions) {

  int           region_size;
  int           tile_size;
  rid_t         r;
  rid_t         *array0_regions;
  rid_t         *array1_regions;
  int           ***array0;
  int           ***array1;
  int           **src1_ptr;
  int           **src2_ptr;
  int           **dst_ptr;
  rid_t         *tmp_regions;
  int           ***tmp;
  rid_t         src1_id;
  rid_t         src2_id;
  rid_t         dst_id;
  int           **tmp_ptr;
  rid_t         tmp_id;
  int           **region_ptr;
  rid_t         region_id;
  unsigned int  seed;
  int           mask;
  int           mask2;
  int           partner_i;
  int           partner_j;
  int           partner_ij;
  int           cur;
  unsigned int  time_start;
  int           i;
  int           j;
  int           k;
  int           ij;


  // Sanity checks
  if (num_regions & (num_regions - 1)) {
    printf("Number of regions should be a power of 2\r\n");
    return;
  }
  if (num_tiles & (num_tiles - 1)) {
    printf("Number of tiles should be a power of 2\r\n");
    return;
  }
  if (num_elements & (num_elements - 1)) {
    printf("Number of elements should be a power of 2\r\n");
    return;
  }

  if (num_tiles % num_regions) {
    printf("%d tiles not divisible by %d regions\r\n", 
           num_tiles, num_regions);
    return;
  }
  region_size = num_tiles / num_regions;
  if (region_size <= 0) {
    printf("Too few tiles\r\n");
    return;
  }

  if (num_elements % num_tiles) {
    printf("%d elements not divisible by %d tiles\r\n", 
           num_elements, num_tiles);
    return;
  }
  tile_size = num_elements / num_tiles;
  if (tile_size <= 0) {
    printf("Too few elements\r\n");
    return;
  }

  // Create all-holding region
  r = sys_ralloc(0, 99); // highest level
  
  // Create regions
  array0_regions = sys_alloc(num_regions * sizeof(rid_t), r);
  array1_regions = sys_alloc(num_regions * sizeof(rid_t), r);
  tmp_regions    = sys_alloc(num_regions * sizeof(rid_t), r);

  for (i = 0; i < num_regions; i++) {
    array0_regions[i] = sys_ralloc(r, 0); // lowest level
  }
  for (i = 0; i < dummy_regions; i++) {
    // Dummy regions to balance leaf schedulers as user wants
    sys_ralloc(r, 0); // lowest level
  }

  for (i = 0; i < num_regions; i++) {
    array1_regions[i] = sys_ralloc(r, 0); // lowest level
  }
  for (i = 0; i < dummy_regions; i++) {
    // Dummy regions to balance leaf schedulers as user wants
    sys_ralloc(r, 0); // lowest level
  }

  for (i = 0; i < num_regions; i++) {
    tmp_regions[i] = sys_ralloc(r, 0); // lowest level
  }


  // Allocate arrays
  array0 = sys_alloc(num_regions * sizeof(int **), r);
  array1 = sys_alloc(num_regions * sizeof(int **), r);
  tmp    = sys_alloc(num_regions * sizeof(int **), r);

  for (i = 0; i < num_regions; i++) {
    array0[i] = sys_alloc(region_size * sizeof(int *), array0_regions[i]);
    array1[i] = sys_alloc(region_size * sizeof(int *), array1_regions[i]);
    tmp[i]    = sys_alloc(region_size * sizeof(int *), tmp_regions[i]);

    if (region_size > 1) {
      sys_balloc(tile_size * sizeof(int), array0_regions[i], 
                 region_size, array0[i]);
      sys_balloc(tile_size * sizeof(int), array1_regions[i], 
                 region_size, array1[i]);
      sys_balloc(tile_size * sizeof(int), tmp_regions[i], 
                 region_size, tmp[i]);
    }
    else {
      array0[i][0] = sys_alloc(tile_size * sizeof(int), array0_regions[i]);
      array1[i][0] = sys_alloc(tile_size * sizeof(int), array1_regions[i]);
      tmp[i][0]    = sys_alloc(tile_size * sizeof(int), tmp_regions[i]);
    }
  }


  // Infomercial
  printf(
    "Bitonic sort of %d elements split into %d total tiles and %d regions\r\n",
    num_elements, num_tiles, num_regions);


  // Initialize first array with random values
  for (i = 0; i < num_regions; i++) {
    for (j = 0; j < region_size; j++) {
      seed = i * 42 + j * 99 + 666;
      for (k = 0; k < tile_size; k++) {
        array0[i][j][k] = (seed = rand(seed)) % 2147483647;
      }
    }
  }


  // Start time
  time_start = sys_free_timer_get_ticks();


  // Sort each tile in parallel before we begin the bitonic phases
  for (i = 0; i < num_regions; i++) {

      region_id = array0_regions[i];
      region_ptr = array0[i];

//printf("region_sort rid %d region ptr 0x%X\r\n", region_id, region_ptr);

      #pragma myrmics task region inout(region_id) \
                           in(region_ptr) safe(region_ptr) \
                           in(region_size, tile_size)
      bitonic_myrmics_region_sort(region_id, region_ptr, 
                                  region_size, tile_size);
  }
  

  // Do the parallel phases
  cur = 0;
  for (mask = 2; mask <= num_tiles; mask <<= 1) {

    for (mask2 = mask >> 1; mask2; mask2 >>= 1) {

      for (i = 0; i < num_regions; i++) {

        // Select regions & buffers
        partner_i = (mask2 < region_size) ? i : i ^ (mask2 / region_size);
        if (cur == 0) {
          src1_ptr = array0[i];
          src2_ptr = array0[partner_i];
          dst_ptr  = array1[i];

          src1_id  = array0_regions[i];
          src2_id  = array0_regions[partner_i];
          dst_id   = array1_regions[i];
        }
        else {
          src1_ptr = array1[i];
          src2_ptr = array1[partner_i];
          dst_ptr  = array0[i];

          src1_id  = array1_regions[i];
          src2_id  = array1_regions[partner_i];
          dst_id   = array0_regions[i];
        }
        tmp_ptr = tmp[i];
        tmp_id  = tmp_regions[i];

//printf("region_merge %d 0x%X + %d 0x%X -> %d 0x%X\r\n", src1_id, src1_ptr, src2_id, src2_ptr, dst_id, dst_ptr);

        if (src1_id != src2_id) {
          #pragma myrmics task region in(src2_id) in(src2_ptr) safe(src1_ptr) \
                               region inout(tmp_id) in(tmp_ptr) safe(tmp_ptr) \
                               in(region_size, tile_size)
          bitonic_myrmics_region_copy(src2_id, src2_ptr, tmp_id, tmp_ptr, 
                                      region_size, tile_size);

          #pragma myrmics task region in(src1_id) in(src1_ptr) safe(src1_ptr) \
                               region in(tmp_id) in(tmp_ptr) safe(tmp_ptr) \
                               region inout(dst_id) in(dst_ptr) safe(dst_ptr) \
                               in(i, partner_i, mask, mask2, region_size, \
                                  tile_size, verify)
          bitonic_myrmics_region_merge(src1_id, src1_ptr, tmp_id, tmp_ptr,
                                       dst_id, dst_ptr, i, partner_i, mask,
                                       mask2, region_size, tile_size, verify);
        }
        else {
          // We must avoid giving the same region id as a dependence twice,
          // the runtime does not check this and it will lead to double the
          // DMA (if it doesn't hang completely). We pass src2_id by value.
          #pragma myrmics task region in(src1_id) in(src1_ptr) safe(src1_ptr) \
                               in(src2_id) safe(src2_id) \
                               in(src2_ptr) safe(src1_ptr) \
                               region inout(dst_id) in(dst_ptr) safe(dst_ptr) \
                               in(i, partner_i, mask, mask2, region_size, \
                                  tile_size, verify)
          bitonic_myrmics_region_merge(src1_id, src1_ptr, src2_id, src2_ptr,
                                       dst_id, dst_ptr, i, partner_i, mask,
                                       mask2, region_size, tile_size, verify);
        }
      }
   
      cur = 1 - cur;
    }
  } 
  
  // Stop time
  #pragma myrmics task region inout(r) in(time_start)
  bitonic_myrmics_time(r, time_start);


  // Verify results
  if (verify == 1) {

    if (cur == 0) {
      #pragma myrmics task region in(r) in(array0) safe(array0) \
                           in(num_regions, region_size, num_tiles, tile_size)
      bitonic_myrmics_verify(r, array0, num_regions, region_size, 
                             num_tiles, tile_size);
    }
    else {
      #pragma myrmics task region in(r) in(array1) safe(array1) \
                           in(num_regions, region_size, num_tiles, tile_size)
      bitonic_myrmics_verify(r, array1, num_regions, region_size,
                             num_tiles, tile_size);
    }

  }

}
