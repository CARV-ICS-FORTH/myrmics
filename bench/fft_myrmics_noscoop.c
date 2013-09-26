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
// Abstract      : 2D-FFT kernel, Myrmics version. Adapted from the MPI 
//                 version.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: fft_myrmics_noscoop.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/12/14 16:06:46 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <myrmics.h>

#define SIGN      1.0
#define MAX_TILES 512


// ===========================================================================
// ===========================================================================
void fft_myrmics_checksum(rid_t r, float ***a, int tile_size, int blk_size) {
  int i, j, k;
  float checksum;

  checksum = 0.0F;
  for (i = 0; i < tile_size; i++) {
    for (j = 0; j < tile_size; j++) {
      for (k = 0; k < 2 * blk_size * blk_size; k++) {
        checksum += a[i][j][k];
      }
    }
  }

  printf("Checksum is %f\r\n", checksum);
}


// ===========================================================================
// ===========================================================================
void fft_myrmics_step(int n, int mj, float *a, float *b, float *c, 
                      float *d, float *w) {

  int j, k, jc, jw, lj, mj2;
  float rp, up;
  float wr[4], wu[4];
  float V0[4], V1[4], V2[4], V3[4], V4[4], V6[4], V7[4];

  mj2 = 2 * mj;
  lj = n / mj2;
  for (j = 0; j < lj; j++) {
    jw = j * mj;
    jc = j * mj2;
    rp = w[2 * jw + 0];
    up = w[2 * jw + 1];
    if (SIGN < 0.0)
      up = -up;
    if (mj < 2) {               /* special case mj=1 */
      d[2*jc+0] = rp * (a[2*jw+0] - b[2*jw+0]) - up * (a[2*jw+1] - b[2*jw+1]);
      d[2*jc+1] = up * (a[2*jw+0] - b[2*jw+0]) + rp * (a[2*jw+1] - b[2*jw+1]);
      c[2*jc+0] = a[2*jw+0] + b[2*jw+0];
      c[2*jc+1] = a[2*jw+1] + b[2*jw+1];
    }
    else {                      /* mj > 1 cases */
      wr[0] = rp;
      wr[1] = rp;
      wr[2] = rp;
      wr[3] = rp;

      wu[0] = -up;
      wu[1] = up;
      wu[2] = -up;
      wu[3] = up;

      //V6 = _mm_load_ps (wr);
      V6[0] = wr[0];
      V6[1] = wr[1];
      V6[2] = wr[2];
      V6[3] = wr[3];

      //V7 = _mm_load_ps (wu);
      V7[0] = wu[0];
      V7[1] = wu[1];
      V7[2] = wu[2];
      V7[3] = wu[3];

      for (k = 0; k < mj; k += 2) {
        //V0 = _mm_load_ps (&a[jw + k][0]);
        V0[0] = a[2 * jw + k + 0];
        V0[1] = a[2 * jw + k + 1];
        V0[2] = a[2 * jw + k + 2];
        V0[3] = a[2 * jw + k + 3];

        //V1 = _mm_load_ps (&b[jw + k][0]);
        V1[0] =b[2 * jw + k + 0];
        V1[1] =b[2 * jw + k + 1];
        V1[2] =b[2 * jw + k + 2];
        V1[3] =b[2 * jw + k + 3];

        //V2 = _mm_add_ps (V0, V1); /* a+b */
        V2[0] = V0[0] + V1[0];
        V2[1] = V0[1] + V1[1];
        V2[2] = V0[2] + V1[2];
        V2[3] = V0[3] + V1[3];

        //_mm_store_ps (&c[jc + k][0], V2); /* c to M */
        c[2 * jc + k + 0] = V2[0];
        c[2 * jc + k + 1] = V2[1];
        c[2 * jc + k + 2] = V2[2];
        c[2 * jc + k + 3] = V2[3];

        //V3 = _mm_sub_ps (V0, V1); /* a-b */
        V3[0] = V0[0] - V1[0];
        V3[1] = V0[1] - V1[1];
        V3[2] = V0[2] - V1[2];
        V3[3] = V0[3] - V1[3];

        //V4 = _mm_shuffle_ps (V3, V3, _MM_SHUFFLE (2, 3, 0, 1));
        V4[0] = V3[1];
        V4[1] = V3[0];
        V4[2] = V3[3];
        V4[3] = V3[2];

        //V0 = _mm_mul_ps (V6, V3);
        V0[0] = V6[0] * V3[0];
        V0[1] = V6[1] * V3[1];
        V0[2] = V6[2] * V3[2];
        V0[3] = V6[3] * V3[3];

        //V1 = _mm_mul_ps (V7, V4);
        V1[0] = V7[0] * V4[0];
        V1[1] = V7[1] * V4[1];
        V1[2] = V7[2] * V4[2];
        V1[3] = V7[3] * V4[3];

        //V2 = _mm_add_ps (V0, V1); /* w*(a-b) */
        V2[0] = V0[0] + V1[0];
        V2[1] = V0[1] + V1[1];
        V2[2] = V0[2] + V1[2];
        V2[3] = V0[3] + V1[3];

        //_mm_store_ps (&d[jc + k][0], V2); /* d to M */
        d[2 * jc + k + 0] = V2[0];
        d[2 * jc + k + 1] = V2[1];
        d[2 * jc + k + 2] = V2[2];
        d[2 * jc + k + 3] = V2[3];
      }
    }
  }
}


// ===========================================================================
// x=in, y=out, w=exp(2*pi*i*k/n), k=0..n/2-1
// ===========================================================================
void fft_myrmics_cfft2(int n, float *x, float *y, float *w) {

  int m, j, mj, tgle, i;

  m = int_log2(n); //(int) (log ((float) n) / log (1.99));
  mj = 1;
  tgle = 1;

  fft_myrmics_step(n, mj, x, x + n, y, y + 2 * mj, w);

  for (j = 0; j < m - 2; j++) {
    mj *= 2;
    if (tgle) {
      fft_myrmics_step(n, mj, y, y + n, x, x + 2 * mj, w);
      tgle = 0;
    }
    else {
      fft_myrmics_step(n, mj, x, x + n, y, y + 2 * mj, w);
      tgle = 1;
    }
  }
  if (tgle) {
    for (i = 0; i < n; i++) {
      y[i] = x[i];
    }
  }
  
  mj = n / 2;
  
  fft_myrmics_step(n, mj, x, x + n, y, y + 2 * mj, w);
}


// ===========================================================================
// ===========================================================================
void fft_myrmics_cfft2_row(rid_t r_single_row, float **a, float *w, int n, 
                           int tile_size, int blk_size) {
  int i, j, k;
  int off;
  float *pa;

//printf("%d: row begin\r\n", sys_get_worker_id());

  // tmp buffer in order to hold the data of a single row
  pa = sys_alloc(2 * n * sizeof(float), r_single_row);

  for (i = 0; i < blk_size; i++) {

    off = 2 * i * blk_size;

    // merge data from all blocks "a" into a single table "pa"
    for (j = 0; j < tile_size; j++) {
      for (k = 0; k < 2*blk_size; k++) {
        pa[j*2*blk_size+k] = a[j][k+off]; 
      }
    }

    fft_myrmics_cfft2(n, pa, pa, w);

    // write results back from "pa" to block "a"
    for (j = 0; j < tile_size; j++) {
      for (k = 0; k < 2*blk_size; k++) {
        a[j][k+off] = pa[j*2*blk_size+k];
      }
    }
  }

  sys_free(pa);

printf("%d: row end\r\n", sys_get_worker_id());
}


// ===========================================================================
// ===========================================================================
void fft_myrmics_cfft2_phase(rid_t *buf_regions, float ***buf, float **w, 
                             int n, int tile_size, int blk_size) {

  rid_t         buf_single_region;
  float         **buf_row;
  float         *w_row;
  int           i;
  void          *args[6];
  unsigned int  deps[6];


  // For each row of tiles, run cfft2
  for (i = 0; i < tile_size; i++) {

    buf_single_region = buf_regions[i];
    buf_row           = buf[i];
    w_row             = w[i]; 

//printf("cfft2_row rid %d ptrs 0x%X 0x%X\r\n", buf_single_region, buf_row, w_row);

    args[0] = (void *) buf_single_region;
    deps[0] = SYS_TYPE_REGION_ARG | SYS_TYPE_INOUT_ARG;
    args[1] = (void *) buf_row;
    deps[1] = SYS_TYPE_BYVALUE_ARG;
    args[2] = (void *) w_row;
    deps[2] = SYS_TYPE_IN_ARG;
    args[3] = (void *) n;
    deps[3] = SYS_TYPE_BYVALUE_ARG;
    args[4] = (void *) tile_size;
    deps[4] = SYS_TYPE_BYVALUE_ARG;
    args[5] = (void *) blk_size;
    deps[5] = SYS_TYPE_BYVALUE_ARG;
    sys_spawn(1, args, deps, 6); // fft_myrmics_cfft2_row()
  }
}


// ===========================================================================
// ===========================================================================
void fft_myrmics_Xpose_blk(float *src, float *dst, int blk_size) {

  int i;
  int j;

  for (i = 0; i < blk_size; i++) {
    for (j = 0; j < blk_size; j++) {
      dst[i * 2 * blk_size + j * 2]     = src[j * 2 * blk_size + 2 * i];
      dst[i * 2 * blk_size + j * 2 + 1] = src[j * 2 * blk_size + 2 * i + 1];
    }
  }
}


// ===========================================================================
// ===========================================================================
void fft_myrmics_Xpose_row(rid_t unused, float **dst_buf_row, int tile_size,
                           int blk_size, ...) {
  va_list       ap;
  float         *src_buf;
  int           i;

//printf("%d: xpose begin\r\n", sys_get_worker_id());

  va_start(ap, blk_size);

  for (i = 0; i < tile_size; i++) {
    src_buf = va_arg(ap, float *);
    fft_myrmics_Xpose_blk(src_buf, dst_buf_row[i], blk_size);
  }

  va_end(ap);

printf("%d: xpose end\r\n", sys_get_worker_id());
}


// ===========================================================================
// ===========================================================================
void fft_myrmics_Xpose(float ***src_buf_copy, rid_t *dst_buf_row, 
                       float ***dst_buf, int tile_size, int blk_size) {

  void          *args[MAX_TILES + 4];
  unsigned int  deps[MAX_TILES + 4];
  int           i;
  int           j;


  sys_assert(tile_size <= MAX_TILES);
  
  for (i = 0; i < tile_size; i++) {
    
    args[0] = (void *) dst_buf_row[i];
    deps[0] = SYS_TYPE_REGION_ARG | SYS_TYPE_INOUT_ARG;

    args[1] = (void *) dst_buf[i];
    deps[1] = SYS_TYPE_BYVALUE_ARG;

    args[2] = (void *) tile_size;
    deps[2] = SYS_TYPE_BYVALUE_ARG;

    args[3] = (void *) blk_size;
    deps[3] = SYS_TYPE_BYVALUE_ARG;

    for (j = 0; j < tile_size; j++) {

      args[j + 4] = (void *) src_buf_copy[j][i];
      deps[j + 4] = SYS_TYPE_IN_ARG;

    }

    sys_spawn(2, args, deps, tile_size + 4); // fft_myrmics_Xpose_row
  }

}


// ===========================================================================
// ===========================================================================
void fft_myrmics_time(rid_t unused, unsigned int time_start) {

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
// fft_myrmics()                Myrmics version of FFT
//                              A n x n input table is partitioned into
//                              tile_size x tile_size sub-blocks where each
//                              block has size (n/tile_size) x (n/tile_size)
// ===========================================================================
// * INPUTS
//   int n                      Number of rows and columns of the input table
//   int tile_size              Number of tile rows and columns in the 2D layout
// ===========================================================================
void fft_myrmics(int tile_size, int n) {

  rid_t                 r;
  rid_t                 *ra_row;
  rid_t                 *rb_row;
  int                   blk_size;
  int                   blk_size_sq;
  int                   i, j, k, l;
  float                 ***a;
  float                 ***a_copy;
  float                 ***b;
  float                 ***b_copy;
  float                 **w;
  unsigned int          time_start;
  void                  *args[4];
  unsigned int          deps[4];



  // Sanity checks
  if (n & (n - 1)) {
    printf("N must be a power of 2\r\n");
    return;
  }
  if (tile_size & (tile_size - 1)) {
    printf("Tile size must be a power of 2\r\n");
    return;
  }


  blk_size = n/tile_size;
  blk_size_sq = blk_size * blk_size;

  // Create holding region
  r = sys_ralloc(0, 99); // highest level
  a = sys_alloc(tile_size * sizeof(float **), r);
  b = sys_alloc(tile_size * sizeof(float **), r);
  w = sys_alloc(tile_size * sizeof(float *), r);

  // Initialize tiles and regions
  ra_row = sys_alloc(tile_size * sizeof(rid_t), r);
  rb_row = sys_alloc(tile_size * sizeof(rid_t), r);

  for (i = 0; i < tile_size; i++) {

    // Create a region for each row of tiles, for buffers a and b
    ra_row[i] = sys_ralloc(r, 0);
    rb_row[i] = sys_ralloc(r, 0);

    // Allocate pointers for each tile in the row
    a[i] = sys_alloc(tile_size * sizeof(float *), ra_row[i]);
    b[i] = sys_alloc(tile_size * sizeof(float *), rb_row[i]);

    // Allocate tiles in the row
    sys_balloc(2* blk_size_sq * sizeof(float), ra_row[i], tile_size, 
               (void *) a[i]);
    sys_balloc(2* blk_size_sq * sizeof(float), rb_row[i], tile_size, 
               (void *) b[i]);

    // Allocate w array
    w[i] = sys_alloc(2 * n * sizeof(float), r);

    // Initialize tables with some values
    for (j = 0; j  < tile_size; j++) {
      for (k = 0; k < blk_size; k++) {
        for (l = 0; l < 2 * blk_size; l++) {
          a[i][j][k * 2 * blk_size + l] = 0.01F;
        }
      }
    }
    for (j = 0; j < 2 * n; j++) {
      w[i][j] = 0.3F;
    }
  }

  // Copy the pointers of the a and b tables, because we'll need them in
  // fft_myrmics_Xpose() to spawn the tasks. We need to do this, because in
  // fft_myrmics_cfft2_phase() we delegate all r_row[*] to children tasks:
  // a[*][*] and b[*][*] are allocated in these regions, so we don't have
  // access there anymore from the master task.
  a_copy = sys_alloc(tile_size * sizeof(float **), 0);
  b_copy = sys_alloc(tile_size * sizeof(float **), 0);

  sys_balloc(tile_size * sizeof(float *), 0, tile_size, (void *) a_copy);
  sys_balloc(tile_size * sizeof(float *), 0, tile_size, (void *) b_copy);

  for (i = 0; i < tile_size; i++) {
    for (j = 0; j < tile_size; j++) {
      a_copy[i][j] = a[i][j];
      b_copy[i][j] = b[i][j];
    }
  }


  // Starting FFT
  printf("FFT 2D-block of %d x %d starting split into %d x %d tiles \r\n",
              n, n, tile_size, tile_size);

  // Start time
  time_start = ar_free_timer_get_ticks();


  // Run first phase on buffer a
  fft_myrmics_cfft2_phase(ra_row, a, w, n, tile_size, blk_size);

  // Transpose a->b
  fft_myrmics_Xpose(a_copy, rb_row, b, tile_size, blk_size);

  // Run second phase on buffer b
  fft_myrmics_cfft2_phase(rb_row, b, w, n, tile_size, blk_size);

  // Transpose b->a
  fft_myrmics_Xpose(b_copy, ra_row, a, tile_size, blk_size);



  // Stop time
  args[0] = (void *) r;
  deps[0] = SYS_TYPE_REGION_ARG | SYS_TYPE_INOUT_ARG;
  args[1] = (void *) time_start;
  deps[1] = SYS_TYPE_BYVALUE_ARG;
  sys_spawn(3, args, deps, 2); // fft_myrmics_time()
  
  // Checksum
  args[0] = (void *) r;
  deps[0] = SYS_TYPE_REGION_ARG | SYS_TYPE_INOUT_ARG;
  args[1] = (void *) a;
  deps[1] = SYS_TYPE_BYVALUE_ARG;
  args[2] = (void *) tile_size;
  deps[2] = SYS_TYPE_BYVALUE_ARG;
  args[3] = (void *) blk_size;
  deps[3] = SYS_TYPE_BYVALUE_ARG;
  sys_spawn(4, args, deps, 4); // fft_myrmics_checksum()

printf("%d: spawns done\r\n", sys_get_worker_id());
}


// ===========================================================================
// ===========================================================================
void (*fft_myrmics_task_table[])() = {
           fft_myrmics,                  // 0
           fft_myrmics_cfft2_row,        // 1
  (func_t) fft_myrmics_Xpose_row,        // 2
           fft_myrmics_time,             // 3
           fft_myrmics_checksum          // 4
};
