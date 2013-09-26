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
// Abstract      : 2D Jacobi kernel, row-wise distribution, bare-metal version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: jacobi_bare.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2012/05/16 13:12:30 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>


// ===========================================================================
// ===========================================================================
static unsigned int buf_adr(int cid, int buf_id, int my_rows, int num_cols) {
  unsigned int buf_adr;

  buf_adr = mm_va_kernel_base(cid);
  if (buf_id == 0) {
    return buf_adr;
  }
  else if (buf_id == 1) {
    return buf_adr + (my_rows + 2) * num_cols * sizeof(float);
  }
  else {
    ar_abort();
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
int jacobi_bare(int num_procs, int num_rows, int num_cols, int num_iter) {

  int           rank;
  int           my_bid;
  int           my_cid;
  int           dst_bid;
  int           dst_cid;
  unsigned int  dst_buf;
  int           my_rows = 0;
  float         *grid[2];
  int           rep;
  int           start_row;
  int           end_row;
  int           cur_buf;
  int           nxt_buf;
#ifdef ARCH_MB
  unsigned int  time_start = 0;
  unsigned int  time_stop;
  unsigned int  time;
#endif
  int           i;
  int           j;
  int           k;


  // Who are we?
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  rank = my_bid * 8 + my_cid;

  // Sanity checks
  if (num_rows % num_procs) {
    if (!rank) {
      kt_printf("%d rows not divisible by %d cores\r\n",
                num_rows, num_procs);
    }
    return 1;
  }
  if (num_cols % 16) {
    // We need this so that a single row is a multiple of 64 bytes (cache line)
    if (!rank) {
      kt_printf("Number of columns must be a multiple of 16\r\n");
    }
    return 1;
  }
  ar_assert(!((num_cols * sizeof(float)) % 64));


  // Allocate buffers for our portion of the table
  if (rank < num_procs) {
    my_rows = num_rows / num_procs;
    grid[0] = (float *) buf_adr(my_cid, 0, my_rows, num_cols);
    grid[1] = (float *) buf_adr(my_cid, 1, my_rows, num_cols);


    // Initialize grid and extra row buffers to 0.0, grid borders to 1.0
    for (i = 0; i < my_rows + 2; i++) {
      for (j = 0; j < num_cols; j++) {
        if (((i == 1) && (rank == 0)) || 
            ((i == my_rows) && (rank == num_procs - 1)) ||
            (j == 0) || 
            (j == num_cols - 1)) {
          grid[0][i * num_cols + j] = 1.0;
          grid[1][i * num_cols + j] = 1.0;
        }
        else {
          grid[0][i * num_cols + j] = 0.0;
          grid[1][i * num_cols + j] = 0.0;
        }
      }
    }
  }


  // Synchronize everyone and print infomercial
  ar_cnt_set(my_cid, 0, -num_cols * sizeof(float));
  ar_cnt_set(my_cid, 1, -num_cols * sizeof(float));
  ar_cnt_set(my_cid, 2, -num_cols * sizeof(float));
  ar_cnt_set(my_cid, 3, -num_cols * sizeof(float));
  if (!rank) {
    for (i = 0; i < 511; i++) {
      ar_assert(ar_mbox_get(my_cid) == 0x1000BABE);
    }
    kt_printf("Jacobi 2D-row of %d x %d starting on %d core(s)\r\n",
              num_rows, num_cols, num_procs);
    for (i = 1; i < 512; i++) {
      ar_mbox_send(my_cid, i / 8, i % 8, 0x2000BABE);
    }
  }
  else {
    ar_mbox_send(my_cid, 0, 0, 0x1000BABE);
    ar_assert(ar_mbox_get(my_cid) == 0x2000BABE);
  }

  // Keep time
  if (!rank) {
#ifdef ARCH_MB
    time_start = ar_glob_timer_read();
#endif
  }

  // This kernel is reentrant for multiple runs. If our core is not
  // part of the current num_procs setup, go to the next barrier directly.
  if (rank >= num_procs) {
    goto skip;
  }


  cur_buf = 0;
  nxt_buf = 1;
  start_row = (rank == 0) ? 2 : 1;
  end_row   = (rank == num_procs - 1) ? my_rows - 1: my_rows;


  for (rep = 0; rep < num_iter; rep++) {

    // Wait until borders arrive for the current buffer
    if (rep > 0) {
      if (rank > 0) {
        while (ar_cnt_get(my_cid, 2 * cur_buf + 1)) {
          ;
        }
      }
      if (rank < num_procs - 1) {
        while (ar_cnt_get(my_cid, 2 * cur_buf + 0)) {
          ;
        }
      }
      ar_cnt_set(my_cid, 2 * cur_buf + 0, -num_cols * sizeof(float));
      ar_cnt_set(my_cid, 2 * cur_buf + 1, -num_cols * sizeof(float));
    }

    // Compute
    for (i = start_row; i <= end_row; i++) {
      k = i * num_cols;
      for (j = 1; j < num_cols - 1; j++) {
        grid[nxt_buf][k + j] = (
                                grid[cur_buf][k - num_cols + j] + // north 
                                grid[cur_buf][k + num_cols + j] + // south 
                                grid[cur_buf][k + j - 1] +        // west
                                grid[cur_buf][k + j + 1]          // east
                               ) * 0.25;
      }
    }


    // Send our top row to our previous rank
    if (rank > 0) {
      dst_bid = (rank - 1) / 8;
      dst_cid = (rank - 1) % 8;
      dst_buf = buf_adr(dst_cid, nxt_buf, my_rows, num_cols) + 
                (my_rows + 1) * num_cols * sizeof(float);

      ar_dma_with_ack(my_cid,
                      my_bid,  my_cid,  (int) &(grid[nxt_buf][1*num_cols]),
                      dst_bid, dst_cid, dst_buf,
                      dst_bid, dst_cid, 2 * nxt_buf + 0,
                      num_cols * sizeof(float), 0, 0, 0);
    }
    // Send our bottom row to our next rank
    if (rank < num_procs - 1) {
      dst_bid = (rank + 1) / 8;
      dst_cid = (rank + 1) % 8;
      dst_buf = buf_adr(dst_cid, nxt_buf, my_rows, num_cols) + 
                0 * num_cols;

      ar_dma_with_ack(my_cid,
                      my_bid,  my_cid, (int) &(grid[nxt_buf][my_rows*num_cols]),
                      dst_bid, dst_cid, dst_buf,
                      dst_bid, dst_cid, 2 * nxt_buf + 1,
                      num_cols * sizeof(float), 0, 0, 0);
    }

    // Switch buffers
    cur_buf = 1 - cur_buf;
    nxt_buf = 1 - nxt_buf;
  }



skip:

  // Barrier
  if (!rank) {
    for (i = 0; i < 511; i++) {
      ar_assert(ar_mbox_get(my_cid) == 0x3000BABE);
    }
#ifdef ARCH_MB
    time_stop = ar_glob_timer_read();
    if (time_stop > time_start) {
      time = time_stop - time_start;
    }
    else {
      time = 0xFFFFFFFF - (time_start - time_stop);
    }
    kt_printf("Time: %10u cycles (%6u msec)\r\n", time, time / 10000);
#endif
    for (i = 1; i < 512; i++) {
      ar_mbox_send(my_cid, i / 8, i % 8, 0x4000BABE);
    }
  }
  else {
    ar_mbox_send(my_cid, 0, 0, 0x3000BABE);
    ar_assert(ar_mbox_get(my_cid) == 0x4000BABE);
  }

  return 0;
}
