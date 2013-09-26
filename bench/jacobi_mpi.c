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
// Abstract      : 2D Jacobi kernel, row-wise distribution, MPI version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: jacobi_mpi.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2012/05/16 13:12:30 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>

#define TAG_TOP    100
#define TAG_BOTTOM 200

// ===========================================================================
// ===========================================================================
void jacobi_mpi_print(int rank, int num_cores, float *grid, int my_rows, 
                      int num_cols) {
  int i;
  int j;
  int r;

  for (r = 0; r < num_cores; r++) {
    if (r == rank) {
      kt_printf("rank = %d:\r\n", rank);
      for (i = 1; i <= my_rows; i++) {
        for (j = 0; j < num_cols; j++) {
          kt_printf("%.1f ", grid[i * num_cols + j]);
        }
        kt_printf("\r\n");
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
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
int jacobi_mpi(int num_procs, int num_rows, int num_cols, int num_iter) {

  int           num_cores;
  int           rank;
  int           my_rows = 0;
  float         *grid[2];
  int           rep;
  int           start_row;
  int           end_row;
  MPI_Request   reqs[2][2];
  MPI_Request   unused[2][2];
  MPI_Status    status[2][2];
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
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // Sanity checks
  if (num_cores < num_procs) {
    if (!rank) {
      kt_printf("Cannot run with %d cores, MPI setup has only %d cores\r\n",
                num_procs, num_cores);
    }
    return 1;
  }
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
    grid[0] = kt_malloc((my_rows + 2) * num_cols * sizeof(float));
    grid[1] = kt_malloc((my_rows + 2) * num_cols * sizeof(float));


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
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    kt_printf("Jacobi 2D-row of %d x %d starting on %d core(s)\r\n",
              num_rows, num_cols, num_procs);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Keep time
  if (!rank) {
#ifdef ARCH_MB
    time_start = ar_glob_timer_read();
#endif
  }

  // This kernel is reentrant for multiple MPI runs. If our core is not
  // part of the current num_procs setup, go to the next barrier directly.
  if (rank >= num_procs) {
    goto skip;
  }


  cur_buf = 0;
  nxt_buf = 1;
  start_row = (rank == 0) ? 2 : 1;
  end_row   = (rank == num_procs - 1) ? my_rows - 1: my_rows;


  for (rep = 0; rep < num_iter; rep++) {

    // Receive top row from our next rank and put it as our bottom border
    if (rank < num_procs - 1) {
      MPI_Irecv(&(grid[nxt_buf][(my_rows + 1) * num_cols]), num_cols, MPI_FLOAT,
                rank + 1, TAG_TOP + nxt_buf, MPI_COMM_WORLD, 
                &(reqs[nxt_buf][0]));
    }
    // Receive bottom row from our previous rank and put it as our top border
    if (rank > 0) {
      MPI_Irecv(&(grid[nxt_buf][0 * num_cols]), num_cols, MPI_FLOAT,
                rank - 1, TAG_BOTTOM + nxt_buf, MPI_COMM_WORLD, 
                &(reqs[nxt_buf][1]));
    }

    // Wait until borders arrive for the current buffer
    if (rep > 0) {
      if ((rank > 0) && (rank < num_procs - 1)) {
        MPI_Waitall(2, &(reqs[cur_buf][0]), &(status[cur_buf][0]));
      }
      else if (rank > 0) {
        MPI_Waitall(1, &(reqs[cur_buf][1]), &(status[cur_buf][1]));
      }
      else if (rank < num_procs - 1) {
        MPI_Waitall(1, &(reqs[cur_buf][0]), &(status[cur_buf][0]));
      }
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
      MPI_Isend(&(grid[nxt_buf][1 * num_cols]), num_cols, MPI_FLOAT,
                rank - 1, TAG_TOP + nxt_buf, MPI_COMM_WORLD, 
                &(unused[nxt_buf][0]));
    }
    // Send our bottom row to our next rank
    if (rank < num_procs - 1) {
      MPI_Isend(&(grid[nxt_buf][my_rows * num_cols]), num_cols, MPI_FLOAT,
                rank + 1, TAG_BOTTOM + nxt_buf, MPI_COMM_WORLD, 
                &(unused[nxt_buf][1]));
    }

    // Switch buffers
    cur_buf = 1 - cur_buf;
    nxt_buf = 1 - nxt_buf;
  }



skip:

  MPI_Barrier(MPI_COMM_WORLD);

  // Keep time
  if (!rank) {
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
  }

  //if (rank < num_procs) {
  //  jacobi_mpi_print(rank, num_procs, grid[cur_buf], my_rows, num_cols);
  //}
  //else {
  //  jacobi_mpi_print(rank, num_procs, NULL, 0, 0);
  //}

  // Free stuff
  if (rank < num_procs) {
    kt_free(grid[0]);
    kt_free(grid[1]);
  }


  return 0;
}
