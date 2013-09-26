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
// Abstract      : Smith-Waterman kernel, MPI version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: swat_mpi.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/05/17 15:37:37 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>


// ===========================================================================
// ===========================================================================
int swat_find_max(int a, int b, int c) {
  int ret;

  ret = (a > b)   ? a : b;
  ret = (c > ret) ? c : ret;

  return ret;
}


// ===========================================================================
// ===========================================================================
void swat_mpi_print(int rank, int num_cores, int *seq_matrix, 
                    int num_stripes, int rows_per_stripe, int seq1_len) {
  int i;
  int j;
  int k;
  int r;

  for (i = 0; i < num_stripes; i++) {
    for (r = 0; r < num_cores; r++) {
      if (r == rank) {
        kt_printf("stripe = %d, rank = %d:\r\n", i, rank);
        for (k = 0; k < rows_per_stripe; k++) {
          for (j = 0; j < seq1_len; j++) {
            kt_printf("%3d ", 
                seq_matrix[(i * rows_per_stripe + k) * seq1_len + j]);
          }
          kt_printf("\r\n");
        }
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
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
int swat_mpi(int num_procs, int seq0_len, int seq1_len, int seq1_step,
             int min_rows_per_stripe) {

  int           num_cores;
  int           rank;
  int           my_rows;
  int           my_stripes;
  int           rows_per_stripe;
  int           col_steps;
  char          base[] = {'a', 'c', 'g', 't'};
  unsigned int  seed;
  char          *seq0 = NULL;
  char          *seq1 = NULL;
  int           *seq_matrix = NULL;
  int           *buf = NULL;
#ifdef ARCH_MB
  unsigned int  time_start = 0;
  unsigned int  time_stop;
  unsigned int  time;
#endif
  int           row_idx;
  int           col_idx;
  int           idx;
  MPI_Status    status;
  MPI_Request   req;
  MPI_Request   unused;
  int           i;
  int           j;
  int           k;
  int           l;


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
  if (seq0_len % num_procs) {
    if (!rank) {
      kt_printf("Sequence0 length %d not divisible by %d cores\r\n",
                seq0_len, num_procs);
    }
    return 1;
  }
  if (seq0_len & (seq0_len - 1)) {
    if (!rank) {
      kt_printf("Sequence0 length must be a power of 2\r\n");
    }
    return 1;
  }
  if (seq1_len % seq1_step) {
    if (!rank) {
      kt_printf("Sequence1 length %d not divisible by step %d\r\n",
                seq1_len, seq1_step);
    }
    return 1;
  }
  if (seq1_step % 16) {
    // We need this so that a step communication is a multiple of 64 bytes 
    // (cache line size)
    if (!rank) {
      kt_printf("Sequence1 step must be a multiple of 16\r\n");
    }
    return 1;
  }

  // Compute row splits
  my_rows = seq0_len / num_procs;
  for (my_stripes = 1, rows_per_stripe = my_rows; 
       rows_per_stripe / 2 >= min_rows_per_stripe;
       my_stripes *= 2, rows_per_stripe /= 2) {
    ;
  }
  if (!rank) {
    kt_printf("my_rows = %d, my_stripes = %d, rows_per_stripe = %d\r\n", my_rows, my_stripes, rows_per_stripe);
    ar_assert(my_stripes * rows_per_stripe == my_rows);
    ar_assert(my_rows * num_procs == seq0_len);
  }

  // Compute column splits
  col_steps = seq1_len / seq1_step;


  // Allocate buffers
  if (rank < num_procs) {
    seq_matrix = kt_malloc(2 * rows_per_stripe * seq1_len * sizeof(int));
    buf = kt_malloc(seq1_len * sizeof(int));
    seq0 = kt_malloc(seq0_len * sizeof(char));
    seq1 = kt_malloc(seq1_len * sizeof(char));

    // Everybody initializes all parts of both sequences
    seed = 0;
    for (i = 0; i < seq0_len; i++) {
      seq0[i] = base[(seed = kt_rand(seed)) % 4];
    }
    for (i = 0; i < seq1_len; i++) {
      seq1[i] = base[(seed = kt_rand(seed)) % 4];
    }

    // Initialize our 2 stripes 
    for (i = 0; i < 2 * rows_per_stripe * seq1_len; i++) {
      seq_matrix[i] = 0;
    }
  }


  // Synchronize everyone and print infomercial
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    kt_printf("Smith-Waterman of %d x %d starting on %d core(s)\r\n",
              seq0_len, seq1_len, num_procs);
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

int snd_msg = 0;
int rcv_msg = 0;
int first_time = 1;

  for (i = 0; i < my_stripes; i++) {
    for (j = 0; j < col_steps; j++) {

      // Do the first col_step of the stripe seperately, so we can use the
      // incoming buffer. Rank 0 must not do that for its first stripe 
      // (i == 0).
      if (rank || i) {

        if (num_procs > 1) {

          // Pipeline preamble
          if (first_time) {
            //ar_assert(j * seq1_step + seq1_step <= seq1_len);
            MPI_Irecv(buf + j * seq1_step, seq1_step, MPI_INT,
                      (num_procs + rank - 1) % num_procs, 
                      rcv_msg, MPI_COMM_WORLD, &req);
            rcv_msg = (rcv_msg + 1) % 256;
            first_time = 0;
          }

          // Wait for the previous recv to arrive
          MPI_Wait(&req, &status);

          // Post a new recv for the next part
          if ((i < my_stripes - 1) || (j < col_steps - 1)) {
            //ar_assert(((j + 1) % col_steps) * seq1_step + seq1_step <= seq1_len);
            MPI_Irecv(buf + ((j + 1) % col_steps) * seq1_step, seq1_step, 
                      MPI_INT, (num_procs + rank - 1) % num_procs, 
                      rcv_msg, MPI_COMM_WORLD, &req);
            rcv_msg = (rcv_msg + 1) % 256;
          }
        }

        k = 0;
        for (l = (!j ? 1 : 0); l < seq1_step; l++) {
          
          // Compute
          row_idx = (i % 2) * rows_per_stripe;
          col_idx = j * seq1_step + l;
          idx     = row_idx * seq1_len + col_idx;

          ar_assert(idx >= 0);
          ar_assert(idx < 2 * rows_per_stripe * seq1_len);
          ar_assert(col_idx > 0);
          ar_assert(col_idx < seq1_len);

          seq_matrix[idx] = swat_find_max(
                                buf[col_idx - 1],       // northwest
                                buf[col_idx],           // north
                                seq_matrix[idx - 1]);   // west

          row_idx = i * rows_per_stripe * num_procs + 
                    rank * rows_per_stripe;

          ar_assert(row_idx > 0);
          ar_assert(col_idx > 0);
          ar_assert(row_idx < seq1_len);
          ar_assert(col_idx < seq1_len);

          if (seq0[row_idx - 1] == seq1[col_idx - 1]) {
            seq_matrix[idx] += 2;
          }
          else {
            seq_matrix[idx] -= 1;
          }

        }
      }

//int dummy;
//for (dummy = 0; dummy < 8; dummy++) {
      // Do the rest of the stripe lines
      for (k = 1; k < rows_per_stripe; k++) {
        for (l = (!j ? 1 : 0); l < seq1_step; l++) {

          // Compute
          row_idx = (i % 2) * rows_per_stripe + k;
          col_idx = j * seq1_step + l;
          idx     = row_idx * seq1_len + col_idx;

          ar_assert(idx - seq1_len - 1 >= 0);
          ar_assert(idx - 1 < 2 * rows_per_stripe * seq1_len);

          if ((i % 2) == 0) {
            seq_matrix[idx] = swat_find_max(
                                  seq_matrix[idx + seq1_len - 1], // northwest
                                  seq_matrix[idx + seq1_len],     // north
                                  seq_matrix[idx - 1]);           // west
          }
          else {
            seq_matrix[idx] = swat_find_max(
                                  seq_matrix[idx - seq1_len - 1], // northwest
                                  seq_matrix[idx - seq1_len],     // north
                                  seq_matrix[idx - 1]);           // west
          }

          row_idx = i * rows_per_stripe * num_procs + 
                    rank * rows_per_stripe +
                    k;
          ar_assert(row_idx > 0);
          ar_assert(col_idx > 0);
          ar_assert(row_idx < seq1_len);
          ar_assert(col_idx < seq1_len);

          if (seq0[row_idx - 1] == seq1[col_idx - 1]) {
            seq_matrix[idx] += 2;
          }
          else {
            seq_matrix[idx] -= 1;
          }

        }
      }
//}

      // Send bottom line of this col_step to the next peer
      if ((rank < num_procs - 1) || (i < my_stripes - 1)) {
        if (num_procs > 1) {
//kt_printf("%d: send i=%d, j=%d\r\n", rank, i, j);
          //ar_assert((i * rows_per_stripe + rows_per_stripe - 1) * seq1_len +
          //          j * seq1_step + seq1_step <= my_rows * seq1_len);
          MPI_Isend(seq_matrix + 
                     ((i % 2) * rows_per_stripe + rows_per_stripe - 1) * 
                        seq1_len + j * seq1_step,
                   seq1_step, MPI_INT, (rank + 1) % num_procs, 
                   snd_msg,
                   MPI_COMM_WORLD, &unused);
          snd_msg = (snd_msg + 1) % 256;
        }
        else {
          for (l = 0; l < seq1_step; l++) {
            buf[j * seq1_step + l] = 
              seq_matrix[((i % 2) * rows_per_stripe + 
                                rows_per_stripe - 1) * seq1_len +
                         j * seq1_step + l];
          }
        }
      }

    }
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
  //  swat_mpi_print(rank, num_procs, seq_matrix, my_stripes, rows_per_stripe, 
  //                 seq1_len);
  //}
  //else {
  //  swat_mpi_print(rank, num_procs, NULL, my_stripes, 0, 0);
  //}

  // Free stuff
  if (rank < num_procs) {
    kt_free(seq0);
    kt_free(seq1);
    kt_free(seq_matrix);
    kt_free(buf);
  }


  return 0;
}
