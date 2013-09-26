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
// Abstract      : Matrix multiplication benchmark
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: matrix_mult_mpi.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/05/16 13:12:30 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>


// ===========================================================================
// ===========================================================================
static inline int tile_to_rank(int tile_size, int tile_row, int tile_col) {
  return (tile_row * tile_size + tile_col);
}

// ===========================================================================
// ===========================================================================
static inline void rank_to_tile(int tile_size, int rank, 
                                int *ret_tile_row, int *ret_tile_col) {
  ar_int_divide(rank, tile_size, ret_tile_row, ret_tile_col);
}

// ===========================================================================
// ===========================================================================
void matrix_mult_mpi_place_partial(float *whole, float *partial, 
                                   int tile_size, int blk_size, int rank) {
  int i;
  int j;
  int tile_row;
  int tile_col;
  int whole_row;
  int whole_col;

  rank_to_tile(tile_size, rank, &tile_row, &tile_col);

  for (i = 0; i < blk_size; i++) {
    for (j = 0; j < blk_size; j++) {
      whole_row = tile_row * blk_size + i;
      whole_col = tile_col * blk_size + j;

      whole[whole_row * tile_size * blk_size + whole_col] = 
        partial[i * blk_size + j];
    }
  }

}

// ===========================================================================
// ===========================================================================
void matrix_mult_mpi_print_matrix(float *matrix, int size) {
  int i, j;
  int integer;
  int comma;

  for (i = 0; i < size; i++) {
    for (j = 0; j < size; j++) {
      integer = matrix[i * size + j];
      comma = (matrix[i * size + j] - integer) * 1000;

      kt_printf("%4d.%01d ", integer, comma);
    }
    kt_printf("\r\n");
  }
}


// ===========================================================================
// matrix_mult_mpi()            MPI version of dense matrix multiplication.
//                              Processor cores are laid out in a 2D space
//                              core_rows * core_cols layout. Each core
//                              has a portion of matrices A and B (each such
//                              portion blk_size * blk_size numbers) and
//                              computes a portion of matrix C = A * B.
// ===========================================================================
// * INPUTS
//   int num_procs              Number of processors to use from the MPI setup
//   int tile_size              Number of core rows and columns in the 2D layout
//   int blk_size               Number of rows and columns for each portion.
//                              The total size of A, B and C is
//                              tile_size ^ 2 * blk_size ^ 2 elements.
//   int verify                 0: don't verify the multiplication
//                              1: core 0 gathers all matrices back and
//                                 does the verification
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int matrix_mult_mpi(int num_procs, int tile_size, int blk_size, int verify) {

  int                   num_cores;
  int                   rank;
  int                   tile_row;
  int                   tile_col;
  int                   blk_size_sq = 0;
  int                   whole_size = 0;
  int                   whole_size_sq = 0;
  int                   phase;
  int                   tag_a;
  int                   tag_b;
  int                   tag_c;
  MPI_Request           reqs[2];
  int                   num_reqs;
  MPI_Status            status;
  float                 *a = NULL;
  float                 *b = NULL;
  float                 *c = NULL;
  int                   peer_rank;
  float                 *peer_a = NULL;
  float                 *peer_b = NULL;
  float                 *work_a = NULL;
  float                 *work_b = NULL;
  float                 *whole_a = NULL;
  float                 *whole_b = NULL;
  float                 *whole_c = NULL;
  float                 *peer_c = NULL;
  float                 verify_val;
#ifdef ARCH_MB
  unsigned int          time_start = 0;
  unsigned int          time_stop;
  unsigned int          time;
#endif
  int                   i;
  int                   j;
  int                   k;


  // Initialize MPI
  //MPI_Init(NULL, NULL);
 
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
  if (tile_size * tile_size != num_procs) {
    if (!rank) {
      kt_printf("%d * %d tiles != %d cores\r\n", 
                tile_size, tile_size, num_procs);
    }
    return 1;
  }
  rank_to_tile(tile_size, rank, &tile_row, &tile_col);


  // Are we running?
  if (rank < num_procs) {

    // Create array portions
    blk_size_sq = blk_size * blk_size;
    whole_size = tile_size * blk_size;
    whole_size_sq = whole_size * whole_size;
    a = kt_malloc(blk_size_sq * sizeof(float));
    b = kt_malloc(blk_size_sq * sizeof(float));
    c = kt_malloc(blk_size_sq * sizeof(float));

    // Initialize our portions of A and B with some values, zero out C
    for (i = 0; i < blk_size_sq; i++) {
      a[i] = rank + i;
      b[i] = rank - i;
      c[i] = 0.0;
    }

    // Create two buffers to receive peer A and B portions
    peer_a = kt_malloc(blk_size_sq * sizeof(float));
    peer_b = kt_malloc(blk_size_sq * sizeof(float));

    // Assign tags for A and B send/recvs
    tag_a = 42;
    tag_b = 666;
    tag_c = 3;
  }

  // Synchronize everyone
  if (!rank) {
    kt_printf("Matrix multiplication of %d x %d starting on %d core(s)\r\n", 
              whole_size, whole_size, num_procs);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  // Keep time
  if (!rank) {
#ifdef ARCH_MB
    time_start = ar_glob_timer_read();
#endif
  }
  
  // If not part of active cores, skip to next barrier
  if (rank >= num_procs) {
    goto skip;
  }


  // For all the phases in the algorithm
  for (phase = 0; phase < tile_size; phase++) {

    // Remember how many waits we'll have to do
    num_reqs = 0;

    // Are we responsible to broadcast our A portion?
    if (tile_col == phase) {
      // Broadcast to others in our tile row
      for (i = 0; i < tile_size; i++) {
        if (i != tile_col) {
          peer_rank = tile_to_rank(tile_size, tile_row, i);
          MPI_Isend(a, blk_size_sq, MPI_FLOAT, peer_rank, tag_a, 
                    MPI_COMM_WORLD, NULL);
        }
      }
      work_a = a;
    }
    else {
      // Receive A portion from someone else from our tile row
      peer_rank = tile_to_rank(tile_size, tile_row, phase);
      MPI_Irecv(peer_a, blk_size_sq, MPI_FLOAT, peer_rank, tag_a, 
                MPI_COMM_WORLD, &reqs[num_reqs++]);
      work_a = peer_a;
    }


    // Are we responsible to broadcast our B portion?
    if (tile_row == phase) {
      // Broadcast to others in our tile col
      for (i = 0; i < tile_size; i++) {
        if (i != tile_row) {
          peer_rank = tile_to_rank(tile_size, i, tile_col);
          MPI_Isend(b, blk_size_sq, MPI_FLOAT, peer_rank, tag_b, 
                    MPI_COMM_WORLD, NULL);
        }
      }
      work_b = b;
    }
    else {
      // Receive B portion from someone else from our tile col
      peer_rank = tile_to_rank(tile_size, phase, tile_col);
      MPI_Irecv(peer_b, blk_size_sq, MPI_FLOAT, peer_rank, tag_b, 
                MPI_COMM_WORLD, &reqs[num_reqs++]);
      work_b = peer_b;
    }


    // Wait for needed receives to arrive
    if (num_reqs) {
      MPI_Waitall(num_reqs, reqs, &status);
    }

    // Add to partial results of C the A * B portion sums
    for (i = 0; i < blk_size; i++) {
      for (j = 0; j < blk_size; j++) {
        for (k = 0; k < blk_size; k++) {
          c[i * blk_size + j] += work_a[i * blk_size + k] * 
                                 work_b[k * blk_size + j];
        }
      }
    }
  }


  // Synchronize
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

  if (!verify) {
    goto finished;
  }

  // Rank 0 gathers all A, B and C's and verifies
  if (rank == 0) {

    kt_printf("Multiplication finished, rank 0 is gathering results...\r\n");

    // Allocate big arrays
    whole_a = kt_malloc(whole_size_sq * sizeof(float));
    whole_b = kt_malloc(whole_size_sq * sizeof(float));
    whole_c = kt_malloc(whole_size_sq * sizeof(float));
    
    // Allocate partial C buffer
    peer_c = kt_malloc(blk_size_sq * sizeof(float));

    // Place my partial arrays
    matrix_mult_mpi_place_partial(whole_a, (float *) a, tile_size, blk_size, 0);
    matrix_mult_mpi_place_partial(whole_b, (float *) b, tile_size, blk_size, 0);
    matrix_mult_mpi_place_partial(whole_c, (float *) c, tile_size, blk_size, 0);

    // Gather from others
    for (peer_rank = 1; peer_rank < num_procs; peer_rank++) {
      MPI_Recv(peer_a, blk_size_sq, MPI_FLOAT, peer_rank, tag_a, 
               MPI_COMM_WORLD, &status);
      MPI_Recv(peer_b, blk_size_sq, MPI_FLOAT, peer_rank, tag_b, 
               MPI_COMM_WORLD, &status);
      MPI_Recv(peer_c, blk_size_sq, MPI_FLOAT, peer_rank, tag_c, 
               MPI_COMM_WORLD, &status);

      matrix_mult_mpi_place_partial(whole_a, (float *) peer_a, tile_size, 
                                    blk_size, peer_rank);
      matrix_mult_mpi_place_partial(whole_b, (float *) peer_b, tile_size, 
                                    blk_size, peer_rank);
      matrix_mult_mpi_place_partial(whole_c, (float *) peer_c, tile_size, 
                                    blk_size, peer_rank);
    }

    // Print
    //kt_printf("A is:\r\n");
    //matrix_mult_mpi_print_matrix(whole_a, whole_size);
    //kt_printf("\r\nB is:\r\n");
    //matrix_mult_mpi_print_matrix(whole_b, whole_size);
    //kt_printf("\r\nC is:\r\n");
    //matrix_mult_mpi_print_matrix(whole_c, whole_size);
    //kt_printf("\r\n");

    // Verify
    for (i = 0; i < whole_size; i++) {
      for (j = 0; j < whole_size; j++) {
        verify_val = 0;
        for (k = 0; k < whole_size; k++) {
          verify_val += whole_a[i * whole_size + k] * 
                        whole_b[k * whole_size + j];
        }
        if (whole_c[i * whole_size + j] != verify_val) {
          kt_printf("Results gathered: "
                    "Verification [31;1mFAILED[0m at C[%d, %d]\r\n", i, j);
          while (1) {
            ;
          }
        }
      }
    }
    kt_printf("Results gathered. Verification [32;1mPASSED[0m\r\n");
  }
  else if (rank < num_procs) {
    // Send partial arrays to rank 0
    MPI_Send(a, blk_size_sq, MPI_FLOAT, 0, tag_a, MPI_COMM_WORLD);
    MPI_Send(b, blk_size_sq, MPI_FLOAT, 0, tag_b, MPI_COMM_WORLD);
    MPI_Send(c, blk_size_sq, MPI_FLOAT, 0, tag_c, MPI_COMM_WORLD);
  }

finished:

  // Free stuff
  if (rank < num_procs) {
    kt_free(a);
    kt_free(b);
    kt_free(c);
    kt_free(peer_a);
    kt_free(peer_b);
    if (verify) {
      kt_free(whole_a);
      kt_free(whole_b);
      kt_free(whole_c);
      kt_free(peer_c);
    }
  }

  return 0;
}
