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
// Abstract      : Parallel bitonic sort, MPI version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: bitonic_mpi.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/05/30 13:14:59 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>


// ===========================================================================
// ===========================================================================
static inline void merge_even(int size, int *src1, int *src2, int *dst) {

  int i;
  int j = 0;
  int k = 0;

  for (i = 0; i < size; i++) {
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
static inline void merge_odd(int size, int *src1, int *src2, int *dst) {

  int i;
  int j = size - 1;
  int k = size - 1;

  for (i = size - 1; i >= 0; i--) {
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
void bitonic_mpi_print(int rank, int num_cores, int *array, int my_elements) {
  int i;
  int r;

  for (r = 0; r < num_cores; r++) {
    if (r == rank) {
      kt_printf("rank = %d:\r\n", rank);
      for (i = 0; i < my_elements; i++) {
        kt_printf("%12d\r\n", array[i]);
      }
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}


// ===========================================================================
// ===========================================================================
void bitonic_mpi_verify(int rank, int num_cores, int *array, int my_elements) {

  int           max = -1;
  MPI_Status    status;
  int           j;
  int           r;


  if (!rank) {
    kt_printf("Sorting finished, rank 0 is verifying...\r\n");
  }

  for (r = 0; r < num_cores; r++) {
    if (!rank) {
      if (r) {
        MPI_Recv(array, my_elements, MPI_INT, r, 42 + r, 
                 MPI_COMM_WORLD, &status);
      }

      for (j = 0; j < my_elements; j++) {
        if (array[j] < max) {
          kt_printf("Verification [31;1mFAILED[0m: rank %d pos %d\r\n",
                    r, j);
          while (1) {
            ;
          }
        }
        else {
          max = array[j];
        }
      }
    }
    else if (r == rank) {
      MPI_Send(array, my_elements, MPI_INT, 0, 42 + rank, MPI_COMM_WORLD);
    }

    MPI_Barrier(MPI_COMM_WORLD);
  }

  if (!rank) {
    kt_printf("Results gathered. Verification [32;1mPASSED[0m\r\n");
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
int bitonic_mpi(int num_procs, int num_elements, int verify) {

  int           num_cores;
  int           rank;
  int           my_elements;
  int           *array = NULL;
  int           *buf = NULL;
  int           *merge_buf = NULL;
  unsigned int  seed;
  unsigned int  tmp;
  int           dim;
  int           mask;
  int           mask2;
  int           partner;
  MPI_Request   reqs[2];
  MPI_Status    status;
#ifdef ARCH_MB
  unsigned int  time_start = 0;
  unsigned int  time_stop;
  unsigned int  time;
#endif
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
  if (num_cores & (num_cores - 1)) {
    if (!rank) {
      kt_printf("Number of cores should be a power of 2\r\n");
    }
    return 1;
  }
  if (num_elements & (num_elements - 1)) {
    if (!rank) {
      kt_printf("Number of elements should be a power of 2\r\n");
    }
    return 1;
  }
  if (num_elements % num_procs) {
    if (!rank) {
      kt_printf("%d elements not divisible by %d cores\r\n",
                num_elements, num_procs);
    }
    return 1;
  }
  my_elements = num_elements / num_procs;
  if (my_elements <= 0) {
    if (!rank) {
      kt_printf("Too few elements\r\n");
    }
    return 1;
  }
  if (my_elements % 16) {
    // We need this so that array slices are multiples of 64 bytes (cache line)
    if (!rank) {
      kt_printf("Array slice not a multiple of 64-B\r\n");
    }
    return 1;
  }


  // Allocate buffers for our portion of the elements
  if (rank < num_procs) {
    array     = kt_malloc(my_elements * sizeof(int)); 
    buf       = kt_malloc(my_elements * sizeof(int)); 
    merge_buf = kt_malloc(my_elements * sizeof(int)); 
  }


  // Synchronize everyone and print infomercial
  MPI_Barrier(MPI_COMM_WORLD);
  if (!rank) {
    kt_printf("Bitonic sort of %d elements starting on %d core(s)\r\n",
              num_elements, num_procs);
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


  // Initialize our array with random values
  seed = rank * 42 + 666;
  for (i = 0; i < my_elements; i++) {
    array[i] = (seed = kt_rand(seed)) % 2147483647;
  }

  // Do a local sort first
  for (k = 2; k <= my_elements; k *= 2) {
    for (j = k / 2; j > 0; j /= 2) {
      for (i = 0; i < my_elements; i++) {
        l = i ^ j;
        if (l < i) {
          i = i + j - 1;
        }
        else if ((array[i] > array[l]) == !(i & k)) {
          tmp = array[i];
          array[i] = array[l];
          array[l] = tmp;
        }
      }
    }
  }

  // Do the parallel phases
  for (i = 2, mask = 2; i <= num_procs; i *= 2, mask <<= 1) {

    dim = kt_int_log2(i);
    mask2 = 1 << (dim - 1);

    for (j = 0; j < dim; j++, mask2 >>= 1) {
      partner = rank ^ mask2;

      // Exchange with partner
      MPI_Irecv(buf, my_elements, MPI_INT, partner, 0, 
                MPI_COMM_WORLD, &(reqs[0]));
      MPI_Isend(array, my_elements, MPI_INT, partner, 0, 
                MPI_COMM_WORLD, &(reqs[1]));
      MPI_Waitall(2, reqs, &status);

      // Merge
      if ((((rank & mask) == 0) && (rank < partner)) ||
          (((rank & mask) != 0) && (rank > partner))) {
        merge_even(my_elements, array, buf, merge_buf);
      }
      else {
        merge_odd(my_elements, array, buf, merge_buf);
      }
      for (k = 0; k < my_elements; k++) {
        array[k] = merge_buf[k];
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


  // Verify results
  if (verify) {
    if (rank < num_procs) {
      bitonic_mpi_verify(rank, num_procs, array, my_elements);
    }
    else {
      bitonic_mpi_verify(rank, num_procs, NULL, 0);
    }
  }

  //if (rank < num_procs) {
  //  bitonic_mpi_print(rank, num_procs, array, my_elements);
  //}
  //else {
  //  bitonic_mpi_print(rank, num_procs, NULL, 0);
  //}

  // Free stuff
  if (rank < num_procs) {
    kt_free(array);
    kt_free(buf);
    kt_free(merge_buf);
  }


  return 0;
}
