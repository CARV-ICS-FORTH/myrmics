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
// Abstract      : MPI communication tests
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: test_mpi.c,v $
// CVS revision  : $Revision: 1.4 $
// Last modified : $Date: 2012/06/05 09:13:27 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <fmpi.h>

#define TEST_PASS "[32;1mPASS[0m"
#define TEST_FAIL "[31;1mFAIL[0m"


// ===========================================================================
// ===========================================================================
void test_mpi_master_multi_receive(int rank, int num_cores, int master) {

  volatile int          *buf;
  int                   peer_rank;
  MPI_Status            st;
  int                   i;
  int                   j;


  if (rank == master) {
    kt_printf("%d: MPI master multi-receive starts\r\n", rank);
  }

  // Allocate communication buffer
  buf = kt_malloc(16 * sizeof(int));

  // Everybody sends to the master
  if (rank != master) {
    for (j = 0; j < 2; j++) {
      for (i = 0; i < 16; i++) {
        buf[i] = rank + i + j;
      }

      MPI_Send(buf, 16, MPI_INT, master, 0, MPI_COMM_WORLD);
    }
  }
  else {
    // Master gathers
    for (peer_rank = 0; peer_rank < num_cores; peer_rank++) {
      if (peer_rank == master) {
        continue;
      }
      for (j = 0; j < 2; j++) {
        MPI_Recv(buf, 16, MPI_INT, peer_rank, MPI_ANY_TAG, MPI_COMM_WORLD, &st);

        // Verify
        for (i = 0; i < 16; i++) {
          if (buf[i] != peer_rank + i + j) {
            kt_printf("%d: peer_rank = %d, buf[%d] = %d [ %s ]\r\n", 
                      rank, peer_rank, i, buf[i], TEST_FAIL);
            while (1) {
              ;
            }
          }
        }
      }
    }
    kt_printf("%d: Received 64 + 64 bytes from %d slaves [ %s ]\r\n", 
              rank, num_cores - 1, TEST_PASS);
  }

  // Free buffer
  kt_free((void *) buf);
}


// ===========================================================================
// ===========================================================================
void test_mpi_pipeline(int rank, int num_cores, int packet_size, 
                       int num_packets) {

  volatile int  *buf;
  MPI_Status    st;
  int           i;
  int           j;


  // Sanity checks
  ar_assert(num_cores > 2);

  // Allocate buffer
  buf = kt_malloc(packet_size);

  // Producer
  if (rank == 0) {

    kt_printf("%d: MPI pipeline starts\r\n", rank);

    for (i = 0; i < num_packets; i++) {
      for (j = 0; j < packet_size / 4; j++) {
        buf[j] = rank + i + j;
      }
      MPI_Send(buf, packet_size / 4, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
  }

  // Middleman
  else if (rank < num_cores - 1) {

    for (i = 0; i < num_packets; i++) {

      MPI_Recv(buf, packet_size / 4, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &st);

      for (j = 0; j < packet_size / 4; j++) {
        buf[j] += rank;
      }

      MPI_Send(buf, packet_size / 4, MPI_INT, rank + 1, 0, MPI_COMM_WORLD);
    }
  }

  // Consumer
  else {

    for (i = 0; i < num_packets; i++) {

      MPI_Recv(buf, packet_size / 4, MPI_INT, rank - 1, 0, MPI_COMM_WORLD, &st);

      for (j = 0; j < packet_size / 4; j++) {

        if (buf[j] != (num_cores - 1) * (num_cores - 2) / 2 + i + j) {
          kt_printf("%d: rep %d, buf[%d] = %d [ %s ]\r\n", 
                    rank, i, j, buf[j], TEST_FAIL);
          while (1) {
            ;
          }
        }
      }
    }
    kt_printf("\n%d: Received %d * %d bytes [ %s ]\r\n", 
              rank, num_packets, packet_size, TEST_PASS);
  }

  // Free buffer
  kt_free((void *) buf);


}


// ===========================================================================
// ===========================================================================
void test_mpi_barrier(int rank, int num_cores) {

  int           i;
  unsigned int  time[10];


  for (i = 0; i < 10; i++) {
    ar_timer_reset();
    MPI_Barrier(MPI_COMM_WORLD);
    time[i] = ar_timer_get_cycles();
  }

  if (!rank) {
    kt_printf("Barrier time = %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d\r\n"
              "               %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }
}


// ===========================================================================
// ===========================================================================
void test_mpi_sendrecv(int rank, int num_cores, int blocking) {

  int           i;
  unsigned int  time[20];
  int           *buf;
  MPI_Status    st;
  MPI_Request   req;


  MPI_Barrier(MPI_COMM_WORLD);

  buf = kt_malloc(64);
  for (i = 0; i < 16; i++) {
    buf[i] = 0;
  }

  // Sender
  if (rank == 0) {
    if (blocking) {
      for (i = 0; i < 10; i++) {
        ar_timer_reset();
        MPI_Send(buf, 16, MPI_INT, 1, 0, MPI_COMM_WORLD);
        time[i] = ar_timer_get_cycles();
      }
    }
    else {
      for (i = 0; i < 20; i++) {
        ar_timer_reset();
        MPI_Isend(buf, 16, MPI_INT, 1, 0, MPI_COMM_WORLD, &req);
        time[i++] = ar_timer_get_cycles();
        ar_timer_reset();
        MPI_Wait(&req, &st);
        time[i] = ar_timer_get_cycles();
      }
    }
  }

  // Receiver
  else if (rank == 1) {
    if (blocking) {
      for (i = 0; i < 10; i++) {
        ar_timer_reset();
        MPI_Recv(buf, 16, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &st);
        time[i] = ar_timer_get_cycles();
      }
    }
    else {
      for (i = 0; i < 20; i++) {
        ar_timer_reset();
        MPI_Irecv(buf, 16, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &req);
        time[i++] = ar_timer_get_cycles();
        ar_timer_reset();
        MPI_Wait(&req, &st);
        time[i] = ar_timer_get_cycles();
      }
    }
  }

  kt_free(buf);


  MPI_Barrier(MPI_COMM_WORLD);

  // Sender print results
  if (rank == 0) {
    if (blocking) {
      kt_printf("MPI_Send time = %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d cycles\r\n",
                time[0], time[1], time[2], time[3], time[4],
                time[5], time[6], time[7], time[8], time[9]);
    }
    else {
      kt_printf("MPI_Isend time = %12d, wait time = %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d cycles\r\n",
                time[0], time[1], time[2], time[3], time[4],
                time[5], time[6], time[7], time[8], time[9],
                time[10], time[11], time[12], time[13], time[14],
                time[15], time[16], time[17], time[18], time[19]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);

  // Receiver print results
  if (rank == 1) {
    if (blocking) {
      kt_printf("MPI_Recv time = %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d\r\n"
                "                %12d cycles\r\n",
                time[0], time[1], time[2], time[3], time[4],
                time[5], time[6], time[7], time[8], time[9]);
    }
    else {
      kt_printf("MPI_Irecv time = %12d, wait time = %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d\r\n"
                "                 %12d,             %12d cycles\r\n",
                time[0], time[1], time[2], time[3], time[4],
                time[5], time[6], time[7], time[8], time[9],
                time[10], time[11], time[12], time[13], time[14],
                time[15], time[16], time[17], time[18], time[19]);
    }
  }

  MPI_Barrier(MPI_COMM_WORLD);
}


// ===========================================================================
// ===========================================================================
void test_mpi_broadcast(int rank, int num_cores) {

  int           i;
  unsigned int  time[10];
  int           *buf;


  buf = kt_malloc(64);
  for (i = 0; i < 16; i++) {
    buf[i] = 0;
  }

  for (i = 0; i < 10; i++) {
    ar_timer_reset();
    MPI_Bcast(buf, 16, MPI_INT, 0, MPI_COMM_WORLD);
    time[i] = ar_timer_get_cycles();
  }

  kt_free(buf);

  if (!rank) {
    kt_printf("Broadcast time = %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }
}


// ===========================================================================
// ===========================================================================
void test_mpi_reduction(int rank, int num_cores) {

  int           i;
  unsigned int  time[10];
  int           *buf_in;
  int           *buf_out;


  buf_in = kt_malloc(64);
  buf_out = kt_malloc(64);
  for (i = 0; i < 16; i++) {
    buf_in[i] = 0;
    buf_out[i] = 0;
  }

  for (i = 0; i < 10; i++) {
    ar_timer_reset();
    MPI_Reduce(buf_out, buf_in, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    time[i] = ar_timer_get_cycles();
  }

  kt_free(buf_in);
  kt_free(buf_out);

  if (!rank) {
    kt_printf("Reduction time = %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d\r\n"
              "                 %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }
}


// ===========================================================================
// ===========================================================================
void test_mpi_alltoall(int rank, int num_cores) {

  int           i;
  unsigned int  time[10];
  int           *buf_in;
  int           *buf_out;


  buf_in = kt_malloc(512 * 4);
  buf_out = kt_malloc(512 * 4);
  for (i = 0; i < 512 * 1; i++) {
    buf_in[i] = 0;
    buf_out[i] = 0;
  }

  for (i = 0; i < 10; i++) {
    if (!rank) {
      kt_printf("entering %d\r\n", i);
    }
    MPI_Barrier(MPI_COMM_WORLD);
    ar_timer_reset();
    MPI_Alltoall(buf_out, 1, MPI_INT, 
                 buf_in,  1, MPI_INT, 
                 MPI_COMM_WORLD);
    time[i] = ar_timer_get_cycles();
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank) {
      kt_printf("done %d\r\n", i);
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }

  kt_free(buf_in);
  kt_free(buf_out);

  if (!rank) {
    kt_printf("Alltoall time = %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d\r\n"
              "                %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }
}


// ===========================================================================
// test_mpi()                   Various MPI communication test patterns
// ===========================================================================
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int test_mpi() {

  int           num_cores;
  int           rank;
  //int           i;


  // Initialize MPI
  //MPI_Init(NULL, NULL);
 
  // Who are we?
  MPI_Comm_size(MPI_COMM_WORLD, &num_cores);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // Multi-receive
#if 0
  for (i = 0; i < 4; i++) {
    test_mpi_master_multi_receive(rank, num_cores, i);
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank) {
      kt_printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif


  // Pipeline
#if 0
  for (i = 64; i <= 1024 * 1024; i *= 2) {
    test_mpi_pipeline(rank, num_cores, i, 42);
    MPI_Barrier(MPI_COMM_WORLD);
    if (!rank) {
      kt_printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#endif

  // Barriers
#if 0
  test_mpi_barrier(rank, num_cores);
#endif

  // Blocking send-recv primitives
#if 0
  test_mpi_sendrecv(rank, num_cores, 1);
#endif
 
  // Non-blocking send-recv primitives
#if 0
  test_mpi_sendrecv(rank, num_cores, 0);
#endif
 
  // Broadcasts
#if 0
  test_mpi_broadcast(rank, num_cores);
#endif

  // Reductions
#if 0
  test_mpi_reduction(rank, num_cores);
#endif

  // Alltoall
#if 1
  test_mpi_alltoall(rank, num_cores);
#endif


  // End
  //MPI_Finalize();

  return 0;
}
