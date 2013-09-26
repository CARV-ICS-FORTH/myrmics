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
// Abstract      : Minimal MPI library: exported MPI functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: mpi.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2012/10/24 13:06:55 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <stdarg.h>

#include <arch.h>
#include <kernel_toolset.h>
#include <memory_management.h>
#include <noc.h>
#include <fmpi.h>


// ==========================================================================
// ==========================================================================
int MPI_Init(int *argc, char ***argv) {

  Context       *context;
  int           cid;
  int           i;


  // Get context
  cid = ar_get_core_id();
  context = mm_get_context(cid);

  // Sanity checks
  ar_assert(sizeof(FMPI_Descriptor) == FMPI_DMA_WORD);

  if (argc || argv) {
    kt_printf("MPI_Init: argc/argv ignored\r\n");
    return MPI_ERR_OTHER;
  }

  // Initialize context fields
  for (i = 0; i < FMPI_NUM_USER_REQUESTS; i++) {
    context->fmpi->user_req[i].status = 0; // free
  }
  context->fmpi->user_req_free = fmpi_list_init(
                                        &context->fmpi->user_req_list[0], 
                                        FMPI_NUM_USER_REQUESTS);
  context->fmpi->user_req_used_first = NULL;
  context->fmpi->user_req_used_last = NULL;

  context->fmpi->cnt_list = kt_malloc(NOC_MAX_COUNTERS * sizeof(FMPI_List));
  context->fmpi->cnt_free = fmpi_list_init(&context->fmpi->cnt_list[0], 
                                           NOC_MAX_COUNTERS);

  context->fmpi->cnt_used_first = NULL;
  context->fmpi->cnt_used_last = NULL;

  context->fmpi->descr_free = fmpi_list_init(&context->fmpi->descr_list[0], 
                                             FMPI_NUM_DESCRIPTORS);
  context->fmpi->descr_used_first = NULL;
  context->fmpi->descr_used_last = NULL;

  context->fmpi->descr_recv_free = fmpi_list_init(
                                        &context->fmpi->descr_recv_list[0], 
                                        FMPI_NUM_DESCR_RECV);
  context->fmpi->descr_recv_used_first = NULL;
  context->fmpi->descr_recv_used_last = NULL;

  context->fmpi->incoming_dmas = ar_ni_status_get(cid) >> 16;

  context->fmpi->mbox_credits = kt_malloc(context->fmpi->num_ranks * 
                                          sizeof(unsigned char));
  for (i = 0; i < context->fmpi->num_ranks; i++) {
    context->fmpi->mbox_credits[i] = FMPI_CREDITS_PER_CORE;
  }

  // create tree of context (required for barrier and bcast)
  fmpi_find_father(context->fmpi);
  fmpi_find_children(context->fmpi);

  // reserve two counters for barrier
  context->fmpi->cnt_breq = context->fmpi->cnt_free->id; // use counter 0 
                                                         // for requests
  context->fmpi->cnt_free = context->fmpi->cnt_free->next;

  context->fmpi->cnt_back = context->fmpi->cnt_free->id; // use counter 1 
                                                         // for acks
  context->fmpi->cnt_free = context->fmpi->cnt_free->next;

  fmpi_barrier_init(context->fmpi);

  context->fmpi->mbx_pending = 0;

  // synchronize everyone
  fmpi_mbx_barrier(context->fmpi);

#ifdef FMPI_INTR_ENABLE
  // enable mailbox interrupts
  context->fmpi->intr_mask = 0xfd;
  context->fmpi->intr_cnt_pos = 1;
  // clear mailbox and counter interrupts
  ar_intr_cpu_set(0x1e00ff);
  // enable interrupts
  fmpi_intr_enable(context->fmpi);
#endif

#ifdef ARCH_MB
  // remember global time. Used in MPI_Wtime
  context->fmpi->time_init = ar_glob_timer_read();
  context->fmpi->time_last = 0;
  context->fmpi->time_ovfl = 0;
#endif

  return MPI_SUCCESS;
}


// ==========================================================================
// ==========================================================================
int MPI_Comm_size(MPI_Comm comm, int *size) {

  Context  *context;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // We only support the world communicator
  if (comm != MPI_COMM_WORLD) {
    return MPI_ERR_COMM;
  }

  // Check that return argument is sane
  if (!size) {
    return MPI_ERR_ARG;
  }

  // Return number of running cores
  *size = context->fmpi->num_ranks;

  // Success
  return MPI_SUCCESS;
}


// ==========================================================================
// ==========================================================================
int MPI_Comm_rank(MPI_Comm comm, int *rank) {

  Context  *context;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // We only support the world communicator
  if (comm != MPI_COMM_WORLD) {
    return MPI_ERR_COMM;
  }

  // Check that return argument is sane
  if (!rank) {
    return MPI_ERR_ARG;
  }

  // Return rank
  *rank = context->fmpi->rank;

  // Success
  return MPI_SUCCESS;
}


// ==========================================================================
// ==========================================================================
int MPI_Finalize() {
  
  // Check for remaining events ?

  // Block forever ?

  return MPI_SUCCESS;
}


// ==========================================================================
// ==========================================================================
int MPI_Send(volatile void *buf, int count, MPI_Datatype datatype, int dest, 
             int tag, MPI_Comm comm) {

  // commit command
  return fmpi_commit_cmd(buf, count, datatype, dest, tag, comm, 
                         FMPI_OPCODE_SEND, NULL, NULL);
}

// ==========================================================================
// ==========================================================================
int MPI_Isend(volatile void *buf, int count, MPI_Datatype datatype, int dest, 
             int tag, MPI_Comm comm, MPI_Request *request) {

  // commit command
  return fmpi_commit_cmd(buf, count, datatype, dest, tag, comm, 
                         FMPI_OPCODE_ISEND, NULL, request);
}

// ==========================================================================
// ==========================================================================
int MPI_Recv(volatile void *buf, int count, MPI_Datatype datatype, int source,
             int tag, MPI_Comm comm, MPI_Status *status) {

  // commit command
  return fmpi_commit_cmd(buf, count, datatype, source, tag, comm, 
                         FMPI_OPCODE_RECV, status, NULL);
}

// ==========================================================================
// ==========================================================================
int MPI_Irecv(volatile void *buf, int count, MPI_Datatype datatype, int source,
             int tag, MPI_Comm comm, MPI_Request *request) {

  // commit command
  return fmpi_commit_cmd(buf, count, datatype, source, tag, comm, 
                         FMPI_OPCODE_IRECV, NULL, request);
}

// ==========================================================================
// ==========================================================================
int MPI_Waitall(int count, MPI_Request *reqs, MPI_Status *status) {

  Context       *context;
  int           i;
  int           wait;

  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // disable interrupts
  fmpi_intr_disable(context->fmpi);

  // wait until all status are done
  wait = 1;
  while (wait) {
    wait = 0;
    for (i = 0; i < count; i++) {
      if (reqs[i].status != FMPI_REQ_DONE) {
        wait = 1;
      }
    }
    if (wait) {
      fmpi_process_pending_events(context->fmpi);
    }
  }

  // enable interrupts
  fmpi_intr_enable(context->fmpi);

  return MPI_SUCCESS;
}


// ==========================================================================
// ==========================================================================
int MPI_Barrier(MPI_Comm comm) {

  Context       *context;
  int           cid;
  int           i;

  
  // Get context
  cid = ar_get_core_id();
  context = mm_get_context(cid);

  // disable interrupts
  fmpi_intr_disable(context->fmpi);

  // wait until all children have sent a barrier request
  while (ar_cnt_get(cid, context->fmpi->cnt_breq) != 0) {
    fmpi_process_pending_events(context->fmpi);
  }

  // if not root
  if (context->fmpi->father != -1) {

    // send barrier request to father
    fmpi_barrier_send(context->fmpi, context->fmpi->father, 
                      context->fmpi->cnt_breq);

    // wait for barrier ack
    while (ar_cnt_get(cid, context->fmpi->cnt_back) != 0) {
      fmpi_process_pending_events(context->fmpi);
    }
  }

  // init for next barrier
  fmpi_barrier_init(context->fmpi);

  // send barrier acks
  for (i = 0; i < context->fmpi->num_children; i++) {
    fmpi_barrier_send(context->fmpi, context->fmpi->children[i], 
                      context->fmpi->cnt_back);
  }
  
  // enable interrupts
  fmpi_intr_enable(context->fmpi);

  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int MPI_Bcast(volatile void *buf, int count, MPI_Datatype datatype, 
              int root, MPI_Comm comm) {

  MPI_Status    status;
  Context       *context;
  MPI_Request   reqs[24]; // total send requests should be less than 24
  int           i;
  int           req_check_results;
  int           rank_recv;
  int           total_reqs;
 

  // Get context
  context = mm_get_context(ar_get_core_id());

  // disable interrupts
  fmpi_intr_disable(context->fmpi);

  // make sure command is ok
  req_check_results = fmpi_req_check(context->fmpi, buf, count, datatype, 
                                     root, 0, comm);
  if (req_check_results != MPI_SUCCESS) {
    ar_abort();
    return req_check_results;
  }

  // receive data from root       
  rank_recv = fmpi_bcast_recv(context->fmpi, root);;
  if (rank_recv != -1) {
    MPI_Recv(buf, count, datatype, rank_recv, MPI_TAG_BCAST, comm, &status);
  }

  // send bcast to father
  total_reqs = 0;
  if (context->fmpi->father != -1 && context->fmpi->father != rank_recv) {
    MPI_Isend(buf, count, datatype, context->fmpi->father, MPI_TAG_BCAST, 
              comm, &reqs[total_reqs++]);
  }

  // send bcast to children
  for (i = 0; i < context->fmpi->num_children; i++) {
    if (rank_recv != context->fmpi->children[i]) {
      MPI_Isend(buf, count, datatype, context->fmpi->children[i], MPI_TAG_BCAST, 
                comm, &reqs[total_reqs++]);
    }
  }
  ar_assert(total_reqs < 24);

  // wait for bcast to finish 
  // (otherwise we may have problem with back-to-back bcast)
  MPI_Waitall(total_reqs, reqs, &status);

  // enable interrupts
  fmpi_intr_enable(context->fmpi);

  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int MPI_Wait(MPI_Request *request, MPI_Status *status) {

  FMPI_Context       *context;

  // Get context
  context = mm_get_context(ar_get_core_id())->fmpi;

  // disable interrupts
  fmpi_intr_disable(context);

  // wait until all status are done
  while (request->status != FMPI_REQ_DONE)
    fmpi_process_pending_events(context);

  // enable interrupts
  fmpi_intr_enable(context);

  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int MPI_Reduce(volatile void *sendbuf, volatile void *recvbuf, int count,
               MPI_Datatype datatype, MPI_Op op, int root, MPI_Comm comm) {

  FMPI_Context  *context;
  MPI_Status status;
  volatile void *gather; // gather data from others
  volatile void *result; // results of reduce operation
  int rank_send;
  int i;

  // Get context
  context = mm_get_context(ar_get_core_id())->fmpi;

  // disable interrupts
  fmpi_intr_disable(context);

  // allocate buffers
  gather = kt_malloc(fmpi_sizeof(datatype, count));
  if (context->rank == root) {
    result = recvbuf;
  } else {
    result = kt_malloc(fmpi_sizeof(datatype, count));
  }

  // init result buffer
  fmpi_reduce(sendbuf, result, count, datatype, FMPI_INIT);

  // reverse operation from broadcast
  // find the path in order to send the results to root
  rank_send = fmpi_bcast_recv(context, root);

  // receive data from parent
  if (context->father != -1 && context->father != rank_send) {
    MPI_Recv(gather, count, datatype, context->father, MPI_TAG_REDUCE,
              comm, &status);

    // fix results
    fmpi_reduce(gather, result, count, datatype, op);
  }

  // receive data from children
  for (i = 0; i < context->num_children; i++) {
    if (rank_send != context->children[i]) {
      MPI_Recv(gather, count, datatype, context->children[i], MPI_TAG_REDUCE,
                comm, &status);
      // fix results
      fmpi_reduce(gather, result, count, datatype, op);
    }
  }

  // send results to root path
  if (rank_send != -1) {
    MPI_Isend(result, count, datatype, rank_send, MPI_TAG_REDUCE, comm, NULL);
  }

  // make sure that all data is sent.
  // (otherwise we may have problem with back-to-back operations)
  MPI_Barrier(comm);

  kt_free((void *)gather);
  if (context->rank != root) kt_free((void *)result);

  // enable interrupts
  fmpi_intr_enable(context);


  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int MPI_Allreduce(volatile void *sendbuf, volatile void *recvbuf, int count,
                  MPI_Datatype datatype, MPI_Op op, MPI_Comm comm) {

  MPI_Reduce(sendbuf, recvbuf, count, datatype, op, 0, comm);
  MPI_Bcast(recvbuf, count, datatype, 0, comm);

  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int MPI_Alltoallv(volatile void *sendbuf, int *sendcounts, int *sdispls,
                  MPI_Datatype sendtype, volatile void *recvbuf, 
                  int *recvcounts, int *rdispls, MPI_Datatype recvtype, 
                  MPI_Comm comm) {

  FMPI_Context  *context;
  int rank;
  MPI_Request reqs[512+FMPI_ALL2ALL_GS];  // maximum 512 cores
  int total_reqs;
  MPI_Status status;
  void *buf_tmp[FMPI_ALL2ALL_GS+1];


  // Get context
  context = mm_get_context(ar_get_core_id())->fmpi;
  int max_rank = context->num_ranks - 1;;

  // disable interrupts
  fmpi_intr_disable(context);

  // if this is Alltoall command, we support only MPI_COMM_WORLD
  if (!sdispls) {
    ar_assert(*sendcounts == *recvcounts);
    ar_assert(sendtype == recvtype);
  }

  ar_assert(max_rank < 512);
  int type_size = fmpi_sizeof(sendtype, 1);

  // In order not to overfllow resources serve 64 ranks at a time
  int rank_start, rank_end;
  int i, j;

  if (sdispls) ar_assert(rdispls);

  int rbuf_size = FMPI_ALL2ALL_RBUF_SZ;
  if (sdispls) buf_tmp[FMPI_ALL2ALL_GS] = kt_malloc(fmpi_sizeof(recvtype, rbuf_size));
  else buf_tmp[FMPI_ALL2ALL_GS] = kt_malloc(fmpi_sizeof(recvtype, *recvcounts));

  for (i = 0; i < (max_rank/FMPI_ALL2ALL_GS)+1; i++) {
    // find group of ranks where all ranks will send data to
    rank_start = FMPI_ALL2ALL_GS*i;
    if (max_rank < FMPI_ALL2ALL_GS*(i+1)-1) rank_end = max_rank;
    else rank_end = FMPI_ALL2ALL_GS*(i+1)-1;

    total_reqs = 0;
    // scatter/send data to the group
    for (rank = rank_start; rank <= rank_end; rank++) {
      // copy data to alligned space before mpi_isend
      if (sdispls) {
        if (sendcounts[rank]) {
          buf_tmp[total_reqs] = kt_malloc(fmpi_sizeof(sendtype, sendcounts[rank]));
          fmpi_reduce(sendbuf+(sdispls[rank]*type_size), buf_tmp[total_reqs], sendcounts[rank], sendtype, FMPI_INIT);
          MPI_Isend(buf_tmp[total_reqs], sendcounts[rank], sendtype, rank, MPI_TAG_ALL, comm, &reqs[total_reqs]);
          total_reqs++;
        }
      } else {
        if (*sendcounts) {
          buf_tmp[total_reqs] = kt_malloc(fmpi_sizeof(sendtype, *sendcounts));
          fmpi_reduce(sendbuf+(rank*(*sendcounts)*type_size), buf_tmp[total_reqs], *sendcounts, sendtype, FMPI_INIT);
          MPI_Isend(buf_tmp[total_reqs], *sendcounts, sendtype, rank, MPI_TAG_ALL, comm, &reqs[total_reqs]);
          total_reqs++;
        }
      }
    }
    ar_assert(total_reqs <= rank_end-rank_start+1);
    ar_assert(total_reqs <= FMPI_ALL2ALL_GS);

    // gather data if we are in the group
    if (context->rank >= rank_start && context->rank <= rank_end) {
      type_size = fmpi_sizeof(recvtype, 1);
      for (rank = 0; rank <= max_rank; rank++) {
        // perform mpi_recv and copy from alligned space to data 
        if (rdispls) {
          // resize if not sufficient space
          if (recvcounts[rank] > rbuf_size) {
            kt_free(buf_tmp[FMPI_ALL2ALL_GS]);
            buf_tmp[FMPI_ALL2ALL_GS] = kt_malloc(fmpi_sizeof(recvtype, recvcounts[rank]+FMPI_DMA_WORD));
            rbuf_size = recvcounts[rank];
          }
          if (recvcounts[rank]) {
            MPI_Recv(buf_tmp[FMPI_ALL2ALL_GS], recvcounts[rank], recvtype, rank, MPI_TAG_ALL, comm, NULL);
            fmpi_reduce(buf_tmp[FMPI_ALL2ALL_GS], recvbuf+(rdispls[rank]*type_size), recvcounts[rank], recvtype, FMPI_INIT);
          }
        } else {
          if (*recvcounts) {
            MPI_Recv(buf_tmp[FMPI_ALL2ALL_GS], *recvcounts, recvtype, rank, MPI_TAG_ALL, comm, NULL);
            fmpi_reduce(buf_tmp[FMPI_ALL2ALL_GS], recvbuf+(rank*(*recvcounts)*type_size), *recvcounts, recvtype, FMPI_INIT);
          }
        }
      }
    }

    // wait for requests to finish in order to free resources
    MPI_Waitall(total_reqs, reqs, &status);

    // free tmp send buffers 
    for (j = 0; j < total_reqs; j++) kt_free(buf_tmp[j]);
  }

  // free tmp recv buffer 
  kt_free(buf_tmp[FMPI_ALL2ALL_GS]);

  // enable interrupts
  fmpi_intr_enable(context);

  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int MPI_Alltoall(volatile void *sendbuf, int sendcount, MPI_Datatype sendtype,
                  volatile void *recvbuf, int recvcount,
                  MPI_Datatype recvtype, MPI_Comm comm) {

  MPI_Alltoallv(sendbuf, &sendcount, NULL, sendtype, recvbuf, &recvcount, NULL, recvtype, comm);

  return MPI_SUCCESS;
}

