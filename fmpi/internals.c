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
// Abstract      : Minimal MPI library: internal functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: internals.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2012/10/24 13:06:55 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <memory_management.h>
#include <noc.h>
#include <fmpi.h>

     
// ===========================================================================
// ===========================================================================
#define fmpi_align_int(num)  ((int) ((num + sizeof(int) - 1) / sizeof(int)))

// ===========================================================================
// ===========================================================================
static inline void fmpi_set_status(FMPI_UserReq *request, int status) {
  request->status = status;
  if (request->ext_request) {
    request->ext_request->status = status;
  }
} 


// ===========================================================================
// ===========================================================================
void fmpi_get_bid_cid(FMPI_Context *context, int rank, 
                                    int *bid, int *cid) {
  ar_assert(rank >= 0);
  ar_assert(rank < context->num_ranks);

  *bid = context->rank_bid_cid[rank] >> 8;
  *cid = context->rank_bid_cid[rank] & 0xFF;
}


// ===========================================================================
// ===========================================================================
int fmpi_req_check(FMPI_Context *context, volatile void *buf, 
                   int count, MPI_Datatype datatype, int dest, 
                   int tag, MPI_Comm comm) {

  // Sanity checks
  // Check for 64-B alignment
  if ((!buf) ||
      ((int) buf & (FMPI_DMA_WORD-1))) {
    return MPI_ERR_BUFFER;
  }

  // MPI_ERR_COMM
  // MPI_ERR_COUNT
  if (count == 0) {
    return MPI_ERR_COUNT;
  }

  // MPI_ERR_TYPE
  int data_size = fmpi_sizeof(datatype, count);
  if (!data_size) {
    return MPI_ERR_TYPE;
  }
  // MPI_ERR_TAG
  if (tag < 0 && tag != MPI_ANY_TAG && tag != MPI_TAG_BCAST && tag != MPI_TAG_REDUCE && tag != MPI_TAG_ALL) {
    return MPI_ERR_TAG;
  }

  // MPI_ERR_RANK
  if (dest < 0 || dest >= context->num_ranks) {
    return MPI_ERR_RANK;
  }

  return MPI_SUCCESS;
}


// ###########################################################################
// ###                                                                     ###
// ###                            FMPI general routines                    ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ===========================================================================
FMPI_List *fmpi_list_init(FMPI_List *first, int nums) {

  int           i;
  FMPI_List     *free = first;


  // create free list
  for (i = 0; i < nums; i++) {
    free[i].id = i;

    if (i == nums-1) {
      free[i].next = NULL;
    } else {
      free[i].next = &free[i+1];
    }
  }
  
  // return first pointer
  return first;
}


// ===========================================================================
// ===========================================================================
FMPI_List *fmpi_list_alloc(FMPI_List **free, FMPI_List **used_first, FMPI_List **used_last) {

  FMPI_List *ret_object = *free;


  ar_assert(*free);  // list is not empty

  // remove from free list
  if (*free) {
    *free = (*free)->next;
  }

  // append to used list
  if (ret_object) {
    ret_object->prev = *used_last;
    ret_object->next = NULL;
    if (*used_last) {
      (*used_last)->next = ret_object;
    }
    *used_last = ret_object;
    ar_assert(*used_last);
  }

  // fix pointer to first element
  if (!*used_first)
    *used_first = ret_object;

  ar_assert(*used_first);

  return ret_object;
}


// ===========================================================================
// ===========================================================================
void fmpi_list_free(FMPI_List **free, FMPI_List **used_first, FMPI_List **used_last, FMPI_List *object) {

  ar_assert(object);
  ar_assert(*used_last);
  ar_assert(*used_first);

  // fix previous node
  if (object->prev) {
    object->prev->next = object->next;
  } else {
    ar_assert(*used_first == object);
    *used_first = object->next;
    if (*used_first) {
      (*used_first)->prev = NULL;
    }
  }
    
  // fix next node
  if (object->next) {
    object->next->prev = object->prev;
  } else {
    ar_assert(*used_last == object);
    *used_last = object->prev;
    if (*used_last) {
      (*used_last)->next = NULL;
    } else {
      ar_assert(!*used_first);
    }
  }

  // add object to free list
  object->next = *free;
  *free = object;
}


// ===========================================================================
// ===========================================================================
#if 0
void fmpi_print_list(FMPI_List *object) {
  FMPI_List *i = object;

  kt_printf("List: ");
  while (i) {
    kt_printf("%d ", (int)i->id);
    i = i->next;
  }
  kt_printf("\r\n");
}
#endif


// ===========================================================================
// ===========================================================================
void fmpi_mbx_barrier(FMPI_Context *context) {
  int i;
  unsigned int word;
  int bid, cid;
  int my_cid;

  my_cid = ar_get_core_id();

  // master 
  if (context->rank == 0) {
    // Gather replies
    for (i = 0; i < context->num_ranks - 1; i++) {
      word = ar_mbox_get(my_cid);
      ar_assert(((word >> 12) & 0xFF) == 0x42);
    }

    // Send acks
    word = (0x45 << 12) | context->rank;
    for (i = 0; i < context->num_ranks; i++) {
      fmpi_get_bid_cid(context, i, &bid, &cid);
      if (i != context->rank) {
        ar_mbox_send(my_cid, bid, cid, word);
      }
    }
  }

  // slaves
  else {
    // Send "I'm alive" message
    word = (0x42 << 12) | context->rank;
    fmpi_get_bid_cid(context, 0, &bid, &cid);
    ar_mbox_send(my_cid, bid, cid, word);

    // Receive ack
    word = ar_mbox_get(my_cid);
    ar_assert(((word >> 12) & 0xFF) == 0x45);
  }
}

// ###########################################################################
// ###                                                                     ###
// ###                            FMPI Interrupt routines                  ###
// ###                                                                     ###
// ###########################################################################


// ===========================================================================
// ===========================================================================
void fmpi_cnt_intr_new(FMPI_Context *context, int cnt) {
#ifdef FMPI_INTR_ENABLE

  int cnt_pos = context->intr_cnt_pos;

  // fix counter interrupt
  unsigned int cnt_intr = ar_intr_cnt_get();
  switch (cnt_pos) {
    case 1: ar_intr_cnt_set((cnt_intr & 0xffffff00) | cnt); break;
    case 2: ar_intr_cnt_set((cnt_intr & 0xffff00ff) | (cnt << 8)); break;
    case 4: ar_intr_cnt_set((cnt_intr & 0xff00ffff) | (cnt << 16)); break;
    case 8: ar_intr_cnt_set((cnt_intr & 0x00ffffff) | (cnt << 24)); break;
    default : ar_abort(); break;
  }

  // fix counter mask
  context->intr_mask &= ~(cnt_pos<<2);

  // if interrupt enable fix interrupt mask
  if (context->intr_enabled) ar_intr_cpu_set(context->intr_mask);
  //kt_printf("cnt=%08X cnt_inr=%08X cpu_intr=%08X\n\r", cnt, ar_intr_cnt_get(), ar_intr_cpu_get());

  // intr_cnt_pos selets counters in RR fashion
  if (context->intr_cnt_pos == 8) context->intr_cnt_pos = 1;
  else context->intr_cnt_pos <<= 1;

#endif
}


// ===========================================================================
// ===========================================================================
#ifdef FMPI_INTR_ENABLE
void fmpi_intr_serve() {
  Context       *context;
  unsigned int  mi_interrupts;

  // Get context
  context = mm_get_context(ar_get_core_id());

  // read interrupts
  mi_interrupts = (ar_intr_status_get() >> 8) & 0xFF;

  // serve pending requests
  fmpi_process_pending_events(context->fmpi);

  // clear old interrupts. Don't clear new interrupts and mailbox!
  context->fmpi->intr_mask = (ar_intr_cpu_get() | mi_interrupts) & (~0x2);
  ar_intr_cpu_set(context->fmpi->intr_mask | (mi_interrupts << 16));
}
#endif


// ###########################################################################
// ###                                                                     ###
// ###                            FMPI send/receive routines               ###
// ###                                                                     ###
// ###########################################################################


// ===========================================================================
// ===========================================================================
void fmpi_serve_recv(FMPI_Context *context, FMPI_UserReq *request) {

  int bid;
  int cid;
  int peer_bid, peer_cid;
  unsigned int src_addr, my_addr;
  int dma_size;
  FMPI_List *descr_recv, *descr_next;
  FMPI_DescrRecv *descr_pending;
  int descr_tag;
  int descr_match;


  bid = ar_get_board_id();
  cid = ar_get_core_id();

  // Search for a matching, ready, recv descriptor. If none, block here 
  // and poll for new events until match is found.
  if (request->status == FMPI_REQ_PENDING) {

    ar_assert(request->descr == NULL);
    ar_assert(request->cnt == NULL);

    //reset match flag
    descr_match = 0;

    // go through all recv descriptors
    descr_recv = context->descr_recv_used_first;

    while (descr_recv) {
      descr_pending = &context->descr_recv[descr_recv->id];

      // maybe descriptor has just arrived?
      if (!descr_pending->arrived) {
        ar_assert(descr_pending->cnt != NULL);
        if (ar_cnt_get(cid, descr_pending->cnt->id) == 0) {
          // mark it's here
          descr_pending->arrived = 1;

          // free counter, so it can be reused for other purposes
          fmpi_list_free(&context->cnt_free, &context->cnt_used_first, 
                         &context->cnt_used_last, descr_pending->cnt);
          descr_pending->cnt = NULL;
        }
      }

      descr_next = descr_recv->next;

      // if descriptor has arrived
      if (!descr_match && descr_pending->arrived) {

        // read tag
        descr_tag = context->descr[descr_pending->descr->id].src_tag;
//        ar_assert(descr_tag < 1000); // delete this

        // check for match
        ar_assert(descr_pending->peer_rank < 1024);

        // if match
        if ((request->peer_rank == MPI_ANY_SOURCE || 
             descr_pending->peer_rank == request->peer_rank) &&
            (request->tag == MPI_ANY_TAG || descr_tag == request->tag)) {

          // move data in user request
          request->descr = descr_pending->descr;
          request->peer_cnt = descr_pending->peer_cnt;

          // fix status
          fmpi_set_status(request, FMPI_REQ_INIT);

          // free recv descriptor
          fmpi_list_free(&context->descr_recv_free, &context->descr_recv_used_first, 
                         &context->descr_recv_used_last, descr_recv);

          // set match flag
          descr_match = 1;
          
          // break; 
        }
      }

      descr_recv = descr_next;
    }
  }


  if (request->status == FMPI_REQ_INIT) {
    ar_assert(request->descr != NULL);
    ar_assert(request->cnt == NULL);
    ar_assert(request->peer_cnt < NOC_MAX_COUNTERS);

    // My DMA engine can support one more outgoing
    if ((ar_ni_status_get(cid) & 0xFF) < FMPI_DMA_THRESHOLD) return;

    // There is a free counter
    if (!context->cnt_free) return;

    // Allocate counter
    request->cnt = fmpi_list_alloc(&context->cnt_free, &context->cnt_used_first, &context->cnt_used_last);

    // Do a DMA from peer buffer to user buffer
    fmpi_get_bid_cid(context, request->peer_rank, &peer_bid, &peer_cid);
        int descr_id = request->descr->id;
    src_addr = (unsigned int) context->descr[descr_id].src_buf;
    my_addr = (unsigned int) request->usr_buf;
    dma_size = (unsigned int) context->descr[descr_id].size;

#ifdef FMPI_EAGER_ENABLE
    if (dma_size <= FMPI_EAGER_SIZE) {
      int i;
      for (i = 0; i < fmpi_align_int(dma_size); i++) {
        ((unsigned int *)my_addr)[i] = context->descr[descr_id].eager_payload[i];
      }
      fmpi_set_status(request, FMPI_REQ_DONE);
    } else {
#endif

    // 64-byte alignment
    dma_size = ((dma_size + (FMPI_DMA_WORD-1)) >> 6) * FMPI_DMA_WORD;

    // send DMA with notification
    ar_cnt_set_notify_cnt(cid, request->cnt->id, -dma_size,
                          peer_bid, peer_cid, request->peer_cnt);

    ar_dma_with_ack(cid, 
                    peer_bid, peer_cid, (unsigned int) src_addr,
                    bid, cid, (unsigned int) my_addr, 
                    bid, cid, request->cnt->id, 
                    dma_size, 0, 0, 0);

    // fix status
    fmpi_set_status(request, FMPI_REQ_PROGRESS);

#ifdef FMPI_EAGER_ENABLE
    }
#endif
  }
  
  if (request->status == FMPI_REQ_PROGRESS) {
    ar_assert(request->descr != NULL);
    ar_assert(request->cnt != NULL);
    
    // if DMA is done
    volatile int cnt_recv = ar_cnt_get(cid, request->cnt->id);
    if (cnt_recv == 0) {
      // done
      fmpi_set_status(request, FMPI_REQ_DONE);
    }
  }
}



// ===========================================================================
// ===========================================================================
void fmpi_serve_send(FMPI_Context *context, FMPI_UserReq *request) {

  int cid;
  volatile FMPI_Descriptor *descr_send;
  unsigned int word_send;
  int bid_dest, cid_dest;
  volatile int cnt_send;


  cid = ar_get_core_id();


  if (request->status == FMPI_REQ_PENDING) {
    // Check if request can proceed:
    // - Incoming DMAs > safe value (context->incoming_dmas)
    if (context->incoming_dmas < FMPI_DMA_THRESHOLD) {
      return;
    }

    // - My DMA engine can support one more outgoing (MNI_STATUS CPU_OPS)
    if ((ar_ni_status_get(cid) & 0xFF) < FMPI_DMA_THRESHOLD) {
      return;
    } 

    // - We have mailbox credits for the peer
    if (!context->mbox_credits[request->peer_rank]) {
      return;
    } 

    // - There is a free descriptor
    if (!context->descr_free) {
      // This should not happen. If we run out of descriptors, it is possible
      // that a deadlock will occur, because a rank cannot receive outstanding
      // sends. Decrease application non-blocking sends/recvs or increase
      // the FMPI number of descriptors if this fails.
      ar_panic("FMPI: Out of descriptors in serve_send()");
      //return;
    }

    // - There is a free send counter
    if (!context->cnt_free) {
      //ar_panic("FMPI: Out of descriptors in serve_send()");
      return;
    } 


    // process request
    // Decrement incoming DMAs
    context->incoming_dmas--;
    ar_assert(context->incoming_dmas > -1);

    // Allocate & prepare descriptor
    request->descr = fmpi_list_alloc(&context->descr_free, 
                                     &context->descr_used_first,
                                     &context->descr_used_last);
    ar_assert(request->descr != NULL);
    descr_send = &context->descr[request->descr->id];
    descr_send->src_tag = request->tag;
    descr_send->src_buf = (void *)request->usr_buf;
    descr_send->size = request->size; 

#ifdef FMPI_EAGER_ENABLE
    if (request->size <= FMPI_EAGER_SIZE) {
      int i;

      // copy data in descriptor. No need for second dma.
      for (i = 0; i < fmpi_align_int(request->size); i++) {
        descr_send->eager_payload[i] = ((unsigned int *)request->usr_buf)[i];
      }
    }
#endif

    // Allocate & initialize counter 
    request->cnt = fmpi_list_alloc(&context->cnt_free, &context->cnt_used_first, &context->cnt_used_last);

    ar_cnt_set(cid, request->cnt->id, -2);


    // Decrease mailbox credit & send mailbox message to peer
    context->mbox_credits[request->peer_rank]--;

    ar_assert(request->cnt->id < NOC_MAX_COUNTERS);
    ar_assert(request->descr->id < FMPI_NUM_DESCRIPTORS);
    ar_assert(context->rank < 1024);

    // pack send word for mailbox
    // FIXME: use FMPI_MboxMessage instead
    word_send = request->cnt->id | (request->descr->id << 8) |
                (context->rank << 20) | (FMPI_MBX_SEND << 30);

    // find destination
    fmpi_get_bid_cid(context, request->peer_rank, &bid_dest, &cid_dest);

    // send mailbox
    ar_mbox_send(cid, bid_dest, cid_dest, word_send);
    
    // fix status
    fmpi_set_status(request, FMPI_REQ_INIT);
  }


  if (request->status == FMPI_REQ_INIT) {
    // Check counter
    cnt_send = (int) ar_cnt_get(cid, request->cnt->id);
    if (cnt_send > -2) {
      // Increase mailbox credit
      context->mbox_credits[request->peer_rank]++;
      
      // fix status
      fmpi_set_status(request, FMPI_REQ_PROGRESS);
    }
  }

  if (request->status == FMPI_REQ_PROGRESS) {
    // Check counter
    cnt_send = ar_cnt_get(cid, request->cnt->id);
#ifdef FMPI_EAGER_ENABLE
    // No dma for eager data so move on.
    if (cnt_send == 0 || request->size <= FMPI_EAGER_SIZE) {
#else
    if (cnt_send == 0) {
#endif
      // Increament incoming DMAs
      context->incoming_dmas++;

      // done
      fmpi_set_status(request, FMPI_REQ_DONE);
    }
  }
}


// ===========================================================================
// ===========================================================================
void fmpi_descr_fetch(FMPI_Context *context, int descr_recv_id) {

  FMPI_DescrRecv        *descr_pending;
  int                   peer_bid;
  int                   peer_cid;
  int                   bid;
  int                   cid;
  unsigned int          word_recv;
  int                   peer_descr_id;
  FMPI_Context          *peer_context;
  FMPI_Descriptor       *peer_descr;
  FMPI_Descriptor       *descr_recv;


  bid = ar_get_board_id();
  cid = ar_get_core_id();
  descr_pending = &context->descr_recv[descr_recv_id];

  // read incoming word
  word_recv = context->mbx_recv;
  context->mbx_pending = 0;

  // unpack word
  // FIXME: use FMPI_MboxMessage instead
  descr_pending->peer_cnt = word_recv & 0xFF;
  peer_descr_id = (word_recv >> 8) & 0xFFF;
  descr_pending->peer_rank = (word_recv >> 20) & 0x3FF;
  ar_assert((descr_pending->peer_rank >= 0) &&
            (descr_pending->peer_rank < context->num_ranks));

  // find bids, cpuids
  fmpi_get_bid_cid(context, descr_pending->peer_rank, &peer_bid, &peer_cid);
  
  // find address of descriptors

  peer_context = context->rank_context[descr_pending->peer_rank];
  peer_descr = &peer_context->descr[peer_descr_id];
  descr_recv = &context->descr[descr_pending->descr->id];

  // Mark it's not here yet
  descr_pending->arrived = 0;

  // start dma in order to fetch descriptor
  ar_cnt_set_notify_cnt(cid, descr_pending->cnt->id, -FMPI_DMA_WORD,
                        peer_bid, peer_cid, descr_pending->peer_cnt);

  // set counter interrupt
  fmpi_cnt_intr_new(context, descr_pending->cnt->id);

  ar_dma_with_ack(cid, 
                  peer_bid, peer_cid, (unsigned int) peer_descr,
                  bid, cid, (unsigned int) descr_recv, 
                  bid, cid, descr_pending->cnt->id, 
                  FMPI_DMA_WORD, 0, 0, 0);
}


// ###########################################################################
// ###                                                                     ###
// ###                            FMPI process commands                    ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ===========================================================================
int fmpi_commit_cmd(volatile void *buf, int count, MPI_Datatype datatype, 
                    int dest, int tag, MPI_Comm comm, FMPI_Opcode opcode, 
                    MPI_Status *status, MPI_Request *request) {

  Context       *context;
  int           req_check_results;
  FMPI_List     *req_new;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // disable interrupts
  int intr_enable = context->fmpi->intr_enabled;
  fmpi_intr_disable(context->fmpi);

  // make sure command is ok
  req_check_results = fmpi_req_check(context->fmpi, buf, count, datatype, 
                                     dest, tag, comm);
  if (req_check_results != MPI_SUCCESS) {
    ar_abort();
    return req_check_results;
  }

  // Transfer to user request. Assign new handle. Mark as pending.
  // If no more space, block here and poll for more progress.

  // block until a free user request is found
  while (!(req_new = fmpi_list_alloc(&context->fmpi->user_req_free, 
                                     &context->fmpi->user_req_used_first,
                                     &context->fmpi->user_req_used_last))) {
    ar_abort();
    fmpi_process_pending_events(context->fmpi);
  }

  // initialize user request
  context->fmpi->user_req[req_new->id].ext_status = status;
  if (status) {
    (context->fmpi->user_req[req_new->id].ext_status)->MPI_SOURCE = dest;
    (context->fmpi->user_req[req_new->id].ext_status)->MPI_TAG = tag;
    (context->fmpi->user_req[req_new->id].ext_status)->MPI_ERROR = 0;
  }
  context->fmpi->user_req[req_new->id].ext_request = request;
  fmpi_set_status(&context->fmpi->user_req[req_new->id], FMPI_REQ_PENDING);
  context->fmpi->user_req[req_new->id].opcode = opcode;
  context->fmpi->user_req[req_new->id].usr_buf = (void *)buf;
  context->fmpi->user_req[req_new->id].size = fmpi_sizeof(datatype, count);
  context->fmpi->user_req[req_new->id].peer_rank = dest;
  context->fmpi->user_req[req_new->id].tag = tag;
  context->fmpi->user_req[req_new->id].descr = NULL;
  context->fmpi->user_req[req_new->id].cnt = NULL;


  // execute all pending commands
  fmpi_process_pending_events(context->fmpi);

  // if blocking command wait until done
  if (opcode == FMPI_OPCODE_SEND || opcode == FMPI_OPCODE_RECV) {
    while (context->fmpi->user_req[req_new->id].status != FMPI_REQ_DONE) {
      fmpi_process_pending_events(context->fmpi);
    }
  }

  // enable interrupts
  if (intr_enable) fmpi_intr_enable(context->fmpi);

  return MPI_SUCCESS;
}


// ===========================================================================
// ===========================================================================
int fmpi_mbx_recv_serve(FMPI_Context *context) {

  FMPI_List *descr_recv;
  int cid;

  cid = ar_get_core_id();

  // Is there a pending mailbox word to process?
  if (!context->mbx_pending) {
    context->mbx_recv = ar_mbox_get(cid);
  }
  context->mbx_pending = 1;

  // get opcode of incoming word
  int mbx_opcode = context->mbx_recv >> 30;

  // serve incoming mail
  if (mbx_opcode == FMPI_MBX_SEND) {

    // Try to allocate receive descriptor

    // Can my DMA engine support one more outgoing?
    if ((ar_ni_status_get(ar_get_core_id()) & 0xFF) < FMPI_DMA_THRESHOLD) {
      ar_panic("FMPI: Out of NI in descr_recv_alloc()");
      goto fail;
    }

    // - There is a free descriptor
    if (!context->descr_free) {
      // This should not happen. If we run out of descriptors, it is possible
      // that a deadlock will occur, because a rank cannot receive outstanding
      // sends. Decrease application non-blocking sends/recvs or increase
      // the FMPI number of descriptors if this fails.
      ar_panic("FMPI: Out of descriptors in descr_recv_alloc()");
    }

    // - There is a free recv counter
    if (!context->descr_recv_free) {
      // See above.
      ar_panic("FMPI: Out of receive descriptors in descr_recv_alloc()");
    }

    // - There is a free recv counter
    if (!context->cnt_free) {
      ar_panic("FMPI: Out of counters in descr_recv_alloc()");
      goto fail;
    }

    // Allocate descr recv
    descr_recv = fmpi_list_alloc(&context->descr_recv_free,    
                                 &context->descr_recv_used_first,
                                 &context->descr_recv_used_last);

    // Allocate descriptor
    context->descr_recv[descr_recv->id].descr = fmpi_list_alloc(
                                  &context->descr_free, &context->descr_used_first, &context->descr_used_last);

    // Allocate counter
    context->descr_recv[descr_recv->id].cnt = fmpi_list_alloc(
                                  &context->cnt_free, &context->cnt_used_first, &context->cnt_used_last);

    // Mark it's not valid
    context->descr_recv[descr_recv->id].arrived = 0;

    if (descr_recv) {

      // fetch descriptor
      fmpi_descr_fetch(context, descr_recv->id);

      return 1;
    }
  }

fail:

  return 0;
}


// ===========================================================================
// ===========================================================================
void fmpi_process_pending_events(FMPI_Context *context) {
  
  int cid;
  FMPI_List *user_req;
  FMPI_List *user_req_next;
  FMPI_UserReq *request;

  cid = ar_get_core_id();

  // If:
  // - mailbox has messages 
  //    => Dequeue mailbox word
  //       Do a DMA from peer descriptor to my desciptor
  // (loop)

  int mbx_done = 0;
  while (!mbx_done) {
    mbx_done = 1;
    // check that there is a word in the mailbox
    if (context->mbx_pending || (ar_mbox_status_get(cid) & 0xFFFF)) {
      // serve incoming mail
      if (fmpi_mbx_recv_serve(context)) {
              mbx_done = 0;
      }
    }
  }

  // try to serve all pending requests
  user_req = context->user_req_used_first;

  while (user_req) {

    request = &context->user_req[user_req->id];
    ar_assert(request->status != FMPI_REQ_DONE);

    switch (context->user_req[user_req->id].opcode) {
      case FMPI_OPCODE_SEND: 
      case FMPI_OPCODE_ISEND: fmpi_serve_send(context, request); break;
      case FMPI_OPCODE_RECV: 
      case FMPI_OPCODE_IRECV: fmpi_serve_recv(context, request); break;
      default : ar_abort(); break;
    }

    user_req_next = user_req->next;

    // free request if done
    if (request->status == FMPI_REQ_DONE) {

      // free counter
      fmpi_list_free(&context->cnt_free, &context->cnt_used_first, 
                     &context->cnt_used_last, request->cnt);
              
      // free descriptor
      fmpi_list_free(&context->descr_free, &context->descr_used_first, 
                     &context->descr_used_last, request->descr);

      // free request
      fmpi_list_free(&context->user_req_free, &context->user_req_used_first, 
                     &context->user_req_used_last, user_req);
    } 

    user_req = user_req_next;
  }
  // fix receive requests. check for pending descriptors.

}


// ###########################################################################
// ###                                                                     ###
// ###                        BARRIER ROUTINES                             ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ===========================================================================
void fmpi_barrier_send(FMPI_Context *context, int dest, int cnt) {

  int cid;
  int bid_dest, cid_dest;

  cid = ar_get_core_id();

  // wait until my DMA engine can support one more outgoing
  while ((ar_ni_status_get(cid) & 0xFF) < FMPI_DMA_THRESHOLD) {
    ;
  }

  // increment counter of dest
  fmpi_get_bid_cid(context, dest, &bid_dest, &cid_dest);
  ar_cnt_incr(cid, bid_dest, cid_dest, cnt, 1);
}



// ===========================================================================
// ===========================================================================
void fmpi_find_children(FMPI_Context *context) {
  int level_step = 1;
  int children_below = 1;
  int level_rank = context->rank;
  int level_max_rank = context->num_ranks - 1;
  int i, child_rank, child_i;
  
  context->num_children = 0;
  child_i = 0;

  while ((level_rank & 0x7) == 0 && (children_below > 0)) {
    if (level_max_rank > level_rank + 7) {
      children_below = 7;
    } else  {
      children_below = level_max_rank - level_rank;
    }
    ar_assert(children_below >= 0);
    
    for (i = 0; i < children_below; i++) {
      child_rank = context->rank + (i+1) * level_step;
      ar_assert(child_i < 24);
      context->children[child_i++] = child_rank;
    }
    level_step *= 8;
    
    
    level_max_rank >>= 3;
    level_rank >>= 3;
  }
  
  context->num_children = child_i;
}


// ===========================================================================
// ===========================================================================
void fmpi_find_father(FMPI_Context *context) {
  int level = 0;

  // root has no father
  if (!context->rank) {
    context->father = -1;
    return;
  }

  while (1) {
    if ((context->rank & (0x7 << (3*level))) != 0) {
      context->father = context->rank & (0xFF8 << (3*level));
      return;
    }
    level++;
  }
}

// ===========================================================================
// ===========================================================================
void fmpi_barrier_init(FMPI_Context *context) {
  int cid;

  cid = ar_get_core_id();

  ar_cnt_set(cid, context->cnt_breq, -context->num_children);
  ar_cnt_set(cid, context->cnt_back, -1);
}


// ==========================================================================
// returns the father or child where the bcast will come from
// ==========================================================================
int fmpi_bcast_recv(FMPI_Context  *context, int src) {

  int dest = context->rank;

  if (dest == src) {
    return -1;
  }

  // find common part
  int mask = 0x7;
  while ((src & ~mask) != (dest & ~mask)) {
    mask = (mask << 3) | 0x7;
  }

  // bcast comes from father
  if ((dest & mask) != 0) {
    return context->father;
  }

  // bcast comes from a child
  return (src & ((~mask) >>3));
}


// ==========================================================================
// reduce operation between two arrays
// ==========================================================================
void fmpi_char_op(char *data, char *result, MPI_Op op) {
    switch (op) {
        case FMPI_INIT : *result = *data; break;
        case MPI_MAX   : if (*data > *result) *result = *data; break;
        case MPI_MIN   : if (*data < *result) *result = *data; break;
        case MPI_SUM   : *result += *data; break;
        // FIXME: add more operations here
        default:                ar_abort(); break;
    }
}

// ==========================================================================
// ==========================================================================
void fmpi_short_op(short *data, short *result, MPI_Op op) {
    switch (op) {
        case FMPI_INIT : *result = *data; break;
        case MPI_MAX   : if (*data > *result) *result = *data; break;
        case MPI_MIN   : if (*data < *result) *result = *data; break;
        case MPI_SUM   : *result += *data; break;
        // FIXME: add more operations here
        default:                ar_abort(); break;
    }
  }

// ==========================================================================
// ==========================================================================
void fmpi_int_op(int *data, int *result, MPI_Op op) {
    switch (op) {
        case FMPI_INIT : *result = *data; break;
        case MPI_MAX   : if (*data > *result) *result = *data; break;
        case MPI_MIN   : if (*data < *result) *result = *data; break;
        case MPI_SUM   : *result += *data; break;
        // FIXME: add more operations here
        default:                ar_abort(); break;
    }
}

// ==========================================================================
// ==========================================================================
void fmpi_float_op(float *data, float *result, MPI_Op op) {
    switch (op) {
        case FMPI_INIT : *result = *data; break;
        case MPI_MAX   : if (*data > *result) *result = *data; break;
        case MPI_MIN   : if (*data < *result) *result = *data; break;
        case MPI_SUM   : *result += *data; break;
        // FIXME: add more operations here
        default:                ar_abort(); break;
    }
}

// ==========================================================================
// ==========================================================================
void fmpi_reduce(volatile void *data, volatile void *result, int count, 
                 MPI_Datatype datatype, MPI_Op op) {
    int i;

    for (i = 0; i< count; i++) {
        switch (datatype) {
            case MPI_BYTE:
            case MPI_CHAR:          fmpi_char_op((char *)data, (char *)result, op); break;
            case MPI_SHORT:         fmpi_short_op((short *)data, (short *)result, op); break;
            case MPI_INT:
            case MPI_LONG:          fmpi_int_op((int *)data, (int *)result, op); break;
            case MPI_FLOAT:         fmpi_float_op((float *)data, (float *)result, op); break;
            case MPI_UNSIGNED_CHAR:
            case MPI_DOUBLE:
            case MPI_LONG_DOUBLE:
            default:                ar_abort(); break;
        }

        data += fmpi_sizeof(datatype, 1);
        result += fmpi_sizeof(datatype, 1);
    }

}


// ===========================================================================
// ===========================================================================
unsigned int FMPI_Wtime(void) {
#ifdef ARCH_MB
  unsigned int now;
  unsigned int time_diff;

  // Get context
  FMPI_Context *context = mm_get_context(ar_get_core_id())->fmpi;

  now = ar_glob_timer_read();
  if (now > context->time_init) 
    time_diff = now - context->time_init;
  else 
    time_diff = 0xFFFFFFFF - (now - context->time_init);
  if (context->time_last > time_diff) context->time_ovfl++;
  context->time_last = time_diff;
  return context->time_last;
#endif
#ifdef ARCH_ARM
  ar_abort(); // fixme
#endif
}

