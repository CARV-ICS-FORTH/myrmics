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
// Author        : Spyros LYBERIS
// Abstract      : Network-on-chip message-based communication functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: message.c,v $
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2012/12/18 16:40:29 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <debug.h>


// ===========================================================================
// noc_msg_send_get_buf()       Returns a buffer of size NOC_MESSAGE_SIZE 
//                              bytes, which will be used for a subsequent
//                              call to noc_msg_send() to send it to
//                              dst_core_id.
// ===========================================================================
// * INPUTS
//   int dst_core_id            Unique core ID of the destination
//
// * RETURN VALUE
//   void *                     The buffer to be filled with data by the
//                              caller
// ===========================================================================
void *noc_msg_send_get_buf(int dst_core_id) {

  Context               *context;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Remember the core ID
  ar_assert(context->noc_msg_core_id == -1);
  ar_assert((dst_core_id >= 0) && (dst_core_id < context->pr_num_cores));
  context->noc_msg_core_id = dst_core_id;

  // Mailbox-only mode
  if (context->noc_mode == NOC_MODE_MAILBOX_ONLY) {

    // The single buffer can be used over and over again, because in this
    // mode noc_msg_send() will wait until the DMA is finished.
    return context->noc_send_buf;
  }

  // Credit-based modes
  else if ((context->noc_mode == NOC_MODE_CREDIT_ONLY) ||
           (context->noc_mode == NOC_MODE_CREDIT_MAILBOX)) {

    // Extra sanity checks
    ar_assert(context->noc_credits);
    ar_assert(context->noc_credits[dst_core_id]);

    // This buffer can be used, because the credit-based system guarantees
    // that we won't be using more buffers than there can be in flight.
    return (void *) (((unsigned int) 
                        context->noc_credits[dst_core_id]->loc_send_bufs) +
                     context->noc_credits[dst_core_id]->cur_send_buf_id * 
                     NOC_MESSAGE_SIZE);
  }

  // Unknown mode
  else {
    ar_abort();
  }

  return NULL;
}


// ===========================================================================
// noc_msg_send()               Send a software message to another core. This
//                              function should be called after a call to
//                              noc_msg_send_get_buf() and after the caller
//                              has filled the returned buffer. The destination
//                              core ID is also specified during the call to
//                              noc_msg_send_get_buf().
//
//                              In mailbox-only mode, the buffer is directly
//                              DMAed to the destination mailbox. The function
//                              waits until the DMA is finished successfully.
//                              Because there are no credits, this may not
//                              succeed because the destination mailbox may
//                              be full. If the transfer fails, the function
//                              retries as many times as needed.
//
//                              In credit-based modes, the function waits until
//                              there are enough credits to guarantee the
//                              transfer. Then it DMAs the buffer to a direct
//                              prearranged buffer on the destination. When
//                              using mailbox polling, without waiting it also
//                              sends a single mailbox word to indicate which
//                              buffer was filled and which remote counter is
//                              used to track the DMA completion.
// ===========================================================================
// * RETURN VALUE
//   int                        0 for success
//                              ERR_OUT_OF_COUNTERS: no available counters
// ===========================================================================
int noc_msg_send() {

  Context               *context;
  int                   my_bid;
  int                   my_cid;
  int                   dst_core_id;
  int                   dst_bid;
  int                   dst_cid;
  unsigned int          src_buf;
  unsigned int          dst_buf;
  ListNode              *cnt_node;
  int                   cnt;
  int                   ret;
  unsigned int          credits;
  unsigned int          word;


  // Get context and architectural board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Verify noc_msg_send_get_buf() has been called before this function
  ar_assert(context->noc_msg_core_id > -1);
  dst_core_id = context->noc_msg_core_id;
  context->noc_msg_core_id = -1; // reset it for the next call

  // Sanity checks
  ar_assert((dst_core_id >= 0) && (dst_core_id < context->pr_num_cores));

  // Discover destination arch board/core ID
  pr_core_arch_bid_cid(dst_core_id, &dst_bid, &dst_cid);

  // Take a note for the statistics gathering
  dbg_scount(context, DBG_STATS_IDX_NUM_MESSAGES, 1);


  // =========================================================================
  // Mailbox-only communication mode
  // =========================================================================
  if (context->noc_mode == NOC_MODE_MAILBOX_ONLY) {

    // Find a free counter, but don't dequeue it; we'll use it momentarily.
    cnt_node = kt_list_head(context->noc_cnt_free);
    if (!cnt_node) {
      return ERR_OUT_OF_COUNTERS;
    }
    cnt = (int) cnt_node->data;


    // Forever
    while (1) {

      // Wait until our DMA engine can support at least one more DMA
      while (!(ar_ni_status_get(my_cid) & 0xFF)) {
        ;
      }

      // Send to dst mailbox, keeping track with a local counter which will
      // notify our mailslot
      ar_cnt_set_notify_mslot(my_bid, my_cid, cnt, -NOC_MESSAGE_SIZE);
      ar_dma_with_ack(my_cid,
                      my_bid,  my_cid,  (unsigned int) context->noc_send_buf,
                      dst_bid, dst_cid, ar_addr_of_mbox(dst_bid, dst_cid),
                      my_bid,  my_cid,  cnt,
                      NOC_MESSAGE_SIZE, 0, 0, 0);


      // Block on mailslot until counter notification arrives
      ret = ar_mslot_get(my_cid);

      if (ret == 1) {
        // Success, counter triggered "ack"
        break;
      }
      else if (ret != 0) {
        // Unkown mailslot contents. Expected a "nack" value of 0.
        ar_abort();
      }

      // Got a "nack", because dst mailbox is full. Retry.
    }
  }
 

  // =========================================================================
  // Credit-based communication modes
  // =========================================================================
  else if ((context->noc_mode == NOC_MODE_CREDIT_MAILBOX) ||
           (context->noc_mode == NOC_MODE_CREDIT_ONLY))  {

    // Extra sanity checks
    ar_assert(context->noc_credits);
    ar_assert(context->noc_credits[dst_core_id]);

    // Wait until we have credits to send to this core
    do {

      // Credit counters wrap around. Their unsigned difference is always the
      // number of available credits, because the hw counter is always ahead.
      credits = ((unsigned int) ar_cnt_get(my_cid, 
                        context->noc_credits[dst_core_id]->cnt_send_credits)) -
                context->noc_credits[dst_core_id]->send_credits;

#ifdef NOC_WARN_OUT_OF_CREDITS
      // Warn once if we're out of credits. This may lead to deadlocks in
      // complex situations under heavy traffic where both sender and receiver
      // send multiple messages to each other and may run out of credits.
      // If encountered, suggested avoidance at the moment is to switch to
      // credit-only mode and user more credits per peer.
      if ((!credits) && (!context->noc_cred_warned)) {
        kt_printf("%d: WARNING: Out of credits encountered\r\n", 
                  context->pr_core_id);
        context->noc_cred_warned = 1;
      }
#endif
    } while (credits == 0);

    // Remember we're using one credit
    context->noc_credits[dst_core_id]->send_credits++;


    // Compute buffer address
    src_buf = ((unsigned int) 
                        context->noc_credits[dst_core_id]->loc_send_bufs) +
              context->noc_credits[dst_core_id]->cur_send_buf_id * 
              NOC_MESSAGE_SIZE;
    dst_buf = ((unsigned int) 
                        context->noc_credits[dst_core_id]->rem_send_bufs) +
              context->noc_credits[dst_core_id]->cur_send_buf_id * 
              NOC_MESSAGE_SIZE;

    // Take a trace
#ifdef DBG_TRC_COMM_ENABLED
    dbg_trace(context, (0x2 << 30) | 
                       (dst_core_id << 20) |
                       ((((PrMsgReq *) src_buf)->req_id) & 0xFFFFF));
#endif

    // Wait until our DMA engine can support at least one more DMA
    while (!(ar_ni_status_get(my_cid) & 0xFF)) {
      ;
    }

    // Send data to destination buffer and send the acknowledgments to the
    // remote hw data counter
    ar_dma_with_ack(my_cid,
                    my_bid,  my_cid,  src_buf,
                    dst_bid, dst_cid, dst_buf,
                    dst_bid, dst_cid, context->noc_credits[dst_core_id]->
                                                                cnt_send_data,
                    NOC_MESSAGE_SIZE, 0, 0, 0);

    // Are we using the mailbox for polling?
    if (context->noc_mode == NOC_MODE_CREDIT_MAILBOX) {

      // Prepare mailbox word
      word = context->pr_core_id << NOC_MBOX_COREID_OFFSET;
      word |= (context->noc_credits[dst_core_id]->cur_send_buf_id << 
                                                        NOC_MBOX_BUFID_OFFSET);

      // Wait until our DMA engine can support at least one more DMA
      while (!(ar_ni_status_get(my_cid) & 0xFF)) {
        ;
      }

      // Send the word to destination mailbox
      ar_mbox_send(my_cid, dst_bid, dst_cid, word);
    }

    // We used this buffer; next time go to the next one
    context->noc_credits[dst_core_id]->cur_send_buf_id++;
    if (context->noc_credits[dst_core_id]->cur_send_buf_id >= 
        context->noc_credits[dst_core_id]->num_send_bufs) {
      context->noc_credits[dst_core_id]->cur_send_buf_id = 0;
    }
  }


  // Unknown noc mode
  else {
    ar_abort();
  }


  // Success
  return 0;
}


// ===========================================================================
// noc_msg_recv()               Get a message from any core that has sent 
//                              one to us. Each message is exactly
//                              NOC_MESSAGE_SIZE bytes.
//
//                              In mailbox-only mode, the buffer is dequeued
//                              from the mailbox and copied to the single
//                              buffer in the context.
//
//                              In credit-based modes, we either learn that a
//                              new message has arrived in two ways: (i) when
//                              using mailbox polling, a single word is 
//                              dequeued which specifies the sender and the
//                              buffer it used, or (ii) when not using the 
//                              mailbox, we scan the local hw counters of all
//                              our peers in round-robin to learn if there's
//                              any new message transfer. In both cases, the
//                              DMA completion is waited upon (this delay
//                              should be trivial, because network guarantees
//                              in-order delivery -- it just guards against
//                              multi-port DRAM race conditions upon L2 cache
//                              rejections) and the buffer address is 
//                              returned.
// ===========================================================================
// * INPUTS
//   int block                  1: Block until there's a new message
//                              0: Return immediately if no new messages
//
// * OUTPUTS
//   void **ret_buf             The address of a buffer containing the
//                              received data will be placed on *ret_buf.
//                              WARNING: the buffer may be overwritten by
//                              subsequent calls to noc_msg_recv(), so don't
//                              rely on it.
//
// * RETURN VALUE
//   int                        0: success, new message in *ret_buf
//                              1: no new messages
// ===========================================================================
int noc_msg_recv(int block, void **ret_buf) {

  Context       *context;
  int           my_bid;
  int           my_cid;
  unsigned int  word;
  int           src_core_id;
  int           src_bid;
  int           src_cid;
  unsigned int  buf_id;
  unsigned int  bytes;
  int           i;
  int           j;
  int           nxt_j;


  // Get context and architectural board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Sanity checks
  ar_assert(ret_buf);

  // Non-blocking? We can return immediately for modes using the mailbox.
  if ((!block) && ((context->noc_mode == NOC_MODE_MAILBOX_ONLY) ||
                   (context->noc_mode == NOC_MODE_CREDIT_MAILBOX))) {

    // If nothing in mailbox, return empty-handed
    if ((ar_mbox_status_get(my_cid) & 0xFFFF) == 0) {
      *ret_buf = NULL;
      return 1;
    }
  }


  // =========================================================================
  // Mailbox-only communication mode
  // =========================================================================
  if (context->noc_mode == NOC_MODE_MAILBOX_ONLY) {

    // Extra checks
    ar_assert(context->noc_recv_buf);

    // Dequeue message and place it in the context buffer
    for (i = 0; i < NOC_MESSAGE_SIZE / 4; i++) {
      ((unsigned int *) context->noc_recv_buf)[i] = ar_mbox_get(my_cid);
    }

    // Return single buffer
    *ret_buf = context->noc_recv_buf;
  }


  // =========================================================================
  // Credit-based communication mode, no mailbox
  // =========================================================================
  else if (context->noc_mode == NOC_MODE_CREDIT_ONLY) {

    ar_assert(context->noc_credits_poll);

    // Polling loop. Loop continuously if we are in blocking mode, poll all
    // peers once and then exit if we don't want to block.
    src_core_id = -1;
    do {
      j = context->noc_poll_rr;
      for (i = 0; i < context->noc_num_peers; i++) {
        // Data counters wrap around. Their unsigned difference is always the
        // number of available data, because the hw counter is always ahead.
        bytes = ((unsigned int) ar_cnt_get(my_cid, 
                          context->noc_credits_poll[j]->cnt_recv_data)) -
                  context->noc_credits_poll[j]->recv_data;

        // Compute next position to poll
        nxt_j = j + 1;
        if (nxt_j >= context->noc_num_peers) {
          nxt_j = 0;
        }

        // Hit at current position?
        if (bytes >= NOC_MESSAGE_SIZE) {

          // Signal the hit outside the loop
          src_core_id = context->noc_credits_poll[j]->rem_core_id;

          // Remember next position to continue the round-robin
          context->noc_poll_rr = nxt_j;

          goto done;
        }

        j = nxt_j;
      }
    } while (block);

done:

    if ((!block) && (src_core_id == -1)) {
      *ret_buf = NULL;
      return 1;
    }


    // Discover source arch board/core ID
    pr_core_arch_bid_cid(src_core_id, &src_bid, &src_cid);

    // We handled this receive, so send a credit back. This is a DMA
    // operation for us, so make sure our DMA engine is not full.
    while (!(ar_ni_status_get(my_cid) & 0xFF)) {
      ;
    }
    ar_cnt_incr(my_cid, src_bid, src_cid, 
                context->noc_credits_poll[j]->cnt_recv_credits, 1);

    // Compute and store the buffer address
    *ret_buf = (void *) (
        ((unsigned int) context->noc_credits_poll[j]->recv_bufs) +
        context->noc_credits_poll[j]->cur_recv_buf_id * NOC_MESSAGE_SIZE);

    // Remember we've read so many bytes and we've used this buffer
    context->noc_credits_poll[j]->recv_data += NOC_MESSAGE_SIZE;
    context->noc_credits_poll[j]->cur_recv_buf_id++;
    if (context->noc_credits_poll[j]->cur_recv_buf_id >= 
        context->noc_credits_poll[j]->num_recv_bufs) {
      context->noc_credits_poll[j]->cur_recv_buf_id = 0;
    }
    
    // Take a trace
#ifdef DBG_TRC_COMM_ENABLED
    dbg_trace(context, (0x3 << 30) | 
                       (src_core_id << 20) |
                       ((((PrMsgReq *) *ret_buf)->req_id) & 0xFFFFF));
#endif
  }


  // =========================================================================
  // Credit-based communication mode, using mailbox for polling
  // =========================================================================
  else if (context->noc_mode == NOC_MODE_CREDIT_MAILBOX) {

    // Dequeue control word
    word = ar_mbox_get(my_cid);

    // Decode it
    src_core_id = (word & NOC_MBOX_COREID_MASK) >> NOC_MBOX_COREID_OFFSET;
    buf_id      = (word & NOC_MBOX_BUFID_MASK)  >> NOC_MBOX_BUFID_OFFSET;

    // Sanity check: verify we should be talking to src_core_id
    ar_assert((src_core_id >= 0) && (src_core_id < context->pr_num_cores));
    ar_assert(context->noc_credits);
    ar_assert(context->noc_credits[src_core_id]);
    ar_assert(buf_id < context->noc_credits[src_core_id]->num_recv_bufs);

    // Discover source arch board/core ID
    pr_core_arch_bid_cid(src_core_id, &src_bid, &src_cid);

    // We dequeued the mailbox word, so send a credit back. This is a DMA
    // operation for us, so make sure our DMA engine is not full.
    while (!(ar_ni_status_get(my_cid) & 0xFF)) {
      ;
    }
    ar_cnt_incr(my_cid, src_bid, src_cid, 
                context->noc_credits[src_core_id]->cnt_recv_credits, 1);

    // Compute and store the buffer address
    *ret_buf = (void *) (((unsigned int) 
                          context->noc_credits[src_core_id]->recv_bufs) +
                         buf_id * NOC_MESSAGE_SIZE);

    // Make sure the buffer DMA is completed
    do {

      // Data counters wrap around. Their unsigned difference is always the
      // number of available data, because the hw counter is always ahead.
      bytes = ((unsigned int) ar_cnt_get(my_cid, 
                        context->noc_credits[src_core_id]->cnt_recv_data)) -
                context->noc_credits[src_core_id]->recv_data;
    } while (bytes < NOC_MESSAGE_SIZE);

    // Remember we've read so many bytes
    context->noc_credits[src_core_id]->recv_data += NOC_MESSAGE_SIZE;
    
    // Take a trace
#ifdef DBG_TRC_COMM_ENABLED
    dbg_trace(context, (0x3 << 30) | 
                       (src_core_id << 20) |
                       ((((PrMsgReq *) *ret_buf)->req_id) & 0xFFFFF));
#endif
  }


  // Unknown noc mode
  else {
    ar_abort();
  }

  // Take a note for the statistics gathering
  dbg_scount(context, DBG_STATS_IDX_NUM_MESSAGES, 1);

  // Success
  return 0;
}

