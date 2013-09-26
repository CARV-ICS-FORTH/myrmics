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
// Abstract      : Network-on-chip initialization
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: init.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/12/18 16:40:29 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>


// Magic command words used during initialization
#define NOC_MAGIC_INIT_TELL             0xDEAD00AA
#define NOC_MAGIC_PEER_TELL             0xDEAD00BB
#define NOC_MAGIC_PEER_LISTEN           0xDEAD00CC
#define NOC_MAGIC_FINISHED              0xDEAD00DD


// ===========================================================================
// noc_init()                   Network-on-chip layer initialization.
//                              Initializes per-peer counters, credits and
//                              memory slots (if applicable).
//
//                              Must be called after pr_init(), because it
//                              depends on parent/children pr layer arrays to
//                              determine the number of peers.
// ===========================================================================
// * INPUTS
//   NocMode mode               NOC_MODE_MAILBOX_ONLY: 64-B messages are
//                                exchanged through mailboxes. Sender waits
//                                until receiver mailbox has received it and
//                                retries if it gets Nacked. Slowest mode,
//                                but can be used for large fanouts of
//                                scheduler-workers, because it does not
//                                dedicate counters to peers.
//
//                              NOC_MODE_CREDIT_ONLY: 64-B messages are
//                                directly written to per-peer memory slots and
//                                DMAed to the other side. 2 counters per peer
//                                are dedicated for incoming data count and
//                                outgoing credits. Receiver also uses them to
//                                poll if there's an incoming message.  Fast
//                                sends, but receives are polling-based using
//                                the hardware counters (polling latency
//                                increases when number of peers increase).
//                                Can use arbitrarily large number of credits,
//                                by allocating send/recv buffers as needed
//                                during initialization.  It cannot scale to
//                                large fanouts of scheduler- workers, because
//                                it dedicates 2 counters per peer.
//
//                              NOC_MODE_CREDIT_MAILBOX: 64-B messages are
//                                directly written to per-peer memory slots and
//                                mailboxes are used for ordering.  2 counters
//                                per peer are dedicated for incoming data
//                                count and outgoing credits. Fastest mode,
//                                because it does not wait almost nowhere,
//                                mailbox space is guaranteed to fit and used
//                                to immediately discover the next buffer to
//                                receive without polling all peer counters. It
//                                cannot scale to large fanouts of scheduler-
//                                workers, because it dedicates 2 counters per
//                                peer. It also has limited credits per peer,
//                                due to the partitioned mailbox space.
//
//   int credit_only_credits    For NOC_MODE_CREDIT_ONLY mode, specify the
//                              credits (and buffers) per peer to allocate.
// ===========================================================================
void noc_init(NocMode mode, int credit_only_credits) {

  Context       *context;
  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  ListNode      *cnt_node;
  int           num_peers;
  unsigned int  max_credits;
  int           peers_done;
  List          *peer_list;
  ListNode      *peer_list_node;
  int           tmp_cmd;
  int           tmp_peers;
  int           tmp_credits;
  int           tmp_core_id;
  int           tmp_bid;
  int           tmp_cid;
  int           tmp_cnt_credit;
  int           tmp_cnt_data;
  int           tmp_buf_base;
  int           tmp_num_bufs;
  int           i;
  int           j;


  // Get context and architectural board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Sanity checks
  ar_assert(context->pr_core_id > -1);
  ar_assert(!(NOC_MESSAGE_SIZE & 0x3F));

  // Remember noc mode
  context->noc_mode = mode;

  // Don't proceed if we are in single-core mode
  if (context->pr_num_cores == 1) {
    return;
  }

  // Create the free list of all available counters. The counter number
  // is used as the generic data pointer.
  context->noc_cnt_free = kt_alloc_list();
  for (i = 0; i < NOC_MAX_COUNTERS; i++) {
    kt_list_insert(context->noc_cnt_free, (void *) i, 
                   kt_list_tail(context->noc_cnt_free), 1);
  }

  // Create the active DMAs list and leave it empty
  context->noc_active_dmas = kt_alloc_list();


  // Mode-dependent init: Mailbox-only mode
  if (mode == NOC_MODE_MAILBOX_ONLY) {
    
    // Allocate a single buffer for reception and a single for sending
    context->noc_send_buf = kt_malloc(NOC_MESSAGE_SIZE);
    context->noc_recv_buf = kt_malloc(NOC_MESSAGE_SIZE);

    // Nothing else is valid in this mode
    ar_assert(!context->noc_credits);
    ar_assert(!context->noc_credits_poll);
  }

  // Mode-dependent init: Credit-based modes
  else if ((mode == NOC_MODE_CREDIT_ONLY) || 
           (mode == NOC_MODE_CREDIT_MAILBOX)) {

    // No single buffers in this mode
    ar_assert(!context->noc_send_buf);
    ar_assert(!context->noc_recv_buf);

    // Allocate a credit structure for every possible core ID
    context->noc_credits = kt_malloc(context->pr_num_cores * 
                                     sizeof(NocCreditMode *));

    // Discover our peers
    num_peers = 0;
    peer_list = kt_alloc_list();
    for (i = 0; i < context->pr_num_cores; i++) {

      // Default is no communication with this core ID
      context->noc_credits[i] = NULL;

      // Maybe this core ID is our parent scheduler?
      if ((context->pr_parent_sched_id != -1) &&
          (pr_scheduler_core_id(context->pr_parent_sched_id) == i)) {
        context->noc_credits[i] = kt_malloc(sizeof(NocCreditMode));
      }

      // Is this core ID one of our children?
      else if (context->pr_role == ROLE_SCHEDULER) {
        for (j = 0; j < context->pr_num_children; j++) {
          if (context->pr_children[j] == i) {
            context->noc_credits[i] = kt_malloc(sizeof(NocCreditMode));
            break;
          }
        }
      }

      // If we don't communicate with this core ID, go to the next one
      if (!context->noc_credits[i]) {
        continue;
      }

      // We found one more
      num_peers++;

      // Fill out its core ID
      context->noc_credits[i]->rem_core_id = i;

      // Assign a local hardware counter for receiving data acks from this peer
      cnt_node = kt_list_head(context->noc_cnt_free);
      ar_assert(cnt_node); // if this fails, there aren't enough hw counters
                           // to use credit-based communication for this
                           // scheduler-to-worker or scheduler-to-scheduler
                           // ratios. Switch to mailbox-only communication.
      context->noc_credits[i]->cnt_recv_data = (int) cnt_node->data;
      ar_cnt_set(my_cid, (int) cnt_node->data, 0);
      context->noc_credits[i]->recv_data = 0;
      kt_list_delete(context->noc_cnt_free, cnt_node, NULL);

      // Assign a local hardware counter for credits when sending to this peer
      cnt_node = kt_list_head(context->noc_cnt_free);
      ar_assert(cnt_node); // if this fails, same problem as above
      context->noc_credits[i]->cnt_send_credits = (int) cnt_node->data;
      ar_cnt_set(my_cid, (int) cnt_node->data, 0); // will be replaced below
      context->noc_credits[i]->send_credits = 0;
      kt_list_delete(context->noc_cnt_free, cnt_node, NULL);

      // No remote counters yet
      context->noc_credits[i]->cnt_send_data = -1;
      context->noc_credits[i]->cnt_recv_credits = -1;

      // No local/remote buffers yet
      context->noc_credits[i]->loc_send_bufs = NULL;
      context->noc_credits[i]->rem_send_bufs = NULL;
      context->noc_credits[i]->num_send_bufs = 0;
      context->noc_credits[i]->cur_send_buf_id = 0;
      context->noc_credits[i]->recv_bufs = NULL;
      context->noc_credits[i]->num_recv_bufs = 0;
      context->noc_credits[i]->cur_recv_buf_id = 0;

      // Add it to a list for future gap filling
      kt_list_insert(peer_list, (void *) i, kt_list_tail(peer_list), 1);
    }


    // Sanity check: number of peers should be our children plus our parent
    ar_assert(num_peers == context->pr_num_children + 
                           ((context->pr_parent_sched_id != -1) ? 1 : 0));
    context->noc_num_peers = num_peers;

    // In credit-only mode, create and fill out a second, condensed array
    // to be used for polling. This is a pointer-linked array, so whatever
    // changes we do below to fill out the first array are visible also to 
    // the second one.
    if (mode == NOC_MODE_CREDIT_ONLY) {
      context->noc_credits_poll = kt_malloc(context->noc_num_peers * 
                                       sizeof(NocCreditMode *));
      j = 0;
      for (i = 0; i < context->pr_num_cores; i++) {
        if (context->noc_credits[i]) {
          context->noc_credits_poll[j++] = context->noc_credits[i];
        }
      }
      ar_assert(j == context->noc_num_peers);
    }

    // Specify how many credits we'll give to each of our peers when he
    // needs to send us messages
    if (mode == NOC_MODE_CREDIT_MAILBOX) {
      ar_uint_divide(AR_MBOX_SIZE / 4, num_peers, &max_credits, NULL);
      ar_assert(max_credits > 0);
    }
    else if (mode == NOC_MODE_CREDIT_ONLY) {
      max_credits = credit_only_credits;
    }
    else {
      ar_abort();
    }
    ar_assert(max_credits > 0);
    

    // Core #0 acts as a hub to exchange information without messing up
    // who relays what to whom
    peers_done = 0;
    peer_list_node = kt_list_head(peer_list);
    if (!context->pr_core_id) {

      // For all cores
      for (i = 0; i < context->pr_num_cores; i++) {

        // If not myself...
        if (i) {

          // Get arch-level board/core IDs
          pr_core_arch_bid_cid(i, &bid, &cid);

          // Command core to tell us about its peers
          ar_mbox_send(my_cid, bid, cid, NOC_MAGIC_INIT_TELL);

          // First words are how many peers he has and how many credits he 
          // gives to each of them
          tmp_peers   = ar_mbox_get(my_cid);
          tmp_credits = ar_mbox_get(my_cid);
        }
        
        // myself
        else {
          tmp_peers = num_peers;
          tmp_credits = max_credits;
        }


        // For all his peers
        for (j = 0; j < tmp_peers; j++) {

          // If not myself...
          if (i) {

            // Command core to tell us about its next peer
            ar_mbox_send(my_cid, bid, cid, NOC_MAGIC_PEER_TELL);

            // Get per-peer words
            tmp_core_id    = ar_mbox_get(my_cid);
            tmp_cnt_credit = ar_mbox_get(my_cid);
            tmp_cnt_data   = ar_mbox_get(my_cid);
            tmp_buf_base   = ar_mbox_get(my_cid);
            tmp_num_bufs   = ar_mbox_get(my_cid);
          }

          // myself
          else {
            ar_assert(peers_done < num_peers);

            // Get next peer from the list
            tmp_core_id = (int) peer_list_node->data;
            ar_assert(tmp_core_id >= 0);
            ar_assert(tmp_core_id < context->pr_num_cores);
            ar_assert(tmp_core_id != context->pr_core_id);
            peer_list_node = kt_list_next(peer_list_node);

            // Assign buffers, now we know how many credits we're giving out.
            // We are using (credits + 2) buffers, because up to (credits)
            // buffers may be sent by a sender, while the receiver has returned
            // another one to the caller of noc_msg_recv() and while the 
            // sender is filling up another one returned by calling
            // noc_msg_get_send_buf().
            context->noc_credits[tmp_core_id]->num_recv_bufs = max_credits + 2;
            context->noc_credits[tmp_core_id]->recv_bufs = kt_malloc(
                            context->noc_credits[tmp_core_id]->num_recv_bufs *
                            NOC_MESSAGE_SIZE);

            // Assign info
            tmp_cnt_credit = 
                         context->noc_credits[tmp_core_id]->cnt_send_credits;
            tmp_cnt_data = 
                         context->noc_credits[tmp_core_id]->cnt_recv_data;
            tmp_buf_base = (int)
                         context->noc_credits[tmp_core_id]->recv_bufs;
            tmp_num_bufs = 
                         context->noc_credits[tmp_core_id]->num_recv_bufs;

            // Remember we did one more
            peers_done++;
          }

          // Find out arch-level board/core IDs of its peer
          pr_core_arch_bid_cid(tmp_core_id, &tmp_bid, &tmp_cid);

          // Relay info to him...
          if (tmp_core_id) {
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, NOC_MAGIC_PEER_LISTEN);
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, i);
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, tmp_credits);
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, tmp_cnt_credit);
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, tmp_cnt_data);
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, tmp_buf_base);
            ar_mbox_send(my_cid, tmp_bid, tmp_cid, tmp_num_bufs);
          }

          // ... or directly to us.
          else {
            context->noc_credits[i]->cnt_recv_credits = tmp_cnt_credit;
            context->noc_credits[i]->cnt_send_data    = tmp_cnt_data;
            context->noc_credits[i]->rem_send_bufs    = (void *) tmp_buf_base;
            context->noc_credits[i]->num_send_bufs    = tmp_num_bufs;

            // Allocate the same number of local send buffers
            context->noc_credits[i]->loc_send_bufs = kt_malloc(tmp_num_bufs * 
                                                NOC_MESSAGE_SIZE);

            // Set correct credit value
            ar_cnt_set(my_cid, 
                       context->noc_credits[i]->cnt_send_credits, tmp_credits);
          }
        }
      }

      // Indicate to everyone we're done. We don't rely only on the mailbox in
      // order to avoid cores leaving this phase and sending normal mailbox
      // messages to others who are still waiting to leave this barrier.
      for (i = 1; i < context->pr_num_cores; i++) {
        pr_core_arch_bid_cid(i, &bid, &cid);
        ar_mbox_send(my_cid, bid, cid, NOC_MAGIC_FINISHED);
      }
      for (i = 1; i < context->pr_num_cores; i++) {
        pr_core_arch_bid_cid(i, &bid, &cid);
        ar_cnt_incr(my_cid, bid, cid, NOC_COUNTER_WAKEUP0, 1);
      }
    }

    // Other cores are slaves
    else {

      // Initialize our barrier counter to -1
      ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP0, -1);

      // Get master arch-level board/core IDs
      pr_core_arch_bid_cid(0, &bid, &cid);

      // Wait for a new command from master
      do {

        // Get command
        tmp_cmd = ar_mbox_get(my_cid);

        switch (tmp_cmd) {

          // Tell number of peers and credits
          case NOC_MAGIC_INIT_TELL:
            ar_mbox_send(my_cid, bid, cid, num_peers);
            ar_mbox_send(my_cid, bid, cid, max_credits);
            break;


          // Tell about next peer
          case NOC_MAGIC_PEER_TELL:
            ar_assert(peers_done < num_peers);

            // Get next peer from the list
            tmp_core_id = (int) peer_list_node->data;
            ar_assert(tmp_core_id >= 0);
            ar_assert(tmp_core_id < context->pr_num_cores);
            ar_assert(tmp_core_id != context->pr_core_id);
            peer_list_node = kt_list_next(peer_list_node);

            // Assign buffers, now we know how many credits we're giving out.
            // See above why max_credits + 2.
            context->noc_credits[tmp_core_id]->num_recv_bufs = max_credits + 2;
            context->noc_credits[tmp_core_id]->recv_bufs = kt_malloc(
                            context->noc_credits[tmp_core_id]->num_recv_bufs *
                            NOC_MESSAGE_SIZE);

            // Send info
            ar_mbox_send(my_cid, bid, cid, tmp_core_id);
            ar_mbox_send(my_cid, bid, cid, 
                         context->noc_credits[tmp_core_id]->cnt_send_credits);
            ar_mbox_send(my_cid, bid, cid,
                         context->noc_credits[tmp_core_id]->cnt_recv_data);
            ar_mbox_send(my_cid, bid, cid, (int)
                         context->noc_credits[tmp_core_id]->recv_bufs);
            ar_mbox_send(my_cid, bid, cid, 
                         context->noc_credits[tmp_core_id]->num_recv_bufs);

            // Remember we did one more
            peers_done++;
            break;


          // Discover info about a peer
          case NOC_MAGIC_PEER_LISTEN:

            // Get info
            i              = ar_mbox_get(my_cid);
            tmp_credits    = ar_mbox_get(my_cid);
            tmp_cnt_credit = ar_mbox_get(my_cid);
            tmp_cnt_data   = ar_mbox_get(my_cid);
            tmp_buf_base   = ar_mbox_get(my_cid);
            tmp_num_bufs   = ar_mbox_get(my_cid);

            // Fill it in our array
            ar_assert(i >= 0);
            ar_assert(i < context->pr_num_cores);

            context->noc_credits[i]->cnt_recv_credits = tmp_cnt_credit;
            context->noc_credits[i]->cnt_send_data    = tmp_cnt_data;
            context->noc_credits[i]->rem_send_bufs    = (void *) tmp_buf_base;
            context->noc_credits[i]->num_send_bufs    = tmp_num_bufs;

            // Allocate the same number of local send buffers
            context->noc_credits[i]->loc_send_bufs = kt_malloc(tmp_num_bufs * 
                                                NOC_MESSAGE_SIZE);

            // Set correct credit value
            ar_cnt_set(my_cid, 
                       context->noc_credits[i]->cnt_send_credits, tmp_credits);

            break;

          // Finished: the while loop will terminate
          case NOC_MAGIC_FINISHED:
            break;


          // Unknown command
          default:
            ar_abort();
        }

      } while (tmp_cmd != NOC_MAGIC_FINISHED);


      // Block until master says we're done
      while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP0)) {
        ;
      }
    }

    // Sanity check: we must have done all our peers
    ar_assert(peers_done == num_peers);

    // Destroy peer list
    kt_free_list(peer_list, NULL);
  }

  // Unknown mode
  else {
    ar_abort();
  }
}

