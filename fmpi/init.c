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
// Abstract      : Minimal MPI library: initialization
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: init.c,v $
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


// Magic command words used during initialization
#define FMPI_MAGIC_TELL                 0xDEAD0AA1
#define FMPI_MAGIC_LISTEN               0xDEAD0BB1
#define FMPI_MAGIC_FINISHED             0xDEAD0CC1


// ===========================================================================
// fmpi_rank_init()             Rank initialization. Finds out which ranks
//                              are mapped to which cores and which cores
//                              are active.
// ===========================================================================
// * INPUTS
//   int num_ranks              Number of total active cores in the system.
//                              Two arguments per core must be supplied to
//                              pr_init(), as specified below. The ranks will
//                              be created in order, i.e. the first two
//                              arguments will refer to rank = 0, the next two
//                              to rank = 1 and the last two to rank =
//                              (num_ranks - 1).
//
//   For each rank, specify two arguments:
//
//   int bid                    Architecture-level board ID of core
//   int cid                    Architecture-level core ID of core
//
// * RETURN VALUE
//   int                        0: success, and this core is part of the
//                                 running setup
//                              1: this core is not part of the setup and
//                                 should not continue running
// ===========================================================================
int fmpi_rank_init(int num_ranks, ...) {

  va_list       ap;
  Context       *context;
  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  FMPI_Context  *tmp_context;
  int           tmp_cmd;
  int           tmp_bid;
  int           tmp_cid;
  int           i;
  int           j;


  // Get global context and board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Allocate FMPI context
  ar_assert(!context->fmpi);
  context->fmpi = kt_zalloc(sizeof(FMPI_Context));

  // Assign the rank-related parts
  context->fmpi->num_ranks = num_ranks;
  context->fmpi->rank = -1; // unknown yet
  context->fmpi->rank_bid_cid = kt_malloc(num_ranks * sizeof(int));
  context->fmpi->rank_context = kt_malloc(num_ranks * sizeof(FMPI_Context *));

  // Init them to dummy values, to catch errors
  for (i = 0; i < num_ranks; i++) {
    context->fmpi->rank_bid_cid[i] = -1;
    context->fmpi->rank_context[i] = NULL;
  }

  // Single-core hack
  if (num_ranks == 1) {
    
    // Just get arch board/core ID
    va_start(ap, num_ranks);
    bid  = va_arg(ap, int);
    cid  = va_arg(ap, int);
    va_end(ap);

    // See if we are this core
    if ((my_bid == bid) && (my_cid == cid)) {

      context->fmpi->rank = 0;
      context->fmpi->rank_bid_cid[0] = (bid << 8) | cid;
      context->fmpi->rank_context[0] = context->fmpi;

      return 0;
    }
    else {
      // Fail cores that don't match
      return 1;
    }
  }


  // Initialize variable argument list
  va_start(ap, num_ranks);

  // Process list to assign ranks
  for (i = 0; i < num_ranks; i++) {

    // Get four arguments for this core ID
    bid  = va_arg(ap, int);
    cid  = va_arg(ap, int);

    // Are we this core ID?
    if ((bid == my_bid) && (cid == my_cid)) {
      context->fmpi->rank = i;
      context->fmpi->rank_context[i] = context->fmpi;
    }

    // Fill out board/core ID translation array
    ar_assert(context->fmpi->rank_bid_cid[i] == -1);
    context->fmpi->rank_bid_cid[i] = (bid << 8) | cid;
  }

  // Close list
  va_end(ap);

  // If we are not in the list, we're not part of the FMPI setup: quit.
  if (context->fmpi->rank == -1) {
    kt_free(context->fmpi->rank_bid_cid);
    kt_free(context->fmpi->rank_context);
    context->fmpi->rank_bid_cid = NULL;
    context->fmpi->rank_context = NULL;
    return 1;
  }


  // =========================================================================
  // Phase 2: All ranks must learn where is the FMPI context of all other
  //          ranks, so they know where to fetch descriptors from.
  // =========================================================================

  // Rank #0 acts as a hub to exchange information without messing up
  // who relays what to whom
  if (!context->fmpi->rank) {

    // For all ranks
    for (i = 0; i < context->fmpi->num_ranks; i++) {

      // If not myself...
      if (i) {

        // Get arch-level board/core IDs
        fmpi_get_bid_cid(context->fmpi, i, &bid, &cid);

        // Command core to tell us its context address
        ar_mbox_send(my_cid, bid, cid, FMPI_MAGIC_TELL);

        // Remember it
        tmp_context = (FMPI_Context *) ar_mbox_get(my_cid);
      }
      
      // myself
      else {
        tmp_context = context->fmpi;
      }


      // For all ranks
      for (j = 0; j < context->fmpi->num_ranks; j++) {

        // Find out arch-level board/core IDs of its peer
        fmpi_get_bid_cid(context->fmpi, j, &tmp_bid, &tmp_cid);

        // Relay context to him...
        if (j) {
          ar_mbox_send(my_cid, tmp_bid, tmp_cid, FMPI_MAGIC_LISTEN);
          ar_mbox_send(my_cid, tmp_bid, tmp_cid, i);
          ar_mbox_send(my_cid, tmp_bid, tmp_cid, (unsigned int) tmp_context);
        }

        // ... or directly to us.
        else {
          context->fmpi->rank_context[i] = tmp_context;
        }
      }
    }

    // Indicate to everyone we're done. We don't rely only on the mailbox in
    // order to avoid cores leaving this phase and sending normal mailbox
    // messages to others who are still waiting to leave this barrier.
    for (i = 1; i < num_ranks; i++) {
      fmpi_get_bid_cid(context->fmpi, i, &bid, &cid);
      ar_mbox_send(my_cid, bid, cid, FMPI_MAGIC_FINISHED);
    }
    for (i = 1; i < num_ranks; i++) {
      fmpi_get_bid_cid(context->fmpi, i, &bid, &cid);
      ar_cnt_incr(my_cid, bid, cid, NOC_COUNTER_WAKEUP0, 1);
    }
  }

  // Other cores are slaves
  else {

    // Initialize our barrier counter to -1
    ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP0, -1);

    // Get master arch-level board/core IDs
    fmpi_get_bid_cid(context->fmpi, 0, &bid, &cid);

    // Wait for a new command from master
    do {

      // Get command
      tmp_cmd = ar_mbox_get(my_cid);

      switch (tmp_cmd) {

        // Tell our context
        case FMPI_MAGIC_TELL:
          ar_mbox_send(my_cid, bid, cid, (unsigned int) context->fmpi);
          break;


        // Discover context of a peer
        case FMPI_MAGIC_LISTEN:

          // Get info
          i           = ar_mbox_get(my_cid);
          tmp_context = (FMPI_Context *) ar_mbox_get(my_cid);

          // Fill it in our array
          ar_assert(i >= 0);
          ar_assert(i < context->fmpi->num_ranks);

          context->fmpi->rank_context[i] = tmp_context;

          break;

        // Finished: the while loop will terminate
        case FMPI_MAGIC_FINISHED:
          break;


        // Unknown command
        default:
          ar_abort();
      }

    } while (tmp_cmd != FMPI_MAGIC_FINISHED);


    // Block until master says we're done
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP0)) {
      ;
    }
  }

  // Sanity check: we must have a valid context address for all ranks
  for (i = 0; i < context->fmpi->num_ranks; i++) {
    ar_assert(context->fmpi->rank_context[i]);
  }

  // Success
  return 0;
}


// ===========================================================================
// fmpi_select_rank_init()      Calls the real initialization function with
//                              the correct core setup
// ===========================================================================
// * INPUTS
//   FMPI_CfgMode select        The core setup
//
// * RETURN VALUE
//   int                        fmpi_rank_init() return value
// ===========================================================================
int fmpi_select_rank_init(FMPI_CfgMode select) {

  switch (select) {

    // ======================================================================
    
    case FMPI_CFG_1_MB:
      return fmpi_rank_init(1,
                            0x00, 0);   // 0

    case FMPI_CFG_1_ARM:
      return fmpi_rank_init(1,
                            0x6B, 0);   // 0

    // ======================================================================

    case FMPI_CFG_2_MB:
      return fmpi_rank_init(2,
                            0x00, 0,    // 0
                            0x00, 1);   // 1

    case FMPI_CFG_2_ARM:
      return fmpi_rank_init(2,
                            0x6B, 0,    // 0
                            0x6B, 1);   // 1

    case FMPI_CFG_2_HET:
      return fmpi_rank_init(2,
                            0x00, 0,    // 0
                            0x6B, 0);   // 1

    // ======================================================================

    case FMPI_CFG_4_MB:
      return fmpi_rank_init(4,
                            0x00, 0,    // 0
                            0x00, 1,    // 1
                            0x00, 2,    // 2
                            0x00, 3);   // 3

    case FMPI_CFG_4_ARM:
      return fmpi_rank_init(4,
                            0x6B, 0,    // 0
                            0x6B, 1,    // 1
                            0x6B, 2,    // 2
                            0x6B, 3);   // 3

    case FMPI_CFG_4_HET:
      return fmpi_rank_init(4,
                            0x00, 0,    // 0
                            0x00, 1,    // 1
                            0x00, 2,    // 2
                            0x6B, 0);   // 3

    // ======================================================================

    case FMPI_CFG_8_MB:
      return fmpi_rank_init(8,
                            0x00, 0,    // 0
                            0x00, 1,    // 1
                            0x00, 2,    // 2
                            0x00, 3,    // 3
                            0x00, 4,    // 4
                            0x00, 5,    // 5
                            0x00, 6,    // 6
                            0x00, 7);   // 7

    case FMPI_CFG_8_HET:
      return fmpi_rank_init(8,
                            0x2B, 0,    // 0
                            0x2B, 1,    // 1
                            0x2B, 2,    // 2
                            0x2B, 3,    // 3
                            0x6B, 0,    // 4
                            0x6B, 1,    // 5
                            0x6B, 2,    // 6
                            0x6B, 3);   // 7

    // ======================================================================

    case FMPI_CFG_12_HET:
      return fmpi_rank_init(12,
                            0x2B, 0,    // 0
                            0x2B, 1,    // 1
                            0x2B, 2,    // 2
                            0x2B, 3,    // 3
                            0x2B, 4,    // 4
                            0x2B, 5,    // 5
                            0x2B, 6,    // 6
                            0x2B, 7,    // 7
                            0x6B, 0,    // 8
                            0x6B, 1,    // 9
                            0x6B, 2,    // 10
                            0x6B, 3);   // 11

    // ======================================================================

    case FMPI_CFG_16_MB:
      return fmpi_rank_init(16,
                            0x00,0, 0x00,1, 0x00,2, 0x00,3, 
                            0x00,4, 0x00,5, 0x00,6, 0x00,7,  // 0-7
                            0x10,0, 0x10,1, 0x10,2, 0x10,3, 
                            0x10,4, 0x10,5, 0x10,6, 0x10,7); // 8-15

    // ======================================================================

    case FMPI_CFG_20_HET:
      return fmpi_rank_init(20,
                            0x2B, 0,    // 0
                            0x2B, 1,    // 1
                            0x2B, 2,    // 2
                            0x2B, 3,    // 3
                            0x2B, 4,    // 4
                            0x2B, 5,    // 5
                            0x2B, 6,    // 6
                            0x2B, 7,    // 7
                            0x1B, 0,    // 8
                            0x1B, 1,    // 9
                            0x1B, 2,    // 10
                            0x1B, 3,    // 11
                            0x1B, 4,    // 12
                            0x1B, 5,    // 13
                            0x1B, 6,    // 14
                            0x1B, 7,    // 15
                            0x6B, 0,    // 16
                            0x6B, 1,    // 17
                            0x6B, 2,    // 18
                            0x6B, 3);   // 19

    // ======================================================================

    case FMPI_CFG_64_MB:
      return fmpi_rank_init(64,
                            0x00,0, 0x00,1, 0x00,2, 0x00,3, 
                            0x00,4, 0x00,5, 0x00,6, 0x00,7,  // 0-7
                            0x10,0, 0x10,1, 0x10,2, 0x10,3, 
                            0x10,4, 0x10,5, 0x10,6, 0x10,7,  // 8-15
                            0x11,0, 0x11,1, 0x11,2, 0x11,3, 
                            0x11,4, 0x11,5, 0x11,6, 0x11,7,  // 16-23
                            0x01,0, 0x01,1, 0x01,2, 0x01,3, 
                            0x01,4, 0x01,5, 0x01,6, 0x01,7,  // 24-31
                            0x05,0, 0x05,1, 0x05,2, 0x05,3, 
                            0x05,4, 0x05,5, 0x05,6, 0x05,7,  // 32-39
                            0x15,0, 0x15,1, 0x15,2, 0x15,3, 
                            0x15,4, 0x15,5, 0x15,6, 0x15,7,  // 40-47
                            0x14,0, 0x14,1, 0x14,2, 0x14,3, 
                            0x14,4, 0x14,5, 0x14,6, 0x14,7,  // 48-55
                            0x04,0, 0x04,1, 0x04,2, 0x04,3, 
                            0x04,4, 0x04,5, 0x04,6, 0x04,7); // 56-63

    // ======================================================================

    case FMPI_CFG_128_MB:
      return fmpi_rank_init(128,
                            0x00,0, 0x00,1, 0x00,2, 0x00,3, 
                            0x00,4, 0x00,5, 0x00,6, 0x00,7,  // 0-7
                            0x10,0, 0x10,1, 0x10,2, 0x10,3, 
                            0x10,4, 0x10,5, 0x10,6, 0x10,7,  // 8-15
                            0x11,0, 0x11,1, 0x11,2, 0x11,3, 
                            0x11,4, 0x11,5, 0x11,6, 0x11,7,  // 16-23
                            0x01,0, 0x01,1, 0x01,2, 0x01,3, 
                            0x01,4, 0x01,5, 0x01,6, 0x01,7,  // 24-31
                            0x05,0, 0x05,1, 0x05,2, 0x05,3, 
                            0x05,4, 0x05,5, 0x05,6, 0x05,7,  // 32-39
                            0x15,0, 0x15,1, 0x15,2, 0x15,3, 
                            0x15,4, 0x15,5, 0x15,6, 0x15,7,  // 40-47
                            0x14,0, 0x14,1, 0x14,2, 0x14,3, 
                            0x14,4, 0x14,5, 0x14,6, 0x14,7,  // 48-55
                            0x04,0, 0x04,1, 0x04,2, 0x04,3, 
                            0x04,4, 0x04,5, 0x04,6, 0x04,7,  // 56-63

                            0x08,0, 0x08,1, 0x08,2, 0x08,3, 
                            0x08,4, 0x08,5, 0x08,6, 0x08,7,  // 64-71
                            0x18,0, 0x18,1, 0x18,2, 0x18,3, 
                            0x18,4, 0x18,5, 0x18,6, 0x18,7,  // 72-79
                            0x19,0, 0x19,1, 0x19,2, 0x19,3, 
                            0x19,4, 0x19,5, 0x19,6, 0x19,7,  // 80-87
                            0x09,0, 0x09,1, 0x09,2, 0x09,3, 
                            0x09,4, 0x09,5, 0x09,6, 0x09,7,  // 88-95
                            0x0D,0, 0x0D,1, 0x0D,2, 0x0D,3, 
                            0x0D,4, 0x0D,5, 0x0D,6, 0x0D,7,  // 96-103
                            0x0C,0, 0x0C,1, 0x0C,2, 0x0C,3, 
                            0x0C,4, 0x0C,5, 0x0C,6, 0x0C,7,  // 104-111
                            0x1C,0, 0x1C,1, 0x1C,2, 0x1C,3, 
                            0x1C,4, 0x1C,5, 0x1C,6, 0x1C,7,  // 112-119
                            0x1D,0, 0x1D,1, 0x1D,2, 0x1D,3, 
                            0x1D,4, 0x1D,5, 0x1D,6, 0x1D,7); // 120-127


    // ======================================================================

    case FMPI_CFG_256_MB:
      return fmpi_rank_init(256,
                            0x00,0, 0x00,1, 0x00,2, 0x00,3, 
                            0x00,4, 0x00,5, 0x00,6, 0x00,7,  // 0-7
                            0x10,0, 0x10,1, 0x10,2, 0x10,3, 
                            0x10,4, 0x10,5, 0x10,6, 0x10,7,  // 8-15
                            0x11,0, 0x11,1, 0x11,2, 0x11,3, 
                            0x11,4, 0x11,5, 0x11,6, 0x11,7,  // 16-23
                            0x01,0, 0x01,1, 0x01,2, 0x01,3, 
                            0x01,4, 0x01,5, 0x01,6, 0x01,7,  // 24-31
                            0x05,0, 0x05,1, 0x05,2, 0x05,3, 
                            0x05,4, 0x05,5, 0x05,6, 0x05,7,  // 32-39
                            0x15,0, 0x15,1, 0x15,2, 0x15,3, 
                            0x15,4, 0x15,5, 0x15,6, 0x15,7,  // 40-47
                            0x14,0, 0x14,1, 0x14,2, 0x14,3, 
                            0x14,4, 0x14,5, 0x14,6, 0x14,7,  // 48-55
                            0x04,0, 0x04,1, 0x04,2, 0x04,3, 
                            0x04,4, 0x04,5, 0x04,6, 0x04,7,  // 56-63

                            0x08,0, 0x08,1, 0x08,2, 0x08,3, 
                            0x08,4, 0x08,5, 0x08,6, 0x08,7,  // 64-71
                            0x18,0, 0x18,1, 0x18,2, 0x18,3, 
                            0x18,4, 0x18,5, 0x18,6, 0x18,7,  // 72-79
                            0x19,0, 0x19,1, 0x19,2, 0x19,3, 
                            0x19,4, 0x19,5, 0x19,6, 0x19,7,  // 80-87
                            0x09,0, 0x09,1, 0x09,2, 0x09,3, 
                            0x09,4, 0x09,5, 0x09,6, 0x09,7,  // 88-95
                            0x0D,0, 0x0D,1, 0x0D,2, 0x0D,3, 
                            0x0D,4, 0x0D,5, 0x0D,6, 0x0D,7,  // 96-103
                            0x0C,0, 0x0C,1, 0x0C,2, 0x0C,3, 
                            0x0C,4, 0x0C,5, 0x0C,6, 0x0C,7,  // 104-111
                            0x1C,0, 0x1C,1, 0x1C,2, 0x1C,3, 
                            0x1C,4, 0x1C,5, 0x1C,6, 0x1C,7,  // 112-119
                            0x1D,0, 0x1D,1, 0x1D,2, 0x1D,3, 
                            0x1D,4, 0x1D,5, 0x1D,6, 0x1D,7,  // 120-127

                            0x2D,0, 0x2D,1, 0x2D,2, 0x2D,3, 
                            0x2D,4, 0x2D,5, 0x2D,6, 0x2D,7,  // 128-135
                            0x3D,0, 0x3D,1, 0x3D,2, 0x3D,3, 
                            0x3D,4, 0x3D,5, 0x3D,6, 0x3D,7,  // 136-143
                            0x3C,0, 0x3C,1, 0x3C,2, 0x3C,3, 
                            0x3C,4, 0x3C,5, 0x3C,6, 0x3C,7,  // 144-151
                            0x2C,0, 0x2C,1, 0x2C,2, 0x2C,3, 
                            0x2C,4, 0x2C,5, 0x2C,6, 0x2C,7,  // 152-159
                            0x28,0, 0x28,1, 0x28,2, 0x28,3, 
                            0x28,4, 0x28,5, 0x28,6, 0x28,7,  // 160-167
                            0x38,0, 0x38,1, 0x38,2, 0x38,3, 
                            0x38,4, 0x38,5, 0x38,6, 0x38,7,  // 168-175
                            0x39,0, 0x39,1, 0x39,2, 0x39,3, 
                            0x39,4, 0x39,5, 0x39,6, 0x39,7,  // 176-183
                            0x29,0, 0x29,1, 0x29,2, 0x29,3, 
                            0x29,4, 0x29,5, 0x29,6, 0x29,7,  // 184-191

                            0x25,0, 0x25,1, 0x25,2, 0x25,3, 
                            0x25,4, 0x25,5, 0x25,6, 0x25,7,  // 192-199
                            0x35,0, 0x35,1, 0x35,2, 0x35,3, 
                            0x35,4, 0x35,5, 0x35,6, 0x35,7,  // 200-207
                            0x34,0, 0x34,1, 0x34,2, 0x34,3, 
                            0x34,4, 0x34,5, 0x34,6, 0x34,7,  // 208-215
                            0x24,0, 0x24,1, 0x24,2, 0x24,3, 
                            0x24,4, 0x24,5, 0x24,6, 0x24,7,  // 216-223
                            0x20,0, 0x20,1, 0x20,2, 0x20,3, 
                            0x20,4, 0x20,5, 0x20,6, 0x20,7,  // 224-231
                            0x30,0, 0x30,1, 0x30,2, 0x30,3, 
                            0x30,4, 0x30,5, 0x30,6, 0x30,7,  // 232-239
                            0x31,0, 0x31,1, 0x31,2, 0x31,3, 
                            0x31,4, 0x31,5, 0x31,6, 0x31,7,  // 240-247
                            0x21,0, 0x21,1, 0x21,2, 0x21,3, 
                            0x21,4, 0x21,5, 0x21,6, 0x21,7); // 248-255


    // ======================================================================

    case FMPI_CFG_512_MB:
      return fmpi_rank_init(512,
                            0x00,0, 0x00,1, 0x00,2, 0x00,3, 
                            0x00,4, 0x00,5, 0x00,6, 0x00,7,  // 0-7
                            0x10,0, 0x10,1, 0x10,2, 0x10,3, 
                            0x10,4, 0x10,5, 0x10,6, 0x10,7,  // 8-15
                            0x11,0, 0x11,1, 0x11,2, 0x11,3, 
                            0x11,4, 0x11,5, 0x11,6, 0x11,7,  // 16-23
                            0x01,0, 0x01,1, 0x01,2, 0x01,3, 
                            0x01,4, 0x01,5, 0x01,6, 0x01,7,  // 24-31
                            0x05,0, 0x05,1, 0x05,2, 0x05,3, 
                            0x05,4, 0x05,5, 0x05,6, 0x05,7,  // 32-39
                            0x15,0, 0x15,1, 0x15,2, 0x15,3, 
                            0x15,4, 0x15,5, 0x15,6, 0x15,7,  // 40-47
                            0x14,0, 0x14,1, 0x14,2, 0x14,3, 
                            0x14,4, 0x14,5, 0x14,6, 0x14,7,  // 48-55
                            0x04,0, 0x04,1, 0x04,2, 0x04,3, 
                            0x04,4, 0x04,5, 0x04,6, 0x04,7,  // 56-63

                            0x08,0, 0x08,1, 0x08,2, 0x08,3, 
                            0x08,4, 0x08,5, 0x08,6, 0x08,7,  // 64-71
                            0x18,0, 0x18,1, 0x18,2, 0x18,3, 
                            0x18,4, 0x18,5, 0x18,6, 0x18,7,  // 72-79
                            0x19,0, 0x19,1, 0x19,2, 0x19,3, 
                            0x19,4, 0x19,5, 0x19,6, 0x19,7,  // 80-87
                            0x09,0, 0x09,1, 0x09,2, 0x09,3, 
                            0x09,4, 0x09,5, 0x09,6, 0x09,7,  // 88-95
                            0x0D,0, 0x0D,1, 0x0D,2, 0x0D,3, 
                            0x0D,4, 0x0D,5, 0x0D,6, 0x0D,7,  // 96-103
                            0x0C,0, 0x0C,1, 0x0C,2, 0x0C,3, 
                            0x0C,4, 0x0C,5, 0x0C,6, 0x0C,7,  // 104-111
                            0x1C,0, 0x1C,1, 0x1C,2, 0x1C,3, 
                            0x1C,4, 0x1C,5, 0x1C,6, 0x1C,7,  // 112-119
                            0x1D,0, 0x1D,1, 0x1D,2, 0x1D,3, 
                            0x1D,4, 0x1D,5, 0x1D,6, 0x1D,7,  // 120-127

                            0x2D,0, 0x2D,1, 0x2D,2, 0x2D,3, 
                            0x2D,4, 0x2D,5, 0x2D,6, 0x2D,7,  // 128-135
                            0x3D,0, 0x3D,1, 0x3D,2, 0x3D,3, 
                            0x3D,4, 0x3D,5, 0x3D,6, 0x3D,7,  // 136-143
                            0x3C,0, 0x3C,1, 0x3C,2, 0x3C,3, 
                            0x3C,4, 0x3C,5, 0x3C,6, 0x3C,7,  // 144-151
                            0x2C,0, 0x2C,1, 0x2C,2, 0x2C,3, 
                            0x2C,4, 0x2C,5, 0x2C,6, 0x2C,7,  // 152-159
                            0x28,0, 0x28,1, 0x28,2, 0x28,3, 
                            0x28,4, 0x28,5, 0x28,6, 0x28,7,  // 160-167
                            0x38,0, 0x38,1, 0x38,2, 0x38,3, 
                            0x38,4, 0x38,5, 0x38,6, 0x38,7,  // 168-175
                            0x39,0, 0x39,1, 0x39,2, 0x39,3, 
                            0x39,4, 0x39,5, 0x39,6, 0x39,7,  // 176-183
                            0x29,0, 0x29,1, 0x29,2, 0x29,3, 
                            0x29,4, 0x29,5, 0x29,6, 0x29,7,  // 184-191

                            0x25,0, 0x25,1, 0x25,2, 0x25,3, 
                            0x25,4, 0x25,5, 0x25,6, 0x25,7,  // 192-199
                            0x35,0, 0x35,1, 0x35,2, 0x35,3, 
                            0x35,4, 0x35,5, 0x35,6, 0x35,7,  // 200-207
                            0x34,0, 0x34,1, 0x34,2, 0x34,3, 
                            0x34,4, 0x34,5, 0x34,6, 0x34,7,  // 208-215
                            0x24,0, 0x24,1, 0x24,2, 0x24,3, 
                            0x24,4, 0x24,5, 0x24,6, 0x24,7,  // 216-223
                            0x20,0, 0x20,1, 0x20,2, 0x20,3, 
                            0x20,4, 0x20,5, 0x20,6, 0x20,7,  // 224-231
                            0x30,0, 0x30,1, 0x30,2, 0x30,3, 
                            0x30,4, 0x30,5, 0x30,6, 0x30,7,  // 232-239
                            0x31,0, 0x31,1, 0x31,2, 0x31,3, 
                            0x31,4, 0x31,5, 0x31,6, 0x31,7,  // 240-247
                            0x21,0, 0x21,1, 0x21,2, 0x21,3, 
                            0x21,4, 0x21,5, 0x21,6, 0x21,7,  // 248-255

                            0x02,0, 0x02,1, 0x02,2, 0x02,3, 
                            0x02,4, 0x02,5, 0x02,6, 0x02,7,  // 256-263
                            0x12,0, 0x12,1, 0x12,2, 0x12,3, 
                            0x12,4, 0x12,5, 0x12,6, 0x12,7,  // 264-271
                            0x13,0, 0x13,1, 0x13,2, 0x13,3, 
                            0x13,4, 0x13,5, 0x13,6, 0x13,7,  // 272-279
                            0x03,0, 0x03,1, 0x03,2, 0x03,3, 
                            0x03,4, 0x03,5, 0x03,6, 0x03,7,  // 280-287
                            0x07,0, 0x07,1, 0x07,2, 0x07,3, 
                            0x07,4, 0x07,5, 0x07,6, 0x07,7,  // 288-295
                            0x17,0, 0x17,1, 0x17,2, 0x17,3, 
                            0x17,4, 0x17,5, 0x17,6, 0x17,7,  // 296-303
                            0x16,0, 0x16,1, 0x16,2, 0x16,3, 
                            0x16,4, 0x16,5, 0x16,6, 0x16,7,  // 304-311
                            0x06,0, 0x06,1, 0x06,2, 0x06,3, 
                            0x06,4, 0x06,5, 0x06,6, 0x06,7,  // 312-319

                            0x0A,0, 0x0A,1, 0x0A,2, 0x0A,3, 
                            0x0A,4, 0x0A,5, 0x0A,6, 0x0A,7,  // 320-327
                            0x1A,0, 0x1A,1, 0x1A,2, 0x1A,3, 
                            0x1A,4, 0x1A,5, 0x1A,6, 0x1A,7,  // 328-335
                            0x1B,0, 0x1B,1, 0x1B,2, 0x1B,3, 
                            0x1B,4, 0x1B,5, 0x1B,6, 0x1B,7,  // 336-343
                            0x0B,0, 0x0B,1, 0x0B,2, 0x0B,3, 
                            0x0B,4, 0x0B,5, 0x0B,6, 0x0B,7,  // 344-351
                            0x0F,0, 0x0F,1, 0x0F,2, 0x0F,3, 
                            0x0F,4, 0x0F,5, 0x0F,6, 0x0F,7,  // 352-359
                            0x0E,0, 0x0E,1, 0x0E,2, 0x0E,3, 
                            0x0E,4, 0x0E,5, 0x0E,6, 0x0E,7,  // 360-367
                            0x1E,0, 0x1E,1, 0x1E,2, 0x1E,3, 
                            0x1E,4, 0x1E,5, 0x1E,6, 0x1E,7,  // 368-375
                            0x1F,0, 0x1F,1, 0x1F,2, 0x1F,3, 
                            0x1F,4, 0x1F,5, 0x1F,6, 0x1F,7,  // 376-383

                            0x2F,0, 0x2F,1, 0x2F,2, 0x2F,3, 
                            0x2F,4, 0x2F,5, 0x2F,6, 0x2F,7,  // 384-391
                            0x3F,0, 0x3F,1, 0x3F,2, 0x3F,3, 
                            0x3F,4, 0x3F,5, 0x3F,6, 0x3F,7,  // 392-399
                            0x3E,0, 0x3E,1, 0x3E,2, 0x3E,3, 
                            0x3E,4, 0x3E,5, 0x3E,6, 0x3E,7,  // 400-407
                            0x2E,0, 0x2E,1, 0x2E,2, 0x2E,3, 
                            0x2E,4, 0x2E,5, 0x2E,6, 0x2E,7,  // 408-415
                            0x2A,0, 0x2A,1, 0x2A,2, 0x2A,3, 
                            0x2A,4, 0x2A,5, 0x2A,6, 0x2A,7,  // 416-423
                            0x3A,0, 0x3A,1, 0x3A,2, 0x3A,3, 
                            0x3A,4, 0x3A,5, 0x3A,6, 0x3A,7,  // 424-431
                            0x3B,0, 0x3B,1, 0x3B,2, 0x3B,3, 
                            0x3B,4, 0x3B,5, 0x3B,6, 0x3B,7,  // 432-439
                            0x2B,0, 0x2B,1, 0x2B,2, 0x2B,3, 
                            0x2B,4, 0x2B,5, 0x2B,6, 0x2B,7,  // 440-447

                            0x27,0, 0x27,1, 0x27,2, 0x27,3, 
                            0x27,4, 0x27,5, 0x27,6, 0x27,7,  // 448-455
                            0x37,0, 0x37,1, 0x37,2, 0x37,3, 
                            0x37,4, 0x37,5, 0x37,6, 0x37,7,  // 456-463
                            0x36,0, 0x36,1, 0x36,2, 0x36,3, 
                            0x36,4, 0x36,5, 0x36,6, 0x36,7,  // 464-471
                            0x26,0, 0x26,1, 0x26,2, 0x26,3, 
                            0x26,4, 0x26,5, 0x26,6, 0x26,7,  // 472-479
                            0x22,0, 0x22,1, 0x22,2, 0x22,3, 
                            0x22,4, 0x22,5, 0x22,6, 0x22,7,  // 480-487
                            0x32,0, 0x32,1, 0x32,2, 0x32,3, 
                            0x32,4, 0x32,5, 0x32,6, 0x32,7,  // 488-495
                            0x33,0, 0x33,1, 0x33,2, 0x33,3, 
                            0x33,4, 0x33,5, 0x33,6, 0x33,7,  // 496-503
                            0x23,0, 0x23,1, 0x23,2, 0x23,3, 
                            0x23,4, 0x23,5, 0x23,6, 0x23,7); // 504-511


    // ======================================================================

    case FMPI_CFG_516_HET:
      return fmpi_rank_init(516,
                            0x00,0, 0x00,1, 0x00,2, 0x00,3, 
                            0x00,4, 0x00,5, 0x00,6, 0x00,7,  // 0-7
                            0x10,0, 0x10,1, 0x10,2, 0x10,3, 
                            0x10,4, 0x10,5, 0x10,6, 0x10,7,  // 8-15
                            0x11,0, 0x11,1, 0x11,2, 0x11,3, 
                            0x11,4, 0x11,5, 0x11,6, 0x11,7,  // 16-23
                            0x01,0, 0x01,1, 0x01,2, 0x01,3, 
                            0x01,4, 0x01,5, 0x01,6, 0x01,7,  // 24-31
                            0x05,0, 0x05,1, 0x05,2, 0x05,3, 
                            0x05,4, 0x05,5, 0x05,6, 0x05,7,  // 32-39
                            0x15,0, 0x15,1, 0x15,2, 0x15,3, 
                            0x15,4, 0x15,5, 0x15,6, 0x15,7,  // 40-47
                            0x14,0, 0x14,1, 0x14,2, 0x14,3, 
                            0x14,4, 0x14,5, 0x14,6, 0x14,7,  // 48-55
                            0x04,0, 0x04,1, 0x04,2, 0x04,3, 
                            0x04,4, 0x04,5, 0x04,6, 0x04,7,  // 56-63

                            0x08,0, 0x08,1, 0x08,2, 0x08,3, 
                            0x08,4, 0x08,5, 0x08,6, 0x08,7,  // 64-71
                            0x18,0, 0x18,1, 0x18,2, 0x18,3, 
                            0x18,4, 0x18,5, 0x18,6, 0x18,7,  // 72-79
                            0x19,0, 0x19,1, 0x19,2, 0x19,3, 
                            0x19,4, 0x19,5, 0x19,6, 0x19,7,  // 80-87
                            0x09,0, 0x09,1, 0x09,2, 0x09,3, 
                            0x09,4, 0x09,5, 0x09,6, 0x09,7,  // 88-95
                            0x0D,0, 0x0D,1, 0x0D,2, 0x0D,3, 
                            0x0D,4, 0x0D,5, 0x0D,6, 0x0D,7,  // 96-103
                            0x0C,0, 0x0C,1, 0x0C,2, 0x0C,3, 
                            0x0C,4, 0x0C,5, 0x0C,6, 0x0C,7,  // 104-111
                            0x1C,0, 0x1C,1, 0x1C,2, 0x1C,3, 
                            0x1C,4, 0x1C,5, 0x1C,6, 0x1C,7,  // 112-119
                            0x1D,0, 0x1D,1, 0x1D,2, 0x1D,3, 
                            0x1D,4, 0x1D,5, 0x1D,6, 0x1D,7,  // 120-127

                            0x2D,0, 0x2D,1, 0x2D,2, 0x2D,3, 
                            0x2D,4, 0x2D,5, 0x2D,6, 0x2D,7,  // 128-135
                            0x3D,0, 0x3D,1, 0x3D,2, 0x3D,3, 
                            0x3D,4, 0x3D,5, 0x3D,6, 0x3D,7,  // 136-143
                            0x3C,0, 0x3C,1, 0x3C,2, 0x3C,3, 
                            0x3C,4, 0x3C,5, 0x3C,6, 0x3C,7,  // 144-151
                            0x2C,0, 0x2C,1, 0x2C,2, 0x2C,3, 
                            0x2C,4, 0x2C,5, 0x2C,6, 0x2C,7,  // 152-159
                            0x28,0, 0x28,1, 0x28,2, 0x28,3, 
                            0x28,4, 0x28,5, 0x28,6, 0x28,7,  // 160-167
                            0x38,0, 0x38,1, 0x38,2, 0x38,3, 
                            0x38,4, 0x38,5, 0x38,6, 0x38,7,  // 168-175
                            0x39,0, 0x39,1, 0x39,2, 0x39,3, 
                            0x39,4, 0x39,5, 0x39,6, 0x39,7,  // 176-183
                            0x29,0, 0x29,1, 0x29,2, 0x29,3, 
                            0x29,4, 0x29,5, 0x29,6, 0x29,7,  // 184-191

                            0x25,0, 0x25,1, 0x25,2, 0x25,3, 
                            0x25,4, 0x25,5, 0x25,6, 0x25,7,  // 192-199
                            0x35,0, 0x35,1, 0x35,2, 0x35,3, 
                            0x35,4, 0x35,5, 0x35,6, 0x35,7,  // 200-207
                            0x34,0, 0x34,1, 0x34,2, 0x34,3, 
                            0x34,4, 0x34,5, 0x34,6, 0x34,7,  // 208-215
                            0x24,0, 0x24,1, 0x24,2, 0x24,3, 
                            0x24,4, 0x24,5, 0x24,6, 0x24,7,  // 216-223
                            0x20,0, 0x20,1, 0x20,2, 0x20,3, 
                            0x20,4, 0x20,5, 0x20,6, 0x20,7,  // 224-231
                            0x30,0, 0x30,1, 0x30,2, 0x30,3, 
                            0x30,4, 0x30,5, 0x30,6, 0x30,7,  // 232-239
                            0x31,0, 0x31,1, 0x31,2, 0x31,3, 
                            0x31,4, 0x31,5, 0x31,6, 0x31,7,  // 240-247
                            0x21,0, 0x21,1, 0x21,2, 0x21,3, 
                            0x21,4, 0x21,5, 0x21,6, 0x21,7,  // 248-255

                            0x02,0, 0x02,1, 0x02,2, 0x02,3, 
                            0x02,4, 0x02,5, 0x02,6, 0x02,7,  // 256-263
                            0x12,0, 0x12,1, 0x12,2, 0x12,3, 
                            0x12,4, 0x12,5, 0x12,6, 0x12,7,  // 264-271
                            0x13,0, 0x13,1, 0x13,2, 0x13,3, 
                            0x13,4, 0x13,5, 0x13,6, 0x13,7,  // 272-279
                            0x03,0, 0x03,1, 0x03,2, 0x03,3, 
                            0x03,4, 0x03,5, 0x03,6, 0x03,7,  // 280-287
                            0x07,0, 0x07,1, 0x07,2, 0x07,3, 
                            0x07,4, 0x07,5, 0x07,6, 0x07,7,  // 288-295
                            0x17,0, 0x17,1, 0x17,2, 0x17,3, 
                            0x17,4, 0x17,5, 0x17,6, 0x17,7,  // 296-303
                            0x16,0, 0x16,1, 0x16,2, 0x16,3, 
                            0x16,4, 0x16,5, 0x16,6, 0x16,7,  // 304-311
                            0x06,0, 0x06,1, 0x06,2, 0x06,3, 
                            0x06,4, 0x06,5, 0x06,6, 0x06,7,  // 312-319

                            0x0A,0, 0x0A,1, 0x0A,2, 0x0A,3, 
                            0x0A,4, 0x0A,5, 0x0A,6, 0x0A,7,  // 320-327
                            0x1A,0, 0x1A,1, 0x1A,2, 0x1A,3, 
                            0x1A,4, 0x1A,5, 0x1A,6, 0x1A,7,  // 328-335
                            0x1B,0, 0x1B,1, 0x1B,2, 0x1B,3, 
                            0x1B,4, 0x1B,5, 0x1B,6, 0x1B,7,  // 336-343
                            0x0B,0, 0x0B,1, 0x0B,2, 0x0B,3, 
                            0x0B,4, 0x0B,5, 0x0B,6, 0x0B,7,  // 344-351
                            0x0F,0, 0x0F,1, 0x0F,2, 0x0F,3, 
                            0x0F,4, 0x0F,5, 0x0F,6, 0x0F,7,  // 352-359
                            0x0E,0, 0x0E,1, 0x0E,2, 0x0E,3, 
                            0x0E,4, 0x0E,5, 0x0E,6, 0x0E,7,  // 360-367
                            0x1E,0, 0x1E,1, 0x1E,2, 0x1E,3, 
                            0x1E,4, 0x1E,5, 0x1E,6, 0x1E,7,  // 368-375
                            0x1F,0, 0x1F,1, 0x1F,2, 0x1F,3, 
                            0x1F,4, 0x1F,5, 0x1F,6, 0x1F,7,  // 376-383

                            0x2F,0, 0x2F,1, 0x2F,2, 0x2F,3, 
                            0x2F,4, 0x2F,5, 0x2F,6, 0x2F,7,  // 384-391
                            0x3F,0, 0x3F,1, 0x3F,2, 0x3F,3, 
                            0x3F,4, 0x3F,5, 0x3F,6, 0x3F,7,  // 392-399
                            0x3E,0, 0x3E,1, 0x3E,2, 0x3E,3, 
                            0x3E,4, 0x3E,5, 0x3E,6, 0x3E,7,  // 400-407
                            0x2E,0, 0x2E,1, 0x2E,2, 0x2E,3, 
                            0x2E,4, 0x2E,5, 0x2E,6, 0x2E,7,  // 408-415
                            0x2A,0, 0x2A,1, 0x2A,2, 0x2A,3, 
                            0x2A,4, 0x2A,5, 0x2A,6, 0x2A,7,  // 416-423
                            0x3A,0, 0x3A,1, 0x3A,2, 0x3A,3, 
                            0x3A,4, 0x3A,5, 0x3A,6, 0x3A,7,  // 424-431
                            0x3B,0, 0x3B,1, 0x3B,2, 0x3B,3, 
                            0x3B,4, 0x3B,5, 0x3B,6, 0x3B,7,  // 432-439
                            0x2B,0, 0x2B,1, 0x2B,2, 0x2B,3, 
                            0x2B,4, 0x2B,5, 0x2B,6, 0x2B,7,  // 440-447

                            0x27,0, 0x27,1, 0x27,2, 0x27,3, 
                            0x27,4, 0x27,5, 0x27,6, 0x27,7,  // 448-455
                            0x37,0, 0x37,1, 0x37,2, 0x37,3, 
                            0x37,4, 0x37,5, 0x37,6, 0x37,7,  // 456-463
                            0x36,0, 0x36,1, 0x36,2, 0x36,3, 
                            0x36,4, 0x36,5, 0x36,6, 0x36,7,  // 464-471
                            0x26,0, 0x26,1, 0x26,2, 0x26,3, 
                            0x26,4, 0x26,5, 0x26,6, 0x26,7,  // 472-479
                            0x22,0, 0x22,1, 0x22,2, 0x22,3, 
                            0x22,4, 0x22,5, 0x22,6, 0x22,7,  // 480-487
                            0x32,0, 0x32,1, 0x32,2, 0x32,3, 
                            0x32,4, 0x32,5, 0x32,6, 0x32,7,  // 488-495
                            0x33,0, 0x33,1, 0x33,2, 0x33,3, 
                            0x33,4, 0x33,5, 0x33,6, 0x33,7,  // 496-503
                            0x23,0, 0x23,1, 0x23,2, 0x23,3, 
                            0x23,4, 0x23,5, 0x23,6, 0x23,7,  // 504-511

                            0x6B,0, 0x6B,1, 0x6B,2, 0x6B,3); // 512-515



    // Undeclared mode
    default: 
      ar_abort();
  }
}

