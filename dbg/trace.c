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
// Abstract      : Trace-collecting functions and output for Paraver
//                 visualization
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: trace.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/12/12 16:37:19 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <debug.h>


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
void dbg_trace_init() {
#ifdef DBG_TRC_ENABLED

  Context       *context;
  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  unsigned int  i;
  unsigned int  my;
  unsigned int  ref;


  // Get context
  context = mm_get_context(ar_get_core_id());
  my_cid = ar_get_core_id();
  my_bid = ar_get_board_id();

  // Init index
  context->dbg_trc_idx = 0;

  // Allocate array
  ar_assert(!context->dbg_trc_data);
  context->dbg_trc_data = kt_malloc(2 * DBG_TRC_MAX_ENTRIES * sizeof(int));

  // Top-level scheduler is the master core for timer offset handshaking
  if (context->pr_parent_sched_id == -1) {

    // For all cores in the setup
    for (i = 0; i < context->pr_num_cores; i++) {

      // Get arch-level board/core ID
      pr_core_arch_bid_cid(i, &bid, &cid);

      // If it's us, we have no offset
      if ((bid == my_bid) && (cid == my_cid)) {
        context->dbg_trc_offset = 0;
        context->dbg_trc_time_start = ar_free_timer_get_ticks();
        continue;
      }

      // Handshake with peer and send to him our timer, so he can record
      // its offset

      // Step 1: send our board/core ID, so he can communicate back
      ar_mbox_send(my_cid, bid, cid, DBG_TRC_MAGIC_HSHAKE1);
      ar_mbox_send(my_cid, bid, cid, (my_bid << 8) | my_cid);
      
      // Step 2: get back acknowledge
      ar_assert(ar_mbox_get(my_cid) == DBG_TRC_MAGIC_HSHAKE2);

      // Step 3: send our timer
      ar_mbox_send(my_cid, bid, cid, ar_free_timer_get_ticks());
    }

    // We're done, inform all others to exit their barriers
    for (i = 0; i < context->pr_num_cores; i++) {
      pr_core_arch_bid_cid(i, &bid, &cid);
      if ((bid != my_bid) || (cid != my_cid)) {
        ar_cnt_incr(my_cid, bid, cid, NOC_COUNTER_WAKEUP2, 1);
      }
    }
  }
  // Other cores are slaves
  else {

    // Initialize our barrier counter to -1
    ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP2, -1);

    // Step 1: receive master board/core ID
    ar_assert(ar_mbox_get(my_cid) == DBG_TRC_MAGIC_HSHAKE1);
    i = ar_mbox_get(my_cid);
    bid = i >> 8;
    cid = i & 0xFF;

    // Step 2: acknowledge
    ar_mbox_send(my_cid, bid, cid, DBG_TRC_MAGIC_HSHAKE2);

    // Step 3: get remote timer and compute offset
    my = ar_free_timer_get_ticks();
    ref = ar_mbox_get(my_cid);
    context->dbg_trc_time_start = my;
    context->dbg_trc_offset = my - ref;

    // Block until master says we're done
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP2)) {
      ;
    }
  }
#endif
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
void dbg_trace_print_all_records(Context *context, unsigned int time_end, 
                                 int offset_end) {
#ifdef DBG_TRC_ENABLED

  int           i;
  unsigned int  elapsed;
  int           drift;
  float         point;
  unsigned int  adjusted;


  // Copute how much time has passed locally, and what's the drift in clock
  // cycles w.r.t. the reference remote clock
  elapsed = time_end - context->dbg_trc_time_start;
  drift = offset_end - context->dbg_trc_offset;

  //kt_printf("core %3d: time_start = %d, time_end = %d, elapsed = %d\r\n"
  //          "          offset_begin = %d, offset_end = %d, drift = %d\r\n",
  //          context->pr_core_id, 
  //          context->dbg_trc_time_start, time_end, elapsed,
  //          context->dbg_trc_offset, offset_end, drift);

  // Adjust timestamps
  for (i = 0; i < context->dbg_trc_idx; i += 2) {

    // Compute where's this trace w.r.t. the local elapsed time, and assign
    // accordingly a portion of the drift to it. Also adjust the time so it
    // starts together with the start of the reference clock.
    point = context->dbg_trc_data[i] - context->dbg_trc_time_start;
    point *= drift;
    point /= elapsed;
    adjusted = context->dbg_trc_data[i] - context->dbg_trc_offset - point;

    context->dbg_trc_data[i] = adjusted;

    // Print trace with adjusted time
    //kt_printf("%X %X %X\r\n", context->pr_core_id, adjusted, 
    //                          context->dbg_trc_data[i + 1]);

    //kt_printf("core %d, adj_time %10d, org_time %10d, evt %X\r\n", 
    //          context->pr_core_id, 
    //          adjusted,
    //          context->dbg_trc_data[i],
    //          context->dbg_trc_data[i + 1]);
  }

  // Print preamble
  kt_printf("=== Core %d begin ===\r\n", context->pr_core_id);

  // Encode buffer and print it
  kt_encode85(context->dbg_trc_data, context->dbg_trc_idx);
  
  // Print epilogue
  kt_printf("=== Core %d end ===\r\n\n", context->pr_core_id);

#endif
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
void dbg_trace_report(char *filename) {
#ifdef DBG_TRC_ENABLED

  Context       *context;
  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  unsigned int  i;
  unsigned int  my;
  unsigned int  ref;


  // Get context
  context = mm_get_context(ar_get_core_id());
  my_cid = ar_get_core_id();
  my_bid = ar_get_board_id();


  // Top-level scheduler is the master core for reporting
  if (context->pr_parent_sched_id == -1) {

    // Print file begin dump header
    ar_uart_flush();
    ar_timer_busy_wait_msec(200);
    kt_printf(DBG_FILE_DUMP_BEGIN_FORMAT, filename);
    ar_uart_flush();
    ar_timer_busy_wait_msec(20);

    // For all cores in the setup
    for (i = 0; i < context->pr_num_cores; i++) {

      // Get arch-level board/core ID
      pr_core_arch_bid_cid(i, &bid, &cid);

      // Is it us? Report with no offset.
      if ((bid == my_bid) && (cid == my_cid)) {
        dbg_trace_print_all_records(context, ar_free_timer_get_ticks(), 0);
        continue;
      }

      // Handshake with peer and send to him our timer, so he can re-estimate
      // its offset with our current clock, as well as where time 0 is set
      // according to our clock.

      // Step 1: send our board/core ID, so he can communicate back
      ar_mbox_send(my_cid, bid, cid, DBG_TRC_MAGIC_HSHAKE3);
      ar_mbox_send(my_cid, bid, cid, (my_bid << 8) | my_cid);
      
      // Step 2: get back acknowledge
      ar_assert(ar_mbox_get(my_cid) == DBG_TRC_MAGIC_HSHAKE4);

      // Step 3: send our timer
      ar_mbox_send(my_cid, bid, cid, ar_free_timer_get_ticks());
      
      // Step 4: get back acknowledge
      ar_assert(ar_mbox_get(my_cid) == DBG_TRC_MAGIC_HSHAKE5);
    }

    // Print file end dump header
    ar_uart_flush();
    ar_timer_busy_wait_msec(200);
    kt_printf(DBG_FILE_DUMP_END_FORMAT, filename);
    ar_uart_flush();
    ar_timer_busy_wait_msec(20);

  }
  // Other cores are slaves
  else {

    // Step 1: receive master board/core ID
    ar_assert(ar_mbox_get(my_cid) == DBG_TRC_MAGIC_HSHAKE3);
    i = ar_mbox_get(my_cid);
    bid = i >> 8;
    cid = i & 0xFF;

    // Step 2: acknowledge
    ar_mbox_send(my_cid, bid, cid, DBG_TRC_MAGIC_HSHAKE4);

    // Step 3: get remote timer
    my = ar_free_timer_get_ticks();
    ref = ar_mbox_get(my_cid);

    // Report with computed offset
    dbg_trace_print_all_records(context, my, my - ref);

    // Flush the UART, so others can follow without scrambling
    ar_uart_flush();
    ar_timer_busy_wait_msec(20);

    // Step 4: say to master we're done
    ar_mbox_send(my_cid, bid, cid, DBG_TRC_MAGIC_HSHAKE5);
  }

#endif
}
