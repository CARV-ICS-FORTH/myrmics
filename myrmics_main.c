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
// Abstract      : Main C language entry point for the Myrmics version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: myrmics_main.c,v $
// CVS revision  : $Revision: 1.67 $
// Last modified : $Date: 2013/03/27 09:55:03 $
// Last author   : $Author: jacob $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <noc.h>
#include <memory_management.h>
#include <processing.h>
#include <video.h>
#include <bench.h>
#include <debug.h>



// ==========================================================================
// main_print_infomercial: Give info about our core setup
// ==========================================================================
void main_print_infomercial(Context *context, int my_bid, int my_cid) {
  if (context->pr_role == ROLE_SCHEDULER) {
    if (context->pr_num_cores > 1) {
      kt_printf(
        "%3d: L%d Scheduler, sch_id = %3d, parent_id = %3d [0x%02X/%d]\r\n",
        context->pr_core_id, context->pr_scheduler_level,
        context->pr_scheduler_id, context->pr_parent_sched_id, my_bid, my_cid);
    }
    else {
      kt_printf(
        "0: Scheduler / Worker (single-core mode) [0x%02X/%d]\r\n", 
        my_bid, my_cid);
    }
  }
  else {
    kt_printf(
      "%3d: Worker,    worker_id = %3d, parent_id = %3d [0x%02X/%d]\r\n",
      context->pr_core_id, context->pr_worker_id, context->pr_parent_sched_id, 
      my_bid, my_cid);
  }
}


// ===========================================================================
// main_run_app()               Initializes a core setup and runs the 
//                              chosen application
// ===========================================================================
// * INPUTS
//   int my_bid                 Arch-level board ID for this core
//   int my_cid                 Arch-level core ID for this core
//   PrCfgMode core_setup       Chosen core setup from include/processing.h
//   int auto_levels            Tree levels for PR_CFG_AUTO_MB setup
//   int *auto_fanouts          Tree fanouts for PR_CFG_AUTO_MB
//   NocMode noc_setup          Chosen network mode from include/noc.h
//   int noc_credits            Credits per peer for NOC_MODE_CREDIT_ONLY
//   int load_vs_locality       If 0-100, changes the default value of
//                              context->pr_load_vs_locality. If -1, leaves
//                              it unchanged.
//   char *stats_filename       If not NULL, host filename to upload the
//                              debug statistics, after app execution
//   char *paraver_filename     If not NULL, host filename to upload the
//                              debug traces for Paraver, after app execution
//   func_t **task_table        Function task pointer table for application.
//                              The main task of the application should be
//                              in task_table[0].
//   ...                        NULL-terminated list of arguments to give
//                              to the function in task_table[0]
// ==========================================================================
void main_run_app(int my_bid, int my_cid, PrCfgMode core_setup, 
                  int auto_levels, int *auto_fanouts, NocMode noc_setup, 
                  int noc_credits, int load_vs_locality, char *stats_filename, 
                  char *paraver_filename, func_t **task_table, ...) {

  Context       *context;
  int           active;
  int           ret;
  int           bid;
  int           cid;
  int           x;
  int           y;
  int           z;
  va_list       ap;
  void          **args;
  void          *a;
  int           barrier;
  int           video_cores_active;
  int           i;


  // Prepare the counter used to handle reentrant main_run_app() calls
  ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP3, -2);

  // Initialize kernel memory allocator
  kt_mem_init();

  // We now have an initialized context
  context = mm_get_context(my_cid);

  // Override pr_load_vs_locality?
  if (load_vs_locality != -1) {
    context->pr_load_vs_locality = load_vs_locality;
  }


  // Configure processing layer with selected core setup
  ret = pr_configure(core_setup, auto_levels, auto_fanouts);

  // Are we part of the Myrmics core setup?
  if (ret == 0) {
    active = 1;
  }
  // If not, skip to the end
  else if (ret == 1) {
    active = 0;
    goto skip;
  }
  // If we're one of the special cores to help the video demo, run the
  // appropriate functions
  else if (ret == 2) {
    demo_myrmics_input_loop(my_bid, my_cid);
  }
  else if (ret == 3) {
    demo_myrmics_output_loop(my_bid, my_cid);
  }
  // Unexpected pr_configure() return value
  else {
    ar_abort();
  }


  // Infomercial: one by one, so the UART doesn't get scrambled and we get
  // a proper listing according to core ID order. We skip this printing when
  // we are in the video demo modes.
  if (context->vid_demo_in_bid == -1) {
    if (context->pr_core_id == 0) {
      for (i = 0; i < context->pr_num_cores; i++) {
        pr_core_arch_bid_cid(i, &bid, &cid);
        if ((bid == my_bid) && (cid == my_cid)) {
          main_print_infomercial(context, my_bid, my_cid);
          continue;
        }
        ar_mbox_send(my_cid, bid, cid, 0x10101010);
        ar_assert(ar_mbox_get(my_cid) == 0x01010101);
      }
      ar_uart_flush();
      ar_timer_busy_wait_msec(20);
      kt_printf("\n");
      ar_uart_flush();
      ar_timer_busy_wait_msec(20);
      for (i = 0; i < context->pr_num_cores; i++) {
        pr_core_arch_bid_cid(i, &bid, &cid);
        if ((bid == my_bid) && (cid == my_cid)) {
          continue;
        }
        ar_mbox_send(my_cid, bid, cid, 0x11111111);
      }
    }
    else {
      ar_assert(ar_mbox_get(my_cid) == 0x10101010);
      main_print_infomercial(context, my_bid, my_cid);
      ar_uart_flush();
      pr_core_arch_bid_cid(0, &bid, &cid);
      ar_mbox_send(my_cid, bid, cid, 0x01010101);
      ar_assert(ar_mbox_get(my_cid) == 0x11111111);
    }
  }


  // Initialize network-on-chip layer
  noc_init(noc_setup, noc_credits);


  // Prepare the final counter-based barrier which will be done after all
  // the cores are done from their main loops
  ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP1, (!context->pr_core_id) ? 
                                              -context->pr_num_cores + 1 : -1);

  // Initialize statistics gathering
  dbg_stats_init();
  
  // Initialize tracing interface
  dbg_trace_init();


  // Initialize user memory allocator
  mm_user_init(context->pr_role);


  // Prepare application argument array. Memory allocated here passes through
  // and is being used and deallocated internally, so don't free it.
  va_start(ap, task_table);
  args = NULL;
  i = 0;
  while (1) {
    a = va_arg(ap, void *);
    if (!a) {
      break;
    }
    args = kt_realloc(args, (i + 1) * sizeof (void *));
    args[i] = a;
    i++;
    ar_assert(i < PR_TASK_MAX_ARGS); // caller forgot the NULL termination?
  }
  va_end(ap);


  // Init application
  pr_init_app(task_table, args, i);


  // Enter event-based loop
  pr_event_main_loop();


  // Do a barrier, so we're sure everyone exited cleanly and there are no
  // pending NoC messages in mailboxes, before we start stats/trace reporting
  ar_uart_flush();
  if (!context->pr_core_id) {
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP1)) {
      ;
    }
    for (i = 1; i < context->pr_num_cores; i++) {
      pr_core_arch_bid_cid(i, &bid, &cid);
      ar_cnt_incr(my_cid, bid, cid, NOC_COUNTER_WAKEUP1, 1);
    }
  }
  else {
    pr_core_arch_bid_cid(0, &bid, &cid);
    ar_cnt_incr(my_cid, bid, cid, NOC_COUNTER_WAKEUP1, 1);
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP1)) {
      ;
    }
  }

  // In order to be prepared for reentry, we must flush all caches. Otherwise,
  // we risk leaving dirty lines that were used by the previous application to
  // some cores, which may be evicted later on and writeback stale data on the
  // new application. While we're on this, reset the cache epoch to 0, so that
  // runs are as identical as possible.
#ifdef ARCH_MB
  ar_cache_flush();
  ar_cache_set_epoch(0);
#endif


  // Report statistics
  if (stats_filename) {
    dbg_stats_report(stats_filename);
  }

  // Report traces
  if (paraver_filename) {
    dbg_trace_report(paraver_filename);
  }


skip:

  // Boot master notifies everybody else, collects answers and notifies
  // them again that we're all done.
  if ((my_bid == AR_BOOT_MASTER_BID) && (!my_cid)) {

    // Reset our barrier counter with the total number of expected answers
    video_cores_active = (context->vid_demo_in_bid > -1) +
                         (context->vid_demo_out_bid > -1);
    ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP3, 
               -(AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1) *
               (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
               (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1) * 
                                             AR_FORMIC_CORES_PER_BOARD -
               (AR_ARM0_BID != -1) * AR_ARM0_CORES_PER_BOARD -
               (AR_ARM1_BID != -1) * AR_ARM1_CORES_PER_BOARD + 1 +
               video_cores_active);

    // For each barrier phase
    for (barrier = 0; barrier < 2; barrier++) {

      // Make Formic slaves reach their barrier
      for (x = AR_FORMIC_MIN_X; x <= AR_FORMIC_MAX_X; x++) {
        for (y = AR_FORMIC_MIN_Y; y <= AR_FORMIC_MAX_Y; y++) {
          for (z = AR_FORMIC_MIN_Z; z <= AR_FORMIC_MAX_Z; z++) {
            for (i = 0; i < AR_FORMIC_CORES_PER_BOARD; i++) {
              bid = (x << 4) | (y << 2) | z;
              if ((bid != my_bid) || (i != my_cid)) {
                ar_cnt_incr(my_cid, bid, i, NOC_COUNTER_WAKEUP3, 1);
              }
            }
          }
        }
      }

      // Make ARM slaves reach their barrier
      if (AR_ARM0_BID != -1) {
        for (i = 0; i < AR_ARM0_CORES_PER_BOARD; i++) {
          if (((AR_ARM0_BID != my_bid) || (i != my_cid)) &&
              ((AR_ARM0_BID != context->vid_demo_in_bid) ||
               (i != context->vid_demo_in_cid)) &&
              ((AR_ARM0_BID != context->vid_demo_out_bid) ||
               (i != context->vid_demo_out_cid))) {
            ar_cnt_incr(my_cid, AR_ARM0_BID, i, NOC_COUNTER_WAKEUP3, 1);
          }
        }
      }
      if (AR_ARM1_BID != -1) {
        for (i = 0; i < AR_ARM1_CORES_PER_BOARD; i++) {
          if (((AR_ARM1_BID != my_bid) || (i != my_cid)) &&
              ((AR_ARM1_BID != context->vid_demo_in_bid) ||
               (i != context->vid_demo_in_cid)) &&
              ((AR_ARM1_BID != context->vid_demo_out_bid) ||
               (i != context->vid_demo_out_cid))) {
            ar_cnt_incr(my_cid, AR_ARM1_BID, i, NOC_COUNTER_WAKEUP3, 1);
          }
        }
      }

      // After we've notified them once, wait until they all answer
      if (barrier == 0) {

        while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP3)) {
          ;
        }

        // Print a separator
        ar_uart_flush();
        ar_timer_busy_wait_msec(200);
        kt_printf("\r\n***** Application run done. *****\r\n\n");
        ar_uart_flush();
        ar_timer_busy_wait_msec(100);
      }
    }
  }

  // Slaves: wait two messages from boot master to continue
  else {
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP3) == -2) {
      ;
    }
    ar_cnt_incr(my_cid, AR_BOOT_MASTER_BID, 0, NOC_COUNTER_WAKEUP3, 1);
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP3) == -1) {
      ;
    }
  }
}


// ==========================================================================
// Main: entry point after boot.s
// ==========================================================================
int main() {

  int           b;
  int           c;
  unsigned int  word;
  int           total_cores;
  int           x;
  int           y;
  int           z;
  int           slave_bid;
  int           i;


  // Get core & board ID
  b = ar_get_board_id();
  c = ar_get_core_id();

  // Arch-specific initialization (FPU/caches/MMU, etc)
  ar_init(b, c);


  // Boot master
  if ((b == AR_BOOT_MASTER_BID) && (!c)) {

    // Greet
    kt_printf("Boot master is 0x%02X/%d\r\n", b, c);

    // Check links, board versions, transfer code and wake up Formic boards
    ar_wake_up_formic_boards(b, c);

    // Compute how many cores should be activated
    total_cores = (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1) *
                  (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
                  (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1) * 
                                                AR_FORMIC_CORES_PER_BOARD +
                  (AR_ARM0_BID != -1) * AR_ARM0_CORES_PER_BOARD +
                  (AR_ARM1_BID != -1) * AR_ARM1_CORES_PER_BOARD;

    // Collect from all cores their "boot ok" message
    for (i = 0; i < total_cores - 1; i++) {
      word = ar_mbox_get(c);
      if ((word >> 16) != 0x1111) {
        kt_printf("Boot master: invalid response 0x%08X\r\n", word);
        ar_abort();
      }
    }
    kt_printf("Collected %d boot answers\r\n\n", i);

    // Make Formic slaves reach their boot barrier
    for (x = AR_FORMIC_MIN_X; x <= AR_FORMIC_MAX_X; x++) {
      for (y = AR_FORMIC_MIN_Y; y <= AR_FORMIC_MAX_Y; y++) {
        for (z = AR_FORMIC_MIN_Z; z <= AR_FORMIC_MAX_Z; z++) {
          for (i = 0; i < AR_FORMIC_CORES_PER_BOARD; i++) {

            slave_bid = (x << 4) | (y << 2) | z;

            if ((slave_bid != b) || (i != c)) {
              ar_cnt_incr(c, slave_bid, i, NOC_COUNTER_WAKEUP1, 1);
            }
          }
        }
      }
    }

    // Make ARM slaves reach their boot barrier
    if (AR_ARM0_BID != -1) {
      for (i = 0; i < AR_ARM0_CORES_PER_BOARD; i++) {
        if ((AR_ARM0_BID != b) || (i != c)) {
          ar_cnt_incr(c, AR_ARM0_BID, i, NOC_COUNTER_WAKEUP1, 1);
        }
      }
    }
    if (AR_ARM1_BID != -1) {
      for (i = 0; i < AR_ARM1_CORES_PER_BOARD; i++) {
        if ((AR_ARM1_BID != b) || (i != c)) {
          ar_cnt_incr(c, AR_ARM1_BID, i, NOC_COUNTER_WAKEUP1, 1);
        }
      }
    }
  }

  // Boot slave
  else {

    // Initialize a boot barrier counter
    ar_cnt_set(c, NOC_COUNTER_WAKEUP1, -1);

    // Signal master we've booted
    ar_mbox_send(c, AR_BOOT_MASTER_BID, 0,
                 0x11110000 | (b << 8) | c);

    // Block on barrier, until master indicates all cores have booted
    while (ar_cnt_get(c, NOC_COUNTER_WAKEUP1)) {
      ;
    }
  }


  // Core setup options:
  //
  // PR_CFG_1_MB              //   single-core setup
  // PR_CFG_1_ARM             //   single-core setup
  // PR_CFG_2_MB_FLAT         //   1 worker
  // PR_CFG_2_ARM_FLAT        //   1 worker
  // PR_CFG_2_HET_FLAT        //   1 worker
  // PR_CFG_3_MB_FLAT         //   2 workers
  // PR_CFG_3_ARM_FLAT        //   2 workers
  // PR_CFG_3_HET_FLAT        //   2 workers
  // PR_CFG_5_MB_FLAT         //   4 workers
  // PR_CFG_5_HET_FLAT        //   4 workers
  // PR_CFG_8_ARM_FLAT        //   7 workers
  // PR_CFG_9_MB_FLAT         //   8 workers
  // PR_CFG_9_HET_FLAT        //   8 workers
  // PR_CFG_11_MB_HIER        //   8 workers (4 per L0 sched)
  // PR_CFG_11_HET_HIER       //   8 workers (4 per L0 sched)
  // PR_CFG_17_MB_FLAT        //  16 workers
  // PR_CFG_17_HET_FLAT       //  16 workers
  // PR_CFG_25_MB_FLAT        //  24 workers
  // PR_CFG_25_HET_FLAT       //  24 workers
  // PR_CFG_28_HET_HIER       //  24 workers (8 per L0 sched)
  // PR_CFG_31_MB_HIER        //  16 workers (2 per L0 sched)
  // PR_CFG_33_MB_FLAT        //  32 workers
  // PR_CFG_33_HET_FLAT       //  32 workers
  // PR_CFG_35_HET_HIER       //  32 workers (16 per L0 sched)
  // PR_CFG_52_HET_HIER       //  48 workers (16 per L0 sched)
  // PR_CFG_65_HET_FLAT       //  64 workers
  // PR_CFG_69_HET_HIER       //  64 workers (16 per L0 sched)
  // PR_CFG_120_HET_HIER      // 112 workers (16 per L0 sched)
  // PR_CFG_129_HET_FLAT      // 128 workers
  // PR_CFG_136_HET_HIER      // 128 workers (18-19 per L0 sched)
  // PR_CFG_232_HET_HIER      // 224 workers (32 per L0 sched)
  // PR_CFG_257_HET_FLAT      // 256 workers
  // PR_CFG_264_HET_HIER      // 256 workers (36-37 per L0 sched)
  // PR_CFG_513_HET_FLAT      // 512 workers
  // PR_CFG_520_HET_HIER3     // 448 workers (7 per L0 sched)
  // PR_CFG_520_HET_HIER2     // 512 workers (73-74 per L0 sched)
  // PR_CFG_AUTO_MB           // auto-generated setup

  // Network setup options:
  //
  // NOC_MODE_MAILBOX_ONLY    // Mailbox-based (slow, any fanout possible)
  // NOC_MODE_CREDIT_ONLY     // Credit-based (both modes fast, up to 
  // NOC_MODE_CREDIT_MAILBOX  // ~50 fanout feasible)

  // Stats and Paraver filenames: specify filenames or NULL for no reporting

  // Finally, select application to run, by specifying its task table. The
  // SCOOP compiler creates a task table by extending the main function name
  // with _task_table. Specify arguments with a NULL-terminated list.
  //
  // WARNING: this means you can't pass a 0 value in the args :)



  // ##########################################################################
  // ###                                                                    ###
  // ###                             TEST RUNS                              ###
  // ###                                                                    ###
  // ##########################################################################

#if 0
  //main_run_app(b, c, 
  //             PR_CFG_11_MB_HIER, 0, 0, NOC_MODE_CREDIT_ONLY, 200, -1,
  //             "stats1.txt", "paraver1.evt",
  //             &test_myrmics_task_table, 
  //             NULL);

  main_run_app(b, c, 
               PR_CFG_11_MB_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats2.txt", "paraver2.evt",
               &test_myrmics_task_table, 
               NULL);

  //main_run_app(b, c, 
  //             PR_CFG_11_MB_HIER, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
  //             "stats3.txt", "paraver3.evt",
  //             &test_myrmics_task_table, 
  //             NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_9_MB_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &matrix_mult_myrmics_task_table, 
               2,       // region size (num of region rows and region cols)
               2,       // tile size (num of tile rows and tile cols per region)
               16,      // block size (num of rows and cols per tile)
               1,       // verify (1 = verify, -1 = don't, 2 = checksum only)
               -1,      // dummy regions (-1 = none)
               NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_25_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &jacobi_myrmics_task_table, 
               3,       // number of regions
               12,      // number of total tiles
               360,     // number of total rows
               1024,    // number of total columns
               2,       // number of iterations
               1,       // verify (1 = verify, -1 = don't)
               NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_25_MB_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &bitonic_myrmics_task_table, 
               4,       // number of regions
               32,      // number of total tiles
               262144,  // number of total elements
               1,       // verify (1 = verify, -1 = don't)
               -1,      // dummy regions (-1 = none)
               NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_17_MB_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &fft_myrmics_task_table, 
               4,       // number of tile rows and columns
               64,      // number of total rows and columns
               NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_11_MB_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &kmeans_myrmics_task_table, 
               16,      // number of output clusters
               1,       // number of regions to use
               3,       // tiles per region
               1024,    // input objects per tile
               4,       // loop repetitions
               NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_11_MB_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &cray_myrmics_task_table, 
               106,     // x resolution (total picture columns)
               80,      // y resolution (total picture rows)
               2,       // regions to split y resolution into
               4,       // tiles per region for further split
               NULL);
#endif

#if 0
  main_run_app(b, c, 
               PR_CFG_11_MB_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats.txt", "paraver.evt",
               &md5_myrmics_task_table, 
               64,      // number of total buffers to do in parallel
               65536,   // buffer size
               2,       // number of regions to use
               NULL);
#endif


  // ##########################################################################
  // ###                                                                    ###
  // ###                            VIDEO DEMOS                             ###
  // ###                                                                    ###
  // ##########################################################################
#if 0

  // Video passthrough
  if ((b == AR_BOOT_MASTER_BID) && (!c)) {
    kt_mem_init();

    Context *context = mm_get_context(c);
    context->pr_core_id = 0;

    demo_myrmics_passthrough();
  }
#endif
#if 0

  Context   *context;
  int       core_config;
  PrCfgMode cfg;
  NocMode   noc;

  // Raytracing demo. We're starting the video in/out cores without running
  // main_run_app(), so that we can use the video in/out cores to select
  // a core setup. We depend on video out core to send a mailslot message
  // to all other active cores, to find out the configuration. For this
  // to work, make sure the video in/out board/core IDs here are the same
  // with the PR_VID_* setups.
  if ((b == AR_ARM0_BID) && ((c == 1) || (c == 2))) {
    kt_mem_init();
    context = mm_get_context(c);
    context->vid_demo_in_bid = AR_ARM0_BID;
    context->vid_demo_in_cid = 1;
    context->vid_demo_out_bid = AR_ARM0_BID;
    context->vid_demo_out_cid = 2;
    if (c == 1) {
      demo_myrmics_input_loop(b, c);
    }
    else {
      demo_myrmics_output_loop(b, c);
    }
  }


  while (1) {
    
    // Find out core configuration for this raytraced image
    core_config = ar_mslot_get(c);
    switch (core_config) {
      
      // Flat scheduling
      case  0: cfg = PR_VID_2_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  1: cfg = PR_VID_3_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  2: cfg = PR_VID_5_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  3: cfg = PR_VID_9_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  4: cfg = PR_VID_17_HET_FLAT;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  5: cfg = PR_VID_33_HET_FLAT;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  6: cfg = PR_VID_65_HET_FLAT;  noc = NOC_MODE_MAILBOX_ONLY;   break;
      case  7: cfg = PR_VID_129_HET_FLAT; noc = NOC_MODE_MAILBOX_ONLY;   break;
      case  8: cfg = PR_VID_257_HET_FLAT; noc = NOC_MODE_MAILBOX_ONLY;   break;
      case  9: cfg = PR_VID_513_HET_FLAT; noc = NOC_MODE_MAILBOX_ONLY;   break;

      // Hier scheduling (setups <16 workers remain flat)
      case 16: cfg = PR_VID_2_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 17: cfg = PR_VID_3_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 18: cfg = PR_VID_5_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 19: cfg = PR_VID_9_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 20: cfg = PR_VID_17_HET_FLAT;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 21: cfg = PR_VID_35_HET_HIER;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 22: cfg = PR_VID_69_HET_HIER;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 23: cfg = PR_VID_133_HET_HIER; noc = NOC_MODE_CREDIT_MAILBOX; break;
      case 24: cfg = PR_VID_261_HET_HIER; noc = NOC_MODE_MAILBOX_ONLY;   break;
      case 25: cfg = PR_VID_517_HET_HIER; noc = NOC_MODE_MAILBOX_ONLY;   break;

      default: ar_abort();
    }

    // Run
    main_run_app(b, c, cfg, 0, 0, noc, 0, -1, NULL, NULL, 
                 &demo_myrmics_task_table, 
                 (cfg == PR_VID_35_HET_HIER) ? 2 :
                 ((cfg == PR_VID_69_HET_HIER) ||
                  (cfg == PR_VID_133_HET_HIER) ||
                  (cfg == PR_VID_261_HET_HIER) ||
                  (cfg == PR_VID_517_HET_HIER)) ? 4 : 1,
                 NULL);
  }
#endif



  // ##########################################################################
  // ###                                                                    ###
  // ###                             BATCH RUNS                             ###
  // ###                                                                    ###
  // ##########################################################################


  // =========================================================================
  // Various microbenchmarks
  // =========================================================================
#if 0
  // ---------------------------------------------------------
  // For the "task spawning" part in test_myrmics.c
  // ---------------------------------------------------------
  main_run_app(b, c, PR_CFG_2_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats2h.txt", "paraver2h.evt", &test_myrmics_task_table, 
               1000, NULL); // spawn: 16,164,824 cycles => 16.1 Kcycles/task
                            // exec:  13,252,721 cycles => 13.2 Kcycles/task

  main_run_app(b, c, PR_CFG_2_MB_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats2m.txt", "paraver2m.evt", &test_myrmics_task_table, 
               1000, NULL); // spawn: 37,406,751 cycles => 37.4 Kcycles/task
                            // exec:  22,529,202 cycles => 22.5 Kcycles/task

  main_run_app(b, c, PR_CFG_2_ARM_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats2a.txt", "paraver2a.evt", &test_myrmics_task_table, 
               1000, NULL); // spawn:  5,394,728 cycles =>  5.4 Kcycles/task
                            // exec:   2,713,152 cycles =>  2.7 Kcycles/task
#endif
#if 0
  // ---------------------------------------------------------
  // For the "delayed task spawning" part in test_myrmics.c.
  //
  // Use the parse_synthetic*.awk scripts to create results.
  // Adjust obvious timer overflows (e.g. for 1 worker, 8M
  // task size) by adding manually 2^32 to the result, before
  // running the script.
  // ---------------------------------------------------------
  int       delay;
  int       workers;
  int       fanouts[4];
  PrCfgMode cfg;
  NocMode   noc;

  // --- Heterogeneous setup ---
  for (workers = 1; workers <= 512; workers <<= 1) {
  
    switch (workers) {
      case   1: cfg = PR_CFG_2_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case   2: cfg = PR_CFG_3_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case   4: cfg = PR_CFG_5_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case   8: cfg = PR_CFG_9_HET_FLAT;   noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  16: cfg = PR_CFG_17_HET_FLAT;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  32: cfg = PR_CFG_33_HET_FLAT;  noc = NOC_MODE_CREDIT_MAILBOX; break;
      case  64: cfg = PR_CFG_65_HET_FLAT;  noc = NOC_MODE_MAILBOX_ONLY;   break;
      case 128: cfg = PR_CFG_129_HET_FLAT; noc = NOC_MODE_MAILBOX_ONLY;   break;
      case 256: cfg = PR_CFG_257_HET_FLAT; noc = NOC_MODE_MAILBOX_ONLY;   break;
      case 512: cfg = PR_CFG_513_HET_FLAT; noc = NOC_MODE_MAILBOX_ONLY;   break;
      default: ar_abort();
    }
  
    for (delay = 8192; delay <= 8388608; delay <<= 1) {
      main_run_app(b, c, cfg, 0, 0, noc, 0, -1, NULL, NULL,
                   &test_myrmics_task_table, 512, delay, NULL);
    }
  }

  // --- Homogeneous setup ---
  for (workers = 1; workers <= 512; workers <<= 1) {

    fanouts[0] = (workers == 512) ? 511 : workers;
    noc        = (workers <= 32)  ? NOC_MODE_CREDIT_MAILBOX : 
                                    NOC_MODE_MAILBOX_ONLY;

    for (delay = 8192; delay <= 8388608; delay <<= 1) {
      main_run_app(b, c, PR_CFG_AUTO_MB, 2, fanouts, noc, 0, -1, NULL, NULL,
                   &test_myrmics_task_table, 
                   (workers == 512) ? 511 : 512, delay, NULL);
    }
  }
#endif

  // =========================================================================
  // Matrix multiplication
  // =========================================================================
#if 0
  // --- Strong scaling, variable task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_005f", NULL, &matrix_mult_myrmics_task_table, 
               1, 2, 384, -1, -1, NULL); // 665.935 sec (236.439 + ovflow)
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_009f", NULL, &matrix_mult_myrmics_task_table, 
               1, 3, 256, -1, -1, NULL); // 278.494 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_017f", NULL, &matrix_mult_myrmics_task_table, 
               1, 4, 192, -1, -1, NULL); // 113.475 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_035h", NULL, &matrix_mult_myrmics_task_table, 
               2, 3, 128, -1, -1, NULL); // 60.847 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_069h", NULL, &matrix_mult_myrmics_task_table, 
               2, 4, 96, -1, -1, NULL); // 26.104 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_136h", NULL, &matrix_mult_myrmics_task_table, 
               4, 3, 64, -1, 5, NULL); // 17.291 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "matrix_strong_264h", NULL, &matrix_mult_myrmics_task_table, 
               4, 4, 48, -1, 5, NULL); // 7.571 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               "matrix_strong_520h", NULL, &matrix_mult_myrmics_task_table, 
               4, 6, 32, -1, 5, NULL); // 5.621 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 6, 128, -1, -1, NULL); // 55.161 sec (vs 60.847 hier)
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 8, 96, -1, -1, NULL); // 26.925 sec (vs 26.104 hier)
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 12, 64, -1, -1, NULL); // 16.015 sec (vs 17.291 hier)
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 16, 48, -1, -1, NULL); // 12.333 sec (vs 7.571 hier)
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 24, 32, -1, 5, NULL); // 43.852 sec (vs 5.621 hier)
#endif
#if 0
  // --- Weak scaling, constant task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 2, 48, -1, -1, NULL); // 0.631 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 4, 48, -1, -1, NULL); // 1.564 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               2, 4, 48, -1, -1, NULL); // 3.153 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               4, 4, 48, -1, 5, NULL); // 7.533 sec

  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 8, 48, -1, -1, NULL); // 3.664 sec
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &matrix_mult_myrmics_task_table, 
               1, 16, 48, -1, -1, NULL); // 12.247 sec
#endif
#if 0
  // --- Impact of load-balancing vs. locality ---

  char buf[64];

  for (i = 0; i <= 100; i += 10) {
    kt_sprintf(buf, "matrix33_balance%03d.stats", i);
    main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, 
                 NOC_MODE_MAILBOX_ONLY, 0, i, buf, NULL,
                 &matrix_mult_myrmics_task_table, 
                 1, 6, 48, -1, -1, NULL);
  }
  //   i        time            DMAed data   Imbalance
  //   0 (loc)  62.953  sec     0.0          100.00%
  //  10        62.950  sec     0.0          100.00%
  //  20        62.950  sec     0.0          100.00%
  //  30        62.952  sec     0.0          100.00%
  //  40        62.951  sec     0.0          100.00%
  //  50         3.151  sec     211.0 K      8.04%
  //  60         2.875  sec     247.2 K      4.66%
  //  70         2.880  sec     246.6 K      3.81%
  //  80         3.153  sec     248.1 K      6.35%
  //  90         2.874  sec     245.5 K      5.08%
  // 100 (lbal)  2.552  sec     329.6 K      1.27%
#endif


  // =========================================================================
  // Jacobi
  // =========================================================================
#if 0
  // --- Strong scaling, variable task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_005f", NULL, &jacobi_myrmics_task_table, 
               2, 8, 4096, 1024, 16, -1, NULL); // 131.627 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_009f", NULL, &jacobi_myrmics_task_table, 
               2, 16, 4096, 1024, 16, -1, NULL); // 66.975 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_017f", NULL, &jacobi_myrmics_task_table,
               2, 32, 4096, 1024, 16, -1, NULL); // 38.914 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_035h", NULL, &jacobi_myrmics_task_table,
               4, 64, 4096, 1024, 16, -1, NULL); // 20.552 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_069h", NULL, &jacobi_myrmics_task_table, 
               8, 128, 4096, 1024, 16, -1, NULL); // 13.210 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_136h", NULL, &jacobi_myrmics_task_table, 
               14, 126, 4032, 1040, 16, -1, NULL); // 8.148 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "jacobi_strong_264h", NULL, &jacobi_myrmics_task_table, 
               14, 252, 4032, 1040, 16, -1, NULL); // 8.047 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               "jacobi_strong_520h", NULL, &jacobi_myrmics_task_table, 
               14, 504, 4032, 1040, 16, -1, NULL); // 12.553 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               2, 64, 4096, 1024, 16, -1, NULL); // 25.839 sec (vs 20.552 hier)
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               2, 128, 4096, 1024, 16, -1, NULL); // 22.300 sec (vs 13.210 hier)
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               1, 128, 4096, 1024, 16, -1, NULL); // 33.439 sec (vs 8.148 hier)
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               1, 256, 4096, 1024, 16, -1, NULL); // 48.029 sec (vs 8.047 hier)
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               1, 512, 4096, 1024, 16, -1, NULL); // 97.609 sec (vs 12.553 hier)
#endif
#if 0
  // --- Weak scaling, constant task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               2, 4, 32, 1024, 16, -1, NULL); // 1.547 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               2, 8, 64, 1024, 16, -1, NULL); // 1.643 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table,
               2, 16, 128, 1024, 16, -1, NULL); // 2.248 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table,
               4, 32, 256, 1024, 16, -1, NULL); // 2.337 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               8, 64, 512, 1024, 16, -1, NULL); // 2.791 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               14, 126, 1008, 1040, 16, -1, NULL); // 3.991 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               14, 252, 2016, 1040, 16, -1, NULL); // 6.413 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               14, 504, 4032, 1040, 16, -1, NULL); // 12.621 sec
               
  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               2, 32, 256, 1024, 16, -1, NULL); // 3.607 sec
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               2, 64, 512, 1024, 16, -1, NULL); // 7.019 sec
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               1, 128, 1024, 1024, 16, -1, NULL); // 20.038 sec
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               1, 256, 2048, 1024, 16, -1, NULL); // 42.417 sec
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &jacobi_myrmics_task_table, 
               1, 512, 4096, 1024, 16, -1, NULL); // 96.866 sec
#endif
#if 0
  // --- Impact of load-balancing vs. locality ---

  char buf[64];

  for (i = 100; i >= 0; i -= 10) {
    kt_sprintf(buf, "jacobi136_balance%03d.stats", i);
    main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, 
                 NOC_MODE_CREDIT_MAILBOX, 0, i, buf, NULL, 
                 &jacobi_myrmics_task_table, 
                 14, 126, 4032, 1040, 16, -1, NULL);
  }
  //   i        time            DMAed data   Imbalance
  //   0 (loc)  54.434 sec      702.2 K      93.55%
  //  10        54.427 sec      702.2 K      93.55%
  //  20        54.444 sec      702.2 K      93.55%
  //  30        54.450 sec      706.1 K      93.55%
  //  40        54.440 sec      5.3 M        64.07%
  //  50         8.700 sec      7.7 M        6.33%
  //  60         8.709 sec      7.8 M        7.23%
  //  70         7.678 sec      7.8 M        5.33%
  //  80         7.832 sec      7.9 M        4.61%
  //  90         7.765 sec      7.9 M        4.82%
  // 100 (lbal)  8.110 sec      8.5 M        1.47%
#endif

  
  // =========================================================================
  // Bitonic Sort
  // =========================================================================
#if 0
  // --- Strong scaling, variable task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_005f", NULL, &bitonic_myrmics_task_table,
               1, 4, 2097152, -1, -1, NULL); // 224.042 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_009f", NULL, &bitonic_myrmics_task_table,
               1, 8, 2097152, -1, -1, NULL); // 106.047 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_017f", NULL, &bitonic_myrmics_task_table,
               1, 16, 2097152, -1, -1, NULL); // 52.143 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_035h", NULL, &bitonic_myrmics_task_table,
               2, 32, 2097152, -1, -1, -1, NULL); // 26.117 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_069h", NULL, &bitonic_myrmics_task_table,
               4, 64, 2097152, -1, -1, NULL); // 17.280 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_136h", NULL, &bitonic_myrmics_task_table,
               8, 128, 2097152, -1, 6, NULL); // 15.986 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "bitonic_strong_264h", NULL, &bitonic_myrmics_task_table,
               8, 256, 2097152, -1, 6, NULL); // 14.217 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               "bitonic_strong_520h", NULL, &bitonic_myrmics_task_table,
               8, 512, 2097152, -1, 6, NULL); // 19.157 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 32, 2097152, -1, -1, NULL); // 28.434 sec (vs 26.117 hier)
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 64, 2097152, -1, -1, NULL); // 23.163 sec (vs 17.280 hier)
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 128, 2097152, -1, -1, NULL); // 26.157 sec (vs 15.986 hier)
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 256, 2097152, -1, -1, NULL); // 44.838 sec (vs 14.217 hier)
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 512, 2097152, -1, -1, NULL); // 115.935 sec (vs 19.157 hier)
#endif
#if 0
  // --- Weak scaling, constant task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 4, 16384, -1, -1, NULL); // 0.686 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 8, 32768, -1, -1, NULL); // 0.897 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 16, 65536, -1, -1, NULL); // 1.422 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               2, 32, 131072, -1, -1, -1, NULL); // 1.904 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               4, 64, 262144, -1, -1, NULL); // 3.020 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               8, 128, 524288, -1, 6, NULL); // 4.522 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               8, 256, 1048576, -1, 6, NULL); // 8.918 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               8, 512, 2097152, -1, 6, NULL); // 19.224 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 32, 131072, -1, -1, NULL); // 2.945 sec
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 64, 262144, -1, -1, NULL); // 6.674 sec
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 128, 524288, -1, -1, NULL); // 16.496 sec
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 256, 1048576, -1, -1, NULL); // 42.368 sec
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &bitonic_myrmics_task_table,
               1, 512, 2097152, -1, -1, NULL); // 115.951 sec
#endif


  // =========================================================================
  // FFT
  // =========================================================================
#if 0
  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats5b.txt", "paraver5.evt", &fft_myrmics_task_table,
               8, 1024, NULL); // 45.365 sec

  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats17a.txt", "paraver17a.evt", &fft_myrmics_task_table,
               16, 1024, NULL); // 14.356 sec

  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "stats17b.txt", "paraver17b.evt", &fft_myrmics_task_table,
               32, 1024, NULL); // 10.712 sec

  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_ONLY, 1024, -1,
               NULL, "paraver69.evt", &fft_myrmics_task_table,
               128, 1024, NULL);  // 128: 17.044 sec, 64: 7.709 sec (w. printfs)
  
  main_run_app(b, c, PR_CFG_120_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, "paraver120.evt", &fft_myrmics_task_table,
               256, 1024, NULL); // untested

  main_run_app(b, c, PR_CFG_232_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, "paraver232.evt", &fft_myrmics_task_table,
               512, 1024, NULL); // untested
#endif


  // =========================================================================
  // k-means
  //
  // NOTE: initialization is serial in k-means, so it takes a long time (few
  //       minutes) to start the parallel phases. TODO: fix this with parallel
  //       init in the application.
  // =========================================================================
#if 0
  // --- Strong scaling, variable task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_005f", NULL, &kmeans_myrmics_task_table, 
               16, 1, 4, 524288, 8, NULL); // 635.801 sec (206.305 + 1*ovflow)
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_009f", NULL, &kmeans_myrmics_task_table, 
               16, 1, 8, 262144, 8, NULL); // 323.131 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_017f", NULL, &kmeans_myrmics_task_table, 
               16, 1, 16, 131072, 8, NULL); // 166.097 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_035h", NULL, &kmeans_myrmics_task_table, 
               16, 2, 16, 65536, 8, NULL); // 83.547 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_069h", NULL, &kmeans_myrmics_task_table, 
               16, 4, 16, 32768, 8, NULL); // 42.408 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_136h", NULL, &kmeans_myrmics_task_table, 
               16, 7, 18, 16644, 8, NULL); // 22.514 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "kmeans_strong_264h", NULL, &kmeans_myrmics_task_table, 
               16, 7, 36, 8322, 8, NULL); // 13.081 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               "kmeans_strong_520h", NULL, &kmeans_myrmics_task_table, 
               16, 7, 73, 4104, 8, NULL); // 9.238 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 32, 65536, 8, NULL); // 88.071 sec (vs 83.547 hier)
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 64, 32768, 8, NULL); // 49.581 sec (vs 42.408 hier)
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 128, 16384, 8, NULL); // 31.350 sec (vs 22.514 hier)
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 256, 8192, 8, NULL); // 26.802 sec (vs 13.081 hier)
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 512, 4096, 8, NULL); // 34.044 sec (vs 9.238 hier)
#endif
#if 0
  // --- Weak scaling, constant task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 4, 4096, 8, NULL); // 5.227 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 8, 4096, 8, NULL); // 5.451 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 16, 4096, 8, NULL); // 5.903 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 2, 16, 4096, 8, NULL); // 5.819 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 4, 16, 4096, 8, NULL); // 5.931 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 7, 18, 4096, 8, NULL); // 6.318 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 7, 36, 4096, 8, NULL); // 7.172 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 7, 73, 4096, 8, NULL); // 9.307 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 32, 4096, 8, NULL); // 6.689 sec
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 64, 4096, 8, NULL); // 8.295 sec
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 128, 4096, 8, NULL); // 11.415 sec
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 256, 4096, 8, NULL); // 18.445 sec
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &kmeans_myrmics_task_table, 
               16, 1, 512, 4096, 8, NULL); // 34.013 sec
#endif
#if 0
  // --- Impact of load-balancing vs. locality ---

  char buf[64];

  for (i = 0; i <= 100; i += 10) {
    kt_sprintf(buf, "kmeans520_balance%03d.stats", i);
    main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, 
                 NOC_MODE_MAILBOX_ONLY, 0, i, buf, NULL,
                 &kmeans_myrmics_task_table, 
                 16, 7, 73, 4096, 8, NULL);
  }
  //   i        time            DMAed data   Imbalance
  //   0 (loc)  363.042 sec     61.5 K       98.83%
  //  10        363.054 sec     61.5 K       98.83%
  //  20        363.048 sec     61.5 K       98.83%
  //  30        363.049 sec     61.5 K       98.83%
  //  40        363.049 sec     61.5 K       98.83%
  //  50          8.789 sec     698.6 K      2.64%
  //  60          8.751 sec     708.0 K      2.64%
  //  70          8.702 sec     708.0 K      2.64%
  //  80          8.849 sec     708.0 K      2.64%
  //  90          8.747 sec     708.0 K      2.64%
  // 100 (lbal)   9.649 sec     1.0 M        2.11%
#endif

  // =========================================================================
  // c-ray
  // =========================================================================
#if 0
  // --- Strong scaling, variable task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_005f", NULL, &cray_myrmics_task_table,
               800, 512, 1, 4, NULL); // 157.712 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_009f", NULL, &cray_myrmics_task_table,
               800, 512, 1, 8, NULL); // 79.788 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_017f", NULL, &cray_myrmics_task_table,
               800, 512, 1, 16, NULL); // 44.016 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_035h", NULL, &cray_myrmics_task_table,
               800, 512, 2, 16, NULL); // 22.589 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_069h", NULL, &cray_myrmics_task_table,
               800, 512, 4, 16, NULL); // 11.428 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_136h", NULL, &cray_myrmics_task_table,
               800, 504, 7, 18, NULL); // 5.912 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               "cray_strong_264h", NULL, &cray_myrmics_task_table,
               800, 504, 7, 36, NULL); // 3.171 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               "cray_strong_520h", NULL, &cray_myrmics_task_table,
               800, 511, 7, 73, NULL); // 2.023 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 512, 1, 32, NULL); // 22.619 sec (vs 22.589 hier)
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 512, 1, 64, NULL); // 11.481 sec (vs 11.428 hier)
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 512, 1, 128, NULL); // 6.037 sec (vs 5.912 hier)
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 512, 1, 256, NULL); // 3.476 sec (vs 3.171 hier)
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 512, 1, 512, NULL); // 2.803 sec (vs 2.023 hier)
#endif
#if 0
  // --- Weak scaling, constant task size ---

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 4, 1, 4, NULL); // 1.273 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 8, 1, 8, NULL); // 1.304 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 16, 1, 16, NULL); // 1.368 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 32, 2, 16, NULL); // 1.462 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 64, 4, 16, NULL); // 1.546 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 126, 7, 18, NULL); // 1.623 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 252, 7, 36, NULL); // 1.763 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 511, 7, 73, NULL); // 2.023 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 32, 1, 32, NULL); // 1.462 sec
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 64, 1, 64, NULL); // 1.539 sec
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 128, 1, 128, NULL); // 1.694 sec
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 256, 1, 256, NULL); // 2.000 sec
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &cray_myrmics_task_table,
               800, 512, 1, 512, NULL); // 2.800 sec
#endif


  // =========================================================================
  // MD5
  // =========================================================================
#if 0
  // NOTE: initialization is serial in MD5, so it takes a long time (few
  //       minutes) to start the parallel phases. This can't be fixed without
  //       affecting the nature of the app -- it's supposed to be file I/O
  //       serial reading.

  main_run_app(b, c, PR_CFG_5_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 1, NULL); // 21.657 sec
  main_run_app(b, c, PR_CFG_9_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 1, NULL); // 11.755 sec
  main_run_app(b, c, PR_CFG_17_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 1, NULL); // 6.799 sec
  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 2, NULL); // 4.509 sec
  main_run_app(b, c, PR_CFG_69_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 4, NULL); // 3.266 sec
  main_run_app(b, c, PR_CFG_136_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               511, 98304, 7, NULL); // 2.608 sec
  main_run_app(b, c, PR_CFG_264_HET_HIER, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               511, 98304, 7, NULL); // 2.450 sec
  main_run_app(b, c, PR_CFG_520_HET_HIER2, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               511, 98304, 7, NULL); // 2.111 sec

  main_run_app(b, c, PR_CFG_33_HET_FLAT, 0, 0, NOC_MODE_CREDIT_MAILBOX, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 1, NULL); // 4.301 sec (vs 4.509 hier)
  main_run_app(b, c, PR_CFG_65_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               512, 98304, 1, NULL); // 3.293 sec (vs 3.266 hier)
  main_run_app(b, c, PR_CFG_129_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               511, 98304, 1, NULL); // 3.814 sec (vs 2.608 hier)
  main_run_app(b, c, PR_CFG_257_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               511, 98304, 1, NULL); // 4.491 sec (vs 2.450 hier)
  main_run_app(b, c, PR_CFG_513_HET_FLAT, 0, 0, NOC_MODE_MAILBOX_ONLY, 0, -1,
               NULL, NULL, &md5_myrmics_task_table,
               511, 98304, 1, NULL); // 3.801 sec (vs 2.111 hier)
#endif


  // =========================================================================
  // Barnes Hut
  // =========================================================================

  #if 1

  main_run_app(b, c, PR_CFG_35_HET_HIER, 0, 0, NOC_MODE_CREDIT_ONLY, 2048, -1,
                 NULL, NULL, &barnes_myrmics_task_table,
                 32, 512, 2, 10, NULL); // 
  #endif

  // =========================================================================
  // Multilevel benchmark. Use the parse_multilevel_*.awk scripts for results.
  //
  // WARNING: This is a synthetic benchmark for the 512-core Formic cube
  //          which uses an unusually high number of schedulers. If 
  //          various out-of-memory errors are encountered (kernel or user), 
  //          do the following in include/memory_management.h:
  //
  //          - Comment out MM_AGRESSIVE_PAGE_REQUESTS
  //          - define MM_LOCAL_FREE_PAGES_MAX  1
  //          - define MM_LOCAL_FREE_PAGES_MIN  1
  //          - define MM_MB_USER_SIZE          (85 * 1024 * 1024)
  //          - define MM_MB_VA_KERNEL_BASE     0x05800000
  //          - define MM_MB_KERNEL_SIZE        (5 * 1024 * 1024)
  // =========================================================================
#if 0

#define DELAY                   5000
#define REPS                    100
#define LEAF_TASKS              1
#define SIZE                    1
#define STRONG_PROBLEM_SIZE     (512 * DELAY)

  int fanouts[4];
  int tasks;

  // --- Strong scaling, variable task size ---
  for (i = 1; i <= 210; i = ((i * 3) + 1) / 2) {
    fanouts[0] = i;
    tasks = fanouts[0] * LEAF_TASKS;
    main_run_app(b, c, PR_CFG_AUTO_MB, 2, fanouts,
                 (i <= 32) ? NOC_MODE_CREDIT_MAILBOX : NOC_MODE_MAILBOX_ONLY, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 -1, -1, -1, LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, STRONG_PROBLEM_SIZE / tasks, NULL);
  }

  if ((b == AR_BOOT_MASTER_BID) && (!c)) kt_printf("\n\n\n@@@\r\n\n\n");

  for (i = 23; i <= 73; i += 10) {
    fanouts[0] = 6;
    fanouts[1] = i;
    tasks = fanouts[1] * fanouts[0] * LEAF_TASKS;
    main_run_app(b, c, PR_CFG_AUTO_MB, 3, fanouts,
                 (i <= 43) ? NOC_MODE_CREDIT_MAILBOX : NOC_MODE_MAILBOX_ONLY, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 -1, -1, fanouts[1], LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, STRONG_PROBLEM_SIZE / tasks, NULL);
  }

  if ((b == AR_BOOT_MASTER_BID) && (!c)) kt_printf("\n\n\n@@@\r\n\n\n");

  for (i = 6; i <= 11; i++) {
    fanouts[0] = 6;
    fanouts[1] = 6;
    fanouts[2] = i;
    tasks = fanouts[2] * fanouts[1] * fanouts[0] * LEAF_TASKS;
    main_run_app(b, c, PR_CFG_AUTO_MB, 4, fanouts,
                 NOC_MODE_CREDIT_MAILBOX, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 -1, fanouts[2], fanouts[1], LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, STRONG_PROBLEM_SIZE / tasks, NULL);
  }

  if ((b == AR_BOOT_MASTER_BID) && (!c)) kt_printf("\n\n\n@@@\r\n\n\n");

  for (i = 2; i <= 2; i++) {
    fanouts[0] = 6;
    fanouts[1] = 6;
    fanouts[2] = 5;
    fanouts[3] = i;
    tasks = fanouts[3] * fanouts[2] * fanouts[1] * fanouts[0] * LEAF_TASKS;
    main_run_app(b, c, PR_CFG_AUTO_MB, 5, fanouts,
                 NOC_MODE_CREDIT_MAILBOX, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 fanouts[3], fanouts[2], fanouts[1], 
                 LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, STRONG_PROBLEM_SIZE / tasks, NULL);
  }


  // --- Weak scaling, constant task size ---

  for (i = 1; i <= 210; i = ((i * 3) + 1) / 2) {
    fanouts[0] = i;
    main_run_app(b, c, PR_CFG_AUTO_MB, 2, fanouts,
                 (i <= 32) ? NOC_MODE_CREDIT_MAILBOX : NOC_MODE_MAILBOX_ONLY, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 -1, -1, -1, LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, DELAY, NULL);
  }

  if ((b == AR_BOOT_MASTER_BID) && (!c)) kt_printf("\n\n\n@@@\r\n\n\n");

  for (i = 23; i <= 73; i += 10) {
    fanouts[0] = 6;
    fanouts[1] = i;
    main_run_app(b, c, PR_CFG_AUTO_MB, 3, fanouts,
                 (i <= 43) ? NOC_MODE_CREDIT_MAILBOX : NOC_MODE_MAILBOX_ONLY, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 -1, -1, fanouts[1], LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, DELAY, NULL);
  }

  if ((b == AR_BOOT_MASTER_BID) && (!c)) kt_printf("\n\n\n@@@\r\n\n\n");

  for (i = 6; i <= 11; i++) {
    fanouts[0] = 6;
    fanouts[1] = 6;
    fanouts[2] = i;
    main_run_app(b, c, PR_CFG_AUTO_MB, 4, fanouts,
                 NOC_MODE_CREDIT_MAILBOX, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 -1, fanouts[2], fanouts[1], LEAF_TASKS * fanouts[0], 
                 SIZE, REPS, DELAY, NULL);
  }

  if ((b == AR_BOOT_MASTER_BID) && (!c)) kt_printf("\n\n\n@@@\r\n\n\n");

  for (i = 2; i <= 2; i++) {
    fanouts[0] = 6;
    fanouts[1] = 6;
    fanouts[2] = 5;
    fanouts[3] = i;
    main_run_app(b, c, PR_CFG_AUTO_MB, 5, fanouts,
                 NOC_MODE_CREDIT_MAILBOX, 0, -1,
                 NULL, NULL, &multilevel_myrmics_task_table,
                 fanouts[3], fanouts[2], fanouts[1], 
                 LEAF_TASKS * fanouts[0], SIZE, REPS, DELAY, NULL);
  }
#endif

  // Instruct server to kill the power
  if ((b == AR_BOOT_MASTER_BID) && (!c)) {
    ar_uart_flush();
    ar_timer_busy_wait_msec(200);
    kt_printf("%s", DBG_KILL_CONNECTION);
  }

  // Stall
  while (1) {
    ;
  }
  
}

