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
// Abstract      : Statistics gathering & reporting functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: stats.c,v $
// CVS revision  : $Revision: 1.3 $
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
void dbg_stats_init() {
#ifdef DBG_STATS_ENABLED

  Context       *context;
  unsigned int  i;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Allocate array
  ar_assert(!context->dbg_stats_data);
  context->dbg_stats_data = kt_malloc(DBG_STATS_NUM_STATS * sizeof(int));

  // Zero them all
  for (i = 0; i < DBG_STATS_NUM_STATS; i++) {
    context->dbg_stats_data[i] = 0;
  }

  // Init the "last timer" value to something meaningful
  context->dbg_stats_last_tmr = ar_free_timer_get_ticks();

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
void dbg_stats_format(unsigned int val, int *ret_dec, int *ret_frac, 
                      char *ret_unit) {

  if (val < 1000) {
    *ret_dec = val;
    *ret_frac = 0;
    *ret_unit = ' ';
  }
  else if (val < 1000000) {
    *ret_dec = val / 1000;
    *ret_frac = (val % 1000) / 100;
    *ret_unit = 'K';
  }
  else if (val < 1000000000) {
    *ret_dec = val / 1000000;
    *ret_frac = (val % 1000000) / 100000;
    *ret_unit = 'M';
  }
  else {
    *ret_dec = val / 1000000000;
    *ret_frac = (val % 1000000000) / 100000000;
    *ret_unit = 'G';
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
void dbg_stats_format_number(char *buf, unsigned int val) {

  char *s;
  int  i;


  s = buf + kt_int_log10(val);
  s += (s - buf) / 3 + 1;
  
  *s-- = 0;

  i = 0;
  do {
    *s-- = (val % 10) + '0';
    if (!(++i % 3)) {
      *s-- = ',';
    }
    val /= 10;
  } while (val);
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
void dbg_stats_report(char *filename) {
#ifdef DBG_STATS_ENABLED

#define SUMMARY_BUF_SIZE 8196

  Context       *context;
  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  unsigned int  **stats = NULL;
  char          *buf = NULL;
  char          fmt_buf[16];
  char          *s;
  unsigned int  i;
  unsigned int  j;

  unsigned int  sch_idle_tot = 0;
  unsigned int  sch_idle_avg = 0;
  int           sch_idle_avg_dec = 0;
  int           sch_idle_avg_frac = 0;
  char          sch_idle_avg_units = 0;
  unsigned int  sch_idle_min = UINT_MAX;
  int           sch_idle_min_idx = -1;
  int           sch_idle_min_dec = 0;
  int           sch_idle_min_frac = 0;
  char          sch_idle_min_units = ' ';
  unsigned int  sch_idle_max = 0;
  int           sch_idle_max_idx = -1;
  int           sch_idle_max_dec = 0;
  int           sch_idle_max_frac = 0;
  char          sch_idle_max_units = 0;

  unsigned int  wrk_idle_tot = 0;
  unsigned int  wrk_idle_avg = 0;
  int           wrk_idle_avg_dec = 0;
  int           wrk_idle_avg_frac = 0;
  char          wrk_idle_avg_units = 0;
  unsigned int  wrk_idle_min = UINT_MAX;
  int           wrk_idle_min_idx = -1;
  int           wrk_idle_min_dec = 0;
  int           wrk_idle_min_frac = 0;
  char          wrk_idle_min_units = ' ';
  unsigned int  wrk_idle_max = 0;
  int           wrk_idle_max_idx = -1;
  int           wrk_idle_max_dec = 0;
  int           wrk_idle_max_frac = 0;
  char          wrk_idle_max_units = 0;

  unsigned int  sch_mem_tot = 0;
  unsigned int  sch_mem_avg = 0;
  int           sch_mem_avg_dec = 0;
  int           sch_mem_avg_frac = 0;
  char          sch_mem_avg_units = 0;
  unsigned int  sch_mem_min = UINT_MAX;
  int           sch_mem_min_idx = -1;
  int           sch_mem_min_dec = 0;
  int           sch_mem_min_frac = 0;
  char          sch_mem_min_units = ' ';
  unsigned int  sch_mem_max = 0;
  int           sch_mem_max_idx = -1;
  int           sch_mem_max_dec = 0;
  int           sch_mem_max_frac = 0;
  char          sch_mem_max_units = 0;

  unsigned int  wrk_work_tot = 0;
  unsigned int  wrk_work_avg = 0;
  int           wrk_work_avg_dec = 0;
  int           wrk_work_avg_frac = 0;
  char          wrk_work_avg_units = 0;
  unsigned int  wrk_work_min = UINT_MAX;
  int           wrk_work_min_idx = -1;
  int           wrk_work_min_dec = 0;
  int           wrk_work_min_frac = 0;
  char          wrk_work_min_units = ' ';
  unsigned int  wrk_work_max = 0;
  int           wrk_work_max_idx = -1;
  int           wrk_work_max_dec = 0;
  int           wrk_work_max_frac = 0;
  char          wrk_work_max_units = 0;

  unsigned int  sch_sched_tot = 0;
  unsigned int  sch_sched_avg = 0;
  int           sch_sched_avg_dec = 0;
  int           sch_sched_avg_frac = 0;
  char          sch_sched_avg_units = 0;
  unsigned int  sch_sched_min = UINT_MAX;
  int           sch_sched_min_idx = -1;
  int           sch_sched_min_dec = 0;
  int           sch_sched_min_frac = 0;
  char          sch_sched_min_units = ' ';
  unsigned int  sch_sched_max = 0;
  int           sch_sched_max_idx = -1;
  int           sch_sched_max_dec = 0;
  int           sch_sched_max_frac = 0;
  char          sch_sched_max_units = 0;

  unsigned int  wrk_wait_tot = 0;
  unsigned int  wrk_wait_avg = 0;
  int           wrk_wait_avg_dec = 0;
  int           wrk_wait_avg_frac = 0;
  char          wrk_wait_avg_units = 0;
  unsigned int  wrk_wait_min = UINT_MAX;
  int           wrk_wait_min_idx = -1;
  int           wrk_wait_min_dec = 0;
  int           wrk_wait_min_frac = 0;
  char          wrk_wait_min_units = ' ';
  unsigned int  wrk_wait_max = 0;
  int           wrk_wait_max_idx = -1;
  int           wrk_wait_max_dec = 0;
  int           wrk_wait_max_frac = 0;
  char          wrk_wait_max_units = 0;

  unsigned int  sch_tasks_tot = 0;
  unsigned int  sch_tasks_avg = 0;
  int           sch_tasks_avg_dec = 0;
  int           sch_tasks_avg_frac = 0;
  char          sch_tasks_avg_units = 0;
  unsigned int  sch_tasks_min = UINT_MAX;
  int           sch_tasks_min_idx = -1;
  int           sch_tasks_min_dec = 0;
  int           sch_tasks_min_frac = 0;
  char          sch_tasks_min_units = ' ';
  unsigned int  sch_tasks_max = 0;
  int           sch_tasks_max_idx = -1;
  int           sch_tasks_max_dec = 0;
  int           sch_tasks_max_frac = 0;
  char          sch_tasks_max_units = 0;

  unsigned int  wrk_tasks_tot = 0;
  unsigned int  wrk_tasks_avg = 0;
  int           wrk_tasks_avg_dec = 0;
  int           wrk_tasks_avg_frac = 0;
  char          wrk_tasks_avg_units = 0;
  unsigned int  wrk_tasks_min = UINT_MAX;
  int           wrk_tasks_min_idx = -1;
  int           wrk_tasks_min_dec = 0;
  int           wrk_tasks_min_frac = 0;
  char          wrk_tasks_min_units = ' ';
  unsigned int  wrk_tasks_max = 0;
  int           wrk_tasks_max_idx = -1;
  int           wrk_tasks_max_dec = 0;
  int           wrk_tasks_max_frac = 0;
  char          wrk_tasks_max_units = 0;

  unsigned int  sch_msg_tot = 0;
  unsigned int  sch_msg_avg = 0;
  int           sch_msg_avg_dec = 0;
  int           sch_msg_avg_frac = 0;
  char          sch_msg_avg_units = 0;
  unsigned int  sch_msg_min = UINT_MAX;
  int           sch_msg_min_idx = -1;
  int           sch_msg_min_dec = 0;
  int           sch_msg_min_frac = 0;
  char          sch_msg_min_units = ' ';
  unsigned int  sch_msg_max = 0;
  int           sch_msg_max_idx = -1;
  int           sch_msg_max_dec = 0;
  int           sch_msg_max_frac = 0;
  char          sch_msg_max_units = 0;

  unsigned int  wrk_msg_tot = 0;
  unsigned int  wrk_msg_avg = 0;
  int           wrk_msg_avg_dec = 0;
  int           wrk_msg_avg_frac = 0;
  char          wrk_msg_avg_units = 0;
  unsigned int  wrk_msg_min = UINT_MAX;
  int           wrk_msg_min_idx = -1;
  int           wrk_msg_min_dec = 0;
  int           wrk_msg_min_frac = 0;
  char          wrk_msg_min_units = ' ';
  unsigned int  wrk_msg_max = 0;
  int           wrk_msg_max_idx = -1;
  int           wrk_msg_max_dec = 0;
  int           wrk_msg_max_frac = 0;
  char          wrk_msg_max_units = 0;

  unsigned int  wrk_ndma_tot = 0;
  unsigned int  wrk_ndma_avg = 0;
  int           wrk_ndma_avg_dec = 0;
  int           wrk_ndma_avg_frac = 0;
  char          wrk_ndma_avg_units = 0;
  unsigned int  wrk_ndma_min = UINT_MAX;
  int           wrk_ndma_min_idx = -1;
  int           wrk_ndma_min_dec = 0;
  int           wrk_ndma_min_frac = 0;
  char          wrk_ndma_min_units = ' ';
  unsigned int  wrk_ndma_max = 0;
  int           wrk_ndma_max_idx = -1;
  int           wrk_ndma_max_dec = 0;
  int           wrk_ndma_max_frac = 0;
  char          wrk_ndma_max_units = 0;

  unsigned int  wrk_sdma_tot = 0;
  unsigned int  wrk_sdma_avg = 0;
  int           wrk_sdma_avg_dec = 0;
  int           wrk_sdma_avg_frac = 0;
  char          wrk_sdma_avg_units = 0;
  unsigned int  wrk_sdma_min = UINT_MAX;
  int           wrk_sdma_min_idx = -1;
  int           wrk_sdma_min_dec = 0;
  int           wrk_sdma_min_frac = 0;
  char          wrk_sdma_min_units = ' ';
  unsigned int  wrk_sdma_max = 0;
  int           wrk_sdma_max_idx = -1;
  int           wrk_sdma_max_dec = 0;
  int           wrk_sdma_max_frac = 0;
  char          wrk_sdma_max_units = 0;


  // Get context
  context = mm_get_context(ar_get_core_id());
  my_cid = ar_get_core_id();
  my_bid = ar_get_board_id();


  // =========================================================================
  // Gather stats from everybody
  // =========================================================================

  // Top-level scheduler is the master core for reporting
  if (context->pr_parent_sched_id == -1) {

    // Allocate space for all
    stats = kt_malloc(context->pr_num_cores * sizeof(unsigned int *));
    for (i = 0; i < context->pr_num_cores; i++) {
      stats[i] = kt_malloc(DBG_STATS_NUM_STATS * sizeof(unsigned int));
    }
    
    // For all cores in the setup
    for (i = 0; i < context->pr_num_cores; i++) {

      // Get arch-level board/core ID
      pr_core_arch_bid_cid(i, &bid, &cid);

      // Is it us?
      if ((bid == my_bid) && (cid == my_cid)) {

        // Copy our own stats
        for (j = 0; j < DBG_STATS_NUM_STATS; j++) {
          stats[i][j] = context->dbg_stats_data[j];
        }

        continue;
      }

      // Handshake with peer and tell him to send us his stats; send our
      // bid/cid so he can communicate back
      ar_mbox_send(my_cid, bid, cid, DBG_STATS_MAGIC_HSHAKE1);
      ar_mbox_send(my_cid, bid, cid, (my_bid << 8) | my_cid);
      
      // Sanity check: get start-of-transmission
      ar_assert(ar_mbox_get(my_cid) == DBG_STATS_MAGIC_HSHAKE2);

      // Copy his stats
      for (j = 0; j < DBG_STATS_NUM_STATS; j++) {
        stats[i][j] = ar_mbox_get(my_cid);
      }
      
      // Sanity check: get end-of-transmission
      ar_assert(ar_mbox_get(my_cid) == DBG_STATS_MAGIC_HSHAKE3);
    }
  }

  // Other cores are slaves
  else {

    // Receive master command and his bid/cid
    ar_assert(ar_mbox_get(my_cid) == DBG_STATS_MAGIC_HSHAKE1);
    i = ar_mbox_get(my_cid);
    bid = i >> 8;
    cid = i & 0xFF;

    // Send start-of-transmission
    ar_mbox_send(my_cid, bid, cid, DBG_STATS_MAGIC_HSHAKE2);

    // Send our stats
    for (i = 0; i < DBG_STATS_NUM_STATS; i++) {
      ar_mbox_send(my_cid, bid, cid, context->dbg_stats_data[i]);
    }

    // Send end-of-transmission
    ar_mbox_send(my_cid, bid, cid, DBG_STATS_MAGIC_HSHAKE3);
  }


  // =========================================================================
  // Create a summary and print it
  // =========================================================================
  
  // Top-level scheduler only, everybody else get out
  if (context->pr_parent_sched_id != -1) {
    return;
  }

  // Scheduler & worker idle time
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_sched_ids[i] > -1) {
      sch_idle_tot += stats[i][DBG_STATS_IDX_TIME_IDLE];
      if (stats[i][DBG_STATS_IDX_TIME_IDLE] < sch_idle_min) {
        sch_idle_min = stats[i][DBG_STATS_IDX_TIME_IDLE];
        sch_idle_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_TIME_IDLE] > sch_idle_max) {
        sch_idle_max = stats[i][DBG_STATS_IDX_TIME_IDLE];
        sch_idle_max_idx = i;
      }
    }
    else {
      wrk_idle_tot += stats[i][DBG_STATS_IDX_TIME_IDLE];
      if (stats[i][DBG_STATS_IDX_TIME_IDLE] < wrk_idle_min) {
        wrk_idle_min = stats[i][DBG_STATS_IDX_TIME_IDLE];
        wrk_idle_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_TIME_IDLE] > wrk_idle_max) {
        wrk_idle_max = stats[i][DBG_STATS_IDX_TIME_IDLE];
        wrk_idle_max_idx = i;
      }
    }
  }
  sch_idle_avg = sch_idle_tot / context->pr_num_schedulers;
  wrk_idle_avg = wrk_idle_tot / context->pr_num_workers;
  
  dbg_stats_format(sch_idle_avg, &sch_idle_avg_dec, &sch_idle_avg_frac, 
                   &sch_idle_avg_units);
  dbg_stats_format(sch_idle_min, &sch_idle_min_dec, &sch_idle_min_frac, 
                   &sch_idle_min_units);
  dbg_stats_format(sch_idle_max, &sch_idle_max_dec, &sch_idle_max_frac, 
                   &sch_idle_max_units);

  dbg_stats_format(wrk_idle_avg, &wrk_idle_avg_dec, &wrk_idle_avg_frac, 
                   &wrk_idle_avg_units);
  dbg_stats_format(wrk_idle_min, &wrk_idle_min_dec, &wrk_idle_min_frac, 
                   &wrk_idle_min_units);
  dbg_stats_format(wrk_idle_max, &wrk_idle_max_dec, &wrk_idle_max_frac, 
                   &wrk_idle_max_units);

  // Scheduler mem time
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_sched_ids[i] > -1) {
      sch_mem_tot += stats[i][DBG_STATS_IDX_TIME_MEM_SERVE];
      if (stats[i][DBG_STATS_IDX_TIME_MEM_SERVE] < sch_mem_min) {
        sch_mem_min = stats[i][DBG_STATS_IDX_TIME_MEM_SERVE];
        sch_mem_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_TIME_MEM_SERVE] > sch_mem_max) {
        sch_mem_max = stats[i][DBG_STATS_IDX_TIME_MEM_SERVE];
        sch_mem_max_idx = i;
      }
    }
  }
  sch_mem_avg = sch_mem_tot / context->pr_num_schedulers;
  
  dbg_stats_format(sch_mem_avg, &sch_mem_avg_dec, &sch_mem_avg_frac, 
                   &sch_mem_avg_units);
  dbg_stats_format(sch_mem_min, &sch_mem_min_dec, &sch_mem_min_frac, 
                   &sch_mem_min_units);
  dbg_stats_format(sch_mem_max, &sch_mem_max_dec, &sch_mem_max_frac, 
                   &sch_mem_max_units);

  // Worker work time
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_work_ids[i] > -1) {
      wrk_work_tot += stats[i][DBG_STATS_IDX_TIME_TASK_EXEC];
      if (stats[i][DBG_STATS_IDX_TIME_TASK_EXEC] < wrk_work_min) {
        wrk_work_min = stats[i][DBG_STATS_IDX_TIME_TASK_EXEC];
        wrk_work_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_TIME_TASK_EXEC] > wrk_work_max) {
        wrk_work_max = stats[i][DBG_STATS_IDX_TIME_TASK_EXEC];
        wrk_work_max_idx = i;
      }
    }
  }
  wrk_work_avg = wrk_work_tot / context->pr_num_workers;
  
  dbg_stats_format(wrk_work_avg, &wrk_work_avg_dec, &wrk_work_avg_frac, 
                   &wrk_work_avg_units);
  dbg_stats_format(wrk_work_min, &wrk_work_min_dec, &wrk_work_min_frac, 
                   &wrk_work_min_units);
  dbg_stats_format(wrk_work_max, &wrk_work_max_dec, &wrk_work_max_frac, 
                   &wrk_work_max_units);

  // Scheduler non-mem time
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_sched_ids[i] > -1) {
      sch_sched_tot += stats[i][DBG_STATS_IDX_TIME_SCH_SERVE];
      if (stats[i][DBG_STATS_IDX_TIME_SCH_SERVE] < sch_sched_min) {
        sch_sched_min = stats[i][DBG_STATS_IDX_TIME_SCH_SERVE];
        sch_sched_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_TIME_SCH_SERVE] > sch_sched_max) {
        sch_sched_max = stats[i][DBG_STATS_IDX_TIME_SCH_SERVE];
        sch_sched_max_idx = i;
      }
    }
  }
  sch_sched_avg = sch_sched_tot / context->pr_num_schedulers;
  
  dbg_stats_format(sch_sched_avg, &sch_sched_avg_dec, &sch_sched_avg_frac, 
                   &sch_sched_avg_units);
  dbg_stats_format(sch_sched_min, &sch_sched_min_dec, &sch_sched_min_frac, 
                   &sch_sched_min_units);
  dbg_stats_format(sch_sched_max, &sch_sched_max_dec, &sch_sched_max_frac, 
                   &sch_sched_max_units);

  // Worker wait time
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_work_ids[i] > -1) {
      wrk_wait_tot += stats[i][DBG_STATS_IDX_TIME_WORKER_WAIT];
      if (stats[i][DBG_STATS_IDX_TIME_WORKER_WAIT] < wrk_wait_min) {
        wrk_wait_min = stats[i][DBG_STATS_IDX_TIME_WORKER_WAIT];
        wrk_wait_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_TIME_WORKER_WAIT] > wrk_wait_max) {
        wrk_wait_max = stats[i][DBG_STATS_IDX_TIME_WORKER_WAIT];
        wrk_wait_max_idx = i;
      }
    }
  }
  wrk_wait_avg = wrk_wait_tot / context->pr_num_workers;
  
  dbg_stats_format(wrk_wait_avg, &wrk_wait_avg_dec, &wrk_wait_avg_frac, 
                   &wrk_wait_avg_units);
  dbg_stats_format(wrk_wait_min, &wrk_wait_min_dec, &wrk_wait_min_frac, 
                   &wrk_wait_min_units);
  dbg_stats_format(wrk_wait_max, &wrk_wait_max_dec, &wrk_wait_max_frac, 
                   &wrk_wait_max_units);

  // Scheduler & worker local tasks
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_sched_ids[i] > -1) {
      sch_tasks_tot += stats[i][DBG_STATS_IDX_NUM_TASKS];
      if (stats[i][DBG_STATS_IDX_NUM_TASKS] < sch_tasks_min) {
        sch_tasks_min = stats[i][DBG_STATS_IDX_NUM_TASKS];
        sch_tasks_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_NUM_TASKS] > sch_tasks_max) {
        sch_tasks_max = stats[i][DBG_STATS_IDX_NUM_TASKS];
        sch_tasks_max_idx = i;
      }
    }
    else {
      wrk_tasks_tot += stats[i][DBG_STATS_IDX_NUM_TASKS];
      if (stats[i][DBG_STATS_IDX_NUM_TASKS] < wrk_tasks_min) {
        wrk_tasks_min = stats[i][DBG_STATS_IDX_NUM_TASKS];
        wrk_tasks_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_NUM_TASKS] > wrk_tasks_max) {
        wrk_tasks_max = stats[i][DBG_STATS_IDX_NUM_TASKS];
        wrk_tasks_max_idx = i;
      }
    }
  }
  sch_tasks_avg = sch_tasks_tot / context->pr_num_schedulers;
  wrk_tasks_avg = wrk_tasks_tot / context->pr_num_workers;
  
  dbg_stats_format(sch_tasks_avg, &sch_tasks_avg_dec, &sch_tasks_avg_frac, 
                   &sch_tasks_avg_units);
  dbg_stats_format(sch_tasks_min, &sch_tasks_min_dec, &sch_tasks_min_frac, 
                   &sch_tasks_min_units);
  dbg_stats_format(sch_tasks_max, &sch_tasks_max_dec, &sch_tasks_max_frac, 
                   &sch_tasks_max_units);

  dbg_stats_format(wrk_tasks_avg, &wrk_tasks_avg_dec, &wrk_tasks_avg_frac, 
                   &wrk_tasks_avg_units);
  dbg_stats_format(wrk_tasks_min, &wrk_tasks_min_dec, &wrk_tasks_min_frac, 
                   &wrk_tasks_min_units);
  dbg_stats_format(wrk_tasks_max, &wrk_tasks_max_dec, &wrk_tasks_max_frac, 
                   &wrk_tasks_max_units);

  // Scheduler & worker number of messages
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_sched_ids[i] > -1) {
      sch_msg_tot += stats[i][DBG_STATS_IDX_NUM_MESSAGES];
      if (stats[i][DBG_STATS_IDX_NUM_MESSAGES] < sch_msg_min) {
        sch_msg_min = stats[i][DBG_STATS_IDX_NUM_MESSAGES];
        sch_msg_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_NUM_MESSAGES] > sch_msg_max) {
        sch_msg_max = stats[i][DBG_STATS_IDX_NUM_MESSAGES];
        sch_msg_max_idx = i;
      }
    }
    else {
      wrk_msg_tot += stats[i][DBG_STATS_IDX_NUM_MESSAGES];
      if (stats[i][DBG_STATS_IDX_NUM_MESSAGES] < wrk_msg_min) {
        wrk_msg_min = stats[i][DBG_STATS_IDX_NUM_MESSAGES];
        wrk_msg_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_NUM_MESSAGES] > wrk_msg_max) {
        wrk_msg_max = stats[i][DBG_STATS_IDX_NUM_MESSAGES];
        wrk_msg_max_idx = i;
      }
    }
  }
  sch_msg_avg = sch_msg_tot / context->pr_num_schedulers;
  wrk_msg_avg = wrk_msg_tot / context->pr_num_workers;
  
  dbg_stats_format(sch_msg_avg, &sch_msg_avg_dec, &sch_msg_avg_frac, 
                   &sch_msg_avg_units);
  dbg_stats_format(sch_msg_min, &sch_msg_min_dec, &sch_msg_min_frac, 
                   &sch_msg_min_units);
  dbg_stats_format(sch_msg_max, &sch_msg_max_dec, &sch_msg_max_frac, 
                   &sch_msg_max_units);

  dbg_stats_format(wrk_msg_avg, &wrk_msg_avg_dec, &wrk_msg_avg_frac, 
                   &wrk_msg_avg_units);
  dbg_stats_format(wrk_msg_min, &wrk_msg_min_dec, &wrk_msg_min_frac, 
                   &wrk_msg_min_units);
  dbg_stats_format(wrk_msg_max, &wrk_msg_max_dec, &wrk_msg_max_frac, 
                   &wrk_msg_max_units);

  // Worker number of DMAs
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_work_ids[i] > -1) {
      wrk_ndma_tot += stats[i][DBG_STATS_IDX_NUM_DMAS];
      if (stats[i][DBG_STATS_IDX_NUM_DMAS] < wrk_ndma_min) {
        wrk_ndma_min = stats[i][DBG_STATS_IDX_NUM_DMAS];
        wrk_ndma_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_NUM_DMAS] > wrk_ndma_max) {
        wrk_ndma_max = stats[i][DBG_STATS_IDX_NUM_DMAS];
        wrk_ndma_max_idx = i;
      }
    }
  }
  wrk_ndma_avg = wrk_ndma_tot / context->pr_num_workers;
  
  dbg_stats_format(wrk_ndma_avg, &wrk_ndma_avg_dec, &wrk_ndma_avg_frac, 
                   &wrk_ndma_avg_units);
  dbg_stats_format(wrk_ndma_min, &wrk_ndma_min_dec, &wrk_ndma_min_frac, 
                   &wrk_ndma_min_units);
  dbg_stats_format(wrk_ndma_max, &wrk_ndma_max_dec, &wrk_ndma_max_frac, 
                   &wrk_ndma_max_units);

  // Worker DMA size
  for (i = 0; i < context->pr_num_cores; i++) {
    if (context->pr_core_work_ids[i] > -1) {
      wrk_sdma_tot += stats[i][DBG_STATS_IDX_DMA_TOTAL_SIZE];
      if (stats[i][DBG_STATS_IDX_DMA_TOTAL_SIZE] < wrk_sdma_min) {
        wrk_sdma_min = stats[i][DBG_STATS_IDX_DMA_TOTAL_SIZE];
        wrk_sdma_min_idx = i;
      }
      if (stats[i][DBG_STATS_IDX_DMA_TOTAL_SIZE] > wrk_sdma_max) {
        wrk_sdma_max = stats[i][DBG_STATS_IDX_DMA_TOTAL_SIZE];
        wrk_sdma_max_idx = i;
      }
    }
  }
  wrk_sdma_avg = wrk_sdma_tot / context->pr_num_workers;
  
  dbg_stats_format(wrk_sdma_avg, &wrk_sdma_avg_dec, &wrk_sdma_avg_frac, 
                   &wrk_sdma_avg_units);
  dbg_stats_format(wrk_sdma_min, &wrk_sdma_min_dec, &wrk_sdma_min_frac, 
                   &wrk_sdma_min_units);
  dbg_stats_format(wrk_sdma_max, &wrk_sdma_max_dec, &wrk_sdma_max_frac, 
                   &wrk_sdma_max_units);



  // Print summary to buffer
  buf = kt_malloc(SUMMARY_BUF_SIZE * sizeof(char));
  s = buf;
  s += kt_sprintf(s, 
"=============================================================================\r\n"
"                            Statistics Summary\r\n"
"=============================================================================\r\n"
"\r\n"
"             Schedulers                               Workers\r\n"
"          ----------------                         -------------\r\n"
"\r\n"
"Idle time:    %3d.%d %c (avg)             Idle time:    %3d.%d %c (avg)\r\n"
"              %3d.%d %c (min, core %3d)                 %3d.%d %c (min, core %3d)\r\n"
"              %3d.%d %c (max, core %3d)                 %3d.%d %c (max, core %3d)\r\n"
"\r\n"
"Mem time:     %3d.%d %c (avg)             Work time:    %3d.%d %c (avg)\r\n"
"              %3d.%d %c (min, core %3d)                 %3d.%d %c (min, core %3d)\r\n"
"              %3d.%d %c (max, core %3d)                 %3d.%d %c (max, core %3d)\r\n"
"\r\n"
"Non-mem time: %3d.%d %c (avg)             Wait time:    %3d.%d %c (avg)\r\n"
"              %3d.%d %c (min, core %3d)                 %3d.%d %c (min, core %3d)\r\n"
"              %3d.%d %c (max, core %3d)                 %3d.%d %c (max, core %3d)\r\n"
"\r\n"
"Local tasks:  %3d.%d %c (avg)             Local tasks:  %3d.%d %c (avg)\r\n"
"              %3d.%d %c (min, core %3d)                 %3d.%d %c (min, core %3d)\r\n"
"              %3d.%d %c (max, core %3d)                 %3d.%d %c (max, core %3d)\r\n"
"\r\n"
"Num messages: %3d.%d %c (avg)             Num messages: %3d.%d %c (avg)\r\n"
"              %3d.%d %c (min, core %3d)                 %3d.%d %c (min, core %3d)\r\n"
"              %3d.%d %c (max, core %3d)                 %3d.%d %c (max, core %3d)\r\n"
"\r\n"
"                                        Num DMAs:     %3d.%d %c (avg)\r\n"
"                                                      %3d.%d %c (min, core %3d)\r\n"
"                                                      %3d.%d %c (max, core %3d)\r\n"
"\r\n"
"                                        DMAed data:   %3d.%d %c (avg)\r\n"
"                                                      %3d.%d %c (min, core %3d)\r\n"
"                                                      %3d.%d %c (max, core %3d)\r\n"
"=============================================================================\r\n",
    sch_idle_avg_dec, sch_idle_avg_frac, sch_idle_avg_units,
    wrk_idle_avg_dec, wrk_idle_avg_frac, wrk_idle_avg_units,
    sch_idle_min_dec, sch_idle_min_frac, sch_idle_min_units, sch_idle_min_idx,
    wrk_idle_min_dec, wrk_idle_min_frac, wrk_idle_min_units, wrk_idle_min_idx,
    sch_idle_max_dec, sch_idle_max_frac, sch_idle_max_units, sch_idle_max_idx,
    wrk_idle_max_dec, wrk_idle_max_frac, wrk_idle_max_units, wrk_idle_max_idx,

    sch_mem_avg_dec, sch_mem_avg_frac, sch_mem_avg_units,
    wrk_work_avg_dec, wrk_work_avg_frac, wrk_work_avg_units,
    sch_mem_min_dec, sch_mem_min_frac, sch_mem_min_units, sch_mem_min_idx,
    wrk_work_min_dec, wrk_work_min_frac, wrk_work_min_units, wrk_work_min_idx,
    sch_mem_max_dec, sch_mem_max_frac, sch_mem_max_units, sch_mem_max_idx,
    wrk_work_max_dec, wrk_work_max_frac, wrk_work_max_units, wrk_work_max_idx,

    sch_sched_avg_dec, sch_sched_avg_frac, sch_sched_avg_units,
    wrk_wait_avg_dec, wrk_wait_avg_frac, wrk_wait_avg_units,
    sch_sched_min_dec, sch_sched_min_frac, sch_sched_min_units, sch_sched_min_idx,
    wrk_wait_min_dec, wrk_wait_min_frac, wrk_wait_min_units, wrk_wait_min_idx,
    sch_sched_max_dec, sch_sched_max_frac, sch_sched_max_units, sch_sched_max_idx,
    wrk_wait_max_dec, wrk_wait_max_frac, wrk_wait_max_units, wrk_wait_max_idx,

    sch_tasks_avg_dec, sch_tasks_avg_frac, sch_tasks_avg_units,
    wrk_tasks_avg_dec, wrk_tasks_avg_frac, wrk_tasks_avg_units,
    sch_tasks_min_dec, sch_tasks_min_frac, sch_tasks_min_units, sch_tasks_min_idx,
    wrk_tasks_min_dec, wrk_tasks_min_frac, wrk_tasks_min_units, wrk_tasks_min_idx,
    sch_tasks_max_dec, sch_tasks_max_frac, sch_tasks_max_units, sch_tasks_max_idx,
    wrk_tasks_max_dec, wrk_tasks_max_frac, wrk_tasks_max_units, wrk_tasks_max_idx,

    sch_msg_avg_dec, sch_msg_avg_frac, sch_msg_avg_units,
    wrk_msg_avg_dec, wrk_msg_avg_frac, wrk_msg_avg_units,
    sch_msg_min_dec, sch_msg_min_frac, sch_msg_min_units, sch_msg_min_idx,
    wrk_msg_min_dec, wrk_msg_min_frac, wrk_msg_min_units, wrk_msg_min_idx,
    sch_msg_max_dec, sch_msg_max_frac, sch_msg_max_units, sch_msg_max_idx,
    wrk_msg_max_dec, wrk_msg_max_frac, wrk_msg_max_units, wrk_msg_max_idx,

    wrk_ndma_avg_dec, wrk_ndma_avg_frac, wrk_ndma_avg_units,
    wrk_ndma_min_dec, wrk_ndma_min_frac, wrk_ndma_min_units, wrk_ndma_min_idx,
    wrk_ndma_max_dec, wrk_ndma_max_frac, wrk_ndma_max_units, wrk_ndma_max_idx,

    wrk_sdma_avg_dec, wrk_sdma_avg_frac, wrk_sdma_avg_units,
    wrk_sdma_min_dec, wrk_sdma_min_frac, wrk_sdma_min_units, wrk_sdma_min_idx,
    wrk_sdma_max_dec, wrk_sdma_max_frac, wrk_sdma_max_units, wrk_sdma_max_idx

  );
  if (s - buf > SUMMARY_BUF_SIZE) {
    ar_panic("Summary buffer overflow");
  }

  // Print it
  ar_uart_flush();
  ar_timer_busy_wait_msec(200);
  kt_printf("\r\n%s\r\n", buf);


  // =========================================================================
  // Dump both the summary and all the per-core analytical stats to a file
  // =========================================================================
  
  // Print file begin dump header
  ar_uart_flush();
  ar_timer_busy_wait_msec(200);
  kt_printf(DBG_FILE_DUMP_BEGIN_FORMAT, filename);
  ar_uart_flush();
  ar_timer_busy_wait_msec(20);

  // Print the summary
  kt_printf("%s\r\n", buf);

  // For all cores in the setup
  for (i = 0; i < context->pr_num_cores; i++) {
    
    // Print the analytical stats
    kt_printf("==================================================\r\n");
    kt_printf("Analytical Statistics for Core %d [%s]\r\n", 
              i,
              (context->pr_core_sched_ids[i] > -1) ? "Scheduler" : 
                                                     "Worker");
    kt_printf("==================================================\r\n");
    for (j = 0; j < DBG_STATS_NUM_STATS; j++) {
      // Print relevant fields
      if ((context->pr_core_sched_ids[i] > -1) &&
           ((j == 1) || (j == 2) || (j == 7) || (j == 8))) {
        continue;
      }
      dbg_stats_format_number(fmt_buf, stats[i][j]);
      kt_printf("%s %15s\r\n",
          (j == 0) ? "Idle time:    " :
          (j == 1) ? "Work time:    " :
          (j == 2) ? "Wait time:    " :
          (j == 3) ? "Mem time:     " :
          (j == 4) ? "Non-mem time: " :
          (j == 5) ? "Local tasks:  " :
          (j == 6) ? "Num messages: " :
          (j == 7) ? "Num DMAs:     " :
          (j == 8) ? "DMAed data:   " : 
                     "!!!ERROR!!!",
          fmt_buf);
    }
    kt_printf("==================================================\r\n\r\n");
  }

  // Print file end dump header
  ar_uart_flush();
  ar_timer_busy_wait_msec(200);
  kt_printf(DBG_FILE_DUMP_END_FORMAT, filename);
  ar_uart_flush();
  ar_timer_busy_wait_msec(20);

  // Free stuff
  kt_free(stats);
  kt_free(buf);


#endif
}
