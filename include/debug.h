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
// Abstract      : Header file for debugging functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: debug.h,v $
// CVS revision  : $Revision: 1.19 $
// Last modified : $Date: 2012/12/12 16:37:19 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _DEBUG_H
#define _DEBUG_H

#include <arch.h>
#include <memory_management.h>


// File dump magic strings for formicarium client/server
#define DBG_FILE_DUMP_BEGIN_FORMAT      "*#!=BEGIN_FILE %s\r\n"
#define DBG_FILE_DUMP_END_FORMAT        "*#!=END_FILE %s\r\n"

// Kill connection magic string
#define DBG_KILL_CONNECTION             "*#!=KILL\r\n"


// ===========================================================================
// Stats gatehring mechanism
// ===========================================================================

// Uncomment to enable stats gathering & reporting
#define DBG_STATS_ENABLED

#define DBG_STATS_MAGIC_HSHAKE1 0xDECA1001      // Handshakes
#define DBG_STATS_MAGIC_HSHAKE2 0xDECA2002
#define DBG_STATS_MAGIC_HSHAKE3 0xDECA3003

#define DBG_STATS_NUM_STATS             9       // Set this equal to the
                                                // number of stat types below
#define DBG_STATS_IDX_TIME_IDLE         0
#define DBG_STATS_IDX_TIME_TASK_EXEC    1
#define DBG_STATS_IDX_TIME_WORKER_WAIT  2
#define DBG_STATS_IDX_TIME_MEM_SERVE    3
#define DBG_STATS_IDX_TIME_SCH_SERVE    4
#define DBG_STATS_IDX_NUM_TASKS         5
#define DBG_STATS_IDX_NUM_MESSAGES      6
#define DBG_STATS_IDX_NUM_DMAS          7
#define DBG_STATS_IDX_DMA_TOTAL_SIZE    8

// Stats time taking macro
#ifdef DBG_STATS_ENABLED
static inline void dbg_stime(Context *context, int idx_to_charge) {
  unsigned int now;

  now = ar_free_timer_get_ticks();

  context->dbg_stats_data[idx_to_charge] += now - context->dbg_stats_last_tmr;
  context->dbg_stats_last_tmr            = now;
}
#else
#define dbg_stime(context, idx_to_charge)
#endif

// Stats count taking macro
#ifdef DBG_STATS_ENABLED
static inline void dbg_scount(Context *context, int idx_to_charge, int val) {
  context->dbg_stats_data[idx_to_charge] += val;
}
#else
#define dbg_scount(context, idx_to_charge, val)
#endif

// Exported functions
extern void dbg_stats_init();
extern void dbg_stats_report(char *filename);



// ===========================================================================
// Paraver tracing mechanism
// ===========================================================================

// Uncomment to enable Paraver tracing & reporting
#define DBG_TRC_ENABLED

// Uncomment to also enable communication traces
//#define DBG_TRC_COMM_ENABLED

#define DBG_TRC_MAX_ENTRIES     10000           // Mem needed: (8 * entries) B

#define DBG_TRC_MAGIC_HSHAKE1   0xBEEF4242      // Handshakes for timer offset
#define DBG_TRC_MAGIC_HSHAKE2   0xBEEF4343
#define DBG_TRC_MAGIC_HSHAKE3   0xBEEF4444
#define DBG_TRC_MAGIC_HSHAKE4   0xBEEF4545
#define DBG_TRC_MAGIC_HSHAKE5   0xBEEF4646

#define DBG_TRC_REQ_EXEC_BEGIN                  ((0x0 << 30) | 1)
#define DBG_TRC_REQ_EXEC_END                    ((0x1 << 30) | 1)
#define DBG_TRC_REQ_ALLOC_BEGIN                 ((0x0 << 30) | 2)
#define DBG_TRC_REQ_ALLOC_END                   ((0x1 << 30) | 2)
#define DBG_TRC_REQ_BALLOC_BEGIN                ((0x0 << 30) | 3)
#define DBG_TRC_REQ_BALLOC_END                  ((0x1 << 30) | 3)
#define DBG_TRC_REQ_FREE_BEGIN                  ((0x0 << 30) | 4)
#define DBG_TRC_REQ_FREE_END                    ((0x1 << 30) | 4)
#define DBG_TRC_REQ_RALLOC_BEGIN                ((0x0 << 30) | 5)
#define DBG_TRC_REQ_RALLOC_END                  ((0x1 << 30) | 5)
#define DBG_TRC_REQ_RALLOC_ORPHAN_BEGIN         ((0x0 << 30) | 6)
#define DBG_TRC_REQ_RALLOC_ORPHAN_END           ((0x1 << 30) | 6)
#define DBG_TRC_REQ_RFREE_BEGIN                 ((0x0 << 30) | 7)
#define DBG_TRC_REQ_RFREE_END                   ((0x1 << 30) | 7)
#define DBG_TRC_REQ_RFREE_UPDATE_PARENT_BEGIN   ((0x0 << 30) | 8)
#define DBG_TRC_REQ_RFREE_UPDATE_PARENT_END     ((0x1 << 30) | 8)
#define DBG_TRC_REQ_PACK_BEGIN                  ((0x0 << 30) | 9)
#define DBG_TRC_REQ_PACK_END                    ((0x1 << 30) | 9)
#define DBG_TRC_REQ_QUERY_POINTER_BEGIN         ((0x0 << 30) | 10)
#define DBG_TRC_REQ_QUERY_POINTER_END           ((0x1 << 30) | 10)
#define DBG_TRC_REQ_GET_PAGES_BEGIN             ((0x0 << 30) | 11)
#define DBG_TRC_REQ_GET_PAGES_END               ((0x1 << 30) | 11)
#define DBG_TRC_REQ_GET_RIDS_BEGIN              ((0x0 << 30) | 12)
#define DBG_TRC_REQ_GET_RIDS_END                ((0x1 << 30) | 12)
#define DBG_TRC_REQ_LOAD_REPORT_BEGIN           ((0x0 << 30) | 13)
#define DBG_TRC_REQ_LOAD_REPORT_END             ((0x1 << 30) | 13)
#define DBG_TRC_REQ_SHUTDOWN_BEGIN              ((0x0 << 30) | 14)
#define DBG_TRC_REQ_SHUTDOWN_END                ((0x1 << 30) | 14)
#define DBG_TRC_REQ_SPAWN_BEGIN                 ((0x0 << 30) | 15)
#define DBG_TRC_REQ_SPAWN_END                   ((0x1 << 30) | 15)
#define DBG_TRC_REQ_DELEGATE_BEGIN              ((0x0 << 30) | 16)
#define DBG_TRC_REQ_DELEGATE_END                ((0x1 << 30) | 16)
#define DBG_TRC_SELF_RALLOC_UPDATE_PARENT_BEGIN ((0x0 << 30) | 17)
#define DBG_TRC_SELF_RALLOC_UPDATE_PARENT_END   ((0x1 << 30) | 17)
#define DBG_TRC_SELF_PACK_MERGE_BEGIN           ((0x0 << 30) | 18)
#define DBG_TRC_SELF_PACK_MERGE_END             ((0x1 << 30) | 18)
#define DBG_TRC_SELF_RFREE_WAIT_CHILDREN_BEGIN  ((0x0 << 30) | 19)
#define DBG_TRC_SELF_RFREE_WAIT_CHILDREN_END    ((0x1 << 30) | 19)
#define DBG_TRC_REQ_DEP_START_BEGIN             ((0x0 << 30) | 20)
#define DBG_TRC_REQ_DEP_START_END               ((0x1 << 30) | 20)
#define DBG_TRC_REPLY_ALLOC_BEGIN               ((0x0 << 30) | 21)
#define DBG_TRC_REPLY_ALLOC_END                 ((0x1 << 30) | 21)
#define DBG_TRC_REPLY_BALLOC_BEGIN              ((0x0 << 30) | 22)
#define DBG_TRC_REPLY_BALLOC_END                ((0x1 << 30) | 22)
#define DBG_TRC_REPLY_FREE_BEGIN                ((0x0 << 30) | 23)
#define DBG_TRC_REPLY_FREE_END                  ((0x1 << 30) | 23)
#define DBG_TRC_REPLY_RALLOC_BEGIN              ((0x0 << 30) | 24)
#define DBG_TRC_REPLY_RALLOC_END                ((0x1 << 30) | 24)
#define DBG_TRC_REPLY_RFREE_BEGIN               ((0x0 << 30) | 25)
#define DBG_TRC_REPLY_RFREE_END                 ((0x1 << 30) | 25)
#define DBG_TRC_REPLY_QUERY_POINTER_BEGIN       ((0x0 << 30) | 26)
#define DBG_TRC_REPLY_QUERY_POINTER_END         ((0x1 << 30) | 26)
#define DBG_TRC_REPLY_PACK_BEGIN                ((0x0 << 30) | 27)
#define DBG_TRC_REPLY_PACK_END                  ((0x1 << 30) | 27)
#define DBG_TRC_REPLY_GET_PAGES_BEGIN           ((0x0 << 30) | 28)
#define DBG_TRC_REPLY_GET_PAGES_END             ((0x1 << 30) | 28)
#define DBG_TRC_REPLY_GET_RIDS_BEGIN            ((0x0 << 30) | 29)
#define DBG_TRC_REPLY_GET_RIDS_END              ((0x1 << 30) | 29)
#define DBG_TRC_REPLY_EXEC_BEGIN                ((0x0 << 30) | 30)
#define DBG_TRC_REPLY_EXEC_END                  ((0x1 << 30) | 30)
#define DBG_TRC_REQ_DEP_ROUTE_BEGIN             ((0x0 << 30) | 31)
#define DBG_TRC_REQ_DEP_ROUTE_END               ((0x1 << 30) | 31)
#define DBG_TRC_REQ_EXEC_DONE_BEGIN             ((0x0 << 30) | 32)
#define DBG_TRC_REQ_EXEC_DONE_END               ((0x1 << 30) | 32)
#define DBG_TRC_WORKER_WAIT_BEGIN               ((0x0 << 30) | 33)
#define DBG_TRC_WORKER_WAIT_END                 ((0x1 << 30) | 33)
#define DBG_TRC_REQ_DEP_ENQUEUE_BEGIN           ((0x0 << 30) | 34)
#define DBG_TRC_REQ_DEP_ENQUEUE_END             ((0x1 << 30) | 34)
#define DBG_TRC_REQ_DEP_OK_BEGIN                ((0x0 << 30) | 35)
#define DBG_TRC_REQ_DEP_OK_END                  ((0x1 << 30) | 35)
#define DBG_TRC_REQ_DEP_STOP_BEGIN              ((0x0 << 30) | 36)
#define DBG_TRC_REQ_DEP_STOP_END                ((0x1 << 30) | 36)
#define DBG_TRC_REQ_DEP_CHILD_FREE_BEGIN        ((0x0 << 30) | 37)
#define DBG_TRC_REQ_DEP_CHILD_FREE_END          ((0x1 << 30) | 37)
#define DBG_TRC_REQ_UPDATE_LOCATION_BEGIN       ((0x0 << 30) | 38)
#define DBG_TRC_REQ_UPDATE_LOCATION_END         ((0x1 << 30) | 38)
#define DBG_TRC_SELF_SCHEDULE_RESULT_BEGIN      ((0x0 << 30) | 39)
#define DBG_TRC_SELF_SCHEDULE_RESULT_END        ((0x1 << 30) | 39)
#define DBG_TRC_REQ_SCHEDULE_BEGIN              ((0x0 << 30) | 40)
#define DBG_TRC_REQ_SCHEDULE_END                ((0x1 << 30) | 40)
#define DBG_TRC_REPLY_SCHEDULE_BEGIN            ((0x0 << 30) | 41)
#define DBG_TRC_REPLY_SCHEDULE_END              ((0x1 << 30) | 41)
#define DBG_TRC_REPLY_SPAWN_BEGIN               ((0x0 << 30) | 42)
#define DBG_TRC_REPLY_SPAWN_END                 ((0x1 << 30) | 42)
#define DBG_TRC_SELF_WAIT_SPAWN_BEGIN           ((0x0 << 30) | 43)
#define DBG_TRC_SELF_WAIT_SPAWN_END             ((0x1 << 30) | 43)

// NoC send trace:    ((0x2 << 30) | (dst_core_id << 20) | (req_id & 0xFFFFF)
// NoC recv trace:    ((0x3 << 30) | (src_core_id << 20) | (req_id & 0xFFFFF)


// Trace taking macro
#ifdef DBG_TRC_ENABLED
static inline void dbg_trace(Context *context, unsigned int trace) {
  if (context->dbg_trc_idx > 2 * DBG_TRC_MAX_ENTRIES - 2) {
    return;
  }
  context->dbg_trc_data[context->dbg_trc_idx++] = ar_free_timer_get_ticks();
  context->dbg_trc_data[context->dbg_trc_idx++] = trace;
}
#else
#define dbg_trace(context, trace)
#endif

// Exported functions
extern void dbg_trace_init();
extern void dbg_trace_report(char *filename);


#endif
