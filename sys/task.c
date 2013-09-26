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
// Abstract      : Put a short summary of the file purpose here
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: task.c,v $
// CVS revision  : $Revision: 1.10 $
// Last modified : $Date: 2012/10/22 13:19:33 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <syscall.h>
#include <debug.h>


// ===========================================================================
// _sys_spawn()                 Spawns a new task for parallel execution.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   unsigned int idx           Which task to run, specified as an index
//                              to a global function pointer table named
//                              void (*_sys_task_table[])
//   void **args                Array of in-order task arguments. All arguments
//                              must fit into a void * typecast.
//   unsigned int *types        Array of in-order argument types. Each type
//                              is a bitmask of SYS_TYPE_* properties as 
//                              defined in syscall.h header.
//   unsigned int num_args      Number of task arguments
// ===========================================================================
void _sys_spawn(char *filename, int line_nr, unsigned int idx, void **args, 
                unsigned int *types, unsigned int num_args) {

  Context       *context;
  PrMsgReq      *req;
  ListNode      *node;
  PrTaskDescr   *ptask;
  PrEvtPending  *event;
  int           i;
  int           j;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // What's the parent task?
  node = kt_list_head(context->pr_ready_queue);
  ar_assert(node);
  ptask = node->data;
  ar_assert(ptask);

  // Sanity checks
  ar_assert(idx < (1 << PR_TASK_IDX_SIZE));
  ar_assert(num_args > 0); // No task should be declared without a memory
                           // footprint. If for some reason this is desirable,
                           // the while loop below must be altered (it won't
                           // work for num_args == 0).

//kt_printf("%d: Parent task = 0x%X (idx = %d), spawning task idx = %d\r\n",
//          context->pr_core_id, ptask->id, ptask->index, idx);


  // Single-core mode?
  if (context->pr_num_cores == 1) {
    // FIXME: don't involve scheduler, ar_exec() directly...
    ar_abort();
  }


  // Make sure we haven't got any previous spawn request pending. We throttle
  // the spawn rate here, in order to avoid race conditions: we are allowed
  // to spawn a new task only when the scheduler has signalled (through a 
  // REPLY_SPAWN) that the old task spawn has progressed up to a safe point
  // which will not create any race.
  if (context->pr_spawn_pending) {
    ar_assert(!pr_event_worker_inner_loop(1, 0, 0, NULL));
  }

  // Mark that a new spawn request is now pending
  ar_assert(!context->pr_spawn_pending);
  context->pr_spawn_pending = 1;


  // Enter multi-part request creation loop
  i = 0;
  j = 0;
  req = NULL;
  while (i < num_args) {

    // New message part?
    if (!j) {

      // Get a request buffer
      req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));

      // Build new message
      req->core_id     = context->pr_core_id;
      req->req_id      = context->pr_message_id; // same ID for all parts
      req->type        = EXT_REQ_SPAWN;
      req->size        = 0;                 // init I/O type bitmap
      req->region      = 0;                 // init region/byvalue type bitmap
      req->ptr         = (void *) ((ptask->id << PR_TASK_IDX_SIZE) | idx);
      req->num_regions = (i != 0);          // new request part?
    }

    // Fill out next argument
    req->data[j] = args[i];

    req->size |= (types[i] & SYS_TYPE_IN_ARG)  ? (1 << (2 * j))     : 0;
    req->size |= (types[i] & SYS_TYPE_OUT_ARG) ? (1 << (2 * j + 1)) : 0;

    req->region |= (types[i] & SYS_TYPE_SAFE_ARG)   ? (1 << (2 * j))     : 0;
    req->region |= (types[i] & SYS_TYPE_REGION_ARG) ? (1 << (2 * j + 1)) : 0;

//kt_printf("Spawn arg %d: %s %s %s %s\r\n", i,
//    (types[i] & SYS_TYPE_IN_ARG) ? "IN" : "",
//    (types[i] & SYS_TYPE_OUT_ARG) ? "OUT" : "",
//    (types[i] & SYS_TYPE_SAFE_ARG) ? "SAFE" : "",
//    (types[i] & SYS_TYPE_REGION_ARG) ? "REGION" : ""
//    );

    i++;
    j++;

    // Finished with this part?
    if ((j == PR_REQ_MAX_SIZE) || (i == num_args)) {

      // More to follow?
      if ((j == PR_REQ_MAX_SIZE) && (i < num_args)) {
        req->num_ptrs = -1;
      }
      else {
        req->num_ptrs = j;
      }

      // Send message to scheduler
      ar_assert(!noc_msg_send());

      // Finished with this request buffer
      j = 0;
      req = NULL;
    }
  }

  // Create a note-to-self event to wait for the reply of this spawn request
  event = kt_malloc(sizeof(PrEvtPending));
  event->req = kt_malloc(sizeof(PrMsgReq));

  event->req->core_id = -1;
  event->req->req_id  = -1;
  event->req->type = SELF_WAIT_SPAWN;
  event->action = PR_ACT_REDO;
  event->prev = NULL;
  event->next = NULL;
  event->data = NULL;
  
  // Store event; we don't expect conflicts on this message ID.
  ar_assert(!kt_trie_insert(context->pr_pending_events, context->pr_message_id, 
                            event));

  // Increase message ID, avoiding value 0 on wrap-arounds
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
}


// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
void _sys_wait_on(char *filename, int line_nr, void **ptrs, 
                  unsigned int num_ptrs) {

  // FIXME
}

