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
// ============================================================================
// The code in this file is derived from the MMA version 1.0 project, which 
// was licensed under the following copyright:
//
// Copyright (c) 2011, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by Spyros Lyberis <lymperis1@llnl.gov>
// LLNL-CODE-637218, OCEC-13-184
// All rights reserved.
// 
// This file is part of MMA, version 1.0. 
// For details, please see http://myrmics.com/download.php
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// Additional BSD Notice
// 
// 1. This notice is required to be provided under our contract with the U.S.
//    Department of Energy (DOE). This work was produced at Lawrence Livermore
//    National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
// 2. Neither the United States Government nor Lawrence Livermore National
//    Security, LLC nor any of their employees, makes any warranty, express or
//    implied, or assumes any liability or responsibility for the accuracy,
//    completeness, or usefulness of any information, apparatus, product, or
//    process disclosed, or represents that its use would not infringe
//    privately-owned rights.
// 3. Also, reference herein to any specific commercial products, process, or
//    services by trade name, trademark, manufacturer or otherwise does not
//    necessarily constitute or imply its endorsement, recommendation, or
//    favoring by the United States Government or Lawrence Livermore National
//    Security, LLC. The views and opinions of authors expressed herein do not
//    necessarily state or reflect those of the United States Government or
//    Lawrence Livermore National Security, LLC, and shall not be used for
//    advertising or product endorsement purposes.
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : Task preparation, dispatching and monitoring functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: task.c,v $
// CVS revision  : $Revision: 1.27 $
// Last modified : $Date: 2013/02/12 11:20:30 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <syscall.h>
#include <debug.h>


// ===========================================================================
// pr_task_create()             Create a new task descriptor for a task,
//                              fill it up and store it in the scheduler's
//                              task Trie.
// ===========================================================================
// * INPUTS
//   unsigned int par_task_id   Task ID of the parent task
//   int task_idx               Index of task function in the task table
//   void **args                Arguments of the task function call. The
//                              array will be hooked to the task descriptor
//                              and will be freed by pr_task_collect().
//   unsigned char *deps        Dependence flags for each arg. The
//                              array will be hooked to the task descriptor
//                              and will be freed by pr_task_collect().
//   int num_args               Number of task function call arguments and
//                              dependencies
//
// * RETURN VALUE
//   PrTaskDescr *              Task descriptor entry of the new task
// ===========================================================================
PrTaskDescr *pr_task_create(unsigned int par_task_id, int task_idx, 
                            void **args, unsigned char *deps, int num_args) {

  Context       *context;
  PrTaskDescr   *task;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(task_idx < (1 << PR_TASK_IDX_SIZE));

  // Take a note for the statistics gathering
  dbg_scount(context, DBG_STATS_IDX_NUM_TASKS, 1);

  // Create task descriptor
  task = kt_malloc(sizeof(PrTaskDescr *));
  task->pid      = par_task_id;
  task->id       = (context->pr_scheduler_id << 
                                 (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS)) | 
                   context->pr_avail_task_id;
  task->index    = task_idx;
  task->args     = args;
  task->deps     = deps;
  task->num_args = num_args;
  if (!num_args) {
    ar_assert(!task_idx);   // special case: main task can be argument-less,
    ar_assert(!task->args); //               but normal tasks can not: they
    ar_assert(!task->deps); //               should own some memory...
  }

  // Count and mark how many dependencies are there. We do that for all tasks,
  // except the main function (who won't wait on anything to start).
  task->deps_waiting = 0;
  task->dmas_started = 0;
  task->dmas_waiting = 0;
  if (task->index) {
    for (i = 0; i < task->num_args; i++) {
      // We only ignore (a) safe, by-value objects and (b) NULL pointers, which
      // we mark as by-value, because the compiler may have not done that (e.g.
      // if they sometimes are NULL, but some other times non-NULL). We don't
      // support safe regions or other variations yet.
      if (!task->args[i]) {
        // Force by-value of NULL pointer. If the assertion fails, it's 
        // probably a user error (he passed the NULL region as a dep).
        ar_assert(!(task->deps[i] & SYS_TYPE_REGION_ARG));
        task->deps[i] = SYS_TYPE_BYVALUE_ARG;
      }
      if (task->deps[i] != SYS_TYPE_BYVALUE_ARG) {
        task->deps[i] |= SYS_TYPE_DEP_PENDING;
        task->deps_waiting++;
      }
    }
  }

  // No packed adr/sizes to DMA yet; they will be filled later when the
  // region/object packing calls return.
  task->dma_addresses = NULL;
  task->dma_sizes = NULL;
  task->num_dmas = 0;

  // Task is not scheduled yet
  task->run_core_id = -1;

  // No active spawns
  task->spawn_msg_id = 0;
  task->spawn_core_id = -1;
  task->spawn_rem_deps = 0;

  // Insert task descriptor into the tasks trie
  ar_assert(!kt_trie_insert(context->pr_tasks, task->id, task));

//kt_printf("%d: created task 0x%X, idx %d, args %d, deps_waiting %d\r\n", context->pr_core_id, task->id, task_idx, task->num_args, task->deps_waiting);

  // Advance next available task ID. We have to look if there's another one
  // with this ID in the trie. This is highly unlikely: it will occur after
  // 2 ^ (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS) task creations, if the 
  // task has not finished. However, there are some cases this may happen,
  // e.g. at the top-level scheduler which will run the main task for the
  // whole program duration.
  do {
    context->pr_avail_task_id = pr_advance_task_id(context->pr_avail_task_id);
  } while (kt_trie_find(context->pr_tasks, 
                        (context->pr_scheduler_id <<
                                   (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS)) | 
                        context->pr_avail_task_id,
                        NULL));

  // Note that we have one more scheduled task to be responsible for, and
  // this adds to our scheduled task load. Possibly report upstream.
  context->pr_cur_sched_load++;
  pr_sched_report_load();

  // Success
  return task;
}


// ===========================================================================
// pr_task_dispatch()           Sends a task that is ready to run (i.e. its
//                              dependencies are resolved and it is scheduled)
//                              towards a worker core for execution.
// ===========================================================================
// * INPUTS
//   PrTaskDescr *task          Task descriptor of task to be dispatched
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_task_dispatch(PrTaskDescr *task) {

  Context       *context;
  PrMsgReq      *req;
  int           child;
  int           msg_id;
  int           cur_dma;
  int           total;
  int           i;
  int           j;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(task);
  ar_assert(task->id < (1 << PR_TASK_ID_SIZE));
  ar_assert((task->run_core_id > -1) && 
            (task->run_core_id < context->pr_num_cores));

  // Find out which subtree of ours we'll be forwarding the task to
  child = pr_core_child_index(context->pr_core_route[task->run_core_id]);

//kt_printf("%d: Dispatching task 0x%X to %d\r\n", context->pr_core_id, task->id, context->pr_children[child]);

  // Single message ID for all cases (and all parts of the multi-request)
  msg_id = context->pr_message_id;

  // Special case: num_args == 0 (should be only for dispatching the main task)
  if (!task->num_args) {
    ar_assert(!task->index);

    // Get a request buffer
    req = noc_msg_send_get_buf(context->pr_children[child]);

    // Build single message
    req->core_id     = context->pr_core_id;
    req->req_id      = msg_id;
    req->type        = EXT_REQ_EXEC;
    req->region      = task->run_core_id;
    req->ptr         = (void *) ((task->id << PR_TASK_IDX_SIZE) | 
                                 task->index);
    req->size = 0;
    req->num_regions = 0;
    req->num_ptrs = 0;

    // Send message to scheduler
    ar_assert(!noc_msg_send());
  }

  // Normal case: multi-part request creation loop
  else {
    i = 0;
    j = 0;
    req = NULL;
    total = task->num_args + 2 * task->num_dmas;

    while (i < total) {

      // New message part?
      if (!j) {

        // Get a request buffer
        req = noc_msg_send_get_buf(context->pr_children[child]);

        // Build new message
        req->core_id     = context->pr_core_id;
        req->req_id      = msg_id;
        req->type        = EXT_REQ_EXEC;
        req->region      = task->run_core_id;
        req->ptr         = (void *) ((task->id << PR_TASK_IDX_SIZE) | 
                                     task->index);
        req->size        = task->num_args;
        req->num_regions = (i != 0);          // new request part?
      }

      // Fill out next argument. We first put all task arguments and then
      // we follow by the DMAs, two args per DMA: one is the address and the
      // other is the size & location & pack_options.
      if (i < task->num_args) {
        req->data[j] = task->args[i];
      }
      else {
        cur_dma = i - task->num_args;
        if (cur_dma % 2 == 0) {
          req->data[j] = (void *) task->dma_addresses[cur_dma / 2];
        }
        else {
          req->data[j] = (void *) task->dma_sizes[cur_dma / 2];
        }
      }
      i++;
      j++;

      // Finished with this part?
      if ((j == PR_REQ_MAX_SIZE) || (i == total)) {

        // More to follow?
        if ((j == PR_REQ_MAX_SIZE) && (i < total)) {
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
  }

  // Increase message ID, avoiding value 0 on wrap-arounds
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // When dispatching directly to a worker child, mark the child task load
  // change, increase our task load and possibly report upstream.
  if (!context->pr_scheduler_level) {
    context->pr_chld_run_load[child]++;
    context->pr_cur_run_load++;
    pr_sched_report_load();
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_task_collect()            Handle a task that it's finished execution
// ===========================================================================
// * INPUTS
//   unsigned int task_id       Task ID of the finished task
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_task_collect(unsigned int task_id) {

  Context       *context;
  PrTaskDescr   *task;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(task_id < (1 << PR_TASK_ID_SIZE));

  // Locate task
  task = NULL;
  kt_trie_find(context->pr_tasks, task_id, (void *) &task);
  ar_assert(task);
  ar_assert(task->id == task_id);

//kt_printf("%d: Collecting task 0x%08X\r\n", context->pr_core_id, task->id);

  // If the main task ended, just take a note
  if (!task->index) {
    ar_assert(context->pr_parent_sched_id == -1);
    ar_assert(!context->pr_main_finished);
    context->pr_main_finished = 1;
  }
  // For normal tasks, notify any next tasks blocked on our args
  else {
    ar_assert(!pr_dep_start_stop(0, 0, task->id, task->num_args, task->args,
                                 task->deps));
  }

  // Free task argument/deps arrays, DMA structures and task descriptor and
  // remove the latter from the trie
  kt_free(task->args);
  kt_free(task->deps);
  kt_free(task->dma_addresses);
  kt_free(task->dma_sizes);
  kt_trie_delete(context->pr_tasks, task_id, kt_free);

  // Note that we have one less scheduled task to be responsible for, and
  // remove it from our scheduled task load. Possibly report upstream.
  ar_assert(context->pr_cur_sched_load > 0);
  context->pr_cur_sched_load--;
  pr_sched_report_load();

  // Success
  return 0;
}


// ===========================================================================
// pr_task_must_delegate()      Decides if a task we're getting from a spawn
//                              request (or a delegation of such a request
//                              which comes from our (grand)parent) must
//                              be further delegated to a child scheduler of
//                              ours, or it must be handled locally.
//
//                              The criterion we're using is to check all task
//                              arguments if they're local or which child
//                              subtree owns them. If ALL TASK ARGUMENTS belong
//                              to the SAME child, we choose to delegate the
//                              task to this scheduler. In this way, we
//                              guarantee that the scheduler which owns a task
//                              will never have to make upcalls: all
//                              objects/regions will belong to either the task
//                              scheduler or to some of its
//                              children/grandchildren.
// ===========================================================================
// * INPUTS
//   PrMultiPartReq *multi      The multi-part request minus the final part,
//                              or NULL if only a final part exists
//   PrMsgReq *final_part       The final part of the multi-part request
//
// * OUTPUTS
//   int *ret_child             Selected child to delegate the task to,
//                              if return value is 1
//
// * RETURN VALUE
//   int                        1: yes, delegate the task to *ret_child
//                              0: no, do not delegate the task
// ===========================================================================
int pr_task_must_delegate(PrMultiPartReq *multi, PrMsgReq *final_part, 
                          int *ret_child) {
  
  Context       *context;
  int           sched_id;
  int           i;
  int           j;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // No candidate child yet
  *ret_child = -1;

  // Do stored arguments first
  if (multi) {
    for (i = 0; i < multi->num_parts; i++) {
      ar_assert(multi->parts[i]->num_ptrs == -1);
      for (j = 0; j < PR_REQ_MAX_SIZE; j++) {
        // Ignore by-value arguments (SYS_TYPE_BYVALUE_ARG: in, !out, safe and
        // !region flags)
        if ( (multi->parts[i]->size   & (1 << (2 * j))) &&
            !(multi->parts[i]->size   & (1 << (2 * j + 1))) &&
             (multi->parts[i]->region & (1 << (2 * j))) &&
            !(multi->parts[i]->region & (1 << (2 * j + 1)))) {
          continue;
        }
        // Also ignore NULL pointers (NULL value, !region flag)
        if (!multi->parts[i]->data[j]  &&
            !(multi->parts[i]->region & (1 << (2 * j + 1)))) {
          continue;
        }
        
        // Region?
        if (multi->parts[i]->region & (1 << (2 * j + 1))) {
          sched_id = mm_distr_region_sched_id(
                                      (rid_t) multi->parts[i]->data[j]);
        }
        // Object
        else {
          sched_id = mm_distr_object_sched_id(multi->parts[i]->data[j]);
        }
        // Object/region should exist and it shouldn't be our parent's. 
        //
        // TODO this is most probably a user error: add error reporting 
        // here instead of aborting.
        ar_assert(sched_id >= 0);
        ar_assert(sched_id != context->pr_parent_sched_id);

        // If it's ours, we can't delegate
        if (sched_id == context->pr_scheduler_id) {
          return 0;
        }
        // If we didn't have a candidate, we just did.
        if (*ret_child == -1) {
          *ret_child = sched_id;
          continue;
        }
        // Did we have another candidate? If yes, we can't delegate.
        if (*ret_child != sched_id) {
          return 0;
        }
      }
    }
  }

  // Repeat the exact same steps for the last part of the request
  for (j = 0; j < final_part->num_ptrs; j++) {
    if ( (final_part->size   & (1 << (2 * j))) &&
        !(final_part->size   & (1 << (2 * j + 1))) &&
         (final_part->region & (1 << (2 * j))) &&
        !(final_part->region & (1 << (2 * j + 1)))) {
      continue;
    }
    if (!final_part->data[j]  &&
        !(final_part->region & (1 << (2 * j + 1)))) {
      continue;
    }
    if (final_part->region & (1 << (2 * j + 1))) {
      sched_id = mm_distr_region_sched_id((rid_t) final_part->data[j]);
    }
    else {
      sched_id = mm_distr_object_sched_id(final_part->data[j]);
    }
    ar_assert(sched_id >= 0);
    ar_assert(sched_id != context->pr_parent_sched_id);
    if (sched_id == context->pr_scheduler_id) {
      return 0;
    }
    if (*ret_child == -1) {
      *ret_child = sched_id;
      continue;
    }
    if (*ret_child != sched_id) {
      return 0;
    }
  }

  // If we reached this point, delegate to the selected candidate
  ar_assert(*ret_child != -1);
  return 1;
}


// ===========================================================================
// pr_task_scheduled()          When a task has been scheduled, this function
//                              updates the location of the arguments and
//                              dispatches it for execution. 
// ===========================================================================
// * INPUTS
//   PrTaskDescr *task          The newly scheduled task
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_task_scheduled(PrTaskDescr *task) {

  Context       *context;
  int           src_core_id;
  int           something_to_move;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(task);
  ar_assert((task->run_core_id >= 0) && 
            (task->run_core_id < context->pr_num_cores));

  // Dispatch the task for execution
  pr_task_dispatch(task);

  // See if anything will be moved from its current location
  something_to_move = 0;
  for (i = 0; i < task->num_dmas; i++) {
    src_core_id = (task->dma_sizes[i] >> MM_PACK_SIZE_BITS) &
                  ((1 << MM_PACK_LOCATION_BITS) - 1);
    if (src_core_id != task->run_core_id) {
      something_to_move = 1;
      break;
    }
  }
//kt_printf("Task 0x%08X scheduled at %d, data move = %d\r\n", task->id, task->run_core_id, something_to_move);

  // Update the locations of all task arguments, possibly generating 
  // messages downstream
  if (something_to_move) {
    ar_assert(!pr_sched_update_location(task->run_core_id, task->num_args,
                                        task->args, task->deps, 0));
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_task_start_dmas()         Starts a new DMA for each of the task argument
//                              packed list, from the last producer of the
//                              argument to the current worker core that the
//                              task will execute
// ===========================================================================
// * INPUTS
//   PrTaskDescr *task          The newly scheduled task
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_task_start_dmas(PrTaskDescr *task) {

  Context       *context;
  int           src_core_id;
  int           size;
  int           options;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(task);
  //ar_assert((task->run_core_id >= 0) && 
  //          (task->run_core_id < context->pr_num_cores));
  ar_assert(task->run_core_id == context->pr_core_id); // All DMAs to here
  ar_assert(!task->dmas_started);


  // Initialize software counter for group DMA notification
  task->dmas_waiting = task->num_dmas;

//kt_printf("%d: Starting %d DMAs for task 0x%08X, var 0x%08X=%d\r\n", context->pr_core_id, task->num_dmas, task->id, &task->dmas_waiting, task->dmas_waiting);

  // Start a new DMA for each packed address/size pair
  for (i = 0; i < task->num_dmas; i++) {

    // Degenerate case? Or, has the scheduler done a good job? :)
    src_core_id = (task->dma_sizes[i] >> MM_PACK_SIZE_BITS) &
                  ((1 << MM_PACK_LOCATION_BITS) - 1);
    if (src_core_id == task->run_core_id) {
      // No DMA
      task->dmas_waiting--;
//kt_printf("> %d: 0x%08X %d->%d CANCELLED\r\n", i, task->dma_addresses[i], src_core_id, task->run_core_id);
    }
    else {

      // Start a DMA
      size = task->dma_sizes[i] & ((1 << MM_PACK_SIZE_BITS) - 1);
      options = (task->dma_sizes[i] >> (MM_PACK_SIZE_BITS + 
                                        MM_PACK_LOCATION_BITS)) & 
                ((1 << MM_PACK_OPTION_BITS) - 1);

      ar_assert((size > 0) && (size <= AR_DMA_MAX_SIZE) && 
                !(size & (AR_DMA_ALIGN - 1)));

      ar_assert((options == MM_PACK_OPTION_RW) || 
                (options == MM_PACK_OPTION_RO));

//kt_printf("%d> %d: 0x%08X %d->%d %dB [%s]\r\n", context->pr_core_id, i, task->dma_addresses[i], src_core_id, task->run_core_id, size, (options == MM_PACK_OPTION_RW) ? "RW" : "RO");

      // For the DMA I/W/C flags, we have two cases:
      //
      // (a) Read/Write dependencies: We set I=1, W=0, C=0.
      //     
      //     We just transfer dirty lines from a last producer with read/write
      //     access to a new producer with read/write access. Since there is
      //     noone else with read-only access, the last producer is always the
      //     one to write-back the dirty line on the DRAM if the cache gets
      //     full. Upon a transfer, the dirty flag must be forgotten by the
      //     source core, so if the line gets evicted and happens to be at the
      //     same board, no stale writeback happens.
      //
      // (b) Read-only dependencies: We set I=0, W=1, C=1.
      //
      //     We allow the dirty lines to remain to the last producer (i.e. the
      //     one that had the last write permission). Every reader will get
      //     clean cache lines, written through to its DRAM so that even if
      //     they are replaced they can be found upon a miss. We always fetch
      //     from the last (write) producer, because we update the locations
      //     only when a new read/write dep is found.
      noc_dma_add_new(src_core_id,
                      (void *) task->dma_addresses[i],
                      task->run_core_id,
                      (void *) task->dma_addresses[i],
                      size,
                      (options == MM_PACK_OPTION_RW) ? 1 : 0,   // I
                      (options == MM_PACK_OPTION_RW) ? 0 : 1,   // W
                      (options == MM_PACK_OPTION_RW) ? 0 : 1,   // C
                      &task->dmas_waiting, 
                      NULL);
      
      // Take notes for the statistics gathering
      dbg_scount(context, DBG_STATS_IDX_NUM_DMAS, 1);
      dbg_scount(context, DBG_STATS_IDX_DMA_TOTAL_SIZE, size);
    }
  }

  // Mark we've started the DMAs for this task
  task->dmas_started = 1;

  // Success
  return 0;
}
