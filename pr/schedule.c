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
// Abstract      : Scheduling functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: schedule.c,v $
// CVS revision  : $Revision: 1.14 $
// Last modified : $Date: 2013/02/26 11:20:29 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <debug.h>
#include <syscall.h>


// ===========================================================================
// pr_sched_schedule()          Main scheduling decision routine. Given a
//                              list of location and weights of a task's
//                              arguments, it decides on which of our children
//                              should the task be run. If we are an L0
//                              scheduler, this decides which of the workers
//                              will get the task. If we are a higher-level
//                              scheduler, this decides which of our subtrees
//                              must get the new task.
//
//                              Two scores are compiled, one based on the
//                              affinity of (most of) the task arguments with
//                              respect to their last producer, and another
//                              based on the system load balance. A total score
//                              based on both of the previous ones is assigned
//                              to each candidate and the one with the best
//                              score wins. A round-robin policy is applied on
//                              draws.
// ===========================================================================
// * INPUTS
//   unsigned int task_id       ID of the task being scheduled
//   int *locations             Array of task argument last producer core IDs
//   int *weights               Array of task argument weights (one unit per
//                              MM_ALLOC_ALIGN bytes of argument size)
//   int num_locations          Number of task arguments and size of two
//                              arrays above
//
// * OUTPUTS
//   int *ret_msg_id            If higher-level scheduler, the message ID
//                              of the message sent to our child scheduler
//                              is returned here
//   int *ret_core_id           If L0 scheduler, the core ID of the finally
//                              selected worker is returned here
//
// * RETURN VALUE
//   int                        0: scheduling completed locally, result
//                                 stored in *ret_core_id
//                              1: scheduling decision delegated to a child
//                                 scheduler, a request with message ID
//                                 stored in *ret_msg_id has been sent
// ===========================================================================
int pr_sched_schedule(unsigned int task_id, int *locations, int *weights, 
                      int num_locations, int *ret_msg_id, int *ret_core_id) {

  Context       *context;
  int           child;
  int           parent_core_id;
  int           hop;
  int           load_min;
  int           load_max;
  int           load_range;
  int           locality_weights[PR_MAX_SCHEDULER_CHILDREN];
  int           total_weight;
  int           discarded_weight;
  int           locality_score;
  int           load_score;
  int           total_score;
  int           best_score;
  PrMsgReq      *req;
  int           i;
  int           j;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->pr_chld_run_load);
  ar_assert(context->pr_num_children > 0);
  ar_assert((num_locations >= 0) && (num_locations <= PR_REQ_MAX_SIZE));

//kt_printf("%d: Task 0x%08X has %d locations:\r\n", context->pr_core_id, task_id, num_locations);
//for (i = 0; i < num_locations; i++) {
//  kt_printf("   Loc %3d, weight %3d\r\n", locations[i], weights[i]);
//}

  // =========================================================================
  // Pass 1: Scan all children to discover min/max load and initialize the
  //         locality score array
  // =========================================================================

  load_min = INT_MAX;
  load_max = -1;

  // For all children
  for (i = 0; i < context->pr_num_children; i++) {

    // Initialize locality array
    locality_weights[i] = 0;

    // Find min/max running task load
    if (context->pr_chld_run_load[i] < load_min) {
      load_min = context->pr_chld_run_load[i];
    }
    if (context->pr_chld_run_load[i] > load_max) {
      load_max = context->pr_chld_run_load[i];
    }
  }


  // =========================================================================
  // Pass 2: Scan all locations. If they are not in our subtree, discard
  //         their weight. Otherwise, see which child of ours represent
  //         and update the child locality score with this weight (this is
  //         still not normalized).
  // =========================================================================

  total_weight = 0;
  discarded_weight = 0;

  // See what's our parent core ID
  if (context->pr_parent_sched_id != -1) {
    parent_core_id = pr_scheduler_core_id(context->pr_parent_sched_id);
  }
  else {
    parent_core_id = -1;
  }

  // For all locations
  for (i = 0; i < num_locations; i++) {
    
    // Discover next hop core ID for current location
    ar_assert((locations[i] >= 0) && (locations[i] < context->pr_num_cores));




    hop = context->pr_core_route[locations[i]];
    ar_assert((hop >= 0) && (hop < context->pr_num_cores));

    // If parent, discard and count how much weight we've discarded
    if (hop == parent_core_id) {
      discarded_weight += weights[i];
    }
    // If child, update location
    else {
      child = context->pr_core_child_idx[hop];
      ar_assert((child >= 0) && (child < context->pr_num_children));
      locality_weights[child] += weights[i];
    }

    // Count total weight
    total_weight += weights[i];
  }


  // =========================================================================
  // Pass 3: Re-scan all children, this time starting from the last round-
  //         robin position so we can favor a different child in case of 
  //         a draw (reminder: the loads are approximate, so draws will 
  //         happen more often than it seems). Now that we know the load
  //         min/max and the total/discarded weight, we can create two
  //         normalized scores: one for locality and one for load. Make
  //         an on-the-fly decision of the total score and select the best
  //         one.
  // =========================================================================

  child = -1;
  best_score = -1;
  load_range = load_max - load_min;

  // For all children
  for (i = 0; i < context->pr_num_children; i++) {

    // Round-robin offset position
    j = (context->pr_sched_rr + i) % context->pr_num_children;

    // Locality score: We normalize the affinity of the child vs. the
    //                 total weight (and ignore the discarded). This policy
    //                 favors a location less if e.g. a 64-B object has a 
    //                 local affinity, but we're going to fetch a 256-KB
    //                 region anyway outside this subtree.
    //
    //                 The locality score is normalized from 0 to 1024, where
    //                 1024 is the perfect score (all DMAs are fetched from
    //                 this child).
    if (total_weight && locality_weights[j]) {
      locality_score = (locality_weights[j] << 10) / total_weight;
    }
    else {
      locality_score = 0;
    }

    // Load score: We see how much balanced is this child in terms of all
    //             children load.
    //
    //             The load score is also normalized from 0 to 1024, where
    //             1024 is the perfect score (child has the least load)
    //             and 0 is the worst one (child has the most load).
    if (load_range) {
      load_score = ((load_max - context->pr_chld_run_load[j]) << 10) /
                   load_range;
    }
    else {
      load_score = 0;
    }

    // Final policy selector: decide total score based on a mix of the 
    // locality and load scores. These are generally conflicting: the best
    // locality would always be achieved by scheduling everything to the same
    // core, which unfortunately produces the worst load balance.

    total_score = context->pr_load_vs_locality * load_score +
                  (100 - context->pr_load_vs_locality) * locality_score;

//kt_printf("%d: core %d[%d] locality %d, load %d, total %d\r\n", context->pr_core_id, context->pr_children[j], j, locality_score, load_score, total_score);

    // See if this is a winner and keep it
    if (total_score > best_score) {
      best_score = total_score;
      child = j;
    }
  }

  // Save the decision to start from the next core the next time
  ar_assert((child > -1) && (child < context->pr_num_children));
  context->pr_sched_rr = (child + 1) % context->pr_num_children;


  // If we're a leaf scheduler, the scheduling has been concluded. Return
  // result and quit.
  if (!context->pr_scheduler_level) {
    ar_assert(ret_core_id);
    *ret_core_id = context->pr_children[child];
//kt_printf("%d: Scheduled task 0x%08X to core %d, score %d [final]\r\n", context->pr_core_id, task_id, *ret_core_id, best_score);
    return 0;
  }
  

  // Continue the scheduling donwstream
//kt_printf("%d: Scheduled task 0x%08X to core %d, score %d [continuing downstream]\r\n", context->pr_core_id, task_id, context->pr_children[child], best_score);
  req = noc_msg_send_get_buf(context->pr_children[child]);

  req->core_id     = context->pr_core_id;
  req->req_id      = context->pr_message_id;
  req->type        = EXT_REQ_SCHEDULE;
  req->size        = task_id;
  req->num_regions = num_locations;
  for (i = 0; i < num_locations; i++) {
    // Note: the EXT_REQ_SCHEDULE request does not have packing options, so
    // it's only the locations and the weights fields.
    req->data[i] = (void *) ((locations[i] << MM_PACK_SIZE_BITS) | weights[i]);
  }
  
  // Send message to the selected child
  ar_assert(!noc_msg_send());

  // Increase message ID
  ar_assert(ret_msg_id);
  *ret_msg_id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  return 1;
}


// ===========================================================================
// pr_sched_start_scheduling()  Gets the final result from all packing
//                              operations for a new task and starts the
//                              scheduler to decide which worker should
//                              run this. If the scheduler decides locally,
//                              the task is dispatched. Otherwise,
//                              a note-to-self is created to wait for the
//                              scheduling decision and then dispatch the
//                              task.
// ===========================================================================
// * INPUTS
//   PrEvtHookPack *state       Region packing final result
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_sched_start_scheduling(PrEvtHookPack *state) {

  Context       *context;
  PrTaskDescr   *task;
  int           locations[PR_REQ_MAX_SIZE];
  int           weights[PR_REQ_MAX_SIZE];
  int           num_locations;
  int           cur_location;
  int           cur_weight;
  int           this_idx;
  int           min_idx;
  int           min_weight;
  int           msg_id;
  int           core_id;
  int           ret;
  PrEvtPending  *event;
  int           i;
  int           j;


  // Sanity checks. TODO: If packing gives an error other than NOTHING_TO_PACK,
  // this must be a user error. Handle this gracefully instead of crashing.
  context = mm_get_context(ar_get_core_id());
  ar_assert((!state->error_status) || 
            (state->error_status == ERR_NOTHING_TO_PACK));

  // Locate task
  task = NULL;
  kt_trie_find(context->pr_tasks, state->native_task_id, (void *) &task);
  ar_assert(task);
  ar_assert(task->id == state->native_task_id);


//kt_printf("%d: Packing set for task 0x%08X is as follows:\r\n", context->pr_core_id, task->id);
//for (i = 0; i < state->num_elements; i++) {
//kt_printf("   [%d] Adr 0x%08X, size %d, location + options = 0x%08X\r\n", 
//i, 
//state->addresses[i], 
//state->sizes[i] & ((1 << MM_PACK_SIZE_BITS) - 1),
//state->sizes[i] >> MM_PACK_SIZE_BITS);
//}


  // Move state->sizes and state->addresses on the task descriptor. Warning:
  // the caller function MUST NOT free them. They will be freed upon task
  // deletion.
  task->dma_addresses = state->addresses;
  task->dma_sizes     = state->sizes;
  task->num_dmas      = state->num_elements;

  // Find the source core locations that are the last producers of the data
  // needed for the new task. We keep only the top PR_REQ_MAX_SIZE of them, so
  // the scheduling communication can fit into a single message. This should be
  // enough to provide accuracy. If not, change the scheduling communication to
  // be in multi-part messages.
  for (i = 0; i < PR_REQ_MAX_SIZE; i++) {
    // Init all weights to invalid
    weights[i] = 0;
  }
  num_locations = 0;
  for (i = 0; i < task->num_dmas; i++) {
    // Get current weight and location. As a weight, we divide the bytes 
    // by the minimum allocation/communication granularity.
    cur_location = (task->dma_sizes[i] >> MM_PACK_SIZE_BITS) &
                   ((1 << MM_PACK_LOCATION_BITS) - 1);
    cur_weight   = (task->dma_sizes[i] & ((1 << MM_PACK_SIZE_BITS) - 1)) /
                   MM_ALLOC_ALIGN;
    ar_assert(cur_weight);

    // Search if we already have this location, while also scanning for the
    // one with the minimum weight
    min_idx = -1;
    this_idx = -1;
    min_weight = INT_MAX;
    for (j = 0; j < PR_REQ_MAX_SIZE; j++) {
      if ((weights[j]) && (locations[j] == cur_location)) {
        this_idx = j;
        break;
      }
      if (weights[j] < min_weight) {
        min_idx = j;
        min_weight = weights[j];
      }
      if (!weights[j]) {
        num_locations++;
        break;
      }
    }

    // Add to weight of existing slot...
    if (this_idx > -1) {
      weights[this_idx] += cur_weight;
    }
    // ... or replace minimum weight/empty slot
    else if (min_idx > -1) {
      locations[min_idx] = cur_location;
      weights[min_idx]   = cur_weight;
    }
    else {
      ar_abort();
    }
  }
//kt_printf("task_id 0x%08X has %d DMAs:\r\n", task->id, task->num_dmas);
//for (j = 0; j < task->num_dmas; j++) {
//kt_printf("  %d: adr 0x%08X size 0x%08X\r\n", j, task->dma_addresses[j], task->dma_sizes[j]);
//}
//kt_printf("task_id 0x%08X has %d locations:\r\n", task->id, num_locations);
//for (j = 0; j < num_locations; j++) {
//kt_printf("  %d: Location %d has weight %d\r\n", j, locations[j], weights[j]);
//}

  ar_assert((num_locations == PR_REQ_MAX_SIZE) || 
            ((num_locations < PR_REQ_MAX_SIZE) && (!weights[num_locations])));

  // Call the actual scheduling routine
  ret = pr_sched_schedule(task->id, locations, weights, num_locations, 
                          &msg_id, &core_id);

  if (ret == 0) {
    // Scheduling done, store core ID of the worker on the task descriptor
    ar_assert((core_id >= 0) && (core_id < context->pr_num_cores));
    task->run_core_id = core_id;

    // Task is scheduled, proceed with location update and dispatch
    pr_task_scheduled(task);
  }
  else if (ret == 1) {
    // Scheduling continues downstream, create a note-to-self to wait the final
    // result and do location update/dispatch when it comes.
    ar_assert(msg_id);
    event = kt_malloc(sizeof(PrEvtPending));
    event->req = kt_malloc(sizeof(PrMsgReq));

    // Fill in the local task descriptor
    event->req->core_id = -1;
    event->req->req_id  = -1;
    event->req->ptr     = task;

    // Create the rest of the note to self
    event->req->type = SELF_SCHEDULE_RESULT;
    event->action = PR_ACT_REDO;
    event->prev = NULL;
    event->next = NULL;
    event->data = NULL;

    // Store event; we don't expect conflicts on this message ID.
    ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));
  }
  else {
    ar_abort();
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_sched_report_load()       Checks if our region and/or task load has
//                              changed enough. If so, sends an upstream
//                              message to our parent scheduler to indicate our
//                              new load.
// ===========================================================================
// * RETURN VALUE
//   int                        0: Not enough load change, no new request
//                              > 0: Message ID of the new load report request
// ===========================================================================
int pr_sched_report_load() {

  Context               *context;
  PrMsgReq              *new_req;
  int                   id;
  int                   diff_reg;
  int                   diff_run;
  int                   diff_sched;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Top-level schedulers are tough and don't report to nobody
  if (context->pr_parent_sched_id == -1) {
    return 0;
  }

  // We have to report when our current loads differ at least a certain
  // percentage from what we had last reported
  if (context->mm_current_load > context->mm_reported_load) {
    diff_reg = context->mm_current_load - context->mm_reported_load;
  }
  else {
    diff_reg = context->mm_reported_load - context->mm_current_load;
  }

  if (context->pr_cur_run_load > context->pr_rep_run_load) {
    diff_run = context->pr_cur_run_load - context->pr_rep_run_load;
  }
  else {
    diff_run = context->pr_rep_run_load - context->pr_cur_run_load;
  }

  if (context->pr_cur_sched_load > context->pr_rep_sched_load) {
    diff_sched = context->pr_cur_sched_load - context->pr_rep_sched_load;
  }
  else {
    diff_sched = context->pr_rep_sched_load - context->pr_cur_sched_load;
  }

  if ((diff_reg   <= MM_LOAD_REPORT_CHANGE * context->mm_reported_load) &&
      (diff_run   <= PR_LOAD_REPORT_CHANGE * context->pr_rep_run_load)  &&
      (diff_sched <= PR_LOAD_REPORT_CHANGE * context->pr_rep_sched_load)) {
    return 0;
  }

  // Remember we're reporting the current load values
  context->mm_reported_load  = context->mm_current_load;
  context->pr_rep_run_load   = context->pr_cur_run_load;
  context->pr_rep_sched_load = context->pr_cur_sched_load;

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_LOAD_REPORT;
  new_req->region  = context->mm_current_load;
  new_req->ptr     = (void *) ((size_t) context->pr_cur_run_load);
  new_req->size    = context->pr_cur_sched_load;
  
  // Send message to parent scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// pr_sched_update_location()   Updates the location of local regions/objects
//                              out of an argument list. For non-local 
//                              regions/objects, it groups them per child
//                              scheduler and forwards the list parts to them.
// ===========================================================================
// * INPUTS
//   int location               Core ID of the location of the regions/objects
//                              to set
//   int num_agrs               Number of dependency arguments to handle
//   void **args                Dependency argument array
//   unsigned char *dep_flags   Dependency flags array (only region information
//                              is required when in request mode; full flags
//                              are required when in task mode)
//   int req_mode               0: request mode (do all updates)
//                              1: task mode (update only the writable deps)
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_sched_update_location(int location, int num_args, void **args, 
                             unsigned char *dep_flags, int req_mode) {

  typedef struct {
    void          **deps;
    unsigned char *flags;
    int           num_deps;
  } per_sched_type;

  Context               *context;
  Trie                  *trie;
  per_sched_type        *per_sched;
  per_sched_type        local_sched;
  PrMsgReq              *req;
  int                   sched_id;
  void                  **queue_args;
  unsigned char         *queue_flags;
  void                  *dep;
  unsigned char         flag;
  int                   head;
  int                   tail;
  int                   cur;
  rid_t                 rid;
  MmRgnTreeNode         *reg_node;
  int                   i;
  int                   j;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert((location >= 0) && (location < context->pr_num_cores));


  // =========================================================================
  // If we have children schedulers, separate arguments into those handled by
  // local scheduler and others handled by our children. Also, examine all
  // regions: if they have children, recurse into them to update their 
  // locations as well.
  // =========================================================================

  // Shortcut for leaf schedulers: everything is local
  if (context->pr_scheduler_level == 0) {
    trie = NULL;
  }

  else {
    
    // Create a trie to store per-scheduler packing needs. We'll be using IDs
    // from 1 to context->pr_num_schedulers, so the MSB is the log2 of
    // context->pr_num_schedulers. The kt_int_log2() function will return the
    // MSB position correctly, even for non-power-of-2 values.
    ar_assert(trie = kt_alloc_trie(kt_int_log2(context->pr_num_schedulers), 0));
  }

  // Nothing so far is local
  local_sched.deps = NULL;
  local_sched.flags = NULL;
  local_sched.num_deps = 0;

  // Create a queue for recursion and initialize it with all arguments
  queue_args  = kt_malloc(num_args * sizeof(void *));
  queue_flags = kt_malloc(num_args * sizeof(unsigned char));

  for (i = 0, j = 0; i < num_args; i++) {

    // Ignore by-value arguments (i.e. EXACTLY objects that are safe and
    // read-only; we don't handle safe regions or other variations yet).
    if (dep_flags[i] == SYS_TYPE_BYVALUE_ARG) {
      continue;
    }

    // In task mode, also ignore arguments which are read-only. This implements
    // the policy of letting the location be the one of the last producer of
    // the data (the last task that had write access).
    if (!req_mode && !(dep_flags[i] & SYS_TYPE_OUT_ARG)) {
      continue;
    }

    queue_args[j] = args[i];
    queue_flags[j] = dep_flags[i];
    j++;
  }

  head = j;
  tail = 0;


  // Iterate on the queue
  while (head > tail) {

    // Get next region/object
    dep = queue_args[tail];
    flag = queue_flags[tail];
    tail++;

    // Find which scheduler should we ask for this region or object
    if (flag & SYS_TYPE_REGION_ARG) {
      sched_id = mm_distr_region_sched_id((rid_t) dep);
    }
    else {
      sched_id = mm_distr_object_sched_id(dep);
    }
    ar_assert(sched_id >= 0);
    ar_assert(sched_id != context->pr_parent_sched_id);

    // Is it ours?
    if (sched_id == context->pr_scheduler_id) {

      // Add it to local pile
      local_sched.deps = kt_realloc(local_sched.deps, 
                                    (local_sched.num_deps + 1) * 
                                    sizeof(void *));
      local_sched.flags = kt_realloc(local_sched.flags, 
                                    (local_sched.num_deps + 1) * 
                                    sizeof(unsigned char));
      local_sched.deps[local_sched.num_deps] = dep;
      local_sched.flags[local_sched.num_deps] = flag;
      local_sched.num_deps++;

      // If it's a region, add its children for recursive inspection
      if (flag & SYS_TYPE_REGION_ARG) {

        // Find local region tree node
        ar_assert(kt_trie_find(context->mm_local_rids, (size_t) dep, 
                               (void *) &reg_node));

        if (reg_node->num_children) {

          // Enlarge queue
          queue_args  = kt_realloc(queue_args, 
                                   (head + reg_node->num_children) * 
                                   sizeof(void *));
          queue_flags = kt_realloc(queue_flags, 
                                   (head + reg_node->num_children) * 
                                   sizeof(unsigned char));

          // Add children, marking they're regions
          for (i = 0; i < reg_node->num_children; i++) {
            queue_args[head + i]  = (void *) reg_node->children_ids[i];
            queue_flags[head + i] = SYS_TYPE_REGION_ARG;
          }

          head += reg_node->num_children;
        }
      }
    }
    
    // It belongs to one of our children
    else {

      ar_assert(trie);

      // We use scheduler ID + 1, because tries cannot handle zero keys
      sched_id++;

      // Do we have anything else for this scheduler ID?
      kt_trie_find(trie, sched_id, (void *) &per_sched);

      // Append this region to existing entry
      if (per_sched) {
        per_sched->deps = kt_realloc(per_sched->deps, 
                                     (per_sched->num_deps + 1) * 
                                     sizeof(void *));
        per_sched->flags = kt_realloc(per_sched->flags, 
                                      (per_sched->num_deps + 1) * 
                                      sizeof(unsigned char));
        per_sched->deps[per_sched->num_deps] = dep;
        per_sched->flags[per_sched->num_deps] = flag;
        per_sched->num_deps++;
      }

      // Create new entry for this scheduler
      else {
        per_sched = kt_malloc(sizeof(per_sched_type));
        per_sched->deps = kt_malloc(sizeof(void *));
        per_sched->flags = kt_malloc(sizeof(unsigned char));
        per_sched->deps[0] = dep;
        per_sched->flags[0] = flag;
        per_sched->num_deps = 1;

        // Insert it into the trie
        ar_assert(!kt_trie_insert(trie, sched_id, per_sched));
      }
    }
  }


  // =========================================================================
  // For each child scheduler we need to communicate with, send to it a 
  // message to update the location of these arguments (or send them 
  // further downstream).
  // =========================================================================

  if (trie) {

    // For all schedulers that we need to communicate with
    for (sched_id = kt_trie_find_minmax(trie, 0, (void *) &per_sched);
         sched_id;
         sched_id = kt_trie_find_next(trie, 1, (void *) &per_sched)) {

      // Restore correct scheduler ID value
      sched_id--;

      // Initialize
      cur = 0;
      req = NULL;

      // Loop until all regions and objects have been sent
      while (cur < per_sched->num_deps) {

        // Start building new message
        req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

        req->type        = EXT_REQ_UPDATE_LOCATION;
        req->core_id     = context->pr_core_id;
        req->req_id      = context->pr_message_id;
        req->num_regions = location;

        // Populate message array and bitmap (only region info here)
        for (req->num_ptrs = 0, req->region = 0;
             (req->num_ptrs < PR_REQ_MAX_SIZE) && (cur < per_sched->num_deps);
             req->num_ptrs++, cur++) {

          req->data[req->num_ptrs] = per_sched->deps[cur];

          req->region |= (per_sched->flags[cur] & SYS_TYPE_REGION_ARG) 
                       ? (1 << (2 * req->num_ptrs + 1)) : 0;
        }

        // Send request to child
        ar_assert(!noc_msg_send());
        req = NULL;
        
        // Increase message ID, avoiding value 0 on wrap-arounds
        context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
      }

      // Free this entry
      kt_free(per_sched->deps);
      kt_free(per_sched->flags);
      kt_free(per_sched);
    }

    // Free the trie
    kt_free_trie(trie, NULL);
  }


  // =========================================================================
  // Finally, for all arguments handled locally, update the location
  // =========================================================================
  for (i = 0; i < local_sched.num_deps; i++) {

    // Region?
    if (local_sched.flags[i] & SYS_TYPE_REGION_ARG) {

//kt_printf("%d: Updating region %d location to %d\r\n", context->pr_core_id, local_sched.deps[i], location);

      // Locate the region node
      ar_assert(kt_trie_find(context->mm_local_rids, 
                             (size_t) local_sched.deps[i], (void *) &reg_node));

      // Update the location
      reg_node->location = location;

      // Are there objects that have exceptions? Remove all of them: the 
      // region has a new location, so all of the objects are taken with
      // it and there are no exceptions.
      kt_reset_trie(reg_node->obj_locations, NULL);
    }

    // Object?
    else {

//kt_printf("%d: Updating object 0x%08X location to %d\r\n", context->pr_core_id, local_sched.deps[i], location);

      // Locate the region node where the object resides
      ar_assert(!mm_region_query_pointer((void *) local_sched.deps[i], NULL,
                                         (void *) &rid));
      ar_assert(kt_trie_find(context->mm_local_rids, rid, (void *) &reg_node));


      // Try to delete (i.e. make sure there isn't) any exception for the
      // object
      kt_trie_delete(reg_node->obj_locations, (size_t) local_sched.deps[i], 
                     NULL);

      // If the current region location is different from the new object
      // location, create an exception for the object
      if (reg_node->location != location) {
        ar_assert(!kt_trie_insert(reg_node->obj_locations, 
                                  (size_t) local_sched.deps[i],
                                  (void *) location));
      }
    }
  }


  // Free stuff
  kt_free(local_sched.deps);
  kt_free(local_sched.flags);
  kt_free(queue_args);
  kt_free(queue_flags);


  // Success
  return 0;
}

