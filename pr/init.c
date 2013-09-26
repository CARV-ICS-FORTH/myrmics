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
// Abstract      : Processing layer initialization and utility functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: init.c,v $
// CVS revision  : $Revision: 1.36 $
// Last modified : $Date: 2013/03/22 12:25:11 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <stdarg.h>

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <syscall.h>


// ===========================================================================
// pr_init_internal()           Processing layer initialization. Assigns
//                              core roles (schedulers/workers) and their
//                              relationships.
//
//                              This is the real function, called by
//                              pr_init().
// ===========================================================================
// * INPUTS
//   int num_cores              Number of total active cores in the system.
//                              See pr_init() below for description.
//
//   int *list                  The variable-argument list in an array form.
//                              Contains 4 arguments per core, see pr_init()
//                              below for description.
//
// * RETURN VALUE
//   int                        0: success, and this core is part of the
//                                 running setup
//                              1: this core is not part of the setup and
//                                 should not continue running
//                              2: for video demo configurations, this core
//                                 should execute the video input loop
//                              3: for video demo configurations, this core
//                                 should execute the video output loop
// ===========================================================================
int pr_init_internal(int num_cores, int *list) {

  Context       *context;
  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  int           pid;
  PrRole        role;
  int           num_workers;
  int           num_schedulers;
  int           sched_level;
  int           target_id;
  int           *parents;
  int           cur_parent;
  int           cur_child;
  int           found;
  int           i;


  // Get context and architectural board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);

  // Size sanity checks. Adjust PR_REQ_MAX_SIZE and PR_RPL_MAX_SIZE to make
  // each request and each reply occupy exactly one NOC message.
  ar_assert(sizeof(PrMsgReq) == NOC_MESSAGE_SIZE);
  ar_assert(sizeof(PrMsgReply) == NOC_MESSAGE_SIZE);

  // Set all related context fields to invalid status
  context->pr_core_id = -1;
  context->pr_role = -1;
  context->pr_core_bid_cid = NULL;
  context->pr_core_work_ids = NULL;
  context->pr_core_sched_ids = NULL;
  context->pr_core_child_idx = NULL;
  context->pr_work_core_ids = NULL;
  context->pr_sched_core_ids = NULL;
  context->pr_num_schedulers = -1;
  context->pr_num_workers = -1;
  context->pr_worker_id = -1;
  context->pr_scheduler_id = -1;
  context->pr_parent_sched_id = -1;
  context->pr_scheduler_level = -1;
  context->pr_children = NULL;
  context->pr_num_children = 0;
  context->pr_chld_sched_load = NULL;
  context->pr_chld_run_load = NULL;
  context->vid_demo_in_bid = -1;
  context->vid_demo_in_cid = -1;
  context->vid_demo_out_bid = -1;
  context->vid_demo_out_cid = -1;

  // First, handle the video demo special cases. If the ROLE_VID_* roles
  // are specified, they must be at the two last core positions.
  if ((list[(num_cores - 1) * 4 + 2] == ROLE_VID_OUTPUT) ||
      (list[(num_cores - 1) * 4 + 2] == ROLE_VID_INPUT)) {
    ar_assert(num_cores >= 4); // 2 video cores, 1 scheduler, 1 worker minimum

    for (i = num_cores - 2; i < num_cores; i++) {
      bid  = list[i * 4];
      cid  = list[i * 4 + 1];
      role = list[i * 4 + 2];
      if (role == ROLE_VID_INPUT) {
        context->vid_demo_in_bid = bid;
        context->vid_demo_in_cid = cid;
        if ((my_bid == bid) && (my_cid == cid)) {
          context->pr_role = role;
        }
      }
      else if (role == ROLE_VID_OUTPUT) {
        context->vid_demo_out_bid = bid;
        context->vid_demo_out_cid = cid;
        if ((my_bid == bid) && (my_cid == cid)) {
          context->pr_role = role;
        }
      }
      else {
        // We expect the last 2 positions to be ROLE_VID_*
        ar_abort();
      }
    }

    // Reduce num_cores, so all the normal Myrmics cores will see a 
    // number being the sum of schedulers and workers. The video cores
    // are special and their arch-level board/core IDs can be found in 
    // the context.
    num_cores -= 2;
  }

  // Video cores return here
  if (context->pr_role == ROLE_VID_INPUT) {
    return 2;
  }
  else if (context->pr_role == ROLE_VID_OUTPUT) {
    return 3;
  }


  // Remember number of Myrmics cores
  context->pr_num_cores = num_cores;


  // Single-core hack
  if (num_cores == 1) {
    
    // Just get arch board/core ID
    bid = list[0];
    cid = list[1];

    // Set single-core mode for the core that matches
    if ((my_bid == bid) && (my_cid == cid)) {

      context->pr_core_id = 0;
      context->pr_role = ROLE_SCHEDULER;
      context->pr_num_schedulers = 1;
      context->pr_num_workers = 1;
      context->pr_worker_id = 0;
      context->pr_scheduler_id = 0;
      context->pr_parent_sched_id = -1;
      context->pr_scheduler_level = 0;
      context->pr_children = kt_malloc(sizeof(int));
      context->pr_children[0] = 0;
      context->pr_num_children = 1;
      context->pr_chld_sched_load = kt_malloc(sizeof(int));
      context->pr_chld_sched_load[0] = 0;
      context->pr_chld_run_load = kt_malloc(sizeof(int));
      context->pr_chld_run_load[0] = 0;

      return 0;
    }
    else {
      // Fail cores that don't match
      return 1;
    }
  }


  // Initialize variables
  num_schedulers = 0;
  num_workers = 0;
  role = -1;

  // Initialize arrays for the first pass
  context->pr_core_bid_cid = kt_malloc(num_cores * sizeof(int));
  context->pr_core_work_ids = kt_malloc(num_cores * sizeof(int));
  context->pr_core_sched_ids = kt_malloc(num_cores * sizeof(int));

  for (i = 0; i < num_cores; i++) {
    context->pr_core_bid_cid[i] = -1;
    context->pr_core_work_ids[i] = -1;
    context->pr_core_sched_ids[i] = -1;
  }



  // First pass of the list, to assign core IDs to worker & scheduler IDs
  for (i = 0; i < num_cores; i++) {

    // Get four arguments for this core ID
    bid  = list[i * 4];
    cid  = list[i * 4 + 1];
    role = list[i * 4 + 2];
    pid  = list[i * 4 + 3];

    // Are we this core ID?
    if ((bid == my_bid) && (cid == my_cid)) {
      context->pr_core_id = i;
      context->pr_role = role;
      if (role == ROLE_WORKER) {
        context->pr_worker_id = num_workers;
      }
      else {
        context->pr_scheduler_id = num_schedulers;
      }
    }

    // Fill out arrays for this core
    context->pr_core_bid_cid[i] = (bid << 8) | cid;
    if (role == ROLE_WORKER) {
      context->pr_core_work_ids[i] = num_workers++;
    }
    else if (role == ROLE_SCHEDULER) {
      context->pr_core_sched_ids[i] = num_schedulers++;
    }
    else {
      ar_abort();
    }
  }

  // Copy number of workers/schedulers
  context->pr_num_workers = num_workers;
  context->pr_num_schedulers = num_schedulers;


  // We must have found ourselves in the core list, otherwise just quit
  if (context->pr_core_id == -1) {
    kt_free(context->pr_core_bid_cid);
    kt_free(context->pr_core_work_ids);
    kt_free(context->pr_core_sched_ids);
    context->pr_core_bid_cid = NULL;
    context->pr_core_work_ids = NULL;
    context->pr_core_sched_ids = NULL;
    return 1;
  }


  // Initialize inverse arrays for the second pass
  context->pr_work_core_ids = kt_malloc(num_workers * sizeof(int));
  context->pr_sched_core_ids = kt_malloc(num_schedulers * sizeof(int));

  for (i = 0; i < num_workers; i++) {
    context->pr_work_core_ids[i] = -1;
  }
  for (i = 0; i < num_schedulers; i++) {
    context->pr_sched_core_ids[i] = -1;
  }

  // Allocate temporary parents array for all cores
  parents = kt_malloc(num_cores * sizeof(int));
  for (i = 0; i < num_cores; i++) {
    parents[i] = -1;
  }



  // Restart counters
  num_workers = 0;
  num_schedulers = 0;

  // Second pass of the list, to assign inverse arrays and discover our
  // father and children
  for (i = 0; i < num_cores; i++) {

    // Get four arguments for this core ID
    bid  = list[i * 4];
    cid  = list[i * 4 + 1];
    role = list[i * 4 + 2];
    pid  = list[i * 4 + 3];

    // Store parent core ID in the temporary array for later processing
    parents[i] = pid;

    // Are we this core ID? Assign parent scheduler ID, now we know from 
    // the first pass all the scheduler IDs.
    if ((bid == my_bid) && (cid == my_cid)) {
      if (pid == -1) {
        context->pr_parent_sched_id = -1; // top-level scheduler
      }
      else {
        ar_assert((pid >= 0) && (pid < num_cores));
        context->pr_parent_sched_id = context->pr_core_sched_ids[pid];
        ar_assert(context->pr_parent_sched_id < context->pr_num_schedulers);
      }
    }

    // Are we a father to this core? Add its core ID to our children array.
    if (pid == context->pr_core_id) {
      ar_assert(context->pr_role == ROLE_SCHEDULER);
      context->pr_children = kt_realloc(context->pr_children,
                                (context->pr_num_children + 1) * sizeof(int));
      context->pr_children[context->pr_num_children++] = i;
    }


    // Fill out inverse arrays for this core
    if (role == ROLE_WORKER) {
      context->pr_work_core_ids[num_workers++] = i;
    }
    else if (role == ROLE_SCHEDULER) {
      context->pr_sched_core_ids[num_schedulers++] = i;
    }
  }


  // We now analyze the parent-children relationships from the temporary array,
  // so we can populate the pr_core_route table.
  context->pr_core_route = kt_malloc(num_cores * sizeof(int));
  for (i = 0; i < num_cores; i++) {

    // Cannot route messages to ourselves
    if (i == context->pr_core_id) {
      context->pr_core_route[i] = -1;
      continue;
    }

    // Assume core i is a (possibly distant) child. Go up until we encounter
    // ourselves. The next to last hop (our direct child) is the core we should
    // forward messages going to i.
    found = 0;
    for (cur_child = i, cur_parent = parents[i]; 
         cur_parent > -1; 
         cur_child = cur_parent, cur_parent = parents[cur_parent]) {
      ar_assert(cur_parent < num_cores);
      if (cur_parent == context->pr_core_id) {
        context->pr_core_route[i] = cur_child;
        found = 1;
        break;
      }
    }

    // If we didn't find it, we assume is either above us in the tree or it's
    // on a different subtree of our parent. In any case, the next hop towards
    // i should be our parent.
    if (!found) {
      ar_assert(context->pr_parent_sched_id > -1);
      context->pr_core_route[i] = 
                    context->pr_sched_core_ids[context->pr_parent_sched_id];
    }
  }


  // Create inverse mapping of core IDs to children index ID, only for
  // schedulers
  if (context->pr_role == ROLE_SCHEDULER) {

    // Allocate child indexing array
    context->pr_core_child_idx = kt_malloc(num_cores * sizeof(int));

    // By the way, create children load arrays
    ar_assert(context->pr_chld_sched_load = kt_malloc(context->pr_num_children *
                                                    sizeof(int)));
    ar_assert(context->pr_chld_run_load = kt_malloc(context->pr_num_children *
                                                    sizeof(int)));
    for (i = 0; i < context->pr_num_children; i++) {
      context->pr_chld_sched_load[i] = 0;
      context->pr_chld_run_load[i] = 0;
    }

    // Initialize everyone as a non-child
    for (i = 0; i < num_cores; i++) {
      context->pr_core_child_idx[i] = -1;
    }

    // Store children indices
    for (i = 0; i < context->pr_num_children; i++) {
      ar_assert((context->pr_children[i] >= 0) && 
                (context->pr_children[i] < num_cores));
      context->pr_core_child_idx[context->pr_children[i]] = i;
    }
  }


  // Finally, if we are a scheduler, we need to know our scheduler level. We do
  // this the easy and slow way, which is do multiple passes on the vararg list
  // to discover that. A faster way would be to have already built a tree from
  // the 1st pass and traverse the tree. Since we won't be needing such a tree
  // for normal operation, that's an overkill for a non-critical initialization
  // routine, so we choose the dumb way to do things.
  if (context->pr_role == ROLE_SCHEDULER) {

    // Start as a lowest level scheduler
    sched_level = 0;

    // Start by searching for ourselves
    target_id = context->pr_core_id;

    // Find first child that has the target as its parent
    while (1) {
      for (i = 0; i < num_cores; i++) {
        bid  = list[i * 4];
        cid  = list[i * 4 + 1];
        role = list[i * 4 + 2];
        pid  = list[i * 4 + 3];

        if (pid == target_id) {
          break;
        }
      }

      // Sanity check: target must be referenced by someone
      ar_assert(i < num_cores);

      // Stop when this child is a worker
      if (role == ROLE_WORKER) {
        break;
      }

      // Increment level and continue search for our son scheduler
      sched_level++;
      target_id = i;
    }

    context->pr_scheduler_level = sched_level;
  }


  // Free temporary stuff
  kt_free(parents);


  // Sanity checks
  ar_assert(context->pr_num_workers);
  ar_assert(context->pr_num_schedulers);
  ar_assert(num_workers == context->pr_num_workers);
  ar_assert(num_schedulers == context->pr_num_schedulers);
  if (context->pr_role == ROLE_SCHEDULER) {
    ar_assert((context->pr_scheduler_id >= 0) && 
              (context->pr_scheduler_id < context->pr_num_schedulers));
    ar_assert((context->pr_parent_sched_id == -1) ||  // top-level sched
              ((context->pr_parent_sched_id >= 0) && // mid/low-level sched
               (context->pr_parent_sched_id < context->pr_num_schedulers)));
    ar_assert(context->pr_scheduler_level >= 0);
    ar_assert(context->pr_children);
    ar_assert(context->pr_num_children);
    ar_assert(context->pr_num_children <= PR_MAX_SCHEDULER_CHILDREN);
  }
  else if (context->pr_role == ROLE_WORKER) {
    ar_assert((context->pr_worker_id >= 0) && 
              (context->pr_worker_id < context->pr_num_workers));
    ar_assert((context->pr_parent_sched_id >= 0) && 
              (context->pr_parent_sched_id < context->pr_num_schedulers));
    ar_assert(context->pr_scheduler_level == -1);
  }
  else {
    ar_abort();
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_init()                    Processing layer initialization. Assigns
//                              core roles (schedulers/workers) and their
//                              relationships.
// ===========================================================================
// * INPUTS
//   int num_cores              Number of total active cores in the system.
//                              This many unique core IDs will be created. 
//                              Four arguments per core must be supplied to
//                              pr_init(), as specified below. The unique
//                              core IDs will be created in order, i.e. the
//                              first four arguments will refer to unique
//                              core ID = 0, the next four to core ID = 1
//                              and the last four to core ID = (num_cores - 1).
//
//                              An exception is if num_cores = 1: there is
//                              no need to define the last two arguments in
//                              this case (role and parent).
//
//   For each active core, specify four arguments:
//
//   int bid                    Architecture-level board ID of core
//   int cid                    Architecture-level core ID of core
//   PrRole role                Core role: ROLE_SCHEDULER or ROLE_WORKER.
//                              Specifically for video demo configurations,
//                              ROLE_VID_INPUT and ROLE_VID_OUTPUT can also
//                              be used for the last 2 core positions.
//   int parent_core_id         Unique core ID of the scheduler responsible for
//                              this core, or -1 to specify that this is the
//                              top-level scheduler
//
// * RETURN VALUE
//   int                        0: success, and this core is part of the
//                                 running setup
//                              1: this core is not part of the setup and
//                                 should not continue running
//                              2: for video demo configurations, this core
//                                 should execute the video input loop
//                              3: for video demo configurations, this core
//                                 should execute the video output loop
// ===========================================================================
int pr_init(int num_cores, ...) {

  va_list       ap;
  int           *list;
  int           ret;
  int           i;


  // Allocate a buffer
  list = kt_malloc(4 * num_cores * sizeof(int));

  // Initialize variable argument list
  va_start(ap, num_cores);

  // Fill the buffer
  for (i = 0; i < 4 * num_cores; i++) {
    list[i] = va_arg(ap, int);
  }

  // Close list
  va_end(ap);

  // Call the real function
  ret = pr_init_internal(num_cores, list);

  // Free the buffer
  kt_free(list);

  // Return
  return ret;

}


// ===========================================================================
// pr_auto_init()               Auto-generates a MB-based tree setup and
//                              calls the real initialization function.
//
//                              For example, for num_levels = 3 and
//                              fanouts[0] = 4 and fanouts[1] = 3 we'd get
//                              a tree with core IDs as follows:
//
//
// L1 sch, 3 children                               15          
// 
//
// 3 x L0 sch, 4 children each       12             13             14
//                                        
//                                         
// 3 x 4 workers                 0  1  2  3     4  5  6  7     8  9 10 11
// ===========================================================================
// * INPUTS
//   int num_levels             The number of tree levels, including the
//                              workers
//   int *fanouts               The fanouts of the tree in an array. The
//                              array should have a size of (num_levels - 1),
//                              with a value for each scheduler fanout.
//                              Specifically, fanouts[0] is the fanout of the
//                              L0 schedulers, fanouts[1] of the L1 schedulers
//                              and so on.
//
// * RETURN VALUE
//   int                        pr_init() return value
// ===========================================================================
int pr_auto_init(int num_levels, int *fanouts) {

  int   max_cores;
  int   num_cores;
  int   *list;
  int   core_id;
  int   x;
  int   y;
  int   z;
  int   c;
  int   level;
  int   *level_cores;
  int   level_rem_cores;
  int   next_level_starts_at;
  int   cur_parent;
  int   cur_parent_children;
  int   ret;
  int   i;


  // Compute per-level and total cores and do sanity checks
  max_cores = (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1) *
              (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
              (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1) *
              AR_FORMIC_CORES_PER_BOARD;
  ar_assert(num_levels >= 2);

  level_cores = kt_malloc(num_levels * sizeof(int));

  level_cores[num_levels - 1] = 1;
  num_cores = 1;
  for (i = num_levels - 2; i >= 0; i--) {
    level_cores[i] = level_cores[i + 1] * fanouts[i];
    num_cores += level_cores[i];
  }
  ar_assert(num_cores <= max_cores);


  // Allocate a buffer
  list = kt_malloc(4 * num_cores * sizeof(int));


  // We do this bottom-up. We start at level 0 and fill out all worker cores.
  // We then go to the upper levels and fill out the schedulers.


  // Initialize counters
  x = AR_FORMIC_MIN_X;
  y = AR_FORMIC_MIN_Y;
  z = AR_FORMIC_MIN_Z;
  c = 0;
  
  level = 0;
  level_rem_cores = level_cores[level];
  next_level_starts_at = level_cores[level];

  cur_parent = next_level_starts_at;
  cur_parent_children = 1;


  // For all core IDs
  for (core_id = 0; core_id < num_cores; core_id++) {

    // Board ID
    list[core_id * 4] = (x << 4) | (y << 2) | z;

    // Core ID
    list[core_id * 4 + 1] = c;

    // Role
    list[core_id * 4 + 2] = (level == 0) ? ROLE_WORKER : ROLE_SCHEDULER;

    // Parent core ID
    list[core_id * 4 + 3] = (cur_parent < num_cores) ? cur_parent : -1;


    // If we're done, break now, otherwise assertions below may crash
    if (core_id == num_cores - 1) {
      break;
    }

    // Prepare next parent
    cur_parent_children++;
    if (cur_parent_children > fanouts[level]) {
      cur_parent_children = 1;
      cur_parent++;
    }

    // Prepare next arch-level board/core IDs
    c++;
    if (c >= AR_FORMIC_CORES_PER_BOARD) {
      // Reset and go to next X board
      c = 0;
      x++;
      if (x > AR_FORMIC_MAX_X) {
        // Reset and go to next Y level
        x = AR_FORMIC_MIN_X;
        y++;
        if (y > AR_FORMIC_MAX_Y) {
          // Reset and go to next Z level
          y = AR_FORMIC_MIN_Y;
          z++;
          if (z > AR_FORMIC_MAX_Z) {
            ar_abort();
          }
        }
      }
    }

    
    // See if we changed level
    level_rem_cores--;
    if (level_rem_cores <= 0) {
      level++;
      next_level_starts_at += level_cores[level];
      level_rem_cores = level_cores[level];
      cur_parent = next_level_starts_at;
      ar_assert(cur_parent_children == 1);
    }
  }


  // Call the real init function
  ret = pr_init_internal(num_cores, list);


  // Free stuff
  kt_free(list);
  kt_free(level_cores);

  // Return
  return ret;
}


// ===========================================================================
// pr_init_app()                Selects which application will run by giving
//                              its task table to all cores. The top-level
//                              scheduler dispatches the first task (index 0)
//                              in the table, using the variadic arguments
//                              given to pr_init_app().
// ===========================================================================
// * INPUTS
//   func_t **task_table        Task table of application to be run
//   void **args                Arguments to be given to task index 0 of the
//                              task table
//   int num_args               Number of arguments
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_init_app(func_t **task_table, void **args, int num_args) {

  Context       *context;
  unsigned char *deps;
  PrTaskDescr   *task;
  int           i;


  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(task_table);

  // Assign task table to everybody
  context->pr_task_table = task_table;

  // Do other task- and dep-related initializations here
  if (context->pr_role == ROLE_SCHEDULER) {

    // Initialize tasks trie
    context->pr_tasks = kt_alloc_trie(PR_TASK_ID_SIZE - 1 , 0);
  }
  else if (context->pr_role == ROLE_WORKER) {

    // Initialize ready queue
    context->pr_ready_queue = kt_alloc_list();
  }
  else {
    ar_abort();
  }

  // Everybody except for the top-level scheduler get out
  if ((context->pr_role != ROLE_SCHEDULER) ||
      (context->pr_parent_sched_id != -1)) {
    return 0;
  }

  // Build argument array
  deps = kt_malloc(num_args * sizeof(unsigned char));
  for (i = 0; i < num_args; i++) {
    deps[i] = SYS_TYPE_BYVALUE_ARG; // no dependencies for main function
  }

  // Create and dispatch main function of application (no dependencies, set
  // to run at the first worker core without scheduling)
  ar_assert(task = pr_task_create(0, 0, args, deps, num_args));
  task->run_core_id = pr_worker_core_id(0);
  ar_assert(!pr_task_dispatch(task));

  // Enqueue main function at the head of the top-level (NULL) region, which
  // has a known region ID of 1
  ar_assert(context->mm_region_tree->id == 1);
  ar_assert(pr_dep_enqueue(context->mm_region_tree->dep_queue, (void *) 1, 
                           SYS_TYPE_INOUT_ARG | SYS_TYPE_REGION_ARG,
                           task->id, 0, -1, 0, 0, NULL, 0, NULL) == 1);


  // Success
  return 0;
}


// ===========================================================================
// pr_get_role()                Find out what our role is as a core
// ===========================================================================
// * RETURN VALUE
//   PrRole                     ROLE_SCHEDULER or ROLE_WORKER
// ===========================================================================
PrRole pr_get_role() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_role;
}


// ===========================================================================
// pr_get_core_id()             Find out what is our unique core ID in the
//                              system. Not to be confused with the
//                              architecture-specific core ID which identifies
//                              a core per board and is returned by
//                              ar_get_core_id().
// ===========================================================================
// * RETURN VALUE
//   int                        Our unique core ID
// ===========================================================================
int pr_get_core_id() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_core_id;
}


// ===========================================================================
// pr_get_worker_id()           Find out what is our worker ID
// ===========================================================================
// * RETURN VALUE
//   int                        Our worker ID
// ===========================================================================
int pr_get_worker_id() {
  Context *context;

  context = mm_get_context(ar_get_core_id());

  if (context->pr_num_cores > 1) {
    ar_assert(context->pr_role == ROLE_WORKER);
    return context->pr_worker_id;
  }
  else {
    ar_assert(context->pr_role == ROLE_SCHEDULER);
    ar_assert(!context->pr_core_id);
    return 0;
  }
}


// ===========================================================================
// pr_get_scheduler_id()        If scheduler: find out which scheduler
//                              instance we are (0 ... num_schedulers - 1).
//                              If worker: find scheduler ID of the scheduler
//                              responsible for us
// ===========================================================================
// * RETURN VALUE
//   int                        scheduler ID, as defined above
// ===========================================================================
int pr_get_scheduler_id() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_scheduler_id;
}


// ===========================================================================
// pr_get_parent_scheduler_id() Find out the scheduler ID we report to.
//                              Returns -1 if we are the top-level scheduler.
// ===========================================================================
// * RETURN VALUE
//   int                        Our parent scheduler ID
// ===========================================================================
int pr_get_parent_scheduler_id() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_parent_sched_id;
}


// ===========================================================================
// pr_get_num_cores()           Find out what is the total number of cores
// ===========================================================================
// * RETURN VALUE
//   int                        Number of cores
// ===========================================================================
int pr_get_num_cores() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_num_cores;
}


// ===========================================================================
// pr_get_num_schedulers()      Find how many schedulers are in the system
// ===========================================================================
// * RETURN VALUE
//   int                        Number of schedulers
// ===========================================================================
int pr_get_num_schedulers() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_num_schedulers;
}


// ===========================================================================
// pr_get_num_workers()         Find how many workers are in the system
// ===========================================================================
// * RETURN VALUE
//   int                        Number of workers
// ===========================================================================
int pr_get_num_workers() {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  return context->pr_num_workers;
}


// ===========================================================================
// pr_worker_core_id()          Translate a worker ID to a unique core ID
// ===========================================================================
// * INPUTS
//   int worker_id              Worker ID that needs translation
//
// * RETURN VALUE
//   int                        Core ID of worker
// ===========================================================================
int pr_worker_core_id(int worker_id) {
  Context *context;

  context = mm_get_context(ar_get_core_id());

  if (context->pr_num_cores > 1) {
    ar_assert((worker_id >= 0) && (worker_id < context->pr_num_workers));
    return context->pr_work_core_ids[worker_id];
  }
  else {
    ar_assert(!worker_id);
    return 0;
  }
}


// ===========================================================================
// pr_scheduler_core_id()       Translate a scheduler ID to a unique core ID
// ===========================================================================
// * INPUTS
//   int sched_id               Scheduler ID that needs translation
//
// * RETURN VALUE
//   int                        Core ID of scheduler
// ===========================================================================
int pr_scheduler_core_id(int sched_id) {
  Context *context;

  context = mm_get_context(ar_get_core_id());

  if (context->pr_num_cores > 1) {
    ar_assert((sched_id >= 0) && (sched_id < context->pr_num_schedulers));
    return context->pr_sched_core_ids[sched_id];
  }
  else {
    ar_assert(!sched_id);
    return 0;
  }
}


// ===========================================================================
// pr_core_worker_id()          Translate a unique core ID to a worker ID
// ===========================================================================
// * INPUTS
//   int core_id                Core ID that needs translation
//
// * RETURN VALUE
//   int                        Worker ID of core
// ===========================================================================
int pr_core_worker_id(int core_id) {
  Context *context;

  context = mm_get_context(ar_get_core_id());

  if (context->pr_num_cores > 1) {
    ar_assert((core_id >= 0) && (core_id < context->pr_num_cores));
    ar_assert(context->pr_core_work_ids[core_id] > -1);
    return context->pr_core_work_ids[core_id];
  }
  else {
    ar_assert(!core_id);
    return 0;
  }
}


// ===========================================================================
// pr_core_scheduler_id()       Translate a unique core ID to a scheduler ID
// ===========================================================================
// * INPUTS
//   int core_id                Core ID that needs translation
//
// * RETURN VALUE
//   int                        Scheduler ID of core
// ===========================================================================
int pr_core_scheduler_id(int core_id) {
  Context *context;

  context = mm_get_context(ar_get_core_id());

  if (context->pr_num_cores > 1) {
    ar_assert((core_id >= 0) && (core_id < context->pr_num_cores));
    ar_assert(context->pr_core_sched_ids[core_id] > -1);
    return context->pr_core_sched_ids[core_id];
  }
  else {
    ar_assert(!core_id);
    return 0;
  }
}


// ===========================================================================
// pr_core_child_index()        Translate a unique core ID to an index of our
//                              children array
// ===========================================================================
// * INPUTS
//   int core_id                Core ID that needs translation
//
// * RETURN VALUE
//   int                        Child index of core
// ===========================================================================
int pr_core_child_index(int core_id) {
  Context *context;

  context = mm_get_context(ar_get_core_id());
  ar_assert(context->pr_role == ROLE_SCHEDULER);

  if (context->pr_num_cores > 1) {
    ar_assert((core_id >= 0) && (core_id < context->pr_num_cores));
    ar_assert(context->pr_core_child_idx[core_id] > -1);
    return context->pr_core_child_idx[core_id];
  }
  else {
    ar_assert(!core_id);
    return 0;
  }
}


// ===========================================================================
// pr_core_arch_bid_cid()       Translate a unique core ID to the arch-level
//                              board ID and core ID
// ===========================================================================
// * INPUTS
//   int core_id                Core ID that needs translation
//
// * OUTPUTS
//   int *ret_bid               Architectural-level board ID of core
//   int *ret_cid               Architectural-level core ID of core
// ===========================================================================
void pr_core_arch_bid_cid(int core_id, int *ret_bid, int *ret_cid) {
  Context *context;

  ar_assert(ret_bid);
  ar_assert(ret_cid);
  context = mm_get_context(ar_get_core_id());

  if (context->pr_num_cores > 1) {
    ar_assert((core_id >= 0) && (core_id < context->pr_num_cores));
    ar_assert(context->pr_core_bid_cid[core_id] > -1);
    *ret_bid = context->pr_core_bid_cid[core_id] >> 8;
    *ret_cid = context->pr_core_bid_cid[core_id] & 0xFF;
  }
  else {
    ar_assert(!core_id);
    *ret_bid = ar_get_board_id();
    *ret_cid = ar_get_core_id();
  }
}

