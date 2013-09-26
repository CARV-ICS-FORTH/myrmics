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
// Abstract      : Dependency analysis subsystem
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: dependency.c,v $
// CVS revision  : $Revision: 1.24 $
// Last modified : $Date: 2012/12/20 12:13:16 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <syscall.h>


// ===========================================================================
// pr_dep_start_stop()           Starts or ends the dependency analysis for 
//                               the argument list given to the function.
//                               Arguments are groupped together per
//                               responsible scheduler. Downstream messages 
//                               are sent for each child scheduler. Any
//                               remaining local arguments are handled after
//                               that.
//
//                               When invoked in start mode, the purpose of
//                               this function is to notify downstream
//                               schedulers (including ourselves) to start the
//                               upstream route-finding process for each
//                               argument, i.e.  discover for each dependency
//                               the path from the dependency up to the 
//                               highest level the parent task is enqueued.
//
//                               When invoked in stop mode, the function 
//                               notifies downstream schedulers (including 
//                               ourselves) to notify the next task enqueued
//                               for each argument, or go up the region tree
//                               if nothing else is enqueued.
// ===========================================================================
// * INPUTS
//   int start_mode             Invocation mode. 1: start mode, 0: stop mode
//   unsigned int par_task_id   ID of the parent task (used only in start mode)
//   unsigned int task_id       ID of the new child task being spawned (start
//                              mode) or the task being collected (stop mode)
//   int num_agrs               Number of dependency arguments to handle
//   void **args                Dependency argument array
//   unsigned char *dep_flags   Dependency flags array
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_dep_start_stop(int start_mode, unsigned int par_task_id, 
                      unsigned int task_id, int num_args, void **args, 
                      unsigned char *dep_flags) {

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
  int                   cur;
  int                   i;


  // Get context
  context = mm_get_context(ar_get_core_id());


  // =========================================================================
  // If we have children schedulers, separate arguments into those handled by
  // local scheduler and others handled by our children.
  // =========================================================================

  // Shortcut for leaf schedulers: everything is local
  if (context->pr_scheduler_level == 0) {
    local_sched.deps = args;
    local_sched.flags = dep_flags;
    local_sched.num_deps = num_args;
    trie = NULL;
  }

  else {

    // Create a trie to store per-scheduler packing needs. We'll be using IDs
    // from 1 to context->pr_num_schedulers, so the MSB is the log2 of
    // context->pr_num_schedulers. The kt_int_log2() function will return the
    // MSB position correctly, even for non-power-of-2 values.
    ar_assert(trie = kt_alloc_trie(kt_int_log2(context->pr_num_schedulers), 0));

    // Nothing so far is local
    local_sched.deps = NULL;
    local_sched.flags = NULL;
    local_sched.num_deps = 0;

    // For all arguments
    for (i = 0; i < num_args; i++) {

      // Ignore by-value arguments (i.e. EXACTLY objects that are safe and
      // read-only; we don't handle safe regions or other variations yet).
      if (dep_flags[i] == SYS_TYPE_BYVALUE_ARG) {
        continue;
      }

      // Find which scheduler should we ask for this region or object
      if (dep_flags[i] & SYS_TYPE_REGION_ARG) {
        sched_id = mm_distr_region_sched_id((rid_t) args[i]);
      }
      else {
        sched_id = mm_distr_object_sched_id(args[i]);
      }

      // TODO these are most probably a user error: add error reporting here
      // instead of aborting.
      ar_assert(sched_id >= 0);
      ar_assert(sched_id != context->pr_parent_sched_id);

      // Is it ours?
      if (sched_id == context->pr_scheduler_id) {
        local_sched.deps = kt_realloc(local_sched.deps, 
                                      (local_sched.num_deps + 1) * 
                                      sizeof(void *));
        local_sched.flags = kt_realloc(local_sched.flags, 
                                      (local_sched.num_deps + 1) * 
                                      sizeof(unsigned char));
        local_sched.deps[local_sched.num_deps] = args[i];
        local_sched.flags[local_sched.num_deps] = dep_flags[i];
        local_sched.num_deps++;
      }
      
      // It belongs to one of our children
      else {

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
          per_sched->deps[per_sched->num_deps] = args[i];
          per_sched->flags[per_sched->num_deps] = dep_flags[i];
          per_sched->num_deps++;
        }

        // Create new entry for this scheduler
        else {
          per_sched = kt_malloc(sizeof(per_sched_type));
          per_sched->deps = kt_malloc(sizeof(void *));
          per_sched->flags = kt_malloc(sizeof(unsigned char));
          per_sched->deps[0] = args[i];
          per_sched->flags[0] = dep_flags[i];
          per_sched->num_deps = 1;

          // Insert it into the trie
          ar_assert(!kt_trie_insert(trie, sched_id, per_sched));
        }
      }
    }
  }


  // =========================================================================
  // For each child scheduler we need to communicate with, send to it a 
  // message to start/stop dependency analysis on these arguments (or send 
  // them further downstream).
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

        if (start_mode) {
          req->type = EXT_REQ_DEP_START;
        }
        else {
          req->type = EXT_REQ_DEP_STOP;
        }
        req->core_id = context->pr_core_id;
        req->req_id  = context->pr_message_id;
        req->ptr = (void *) par_task_id;
        req->num_regions = task_id;

        // Populate message array and bitmaps
        for (req->num_ptrs = 0, req->size = 0, req->region = 0;
             (req->num_ptrs < PR_REQ_MAX_SIZE) && (cur < per_sched->num_deps);
             req->num_ptrs++, cur++) {

          req->data[req->num_ptrs] = per_sched->deps[cur];

          req->size |= (per_sched->flags[cur] & SYS_TYPE_IN_ARG) 
                     ? (1 << (2 * req->num_ptrs))     : 0;
          req->size |= (per_sched->flags[cur] & SYS_TYPE_OUT_ARG) 
                     ? (1 << (2 * req->num_ptrs + 1)) : 0;

          req->region |= (per_sched->flags[cur] & SYS_TYPE_SAFE_ARG)
                       ? (1 << (2 * req->num_ptrs))     : 0;
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
  // Finally, for all arguments handled locally, begin the routing upstream
  // towards the parent scheduler (in start mode) or the next-in-line or 
  // all-free upstream notification (in stop mode)
  // =========================================================================
  for (i = 0; i < local_sched.num_deps; i++) {
    
    // Ignore by-value arguments once more (in case we're in the local-only mode
    // and we haven't scanned the array and discarded them already)
    if (local_sched.flags[i] == SYS_TYPE_BYVALUE_ARG) {
      continue;
    }

    if (start_mode) {
      pr_dep_route(par_task_id, task_id, 
                   local_sched.deps[i], local_sched.flags[i], 1,
                   0, NULL, 0);
    }
    else {
      pr_dep_next(task_id, local_sched.deps[i], 
                  (local_sched.flags[i] & SYS_TYPE_REGION_ARG) ? 1 : 0,
                  0, 0, 0);
    }
  }


  // Free local deps arrays (only if not in shortcut mode, otherwise these
  // are just pointers to the function arrays which will be freed by the 
  // caller)
  if (context->pr_scheduler_level > 0) {
    kt_free(local_sched.deps);
    kt_free(local_sched.flags);
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_dep_route()                Starts or continues the upstream 
//                               route-finding from a dependency towards the
//                               parent task scheduler. Walks the region
//                               tree towards the root until the parent task
//                               scheduler is reached and notes the highest
//                               level that the parent task is enqueued in.
// ===========================================================================
// * INPUTS
//   unsigned int parent_task_id ID of the parent task
//   unsigned int new_task_id    ID of the new child task being spawned
//   void *dep                   The dependency we're interested in
//   unsigned char dep_flags     Flags of the dependency we're interested in
//   int dep_is_local            0: dependency belongs to one of our children
//                                  schedulers, who has started the upstream
//                                  route-finding and has already found
//                                  parent region IDs
//                               1: dependency belongs to us, and so there
//                                  are no existing parent region IDs
//   int num_remote_parents      Number of remote parent regions already found
//   rid_t *remote_parents       Remote parent region IDs array
//   rid_t highest_remote_parent Highest-level region ID so far found to 
//                               contain the parent task ID, or 0 if no such
//                               region has been found, or -1 to denote that
//                               the parent task is enqueued on the dependency
//                               object by itself (and that's an object, not
//                               a region -- if it's enqueued on the dep region,
//                               then its normal rid will be here)
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_dep_route(unsigned int parent_task_id, unsigned int new_task_id,
                 void *dep, unsigned char dep_flags, int dep_is_local,
                 int num_remote_parents, rid_t *remote_parents,
                 rid_t highest_remote_parent) {

  Context               *context;
  rid_t                 par_rids[PR_REQ_MAX_SIZE];
  int                   num_parents;
  int                   ret;
  int                   last_parent_local;
  int                   parent_sched_id;
  PrTaskDescr           *parent_task;
  rid_t                 highest_rid;
  PrObjDep              *obj_dep;
  MmRgnTreeNode         *reg_node;
  List                  *dep_queue;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  int                   i;

  
  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(dep);

  // Who is the responsible scheduler for the parent task?
  parent_sched_id = parent_task_id >> (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS);

  // What's our starting point? It may either be a local object/region (dep) if
  // the routing begins from a local object/region, or the last remote parent
  // (remote_parents[num_remote_parents-1]) that must be a local region, if the
  // routing has begun from a child scheduler.
  if (dep_is_local) {
    ar_assert(!num_remote_parents);
    ar_assert(!remote_parents);

    // Discover and put in the array the parent region of the local argument
    if (dep_flags & SYS_TYPE_REGION_ARG) {

      // Find parent
      ar_assert(!mm_region_query_region((rid_t) dep, &par_rids[0]));

      // See if the parent task is enqueued in the dep region queue itself
      ar_assert(kt_trie_find(context->mm_local_rids, (rid_t) dep,
                             (void *) &reg_node));
      if (pr_dep_search_queue(reg_node->dep_queue, parent_task_id)) {
        highest_rid = (rid_t) dep;
      }
      else {
        highest_rid = 0;
      }
    }
    else {
      
      // Find parent. TODO: mm_region_query_pointer() will fail and this
      // assertion will crash the system if an invalid user pointer reaches
      // this step. Fix this to terminate gracefully...
      if(mm_region_query_pointer(dep, NULL, &par_rids[0])) {
        kt_printf("%d: Failed query_pointer 0x%08X\r\n", context->pr_core_id, dep);
        ar_abort();
      }
      // FIXME 
      // ar_assert(!mm_region_query_pointer(dep, NULL, &par_rids[0]));
      ar_assert(kt_trie_find(context->mm_local_rids, par_rids[0],
                             (void *) &reg_node));

      // See if the parent task is enqueued in the dep object queue, if
      // it has such a queue
      highest_rid = 0;
      if (kt_trie_find(reg_node->obj_dep_queues, (size_t) dep,
                        (void *) &obj_dep)) {
        dep_queue = obj_dep->dep_queue;
        ar_assert(dep_queue);
        if (pr_dep_search_queue(dep_queue, parent_task_id)) {
          highest_rid = -1; // Special number, not used for region IDs, which
                            // use only positive numbers (up to 31 bits in
                            // 32-bit archs)
        }
      }
    }
    
    num_parents = 1;
  }
  else {
    ar_assert(remote_parents);
    ar_assert(num_remote_parents);

    // Copy existing remote parents into the local array
    for (i = 0; i < num_remote_parents; i++) {
      par_rids[i] = remote_parents[i];
    }

    ar_assert(num_remote_parents + 1 < PR_REQ_MAX_SIZE); // If this fails,
                                                         // expand
                                                         // EXT_REQ_DEP_ROUTE
                                                         // to be in multiple
                                                         // parts.
    num_parents = num_remote_parents;
    highest_rid = highest_remote_parent;
  }



  // We've got our first local parent region, now go upwards the region tree as
  // much as possible locally
  last_parent_local = 1;
  ar_assert(num_parents > 0);
  while (1) {
    ar_assert(num_parents + 1 < PR_REQ_MAX_SIZE); // If this fails, expand
                                                  // EXT_REQ_DEP_ROUTE to be
                                                  // in multiple parts.

    // See if the last parent is actually local and discover the new parent
    ret = mm_region_query_region(par_rids[num_parents - 1], 
                                 &par_rids[num_parents]);

    // If last parent not local...
    if (ret == ERR_NO_SUCH_REGION) {
      // ... then its parent belongs to our parent scheduler. Make sure 
      // we have a parent scheduler and stop here.
      ar_assert(context->pr_parent_sched_id != -1);
      last_parent_local = 0;
      break;
    }
    ar_assert(ret == 0);

    // Find if parent task ID is enqueued in the last parent and remember the
    // highest region that fulfills this condition
    ar_assert(kt_trie_find(context->mm_local_rids, par_rids[num_parents - 1],
                           (void *) &reg_node));
    if (pr_dep_search_queue(reg_node->dep_queue, parent_task_id)) {
      highest_rid = par_rids[num_parents - 1];
    }

    // Is the last parent the NULL top-level region (which does not have
    // a parent)?
    if (par_rids[num_parents] == 0) {

      ar_assert(par_rids[num_parents - 1] == 1); // ID of the NULL region

      // Stop and check that we are the top-level scheduler in this case
      ar_assert(context->pr_parent_sched_id == -1);

      // If we reached the top-level region, we must also be the responsible
      // scheduler for this task
      ar_assert(context->pr_scheduler_id == parent_sched_id);
      break;
    }

    // Increase parents counter
    num_parents++;
  }


  // If we are the responsible scheduler of the parent task, start descend &
  // enqueue of the dependency. This happens because (i) if we are the parent
  // task scheduler, then by the programming model any child task spawned must
  // request only sub-regions or objects that are under the ones owned by the
  // parent task, and (ii) a scheduler is responsible for a task only if all
  // the regions/objects the task owns are his or belong to its children.
  if (parent_sched_id == context->pr_scheduler_id) {

    // For reasons explained above, we must have found the parent task enqueued
    // somewhere, otherwise this is an error. TODO: this should be a user
    // error, handle the error gracefully instead of crashing.
    ar_assert((highest_rid > 0) || (highest_rid == -1));

    // Start descend & enqueue
    pr_dep_descend_enqueue(parent_task_id, new_task_id, dep, dep_flags, 
                           num_parents, par_rids, highest_rid, 0, 0);

    // Locate parent task
    parent_task = NULL;
    kt_trie_find(context->pr_tasks, parent_task_id, (void *) &parent_task);
    ar_assert(parent_task);
    ar_assert(parent_task->id == parent_task_id);

    // Note that one more dependency has been sent to its final downstream trip
    ar_assert(parent_task->spawn_rem_deps);
    parent_task->spawn_rem_deps--;

    // If it was the final one, we must send a reply to the task spawn request
    if (!parent_task->spawn_rem_deps) {

      // Build reply
      reply = noc_msg_send_get_buf(parent_task->spawn_core_id);

      reply->type   = REPLY_SPAWN;
      reply->req_id = parent_task->spawn_msg_id;
      reply->status = 0;

      // Send it to the one who had sent the request
      ar_assert(!noc_msg_send());

      // Remember that no spawn is now active for the parent task
      parent_task->spawn_msg_id = 0;
      parent_task->spawn_core_id = -1;
    }

  }

  // Otherwise, we must send a message to our parent scheduler with the
  // so-far found dependency routing so he can continue the routing upstream
  // towards the parent task scheduler.
  else {
    ar_assert(!last_parent_local);

    // Start building new message
    req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    req->type        = EXT_REQ_DEP_ROUTE;
    req->core_id     = context->pr_core_id;
    req->req_id      = context->pr_message_id;
    req->size        = parent_task_id;
    req->region      = highest_rid;
    req->ptr         = dep;
    req->num_regions = num_parents;
    req->num_ptrs    = (new_task_id << PR_TASK_IDX_SIZE) | dep_flags;

    for (i = 0; i < num_parents; i++) {
      req->data[i] = (void *) par_rids[i];
    }

    // Send request to parent
    ar_assert(!noc_msg_send());
    
    // Increase message ID, avoiding value 0 on wrap-arounds
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_dep_descend_enqueue()      Descends the region tree using the route
//                               of parents, until the given dependency is
//                               found. If non-local parent regions are 
//                               encountered, a message is sent to the
//                               child scheduler to continue the traversal.
//
//                               The dependency is enqueued at the first 
//                               blocked level between the highest_rid where
//                               the parent task is enqueued down to the 
//                               requested obj/region dependency itself.
// ===========================================================================
// * INPUTS
//   unsigned int parent_task_id ID of the parent task
//   unsigned int new_task_id    ID of the new child task being spawned
//   void *dep                   The dependency to be enqueued
//   unsigned char dep_flags     Flags of the dependency 
//   int num_parents             Number of parent regions of dep
//   rid_t *par_rids             Remote parent region IDs array
//   rid_t highest_rid           Highest-level region ID that contains the 
//                               parent task ID, or -1 to denote that the
//                               parent task is enqueued on the dependency
//                               object by itself (and that's an object, not a
//                               region -- if it's enqueued on the dep region,
//                               then its normal rid will be here)
//   rid_t dont_incr_parent_enq  If > 0, do not increase the parent_enqs on 
//                               this region ID if we pass through it
//   int resume                  If > 0, it means we are resuming the
//                               operation after being unblocked in an 
//                               intermediate queue. This causes not to 
//                               (re)enqueue the dep into the queue it was
//                               just unblocked from.
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_dep_descend_enqueue(unsigned int parent_task_id, 
                           unsigned int new_task_id, void *dep, 
                           unsigned char dep_flags, int num_parents, 
                           rid_t *par_rids, rid_t highest_rid, 
                           rid_t dont_incr_parent_enq, int resume) {

  Context       *context;
  int           cur_route;
  int           enqueue_started;
  int           not_ours;
  int           sched_id;
  MmRgnTreeNode *reg_node;
  PrMsgReq      *req;
  List          *dep_queue;
  PrObjDep      *obj_dep;
  ListNode      *list_node;
  PrDep         *list_dep;
  PrDep         *new_pr_dep;
  rid_t         tmp_rid;
  int           ret;
  int           ret_override;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert((highest_rid > 0) || (highest_rid == -1));
  ar_assert(num_parents >= 0);

//kt_printf(
//   "%d: Desc/enq 0x%08X [%c%c%c%c%c] par 0x%X, new 0x%X, par enq rid %d:\r\n",
//   context->pr_core_id, dep,
//   (dep_flags & SYS_TYPE_IN_ARG)      ? 'I' : '-',
//   (dep_flags & SYS_TYPE_OUT_ARG)     ? 'O' : '-',
//   (dep_flags & SYS_TYPE_SAFE_ARG)    ? 'S' : '-',
//   (dep_flags & SYS_TYPE_REGION_ARG)  ? 'R' : '-',
//   (dep_flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-',
//   parent_task_id, new_task_id, highest_rid);
//for (i = 0; i < num_parents; i++) {
//  kt_printf("    Parent[%d] = %d\r\n", i, par_rids[i]);
//}

  // We need to enqueue the new task at a point between (a) the highest point
  // the parent task is enqueued and (b) the target dependency itself. We need
  // to traverse the path from (a) down to (b) and enqueue at the highest 
  // blocked point in between.

  // Since this function can be called from any scheduler during the region
  // tree traversal, first we need to determine if we already are between
  // (a) and (b) or if we still descending towards the starting point (a).
  if (((!(dep_flags & SYS_TYPE_REGION_ARG)) && (highest_rid == -1)) ||
      ((dep_flags & SYS_TYPE_REGION_ARG) && (highest_rid == (rid_t) dep))) {
    enqueue_started = 0;
  }
  else {
    // Find which parent is the highest rid
    cur_route = -1;
    for (i = num_parents - 1; i >= 0; i--) {
      if (par_rids[i] == highest_rid) {
        cur_route = i;
        break;
      }
    }

    // If highest_rid is one of the remaining parents, we are still outside
    // (a)-(b).  Otherwise, we assume we have reached highest_rid in a past
    // call of this function and we are inside (a)-(b).
    if (cur_route > -1) {
      enqueue_started = 0;
    }
    else {
      enqueue_started = 1;
    }
  }

  // For all remaining parents
  not_ours = 0;
  for (i = num_parents - 1; i >= 0; i--) {

    // Who manages this region?
    sched_id = mm_distr_region_sched_id(par_rids[i]);
    ar_assert(sched_id >= 0);

    if (sched_id != context->pr_scheduler_id) {
      not_ours = 1;
      break;
    }

    // If we are inside (a)-(b), enqueue at the first place which is blocked,
    // ignoring (passing through) all other regions that are unblocked.
    reg_node = NULL;
    if (enqueue_started) {

      // Find region tree node
      ar_assert(kt_trie_find(context->mm_local_rids, par_rids[i],
                             (void *) &reg_node));

      // Mark that we came here because our parent requested it, if this is so
      if ((par_rids[i] != highest_rid) && 
          (par_rids[i] != dont_incr_parent_enq)) {
        if (dep_flags & SYS_TYPE_OUT_ARG) {
          reg_node->parent_enqs_rw++;
        }
        else {
          reg_node->parent_enqs_ro++;
        }
      }

      // Blocked here? Enqueue and we're done for now. For the first time
      // we encounter this, we have to ignore the enqueue if we're resuming
      // (because this is the region we're resuming from).
      if (!resume) {
        if (kt_list_size(reg_node->dep_queue)) {

          // One exception to blocking the dep here, is if it has a read-only
          // and here is a region where everybody ahead of us already has
          // read-only access. In this case, we let the dependency pass
          // through.
//kt_printf("region %d [queue size %d], examining task 0x%08X, real dep 0x%08X\r\n", reg_node->id, kt_list_size(reg_node->dep_queue), new_task_id, dep);
          if (!(dep_flags & SYS_TYPE_OUT_ARG)) {
            for (list_node = kt_list_head(reg_node->dep_queue);
                 list_node;
                 list_node = kt_list_next(list_node)) {

              ar_assert(list_node->data);
              list_dep = (PrDep *) list_node->data;

//kt_printf("> task 0x%08X, access [%c%c%c]\r\n", list_dep->task_id, (list_dep->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (list_dep->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (list_dep->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-');

              // Stop if we find a non-compatible access pattern (read-write
              // or pending (i.e. possibly waiting for read-write children))
              if ((list_dep->flags & SYS_TYPE_OUT_ARG) ||
                  (list_dep->flags & SYS_TYPE_DEP_PENDING)) {
//kt_printf("> stopping\r\n");
                break;
              }
            }

            // If we managed to get through the whole list, everyone is
            // compatible, so don't block here
            if (!list_node) {
//kt_printf("Averted!\r\n");
              goto blocking_averted;
            }
          }

//kt_printf("*1* region %d, task 0x%08X blocked, real dep 0x%08X\r\n", reg_node->id, new_task_id, dep); ListNode *foo; PrDep *bar; for (foo = kt_list_head(reg_node->dep_queue); foo; foo = kt_list_next(foo)) {ar_assert(foo->data); bar = foo->data; kt_printf("0x%08X [%c%c%c] ", bar->task_id, (bar->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (bar->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (bar->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-'); } kt_printf("\r\n");
          
          // Block here at the end of the queue
          ret = pr_dep_enqueue(reg_node->dep_queue, dep, dep_flags,
                               new_task_id, parent_task_id, highest_rid, 
                               0, 0, par_rids, i, NULL);
          ar_assert(ret != 1);
          return 0;
        }
      }
      else {
        resume = 0;
      }
    }

blocking_averted:

    // Maybe this is point (a)? If so, remember that we are inside (a)-(b)
    // for the next regions we traverse.
    if (!enqueue_started) {
      if (par_rids[i] == highest_rid) {
        enqueue_started = 1;
        ar_assert(kt_trie_find(context->mm_local_rids, par_rids[i],
                               (void *) &reg_node));
      }
    }

    // Remember that we are sending one more task to a child region of ours.
    // Note that we objects are tracked by different counters.
    if (enqueue_started) {
      ar_assert(reg_node);
      if ((i > 0) || (dep_flags & SYS_TYPE_REGION_ARG)) {
        if (dep_flags & SYS_TYPE_OUT_ARG) {
          reg_node->children_enqs_rw++;
        }
        else {
          reg_node->children_enqs_ro++;
        }
//kt_printf("%d: Region %d at %d/%d children_enqs\r\n", context->pr_core_id, reg_node->id, reg_node->children_enqs_ro, reg_node->children_enqs_rw);
      }
      else {
        if (dep_flags & SYS_TYPE_OUT_ARG) {
          reg_node->obj_enqs_rw++;
        }
        else {
          reg_node->obj_enqs_ro++;
        }
//kt_printf("%d: Region %d at %d/%d obj_enqs\r\n", context->pr_core_id, reg_node->id, reg_node->obj_enqs_ro, reg_node->obj_enqs_rw);
      }
    }
  }

  // If we aborted the loop above, we found a parent region that's not ours.
  // Even if we did not abort, the final object or region may not be ours.
  // Discover this case here.
  if (!not_ours) {
    // Region?
    if (dep_flags & SYS_TYPE_REGION_ARG) {
      sched_id = mm_distr_region_sched_id((rid_t) dep);
      if (sched_id != context->pr_scheduler_id) {
        not_ours = 1;
      }
    }
    // Object?
    else {
      // If it's an object, its parent _should_ be ours and therefore must
      // have been discovered as the last parent rid, so this is only a 
      // sanity check.
      ar_assert(!mm_region_query_pointer((void *) dep, NULL, NULL));
    }
  }

  // If remaining parents and/or final obj/region are not ours, send a messsage
  // to our child scheduler to continue the descension.
  if (not_ours) {

    // Start building new message
    req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

    req->type        = EXT_REQ_DEP_ENQUEUE;
    req->core_id     = context->pr_core_id;
    req->req_id      = context->pr_message_id;
    req->size        = parent_task_id;
    req->region      = highest_rid;
    req->ptr         = dep;
    req->num_regions = i + 1;
    req->num_ptrs    = (new_task_id << PR_TASK_IDX_SIZE) | dep_flags;

    for (; i >= 0; i--) {
      req->data[i] = (void *) par_rids[i];
    }

    // Send request to child
    ar_assert(!noc_msg_send());
    
    // Increase message ID, avoiding value 0 on wrap-arounds
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // We're done
    return 0;
  }

  // Finally, if we reached this point, the dep is ours and we must enqueue
  // the new task at this queue. Locate the queue.
  if (dep_flags & SYS_TYPE_REGION_ARG) {
    ar_assert(kt_trie_find(context->mm_local_rids, (size_t) dep, 
                           (void *) &reg_node));
    dep_queue = reg_node->dep_queue;

    // If we came here by traversing our parent, remember it
    if (enqueue_started && (reg_node->id != dont_incr_parent_enq)) {
      if (dep_flags & SYS_TYPE_OUT_ARG) {
        reg_node->parent_enqs_rw++;
      }
      else {
        reg_node->parent_enqs_ro++;
      }
    }
  }
  else {
    // Object: if it has other dependencies, it will be in our trie. If it's
    // not, we create its dep queue now.
    ar_assert(!mm_region_query_pointer(dep, NULL, &tmp_rid)); // sanity check
    ar_assert(tmp_rid == par_rids[0]);
    ar_assert(kt_trie_find(context->mm_local_rids, par_rids[0],
                           (void *) &reg_node));

    if (!kt_trie_find(reg_node->obj_dep_queues, (size_t) dep, 
                      (void *) &obj_dep)) {
      // Create a new obj dep structure and a dependency queue
      obj_dep = kt_malloc(sizeof(PrObjDep));
      obj_dep->parent_enqs_ro = 0;
      obj_dep->parent_enqs_rw = 0;
      obj_dep->dep_queue = kt_alloc_list();

      // Insert it into the trie
      ar_assert(!kt_trie_insert(reg_node->obj_dep_queues, (size_t) dep,
                                obj_dep));

    }
    dep_queue = obj_dep->dep_queue;

    // If we came here by traversing our parent, remember it
    if (enqueue_started) {
      if (dep_flags & SYS_TYPE_OUT_ARG) {
        obj_dep->parent_enqs_rw++;
      }
      else {
        obj_dep->parent_enqs_ro++;
      }
    }

  }
  ar_assert(dep_queue);

  // Enqueue. If the highest_rid refers to the same obj/region as dep, then we
  // must be inserted ahead of the parent task. Otherwise, insert at the tail.
  if ((highest_rid == -1) || 
      ((dep_flags & SYS_TYPE_REGION_ARG) && (highest_rid == (size_t) dep))) {
//kt_printf("*2* dep 0x%08X enq task 0x%08X\r\n", dep, new_task_id); ListNode *foo; PrDep *bar; for (foo = kt_list_head(dep_queue); foo; foo = kt_list_next(foo)) {ar_assert(foo->data); bar = foo->data; kt_printf("0x%08X [%c%c%c] ", bar->task_id, (bar->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (bar->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (bar->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-'); } kt_printf("\r\n");
    ret = pr_dep_enqueue(dep_queue, dep, dep_flags, new_task_id, 
                         parent_task_id, highest_rid, parent_task_id,
                         1, NULL, 0, &new_pr_dep);
  }
  else {
//kt_printf("*3* dep 0x%08X enq task 0x%08X\r\n", dep, new_task_id); ListNode *foo; PrDep *bar; for (foo = kt_list_head(dep_queue); foo; foo = kt_list_next(foo)) {ar_assert(foo->data); bar = foo->data; kt_printf("0x%08X [%c%c%c] ", bar->task_id, (bar->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (bar->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (bar->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-'); } kt_printf("\r\n");
    ret = pr_dep_enqueue(dep_queue, dep, dep_flags, new_task_id, 
                         parent_task_id, highest_rid, 0, 0, NULL, 0, 
                         &new_pr_dep);
  }

  // If we are at the head of this final queue, we're going to proceed below to
  // inform the new task that this dependency is ready. Before we do that, we
  // have to check one more case: even if we are NOT at the head of the queue,
  // it may happen that everybody before us is a read-only dependency and we
  // want also to have a read-only access. Discover this here and override the
  // ret=0 value so we can include this case.
  if (!ret && !(dep_flags & SYS_TYPE_OUT_ARG)) {

    // Assume the best
    ret_override = 1;

    // For all the deps in the queue
    for (list_node = kt_list_head(dep_queue);
         list_node;
         list_node = kt_list_next(list_node)) {

      ar_assert(list_node->data);
      list_dep = (PrDep *) list_node->data;

      // Stop when we find ourselves in the queue
      if (list_dep->task_id == new_task_id) {
        break;
      }

      // Stop when we find a non-pending or a non-read-only dep and revert
      // ret_override value
      if ((list_dep->flags & SYS_TYPE_DEP_PENDING) ||
          (list_dep->flags & SYS_TYPE_OUT_ARG)) {
        ret_override = 0;
        break;
      }
    }

    ar_assert(list_node); // sanity check: we must have found ourselves
  }
  else {
    ret_override = 0;
  }
  
//kt_printf("%d: *** descend task 0x%08X dep 0x%08X [%d %d]\r\n   ", context->pr_core_id, new_task_id, dep, ret, ret_override); ListNode *foo; PrDep *bar; for (foo = kt_list_head(dep_queue); foo; foo = kt_list_next(foo)) {ar_assert(foo->data); bar = foo->data; kt_printf("0x%08X [%c%c%c] ", bar->task_id, (bar->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (bar->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (bar->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-'); } kt_printf("\r\n");

  // If we are at the head (or we have head privileges, as fixed above), tell
  // the responsible scheduler of the new task that this dependency is ready.
  // For read/write regions, we also have to check if we have delegated any
  // part to children regions or objects; if so, we can't yet grant the whole
  // region to the task at the head of the queue: this will be done later on
  // when all children have finished. For read-only regions, we must only check
  // that we have not delegated any part to children regions or objects in
  // read/write mode -- read-only delegations are compatible (everybody just
  // reads). Finally, for regions that we have already found others
  // be in read-only access, we don't need any check at all: any children must
  // be in read-only access as well, otherwise it's a progammer error.
  if (((ret == 1) && 
       (!(dep_flags & SYS_TYPE_REGION_ARG) || 
        ((dep_flags & SYS_TYPE_OUT_ARG) &&
         (reg_node->children_enqs_rw + reg_node->children_enqs_ro == 0) && 
         (reg_node->obj_enqs_rw + reg_node->obj_enqs_ro == 0)) ||
        (!(dep_flags & SYS_TYPE_OUT_ARG) &&
         (reg_node->children_enqs_rw == 0) && 
         (reg_node->obj_enqs_rw == 0)))) ||
      (ret_override == 1)) {

    // Whose is it?
    sched_id = new_task_id >> (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS);

    if (sched_id == context->pr_scheduler_id) {
      // Local call
//kt_printf("%d: descend->at_head dep 0x%08X task_id 0x%08X\r\n", context->pr_core_id, dep, new_task_id);
      ar_assert(!pr_dep_at_head_of_queue(new_task_id, dep));
    }
    else {
      // Start building new message
      ar_assert(context->pr_parent_sched_id != -1);

      req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));

      req->type    = REQ_DEP_OK;
      req->core_id = context->pr_core_id;
      req->req_id  = context->pr_message_id;
      req->size    = new_task_id;
      req->ptr     = dep;

      // Send request to parent
      ar_assert(!noc_msg_send());
      
      // Increase message ID, avoiding value 0 on wrap-arounds
      context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
    }

    // Also remember that this dependency is non-pending:
    // pr_dep_at_head_of_queue() will do this for the task structure, but
    // we also need it in the dependency queue.
    ar_assert(new_pr_dep->flags & SYS_TYPE_DEP_PENDING);
    new_pr_dep->flags ^= SYS_TYPE_DEP_PENDING;
  }

//if (dep == 0x00000006) {kt_printf("%d: +++ descend task 0x%08X dep 0x%08X [%d %d]\r\n   ", context->pr_core_id, new_task_id, dep, ret, ret_override); ListNode *foo; PrDep *bar; for (foo = kt_list_head(dep_queue); foo; foo = kt_list_next(foo)) {ar_assert(foo->data); bar = foo->data; kt_printf("0x%08X [%c%c%c] ", bar->task_id, (bar->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (bar->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (bar->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-'); } kt_printf("\r\n");}

  // Success
  return 0;
}


// ===========================================================================
// pr_dep_next()                Removes a task ID from a dependency queue
//                              and notifies the next task waiting in this
//                              queue; if none are waiting, it notifies one
//                              level up the region tree that this queue is
//                              empty.
// ===========================================================================
// * INPUTS
//   unsigned int task_id       Task ID to locate in the queue and remove
//                              when child_decr_mode == 0
//   void *dep                  The dependency whose queue we are searching
//   int dep_is_region          0: dep is an object pointer
//                              1: dep is a region ID
//   int child_decr_mode        0: remove task_id from dep and go up from there
//                              1: just decrement child_decr from dep and go up 
//                                 from there
//   int child_decr_ro          Children enqueues to decrement from dep when
//                              child_decr_mode == 1, for read-only enqueues
//   int child_decr_rw          As above, for read/write enqueues
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_dep_next(unsigned int task_id, void *dep, int dep_is_region,
                int child_decr_mode, int child_decr_ro, int child_decr_rw) {

  Context       *context;
  MmRgnTreeNode *reg_node;
  List          *dep_queue;
  ListNode      *pos;
  ListNode      *head;
  ListNode      *cur;
  ListNode      *nxt;
  PrObjDep      *obj_dep;
  PrDep         *pr_dep;
  rid_t         par_rid;
  int           wake_others;
  int           sched_id;
  PrMsgReq      *req;
  MmRgnTreeNode *par_reg_node;
  int           read_only_mode;


  // Ok, so this function uses "harmful" gotos. Unfortunately, it's a loop that
  // we need to enter and reenter from various points. It can be written more
  // clearly using recursion, but I'd prefer to avoid it, because it's possible
  // for a large amount of tasks to be blocked after a bottleneck task. It can
  // also be written more clearly if it's broken into small functions, but it's
  // a shame to make it so inefficient...


  // Get context
  context = mm_get_context(ar_get_core_id());

  if (child_decr_mode) {
    sched_id = context->pr_scheduler_id;
    par_rid = (rid_t) dep;
    goto decrement_parent;
  }

  // Find the dependency queue
  if (dep_is_region) {

    // Locate the region itself
    ar_assert(kt_trie_find(context->mm_local_rids, (rid_t) dep,
                           (void *) &reg_node));
    dep_queue = reg_node->dep_queue;

    obj_dep = NULL;
  }
  else {
    
    // Find the object's parent region
    ar_assert(!mm_region_query_pointer(dep, NULL, &par_rid));
    ar_assert(kt_trie_find(context->mm_local_rids, par_rid,
                           (void *) &reg_node));

    // Locate the dependency queue
    ar_assert(kt_trie_find(reg_node->obj_dep_queues, (size_t) dep,
                           (void *) &obj_dep));
    dep_queue = obj_dep->dep_queue;
  }

//kt_printf("%d: *** dep_next dep 0x%08X\r\n   ", context->pr_core_id, dep); ListNode *foo; PrDep *bar; for (foo = kt_list_head(dep_queue); foo; foo = kt_list_next(foo)) {ar_assert(foo->data); bar = foo->data; kt_printf("0x%08X [%c%c%c] ", bar->task_id, (bar->flags & SYS_TYPE_IN_ARG) ? 'I' : '-', (bar->flags & SYS_TYPE_OUT_ARG) ? 'O' : '-', (bar->flags & SYS_TYPE_DEP_PENDING) ? 'P' : '-'); } kt_printf("\r\n");

  // Locate task ID in the queue
  pos = pr_dep_search_queue(dep_queue, task_id);
  ar_assert(pos);

  // See if we must wake anyone after we delete ourselves from the queue. There
  // are two cases where we must NOT wake anyone else: (a) we may be a function
  // that recurses, so its child will be inserted before us and therefore will
  // have taken the head of the queue already and marked non-pending, or (b)
  // we may be a read-only dep which is non-pending along with others, starting
  // from the head of the queue, including us, and possibly even further beyond
  // us. In both cases, it means at least one more task is already active and
  // has this dependency marked non-pending. To check this, it suffices to 
  // see if after our deletion there's at least one task left at the head of
  // the queue and marked non-pending. 
  //
  // Reminder: in the case that the whole queue is read-only and all have
  // access, any new task coming which is also read-only will be automatically
  // be marked non-pending by the code at pr_dep_descend_enqueue().
  if (pos != kt_list_head(dep_queue)) {
    wake_others = 0;
  }
  else {
    nxt = kt_list_next(kt_list_head(dep_queue));
    if (nxt) {
      pr_dep = (PrDep *) nxt->data;
      ar_assert(pr_dep);
      if (pr_dep->flags & SYS_TYPE_DEP_PENDING) {
        wake_others = 1;
      }
      else {
        wake_others = 0;
      }
    }
    else {
      wake_others = 1;
    }
  }

  // Delete this dependency from the list and free the PrDep
  kt_list_delete(dep_queue, pos, kt_free);

  // If other task(s) already have privileges, we're done
  if (!wake_others) {
    return 0;
  }


process_next:

  // Is there anyone else at the head of the list?
  head = kt_list_head(dep_queue);
  if (head) {

    pr_dep = head->data;
    ar_assert(pr_dep);

//kt_printf("> process_next dep 0x%08X, task 0x%08X\r\n", dep, pr_dep->task_id);

    // Does it really wait for this queue, or was it blocked here and must
    // continue downstream?
    if (pr_dep->final_target != dep) {

      // Delete this dependency from the list, but keep the PrDep
      kt_list_delete(dep_queue, head, NULL);

      // Proceed downstream
//kt_printf("> dep 0x%08X, task 0x%08X descends for real dep 0x%08X\r\n", dep, pr_dep->task_id, pr_dep->final_target);
      pr_dep_descend_enqueue(pr_dep->par_task_id, pr_dep->task_id,
                             pr_dep->final_target, pr_dep->flags,
                             pr_dep->cur_route + 1, pr_dep->route, 
                             pr_dep->par_highest_rid, (rid_t) dep, 1);

      // This temporarily-enqueued dependency has been forwarded downstream.
      // Free the PrDep and proceed with the next one.
      kt_free(pr_dep);
      goto process_next;
    }

    // Tell the responsible scheduler of this task that this dependency is
    // ready. For regions, we also have to check if we have delegated any part
    // to children regions or objects; if so, we can't yet grant the whole
    // region to the task at the head of the queue: this will be done later on
    // when all children have finished. See the comments in the related part in
    // pr_dep_descend_enqueue() for the intricacies of read/write vs. read-only
    // modes.
    if (!dep_is_region || 
        ((pr_dep->flags & SYS_TYPE_OUT_ARG) &&
         (reg_node->children_enqs_rw + reg_node->children_enqs_ro == 0) && 
         (reg_node->obj_enqs_rw + reg_node->obj_enqs_ro == 0)) ||
        (!(pr_dep->flags & SYS_TYPE_OUT_ARG) &&
         (reg_node->children_enqs_rw == 0) && 
         (reg_node->obj_enqs_rw == 0))) {

      cur = head;
      nxt = kt_list_next(head);

      if (!(pr_dep->flags & SYS_TYPE_OUT_ARG)) {
        read_only_mode = 1;
      }
      else {
        read_only_mode = 0;
      }

      while (1) {

        // pr_dep always exists the first time we enter the loop, but may
        // not exist if we're skipping tasks that are not really waiting
        // for this dep, but continue descending instead.
        if (pr_dep) {

          // Determine responsible scheduler
          sched_id = pr_dep->task_id >> (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS);

          if (sched_id == context->pr_scheduler_id) {
            // Local call
//kt_printf(">* at_head dep 0x%08X, task 0x%08X\r\n", dep, pr_dep->task_id, pr_dep->final_target);
            ar_assert(!pr_dep_at_head_of_queue(pr_dep->task_id, dep));
          }
          else {
            // Start building new message
            ar_assert(context->pr_parent_sched_id != -1);

            req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                  context->pr_parent_sched_id));

            req->type    = REQ_DEP_OK;
            req->core_id = context->pr_core_id;
            req->req_id  = context->pr_message_id;
            req->size    = pr_dep->task_id;
            req->ptr     = dep;

            // Send request to parent
            ar_assert(!noc_msg_send());
            
            // Increase message ID, avoiding value 0 on wrap-arounds
            context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
          }

          // Also remember that this dependency is non-pending:
          // pr_dep_at_head_of_queue() will do this for the task structure, but
          // we also need it in the dependency queue.
          ar_assert(pr_dep->flags & SYS_TYPE_DEP_PENDING);
          pr_dep->flags ^= SYS_TYPE_DEP_PENDING;
        }


        // Does the original task at the head for this dependency ask for
        // read-only access? If so, is there a next one waiting?
        // And if yes, does it also want the dependency read-only? Then,
        // proceed to also wake this one up.
        if (read_only_mode && nxt) {
          cur = nxt;
          pr_dep = (PrDep *) nxt->data;
          ar_assert(pr_dep);
          if (pr_dep->flags & SYS_TYPE_OUT_ARG) {
            break;
          }
          nxt = kt_list_next(cur);
        }
        else {
          break;
        }

        // ... but before we do that, check (as we did above for the list head)
        // that the next one also really waits for this queue.
        if (pr_dep->final_target != dep) {

          // Delete this dependency from the list, but keep the PrDep
          kt_list_delete(dep_queue, cur, NULL);

          // Proceed downstream
//kt_printf(">! dep 0x%08X, task 0x%08X descends for real dep 0x%08X\r\n", dep, pr_dep->task_id, pr_dep->final_target);
          pr_dep_descend_enqueue(pr_dep->par_task_id, pr_dep->task_id,
                                 pr_dep->final_target, pr_dep->flags,
                                 pr_dep->cur_route + 1, pr_dep->route, 
                                 pr_dep->par_highest_rid, (rid_t) dep, 1);

          // This temporarily-enqueued dependency has been forwarded downstream.
          // Free the PrDep and proceed with the next one.
          kt_free(pr_dep);

          cur = NULL;
          pr_dep = NULL;
        }

      }
    }

    // Either if we notified the new task(s) or we're still waiting for
    // children regions, there's nothing more to do here...
    return 0;
  }


  // =========================================================================
  // If we reached this point, we should go up the region tree and notify
  // that we are empty, provided we have no children pending (this can happen
  // if we unblocked deps that were waiting here but their final target was
  // deeper in the tree).
  // =========================================================================

  // If we're an object, we must additionally destroy the empty list and remove
  // it and its holder structure from the parent region trie
  if (!dep_is_region) {
    
    // Remember how much we have to decrement the parent region, before we
    // destroy the object structure
    ar_assert(obj_dep);
    ar_assert(obj_dep->dep_queue == dep_queue);
    child_decr_ro = obj_dep->parent_enqs_ro;
    child_decr_rw = obj_dep->parent_enqs_rw;
    
    // Free queue and holder struct and delete from trie
    kt_free_list(dep_queue, NULL);
    ar_assert(!kt_trie_delete(reg_node->obj_dep_queues, (size_t) dep, kt_free));
    
    // Region to go up is the object parent, which should be local
    ar_assert(!mm_region_query_pointer(dep, NULL, &par_rid));
    sched_id = context->pr_scheduler_id;

  }
  else {

    // Make sure we don't have any children pending
    //
    // FIXME: does this hurt parallelism for read-only accesses? Probably 
    //        not, but needs more careful thinking...
    if ((reg_node->children_enqs_ro + reg_node->children_enqs_rw) || 
        (reg_node->obj_enqs_ro + reg_node->obj_enqs_rw)) {
      return 0;
    }

    // Region to go up is the parent of this region
    ar_assert(!mm_region_query_region((rid_t) dep, &par_rid));

    // Find out who's the scheduler to do this
    sched_id = mm_distr_region_sched_id(par_rid);
    ar_assert(sched_id >= 0);

    // We should notify the parent we done all the tasks he had delegated
    // to us by the last time we reported. Reset our parent counter.
    child_decr_ro = reg_node->parent_enqs_ro;
    child_decr_rw = reg_node->parent_enqs_rw;
    reg_node->parent_enqs_ro = 0;
    reg_node->parent_enqs_rw = 0;
  }

decrement_parent:

  // Parent region local?
  if (sched_id == context->pr_scheduler_id) {

    // Locate parent region in the region tree
    ar_assert(kt_trie_find(context->mm_local_rids, par_rid, 
                           (void *) &par_reg_node));

    // Decrement the parent's children counter
    if (!dep_is_region) {
//kt_printf("%d: Region %d minus %d/%d, now at %d/%d obj_enqs\r\n", context->pr_core_id, par_reg_node->id, child_decr_ro, child_decr_rw, par_reg_node->obj_enqs_ro - child_decr_ro, par_reg_node->obj_enqs_rw - child_decr_rw);
      ar_assert(par_reg_node->obj_enqs_ro >= child_decr_ro);
      ar_assert(par_reg_node->obj_enqs_rw >= child_decr_rw);
      par_reg_node->obj_enqs_ro -= child_decr_ro;
      par_reg_node->obj_enqs_rw -= child_decr_rw;
    }
    else {
//kt_printf("%d: Region %d minus %d/%d, now at %d/%d children_enqs\r\n", context->pr_core_id, par_reg_node->id, child_decr_ro, child_decr_rw, par_reg_node->children_enqs_ro - child_decr_ro, par_reg_node->children_enqs_rw - child_decr_rw);
      ar_assert(par_reg_node->children_enqs_ro >= child_decr_ro);
      ar_assert(par_reg_node->children_enqs_rw >= child_decr_rw);
      par_reg_node->children_enqs_ro -= child_decr_ro;
      par_reg_node->children_enqs_rw -= child_decr_rw;
    }

    // Stop going up if we're waiting for more children. Also stop going up if
    // the parent is the NULL region. By convention, the NULL region cannot be
    // delegated to other tasks than the main task, which cannot recurse. So,
    // it can't have anything waiting for it.
    //
    // FIXME: also check if this hurt parallelism for read-only accesses?
    if ((par_reg_node->children_enqs_rw + par_reg_node->children_enqs_ro) || 
        (par_reg_node->obj_enqs_rw + par_reg_node->obj_enqs_ro) || 
        (par_rid == 1)) {
      return 0;
    }

    // Proceed with the parent region to see if it has anything enqueued
    dep = (void *) par_rid;
    dep_is_region = 1;
    reg_node = par_reg_node;
    dep_queue = reg_node->dep_queue;

    // One last exception before we proceed: it might be the parent region
    // already has a non-pending task working on it. This might be (a) because
    // it's a parent task that spawned some children and we're the last one
    // to finish, but the parent task is still running, or (b) because it's
    // a read-only dep and there's multiple access to it by other tasks
    // (which extends to all their children). In any case, we must not 
    // continue to wake up something that's already awake.
    if (kt_list_size(dep_queue) && 
        !(((PrDep *) kt_list_head(dep_queue)->data)->flags & 
                                                SYS_TYPE_DEP_PENDING)) {
      return 0;
    }


    goto process_next;
  }
  else {
    // Send a message to our parent scheduler
    ar_assert(dep_is_region);
    ar_assert(sched_id == context->pr_parent_sched_id);
    ar_assert(context->pr_parent_sched_id != -1);

    req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

    req->type    = REQ_DEP_CHILD_FREE;
    req->core_id = context->pr_core_id;
    req->req_id  = context->pr_message_id;
    req->region  = par_rid;
    req->size    = child_decr_ro;
    req->ptr     = (void *) child_decr_rw;

    // Send request to parent
    ar_assert(!noc_msg_send());
    
    // Increase message ID, avoiding value 0 on wrap-arounds
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
  }

  // Success
  return 0;
}


// ===========================================================================
// pr_dep_enqueue()              Enqueues the task ID into an object/region
//                               dependency queue and sets it as pending
// ===========================================================================
// * INPUTS
//   List *dep_queue             The dependency queue of the object or region
//                               to enqueue the task at this time
//   void *final_dep             The final target dependency object or region
//                               where the task needs to be enqueued
//   unsigned char dep_flags     Dependency flags of final_dep, as defined in
//                               syscall.h. Note that SYS_TYPE_DEP_PENDING
//                               is forced to 1, regardless of its value
//                               in dep_flags.
//   unsigned int task_id        ID of the task to be added to the queue
//   unsigned int par_task_id    ID of its parent task
//   rid_t par_highest_rid       Region ID of the highest region the parent
//                               task is enqueued, or -1 for recursions
//   unsigned int ref_task_id    ID of the reference task that must exist
//                               in the queue for insertion modes 1 and 2
//   int insert_mode             0: insert at the tail of the queue
//                               1: insert before ref_task_id (i.e. closer to
//                                  the head of the queue than ref_task_id)
//                               2: insert after ref_task_id (i.e. closer to
//                                  the tail of the queue than ref_task_id)
//   rid_t *route                Array of parent regions of the final_dep
//                               (entry 0 is the immediate parent). Can be
//                               NULL only for the top-level region queue.
//   int cur_route               Current array entry which represents the
//                               region of dep
// * OUTPUTS
//   PrDep **ret_pr_dep          If not NULL, the newly allocated PrDep *
//                               is returned here
//
// * RETURN VALUE
//   int                         1: task_id is at the head of the queue
//                               0: task_id not at the head of the queue
// ===========================================================================
int pr_dep_enqueue(List *dep_queue, void *final_dep, unsigned char dep_flags,
                   unsigned int task_id, unsigned int par_task_id,
                   rid_t par_highest_rid, unsigned int ref_task_id, 
                   int insert_mode, rid_t *route, int cur_route, 
                   PrDep **ret_pr_dep) {

  Context       *context;
  PrDep         *pr_dep;
  ListNode      *list_node;
  ListNode      *new_list_node;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(dep_queue);
  ar_assert(final_dep);
  ar_assert(task_id);

  // Create new structure
  pr_dep = kt_malloc(sizeof(PrDep));
  pr_dep->task_id = task_id;
  pr_dep->par_task_id = par_task_id;
  pr_dep->par_highest_rid = par_highest_rid;
  pr_dep->final_target = final_dep;
  pr_dep->flags = dep_flags;
  pr_dep->flags |= SYS_TYPE_DEP_PENDING; // We force it to 1, because if
                                         // we are called by a request, this
                                         // flag is not encoded in it (although
                                         // it is set in the scheduler
                                         // responsible for the task).
  if (route) {
    ar_assert(cur_route + 1 < PR_REQ_MAX_SIZE);
    for (i = 0; i < cur_route + 1; i++) {
      pr_dep->route[i] = route[i];
    }
    pr_dep->cur_route = cur_route;
  }
  else {
    pr_dep->route[0] = 0;
    pr_dep->cur_route = 0;
  }
  if (ret_pr_dep) {
    *ret_pr_dep = pr_dep;
  }


  // How do we insert the new task?
  switch (insert_mode) {
    
    // At the tail
    case 0:
      new_list_node = kt_list_insert(dep_queue, pr_dep, 
                                     kt_list_tail(dep_queue), 1);
      break;

    // Before or after ref_task_id
    case 1:
    case 2:

      // Locate ref_task_id
      ar_assert(ref_task_id);
      for (list_node = kt_list_head(dep_queue);
           list_node;
           list_node = kt_list_next(list_node)) {
        ar_assert(list_node->data);
        if (((PrDep *) list_node->data)->task_id == ref_task_id) {
          break;
        }
      }
      ar_assert(list_node); // Sanity check: we must have found ref_task_id

      if (insert_mode == 1) {
        // Insert before ref_task_id
        new_list_node = kt_list_insert(dep_queue, pr_dep, list_node, 0);
      }
      else {
        // Insert after ref_task_id
        new_list_node = kt_list_insert(dep_queue, pr_dep, list_node, 1);
      }
      break;

    // Unknown mode
    default:
      ar_abort();
  }

  // Are we now at the head?
  if (new_list_node == kt_list_head(dep_queue)) {
    return 1;
  }
  else {
    return 0;
  }
}


// ===========================================================================
// pr_dep_search_queue()        Searches a dependency queue to find if a 
//                              given task ID is in there somewhere
// ===========================================================================
// * INPUTS
//   List *dep_queue            The dependency queue
//   unsigned int task_id       The task ID to locate
//
// * RETURN VALUE
//   ListNode *                 The located task ID dep node, or NULL if the
//                              task ID is not enqueued
// ===========================================================================
ListNode *pr_dep_search_queue(List *dep_queue, unsigned int task_id) {

  ListNode      *list_node;


  // Sanity checks
  ar_assert(dep_queue);
  ar_assert(task_id);

  // Locate task ID
  for (list_node = kt_list_head(dep_queue);
       list_node;
       list_node = kt_list_next(list_node)) {
    ar_assert(list_node->data);
    if (((PrDep *) list_node->data)->task_id == task_id) {
      return list_node;
    }
  }

  return NULL;
}


// ===========================================================================
// pr_dep_at_head_of_queue()    Sets a dependency to non-pending state and
//                              when all task dependencies are non-pending
//                              it makes the task ready to run. Note that
//                              this function may also be called for deps
//                              that are not exactly at the head of the queue,
//                              e.g. when a read-only dep is at the head and
//                              we are also looking ahead for more read-only
//                              deps.
// ===========================================================================
// * INPUTS
//   unsigned int task_id       The task ID of the task
//   void *dep                  The dependency to set non-pending
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int pr_dep_at_head_of_queue(unsigned int task_id, void *dep) {
  
  PrTaskDescr   *task;
  Context       *context;
  PrEvtHookPack *state;
  int           status;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(task_id);
  ar_assert(dep);

  // Locate task
  ar_assert(kt_trie_find(context->pr_tasks, task_id, (void *) &task));
  ar_assert(task->id == task_id);

  // Locate dependency and make it non-pending
  for (i = 0; i < task->num_args; i++) {
    if ((task->deps[i] != SYS_TYPE_BYVALUE_ARG) && (task->args[i] == dep)) {
//if (!(task->deps[i] & SYS_TYPE_DEP_PENDING)) {
//kt_printf("%d: at_head task 0x%08X index %d already non-pending arg %d\r\n", context->pr_core_id, task->id, task->index, i);
//for (i = 0; i < task->num_args; i++) {
//  kt_printf("arg %d = 0x%08X flags 0x%02X\r\n", i, task->args[i], task->deps[i]);
//}
//ar_abort();
//}
      ar_assert(task->deps[i] & SYS_TYPE_DEP_PENDING);
      task->deps[i] ^= SYS_TYPE_DEP_PENDING;
      break;
    }
  }
  ar_assert(i < task->num_args);

  ar_assert(task->deps_waiting);
  task->deps_waiting--;

//kt_printf("%d: Task 0x%08X has OK dependency 0x%08X; Remaining %d more\r\n", context->pr_core_id, task->id, dep, task->deps_waiting);

  // Was this the last dependency we were waiting for?
  if (!task->deps_waiting) {

//kt_printf("%d: task 0x%08X idx %d dep-free\r\n", context->pr_core_id, task->id, task->index);

    // Pack the task arguments
    state = NULL;

    status = mm_pack(NULL, task, &state);

    // If the packing involves other schedulers, nothing more to do here. 
    // Event(s) will have been set up to call the scheduler back when 
    // all reply(ies) are gathered.
    if (status == ERR_REPLY_POSTPONED) {
      return 0;
    }

    // TODO: If packing gives an error, this is most probably a user error.
    //       Handle this gracefully instead of crashing.
    ar_assert(!status);

    // Start scheduling this task
    pr_sched_start_scheduling(state);

    // Free the state, but don't free state->sizes and state->addresses:
    // pr_sched_start_scheduling() will prune & graft them on the task
    // descriptor.
    kt_free(state);
  }


  // Success
  return 0;
}

