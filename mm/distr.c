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
// Author        : Spyros LYBERIS
// Abstract      : Distributed nature of memory allocation; free pages/region
//                 IDs trading, as well as remote regions/objects handling
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: distr.c,v $
// CVS revision  : $Revision: 1.7 $
// Last modified : $Date: 2012/11/07 16:17:03 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <noc.h>


// ===========================================================================
// mm_distr_region_sched_id()   Finds out to which scheduler is responsible
//                              for a given region. This is approximate: it 
//                              indicates the direction to be searched
//                              (one level up or one level down). A query
//                              must be made at that level to learn more.
// ===========================================================================
// * INPUTS
//   rid_t region               The region ID we're searching for
//
// * RETURN VALUE
//   int                        scheduler_id on success
//                              ERR_NO_SUCH_REGION: region does not exist
//                                 anywhere in the system
// ===========================================================================
int mm_distr_region_sched_id(rid_t region) {

  Context       *context;
  MmRgnRidRange *rid_range;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(region);
  ar_assert(context->mm_used_rids);

  // See if region belongs to our scheduler subtree. We have to check that,
  // because direction of search initially begins upwards (from child
  // scheduler to parent scheduler) but later on, if it passes the top-level,
  // it becomes downwards (from parent scheduler to child scheduler) until it
  // reaches the correct scheduler.
  kt_trie_find_approx(context->mm_used_rids, 0, region, (void *) &rid_range);
  if (rid_range) {
    // We found something, but does it contain this region ID?
    ar_assert(region >= rid_range->rid);
    if (region >= rid_range->rid + rid_range->num_rids) {
      rid_range = NULL;
    }
  }

  // We know nothing about it
  if (!rid_range) {

    // Are we the top-level scheduler?
    if (context->pr_parent_sched_id == -1) {
      // No such region exists anywhere in the system
      return ERR_NO_SUCH_REGION;
    }

    // Otherwise, our parent scheduler should know more
    else {
      return context->pr_parent_sched_id;
    }
  }

  // Region belongs to our subtree
  else {

    // Return appropriate child scheduler
    return rid_range->sched_id;
  }
}


// ===========================================================================
// mm_distr_object_sched_id()   Finds out which scheduler is responsible for
//                              a given object. This is approximate: it 
//                              indicates the direction to be searched
//                              (one level up or one level down). A query
//                              must be made at that level to learn more.
// ===========================================================================
// * INPUTS
//   void *ptr                  The object we're searching for
//
// * RETURN VALUE
//   int                        scheduler_id on success
//                              ERR_OUT_OF_RANGE: pointer does not exist
//                                 anywhere in the system
// ===========================================================================
int mm_distr_object_sched_id(void *ptr) {

  Context       *context;
  MmRgnAdrRange *range;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(ptr);
  ar_assert(context->mm_used_ranges);

  // See if pointer belongs to our scheduler subtree. We have to check that,
  // because direction of search initially begins upwards (from child
  // scheduler to parent scheduler) but later on, if it passes the top-level,
  // it becomes downwards (from parent scheduler to child scheduler) until it
  // reaches the correct scheduler.
  kt_trie_find_approx(context->mm_used_ranges, 0, (size_t) ptr, 
                      (void *) &range);
  if (range) {
    // We found something, but does it contain the pointer?
    ar_assert(range->address <= (size_t) ptr);
    if ((size_t) ptr >= range->address + range->num_slabs * MM_SLAB_SIZE) {
      range = NULL;
    }
  }

  // We know nothing about it
  if (!range) {

    // Are we the top-level scheduler?
    if (context->pr_parent_sched_id == -1) {
      // No such pointer exists anywhere in the system
      return ERR_OUT_OF_RANGE;
    }

    // Otherwise, our parent scheduler should know
    else {
      return context->pr_parent_sched_id;
    }
  }

  // Region belongs to our subtree
  else {

    // Return appropriate child scheduler
    return range->sched_id;
  }
}


// ===========================================================================
// mm_distr_get_pages()         Requests more free pages from our parent 
//                              scheduler
// ===========================================================================
// * INPUTS
//   int num_pages              Number of pages we need
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_OUT_OF_MEMORY: failure, no more memory
//                                 in the whole system
// ===========================================================================
int mm_distr_get_pages(int num_pages) {

  Context               *context;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(num_pages);

  // Are we the top-level scheduler?
  if (context->pr_parent_sched_id == -1) {

    // No more memory in the whole system
    return ERR_OUT_OF_MEMORY;
  }

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_GET_PAGES;
  new_req->size    = num_pages;
  
  // Send message to parent
  ar_assert(!noc_msg_send());
  
  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_get_rids()          Requests more free region IDs from our parent 
//                              scheduler
// ===========================================================================
// * INPUTS
//   int num_rids               Number of region IDs we need
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_OUT_OF_MEMORY: failure, no more memory
//                                 in the whole system
// ===========================================================================
int mm_distr_get_rids(int num_rids) {

  Context               *context;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(num_rids);

  // Are we the top-level scheduler?
  if (context->pr_parent_sched_id == -1) {

    // No more region IDs in the whole system
    return ERR_OUT_OF_RIDS;
  }

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_GET_RIDS;
  new_req->size    = num_rids;
  
  // Send message to parent
  ar_assert(!noc_msg_send());
  
  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_alloc()             Requests a new object allocation in a remote
//                              region from the appropriate scheduler
// ===========================================================================
// * INPUTS
//   size_t size                Size of the new object
//   rid_t region               Remote region ID of the new object
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_NO_SUCH_REGION: failure, we are sure
//                                 that this region ID does not exist in the
//                                 whole system
// ===========================================================================
int mm_distr_alloc(size_t size, rid_t region) {

  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(size);
  ar_assert(region);

  // Find which scheduler should we ask for this region
  sched_id = mm_distr_region_sched_id(region);
  if (sched_id == ERR_NO_SUCH_REGION) {

    // No such region exists anywhere in the system
    return ERR_NO_SUCH_REGION;
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_ALLOC;
  new_req->size    = size;
  new_req->region  = region;
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_balloc()            Requests a new bulk-objects allocation in a 
//                              remote region from the appropriate scheduler
// ===========================================================================
// * INPUTS
//   size_t size                Size of the new object
//   rid_t region               Remote region ID of the new object
//   int num_elements           Number of objects to bulk-allocate
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_NO_SUCH_REGION: failure, we are sure
//                                 that this region ID does not exist in the
//                                 whole system
// ===========================================================================
int mm_distr_balloc(size_t size, rid_t region, int num_elements) {
                    
  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(size);
  ar_assert(region);
  ar_assert(num_elements);

  // Find which scheduler should we ask for this region
  sched_id = mm_distr_region_sched_id(region);
  if (sched_id == ERR_NO_SUCH_REGION) {

    // No such region exists anywhere in the system
    return ERR_NO_SUCH_REGION;
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_BALLOC;
  new_req->size    = size;
  new_req->region  = region;
  new_req->ptr     = (void *) ((size_t) num_elements);
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_free()              Requests a new bulk-objects allocation in a 
//                              remote region from the appropriate scheduler
// ===========================================================================
// * INPUTS
//   void *ptr                  The remote object that must be freed
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_OUT_OF_RANGE: failure, we are sure
//                                 that this pointer is not handled by anyone
//                                 in the whole system
// ===========================================================================
int mm_distr_free(void *ptr) {

  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(ptr);

  // Find out which scheduler should we ask for this object
  sched_id = mm_distr_object_sched_id(ptr);
  if (sched_id == ERR_OUT_OF_RANGE) {

    // No such object exists anywhere in the system
    return ERR_OUT_OF_RANGE;
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_FREE;
  new_req->ptr     = ptr;
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_ralloc()            Requests a new region creation under a 
//                              remote parent region from the appropriate
//                              scheduler
// ===========================================================================
// * INPUTS
//   rid_t parent               Region ID of the remote parent region
//   int level_hint             Hint on the scheduler level that the user
//                              wants this region to be handled
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_NO_SUCH_REGION: failure, we are sure
//                                 that this parent region is not valid in
//                                 the whole system
// ===========================================================================
int mm_distr_ralloc(rid_t parent, int level_hint) {

  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(parent);

  // Find which scheduler should we ask for the parent region
  sched_id = mm_distr_region_sched_id(parent);
  if (sched_id == ERR_NO_SUCH_REGION) {

    // No such region exists anywhere in the system
    return ERR_NO_SUCH_REGION;
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_RALLOC;
  new_req->region  = parent;
  new_req->ptr     = (void *) ((size_t) level_hint);
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_ralloc_orphan()     Requests a new region creation under a 
//                              remote parent region from the appropriate
//                              scheduler, and tells him that it's an orphan
//                              for him, i.e. with a remote region ID.
// ===========================================================================
// * INPUTS
//   rid_t parent               Region ID of the remote parent region
//   int src_core_id            Core ID of the scheduler that did the request,
//                              so we can take it into account when deciding
//                              which scheduler should own the new orphan
//
// * RETURN VALUE
//   int                        Message ID of the new request
// ===========================================================================
int mm_distr_ralloc_orphan(rid_t parent, int src_core_id) {

  Context               *context;
  int                   src_sched_id;
  MmRgnTreeNode         *parent_node;
  int                   child_id;
  PrMsgReq              *new_req;
  int                   id;
  int                   i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(parent);
  ar_assert(context->pr_scheduler_level);
  ar_assert(context->pr_num_children);
  src_sched_id = pr_core_scheduler_id(src_core_id);
  ar_assert((src_sched_id >= 0) && 
            (src_sched_id < context->pr_num_schedulers));

  // The parent region should belong to us. Find the region node.
  ar_assert(kt_trie_find(context->mm_local_rids, parent, 
                         (void *) &parent_node));


  // In the stand-alone allocator version, we only load-balance if requests
  // come from "above" (i.e. the request from a worker reaches the top-level
  // scheduler and starts to descend). In the Myrmics version, there is 
  // no such case: a task is allowed to see only subtrees of objects/regions
  // it owns, so requests always come from "below". Thus, we always
  // load-balance.

#if 0
  // Request comes from one of our children
  if (src_sched_id != context->pr_parent_sched_id) {
    // Verify he's one of our children
    child_id = -1;
    for (i = 0; i < context->pr_num_children; i++) {
      if (context->pr_children[i] == src_core_id) {
        child_id = i;
        break;
      }
    }
    ar_assert(child_id > -1);
  }
  
  // Request comes from our parent
  else {
#endif  
    // Search for least load, favoring the lastly decided round robin
    // position, in case there's a draw
    child_id = context->mm_load_rrobin;

    for (i = 0; i < context->pr_num_children; i++) {
      if (context->mm_children_load[i] < context->mm_children_load[child_id]) {
        child_id = i;
      }
    }

    // New round robin position is the next child from the one we decided
    context->mm_load_rrobin = child_id + 1;
    if (context->mm_load_rrobin >= context->pr_num_children) {
      context->mm_load_rrobin = 0;
    }
#if 0
  }
#endif

  // Build message
  new_req = noc_msg_send_get_buf(context->pr_children[child_id]);

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_RALLOC_ORPHAN;
  new_req->region  = parent;
  new_req->size    = parent_node->location;
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_rfree()             Requests the freeing of a remote region
//                              from the appropriate scheduler
// ===========================================================================
// * INPUTS
//   rid_t region               Region ID of the remote region to be freed
//   int update_remote_parent   1: Tell the other scheduler to update the
//                                 deleted region's parent, if this parent
//                                 is remote to the other scheduler
//                              0: Tell the other scheduler not to update the
//                                 deleted region's parent if it's remote
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_NO_SUCH_REGION: failure, we are sure
//                                 that this region is not valid in the 
//                                 whole system
// ===========================================================================
int mm_distr_rfree(rid_t region, int update_remote_parent) {

  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(region);

  // Find which scheduler should we ask for this region
  sched_id = mm_distr_region_sched_id(region);
  if (sched_id == ERR_NO_SUCH_REGION) {

    // No such region exists anywhere in the system
    return ERR_NO_SUCH_REGION;
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_RFREE;
  new_req->region  = region;
  new_req->ptr     = (void *) ((size_t) update_remote_parent);
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_rfree_update_parent()  Requests the updating of a remote parent
//                                 as to the loss of one of its children
// ===========================================================================
// * INPUTS
//   rid_t parent               Remote region ID of parent
//   rid_t child                Region ID of the child that is deleted
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
// ===========================================================================
int mm_distr_rfree_update_parent(rid_t parent, rid_t child) {

  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(parent);
  ar_assert(child);

  // Find which scheduler should we ask for the parent
  sched_id = mm_distr_region_sched_id(parent);
  if (sched_id == ERR_NO_SUCH_REGION) {
    // Corruption, we should have a valid parent ID
    ar_abort();
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_RFREE_UPDATE_PARENT;
  new_req->region  = parent;
  new_req->ptr     = (void *) child;
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


// ===========================================================================
// mm_distr_pack()              Create multiple packing requests to
//                              appropriate schedulers for packing multiple
//                              objects and regions. 
//
//                              The regions and objects given to this function
//                              are checked one by one as to which scheduler
//                              should be responsible for them. Before any
//                              messages are sent, per-scheduler arrays are
//                              created to gather all these decisions. When
//                              this sorting process is finished, request(s)
//                              to each scheduler are made, as many as needed
//                              depending on the number of regions/objects
//                              that are needed by each scheduler.
//
//                              Because many messages to many schedulers are
//                              sent, an array of new message IDs is retured.
// ===========================================================================
// * INPUTS
//   rid_t *regions             Array of remote region IDs to be packed
//   int *region_options        Array of remote region packing options
//   int num_regions            Number of regions in array
//   void **objects             Array of remote objects to be packed
//   int *object_options        Array of remote object packing options
//   int num_objects            Number of objects in array
//   
// * OUTPUTS
//   int **ret_msg_ids          Array of returned message IDs, one for each
//                              message sent
//   int *ret_num_messages      Number of messages generated
//   void **ret_error_ptr       In case of error, region ID or object that
//                              triggered it
//
// * RETURN VALUE
//   int                        0: success, *ret_num_messages generated
//                              ERR_NO_SUCH_REGION: failure, *ret_error_ptr
//                                 region ID is invalid in the whole system.
//                                 No messages generated.
//                              ERR_OUT_OF_RANGE: failure, *ret_error_ptr
//                                 object is invalid in the whole system.
//                                 No messages generated.
// ===========================================================================
int mm_distr_pack(rid_t *regions, int *region_options, int num_regions, 
                  void **objects, int *object_options, int num_objects, 
                  int **ret_msg_ids, int *ret_num_messages,
                  void **ret_error_ptr) {

  typedef struct {
    rid_t *regions;
    int   *region_options;
    int   num_regions;
    void  **objects;
    int   *object_options;
    int   num_objects;
  } pack_per_sched_type;

  Context               *context;
  PrMsgReq              *req;
  Trie                  *trie;
  pack_per_sched_type   *per_sched;
  int                   sched_id;
  int                   error_status;
  int                   cur_reg;
  int                   cur_obj;
  int                   i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(ret_msg_ids);
  ar_assert(ret_num_messages);
  ar_assert(ret_error_ptr);

  // Initialize return variables
  *ret_msg_ids = NULL;
  *ret_num_messages = 0;
  *ret_error_ptr = NULL;

  // Create a trie to store per-scheduler packing needs. We'll be using IDs
  // from 1 to context->pr_num_schedulers, so the MSB is the log2 of
  // context->pr_num_schedulers. The kt_int_log2() function will return the
  // MSB position correctly, even for non-power-of-2 values.
  ar_assert(trie = kt_alloc_trie(kt_int_log2(context->pr_num_schedulers), 0));


  // For all regions
  for (i = 0; i < num_regions; i++) {

    // Find which scheduler should we ask for this region
    sched_id = mm_distr_region_sched_id(regions[i]);
    if (sched_id == ERR_NO_SUCH_REGION) {

      // No such region exists anywhere in the system. Abort.
      *ret_error_ptr = (void *) regions[i];
      error_status = ERR_NO_SUCH_REGION;
      goto error;
    }
    // It can't be ours, local functions should have picked it up
    ar_assert(sched_id != context->pr_scheduler_id);

    // We use scheduler ID + 1, because tries cannot handle zero keys
    sched_id++;

    // Do we have anything else for this scheduler ID?
    kt_trie_find(trie, sched_id, (void *) &per_sched);

    // Append this region to existing entry
    if (per_sched) {
      per_sched->regions = kt_realloc(per_sched->regions,
                                      (per_sched->num_regions + 1) *
                                                sizeof(rid_t));
      per_sched->region_options = kt_realloc(per_sched->region_options,
                                             (per_sched->num_regions + 1) *
                                                       sizeof(int));
      per_sched->regions[per_sched->num_regions] = regions[i];
      per_sched->region_options[per_sched->num_regions] = region_options[i];
      per_sched->num_regions++;
    }

    // Create new entry for this scheduler
    else {
      per_sched = kt_malloc(sizeof(pack_per_sched_type));
      per_sched->regions = kt_malloc(sizeof(rid_t));
      per_sched->region_options = kt_malloc(sizeof(int));
      per_sched->regions[0] = regions[i];
      per_sched->region_options[0] = region_options[i];
      per_sched->num_regions = 1;
      per_sched->objects = NULL;
      per_sched->object_options = NULL;
      per_sched->num_objects = 0;

      // Insert it into the trie
      ar_assert(!kt_trie_insert(trie, sched_id, per_sched));
    }
  }


  // For all objects
  for (i = 0; i < num_objects; i++) {

    // Find which scheduler should we ask for this object
    sched_id = mm_distr_object_sched_id(objects[i]);
    if (sched_id == ERR_OUT_OF_RANGE) {

      // No such object exists anywhere in the system. Abort.
      *ret_error_ptr = objects[i];
      error_status = ERR_OUT_OF_RANGE;
      goto error;
    }
    // It can't be ours, local functions should have picked it up
    ar_assert(sched_id != context->pr_scheduler_id);

    // We use scheduler ID + 1, because tries cannot handle zero keys
    sched_id++;

    // Do we have anything else for this scheduler ID?
    kt_trie_find(trie, sched_id, (void *) &per_sched);

    // Append this object to existing entry
    if (per_sched) {
      per_sched->objects = kt_realloc(per_sched->objects,
                                      (per_sched->num_objects + 1) *
                                                sizeof(void *));
      per_sched->object_options = kt_realloc(per_sched->object_options,
                                      (per_sched->num_objects + 1) *
                                                sizeof(int));
      per_sched->objects[per_sched->num_objects] = objects[i];
      per_sched->object_options[per_sched->num_objects] = object_options[i];
      per_sched->num_objects++;
    }

    // Create new entry for this scheduler
    else {
      per_sched = kt_malloc(sizeof(pack_per_sched_type));
      per_sched->objects = kt_malloc(sizeof(void *));
      per_sched->object_options = kt_malloc(sizeof(int));
      per_sched->objects[0] = objects[i];
      per_sched->object_options[0] = object_options[i];
      per_sched->num_objects = 1;
      per_sched->regions = NULL;
      per_sched->region_options = NULL;
      per_sched->num_regions = 0;

      // Insert it into the trie
      ar_assert(!kt_trie_insert(trie, sched_id, per_sched));
    }
  }


  // For all schedulers that we need to communicate with
  for (sched_id = kt_trie_find_minmax(trie, 0, (void *) &per_sched);
       sched_id;
       sched_id = kt_trie_find_next(trie, 1, (void *) &per_sched)) {

    // Restore correct scheduler ID value
    sched_id--;

    // Initialize
    cur_reg = 0;
    cur_obj = 0;
    req = NULL;

    // Loop until all regions and objects have been sent
    while ((cur_reg < per_sched->num_regions) || 
           (cur_obj < per_sched->num_objects)) {

      // Start building new message
      req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

      req->core_id = context->pr_core_id;
      req->req_id  = context->pr_message_id;
      req->size    = 0;

      // Append this message ID to the array of IDs we'll return
      ar_assert(*ret_msg_ids = kt_realloc(*ret_msg_ids, 
                                            (*ret_num_messages + 1) * 
                                            sizeof(int)));
      (*ret_msg_ids)[*ret_num_messages] = req->req_id;
      (*ret_num_messages)++;

      // Increase message ID, avoiding value 0 on wrap-arounds
      context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

      // Embed first remaining region and first remaining object to basic
      // request
      if (cur_reg < per_sched->num_regions) {
        ar_assert(per_sched->regions[cur_reg]);
        req->region = per_sched->regions[cur_reg];
        req->size |= per_sched->region_options[cur_reg] << 0;
        cur_reg++;
      }
      else {
        req->region = 0;
      }
      if (cur_obj < per_sched->num_objects) {
        ar_assert(per_sched->objects[cur_obj]);
        req->ptr = per_sched->objects[cur_obj];
        req->size |= per_sched->object_options[cur_obj] << MM_PACK_OPTION_BITS;
        cur_obj++;
      }
      else {
        req->ptr = NULL;
      }

      // Did we fit or do we need an extended request?
      if ((cur_reg >= per_sched->num_regions) && 
          (cur_obj >= per_sched->num_objects)) {
        req->type = REQ_PACK;

        // Send basic request only to scheduler
        ar_assert(!noc_msg_send());
        req = NULL;
      }
      else {
        req->type = EXT_REQ_PACK;

        // Put regions in first slots of extended array
        for (req->num_regions = 0;
             (req->num_regions < PR_REQ_MAX_SIZE) && 
             (cur_reg < per_sched->num_regions);
             req->num_regions++, cur_reg++) {
          ar_assert(per_sched->regions[cur_reg]);
          req->data[req->num_regions] = (void *) per_sched->regions[cur_reg];
          req->size |= per_sched->region_options[cur_reg] <<
                              ((2 + req->num_regions) * MM_PACK_OPTION_BITS);
        }

        // Put objects on following slots, if regions are finished
        for (req->num_ptrs = 0;
             (req->num_regions + req->num_ptrs < PR_REQ_MAX_SIZE) &&
             (cur_obj < per_sched->num_objects);
             req->num_ptrs++, cur_obj++) {
          ar_assert(per_sched->objects[cur_obj]);
          req->data[req->num_regions + req->num_ptrs] = 
                                                  per_sched->objects[cur_obj];
          req->size |= per_sched->object_options[cur_obj] <<
                ((2 + req->num_regions + req->num_ptrs) * MM_PACK_OPTION_BITS);
        }

        // Send extended request to scheduler
        ar_assert(!noc_msg_send());
        req = NULL;
      }
    }

    // Free this entry
    kt_free(per_sched->regions);
    kt_free(per_sched->region_options);
    kt_free(per_sched->objects);
    kt_free(per_sched->object_options);
    kt_free(per_sched);
  }

  // Free the trie
  kt_free_trie(trie, NULL);

  // Success
  return 0;


error:

  // Free all allocated per-scheduler entries
  for (sched_id = kt_trie_find_minmax(trie, 0, (void *) &per_sched);
       sched_id;
       sched_id = kt_trie_find_next(trie, 1, (void *) &per_sched)) {
    kt_free(per_sched->regions);
    kt_free(per_sched->region_options);
    kt_free(per_sched->objects);
    kt_free(per_sched->object_options);
    kt_free(per_sched);
  }

  // Free the trie
  kt_free_trie(trie, NULL);

  // Return error code; error_ptr is set above.
  return error_status;

}


// ===========================================================================
// mm_distr_query_pointer()     Requests a pointer query operation from a 
//                              remote region from the appropriate scheduler
// ===========================================================================
// * INPUTS
//   void *ptr                  The remote object that must be queried
//
// * RETURN VALUE
//   int                        Message ID of the new request on success
//                              ERR_OUT_OF_RANGE: failure, we are sure
//                                 that this pointer is not handled by anyone
//                                 in the whole system
// ===========================================================================
int mm_distr_query_pointer(void *ptr) {

  Context               *context;
  int                   sched_id;
  PrMsgReq              *new_req;
  int                   id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(ptr);

  // Find out which scheduler should we ask for this object
  sched_id = mm_distr_object_sched_id(ptr);
  if (sched_id == ERR_OUT_OF_RANGE) {

    // No such object exists anywhere in the system
    return ERR_OUT_OF_RANGE;
  }
  // It can't be ours, local functions should have picked it up
  ar_assert(sched_id != context->pr_scheduler_id);

  // Build message
  new_req = noc_msg_send_get_buf(pr_scheduler_core_id(sched_id));

  new_req->core_id = context->pr_core_id;
  new_req->req_id  = context->pr_message_id;
  new_req->type    = REQ_QUERY_POINTER;
  new_req->ptr     = ptr;
  
  // Send message to the selected scheduler
  ar_assert(!noc_msg_send());

  // Increase message ID
  id = context->pr_message_id;
  context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

  // Success
  return id;
}


