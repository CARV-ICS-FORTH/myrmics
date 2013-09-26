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
// Abstract      : Memory allocator top-level interface. Manages all the
//                 initialization and provides the entry points for handling
//                 user requests.
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: mm.c,v $
// CVS revision  : $Revision: 1.21 $
// Last modified : $Date: 2013/04/09 14:54:13 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <syscall.h>


// ###########################################################################
// ###                                                                     ###
// ###                          INTERNAL FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// mm_sched_level_multiplier()  Returns the multiplier which must be applied
//                              to pages and region IDs requests, as indicated
//                              by our scheduler level
// ===========================================================================
// * RETURN VALUE
//   int                        Value is: pr_num_children ^ scheduler_level
// ===========================================================================
static inline int mm_sched_level_multiplier() {
  Context *context;
  int     mult;
  int     i;

  context = mm_get_context(ar_get_core_id());
  mult = 1;

  for (i = 0; i < context->pr_scheduler_level; i++) {
#ifdef MM_AGGRESSIVE_PAGE_REQUESTS
    mult *= context->pr_num_children;
#endif
  }

  return mult;
}


// ===========================================================================
// mm_ask_pages_and_postpone()  Asks for pages using distributed allocation,
//                              and puts the current request on hold until
//                              pages arrive. If a page mesage is already
//                              pending, the request is hooked up on the
//                              existing mesage.
// ===========================================================================
// * INPUTS
//   size_t size                Wanted size in bytes
//   PrMsgReq *req              Request we're trying to serve. Can be
//                              basic or extended.
//   int action                 Action to do after pages arrive
//   void *action_data          Event-specific data hook to attach to action
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_OUT_OF_MEMORY: distributed mem allocation
//                                 failed
// ===========================================================================
int mm_ask_pages_and_postpone(size_t size, PrMsgReq *req, int action,
                              void *action_data) {

  Context               *context;
  int                   msg_id;
  PrEvtPending          *event;
  PrEvtPending          *list;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(size);
  ar_assert(req);

  // Single-core hack mode?
  if (context->pr_num_cores == 1) {
    return ERR_OUT_OF_MEMORY;
  }

  // Do a new request for pages only when there's no pending one
  if (!context->pr_pages_msg_id) {

    // Get more pages. Piggyback the wanted size, in pages, along with
    // what's missing from our local free pool to get back to the maximum
    // size.
    msg_id = mm_distr_get_pages((size +
                                 MM_LOCAL_FREE_PAGES_MAX *
                                   mm_sched_level_multiplier() *
                                   MM_PAGE_SIZE -
                                 context->mm_free_num_slabs * MM_SLAB_SIZE +
                                 // Add MM_PAGE_SIZE - 1 to round up division
                                 MM_PAGE_SIZE - 1) / MM_PAGE_SIZE);
    // Success
    if (msg_id > 0) {
      context->pr_pages_msg_id = msg_id;
    }

    // If distributed page allocation fails, we're definitely out of memory.
    else if (msg_id == ERR_OUT_OF_MEMORY) {
      return ERR_OUT_OF_MEMORY;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }


  // Defer this call to later on, when the pages have arrived. It'll have
  // to be redone from scratch.
  ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
  ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

  // Copy the whole request
  kt_memcpy(event->req, req, sizeof(PrMsgReq));
  event->action = action;
  event->prev = NULL;
  event->next = NULL;
  event->data = action_data;


  // Try to store event, to be woken up when more pages have arrived
  if (kt_trie_insert(context->pr_pending_events, context->pr_pages_msg_id,
                     event)) {

    // If we failed, other events are there. Get the first one in the list.
    kt_trie_find(context->pr_pending_events,
                 context->pr_pages_msg_id, (void *) &list);
    ar_assert(list);

    // Walk list till the tail
    ar_assert(!list->prev);
    while (list->next) {
      list = list->next;
    }

    // Hook us up on the tail
    list->next = event;
    event->prev = list;
  }

  // Success
  return 0;
}


// ===========================================================================
// mm_ask_rids_and_postpone()   Asks for region IDs using distributed
//                              allocation, and puts the current request on
//                              hold until they arrive. If a region IDs
//                              allocation mesage is already pending, the
//                              request is hooked up on the existing mesage.
// ===========================================================================
// * INPUTS
//   int num_rids               Wanted number of region IDs
//   PrMsgReq *req              Request we're trying to serve. Can be
//                              basic or extended.
//   int action                 Action to do after region IDs arrive
//   void *action_data          Event-specific data hook to attach to action
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_OUT_OF_RIDS: distributed allocation failed
// ===========================================================================
int mm_ask_rids_and_postpone(int num_rids, PrMsgReq *req, int action,
                             void *action_data) {

  Context               *context;
  int                   msg_id;
  PrEvtPending          *event;
  PrEvtPending          *list;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(num_rids);
  ar_assert(req);

  // Single-core hack mode?
  if (context->pr_num_cores == 1) {
    return ERR_OUT_OF_RIDS;
  }

  // Do a new request for region IDs only when there's no pending one
  if (!context->pr_rids_msg_id) {


    // Get more IDs. Piggyback the number of IDs we were asked to give,
    // plus what's missing from our local pool to get back to the maximum
    // amount.
    msg_id = mm_distr_get_rids(num_rids +
                                MM_LOCAL_FREE_RIDS_MAX *
                                  mm_sched_level_multiplier() -
                                context->mm_free_num_rids);

    // Success
    if (msg_id > 0) {
      context->pr_rids_msg_id = msg_id;
    }

    // If distributed region IDs fails, we're definitely out of IDs
    else if (msg_id == ERR_OUT_OF_RIDS) {
      return ERR_OUT_OF_RIDS;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }


  // Defer this call to later on, when the IDs have arrived. It'll have
  // to be redone from scratch.
  ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
  ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

  // Copy the whole request
  kt_memcpy(event->req, req, sizeof(PrMsgReq));
  event->action = action;
  event->prev = NULL;
  event->next = NULL;
  event->data = action_data;


  // Try to store event, to be woken up when more region IDs have arrived
  if (kt_trie_insert(context->pr_pending_events, context->pr_rids_msg_id,
                     event)) {

    // If we failed, other events are there. Get the first one in the list.
    kt_trie_find(context->pr_pending_events,
                 context->pr_rids_msg_id, (void *) &list);
    ar_assert(list);

    // Walk list till the tail
    ar_assert(!list->prev);
    while (list->next) {
      list = list->next;
    }

    // Hook us up on the tail
    list->next = event;
    event->prev = list;
  }

  // Success
  return 0;
}


// ###########################################################################
// ###                                                                     ###
// ###                          EXPORTED FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################


// ===========================================================================
// mm_user_init()               Initializes the user memory pool
// ===========================================================================
// * INPUTS
//   PrRole role                Core role (ROLE_SCHEDULER or ROLE_WORKER)
// ===========================================================================
void mm_user_init(PrRole role) {

  Context *context;
  int     my;
  int     theirs;
  int     i;

  context = mm_get_context(ar_get_core_id());

  // Worker
  if (role == ROLE_WORKER) {
    // In the future, request from the hardware access to all user pages. For
    // the moment there are no virtual/physical addresses complications, so
    // there is nothing to be done here. Both ARM and Formic arch-specific init
    // has already installed page tables/ART entries that enable access to the
    // user space.
  }


  // Scheduler: Initialize user memory regions metadata (we keep these in
  //            the kernel data pool).
  else if (role == ROLE_SCHEDULER) {

    // Top-level scheduler?
    if (context->pr_parent_sched_id == -1) {

      // Initialize regions and get all memory and region IDs
      ar_assert(!mm_region_init(1));
    }
    else {

      // Initialize regions, but don't claim memory or region IDs
      ar_assert(!mm_region_init(0));

      // Request some free pages from our parent, depending on our level and
      // number of children schedulers. We'll ask enough so that after all
      // schedulers have asked their parents, each one remains with its ideal
      // number of pages according to its level.
      my = 1;
      theirs = 0;
#ifdef MM_AGGRESSIVE_PAGE_REQUESTS
      for (i = 0; i < context->pr_scheduler_level; i++) {
        // How many we'll need to give immediately to our children schedulers
        theirs = context->pr_num_children * (my + theirs);
        // How many we'll need to keep for ourselves
        my *= context->pr_num_children;
      }
#endif

      // Ask for pages and keep the message ID as pending, so that further
      // memory requests are not processed
      context->pr_pages_msg_id = mm_distr_get_pages(
                                     MM_LOCAL_FREE_PAGES_MAX * (my + theirs));
      ar_assert(context->pr_pages_msg_id > 0);

      // Ask for region IDs and keep the message ID as pending, so that further
      // region alloc requests are not processed
      context->pr_rids_msg_id = mm_distr_get_rids(
                                     MM_LOCAL_FREE_RIDS_MAX * (my + theirs));
      ar_assert(context->pr_rids_msg_id > 0);
    }

    // Create children load array
    ar_assert(context->mm_children_load = kt_malloc(context->pr_num_children *
                                                    sizeof(int)));
    for (i = 0; i < context->pr_num_children; i++) {
      context->mm_children_load[i] = 0;
    }
  }


  // Unknown role
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_get_pages()               Tries finding some free pages locally. If it
//                              fails, it tries to forward the request to an
//                              appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_GET_PAGES request. Useful fields:
//                                req->size: number of needed pages
//
// * OUTPUTS
//   int *ret_adr               If successful, base address of the requested
//                              number of free pages
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_OUT_OF_MEMORY: failure, no more memory in
//                                 the whole system
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_get_pages(PrMsgReq *req, size_t *ret_adr) {

  Context               *context;
  int                   num_pages;
  size_t                free_adr;
  int                   free_num_slabs;
  int                   ret;


  // Get needed fields
  num_pages = req->size;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(num_pages);
  ar_assert(ret_adr);


  // Try to find some local free pages to give to requestor
  //
  // FIXME: this does not ensure that range is page-aligned. We'll need that
  //        for the Formic hw. Also fix page-align check at insert_pages.
  ret = mm_region_get_adr_range(num_pages * MM_PAGE_SIZE / MM_SLAB_SIZE,
                                &free_adr, &free_num_slabs);

  // Success
  if (!ret) {

//kt_printf("%d: got %d pages locally\r\n", context->pr_core_id, num_pages);

    // We don't expect more pages than we asked to, since page size should
    // be aligned to MM_ADR_RANGE_CHUNK_{MIN,MAX}...
    ar_assert(free_num_slabs == num_pages * MM_PAGE_SIZE / MM_SLAB_SIZE);

    // Update used ranges to show that we gave this range to our next level
    // scheduler. Note that he may give it to his child, and we won't know it:
    // we only know enough to forward messages up to the next scheduler.
    mm_region_add_adr_range(context->mm_used_ranges, NULL,
                            pr_core_scheduler_id(req->core_id),
                            free_adr, free_num_slabs, NULL);
    context->mm_free_num_slabs -= free_num_slabs;
    ar_assert(context->mm_free_num_slabs >= 0);

    // Success
    *ret_adr = free_adr;
    return 0;
  }

  // No more local pages
  else if (ret == ERR_OUT_OF_MEMORY) {

    // Revert to distributed version
    ret = mm_ask_pages_and_postpone(num_pages * MM_PAGE_SIZE, req, PR_ACT_REDO,
                                    NULL);

    // Success?
    if (!ret) {
//kt_printf("%d: asked for %d pages from parent\r\n", context->pr_core_id, num_pages);
      return ERR_REPLY_POSTPONED;
    }

    // If distributed page allocation fails, we're definitely out of memory.
    else if (ret == ERR_OUT_OF_MEMORY) {
      return ERR_OUT_OF_MEMORY;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Unexpected return code
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_get_rids()                Tries finding a chunk of free region IDs
//                              locally. If it fails, it tries to forward the
//                              request to an appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_GET_RIDS request. Useful fields:
//                                req->size: number of needed region IDs
//
// * OUTPUTS
//   size_t *ret_base           If successful, starting reigon ID of the
//                              requested chunk of IDs
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_OUT_OF_RIDS: failure, no more region IDs in
//                                 the whole system
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_get_rids(PrMsgReq *req, size_t *ret_base) {

  Context               *context;
  int                   num_rids;
  size_t                free_base;
  int                   ret;


  // Get needed fields
  num_rids = req->size;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(num_rids);
  ar_assert(ret_base);

  // Try to find some local free IDs to give to requestor
  ret = mm_region_get_rid_range(num_rids, &free_base);

  // Success
  if (!ret) {

    // Update used ranges to show that we gave this range to our next level
    // scheduler. Note that he may give it to his child, and we won't know it:
    // we only know enough to forward messages up to the next scheduler.
    mm_region_add_rid_range(context->mm_used_rids,
                            pr_core_scheduler_id(req->core_id),
                            free_base, num_rids, NULL);
    context->mm_free_num_rids -= num_rids;
    ar_assert(context->mm_free_num_rids >= 0);

    // Success
    *ret_base = free_base;
    return 0;
  }

  // No more local IDs
  else if (ret == ERR_OUT_OF_RIDS) {

    // Revert to distributed version
    ret = mm_ask_rids_and_postpone(num_rids, req, PR_ACT_REDO, NULL);

    // Success?
    if (!ret) {
      return ERR_REPLY_POSTPONED;
    }

    // If distributed IDs allocation fails, we're definitely out of them.
    else if (ret == ERR_OUT_OF_RIDS) {
      return ERR_OUT_OF_RIDS;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Unexpected return code
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_insert_free_pages()       Inserts a number of incoming free pages to our
//                              free pool
// ===========================================================================
// * INPUTS
//   size_t adr                 Free chunk starting address
//   int num_pages              Number of pages in free chunk
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_insert_free_pages(size_t adr, int num_pages) {

  Context *context;
  int     num_slabs;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(adr);
  //ar_assert(!(adr & (MM_PAGE_SIZE - 1))); // FIXME - fix get_pages
  ar_assert(num_pages);

  // Convert to slabs
  num_slabs = num_pages * MM_PAGE_SIZE / MM_SLAB_SIZE;

  // Put them in the free pool
  mm_region_add_adr_range(context->mm_free_ranges, NULL, -1,
                          adr, num_slabs, NULL);
  ar_assert(context->mm_free_num_slabs >= 0);
  context->mm_free_num_slabs += num_slabs;

  // Success
  return 0;
}


// ===========================================================================
// mm_insert_free_rids()        Inserts a number of incoming free region IDs
//                              to our free pool
// ===========================================================================
// * INPUTS
//   size_t base                Free chunk starting address
//   int num_rids               Number of region IDs in free chunk
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_insert_free_rids(size_t base, int num_rids) {

  Context *context;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(base);
  ar_assert(num_rids);

  // Put them in the free pool
  mm_region_add_rid_range(context->mm_free_rids, -1,
                          base, num_rids, NULL);
  ar_assert(context->mm_free_num_rids >= 0);
  context->mm_free_num_rids += num_rids;

  // Success
  return 0;
}


// ===========================================================================
// mm_alloc()                   Tries allocating a new object locally. If it
//                              fails, it tries to forward the request to an
//                              appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_ALLOC request. Useful fields:
//                                req->size:   size of object in bytes
//                                req->region: region ID to allocate object in
//
// * OUTPUTS
//   void **ret_ptr             If successful, new allocated object
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_NO_SUCH_REGION: failure, this region ID
//                                 is invalid in the whole system
//                              ERR_OUT_OF_MEMORY: failure, no more memory in
//                                 the whole system
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_alloc(PrMsgReq *req, void **ret_ptr) {

  Context               *context;
  size_t                size;
  rid_t                 region;
  PrEvtPending          *event;
  int                   ret;
  int                   msg_id;


  // Get needed fields
  size = req->size;
  region = req->region;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(size);
  ar_assert(size < (1 << 31));
  ar_assert(!(size & (MM_ALLOC_ALIGN - 1)));
  ar_assert(region);
  ar_assert(ret_ptr);

  // Try to allocate object locally
  ret = mm_region_alloc_object(region, size, ret_ptr);

  // Success
  if (!ret) {
    return 0;
  }

  // Region not found locally
  else if (ret == ERR_NO_SUCH_REGION) {

    // Single-core hack mode?
    if (context->pr_num_cores == 1) {
      return ERR_NO_SUCH_REGION;
    }

    // Revert to distributed allocation
    msg_id = mm_distr_alloc(size, region);

    // Success
    if (msg_id > 0) {

      // Ok, but reply impossible at this point. Generate a pending event
      // to be handled when the reply for the distributed alloc is received.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Now it's definite, there's no such region.
    else if (msg_id == ERR_NO_SUCH_REGION) {
      return ERR_NO_SUCH_REGION;
    }

    // Unexpected return code
    else {
      ar_abort();
    }

    return ERR_REPLY_POSTPONED;
  }

  // Out of memory
  else if (ret == ERR_OUT_OF_MEMORY) {

    // Ask for more pages and put the request to wait on them
    ret = mm_ask_pages_and_postpone(size, req, PR_ACT_REDO, NULL);

    // Success?
    if (!ret) {
      return ERR_REPLY_POSTPONED;
    }

    // If distributed page allocation fails, we're definitely out of memory.
    else if (ret == ERR_OUT_OF_MEMORY) {
      return ERR_OUT_OF_MEMORY;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Unexpected return code
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_balloc()                  Tries allocating multiple new objects locally.
//                              If it fails, it tries to forward the request
//                              to an appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_BALLOC request. Useful fields:
//                                req->size:   size of each object in bytes
//                                req->region: region ID to allocate objects
//                                req->ptr:    number of objects
//
// * INOUTS
//   PrEvtHookBalloc **state    Reentrant state which keeps track of our
//                              progress. If *state is NULL, we allocate a new
//                              state. Final caller should free this when the
//                              final reply is created out of it.
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_NO_SUCH_REGION: failure, this region ID
//                                 is invalid in the whole system
//                              ERR_OUT_OF_MEMORY: failure, no more memory in
//                                 the whole system
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_balloc(PrMsgReq *req, PrEvtHookBalloc **state) {

  Context               *context;
  size_t                size;
  rid_t                 region;
  int                   num_elements;
  PrEvtPending          *event;
  int                   ret;
  int                   msg_id;
  int                   i;


  // Get needed fields
  size = req->size;
  region = req->region;
  num_elements = (size_t) req->ptr;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(size);
  ar_assert(size < (1 << 31));
  ar_assert(!(size & (MM_ALLOC_ALIGN - 1)));
  ar_assert(region);
  ar_assert(num_elements);
  ar_assert(state);


  // Is the region local?
  ret = mm_region_query_region(region, NULL);

  // Success, continue handling balloc locally
  if (!ret) {

    // If this is the first time we start handling this balloc, allocate state
    // to hold our results and mark our progress
    if (!(*state)) {
      ar_assert(*state = kt_malloc(sizeof(PrEvtHookBalloc)));
      ar_assert((*state)->objects = kt_malloc(num_elements * sizeof(void *)));
      (*state)->num_objects = 0;
    }

    // For all (remaining) elements
    for (i = 0; i < num_elements; i++) {

      // Allocate object
      ret = mm_region_alloc_object(region, size,
                                   (*state)->objects + (*state)->num_objects);

      if (!ret) {
        // Success; continue.
        (*state)->num_objects++;
      }
      else if (ret == ERR_NO_SUCH_REGION) {
        // Impossible: we checked that region is local.
        ar_abort();
      }
      else if (ret == ERR_OUT_OF_MEMORY) {

        // We ran out of memory. We need to stop here, ask for more memory and
        // resume allocating objects when more pages arrive. Fiddle original
        // request to be a balloc() for the remaining objects.
        req->ptr = (void *) ((size_t) (num_elements - (*state)->num_objects));

        // Ask for more pages and put the modified request to wait on them.
        // It's a reentrant request, with our current state attached.
        ret = mm_ask_pages_and_postpone(size * (size_t) req->ptr, req,
                                        PR_ACT_REENTRANT, *state);

        // Success?
        if (!ret) {
          return ERR_REPLY_POSTPONED;
        }

        // If distributed page allocation fails, we're definitely out of memory.
        else if (ret == ERR_OUT_OF_MEMORY) {
          return ERR_OUT_OF_MEMORY;
        }

        // Unexpected return code
        else {
          ar_abort();
        }

      }
      else {
        // Unknown return code
        ar_abort();
      }
    }

    // Success
    return 0;
  }

  // Region not local. Do not do multiple allocs() from here; instead,
  // pass a single balloc request to be done by the responsible scheduler.
  else if (ret == ERR_NO_SUCH_REGION) {
    ar_assert(!(*state));

    // Single-core hack mode? Non-local means non-existent.
    if (context->pr_num_cores == 1) {
      return ERR_NO_SUCH_REGION;
    }

    // Revert to distributed bulk allocation
    msg_id = mm_distr_balloc(size, region, num_elements);

    // Success
    if (msg_id > 0) {

      // Ok, but reply impossible at this point. Generate a pending event
      // to be handled when the reply for the distributed balloc is received.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Now it's definite, there's no such region.
    else if (msg_id == ERR_NO_SUCH_REGION) {
      return ERR_NO_SUCH_REGION;
    }

    // Unexpected return code
    else {
      ar_abort();
    }

    return ERR_REPLY_POSTPONED;
  }

  // Unexpected return code
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_free()                    Tries freeing an object locally. If it fails,
//                              it tries to forward the request to an
//                              appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_FREE request. Useful fields:
//                                req->ptr: object to be freed
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_MISALIGNED: failure, this pointer is bad
//                                 (not aligned properly to LSBs or slab)
//                              ERR_NOT_ALLOCED: failure, this object is
//                                 nowhere allocated in the whole system (and
//                                 specifically, its coarse address range is
//                                 correct, but the object is not claimed by
//                                 the slab system)
//                              ERR_OUT_OF_RANGE: failure, this object is
//                                 nowhere allocated in the whole system (and
//                                 specifically, its coarse address range is
//                                 invalid)
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_free(PrMsgReq *req) {

  Context               *context;
  void                  *ptr;
  PrEvtPending          *event;
  int                   ret;
  int                   msg_id;


  // Get needed fields
  ptr = req->ptr;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(ptr);

  // Try to free it
  ret = mm_region_free_object(ptr);

  if (!ret) {
    // Success
    return 0;
  }

  // Object not found locally
  else if (ret == ERR_OUT_OF_RANGE) {

    // Single-core hack mode?
    if (context->pr_num_cores == 1) {
      return ERR_OUT_OF_RANGE;
    }

    // Revert to distributed allocation
    msg_id = mm_distr_free(ptr);

    // Success
    if (msg_id > 0) {

      // Ok, but reply impossible at this point. Generate a pending event
      // to be handled when the reply for the distributed alloc is received.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Now it's definite, there's no such pointer
    else if (msg_id == ERR_OUT_OF_RANGE) {
      return ERR_OUT_OF_RANGE;
    }

    // Unexpected return code
    else {
      ar_abort();
    }

    return ERR_REPLY_POSTPONED;
  }

  else if (ret == ERR_MISALIGNED) {
    // Bad pointer
    return ERR_MISALIGNED;
  }
  else if (ret == ERR_NOT_ALLOCED) {
    // Local range, but non-allocated pointer
    return ERR_NOT_ALLOCED;
  }

  // Unexpected return code
  ar_abort();
}


// ===========================================================================
// mm_ralloc()                  Tries allocating a new region locally. If it
//                              fails, it tries to forward the request to an
//                              appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_RALLOC request. Useful fields:
//                                req->region: parent region ID for new region
//                                req->ptr:    level hint for new region
//
// * OUTPUTS
//   rid_t *ret_region          If successful, new region ID
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_NO_SUCH_REGION: failure, this parent
//                                 region ID is invalid in the whole system
//                              ERR_OUT_OF_MEMORY: failure, no more memory in
//                                 the whole system
//                              ERR_OUT_OF_RIDS: failure, no more region IDs
//                                 in the whole system
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_ralloc(PrMsgReq *req, rid_t *ret_region) {

  Context               *context;
  rid_t                 parent;
  int                   src_core_id;
  int                   level_hint;
  PrEvtPending          *event;
  int                   ret;
  int                   msg_id;


  // Get needed fields
  parent = req->region;
  src_core_id = req->core_id;
  level_hint = (size_t) req->ptr;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(parent);
  ar_assert(ret_region);

  // Try to create region locally
  ret = mm_region_create_region(parent, 0, -1, level_hint, ret_region);

  // Success
  if (!ret) {

    // Update our region load, and check whether we need to report upstream
    context->mm_current_load++;
    pr_sched_report_load();

    // Success
    return 0;
  }

  // Parent region not found locally
  else if (ret == ERR_NO_SUCH_REGION) {

    // Single-core hack mode?
    if (context->pr_num_cores == 1) {
      return ERR_NO_SUCH_REGION;
    }

    // Revert to distributed allocation
    msg_id = mm_distr_ralloc(parent, level_hint);

    // Success
    if (msg_id > 0) {

      // Ok, but reply impossible at this point. Generate a pending event
      // to be handled when the reply for the distributed ralloc is received.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Now it's definite, there's no such region.
    else if (msg_id == ERR_NO_SUCH_REGION) {
      return ERR_NO_SUCH_REGION;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Child creation rejected; child is on the scheduler boundary
  else if (ret == ERR_BOUNDARY_REGION) {

    // Select one of our child schedulers and tell it to create the child
    // region there, using the parent ID. Note that this scheduler handles
    // this as an "orphan" child, i.e. a child with a remote parent ID.
    msg_id = mm_distr_ralloc_orphan(parent, src_core_id);

    // Success
    if (msg_id > 0) {

      // When the child has been created, we need to update the parent locally
      // with the newly arrived child region ID. Create a note-to-self
      // request here, to be woken up when the reply arrives.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event processor knows where to send
      // the final reply to, as well as the parent region ID, so we know
      // which region to update.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->req->region  = parent;

      // Create the rest of the note to self
      event->req->type = SELF_RALLOC_UPDATE_PARENT;
      event->action = PR_ACT_REDO;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Unexpected return code
    else {
      ar_abort();
    }

  }

  // Out of memory
  else if (ret == ERR_OUT_OF_MEMORY) {

    // Revert to distributed version
    ret = mm_ask_pages_and_postpone(MM_PAGE_SIZE, req,
                                    PR_ACT_REDO, NULL);

    // Success?
    if (!ret) {
      return ERR_REPLY_POSTPONED;
    }

    // If distributed page allocation fails, we're definitely out of memory.
    else if (ret == ERR_OUT_OF_MEMORY) {
      return ERR_OUT_OF_MEMORY;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Out of region IDs
  else if (ret == ERR_OUT_OF_RIDS) {

    // Revert to distributed version
    ret = mm_ask_rids_and_postpone(1, req, PR_ACT_REDO, NULL);

    // Success?
    if (!ret) {
      return ERR_REPLY_POSTPONED;
    }

    // If distributed IDs allocation fails, we're definitely out of them.
    else if (ret == ERR_OUT_OF_RIDS) {
      return ERR_OUT_OF_RIDS;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Unexpected return code
  ar_abort();
}


// ===========================================================================
// mm_ralloc_orphan()           Allocates a new region locally, with a
//                              remote parent region ID
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_RALLOC_ORPHAN request. Useful fields:
//                                req->region: parent region ID
//                                req->size:   parent location core ID
//
// * OUTPUTS
//   rid_t *ret_region          If successful, new region ID
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_OUT_OF_MEMORY: failure, no more memory in
//                                 the whole system
//                              ERR_OUT_OF_RIDS: failure, no more region IDs
//                                 in the whole system
//                              ERR_REPLY_POSTPONED: request paused while
//                                 waiting for more pages/region IDS. Pending
//                                 event created to handle replies.
// ===========================================================================
int mm_ralloc_orphan(PrMsgReq *req, rid_t *ret_region) {

  Context *context;
  rid_t   parent;
  int     parent_location;
  int     ret;


  // Get needed fields
  context = mm_get_context(ar_get_core_id());
  parent = req->region;
  parent_location = req->size;

  // Sanity checks
  ar_assert(parent);
  ar_assert(ret_region);

  // Allocate region locally with remote parent
  ret = mm_region_create_region(parent, 1, parent_location, -1, ret_region);

  // Success
  if (!ret) {

    // Update our region load, and check whether we need to report upstream
    context->mm_current_load++;
    pr_sched_report_load();

    // Success
    return 0;
  }

  // Parent region not found locally
  else if (ret == ERR_NO_SUCH_REGION) {
    // Impossible: the scheduler that told us to do so should be our parent,
    // and he should know that we have this region ID. This is a bug.
    ar_abort();
  }

  // Child creation rejected
  else if (ret == ERR_BOUNDARY_REGION) {
    // Also impossible: this should not happen for remote parents.
    ar_abort();
  }

  // Out of memory
  else if (ret == ERR_OUT_OF_MEMORY) {

    // Revert to distributed version
    ret = mm_ask_pages_and_postpone(MM_PAGE_SIZE, req, PR_ACT_REDO, NULL);

    // Success?
    if (!ret) {
      return ERR_REPLY_POSTPONED;
    }

    // If distributed page allocation fails, we're definitely out of memory.
    else if (ret == ERR_OUT_OF_MEMORY) {
      return ERR_OUT_OF_MEMORY;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Out of region IDs
  else if (ret == ERR_OUT_OF_RIDS) {

    // Revert to distributed version
    ret = mm_ask_rids_and_postpone(1, req, PR_ACT_REDO, NULL);

    // Success?
    if (!ret) {
      return ERR_REPLY_POSTPONED;
    }

    // If distributed IDs allocation fails, we're definitely out of them.
    else if (ret == ERR_OUT_OF_RIDS) {
      return ERR_OUT_OF_RIDS;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Unexpected return code
  ar_abort();
}


// ===========================================================================
// mm_ralloc_update_parent()    Updates a local parent that he has a new
//                              remote child region
// ===========================================================================
// * INPUTS
//   PrMsgReq *stored_req       SELF_RALLOC_UPDATE_PARENT request. Useful
//                              fields:
//                                stored_req->region: parent region ID
//   PrMsgReply *child_reply    REPLY_RALLOC reply from child scheduler.
//                              Useful fields:
//                                child_reply->result: new child region ID
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_NO_SUCH_REGION: failure, this parent
//                                 region ID is invalid locally
// ===========================================================================
int mm_ralloc_update_parent(PrMsgReq *stored_req,
                            PrMsgReply *child_reply) {

  rid_t    parent;
  rid_t    child;


  // Get needed fields
  parent = stored_req->region;
  child = (rid_t) child_reply->result;

  // Just wrap to region call
  return mm_region_update_parent(parent, child, 1);
}


// ===========================================================================
// mm_rfree()                   Tries freeing a region locally. If it fails,
//                              it tries to forward the request to an
//                              appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_RFREE request. Useful fields:
//                                req->region: region ID to be freed
//                                req->ptr:    if 1, create request to update
//                                             parent (if remote)
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_NO_SUCH_REGION: failure, this region ID
//                                 is invalid in the whole system
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_rfree(PrMsgReq *req) {

  Context               *context;
  rid_t                 region;
  int                   update_remote_parent;
  rid_t                 *err_regions;
  int                   num_err_regions;
  int                   num_frees;
  PrEvtHookRfree        *state;
  PrEvtPending          *event;
  int                   ret;
  int                   msg_id;
  int                   i;


  // Get needed fields
  region = req->region;
  update_remote_parent = (size_t) req->ptr;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(region);

  // Try to destroy region locally
  ret = mm_region_destroy_region(region, &err_regions, &num_err_regions,
                                 &num_frees);

  // Success?
  if (!ret) {

    // Update our region load, and check whether we need to report upstream
    context->mm_current_load -= num_frees;
    pr_sched_report_load();

    // Success
    return 0;
  }

  // Region not found locally
  else if (ret == ERR_NO_SUCH_REGION) {

    // Single-core hack mode?
    if (context->pr_num_cores == 1) {
      return ERR_NO_SUCH_REGION;
    }

    // Revert to distributed region freeing
    msg_id = mm_distr_rfree(region, update_remote_parent);

    // Success
    if (msg_id > 0) {

      // Ok, but reply impossible at this point. Generate a pending event
      // to be handled when the reply for the distributed ralloc is received.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Now it's definite, there's no such region.
    else if (msg_id == ERR_NO_SUCH_REGION) {
      return ERR_NO_SUCH_REGION;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  // Region has remote children, we have to free them first
  else if (ret == ERR_REMOTE_CHILDREN) {

    ar_assert(err_regions);
    ar_assert(num_err_regions);

    // Create a state variable to hold the number of pending children replies
    ar_assert(state = kt_malloc(sizeof(PrEvtHookRfree)));

    // We know how many messages we need to wait for
    state->wait_replies = num_err_regions;

    // For all children
    for (i = 0; i < num_err_regions; i++) {

      // Call appropriate scheduler to free this child. Since we own their
      // parent region, and we already have deleted these children from the
      // parent array, ask that we do not want to receive update_parent
      // messages. Also, we shouldn't encounter any errors: these region IDs
      // are stored by us, so they must definitely exist somewhere.
      msg_id = mm_distr_rfree(err_regions[i], 0);
      ar_assert(msg_id > 0);

      // When each reply arrives, we need to monitor it, so we can repeat the
      // local region freeing when _all_ of them have arrived. Create a
      // note-to-self request for each mesage to do that.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // Things that matter from the original request are the requestor core
      // ID and message ID, so the event processor knows where to send the
      // final reply to, as well as what's the region we are trying to free,
      // so the attempt can be repeated.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->req->region  = region;
      event->req->ptr     = (void *) ((size_t) update_remote_parent);

      // Create the rest of the note to self, attaching the current state.
      // We also want to handle errors ourselves, because we don't expect any
      // and we need to sanity-check this.
      event->req->type = SELF_RFREE_WAIT_CHILDREN;
      event->action = PR_ACT_REENTRANT_ERRORS;
      event->prev = NULL;
      event->next = NULL;
      event->data = state;

      // Store event; we don't expect conflicts on any of the message IDs
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));
    }

    // Free remote errors array
    kt_free(err_regions);

    // Success
    return ERR_REPLY_POSTPONED;
  }

  // Region freed successfully, but has a remote parent that needs to be
  // updated
  else if (ret == ERR_REMOTE_PARENT) {

    // Update our region load, and check whether we need to report upstream
    context->mm_current_load -= num_frees;
    pr_sched_report_load();

    // Sanity checks
    ar_assert(err_regions);
    ar_assert(num_err_regions == 1);

    // Should we really send a message to update the parent?
    if (update_remote_parent) {

      // Send request to the parent to do that
      msg_id = mm_distr_rfree_update_parent(*err_regions, region);
      ar_assert(msg_id > 0);

      // Free remote errors array
      kt_free(err_regions);

      // When parent is freed, just forward the reply
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // The request told us not to update the remote parent (because it's
    // coming from the scheduler that owns this parent and it's already taken
    // care of)
    else {

      // Free remote errors array
      kt_free(err_regions);

      // We're done
      return 0;
    }
  }

  // Unexpected return code
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_rfree_update_parent()     Updates a local region that it has lost a
//                              remote child
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_RFREE_UPDATE_PARENT request. Useful fields:
//                                req->region: parent region ID
//                                req->ptr:    remote deleted child ID
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_rfree_update_parent(PrMsgReq *req) {

  Context       *context;
  rid_t         parent;
  rid_t         child;


  // Get needed fields
  parent = req->region;
  child  = (rid_t) req->ptr;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(parent);
  ar_assert(child);

  // Update parent. It's our internal call, so this cannot fail.
  ar_assert(!mm_region_update_parent(parent, child, 0));

  // Success
  return 0;
}


// ===========================================================================
// mm_pack()                    Tries packing multiple regions and objects
//                              locally (i.e. create lists of starting
//                              addresses and sizes for them). For those
//                              regions/object that are not local, it tries to
//                              forward requests to appropriate schedulers.
//
//                              The function receives the regions/objects to
//                              pack in one of two ways:
//
//                              - If req != NULL, from the request structure
//                                (basic or extended). The task should be
//                                NULL in this case and is irrelevant. This
//                                mode is used when a packing is requested
//                                by an incoming message.
//
//                              - If req == NULL, then task != NULL and the
//                                task arguments are the regions/objects to
//                                be packed. This mode is used when a packing
//                                is requested natively.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_PACK or EXT_REQ_PACK request. For more
//                              details, see the .h file on how these are
//                              organized. Can be NULL only if task != NULL.
//   PrTaskDescr *task          A native task descriptor, used only when
//                              req == NULL.
//
// * OUTPUTS
//   PrEvtHookPack **state      Reentrant state which will keep track of our
//                              progress for the self pack merging events.
//                              Final caller should free this when the final
//                              replies are created out of it.
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_MISALIGNED: failure, an object pointer was
//                                 bad (not aligned properly to LSBs or slab)
//                              ERR_NOT_ALLOCED: failure, an object was
//                                 nowhere allocated in the whole system (and
//                                 specifically, its coarse address range is
//                                 correct, but the object is not claimed by
//                                 the slab system)
//                              ERR_OUT_OF_RANGE: failure, an object was
//                                 nowhere allocated in the whole system (and
//                                 specifically, its coarse address range is
//                                 invalid)
//                              ERR_NO_SUCH_REGION: failure, a region was
//                                 nowhere allocated in the whole system
//                              ERR_REPLY_POSTPONED: request(s) to other(s)
//                                 scheduler(s) were needed, pending event(s)
//                                 created to handle replies
// ===========================================================================
int mm_pack(PrMsgReq *req, PrTaskDescr *task, PrEvtHookPack **state) {

  Context               *context;
  PrEvtPending          *event;
  size_t                size;
  int                   ret;
  rid_t                 *err_children;
  int                   num_err_children;
  rid_t                 *remote_regions;
  int                   *remote_region_options;
  int                   num_remote_regions;
  void                  **remote_objects;
  int                   *remote_object_options;
  int                   num_remote_objects;
  rid_t                 cur_region;
  void                  *cur_object;
  int                   cur_object_location;
  rid_t                 cur_object_parent;
  MmRgnTreeNode         *cur_object_parent_node;
  int                   cur_options;
  size_t                new_object;
  size_t                new_size;
  int                   total;
  int                   *msg_ids;
  int                   num_messages;
  void                  *error_ptr;
  int                   i;
  int                   j;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(state);
  ar_assert(!(*state)); // Pack is not reentrant; the replies merge version is
  ar_assert((req && !task) || (!req && task)); // Req mode or task mode

  // Allocate state to hold our results
  ar_assert(*state = kt_malloc(sizeof(PrEvtHookPack)));
  (*state)->alloc_elements = 32;
  ar_assert((*state)->sizes = kt_malloc((*state)->alloc_elements *
                                        sizeof(int)));
  ar_assert((*state)->addresses = kt_malloc((*state)->alloc_elements *
                                            sizeof(size_t)));
  (*state)->num_elements = 0;
  (*state)->wait_replies = 0;
  (*state)->native_task_id = task ? task->id : 0;
  (*state)->error_ptr = NULL;
  (*state)->error_status = 0;

  // Initialize arrays for non-local regions/objects
  remote_regions = NULL;
  remote_region_options = NULL;
  num_remote_regions = 0;
  remote_objects = NULL;
  remote_object_options = NULL;
  num_remote_objects = 0;


  // For all regions we have to pack
  if (req) {
    total = ((req->region) ? 1 : 0) +
            ((req->type == EXT_REQ_PACK) ? req->num_regions : 0);
  }
  else {
    total = task->num_args; // that's the maximum, non-regions will be ignored
  }
  for (i = 0; i < total; i++) {

    // Determine which region we're processing
    if (req) {
      if (req->region) {
        if (i == 0) {
          cur_region  = req->region;
          cur_options = req->size & ((1 << MM_PACK_OPTION_BITS) - 1);
        }
        else {
          cur_region  = (rid_t) req->data[i - 1];
          cur_options = (req->size >> ((2 + i - 1) * MM_PACK_OPTION_BITS)) &
                        ((1 << MM_PACK_OPTION_BITS) - 1);
        }
      }
      else {
        cur_region  = (rid_t) req->data[i];
        cur_options = (req->size >> ((2 + i) * MM_PACK_OPTION_BITS)) &
                      ((1 << MM_PACK_OPTION_BITS) - 1);
      }
    }
    else {
      if (task->deps[i] & SYS_TYPE_REGION_ARG) {
        cur_region  = (rid_t) task->args[i];
        cur_options = (task->deps[i] & SYS_TYPE_OUT_ARG) ? MM_PACK_OPTION_RW :
                                                           MM_PACK_OPTION_RO;
      }
      else {
        // Object, ignore
        continue;
      }
    }

    // Pack region
    ret = mm_region_pack_region(cur_region, cur_options, 1,
                                &((*state)->alloc_elements),
                                &((*state)->num_elements),
                                &((*state)->addresses),
                                &((*state)->sizes),
                                &err_children, &num_err_children);

    // Was region itself remote?
    if (ret == ERR_NO_SUCH_REGION) {
      remote_regions = kt_realloc(remote_regions,
                                  (num_remote_regions + 1) * sizeof(rid_t));
      remote_region_options = kt_realloc(remote_region_options,
                                  (num_remote_regions + 1) * sizeof(int));
      remote_regions[num_remote_regions] = cur_region;
      remote_region_options[num_remote_regions] = cur_options;
      num_remote_regions++;
    }

    // Did region have remote children?
    else if (ret == ERR_REMOTE_CHILDREN) {
      ar_assert(err_children);
      ar_assert(num_err_children);
      remote_regions = kt_realloc(remote_regions,
                                  (num_remote_regions +
                                   num_err_children) * sizeof(rid_t));
      remote_region_options = kt_realloc(remote_region_options,
                                         (num_remote_regions +
                                          num_err_children) * sizeof(int));
      for (j = 0; j < num_err_children; j++) {
        remote_regions[num_remote_regions] = err_children[j];
        // children regions inherit the packing options of the parent
        remote_region_options[num_remote_regions] = cur_options;
        num_remote_regions++;
      }
      kt_free(err_children);
    }

    // Unknown error code
    else if (ret) {
      ar_abort();
    }
  }


  // For all objects we have to pack
  if (req) {
    total = ((req->ptr) ? 1 : 0) +
            ((req->type == EXT_REQ_PACK) ? req->num_ptrs : 0);
  }
  else {
    total = task->num_args; // that's the maximum, non-objects as well as safe,
                            // by-value objects will be ignored
  }
  for (i = 0; i < total; i++) {

    // Determine which object we're processing
    if (req) {
      if (req->ptr) {
        if (i == 0) {
          cur_object  = req->ptr;
          cur_options = (req->size >> MM_PACK_OPTION_BITS) &
                        ((1 << MM_PACK_OPTION_BITS) - 1);
        }
        else {
          cur_object  = req->data[req->num_regions + i - 1];
          cur_options = (req->size >> ((2 + req->num_regions + i - 1) *
                                                      MM_PACK_OPTION_BITS)) &
                        ((1 << MM_PACK_OPTION_BITS) - 1);
        }
      }
      else {
        cur_object  = req->data[req->num_regions + i];
        cur_options = (req->size >> ((2 + req->num_regions + i) *
                                                    MM_PACK_OPTION_BITS)) &
                      ((1 << MM_PACK_OPTION_BITS) - 1);
      }
    }
    else {
      if ((!(task->deps[i] & SYS_TYPE_REGION_ARG)) &&
          (task->deps[i] != SYS_TYPE_BYVALUE_ARG)) {
        cur_object  = task->args[i];
        cur_options = (task->deps[i] & SYS_TYPE_OUT_ARG) ? MM_PACK_OPTION_RW :
                                                           MM_PACK_OPTION_RO;
      }
      else {
        // Region or by-value object, ignore
        continue;
      }
    }

    // Query object
    ret = mm_region_query_pointer(cur_object, &size, &cur_object_parent);

    // Success
    if (!ret) {

      // Find local parent node
      ar_assert(kt_trie_find(context->mm_local_rids, cur_object_parent,
                             (void *) &cur_object_parent_node));

      // See if the object has an exceptional location, otherwise assign the
      // default location of the region
      if (cur_object_parent_node->obj_locations) {
        if (!kt_trie_find(cur_object_parent_node->obj_locations,
                          (size_t) cur_object,
                          (void *) &cur_object_location)) {
          cur_object_location = cur_object_parent_node->location;
        }
      }
      else {
        cur_object_location = cur_object_parent_node->location;
      }
      ar_assert(cur_object_location < (1 << MM_PACK_LOCATION_BITS));
      ar_assert(cur_options < (1 << MM_PACK_OPTION_BITS));

      // Store object, breaking it in chunks of maximum DMA size
      new_object = (size_t) cur_object;
      new_size = size;
      while (new_size) {

        if (new_size >= (1 << MM_PACK_SIZE_BITS)) {
          new_size = (1 << MM_PACK_SIZE_BITS) - AR_DMA_ALIGN;
        }

        // Reallocate result arrays, if needed, and store the new adr and size
        if ((*state)->num_elements >= (*state)->alloc_elements) {
          (*state)->alloc_elements *= 2;
          ar_assert((*state)->sizes = kt_realloc((*state)->sizes,
                                   (*state)->alloc_elements * sizeof(int)));
          ar_assert((*state)->addresses = kt_realloc((*state)->addresses,
                                   (*state)->alloc_elements * sizeof(size_t)));
        }
        ((*state)->addresses)[(*state)->num_elements] = new_object;
        ((*state)->sizes)[(*state)->num_elements] =
                  new_size |
                  (cur_object_location << MM_PACK_SIZE_BITS) |
                  (cur_options << (MM_PACK_SIZE_BITS + MM_PACK_LOCATION_BITS));
        ((*state)->num_elements)++;

        // Proceed with next chunk
        new_object += new_size;
        size -= new_size;
        new_size = size;
      }
    }

    // Object is non-local
    else if (ret == ERR_OUT_OF_RANGE) {
      remote_objects = kt_realloc(remote_objects,
                                  (num_remote_objects + 1) * sizeof(void *));
      remote_object_options = kt_realloc(remote_object_options,
                                  (num_remote_objects + 1) * sizeof(int));
      remote_objects[num_remote_objects] = cur_object;
      remote_object_options[num_remote_objects] = cur_options;
      num_remote_objects++;
    }

    // Bad pointer
    else if (ret == ERR_MISALIGNED) {
      (*state)->error_ptr = cur_object;
      (*state)->error_status = ERR_MISALIGNED;
      return ERR_MISALIGNED;
    }

    // Local range, but non-allocated pointer
    else if (ret == ERR_NOT_ALLOCED) {
      (*state)->error_ptr = cur_object;
      (*state)->error_status = ERR_NOT_ALLOCED;
      return ERR_NOT_ALLOCED;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }


  // Did we manage to finish without any pending remote objects/regions?
  if (!num_remote_objects && !num_remote_regions) {
    // Success
    return 0;
  }


  // Switch into distributed version: ask other scheduler(s) to pack any
  // remote objects & regions we gathered. Because multiple schedulers might
  // be involved, the function returns multiple message IDs.
  ret = mm_distr_pack(remote_regions, remote_region_options, num_remote_regions,
                      remote_objects, remote_object_options, num_remote_objects,
                      &msg_ids, &num_messages, &error_ptr);

  // Free helper arrays
  kt_free(remote_regions);
  kt_free(remote_region_options);
  kt_free(remote_objects);
  kt_free(remote_object_options);

  // Success
  if (!ret) {

    // For all messages that were sent to other schedulers
    ar_assert(num_messages);
    ar_assert(msg_ids);
    for (i = 0; i < num_messages; i++) {

      // When each reply arrives, we need to merge it with all collected
      // addresses/sizes we have found locally. Create a note-to-self request
      // for each message to do that.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event processor knows where to send
      // the final reply to.
      if (req) {
        event->req->core_id = req->core_id;
        event->req->req_id  = req->req_id;
      }
      else {
        event->req->core_id = -1;
        event->req->req_id  = -1;
      }

      // Create the rest of the note to self, attaching the current state.
      // We also want to handle errors ourselves, because we'll need to wait
      // replies from multiple sources.
      event->req->type = SELF_PACK_MERGE;
      event->action = PR_ACT_REENTRANT_ERRORS;
      event->prev = NULL;
      event->next = NULL;
      event->data = *state;

      // Store event; we don't expect conflicts on any of the message IDs
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_ids[i], event));
    }

    // Mark that this many replies must be received before we can finally
    // reply to original requestor
    (*state)->wait_replies = num_messages;

    // Free the message array
    kt_free(msg_ids);

    // Success
    return ERR_REPLY_POSTPONED;
  }

  // One of the regions failed; don't create any pending event. Any packing
  // replies that may come will be ignored.
  else if (ret == ERR_NO_SUCH_REGION) {
    (*state)->wait_replies = 0;
    (*state)->error_ptr = error_ptr;
    (*state)->error_status = ERR_NO_SUCH_REGION;
    kt_free(msg_ids);
    return ERR_NO_SUCH_REGION;
  }

  // One of the objects failed; don't create any pending event. Any packing
  // replies that may come will be ignored.
  else if (ret == ERR_OUT_OF_RANGE) {
    (*state)->wait_replies = 0;
    (*state)->error_ptr = error_ptr;
    (*state)->error_status = ERR_OUT_OF_RANGE;
    kt_free(msg_ids);
    return ERR_OUT_OF_RANGE;
  }

  // Unexpected return code
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_query_pointer()           Tries to get an object region ID and size. If
//                              it fails, it tries to forward the request to
//                              an appropriate scheduler.
// ===========================================================================
// * INPUTS
//   PrMsgReq *req              REQ_QUERY_POINTER request. Useful fields:
//                                req->ptr: object to be queried
//
// * OUTPUTS
//   size_t *ret_size           Returned size of object
//   rid_t *ret_rid             Returned region ID of object
//
// * RETURN VALUE
//   int                        0 for success
//                              ERR_MISALIGNED: failure, this pointer is bad
//                                 (not aligned properly to LSBs or slab)
//                              ERR_NOT_ALLOCED: failure, this object is
//                                 nowhere allocated in the whole system (and
//                                 specifically, its coarse address range is
//                                 correct, but the object is not claimed by
//                                 the slab system)
//                              ERR_OUT_OF_RANGE: failure, this object is
//                                 nowhere allocated in the whole system (and
//                                 specifically, its coarse address range is
//                                 invalid)
//                              ERR_REPLY_POSTPONED: request to another
//                                 scheduler was needed, pending event created
//                                 to handle its reply
// ===========================================================================
int mm_query_pointer(PrMsgReq *req, size_t *ret_size, rid_t *ret_rid) {

  Context               *context;
  void                  *ptr;
  PrEvtPending          *event;
  int                   ret;
  int                   msg_id;


  // Get needed fields
  ptr = req->ptr;

  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(ptr);
  ar_assert(ret_size);
  ar_assert(ret_rid);

  // Find out about this pointer
  ret = mm_region_query_pointer(ptr, ret_size, ret_rid);

  if (!ret) {
    // Success
    return 0;
  }

  // Object not found locally
  else if (ret == ERR_OUT_OF_RANGE) {

    // Single-core hack mode?
    if (context->pr_num_cores == 1) {
      return ERR_OUT_OF_RANGE;
    }

    // Revert to distributed allocation
    msg_id = mm_distr_query_pointer(ptr);

    // Success
    if (msg_id > 0) {

      // Ok, but reply impossible at this point. Generate a pending event
      // to be handled when the reply for the distributed query is received.
      ar_assert(event = kt_malloc(sizeof(PrEvtPending)));
      ar_assert(event->req = kt_malloc(sizeof(PrMsgReq)));

      // The only thing that matters from the original request is the requestor
      // core ID and message ID, so the event handler knows where to forward
      // the reply to.
      event->req->core_id = req->core_id;
      event->req->req_id  = req->req_id;
      event->action = PR_ACT_FORWARD;
      event->prev = NULL;
      event->next = NULL;
      event->data = NULL;

      // Store event; we don't expect conflicts on this message ID.
      ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));

      // Success
      return ERR_REPLY_POSTPONED;
    }

    // Now it's definite, there's no such pointer
    else if (msg_id == ERR_OUT_OF_RANGE) {
      return ERR_OUT_OF_RANGE;
    }

    // Unexpected return code
    else {
      ar_abort();
    }
  }

  else if (ret == ERR_MISALIGNED) {
    // Bad pointer
    return ERR_MISALIGNED;
  }
  else if (ret == ERR_NOT_ALLOCED) {
    // Local range, but non-allocated pointer
    return ERR_NOT_ALLOCED;
  }

  // Unexpected return code
  ar_abort();
}
