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
// Abstract      : Memory allocation system calls
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: alloc.c,v $
// CVS revision  : $Revision: 1.4 $
// Last modified : $Date: 2012/10/22 13:19:33 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <noc.h>
#include <syscall.h>
#include <debug.h>


// ===========================================================================
// _sys_alloc()                 Allocates a new object in a region.
//
//                              The system call aborts the execution upon any 
//                              error and prints a message.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   size_t size                The requested size, in bytes
//   rid_t region               Region ID of the region, or 0 for the NULL
//                              region
//
// * RETURN VALUE
//   void *                     Pointer to the new allocated object, or NULL
//                              if the size is 0.
// ===========================================================================
void *_sys_alloc(char *filename, int line_nr, size_t size, rid_t region) {

  Context               *context;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  PrMsgReq              single_core_req;
  PrMsgReply            single_core_reply;
  int                   wait_msg_id;
#if 0
  unsigned int          sched_sec;
  unsigned int          sched_usec;
#endif


  // No size?
  if (!size) {
    return NULL;
  }

  // Get context
  context = mm_get_context(ar_get_core_id());

  // Clamp requests up to 2-GB size. We won't support more, even for the
  // x86_64 port. It makes it easier to work internally with signed integers,
  // because error checking and assertions work way better.
  if (size >= (1 << 31)) {
    _sys_error(filename, line_nr, "Huge size (%lu) for sys_alloc", size);
  }

  // Align size request to nearest allowed size
  if (size & (MM_ALLOC_ALIGN - 1)) {
    size = (size & ~(MM_ALLOC_ALIGN - 1)) + MM_ALLOC_ALIGN;
  }

  // Fix NULL region: internally, it has id 1. User should not directly use
  // region ID 1, it's illegal (it's never returned by any of the system
  // calls).
  if (region == 1) {
    _sys_error(filename, line_nr, "Illegal region id 1");
  }
  if (!region) {
    region = 1;
  }

  // Select request/reply buffers
  if (context->pr_num_cores == 1) {
    req   = &single_core_req;
    reply = &single_core_reply;
  }
  else {
    req   = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    reply = NULL;
  }

  // Build message
  req->core_id = context->pr_core_id;
  req->req_id  = context->pr_message_id;
  req->type    = REQ_ALLOC;
  req->size    = size;
  req->region  = region;

#if 0
  // Track time
  ar_get_time(&sched_sec, &sched_usec);
#endif

  // Single-core fast tracking?
  if (context->pr_num_cores == 1) {
    reply->status = mm_alloc(req, &(reply->result));
  }
  // Multi-core normal path
  else {

    // Remember which message we're waiting for, send message to scheduler
    // and increase system message ID
    wait_msg_id = req->req_id;
    ar_assert(!noc_msg_send());
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Wait until we get the reply we want
    reply = pr_event_worker_inner_loop(0, REPLY_ALLOC, wait_msg_id, NULL);
    ar_assert(reply);
  }

#if 0
  // Track time
  _sys_track_sched_time(sched_sec, sched_usec);
  context->sys_sched_calls++;
#endif

  // See what happened
  switch (reply->status) {
    case 0:
      return reply->result;
    case ERR_NO_SUCH_REGION:
      _sys_error(filename, line_nr, "Invalid region id %lu", region);
    case ERR_OUT_OF_MEMORY:
      _sys_error(filename, line_nr, "Out of memory");
    default:
      // Unknown return code
      ar_abort();
  }
}


// ===========================================================================
// _sys_balloc()                Allocates multiple objects of the same size
//                              in a single call and returns multiple pointers
//                              back to the user.
//
//                              The system call aborts the execution upon any 
//                              error and prints a message.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   size_t size                The requested size, in bytes
//   rid_t region               Region ID of the region, or 0 for the NULL
//                              region
//   int num_objects            Number of objects to be allocated
//
// * OUTPUTS
//   void **objects             Array of pointers that were allocated. The 
//                              user must allocate this array before the call.
// ===========================================================================
void _sys_balloc(char *filename, int line_nr, size_t size, rid_t region, 
                 int num_objects, void **objects) {

  Context               *context;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  PrMsgReq              single_core_req;
  PrMsgReply            single_core_reply;
  int                   wait_msg_id;
#if 0
  unsigned int          sched_sec;
  unsigned int          sched_usec;
#endif
  int                   i;
  int                   j;


  // Sanity checks
  if (num_objects <= 1) {
    _sys_error(filename, line_nr,
               "Invalid number of objects (%d) for sys_balloc", num_objects);
  }
  if (!size) {
    _sys_error(filename, line_nr, "Invalid size (0) for sys_balloc");
  }
  if (!objects) {
    _sys_error(filename, line_nr, "NULL objects array in sys_balloc");
  }

  // Get context
  context = mm_get_context(ar_get_core_id());

  // Clamp requests up to 2-GB size. We won't support more, even for the
  // x86_64 port. It makes it easier to work internally with signed integers,
  // because error checking and assertions work way better.
  if (size >= (1 << 31)) {
    _sys_error(filename, line_nr, "Huge size (%lu) for sys_balloc", size);
  }

  // Align size request to nearest allowed size
  if (size & (MM_ALLOC_ALIGN - 1)) {
    size = (size & ~(MM_ALLOC_ALIGN - 1)) + MM_ALLOC_ALIGN;
  }

  // Fix NULL region: internally, it has id 1. User should not directly use
  // region ID 1, it's illegal (it's never returned by any of the system
  // calls).
  if (region == 1) {
    _sys_error(filename, line_nr, "Illegal region id 1");
  }
  if (!region) {
    region = 1;
  }

  
  // Select request/reply buffers
  if (context->pr_num_cores == 1) {
    req   = &single_core_req;
    reply = &single_core_reply;
  }
  else {
    req   = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    reply = NULL;
  }

  // Build message
  req->core_id = context->pr_core_id;
  req->req_id  = context->pr_message_id;
  req->type    = REQ_BALLOC;
  req->size    = size;
  req->region  = region;
  req->ptr     = (void *) ((size_t) num_objects);

#if 0
  // Track time
  ar_get_time(&sched_sec, &sched_usec);
#endif

  // Single-core fast tracking?
  if (context->pr_num_cores == 1) {
    for (i = 0; i < num_objects; i++) {
      reply->status = mm_alloc(req, objects + i);
      if (reply->status) {
        break;
      }
    }
  }
  // Multi-core normal path
  else {

    // Remember which message we're waiting for, send message to scheduler
    // and increase system message ID
    wait_msg_id = req->req_id;
    ar_assert(!noc_msg_send());
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Wait until we get all the replies we want
    i = 0;
    while (1) {

      // Get a new reply
      reply = pr_event_worker_inner_loop(0, EXT_REPLY_BALLOC, wait_msg_id, 
                                         NULL);
      ar_assert(reply);

      // Is there an error?
      if (reply->status) {
        break;
      }

      // Get all the addresses
      for (j = 0; 
           j < ((reply->num_elements == -1) ? PR_RPL_MAX_SIZE
                                            : reply->num_elements); 
           i++, j++) {
        ar_assert(i < num_objects);
        objects[i] = reply->adr_ar[j];
      }

      // If no more replies must come, stop
      if (reply->num_elements > -1) {
        break;
      }
    }
  }

#if 0
  // Track time
  _sys_track_sched_time(sched_sec, sched_usec);
  context->sys_sched_calls++;
#endif

  // See what happened
  switch (reply->status) {
    case 0:
      ar_assert(i == num_objects);
      break;
    case ERR_NO_SUCH_REGION:
      _sys_error(filename, line_nr, "Invalid region id %lu", region);
    case ERR_OUT_OF_MEMORY:
      _sys_error(filename, line_nr, "Out of memory");
    default:
      // Unknown return code
      ar_abort();
  }
}


// ===========================================================================
// _sys_free()                  Frees a previously allocated object.
//
//                              The system call aborts the execution upon any 
//                              error and prints a message.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   void *ptr                  Pointer to the object to be freed
// ===========================================================================
void _sys_free(char *filename, int line_nr, void *ptr) {

  Context               *context;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  PrMsgReq              single_core_req;
  PrMsgReply            single_core_reply;
  int                   wait_msg_id;
#if 0
  unsigned int          sched_sec;
  unsigned int          sched_usec;
#endif

  // Allow NULL frees
  if (!ptr) {
    return;
  }

  // Get context
  context = mm_get_context(ar_get_core_id());

  // Pointer must be aligned
  if ((size_t) ptr & (MM_ALLOC_ALIGN - 1)) {
    _sys_error(filename, line_nr, "Tried to free misaligned pointer %p", ptr);
  }

  // Select request/reply buffers
  if (context->pr_num_cores == 1) {
    req   = &single_core_req;
    reply = &single_core_reply;
  }
  else {
    req   = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    reply = NULL;
  }

  // Build message
  req->core_id = context->pr_core_id;
  req->req_id  = context->pr_message_id;
  req->type    = REQ_FREE;
  req->ptr     = ptr;

#if 0
  // Track time
  ar_get_time(&sched_sec, &sched_usec);
#endif

  // Single-core fast tracking?
  if (context->pr_num_cores == 1) {
    reply->status = mm_free(req);
  }
  // Multi-core normal path
  else {

    // Remember which message we're waiting for, send message to scheduler
    // and increase system message ID
    wait_msg_id = req->req_id;
    ar_assert(!noc_msg_send());
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Wait until we get the reply we want
    reply = pr_event_worker_inner_loop(0, REPLY_FREE, wait_msg_id, NULL);
    ar_assert(reply);
  }

#if 0
  // Track time
  _sys_track_sched_time(sched_sec, sched_usec);
  context->_sys_sched_calls++;
#endif

  // See what happened
  switch (reply->status) {
    case 0:
      // Success
      return;
    case ERR_OUT_OF_RANGE:
      _sys_error(filename, line_nr,
                 "Tried to free out of range pointer %p", ptr);
    case ERR_MISALIGNED:
      _sys_error(filename, line_nr,
                 "Tried to free misaligned pointer %p", ptr);
    case ERR_NOT_ALLOCED:
      _sys_error(filename, line_nr,
                 "Tried to free non-allocated pointer %p", ptr);
    default:
      // Unknown return code
      ar_abort();
  }
}


// ===========================================================================
// _sys_realloc()               Reallocates an existing object to a new size
//                              and/or to a new region.
//
//                              The system call aborts the execution upon any 
//                              error and prints a message.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   void *old_ptr              The old object to be reallocated
//   size_t new_size            The new requested size, in bytes. The size
//                              can be the same as the old size (only region
//                              reallocation is done).
//   rid_t new_region           Region ID of the new region, or 0 for the NULL
//                              region. The region ID can be the same as the
//                              old region ID (only size adjustment is done).
//
// * RETURN VALUE
//   void *                     Pointer to the reallocated object, or NULL
//                              if the size is 0.
// ===========================================================================
void *_sys_realloc(char *filename, int line_nr, void *old_ptr, size_t new_size, 
                   rid_t new_region) {

  Context               *context;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  PrMsgReq              single_core_req;
  PrMsgReply            single_core_reply;
  int                   wait_msg_id;
#if 0
  unsigned int          sched_sec;
  unsigned int          sched_usec;
#endif
  void                  *new_ptr;
  size_t                old_size;
  rid_t                 old_region;


  // Handle degenerate cases here
  if (!old_ptr) {
    return _sys_alloc(filename, line_nr, new_size, new_region);
  }
  if (!new_size) {
    _sys_free(filename, line_nr, old_ptr);
    return NULL;
  }

  // NOTE: This implementation is non-optimal: it first does a new alloc
  //       call to scheduler, than a second call to find out the old pointer
  //       size, then copies the data and finally a third call to free the old
  //       pointer. This seems bad at first glance, but remember that we use
  //       an underlying slab allocator; it's very rare that it'll return the
  //       same pointer for a different slot size (it must be a > slab size
  //       object and manage to find adjacent free slabs to extend it). So,
  //       the really redundant stuff is the packing call.
  //
  //       Trying to optimize it by doing a single REALLOC call to the
  //       scheduler will involve several challenges. First and foremost, the
  //       scheduler itself is in another core and address space: so, in the
  //       common case that the old/new pointers are different, this will need
  //       to revert mostly to this version, where the data is copied by the
  //       requesting worker kernel space and a second message will be sent
  //       back to the scheduler to free the old pointer.
  //
  //       Note that freeing the old pointer must happen AFTER allocating a
  //       new one and copying the data, because otherwise the old pointer may
  //       be given to somebody else who could corrupt the data before we
  //       finish the copying. Even in a single address space this ordering
  //       ensures that a number of weird scenarios will never happen. Bad
  //       examples if the alloc/free order changes could include:
  //
  //       - Metadata allocation for single/both pools manage to get the old
  //         pointer position and corrupt it, before copy is made, because
  //         their pool is exhausted and they get stuff from the common free
  //         space
  //
  //       - Data allocation between pools steal slabs from the original
  //         pool, because after the free it has lots of free space
  //
  //       - New slot range is similar (but not same), so during copy we
  //         screw up our own data if we choose the wrong direction to copy
  //
  //       So, think very, very carefully before optimizing this realloc to
  //       avoid copies about what could possibly go wrong here. 


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Select request/reply buffers
  if (context->pr_num_cores == 1) {
    req   = &single_core_req;
    reply = &single_core_reply;
  }
  else {
    req   = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    reply = NULL;
  }

  // Build message for asking the old pointer size & region ID
  req->core_id = context->pr_core_id;
  req->req_id  = context->pr_message_id;
  req->type    = REQ_QUERY_POINTER;
  req->size    = 0;
  req->region  = 0;
  req->ptr     = old_ptr;

#if 0
  // Track time
  ar_get_time(&sched_sec, &sched_usec);
#endif

  // Single-core fast tracking?
  if (context->pr_num_cores == 1) {
    reply->status = mm_region_query_pointer(req->ptr, &old_size, &old_region);
    reply->result = (void *) old_region;
    reply->size = old_size;
  }
  // Multi-core normal path
  else {

    // Remember which message we're waiting for, send message to scheduler
    // and increase system message ID
    wait_msg_id = req->req_id;
    ar_assert(!noc_msg_send());
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Wait until we get the reply we want
    reply = pr_event_worker_inner_loop(0, REPLY_QUERY_POINTER, wait_msg_id, 
                                       NULL);
    ar_assert(reply);
  }

#if 0
  // Track time
  _sys_track_sched_time(sched_sec, sched_usec);
  context->_sys_sched_calls++;
#endif

  // See what happened
  switch (reply->status) {
    case 0:
      old_size = reply->size;
      old_region = (rid_t) reply->result;
      break;
    case ERR_OUT_OF_RANGE:
      _sys_error(filename, line_nr,
                 "Tried to realloc out of range pointer %p", old_ptr);
    case ERR_MISALIGNED:
      _sys_error(filename, line_nr,
                 "Tried to realloc misaligned pointer %p", old_ptr);
    case ERR_NOT_ALLOCED:
      _sys_error(filename, line_nr,
                 "Tried to realloc non-allocated pointer %p", old_ptr);
    default:
      // Unknown return code
      ar_abort();
  }

  // It may happen that the user is asking for the same quantity in the same
  // region (especially if he's asking for small increases, e.g. 4 bytes, but
  // we are forcing all allocations in MM_ALLOC_ALIGN objects). If so, just
  // return the same pointer.
  if ((old_size == new_size) && (old_region == new_region)) {
    return old_ptr;
  }

  // Allocate new pointer
  new_ptr = _sys_alloc(filename, line_nr, new_size, new_region);

  // It should be valid (zero-size has been handled above)
  ar_assert(new_ptr);

  // Copy contents
  kt_memcpy(new_ptr, old_ptr, (old_size < new_size) ? old_size : new_size);

  // Free old pointer
  _sys_free(filename, line_nr, old_ptr);

  // Success
  return new_ptr;
}


// ===========================================================================
// _sys_ralloc()                Allocates a new region.
//
//                              The system call aborts the execution upon any 
//                              error and prints a message.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   rid_t parent               Region ID of the existing region into which
//                              the new one will be created, or 0 for the 
//                              NULL region.
//   int level_hint             Scheduler level preferred to handle the new
//                              region. May be in the range of 0 (lowest level
//                              scheduler) to L-1, where L is the highest
//                              level of scheduler in the system. The hint
//                              may be silently overridden, depending on
//                              where the parent region is located, if the
//                              level is sane, etc.
//
// * RETURN VALUE
//   rid_t                      Region ID of the newly allocated region
// ===========================================================================
rid_t _sys_ralloc(char *filename, int line_nr, rid_t parent, int level_hint) {

  Context               *context;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  PrMsgReq              single_core_req;
  PrMsgReply            single_core_reply;
  int                   wait_msg_id;
#if 0
  unsigned int          sched_sec;
  unsigned int          sched_usec;
#endif

  // Get context
  context = mm_get_context(ar_get_core_id());

  // Fix NULL region: internally, it has id 1. User should not directly use
  // region ID 1, it's illegal (it's never returned by any of the system
  // calls).
  if (parent == 1) {
    _sys_error(filename, line_nr, "Illegal region id 1");
  }
  if (!parent) {
    parent = 1;
  }

  // Don't accept negatives for level hint. We cast it back and forth
  // with unsigned types during communication, and it doesn't make sense
  // anyway: the lower level of scheduler is always 0.
  if (level_hint < 0) {
    _sys_error(filename, line_nr,
               "Illegal negative level hint %d", level_hint);
  }

  // Select request/reply buffers
  if (context->pr_num_cores == 1) {
    req   = &single_core_req;
    reply = &single_core_reply;
  }
  else {
    req   = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    reply = NULL;
  }

  // Build message
  req->core_id = context->pr_core_id;
  req->req_id  = context->pr_message_id;
  req->type    = REQ_RALLOC;
  req->region  = parent;
  req->ptr     = (void *) ((size_t) level_hint);

#if 0
  // Track time
  ar_get_time(&sched_sec, &sched_usec);
#endif

  // Single-core fast tracking?
  if (context->pr_num_cores == 1) {
    reply->status = mm_ralloc(req, (void *) &(reply->result));
  }
  // Multi-core normal path
  else {

    // Remember which message we're waiting for, send message to scheduler
    // and increase system message ID
    wait_msg_id = req->req_id;
    ar_assert(!noc_msg_send());
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Wait until we get the reply we want
    reply = pr_event_worker_inner_loop(0, REPLY_RALLOC, wait_msg_id, NULL);
    ar_assert(reply);
  }

#if 0
  // Track time
  _sys_track_sched_time(sched_sec, sched_usec);
  context->sys_sched_calls++;
#endif

  // See what happened
  switch (reply->status) {
    case 0:
      return (rid_t) reply->result;
    case ERR_NO_SUCH_REGION:
      _sys_error(filename, line_nr, "Invalid region id %lu", parent);
    case ERR_OUT_OF_MEMORY:
      _sys_error(filename, line_nr, "Out of memory");
    case ERR_OUT_OF_RIDS:
      _sys_error(filename, line_nr, "Out of region IDs");
    default:
      // Unknown return code
      ar_abort();
  }
}


// ===========================================================================
// _sys_rfree()                 Destroys an allocated region, along with all
//                              its objects and children regions.
//
//                              The system call aborts the execution upon any 
//                              error and prints a message.
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//   rid_t region               Region ID of the allocated region to be freed
// ===========================================================================
void _sys_rfree(char *filename, int line_nr, rid_t region) {

  Context               *context;
  PrMsgReq              *req;
  PrMsgReply            *reply;
  PrMsgReq              single_core_req;
  PrMsgReply            single_core_reply;
  int                   wait_msg_id;
#if 0
  unsigned int          sched_sec;
  unsigned int          sched_usec;
#endif


  // No NULL region freeing
  if (!region) {
    _sys_error(filename, line_nr, "Cannot free NULL region");
  }
  // Also prevent from illegal use of region ID 1 (which we interpret
  // internally as the NULL region)
  if (region == 1) {
    _sys_error(filename, line_nr, "Illegal region id 1");
  }

  // Get context
  context = mm_get_context(ar_get_core_id());

  // Select request/reply buffers
  if (context->pr_num_cores == 1) {
    req   = &single_core_req;
    reply = &single_core_reply;
  }
  else {
    req   = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    reply = NULL;
  }

  // Build message
  req->core_id = context->pr_core_id;
  req->req_id  = context->pr_message_id;
  req->type    = REQ_RFREE;
  req->region  = region;
  req->ptr     = (void *) 1;

#if 0
  // Track time
  ar_get_time(&sched_sec, &sched_usec);
#endif

  // Single-core fast tracking?
  if (context->pr_num_cores == 1) {
    reply->status = mm_rfree(req);
  }
  // Multi-core normal path
  else {

    // Remember which message we're waiting for, send message to scheduler
    // and increase system message ID
    wait_msg_id = req->req_id;
    ar_assert(!noc_msg_send());
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Wait until we get the reply we want
    reply = pr_event_worker_inner_loop(0, REPLY_RFREE, wait_msg_id, NULL);
    ar_assert(reply);
  }

#if 0
  // Track time
  _sys_track_sched_time(sched_sec, sched_usec);
  context->sys_sched_calls++;
#endif

  // See what happened
  switch (reply->status) {
    case 0:
      return;
    case ERR_NO_SUCH_REGION:
      _sys_error(filename, line_nr, "Invalid region id %lu", region);
    default:
      // Unknown return code
      ar_abort();
  }
}

