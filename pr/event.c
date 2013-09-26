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
// Abstract      : Scheduler/memory allocator event-based driven processing
//                 functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: event.c,v $
// CVS revision  : $Revision: 1.35 $
// Last modified : $Date: 2013/02/20 08:45:17 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <processing.h>
#include <noc.h>
#include <debug.h>
#include <syscall.h>


// ###########################################################################
// ###                                                                     ###
// ###                          INTERNAL FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################



// ===========================================================================
// pr_event_reply_type()        Returns the appropriate type of reply for
//                              the given request type
// ===========================================================================
// * INPUTS
//   PrMsgType                  The request message type
//
// * RETURN VALUE
//   PrMsgType                  The appropriate reply for the given request
// ===========================================================================
PrMsgType pr_event_reply_type(PrMsgType type) {

  switch (type) {

    case REQ_ALLOC:                     return REPLY_ALLOC;
    case REQ_BALLOC:                    return EXT_REPLY_BALLOC;
    case REQ_FREE:                      return REPLY_FREE;
    case REQ_RALLOC:                    return REPLY_RALLOC;
    case REQ_RALLOC_ORPHAN:             return REPLY_RALLOC;
    case REQ_RFREE:                     return REPLY_RFREE;
    case REQ_RFREE_UPDATE_PARENT:       return REPLY_RFREE;
    case REQ_PACK:                      return EXT_REPLY_PACK;
    case EXT_REQ_PACK:                  return EXT_REPLY_PACK;
    case REQ_QUERY_POINTER:             return REPLY_QUERY_POINTER;
    case REQ_GET_PAGES:                 return REPLY_GET_PAGES;
    case REQ_RETURN_PAGES:              return REPLY_RETURN_PAGES;
    case REQ_GET_RIDS:                  return REPLY_GET_RIDS;
    case REQ_RETURN_RIDS:               return REPLY_RETURN_RIDS;
    case REQ_LOAD_REPORT:               ar_abort(); // No reply
    case REQ_SHUTDOWN:                  ar_abort(); // No reply
    case EXT_REQ_SPAWN:                 return REPLY_SPAWN;
    case EXT_REQ_DELEGATE:              ar_abort(); // No reply
    case EXT_REQ_DEP_START:             ar_abort(); // No reply
    case EXT_REQ_DEP_ROUTE:             ar_abort(); // No reply
    case EXT_REQ_DEP_ENQUEUE:           ar_abort(); // No reply
    case REQ_DEP_OK:                    ar_abort(); // No reply
    case EXT_REQ_DEP_STOP:              ar_abort(); // No reply
    case REQ_DEP_CHILD_FREE:            ar_abort(); // No reply
    case EXT_REQ_UPDATE_LOCATION:       ar_abort(); // No reply
    case EXT_REQ_EXEC:                  ar_abort(); // No reply
    case REQ_EXEC_DONE:                 ar_abort(); // No reply
    case EXT_REQ_SCHEDULE:              return REPLY_SCHEDULE;
    case SELF_RALLOC_UPDATE_PARENT:     return REPLY_RALLOC;
    case SELF_PACK_MERGE:               ar_abort(); // No reply
    case SELF_RFREE_WAIT_CHILDREN:      ar_abort(); // No reply
    case SELF_SCHEDULE_RESULT:          ar_abort(); // No reply
    case SELF_WAIT_SPAWN:               ar_abort(); // No reply

    // Unknown request
    default:                            ar_abort();
  }
}


// ===========================================================================
// pr_event_is_extended()       Determines if the given message type is basic 
//                              or extended
// ===========================================================================
// * INPUTS
//   PrMsgType                  The request message type
//
// * RETURN VALUE
//   int                        0: basic size
//                              1: extended size
// ===========================================================================
int pr_event_is_extended(PrMsgType type) {

  switch (type) {

    case REQ_ALLOC:                     return 0;
    case REQ_BALLOC:                    return 0;
    case REQ_FREE:                      return 0;
    case REQ_RALLOC:                    return 0;
    case REQ_RALLOC_ORPHAN:             return 0;
    case REQ_RFREE:                     return 0;
    case REQ_RFREE_UPDATE_PARENT:       return 0;
    case REQ_PACK:                      return 0;
    case EXT_REQ_PACK:                  return 1;
    case REQ_QUERY_POINTER:             return 0;
    case REQ_GET_PAGES:                 return 0;
    case REQ_RETURN_PAGES:              return 0;
    case REQ_GET_RIDS:                  return 0;
    case REQ_RETURN_RIDS:               return 0;
    case REQ_LOAD_REPORT:               return 0;
    case REQ_SHUTDOWN:                  return 0;
    case EXT_REQ_SPAWN:                 return 1;
    case EXT_REQ_DELEGATE:              return 1;
    case EXT_REQ_DEP_START:             return 1;
    case EXT_REQ_DEP_ROUTE:             return 1;
    case EXT_REQ_DEP_ENQUEUE:           return 1;
    case REQ_DEP_OK:                    return 0;
    case EXT_REQ_DEP_STOP:              return 1;
    case REQ_DEP_CHILD_FREE:            return 0;
    case EXT_REQ_UPDATE_LOCATION:       return 1;
    case EXT_REQ_EXEC:                  return 1;
    case REQ_EXEC_DONE:                 return 0;
    case EXT_REQ_SCHEDULE:              return 1;
    case SELF_RALLOC_UPDATE_PARENT:     return 0;
    case SELF_PACK_MERGE:               return 0;
    case SELF_RFREE_WAIT_CHILDREN:      return 0;
    case SELF_SCHEDULE_RESULT:          return 0;
    case SELF_WAIT_SPAWN:               return 0;
    case REPLY_ALLOC:                   return 0;
    case EXT_REPLY_BALLOC:              return 1;
    case REPLY_FREE:                    return 0;
    case REPLY_RALLOC:                  return 0;
    case REPLY_RFREE:                   return 0;
    case EXT_REPLY_PACK:                return 1;
    case REPLY_QUERY_POINTER:           return 0;
    case REPLY_GET_PAGES:               return 0;
    case REPLY_RETURN_PAGES:            return 0;
    case REPLY_GET_RIDS:                return 0;
    case REPLY_RETURN_RIDS:             return 0;
    case REPLY_SCHEDULE:                return 0;
    case REPLY_SPAWN:                   return 0;

    // Unknown type
    default:                            ar_abort();
  }
}


// ===========================================================================
// pr_event_find_multi()        Finds or creates a node in the incomplete
//                              partial requests list, if an incoming request 
//                              continues further into more messages. Depending
//                              on the mode called, the incoming message can
//                              be stored as a part in the structure.
// ===========================================================================
// * INPUTS
//   Context *context           The context
//   PrMsgReq *in_req           The incoming message
//   int store_part             0: do not store in_req on the structure
//                              1: store in_req on the structure
//
// * OUTPUTS
//   ListNode **ret_list_pos    Returns the list node of the returned
//                              multi-part structure in the incomplete
//                              partial requests list
//
// * RETURN VALUE
//   PrMultiPartReq *           The located multi-part structure, or NULL
//                              if the request is in a single message. Note
//                              that for a normal multi-part message with 
//                              stored parts, the returned structure will have
//                              all the parts except the last one (which is
//                              the in_req itself, not copied to save the
//                              last copy).
// ===========================================================================
PrMultiPartReq *pr_event_find_multi(Context *context, PrMsgReq *in_req, 
                                    int store_part, ListNode **ret_list_pos) {

  PrMultiPartReq        *multi;
  ListNode              *node;


  // Sanity checks
  ar_assert(ret_list_pos);

  // Is this a continuation of a multi-part request?
  multi = NULL;
  node = NULL;
  if (in_req->num_regions) {

    // Find the previous part(s)
    for (node = context->pr_incomplete_req->head; node; node = node->next) {
      if ((((PrMultiPartReq *) node->data)->req_id  == in_req->req_id) &&
          (((PrMultiPartReq *) node->data)->core_id == in_req->core_id)) {
        // Match
        multi = node->data;
        break;
      }
    }
    ar_assert(node);
    ar_assert(multi);
  }

  // Is this request going to be continued?
  if (in_req->num_ptrs == -1) {

    // Did we find any previous parts?
    if (multi) {

      // Store the new part?
      if (store_part) {

        // Attach the new request to the multi-part stored list node
        multi->parts = kt_realloc(multi->parts, 
                                  (multi->num_parts + 1) * sizeof(PrMsgReq *));
        multi->parts[multi->num_parts] = kt_malloc(sizeof(PrMsgReq));
        kt_memcpy(multi->parts[multi->num_parts], in_req, sizeof(PrMsgReq));
        multi->num_parts++;
      }
    }

    // No previous multi, create one now
    else {
      
      // Create new multi-part node
      multi = kt_malloc(sizeof(PrMultiPartReq));
      multi->req_id = in_req->req_id;
      multi->core_id = in_req->core_id;
      multi->forward_idx = -1;
      multi->forward_msg_id = -1;

      // Store the new part?
      if (store_part) {
        multi->parts = kt_malloc(1 * sizeof (PrMsgReq *));
        multi->parts[0] = kt_malloc(sizeof(PrMsgReq));
        kt_memcpy(multi->parts[0], in_req, sizeof(PrMsgReq));
        multi->num_parts = 1;
      }
      else {
        multi->parts = NULL;
        multi->num_parts = 0;
      }

      // Append it to the list
      kt_list_insert(context->pr_incomplete_req, multi,
                     context->pr_incomplete_req->tail, 1);
    }
  }

  // Return the multi node, if we found one
  *ret_list_pos = node;
  return multi;
}


// ===========================================================================
// pr_event_free_multi()        Frees a multi-part structure node, its stored
//                              parts and removes it from the incomplete
//                              requests list
// ===========================================================================
// * INPUTS
//   Context *context           The context
//   ListNode *multi_node       The list node of the multi-part structure in 
//                              the incomplete requests list
// ===========================================================================
void pr_event_free_multi(Context *context, ListNode *multi_node) {

  PrMultiPartReq        *multi;
  int                   i;


  // Assign multi node
  multi = multi_node->data;
  ar_assert(multi);

  // Delete it from the incomplete partial request list
  kt_list_delete(context->pr_incomplete_req, multi_node, NULL);
  
  // Free the stored parts, if any
  if (multi->num_parts) {
    for (i = 0; i < multi->num_parts; i++) {
      kt_free(multi->parts[i]);
    }
    kt_free(multi->parts);
  }
  else {
    ar_assert(!multi->parts);
  }

  // Free the multi structure
  kt_free(multi);
}


// ===========================================================================
// pr_event_copy_multi_arrays() Allocates arrays and populates them with the
//                              task arguments and dependencies from multi-
//                              part requests
// ===========================================================================
// * INPUTS
//   PrMultiPartReq *multi      The multi-part request minus the final part,
//                              or NULL if only a final part exists
//   PrMsgReq *final_part       The final part of the multi-part request
//
// * OUTPUTS
//   void ***ret_args           Returned array of arguments (must exist)
//   unsigned char **ret_deps   If not NULL, returned array of dependencies
//
// * RETURN VALUE
//   int                        Number of arguments copied
// ===========================================================================
int pr_event_create_multi_arrays(PrMultiPartReq *multi, PrMsgReq *final_part,
                                 void ***ret_args, unsigned char **ret_deps) {

  int   num_args;
  int   i;
  int   j;
  int   k;


  // Calculate total number of arguments
  if (final_part->type != EXT_REQ_EXEC) {
    num_args = final_part->num_ptrs;
    if (multi) {
      num_args += multi->num_parts * PR_REQ_MAX_SIZE;
    }
  }
  // EXT_REQ_EXEC is an exception, because the DMA list follows the arguments.
  // The number of arguments is explicitly specified.
  else {
    ar_assert(!ret_deps); // sanity check, size is not a bitmap
    num_args = final_part->size;
  }

  // No arguments?
  if (!num_args) {
    *ret_args = NULL;
    if (ret_deps) {
      *ret_deps = NULL;
    }
    return 0;
  }

  // Allocate array(s)
  *ret_args = kt_malloc(num_args * sizeof(void *));
  if (ret_deps) {
    *ret_deps = kt_malloc(num_args * sizeof(unsigned char));
  }

  // Copy arguments/dependencies for all parts except the final
  k = 0;
  if (multi) {
    for (i = 0; i < multi->num_parts; i++) {
      ar_assert(multi->parts[i]->num_ptrs == -1);
      for (j = 0; j < PR_REQ_MAX_SIZE; j++) {
        // In EXT_REQ_EXEC, finish the loop without going into the DMA list
        if ((final_part->type == EXT_REQ_EXEC) && (k >= num_args)) {
          break;
        }
        (*ret_args)[k] = multi->parts[i]->data[j];
        if (ret_deps) {
          (*ret_deps)[k] = 0;
          (*ret_deps)[k] |= (multi->parts[i]->size & (1 << (2 * j))) 
                          ? SYS_TYPE_IN_ARG : 0;
          (*ret_deps)[k] |= (multi->parts[i]->size & (1 << (2 * j + 1))) 
                          ? SYS_TYPE_OUT_ARG : 0;
          (*ret_deps)[k] |= (multi->parts[i]->region & (1 << (2 * j)))
                          ? SYS_TYPE_SAFE_ARG : 0;
          (*ret_deps)[k] |= (multi->parts[i]->region & (1 << (2 * j + 1)))
                          ? SYS_TYPE_REGION_ARG : 0;
        }
        k++;
      }
    }
  }

  // Copy arguments/dependencies for the final part
  for (j = 0; j < final_part->num_ptrs; j++) {
    // In EXT_REQ_EXEC, finish the loop without going into the DMA list
    if ((final_part->type == EXT_REQ_EXEC) && (k >= num_args)) {
      break;
    }
    (*ret_args)[k] = final_part->data[j];
    if (ret_deps) {
      (*ret_deps)[k] = 0;
      (*ret_deps)[k] |= (final_part->size & (1 << (2 * j))) 
                      ? SYS_TYPE_IN_ARG : 0;
      (*ret_deps)[k] |= (final_part->size & (1 << (2 * j + 1))) 
                      ? SYS_TYPE_OUT_ARG : 0;
      (*ret_deps)[k] |= (final_part->region & (1 << (2 * j)))
                      ? SYS_TYPE_SAFE_ARG : 0;
      (*ret_deps)[k] |= (final_part->region & (1 << (2 * j + 1)))
                      ? SYS_TYPE_REGION_ARG : 0;
    }
    k++;
  }
  ar_assert(k == num_args);

  // Success
  return num_args;
}


// ===========================================================================
// pr_event_req_alloc()         Handler for REQ_ALLOC events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_alloc(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;
  void                  *result;


  // Allocate object
  status = mm_alloc(in_req, &result);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_ALLOC;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;
  out_reply->result = result;

  // Send reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_balloc()        Handler for REQ_BALLOC events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
//   PrEvtHookBalloc *state     Reentrant state structure
// ===========================================================================
void pr_event_req_balloc(PrMsgReq *in_req, PrEvtHookBalloc *state) {

  PrMsgReply            *out_reply;
  int                   status;
  int                   i;
  int                   j;


  // Sanity check
  ar_assert(in_req->ptr);
  
  // Do the bulk alloc
  status = mm_balloc(in_req, &state);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Did the call otherwise fail?
  else if (status) {

    // Send a failure extended reply
    out_reply = noc_msg_send_get_buf(in_req->core_id);

    out_reply->type = EXT_REPLY_BALLOC;
    out_reply->req_id = in_req->req_id;
    out_reply->status = status;

    ar_assert(!noc_msg_send());

    // Free stuff & return
    if (state) {
      kt_free(state->objects);
      kt_free(state);
    }
    return;
  }

  // Sanity check: results are in saved state structure
  ar_assert(state);
  ar_assert(state->objects);
  ar_assert(state->num_objects);

  // Create as many extended replies as needed
  out_reply = NULL;
  j = 0;
  for (i = 0; i < state->num_objects; i++) {

    // New packet?
    if (!j) {

      // Get a new reply buffer
      out_reply = noc_msg_send_get_buf(in_req->core_id);

      // Fill standards
      out_reply->type = EXT_REPLY_BALLOC;
      out_reply->req_id = in_req->req_id;
      out_reply->status = 0;
      out_reply->result = NULL;
    }

    // Fill current field
    out_reply->adr_ar[j] = state->objects[i];
    out_reply->size_ar[j] = in_req->size;
    j++;

    // Last field for this packet?
    if ((j == PR_RPL_MAX_SIZE) || (i == state->num_objects - 1)) {

      if (i == state->num_objects - 1) {
        // Last out_reply, number of elements is exact
        out_reply->num_elements = j;
      }
      else {
        // More replies will follow
        out_reply->num_elements = -1;
      }

      // Send current reply
      ar_assert(!noc_msg_send());

      out_reply = NULL;
      j = 0;
    }
  }

  // Free stuff
  kt_free(state->objects);
  kt_free(state);
}


// ===========================================================================
// pr_event_req_free()          Handler for REQ_FREE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_free(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;


  // Free object
  status = mm_free(in_req);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_FREE;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;

  // Send reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_ralloc()        Handler for REQ_RALLOC events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_ralloc(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;
  void                  *result;


  // Allocate the region
  status = mm_ralloc(in_req, (void *) &result);
  
  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_RALLOC;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;
  out_reply->result = result;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_ralloc_orphan() Handler for REQ_RALLOC_ORPHAN events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_ralloc_orphan(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;
  void                  *result;


  // Allocate the orphan
  status = mm_ralloc_orphan(in_req, (void *) &result);
  
  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_RALLOC;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;
  out_reply->result = result;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_rfree()         Handler for REQ_RFREE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_rfree(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;


  // Free the region
  status = mm_rfree(in_req);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_RFREE;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_rfree_update_parent()  Handler for REQ_RFREE_UPDATE_PARENT
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_rfree_update_parent(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;


  // Update the parent
  status = mm_rfree_update_parent(in_req);

  // This can't fail
  ar_assert(!status);

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_RFREE;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_reply_pack()        Helper function that creates replies for the
//                              accumulated state of packing operations
// ===========================================================================
// * INPUTS
//   int core_id                Core ID to reply to
//   int msg_id                 Message ID of request we're replying to
//   PrEvtHookPack *state       Reentrant state structure
// ===========================================================================
void pr_event_reply_pack(int core_id, int req_id, PrEvtHookPack *state) {

  PrMsgReply            *out_reply;
  int                   i;
  int                   j;


  // Sanity checks
  ar_assert(state);
  ar_assert(!state->wait_replies);
  ar_assert(!state->native_task_id);
  ar_assert(state->addresses); // these two are always allocated, even if
  ar_assert(state->sizes);     // num_elements is 0

  // Corner case: it may be that the user's querying empty regions, so nothing
  // is returned. Generate a distinct error code (which will be squelched
  // later on internally) to handle communication on this case.
  if (!state->error_status && !state->num_elements) {
    
    // Just send a single failure extended reply
    out_reply = noc_msg_send_get_buf(core_id);

    out_reply->type = EXT_REPLY_PACK;
    out_reply->req_id = req_id;
    out_reply->status = ERR_NOTHING_TO_PACK;
    out_reply->result = NULL;
    ar_assert(!noc_msg_send());

    return;
  }

  // Other failure?
  if (state->error_status) {

    // Just send a single failure extended reply
    out_reply = noc_msg_send_get_buf(core_id);

    out_reply->type = EXT_REPLY_PACK;
    out_reply->req_id = req_id;
    out_reply->status = state->error_status;
    out_reply->result = state->error_ptr;
    ar_assert(!noc_msg_send());

    return;
  }

  // Create as many extended replies as needed
  out_reply = NULL;
  j = 0;
  for (i = 0; i < state->num_elements; i++) {

    // New packet?
    if (!j) {
      out_reply = noc_msg_send_get_buf(core_id);

      out_reply->type = EXT_REPLY_PACK;
      out_reply->req_id = req_id;
      out_reply->status = 0;
      out_reply->result = NULL;
    }

    // Add new field
    out_reply->adr_ar[j] = (void *) state->addresses[i];
    out_reply->size_ar[j] = state->sizes[i];
    j++;

    if ((j == PR_RPL_MAX_SIZE) || (i == state->num_elements - 1)) {

      if (i == state->num_elements - 1) {
        // Last out_reply, number of elements is exact
        out_reply->num_elements = j;
      }
      else {
        // More replies will follow
        out_reply->num_elements = -1;
      }

      // Send current reply
      ar_assert(!noc_msg_send());

      out_reply = NULL;
      j = 0;
    }
  }
}


// ===========================================================================
// pr_event_req_pack()          Handler for REQ_PACK events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_pack(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  PrEvtHookPack         *state;
  int                   status;


  // Do the packing
  state = NULL;
  status = mm_pack(in_req, NULL, &state);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Sanity check: if we came from an event, this is not a native call and
  // we must always reply with a message (error or success)
  ar_assert(!state->native_task_id);

  // Failure?
  if (status) {

    // Send a failure extended reply
    out_reply = noc_msg_send_get_buf(in_req->core_id);

    out_reply->type = EXT_REPLY_PACK;
    out_reply->req_id = in_req->req_id;
    out_reply->status = status;
    out_reply->result = state->error_ptr;

    ar_assert(!noc_msg_send());

    // Free stuff and return
    kt_free(state->sizes);
    kt_free(state->addresses);
    kt_free(state);
    return;
  }

  // Send reply based on what's gathered in the state
  pr_event_reply_pack(in_req->core_id, in_req->req_id, state);

  // Free stuff
  kt_free(state->sizes);
  kt_free(state->addresses);
  kt_free(state);
}


// ===========================================================================
// pr_event_req_query_pointer() Handler for REQ_QUERY_POINTER events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_query_pointer(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;
  size_t                size;
  void                  *result;


  // Query
  status = mm_query_pointer(in_req, &size, (void *) &result);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_QUERY_POINTER;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;
  out_reply->size = size;
  out_reply->result = result;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_get_pages()     Handler for REQ_GET_PAGES events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_get_pages(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;
  void                  *result;

  
  // Get the pages
  status = mm_get_pages(in_req, (void *) &result);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_GET_PAGES;
  out_reply->req_id = in_req->req_id;
  out_reply->status = status;
  out_reply->size = in_req->size;
  out_reply->result = result;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_get_rids()      Handler for REQ_GET_RIDS events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_get_rids(PrMsgReq *in_req) {

  PrMsgReply            *out_reply;
  int                   status;
  void                  *result;

  
  // Get the region IDs
  status = mm_get_rids(in_req, (void *) &result);

  // Should we send a reply now?
  if (status == ERR_REPLY_POSTPONED) {
    return;
  }

  // Prepare reply
  out_reply = noc_msg_send_get_buf(in_req->core_id);

  out_reply->type = REPLY_GET_RIDS;
  out_reply->req_id = in_req->req_id;
  out_reply->size = in_req->size;
  out_reply->status = status;
  out_reply->result = result;

  // Reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_req_load_report()   Handler for REQ_LOAD_REPORT events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
//
// * RETURN VALUE
//   int                        0: normal exit
//                              1: shutdown was sent to our children
// ===========================================================================
int pr_event_req_load_report(PrMsgReq *in_req) {

  Context       *context;
  int           region_load;
  int           sched_task_load;
  int           run_task_load;
  int           child_id;
  int           status;
  PrMsgReq      *req;
  int           i;


  // Get needed fields
  context = mm_get_context(ar_get_core_id());
  region_load = in_req->region;
  run_task_load = (size_t) in_req->ptr;
  sched_task_load = in_req->size;

  // Set default exit status to normal
  status = 0;

  // Find out which child of ours is this
  child_id = pr_core_child_index(in_req->core_id);
  ar_assert((child_id >= 0) && (child_id < context->pr_num_children));

  // Update our totals
  context->mm_current_load   += region_load - 
                                        context->mm_children_load[child_id];
  context->pr_cur_run_load   += run_task_load - 
                                        context->pr_chld_run_load[child_id];
  context->pr_cur_sched_load += sched_task_load - 
                                        context->pr_chld_sched_load[child_id];

  // Update child's value
  context->mm_children_load[child_id]   = region_load;
  context->pr_chld_run_load[child_id]   = run_task_load;
  context->pr_chld_sched_load[child_id] = sched_task_load;

//kt_printf("%d: New load %d/%d/%d (core %d reported %d/%d/%d)\r\n", context->pr_core_id, context->mm_current_load, context->pr_cur_run_load, context->pr_cur_sched_load, in_req->core_id, region_load, run_task_load, sched_task_load);

  // Are we the top-level scheduler?
  if (context->pr_parent_sched_id == -1) {

    // Check if the main task has finished and no task load exist anywhere.
    if ((context->pr_main_finished) && 
        (!context->pr_cur_run_load) && (!context->pr_cur_sched_load)) {

//kt_printf("%d: Begin shutdown\r\n", context->pr_core_id);

      // Begin the system shutdown
      for (i = 0; i < context->pr_num_children; i++) {

        // Send request
        req = noc_msg_send_get_buf(context->pr_children[i]);
        req->type = REQ_SHUTDOWN;
        req->req_id = context->pr_message_id;
        ar_assert(!noc_msg_send());

        // Increase message ID
        context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
      }

      // Indicate shutdown to the event handler
      status = 1;
    }

  }
  else {
    // Check if we need to report it upstream
    pr_sched_report_load();
  }


  return status;
}


// ===========================================================================
// pr_event_req_shutdown()      Handler for REQ_SHUTDOWN events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_shutdown(PrMsgReq *in_req) {

  Context               *context;
  PrMsgReq              *out_req;
  int                   i;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Notify our children
  if (context->pr_children) {
    for (i = 0; i < context->pr_num_children; i++) {

      // Prepare message
      out_req = noc_msg_send_get_buf(context->pr_children[i]);

      out_req->type = REQ_SHUTDOWN;
      out_req->req_id = context->pr_message_id;

      // Send it
      ar_assert(!noc_msg_send());

      // Increase message ID
      context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
    }
  }

}


// ===========================================================================
// pr_event_req_spawn_forward() Handler for EXT_REQ_SPAWN events, in case 
//                              we're just an intermediate scheduler
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_spawn_forward(PrMsgReq *in_req) {

  Context               *context;
  ListNode              *multi_node;
  PrMultiPartReq        *multi;
  PrMsgReq              *fwd_req;
  PrEvtPending          *event;
  int                   msg_id;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->pr_parent_sched_id != -1);

  // If it's a multi-part request, find or create the relevant structure, 
  // but don't store any parts: we'll forward each part on the fly.
  multi = pr_event_find_multi(context, in_req, 0, &multi_node);

  // Have we decided a message ID for this forwarding before?
  if ((multi) && (multi->forward_msg_id != -1)) {
    msg_id = multi->forward_msg_id;
  }
  else {
    // Use a new message ID
    msg_id = context->pr_message_id;
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
    if (multi) {
      multi->forward_msg_id = msg_id;
    }

    // Create an event to wait for the reply
    event = kt_malloc(sizeof(PrEvtPending));
    event->req = kt_malloc(sizeof(PrMsgReq));

    // The only thing that matters from the original request is the requestor
    // core ID and message ID, so the event handler knows where to forward
    // the reply to.
    event->req->core_id = in_req->core_id;
    event->req->req_id  = in_req->req_id;
    event->action = PR_ACT_FORWARD;
    event->prev = NULL;
    event->next = NULL;
    event->data = NULL;

    // Store event; we don't expect conflicts on this message ID.
    ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));
  }

  // Forward to our parent scheduler
  fwd_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));

  kt_memcpy(fwd_req, in_req, sizeof(PrMsgReq));
  fwd_req->core_id = context->pr_core_id; // Overwrite with our core ID
  fwd_req->req_id  = msg_id;              // Overwrite with determined msg ID
  
  ar_assert(!noc_msg_send());


  // If it was the last part, remove multi-part node from the list and
  // free the multi-part node
  if ((multi) && (in_req->num_ptrs != -1)) {
    pr_event_free_multi(context, multi_node);
  }
}


// ===========================================================================
// pr_event_req_spawn_final()   Handler for EXT_REQ_SPAWN events, if we are
//                              the responsible scheduler of the parent
//                              task that's doing the spawning. 
//
//                              We use the same function to also handle the
//                              EXT_REQ_DELEGATE events, since the code
//                              path is almost identical.
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_spawn_final(PrMsgReq *in_req) {

  Context               *context;
  ListNode              *multi_node;
  PrMultiPartReq        *multi;
  int                   child;
  void                  **args;
  unsigned char         *deps;
  int                   num_args;
  PrTaskDescr           *task;
  unsigned int          parent_task_id;
  PrTaskDescr           *parent_task;
  int                   child_core_id;
  PrMsgReq              *fwd_req;
  int                   msg_id;
  int                   i;
  int                   j;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // If it's a multi-part request, find or create the relevant structure
  // and store the new part in it if the request will continue further
  multi = pr_event_find_multi(context, in_req, 1, &multi_node);

  // If the request continues, wait for all the parts. We're done for now.
  if (in_req->num_ptrs == -1) {
    ar_assert(multi);
    return;
  }


  // Determine parent task ID (ID of the task that's doing the spawning)
  parent_task_id = ((size_t) in_req->ptr) >> PR_TASK_IDX_SIZE;


  // Are we the parent task scheduler?
  if ((parent_task_id >> (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS)) == 
      context->pr_scheduler_id) {

    // Locate parent task
    parent_task = NULL;
    kt_trie_find(context->pr_tasks, parent_task_id, (void *) &parent_task);
    ar_assert(parent_task);
    ar_assert(parent_task->id == parent_task_id);

    // Store the number of arguments that will be required for dependency
    // analysis, so we can know when to send the reply for the spawn.
    // Ignore by-value arguments (SYS_TYPE_BYVALUE_ARG: in, !out, safe and
    // !region flags) and NULL pointers (NULL data, !region flag)
    parent_task->spawn_rem_deps = 0;
    for (i = 0; i < in_req->num_ptrs; i++) {
      if (!(( (in_req->size   & (1 << (2 * i))) &&
             !(in_req->size   & (1 << (2 * i + 1))) &&
              (in_req->region & (1 << (2 * i))) &&
             !(in_req->region & (1 << (2 * i + 1)))) ||
            (!in_req->data[i] &&  
             !(in_req->region & (1 << (2 * i + 1)))))) {
        parent_task->spawn_rem_deps++;
      }
    }
    if (multi) {
      for (i = 0; i < multi->num_parts; i++) {
        ar_assert(multi->parts[i]->num_ptrs == -1);
        for (j = 0; j < PR_REQ_MAX_SIZE; j++) {
          if (!(( (multi->parts[i]->size   & (1 << (2 * j))) &&
                 !(multi->parts[i]->size   & (1 << (2 * j + 1))) &&
                  (multi->parts[i]->region & (1 << (2 * j))) &&
                 !(multi->parts[i]->region & (1 << (2 * j + 1)))) ||
                (!multi->parts[i]->data[j] &&
                 !(multi->parts[i]->region & (1 << (2 * j + 1)))))) { 
            parent_task->spawn_rem_deps++;
          }
        }
      }
    }

    // Store the request message ID and core ID, so we can know how to reply
    ar_assert(!parent_task->spawn_msg_id);
    ar_assert(parent_task->spawn_core_id == -1);
    parent_task->spawn_msg_id = in_req->req_id;
    parent_task->spawn_core_id = in_req->core_id;
  }


  // Are we a non-leaf scheduler and if so, should we delegate this task
  // to one of our children schedulers?
  if ((context->pr_scheduler_level > 0) && 
      (pr_task_must_delegate(multi, in_req, &child))) {

    // For both cases (EXT_REQ_SPAWN and EXT_REQ_DELEGATE), forward all message
    // parts using a EXT_REQ_DELEGATE message type and remember the message ID
    // we're using.
    msg_id = context->pr_message_id;
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
    child_core_id = pr_scheduler_core_id(child);

//kt_printf("%d: Delegating task to core %d\r\n", context->pr_core_id, child_core_id);

    // Forward all parts except the last one
    if (multi) {
      for (i = 0; i < multi->num_parts; i++) {
        fwd_req = noc_msg_send_get_buf(child_core_id);
        kt_memcpy(fwd_req, multi->parts[i], sizeof(PrMsgReq));
        fwd_req->type = EXT_REQ_DELEGATE;       // Overwrite request type
        fwd_req->core_id = context->pr_core_id; // Overwrite with our core ID
        fwd_req->req_id  = msg_id;        // Overwrite with determined msg ID
        ar_assert(!noc_msg_send());
      }
    }

    // Forward last part
    fwd_req = noc_msg_send_get_buf(child_core_id);
    kt_memcpy(fwd_req, in_req, sizeof(PrMsgReq));
    fwd_req->type = EXT_REQ_DELEGATE;       // Overwrite request type
    fwd_req->core_id = context->pr_core_id; // Overwrite with our core ID
    fwd_req->req_id  = msg_id;        // Overwrite with determined msg ID
    ar_assert(!noc_msg_send());

    // Remove multi-part node from the list and free the parts
    if (multi) {
      pr_event_free_multi(context, multi_node);
    }

    // That's all
    return;
  }


//kt_printf("%d: Chose not to delegate spawn of idx = %d\r\n", context->pr_core_id, ((size_t) in_req->ptr) & ((1 << PR_TASK_IDX_SIZE) - 1));


  // Allocate arguments and dependencies arrays and populate them from the
  // (possibly multi-part) request
  num_args = pr_event_create_multi_arrays(multi, in_req, &args, &deps);
  ar_assert(num_args);
  
  // Create task
  task = pr_task_create(parent_task_id, 
                        ((size_t) in_req->ptr) & ((1 << PR_TASK_IDX_SIZE) - 1),
                        args, deps, num_args);
  ar_assert(task);


  // Start dependency analysis for the arguments
  pr_dep_start_stop(1, parent_task_id, task->id, num_args, args, deps);

  // Remove multi-part node from the list and free the parts
  if (multi) {
    pr_event_free_multi(context, multi_node);
  }
}


// ===========================================================================
// pr_event_req_dep_start()     Handler for EXT_REQ_DEP_START events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_dep_start(PrMsgReq *in_req) {

  Context               *context;
  void                  **args;
  unsigned char         *deps;
  int                   num_args;
  unsigned int          parent_task_id;
  unsigned int          new_task_id;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  parent_task_id = (unsigned int) in_req->ptr;
  new_task_id = in_req->num_regions;

  // Allocate and populate arguments and dependencies arrays
  num_args = pr_event_create_multi_arrays(NULL, in_req, &args, &deps);
  ar_assert(num_args);
  
  // Start dependency analysis for the arguments
  pr_dep_start_stop(1, parent_task_id, new_task_id, num_args, args, deps);

  // Free args/deps arrays
  kt_free(args);
  kt_free(deps);
}


// ===========================================================================
// pr_event_req_dep_stop()      Handler for EXT_REQ_DEP_STOP events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_dep_stop(PrMsgReq *in_req) {

  Context       *context;
  void          **args;
  unsigned char *deps;
  int           num_args;
  unsigned int  task_id;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  task_id = in_req->num_regions;

  // Allocate and populate arguments and dependencies arrays
  num_args = pr_event_create_multi_arrays(NULL, in_req, &args, &deps);
  ar_assert(num_args);
  
  // Go to next dependency for the arguments
  pr_dep_start_stop(0, 0, task_id, num_args, args, deps);

  // Free args/deps arrays
  kt_free(args);
  kt_free(deps);
}


// ===========================================================================
// pr_event_req_dep_route()     Handler for EXT_REQ_DEP_ROUTE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_dep_route(PrMsgReq *in_req) {

  Context               *context;
  unsigned int          parent_task_id;
  unsigned int          new_task_id;
  unsigned char         dep_flags;
  void                  *dep;
  rid_t                 highest_rid;
  rid_t                 *parents;
  int                   num_parents;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  parent_task_id = in_req->size;
  highest_rid    = in_req->region;
  dep            = in_req->ptr;
  num_parents    = in_req->num_regions;
  new_task_id    = in_req->num_ptrs >> PR_TASK_IDX_SIZE;
  dep_flags      = in_req->num_ptrs & ((1 << PR_TASK_IDX_SIZE) - 1);
  parents        = (rid_t *) in_req->data;

  // Process message
  ar_assert(!pr_dep_route(parent_task_id, new_task_id, dep, dep_flags, 
                          0, num_parents, parents, highest_rid));
}


// ===========================================================================
// pr_event_req_dep_enqueue()   Handler for EXT_REQ_DEP_ENQUEUE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_dep_enqueue(PrMsgReq *in_req) {

  Context               *context;
  unsigned int          parent_task_id;
  unsigned int          new_task_id;
  unsigned char         dep_flags;
  void                  *dep;
  rid_t                 highest_rid;
  rid_t                 *parents;
  int                   num_parents;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  parent_task_id = in_req->size;
  highest_rid    = in_req->region;
  dep            = in_req->ptr;
  num_parents    = in_req->num_regions;
  new_task_id    = in_req->num_ptrs >> PR_TASK_IDX_SIZE;
  dep_flags      = in_req->num_ptrs & ((1 << PR_TASK_IDX_SIZE) - 1);
  parents        = (rid_t *) in_req->data;

  // Process message
  ar_assert(!pr_dep_descend_enqueue(parent_task_id, new_task_id, dep, dep_flags,
                                    num_parents, parents, highest_rid, 0, 0));
}


// ===========================================================================
// pr_event_req_dep_ok()        Handler for REQ_DEP_OK events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_dep_ok(PrMsgReq *in_req) {

  Context               *context;
  unsigned int          new_task_id;
  void                  *dep;
  PrMsgReq              *fwd_req;
  int                   sched_id;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  new_task_id   = in_req->size;
  dep           = in_req->ptr;

  // Is it ours?
  sched_id = new_task_id >> (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS);
  if (sched_id == context->pr_scheduler_id) {
    // Handle it locally
    ar_assert(!pr_dep_at_head_of_queue(new_task_id, dep));
  }
  else {
    // Forward message to parent scheduler
    ar_assert(context->pr_parent_sched_id != -1);

    fwd_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                context->pr_parent_sched_id));
    fwd_req->type    = REQ_DEP_OK;
    fwd_req->core_id = context->pr_core_id;    // Overwrite with our core ID
    fwd_req->req_id  = context->pr_message_id; // Overwrite with our message ID
    fwd_req->size    = in_req->size;
    fwd_req->ptr     = in_req->ptr;
    
    ar_assert(!noc_msg_send());

    // Increase message ID, avoiding value 0 on wrap-arounds
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
  }
}


// ===========================================================================
// pr_event_req_dep_child_free() Handler for REQ_DEP_CHILD_FREE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req            Input request buffer
// ===========================================================================
void pr_event_req_dep_child_free(PrMsgReq *in_req) {

  Context       *context;
  void          *dep;
  int           child_decr_ro;
  int           child_decr_rw;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  dep           = (void *) in_req->region;
  child_decr_ro = in_req->size;
  child_decr_rw = (int) in_req->ptr;

  // Handle it locally
  ar_assert(!pr_dep_next(0, dep, 1, 1, child_decr_ro, child_decr_rw));
}


// ===========================================================================
// pr_event_req_dep_update_loc() Handler for EXT_REQ_UPDATE_LOCATION events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req            Input request buffer
// ===========================================================================
void pr_event_req_dep_update_loc(PrMsgReq *in_req) {

  Context               *context;
  void                  **args;
  unsigned char         *deps;
  int                   num_args;
  int                   location;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  location = in_req->num_regions;

  // Allocate and populate arguments and dependencies arrays
  num_args = pr_event_create_multi_arrays(NULL, in_req, &args, &deps);
  ar_assert(num_args);
  
  // Update location of the arguments
  pr_sched_update_location(location, num_args, args, deps, 1);

  // Free args/deps arrays
  kt_free(args);
  kt_free(deps);
}


// ===========================================================================
// pr_event_req_exec_scheduler() Handler for EXT_REQ_EXEC events, in case 
//                               we're a scheduler core
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_exec_scheduler(PrMsgReq *in_req) {

  Context               *context;
  ListNode              *multi_node;
  PrMultiPartReq        *multi;
  PrMsgReq              *fwd_req;
  int                   child;
  int                   msg_id;


  // Get context
  context = mm_get_context(ar_get_core_id());


  // If it's a multi-part request, find or create the relevant structure, 
  // but don't store any parts; we'll forward parts on the fly.
  multi = pr_event_find_multi(context, in_req, 0, &multi_node);


  // Have we decided where are we forwarding this before?
  if ((multi) && (multi->forward_idx != -1)) {
    child = multi->forward_idx;
  }
  else {
    // The final destination unique core ID is in the message. Find which
    // child subtree we must forward the message based on this.
    ar_assert(in_req->region < context->pr_num_cores);
    child = pr_core_child_index(context->pr_core_route[in_req->region]);
    if (multi) {
      multi->forward_idx = child;
    }
  }

  // Have we decided a message ID for this forwarding before?
  if ((multi) && (multi->forward_msg_id != -1)) {
    msg_id = multi->forward_msg_id;
  }
  else {
    // Use a new message ID
    msg_id = context->pr_message_id;
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
    if (multi) {
      multi->forward_msg_id = msg_id;
    }
  }

  // Do the actual forwarding
  fwd_req = noc_msg_send_get_buf(context->pr_children[child]); 

  kt_memcpy(fwd_req, in_req, sizeof(PrMsgReq));
  fwd_req->core_id = context->pr_core_id; // Overwrite with our core ID
  fwd_req->req_id  = msg_id;              // Overwrite with determined msg ID
  
  ar_assert(!noc_msg_send());

  
  // If it was the last part, remove multi-part node from the list and
  // free the multi-part node
  if ((multi) && (in_req->num_ptrs != -1)) {
    pr_event_free_multi(context, multi_node);
  }
  
  // When forwarding of last (or only) part is finished, if we're the leaf
  // scheduler responsible for a worker child we're forwarding the request to,
  // mark the child run task load change, increase our run task load and
  // possibly report upstream.
  if ((!context->pr_scheduler_level) && (in_req->num_ptrs != -1)) {
    context->pr_chld_run_load[child]++;
    context->pr_cur_run_load++;
    pr_sched_report_load();
  }
}


// ===========================================================================
// pr_event_req_exec_worker()   Handler for EXT_REQ_EXEC events, in case 
//                              we're a worker core
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_exec_worker(PrMsgReq *in_req) {

  Context               *context;
  ListNode              *multi_node;
  PrMultiPartReq        *multi;
  PrMsgReq              *out_req;
  PrTaskDescr           *task;
  ListNode              *node;
  int                   total;
  int                   idx;
  ListNode              *nxt_node;
  PrTaskDescr           *nxt_task;
  int                   i;
  int                   j;
  int                   k;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // If it's a multi-part request, find or create the relevant structure
  // and store the new part in it if the request will continue further
  multi = pr_event_find_multi(context, in_req, 1, &multi_node);

  // If the request continues, wait for all the parts. We're done for now.
  if (in_req->num_ptrs == -1) {
    ar_assert(multi);
    return;
  }

  // Take a note for the statistics gathering
  dbg_scount(context, DBG_STATS_IDX_NUM_TASKS, 1);

  // Create a new task descriptor
  task = kt_malloc(sizeof(PrTaskDescr));
  task->id    = ((size_t) in_req->ptr) >> PR_TASK_IDX_SIZE;
  task->index = ((size_t) in_req->ptr) & ((1 << PR_TASK_IDX_SIZE) - 1);
  task->run_core_id = in_req->region;

  // Allocate argument array and copy arguments
  task->num_args = pr_event_create_multi_arrays(multi, in_req, &(task->args), 
                                                NULL);
  if (!task->num_args) {
    ar_assert(!task->index); // Sanity check: Only main task can have no args
  }

  // Unused fields in worker mode
  task->pid            = -1;
  task->deps           = NULL;
  task->deps_waiting   = NULL;
  task->spawn_msg_id   = 0;
  task->spawn_core_id  = -1;
  task->spawn_rem_deps = 0;

  // Extract DMA address/sizes/location array from the request, they are
  // located after the (already extracted) task arguments.
  task->dmas_started = 0;
  task->dmas_waiting = 0;

  total = in_req->num_ptrs;
  if (multi) {
    total += multi->num_parts * PR_REQ_MAX_SIZE;
  }
  total -= task->num_args;
  ar_assert(total >= 0);
  ar_assert(total % 2 == 0);

  // Allocate task DMA arrays
  task->num_dmas = total / 2;
  task->dma_addresses = kt_malloc(task->num_dmas * sizeof(size_t));
  task->dma_sizes     = kt_malloc(task->num_dmas * sizeof(int));

  // Do the multi parts, if available
  k = 0;
  if (multi) {
    for (i = 0; i < multi->num_parts; i++) {
      ar_assert(multi->parts[i]->num_ptrs == -1);
      for (j = 0; j < PR_REQ_MAX_SIZE; j++) {
        if (k >= task->num_args) {
          idx = k - task->num_args;
          if (idx % 2 == 0) {
            task->dma_addresses[idx / 2] = (size_t) multi->parts[i]->data[j];
          }
          else {
            task->dma_sizes[idx / 2] = (int) multi->parts[i]->data[j];
          }
        }
        k++;
      }
    }
  }

  // Continue with the final part
  for (j = 0; j < in_req->num_ptrs; j++) {
    if (k >= task->num_args) {
      idx = k - task->num_args;
      if (idx % 2 == 0) {
        task->dma_addresses[idx / 2] = (size_t) in_req->data[j];
      }
      else {
        task->dma_sizes[idx / 2] = (int) in_req->data[j];
      }
    }
    k++;
  }
  ar_assert(k == total + task->num_args);



  // Remove multi-part node from the list and free the parts
  if (multi) {
    pr_event_free_multi(context, multi_node);
  }

  // Append task at the end of our ready queue
  node = kt_list_insert(context->pr_ready_queue, task, 
                        kt_list_tail(context->pr_ready_queue), 1);
  ar_assert(node);


  // If we are at the head of the queue (i.e. the task that will currently run
  // because the queue was empty), or if we are the second task from the
  // the head (i.e. we are enqueued as the next task behind the one that is
  // running currently), start the DMAs for the task.
  if ((kt_list_head(context->pr_ready_queue) == node) || 
      (kt_list_next(kt_list_head(context->pr_ready_queue)) == node)) {
    ar_assert(!pr_task_start_dmas(task));
  }

  // If we are not at the top of the queue, it means another task is already
  // running and we're called here to handle a new execution request while
  // the running task is communicating with the scheduler and waits for some
  // reply. We're done in this case.
  if (kt_list_head(context->pr_ready_queue) != node) {
    return;
  }

  
  // Enter a loop of servicing ready queue tasks. Although when we enter this
  // loop we're the only task at the head of the queue, this function can
  // be called to handle new execution requests while the first task is
  // running. When the first task has finished, there may be others queued up
  // by the code above.
  while (kt_list_size(context->pr_ready_queue)) {

    // Current task is the head of the queue
    node = kt_list_head(context->pr_ready_queue);
    ar_assert(node);
    task = node->data;
    ar_assert(task);

    // Are the current task DMAs complete? If not, wait until they are.
    ar_assert(task->dmas_started);
    if (task->dmas_waiting) {
//kt_printf("%d: Waiting DMAs for task 0x%08X...\r\n", context->pr_core_id, task->id);
      ar_assert(!pr_event_worker_inner_loop(2, 0, 0, task));
//kt_printf("%d: DMAs for task 0x%08X done\r\n", context->pr_core_id, task->id);
    }
    ar_assert(!task->dmas_waiting);

    // Advance epoch
#ifdef ARCH_MB
    context->pr_cur_epoch = (context->pr_cur_epoch + 1) % AR_CACHE_NUM_EPOCHS;
    ar_cache_set_epoch(context->pr_cur_epoch);
#endif

    // We want to maintain the invariant that we always have started the DMAs
    // of both the "current" (this one, head of the queue) and the "next" task
    // (2nd in the queue). The insertion code above does that; we must also do
    // it inside this loop, as we are dequeueing executed tasks.
    nxt_node = kt_list_next(node);
    if (nxt_node) {
      nxt_task = nxt_node->data;
      ar_assert(nxt_task);
      if (!nxt_task->dmas_started) {
        ar_assert(!pr_task_start_dmas(nxt_task));
      }
    }

//kt_printf("%d: Running task 0x%08X\r\n", context->pr_core_id, task->id);
//for (i = 0; i < task->num_dmas; i++) {
//  kt_printf("   DMA %d: 0x%08X, %d bytes, %d -> %d\r\n", 
//      i, task->dma_addresses[i],
//      task->dma_sizes[i] & ((1 << MM_PACK_SIZE_BITS) - 1),
//      task->dma_sizes[i] >> MM_PACK_SIZE_BITS, context->pr_core_id);
//}

    // Execute it
    ar_exec(context->pr_task_table[task->index], task->num_args, task->args, 
            NULL);

#ifdef ARCH_MB
//*MBS_CHE_MAINT = 0x20000; // flush L2 only
//*MBS_CHE_MAINT = 0x20001; // flush L2, clear data L1
#endif

//kt_printf("%d: Finished task 0x%08X\r\n", context->pr_core_id, task->id);

    // Being at this point means the task has completed (any intermediate
    // communication will be handled by the system calls which will handle
    // their own events). 


    // If a task spawn is pending, we must wait for it. Otherwise, we risk
    // this task being collected while the task scheduler tries to look up
    // the parent task at some point during the pending spawn and not finding
    // it there...
    if (context->pr_spawn_pending) {
      ar_assert(!pr_event_worker_inner_loop(1, 0, 0, NULL));
    }
    ar_assert(!context->pr_spawn_pending);

    // Send a message to parent scheduler that the task is finished
    out_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                  context->pr_parent_sched_id));

    out_req->type    = REQ_EXEC_DONE;
    out_req->core_id = context->pr_core_id;
    out_req->req_id  = context->pr_message_id;
    out_req->ptr     = (void *) ((task->id << PR_TASK_IDX_SIZE) | task->index);

    ar_assert(!noc_msg_send());

    // Increase message ID
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);

    // Free task structure and delete it from the list
    ar_assert(kt_list_head(context->pr_ready_queue) == node);
    kt_free(task->args);
    ar_assert(!task->deps);
    kt_free(task->dma_addresses);
    kt_free(task->dma_sizes);
    kt_list_delete(context->pr_ready_queue, node, kt_free);
  }

}


// ===========================================================================
// pr_event_req_exec_done()     Handler for REQ_EXEC_DONE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
//
// * RETURN VALUE
//   int                        0: normal exit
//                              1: shutdown was sent to our children
// ===========================================================================
int pr_event_req_exec_done(PrMsgReq *in_req) {

  Context               *context;
  int                   child;
  unsigned int          task_id;
  PrMsgReq              *out_req;
  int                   status;
  int                   i;


  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->pr_role == ROLE_SCHEDULER);

  // Set default status to normal
  status = 0;

  // Are we the responsible scheduler of the finished task?
  task_id = ((size_t) in_req->ptr) >> PR_TASK_IDX_SIZE;

  if ((task_id >> (PR_TASK_ID_SIZE - PR_TASK_ID_SCH_BITS)) ==
      context->pr_scheduler_id) {
    ar_assert(!pr_task_collect(task_id));
  }

  else {

    // If we're not the responsible scheduler, forward this upstream
    out_req = noc_msg_send_get_buf(pr_scheduler_core_id(
                                                  context->pr_parent_sched_id));
    out_req->type    = REQ_EXEC_DONE;
    out_req->core_id = context->pr_core_id;
    out_req->req_id  = context->pr_message_id;
    out_req->ptr = in_req->ptr; // copy only field valid for REQ_EXEC_DONE
    ar_assert(!noc_msg_send());

    // Increase message ID
    context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
  }

  // If we are the leaf scheduler responsible for a worker, decrement
  // the child run task load, our total run task load and possibly report
  // upstream.
  if (!context->pr_scheduler_level) {
    child = context->pr_core_child_idx[in_req->core_id];
    ar_assert((child >= 0) && (child < context->pr_num_children));
    ar_assert(context->pr_chld_run_load);
    ar_assert(context->pr_chld_run_load[child] >= 1);
    ar_assert(context->pr_cur_run_load >= 1);
    context->pr_chld_run_load[child]--;
    context->pr_cur_run_load--;
    pr_sched_report_load();
//kt_printf("%d: child %d load decremented to %d\r\n", context->pr_core_id, child, context->pr_chld_run_load[child]);

    // Special case: if we're the ONLY scheduler in the system, it may happen
    // the main task has been collected. We have to begin the shutdown sequence
    // here, because no other schedulers exist to receive a load report change
    // and trigger the shutdown via the event handler.
    if ((context->pr_num_schedulers == 1) && 
        (context->pr_main_finished) &&
        (!context->pr_cur_run_load) && (!context->pr_cur_sched_load)) {

      for (i = 0; i < context->pr_num_children; i++) {

        // Send shutdown request
        out_req = noc_msg_send_get_buf(context->pr_children[i]);
        out_req->type = REQ_SHUTDOWN;
        out_req->req_id = context->pr_message_id;
        ar_assert(!noc_msg_send());

        // Increase message ID
        context->pr_message_id = pr_advance_msg_id(context->pr_message_id);
      }

      // Indicate shutdown to the event handler
      status = 1;
    }
  }

  return status;
}


// ===========================================================================
// pr_event_req_schedule()      Handler for EXT_REQ_SCHEDULE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *in_req           Input request buffer
// ===========================================================================
void pr_event_req_schedule(PrMsgReq *in_req) {

  Context       *context;
  int           locations[PR_REQ_MAX_SIZE];
  int           weights[PR_REQ_MAX_SIZE];
  int           num_locations;
  unsigned int  task_id;
  int           ret;
  int           msg_id;
  int           core_id;
  PrEvtPending  *event;
  PrMsgReply    *out_reply;
  int           i;


  // Get context
  context = mm_get_context(ar_get_core_id());

  // Get request fields
  task_id       = in_req->size;
  num_locations = in_req->num_regions;
  for (i = 0; i < num_locations; i++) {
    locations[i] = ((size_t) in_req->data[i]) >> MM_PACK_SIZE_BITS;
    weights[i]   = ((size_t) in_req->data[i]) & ((1 << MM_PACK_SIZE_BITS) - 1);
  }

  // Call the scheduler
  ret = pr_sched_schedule(task_id, locations, weights, num_locations,
                          &msg_id, &core_id);

  if (ret == 0) {

    // Scheduling done, create a reply with the resulting core ID
    ar_assert((core_id >= 0) && (core_id < context->pr_num_cores));

    out_reply = noc_msg_send_get_buf(in_req->core_id);

    out_reply->type = REPLY_SCHEDULE;
    out_reply->req_id = in_req->req_id;
    out_reply->status = 0;
    out_reply->result = (void *) core_id;

    // Send reply
    ar_assert(!noc_msg_send());
  }

  else if (ret == 1) {

    // Scheduling continues downstream, create an event to wait for the reply
    // to the message the scheduler has sent and forward it back to the
    // requesting scheduler.
    ar_assert(msg_id);
    event = kt_malloc(sizeof(PrEvtPending));
    event->req = kt_malloc(sizeof(PrMsgReq));

    // The only thing that matters from the original request is the requestor
    // core ID and message ID, so the event handler knows where to forward
    // the reply to.
    event->req->core_id = in_req->core_id;
    event->req->req_id  = in_req->req_id;
    event->action = PR_ACT_FORWARD;
    event->prev = NULL;
    event->next = NULL;
    event->data = NULL;

    // Store event; we don't expect conflicts on this message ID.
    ar_assert(!kt_trie_insert(context->pr_pending_events, msg_id, event));
  }
  else {
    ar_abort();
  }


}


// ===========================================================================
// pr_event_reply_get_pages()   Handler for REPLY_GET_PAGES events
// ===========================================================================
// * INPUTS
//   PrMsgReply *in_reply       Input reply buffer
// ===========================================================================
void pr_event_reply_get_pages(PrMsgReply *in_reply) {

  Context               *context;


  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // Insert pages to free pool
  ar_assert(!mm_insert_free_pages((size_t) in_reply->result, in_reply->size));
}


// ===========================================================================
// pr_event_reply_get_rids()    Handler for REPLY_GET_RIDS events
// ===========================================================================
// * INPUTS
//   PrMsgReply *in_reply       Input reply buffer
// ===========================================================================
void pr_event_reply_get_rids(PrMsgReply *in_reply) {

  Context               *context;


  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // Insert region IDs to free pool
  ar_assert(!mm_insert_free_rids((size_t) in_reply->result, in_reply->size));
}


// ===========================================================================
// pr_event_self_ralloc_update_parent()  Handler for SELF_RALLOC_UPDATE_PARENT
//                                       events
// ===========================================================================
// * INPUTS
//   PrMsgReq *stored_req       Stored request buffer
//   PrMsgReply *child_reply    Child that did the RALLOC reply buffer
// ===========================================================================
void pr_event_self_ralloc_update_parent(PrMsgReq *stored_req, 
                                        PrMsgReply *child_reply) {

  PrMsgReply            *out_reply;
  Context               *context;


  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // Update the parent
  ar_assert(!mm_ralloc_update_parent(stored_req, child_reply));

  // Prepare reply. We don't expect any failures of any kind during this call.
  out_reply = noc_msg_send_get_buf(stored_req->core_id);

  out_reply->type = REPLY_RALLOC;
  out_reply->req_id = stored_req->req_id;
  out_reply->status = 0;
  out_reply->result = child_reply->result;

  // Send reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_self_pack_merge()   Handler for SELF_PACK_MERGE events
// ===========================================================================
// * INPUTS
//   PrMsgReq *stored_req       Stored request buffer
//   PrMsgReply *pack_reply     Other scheduler PACK reply buffer
//   PrEvtHookPack *state       Reentrant state structure
// ===========================================================================
void pr_event_self_pack_merge(PrMsgReq *stored_req, PrMsgReply *pack_reply, 
                              PrEvtHookPack *state) {

  Context               *context;
  int                   total;
  int                   i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(state);
  ar_assert(state->wait_replies);
  

  // Was there an error?
  if (pack_reply->status) {

    // A NOTHING_TO_PACK error reply can be normally merged with other stuff
    // we have here, so it's not an error as far as we are concerned. Note
    // that if all replies are empty, then the final reply we'll send to the
    // requestor will regenerate this. So, it's safe to ignore this here.
    if (pack_reply->status == ERR_NOTHING_TO_PACK) {
      total = 0;
    }

    // We may receive multiple errors from multiple scheduler replies.
    // Only keep the first one of them, since we send back a single error
    // object/region.
    else if (!state->error_status) {
      state->error_status = pack_reply->status;
      state->error_ptr = pack_reply->result;
    }
  }

  // If no error (and we have something to pack), merge results with the
  // existing ones in the state
  else {
    total = (pack_reply->num_elements == -1) ? PR_RPL_MAX_SIZE :
                                               pack_reply->num_elements;

    // Reallocate arrays, if needed
    if (state->num_elements + total >= state->alloc_elements) {
      state->alloc_elements = state->num_elements + total;
      ar_assert(state->sizes = kt_realloc(state->sizes,
                                     state->alloc_elements * sizeof(int)));
      ar_assert(state->addresses = kt_realloc(state->addresses,
                                     state->alloc_elements * sizeof(size_t)));
    }

    // Copy stuff
    for (i = 0; i < total; i++) {
      (state->addresses)[state->num_elements + i] = (size_t) 
                                                        pack_reply->adr_ar[i];
      (state->sizes)[state->num_elements + i] = pack_reply->size_ar[i];
    }
    state->num_elements += total;
  }

  // Extended replies may each one be in multiple pieces. If this is the last
  // part of the reply from the given scheduler, note that we wait for one
  // less scheduler reply.
  if (pack_reply->num_elements != -1) {
    state->wait_replies--;
  }

  // If all replies were gathered at last, reply to original requestor
  if (!state->wait_replies) {

    // Is this a native call?
    if (state->native_task_id) {

      // Start scheduling, the packing for this task is complete
      pr_sched_start_scheduling(state);

      // Free the state, but don't free state->sizes and state->addresses:
      // pr_sched_start_scheduling() will prune & graft them on the task
      // descriptor.
      kt_free(state);
    }
    else {

      // Send reply based on what's gathered in the state
      pr_event_reply_pack(stored_req->core_id, stored_req->req_id, state);
      
      // Free the state, including state->sizes and state->addresses
      kt_free(state->sizes);
      kt_free(state->addresses);
      kt_free(state);
    }

  }
}


// ===========================================================================
// pr_event_self_rfree_wait_children()  Handler for SELF_RFREE_WAIT_CHILDREN
//                                      events
// ===========================================================================
// * INPUTS
//   PrMsgReq *stored_req       Stored request buffer
//   PrMsgReply *child_reply    Child that did the RFREE reply buffer
//   PrEvtHookRfree *state      Reentrant state structure
// ===========================================================================
void pr_event_self_rfree_wait_children(PrMsgReq *stored_req, 
                                       PrMsgReply *child_reply, 
                                       PrEvtHookRfree *state) {

  Context               *context;
  PrMsgReply            *out_reply;
  int                   ret;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(state);
  ar_assert(state->wait_replies);
  
  // Was there an error?
  if (child_reply->status) {
    // There shouldn't be. It's a child region ID we attempted to free, so
    // it's supposed to be sane.
    ar_abort();
  }

  // One more down
  state->wait_replies--;

  // Do we wait for others?
  if (state->wait_replies) {
    return;
  }

  // We're done. Free state here.
  kt_free(state);

  // Redo the region free. This may succeed or it may again be postponed. The
  // first time we did the region free (and got the ERR_REMOTE_CHILDREN
  // reply), the function disassociated these children IDs from the local
  // parent, so all that's left to do now is actually free the parent (and any
  // local children). However, a remote parent may still exist and will need
  // to be updated, which will create a new postponed reply...
  ret = mm_rfree(stored_req);
  if (ret == ERR_REPLY_POSTPONED) {
    return;
  }

  // ... but we don't expect any other error here.
  ar_assert(!ret);

  // Prepare reply
  out_reply = noc_msg_send_get_buf(stored_req->core_id);

  out_reply->type = REPLY_RFREE;
  out_reply->req_id = stored_req->req_id;
  out_reply->status = 0;

  // Send reply
  ar_assert(!noc_msg_send());
}


// ===========================================================================
// pr_event_self_schedule_result()      Handler for SELF_SCHEDULE_RESULT
//                                      events
// ===========================================================================
// * INPUTS
//   PrMsgReq *stored_req       Stored request buffer
//   PrMsgReply *child_reply    Child that did the RALLOC reply buffer
// ===========================================================================
void pr_event_self_schedule_result(PrMsgReq *stored_req, 
                                   PrMsgReply *child_reply) {

  Context               *context;
  PrTaskDescr           *task;

  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // Update the task descriptor
  task = stored_req->ptr;
  ar_assert(task->run_core_id == -1);
  ar_assert(child_reply->type == REPLY_SCHEDULE);
  task->run_core_id = (size_t) child_reply->result;

  // Task is scheduled, proceed with location update and dispatch
  pr_task_scheduled(task);
}


// ===========================================================================
// pr_event_self_wait_spawn()   Handler for SELF_WAIT_SPAWN events
// ===========================================================================
void pr_event_self_wait_spawn() {

  Context               *context;

  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // Update the spawn status
  ar_assert(context->pr_spawn_pending);
  context->pr_spawn_pending = 0;
}



// ###########################################################################
// ###                                                                     ###
// ###                          EXPORTED FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// pr_event_wake_up_pending()   Looks at an incoming reply message ID number
//                              and wakes up all events that wait for this
//                              reply. For each of them, does the required
//                              action by calling the process handler
//                              accordingly.
// ===========================================================================
// * INPUTS
//   PrMsgReply *in_reply       Input reply buffer
// ===========================================================================
void pr_event_wake_up_pending(PrMsgReply *in_reply) {

  Context               *context;
  PrMsgReply            *out_reply;
  PrEvtPending          *event;
  PrEvtPending          *del_event;
  PrEvtPending          *keep_head;
  PrEvtPending          *keep_tail;

  
  // Get context
  context = mm_get_context(ar_get_core_id());

  // Look for pending events that wait reply for this message ID
  kt_trie_find(context->pr_pending_events, in_reply->req_id, (void *) &event);
  if (!event) {
    return;
  }

  // For all events in the linked list
  keep_head = NULL;
  keep_tail = NULL;
  while (event) {

    ar_assert(event->req);

    // What should we do?
    switch (event->action) {
      
      // =====================================================================
      // Forward the reply to original requestor
      // =====================================================================
      case PR_ACT_FORWARD:

        // Prepare reply
        out_reply = noc_msg_send_get_buf(event->req->core_id);

        // Extended reply?
        if (pr_event_is_extended(in_reply->type)) {

          // Copy all the fields from incoming reply
          kt_memcpy(out_reply, in_reply, NOC_MESSAGE_SIZE);
        }
        else {

          // Copy only basic fields from incoming reply
          out_reply->type   = in_reply->type;
          out_reply->result = in_reply->result;
          out_reply->status = in_reply->status;
          out_reply->size   = in_reply->size;
        }

        // Copy (basic) or overwrite (extended) reply msg ID with the original
        // request msg ID
        out_reply->req_id = event->req->req_id;
        
        // Reply to original requestor
        ar_assert(!noc_msg_send());
        break;


      // =====================================================================
      // Reprocess the pending event from the beginning, or using saved state
      // =====================================================================
      case PR_ACT_REDO:
      case PR_ACT_REENTRANT:
      case PR_ACT_REDO_ERRORS:
      case PR_ACT_REENTRANT_ERRORS:

        // Is reply successful, or do we also want to handle errors?
        if ((!in_reply->status) ||
            (event->action == PR_ACT_REDO_ERRORS) ||
            (event->action == PR_ACT_REENTRANT_ERRORS)) {
          
          if ((event->action == PR_ACT_REDO) ||
              (event->action == PR_ACT_REDO_ERRORS)) {

            ar_assert(!event->data);


            // Redo the pending event. The reply itself is also given to the
            // event dispatcher as the auxiliary input buffer, because some
            // events will need to get fields from it.
            pr_event_process_message(event->req, in_reply, NULL);
          }

          else if ((event->action == PR_ACT_REENTRANT) ||
                   (event->action == PR_ACT_REENTRANT_ERRORS)) {

            ar_assert(event->data);

            // Continue the pending event from its saved state. The reply
            // itself is also given to the event dispatcher as the auxiliary
            // input buffer, because some events will need to get fields from
            // it.
            pr_event_process_message(event->req, in_reply, event->data);
          }
        }

        // If unsuccessful, just create a failed reply of the correct type
        // and forward the failed status in it.
        else {

          // Prepare reply
          out_reply = noc_msg_send_get_buf(event->req->core_id);

          // Outgoing reply msg ID is the original request msg ID
          out_reply->req_id = event->req->req_id;

          // Determine correct reply type
          out_reply->type = pr_event_reply_type(event->req->type);
          
          // Copy failed status from incoming reply
          out_reply->status = in_reply->status;
          out_reply->result = NULL;
          out_reply->size   = 0;
          
          // Reply to original requestor
          ar_assert(!noc_msg_send());
        }
        break;


      // Unknown action
      default:
        ar_abort();
    }


    // We must not delete the event, if it's a successful (or we want to
    // handle errors), extended reply which continues in more than 1 pieces
    // and this is not the last piece.  Maintain a "keep list" of events that
    // will remain undeleted.
    if (((!in_reply->status) ||
         (event->action == PR_ACT_REDO_ERRORS) ||
         (event->action == PR_ACT_REENTRANT_ERRORS)) &&
        (pr_event_is_extended(in_reply->type)) && 
        (in_reply->num_elements == -1)) {
      if (!keep_head) {
        ar_assert(!keep_tail);
        keep_head = event;
        event->prev = NULL;
      }
      else {
        ar_assert(keep_tail);
        keep_tail->next = event;
        event->prev = keep_tail;
      }
      keep_tail = event;
      event = event->next;
      keep_tail->next = NULL;
    }

    // Otherwise, delete it and move on to the next in list. We don't touch
    // the event's data hook; the process handler must have dealt with it
    // (it may have migrated it to a new event, or freed it, or whatever).
    else {
      del_event = event;
      event = event->next;
      kt_free(del_event->req);
      kt_free(del_event);
    }
  }

  // Delete/replace msg ID entry, depending on if there are any events left
  ar_assert(!kt_trie_delete(context->pr_pending_events, in_reply->req_id,
                            NULL));
  if (keep_head) {
    ar_assert(!kt_trie_insert(context->pr_pending_events, in_reply->req_id,
                              keep_head));
  }
}


// ===========================================================================
// pr_event_process_message()   Event dispatcher. Takes a generic incoming
//                              (or stored) message buffer, finds out what
//                              its event type is and calls the appropriate
//                              handler to handle it.
// ===========================================================================
// * INPUTS
//   void *main_in_buf          Main input request/reply buffer, which
//                              contains the event we must process
//   void *aux_in_buf           Secondary input request/reply buffer, which
//                              (if not NULL) contains a secondary message
//                              that will be needed by some process handlers
//                              (e.g. like intermediate replies to reentrant
//                              events)
//   void *saved_state          If not NULL, saved state for reentrant events
//                              that will be passed to their handlers
//
// * RETURN VALUE
//   int                        0: success
//                              1: message was REQ_SHUTDOWN, need to exit
// ===========================================================================
int pr_event_process_message(void *main_in_buf, void *aux_in_buf,
                             void *saved_state) {

  Context               *context;
  PrMsgReq              *main_in_req;
  PrMsgReply            *main_in_reply;
  PrMsgReq              *aux_in_req;
  PrMsgReply            *aux_in_reply;
  int                   tmp;


  // Get context
  context = mm_get_context(ar_get_core_id());
  
  // Assign basic req/reply pointers to input buffers
  main_in_req = (PrMsgReq *) main_in_buf;
  main_in_reply = (PrMsgReply *) main_in_buf;
  aux_in_req = (PrMsgReq *) aux_in_buf;
  aux_in_reply = (PrMsgReply *) aux_in_buf;
  
  // Call appropriate handler. We find the type as if it were a request:
  // requests and replies share the same first field.
  switch (main_in_req->type) {


    // ====================================================================
    // Part 1: Requests, from other schedulers or workers. Processing
    //         functions may create pending events for later processing,
    //         when replies to outgoing requests are received.
    // ====================================================================

    // A worker or scheduler tells us to allocate an object
    case REQ_ALLOC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_ALLOC_BEGIN);
      pr_event_req_alloc(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_ALLOC_END);
      break;


    // A worker or scheduler tells us to allocate multiple objects
    case REQ_BALLOC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_BALLOC_BEGIN);
      pr_event_req_balloc(main_in_req, saved_state);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_BALLOC_END);
      break;


    // A worker or scheduler tells us to free an object
    case REQ_FREE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_FREE_BEGIN);
      pr_event_req_free(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_FREE_END);
      break;


    // A worker or scheduler tells us to allocate a region
    case REQ_RALLOC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_RALLOC_BEGIN);
      pr_event_req_ralloc(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_RALLOC_END);
      break;


    // A scheduler tells us to allocate a region without a parent, but
    // with a remote parent ID
    case REQ_RALLOC_ORPHAN:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_RALLOC_ORPHAN_BEGIN);
      pr_event_req_ralloc_orphan(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_RALLOC_ORPHAN_END);
      break;


    // A worker or scheduler tells us to free a region
    case REQ_RFREE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_RFREE_BEGIN);
      pr_event_req_rfree(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_RFREE_END);
      break;


    // A scheduler tells us to update one of our regions, because it's a
    // parent to a just deleted region
    case REQ_RFREE_UPDATE_PARENT:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_RFREE_UPDATE_PARENT_BEGIN);
      pr_event_req_rfree_update_parent(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_RFREE_UPDATE_PARENT_END);
      break;


    // A worker or scheduler asks us to pack objects and/or regions
    case REQ_PACK:
    case EXT_REQ_PACK:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_PACK_BEGIN);
      pr_event_req_pack(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_PACK_END);
      break;


    // A worker or scheduler asks us about an object
    case REQ_QUERY_POINTER:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_QUERY_POINTER_BEGIN);
      pr_event_req_query_pointer(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_QUERY_POINTER_END);
      break;


    // A scheduler asks us to provide free pages
    case REQ_GET_PAGES:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_GET_PAGES_BEGIN);
      pr_event_req_get_pages(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_GET_PAGES_END);
      break;


    // A scheduler gives us free pages, because he has too many
    case REQ_RETURN_PAGES:
      // FIXME
      ar_abort();
      break;


    // A scheduler asks us to provide free region IDs
    case REQ_GET_RIDS:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_GET_RIDS_BEGIN);
      pr_event_req_get_rids(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REQ_GET_RIDS_END);
      break;


    // A scheduler gives us free region IDs, because he has too many
    case REQ_RETURN_RIDS:
      // FIXME
      ar_abort();
      break;


    // A child scheduler reports his load to us
    case REQ_LOAD_REPORT:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_LOAD_REPORT_BEGIN);
      tmp = pr_event_req_load_report(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_LOAD_REPORT_END);

      // Should we exit?
      if (tmp) {
        return 1;
      }
      break;


    // A (higher-level) scheduler tells us to shut down
    case REQ_SHUTDOWN:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_SHUTDOWN_BEGIN);
      pr_event_req_shutdown(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_SHUTDOWN_END);

//if (context->pr_role == ROLE_SCHEDULER) kt_printf("%d: Got shutdown, my load = %d/%d/%d\r\n", context->pr_core_id, context->mm_current_load, context->pr_cur_run_load, context->pr_cur_sched_load);

      // Return with an indication to exit
      return 1;


    // A worker or scheduler requests from us to create a new task
    case EXT_REQ_SPAWN:
      // Are we the final destination? Call the appropriate handler.
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_SPAWN_BEGIN);
      if ((((size_t) main_in_req->ptr) >> (PR_TASK_IDX_SIZE + PR_TASK_ID_SIZE - 
                                           PR_TASK_ID_SCH_BITS)) == 
          context->pr_scheduler_id) {
        pr_event_req_spawn_final(main_in_req);
      }
      else {
        pr_event_req_spawn_forward(main_in_req);
      }
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_SPAWN_END);
      break;


    // A parent scheduler delegates to us the creation of a task
    case EXT_REQ_DELEGATE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DELEGATE_BEGIN);
      pr_event_req_spawn_final(main_in_req); // same handler as spawn_final
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DELEGATE_END);
      break;


    // A parent scheduler tells us to start analysis on some task arguments
    case EXT_REQ_DEP_START:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DEP_START_BEGIN);
      pr_event_req_dep_start(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DEP_START_END);
      break;


    // A child scheduler informs us what's the region IDs we have to traverse
    // down the region tree in order to reach a target object/region
    case EXT_REQ_DEP_ROUTE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DEP_ROUTE_BEGIN);
      pr_event_req_dep_route(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DEP_ROUTE_END);
      break;


    // A parent scheduler tells us to descend the region tree to enqueue a
    // dependency for a new task
    case EXT_REQ_DEP_ENQUEUE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DEP_ENQUEUE_BEGIN);
      pr_event_req_dep_enqueue(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DEP_ENQUEUE_END);
      break;


    // A child scheduler informs us that a dependency we're waiting for a 
    // new task is now at the head of its queue
    case REQ_DEP_OK:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DEP_OK_BEGIN);
      pr_event_req_dep_ok(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DEP_OK_END);
      break;


    // A parent scheduler tells us to go to the next task waiting for a number
    // of task arguments
    case EXT_REQ_DEP_STOP:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DEP_STOP_BEGIN);
      pr_event_req_dep_stop(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DEP_STOP_END);
      break;


    // A parent scheduler tells us to go to the next task waiting for a number
    // of task arguments
    case REQ_DEP_CHILD_FREE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_DEP_CHILD_FREE_BEGIN);
      pr_event_req_dep_child_free(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_DEP_CHILD_FREE_END);
      break;


    // A parent scheduler tells us to update the location of some arguments
    case EXT_REQ_UPDATE_LOCATION:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_UPDATE_LOCATION_BEGIN);
      pr_event_req_dep_update_loc(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_UPDATE_LOCATION_END);
      break;


    // A scheduler requests from us to execute a new task
    case EXT_REQ_EXEC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_EXEC_BEGIN);
      if (context->pr_role == ROLE_SCHEDULER) {
        pr_event_req_exec_scheduler(main_in_req);
        dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      }
      else {
        pr_event_req_exec_worker(main_in_req);
        dbg_stime(context, DBG_STATS_IDX_TIME_TASK_EXEC);
      }
      dbg_trace(context, DBG_TRC_REQ_EXEC_END);
      break;

    // A worker finished a task it was assigned to it for execution
    case REQ_EXEC_DONE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_EXEC_DONE_BEGIN);
      tmp = pr_event_req_exec_done(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_EXEC_DONE_END);

      // Should we exit?
      if (tmp) {
        return 1;
      }
      break;

    // A parent scheduler tells us to continue the scheduling of a task
    case EXT_REQ_SCHEDULE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REQ_SCHEDULE_BEGIN);
      pr_event_req_schedule(main_in_req);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REQ_SCHEDULE_END);
      break;



    // ====================================================================
    // Part 2: Requests, from ourselves. These are actually received via
    //         the event wake-up system, and are like "notes to self" to
    //         do actions that depend on replies from other schedulers.
    // ====================================================================

    // A reply for an orphan region has arrived. Update the local parent with
    // its new remote child ID.
    case SELF_RALLOC_UPDATE_PARENT:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_SELF_RALLOC_UPDATE_PARENT_BEGIN);
      pr_event_self_ralloc_update_parent(main_in_req, aux_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_SELF_RALLOC_UPDATE_PARENT_END);
      break;


    // A reply from a remote pack operation has arrived. Merge its results
    // with the existing ones.
    case SELF_PACK_MERGE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_SELF_PACK_MERGE_BEGIN);
      pr_event_self_pack_merge(main_in_req, aux_in_reply, saved_state);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_SELF_PACK_MERGE_END);
      break;


    // A reply from a region free has arrived. Decrement pending counter,
    // and if no others are pending redo the local region free.
    case SELF_RFREE_WAIT_CHILDREN:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_SELF_RFREE_WAIT_CHILDREN_BEGIN);
      pr_event_self_rfree_wait_children(main_in_req, aux_in_reply, 
                                        saved_state);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_SELF_RFREE_WAIT_CHILDREN_END);
      break;


    // A final reply from a scheduling decision has arrived. Update the 
    // local task descriptor and proceed to dispatch the task.
    case SELF_SCHEDULE_RESULT:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_SELF_SCHEDULE_RESULT_BEGIN);
      pr_event_self_schedule_result(main_in_req, aux_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_SELF_SCHEDULE_RESULT_END);
      break;


    // A reply from the parent task scheduler arrived that our last spawn
    // request has been completed.
    case SELF_WAIT_SPAWN:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_SELF_WAIT_SPAWN_BEGIN);
      pr_event_self_wait_spawn();
      dbg_stime(context, DBG_STATS_IDX_TIME_WORKER_WAIT);
      dbg_trace(context, DBG_TRC_SELF_WAIT_SPAWN_END);
      break;




    // ====================================================================
    // Part 3: Replies, from other schedulers or workers. These may trigger 
    //         pending events that wait for them to wake up.
    // ====================================================================

    // Another scheduler allocated an object for us
    case REPLY_ALLOC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_ALLOC_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_ALLOC_END);
      break;


    // Another scheduler allocated multiple objects for us
    case EXT_REPLY_BALLOC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_BALLOC_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_BALLOC_END);
      break;
   

    // Another scheduler freed an object for us
    case REPLY_FREE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_FREE_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_FREE_END);
      break;


    // Another scheduler allocated a region for us
    case REPLY_RALLOC:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_RALLOC_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_RALLOC_END);
      break;


    // Another scheduler freed a region for us, or he just updated his parent
    // region to reflect our child region deletion
    case REPLY_RFREE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_RFREE_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_RFREE_END);
      break;


    // Another scheduler queried an object for us
    case REPLY_QUERY_POINTER:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_QUERY_POINTER_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_QUERY_POINTER_END);
      break;


    // Another scheduler packed region(s) for us
    case EXT_REPLY_PACK:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_PACK_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_PACK_END);
      break;


    // Another scheduler replied to our request for free pages
    case REPLY_GET_PAGES:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_GET_PAGES_BEGIN);

      // Mark that we no longer wait for pages
      ar_assert(context->pr_pages_msg_id == main_in_reply->req_id);
      context->pr_pages_msg_id = 0;
      
      // See what happened
      if (!main_in_reply->status) {
        pr_event_reply_get_pages(main_in_reply);
      }
      pr_event_wake_up_pending(main_in_reply);

      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_GET_PAGES_END);
      break;


    // Another scheduler replied to our request, when we offered free pages
    case REPLY_RETURN_PAGES:
      // FIXME
      ar_abort();
      break;


    // Another scheduler replied to our request for free region IDs
    case REPLY_GET_RIDS:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_GET_RIDS_BEGIN);

      // Mark that we no longer wait for IDs
      ar_assert(context->pr_rids_msg_id == main_in_reply->req_id);
      context->pr_rids_msg_id = 0;
      
      // See what happened
      if (!main_in_reply->status) {
        pr_event_reply_get_rids(main_in_reply);
      }
      pr_event_wake_up_pending(main_in_reply);
      
      dbg_stime(context, DBG_STATS_IDX_TIME_MEM_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_GET_RIDS_END);
      break;


    // Another scheduler replied to our request, when we offered free region
    // IDs
    case REPLY_RETURN_RIDS:
      // FIXME
      ar_abort();
      break;

    // A child scheduler completed the scheduling we had requested
    case REPLY_SCHEDULE:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_SCHEDULE_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_SCHEDULE_END);
      break;

    // A parent scheduler completed the spawn we had requested
    case REPLY_SPAWN:
      dbg_stime(context, DBG_STATS_IDX_TIME_IDLE);
      dbg_trace(context, DBG_TRC_REPLY_SPAWN_BEGIN);
      pr_event_wake_up_pending(main_in_reply);
      dbg_stime(context, DBG_STATS_IDX_TIME_SCH_SERVE);
      dbg_trace(context, DBG_TRC_REPLY_SPAWN_END);
      break;



    // ====================================================================
    // Unknown request or reply type
    // ====================================================================
    default:
      ar_abort();
  }

  // Return with normal status
  return 0;
}


// ===========================================================================
// pr_event_main_loop()         Scheduler and worker main loop. It blocks
//                              waiting for requests from any source. When it
//                              finds a new one, it dequeues, processes it and
//                              replies to the requesting core.
// ===========================================================================
void pr_event_main_loop() {

  Context       *context;
  void          *in_buf;
  int           shutdown;


  // Type-related sanity checks
  ar_assert(8 * sizeof(size_t) > 2 * PR_REQ_MAX_SIZE); // EXT_REQ_SPAWN bitmaps
  ar_assert(8 * sizeof(rid_t) > 2 * PR_REQ_MAX_SIZE);  // EXT_REQ_SPAWN bitmaps
  ar_assert((2 + PR_REQ_MAX_SIZE) * MM_PACK_OPTION_BITS <= 32); // REQ_PACK and
                                                        // EXT_REQ_PACK bitmaps
  ar_assert(PR_TASK_ID_SIZE + PR_TASK_IDX_SIZE <= 32); // task ID + index
  ar_assert(MM_PACK_SIZE_BITS + MM_PACK_LOCATION_BITS + 
            MM_PACK_OPTION_BITS == 32); // pack size field

  // Initialize pending events trie
  context = mm_get_context(ar_get_core_id());
  ar_assert(!context->pr_pending_events);
  ar_assert(context->pr_pending_events = kt_alloc_trie(PR_TRIES_MSG_ID_MSB, 0));

  // Initialize multi-part requests incomplete list
  ar_assert(!context->pr_incomplete_req);
  ar_assert(context->pr_incomplete_req = kt_alloc_list());


  // Loop until REQ_SHUTDOWN is received
  shutdown = 0;
  while (!shutdown) {

    // Block until any request/reply message arrives
    ar_assert(!noc_msg_recv(1, &in_buf));
    
    // Process it
    shutdown = pr_event_process_message(in_buf, NULL, NULL);
  }
}


// ===========================================================================
// pr_event_worker_inner_loop()  Worker secondary loop. While a worker is
//                               executing a task, it may want to receive
//                               a certain reply. We do that in a loop,
//                               processing other events normally. In the 
//                               meantime, if active DMAs exist, their 
//                               progress is monitored so that failed DMAs
//                               can be restarted.
// ===========================================================================
// * INPUTS
//   int wait_mode               0: Wait for a message matching type and 
//                                  msg_id, return it and process all else 
//                               1: Wait until previous spawn requests have
//                                  been acknowledged
//                               2: Wait until DMAs for task are completed
//
// * RETURN VALUE
//   void *                      The received matching message for wait_mode 
//                               0, NULL in other cases
// ===========================================================================
void *pr_event_worker_inner_loop(int wait_mode, PrMsgType type, int msg_id,
                                 PrTaskDescr *task) {
  
  Context       *context;
  PrMsgReq      *in_buf;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());

  dbg_stime(context, DBG_STATS_IDX_TIME_TASK_EXEC);
  dbg_trace(context, DBG_TRC_REQ_EXEC_END);
  dbg_trace(context, DBG_TRC_WORKER_WAIT_BEGIN);

  while (1) {


    // Do we have any pending DMAs?
    if (kt_list_head(context->noc_active_dmas)) {

      // See if there's any progress with them and restart failed ones
      noc_dma_check_progress(NULL);

      // Were the DMAs we're waiting for completed?
      if ((wait_mode == 2) && (!task->dmas_waiting)) {
        dbg_stime(context, DBG_STATS_IDX_TIME_WORKER_WAIT);
        dbg_trace(context, DBG_TRC_WORKER_WAIT_END);
        dbg_trace(context, DBG_TRC_REQ_EXEC_BEGIN);
        return NULL;
      }

//if(wait_mode==2){kt_printf("%d: DMA wait task 0x%08X var 0x%08X=%d\r\n", context->pr_core_id, task->id, &task->dmas_waiting, task->dmas_waiting);}

      // See if there's any new request/reply message, but don't block
      noc_msg_recv(0, (void *) &in_buf);
    }
    else {
      // Block until any request/reply message arrives
      ar_assert(!noc_msg_recv(1, (void *) &in_buf));
    }

    // Found the one we're waiting for?
    if ((wait_mode == 0) && (in_buf) && 
        (in_buf->type == type) && (in_buf->req_id == msg_id)) {
      dbg_stime(context, DBG_STATS_IDX_TIME_WORKER_WAIT);
      dbg_trace(context, DBG_TRC_WORKER_WAIT_END);
      dbg_trace(context, DBG_TRC_REQ_EXEC_BEGIN);
      return in_buf;
    }

    // Process it
    if (in_buf) {
      dbg_stime(context, DBG_STATS_IDX_TIME_WORKER_WAIT);
      dbg_trace(context, DBG_TRC_WORKER_WAIT_END);

      ar_assert(!pr_event_process_message(in_buf, NULL, NULL));

      // Was it a spawn reply?
      if ((wait_mode == 1) && (!context->pr_spawn_pending)) {
        dbg_trace(context, DBG_TRC_REQ_EXEC_BEGIN);
        return NULL;
      }

      dbg_trace(context, DBG_TRC_WORKER_WAIT_BEGIN);
    }
  }
}


// ###########################################################################
// ###                                                                     ###
// ###                         DEBUGGING FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// pr_event_dbg_type_string()   Enumeration-to-string for message type
// ===========================================================================
// * INPUTS
//   PrMsgType type             The message type
//
// * RETURN VALUE
//   char *                     Its string representation
// ===========================================================================
char *pr_event_dbg_type_string(PrMsgType type) {

  switch (type) {
    case REQ_ALLOC:                     return "REQ_ALLOC";
    case REQ_BALLOC:                    return "REQ_BALLOC";
    case REQ_FREE:                      return "REQ_FREE";
    case REQ_RALLOC:                    return "REQ_RALLOC";
    case REQ_RALLOC_ORPHAN:             return "REQ_RALLOC_ORPHAN";
    case REQ_RFREE:                     return "REQ_RFREE";
    case REQ_RFREE_UPDATE_PARENT:       return "REQ_RFREE_UPDATE_PARENT";
    case REQ_PACK:                      return "REQ_PACK";
    case EXT_REQ_PACK:                  return "EXT_REQ_PACK";
    case REQ_QUERY_POINTER:             return "REQ_QUERY_POINTER";
    case REQ_GET_PAGES:                 return "REQ_GET_PAGES";
    case REQ_RETURN_PAGES:              return "REQ_RETURN_PAGES";
    case REQ_GET_RIDS:                  return "REQ_GET_RIDS";
    case REQ_RETURN_RIDS:               return "REQ_RETURN_RIDS";
    case REQ_LOAD_REPORT:               return "REQ_LOAD_REPORT";
    case REQ_SHUTDOWN:                  return "REQ_SHUTDOWN";
    case EXT_REQ_SPAWN:                 return "EXT_REQ_SPAWN";
    case EXT_REQ_DELEGATE:              return "EXT_REQ_DELEGATE";
    case EXT_REQ_DEP_START:             return "EXT_REQ_DEP_START";
    case EXT_REQ_DEP_ROUTE:             return "EXT_REQ_DEP_ROUTE";
    case EXT_REQ_DEP_ENQUEUE:           return "EXT_REQ_DEP_ENQUEUE";
    case REQ_DEP_OK:                    return "REQ_DEP_OK";
    case EXT_REQ_DEP_STOP:              return "EXT_REQ_DEP_STOP";
    case REQ_DEP_CHILD_FREE:            return "REQ_DEP_CHILD_FREE";
    case EXT_REQ_UPDATE_LOCATION:       return "EXT_REQ_UPDATE_LOCATION";
    case EXT_REQ_EXEC:                  return "EXT_REQ_EXEC";
    case REQ_EXEC_DONE:                 return "REQ_EXEC_DONE";
    case EXT_REQ_SCHEDULE:              return "EXT_REQ_SCHEDULE";
    case SELF_RALLOC_UPDATE_PARENT:     return "SELF_RALLOC_UPDATE_PARENT";
    case SELF_PACK_MERGE:               return "SELF_PACK_MERGE";
    case SELF_RFREE_WAIT_CHILDREN:      return "SELF_RFREE_WAIT_CHILDREN";
    case SELF_SCHEDULE_RESULT:          return "SELF_SCHEDULE_RESULT";
    case SELF_WAIT_SPAWN:               return "SELF_WAIT_SPAWN";
    case REPLY_ALLOC:                   return "REPLY_ALLOC";
    case EXT_REPLY_BALLOC:              return "EXT_REPLY_BALLOC";
    case REPLY_FREE:                    return "REPLY_FREE";
    case REPLY_RALLOC:                  return "REPLY_RALLOC";
    case REPLY_RFREE:                   return "REPLY_RFREE";
    case EXT_REPLY_PACK:                return "EXT_REPLY_PACK";
    case REPLY_QUERY_POINTER:           return "REPLY_QUERY_POINTER";
    case REPLY_GET_PAGES:               return "REPLY_GET_PAGES";
    case REPLY_RETURN_PAGES:            return "REPLY_RETURN_PAGES";
    case REPLY_GET_RIDS:                return "REPLY_GET_RIDS";
    case REPLY_RETURN_RIDS:             return "REPLY_RETURN_RIDS";
    case REPLY_SCHEDULE:                return "REPLY_SCHEDULE";
    case REPLY_SPAWN:                   return "REPLY_SPAWN";
    default:                            return "!! ERROR !!";
  }
}
