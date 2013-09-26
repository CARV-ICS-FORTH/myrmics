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
// Abstract      : Functionality related to local memory regions handling.
//                 Local region tree, local user free pool and region/object
//                 allocation/destruction are done here.
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: region.c,v $
// CVS revision  : $Revision: 1.13 $
// Last modified : $Date: 2013/04/09 14:54:35 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>


// ###########################################################################
// ###                                                                     ###
// ###                          INTERNAL FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################


// ===========================================================================
// mm_region_create_node()      Creates a new region tree node, links it
//                              to parent (if parent is local) and updates the
//                              local tries.
// ===========================================================================
// * INPUTS
//   rid_t parent_id            The parent region ID, or 0 if the region tree
//                              root (NULL region) is being created
//   rid_t rid                  Region ID to be given to new node
//   size_t start_adr           Start of free chunk to be given to region
//   int num_slabs              Number of free slabs to be given to region
//   int location               Core ID of the initial region owner
//
// * OUTPUTS
//   MmRgnTreeNode *ret_node    The newly created region tree node
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_region_create_node(rid_t parent_id, rid_t rid, size_t start_adr,
                          int num_slabs, int location,
                          MmRgnTreeNode **ret_node) {

  Context       *context;
  MmRgnTreeNode *node;
  MmRgnAdrRange *adr_range;
  MmRgnTreeNode *parent;


  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(ret_node);

  // Create region tree node
  ar_assert(node = kt_zalloc(sizeof(MmRgnTreeNode)));

  // Create slab pool with these addresses. Currently, all metadata for new
  // regions are kept in the kernel. If region migration is needed in the
  // future, a secondary slab pool must be created and kept in the node.
  ar_assert(node->pool = mm_slab_create_pool(context->mm_kernel_pool,
                                             start_adr, num_slabs, 0));

  // Assign region ID
  node->id = rid;
  ar_assert(!kt_trie_insert(context->mm_local_rids, rid, node));

  // Link to parent
  node->parent_id = parent_id;

  // If parent is local, be a part of his children
  if (parent_id) {
    kt_trie_find(context->mm_local_rids, parent_id, (void *) &parent);
    if (parent) {
      ar_assert(parent->children_ids = kt_realloc(parent->children_ids,
                                              (parent->num_children + 1) *
                                                  sizeof(rid_t)));
      parent->children_ids[parent->num_children++] = rid;
    }
  }

  // Update used ranges to point to the new region. It's a new region, so the
  // address range must not be merged with anything else.
  ar_assert(mm_region_add_adr_range(context->mm_used_ranges, node,
                                    context->pr_scheduler_id,
                                    start_adr, num_slabs, &adr_range) == 2);
  context->mm_free_num_slabs -= num_slabs;
  ar_assert(context->mm_free_num_slabs >= 0);

  // Update used region IDs to point to the new region
  mm_region_add_rid_range(context->mm_used_rids, context->pr_scheduler_id,
                          rid, 1, NULL);
  context->mm_free_num_rids--;
  ar_assert(context->mm_free_num_rids >= 0);

  // Initialize region tree node array of used ranges
  ar_assert(node->ranges = kt_malloc(sizeof(MmRgnAdrRange *)));
  node->ranges[0] = adr_range;
  node->num_ranges = 1;

  // Initialize dependency structures
  node->dep_queue = kt_alloc_list();
  node->obj_dep_queues = kt_alloc_trie(MM_TRIES_SLOT_MSB, MM_TRIES_SLOT_LSB);
  ar_assert(!node->children_enqs_ro);
  ar_assert(!node->children_enqs_rw);
  ar_assert(!node->obj_enqs_rw);
  ar_assert(!node->obj_enqs_ro);
  ar_assert(!node->parent_enqs_ro);
  ar_assert(!node->parent_enqs_rw);

  // Initialize location structures
  node->location = location;
  node->obj_locations = kt_alloc_trie(MM_TRIES_SLOT_MSB, MM_TRIES_SLOT_LSB);

  // Success
  *ret_node = node;
  return 0;
}


// ===========================================================================
// mm_region_free_node()        Destroys a childless region tree node, after
//                              massively freeing all its slabs and reclaiming
//                              them in the local free pool.
//
//                              It is possible that various pieces of the
//                              address ranges used by the region are left
//                              in the caches of the last producers who
//                              touched these ranges. It is also possible that
//                              after we free this region, the address ranges
//                              will be reused for other allocations. Thus,
//                              before desctruction, we pack this region and
//                              create DMAs to all known producers to transfer
//                              their ranges from the caches to the DRAM. We
//                              have to wait these DMAs before returning,
//                              otherwise we introduce a race condition with
//                              new allocation calls.
// ===========================================================================
// * INOUTS
//   MmRgnTreeNode *node        The region tree node to be destroyed. It must
//                              be already childless.
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_region_free_node(MmRgnTreeNode *node) {

  Context       *context;
  MmRgnRidRange *rid_range;
  rid_t         *err_children;
  int           num_err_children;
  int           alloc_num_elements;
  int           num_elements;
  size_t        *addresses;
  int           *sizes;
  int           dmas_waiting;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(node);
  ar_assert(node->pool);
  ar_assert(!node->children_ids);
  ar_assert(!node->num_children);


  // Pack the region
  alloc_num_elements = 32;
  sizes = kt_malloc(alloc_num_elements * sizeof(int));
  addresses = kt_malloc(alloc_num_elements * sizeof(size_t));
  num_elements = 0;

  ar_assert(!mm_region_pack_region(node->id, 0, 0, &alloc_num_elements,
                                   &num_elements, &addresses, &sizes,
                                   &err_children, &num_err_children));
  ar_assert(!err_children);
  ar_assert(!num_err_children);


  // Start the DMAs
  dmas_waiting = num_elements;
  for (i = 0; i < num_elements; i++) {
    noc_dma_add_new((sizes[i] >> MM_PACK_SIZE_BITS) &
                    ((1 << MM_PACK_LOCATION_BITS) - 1),
                    (void *) addresses[i],
                    -1, // dump to DRAM memory
                    (void *) addresses[i],
                    sizes[i] & ((1 << MM_PACK_SIZE_BITS) - 1),
                    1,  // I
                    0,  // W
                    0,  // C
                    &dmas_waiting,
                    NULL);
  }


  // Free size/address arrays used for packing
  kt_free(sizes);
  kt_free(addresses);


  // Destroy pool node
  ar_assert(!mm_slab_destroy_pool(node->pool, 1));

  // Reclaim all used ranges of this node and move them to the free pool
  for (i = 0; i < node->num_ranges; i++) {
    mm_region_add_adr_range(context->mm_free_ranges, NULL, -1,
                            node->ranges[i]->address,
                            node->ranges[i]->num_slabs, NULL);
    context->mm_free_num_slabs += node->ranges[i]->num_slabs;

    // The node->ranges contains exactly all separate ranges that belong to
    // this node (already maximally groupped by the node pointer). So, we
    // don't need to use reduce_adr_range here, we can delete them as they are.
    ar_assert(!kt_trie_delete(context->mm_used_ranges,
                              node->ranges[i]->address, kt_free));
  }

  // Free the ranges array
  kt_free(node->ranges);

  // Remove region ID from the local mappings
  ar_assert(!kt_trie_delete(context->mm_local_rids, node->id, NULL));

  // Reclaim the region ID and move it to the free pool
  mm_region_add_rid_range(context->mm_free_rids, context->pr_scheduler_id,
                          node->id, 1, NULL);
  context->mm_free_num_rids++;
  kt_trie_find_approx(context->mm_used_rids, 0, node->id, (void *) &rid_range);
  ar_assert(rid_range);
  ar_assert(node->id >= rid_range->rid);
  ar_assert(node->id < rid_range->rid + rid_range->num_rids);
  mm_region_reduce_rid_range(context->mm_used_rids, rid_range,
                             node->id, 1, NULL);

  // Dependency queues should be empty; free their structures.
  ar_assert(node->dep_queue);
  ar_assert(!kt_list_size(node->dep_queue));
  kt_free_list(node->dep_queue, NULL);

  ar_assert(node->obj_dep_queues);
  ar_assert(!kt_trie_size(node->obj_dep_queues));
  kt_free_trie(node->obj_dep_queues, NULL);

  // Free object location structure
  ar_assert(node->obj_locations);
  kt_free_trie(node->obj_locations, NULL);

  // Free the region tree node itself. We only destruct childless nodes,
  // so node->children_ids is supposed to be NULL (asserted above).
  kt_free(node);

  // Wait for any pending DMAs before we return
  while (dmas_waiting) {
    noc_dma_check_progress(NULL);
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
// mm_region_init()             Initializes the local region structures from
//                              an initial address set and creates the NULL
//                              pool region as the region tree root.
// ===========================================================================
// * INPUTS
//   int top_level              If 1, we are the top-level scheduler. Apart
//                              from initializing, claim all user memory and
//                              all possible region IDs.
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_region_init(int top_level) {

  Context       *context;
  size_t        free_adr;
  int           free_num_slabs;
  rid_t         free_rid;
  int           sched_id;


  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(!context->mm_region_tree);
  ar_assert(!context->mm_used_rids);
  ar_assert(!context->mm_free_rids);
  ar_assert(!context->mm_used_ranges);
  ar_assert(!context->mm_free_ranges);
  ar_assert(!context->mm_local_rids);


  // Create tries
  ar_assert(context->mm_used_rids = kt_alloc_trie(MM_TRIES_RID_MSB, 0));
  ar_assert(context->mm_free_rids = kt_alloc_trie(MM_TRIES_RID_MSB, 0));
  ar_assert(context->mm_used_ranges = kt_alloc_trie(MM_TRIES_POINTER_MSB,
                                                    MM_TRIES_POINTER_LSB));
  ar_assert(context->mm_free_ranges = kt_alloc_trie(MM_TRIES_POINTER_MSB,
                                                    MM_TRIES_POINTER_LSB));
  ar_assert(context->mm_local_rids = kt_alloc_trie(MM_TRIES_RID_MSB, 0));


  // Top-level scheduler?
  if (top_level) {

    // Initialize address ranges free pool
    mm_region_add_adr_range(context->mm_free_ranges, NULL, -1,
                            MM_VA_USER_BASE,
                            MM_USER_SIZE / MM_SLAB_SIZE,
                            NULL);
    ar_assert(!context->mm_free_num_slabs);
    context->mm_free_num_slabs = MM_USER_SIZE / MM_SLAB_SIZE;

    // Initialize region ID ranges free pool
    sched_id = context->pr_scheduler_id;
    mm_region_add_rid_range(context->mm_free_rids, sched_id,
                            1, MM_MAX_RID - 1, NULL);
    context->mm_free_num_rids = MM_MAX_RID - 1;

    // Create NULL region, and set it to be the region tree root. Set the
    // location that of the first worker core (where the main task runs
    // by definition).
    ar_assert(!mm_region_get_adr_range(0, &free_adr, &free_num_slabs));
    ar_assert(!mm_region_get_rid_range(1, &free_rid));
    ar_assert(!mm_region_create_node(0, free_rid, free_adr, free_num_slabs,
                                     pr_worker_core_id(0),
                                     &context->mm_region_tree));

    // NULL region has id 1. We depend on that later on.
    ar_assert(context->mm_region_tree->id == 1);
  }

  // Success
  return 0;
}


// ===========================================================================
// mm_region_create_region()    Creates a new region, parent of an existing
//                              one, and allocates memory for it.
// ===========================================================================
// * INPUTS
//   rid_t parent               Parent region id
//   int remote_mode            If 0, try to create a child region of a
//                              local parent, and fail if this can't happen.
//                              If 1, try to create a child region of a
//                              remote parent.
//   int remote_location        When remote_mode == 1, the initial location
//                              of the new node (i.e. the current core ID
//                              who has write permissions on the parent region)
//   int level_hint             Hint as to which scheduler level should
//                              maintain the new child. May be overridden.
//
// * OUTPUTS
//   rid_t *ret_rid             The newly created region's id
//
// * RETURN VALUE
//   int                        0: success, ret_rid is valid
//                              ERR_NO_SUCH_REGION: parent rid is not local,
//                                                  and remote_mode == 0
//                              ERR_OUT_OF_MEMORY: out of memory
//                              ERR_OUT_OF_RIDS: out of region IDs
//                              ERR_BOUNDARY_REGION: child region must be
//                                                   created by child
//                                                   scheduler
// ===========================================================================
int mm_region_create_region(rid_t parent, int remote_mode, int remote_location,
                            int level_hint, rid_t *ret_rid) {

  Context       *context;
  MmRgnTreeNode *parent_node;
  MmRgnTreeNode *new_node;
  size_t        free_adr;
  int           free_num_slabs;
  rid_t         free_rid;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_used_rids);
  ar_assert(context->mm_free_rids);
  ar_assert(context->mm_used_ranges);
  ar_assert(context->mm_free_ranges);
  ar_assert(context->mm_local_rids);
  ar_assert(ret_rid);
  ar_assert(parent);

  // Locate the parent, if needed
  if (!remote_mode) {

    kt_trie_find(context->mm_local_rids, parent, (void *) &parent_node);

    // If they told us to create a child of a local region, fail if we can't
    // find the parent.
    if (!parent_node) {
      return ERR_NO_SUCH_REGION;
    }
  }
  else {
    parent_node = NULL;
  }


  // If in local mode, which means we found the parent on this level, check if
  // the child must be ours. The only choices for the child is to belong at
  // the same level of its parent (i.e. our scheduler level) or be one level
  // lower (closer to workers). We depend on the hint for this decision.
  // Obviously, if we're at the lowest possible level, we can't go lower.
  if (!remote_mode &&
      (context->pr_scheduler_level > 0) &&
      (level_hint < context->pr_scheduler_level)) {
    return ERR_BOUNDARY_REGION;
  }

  // Get us a free region ID
  if (mm_region_get_rid_range(1, &free_rid)) {
    return ERR_OUT_OF_RIDS;
  }

  // Get us some free slabs
  if (mm_region_get_adr_range(0, &free_adr, &free_num_slabs)) {

    // Out of local memory. Give back the region ID we got before...
    mm_region_add_rid_range(context->mm_free_rids, context->pr_scheduler_id,
                            free_rid, 1, NULL);

    // ... and complain, so we can get more memory.
    return ERR_OUT_OF_MEMORY;
  }

  // Create the new region
  ar_assert(!mm_region_create_node(parent, free_rid, free_adr, free_num_slabs,
                                   remote_mode ? remote_location :
                                                 parent_node->location,
                                   &new_node));

  // Success
  *ret_rid = new_node->id;
//kt_printf("%d: Created region %d locally [location: %d]\r\n", context->pr_core_id, *ret_rid, new_node->location);
  return 0;
}


// ===========================================================================
// mm_region_update_parent()    Updates a local parent region that he either
//                              has a new, remote child or that he has to
//                              delete one of his existing, remote children.
// ===========================================================================
// * INPUTS
//   rid_t parent               Local parent region ID
//   rid_t child                Remote child region ID
//   int add_mode               1: Insert new child into parent
//                              0: Delete existing child from parent
//
// * RETURN VALUE
//   int                        0: success, ret_rid is valid
//                              ERR_NO_SUCH_REGION: parent rid is not local
// ===========================================================================
int mm_region_update_parent(rid_t parent, rid_t child, int add_mode) {

  Context       *context;
  MmRgnTreeNode *parent_node;
  int           found;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_local_rids);
  ar_assert(parent);
  ar_assert(child);

  // Locate the parent node
  kt_trie_find(context->mm_local_rids, parent, (void *) &parent_node);
  if (!parent_node) {
    return ERR_NO_SUCH_REGION;
  }

  // Update it
  if (add_mode) {
    ar_assert(parent_node->children_ids =
                  kt_realloc(parent_node->children_ids,
                             (parent_node->num_children + 1) * sizeof(rid_t)));
    parent_node->children_ids[parent_node->num_children++] = child;
  }
  else {
    // Locate region ID in parent's children
    found = 0;
    for (i = 0; i < parent_node->num_children; i++) {
      if (parent_node->children_ids[i] == child) {
        // Swap it with the last child
        found = 1;
        parent_node->children_ids[i] =
                      parent_node->children_ids[parent_node->num_children - 1];
        break;
      }
    }
    ar_assert(found);

    // Reduce array size
    parent_node->num_children--;
    if (parent_node->num_children) {
      ar_assert(parent_node->children_ids =
          kt_realloc(parent_node->children_ids,
                     parent_node->num_children * sizeof(rid_t)));
    }
    else {
      kt_free(parent_node->children_ids);
      parent_node->children_ids = NULL;
    }
  }

  // Success
  return 0;
}


// ===========================================================================
// mm_region_destroy_region()   Destroys a region along with any children
//                              it may have. All objects are massively
//                              freed (actually, all slabs and metadata
//                              objects are massively freed).
//
//                              Function fails if non-local regions exist,
//                              and gives back all the non-local region IDs
//                              that must be destroyed first. This is by
//                              choice, because otherwise we'd have zombie
//                              remote children (on a non-existent parent
//                              region ID, or (worse) on a parent ID which
//                              could have been given to a new node).
// ===========================================================================
// * INPUTS
//   rid_t region               The region ID of the region to be destroyed
//
// * OUTPUTS
//   rid_t **ret_err_regions    Array of remote region IDs identified as
//                              children and grandchildren of this region
//                              (valid when ERR_REMOTE_CHILDREN is returned)
//                              or single entry of the remote parent region
//                              (valid when ERR_REMOTE_PARENT is returned).
//                              Caller must free this when done with it.
//   int *ret_num_err_regions   Number of remote region IDs above
//   int *ret_num_frees         Number of successful region frees performed
//
// * RETURN VALUE
//   int                        0: success
//                              ERR_NO_SUCH_REGION: region id does not exist
//                              ERR_REMOTE_CHILDREN: remote children found.
//                                All their IDs are returned in
//                                ret_err_regions and local region (and its
//                                local children) are not yet freed, but
//                                all remote IDs are removed from the children
//                                arrays. Redo this function when the remote
//                                children have been freed, to complete the
//                                local destructions.
//                              ERR_REMOTE_PARENT: remote parent found. The
//                                local region is destroyed successfully, but
//                                the parent (whose ID is returned in
//                                ret_err_regions) must be updated.
// ===========================================================================
int mm_region_destroy_region(rid_t region, rid_t **ret_err_regions,
                             int *ret_num_err_regions, int *ret_num_frees) {

  Context       *context;
  rid_t         parent_id;
  MmRgnTreeNode **stack;
  int           stack_size;
  MmRgnTreeNode **queue;
  int           head;
  int           tail;
  MmRgnTreeNode *region_node;
  MmRgnTreeNode *node;
  MmRgnTreeNode *child_node;
  int           ret;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_used_rids);
  ar_assert(context->mm_free_rids);
  ar_assert(context->mm_used_ranges);
  ar_assert(context->mm_free_ranges);
  ar_assert(context->mm_local_rids);
  ar_assert(region);
  ar_assert(ret_err_regions);
  ar_assert(ret_num_err_regions);
  ar_assert(ret_num_frees);

  // Locate the region
  kt_trie_find(context->mm_local_rids, region, (void *) &region_node);
  if (!region_node) {
    return ERR_NO_SUCH_REGION;
  }

  // Remember its parent ID
  parent_id = region_node->parent_id;
  ar_assert(parent_id);

  // Nothing freed yet
  *ret_num_frees = 0;


  // =========================================================================
  // Phase 1: Use a queue to walk all region children, to make sure that no
  //          remote region IDs are there. If there are, gather them all.
  // =========================================================================

  // Create a queue with the region node
  ar_assert(queue = kt_malloc(sizeof(MmRgnTreeNode *)));
  queue[0] = region_node;
  head = 1;
  tail = 0;
  *ret_err_regions = NULL;
  *ret_num_err_regions = 0;

  // Iterate on the queue
  while (head > tail) {

    // Get next node
    node = queue[tail];
    tail++;

    // For all children IDs
    for (i = 0; i < node->num_children; i++) {
      kt_trie_find(context->mm_local_rids, node->children_ids[i],
                   (void *) &child_node);

      // Enqueue all local children...
      if (child_node) {
        ar_assert(queue = kt_realloc(queue,
                                     (head + 1) * sizeof(MmRgnTreeNode *)));
        queue[head] = child_node;
        head++;
      }

      // ... and return all the remote ones, removing them from the children
      // array.
      else {
        ar_assert(*ret_err_regions = kt_realloc(*ret_err_regions,
                                                 (*ret_num_err_regions + 1) *
                                                   sizeof(rid_t)));
        (*ret_err_regions)[(*ret_num_err_regions)] = node->children_ids[i];
        (*ret_num_err_regions)++;

        // Swap this child with the last one, and redo this child position
        node->children_ids[i] = node->children_ids[node->num_children - 1];
        node->num_children--;
        i--;
      }
    }

    // If we exhausted all the children, free the array
    if (!node->num_children) {
      kt_free(node->children_ids);
      node->children_ids = NULL;
    }
  }

  // Free queue
  kt_free(queue);

  // If we found any remotes, quit here
  if (*ret_num_err_regions) {
    return ERR_REMOTE_CHILDREN;
  }


  // =========================================================================
  // Phase 2: Now that we're sure no remote regions exist, use a stack to
  //          destroy the region subtree. Children are destroyed first,
  //          parents later.
  // =========================================================================

  // Create a stack with the region node
  ar_assert(stack = kt_malloc(sizeof(MmRgnTreeNode *)));
  stack[0] = region_node;
  stack_size = 1;

  // Iterate on the stack
  while (stack_size) {

    // Pop next node
    node = stack[stack_size - 1];
    stack_size--;
    stack = kt_realloc(stack, stack_size * sizeof(MmRgnTreeNode *));

    // Does it have children?
    if (node->num_children) {

      // Add node back to the stack, followed by all its children
      stack = kt_realloc(stack, (stack_size + 1 + node->num_children) *
                                sizeof(MmRgnTreeNode *));
      stack[stack_size] = node;
      for (i = 0; i < node->num_children; i++) {
        // All children should be local now
        kt_trie_find(context->mm_local_rids, node->children_ids[i],
                   (void *) &child_node);
        ar_assert(child_node);
        stack[stack_size + 1 + i] = child_node;
      }
      stack_size += 1 + node->num_children;

      // Node's children are now in the stack. Free node's array here, so that
      // next time we encounter the node we can delete it; by then it will be
      // really childless.
      kt_free(node->children_ids);
      node->children_ids = NULL;
      node->num_children = 0;
    }
    else {
      // Free this node
      ar_assert(!mm_region_free_node(node));
      (*ret_num_frees)++;
    }
  }


  // =========================================================================
  // Phase 3: Region is destroyed, but its parent is not updated. Try to do
  //          this here, and quit with an error if the parent is remote.
  // =========================================================================

  // Try updating the parent locally
  ret =  mm_region_update_parent(parent_id, region, 0);

  if (!ret) {
    // Success
    return 0;
  }
  else if (ret == ERR_NO_SUCH_REGION) {
    // The parent region is not handled by us. Store parent ID and fail.
    ar_assert(*ret_err_regions = kt_malloc(sizeof(rid_t)));
    **ret_err_regions = parent_id;
    *ret_num_err_regions = 1;
    return ERR_REMOTE_PARENT;
  }
  else {
    // Unexpected return code
    ar_abort();
  }

  // Success
  return 0;
}


// ===========================================================================
// mm_region_alloc_object()     Allocates a new object into an existing
//                              region. If the region memory is depleted,
//                              it tries to find some additional local memory
//                              to expand it and then retries the allocation.
// ===========================================================================
// * INPUTS
//   rid_t region               Region ID of the region to allocate object
//   size_t size                Allocation size
//
// * OUTPUTS
//   void **ret_ptr             The newly allocated pointer
//
// * RETURN VALUE
//   int                        0: success, ret_ptr is valid
//                              ERR_NO_SUCH_REGION: region id does not exist
//                              ERR_OUT_OF_MEMORY: out of memory
// ===========================================================================
int mm_region_alloc_object(rid_t region, size_t size, void **ret_ptr) {

  Context       *context;
  MmRgnTreeNode *node;
  size_t        free_adr;
  int           free_num_slabs;
  MmRgnAdrRange *range;
  int           ret;
  int           found;
  int           i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(region);
  ar_assert(size);
  ar_assert(context->mm_used_rids);
  ar_assert(context->mm_free_rids);
  ar_assert(context->mm_used_ranges);
  ar_assert(context->mm_free_ranges);
  ar_assert(context->mm_local_rids);
  ar_assert(ret_ptr);

  // Locate the region
  kt_trie_find(context->mm_local_rids, region, (void *) &node);
  if (!node) {
    return ERR_NO_SUCH_REGION;
  }

  // Try allocating the object
  if (mm_slab_alloc_slot(node->pool, size, (size_t *) ret_ptr)) {

    // Out of memory. Try to find some more to give to this region.
    if (mm_region_get_adr_range((size + MM_SLAB_SIZE - MM_ALLOC_ALIGN) /
                                MM_SLAB_SIZE, // size/slab_size, rounded up
                                &free_adr, &free_num_slabs)) {

      // No more local memory; fail miserably. We'll be called again when
      // our local memory pool has been extended.
      return ERR_OUT_OF_MEMORY;
    }

    // Extend region
    ar_assert(!mm_slab_expand_pool(node->pool, free_adr, free_num_slabs));

    // Update used ranges
    ret = mm_region_add_adr_range(context->mm_used_ranges, node,
                                  context->pr_scheduler_id,
                                  free_adr, free_num_slabs, &range);
    context->mm_free_num_slabs -= free_num_slabs;
    ar_assert(context->mm_free_num_slabs >= 0);

    // A new range may have been generated (if the new slabs were not
    // consecutive and could not be merged with one of the existing ranges we
    // have for this node). In this case, add the new range to the region tree
    // node array.
    if (ret == 2) {
      ar_assert(node->ranges = kt_realloc(node->ranges,
                                          (node->num_ranges + 1) *
                                          sizeof(MmRgnAdrRange *)));
      node->ranges[node->num_ranges++] = range;
    }
    // Also, the range may have filled a hole, merged with previous/next
    // neighbors and deleted the right neighbor. If so, we need to find the
    // deleted range and update the region tree node array.
    else if (ret == 1) {
      found = 0;
      for (i = 0; i < node->num_ranges; i++) {
        if (node->ranges[i] == range) {
          // Got it. Swap it with the last element.
          found = 1;
          node->ranges[i] = node->ranges[node->num_ranges - 1];
          break;
        }
      }
      ar_assert(found);
      ar_assert(node->ranges = kt_realloc(node->ranges,
                        --node->num_ranges * sizeof(MmRgnAdrRange *)));
    }

    // Try again; it should not fail, since we asked for minimum slabs
    // that can accomodate the allocation size.
    ar_assert(!mm_slab_alloc_slot(node->pool, size, (size_t *) ret_ptr));
  }

  // Success
  ar_assert(*ret_ptr);
  return 0;
}


// ===========================================================================
// mm_region_free_object()      Frees an object from an existing region.
//                              The pointer is mapped to the allegedly
//                              responsible region (based on the ranges we
//                              have distributed) and the slab free function
//                              takes care of the rest.
// ===========================================================================
// * INPUTS
//   void *ptr                  Pointer of the object to be freed
//
// * RETURN VALUE
//   int                        0: success
//                              ERR_OUT_OF_RANGE: pointer's range is not
//                                 handled locally
//                              ERR_NOT_ALLOCED: pointer is in our range, but
//                                 doesn't exist
//                              ERR_MISALIGNED: pointer not aligned to
//                                 MM_ALLOC_ALIGN or to slab slot
// ===========================================================================
int mm_region_free_object(void *ptr) {

  Context       *context;
  MmRgnAdrRange *range;
  int           size;
  int           core_id;
  int           dma_waiting;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_used_rids);
  ar_assert(context->mm_free_rids);
  ar_assert(context->mm_used_ranges);
  ar_assert(context->mm_free_ranges);
  ar_assert(context->mm_local_rids);
  ar_assert(ptr);

  // Locate the region responsible for this pointer
  kt_trie_find_approx(context->mm_used_ranges, 0, (size_t) ptr,
                      (void *) &range);
  if (!range) {
    return ERR_OUT_OF_RANGE;
  }
  ar_assert(range->address <= (size_t) ptr);
  if ((size_t) ptr >= range->address + range->num_slabs * MM_SLAB_SIZE) {
    return ERR_OUT_OF_RANGE;
  }

  // We found a region, but is it local?
  if (!range->node) {
    ar_assert(range->sched_id != context->pr_scheduler_id);
    return ERR_OUT_OF_RANGE;
  }
  ar_assert(range->node->pool);

  // Find out its size
  size = mm_slab_query_pointer(range->node->pool, (size_t) ptr);
  ar_assert(size != ERR_OUT_OF_RANGE); // used_ranges says it'll be there
  if (size < 0) {
    return size; // other error found
  }
  ar_assert(size);

  // Find out who is the last writer
  if (range->node->obj_locations) {
    if (!kt_trie_find(range->node->obj_locations, (size_t) ptr,
                      (void *) &core_id)) {
      core_id = range->node->location;
    }
  }
  else {
    core_id = range->node->location;
  }

  // Start a DMA to flush the object from the last writer's cache. Otherwise,
  // we risk the cache line to remain there, while the address is reallocated
  // to some other core. This can lead to stale data from the old writer
  // getting a writeback to the cache and corrupting the new writer location.
  dma_waiting = 1;
  noc_dma_add_new(core_id,
                  ptr,
                  -1, // dump to DRAM memory
                  ptr,
                  size,
                  1,  // I
                  0,  // W
                  0,  // C
                  &dma_waiting,
                  NULL);

  // Wait DMA to finish
  while (dma_waiting) {
    noc_dma_check_progress(NULL);
  }


  // Free it
  ar_assert(!mm_slab_free_slot(range->node->pool, (size_t) ptr));

  // Success
  return 0;
}


// ===========================================================================
// mm_region_query_region()     Tests if a region is local or not.
// ===========================================================================
// * INPUTS
//   rid_t region               Region to be tested
//
// * OUTPUTS
//   rid_t *ret_par_rid         If not NULL, the parent region ID of the
//                              given region
//
// * RETURN VALUE
//   int                        0: success, region is local
//                              ERR_NO_SUCH_REGION: region is not local
// ===========================================================================
int mm_region_query_region(rid_t region, rid_t *ret_par_rid) {

  Context               *context;
  MmRgnTreeNode         *r;


  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(region);

  // Locate the region
  if (kt_trie_find(context->mm_local_rids, region, (void *) &r)) {
    if (ret_par_rid) {
      *ret_par_rid = r->parent_id;
    }
    return 0;
  }
  else {
    return ERR_NO_SUCH_REGION;
  }
}


// ===========================================================================
// mm_region_query_pointer()    Locates a pointer inside the region and,
//                              if found, returns its size.
// ===========================================================================
// * INPUTS
//   void *ptr                  Pointer of the object to be queried
//
// * OUTPUTS
//   size_t *ret_size           If not NULL, returned size of ptr
//   rid_t *ret_rid             If not NULL, returned region ID of ptr
//
// * RETURN VALUE
//   int                        0: success
//                              ERR_OUT_OF_RANGE: ptr range not handled locally
//                              ERR_NOT_ALLOCED: ptr in range, but doesn't exist
//                              ERR_MISALIGNED: ptr not aligned properly
// ===========================================================================
int mm_region_query_pointer(void *ptr, size_t *ret_size, rid_t *ret_rid) {

  Context       *context;
  MmRgnAdrRange *range;
  int           ret;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_used_rids);
  ar_assert(context->mm_free_rids);
  ar_assert(context->mm_used_ranges);
  ar_assert(context->mm_free_ranges);
  ar_assert(context->mm_local_rids);
  ar_assert(ptr);
  if (ret_size) {
    *ret_size = 0;
  }
  if (ret_rid) {
    *ret_rid = 0;
  }

  // Locate the region responsible for this pointer
  kt_trie_find_approx(context->mm_used_ranges, 0, (size_t) ptr,
                      (void *) &range);
  if (!range) {
    return ERR_OUT_OF_RANGE;
  }
  ar_assert(range->address <= (size_t) ptr);
  if ((size_t) ptr >= range->address + range->num_slabs * MM_SLAB_SIZE) {
    return ERR_OUT_OF_RANGE;
  }
  if (!range->node) {
    return ERR_OUT_OF_RANGE;
  }
  ar_assert(range->node->pool);

  // Query it
  ret = mm_slab_query_pointer(range->node->pool, (size_t) ptr);
  ar_assert(ret != ERR_OUT_OF_RANGE); // used_ranges says it'll be here
  if (ret < 0) {
    return ret;
  }

  // Success
  if (ret_size) {
    *ret_size = ret;
  }
  if (ret_rid) {
    *ret_rid = range->node->id;
  }
  return 0;
}


// ===========================================================================
// mm_region_pack_region()      Creates a list of address ranges that
//                              represent the used slabs of the given region.
//                              If the region has children, they are also
//                              included in the list, if caller needs it.
// ===========================================================================
// * INPUTS
//   rid_t region               The region to be packed
//   int pack_options           Packing option flags (MM_PACK_OPTION_*)
//   int recurse                0: pack only this region, ignore its children
//                              1: if region has children, pack them as well
//
// * INOUTS
//   int *alloc_num_elements    Allocated size of the two arrays below
//                              (>= num_elements)
//   int *num_elements          Actual number of packed ranges in the
//                              two arrays below (<= alloc_num_elements)
//   int **addresses            Allocated array of addresses, one for each
//                              start of a packed range. Caller must
//                              initialize and free it.
//   int **sizes                Allocated array of sizes for each address
//                              range above. Caller must init/free it.
//
// * OUTPUTS
//   rid_t **ret_err_children   Array of remote region IDs identified as
//                              children and grandchildren of this region
//                              (valid when ERR_REMOTE_CHILDREN is returned).
//                              Caller must free this when done with it.
//   int *ret_num_err_children  Number of remote region IDs above
//
// * RETURN VALUE
//   int                        0: success
//                              ERR_NO_SUCH_REGION: region id not found
//                              ERR_REMOTE_CHILDREN: all local children are
//                                packed, but some remote children IDs were
//                                also found. In order to be complete, all
//                                the IDs in ret_err_children must be also
//                                packed.
// ===========================================================================
int mm_region_pack_region(rid_t region, int pack_options, int recurse,
                          int *alloc_num_elements, int *num_elements,
                          size_t **addresses, int **sizes,
                          rid_t **ret_err_children,
                          int *ret_num_err_children) {

  Context               *context;
  MmRgnTreeNode         **queue;
  int                   head;
  int                   tail;
  MmRgnTreeNode         *node;
  MmRgnTreeNode         *child_node;
  int                   i;


  // Sanity checks
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_used_rids);
  ar_assert(context->mm_free_rids);
  ar_assert(context->mm_used_ranges);
  ar_assert(context->mm_free_ranges);
  ar_assert(context->mm_local_rids);
  ar_assert(region);
  ar_assert(alloc_num_elements);
  ar_assert(num_elements);
  ar_assert(addresses);
  ar_assert(sizes);
  ar_assert(ret_err_children);
  ar_assert(ret_num_err_children);

  // Locate the region
  kt_trie_find(context->mm_local_rids, region, (void *) &node);
  if (!node) {
    return ERR_NO_SUCH_REGION;
  }

  // Create a queue with the region node
  ar_assert(queue = kt_malloc(sizeof(MmRgnTreeNode *)));
  queue[0] = node;
  head = 1;
  tail = 0;
  *ret_err_children = NULL;
  *ret_num_err_children = 0;

  // Iterate on the queue
  while (head > tail) {

    // Get next node
    node = queue[tail];
    tail++;

    // Pack it
    ar_assert(!mm_slab_query_pool(node->pool, pack_options, node->location,
                                  node->obj_locations, alloc_num_elements,
                                  num_elements, sizes, addresses));

    // Does it have children (and do we care)?
    if (recurse && node->num_children) {

      for (i = 0; i < node->num_children; i++) {

        // Is child local?
        kt_trie_find(context->mm_local_rids, node->children_ids[i],
                     (void *) &child_node);

        // Enqueue all local children...
        if (child_node) {
          ar_assert(queue = kt_realloc(queue,
                                       (head + 1) * sizeof(MmRgnTreeNode *)));
          queue[head] = child_node;
          head++;
        }

        // ... and return all the remote ones, removing them from the children
        // array.
        else {
          ar_assert(*ret_err_children = kt_realloc(*ret_err_children,
                                                   (*ret_num_err_children + 1) *
                                                     sizeof(rid_t)));
          (*ret_err_children)[(*ret_num_err_children)] = node->children_ids[i];
          (*ret_num_err_children)++;
        }
      }
    }
  }

  // Free queue
  kt_free(queue);

  // Found any remotes?
  if (*ret_num_err_children) {
    return ERR_REMOTE_CHILDREN;
  }
  else {
    // Success, no remotes
    return 0;
  }
}


// ===========================================================================
// mm_region_add_rid_range()    Adds a new region ID range into a trie of
//                              ranges. If contiguous (and scheduler
//                              compatible) with left- or right-neighbors, it
//                              merges the new IDs with them. If not, creates
//                              a new range into the trie.
// ===========================================================================
// * INPUTS
//   int sched_id               The scheduler ID, to be checked for
//                              neighbor compatibility (or to be set to the
//                              new range)
//   rid_t rid                  Starting region ID of the new range
//   int num_rids               Number of rids of the new range
//
// * OUTPUTS
//   MmRgnRidRange **ret_range  If not NULL, when a range is deleted or
//                              inserted (as per the return status) it is
//                              returned here.
//
// * INOUTS
//   Trie *t                    The trie to be updated
//
// * RETURN VALUE
//   MmRgnRidRange *            The neighbor that was merged/updated or the
//                              new range that was created
//   int                        0: range merged, no insertions/deletions
//                              1: range merged and filled a hole: ret_range
//                                 contains the dead pointer of the deleted
//                                 range (warning: its contents are invalid).
//                              2: range could not be merged; ret_range
//                                 contains the new range that was inserted
// ===========================================================================
int mm_region_add_rid_range(Trie *t, int sched_id, rid_t rid, int num_rids,
                            MmRgnRidRange **ret_range) {

  MmRgnRidRange *left;
  MmRgnRidRange *right;
  MmRgnRidRange *range;


  // Sanity checks
  ar_assert(t);
  ar_assert(rid);
  ar_assert(num_rids > 0);


  // Search ranges for a compatible left neighbor
  left = NULL;
  if (kt_trie_find_approx(t, 0, rid, (void *) &left)) {
    ar_assert(left->rid < rid);
    if (left->rid + left->num_rids != rid) {
      // Left neighbor is far away, not suitable for merging
      left = NULL;
    }
    else if (left->sched_id != sched_id) {
      // Left neighbor belongs to another scheduler, so it's incompatible to
      // merge
      left = NULL;
    }
  }

  // Search ranges for a right neighbor
  right = NULL;
  kt_trie_find(t, rid + num_rids, (void *) &right);
  if (right && (right->sched_id != sched_id)) {
    // Right neighbor belongs to another scheduler, so it's incompatible to
    // merge
    right = NULL;
  }

  // Four combinations of left/right
  if (left && !right) {
    // Expand left neighbor to include the new range
    left->num_rids += num_rids;
    if (ret_range) {
      *ret_range = NULL;
    }
    return 0;
  }
  else if (!left && right) {
    // Expand right neighbor to include the new range
    ar_assert(!kt_trie_delete(t, right->rid, NULL));
    right->rid -= num_rids;
    right->num_rids += num_rids;
    ar_assert(!kt_trie_insert(t, right->rid, right));
    if (ret_range) {
      *ret_range = NULL;
    }
    return 0;
  }
  else if (left && right) {
    // Expand left neighbor to include both the new range and the right
    // neighbor. Delete right neighbor, after returning it.
    left->num_rids += num_rids + right->num_rids;
    if (ret_range) {
      *ret_range = right;
    }
    ar_assert(!kt_trie_delete(t, right->rid, kt_free));
    return 1;
  }
  else {
    // No neighbors; create new range
    range = kt_zalloc(sizeof(MmRgnRidRange));
    range->rid = rid;
    range->num_rids = num_rids;
    range->sched_id = sched_id;
    ar_assert(!kt_trie_insert(t, rid, range));
    if (ret_range) {
      *ret_range = range;
    }
    return 2;
  }
}


// ===========================================================================
// mm_region_reduce_rid_range() Takes a specific part out of a given region ID
//                              range. If the range is exhausted, it is
//                              deleted. If not, it is modified accordingly.
//                              It may also happen that a hole is created in
//                              its middle; in this case, a second range is
//                              created to handle the remainder.
// ===========================================================================
// * INPUTS
//   rid_t rid                  Starting region ID of the taken range
//   int num_rids               Number of region IDs to be taken
//
// * OUTPUTS
//   MmRgnRidRange **ret_new    If not NULL and only when a new range is
//                              created (because range is split in two), the
//                              new range is returned.
//
// * INOUTS
//   Trie *t                    The trie to be updated
//   MmRgnRidRange *range       The range from which the IDs are taken
//
// * RETURN VALUE
//   int                        0: range modified, but still exists
//                              1: range exhausted and deleted
//                              2: range modified and split in two; ret_new
//                                 contains the new range
// ===========================================================================
int mm_region_reduce_rid_range(Trie *t, MmRgnRidRange *range, rid_t rid,
                               int num_rids, MmRgnRidRange **ret_new) {
  MmRgnRidRange *range2;


  // Sanity checks
  ar_assert(t);
  ar_assert(range);
  ar_assert(range->rid <= rid);
  ar_assert(range->rid + range->num_rids >= rid + num_rids);

  // Take IDs from the bottom part of range?
  if (range->rid == rid) {
    if (range->num_rids == num_rids) {
      // Range exhausted; delete it
      ar_assert(!kt_trie_delete(t, rid, kt_free));
      if (ret_new) {
        *ret_new = NULL;
      }
      return 1;
    }
    else {
      // Reinsert it with changed base address
      ar_assert(!kt_trie_delete(t, rid, NULL));
      range->rid += num_rids;
      range->num_rids -= num_rids;
      ar_assert(!kt_trie_insert(t, range->rid, range));
      if (ret_new) {
        *ret_new = NULL;
      }
      return 0;
    }
  }
  // Take IDs from the mid or high part of range?
  else if (range->rid < rid) {
    // High part?
    if (range->rid + range->num_rids == rid + num_rids) {
      // Just decrease the range
      range->num_rids -= num_rids;
      if (ret_new) {
        *ret_new = NULL;
      }
      return 0;
    }
    // Mid part?
    else {
      // Range must be split into two. Create a new one for the high
      // part...
      ar_assert(range2 = kt_malloc(sizeof(MmRgnRidRange)));
      range2->rid = rid + num_rids;
      range2->num_rids = (range->rid + range->num_rids) -
                         (rid + num_rids);
      range2->sched_id = range->sched_id;
      ar_assert(!kt_trie_insert(t, range2->rid, range2));

      // ... and decrease the range of the lower part.
      range->num_rids -= num_rids + range2->num_rids;
      ar_assert(range->num_rids > 0);
      ar_assert(range2->num_rids > 0);

      if (ret_new) {
        *ret_new = range2;
      }
      return 2;
    }
  }
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_region_get_rid_range()    Tries to find a number of free region IDs from
//                              the local free ID pool.
// ===========================================================================
// * INPUTS
//   int num_rids               Number of IDs requested
//
// * OUTPUTS
//   rid_t *ret_rid             Starting region ID of the free chunk found
//
// * RETURN VALUE
//   int                        0: success, ret_rid valid
//                              ERR_OUT_OF_RIDS: out of region IDs
// ===========================================================================
int mm_region_get_rid_range(int num_rids, size_t *ret_rid) {

  Context       *context;
  MmRgnRidRange *range;


  // Sanity checks
  ar_assert(ret_rid);
  ar_assert(num_rids > 0);

  // Get context
  context = mm_get_context(ar_get_core_id());

  // Search all free pool ranges
  for (kt_trie_find_minmax(context->mm_free_rids, 0, (void *) &range);
       range;
       kt_trie_find_next(context->mm_free_rids, 1, (void *) &range)) {

    // Does it cover the number of IDs we're looking for?
    if (range->num_rids >= num_rids) {

      // Take the IDs from the lowest address
      *ret_rid = range->rid;
      mm_region_reduce_rid_range(context->mm_free_rids, range, range->rid,
                                 num_rids, NULL);

      // Success
      return 0;
    }
  }

  // Nothing found?
  return ERR_OUT_OF_RIDS;
}


// ===========================================================================
// mm_region_add_adr_range()    Adds a new range into a trie of ranges. If
//                              contiguous (and node compatible) with left- or
//                              right-neighbors, it merges the new slabs with
//                              them. If not, creates a new range into the
//                              trie.
// ===========================================================================
// * INPUTS
//   MmRgnTreeNode *node        The region tree node, to be checked for
//                              neighbor compatibility (or to be set to the
//                              new range)
//   int sched_id               The scheduler ID, for the same purposes
//   size_t adr                 Starting address of the new range
//   int num_slabs              Number of slabs of the new range
//
// * OUTPUTS
//   MmRgnAdrRange **ret_range  If not NULL, contains the deleted or inserted
//                              ranges, if they happened (as indicated by the
//                              return value)
//
// * INOUTS
//   Trie *t                    The trie to be updated
//
// * RETURN VALUE
//   int                        0: range merged, no insertions/deletions
//                              1: range merged and filled a hole: ret_range
//                                 contains the dead pointer of the deleted
//                                 range (warning: its contents are invalid).
//                              2: range could not be merged; ret_range
//                                 contains the new range that was inserted
// ===========================================================================
int mm_region_add_adr_range(Trie *t, MmRgnTreeNode *node, int sched_id,
                            size_t adr, int num_slabs,
                            MmRgnAdrRange **ret_range) {

  MmRgnAdrRange *left;
  MmRgnAdrRange *right;
  MmRgnAdrRange *range;


  // Sanity checks
  ar_assert(t);
  ar_assert(adr);
  ar_assert(num_slabs > 0);


  // Search ranges for a compatible left neighbor
  left = NULL;
  if (kt_trie_find_approx(t, 0, adr, (void *) &left)) {
    ar_assert(left->address < adr);
    if (left->address + left->num_slabs * MM_SLAB_SIZE != adr) {
      // Left neighbor is far away, not suitable for merging
      left = NULL;
    }
    else if ((left->node != node) || (left->sched_id != sched_id)) {
      // Left neighbor belongs to another node and/or scheduler, so it's
      // incompatible to merge
      left = NULL;
    }
  }

  // Search ranges for a right neighbor
  right = NULL;
  kt_trie_find(t, adr + num_slabs * MM_SLAB_SIZE, (void *) &right);
  if (right && ((right->node != node) || (right->sched_id != sched_id))) {
    // Right neighbor belongs to another node and/or scheduler, so it's
    // incompatible to merge
    right = NULL;
  }

  // Four combinations of left/right
  if (left && !right) {
    // Expand left neighbor to include the new range
    left->num_slabs += num_slabs;
    if (ret_range) {
      *ret_range = NULL;
    }
    return 0;
  }
  else if (!left && right) {
    // Expand right neighbor to include the new range
    ar_assert(!kt_trie_delete(t, right->address, NULL));
    right->address -= num_slabs * MM_SLAB_SIZE;
    right->num_slabs += num_slabs;
    ar_assert(!kt_trie_insert(t, right->address, right));
    if (ret_range) {
      *ret_range = NULL;
    }
    return 0;
  }
  else if (left && right) {
    // Expand left neighbor to include both the new range and the right
    // neighbor. Delete the right one after returning its address.
    left->num_slabs += num_slabs + right->num_slabs;
    if (ret_range) {
      *ret_range = right;
    }
    ar_assert(!kt_trie_delete(t, right->address, kt_free));
    return 1;
  }
  else {
    // No neighbors; create new range
    range = kt_zalloc(sizeof(MmRgnAdrRange));
    range->address = adr;
    range->num_slabs = num_slabs;
    range->node = node;
    range->sched_id = sched_id;
    ar_assert(!kt_trie_insert(t, adr, range));
    if (ret_range) {
      *ret_range = range;
    }
    return 2;
  }
}


// ===========================================================================
// mm_region_reduce_adr_range() Takes a specific part out of a given range. If
//                              the range is exhausted, it is deleted. If not,
//                              it is modified accordingly. It may also happen
//                              that a hole is created in its middle; in this
//                              case, a second range is created to handle the
//                              remainder.
// ===========================================================================
// * INPUTS
//   size_t adr                 Starting address of the taken slabs
//   int num_slabs              Number of slabs to be taken
//
// * OUTPUTS
//   MmRgnAdrRange **ret_new    If not NULL and only when a new range is
//                              created (because range is split in two), the
//                              new range is returned.
//
// * INOUTS
//   Trie *t                    The trie to be updated
//   MmRgnAdrRange *range       The range from which the slabs are taken
//
// * RETURN VALUE
//   int                        0: range modified, but still exists
//                              1: range exhausted and deleted
//                              2: range modified and split in two; ret_new
//                                 contains the new range
// ===========================================================================
int mm_region_reduce_adr_range(Trie *t, MmRgnAdrRange *range, size_t adr,
                               int num_slabs, MmRgnAdrRange **ret_new) {
  MmRgnAdrRange *range2;


  // Sanity checks
  ar_assert(t);
  ar_assert(range);
  ar_assert(range->address <= adr);
  ar_assert(range->address + range->num_slabs * MM_SLAB_SIZE >=
            adr + num_slabs * MM_SLAB_SIZE);

  // Take slabs from the bottom part of range?
  if (range->address == adr) {
    if (range->num_slabs == num_slabs) {
      // Range exhausted; delete it
      ar_assert(!kt_trie_delete(t, adr, kt_free));
      if (ret_new) {
        *ret_new = NULL;
      }
      return 1;
    }
    else {
      // Reinsert it with changed base address
      ar_assert(!kt_trie_delete(t, adr, NULL));
      range->address += num_slabs * MM_SLAB_SIZE;
      range->num_slabs -= num_slabs;
      ar_assert(!kt_trie_insert(t, range->address, range));
      if (ret_new) {
        *ret_new = NULL;
      }
      return 0;
    }
  }
  // Take slabs from the mid or high part of range?
  else if (range->address < adr) {
    // High part?
    if (range->address + range->num_slabs * MM_SLAB_SIZE ==
        adr + num_slabs * MM_SLAB_SIZE) {
      // Just decrease the range
      range->num_slabs -= num_slabs;
      if (ret_new) {
        *ret_new = NULL;
      }
      return 0;
    }
    // Mid part?
    else {
      // Range must be split into two. Create a new one for the high
      // part...
      ar_assert(range2 = kt_malloc(sizeof(MmRgnAdrRange)));
      range2->address = adr + num_slabs * MM_SLAB_SIZE;
      range2->num_slabs = ((range->address + range->num_slabs * MM_SLAB_SIZE) -
                           (adr + num_slabs * MM_SLAB_SIZE)) /
                          MM_SLAB_SIZE;
      range2->node = range->node;
      range2->sched_id = range->sched_id;
      ar_assert(!kt_trie_insert(t, range2->address, range2));

      // ... and decrease the range of the lower part.
      range->num_slabs -= num_slabs + range2->num_slabs;
      ar_assert(range->num_slabs > 0);
      ar_assert(range2->num_slabs > 0);

      if (ret_new) {
        *ret_new = range2;
      }
      return 2;
    }
  }
  else {
    ar_abort();
  }
}


// ===========================================================================
// mm_region_get_adr_range()    Tries to find a number of free slabs from the
//                              local address pool. It first attempts to
//                              find enough slabs in the free pool. If this
//                              fails, it tries to take back slabs from
//                              existing regions. This happens progressively.
//                              First, all regions are asked to give back
//                              MM_ADR_RANGE_CHUNK_MAX slabs, if they can.
//                              This continues down to MM_ADR_RANGE_CHUNK_MIN.
//                              Whenever a region complies, the process is
//                              stopped and the point is stored. New calls to
//                              this function will continue from this last
//                              point.  When all regions are fully harvested,
//                              this process never repeats: only free pool
//                              requests are served, even if the local address
//                              pool receives more pages.
// ===========================================================================
// * INPUTS
//   int min_slabs              If >0, the minimum number of slabs that must
//                              be returned (otherwise, fail with "out of
//                              memory response). If 0, choose freely how
//                              many slabs to return based on free pool and/or
//                              existing regions free slabs capacities.
//
// * OUTPUTS
//   size_t *ret_adr            Starting address of the free chunk found
//   int *ret_num_slabs         Number of slabs in the free chunk
//
// * RETURN VALUE
//   int                        0: success, ret_adr and ret_num_slabs valid
//                              ERR_OUT_OF_MEMORY: out of memory
// ===========================================================================
int mm_region_get_adr_range(int min_slabs, size_t *ret_adr,
                            int *ret_num_slabs) {

  Context       *context;
  MmRgnAdrRange *range;
  size_t        used_adr;
  MmRgnTreeNode *node;
  size_t        cur_id;
  int           cur_chunk;
  int           steal_status;
  MmRgnAdrRange *steal_new_range;
  int           i;


  // Sanity checks
  ar_assert(ret_adr);
  ar_assert(ret_num_slabs);

  // If we got a minimum slabs request, align it to minimum chunk size
  if (min_slabs && (min_slabs & (MM_ADR_RANGE_CHUNK_MIN - 1))) {
    min_slabs = (min_slabs & ~(MM_ADR_RANGE_CHUNK_MIN - 1)) +
                MM_ADR_RANGE_CHUNK_MIN;
  }

  // Get context
  context = mm_get_context(ar_get_core_id());
  ar_assert(context->mm_range_chunk >= MM_ADR_RANGE_CHUNK_MIN);
  ar_assert(context->mm_range_chunk <= MM_ADR_RANGE_CHUNK_MAX);

  // Try to find a number of free slabs equal to the current free allocation
  // chunk; if it fails, try to take as many slabs are possible, down to the
  // minimum requested (min_slabs) or minimum allowed (MM_ADR_RANGE_CHUNK_MIN).
  // We're not lowering context->mm_range_chunk here; we simply try to
  // accomodate as many regions as possible from the free pool, before taking
  // back slabs from them -- even if this makes them unbalanced.
  for (*ret_num_slabs = ((min_slabs > context->mm_range_chunk) ?
                                                    min_slabs :
                                                    context->mm_range_chunk);
       *ret_num_slabs >= ((min_slabs) ? min_slabs : MM_ADR_RANGE_CHUNK_MIN);
       *ret_num_slabs /= 2) {

    // Search all free pool ranges
    for (kt_trie_find_minmax(context->mm_free_ranges, 0, (void *) &range);
         range;
         kt_trie_find_next(context->mm_free_ranges, 1, (void *) &range)) {

      // We never should have to deal in too fine-grained quantities
      ar_assert(range->num_slabs >= MM_ADR_RANGE_CHUNK_MIN);

      // Does it cover the number of slabs we're looking for?
      if (range->num_slabs >= *ret_num_slabs) {

        // Take the slabs from the lowest address
        *ret_adr = range->address;
        mm_region_reduce_adr_range(context->mm_free_ranges, range,
                                   range->address, *ret_num_slabs, NULL);

        // Success
        return 0;
      }
    }
  }


  // Find first region after the last one we harvested
  if (!context->mm_last_harvest) {
    // Never harvested? Grab the first one.
    cur_id = kt_trie_find_minmax(context->mm_local_rids, 0, (void *) &node);
  }
  else {
    // Search for the last one or its immediately next (last one may have
    // been freed by now)
    cur_id = kt_trie_find_approx(context->mm_local_rids, 1,
                                context->mm_last_harvest, (void *) &node);
    if (cur_id == context->mm_last_harvest) {
      // If last one still exists, get the next
      cur_id = kt_trie_find_next(context->mm_local_rids, 1, (void *) &node);
    }
  }

  // If we are already at the minimum allocation chunk point and have searched
  // in the past among all regions, immediately ask for more pages. This is a
  // way to get better performance, sacrificing best memory fit (e.g. regions
  // could have freed lots of objects from the time we last harvested slabs).
  // In the common case, though, it should be ok. Note that if a region was
  // fully freed, we'd find its slabs in the free pool.
  if (!cur_id && (context->mm_range_chunk == MM_ADR_RANGE_CHUNK_MIN)) {
    return ERR_OUT_OF_MEMORY;
  }


  // Search all region pools for free slabs. If none of them can give what
  // we ask for, halve the quantity we ask for until we drop to the minimum
  // chunk we're willing to work with.
  cur_chunk = context->mm_range_chunk;
  while (cur_chunk >= MM_ADR_RANGE_CHUNK_MIN) {

    // Remember chunk size for future harvests
    context->mm_range_chunk = cur_chunk;

    // If a minimum number is requested, don't get past this point
    if (min_slabs && (min_slabs > cur_chunk)) {
      return ERR_OUT_OF_MEMORY;
    }

    // If no region id, restart from the beginning
    if (!cur_id) {
      cur_id = kt_trie_find_minmax(context->mm_local_rids, 0, (void *) &node);
    }
    while (cur_id) {
      // Remember point for future harvests
      context->mm_last_harvest = cur_id;

      // Never empty regions completely, even if they are momentarily unused.
      // User may want to create a whole region tree first, and then do
      // allocations.
      if (node->pool->free_space <= cur_chunk * MM_SLAB_SIZE) {
        goto cont;
      }

      // Can we get slabs from this region's pool?
      if (!mm_slab_reduce_pool(node->pool, cur_chunk, ret_adr)) {

        // Store return number of slabs. Return address is already stored.
        *ret_num_slabs = cur_chunk;

        // Look up the address in the used ranges
        ar_assert(used_adr = kt_trie_find_approx(context->mm_used_ranges, 0,
                                                 *ret_adr, (void *) &range));

        // Update used ranges to reflect we took away these slabs.
        ar_assert(range);
        ar_assert(range->node == node);
        ar_assert(range->sched_id == context->pr_scheduler_id);
        ar_assert(used_adr == range->address);
        ar_assert(*ret_num_slabs <= range->num_slabs);
        steal_status = mm_region_reduce_adr_range(context->mm_used_ranges,
                                                  range, *ret_adr,
                                                  *ret_num_slabs,
                                                  &steal_new_range);
        context->mm_free_num_slabs += *ret_num_slabs;

        // Update region tree node array of used ranges
        if (steal_status == 1) {

          // Range was exhausted and deleted. Find it and swap it with the
          // last element of the array.
          ar_assert(node->num_ranges > 1); // pools never remain with 0 mem
          for (i = 0; i < node->num_ranges; i++) {
            if (node->ranges[i] == range) {
              if (i != node->num_ranges - 1) {
                // It's not the last entry, so swap with it
                node->ranges[i] = node->ranges[node->num_ranges - 1];
              }
              ar_assert(node->ranges = kt_realloc(node->ranges,
                                --node->num_ranges * sizeof(MmRgnAdrRange *)));
              break;
            }
          }

        }
        else if (steal_status == 2) {

          // Range was modified and split in two. Add the new range.
          ar_assert(node->ranges = kt_realloc(node->ranges,
                            (node->num_ranges + 1) * sizeof(MmRgnAdrRange *)));
          node->ranges[node->num_ranges++] = steal_new_range;
        }

        // Success
        return 0;
      }
cont:
      // Go to the next one
      cur_id = kt_trie_find_next(context->mm_local_rids, 1, (void *) &node);
    }

    // Halve chunk. Note that when we reach the minimum, we never increase
    // it again. So, don't set MM_ADR_RANGE_CHUNK_MIN ridiculously low.
    cur_chunk /= 2;
  }

  // Nothing could be found. Request more memory.
  return ERR_OUT_OF_MEMORY;
}



// ###########################################################################
// ###                                                                     ###
// ###                         DEBUGGING FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################

#if 0
// ===========================================================================
// mm_region_dbg_print_all()    Debug function to print all local info
// ===========================================================================
void mm_region_dbg_print_all() {

  size_t key;
  MmRgnTreeNode *node;
  MmRgnAdrRange *range;
  Context *context = mm_get_context();


  printf("===============================================================\n");

  printf("Region IDs:\n");
  for (key = kt_trie_find_minmax(context->mm_local_rids, 0, (void **) &node);
       key;
       key = kt_trie_find_next(context->mm_local_rids, 1, (void **) &node)) {
    ar_assert(node->pool);
    printf("  id %8lu -> metadata pool %p, used %lu, free %lu bytes\n",
        key,
        node->pool->metadata_pool,
        node->pool->used_space,
        node->pool->free_space);
  }
  printf("\n");

  printf("Used ranges:\n");
  for (key = kt_trie_find_minmax(context->mm_used_ranges, 0, (void **) &range);
       key;
       key = kt_trie_find_next(context->mm_used_ranges, 1, (void **) &range)) {
    ar_assert(range->address == key);
    ar_assert(range->node);
    printf("  [0x%08lx - 0x%08lx] -> region id %8lu\n",
        key,
        range->address + range->num_slabs * MM_SLAB_SIZE,
        range->node->id);
  }
  printf("\n");

  printf("Free ranges:\n");
  for (key = kt_trie_find_minmax(context->mm_free_ranges, 0, (void **) &range);
       key;
       key = kt_trie_find_next(context->mm_free_ranges, 1, (void **) &range)) {
    ar_assert(range->address == key);
    ar_assert(!range->node);
    printf("  [0x%08lx - 0x%08lx]\n",
        key, range->address + range->num_slabs * MM_SLAB_SIZE);
  }

  printf("===============================================================\n");
}
#endif
