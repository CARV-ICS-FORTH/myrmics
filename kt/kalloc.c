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
// Abstract      : Kernel memory allocation functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: kalloc.c,v $
// CVS revision  : $Revision: 1.21 $
// Last modified : $Date: 2013/03/06 16:06:29 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <memory_management.h>
#include <kernel_toolset.h>
#include <arch.h>


// Memory map during bootstrap:
//
// KERNEL_HEAP_START                    MM_ALLOC_ALIGN objects 
// KERNEL_HEAP_START + 
//    MM_BOOTSTRAP_SLABS_STEP *
//    MM_SLAB_SIZE                      2 * MM_ALLOC_ALIGN objects
// KERNEL_HEAP_START + 
//    2 * MM_BOOTSTRAP_SLABS_STEP *
//    MM_SLAB_SIZE                      3 * MM_ALLOC_ALIGN objects
// ...
// KERNEL_HEAP_END -    Counter of how many pointers to be freed follow
//   MM_PAGE_SIZE       
// KERNEL_HEAP_END -    
//   MM_PAGE_SIZE +     
//   sizeof(int)        First pointer to be freed
// KERNEL_HEAP_END -    
//   MM_PAGE_SIZE +     
//   sizeof(int) +
//   sizeof(void *)     Second pointer to be freed
// KERNEL_HEAP_END -    
//   MM_PAGE_SIZE +     
//   sizeof(int) +
//   2 * sizeof(void *) Third pointer to be freed
// ...
// KERNEL_HEAP_END      Counter for alloc() calls of MM_ALLOC_ALIGN objects
// KERNEL_HEAP_END + 4  Counter for alloc() calls of MM_ALLOC_ALIGN * 2 objects
// KERNEL_HEAP_END + 8  Counter for alloc() calls of MM_ALLOC_ALIGN * 3 objects
// ...


// ===========================================================================
// kt_mem_init()                Initializes the kernel memory allocator,
//                              bootstraps it, creates the kernel slab pool
//                              and allocates/initalizes the global Context.
// ===========================================================================
void kt_mem_init() {

  Context       *context;
  int           my_cid;
  size_t        kernel_base;
  size_t        kernel_end;
  int           i;


  // Get core ID and kernel limits
  my_cid = ar_get_core_id();
  kernel_base = mm_va_kernel_base(my_cid);
  kernel_end = kernel_base + MM_KERNEL_SIZE - 1024 * 1024;
  

  // Prepare bootstrapping structures on 2 last kernel heap pages. We are
  // going to keep track of how many objects we allocate per slab slot size,
  // to a maximum of MM_BOOTSTRAP_MAX_SLOT, as well as which of these
  // objects are freed during bootstrap.
  ar_assert(MM_BOOTSTRAP_MAX_SLOT % MM_ALLOC_ALIGN == 0);
  for (i = 0; i < (MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN); i++) {
    // A counter for alloc() calls
    ((int *) kernel_end)[i] = 0;
  }
  // Single counter for free() calls
  *((int *) (kernel_end - MM_PAGE_SIZE)) = 0;


  // Make sure we have enough space from the maximum number of objects during
  // bootstrap, so that they do not reach the two last kernel space pages
  // which are used for the counters above
  ar_assert((MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN) *
              MM_BOOTSTRAP_SLABS_STEP * MM_SLAB_SIZE <  // max bootstrapped adr
             kernel_end - MM_PAGE_SIZE);                // free() counter


  // Prepare context
  context = mm_get_context(my_cid);

  context->dbg_trc_data = NULL;
  context->dbg_trc_idx = 0;
  context->dbg_trc_time_start = 0;
  context->dbg_trc_offset = 0;
  
  context->dbg_stats_data = NULL;
  context->dbg_stats_last_tmr = 0;

  context->noc_mode = -1;
  context->noc_cnt_free = NULL;
  context->noc_send_buf = NULL;
  context->noc_recv_buf = NULL;
  context->noc_msg_core_id = -1;
  context->noc_credits = NULL;
  context->noc_credits_poll = NULL;
  context->noc_num_peers = 0;
  context->noc_poll_rr = 0;
  context->noc_active_dmas = NULL;
#ifdef NOC_WARN_OUT_OF_CREDITS
  context->noc_cred_warned = 0;
#endif

  context->mm_alloc_bootstrap = 1;
  context->mm_frees_bootstrap = 1;
  context->mm_kernel_pool = NULL;
  context->mm_recursion_depth = 0;
  for (i = 0; i < MM_SLAB_SIZE / MM_ALLOC_ALIGN; i++) {
    context->mm_prealloc_flags[i] = 0;
  }
  context->mm_prealloc_needed = 0;
  context->mm_busy_freeing = 0;
  context->mm_defer_frees = NULL;
  context->mm_num_defer_frees = 0;

  context->mm_region_tree = NULL;
  context->mm_used_rids = NULL;
  context->mm_free_rids = NULL;
  context->mm_used_ranges = NULL;
  context->mm_free_ranges = NULL;
  context->mm_free_num_slabs = 0;
  context->mm_free_num_rids = 0;
  context->mm_local_rids = NULL;
  context->mm_range_chunk = MM_ADR_RANGE_CHUNK_MAX;
  context->mm_last_harvest = 0;
  context->mm_load_rrobin = 0;
  context->mm_current_load = 0;
  context->mm_reported_load = 0;
  context->mm_children_load = NULL;

  context->pr_core_id = -1;
  context->pr_num_cores = -1;
  context->pr_role = -1;
  context->pr_core_bid_cid = NULL;
  context->pr_core_work_ids = NULL;
  context->pr_core_sched_ids = NULL;
  context->pr_core_child_idx = NULL;
  context->pr_work_core_ids = NULL;
  context->pr_sched_core_ids = NULL;
  context->pr_core_route = NULL;
  context->pr_num_schedulers = -1;
  context->pr_num_workers = -1;
  context->pr_worker_id = -1;
  context->pr_scheduler_id = -1;
  context->pr_parent_sched_id = -1;
  context->pr_scheduler_level = -1;
  context->pr_children = NULL;
  context->pr_num_children = -1;
  context->pr_cur_epoch = 0;
  context->pr_task_table = NULL;
  context->pr_ready_queue = NULL;
  context->pr_spawn_pending = 0;
  context->pr_tasks = NULL;
  context->pr_avail_task_id = 1;
  context->pr_sched_rr = 0;
  context->pr_load_vs_locality = PR_LOAD_VS_LOCALITY_DEFAULT;
  context->pr_cur_sched_load = 0;
  context->pr_cur_run_load = 0;
  context->pr_rep_sched_load = 0;
  context->pr_rep_run_load = 0;
  context->pr_chld_sched_load = NULL;
  context->pr_chld_run_load = NULL;
  context->pr_main_finished = 0;

  context->pr_message_id = 1;
  context->pr_pending_events = NULL;
  context->pr_incomplete_req = NULL;
  context->pr_pages_msg_id = 0;
  context->pr_rids_msg_id = 0;

  context->vid_demo_in_bid = -1;
  context->vid_demo_in_cid = -1;
  context->vid_demo_out_bid = -1;
  context->vid_demo_out_cid = -1;

#if 0
  context->sys_sched_sec = 0;
  context->sys_sched_usec = 0;
  context->sys_sched_calls = 0;
  context->sys_worker_sec = 0;
  context->sys_worker_usec = 0;
  context->sys_worker_sent = 0;
  context->sys_worker_recv = 0;
  context->sys_barrier_sec = 0;
  context->sys_barrier_usec = 0;
#endif

#ifdef ARCH_ARM
  context->fs_flash_num_sectors = 0;
  context->fs_num_blocks = 0;
  context->fs_max_inodes = 0;
  context->fs_state = NULL;
  context->fs_fdesc = NULL;
#endif

#ifdef FMPI
  context->fmpi = NULL;
#endif


  // Context is at a predefined location. Make sure the bootstrapping structs
  // know about it, so the space can be accounted for later. This is the 
  // first malloc we ever call, so no other object can take this predefined
  // location (which mm_get_context() has already returned above).
  ar_assert(kt_malloc(sizeof(Context)) == context);

  // Bootstrap memory system and create kernel memory pool
  ar_assert((MM_ALLOC_ALIGN & (MM_ALLOC_ALIGN - 1)) == 0);
  ar_assert(MM_PAGE_SIZE % MM_SLAB_SIZE == 0);
  ar_assert(MM_SLAB_SIZE % MM_ALLOC_ALIGN == 0);
  ar_assert(kernel_base % MM_SLAB_SIZE == 0);
  ar_assert(kernel_end % MM_SLAB_SIZE == 0);

  // Create kernel pool. The function will handle the rest of the
  // bootstrapping process, end the bootstrap mode and fill
  // context->mm_kernel_pool with the kernel pool.
  mm_slab_create_pool(NULL, kernel_base, MM_KERNEL_SIZE / MM_SLAB_SIZE, 1);
  ar_assert(!context->mm_alloc_bootstrap);
  ar_assert(!context->mm_frees_bootstrap);
}


// ===========================================================================
// kt_malloc()                  Kernel basic allocation function. Allocates
//                              serially in predefined places (see above) 
//                              during bootstrap, or calls the slab allocator
//                              out of bootstrap.
// ===========================================================================
// * INPUTS
//   size_t size                Number of bytes to be allocated (can be 0)
//
// * RETURN VALUE
//   void *                     Pointer to new allocated object. Note that
//                              NULL will not be returned for a non-zero size
//                              -- out of memory in kernel space will trigger
//                              an abort.
// ===========================================================================
void *kt_malloc(size_t size) {
  
  Context       *context;
  int           my_cid;
  size_t        kernel_base;
  size_t        kernel_end;
  int           *bootstrap_slots;
  size_t        ptr;
  int           i;


  // Allow dummy mallocs
  if (!size) {
    return NULL;
  }
  
  // Clamp requests up to 2-GB size. We won't support more, even for the
  // x86_64 port. It makes it easier to work internally with signed integers,
  // because error checking and assertions work way better.
  ar_assert (size < (1 << 31));

  // Get global context and boundaries
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);
  kernel_base = mm_va_kernel_base(my_cid);
  kernel_end = kernel_base + MM_KERNEL_SIZE - 1024 * 1024;

  // Align size request to nearest allowed size
  if (size & (MM_ALLOC_ALIGN - 1)) {
    size = (size & ~(MM_ALLOC_ALIGN - 1)) + MM_ALLOC_ALIGN;
  }

  // Bootstrapping code
  if (context->mm_alloc_bootstrap) {
    bootstrap_slots = (int *) kernel_end;

    // We support only up to MM_BOOTSTRAP_MAX_SLOT sized requests; increase
    // that if this assertion fails (it means bigger objects are needed and
    // bootstrap has to allow this)
    ar_assert(size <= MM_BOOTSTRAP_MAX_SLOT);

    // Allocate directly from the beginning of the kernel heap, keeping track
    // on the last heap page
    i = size / MM_ALLOC_ALIGN - 1; // find slot
    ptr = kernel_base +                                // base address
          i * MM_SLAB_SIZE * MM_BOOTSTRAP_SLABS_STEP + // slabs per slot
          bootstrap_slots[i] * size;                   // slot address
    
    // Remember how many objects we've allocated
    bootstrap_slots[i]++;

    // Make sure enough slots can fit into MM_BOOTSTRAP_SLABS_STEP; otherwise,
    // the define must be increased
    ar_assert(bootstrap_slots[i] * size <= 
           MM_BOOTSTRAP_SLABS_STEP * MM_SLAB_SIZE);

    return (void *) ptr;
  }


  // Do a normal allocation from the kernel pool
  if (mm_slab_alloc_slot(context->mm_kernel_pool, size, &ptr)) {
    // Kernel memory should never get full
    ar_abort();
  }
  ar_assert(ptr >= kernel_base);
  ar_assert(ptr <  kernel_end + MM_PAGE_SIZE);

  return (void *) ptr;
}


// ===========================================================================
// kt_zalloc()                  Wrapper of kt_malloc which also zeroes out
//                              the new object.
// ===========================================================================
// * INPUTS
//   size_t size                Number of bytes to be allocated
//
// * RETURN VALUE
//   void *                     Pointer to new allocated object. Note that
//                              NULL must not be returned -- out of memory
//                              in kernel space should trigger an abort.
// ===========================================================================
void *kt_zalloc(size_t size) {
  void *ptr;
  ptr = kt_malloc(size);
  kt_bzero(ptr, size);
  return ptr;
}


// ===========================================================================
// kt_free()                    Frees an object. During bootstrap, it simply
//                              tracks it for freeing later on. Out of
//                              bootstrap, calls the slab allocator to free
//                              it.
// ===========================================================================
// * INPUTS
//   void *ptr                  Object to be freed
// ===========================================================================
void kt_free(void *ptr) {

  Context       *context;
  int           my_cid;
  size_t        kernel_base;
  size_t        kernel_end;
  int           *bootstrap_slots;
  int           *counter;
  void          **free_slots;
  int           slot_id;
  int           slot_offset;


  // Dummy free?
  if (!ptr) {
    return;
  }
  
  // Get global context and boundaries
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);
  kernel_base = mm_va_kernel_base(my_cid);
  kernel_end = kernel_base + MM_KERNEL_SIZE - 1024 * 1024;

  // Sanity checks
  ar_assert (!((size_t) ptr & (MM_ALLOC_ALIGN - 1)));
  ar_assert((size_t) ptr >= kernel_base);
  ar_assert((size_t) ptr <  kernel_end + MM_PAGE_SIZE);

  // Get global context
  context = mm_get_context(ar_get_core_id());

  // Bootstrapping?
  if (context->mm_frees_bootstrap) {

    // Verify it's about a slot we actually gave
    bootstrap_slots = (int *) kernel_end;
    slot_id = ((size_t) ptr - kernel_base) / 
        (MM_BOOTSTRAP_SLABS_STEP * MM_SLAB_SIZE);
    ar_assert(slot_id * MM_ALLOC_ALIGN <= MM_BOOTSTRAP_MAX_SLOT);
    ar_uint_divide((size_t) ptr - kernel_base -
                       (slot_id * MM_BOOTSTRAP_SLABS_STEP * MM_SLAB_SIZE),
                   (slot_id + 1) * MM_ALLOC_ALIGN,
                   (unsigned int *) &slot_offset, NULL);
    ar_assert(slot_offset < bootstrap_slots[slot_id]);
    
    // Record the free request
    counter = (int *) ((size_t) kernel_end - MM_PAGE_SIZE);
    free_slots = (void **) ((size_t) kernel_end - MM_PAGE_SIZE + sizeof(int));
    free_slots[*counter] = ptr;
    (*counter)++;
    
    // Make sure we don't overflow the array
    ar_assert((MM_PAGE_SIZE - sizeof(int)) / sizeof(void *) > *counter);

    return;
  }

  // Do a normal free from the kernel pool
  ar_assert(!mm_slab_free_slot(context->mm_kernel_pool, (size_t) ptr));
}


// ===========================================================================
// kt_realloc()                 Reallocates an object to a new size
// ===========================================================================
// * INPUTS
//   void *old_ptr              The old object to be reallocated
//   size_t new_size            Number of bytes for the new object
//
// * RETURN VALUE
//   void *                     Pointer to new allocated object. Note that
//                              NULL must not be returned -- out of memory
//                              in kernel space should trigger an abort.
// ===========================================================================
void *kt_realloc(void *old_ptr, size_t new_size) {
  
  int           old_size;
  void          *new_ptr;
  Context       *context;


  // Query old pointer. We accept NULL pointers in realloc.
  context = mm_get_context(ar_get_core_id());
  if (old_ptr) {
    old_size = mm_slab_query_pointer(context->mm_kernel_pool, (size_t) old_ptr);
    ar_assert(old_size > 0);
  }
  else {
    old_size = 0;
  }

  // Clamp requests up to 2-GB size (see kt_malloc above).
  ar_assert (new_size < (1 << 31));

  // Align size request to nearest allowed size
  if (new_size & (MM_ALLOC_ALIGN - 1)) {
    new_size = (new_size & ~(MM_ALLOC_ALIGN - 1)) + MM_ALLOC_ALIGN;
  }

  // It may happen that the same size is requested, e.g. by asking for 
  // 4, then 8, then 12, ..., but alignment forces all of them to be 64 bytes.
  // If so, do nothing.
  if (old_size == new_size) {
    return old_ptr;
  }

  // Underlying slab allocator will never* return the same pointer for a
  // different slot size. Allocate a new pointer before freeing the old,
  // so that we avoid Tries or other allocations taking its position before
  // we manage to copy its contents.
  //
  // *: correct is "almost never". Definitely true for < slab_size requests 
  //    where not enough empty slabs are preallocated. For other cases, it can
  //    happen that if we freed the pointer before the new malloc, there is a
  //    chance that will happen. For kernel objects, which are small and
  //    dominated by mallocs, not frees, it is nearly impossible as a case. We
  //    choose not to follow it, because it can seriously mess stuff up if an
  //    intermediate malloc (e.g. a Trie kt_malloc) takes up this space and
  //    overwrites the data before we finish here -- and there's no easy way
  //    of knowing that.
  if (new_size) {
    ar_assert(new_ptr = kt_malloc(new_size));
  }
  else {
    new_ptr = NULL;
  }

  // Copy contents
  if (new_ptr && old_size) {
    kt_memcpy(new_ptr, old_ptr, (old_size < new_size) ? old_size : new_size);
  }

  // Free old pointer
  if (old_ptr) {
    kt_free(old_ptr);
  }

  // Success
  return new_ptr;
}
