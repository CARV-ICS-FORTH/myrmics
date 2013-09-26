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
// Abstract      : Slab allocator functions. Management of abstract slab 
//                 pools, which handle allocation/freeing of slots onto
//                 largish slabs. Same slot sizes are allocated from the same
//                 slabs, so that the pool can be kept relatively packed even
//                 after large requests of allocs/frees.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: slab.c,v $
// CVS revision  : $Revision: 1.5 $
// Last modified : $Date: 2013/02/21 07:57:48 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <memory_management.h>
#include <kernel_toolset.h>
#include <arch.h>



// ###########################################################################
// ###                                                                     ###
// ###                          INTERNAL FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################

#define mm_slab_is_full(slab) \
               ((!slab->free_slots) && \
                (!slab->offset || !slab->trailer_bit))

#define mm_slab_is_empty(slab) \
               ((slab->free_slots == slab->bitmap_size) && \
                (!slab->offset || (slab->trailer_bit == 1)))


// ===========================================================================
// mm_slab_init_bitmap()        Initializes a slab bitmap
// ===========================================================================
// * INPUTS
//   int bitmap size            Actual size of the bitmap
//
// * INOUTS
//   unsigned int *bitmap       The bitmap
// ===========================================================================
static inline void mm_slab_init_bitmap(unsigned int *bitmap, 
                                       int bitmap_size) {
  int i;

  for (i = 0; i < bitmap_size; i += 32) {
    bitmap[i/32] = 0xFFFFFFFF;
  }
}


// ===========================================================================
// mm_slab_find_free_bit()      Searches for the first free bit (value 1) in
//                              a slab bitmap. Search is done from LSB to MSB.
// ===========================================================================
// * INPUTS
//   int bitmap size            Actual size of the bitmap
//   unsigned int *bitmap       The bitmap
//
// * RETURN VALUE
//   int                        -1: no free bit found
//                              otherwise: first free bit index (LSB = 0)
// ===========================================================================
static inline int mm_slab_find_free_bit(unsigned int *bitmap, 
                                        int bitmap_size) {
  int i;
  int bit;

  for (i = 0; i < bitmap_size; i += 32) {
    bit = ar_ffs(bitmap[i / 32]);

    if (bit) {
      if (i + bit > bitmap_size) {
        return -1;
      }
      else {
        return (i + bit - 1);
      }
    }
  }
  return -1;
}


// ===========================================================================
// mm_slab_set_used_bit()       Marks a bit as "used" in a slab bitmap
// ===========================================================================
// * INPUTS
//   int bitmap size            Actual size of the bitmap
//
// * INOUTS
//   unsigned int *bitmap       The bitmap
// ===========================================================================
static inline void mm_slab_set_used_bit(unsigned int *bitmap, 
                                          int bitmap_size, int bit) {
  ar_assert(bit >= 0);
  ar_assert(bit < bitmap_size);
  bitmap[bit >> 5] &= ~(1 << (bit - (bit & 0xFFFFFFE0)));
}


// ===========================================================================
// mm_slab_set_free_bit()       Marks a bit as "free" in a slab bitmap
// ===========================================================================
// * INPUTS
//   int bitmap size            Actual size of the bitmap
//
// * INOUTS
//   unsigned int *bitmap       The bitmap
// ===========================================================================
static inline void mm_slab_set_free_bit(unsigned int *bitmap, 
                                        int bitmap_size, int bit) {
  ar_assert(bit >= 0);
  ar_assert(bit < bitmap_size);
  bitmap[bit >> 5] |= (1 << (bit - (bit & 0xFFFFFFE0)));
}


// ===========================================================================
// mm_slab_get_bit()            Reads a bit in a slab bitmap
// ===========================================================================
// * INPUTS
//   unsigned int *bitmap       The bitmap
//   int bitmap size            Actual size of the bitmap
//   int bit                    The bit position
//
// * RETURN VALUE
//   int                        0:   bit is used
//                              > 0: bit is free
// ===========================================================================
static inline int mm_slab_get_bit(unsigned int *bitmap, int bitmap_size, 
                                  int bit) {
  ar_assert(bit >= 0);
  ar_assert(bit < bitmap_size);
  return (bitmap[bit >> 5] & (1 << (bit - (bit & 0xFFFFFFE0))));
}


// ===========================================================================
// mm_slab_swap_list_nodes()    Swaps two nodes in a doubly-linked list.
//                              This will work even if n1 and n2 are
//                              neighbors. Taken from ptspts.blogspot.com/2010/
//                              01/how-to-swap-two-nodes-in-doubly-linked.html
//                              (swapping in doubly linked lists is
//                              surprisingly difficult to do it elegantly...)
// ===========================================================================
// * INOUTS
//   MmSlab *n1                 The first node
//   MmSlab *n2                 The second node
// ===========================================================================
static inline void mm_slab_swap_list_nodes(MmSlab *n1, MmSlab *n2) {
  MmSlab *tmp;

  tmp = n1->next;
  n1->next = n2->next;
  n2->next = tmp;
  if (n1->next != NULL) {
    n1->next->prev = n1;
  }
  if (n2->next != NULL) {
    n2->next->prev = n2;
  }
  tmp = n1->prev;
  n1->prev = n2->prev;
  n2->prev = tmp;
  if (n1->prev != NULL) {
    n1->prev->next = n1;
  }
  if (n2->prev == NULL) {
    return;
  }
  n2->prev->next = n2;
}


// ===========================================================================
// mm_slab_locate_pointer()     Locates a given pointer in the slab pool and
//                              returns the slab and bitmap position that 
//                              it is represented by.
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool
//   size_t ptr                 The pointer to be searched
//
// * OUTPUTS
//   MmSlab *ret_slab           If not NULL, the slab the pointer belongs to
//   int *ret_bit               If not NULL, the bit position in the bitmap
//                              of the slab the pointer belongs to
//
// * RETURN VALUE
//   int                        0: success, ret_slab and bit are valid
//                              ERR_MISALIGNED: pointer not aligned to 
//                                 MM_ALLOC_ALIGN or to slab slot
//                              ERR_OUT_OF_RANGE: pointer not inside this pool
//                                 used/partials slabs
//                              ERR_NOT_ALLOCED: pointer in range, but mapped 
//                                 to empty slot or free slab
// ===========================================================================
int mm_slab_locate_pointer(MmSlabPool *pool, size_t ptr, MmSlab **ret_slab,
                           int *ret_bit) {
  
  MmSlab        *slab;
  size_t        key;
  unsigned int  bit;
  unsigned int  rem;
 

  // Check the pointer for alignment
  if ((ptr & (MM_ALLOC_ALIGN - 1))) {
    return ERR_MISALIGNED;
  }

  // Find slab. The Trie will ignore least significant bits, so don't bother
  // to isolate the slab-aligned address.
  key = kt_trie_find(pool->used_slabs, ptr, (void *) &slab);
  if (!key) {
    return ERR_OUT_OF_RANGE;
  }
  ar_assert(slab);
  ar_assert(key == slab->address);

  // Slab must not be empty. We may find empty slabs in the used_slabs trie,
  // because they can be part of partial lists.
  if (slab->state == SLAB_EMPTY) {
    return ERR_NOT_ALLOCED;
  }

  // Map pointer LSBs to a bit position
  ptr &= (1UL << MM_TRIES_POINTER_LSB) - 1;
  if (ptr < slab->offset) {
    return ERR_MISALIGNED;
  }
  ptr -= slab->offset;
  ar_uint_divide(ptr, slab->slot_size, &bit, &rem);
  if (rem) {
    return ERR_MISALIGNED;
  }
  ar_assert(bit < slab->bitmap_size);

  // Is the bit used?
  if (!mm_slab_get_bit(slab->bitmap, slab->bitmap_size, bit)) {
    if (ret_slab) {
      *ret_slab = slab;
    }
    if (ret_bit) {
      *ret_bit = bit;
    }
    return 0;
  }
  else {
    return ERR_NOT_ALLOCED;
  }
}


// ===========================================================================
// mm_slab_promote_to_full()    Internal function that promotes an empty or
//                              partial slab to full state, updating the
//                              partial linked list 
// ===========================================================================
// * INOUTS
//   MmSlabPool *pool           The slab pool
//   MmSlabPartialHead *head    The partial linked list head for this slot
//                              (cannot be NULL: slab should be part of it)
//   MmSlab *slab               The slab to be promoted
// ===========================================================================
void mm_slab_promote_to_full(MmSlabPool *pool, MmSlabPartialHead *head, 
                             MmSlab *slab) {

  // Sanity checks
  ar_assert(head);
  ar_assert(head->head);

  // Head bookkeeping
  head->num_members--;
  if (slab->state == SLAB_EMPTY) {
    head->num_empty_members--;
  }
  head->used_slots -= slab->bitmap_size;

  // State is now full
  ar_assert(!slab->free_slots);
  slab->state = SLAB_FULL;
  
  // Remove slab from the partial list
  if (slab->next) {
    slab->next->prev = slab->prev;
  }
  if (slab->prev) {
    slab->prev->next = slab->next;
  }
  
  // Maybe we have to select a new head, if we were the head
  if (head->head == slab) {
    ar_assert(!slab->prev);
    if (slab->next) {
      head->head = slab->next;
    }
    else {
      // No more partial slabs left in the linked list; delete the head
      ar_assert(!head->num_members);
      ar_assert(!head->num_empty_members);
      ar_assert(!head->used_slots);
      ar_assert(!head->free_slots);
      ar_assert(!kt_trie_delete(pool->partial_slabs, slab->slot_size, kt_free));
    }
  }

  // These pointers are useless for a full slab
  slab->prev = slab->next = NULL;
}


// ===========================================================================
// mm_slab_demote_from_full()   Internal function that demotes a full slab to
//                              a partial or empty state, updating the
//                              partial linked list 
// ===========================================================================
// * INOUTS
//   MmSlabPool *pool           The slab pool
//   MmSlabPartialHead *head    The partial linked list head for this slot
//                              (can be NULL if no partial/empty slabs exist)
//   MmSlab *slab               The slab to be demoted
//
// * RETURN VALUE
//   MmSlabPartialHead *        The partial linked list head for this slot,
//                              which is either the argument the user passed,
//                              or a new head if NULL was given.
// ===========================================================================
MmSlabPartialHead *mm_slab_demote_from_full(MmSlabPool *pool, 
                                            MmSlabPartialHead *head, 
                                            MmSlab *slab) {

  // Verify head, or create a new one
  if (head) {
    ar_assert(head->head);
  }

  // Update slab state
  if (mm_slab_is_empty(slab)) {
    slab->state = SLAB_EMPTY;
  }
  else {
    slab->state = SLAB_PARTIAL;
  }

  // Do not continue with partial list arrangements if slab got empty, but
  // the slot size is big. It won't be part of the partial list, it's going
  // to be transfered to the free pool.
  if ((slab->state == SLAB_EMPTY) && (slab->slot_size >= MM_SLAB_SIZE)) {
    ar_assert(!slab->prev);
    ar_assert(!slab->next);
    return head;
  }

  // We'll definitely need a head if there is none. This can happen if all
  // slabs of the slot size are full and we are demoting the first one of them
  // to a partial/empty status.
  if (!head) {
    ar_assert(head = kt_zalloc(sizeof(MmSlabPartialHead)));
    ar_assert(!kt_trie_insert(pool->partial_slabs, slab->slot_size, head));
  }
 
  // Rest of the code is for head arrangements
  if (!head) {
    ar_assert(!slab->prev);
    ar_assert(!slab->next);
    return NULL;
  }

  // Head bookkeeping
  head->num_members++;
  if (slab->state == SLAB_EMPTY) {
    head->num_empty_members++;
  }
  head->free_slots += slab->free_slots;
  head->used_slots += slab->bitmap_size - slab->free_slots;

  // Add slab to the head of the partial list
  slab->prev = NULL;
  if (head->head) {
    ar_assert(!head->head->prev);
    head->head->prev = slab;
    slab->next = head->head;
  }
  else {
    slab->next = NULL;
  }
  head->head = slab;

  return head;
}


// ===========================================================================
// mm_slab_alloc_new()          Internal function that allocates a new slab,
//                              takes care of a few fields, connects it to
//                              lists and tries, etc.
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool for the new slab
//   int claim_slot             When >0: mark a slot as "used"
//                              (bitmap[claim_slot - 1] is claimed)
//   int claim_trailer          When 1: mark the trailer_bit as "used"
//   size_t address             Slab address
//   int slot_size              Slot size
//   int offset                 Offset of first object start
//   int num_slabs              Number of slabs claimed by last slot
//                              (see .h file for field details)
//   int insert_pos             0: replace head in partial linked list
//                              1: insert after head in partial linked list
//
// * RETURN VALUE
//   MmSlab *                   The newly allocated slab
// ===========================================================================
MmSlab *mm_slab_alloc_new(MmSlabPool *pool, int claim_slot, int claim_trailer,
                          size_t address, int slot_size, int offset, 
                          int num_slabs, int insert_pos) {
  
  MmSlab                *new_slab;
  MmSlabPartialHead     *head;
  size_t                head_adr;


  // Allocate a new slab
  ar_assert(new_slab = kt_zalloc(sizeof(MmSlab)));

  // Create the bitmap, all slots set free
  ar_assert(slot_size >= MM_ALLOC_ALIGN);
  // that's (MM_SLAB_SIZE - offset) / slot_size, rounded up
  ar_uint_divide(MM_SLAB_SIZE - offset + slot_size - MM_ALLOC_ALIGN,
                 slot_size, (unsigned int *) &(new_slab->bitmap_size), NULL);
  mm_slab_init_bitmap(new_slab->bitmap, new_slab->bitmap_size);
  new_slab->free_slots = new_slab->bitmap_size;
  new_slab->trailer_bit = (offset) ? 1 : 42; // 42 is don't care value

  // Claim anything?
  if (claim_trailer) {
    ar_assert(offset);
    new_slab->trailer_bit = 0;
  }
  if (claim_slot) {
    mm_slab_set_used_bit(new_slab->bitmap, new_slab->bitmap_size, 
                         claim_slot - 1);
    new_slab->free_slots--;
  }

  // Fill out mundane entries
  new_slab->slot_size = slot_size;
  new_slab->offset = offset;
  new_slab->num_slabs = num_slabs;
  new_slab->address = address;

  // Mark its state
  if (mm_slab_is_full(new_slab)) {
    new_slab->state = SLAB_FULL;
  }
  else if (mm_slab_is_empty(new_slab)) {
    new_slab->state = SLAB_EMPTY;
  }
  else {
    new_slab->state = SLAB_PARTIAL;
  }


  // If it's not a full slab, hook it up on the partial slabs list
  if (new_slab->state != SLAB_FULL) {
    head_adr = kt_trie_find(pool->partial_slabs, slot_size, (void *) &head);
    // Is there a list at all for this slot size?
    if (head_adr) {
      ar_assert(head->head);
      ar_assert(!head->head->prev);

      // Update head bookkeeping 
      head->num_members++;
      if (new_slab->state == SLAB_EMPTY) {
        head->num_empty_members++;
      }
      head->free_slots += new_slab->free_slots;
      head->used_slots += new_slab->bitmap_size - new_slab->free_slots;

      // Insert new slab before the head
      if (insert_pos == 0) {
        new_slab->next = head->head;
        head->head->prev = new_slab;
        head->head = new_slab;
      }
      // Insert after the head
      else if (insert_pos == 1) {
        if (head->head->next) {
          head->head->next->prev = new_slab;
        }
        new_slab->next = head->head->next;
        head->head->next = new_slab;
        new_slab->prev = head->head;
      }
      else {
        ar_abort();
      }
      
    }
    // No list: create one with us as a head
    else {
      ar_assert(head = kt_malloc(sizeof(MmSlabPartialHead)));
      head->head = new_slab;
      head->num_members = 1;
      head->num_empty_members = (new_slab->state == SLAB_EMPTY) ? 1 : 0;
      head->free_slots = new_slab->free_slots;
      head->used_slots = new_slab->bitmap_size - new_slab->free_slots;
      ar_assert(!kt_trie_insert(pool->partial_slabs, slot_size, head));
    }
  }

  // Empty, partial or full, also insert it to the used slabs trie. All slabs
  // must be addressable. We also insert empty slabs here, so we can easily
  // locate them by address in alloc() calls.
  ar_assert(!kt_trie_insert(pool->used_slabs, new_slab->address, new_slab));

  return new_slab;
}


// ===========================================================================
// mm_slab_destroy_empty()      Internal function that deallocates an empty
//                              slab and removes it from the partial list
//                              and pool tries
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool for the new slab
//   MmSlab *slab               The slab to be destroyed
//   MmSlabPartialHead *head    The partial list head for the slab slot size.
//                              Can be NULL if slot size is big enough to
//                              not support preallocated empty slabs and no
//                              partial slabs exist for this slot size.
//   int is_in_list             1 if slab is in the list
// ===========================================================================
void mm_slab_destroy_empty(MmSlabPool *pool, MmSlab *slab, 
                           MmSlabPartialHead *head, int is_in_list) {
  

  // Sanity checks
  ar_assert(slab);
  ar_assert(slab->state == SLAB_EMPTY);
  ar_assert(slab->free_slots == slab->bitmap_size);
  if (slab->offset) {
    ar_assert(slab->trailer_bit);
  }
  if (head) {
    ar_assert(head->head);
    ar_assert(!head->head->prev);
  }
  else {
    ar_assert(slab->slot_size >= MM_SLAB_SIZE);
  }

  // Do we deal with the list?
  if (is_in_list) {

    // Head bookkeeping
    ar_assert(head);
    head->num_members--;
    head->num_empty_members--;
    head->free_slots -= slab->free_slots;
    ar_assert(head->num_members >= 0);
    ar_assert(head->num_empty_members >= 0);
    ar_assert(head->num_empty_members <= head->num_members);

    // Remove slab from the partial list
    if (slab->next) {
      slab->next->prev = slab->prev;
    }
    if (slab->prev) {
      slab->prev->next = slab->next;
    }
    
    // Maybe we have to select a new head, if we were the head
    if (head->head == slab) {
      ar_assert(!slab->prev);
      if (slab->next) {
        head->head = slab->next;
      }
      else {
        // No more partial slabs left in the linked list; delete the head
        ar_assert(!head->num_members);
        ar_assert(!head->num_empty_members);
        ar_assert(!head->used_slots);
        ar_assert(!head->free_slots);
        ar_assert(!kt_trie_delete(pool->partial_slabs, slab->slot_size, 
                                  kt_free));
      }
    }
  }
  
  // Make sure nothing very stupid happened
  else {
    ar_assert(!slab->prev);
    ar_assert(!slab->next);
    if (head) {
      ar_assert(head->head != slab);
    }
  }

  // Delete it from the used slabs trie and free the slab
  ar_assert(!kt_trie_delete(pool->used_slabs, slab->address, kt_free));
}


// ===========================================================================
// mm_slab_claim_trailer()      Searches for a partial or empty slab at a
//                              specific address, which has a free trailer
//                              that can hold a given quantity. If it finds 
//                              it, it claims the trailer. 
//
//                              A variation happens if slab is actually given
//                              by the caller. In this case, it is simply 
//                              claimed (and not searched for).
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool
//   MmSlabPartialHead *head    The head of the partial linked list for this
//                              slot size
//   MmSlab *slab               If NULL, the function searches for a partial
//                              or empty slab that matches the specs. If not
//                              NULL, this slab is used directly.
//   size_t adr                 The address of the slab we're looking for
//                              if slab is NULL. Otherwise, don't care.
//   int slot_size              The slot size it has to have to be compatible
//   int claim_size             Size of the trailer it has to have to be 
//                              compatible
//
// * RETURN VALUE
//   int                        0 for success, ERR_GENERIC for failure
// ===========================================================================
int mm_slab_claim_trailer(MmSlabPool *pool, MmSlabPartialHead *head,
                          MmSlab *slab, size_t adr, int slot_size, 
                          int claim_size) {

  // Sanity checks
  ar_assert(head);
  ar_assert(head->head);
  ar_assert(claim_size < slot_size);

  // Do we need to search?
  if (!slab) {
    // Search among all used slabs, looking for a partial or empty
    // of compatible size
    if (!kt_trie_find(pool->used_slabs, adr, (void *) &slab)) {
      // No such address
      return ERR_GENERIC;
    }
  }
    
  // Verify its compatibility. It might not belong to this slot size,
  // since we searched all used slabs of the pool. We don't steal from
  // others, even if it's empty (it would violate the reasons we keep
  // empty slabs in partial lists).
  if (slab->slot_size != slot_size) {
    return ERR_GENERIC;
  }

  // Is it partial?
  if (slab->state == SLAB_PARTIAL) {
    // Offset must be compatible
    if (slab->offset != claim_size) {
      return ERR_GENERIC;
    }

    // If it is, the trailer_bit must be free so we can use the
    // first slot. 
    if (!slab->trailer_bit) {
      return ERR_GENERIC;
    }
  }
  // Is it empty?
  else if (slab->state == SLAB_EMPTY) {
    // An empty slab is definitely a hit. Offset is of no importance, we can
    // change that. Note that the offset may be "wrong", because slab might
    // have been brought by a free slab. Also, bitmap size and free slots may
    // be accordingly wrong. Fix these here.
    slab->offset = claim_size;
    // that's (MM_SLAB_SIZE - claim_size) / slot_size, rounded up
    ar_uint_divide(MM_SLAB_SIZE - claim_size + slot_size - MM_ALLOC_ALIGN,
                   slot_size, (unsigned int *) &(slab->bitmap_size), NULL);
    if (slab->free_slots != slab->bitmap_size) {
      head->free_slots += slab->bitmap_size - slab->free_slots;
      slab->free_slots = slab->bitmap_size;
    }
  }
  else {
    // Full slabs are of no interest
    return ERR_GENERIC;
  }

  // Ok, it passed all the tests: slab is valid to be used for the
  // trailer of slab. Claim the trailer.
  slab->trailer_bit = 0;

  // If no other slots are available, it just become full. Change status and
  // remove from the partial linked list.
  if (!slab->free_slots) {
    mm_slab_promote_to_full(pool, head, slab);
  }
  // Still has free slots? Leave it in the partial list as it is.
  // Make sure it's marked PARTIAL, because it may have been EMPTY.
  else if (slab->state == SLAB_EMPTY) {
    slab->state = SLAB_PARTIAL;
    head->num_empty_members--;
  }

  // Success
  return 0;
}


// ===========================================================================
// mm_slab_take_from_free()     Internal function that takes some slabs out
//                              of a selected free slab. It creates up to two
//                              new slabs to store the new slot, depending
//                              on the slot size and if the slot has started
//                              or ended already (externally to this function)
//                              on a partial slab. The free slab is modified
//                              to depict the reduced size, or deleted 
//                              altogether if it has no more slabs.
//
//                              Preallocation also happens here: for small
//                              slot sizes, when the free slab is reduced a 
//                              number of slabs are possibly moved to the
//                              partial linked list, ready to be allocated
//                              later by setting bits on their bitmaps.
// ===========================================================================
// * INPUTS
//   MmSlab *free_slab          The free slab to take slabs from
//   int num_slabs              How many slabs to take
//   int slot_size              The slot size we're working on
//   int size_on_partial        How many bytes we've already externally 
//                              allocated on a partial slab (0 if no partial
//                              allocation has been done)
//   int partial_is_trailer     0: size_on_partial refers to start of object
//                              1: size_on_partial refers to object trailer
//   int take_from_bottom       0: take slabs from the high part of free_slab
//                              1: take slabs from the low part of free_slab
//
// * INOUTS
//   MmSlabPool *pool           The slab pool where it all takes place
//
// * RETURN VALUE
//   size_t                     Address of the start of the new slot. Only 
//                              meaningful if size_on_partial == 0 or 
//                              partial_is_trailer == 1.
// ===========================================================================
size_t mm_slab_take_from_free(MmSlabPool *pool, MmSlab *free_slab, 
                              int num_slabs, int slot_size, 
                              int size_on_partial, int partial_is_trailer,
                              int take_from_bottom) {
  MmSlab        *new_slab;
  int           rem_size;
  int           extra_slabs;
  size_t        adr;
  int           i;
  int           slot_start;
  int           slot_id;


  // Sanity checks
  ar_assert(pool);
  ar_assert(free_slab);
  ar_assert(num_slabs >= 1);
  ar_assert(!(slot_size & (MM_ALLOC_ALIGN - 1)));
  if (partial_is_trailer) {
    // If a trailer is allocated, highest part of free_slab is needed to
    // prepend it
    ar_assert(!take_from_bottom);
  }

  // Size that must be claimed is the slot size, minus whatever space
  // may have already been externally allocated on a partial slab
  rem_size = slot_size - size_on_partial;
  ar_assert(!(rem_size & (MM_ALLOC_ALIGN - 1)));

  // If the chunk we take from the free slab does not begin on a partial
  // slab, or if it ends on a partial slab, then a new slab must be allocated
  // for the start of the object.
  if (!size_on_partial) {
    // Subcase 1: remaining size starts from offset 0, slot 0
    new_slab = mm_slab_alloc_new(
                  pool, 
                  1,
                  0,
                  (take_from_bottom) ? 
                          free_slab->address :
                          free_slab->address + (free_slab->num_slabs - 
                                               num_slabs) * MM_SLAB_SIZE,
                  slot_size,
                  0,
                  num_slabs,
                  0);
    // Start of slot is start of slab
    adr = new_slab->address;
  }
  else if (partial_is_trailer) {
    // Subcase 2: take_bottom=0. Remaining size starts from a variable
    //            offset and slot, depending on how much it is.
    ar_assert(rem_size > 0);
    slot_start = (MM_SLAB_SIZE - (rem_size % MM_SLAB_SIZE)) % MM_SLAB_SIZE;
    ar_uint_divide(slot_start, slot_size, (unsigned int *) &slot_id, NULL);
    new_slab = mm_slab_alloc_new(
                  pool, 
                  slot_id + 1,
                  0,
                  free_slab->address + (free_slab->num_slabs - 
                                       num_slabs) * MM_SLAB_SIZE,
                  slot_size,
                  slot_start - slot_id * slot_size,
                  num_slabs + 1, // + 1: slab extends to external trailer!
                  0);
    // Start of slot is start of slab plus slot offset
    adr = new_slab->address + slot_start;
  }
  else {
    // Adr is don't care; we won't be starting an object during this function.
    adr = 0;
  }
  
  // If the taken size from the slabs is not fully aligned to the end of the
  // last slab, we definitely need a new slab to store the trailer of the
  // slot. 
  //
  // Things are a little bit different if we already allocated the object
  // start on a new slab above: in this case, if the slot is small enough to
  // fit in a single slab, we don't need a second slab.
  if (((!size_on_partial) && 
       (slot_size > MM_SLAB_SIZE) && (rem_size % MM_SLAB_SIZE)) ||
      ((size_on_partial) && 
       (rem_size % MM_SLAB_SIZE) && (!partial_is_trailer))) {
    new_slab = mm_slab_alloc_new(
                  pool, 
                  0,
                  1,
                  (take_from_bottom) ? 
                          free_slab->address + (num_slabs - 1) * MM_SLAB_SIZE :
                          free_slab->address + (free_slab->num_slabs - 1) * 
                                               MM_SLAB_SIZE,
                  slot_size,
                  rem_size % MM_SLAB_SIZE,
                  1,
                  0);
  }

  // We already took num_slabs from the free pool. If the slot size is small
  // (< slab size), we prefer to take in total MM_BOOTSTRAP_SLABS_STEP slabs
  // from the free pool and hook the extras onto the partial list. We do that
  // for two reasons:
  // - Reduce frequency of changes to the free slab Trie, which may result in
  //   freeing a tree branch and allocating a new one (if we take addresses
  //   from the bottom side)
  // - Bootstrapping code depends on it. E.g. to allocate 512-B slots mingled
  //   with 64-B slots, bootstrap assumes that 512-B slots will be at a given
  //   offset. For each allowed bootstrap slot size, MM_BOOTSTRAP_SLABS_STEP
  //   slabs are assumed.
  if ((num_slabs < MM_BOOTSTRAP_SLABS_STEP) && (slot_size < MM_SLAB_SIZE)) {
    
    // How many are we missing?
    ar_assert(free_slab->num_slabs >= num_slabs);
    if (MM_BOOTSTRAP_SLABS_STEP <= free_slab->num_slabs) {
      extra_slabs = MM_BOOTSTRAP_SLABS_STEP - num_slabs;
    }
    else {
      extra_slabs = free_slab->num_slabs - num_slabs;
    }

    // Create a new slab for each of them. We append them after the current
    // head in the partial list, in order to fill the partial one first.
    // We create them inversely addressed (highest adr first) so that they
    // are inserted in ascending order in the list.
    for (i = extra_slabs - 1; i >= 0; i--) {
      new_slab = mm_slab_alloc_new(
                  pool, 
                  0,
                  0,
                  (take_from_bottom) ? 
                          free_slab->address + (num_slabs + i) * MM_SLAB_SIZE :
                          free_slab->address + (free_slab->num_slabs -
                                                num_slabs - extra_slabs + i) * 
                                               MM_SLAB_SIZE,
                  slot_size,
                  0,
                  1,
                  1);
    }
  }
  else {
    extra_slabs = 0;
  }

  // Update free slab with the numner of slabs it's lost
  free_slab->num_slabs -= num_slabs + extra_slabs;

  // It may happen the free slab is exhausted...
  if (!free_slab->num_slabs) {
    // ... in that case, delete it from the free slabs pool
    ar_assert(!kt_trie_delete(pool->free_slabs, free_slab->address, kt_free));
  }
  else {
    if (take_from_bottom) {

      // Slab start address itself is also changing, so we also need to update
      // the trie with the new key
      ar_assert(!kt_trie_delete(pool->free_slabs, free_slab->address, NULL));
      free_slab->address += (num_slabs + extra_slabs) * MM_SLAB_SIZE;
      ar_assert(!kt_trie_insert(pool->free_slabs, free_slab->address, 
                free_slab));
    }
    // else: nothing to be done, free_slab->num_slabs is updated, it still has
    //       some slabs left and it's correctly addressed in the free slabs
    //       pool (its base address has not changed).
  }

  return adr;
}


// ===========================================================================
// mm_slab_add_to_free()        Internal function that adds a range of slabs
//                              to the free pool. It tries to merge the slabs
//                              with the left and right neighbors, if
//                              possible. If impossible, creates a new free
//                              slab into the free pool that represents this
//                              range.
// ===========================================================================
// * INPUTS
//   size_t adr                 Starting address of slabs to be added
//   int num_slabs              Number of slabs to be added
//
// * INOUTS
//   MmSlabPool *pool           The slab pool
// ===========================================================================
void mm_slab_add_to_free(MmSlabPool *pool, size_t adr, int num_slabs) {

  MmSlab        *left;
  MmSlab        *right;
  MmSlab        *slab;


  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->free_slabs);
  ar_assert(adr);
  ar_assert(num_slabs > 0);


  // Search free slabs for a neighbor on our left side (address below ours).
  // We can't search for an exact match, because left neighbor base can be
  // anywhere. Instead, we search for the immediately smaller address in the
  // free pool from the compaction address.
  left = NULL;
  if (kt_trie_find_approx(pool->free_slabs, 0, adr, (void *) &left)) {
    // Exact match not possible: we're not part of the free pool
    ar_assert(left->address < adr);
    if (left->address + left->num_slabs * MM_SLAB_SIZE != adr) {
      // Left neighbor is far away, not suitable for compaction
      left = NULL;
    }
  }

  // Search free slabs for a neighbor on our right side (address above ours).
  // To be suitable, this has to be an exact match: it has to start exactly
  // on where our compaction range ends.
  right = NULL;
  kt_trie_find(pool->free_slabs, adr + num_slabs * MM_SLAB_SIZE, 
               (void *) &right);

  // Four combinations of left/right
  if (left && !right) {
    // Expand left neighbor to include the new range
    left->num_slabs += num_slabs;
  }
  else if (!left && right) {
    // Expand right neighbor to include the new range
    ar_assert(!kt_trie_delete(pool->free_slabs, right->address, NULL));
    right->address -= num_slabs * MM_SLAB_SIZE;
    right->num_slabs += num_slabs;
    ar_assert(!kt_trie_insert(pool->free_slabs, right->address, right));
  }
  else if (left && right) {
    // Expand left neighbor to include both the new range and the right
    // neighbor
    left->num_slabs += num_slabs + right->num_slabs;
    ar_assert(!kt_trie_delete(pool->free_slabs, right->address, kt_free));
  }
  else {
    // No neighbors; create new free slab for the new range
    slab = kt_zalloc(sizeof(MmSlab));
    slab->state = SLAB_EMPTY;
    slab->address = adr;
    slab->num_slabs = num_slabs;
    ar_assert(!kt_trie_insert(pool->free_slabs, adr, slab));
  }
}


// ===========================================================================
// mm_slab_store_chunk()        FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
void mm_slab_store_chunk(MmSlabPool *pool, int pack_options, 
                         int default_location, Trie *exception_locations, 
                         int *alloc_num_elements, int *num_elements, 
                         int **sizes, size_t **addresses, size_t *cur_exc_adr, 
                         int *cur_exc_location, size_t *chunk_adr, 
                         int *chunk_size) {

  size_t        new_adr;
  int           new_size;
  int           new_location;
  int           piece_size;


  do {

    // See if the old chunk we're storing has any exceptional objects
    // that have a different location. Cut the chunk into pieces in
    // this case, storing the exceptional objects and the chunk pieces
    // in order.
    if (*cur_exc_adr == *chunk_adr) {
      new_adr = *cur_exc_adr;
      new_size = mm_slab_query_pointer(pool, *cur_exc_adr);
      ar_assert(new_size > 0);
      new_location = *cur_exc_location;

      *cur_exc_adr = kt_trie_find_next(exception_locations, 1, 
                                       (void *) cur_exc_location);

      *chunk_adr += new_size;
      *chunk_size -= new_size;
    }
    else if ((*cur_exc_adr > *chunk_adr) && 
             (*cur_exc_adr < *chunk_adr + *chunk_size)) {
      new_adr = *chunk_adr;
      new_size = *cur_exc_adr - *chunk_adr;
      new_location = default_location;

      *chunk_adr += new_size;
      *chunk_size -= new_size;
    }
    else {
      new_adr = *chunk_adr;
      new_size = *chunk_size;
      new_location = default_location;

      *chunk_size = 0;
    }
//kt_printf("+ 0x%08X size %d loc %d || 0x%08X size %d\r\n", new_adr, new_size, new_location, *chunk_adr, *chunk_size);

    ar_assert(*chunk_size >= 0);
    ar_assert(new_size > 0);

    // Loop on this smaller piece, which may be needed to split again if
    // it's too big for the DMA engine.
    do {

      // Possibly break into pieces
      if (new_size < (1 << MM_PACK_SIZE_BITS)) {
        piece_size = new_size;
      }
      else {
        piece_size = (1 << MM_PACK_SIZE_BITS) - AR_DMA_ALIGN;
      }

      // Reallocate, if needed
      if (*num_elements >= *alloc_num_elements) {
        *alloc_num_elements *= 2;
        ar_assert(*sizes = kt_realloc(*sizes, 
                                      *alloc_num_elements * sizeof(int)));
        ar_assert(*addresses = kt_realloc(*addresses,
                                      *alloc_num_elements * sizeof(size_t)));
      }
      
      // Store
      ar_assert(piece_size < (1 << MM_PACK_SIZE_BITS));
      ar_assert(new_location < (1 << MM_PACK_LOCATION_BITS));
         
      (*addresses)[*num_elements] = new_adr;
      (*sizes)[*num_elements] = 
              piece_size | 
              (new_location << MM_PACK_SIZE_BITS) |
              (pack_options << (MM_PACK_SIZE_BITS + MM_PACK_LOCATION_BITS));
      (*num_elements)++;

      // Prepare for next small piece
      new_adr += piece_size;
      new_size -= piece_size;

    } while (new_size);

  } while (*chunk_size);
}

 
// ###########################################################################
// ###                                                                     ###
// ###                          EXPORTED FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################


// ===========================================================================
// mm_slab_create_pool()        Creates a new slab pool. Also, takes care
//                              of allocating bootstrapped objects for the
//                              special kernel slab pool.
// ===========================================================================
// * INPUTS
//   MmSlabPool *metadata_pool  Slab pool to allocate this pool's metadata.
//                              Can be NULL only for the kernel pool, when
//                              bootstrapping.
//   size_t adr                 Start address of slab pool initial free slabs
//   int num_slabs              Number of slab pool initial free slabs
//   int is_kernel              1 if this is the kernel slab pool
//
// * RETURN VALUE
//   MmSlabPool *               The newly allocated slab pool
// ===========================================================================
MmSlabPool *mm_slab_create_pool(MmSlabPool *metadata_pool, size_t adr, 
                                int num_slabs, int is_kernel) {

  Context       *context;
  int           my_cid;
  size_t        kernel_base;
  size_t        kernel_end;
  MmSlabPool    *pool;
  MmSlab        *slab;
  int           *bootstrap_slots;
  void          **free_slots;
  int           bootstrap_malloc_done[MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN];
  int           bootstrap_free_done[MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN];
  size_t        ret_adr;
  int           i;
  int           j;


  // Get global context
  my_cid = ar_get_core_id();
  context = mm_get_context(my_cid);
  
  // Allow non-existent metadata pool only when bootstrapping for the kernel
  // pool creation
  if (!metadata_pool) {
    ar_assert(is_kernel);
    if (!context->mm_alloc_bootstrap) {
      ar_abort();
    }
  }

  // Allocate the new pool
  pool = kt_malloc(sizeof(MmSlabPool));
  pool->is_kernel = is_kernel;
  pool->metadata_pool = metadata_pool;
  pool->used_slabs = kt_alloc_trie(MM_TRIES_POINTER_MSB, MM_TRIES_POINTER_LSB);
  pool->partial_slabs = kt_alloc_trie(MM_TRIES_SLOT_MSB, MM_TRIES_SLOT_LSB);
  pool->free_slabs = kt_alloc_trie(MM_TRIES_POINTER_MSB, MM_TRIES_POINTER_LSB);
  pool->used_space = 0;
  pool->free_space = num_slabs * MM_SLAB_SIZE;

  // Push all given slabs to the free slabs trie
  ar_assert(adr % MM_SLAB_SIZE == 0);
  slab = kt_zalloc(sizeof(MmSlab));
  slab->state = SLAB_EMPTY;
  slab->address = adr;
  slab->num_slabs = num_slabs;
  ar_assert(!kt_trie_insert(pool->free_slabs, adr, slab));


  // If bootstrapping, we need to normally allocate all the stuff we're using
  if (!metadata_pool) {

    // Find kernel limits
    kernel_base = mm_va_kernel_base(my_cid);
    kernel_end = kernel_base + MM_KERNEL_SIZE - 1024 * 1024;

    // First of all, do a dummy malloc and free for each possible bootstrap
    // size we are supporting. This will initialize the
    // MM_BOOTSTRAP_SLABS_STEP slabs per slot size, so the addresses that
    // kt_malloc() returns during bootstrap will be matching the ones we
    // expect here.
    for (i = 0; i < MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN; i++) {
      kt_free(kt_malloc((i + 1) * MM_ALLOC_ALIGN));
    }

    // Initialize our own bootstrapping completion arrays. These are needed
    // because we're still at bootstrap mode, so any call to kt_malloc() and
    // kt_free (by slab or trie functions) will allocate/free even more stuff.
    // We'll track how much we have compensated for by comparing this array to
    // the ones kt_malloc() and kt_free() keep on the last 2 kernel heap pages.
    for (i = 0; i < MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN; i++) {
      bootstrap_malloc_done[i] = 0;
      bootstrap_free_done[i] = 0;
    }


    // For each allocated object during bootstrap, call slab_alloc_slot() to
    // track the address we already have given it.
    // 
    // The loop will converge eventually, because slab allocation happens
    // for multi-KB slabs, whereas Trie nodes are created per slab. So, most
    // of the time the slab_alloc_slot() call will manage to return without
    // calling Trie node mallocs() -- it'll simply toggle bits on partial
    // slabs.
    bootstrap_slots = (int *) kernel_end;
    do {
      for (i = 0; i < MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN; i++) {
        for (j = bootstrap_malloc_done[i]; j < bootstrap_slots[i]; j++) {

          // Allocate the object using the same address
          ar_assert(!mm_slab_alloc_slot(pool, 
                                 (i + 1) * MM_ALLOC_ALIGN, 
                                 &ret_adr));
          ar_assert(ret_adr == kernel_base + 
                               i * MM_SLAB_SIZE * MM_BOOTSTRAP_SLABS_STEP +
                               j * (i + 1) * MM_ALLOC_ALIGN);

          // Remember that this object is done
          bootstrap_malloc_done[i]++;
        }
      }

      // Check that no new objects appeared on small slots while we were
      // doing the bigger slots
      j = 0;
      for (i = 0; i < MM_BOOTSTRAP_MAX_SLOT / MM_ALLOC_ALIGN; i++) {
        if (bootstrap_malloc_done[i] != bootstrap_slots[i]) {
          j = 1;
          break;
        }
      }
    } while (j);

      
    // All mallocs are reconciled; Assume bootstrapping finished. Further
    // calls to malloc() will be done properly from now on.
    context->mm_kernel_pool = pool;
    context->mm_alloc_bootstrap = 0;

    // Finally, allocate the deferred frees array (may require some frees()
    // which are still bootstrapped)
    ar_assert(context->mm_defer_frees = kt_malloc(MM_MAX_DEFERRED_FREES *
                                                  sizeof(size_t)));

    // With the deferred array allocated, we're ready to switch to normal free
    // mode
    context->mm_frees_bootstrap = 0;

    // We still owe to free objects that were freed during bootstrap. We don't
    // have to enter a loop here; frees() may call further mallocs(), but
    // we don't care, they will be handled properly.
    free_slots = (void **) ((size_t) kernel_end - MM_PAGE_SIZE +
                                     sizeof(int));
    for (i = 0; 
         i < *((int *) ((size_t) kernel_end - MM_PAGE_SIZE)); 
         i++) {

      // Free this pointer
      kt_free(free_slots[i]);
    }
  }

  // Success
  return pool;
}


// ===========================================================================
// mm_slab_destroy_pool()       Frees all objects in a slab pool and posisbly
//                              also frees the pool itself
// ===========================================================================
// * INPUTS
//   int complete_destruction   0: free all objects, destroy full, partial and
//                                 preallocated empty slabs but move them to
//                                 the free_slabs trie. This frees the pool
//                                 memory but leaves the pool intact to reuse.
//                              1: completely free everything, including the
//                                 pool free_slabs and the pool object itself.
//
// * INOUTS
//   MmSlabPool *pool           The pool to be emptied or destroyed completely
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_slab_destroy_pool(MmSlabPool *pool, int complete_destruction) {

  size_t        chunk_adr;
  int           chunk_num_slabs;
  MmSlab        *slab;
  size_t        slab_adr;
  size_t        prev_slab_adr;
  int           num_slabs;
  size_t        size;
  size_t        prev_size;


  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->used_slabs);
  ar_assert(pool->partial_slabs);
  ar_assert(pool->free_slabs);

  // Complete destruction is the easiest
  if (complete_destruction) {

    // Remove partial list
    kt_free_trie(pool->partial_slabs, kt_free);

    // Remove all slabs and the used trie
    kt_free_trie(pool->used_slabs, kt_free);

    // Remove any free slabs and the free trie
    kt_free_trie(pool->free_slabs, kt_free);

    // Free ourselves
    kt_free(pool);

    return 0;
  }

  // Remove partial list heads, we won't be needing them
  prev_size = 0;
  for (size = kt_trie_find_minmax(pool->partial_slabs, 0, NULL);
       (size || prev_size);
       size = kt_trie_find_next(pool->partial_slabs, 1, NULL)) {

    // Delete previous node
    if (prev_size) {
      ar_assert(!kt_trie_delete(pool->partial_slabs, prev_size, kt_free));
    }

    // Finished?
    if (!size) {
      break;
    }

    // Mark current node to be deleted after we walk to the next one
    // (we can't delete it now, trie won't know how to walk to the next one)
    prev_size = size;
  }


  // Walk all used slabs nodes in ascending order
  chunk_adr = 0;
  chunk_num_slabs = 0;
  prev_slab_adr = 0;
  for (slab_adr = kt_trie_find_minmax(pool->used_slabs, 0, (void *) &slab);
       (slab_adr || prev_slab_adr);
       slab_adr = kt_trie_find_next(pool->used_slabs, 1, (void *) &slab)) {

    // Delete previous node
    if (prev_slab_adr) {
      ar_assert(!kt_trie_delete(pool->used_slabs, prev_slab_adr, kt_free));
    }

    // Finished?
    if (!slab_adr) {
      // Add remaining chunk to free pool
      ar_assert(chunk_adr);
      mm_slab_add_to_free(pool, chunk_adr, chunk_num_slabs);
      break;
    }

    // How many slabs does this one represent?  If slab->num_slabs is 1, it
    // means the last slab slot does not extend into any other slab (so we
    // count it for 1 here). If it's > 1, it means it uses slab->num_slabs - 1
    // adjacent slabs.
    //
    // The question is: are we going to find the last one as a new used_slabs
    // entry? If the last slot is aligned to end on the end of the last
    // adjacent slab, we won't: we have to count it here.  Otherwise, we'll
    // count it when we find the next entry.
    ar_assert(slab->num_slabs > 0);
    if (slab->num_slabs == 1) {
      num_slabs = 1;
    }
    else if ((slab->slot_size - MM_SLAB_SIZE - slab->offset -
              (slab->bitmap_size - 1) * slab->slot_size) % MM_SLAB_SIZE) {
      num_slabs = slab->num_slabs - 1;
    }
    else {
      num_slabs = slab->num_slabs;
    }

    // Try to create a chunk of many slabs, in order not to do as few
    // free_slabs additions as possible
    if (!chunk_adr) {
      ar_assert(slab_adr == slab->address);
      chunk_adr = slab_adr;
      chunk_num_slabs = num_slabs;

    }
    else {
      // slab is a continuation of existing chunk?
      if (slab_adr == chunk_adr + chunk_num_slabs * MM_SLAB_SIZE) {
        ar_assert(slab_adr == slab->address);
        chunk_num_slabs += num_slabs;
      }
      // chunk ended; process it and begin a new one
      else {
        
        // Add old chunk to free pool
        mm_slab_add_to_free(pool, chunk_adr, chunk_num_slabs);

        // Begin new chunk with this slab
        ar_assert(slab_adr == slab->address);
        chunk_adr = slab_adr;
        chunk_num_slabs = num_slabs;
      }
    }

    // Mark current node to be deleted after we walk to the next one
    // (we can't delete it now, trie won't know how to walk to the next one)
    prev_slab_adr = slab_adr;
  }

  // All is well, but the space variables in the pool are now inconsistent.
  // Fix them here.
  pool->used_space = 0;
  pool->free_space = 0;
  for (slab_adr = kt_trie_find_minmax(pool->free_slabs, 0, (void *) &slab);
       slab_adr;
       slab_adr = kt_trie_find_next(pool->free_slabs, 1, (void *) &slab)) {
    pool->free_space += slab->num_slabs * MM_SLAB_SIZE;
  }

  // Success
  return 0;
}


// ===========================================================================
// mm_slab_preallocate()        Handles the kernel slab pool preallocation,
//                              i.e. for all needed slot sizes that are low on
//                              slabs and have been marked, it fetches from
//                              the free pool and assigns them as empty ones.
// ===========================================================================
void mm_slab_preallocate() {

  Context               *context;
  MmSlab                *free_slab;
  size_t                free_slab_adr;
  MmSlabPool            *pool;
  int                   s;
  int                   i;


  // Get context
  context = mm_get_context(ar_get_core_id());
  pool = context->mm_kernel_pool;

  // For all slot sizes...
  for (s = MM_ALLOC_ALIGN; s < MM_SLAB_SIZE; s += MM_ALLOC_ALIGN) {

    // ... that actually need preallocation
    if (!context->mm_prealloc_flags[s / MM_ALLOC_ALIGN]) {
      continue;
    }

    // Find MM_BOOTSTRAP_SLABS_STEP + 1 more
    for (free_slab_adr = kt_trie_find_minmax(pool->free_slabs, 0, 
                                             (void *) &free_slab);
         free_slab_adr;
         free_slab_adr = kt_trie_find_next(pool->free_slabs, 1,
                                           (void *) &free_slab)) {
      if (free_slab->num_slabs > MM_BOOTSTRAP_SLABS_STEP) {
        break;
      }
    }
    if (!free_slab_adr) {
      // Not enough kernel memory: It's a fine time to panic.
      ar_panic("Out of kernel memory");
    }

    // Create a new slab for each of them
    for (i = MM_BOOTSTRAP_SLABS_STEP - 1; i >= 0; i--) {
      mm_slab_alloc_new(pool, 0, 0, free_slab->address + i * MM_SLAB_SIZE,
                        s, 0, 1, 1);
    }

    // Update free slab with the number of slabs it's lost. It can't be
    // exhausted (we asked MM_BOOTSTRAP_SLABS_STEP + 1)
    free_slab->num_slabs -= MM_BOOTSTRAP_SLABS_STEP;
    ar_assert(free_slab->num_slabs > 0);

    // Slab start address itself is also changing, so we also need to
    // update the trie with the new key
    ar_assert(!kt_trie_delete(pool->free_slabs, free_slab->address, 
                              NULL));
    free_slab->address += MM_BOOTSTRAP_SLABS_STEP * MM_SLAB_SIZE;
    ar_assert(!kt_trie_insert(pool->free_slabs, free_slab->address, 
                              free_slab));

    // Mark that preallocation is finished for this one
    context->mm_prealloc_flags[s / MM_ALLOC_ALIGN] = 0;
  }

  // Finished preallocation in general
  context->mm_prealloc_needed = 0;
}


// ===========================================================================
// mm_slab_alloc_slot()         Allocate a slot in a slab (basic allocation
//                              function). 
//
//                              Finds an appropriate slab in the slab pool
//                              and returns the address of a free slot. 
//                              Searches first among partial slabs and if
//                              nothing comes up, searches among free slabs.
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool for the new slot
//   int size                   Slot size for the allocation
//
// * OUTPUTS
//   size_t *ret_adr            Return address of the allocated slot
//
// * RETURN VALUE
//   int                        0: success (ret_adr is valid)
//                              ERR_OUT_OF_MEMORY: out of pool memory (ret_adr 
//                                 invalid).No suitable partial or free slabs
//                                 could be found.
// ===========================================================================
int mm_slab_alloc_slot(MmSlabPool *pool, int size, size_t *ret_adr) {

  Context               *context;
  MmSlabPool            *pool_save;
  size_t                min_partial_adr;
  MmSlabPartialHead     *head;
  size_t                head_adr;
  MmSlab                *slab;
  MmSlab                *free_slab;
  MmSlab                *empty_slab;
  size_t                free_slab_adr;
  size_t                new_slab_adr;
  int                   bit;
  unsigned int          slot_start;
  int                   num_slabs;
  int                   free_num_slabs;
  int                   take_from_bottom;
  int                   rem_size;
  size_t                adr;


  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->used_slabs);
  ar_assert(pool->partial_slabs);
  ar_assert(pool->free_slabs);
  ar_assert(!(size & (MM_ALLOC_ALIGN - 1)));

  // For kernel pool only, mark the recursion depth of kernel mallocs. We
  // depend on having enough preallocated empty slabs of small sizes so that
  // we can fullfil request without calling any more kernel mallocs which have
  // the ability to steal the pointer we're trying to acquire.
  context = mm_get_context(ar_get_core_id());
  if (pool->is_kernel) {
    ar_assert(context->mm_recursion_depth >= 0);
    context->mm_recursion_depth++;
  }

  // All allocations/frees done here should be done on pool's metadata pool.
  // Push it on the "pool stack" (single pool_save per function recursion).
  pool_save = context->mm_kernel_pool;
  context->mm_kernel_pool = pool->metadata_pool;
  if (!context->mm_kernel_pool) {
    // Exception: stop at the lowest pool recursion level
    ar_assert(pool->is_kernel);
    context->mm_kernel_pool = pool_save;
  }


  // =========================================================================
  // Attempt 1: Find a partial slab that can allocate this request
  // =========================================================================

  // Are there partial slabs with this slot size?
  min_partial_adr = 0;
  head_adr = kt_trie_find(pool->partial_slabs, size, (void *) &head);
  if (head_adr) {

    // Explore the whole linked list
    ar_assert(head->head);
    ar_assert(!head->head->prev);
    for (slab = head->head; slab; slab = slab->next) {

      ar_assert(slab->state != SLAB_FULL);
      ar_assert(slab->slot_size == size);

      // Gather statistics for free slab allocation, if partial search fails
      if ((!min_partial_adr) || (min_partial_adr > slab->address)) {
        min_partial_adr = slab->address;
      }

      // It's a partial or empty slab, so there is a slot available. Where?
      bit = mm_slab_find_free_bit(slab->bitmap, slab->bitmap_size);
      if (bit == -1) {
        // We gave priority to finding a "real" empty slot, but we failed.
        // Verify that this is the case that only the trailer bit is left free.
        ar_assert(!slab->free_slots);
        ar_assert(slab->offset);
        ar_assert(slab->trailer_bit);
        
        // We're dealing here with yet another malevolent case. We should see
        // if previously adjacent free slab(s) can cover. In essence, this
        // guarantees that we can handle cases where mallocs can happen in
        // reverse order (start by something that has a trailer free, and
        // allocate backwards). Note that if a partial slab can accomodate, it
        // will be covered normally while searching (outside this if scope
        // below) on other partial slabs of the list which will find this
        // trailer because they'll look for a forward empty trailer.
        rem_size = size - slab->offset;
        ar_assert(rem_size > 0);
        
        // We only allow to grab free slabs only when we're not recursing 
        // for kernel mallocs
        if (context->mm_recursion_depth > 1) {
          continue;
        }

        // We need rem_size / MM_SLAB_SIZE free slabs, rounded up
        num_slabs = (rem_size + MM_SLAB_SIZE - MM_ALLOC_ALIGN) / MM_SLAB_SIZE;
        ar_assert(slab->address >= num_slabs * MM_SLAB_SIZE);

        // Find them: search for the immediately previous free slab to our
        // address
        free_slab_adr = kt_trie_find_approx(pool->free_slabs, 0, slab->address,
                                            (void *) &free_slab);

        // Did we fail?
        if (!free_slab_adr) {
          // No free previous slab
          continue;
        }
        ar_assert(slab->address != free_slab_adr); // slab can't be in free pool
        ar_assert(free_slab->address == free_slab_adr);
        if ((free_slab->num_slabs < num_slabs) ||  // has enough slabs?
            (free_slab_adr + 
              free_slab->num_slabs * MM_SLAB_SIZE !=
             slab->address)) {                     // are they adjacent?
          // Free slab exists, but either it's not
          // big enough, or it's not adjacent
          continue;
        }

        // Claim the trailer on the partial slab
        ar_assert(!mm_slab_claim_trailer(pool, head, slab, 0, size, 
                                         slab->offset));

        // Take the rest from the high part of the free slabs we found
        new_slab_adr =  mm_slab_take_from_free(pool, free_slab, num_slabs, size,
                                               slab->offset, 1, 0);

        // That's it
        if (ret_adr) {
          *ret_adr = new_slab_adr;
        }
        goto success;
      }

      // Does it fit? If it's the last slot (and the slab size is not a
      // multiple of this slot size, adjusted also for slab offset), it may
      // need adjacent slab(s) to fit.
      slot_start = slab->offset + bit * size;
      ar_assert(slot_start < MM_SLAB_SIZE);
      free_slab = NULL;
      empty_slab = NULL;
      num_slabs = 0;
      rem_size = 0;
      if (slot_start + size > MM_SLAB_SIZE) {
        ar_assert(bit == slab->bitmap_size - 1);

        // We need adjacent slabs. Specifically, we need 
        // (remaining size / MM_SLAB_SIZE), rounded up, adjacent slabs.
        // Note that the order we search below (first free and then partials
        // or empty) is not important, since we're looking for exact
        // addresses: the needed slabs will be either in free_slabs or in
        // used_slabs, but not both.
        rem_size = size - (MM_SLAB_SIZE - slot_start);
        ar_assert(rem_size > 0);
        num_slabs = (rem_size + MM_SLAB_SIZE - MM_ALLOC_ALIGN) / MM_SLAB_SIZE;

        // Search first among the free slabs. If they exist, they will be at a
        // single chunk addressed exactly at the point we say, since we are
        // currently at a partial slab which ends at this address.
        //
        // Again, we take from free slabs only if we're not recursing for
        // kernel mallocs.
        if (context->mm_recursion_depth <= 1) {
          free_slab_adr = kt_trie_find(pool->free_slabs, 
                                       slab->address + MM_SLAB_SIZE,
                                       (void *) &free_slab);
        }
        else {
          free_slab_adr = 0;
        }

        if (free_slab_adr) {
          // Free slab exists at this address, but how many slabs can it 
          // provide? there are 3 subcases.
          ar_assert(num_slabs >= 1);
          ar_assert(free_slab->num_slabs >= 1);

          // Subcase 1: free slab has enough slabs to cover all our trailer.
          if (free_slab->num_slabs >= num_slabs) {
            // Allocate the remainder of the slot using space from the free 
            // pool. Create a new slab for the trailer, if the remainder of
            // the slot does not consume the full slab.
            ar_assert(!mm_slab_take_from_free(pool, free_slab, num_slabs, size, 
                                              size - rem_size, 0, 1));
            
            // The original slab just increased its num_slabs
            ar_assert(slab->num_slabs == 1);
            slab->num_slabs += num_slabs;
          }

          // Subcase 2: free slab has enough slabs to cover all but the last
          //            slab of our trailer. This smells like there is a 
          //            partial slab somewhere with a compatible offset that
          //            can hold this last slab.
          else if (free_slab->num_slabs == num_slabs - 1) {
            ar_assert(num_slabs >= 2);

            // Search for this mysterious slab with compatible slot_size and
            // offset to hold this last part. If such a slab is found, claim
            // its trailer.
            if (mm_slab_claim_trailer(pool, 
                                      head, 
                                      NULL,
                                      free_slab->address + 
                                          (num_slabs - 1) * MM_SLAB_SIZE,
                                      size, 
                                      rem_size - 
                                          (num_slabs - 1) * MM_SLAB_SIZE)) {
              // Failure; search for other partials
              continue;
            }

            // The free slabs are now intermediates between slab (whose slot
            // will be claimed later on) and a trailer (already claimed on a
            // previously partial/empty slab). These intermediate slabs must
            // not appear anywhere. Delete them from the free pool and forget
            // about them.
            ar_assert(!kt_trie_delete(pool->free_slabs, free_slab->address, 
                                      kt_free));
              
            // The original slab just increased its num_slabs
            ar_assert(slab->num_slabs == 1);
            slab->num_slabs += num_slabs;
          }

          // Subcase 3: free pool can cover 2 slabs or even less than what we
          // need. We can't do anything; it's a large slot size (since we
          // need >= 2 slabs) and we don't keep empty preallocated to cover
          // the difference.
          else {
            ar_assert(free_slab->num_slabs < num_slabs - 1);
            // Failure; search for other partials
            continue;
          }
        }
        
        // We need to search the partial slabs for the specific adjacent
        // slabs. They may be partial or even empty (preallocated or emptied
        // due to frees that have not compacted the free space).
        else {

          // We do that only for small slot sizes; larger slot sizes do not
          // keep empty slabs in the partial list.
          if (rem_size >= MM_SLAB_SIZE) {
            // Failure; search for other partials
            continue;
          }
          ar_assert(num_slabs == 1);
          
          // Search now among all used slabs, looking for a partial or empty
          // of compatible slot_size and offset to hold rem_size bytes.
          // If a slab is found, claim its trailer.
          if (mm_slab_claim_trailer(pool, head, NULL, 
                                    slab->address + MM_SLAB_SIZE,
                                    size, rem_size)) {
            // Failure; search for other partials
            continue;
          }
            
          // The original slab just increased its num_slabs
          ar_assert(slab->num_slabs == 1);
          slab->num_slabs += num_slabs;
        }
      }

      // Ok, slab fits and its trailer (if any) has been taken care of (i.e.
      // we either claimed it on a partial slab, or we took space from a free
      // slab and claimed it there, or both of them). Now, claim the slot in
      // the original slab.
      mm_slab_set_used_bit(slab->bitmap, slab->bitmap_size, bit);
      slab->free_slots--;
      head->used_slots++;
      head->free_slots--;
      ar_assert(slab->free_slots >= 0);
      ar_assert(head->free_slots >= 0);

      // If any slots (or simply the trailer) are still available
      if (!mm_slab_is_full(slab)) {

        // Mark it PARTIAL (it may have been EMPTY)
        if (slab->state == SLAB_EMPTY) {
          slab->state = SLAB_PARTIAL;
          head->num_empty_members--;
        }

        // Since we found space on this slab, make it the head of the 
        // linked list. Future searches may benefit from this.
        if (head->head != slab) {
          mm_slab_swap_list_nodes(head->head, slab);
          head->head = slab;
        }
      }

      // If nothing more can be allocated on this slab...
      else {

        // ... it's not partial anymore. Change status to full and remove
        // from the partial linked list.
        mm_slab_promote_to_full(pool, head, slab);
      }
        
      // We're done
      if (ret_adr) {
        *ret_adr = slab->address + slot_start;
      }
      goto success;
    }
  }

  // In kernel malloc recursions, if we didn't find anything we have a
  // problem: not enough preallocated empty slabs were there to cover the
  // requests. Increase MM_KERNEL_MIN_EMPTY, MM_BOOTSTRAP_SLABS_STEP or
  // MM_BOOTSTRAP_MAX_SLOT depending on what failed.
  if (context->mm_recursion_depth > 1) {
    ar_abort();
  }


  // =========================================================================
  // Attempt 2: Find a free slab that can allocate this request
  // =========================================================================

  // Generally, we need (size / MM_SLAB_SIZE), rounded up, adjacent slabs 
  num_slabs = (size + MM_SLAB_SIZE - MM_ALLOC_ALIGN) / MM_SLAB_SIZE;

  // If we're the kernel pool, we need at least MM_BOOTSTRAP_SLABS_STEP
  // slabs for small slot sizes, so we can always preallocate enough mem
  if ((pool->is_kernel) && (size < MM_SLAB_SIZE) && 
      (num_slabs < MM_BOOTSTRAP_SLABS_STEP)) {
    free_num_slabs = MM_BOOTSTRAP_SLABS_STEP;
  }
  else {
    free_num_slabs = num_slabs;
  }


  // There are no suitable partial slabs, so search for free_num_slabs adjacent
  // free slabs. In order to keep the pool addresses as packed as possible, if
  // we encountered so far partial (but unsuitable) slabs, prefer free ones
  // near them. If there are not, just get the minimum-addressed free slabs we
  // own.
  take_from_bottom = 1;
  if (min_partial_adr) {
    free_slab_adr = kt_trie_find_approx(pool->free_slabs, 1, 
                                   min_partial_adr, (void *) &free_slab);
    if ((!free_slab_adr) || (free_slab->num_slabs < free_num_slabs)) {
      free_slab_adr = kt_trie_find_approx(pool->free_slabs, 0, 
                                     min_partial_adr, (void *) &free_slab);
      take_from_bottom = 0;
      if ((!free_slab_adr) || (free_slab->num_slabs < free_num_slabs)) {
        goto find_any;
      }
    }
  }
  else {

find_any:

    // Near-matches (if any) failed us above. Search every possible chunk.
    take_from_bottom = 1;
    for (free_slab_adr = kt_trie_find_minmax(pool->free_slabs, 0, 
                                             (void *) &free_slab);
         free_slab_adr;
         free_slab_adr = kt_trie_find_next(pool->free_slabs, 1,
                                           (void *) &free_slab)) {
      if (free_slab->num_slabs >= free_num_slabs) {
        break;
      }
    }
    if (!free_slab_adr) {
      // No suitable slabs could be found anywhere
      if (pool->is_kernel) {
        // Not enough kernel memory: It's also a good time to panic.
        ar_panic("Out of kernel memory");
      }
      context->mm_kernel_pool = pool_save;
      return ERR_OUT_OF_MEMORY;
    }
  }

  // Sanity checks
  ar_assert(free_slab);
  ar_assert(free_slab->state == SLAB_EMPTY);

  // Take a chunk from the free pool to allocate a slot starting at offset 0.
  // This will create 1 or 2 new slabs: one (always) for the start of the slot
  // and a second one (sometimes) for its trailer, if the number of slabs is
  // more than 1 and the slot size is not a multiple of the slab size.
  new_slab_adr = mm_slab_take_from_free(pool, free_slab, num_slabs, size, 
                                        0, 0, take_from_bottom);

  // Verify and store the return address (of the 1st new slab)
  ar_assert(new_slab_adr);
  if (ret_adr) {
    *ret_adr = new_slab_adr;
  }

success:

  // Just before we return, make sure that we have kept enough empty
  // preallocated slabs (only for the kernel mem pool). This is a
  // safeguard to avoid infinite recursion and/or weird inconsistencies,
  // where we select a pointer but a Trie malloc() recurses and steals it.
  // 
  // WARNING: Correct behavior is guaranteed only when:
  // - Enough empty slabs are preallocated to support any single allocation
  // - All the code is written in a way that when enough empty slabs are
  //   kept, then no malloc() is needed for the same slot size. This is
  //   true for the code here, because a partial hit turns a bit on a slab
  //   and maybe swaps the partial list nodes or even promotes a slab
  //   to be full. None of these should ever require further mallocs.
  if ((pool->is_kernel) && (size < MM_SLAB_SIZE)) {

    // Even after allocations, we must not be drained
    head_adr = kt_trie_find(pool->partial_slabs, size, (void *) &head);
    ar_assert(head);
    ar_assert(head->num_empty_members > 0);
    
    // Do we have enough?
    if (head->num_empty_members < MM_KERNEL_MIN_EMPTY) {

      // This should not happen during bootstrapping. If this assertion
      // fails, just increase MM_BOOTSTRAP_SLABS_STEP. 
      // Actually, what happens is that although the bootstrap slabs
      // are marginally okay (margin is less than MM_KERNEL_MIN_EMPTY,
      // however), the preallocation we do here inserts the empties after
      // the partial head. When bootstrap reconcilation happens, these
      // new empties are grabbed and the pointer we return from them 
      // does not match the one supposed during bootstrap.
      ar_assert(!context->mm_alloc_bootstrap);

      // Mark that we need preallocation (generally), on this slot size
      // (specifically)
      context->mm_prealloc_needed = 1;
      context->mm_prealloc_flags[size / MM_ALLOC_ALIGN] = 1;
    }
  }

  // The actual preallocation must be done only when all the recursions have
  // been finished, as well as we're not in the middle of freeing anything.
  // This will guarantee that no pointers can be stolen from intermediate Trie
  // operations.
  if ((pool->is_kernel) && (context->mm_recursion_depth == 1) &&
      (!context->mm_busy_freeing) && (context->mm_prealloc_needed)) {
    mm_slab_preallocate();
  }

  // Pool bookkeeping
  ar_assert(pool->free_space >= size);
  pool->used_space += size;
  pool->free_space -= size;
      
  // Decrease kernel recursion depth counter
  if (pool->is_kernel) {
    context->mm_recursion_depth--;
    ar_assert(context->mm_recursion_depth >= 0);

    // If we finished with all the nested mallocs, maybe we owe some frees to
    // do? This may happen because mm_slab_free_slot() postpones all
    // operations while we are active.
    if (!context->mm_recursion_depth && context->mm_num_defer_frees &&
        !context->mm_busy_freeing) {

      // Pop the first one from the stack and call mm_slab_free_slot() to
      // process them all
      adr = context->mm_defer_frees[--context->mm_num_defer_frees];
      ar_assert(!mm_slab_free_slot(pool, adr));
    }
  }

  // Pop pool stack
  context->mm_kernel_pool = pool_save;

  // Success
  return 0;
}


// ===========================================================================
// mm_slab_free_slot()          Free a slot in a slab (basic free function). 
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool for the new slot
//   size_t *adr                The slot pointer to be freed
//
// * RETURN VALUE
//   int                        0: success
//                              ERR_MISALIGNED: pointer not aligned to 
//                                 MM_ALLOC_ALIGN or to slab slot
//                              ERR_OUT_OF_RANGE: pointer not inside this pool
//                                 used/partials slabs
//                              ERR_NOT_ALLOCED: pointer in range, but mapped 
//                                 to empty slot or free slab
// ===========================================================================
int mm_slab_free_slot(MmSlabPool *pool, size_t adr) {

  Context               *context;
  MmSlabPool            *pool_save;
  MmSlab                *slab;
  int                   slab_previously_full;
  MmSlab                *slab2;
  int                   slab2_previously_full;
  int                   bit;
  MmSlabPartialHead     *head;
  int                   num_slabs;
  size_t                compact_start_adr;
  int                   compact_num_slabs;
  int                   ret;


  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->used_slabs);
  ar_assert(pool->partial_slabs);
  ar_assert(pool->free_slabs);
  ar_assert(adr);

  // All allocations/frees done here should be done on pool's metadata pool
  context = mm_get_context(ar_get_core_id());
  pool_save = context->mm_kernel_pool;
  context->mm_kernel_pool = pool->metadata_pool;
  if (!context->mm_kernel_pool) {
    // Exception: stop at the lowest pool recusrion level
    ar_assert(pool->is_kernel);
    context->mm_kernel_pool = pool_save;
  }

  // Are we already in the middle of freeing something else in the kernel?
  // Or, more weirdly, are we in the middle of allocating somethign in the
  // kernel? In both cases, we must not proceed with the actual freeing, since
  // it may mess up Tries that are in the middle of updating.
  if (pool->is_kernel) {
    if (context->mm_busy_freeing || context->mm_recursion_depth) {
      // Store this free for later
      ar_assert(context->mm_defer_frees);
      ar_assert(context->mm_num_defer_frees < MM_MAX_DEFERRED_FREES);
      context->mm_defer_frees[context->mm_num_defer_frees++] = adr;
      return 0;
    }
    else {
      context->mm_busy_freeing = 1;
    }
  }

begin:

  // =========================================================================
  // Free slot: locate the slot, free it and also free the possible trailer
  //            on a neighbouring slab.
  // =========================================================================

  // Find the slab and bit position this pointer refers to
  ret = mm_slab_locate_pointer(pool, adr, &slab, &bit);
  if (ret) {
    context->mm_kernel_pool = pool_save;
    return ret;
  }

  // Locate head of the partial list for this slot size
  if (!kt_trie_find(pool->partial_slabs, slab->slot_size, (void *) &head)) {
    head = NULL;
    // No head implies that all slabs (including the one we found) are full
    ar_assert(slab->state == SLAB_FULL);
  }

  // Set slot free
  mm_slab_set_free_bit(slab->bitmap, slab->bitmap_size, bit);
  ar_assert(slab->free_slots < slab->bitmap_size);
  slab->free_slots++;

  // Pool bookkeeping
  pool->used_space -= slab->slot_size;
  pool->free_space += slab->slot_size;

  // Change state
  if (slab->state == SLAB_FULL) {
    // This takes care of the state change and also manages the partial list
    // insertion (head is ok to be NULL; the function will create it, if
    // needed).
    head = mm_slab_demote_from_full(pool, head, slab);
    slab_previously_full = 1;
  }
  else if (slab->state == SLAB_PARTIAL) {
    ar_assert(head);
    head->used_slots--;
    head->free_slots++;

    if (mm_slab_is_empty(slab)) {
      // Just change state, leave partial list as it is
      slab->state = SLAB_EMPTY;
      head->num_empty_members++;
      ar_assert(head->num_empty_members <= head->num_members);
    }
    slab_previously_full = 0;
  }
  else {
    // We just freed a slot in an empty slab...
    ar_abort();
    slab_previously_full = 0;
  }

  // Are we're freeing a last slot, which has implications beyond this slab?
  ar_assert(slab->bitmap_size > 0);
  ar_assert(slab->num_slabs > 0);
  slab2 = NULL;
  slab2_previously_full = 0;
  if ((bit == slab->bitmap_size - 1) && (slab->num_slabs > 1)) {

    // Remember this number of slabs
    num_slabs = slab->num_slabs;

    // This slab does not anymore require adjacent slabs
    slab->num_slabs = 1;

    // There might be a trailer some slabs ahead which must be freed as well
    if ((slab->slot_size - 
        (MM_SLAB_SIZE - slab->offset -  // Remainder of the slot not 
         bit * slab->slot_size)) %      // aligned to slab size
       MM_SLAB_SIZE) {                  

      // Find the slab with the trailer
      ar_assert(kt_trie_find(pool->used_slabs, 
                             slab->address + (num_slabs - 1) * MM_SLAB_SIZE,
                             (void *) &slab2));
      ar_assert(slab2->state != SLAB_EMPTY);
      ar_assert(slab2->offset);
      ar_assert(!slab2->trailer_bit);

      // Set the trailer free
      slab2->trailer_bit = 1;

      // Possibly change state, as we did above for slab
      if (slab2->state == SLAB_FULL) {
        head = mm_slab_demote_from_full(pool, head, slab2);
        slab2_previously_full = 1;
      }
      else if (slab2->state == SLAB_PARTIAL) {
        ar_assert(head);
        if (mm_slab_is_empty(slab2)) {
          slab2->state = SLAB_EMPTY;
          head->num_empty_members++;
          ar_assert(head->num_empty_members <= head->num_members);
        }
        slab2_previously_full = 0;
      }
      else {
        ar_abort();
      }

      // Slot began at slab and ended at slab2. This means that num_slabs - 2 
      // adjacent previously full slabs are now empty, but unaccounted for.
      // Remember this number of "intermediate" (i.e. after slab and before
      // slab2) number of slabs that must be freed.
      num_slabs -= 2;
    }
    else {
      // Slot began at slab and neatly ended at the end of an adjacent
      // MM_SLAB_SIZE boundary. This means that num_slabs - 1 adjacent
      // previously full slabs are now unaccounted for.
      num_slabs -= 1;
    }
  }

  // We freed a non-last slot, or a last one which does not extend into
  // adjacent slabs. We have no intermediate slabs to be freed.
  else {
    num_slabs = 0;
  }


  // =========================================================================
  // Compacting phase: think what we must do with the up to two empty slabs
  //                   (the slot one and the trailer one) and/or the emptied,
  //                   unaccounted, intermediate slabs.
  // =========================================================================

  compact_start_adr = 0;
  compact_num_slabs = 0;

  // Two big cases: for small slot sizes we compact only if more than
  // MM_BOOTSTRAP_SLABS_STEP empty slabs are in the partial list.
  if (slab->slot_size < MM_SLAB_SIZE) {

    // No intermediate slabs are possible here: slot size is small.
    ar_assert(!num_slabs);

    // Compact slab?
    if ((slab->state == SLAB_EMPTY) &&
        (head->num_empty_members > MM_BOOTSTRAP_SLABS_STEP)) {
      compact_start_adr = slab->address;
      compact_num_slabs = 1;
      mm_slab_destroy_empty(pool, slab, head, 1);
    }

    // Compact slab2?
    if ((slab2) && (slab2->state == SLAB_EMPTY) &&
        (head->num_empty_members > MM_BOOTSTRAP_SLABS_STEP)) {
      if (compact_start_adr) {
        // slab and slab2 are adjacent
        ar_assert(slab2->address == compact_start_adr + MM_SLAB_SIZE);
        compact_num_slabs = 2;
      }
      else {
        compact_start_adr = slab2->address;
        compact_num_slabs = 1;
      }
      mm_slab_destroy_empty(pool, slab2, head, 1);
    }
  }

  // For larger slot sizes we don't keep empty slabs at all.
  else {

    // Compact slab?
    if (slab->state == SLAB_EMPTY) {
      compact_start_adr = slab->address;
      compact_num_slabs = 1;
      mm_slab_destroy_empty(pool, slab, head, !slab_previously_full);
    }

    // Attach any intermediates
    if (num_slabs) {
      if (compact_start_adr) {
        compact_num_slabs += num_slabs;
      }
      else {
        // slab has not been destroyed above; slab->address is valid
        compact_start_adr = slab->address + MM_SLAB_SIZE;
        compact_num_slabs = num_slabs;
      }
    }

    // Compact slab2?
    if ((slab2) && (slab2->state == SLAB_EMPTY)) {
      if (compact_start_adr) {
        // slab2 is adjacent to any intermediate slabs (or slab if no
        // intermediates exist)
        ar_assert(compact_start_adr + compact_num_slabs * MM_SLAB_SIZE ==
                  slab2->address);
        compact_num_slabs++;
      }
      else {
        compact_start_adr = slab2->address;
        compact_num_slabs = 1;
      }
      mm_slab_destroy_empty(pool, slab2, head, !slab2_previously_full);
    }
  }

  // After all this mess, do we have anything left to compact?
  if (!compact_start_adr) {
    context->mm_kernel_pool = pool_save;
    goto success;
  }
  ar_assert(compact_num_slabs > 0);

  // Add slabs to free pool
  mm_slab_add_to_free(pool, compact_start_adr, compact_num_slabs);


success:

  // Kernel pool receives special attention
  if (pool->is_kernel) {
    ar_assert(context->mm_busy_freeing);
    ar_assert(context->mm_defer_frees);

    // Finished with this free. Any pointers that we owe?
    if (context->mm_num_defer_frees) {
      adr = context->mm_defer_frees[--context->mm_num_defer_frees];
      goto begin;
    }
    else {
      context->mm_busy_freeing = 0;
    }

    // Finished with all pending frees. See if alloc_slot has left over
    // preallocations to do -- it may be so, because it doesn't enter
    // preallocation while we're busy freeing stuff.
    ar_assert(!context->mm_recursion_depth);
    if (context->mm_prealloc_needed) {
      context->mm_recursion_depth = 1;
      mm_slab_preallocate();
      context->mm_recursion_depth = 0;
      ar_assert(!context->mm_prealloc_needed);
    }
  }

  // Success
  context->mm_kernel_pool = pool_save;
  return 0;
}


// ===========================================================================
// mm_slab_expand_pool()        Increases slab pool's free pool by a given
//                              chunk of free slabs. Basically a wrapper
//                              to mm_slab_add_to_free() which does the 
//                              job and handles merging with neighbors.
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_slab_expand_pool(MmSlabPool *pool, size_t adr, int num_slabs) {

  // Add slabs to free pool, merging them as needed
  mm_slab_add_to_free(pool, adr, num_slabs);

  // Increase free space
  pool->free_space += num_slabs * MM_SLAB_SIZE;

  // Success
  return 0;
}


// ===========================================================================
// mm_slab_reduce_pool()        Decreases slab pool's free pool by a specific
//                              number of slabs. The function will do the
//                              decrement only if exactly num_slabs can be
//                              removed in a single chunk from the free pool.
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The pool to give back free slabs
//   int num_slabs              Number of slabs to be given back.
//
// * OUTPUTS
//   size_t *ret_adr            Address of first slab to be given back
//
// * RETURN VALUE
//   int                        0: success; ret_adr valid
//                              ERR_OUT_OF_MEMORY: could not give back num_slabs
//                                 (ret_adr invalid) because no num_slabs
//                                 consecutive slabs could be located while
//                                 maintaining the integrity of the used
//                                 pool's memory.
// ===========================================================================
int mm_slab_reduce_pool(MmSlabPool *pool, int num_slabs, size_t *ret_adr) {

  size_t        free_slab_adr;
  MmSlab        *free_slab;


  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->free_slabs);
  ar_assert(ret_adr);

  // Search free pool from the top down. Higher addresses are more probable to
  // be contiguous.
  for (free_slab_adr = kt_trie_find_minmax(pool->free_slabs, 1, 
                                           (void *) &free_slab);
       free_slab_adr;
       free_slab_adr = kt_trie_find_next(pool->free_slabs, 0,
                                         (void *) &free_slab)) {
    if (free_slab->num_slabs >= num_slabs) {
      break;
    }
  }

  // No suitable slabs?
  if (!free_slab_adr) {
    return ERR_OUT_OF_MEMORY;
  }

  // Modify the free slab
  ar_assert(free_slab);
  ar_assert(free_slab->state == SLAB_EMPTY);
  free_slab->num_slabs -= num_slabs;

  // Store return address
  *ret_adr = free_slab_adr + free_slab->num_slabs * MM_SLAB_SIZE;

  // Is it exhausted?  Delete it.
  if (!free_slab->num_slabs) {
    ar_assert(!kt_trie_delete(pool->free_slabs, free_slab_adr, kt_free));
  }

  // Update free space
  ar_assert(pool->free_space >= num_slabs * MM_SLAB_SIZE);
  pool->free_space -= num_slabs * MM_SLAB_SIZE;

  // Success
  return 0;
}


// ===========================================================================
// mm_slab_query_pointer()      Looks up for a previously allocated pointer
//                              in the slab pool and returns its size, if
//                              the pointer exists.
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The slab pool for the new slot
//   size_t *adr                The slot pointer to be queried
//
// * RETURN VALUE
//   size_t                     > 0: success, pointer size
//                              ERR_MISALIGNED: pointer not aligned to 
//                                 MM_ALLOC_ALIGN or to slab slot
//                              ERR_OUT_OF_RANGE: pointer not inside this pool
//                                 used/partials slabs
//                              ERR_NOT_ALLOCED: pointer in range, but mapped 
//                                 to empty slot or free slab
// ===========================================================================
int mm_slab_query_pointer(MmSlabPool *pool, size_t adr) {

  MmSlab                *slab;
  int                   ret;

  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->used_slabs);
  ar_assert(adr);

  // Find the slab this pointer refers to
  ret = mm_slab_locate_pointer(pool, adr, &slab, NULL);
  if (ret) {
    return ret;
  }

  // Return size
  return slab->slot_size;
}


// ===========================================================================
// mm_slab_query_pool()         Creates a list of address ranges that
//                              represent the used slabs of the given pool.
// ===========================================================================
// * INPUTS
//   MmSlabPool *pool           The pool to be queried
//   int pack_options           Packing option flags (MM_PACK_OPTION_*)
//   int default_location       The core ID which currently "owns" this pool,
//                              except for the exception_locations below.
//   Trie *exception_locations  Object slots that are at a different location
//                              than default_location (key: slot adr, data:
//                              core ID). If no exceptions are valid, NULL
//                              can be specified.
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
//                              The MM_PACK_LOCATION_BITS bits will contain
//                              the location of the addess range entry and
//                              MM_PACK_OPTION_BITS will contain the 
//                              packing option flags.
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int mm_slab_query_pool(MmSlabPool *pool, int pack_options, 
                       int default_location, Trie *exception_locations, 
                       int *alloc_num_elements, int *num_elements, 
                       int **sizes, size_t **addresses) {

  size_t        chunk_adr;
  int           chunk_size;
  MmSlab        *slab;
  size_t        slab_adr;
  size_t        cur_exc_adr;
  int           cur_exc_location;
  int           num_slabs;


  // Sanity checks
  ar_assert(pool);
  ar_assert(pool->used_slabs);
  ar_assert(pool->partial_slabs);
  ar_assert(pool->free_slabs);
  ar_assert(alloc_num_elements);
  ar_assert(num_elements);
  ar_assert(sizes);
  ar_assert(addresses);
  ar_assert(pack_options < (1 << MM_PACK_OPTION_BITS));

  // Find first exception address, if it exists
  if (exception_locations) {
    cur_exc_adr = kt_trie_find_minmax(exception_locations, 0, 
                                      (void *) &cur_exc_location);
  }
  else {
    cur_exc_adr = 0;
  }

  // Walk all used slabs nodes in ascending order
  chunk_adr = 0;
  chunk_size = 0;
  for (slab_adr = kt_trie_find_minmax(pool->used_slabs, 0, (void *) &slab);
       slab_adr;
       slab_adr = kt_trie_find_next(pool->used_slabs, 1, (void *) &slab)) {
    
    // Ignore preallocated empty slabs
    if (slab->state == SLAB_EMPTY) {
      continue;
    }

    // How many slabs does this one represent?  If slab->num_slabs is 1, it
    // means the last slab slot does not extend into any other slab (so we
    // count it for 1 here). If it's > 1, it means it uses slab->num_slabs - 1
    // adjacent slabs.
    //
    // The question is: are we going to find the last one as a new used_slabs
    // entry? If the last slot is aligned to end on the end of the last
    // adjacent slab, we won't: we have to count it here.  Otherwise, we'll
    // count it when we find the next entry.
    ar_assert(slab->num_slabs > 0);
    ar_assert(slab_adr == slab->address);
    if (slab->num_slabs == 1) {
      num_slabs = 1;
    }
    else if ((slab->offset + slab->bitmap_size * slab->slot_size) % 
             MM_SLAB_SIZE) {
      num_slabs = slab->num_slabs - 1;
    }
    else {
      num_slabs = slab->num_slabs;
    }

    // Try to create a chunk of as many slabs as possible
    if (!chunk_adr) {
      chunk_adr = slab_adr;
      chunk_size = num_slabs * MM_SLAB_SIZE;
    }
    else {
      // slab is a continuation of existing chunk? Add it to the current
      // chunk and continue with the next slab.
      if (slab_adr == chunk_adr + chunk_size) {
        chunk_size += num_slabs * MM_SLAB_SIZE;
      }

      // Existing chunk ends here; store it.
      else {
//kt_printf("Store 0x%08X size 0x%08X loc %d [exc: 0x%08X loc %d]\r\n", chunk_adr, chunk_size, default_location, cur_exc_adr, cur_exc_location);
        mm_slab_store_chunk(pool, pack_options, default_location, 
                            exception_locations, alloc_num_elements, 
                            num_elements, sizes, addresses, 
                            &cur_exc_adr, &cur_exc_location, &chunk_adr, 
                            &chunk_size);

        // Begin new chunk with this slab
        chunk_adr = slab_adr;
        chunk_size = num_slabs * MM_SLAB_SIZE;
      }
    }
  }

  // Any leftover?
  if (chunk_adr) {

    // Store it
//kt_printf("* Store 0x%08X size 0x%08X loc %d [exc: 0x%08X loc %d]\r\n", chunk_adr, chunk_size, default_location, cur_exc_adr, cur_exc_location);
    mm_slab_store_chunk(pool, pack_options, default_location, 
                        exception_locations, alloc_num_elements, 
                        num_elements, sizes, addresses, 
                        &cur_exc_adr, &cur_exc_location, &chunk_adr, 
                        &chunk_size);
  }

  // There must be no exceptions left: every object should be part of a slab,
  // which we must have processed by now.
  ar_assert(!cur_exc_adr);

  // Success
  return 0;
}


// ###########################################################################
// ###                                                                     ###
// ###                         DEBUGGING FUNCTIONS                         ###
// ###                                                                     ###
// ###########################################################################
#if 0
// ===========================================================================
// mm_slab_dbg_slab_print()     Debug function to print a slab structure
// ===========================================================================
// * INPUTS
//   MmSlab *s                  The slab
// ===========================================================================
void mm_slab_dbg_slab_print(MmSlab *s) {
  int i;

  if (s->state == SLAB_EMPTY) {
    printf("    num_slabs:   %d\n", s->num_slabs);
    printf("    slot_size:   %d\n", s->slot_size);
  }
  else if ((s->state == SLAB_PARTIAL) || (s->state == SLAB_FULL)) {
    printf("    bitmap:      [");
    for (i = 0; i < s->bitmap_size; i += 32) {
      printf("0x%08X ", s->bitmap[i/32]);
    }
    printf("\b]\n");
    printf("    bitmap_size: %d\n", s->bitmap_size);
    printf("    trailer_bit: %d\n", s->trailer_bit);
    printf("    free_slots:  %d\n", s->free_slots);
    printf("    slot_size:   %d\n", s->slot_size);
    printf("    offset:      %d\n", s->offset);
    printf("    num_slabs:   %d\n", s->num_slabs);
  }
  else {
    printf("    Invalid state!\n");
    ar_abort();
  }
}


// ===========================================================================
// mm_slab_dbg_pool_print()     Debug function to print a slab pool structure
// ===========================================================================
// * INPUTS
//   MmSlabPool *p              The slab pool
// ===========================================================================
void mm_slab_dbg_pool_print(MmSlabPool *p) {
  size_t key;
  MmSlab *data;
  MmSlabPartialHead *head;

  printf("\n==============================================================\n");
  printf("Pool %p, metadata pool %p\n", p, p->metadata_pool);
  printf("Used %lu, free %lu\n", p->used_space, p->free_space);
  printf("==============================================================\n");

  printf("Used slabs:\n");
  for (key = kt_trie_find_minmax(p->used_slabs, 0, (void **) &data); 
       key; 
       key = kt_trie_find_next(p->used_slabs, 1, (void **) &data)) {
    ar_assert(key == data->address);
    printf("  Slab 0x%lx, state %s:\n", key, 
        (data->state == SLAB_EMPTY)   ? "EMPTY" :
        (data->state == SLAB_PARTIAL) ? "PARTIAL" :
        (data->state == SLAB_FULL)    ? "FULL" : 
                                        "XXXXXXX ERROR! XXXXXXX");
    mm_slab_dbg_slab_print(data);
  }
  printf("\n");

  printf("Partial slabs:\n");
  for (key = kt_trie_find_minmax(p->partial_slabs, 0, (void **) &head); 
       key; 
       key = kt_trie_find_next(p->partial_slabs, 1, (void **) &head)) {
    printf(
        "  Slot size %lu: Members %d (%d empty), Slots: used %lu, free %lu\n", 
        key, head->num_members, head->num_empty_members, head->used_slots,
        head->free_slots);
    ar_assert(head->head);
    for (data = head->head; data; data = data->next) {
      ar_assert(data->slot_size == key);
      printf("    Slab 0x%lx (%s)\n", data->address,
        (data->state == SLAB_EMPTY)   ? "EMPTY" :
        (data->state == SLAB_PARTIAL) ? "PARTIAL" : 
                                        "XXXXXXX ERROR! XXXXXXX");
    }
  }
  printf("\n");

  printf("Free slabs:\n");
  for (key = kt_trie_find_minmax(p->free_slabs, 0, (void **) &data); 
       key; 
       key = kt_trie_find_next(p->free_slabs, 1, (void **) &data)) {
    ar_assert(key == data->address);
    printf("  Slab 0x%lx, state %s:\n", key, 
        (data->state == SLAB_EMPTY) ? "EMPTY" :
        (data->state == SLAB_PARTIAL) ? "PARTIAL" :
        (data->state == SLAB_FULL) ? "FULL" : "XXXXXXXXX ERROR! XXXXXXXXX");
    mm_slab_dbg_slab_print(data);
  }
  printf("\n");
  
  printf("==============================================================\n\n");
}


// ===========================================================================
// mm_slab_dbg_pool_dotty()     Debug function to create a DOT graph file for 
//                              a slab pool
// ===========================================================================
// * INPUTS
//   MmSlabPool *p              The slab pool
//   char *filename             Filename for the graph output
// ===========================================================================
void mm_slab_dbg_pool_dotty(MmSlabPool *p, char *filename) {
  size_t key;
  MmSlab *data;
  MmSlabPartialHead *head;
  FILE *f;
  char buf[128], prev_buf[128];
  int max;
  int i;
  int header;


  // Find max slot size
  max = 0;
  for (key = kt_trie_find_minmax(p->used_slabs, 0, (void **) &data); 
       key; 
       key = kt_trie_find_next(p->used_slabs, 1, (void **) &data)) {
    if (data->slot_size > max) {
      max = data->slot_size;
    }
  }


  // Open file
  ar_assert(f = fopen(filename, "w"));

  // Print out header
  fprintf(f, 
"digraph \"Pool %p, metadata pool %p\" {\n\n", 
          p, p->metadata_pool);
  fprintf(f, 
"  rankdir = LR;\n"
  );
  

  // Create a per-slot subgraph. This is slow, but it doesn't alloc any memory
  // (which could screw with the kernel pool, if we're dumping it) and it's a
  // debugging mode anyway...
  for (i = max; i >= 64; i -= 64) {
  
    prev_buf[0] = 0;
    header = 0;
    for (key = kt_trie_find_minmax(p->used_slabs, 0, (void **) &data); 
         key; 
         key = kt_trie_find_next(p->used_slabs, 1, (void **) &data)) {
      ar_assert(key == data->address);
      if (data->slot_size != i) {
        continue;
      }
      if (!header) {
        header = 1;
        fprintf(f, 
"\n  subgraph cluster_%d { label=\"Slot size %d\"; \n", 
                i, i);
      }

      if (data->state == SLAB_EMPTY) {
        fprintf(f,
"    slab0x%lx [shape=record, style=filled, fillcolor=magenta, label=\"adr: 0x%lx | state: Empty | num_slabs: %d\"];\n", 
                key, key, data->num_slabs);
      }
      else if (data->state == SLAB_PARTIAL) {
        fprintf(f,
"    slab0x%lx [shape=record, style=filled, fillcolor=green, label=\"adr: 0x%lx | state: Partial | slot_size: %d | free_slots: %d/%d | offset: %d | trailer_bit: %d | num_slabs: %d\"];\n", 
                key, key, data->slot_size, data->free_slots, data->bitmap_size,
                data->offset, data->trailer_bit, data->num_slabs);
      }
      else if (data->state == SLAB_FULL) {
        fprintf(f,
"    slab0x%lx [shape=record, style=filled, fillcolor=royalblue, label=\"adr: 0x%lx | state: Full | slot_size: %d | free_slots: %d/%d | offset: %d | trailer_bit: %d | num_slabs: %d\"];\n", 
                key, key, data->slot_size, data->free_slots, data->bitmap_size,
                data->offset, data->trailer_bit, data->num_slabs);
      }
      else {
        ar_abort();
      }
      sprintf(buf, "slab0x%lx", key);
      if (strlen(prev_buf)) {
        fprintf(f,
"  %s -> %s [style=invis, weight=42]\n", prev_buf, buf);
      }
      strcpy(prev_buf, buf);
    }
    if (header) {
      fprintf(f,
"  }\n\n");
    }

  }

  // Do the free pool subgraph
  fprintf(f, "\n  subgraph \"cluster_free\" { label=\"Free pool\";\n");
  prev_buf[0] = 0;
  for (key = kt_trie_find_minmax(p->free_slabs, 0, (void **) &data); 
       key; 
       key = kt_trie_find_next(p->free_slabs, 1, (void **) &data)) {
    ar_assert(key == data->address);
    fprintf(f,
"    slab0x%lx [shape=record, style=filled, fillcolor=red, label=\"adr: 0x%lx | state: Empty | num_slabs: %d\"];\n", 
            key, key, data->num_slabs);
    sprintf(buf, "slab0x%lx", key);
    if (strlen(prev_buf)) {
      fprintf(f,
"  %s -> %s [style=invis, weight=42]\n", prev_buf, buf);
    }
    strcpy(prev_buf, buf);
  }
  fprintf(f, "  }\n\n");
  

  // Link the partial nodes
  for (key = kt_trie_find_minmax(p->partial_slabs, 0, (void **) &head); 
       key; 
       key = kt_trie_find_next(p->partial_slabs, 1, (void **) &head)) {
    ar_assert(head->head);
    prev_buf[0] = 0;
    for (data = head->head; data; data = data->next) {
      ar_assert(data->slot_size == key);
      sprintf(buf, "slab0x%lx", data->address);
      if (strlen(prev_buf)) {
        fprintf(f,
"  %s -> %s;\n",
                prev_buf, buf);
      }
      strcpy(prev_buf, buf);
    }
  }

  fprintf(f, "\n}\n");
  fflush(f);
  fclose(f);
}
#endif
