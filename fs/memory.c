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
// Abstract      : AntFS utility functions that manage the in-memory-only
//                 structures of the filesystem
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: memory.c,v $
// CVS revision  : $Revision: 1.9 $
// Last modified : $Date: 2012/03/15 15:31:47 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <types.h>
#include <kernel_toolset.h>
#include <memory_management.h>
#include <arch.h>
#include <filesystem.h>
#include <errno.h>


// ===========================================================================
// fs_mem_set_free_block()      Updates the memory state to specify that
//                              the given block is free. The block is freed
//                              only in the future_free_blocks structure;
//                              free_blocks remains consistent with the last
//                              checkpoint's view of free blocks, which is
//                              used during log replay.
// ===========================================================================
// * INPUTS
//   unsigned int blk_num       Block number
//
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_set_free_block(Context *context, unsigned int blk_num) {
  unsigned int index;
  unsigned int entry_index;
  unsigned char mask;
  FSCheckpointFreeBitmap *bitmap, *future_bitmap;

  // Sanity check
  if (blk_num >= context->fs_num_blocks) {
    return -EINVAL;
  }

  // Find blk_num position in checkpoint free bitmaps
  ar_uint_divide(blk_num / 8, FS_CHK_FREE_BITMAP_BYTES, &index, &entry_index);
  if (index >= context->fs_state->num_chk_free) {
    ar_panic(
        "AntFS: set_free_block needs more checkpoint free blocks than normal");
  }
  bitmap = context->fs_state->chk_free[index];
  future_bitmap = context->fs_state->chk_future_free[index];
  mask = 1 << (blk_num % 8);

  // Check if it's already free
  if ((bitmap->bitmap[entry_index] & mask) ||
      (future_bitmap->bitmap[entry_index] & mask)) {
    return -EINVAL;
  }

  // Make the block free (future version only)
  future_bitmap->bitmap[entry_index] |= mask;
  context->fs_state->num_future_free_blocks++;

  return 0;
}


// ===========================================================================
// fs_mem_unset_free_block()    Updates the memory state to specify that
//                              the given block is not free anymore. Both
//                              free_blocks and future_free_blocks are 
//                              modified.
// ===========================================================================
// * INPUTS
//   unsigned int blk_num       Block number
//
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_unset_free_block(Context *context, unsigned int blk_num) {
  unsigned int index;
  unsigned int entry_index;
  unsigned char mask;
  FSCheckpointFreeBitmap *bitmap, *future_bitmap;

  // Sanity check
  if (blk_num >= context->fs_num_blocks) {
    return -EINVAL;
  }

  // Find blk_num position in checkpoint free bitmaps
  ar_uint_divide(blk_num / 8, FS_CHK_FREE_BITMAP_BYTES, &index, &entry_index);
  if (index >= context->fs_state->num_chk_free) {
    ar_panic(
      "AntFS: unset_free_block needs more checkpoint free blocks than normal");
  }
  bitmap = context->fs_state->chk_free[index];
  future_bitmap = context->fs_state->chk_future_free[index];
  mask = 1 << (blk_num % 8);


  // Check if it's already non-free
  if ((!(bitmap->bitmap[entry_index] & mask)) ||
      (!(future_bitmap->bitmap[entry_index] & mask))) {
    return -EINVAL;
  }

  // Make the block non-free (both current and future versions)
  bitmap->bitmap[entry_index] &= ~mask;
  context->fs_state->num_free_blocks--;
  future_bitmap->bitmap[entry_index] &= ~mask;
  context->fs_state->num_future_free_blocks--;

  return 0;
}


// ===========================================================================
// fs_mem_free_block_status()   Returns either the conservative or the future
//                              view of whether a given block number is free 
//                              or not.
// ===========================================================================
// * INPUTS
//   Context *context           Current context
//   unsigned int blk_num       Block number
//   int future                 1: the future view is returned
//                              0: the conservative view is returned
//
// * OUTPUTS
//   unsigned int *ret_status   1 if the block is free, 0 if it's not
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_free_block_status(Context *context, unsigned int blk_num,
                             int future, int *ret_status) {

  unsigned int index;
  unsigned int entry_index;
  unsigned char mask;
  FSCheckpointFreeBitmap *bitmap;

  // Sanity check
  if (blk_num >= context->fs_num_blocks) {
    return -EINVAL;
  }
  if (!ret_status) {
    ar_panic("AntFS: free_block_status called with no return");
  }

  // Find blk_num position in checkpoint free bitmaps
  ar_uint_divide(blk_num / 8, FS_CHK_FREE_BITMAP_BYTES, &index, &entry_index);
  if (index >= context->fs_state->num_chk_free) {
    ar_panic(
      "AntFS: free_block_status needs more checkpoint free blocks than normal");
  }
  if (future) {
    bitmap = context->fs_state->chk_future_free[index];
  }
  else {
    bitmap = context->fs_state->chk_free[index];
  }
  mask = 1 << (blk_num % 8);

  // Return the status
  if (bitmap->bitmap[entry_index] & mask) {
    *ret_status = 1;
  }
  else {
    *ret_status = 0;
  }

  return 0;
}


// ===========================================================================
// fs_mem_find_free_block()     Returns the next available free block, 
//                              starting from the current log tail position.
//                              We use the conservative version of free_blocks,
//                              which is consistent with the last checkpoint's
//                              view, used during log replay.
// ===========================================================================
// * OUTPUTS
//   unsigned int *ret_blk_num  Returned block number of the next free block
//
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_find_free_block(Context *context, unsigned int *ret_blk_num) {
  
  FSCurrentState *state = context->fs_state;
  unsigned int blk_num;
  int status;


  // Is there a free next block at all?
  if (!state->num_free_blocks) {
    return -ENOSPC;
  }

  // Start from the tail of the log
  blk_num = state->log_tail;
  while (1) {
    
    // Is the block (conservatively) free?
    if (fs_mem_free_block_status(context, blk_num, 0, &status)) {
      ar_panic("AntFS: find_free_block could not read free block status");
    }
    
    if (status) {
      // We found it
      *ret_blk_num = blk_num;
      return 0;
    }

    // Advance block number and wrap around. Remember that the very first and
    // last blocks are always checkpoints and thus never free.
    if (blk_num >= context->fs_num_blocks - 2) {
      blk_num = 1;
    }
    else {
      blk_num++;
    }

    // If we wrapped around to reach the log tail again, it's a consistency
    // problem: state->num_free_blocks should have told us in the first 
    // place that we had no free space. Abort.
    if (blk_num == state->log_tail) {
      ar_panic("AntFS: No free block in bitmap, inconsistent num_free_blocks");
    }
  }
 
  ar_panic("AntFS: fs_mem_find_free_block loop inconsistency");

  // Never reached
  return 0;
}


// ===========================================================================
// fs_mem_sync_free_blocks()    Copies the real free blocks bitmap and count
//                              onto the conservative view.
// ===========================================================================
// * INOUTS
//   Context *context           Current context
// ===========================================================================
void fs_mem_sync_free_blocks(Context *context) {
  
  FSCheckpointFreeBitmap *bitmap, *future_bitmap;
  int i, j;


  // For all blocks
  for (i = 0; i < context->fs_state->num_chk_free; i++) {

    bitmap = context->fs_state->chk_free[i];
    future_bitmap = context->fs_state->chk_future_free[i];
    
    // Synchronize the bitmaps
    for (j = 0; j < FS_CHK_FREE_BITMAP_BYTES; j++) {
      bitmap->bitmap[j] = future_bitmap->bitmap[j];
    }
  }


  // Synchronize the counts
  context->fs_state->num_free_blocks = 
    context->fs_state->num_future_free_blocks;
}


// ===========================================================================
// fs_mem_init_free_blocks()    Initializes both the conservative and the
//                              future views of the free block bitmaps
//                              to "free"
// ===========================================================================
// * INOUTS
//   Context *context           Current context
// ===========================================================================
void fs_mem_init_free_blocks(Context *context) {
  int i;

  // For all blocks
  for (i = 0; i < context->fs_state->num_chk_free; i++) {

    // Initialize the bitmaps
    kt_memset(context->fs_state->chk_free[i]->bitmap, 
              0xFF, FS_CHK_FREE_BITMAP_BYTES);
    kt_memset(context->fs_state->chk_future_free[i]->bitmap, 
              0xFF, FS_CHK_FREE_BITMAP_BYTES);
  }

  // Initialize free block counters
  context->fs_state->num_free_blocks = context->fs_num_blocks;
  context->fs_state->num_future_free_blocks = context->fs_num_blocks;
}


// ===========================================================================
// fs_mem_set_data_block_ptr()  Updates a data block structure (i.e. a 
//                              collection of data block pointers of 3 bytes)
//                              with a given pointer
// ===========================================================================
// * INPUTS
//   Context *context           Current context
//   int index                  Index of data_blocks to update
//   unsigned int new_ptr       New block pointer to write
//
// * INOUTS
//   void *block                An Inode, Indirect or Checkpoint block whose
//                              data_blocks structure will be updated.
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_set_data_block_ptr(Context *context, void *block, int index, 
                          unsigned int new_ptr) {
  unsigned int header;
  unsigned char *array;

  // Sanity check based on block header
  if (!block) {
    return -EINVAL;
  }
  kt_memcpy(&header, block, sizeof(unsigned int));
  if ((index < 0) ||                           // Invalid index
      (new_ptr >= context->fs_num_blocks) ||   // New ptr too big
      (header & FS_HDR_TYPE_DATA) ||           // Data blocks have no ptrs
      (header & FS_HDR_TYPE_CHK_FREE_BITMAP) || // Free bitmaps have no ptrs
      !((header & FS_HDR_TYPE_INODE) || 
        (header & FS_HDR_TYPE_INDIRECT) ||
        (header & FS_HDR_TYPE_DELETE) ||
        (header & FS_HDR_TYPE_CHK_INDIRECT) ||
        (header & FS_HDR_TYPE_CHK_ROOT)) ||    // Header has no type
      ((header & FS_HDR_TYPE_INODE) && 
       (index >= FS_INODE_NUM_DATA_BLOCKS)) || // Inode #blocks limit error
      (((header & FS_HDR_TYPE_INDIRECT) || 
        (header & FS_HDR_TYPE_CHK_INDIRECT)) && 
       (index >= FS_IND_NUM_DATA_BLOCKS)) ||   // Indirect #blocks limit error
      ((header & FS_HDR_TYPE_DELETE) && 
       (index >= FS_DELETE_NUM_BLOCKS)) ||     // Delete #blocks limit error
      ((header & FS_HDR_TYPE_CHK_ROOT) && 
       (index >= FS_CHK_NUM_DATA_BLOCKS))      // Chkpoint #blocks limit error
     ) {
    return -EINVAL;
  }

  // Find data blocks array. We checked above that header is one of the
  // three expected types.
  array = NULL;
  if (header & FS_HDR_TYPE_INODE) {
    array = ((FSInode *) block)->data_blocks;
  }
  if ((header & FS_HDR_TYPE_INDIRECT) || (header & FS_HDR_TYPE_CHK_INDIRECT)) {
    array = ((FSIndirect *) block)->data_blocks;
  }
  if (header & FS_HDR_TYPE_CHK_ROOT) {
    array = ((FSCheckpointRoot *) block)->data_blocks;
  }
  if (header & FS_HDR_TYPE_DELETE) {
    array = ((FSDelete *) block)->free_blocks;
  }

  // Write new pointer
  array[3 * index    ] = (new_ptr >> 16) & 0xFF;
  array[3 * index + 1] = (new_ptr >>  8) & 0xFF;
  array[3 * index + 2] =  new_ptr        & 0xFF;

  return 0;
}


// ===========================================================================
// fs_mem_get_data_block_ptr()  Reads a data block pointer from a data 
//                              block structure (i.e. a collection of data 
//                              block pointers of 3 bytes)
// ===========================================================================
// * INPUTS
//   Context *context           Current context
//   int index                  Index of data_blocks to update
//
// * INOUTS
//   void *block                An Inode, Indirect or Checkpoint block whose
//                              data_blocks structure will be updated.
// * OUTPUTS
//   unsigned int *ret_ptr      Returned pointer read from the data_blocks
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_get_data_block_ptr(Context *context, void *block, int index, 
                              unsigned int *ret_ptr) {
  unsigned int header;
  unsigned char *array;

  // Sanity check based on block header
  if (!block) {
    return -EINVAL;
  }
  kt_memcpy(&header, block, sizeof(unsigned int));
  if ((index < 0) ||                           // Invalid index
      (!ret_ptr) ||                            // No pointer to return
      (header & FS_HDR_TYPE_DATA) ||           // Data blocks have no ptrs
      (header & FS_HDR_TYPE_CHK_FREE_BITMAP) || // Free bitmaps have no ptrs
      !((header & FS_HDR_TYPE_INODE) || 
        (header & FS_HDR_TYPE_INDIRECT) ||
        (header & FS_HDR_TYPE_DELETE) ||
        (header & FS_HDR_TYPE_CHK_INDIRECT) ||
        (header & FS_HDR_TYPE_CHK_ROOT)) ||    // Header has no type
      ((header & FS_HDR_TYPE_INODE) && 
       (index >= FS_INODE_NUM_DATA_BLOCKS)) || // Inode #blocks limit error
      (((header & FS_HDR_TYPE_INDIRECT) || 
        (header & FS_HDR_TYPE_CHK_INDIRECT)) && 
       (index >= FS_IND_NUM_DATA_BLOCKS)) ||   // Indirect #blocks limit error
      ((header & FS_HDR_TYPE_DELETE) && 
       (index >= FS_DELETE_NUM_BLOCKS)) ||     // Delete #blocks limit error
      ((header & FS_HDR_TYPE_CHK_ROOT) && 
       (index >= FS_CHK_NUM_DATA_BLOCKS))      // Chkpoint #blocks limit error
     ) {
    return -EINVAL;
  }

  // Find data blocks array. We checked above that header is one of the
  // three expected types.
  array = NULL;
  if (header & FS_HDR_TYPE_INODE) {
    array = ((FSInode *) block)->data_blocks;
  }
  if ((header & FS_HDR_TYPE_INDIRECT) || (header & FS_HDR_TYPE_CHK_INDIRECT)) {
    array = ((FSIndirect *) block)->data_blocks;
  }
  if (header & FS_HDR_TYPE_CHK_ROOT) {
    array = ((FSCheckpointRoot *) block)->data_blocks;
  }
  if (header & FS_HDR_TYPE_DELETE) {
    array = ((FSDelete *) block)->free_blocks;
  }

  // Read the pointer
  *ret_ptr = 0;
  *ret_ptr |= array[3 * index    ] << 16;
  *ret_ptr |= array[3 * index + 1] << 8;
  *ret_ptr |= array[3 * index + 2];

  return 0;
}


// ===========================================================================
// fs_mem_read_inode_map()      Reads the inode map and returns the given
//                              entry
// ===========================================================================
// * INPUTS
//   unsigned int inode_num     The inode number that is being updated
//   Context *context           The current context, whose inode map
//                              (live checkpoint block or indirect) is 
//                              updated
// * OUTPUTS 
//   unsigned int               The returned inode block number
//      *ret_inode_blk_num 
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int fs_mem_read_inode_map(Context *context, unsigned int inode_num,
                          unsigned int *ret_inode_blk_num) {
  unsigned char *block;
  int ind_index;
  int entry_index;
  
  // Sanity check
  if (inode_num >= FS_MAX_INODES) {
    ar_panic("AntFS: fs_mem_read_inode_map called with invalid inode_num");
  }
  if (!ret_inode_blk_num) {
    ar_panic("AntFS: fs_mem_read_inode_map called with no return ptr");
  }
  
  // See which block is responsible for inode_num
  if (inode_num < FS_CHK_NUM_DATA_BLOCKS) {
    block = (unsigned char *) context->fs_state->chk;
    entry_index = inode_num;
  }
  else {
    ar_int_divide(inode_num - FS_CHK_NUM_DATA_BLOCKS, 
                  FS_IND_NUM_DATA_BLOCKS, 
                  &ind_index, &entry_index);
    if (ind_index >= context->fs_state->num_chk_ind) {
      ar_panic(
          "AntFS: fs_mem_read_inode_map needs more indirects than normal");
    }
    block = (unsigned char *) context->fs_state->chk_ind[ind_index];
  }

  // Read inode map
  if (fs_mem_get_data_block_ptr(context, block, entry_index, 
                                ret_inode_blk_num)) {
    ar_panic("AntFS: fs_mem_read_inode_map: Could not get inode in inode map");
  }

  return 0;
}


// ===========================================================================
// fs_mem_update_inode_map()    Links an inode with a block number in the
//                              inode map
// ===========================================================================
// * INPUTS
//   unsigned int inode_num     The inode number that is being updated
//   unsigned int inode_blk_num The new block number that inode_num is
//                              associated to
//
// * INOUTS 
//   Context *context           The current context, whose inode map
//                              (live checkpoint block or indirect) is 
//                              updated
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int fs_mem_update_inode_map(Context *context, unsigned int inode_num,
                            unsigned int inode_blk_num) {
  unsigned char *block;
  int ind_index;
  int entry_index;
  
  // Sanity check
  if (inode_num >= FS_MAX_INODES) {
    ar_panic("AntFS: fs_mem_update_inode_map called with invalid inode_num");
  }
  
  // See which block is responsible for inode_num
  if (inode_num < FS_CHK_NUM_DATA_BLOCKS) {
    block = (unsigned char *) context->fs_state->chk;
    entry_index = inode_num;
  }
  else {
    ar_int_divide(inode_num - FS_CHK_NUM_DATA_BLOCKS, 
                  FS_IND_NUM_DATA_BLOCKS, 
                  &ind_index, &entry_index);
    if (ind_index >= context->fs_state->num_chk_ind) {
      ar_panic(
          "AntFS: fs_mem_update_inode_map needs more indirects than normal");
    }
    block = (unsigned char *) context->fs_state->chk_ind[ind_index];
  }

  // Update inode map
  if (fs_mem_set_data_block_ptr(context, block, entry_index, inode_blk_num)) {
    ar_panic("AntFS: create_inode: Could not set inode in inode map");
  }

  return 0;
}


// ===========================================================================
// fs_mem_find_free_inode()     Returns the next available free inode, 
//                              starting always from 0
// ===========================================================================
// * OUTPUTS
//   unsigned int *ret_inode_num Returned inode number of the next free inode
//
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_find_free_inode(Context *context, unsigned int *ret_inode_num) {
  unsigned int i;
  unsigned int blk_num;


  // Sanity check
  if (!ret_inode_num) {
    ar_panic("AntFS: fs_mem_find_free_inode called with no return");
  }

  // Start from scratch
  for (i = 0; i < FS_MAX_INODES; i++) {
    
    // Check this inode
    if (fs_mem_read_inode_map(context, i, &blk_num)) {
      ar_panic("AntFS: fs_mem_find_free_inode could not read inode map");
    }

    // Did we find it?
    if (!blk_num) {
      *ret_inode_num = i;
      return 0;
    }
  }
  
  // Out of inodes
  return -ENOSPC;
}


// ===========================================================================
// fs_mem_init()                Initializes the memory state that we need
//                              to hold in the context
// ===========================================================================
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_mem_init(Context *context) {
  FSCurrentState *state;
  int i;

  state = context->fs_state;

  // Decide how many indirect blocks we'll need for the inode map.
  // The checkpoint block iself holds FS_CHK_NUM_DATA_BLOCKS inode pointers, 
  // each indirect block holds FS_IND_NUM_DATA_BLOCKS.
  ar_int_divide(FS_MAX_INODES - FS_CHK_NUM_DATA_BLOCKS, 
                FS_IND_NUM_DATA_BLOCKS, 
                &(state->num_chk_ind), &i);
  if (i) {
    state->num_chk_ind++;
  }

  // Allocate live checkpoint and its indirect blocks
  if (!(state->chk = kt_zalloc(sizeof(FSCheckpointRoot)))) {
    goto out_of_memory;
  }
  if (!(state->chk_ind = 
        kt_zalloc(state->num_chk_ind * sizeof(FSIndirect *)))) {
    goto out_of_memory;
  }
  for (i = 0; i < state->num_chk_ind; i++) {
    if (!(state->chk_ind[i] = kt_zalloc(sizeof(FSIndirect)))) {
      goto out_of_memory;
    }
  }

  // Decide how many free bitmap blocks we'll need. Allocate twice that 
  // amount, one for the conservative and one for the real estimation.
  ar_int_divide((context->fs_num_blocks + 7) / 8,
                FS_CHK_FREE_BITMAP_BYTES, 
                &(state->num_chk_free), &i);
  if (i) {
    state->num_chk_free++;
  }
  
  // Allocate free blocks bitmaps
  if (!(state->chk_free = 
        kt_zalloc(state->num_chk_free * sizeof(FSCheckpointFreeBitmap *)))) {
    goto out_of_memory;
  }
  if (!(state->chk_future_free = 
        kt_zalloc(state->num_chk_free * sizeof(FSCheckpointFreeBitmap *)))) {
    goto out_of_memory;
  }
  for (i = 0; i < state->num_chk_free; i++) {
    if (!(state->chk_free[i] = kt_zalloc(sizeof(FSCheckpointFreeBitmap)))) {
      goto out_of_memory;
    }
    if (!(state->chk_future_free[i] = kt_zalloc(
                                        sizeof(FSCheckpointFreeBitmap)))) {
      goto out_of_memory;
    }
  }

  // Some info
  kt_printf("Filesystem has %d blocks, %d max inodes, \r\n"
            "               %d chk indirects, %d chk free blocks\r\n",
        context->fs_num_blocks, FS_MAX_INODES,
        context->fs_state->num_chk_ind, context->fs_state->num_chk_free);

  // Success
  return 0;

out_of_memory:  
  return -ENOMEM;
}

