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
// Abstract      : AntFS low-level functions that deal with the physical
//                 on-disk log
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: log.c,v $
// CVS revision  : $Revision: 1.7 $
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
// fs_log_write_block()         Writes a 4,088 bytes block to the disk, 
//                              appending 8 bytes of its CRC-64 to complete
//                              a 4,096 disk block.
// ===========================================================================
// * INOUTS
//   unsigned char *buf         4,096 bytes buffer. First 4,088 bytes contain
//                              data, last 8 bytes will be filled here with
//                              their CRC-64.
// * INPUTS
//   unsigned int adr           Block address (0, 1, 2, ...) 
//                              
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_log_write_block(unsigned char *buf, unsigned int adr) {
  unsigned int crc_hi, crc_lo;
  int i;
  Context *context;

  // Sanity check
  context = mm_get_context(ar_get_core_id());
  if (adr >= context->fs_num_blocks) {
    return -EINVAL;
  }
  if (buf == NULL) {
    return -EINVAL;
  }

  // Compute CRC-64
  crc_hi = crc_lo = 0;
  for (i = 0; i < 4088; i++) {
    fs_crc64(buf[i], &crc_hi, &crc_lo);
  }
  kt_memcpy(buf + 4088, &crc_hi, 4);
  kt_memcpy(buf + 4092, &crc_lo, 4);

  // Write the whole block
  if (ar_flash_write_sectors(adr * 8, 8, buf)) {
    return -EIO;
  }

  return 0;
}


// ===========================================================================
// fs_log_read_block()          Reads a 4,096 bytes block from the disk and 
//                              checks that its CRC-64 matches its 4,088 bytes
//                              of data matches it.
// ===========================================================================
// * OUTPUTS
//   unsigned char *buf         4,096 bytes buffer. First 4,088 bytes will 
//                              be filled with the daat and last 8 with the
//                              CRC-64 value.
// * INPUTS
//   unsigned int adr           Block address (0, 1, 2, ...) 
//                              
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_log_read_block(unsigned char *buf, unsigned int adr) {
  unsigned int crc_hi, crc_lo;
  int i;
  Context *context;

  // Sanity check
  context = mm_get_context(ar_get_core_id());
  if (adr >= context->fs_num_blocks) {
    ar_panic("AntFS: fs_log_read_block called with invalid block number");
  }
  if (buf == NULL) {
    ar_panic("AntFS: fs_log_read_block called with no return buffer");
  }

  // Read the block
  if (ar_flash_read_sectors(adr * 8, 8, buf)) {
    return -EIO;
  }

  // Compute data CRC-64
  crc_hi = crc_lo = 0;
  for (i = 0; i < 4088; i++) {
    fs_crc64(buf[i], &crc_hi, &crc_lo);
  }

  // Compare to stored CRC-64
  if (kt_memcmp(&crc_hi, buf + 4088, 4)) {
    return -ECRC;
  }
  if (kt_memcmp(&crc_lo, buf + 4092, 4)) {
    return -ECRC;
  }

  return 0;
}


// ===========================================================================
// fs_log_append()              Appends a block to the log
// ===========================================================================
// * INPUTS
//   unsigned char *blk         Block to be written
//
// * INOUTS
//   Context *context           The current context
//
// * OUTPUTS
//   unsigned int *ret_blk_num  Returned block number
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_log_append(Context *context, unsigned char *blk, 
                  unsigned int *ret_blk_num) {
  unsigned int blk_num;

  // Sanity check
  if (!ret_blk_num) {
    ar_panic("AntFS: log_append called with no return element");
  }

  // Find next free block and take it
  if (fs_mem_find_free_block(context, &blk_num)) {
    return -ENOSPC;
  }
  if (fs_mem_unset_free_block(context, blk_num)) {
    ar_panic("AntFS: log_append free block inconsistency");
  }

  // Tag the block with the incrementing log ID
  kt_memcpy(blk + sizeof(unsigned int), 
            &(context->fs_state->log_seq_id), sizeof(unsigned int));
  context->fs_state->log_seq_id++;

  // Write the block to disk
  if (fs_log_write_block(blk, blk_num)) {
    return -EIO;
  }

  // Update log tail. 
  //
  // Log tail is always the next candidate of free space. It equals the
  // last written block + 1, wrapped around, omitting the very first and last
  // blocks that are used for the checkpoints.
  if (blk_num >= context->fs_num_blocks - 2) {
    context->fs_state->log_tail = 1;
  }
  else {
    context->fs_state->log_tail = blk_num + 1;
  }
  
  // Return location of written block
  *ret_blk_num = blk_num;

  return 0;
}


// ===========================================================================
// fs_log_checkpoint()          Writes a new checkpoint to disk, either in 
//                              the first or the last block of the device
//                              (alternatively). All the live inode map
//                              indirect blocks and the live free bitmap
//                              blocks are appended to the log. The free 
//                              block lists are synchronized.
// ===========================================================================
// * INOUTS
//   int free_previous          1: free previous checkpoint's indirect blocks
//                                 and free bitmap blocks
//                              0: don't free anything
//
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_log_checkpoint(Context *context, int free_previous) {
  FSCurrentState *state;
  unsigned int blk_num;
  int ret;
  int i;

  // Get filesystem state variable
  state = context->fs_state;
  
  // Do we need to free the previous checkpoint stuff?
  if (free_previous) {

    // Free in memory all previous checkpoint indirect blocks
    blk_num = state->chk->indirect;
    for (i = 0; i < state->num_chk_ind; i++) {
      if (fs_mem_set_free_block(context, blk_num)) {
        ar_panic("fs_log_checkpoint: could not free checkpoint indirect");
      }
      blk_num = state->chk_ind[i]->indirect;
    }

    // Free in memory all previous free bitmap blocks (future version)
    blk_num = state->chk->free_bitmap;
    for (i = 0; i < state->num_chk_free; i++) {
      if (fs_mem_set_free_block(context, blk_num)) {
        ar_panic("fs_log_checkpoint: could not free checkpoint free bitmap");
      }
      blk_num = state->chk_future_free[i]->next;
    }
  }

  // Sanity check: in-memory live checkpoint should be ok
  if (state->chk->header != FS_HDR_TYPE_CHK_ROOT) {
    ar_panic("fs_log_checkpoint: live checkpoint failure");
  }

  // Append to the log the current checkpoint indirect blocks, last to first. 
  blk_num = 0;
  for (i = state->num_chk_ind - 1; i >= 0; i--) {
    if (state->chk_ind[i]->header != FS_HDR_TYPE_CHK_INDIRECT) {
      ar_panic("fs_log_checkpoint: live inode map indirect failure");
    }
    state->chk_ind[i]->indirect = blk_num;

    // Write the indirect block
    if ((ret = fs_log_append(context, (unsigned char *) state->chk_ind[i], 
                             &blk_num))) {
      return ret;
    }
  }

  // Stitch this list to the live checkpoint block
  state->chk->indirect = blk_num;

  // Append to the log the current free bitmap blocks, last to first. We are
  // dealing with the future version of the free bitmaps here, because the
  // checkpoint must be consistent with reality. 
  //
  // Note that here a race condition is introduced: we are appending to the
  // log the real free block list, but we are also changing it for each block
  // we are writing. This may mean that the free block bitmap which will be
  // read from the checkpoint at mount may still mark some of these blocks
  // as free. This is taken care at mount, by making sure that all the free
  // bitmap blocks are actually marked as used.
  blk_num = 0;
  for (i = state->num_chk_free - 1; i >= 0; i--) {
    if (state->chk_future_free[i]->header != FS_HDR_TYPE_CHK_FREE_BITMAP) {
      ar_panic("fs_log_checkpoint: live free bitmap failure");
    }
    state->chk_future_free[i]->next = blk_num;

    // Write the free bitmap block
    if ((ret = fs_log_append(context, 
                             (unsigned char *) state->chk_future_free[i], 
                             &blk_num))) {
      return ret;
    }
  }

  // Stitch this list to the live checkpoint block
  state->chk->free_bitmap = blk_num;

  // Assign current log seq ID to the checkpoint
  state->chk->log_seq_id = context->fs_state->log_seq_id++;

  // Checkpoint head is where our current tail is
  state->chk->log_head = context->fs_state->log_tail;

  // Write the checkpoint to the other end of the device than previously
  if (state->chk_pos) {
    if (fs_log_write_block((unsigned char *) state->chk, 0)) {
      return -EIO;
    }
  }
  else {
    if (fs_log_write_block((unsigned char *) state->chk, 
                           context->fs_num_blocks - 1)) {
      return -EIO;
    }
  }
  state->chk_pos = 1 - state->chk_pos;

  // Synchronize the free lists: the conservative view is now a 
  // snapshot of the real view
  fs_mem_sync_free_blocks(context);

  // Success
  return 0;
}


// ===========================================================================
// fs_log_format()              Initializes the disk structure and the current
//                              memory state to hold an empty filesystem
//
//                              WARNING: All data is lost when you call this.
// ===========================================================================
// * INPUTS
//   unsigned int log_seed      Initial value of the log sequential ID.
//                              Assign this randomly, so that upon repeated
//                              formats/mounts at the same block positions
//                              the log is not confused.
//
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_log_format(Context *context, unsigned int log_seed) {
  FSCurrentState *state;
  int ret;
  int i;

  // Get filesystem state variable
  state = context->fs_state;

  // All blocks are free, minus the two special checkpoint blocks
  fs_mem_init_free_blocks(context);
  if (fs_mem_unset_free_block(context, 0)) {
    ar_panic("fs_log_format: could not claim first block");
  }
  if (fs_mem_unset_free_block(context, context->fs_num_blocks - 1)) {
    ar_panic("fs_log_format: could not claim last block");
  }

  // Start preparing the live, in-memory checkpoint
  state->chk->header = FS_HDR_TYPE_CHK_ROOT;

  // Initialize log from first viable block that can be written
  // (blocks 0 and context->fs_num_blocks - 1 are reserved for checkpoints)
  state->log_seq_id = log_seed;
  state->log_tail = 1;

  // Create root inode. Root inode is inode #0.
  if ((ret = fs_lay_create_inode(context, 0, 
                FS_ATTR_DIRECTORY | 
                FS_ATTR_U_R | FS_ATTR_U_W | FS_ATTR_U_X |
                FS_ATTR_G_R | FS_ATTR_G_X | FS_ATTR_O_R | FS_ATTR_O_X,
                NULL, 0))) {
    return ret;
  }

  // Prepare checkpoint indirect blocks headers. They are already allocated
  // and zeroed by fs_mem_init().
  for (i = 0; i < state->num_chk_ind; i++) {
    state->chk_ind[i]->header = FS_HDR_TYPE_CHK_INDIRECT;
  }

  // Prepare checkpoint free bitmap headers. They are already allocated by
  // fs_mem_init() and initialized by fs_mem_init_free_blocks() above.
  for (i = 0; i < state->num_chk_free; i++) {
    
    // Future view is used normally
    state->chk_future_free[i]->header = FS_HDR_TYPE_CHK_FREE_BITMAP;

    // Conservative view is not used for on-disk purposes
    state->chk_free[i]->header = 0; 
  }

  // Take a checkpoint on the first block of device
  state->chk_pos = 1;
  if ((ret = fs_log_checkpoint(context, 0))) {
    return ret;
  }

  // Take a second checkpoint on the last block of device
  if ((ret = fs_log_checkpoint(context, 1))) {
    return ret;
  }

  // Success
  return 0;
}


// ===========================================================================
// fs_log_mount()               Reads an AntFS filesystem and initializes the
//                              memory state to bring it up to speed with it.
//                              The log is replayed from the latest valid
//                              checkpoint until an invalid log entry is
//                              encountered.
// ===========================================================================
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_log_mount(Context *context) {
  FSCurrentState *state;
  unsigned char buf0[4096], buf1[4096];
  FSCheckpointRoot *cand0;
  FSCheckpointRoot *cand1;
  int cand0_valid;
  int cand1_valid;
  int ret;
  int status;
  int i, j;
  unsigned int blk_num;
  unsigned int old_blk_num;
  unsigned int trans_blk_num;
  unsigned int trans_new_inode_num[FS_MAX_TRANSACTION_NEW_INODES];
  unsigned int trans_new_inode_blk_num[FS_MAX_TRANSACTION_NEW_INODES];
  unsigned int trans_free_inode_num[FS_MAX_TRANSACTION_FREE_INODES];
  unsigned int trans_free_blk_num[FS_MAX_TRANSACTION_FREE_BLOCKS];
  unsigned int trans_used_blk_num[FS_MAX_TRANSACTION_USED_BLOCKS];
  unsigned int trans_log_seq_id;
  unsigned int trans_log_tail;
  int cnt_trans_new_inode;
  int cnt_trans_free_inode;
  int cnt_trans_free_blk_num;
  int cnt_trans_used_blk_num;
  int keep;
  unsigned int header;
  unsigned int log_seq_id;
  int replay_cnt;
  FSInode *inode;
  FSIndirect *indirect;
  FSData *data;
  FSDelete *delete;


  // Get filesystem state variable
  state = context->fs_state;

  // Read both candidate checkpoints
  cand0_valid = 0;
  cand1_valid = 0;
  cand0 = NULL;
  cand1 = NULL;
  if (!fs_log_read_block(buf0, 0)) {
    cand0_valid = 1;
    cand0 = (FSCheckpointRoot *) buf0;
  }
  if (!fs_log_read_block(buf1, context->fs_num_blocks - 1)) {
    cand1_valid = 1;
    cand1 = (FSCheckpointRoot *) buf1;
  }

  // Check their headers
  if (cand0_valid) {
    if (cand0->header != FS_HDR_TYPE_CHK_ROOT) {
      cand0_valid = 0;
    }
  }
  if (cand1_valid) {
    if (cand1->header != FS_HDR_TYPE_CHK_ROOT) {
      cand1_valid = 0;
    }
  }

  // Did any of them survive?
  if ((!cand0_valid) && (!cand1_valid)) {
    return -EIO;
  }

  // If both survived, pick a winner based on the newest log seq ID
  if (cand0_valid && cand1_valid) {

    // Wrap-around scenario #1: cand0 is in first quadrant, cand1 is in last
    if ((cand0->log_seq_id < (1 << 30)) &&
        (cand1->log_seq_id > 3 * (unsigned int) (1 << 30))) {
      cand1_valid = 0;
    }

    // Wrap-around scenario #2: cand1 is in first quadrant, cand0 is in last
    else if ((cand1->log_seq_id < (1 << 30)) &&
        (cand0->log_seq_id > 3 * (unsigned int) (1 << 30))) {
      cand0_valid = 0;
    }

    // Default: take the highest log seq ID
    else {
      if (cand0->log_seq_id < cand1->log_seq_id) {
        cand0_valid = 0;
      }
      else {
        cand1_valid = 0;
      }
    }
  }

  // Copy the winner to the live checkpoint
  if (cand0_valid) {
    if (cand1_valid) {
      ar_panic("fs_log_mount: candidate selection inconsistency");
    }
    kt_memcpy(state->chk, cand0, 4096);
    state->chk_pos = 0;
  }
  else if (cand1_valid) {
    kt_memcpy(state->chk, cand1, 4096);
    state->chk_pos = 1;
  }
  else {
    ar_panic("fs_log_mount: candidate selection inconsistency");
  }

  // Set log tail and seq ID
  state->log_tail = state->chk->log_head;
  state->log_seq_id = state->chk->log_seq_id + 1;


  // Read all checkpoint free bitmap blocks
  blk_num = state->chk->free_bitmap;
  for (i = 0; i < state->num_chk_free; i++) {

    // Sanity check the block number
    if ((blk_num < 1) || (blk_num >= context->fs_num_blocks - 1)) {
      return -EIO;
    }

    // Read the free bitmap block
    if ((ret = fs_log_read_block((unsigned char *) state->chk_future_free[i], 
                                 blk_num))) {
      return ret;
    }

    // Get next free bitmap block number
    blk_num = state->chk_future_free[i]->next;
  }

  // It may have happened that the free bitmap blocks themselves are
  // still marked free, because we write them and updating the free list
  // for them at the same time. Fix this here.
  blk_num = state->chk->free_bitmap;
  for (i = 0; i < state->num_chk_free; i++) {

    // Check free block status
    if (fs_mem_free_block_status(context, blk_num, 1, &status)) {
      ar_panic("AntFS: fs_log_mount could not read free block status");
    }
    if (status) {
      
      // Mark it as used
      if (fs_mem_unset_free_block(context, blk_num)) {
        ar_panic("fs_log_format: could not claim last block");
      }
    }

    // Get next free bitmap block number
    blk_num = state->chk_future_free[i]->next;
  }

  // Synchronize the free lists: the conservative view is now a 
  // snapshot of the real view
  fs_mem_sync_free_blocks(context);

  // Count the free blocks in the bitmap, so we can synchronize the
  // free block counters
  state->num_free_blocks = 0;
  for (i = 0; i < context->fs_num_blocks; i++) {
    if (fs_mem_free_block_status(context, i, 1, &status)) {
      ar_panic("AntFS: fs_log_mount could not read initial free block status");
    }
    state->num_free_blocks += status;
  }
  state->num_future_free_blocks = state->num_free_blocks;

  
  // Read all checkpoint indirect blocks
  blk_num = state->chk->indirect;
  for (i = 0; i < state->num_chk_ind; i++) {

    // Sanity check the block number
    if ((blk_num < 1) || (blk_num >= context->fs_num_blocks - 1)) {
      return -EIO;
    }

    // Read the indirect
    if ((ret = fs_log_read_block((unsigned char *) state->chk_ind[i], 
                                 blk_num))) {
      return ret;
    }

    // Verify that this indirect block is not free. This must always be
    // correct, because the indirects are written before the free bitmap.
    if (fs_mem_free_block_status(context, blk_num, 1, &status)) {
      ar_panic("AntFS: fs_log_mount could not read free block status");
    }
    if (status) {
      return -EIO;
    }

    // Get next indirect block number
    blk_num = state->chk_ind[i]->indirect;
  }

  // Say something
  kt_printf("Valid AntFS checkpoint loaded. Replaying log...\r\n");

  // For all transactions
  replay_cnt = 0;
  while (1) {

    // Initialize transaction counters and valid bits
    cnt_trans_new_inode = 0;
    cnt_trans_free_inode = 0;
    cnt_trans_free_blk_num = 0;
    cnt_trans_used_blk_num = 0;
    trans_log_seq_id = state->log_seq_id;
    trans_log_tail = state->log_tail;
    trans_blk_num = 0;
    i = 0;

    // For all blocks of the current transaction
    while (1) {
      
      // Find next free block according to conservative view
      if (fs_mem_find_free_block(context, &blk_num)) {
        ar_panic("AntFS: fs_log_mount could not find next free block");
      }
      
      // Remember first block number it for the informational message
      if (i == 0) {
        trans_blk_num = blk_num;
      }

      // Try reading the block
      if ((ret = fs_log_read_block(buf0, blk_num))) {

        // I/O errors must be reported to top level
        if (ret == -EIO) {
          return ret;
        }

        // CRC errors may mean that a block was half-written or that
        // the log had never reached this point before. Simply consider
        // that transactions end here.
        else if (ret == -ECRC) {
          keep = 0;
          break;
        }

        // Other unknown errors
        else {
          return ret;
        }
      }

      // Parse vitals
      kt_memcpy(&header, buf0, 4);
      kt_memcpy(&log_seq_id, buf0 + 4, 4);

      // Check for expected log seq ID
      if (log_seq_id != state->log_seq_id) {
        keep = 0;
        break;
      }

      // The first block should have ATS set
      if ((i == 0) && !(header & FS_HDR_ATOMIC_START)) {
        keep = 0;
        break;
      }

      // Check for valid header. Note that all checkpoint-related structures
      // should not be encountered: the current checkpoint will have written
      // them before the log head. Older structures (by previous checkpoint) 
      // will be marked as free, and their blocks should be reused if they
      // appear in the log order with a valid log seq ID.
      inode = NULL;
      indirect = NULL;
      data = NULL;
      delete = NULL;
      if (header & FS_HDR_TYPE_INODE) {
        inode = (FSInode *) buf0;
      }
      else if (header & FS_HDR_TYPE_INDIRECT) {
        indirect = (FSIndirect *) buf0;
      }
      else if (header & FS_HDR_TYPE_DATA) {
        data = (FSData *) buf0;
      }
      else if (header & FS_HDR_TYPE_DELETE) {
        delete = (FSDelete *) buf0;
      }
      else {
        keep = 0;
        break;
      }

      // Remember that this block number should be non-free. 
      if (cnt_trans_used_blk_num >= FS_MAX_TRANSACTION_USED_BLOCKS) {
        ar_panic("AntFS: fs_log_mount used blocks static buffer overflow");
      }
      trans_used_blk_num[cnt_trans_used_blk_num++] = blk_num;

      // Note: if it's a Delete block, we'll also enter it to be freed below.
      // This is consistent with how it's done originally, because it must lead
      // to the situation of the block being freed only in the original view
      // (and not the future one).


      // Is it an inode?
      if (inode) {

        // Remember new inode numbers and block numbers
        if (cnt_trans_new_inode >= FS_MAX_TRANSACTION_NEW_INODES) {
          ar_panic("AntFS: fs_log_mount new inodes static buffer overflow");
        }
        trans_new_inode_num[cnt_trans_new_inode] = inode->inode_num;
        trans_new_inode_blk_num[cnt_trans_new_inode++] = blk_num;

        // Remember old block numbers to free
        if (inode->free_block0) {
          if (cnt_trans_free_blk_num >= FS_MAX_TRANSACTION_FREE_BLOCKS) {
            ar_panic("AntFS: fs_log_mount free blocks static buffer overflow");
          }
          trans_free_blk_num[cnt_trans_free_blk_num++] = inode->free_block0;
        }
        if (inode->free_block1) {
          if (cnt_trans_free_blk_num >= FS_MAX_TRANSACTION_FREE_BLOCKS) {
            ar_panic("AntFS: fs_log_mount free blocks static buffer overflow");
          }
          trans_free_blk_num[cnt_trans_free_blk_num++] = inode->free_block1;
        }
      }

      // Is it an indirect block?
      if (indirect) {

        // Remember old block number to free
        if (indirect->free_block) {
          if (cnt_trans_free_blk_num >= FS_MAX_TRANSACTION_FREE_BLOCKS) {
            ar_panic("AntFS: fs_log_mount free blocks static buffer overflow");
          }
          trans_free_blk_num[cnt_trans_free_blk_num++] = indirect->free_block;
        }
      }

      // Is it a data block?
      if (data) {

        // Remember old block number to free
        if (data->free_block) {
          if (cnt_trans_free_blk_num >= FS_MAX_TRANSACTION_FREE_BLOCKS) {
            ar_panic("AntFS: fs_log_mount free blocks static buffer overflow");
          }
          trans_free_blk_num[cnt_trans_free_blk_num++] = data->free_block;
        }
      }

      // Is it a delete block?
      if (delete) {
        
        // The delete block itself should be freed
        if (cnt_trans_free_blk_num >= FS_MAX_TRANSACTION_FREE_BLOCKS) {
          ar_panic("AntFS: fs_log_mount free blocks static buffer overflow");
        }
        trans_free_blk_num[cnt_trans_free_blk_num++] = blk_num;

        // Remember old inode number to free
        if (delete->free_inode) {
          if (cnt_trans_free_inode >= FS_MAX_TRANSACTION_FREE_INODES) {
            ar_panic("AntFS: fs_log_mount free inodes static buffer overflow");
          }
          trans_free_inode_num[cnt_trans_free_inode++] = delete->free_inode;
        }

        // Remember old block numbers to free
        for (j = 0; j < FS_DELETE_NUM_BLOCKS; j++) {

          // Read next block to free
          if (fs_mem_get_data_block_ptr(context, delete, j, &old_blk_num)) {
            ar_panic("AntFS: fs_log_mount could not read delete block ptrs");
          }

          // A zero entry means we're done
          if (!old_blk_num) {
            break;
          }

          // Remember block number
          if (cnt_trans_free_blk_num >= FS_MAX_TRANSACTION_FREE_BLOCKS) {
            ar_panic("AntFS: fs_log_mount free blocks static buffer overflow");
          }
          trans_free_blk_num[cnt_trans_free_blk_num++] = old_blk_num;
        }
      }

      // Advance in-transaction index
      i++;

      // Advance log seq ID
      state->log_seq_id++;

      // Advance log tail
      if (blk_num >= context->fs_num_blocks - 2) {
        state->log_tail = 1;
      }
      else {
        state->log_tail = blk_num + 1;
      }

      // If ATE is set, we're done
      if (header & FS_HDR_ATOMIC_END) {
        keep = 1;
        break;
      }
    }

    // Was the transaction complete?
    if (keep) {

      // Inform about what we're doing
      kt_printf(
  "  Replaying transaction %d [%d blocks, log entry 0x%08X, blk 0x%05X]\r\n",
                replay_cnt, i, trans_log_seq_id, trans_blk_num);

      // Replay the transaction's used blocks
      for (j = 0; j < cnt_trans_used_blk_num; j++) {
        if (fs_mem_unset_free_block(context, trans_used_blk_num[j])) {
          ar_panic("AntFS: fs_log_mount could not replay used blocks");
        }
      }

      // Replay the transaction's free blocks
      for (j = 0; j < cnt_trans_free_blk_num; j++) {
        if (fs_mem_set_free_block(context, trans_free_blk_num[j])) {
          ar_panic("AntFS: fs_log_mount could not replay free blocks");
        }
      }

      // Replay the transaction's new inodes
      for (j = 0; j < cnt_trans_new_inode; j++) {
        if (fs_mem_update_inode_map(context, trans_new_inode_num[j],
                                    trans_new_inode_blk_num[j])) {
          ar_panic("AntFS: fs_log_mount could not replay new inodes");
        }
      }

      // Replay the transaction's free inodes
      for (j = 0; j < cnt_trans_free_inode; j++) {
        if (fs_mem_update_inode_map(context, trans_free_inode_num[j], 0)) {
          ar_panic("AntFS: fs_log_mount could not replay free inodes");
        }
      }

      // Remember how many transactions we've completed
      replay_cnt++;
    }

    // Abort transaction and exit log replay
    else {
      
      // Restore anything we touched as it was at the start of the latest 
      // transaction
      state->log_seq_id = trans_log_seq_id;
      state->log_tail = trans_log_tail;

      // No more transactions
      break;
    }
  }

  // Say something more
  if (replay_cnt) {
    kt_printf("%d transactions replayed successfully.\r\n", replay_cnt);
  }
  else {
    kt_printf("Filesystem was clean.\r\n");
  }

  // Success
  return 0;
}


