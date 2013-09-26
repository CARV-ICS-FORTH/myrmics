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
// Abstract      : AntFS functions that handle the on-disk layout of data
//                 and metadata structures
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: layout.c,v $
// CVS revision  : $Revision: 1.8 $
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
// fs_lay_create_inode()        Creates a new, empty inode for either a 
//                              file or a directory. This function will create:
//
//                              1. A new block for the inode
// ===========================================================================
// * INPUTS
//   unsigned int inode_num     Inode number to be created
//   unsigned int attr          Attributes of the inode to be created
//   unsigned int atomic_flags  ATS/ATE flags to be used for the transaction
//
// * OUTPUTS
//   unsigned int *ret_blk_num  If not NULL, returned block number of the 
//                              created inode
// * INOUTS
//   Context *context           Current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_create_inode(Context *context, unsigned int inode_num, 
                        unsigned int attr, unsigned int *ret_blk_num, 
                        unsigned int atomic_flags) {

  unsigned char buf[4096];
  FSInode *inode_blk;
  unsigned int inode_blk_num;
  int ret;

  // Create the new inode block
  kt_memset(buf, 0, 4096);
  inode_blk = (FSInode *) buf;
  inode_blk->header |= atomic_flags;
  inode_blk->header |= FS_HDR_TYPE_INODE;

  inode_blk->inode_num = inode_num;
  inode_blk->attr = attr;

  // Write the inode block to disk
  if ((ret = fs_log_append(context, buf, &inode_blk_num))) {
    return ret;
  }

  // Update inode map
  if (fs_mem_update_inode_map(context, inode_num, inode_blk_num)) {
    ar_panic("AntFS: create_inode: Could not set inode in inode map");
  }

  // Return new inode block number
  if (ret_blk_num) {
    *ret_blk_num = inode_blk_num;
  }

  return 0;
}


// ===========================================================================
// fs_lay_find_last_entry()     Scans a directory, by reading its inode's
//                              data blocks containing the directory entries,
//                              and finds the last directory entry.
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned int inode_blk_num The directory inode number to be searched
//
// * OUTPUTS
//   unsigned char              If not NULL, the directory inode contents
//        *ret_inode_blk        are stored in this caller-allocated buffer
//   unsigned char              If not NULL, the dirdata block that contains
//        *ret_dirdata_blk      the found dir entry is stored in this 
//                              caller-allocated buffer
//   unsigned int               If not NULL, the dirdata block number that 
//        *ret_dirdata_blk_num  contains the found dir entry is stored here
//   int *ret_dirdata_index     If not NULL, the index of the dirdata block
//                              containing the direntry is stored here. If
//                              directory is empty, -1 will be returned here.
//   int *ret_direntry_index    If not NULL, the index of the direntry in its
//                              dirdata block is stored here
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_find_last_entry(Context *context, unsigned int inode_blk_num,
                           unsigned char *ret_inode_blk,
                           unsigned char *ret_dirdata_blk,
                           unsigned int *ret_dirdata_blk_num,
                           int *ret_dirdata_index,
                           int *ret_direntry_index) {

  unsigned char *inode_blk;
  FSInode *inode;
  unsigned char *dirdata_blk;
  unsigned int dirdata_blk_num;
  unsigned int last_dirdata_blk_num;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  int ret;
  int i;
  int dirdata_index;
  int direntry_index;

  // Buffer management
  if (ret_inode_blk) {
    inode_blk = ret_inode_blk;
  }
  else {
    inode_blk = kt_malloc(4096);
  }
  if (ret_dirdata_blk) {
    dirdata_blk = ret_dirdata_blk;
  }
  else {
    dirdata_blk = kt_malloc(4096);
  }

  // Read directory inode block
  if ((ret = fs_log_read_block(inode_blk, inode_blk_num))) {
    if (!ret_inode_blk) {
      kt_free(inode_blk);
    }
    if (!ret_dirdata_blk) {
      kt_free(dirdata_blk);
    }
    return ret;
  }

  // Sanity checks
  inode = (FSInode *) inode_blk;
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic("AntFS: find_last_entry called on a non-inode");
  }
  if (!(inode->attr & FS_ATTR_DIRECTORY)) {
    ar_panic("AntFS: find_last_entry called on a non-directory");
  }

  // Initialize all return fields to 0
  dirdata_blk_num = 0;
  dirdata_index = 0;
  direntry_index = 0;
  last_dirdata_blk_num = 0;

  // For all its data blocks
  for (i = 0; i < FS_INODE_NUM_DATA_BLOCKS; i++) {

    // Get the next block number
    if (fs_mem_get_data_block_ptr(context, inode_blk, i, &dirdata_blk_num)) {
      ar_panic("AntFS: find_last_entry could not get data block pointer");
    }

    // Are the blocks finished?
    if (!dirdata_blk_num) {
      dirdata_index = i - 1;
      break;
    }

    // Remember this block number
    last_dirdata_blk_num = dirdata_blk_num;
  }

  // Was there a dirdata block at all?
  if (dirdata_index < 0) {
    // All return fields will be set to 0, dirdata_index will be -1
    goto finished;
  }

  // Read the last dirdata block
  dirdata_blk_num = last_dirdata_blk_num;
  if ((ret = fs_log_read_block(dirdata_blk, dirdata_blk_num))) {
    return ret;
  }

  // Sanity check
  dirdata = (FSData *) dirdata_blk;
  if (!(dirdata->header & FS_HDR_TYPE_DATA)) {
    ar_panic("AntFS: find_last_entry: inode pointer to to non-data block");
  }
  dirdata_entries = (FSDirData *) dirdata->data;

  // For all the directory entries of this dirdata block
  for (i = 0; i < FS_MAX_DIRENTRIES; i++) {
    if (!dirdata_entries->dir_entry[i].inode_num) {
      // inode_num = 0 is the NULL entry; last one is the one before that.
      // Note that i cannot be 0, because then it's an empty dirdata block,
      // which is not allowed. The previous check would have already returned
      // from this function.
      if (!i) {
        ar_panic("AntFS: find_last_entry: empty dirdata block");
      }
      break;
    }
  }
  // Either we reached the i-th NULL entry and we breaked (so we need i-1)
  // or loop was finished with all FS_MAX_DIRENTRIES non-NULL (so we still
  // need i-1)
  direntry_index = i - 1;

finished:
  // Return stuff and free related buffers 
  if (!ret_inode_blk) {
    kt_free(inode_blk);
  }
  if (!ret_dirdata_blk) {
    kt_free(dirdata_blk);
  }
  if (ret_dirdata_blk_num) {
    *ret_dirdata_blk_num = dirdata_blk_num;
  }
  if (ret_dirdata_index) {
    *ret_dirdata_index = dirdata_index;
  }
  if (ret_direntry_index) {
    *ret_direntry_index = direntry_index;
  }

  return 0;
}

// ===========================================================================
// fs_lay_find_entry()          Scans a directory, by following its inode's
//                              data blocks containing the directory entries,
//                              trying to find the given file/dir name.
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned int inode_blk_num The directory inode number to be searched
//   char *name                 Name of file/dir to be searched for
//
// * OUTPUTS
//   unsigned char              If not NULL, the directory inode contents
//        *ret_inode_blk        are stored in this caller-allocated buffer
//   unsigned char              If not NULL, the dirdata block that contains
//        *ret_dirdata_blk      the found dir entry is stored in this 
//                              caller-allocated buffer
//   unsigned int               If not NULL, the dirdata block number that 
//        *ret_dirdata_blk_num  contains the found dir entry is stored here
//   int *ret_dirdata_index     If not NULL, the index of the dirdata block
//                              containing the direntry is stored here
//   int *ret_direntry_index    If not NULL, the index of the direntry in its
//                              dirdata block is stored here
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_find_entry(Context *context, unsigned int inode_blk_num, 
                      char *name, unsigned char *ret_inode_blk,
                      unsigned char *ret_dirdata_blk,
                      unsigned int *ret_dirdata_blk_num,
                      int *ret_dirdata_index,
                      int *ret_direntry_index) {
  
  unsigned char *inode_blk;
  FSInode *inode;
  unsigned char *dirdata_blk;
  unsigned int dirdata_blk_num;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  int ret;
  int i, j;
  int found;

  // Buffer management
  if (ret_inode_blk) {
    inode_blk = ret_inode_blk;
  }
  else {
    inode_blk = kt_malloc(4096);
  }
  if (ret_dirdata_blk) {
    dirdata_blk = ret_dirdata_blk;
  }
  else {
    dirdata_blk = kt_malloc(4096);
  }

  // Read directory inode block
  if ((ret = fs_log_read_block(inode_blk, inode_blk_num))) {
    if (!ret_inode_blk) {
      kt_free(inode_blk);
    }
    if (!ret_dirdata_blk) {
      kt_free(dirdata_blk);
    }
    return ret;
  }

  // Sanity checks
  inode = (FSInode *) inode_blk;
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic("AntFS: find_entry called on a non-inode");
  }
  if (!(inode->attr & FS_ATTR_DIRECTORY)) {
    ar_panic("AntFS: find_entry called on a non-directory");
  }

  // Initialize all return fields to 0
  dirdata_blk_num = 0;
  j = 0;
  found = 0;

  // For all its data blocks
  for (i = 0; i < FS_INODE_NUM_DATA_BLOCKS; i++) {

    // Get the next block number
    if (fs_mem_get_data_block_ptr(context, inode_blk, i, &dirdata_blk_num)) {
      ar_panic("AntFS: find_entry could not get data block pointer");
    }

    // Are the blocks finished?
    if (!dirdata_blk_num) {
      goto finished;
    }

    // Read the dirdata block
    if ((ret = fs_log_read_block(dirdata_blk, dirdata_blk_num))) {
      return ret;
    }

    // Sanity check
    dirdata = (FSData *) dirdata_blk;
    if (!(dirdata->header & FS_HDR_TYPE_DATA)) {
      ar_panic("AntFS: find_entry: inode pointer to to non-data block");
    }
    dirdata_entries = (FSDirData *) dirdata->data;

    // For all the directory entries of this dirdata block
    for (j = 0; j < FS_MAX_DIRENTRIES; j++) {

      // Are we finished?
      if (!dirdata_entries->dir_entry[j].inode_num) {
        // Sanity check
        if (!j) {
          ar_panic("AntFS: find_entry: empty dirdata block");
        }
        goto finished;
      }

      // Check the name
      if (!kt_strcmp(dirdata_entries->dir_entry[j].name, name)) {
        found = 1;
        goto finished;
      }
    }
  }

finished:

  // Return stuff and free related buffers 
  if (!ret_inode_blk) {
    kt_free(inode_blk);
  }
  if (!ret_dirdata_blk) {
    kt_free(dirdata_blk);
  }
  if (ret_dirdata_blk_num) {
    *ret_dirdata_blk_num = dirdata_blk_num;
  }
  if (ret_dirdata_index) {
    *ret_dirdata_index = i;
  }
  if (ret_direntry_index) {
    *ret_direntry_index = j;
  }

  // Return error or success
  if (!found) {
    return -ENOENT;
  }

  return 0;
}


// ===========================================================================
// fs_lay_add_dir_entry()       Adds a directory entry to a given inode.
//                              This function will modify blocks and replace
//                              them with new in the log as follows:
//
//                              (a) If last dirdata block has space:
//
//                                  1. Modified dirdata block
//                                  2. Modified inode block
//
//                              (b) If last dirdata doesn't exist or if it
//                                  doesn't have space:
//
//                                  1. New dirdata block
//                                  2. Modified inode block
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned char *inode_blk   The parent inode block
//   unsigned int inode_blk_num The parent inode block number
//   unsigned char              The dirdata block that contains the last
//      *last_dirdata_blk       direntry of the inode
//   unsigned int               The block number of last_dirdata_blk
//      *last_dirdata_blk_num   
//   int last_dirdata_index     The index of last_dirdata block in the inode
//   int last_direntry_index    The index of the direntry in the dirdata block
//   char *new_name             New entry name
//   unsigned int new_inode_num New entry inode number
//   unsigned int atomic_flags  ATS/ATE flags to be used for the transaction
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_add_dir_entry(Context *context, unsigned char *inode_blk,
                         unsigned int inode_blk_num,
                         unsigned char *last_dirdata_blk, 
                         unsigned int last_dirdata_blk_num,
                         int last_dirdata_index, 
                         int last_direntry_index, 
                         char *new_name, unsigned int new_inode_num,
                         unsigned int atomic_flags) {
  FSInode *inode;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  int ret;
  unsigned char buf[4096];

  // Sanity checks
  inode = (FSInode *) inode_blk;

  if (!inode_blk) {
    ar_panic("AntFs: fs_lay_add_dir_entry called with no inode_blk");
  }
  if ((inode_blk_num < 1) || (inode_blk_num >= context->fs_num_blocks)) {
    ar_panic("AntFs: fs_lay_add_dir_entry called with invalid inode_blk_num");
  }
  if ((last_dirdata_index < -1) || 
      (last_dirdata_index >= FS_INODE_NUM_DATA_BLOCKS)) {
    ar_panic("AntFs: fs_lay_add_dir_entry called with bad last_dirdata_index");
  }
  if ((last_direntry_index < 0) || (last_direntry_index >= FS_MAX_DIRENTRIES)) {
    ar_panic("AntFs: fs_lay_add_dir_entry called with bad last_direntry_index");
  }
  if ((!new_name) || (kt_strlen(new_name) >= FS_MAX_DIRENTRY_NAME)) {
    ar_panic("AntFs: fs_lay_add_dir_entry called with bad new_name");
  }
  if ((new_inode_num < 1) || (new_inode_num >= FS_MAX_INODES)) {
    ar_panic("AntFs: fs_lay_add_dir_entry called with bad new_inode_num");
  }
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic("AntFS: fs_lay_add_dir_entry called on a non-inode");
  }
  if (!(inode->attr & FS_ATTR_DIRECTORY)) {
    ar_panic("AntFS: fs_lay_add_dir_entry called on a non-directory");
  }
  if (inode->size >= FS_INODE_NUM_DATA_BLOCKS * FS_MAX_DIRENTRIES - 1) {
    return -EOVERFLOW;
  }
  

  // Does last dirdata block even exist? Does it also have space?
  if ((last_dirdata_index > -1) && 
      (last_direntry_index < FS_MAX_DIRENTRIES - 1)) {

    // We'll modify the existing last_dirdata_blk
    last_direntry_index++;
    dirdata = (FSData *) last_dirdata_blk;
    dirdata_entries = (FSDirData *) dirdata->data;

    // Add the new entry
    kt_strcpy(dirdata_entries->dir_entry[last_direntry_index].name, new_name);
    dirdata_entries->dir_entry[last_direntry_index].inode_num = new_inode_num;

    // The old block is now free
    if (fs_mem_set_free_block(context, last_dirdata_blk_num)) {
      ar_panic(
          "AntFS: fs_lay_add_dir_entry could not free last_dirdata_blk_num");
    }
    dirdata->free_block = last_dirdata_blk_num;

    // Clear leftover atomic flags and put our of flags to it
    dirdata->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
    if (atomic_flags & FS_HDR_ATOMIC_START) {
      dirdata->header |= FS_HDR_ATOMIC_START;
    }
    
    // Write modified block to disk
    if ((ret = fs_log_append(context, last_dirdata_blk, 
                             &last_dirdata_blk_num))) {
      return ret;
    }

  }
  else {
    // Old block is full, or does not exist at all. We don't touch it, but
    // we'll add a new dirdata block.
    last_dirdata_index++;
    kt_memset(buf, 0, 4096);
    dirdata = (FSData *) buf;
    dirdata_entries = (FSDirData *) dirdata->data;

    // Set header and part of our atomic flags
    dirdata->header |= FS_HDR_TYPE_DATA;
    if (atomic_flags & FS_HDR_ATOMIC_START) {
      dirdata->header |= FS_HDR_ATOMIC_START;
    }
    
    // Add the new entry
    kt_strcpy(dirdata_entries->dir_entry[0].name, new_name);
    dirdata_entries->dir_entry[0].inode_num = new_inode_num;

    // Write new block to disk
    if ((ret = fs_log_append(context, buf, &last_dirdata_blk_num))) {
      return ret;
    }
  }


  // Modify inode block
  inode->mtime++;       // If we had a real-time clock, this would be set to it
  inode->size++;        // One more directory entry
  
  // Add the (updated or new) pointer to last dirdata block
  if (fs_mem_set_data_block_ptr(context, inode, 
                            last_dirdata_index, last_dirdata_blk_num)) {
    ar_panic("AntFS: fs_lay_add_dir_entry: Could not set dirdata prt in inode");
  }

  // The old inode block is now free
  if (fs_mem_set_free_block(context, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_add_dir_entry could not free inode_blk_num");
  }
  inode->free_block0 = inode_blk_num;
  inode->free_block1 = 0;

  // Clear leftover atomic flags and put our of flags to it
  inode->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
  if (atomic_flags & FS_HDR_ATOMIC_END) {
    inode->header |= FS_HDR_ATOMIC_END;
  }

  // Write modified inode block to disk
  if ((ret = fs_log_append(context, inode_blk, &inode_blk_num))) {
    return ret;
  }

  // Update inode map with new inode block number
  if (fs_mem_update_inode_map(context, inode->inode_num, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_add_dir_entry could not set inode in inode map");
  }

  // Success
  return 0;
}
                     

// ===========================================================================
// fs_lay_delete_dir_entry()    Deletes a directory entry from a directory.
//                              This function will modify blocks and replace
//                              them with new in the log as follows:
//
//                              (a) If entry being deleted and the last entry
//                                  are not in the same data block, and the
//                                  deletion from the last block is not
//                                  emptying it:
//
//                                  1. Modified dirdata block
//                                  2. Modified last_dirdata block
//                                  3. Modified inode block
//
//                              (b) If entry being deleted and the last entry
//                                  are in the same block, and there
//                                  are more entries left:
//
//                                  1. Modified dirdata block
//                                  2. Modified inode block
//
//                              (c) If we're deleting the only entry left 
//                                  in the directory:
//
//                                  1. Modified inode block
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned char *inode_blk   The parent inode block
//   unsigned int inode_blk_num The parent inode block number
//   unsigned char              The dirdata block that contains the direntry
//      *dirdata_blk            which we're deleting from the inode
//   unsigned int               The block number of dirdata_blk
//      *dirdata_blk_num   
//   int dirdata_index          The index of dirdata block in the inode
//   int direntry_index         The index of the direntry in the dirdata block
//   unsigned char              The dirdata block that contains the last
//      *last_dirdata_blk       direntry of the inode
//   unsigned int               The block number of last_dirdata_blk
//      *last_dirdata_blk_num   
//   int last_dirdata_index     The index of last_dirdata block in the inode
//   int last_direntry_index    The index of the direntry in the dirdata block
//   unsigned int atomic_flags  ATS/ATE flags to be used for the transaction
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_delete_dir_entry(Context *context, unsigned char *inode_blk,
                            unsigned int inode_blk_num, 
                            unsigned char *dirdata_blk, 
                            unsigned int dirdata_blk_num,
                            int dirdata_index, int direntry_index, 
                            unsigned char *last_dirdata_blk, 
                            unsigned int last_dirdata_blk_num,
                            int last_dirdata_index, int last_direntry_index, 
                            unsigned int atomic_flags) {
  FSInode *inode;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  FSData *last_dirdata;
  FSDirData *last_dirdata_entries;
  int ret;
  int start_written;

  // Sanity checks
  inode = (FSInode *) inode_blk;

  if (!inode_blk) {
    ar_panic("AntFs: fs_lay_delete_dir_entry called with no inode_blk");
  }
  if ((inode_blk_num < 1) || (inode_blk_num >= context->fs_num_blocks)) {
    ar_panic(
        "AntFs: fs_lay_delete_dir_entry called with invalid inode_blk_num");
  }
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic("AntFS: fs_lay_delete_dir_entry called on a non-inode");
  }
  if (!(inode->attr & FS_ATTR_DIRECTORY)) {
    ar_panic("AntFS: fs_lay_delete_dir_entry called on a non-directory");
  }
  if (!inode->size) {
    ar_panic("AntFS: fs_lay_delete_dir_entry called on empty dir");
  }
  if ((last_dirdata_index < -1) || 
      (last_dirdata_index >= FS_INODE_NUM_DATA_BLOCKS)) {
    ar_panic(
        "AntFs: fs_lay_delete_dir_entry called with bad last_dirdata_index");
  }
  if ((last_direntry_index < 0) || (last_direntry_index >= FS_MAX_DIRENTRIES)) {
    ar_panic(
        "AntFs: fs_lay_delete_dir_entry called with bad last_direntry_index");
  }
  if ((dirdata_index < -1) || (dirdata_index >= FS_INODE_NUM_DATA_BLOCKS)) {
    ar_panic(
        "AntFs: fs_lay_delete_dir_entry called with bad dirdata_index");
  }
  if ((direntry_index < 0) || (direntry_index >= FS_MAX_DIRENTRIES)) {
    ar_panic(
        "AntFs: fs_lay_delete_dir_entry called with bad direntry_index");
  }
  if ((dirdata_index * FS_MAX_DIRENTRIES + direntry_index) > 
      (last_dirdata_index * FS_MAX_DIRENTRIES + last_direntry_index)) {
    ar_panic(
        "AntFs: fs_lay_delete_dir_entry called with entry > last");
  }

  // Initialize convenience variables
  start_written = 0;
  dirdata = (FSData *) dirdata_blk;
  dirdata_entries = (FSDirData *) dirdata->data;
  
  // If last_dirdata and dirdata are the same block, alias the pointers
  // so that we can refer to one of the two buffers
  if (last_dirdata_index != dirdata_index) {
    last_dirdata = (FSData *) last_dirdata_blk;
    last_dirdata_entries = (FSDirData *) last_dirdata->data;
  }
  else {
    last_dirdata_blk = dirdata_blk;
    last_dirdata = dirdata;
    last_dirdata_entries = dirdata_entries;
  }

  // We can do the switch if last direntry is not the same with the one
  // we are deleting. This can happen either because it's a directory with
  // one entry (inode->size == 1) or simply because we're deleting the
  // last entry, or both of the above.
  if ((last_dirdata_index != dirdata_index) ||
      (last_direntry_index != direntry_index)) {

    // Copy last to current
    kt_memcpy(&(dirdata_entries->dir_entry[direntry_index]),
              &(last_dirdata_entries->dir_entry[last_direntry_index]),
              sizeof(FSDirEntry));
  }

  // We are going to write the modified dirdata block only if it's different
  // from the last_dirdata block. Otherwise, the modification will be 
  // written anyway when we're dealing with the last_dirdata below.
  if (last_dirdata_blk_num != dirdata_blk_num) {

    // The old dirdata block is now free
    if (fs_mem_set_free_block(context, dirdata_blk_num)) {
      ar_panic(
          "AntFS: fs_lay_delete_dir_entry could not free dirdata_blk_num");
    }
    dirdata->free_block = dirdata_blk_num;

    // Clear leftover atomic flags and put our of flags to it
    dirdata->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
    if (atomic_flags & FS_HDR_ATOMIC_START) {
      dirdata->header |= FS_HDR_ATOMIC_START;
    }
    start_written = 1;
    
    // Write modified dirdata block to disk
    if ((ret = fs_log_append(context, dirdata_blk, &dirdata_blk_num))) {
      return ret;
    }

    // Update inode with the new block number for this dirdata
    if (fs_mem_set_data_block_ptr(context, inode, 
                                  dirdata_index, dirdata_blk_num)) {
      ar_panic(
        "AntFS: fs_lay_delete_dir_entry: Could not set dirdata prt in inode");
    }
  }

  // Invalidate inode free blocks. They can be filled below or not, depending
  // on the case.
  inode->free_block0 = 0;
  inode->free_block1 = 0;
  
  // Is the last dirdata block being emptied?
  if (!last_direntry_index) {

    // Free old last_dirdata block
    if (fs_mem_set_free_block(context, last_dirdata_blk_num)) {
      ar_panic(
        "AntFS: fs_lay_delete_dir_entry could not free last_dirdata_blk_num");
    }
    
    // Mark it on the second inode entry for deletion
    inode->free_block1 = last_dirdata_blk_num;

    // Invalidate this data pointer on the inode
    if (fs_mem_set_data_block_ptr(context, inode, last_dirdata_index, 0)) {
      ar_panic(
       "AntFS: fs_lay_delete_dir_entry: Could not set last_dirdata in inode");
    }
  }
  else {
    
    // Delete last entry
    kt_memset(&(last_dirdata_entries->dir_entry[last_direntry_index]),
              0, sizeof(FSDirEntry));
    
    // The old last_dirdata block is now free
    if (fs_mem_set_free_block(context, last_dirdata_blk_num)) {
      ar_panic(
        "AntFS: fs_lay_delete_dir_entry could not free last_dirdata_blk_num");
    }
    last_dirdata->free_block = last_dirdata_blk_num;

    // Clear leftover atomic flags and put our of flags to it
    last_dirdata->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
    if ((!start_written) && (atomic_flags & FS_HDR_ATOMIC_START)) {
      last_dirdata->header |= FS_HDR_ATOMIC_START;
      start_written = 1;
    }
    
    // Write modified last_dirdata block to disk
    if ((ret = fs_log_append(context, last_dirdata_blk, 
                             &last_dirdata_blk_num))) {
      return ret;
    }

    // Update inode with the new block number for the last_dirdata
    if (fs_mem_set_data_block_ptr(context, inode, 
                                  last_dirdata_index, last_dirdata_blk_num)) {
      ar_panic(
        "AntFS: fs_lay_delete_dir_entry: Could not set last_dirdata in inode");
    }
  }

  // Modify inode block
  inode->mtime++;       // If we had a real-time clock, this would be set to it
  inode->size--;        // One less directory entry
  
  // The old inode block is now free
  if (fs_mem_set_free_block(context, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_delete_dir_entry could not free inode_blk_num");
  }
  inode->free_block0 = inode_blk_num;

  // Clear leftover atomic flags and put our of flags to it
  inode->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
  if ((!start_written) && (atomic_flags & FS_HDR_ATOMIC_START)) {
    inode->header |= FS_HDR_ATOMIC_START;
  }
  if (atomic_flags & FS_HDR_ATOMIC_END) {
    inode->header |= FS_HDR_ATOMIC_END;
  }

  // Write modified inode block to disk
  if ((ret = fs_log_append(context, inode_blk, &inode_blk_num))) {
    return ret;
  }

  // Update inode map with new inode block number
  if (fs_mem_update_inode_map(context, inode->inode_num, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_delete_dir_entry could not set inode in inode map");
  }

  // Success
  return 0;
}
                     

// ===========================================================================
// fs_lay_file_seek()           Starts from a file inode and tries to find
//                              the given data block index. It will read
//                              the inode and as many indirect blocks as it
//                              needs until the data block number is found.
//                              Finally, it will read the data block.
//
//                              The contents of the read inode, indirect 
//                              blocks and the data block are all returned
//                              to the caller.
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned int inode_num     Inode number of the file to be seeked
//   int data_blk_index         Which data block to read (0 ... last) or
//                              -1 to indicate "last block of file"
// * OUTPUTS
//   unsigned char              Caller buffer that will hold the file inode
//      *ret_inode_blk
//   unsigned int               If not NULL, the inode block number is
//      *ret_inode_blk_num      returned here
//   unsigned int               If not NULL, the actual last data block index
//      *ret_data_blk_index     will be returned here (equal to data_blk_index,
//                              unless the latter is -1, where the real one
//                              will be returned here).
//   unsigned char              Caller buffer that will hold the requested data
//      *ret_data_blk           block contents
//   unsigned int               If not NULL, the data block number is
//      *ret_inode_blk_num      returned here
//   unsigned char              Caller array of buffers that will hold the
//      *ret_indirect_blks      encountered indirect nodes. It must be big
//                              enough to hold the maximum file size in bytes,
//                              i.e. FS_MAX_FILE_SIZE * FS_MAX_INDIRECTS
//   unsigned int               Caller array that will hold the encountered
//      *ret_indirect_blk_nums  indirect nodes block numbers. See above for
//                              needed size.
//   int *ret_num_indirects     If not NULL, the number of actually 
//                              encountered indirect blocks will be returned
//                              here
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_file_seek(Context *context, unsigned int inode_num,
                     int data_blk_index, unsigned char *ret_inode_blk,
                     unsigned int *ret_inode_blk_num,
                     unsigned int *ret_data_blk_index,
                     unsigned char *ret_data_blk,
                     unsigned int *ret_data_blk_num,
                     unsigned char *ret_indirect_blks,
                     unsigned int *ret_indirect_blk_nums,
                     int *ret_num_indirects) {
  
  unsigned int inode_blk_num;
  unsigned int data_blk_num;
  int num_indirects;
  int data_blk_rem;
  int ret;
  FSInode *inode;
  int i;


  // Sanity checks
  if (inode_num >= FS_MAX_INODES) {
    ar_panic("AntFS: fs_lay_file_seek called with invalid inode num");
  }
  if (!ret_inode_blk) {
    ar_panic("AntFS: fs_lay_file_seek called with no ret_inode_blk");
  }
  if (!ret_data_blk) {
    ar_panic("AntFS: fs_lay_file_seek called with no ret_data_blk");
  }
  if (!ret_indirect_blks) {
    ar_panic("AntFS: fs_lay_file_seek called with no ret_indirect_blks");
  }
  if (!ret_indirect_blk_nums) {
    ar_panic("AntFS: fs_lay_file_seek called with no ret_indirect_blk_nums");
  }

  // Get file inode block number
  if (fs_mem_read_inode_map(context, inode_num, &inode_blk_num)) {
    ar_panic("AntFS: fs_lay_file_seek could not get inode block num");
  }

  // Read the file inode
  if ((ret = fs_log_read_block(ret_inode_blk, inode_blk_num))) {
    return ret;
  }
  inode = (FSInode *) ret_inode_blk;
  // Check it's a file
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic("AntFS: fs_lay_file_seek found inconsistent inode num");
  }
  if (!(inode->attr & FS_ATTR_FILE)) {
    ar_panic("AntFS: fs_lay_file_seek got inode num of a directory");
  }

  // Check file size
  if (!inode->size) {
    // No data block exists; return appropriate values
    if (ret_data_blk_num) {
      *ret_data_blk_num = 0;
    }
    if (ret_inode_blk_num) {
      *ret_inode_blk_num = inode_blk_num;
    }
    if (ret_num_indirects) {
      *ret_num_indirects = 0;
    }
    if (ret_data_blk_index) {
      *ret_data_blk_index = -1;
    }
    return 0;
  }

  // Handle end-of-file seeks
  if (data_blk_index == -1) {
    ar_int_divide(inode->size - 1, FS_DATA_BLOCK_PAYLOAD, 
                  &data_blk_index, NULL);
  }

  // Where is our data block? Do we need indirects?
  if (data_blk_index < FS_INODE_NUM_DATA_BLOCKS) {
    num_indirects = 0;
  }
  else {
    ar_int_divide(data_blk_index - FS_INODE_NUM_DATA_BLOCKS,
                  FS_IND_NUM_DATA_BLOCKS,
                  &num_indirects, &data_blk_rem);
    num_indirects++;

    if (num_indirects > FS_MAX_INDIRECTS) {
      ar_panic("AntFS: fs_lay_file_seek needed more indirects than expected");
    }
  }

  // Read all needed indirect blocks on the return buffers
  ret_indirect_blk_nums[0] = inode->indirect;
  for (i = 0; i < num_indirects; i++) {

    // Check that indirect block number exists
    if (!ret_indirect_blk_nums[i]) {
      ar_panic("AntFS: fs_lay_file_seek found NULL indirect");
    }

    // Read the indirect node
    if ((ret = fs_log_read_block(ret_indirect_blks + 4096 * i, 
                                 ret_indirect_blk_nums[i]))) {
      return ret;
    }

    // Write the next indirect block number
    if (i + 1 < FS_MAX_INDIRECTS) {
      ret_indirect_blk_nums[i+1] = 
          ((FSIndirect *) (ret_indirect_blks + 4096 * i))->indirect;
    }
  }

  // Get the data block number
  if (!num_indirects) {
    if (fs_mem_get_data_block_ptr(context, ret_inode_blk, 
                                  data_blk_index, &data_blk_num)) {
      ar_panic("AntFS: fs_lay_file_seek could not get inode data block ptr");
    }
  }
  else {
    if (fs_mem_get_data_block_ptr(context, 
                                  ret_indirect_blks + 4096 * (num_indirects-1),
                                  data_blk_rem, &data_blk_num)) {
      ar_panic("AntFS: fs_lay_file_seek could not get inode data block ptr");
    }
  }

  // Read the data block
  if ((ret = fs_log_read_block(ret_data_blk, data_blk_num))) {
    return ret;
  }

  // Return stuff
  if (ret_data_blk_num) {
    *ret_data_blk_num = data_blk_num;
  }
  if (ret_inode_blk_num) {
    *ret_inode_blk_num = inode_blk_num;
  }
  if (ret_num_indirects) {
    *ret_num_indirects = num_indirects;
  }
  if (ret_data_blk_index) {
    *ret_data_blk_index = data_blk_index;
  }

  return 0;
}


// ===========================================================================
// fs_lay_file_read()           Reads any number of bytes from a file. This 
//                              function expects that the file will be already 
//                              analyzed up to the start of the wanted 
//                              position. The inode, the correct first data 
//                              block and all the indirect node buffers up to 
//                              there must be given.
// ===========================================================================
// * INPUTS
//   unsigned int num_bytes     Bytes to read from the file
//   unsigned char *inode_blk   The inode block itself
//   unsigned int inode_blk_num The inode block number
//   int data_blk_index         Index of the last data block of file
//   int data_blk_offset        Offset in the last data block to start reading
//   unsigned int data_blk_num  The last data block number
//   int num_indirects          The number of actually used indirect blocks
//                              up to the currently seeked position.
//
// * INOUTS
//   unsigned char *data_blk    The last data block itself. If more data
//                              blocks are needed, this will be overwritten
//                              by this function.
//   unsigned char              Caller array of buffers that hold the
//      *indirect_blks          encountered indirect nodes. It must be big
//                              enough to hold the maximum file size in bytes,
//                              i.e. FS_MAX_FILE_SIZE * FS_MAX_INDIRECTS. Only
//                              num_indirects * 4096 bytes must be filled with
//                              actual indirects; the rest may be garbage,
//                              but still need to exist. If more data blocks
//                              are read by the function, some of them may
//                              be filled (*ret_num_indirects will be valid
//                              in total).
//   unsigned int               Caller array that holds the file indirect
//      *indirect_blk_nums      nodes block numbers. See above for its size.
//   Context *context           The current context
//
//
// * OUTPUTS
//   unsigned char *ret_buf     Buffer of data read from the file
//   unsigned int               Actual bytes read from the file. This may
//        *ret_num_bytes        be < num_bytes if we encountered end-of-file. 
//   int *ret_num_indirects     The new number of indirect blocks. This may
//                              be > num_indirects if more data blocks had
//                              to be read which required new indirects.
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_file_read(Context *context, unsigned int num_bytes, 
                     unsigned char *inode_blk, 
                     unsigned int inode_blk_num, int data_blk_index, 
                     int data_blk_offset, unsigned char *data_blk, 
                     unsigned int data_blk_num, unsigned char *indirect_blks, 
                     unsigned int *indirect_blk_nums, int num_indirects,
                     unsigned char *ret_buf, unsigned int *ret_num_bytes,
                     int *ret_num_indirects) {

  int bytes_to_read;
  int bytes_read;
  int have_block;
  FSInode *inode;
  FSData *data;
  FSIndirect *indirect;
  int qnt;
  int ind_index;
  int ind_rem;
  int ret;

  // Sanity checks
  if (!num_bytes) {
    ar_panic("AntFS: fs_lay_file_read called to read 0 bytes");
  }
  if (!inode_blk) {
    ar_panic("AntFS: fs_lay_file_read called with no inode_blk");
  }
  if ((inode_blk_num < 1) || (inode_blk_num >= context->fs_num_blocks - 1)) {
    ar_panic("AntFS: fs_lay_file_read called with invalid inode block num");
  }
  if ((data_blk_index < 0) || 
      (data_blk_index >= FS_INODE_NUM_DATA_BLOCKS + 
                         FS_MAX_INDIRECTS * FS_IND_NUM_DATA_BLOCKS)) {
    ar_panic("AntFS: fs_lay_file_read called with invalid data block index");
  }
  if ((data_blk_offset < 0) || (data_blk_offset >= FS_DATA_BLOCK_PAYLOAD)) {
    ar_panic("AntFS: fs_lay_file_read called with invalid data block offset");
  }
  if (!data_blk) {
    ar_panic("AntFS: fs_lay_file_read called with no data_blk");
  }
  if ((data_blk_num < 1) || (data_blk_num >= context->fs_num_blocks - 1)) {
    ar_panic("AntFS: fs_lay_file_read called with invalid data block num");
  }
  if (!indirect_blks) {
    ar_panic("AntFS: fs_lay_file_read called with no indirect_blks");
  }
  if (!indirect_blk_nums) {
    ar_panic("AntFS: fs_lay_file_read called with no indirect_blk_nums");
  }
  if (num_indirects < 0) {
    ar_panic("AntFS: fs_lay_file_read called with negative num_indirects");
  }
  if (!ret_buf) {
    ar_panic("AntFS: fs_lay_file_read called with no return buf");
  }

  // Initialize
  inode = (FSInode *) inode_blk;
  bytes_read = 0;
  have_block = 1;

  // Compute max number of bytes that we can read from the file
  bytes_to_read = inode->size - 
                  (data_blk_index * FS_DATA_BLOCK_PAYLOAD) - 
                  data_blk_offset;

  // Adapt it to what is requested
  if (num_bytes < bytes_to_read) {
    bytes_to_read = num_bytes;
  }

  // Loop until we read everything
  while (bytes_read < bytes_to_read) {

    // Do we need to read another data block?
    if (!have_block) {
      data_blk_offset = 0;
      data_blk_index++;
    
      // Do we need indirects?
      if (data_blk_index < FS_INODE_NUM_DATA_BLOCKS) {

        // Read data block number directly from the inode
        if (fs_mem_get_data_block_ptr(context, inode, 
                                      data_blk_index, &data_blk_num)) {
          ar_panic(
              "AntFS: fs_lay_file_read: could not get data prt from inode");
        }
      }
      else {

        // Which indirect is for us?
        ar_int_divide(data_blk_index - FS_INODE_NUM_DATA_BLOCKS,
                      FS_IND_NUM_DATA_BLOCKS, &ind_index, &ind_rem);

        // Do we have so many indirects already?
        if (ind_index < num_indirects) {
          indirect = (FSIndirect *) (indirect_blks + 4096 * ind_index);
        }
        // Read a new indirect. Caller arrays will have space for
        // as many as needed for the maximum file size.
        else if (ind_index == num_indirects) {

          // Sanity check
          if (num_indirects + 1 > FS_MAX_INDIRECTS) {
            ar_panic(
                "AntFS: fs_lay_file_read needed more indirects than expected");
          }
          
          // Check that the next indirect block number exists
          if (!ind_index) {
            indirect_blk_nums[0] = inode->indirect;
          }
          if (!indirect_blk_nums[ind_index]) {
            ar_panic("AntFS: fs_lay_file_read found NULL indirect");
          }

          // Read the new indirect node to caller array
          if ((ret = fs_log_read_block(indirect_blks + 4096 * ind_index, 
                                       indirect_blk_nums[ind_index]))) {
            goto finished;
          }

          // Write the next indirect block number to caller array
          indirect = (FSIndirect *) (indirect_blks + 4096 * ind_index);
          if (ind_index + 1 < FS_MAX_INDIRECTS) {
            indirect_blk_nums[ind_index + 1] = indirect->indirect;
          }

          // One more indirect was just read
          num_indirects++;
        }
        else {
          indirect = NULL;
          ar_panic("AntFS: fs_lay_file_read: missing indirect block");
        }

        // Read data block number from the indirect block
        if (fs_mem_get_data_block_ptr(context, indirect,
                                      ind_rem, &data_blk_num)) {
          ar_panic(
              "AntFS: fs_lay_file_read could not get data ptr from indirect");
        }
      }
    
      // Now that we finally got the data block number, read the data block
      if ((ret = fs_log_read_block(data_blk, data_blk_num))) {
        goto finished;
      }
    
    }

    // Check that the block is valid
    data = (FSData *) data_blk;
    if (!(data->header & FS_HDR_TYPE_DATA)) {
      ar_panic("AntFS: fs_lay_file_read encountered non-data block");
    }

    // How much will we read from this block?
    qnt = FS_DATA_BLOCK_PAYLOAD - data_blk_offset;
    if (qnt > (bytes_to_read - bytes_read)) {
      qnt = bytes_to_read - bytes_read;
    }

    // Copy to return buffer
    kt_memcpy(ret_buf + bytes_read, data->data + data_blk_offset, qnt);
    bytes_read += qnt;

    // This block is used up
    have_block = 0;
  }


  // If we reached this, return success
  ret = 0;

finished:

  // Return stuff
  if (ret_num_bytes) {
    *ret_num_bytes = bytes_read;
  }
  if (ret_num_indirects) {
    *ret_num_indirects = num_indirects;
  }

  return ret;
}


// ===========================================================================
// fs_lay_file_append()         Appends any number of bytes at the end
//                              of a file. This function expects that the
//                              file will be already analyzed: the inode,
//                              the last data block and all the indirect
//                              node buffers must be given.
//
//                              This function writes on the log:
//
//                              1. The modified, partially filled last 
//                                 data block of the file (if possible)
//                              2. As many new data blocks are needed
//                              3. All indirect nodes (modified old ones
//                                 plus new ones), last to first
//                              4. The modified inode block
// ===========================================================================
// * INPUTS
//   unsigned char *buf         Buffer of data to be appended to file
//   int num_buf_bytes          Bytes in the data buffer above
//   unsigned char *inode_blk   The inode block itself
//   unsigned int inode_blk_num The inode block number
//   int data_blk_index         Index of the last data block of file, or 
//                              -1 to indicate that the file didn't have
//                              any data blocks (0 length)
//   unsigned char *data_blk    The last data block itself, or a garbage
//                              block buffer that will be overwritten if
//                              the file has no data blocks
//   unsigned int data_blk_num  The last data block number, if it exists
//   unsigned char              Caller array of buffers that hold the
//      *indirect_blks          encountered indirect nodes. It must be big
//                              enough to hold the maximum file size in bytes,
//                              i.e. FS_MAX_FILE_SIZE * FS_MAX_INDIRECTS. Only
//                              num_indirects * 4096 bytes must be filled with
//                              actual indirects; the rest may be garbage,
//                              but still need to exist.
//   unsigned int               Caller array that holds the file indirect
//      *indirect_blk_nums      nodes block numbers. See above for its size.
//   int num_indirects          The number of actually used indirect blocks
//                              for the current file size. 
//   unsigned int atomic_flags  ATS/ATE flags to be used for the transaction
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_file_append(Context *context, unsigned char *buf, int num_buf_bytes,
                       unsigned char *inode_blk, unsigned int inode_blk_num, 
                       int data_blk_index, unsigned char *data_blk, 
                       unsigned int data_blk_num, unsigned char *indirect_blks,
                       unsigned int *indirect_blk_nums, int num_indirects,
                       unsigned int atomic_flags) {
  
  unsigned char *buf_pos;
  int bytes_written;
  FSInode *inode;
  FSData *data;
  FSIndirect *indirect;
  int rem;
  int have_data_block;
  int ret;
  int ind_index;
  int ind_rem;

  // Sanity checks
  if (!buf) {
    ar_panic("AntFS: fs_lay_file_append called with no buf");
  }
  if (!num_buf_bytes) {
    ar_panic("AntFS: fs_lay_file_append called to append 0 bytes");
  }
  if (!inode_blk) {
    ar_panic("AntFS: fs_lay_file_append called with no inode_blk");
  }
  if ((inode_blk_num < 1) || (inode_blk_num >= context->fs_num_blocks - 1)) {
    ar_panic("AntFS: fs_lay_file_append called with invalid inode block num");
  }
  if (!data_blk) {
    ar_panic("AntFS: fs_lay_file_append called with no data_blk");
  }
  if ((data_blk_index > -1) && 
      ((data_blk_num < 1) || (data_blk_num >= context->fs_num_blocks - 1))) {
    ar_panic("AntFS: fs_lay_file_append called with invalid data block num");
  }
  if (!indirect_blks) {
    ar_panic("AntFS: fs_lay_file_append called with no indirect_blks");
  }
  if (!indirect_blk_nums) {
    ar_panic("AntFS: fs_lay_file_append called with no indirect_blk_nums");
  }
  if (num_indirects < 0) {
    ar_panic("AntFS: fs_lay_file_append called with negative num_indirects");
  }

  // Initialize
  bytes_written = 0;
  inode = (FSInode *) inode_blk;
  data = (FSData *) data_blk;
  have_data_block = 0;
  buf_pos = buf;

  // Check that file size is acceptable (otherwise our indirect array will
  // be too small)
  if (inode->size + num_buf_bytes > FS_MAX_FILE_SIZE) {
    return -EOVERFLOW;
  }

  // No matter if we partially fill last data block or we create a new one,
  // its header should be data type and set the start of our atomic flags.
  data->header = FS_HDR_TYPE_DATA;
  if (atomic_flags & FS_HDR_ATOMIC_START) {
    data->header |= FS_HDR_ATOMIC_START;
  }

  // Can we start by partially filling last block?
  if (data_blk_index > -1) {
    ar_int_divide(inode->size, FS_DATA_BLOCK_PAYLOAD, NULL, &rem);
    if (rem) {

      // Yes, we can. Compute how much we can fill.
      bytes_written = FS_DATA_BLOCK_PAYLOAD - rem;
      if (bytes_written > num_buf_bytes) {
        bytes_written = num_buf_bytes;
      }
      
      // Fill last part of last data block
      kt_memcpy(data->data + rem, buf_pos, bytes_written);
      buf_pos += bytes_written;

      // Last block is modified, write it to the disk and free the old block
      if (fs_mem_set_free_block(context, data_blk_num)) {
        ar_panic("AntFS: fs_lay_file_append could not free data_blk_num");
      }
      data->free_block = data_blk_num;
      if ((ret = fs_log_append(context, data_blk, &data_blk_num))) {
        return ret;
      }

      // Mark that we have a data block to be linked with inode/indirects,
      // so that the loop below will start working with it 
      have_data_block = 1;
    }
  }

  // Loop while we have a data block already to link, or if we don't but
  // we still have bytes to write
  while ((have_data_block) || (bytes_written < num_buf_bytes)) {

    // If we don't have a data block, create a new one
    if (!have_data_block) {

      // It's a new block, with a new index (if it was -1, it gets 0)
      data_blk_index++;

      // Overwrite caller data block buffer with as much as possible
      rem = num_buf_bytes - bytes_written;
      if (rem < FS_DATA_BLOCK_PAYLOAD) {
        kt_memcpy(data->data, buf_pos, rem);
        buf_pos += rem;
        bytes_written += rem;
        
        // Zero out remainder of block
        data->free_block = 0;
        kt_memset(data->data + rem, 0, FS_DATA_BLOCK_PAYLOAD - rem);
      }
      else {
        kt_memcpy(data->data, buf_pos, FS_DATA_BLOCK_PAYLOAD);
        buf_pos += FS_DATA_BLOCK_PAYLOAD;
        bytes_written += FS_DATA_BLOCK_PAYLOAD;
      }

      // Write new block to disk
      if ((ret = fs_log_append(context, data_blk, &data_blk_num))) {
        return ret;
      }

      // Clear atomic flags for next loop
      data->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
    }

    // Now link the new block to the inode or indirect block it belongs to.
    // This happens only in the memory buffers. All indirects and inode will
    // be written later to disk.
    if (data_blk_index < FS_INODE_NUM_DATA_BLOCKS) {

      // Link it directly on the inode
      if (fs_mem_set_data_block_ptr(context, inode, 
                                    data_blk_index, data_blk_num)) {
        ar_panic("AntFS: fs_lay_file_append: could not set data prt in inode");
      }
    }
    else {

      // Which indirect is for us?
      ar_int_divide(data_blk_index - FS_INODE_NUM_DATA_BLOCKS,
                    FS_IND_NUM_DATA_BLOCKS, &ind_index, &ind_rem);

      // Do we have so many indirects already?
      if (ind_index == num_indirects) {

        // Sanity check
        if (num_indirects + 1 > FS_MAX_INDIRECTS) {
          ar_panic(
              "AntFS: fs_lay_file_append needed more indirects than expected");
        }

        // Create a new indirect. Caller array will have space for
        // as many as needed for the maximum file size.
        kt_memset(indirect_blks + 4096 * ind_index, 0, 4096);
        indirect_blk_nums[ind_index] = 0;
        indirect = (FSIndirect *) (indirect_blks + 4096 * ind_index);
        indirect->header = FS_HDR_TYPE_INDIRECT;
        num_indirects++;
      }
      else if (ind_index > num_indirects) {
        ar_panic("AntFS: fs_lay_file_append: missing indirect block");
      }

      // Link data with its indirect block
      if (fs_mem_set_data_block_ptr(context, indirect_blks + 4096 * ind_index,
                                    ind_rem, data_blk_num)) {
        ar_panic(
            "AntFS: fs_lay_file_append: could not set data prt in indirect");
      }
    }

    // Start next loop with a new data block
    have_data_block = 0;
  }

  
  // Now write to disk all indirect blocks, linking them in a chain from
  // last to first. We always need to write all the indirects, because
  // modifying the last of them (or adding a new last) needs a change of
  // pointer to its previous, which means a change of pointer to its
  // previous, etc.
  for (ind_index = num_indirects - 1; ind_index >= 0; ind_index--) {

    // Link it with next indirect block (if it exists)
    indirect = (FSIndirect *) (indirect_blks + 4096 * ind_index);
    if (ind_index < num_indirects - 1) {
      indirect->indirect = indirect_blk_nums[ind_index + 1];
    }
    else {
      indirect->indirect = 0;
    }

    // Free up old block, if it existed. New indirect blocks are marked
    // above with blk num 0.
    if (indirect_blk_nums[ind_index]) {
      if (fs_mem_set_free_block(context, indirect_blk_nums[ind_index])) {
        ar_panic("AntFS: fs_lay_file_append could not free old indirect");
      }
      indirect->free_block = indirect_blk_nums[ind_index];
    }
    else {
      indirect->free_block = 0;
    }

    // Clear any atomic flags
    indirect->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);

    // Write it on disk
    if ((ret = fs_log_append(context, indirect_blks + 4096 * ind_index, 
                             &(indirect_blk_nums[ind_index])))) {
      return ret;
    }
  }

  
  // Finally, modify the inode as needed and write it back to disk.
  inode->size += num_buf_bytes;
  inode->mtime++;
  if (num_indirects) {
    inode->indirect = indirect_blk_nums[0];
  }
  else {
    inode->indirect = 0;
  }
  
  // Clear old atomic flags and set our own
  inode->header &= ~(FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END);
  if (atomic_flags & FS_HDR_ATOMIC_END) {
    inode->header |= FS_HDR_ATOMIC_END;
  }
 
  // Free old inode block
  if (fs_mem_set_free_block(context, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_file_append could not free inode_blk_num");
  }
  inode->free_block0 = inode_blk_num;
  inode->free_block1 = 0;

  // Write the modified one to disk
  if ((ret = fs_log_append(context, inode_blk, &inode_blk_num))) {
    return ret;
  }

  // Update inode map with new inode block number
  if (fs_mem_update_inode_map(context, inode->inode_num, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_file_append could not set inode in inode map");
  }

  return 0;
}


// ===========================================================================
// fs_lay_delete_inode()        Deletes a file or directory inode, along with
//                              all its data blocks and indirect blocks. The 
//                              inode can be a file or directory, but only
//                              empty directories can be deleted.
//
//                              This function frees blocks and a single inode
//                              from the memory structures. It writes to the
//                              on-disk log one or more Delete blocks, which
//                              contain the block numbers of all the blocks
//                              that were freed, along with the inode number
//                              that was freed.
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned int inode_num     Inode number to be deleted
//   int is_file                1: inode should be a file
//                              0: inode should be a directory
//   unsigned int atomic_flags  ATS/ATE flags to be used for the transaction
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_lay_delete_inode(Context *context, unsigned int inode_num, int is_file,
                        unsigned int atomic_flags) {

  unsigned int inode_blk_num;
  unsigned char inode_blk[4096];
  unsigned int indirect_blk_num;
  unsigned char indirect_blk[4096];
  unsigned int delete_blk_num;
  unsigned char delete_blk[4096];
  FSInode *inode;
  FSIndirect *indirect;
  FSDelete *delete;
  int delete_counter;
  unsigned int data_blk_num;
  int ret;
  int i, j;


  // Sanity checks
  if (inode_num >= FS_MAX_INODES) {
    ar_panic("AntFS: fs_lay_delete_inode called with invalid inode num");
  }

  // Get inode block number
  if (fs_mem_read_inode_map(context, inode_num, &inode_blk_num)) {
    ar_panic("AntFS: fs_lay_delete_inode could not get inode block num");
  }

  // Read the inode
  if ((ret = fs_log_read_block(inode_blk, inode_blk_num))) {
    return ret;
  }
  inode = (FSInode *) inode_blk;


  // Check it's an inode
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic("AntFS: fs_lay_delete_inode found inconsistent inode num");
  }

  // If a directory, check it's empty and sane
  if (inode->attr & FS_ATTR_DIRECTORY) {
    // Are we trying to delete a direcory?
    if (is_file) {
      return -EISDIR;
    }
    // Is it  empty?
    if (inode->size) {
      return -ENOTEMPTY;
    }
    // Sanity-check this with empty dirdata blocks and no indirect
    if (fs_mem_get_data_block_ptr(context, inode, 0, &data_blk_num)) {
      ar_panic("AntFS: fs_lay_delete_inode could not get data ptr from inode");
    }
    if (data_blk_num) {
      ar_panic(
      "AntFS: fs_lay_delete_inode found dir with size 0 but existing dirdata");
    }
    if (inode->indirect) {
      ar_panic("AntFS: fs_lay_delete_inode found dir with indirect block");
    }
  }
  else {
    // Are we trying to delete a file?
    if (!is_file) {
      return -ENOTDIR;
    }
  }

  // Prepare the initial delete block
  kt_memset(delete_blk, 0, 4096);
  delete = (FSDelete *) delete_blk;
  delete->header = FS_HDR_TYPE_DELETE;
  delete_counter = 0;

  // We're deleting this inode
  delete->free_inode = inode_num;

  // Set atomic start, if we have it. We don't know yet if it will be the 
  // last block, so we don't set end here.
  if (atomic_flags & FS_HDR_ATOMIC_START) {
    delete->header |= FS_HDR_ATOMIC_START;
  }


  // Delete all data blocks from the inode
  for (i = 0; i < FS_INODE_NUM_DATA_BLOCKS; i++) {

    // Get next data block number
    if (fs_mem_get_data_block_ptr(context, inode, i, &data_blk_num)) {
      ar_panic("AntFS: fs_lay_delete_inode could not get data ptr from inode");
    }

    // Did we reached the end of the inode data blocks?
    if (!data_blk_num) {
      break;
    }

    // Free this data block
    if (fs_mem_set_free_block(context, data_blk_num)) {
      ar_panic("AntFS: fs_lay_delete_inode could not free data block");
    }

    // Add entry to the delete block
    if (fs_mem_set_data_block_ptr(context, delete_blk, delete_counter++, 
                                  data_blk_num)) {
      ar_panic("AntFS: fs_lay_delete_inode could not add ptr to delete block");
    }

    // Write delete block if it overflowed and start a new one
    if (delete_counter >= FS_DELETE_NUM_BLOCKS) {

      // Write it
      if ((ret = fs_log_append(context, delete_blk, &delete_blk_num))) {
        return ret;
      }

      // This block is immediately free after the next checkpoint
      if (fs_mem_set_free_block(context, delete_blk_num)) {
        ar_panic("AntFS: fs_lay_delete_inode could not free delete block");
      }

      // Start a new one on the same buffer. No atomic flags, no inode num.
      kt_memset(delete_blk, 0, 4096);
      delete->header = FS_HDR_TYPE_DELETE;
      delete_counter = 0;
    }
  }

  // For all the indirect blocks
  indirect_blk_num = inode->indirect;
  for (i = 0; (indirect_blk_num) && (i < FS_MAX_INDIRECTS); i++) {

    // Read the indirect block
    if ((ret = fs_log_read_block(indirect_blk, indirect_blk_num))) {
      return ret;
    }
    indirect = (FSIndirect *) indirect_blk;

    // Check it's an indrect
    if (!(indirect->header & FS_HDR_TYPE_INDIRECT)) {
      ar_panic(
          "AntFS: fs_lay_delete_inode found inconsistent indirect blk num");
    }

    // Delete all data blocks from the indirect
    for (j = 0; j < FS_IND_NUM_DATA_BLOCKS; j++) {

      // Get next data block number
      if (fs_mem_get_data_block_ptr(context, indirect, j, &data_blk_num)) {
        ar_panic(
            "AntFS: fs_lay_delete_inode could not get data ptr from indirect");
      }

      // Did we reached the end of the indirect data blocks?
      if (!data_blk_num) {
        break;
      }

      // Free this data block
      if (fs_mem_set_free_block(context, data_blk_num)) {
        ar_panic("AntFS: fs_lay_delete_inode could not free data block");
      }
      
      // Add entry to the delete block
      if (fs_mem_set_data_block_ptr(context, delete_blk, delete_counter++, 
                                    data_blk_num)) {
        ar_panic(
            "AntFS: fs_lay_delete_inode could not add ptr to delete block");
      }

      // Write delete block if it overflowed and start a new one
      if (delete_counter >= FS_DELETE_NUM_BLOCKS) {

        // Write it
        if ((ret = fs_log_append(context, delete_blk, &delete_blk_num))) {
          return ret;
        }

        // This block is immediately free after the next checkpoint
        if (fs_mem_set_free_block(context, delete_blk_num)) {
          ar_panic("AntFS: fs_lay_delete_inode could not free delete block");
        }

        // Start a new one on the same buffer. No atomic flags, no inode num.
        kt_memset(delete_blk, 0, 4096);
        delete->header = FS_HDR_TYPE_DELETE;
        delete_counter = 0;
      }
    }

    // Free the indirect itself
    if (fs_mem_set_free_block(context, indirect_blk_num)) {
      ar_panic("AntFS: fs_lay_delete_inode could not free indirect block");
    }

    // Add entry to the delete block
    if (fs_mem_set_data_block_ptr(context, delete_blk, delete_counter++, 
                                  indirect_blk_num)) {
      ar_panic(
          "AntFS: fs_lay_delete_inode could not add ptr to delete block");
    }

    // Write delete block if it overflowed and start a new one
    if (delete_counter >= FS_DELETE_NUM_BLOCKS) {

      // Write it
      if ((ret = fs_log_append(context, delete_blk, &delete_blk_num))) {
        return ret;
      }

      // This block is immediately free after the next checkpoint
      if (fs_mem_set_free_block(context, delete_blk_num)) {
        ar_panic("AntFS: fs_lay_delete_inode could not free delete block");
      }

      // Start a new one on the same buffer. No atomic flags, no inode num.
      kt_memset(delete_blk, 0, 4096);
      delete->header = FS_HDR_TYPE_DELETE;
      delete_counter = 0;
    }

    // Go to next indirect
    indirect_blk_num = indirect->indirect;
    if ((i == FS_MAX_INDIRECTS - 1) && (indirect_blk_num)) {
      ar_panic("AntFS: fs_lay_delete_inode found more indirects than allowed");
    }
  }

  // Free the inode itself
  if (fs_mem_set_free_block(context, inode_blk_num)) {
    ar_panic("AntFS: fs_lay_delete_inode could not free inode block");
  }

  // Add entry to the delete block
  if (fs_mem_set_data_block_ptr(context, delete_blk, delete_counter++, 
                                inode_blk_num)) {
    ar_panic("AntFS: fs_lay_delete_inode could not add ptr to delete block");
  }

  // This is the last (or the only) delete block. Set atomic end, if we have it.
  if (atomic_flags & FS_HDR_ATOMIC_END) {
    delete->header |= FS_HDR_ATOMIC_END;
  }

  // Write it
  if ((ret = fs_log_append(context, delete_blk, &delete_blk_num))) {
    return ret;
  }

  // This block is immediately free after the next checkpoint
  if (fs_mem_set_free_block(context, delete_blk_num)) {
    ar_panic("AntFS: fs_lay_delete_inode could not free delete block");
  }

  // Free the inode in the inode map
  if (fs_mem_update_inode_map(context, inode_num, 0)) {
    ar_panic("AntFS: fs_lay_delete_inode could not update inode map");
  }

  // Success
  return 0;
}

