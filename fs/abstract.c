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
// Abstract      : AntFS functions that handle the abstract, logical
//                 filesystem aspects
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: abstract.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/01/27 12:30:15 $
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
// fs_abs_create_entity()       Creates a new, empty file or directory. It 
//                              creates a new inode for it and then goes
//                              into the parent directory inode and modifies
//                              it to contain the new inode.
//
//                              This function will append to the log:
//
//                              1. The new inode for the file/directory
//                              2. The (modified or new) data block of the
//                                 parent directory inode which contains the
//                                 new direntry
//                              3. The modified inode block of the parent
//                                 directory
// ===========================================================================
// * INPUTS
//   unsigned int parent_inode_num  The parent directory inode number
//   char *new_name                 The name of the new file or directory
//   unsigned int new_inode_num     The inode number of the new file/directory
//   unsigned int attr              Attributes of the new file/dir inode
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_abs_create_entity(Context *context, unsigned int parent_inode_num,
                         char *new_name, unsigned int new_inode_num,
                         unsigned int attr) {

  unsigned char parent_inode_blk[4096];
  unsigned int parent_inode_blk_num;
  unsigned char last_dirdata_blk[4096];
  unsigned int last_dirdata_blk_num;
  int last_dirdata_index;
  int last_direntry_index;
  int ret;

  // Sanity checks
  if ((!new_name) || (kt_strlen(new_name) >= FS_MAX_DIRENTRY_NAME)) {
    return -EINVAL;
  }

  // Get parent inode block number
  if (fs_mem_read_inode_map(context, parent_inode_num, 
                            &parent_inode_blk_num)) {
    ar_panic(
        "AntFS: fs_abs_create_entity could not get parent inode block num");
  }

  // Find last directory entry of parent dir
  if ((ret = fs_lay_find_last_entry(context, parent_inode_blk_num, 
                                parent_inode_blk, last_dirdata_blk, 
                                &last_dirdata_blk_num, &last_dirdata_index, 
                                &last_direntry_index))) {
    return ret;
  }

  // Create new inode
  if ((ret = fs_lay_create_inode(context, new_inode_num, attr,
                NULL, FS_HDR_ATOMIC_START))) {
    return ret;
  }

  // Modify parent inode to contain the new entry
  if ((ret = fs_lay_add_dir_entry(context, 
                              parent_inode_blk, parent_inode_blk_num,
                              last_dirdata_blk, last_dirdata_blk_num,
                              last_dirdata_index, last_direntry_index, 
                              new_name, new_inode_num,
                              FS_HDR_ATOMIC_END))) {
    return ret;
  }

  return 0;
}

                     
// ===========================================================================
// fs_abs_delete_entity()       Deletes a file or directory from a parent
//                              directory. It will delete its inode block
//                              and then update the parent directory entries
//                              accordingly. In the case of a directory 
//                              deletion, the victim directory needs to 
//                              be empty.
//
//                              This function will free the victim inode, and
//                              if it's a file it will also free all its
//                              data blocks and/or indirect blocks. It will
//                              append to the log:
//
//                              1. The modified data block of the parent entry
//                              2. The modified last data block of the parent
//                              3. The modified inode block of the parent
//                                 directory
//
//                              Note that this is a worst-case estimate. If
//                              parent directory can hold its entries in a 
//                              single block, only (2) and (3) will be written.
//                              If the parent directory is being emptied (i.e.
//                              its last entry is being deleted) only (3)
//                              will be appended.
// ===========================================================================
// * INPUTS
//   unsigned int parent_inode_num  The parent directory inode number
//   char *name                     The name of the file or directory to be
//                                  deleted
//   int is_file                    1: name should be a file
//                                  0: name should be a directory
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_abs_delete_entity(Context *context, unsigned int parent_inode_num,
                         char *name, int is_file) {

  unsigned char parent_inode_blk[4096];
  unsigned int parent_inode_blk_num;
  unsigned char dirdata_blk[4096];
  unsigned int dirdata_blk_num;
  int dirdata_index;
  int direntry_index;
  unsigned char last_dirdata_blk[4096];
  unsigned int last_dirdata_blk_num;
  int last_dirdata_index;
  int last_direntry_index;
  int ret;
  FSInode *inode;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  unsigned int child_inode_num;

  // Sanity checks
  if ((!name) || (kt_strlen(name) >= FS_MAX_DIRENTRY_NAME)) {
    return -EINVAL;
  }

  // Get parent inode block number
  if (fs_mem_read_inode_map(context, parent_inode_num, 
                            &parent_inode_blk_num)) {
    ar_panic(
        "AntFS: fs_abs_delete_entity could not get parent inode block num");
  }

  // Find the entry in the parent directory
  if ((ret = fs_lay_find_entry(context, parent_inode_blk_num, name, 
                               parent_inode_blk, dirdata_blk, &dirdata_blk_num,
                               &dirdata_index, &direntry_index))) {
    return ret;
  }

  // Sanity check the parent inode
  inode = (FSInode *) parent_inode_blk;
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic(
     "AntFS: fs_delete_entity got a non-inode return from fs_lay_find_entry");
  }

  // Get the child block number
  dirdata = (FSData *) dirdata_blk;
  dirdata_entries = (FSDirData *) dirdata->data;
  if ((direntry_index < 0) || (direntry_index >= FS_MAX_DIRENTRIES)) {
    ar_panic("AntFS: fs_delete_entity got an invalid direntry_index");
  }
  child_inode_num = dirdata_entries->dir_entry[direntry_index].inode_num;
  if ((child_inode_num < 1) || (child_inode_num >= FS_MAX_INODES)) {
    ar_panic("AntFS: fs_delete_entity found an invalid child inode num");
  }

  // Also find the last entry in the parent directory
  // FIXME: this re-reads the parent_inode, which was read above. Pass
  //        a switch to avoid that.
  if ((ret = fs_lay_find_last_entry(context, parent_inode_blk_num, 
                                parent_inode_blk, last_dirdata_blk, 
                                &last_dirdata_blk_num, &last_dirdata_index, 
                                &last_direntry_index))) {
    return ret;
  }

  // Delete the victim inode
  if ((ret = fs_lay_delete_inode(context, child_inode_num, is_file,
                                 FS_HDR_ATOMIC_START))) {
    return ret;
  }

  // Delete the entry from the parent directory
  if ((ret = fs_lay_delete_dir_entry(context, 
                                     parent_inode_blk, parent_inode_blk_num,
                                     dirdata_blk, dirdata_blk_num,
                                     dirdata_index, direntry_index, 
                                     last_dirdata_blk, last_dirdata_blk_num,
                                     last_dirdata_index, last_direntry_index, 
                                     FS_HDR_ATOMIC_END))) {
    return ret;
  }

  return 0;
}

                     
// ===========================================================================
// fs_abs_file_read()           Reads bytes from a file.
// ===========================================================================
// * INPUTS
//   Context *context           The current context
//   unsigned int inode_num     The inode number of the file to be read
//   unsigned int num_bytes     Bytes to be read from the file
//   unsigned int seek_pos      Byte offset from start of file to start 
//                              reading from
//
// * OUTPUTS
//   unsigned char *ret_buf     Buffer to write the data that was read
//   unsigned int               Actual number of bytes read
//      *ret_num_bytes         
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_abs_file_read(Context *context, unsigned int inode_num,
                     unsigned int num_bytes, unsigned int seek_pos,
                     unsigned char *ret_buf, unsigned int *ret_num_bytes) {

  unsigned char inode_blk[4096];
  unsigned char data_blk[4096];
  unsigned char indirect_blks[FS_MAX_INDIRECTS * 4096];
  unsigned int  indirect_blk_nums[FS_MAX_INDIRECTS];
  unsigned int inode_blk_num;
  unsigned int data_blk_index;
  unsigned int data_blk_offset;
  unsigned int data_blk_num;
  int num_indirects;
  int ret;

  // Make sure seek_pos is sane
  if (seek_pos > FS_MAX_FILE_SIZE) {
    return -EINVAL;
  }

  // Compute the data block that the seek position belongs
  ar_uint_divide(seek_pos, FS_DATA_BLOCK_PAYLOAD, 
                 &data_blk_index, &data_blk_offset);

  // Seek to this block
  if ((ret = fs_lay_file_seek(context, inode_num, data_blk_index,
                              inode_blk, &inode_blk_num,
                              NULL, data_blk, &data_blk_num,
                              indirect_blks, indirect_blk_nums,
                              &num_indirects))) {
    return ret;
  }

  // Read data from file
  if ((ret = fs_lay_file_read(context, num_bytes, inode_blk, inode_blk_num, 
                              data_blk_index, data_blk_offset, data_blk,
                              data_blk_num, indirect_blks, indirect_blk_nums, 
                              num_indirects, ret_buf, ret_num_bytes, NULL))) {
    return ret;
  }

  return 0;
}


// ===========================================================================
// fs_abs_file_append()         Appends a chunk of data to the end of a file.
//                              This function will append to the log:
//
//                              1. The modified, partially filled last 
//                                 data block of the file (if possible)
//                              2. As many new data blocks are needed
//                              3. All indirect nodes (modified old ones
//                                 plus new ones), last to first
//                              4. The modified inode block
// ===========================================================================
// * INPUTS
//   unsigned int inode_num     The inode number of the file to be appended
//   unsigned char *buf         Data to be appended to the file
//   int num_buf_bytes          Bytes to be written from the data buffer
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h for values
// ===========================================================================
int fs_abs_file_append(Context *context, unsigned int inode_num,
                       unsigned char *buf, int num_buf_bytes,
                       unsigned int *ret_new_size) {
  unsigned char inode_blk[4096];
  unsigned char data_blk[4096];
  unsigned char indirect_blks[FS_MAX_INDIRECTS * 4096];
  unsigned int  indirect_blk_nums[FS_MAX_INDIRECTS];
  unsigned int inode_blk_num;
  unsigned int data_blk_index;
  unsigned int data_blk_num;
  int num_indirects;
  int ret;
  FSInode *inode;

  // Seek to the end of file
  if ((ret = fs_lay_file_seek(context, inode_num, -1, 
                              inode_blk, &inode_blk_num,
                              &data_blk_index, data_blk, &data_blk_num,
                              indirect_blks, indirect_blk_nums,
                              &num_indirects))) {
    return ret;
  }

  // Append data to file
  if ((ret = fs_lay_file_append(context, buf, num_buf_bytes, inode_blk, 
                                inode_blk_num, data_blk_index, data_blk, 
                                data_blk_num, indirect_blks, indirect_blk_nums,
                                num_indirects, 
                                FS_HDR_ATOMIC_START | FS_HDR_ATOMIC_END))) {
    return ret;
  }

  // Return new size
  inode = (FSInode *) inode_blk;
  if (ret_new_size) {
    *ret_new_size = inode->size;
  }

  return 0;
}


// ===========================================================================
// fs_abs_search_dir()            Searches a directory to find a named entry
// ===========================================================================
// * INPUTS
//   Context *context             The current context
//   unsigned int dir_inode_num   The directory inode number
//   char *name                   The child name to search for
//
// * OUTPUTS
//   unsigned int *ret_inode_num  Child inode number, if found
//   unsigned int *ret_inode_attr Child inode attributes, if found
//   unsigned int *ret_size       Child inode size, if found
//
// * RETURN VALUE
//   int                          0 for success
//                                -ENOENT if child is not found
//                                other errors according to include/errno.h
// ===========================================================================
int fs_abs_search_dir(Context *context, unsigned int dir_inode_num, char *name,
                      unsigned int *ret_inode_num, 
                      unsigned int *ret_inode_attr,
                      unsigned int *ret_size) {
  int ret;
  unsigned int dir_inode_blk_num;
  unsigned char inode_blk[4096];
  unsigned char dirdata_blk[4096];
  int direntry_index;
  FSInode *inode;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  unsigned int child_inode_num;
  unsigned int child_inode_blk_num;

  // Get directory inode block number
  if (fs_mem_read_inode_map(context, dir_inode_num, &dir_inode_blk_num)) {
    ar_panic(
        "AntFS: fs_abs_search_dir could not get directory inode block num");
  }

  // Search for the entry
  if ((ret = fs_lay_find_entry(context, dir_inode_blk_num, name, inode_blk,
                               dirdata_blk, NULL, NULL, &direntry_index))) {
    return ret;
  }

  // Sanity check the parent inode
  inode = (FSInode *) inode_blk;
  if (!(inode->header & FS_HDR_TYPE_INODE)) {
    ar_panic(
     "AntFS: fs_abs_search_dir got a non-inode return from fs_lay_find_entry");
  }

  // Get the child block number
  dirdata = (FSData *) dirdata_blk;
  dirdata_entries = (FSDirData *) dirdata->data;
  if ((direntry_index < 0) || (direntry_index >= FS_MAX_DIRENTRIES)) {
    ar_panic("AntFS: fs_abs_search_dir got an invalid direntry_index");
  }
  child_inode_num = dirdata_entries->dir_entry[direntry_index].inode_num;
  if ((child_inode_num < 1) || (child_inode_num >= FS_MAX_INODES)) {
    ar_panic("AntFS: fs_abs_search_dir found an invalid child inode num");
  }

  // Do we need to read the child inode?
  if (ret_inode_attr || ret_size) {
    
    // Get child inode block number
    if (fs_mem_read_inode_map(context, child_inode_num, &child_inode_blk_num)) {
      ar_panic(
          "AntFS: fs_abs_search_dir could not get child inode block num");
    }

    // Read the child inode
    if ((ret = fs_log_read_block(inode_blk, child_inode_blk_num))) {
      return ret;
    }

    // Sanity check the child inode
    inode = (FSInode *) inode_blk;
    if (!(inode->header & FS_HDR_TYPE_INODE)) {
      ar_panic("AntFS: fs_abs_search_dir found a non-inode child");
    }
  }


  // Return the related fields
  if (ret_inode_num) {
    *ret_inode_num = child_inode_num;
  }
  if (ret_inode_attr) {
    *ret_inode_attr = inode->attr;
  }
  if (ret_size) {
    *ret_size = inode->size;
  }

  return 0;
}
