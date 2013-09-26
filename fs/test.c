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
// Abstract      : AntFS test routines
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: test.c,v $
// CVS revision  : $Revision: 1.17 $
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


    
int fs_ls_minus_lr(Context *context, unsigned int inode_num, char *name) {
  unsigned int inode_blk_num;
  unsigned char buf[4096];
  unsigned char buf2[4096];
  char buf3[1024];
  FSInode *inode;
  FSData *dirdata;
  FSDirData *dirdata_entries;
  unsigned int dirdata_blk_num;
  int ret;
  int i, j;

  // Get block number
  if (fs_mem_read_inode_map(context, inode_num, &inode_blk_num)) {
    ar_panic("AntFS: fs_ls_inode could not get inode block number");
  }

  // Sanity check
  if ((inode_blk_num < 1) || (inode_blk_num >= context->fs_num_blocks)) {
    ar_panic("AntFS: fs_ls_inode called with invalid inode blk num");
  }

  // Read inode
  if ((ret = fs_log_read_block(buf, inode_blk_num))) {
    return ret;
  }
  
  inode = (FSInode *) buf;

  // Print file or directory
  kt_printf("%5u %c%c%c%c%c%c%c%c%c%c root  root  %8u %3d %s\r\n",
      inode_num, 
      (inode->attr & FS_ATTR_DIRECTORY) ? 'd' : '-',
      (inode->attr & FS_ATTR_U_R)       ? 'r' : '-',
      (inode->attr & FS_ATTR_U_W)       ? 'w' : '-',
      (inode->attr & FS_ATTR_U_X)       ? 'x' : '-',
      (inode->attr & FS_ATTR_G_R)       ? 'r' : '-',
      (inode->attr & FS_ATTR_G_W)       ? 'w' : '-',
      (inode->attr & FS_ATTR_G_X)       ? 'x' : '-',
      (inode->attr & FS_ATTR_O_R)       ? 'r' : '-',
      (inode->attr & FS_ATTR_O_W)       ? 'w' : '-',
      (inode->attr & FS_ATTR_O_X)       ? 'x' : '-',
      inode->size,
      inode->mtime,
      (kt_strcmp(name, "")) ? name : "/");
  
  // Is it a file?
  if (inode->attr & FS_ATTR_FILE) {
    return 0;
  }

  // If it's a directory, recurse on its children
  for (i = 0; i < FS_INODE_NUM_DATA_BLOCKS; i++) {

    // Get the next block number
    if (fs_mem_get_data_block_ptr(context, buf, i, &dirdata_blk_num)) {
      ar_panic("AntFS: ls_minus_lr could not get data block pointer");
    }

    // Are the blocks finished?
    if (!dirdata_blk_num) {
      break;
    }

    // Read this dirdata block
    if ((ret = fs_log_read_block(buf2, dirdata_blk_num))) {
      return ret;
    }

    // Sanity check
    dirdata = (FSData *) buf2;
    if (!(dirdata->header & FS_HDR_TYPE_DATA)) {
      ar_panic("AntFS: ls_minus_lr found inode pointer to to non-data block");
    }
    dirdata_entries = (FSDirData *) dirdata->data;

    // For all the directory entries of this dirdata block
    for (j = 0; j < FS_MAX_DIRENTRIES; j++) {
      
      // Are the direntries finished?
      if (!dirdata_entries->dir_entry[j].inode_num) {
        break;
      }

      // List this entry recursively
      if (kt_strlen(dirdata_entries->dir_entry[j].name) + 
          kt_strlen(name) > 1020) {
        ar_panic("ls_minus_lr static name buffer overflow");
      }
      kt_sprintf(buf3, "%s/%s", name, dirdata_entries->dir_entry[j].name);

      if ((ret = fs_ls_minus_lr(context, 
                                dirdata_entries->dir_entry[j].inode_num, 
                                buf3))) {
        return ret;
      }
    }
  }


  return 0;
}



void fs_print_block(unsigned char *buf) {
  int i, j;

  for (i = 0; i < 128; i++) {
    kt_printf("%3d. ", i);
    for (j = 0; j < 32; j++) {
      kt_printf("%02X", buf[32*i + j]);
    }
    kt_printf("\r\n");
  }
  kt_printf("\r\n");
}


void fs_peek_all_blocks() {
  Context *context;
  unsigned int blk;
  int status;
  unsigned int header;
  unsigned int log_seq_id;
  unsigned char buf[4096];
  int ret;
  FSInode *inode;
  FSIndirect *indirect;
  FSCheckpointRoot *checkpoint;
  FSCheckpointFreeBitmap *bitmap;
  FSDelete *delete;
  FSData *data;
  unsigned int ptr;
  int i;

  kt_printf("\r\nFilesystem status:\r\n");

  context = mm_get_context(ar_get_core_id());

  for (blk = 0; blk < context->fs_num_blocks; blk++) {
    if (fs_mem_free_block_status(context, blk, 0, &status)) {
      ar_panic("peek: free block error");
    }
    if (!status) {

      if ((ret = fs_log_read_block(buf, blk))) {
        switch (ret) {
          case -EINVAL: kt_printf("Read 0x%05X: Invalid arg\r\n", blk); break;
          case -EIO: kt_printf("Read 0x%05X: I/O error\r\n", blk); break;
          case -ECRC: kt_printf("Read 0x%05X: CRC error\r\n", blk); break;
          default: kt_printf("Read 0x%05X: Unknown error\r\n", blk);
        }
        continue;
      }

      kt_memcpy(&header, buf, 4);

      // Show only Delete blocks from conservative version
      if (!(header & FS_HDR_TYPE_DELETE)) {
        if (fs_mem_free_block_status(context, blk, 1, &status)) {
          ar_panic("peek: future free block error");
        }
        if (status) {
          continue;
        }
      }

      kt_printf("  Block 0x%05X: ", blk);

      kt_memcpy(&log_seq_id, buf + 4, 4);
      kt_printf("LogSeqID = 0x%08X, ", log_seq_id);


      if (header & FS_HDR_ATOMIC_START) {
        kt_printf("ATS, ");
      }
      if (header & FS_HDR_ATOMIC_END) {
        kt_printf("ATE, ");
      }

      if (header & FS_HDR_TYPE_INODE) {
        kt_printf("Inode.\r\n");
        inode = (FSInode *) buf;
        kt_printf(
            "                 Inode num = %u, size = %u, indirect = 0x%05X\r\n",
            inode->inode_num, inode->size, inode->indirect);
        for (i = 0; i < FS_INODE_NUM_DATA_BLOCKS; i++) {
          if (fs_mem_get_data_block_ptr(context, buf, i, &ptr)) {
            ar_panic("Peek: could not get inode ptr");
          }
          kt_printf(
            "                 data blk ptr %d = 0x%05X\r\n", i, ptr);
          if (!ptr) {
            break;
          }
        }
        if (inode->free_block0) {
          kt_printf(
            "                 free blk ptr 0 = 0x%05X\r\n", inode->free_block0);
        }
        if (inode->free_block1) {
          kt_printf(
            "                 free blk ptr 1 = 0x%05X\r\n", inode->free_block1);
        }
      }
      if ((header & FS_HDR_TYPE_INDIRECT) || 
          (header & FS_HDR_TYPE_CHK_INDIRECT)) {
        if (header & FS_HDR_TYPE_INDIRECT) {
          kt_printf("Indirect.\r\n");
        }
        else {
          kt_printf("Checkpoint Indirect.\r\n");
        }
        indirect = (FSIndirect *) buf;
        kt_printf(
            "                 Indirect = 0x%05X\r\n", 
            indirect->indirect);
        for (i = 0; i < FS_IND_NUM_DATA_BLOCKS; i++) {
          if (fs_mem_get_data_block_ptr(context, buf, i, &ptr)) {
            ar_panic("Peek: could not get indirect ptr");
          }
          kt_printf(
            "                 data blk ptr %d = 0x%05X\r\n", i, ptr);
          if (!ptr) {
            break;
          }
        }
        if (indirect->free_block) {
          kt_printf(
            "                 free blk ptr = 0x%05X\r\n", indirect->free_block);
        }
      }
      if (header & FS_HDR_TYPE_CHK_FREE_BITMAP) {
        kt_printf("Checkpoint Free Bitmap.\r\n");
        bitmap = (FSCheckpointFreeBitmap *) buf;
        kt_printf(
            "                 Next = 0x%05X\r\n", 
            bitmap->next);
      }
      if (header & FS_HDR_TYPE_CHK_ROOT) {
        kt_printf("Checkpoint Root.\r\n");
        checkpoint = (FSCheckpointRoot *) buf;
        kt_printf(
"                 Log head = 0x%05X, indirect = 0x%05X, bitmap = 0x%05X\r\n",
            checkpoint->log_head, checkpoint->indirect, 
            checkpoint->free_bitmap);
        for (i = 0; i < FS_CHK_NUM_DATA_BLOCKS; i++) {
          if (fs_mem_get_data_block_ptr(context, buf, i, &ptr)) {
            ar_panic("Peek: could not get checkpoint ptr");
          }
          kt_printf(
            "                 data blk ptr %d = 0x%05X\r\n", i, ptr);
          if (!ptr) {
            break;
          }
        }
      }
      if (header & FS_HDR_TYPE_DATA) {
        kt_printf("Data.\r\n");
        data = (FSData *) buf;
        kt_printf("                 > [");
        for (i = 0; i < 20; i++) {
          kt_printf("%c", data->data[i]);
        }
        kt_printf("]\r\n");
        if (data->free_block) {
          kt_printf(
            "                 free blk ptr = 0x%05X\r\n", data->free_block);
        }
      }
      if (header & FS_HDR_TYPE_DELETE) {
        kt_printf("Delete.\r\n");
        delete = (FSDelete *) buf;
        for (i = 0; i < FS_DELETE_NUM_BLOCKS; i++) {
          if (fs_mem_get_data_block_ptr(context, buf, i, &ptr)) {
            ar_panic("Peek: could not get delete ptr");
          }
          kt_printf(
            "                 blk ptr %d = 0x%05X\r\n", i, ptr);
          if (!ptr) {
            break;
          }
        }
        if (delete->free_inode) {
          kt_printf(
            "                 free inode = %u\r\n", delete->free_inode);
        }
      }
    }
  }

  kt_printf("  >> Free blocks: %u (conservative), %u (real)\r\n\n", 
      context->fs_state->num_free_blocks,
      context->fs_state->num_future_free_blocks);
}



void fs_test() {
  int sectors;
  Context *context;
  //unsigned int new_inode_num;
  //char buf[4*4080];


  kt_printf("\r\n"
            "============================================================\r\n"
            "= Entering FileSystem tester\r\n"
            "============================================================\r\n");

  // Init & query device
  sectors = ar_flash_init_device();
  if (sectors <= 0) {
    ar_panic("Invalid CompactFlash status");
  }

  // We operate in 4,096 byte blocks. Update context.
  context = mm_get_context(ar_get_core_id());
  context->fs_num_blocks = context->fs_flash_num_sectors / 8;

  // Create the virtual file descriptors table
  if (!(context->fs_fdesc = kt_zalloc(FS_MAX_DESCRIPTORS * 
                                      sizeof(FSDescriptor)))) {
    goto out_of_memory;
  }

  // Create the 3 non-physical file descriptors
  context->fs_fdesc[0].valid = 1;
  context->fs_fdesc[0].flags = FS_VFS_FLAG_RDONLY; // stdin
  context->fs_fdesc[1].valid = 1;
  context->fs_fdesc[1].flags = FS_VFS_FLAG_WRONLY; // stdout
  context->fs_fdesc[2].valid = 1;
  context->fs_fdesc[2].flags = FS_VFS_FLAG_WRONLY; // stderr


  // Allocate filesystem state
  if (!(context->fs_state = kt_zalloc(sizeof(FSCurrentState)))) {
    goto out_of_memory;
  }

  // Initialize filesystem state
  if (fs_mem_init(context)) {
    ar_panic("Filesystem memory initialization error");
  }

  // Format. WARNING: All data is lost here!
  if (fs_log_format(context, 0xABCD1000)) {
    ar_panic("Filesystem format error");
  }
  fs_peek_all_blocks();


/*

  char buf[128];
  int i;

  for (i = 0; i < 80; i++) {
    kt_sprintf(buf, "/dir%02d", i);
    if (fs_vfs_mkdir(context, buf)) {
      ar_panic("mkdir returned an error");
    }
  }

  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  fs_peek_all_blocks();


  int ret;
  
  for (i = 0; i < 80; i++) {
    kt_sprintf(buf, "/dir%02d", i);
    kt_printf("Deleting %s\r\n", buf);
    if (fs_vfs_rmdir(context, buf)) {
      ar_panic("rmdir returned an error");
    }
  }
  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  fs_peek_all_blocks();
*/

  int fd;
  unsigned int bytes;
  int ret;
  char buf[4096];
  long offset;


  if ((ret = fs_vfs_open(context, "/file", 
                         FS_VFS_FLAG_CREAT | FS_VFS_FLAG_APPEND, 
                         &fd))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, fd = %d\r\n",fd);
  }

  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  //fs_peek_all_blocks();


  kt_sprintf(buf, "Hallo file!");
  kt_printf("Writing %d bytes\r\n", kt_strlen(buf));

  if ((ret = fs_vfs_write(context, fd, (unsigned char *) buf, 
                          kt_strlen(buf), &bytes))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, bytes written = %u\r\n", bytes);
  }

  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  //fs_peek_all_blocks();



  if ((ret = fs_vfs_seek(context, fd, 0, FS_VFS_SEEK_SET, &offset))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, offset = %ld\r\n", offset);
  }


  kt_memset(buf, 0, sizeof(buf));

  if ((ret = fs_vfs_read(context, fd, sizeof(buf), (unsigned char *) buf, 
                         &bytes))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, bytes read = %u\r\n", bytes);
    kt_printf("Data read: [%s]\r\n", buf);
  }

  kt_memset(buf, 0, sizeof(buf));

  kt_sprintf(buf, "append...");
  kt_printf("Appending %d bytes\r\n", kt_strlen(buf));

  if ((ret = fs_vfs_write(context, fd, (unsigned char *) buf, 
                          kt_strlen(buf), &bytes))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, bytes written = %u\r\n", bytes);
  }

  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  //fs_peek_all_blocks();



  if ((ret = fs_vfs_seek(context, fd, 0, FS_VFS_SEEK_SET, &offset))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, offset = %ld\r\n", offset);
  }
  
  
  kt_memset(buf, 0, sizeof(buf));

  if ((ret = fs_vfs_read(context, fd, sizeof(buf), (unsigned char *) buf, 
                         &bytes))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("Ok, bytes read = %u\r\n", bytes);
    kt_printf("Data read: [%s]\r\n", buf);
  }




context->fs_fdesc[3].valid = 0;

  if ((ret = fs_vfs_unlink(context, "/file"))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  else {
    kt_printf("ok, deleted\r\n");
  }
  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  fs_peek_all_blocks();

/*
  kt_printf("\nTaking checkpoint.\n\n");
  if ((ret = fs_log_checkpoint(context, 1))) {
    kt_printf("Error = %d\r\n",ret);
    ar_panic("");
  }
  fs_peek_all_blocks();
*/

  kt_printf("\r\n"
            "============================================================\r\n"
            "= Mount\r\n"
            "============================================================\r\n");


  // Re-init filesystem state
  if (fs_mem_init(context)) {
    ar_panic("Filesystem memory initialization error");
  }

  // mount
  if (fs_log_mount(context)) {
    ar_panic("Mount error");
  }
  if (fs_ls_minus_lr(context, 0, "")) {
    ar_panic("fs_ls_minus_lr returned an error");
  }
  fs_peek_all_blocks();


  ar_halt();



out_of_memory:  
  ar_panic("Out of memory");
}


