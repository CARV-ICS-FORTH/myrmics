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
// Abstract      : AntFS: a CompactFlash filesystem
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: filesystem.h,v $
// CVS revision  : $Revision: 1.16 $
// Last modified : $Date: 2012/03/15 15:31:47 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _FILESYSTEM_H
#define _FILESYSTEM_H

#include <memory_management.h>
#include <types.h>

// Virtual file descriptors definitions
#define FS_MAX_DESCRIPTORS              16

#define FS_VFS_FLAG_CREAT               (1 << 0)
#define FS_VFS_FLAG_APPEND              (1 << 1)
#define FS_VFS_FLAG_TRUNC               (1 << 2)
#define FS_VFS_FLAG_RDONLY              (1 << 3)
#define FS_VFS_FLAG_WRONLY              (1 << 4)

#define FS_VFS_SEEK_SET                 1
#define FS_VFS_SEEK_CUR                 2
#define FS_VFS_SEEK_END                 3

// AntFS standard definitions. Note that changing any of those generally 
// makes existing filesystems on CompactFlash incompatible.
#define FS_MAX_INODES                   10000              // 10000 total files
#define FS_MAX_FILE_SIZE                (64 * 1024 * 1024) // 64MB max filesize
#define FS_MAX_DIRENTRY_NAME            64                 // 63 chars filename
#define FS_MAX_DIRENTRIES               60

#define FS_CHK_NUM_DATA_BLOCKS          1356
#define FS_IND_NUM_DATA_BLOCKS          1357
#define FS_INODE_NUM_DATA_BLOCKS        1345
#define FS_DELETE_NUM_BLOCKS            1358

#define FS_DATA_BLOCK_PAYLOAD           4076

#define FS_CHK_FREE_BITMAP_BYTES        4076

#define FS_MAX_INDIRECTS  (\
               (FS_MAX_FILE_SIZE - \
                    (FS_INODE_NUM_DATA_BLOCKS * FS_DATA_BLOCK_PAYLOAD) \
                    + (FS_IND_NUM_DATA_BLOCKS * FS_DATA_BLOCK_PAYLOAD) - 1) / \
               (FS_IND_NUM_DATA_BLOCKS * FS_DATA_BLOCK_PAYLOAD) + 1)

#define FS_MAX_TRANSACTION_NEW_INODES   2
#define FS_MAX_TRANSACTION_FREE_INODES  2
#define FS_MAX_TRANSACTION_FREE_BLOCKS  2000 // should be max file size / 4K
#define FS_MAX_TRANSACTION_USED_BLOCKS  512

// Header flags, common for all types of blocks
#define FS_HDR_TYPE_INODE            (1 <<  0)
#define FS_HDR_TYPE_INDIRECT         (1 <<  1)
#define FS_HDR_TYPE_DATA             (1 <<  2)
#define FS_HDR_TYPE_DELETE           (1 <<  3)
#define FS_HDR_TYPE_CHK_ROOT         (1 <<  4)
#define FS_HDR_TYPE_CHK_INDIRECT     (1 <<  5)
#define FS_HDR_TYPE_CHK_FREE_BITMAP  (1 <<  6)

#define FS_HDR_ATOMIC_START          (1 <<  8)
#define FS_HDR_ATOMIC_END            (1 <<  9)


// FSInode attr bit flags
#define FS_ATTR_DIRECTORY               (1 <<  0)
#define FS_ATTR_FILE                    (1 <<  1)

#define FS_ATTR_U_R                     (1 <<  8)
#define FS_ATTR_U_W                     (1 <<  9)
#define FS_ATTR_U_X                     (1 << 10)

#define FS_ATTR_G_R                     (1 << 16)
#define FS_ATTR_G_W                     (1 << 17)
#define FS_ATTR_G_X                     (1 << 18)

#define FS_ATTR_O_R                     (1 << 24)
#define FS_ATTR_O_W                     (1 << 25)
#define FS_ATTR_O_X                     (1 << 26)


// ============================================================================
// The four low-level block types in AntFS
// ============================================================================

// Inode block structure [4096 bytes]
typedef struct {

  // Header [8 bytes]
  unsigned int  header;         // Header (bitmap of flags defined above)
  unsigned int  log_seq_id;     // Sequential ID of log entry
  
  // Metadata [32 bytes]
  unsigned int  inode_num;      // Inode number 

  unsigned int  ctime;          // Creation time
  unsigned int  mtime;          // Last modified time
  unsigned int  atime;          // Last accessed time

  unsigned int  size;           // File size (or #children for directories)
  unsigned int  uid;            // User id
  unsigned int  gid;            // Group id
  unsigned int  attr;           // Attributes (bitmap of FS_ATTR_* flags 
                                // defined above)

  // Data block pointers [4036 bytes]
  // This is an array of 1345 data block pointers, each 3 bytes.
  // 1345 * 3 = 4035, so last byte is unused.
  // Map this value to FS_INODE_NUM_DATA_BLOCKS
  unsigned char data_blocks[4036];

  // Pointer to Indirect block [4 bytes]. MS byte always 0.
  unsigned int indirect;

  // Block numbers freed by this inode log entry (0 means not used)
  unsigned int free_block0;
  unsigned int free_block1;

  // CRC-64 [8 bytes]
  unsigned char crc_64[8];

} FSInode;


// Indirect block structure [4096 bytes]. Also used for Checkpoint Indirect 
// blocks.
typedef struct {

  // Header [8 bytes]
  unsigned int  header;         // Header (bitmap of flags defined above)
  unsigned int  log_seq_id;     // Sequential ID of log entry
  
  // Data block pointers [4072 bytes].
  // This is an array of 1357 data block pointers, each 3 bytes.
  // 1357 * 3 = 4071, so last byte is unused.
  // Map this value to FS_IND_NUM_DATA_BLOCKS
  unsigned char data_blocks[4072];

  // Next pointer to Indirect block pointer [4 bytes]. MS byte always 0.
  unsigned int indirect;

  // Block number freed by this indirect log entry (0 means not used)
  unsigned int free_block;

  // CRC-64 [8 bytes]
  unsigned char crc_64[8];

} FSIndirect;


// Data block structure [4096 bytes]
typedef struct {

  // Header [8 bytes]
  unsigned int  header;         // Header (bitmap of flags defined above)
  unsigned int  log_seq_id;     // Sequential ID of log entry
  
  // Data [4076 bytes]
  // Map this value to FS_DATA_BLOCK_PAYLOAD
  unsigned char data[4076];

  // Block number freed by this indirect log entry (0 means not used)
  unsigned int free_block;

  // CRC-64 [8 bytes]
  unsigned char crc_64[8];

} FSData;


// Delete block structure [4096 bytes]
typedef struct {

  // Header [8 bytes]
  unsigned int  header;         // Header (bitmap of flags defined above)
  unsigned int  log_seq_id;     // Sequential ID of log entry
  
  // Free block pointers [4076 bytes]
  // This is an array of 1358 block pointers, each 3 bytes.
  // 1358 * 3 = 4074, so last two bytes are unused.
  // Map this value to FS_DELETE_NUM_BLOCKS
  unsigned char free_blocks[4076];

  // Inode number freed by this delete log entry (0 means not used)
  unsigned int free_inode;

  // CRC-64 [8 bytes]
  unsigned char crc_64[8];

} FSDelete;


// Checkpoint Free Bitmap structure [4096 bytes]
typedef struct {

  // Header [8 bytes]
  unsigned int  header;         // Header (bitmap of flags defined above)
  unsigned int  log_seq_id;     // Sequential ID of log entry
  
  // Free blocks bitmap [4076 bytes]
  // Map this value to FS_CHK_FREE_BITMAP_BYTES
  unsigned char bitmap[4076];

  // Pointer to next Checkpoint Free bitmap block [4 bytes]. MS byte always 0.
  unsigned int next;

  // CRC-64 [8 bytes]
  unsigned char crc_64[8];

} FSCheckpointFreeBitmap;


// Checkpoint Root block structure [4096 bytes]
typedef struct {

  // Header [8 bytes]
  unsigned int  header;         // Header (bitmap of flags defined above)
  unsigned int  log_seq_id;     // Sequential ID of log entry (also serves
                                // to distinguish newest of 2 checkpoints)
  
  // Necessary pointers [4 bytes]
  unsigned int  log_head;       // Log head block pointer

  // Inode map inode block pointers [4068 bytes]
  // This is an array of 1356 inode block pointers, each 3 bytes.
  // Map this value to FS_CHK_NUM_DATA_BLOCKS
  unsigned char data_blocks[4068];

  // Pointer to Checkpoint Indirect block [4 bytes]. MS byte always 0.
  unsigned int indirect;

  // Pointer to Checkpoint Free bitmap block [4 bytes]. MS byte always 0.
  unsigned int free_bitmap;

  // CRC-64 [8 bytes]
  unsigned char crc_64[8];

} FSCheckpointRoot;


// ============================================================================
// Directory-related structures
// ============================================================================

// Directory entry type
typedef struct {
  char          name[FS_MAX_DIRENTRY_NAME];
  unsigned int  inode_num;
} FSDirEntry;

// Directory data [4080 bytes, fits in the data of a Data block]
typedef struct {

  // Array of 60 directory entries
  FSDirEntry    dir_entry[FS_MAX_DIRENTRIES];

} FSDirData;


// ============================================================================
// Current state
// ============================================================================
typedef struct fs_current_state_struct {
  
  // Current Inode map in a live-checkpoint block collection
  FSCheckpointRoot *chk;         // Live checkpoint block
  int              chk_pos;      // Last checkpoint position: 
                                 // 0=first device block, 1=last device block
  FSIndirect       **chk_ind;    // Array of checkpoint indirect blocks
  int              num_chk_ind;  // Number of checkpoint indirect blocks

  // Free blocks status in a live-checkpoint block collection
  // (conservative version, consistent with last checkpoint)
  FSCheckpointFreeBitmap 
                   **chk_free;   // Array of checkpoint free blocks
  int              num_chk_free; // Number of checkpoint free blocks
  unsigned int  num_free_blocks; // How many bits in above bitmap are free
  
  // Free blocks status in a live-checkpoint block collection
  // (future version, consistent with reality)
  FSCheckpointFreeBitmap 
                   **chk_future_free;   // Array of future free blocks
                                        // (size is same as num_chk_free)
  unsigned int  num_future_free_blocks; // How many bits in the bitmap are free

  // Log stuff
  unsigned int  log_seq_id;      // Log sequential ID
  unsigned int  log_tail;        // Log Tail block number

} FSCurrentState;


// ============================================================================
// Virtual file system descriptors
// ============================================================================
typedef struct fs_fdesc_struct {
  char          valid;           // Valid bit
  char          physical;        // 0: descriptor refers to UART-mapped I/O
                                 // 1: descriptor refers to an AntFS file
  unsigned int  inode_num;       // Inode number (for physical files)
  unsigned int  size;            // Current file size
  unsigned int  pos;             // Current file seek position
  unsigned int  flags;           // open() bitmap of flags for this descriptor
                                 // (see FS_VFS_FLAG_* definitions above)
} FSDescriptor;



// ============================================================================
// AntFS routines
// ============================================================================


// CRC
extern void fs_crc64(unsigned char msg, 
                     unsigned int *crc_hi, unsigned int *crc_lo);

// In-memory / util
extern int fs_mem_set_free_block(Context *context, unsigned int blk_num);
extern int fs_mem_unset_free_block(Context *context, unsigned int blk_num);
extern int fs_mem_find_free_block(Context *context, unsigned int *ret_blk_num);
extern int fs_mem_free_block_status(Context *context, unsigned int blk_num,
                                    int future, int *ret_status);
extern void fs_mem_sync_free_blocks(Context *context);
extern void fs_mem_init_free_blocks(Context *context);

extern int fs_mem_read_inode_map(Context *context, unsigned int inode_num,
                                 unsigned int *ret_inode_blk_num);
extern int fs_mem_update_inode_map(Context *context, unsigned int inode_num,
                                   unsigned int inode_blk_num);
extern int fs_mem_find_free_inode(Context *context, 
                                  unsigned int *ret_inode_num);

extern int fs_mem_set_data_block_ptr(Context *context, void *block, int index,
                                     unsigned int new_ptr);
extern int fs_mem_get_data_block_ptr(Context *context, void *block, int index,
                                     unsigned int *ret_ptr);

extern int fs_mem_init(Context *context);


// Log
extern int fs_log_write_block(unsigned char *buf, unsigned int adr);
extern int fs_log_read_block(unsigned char *buf, unsigned int adr);
extern int fs_log_append(Context *context, unsigned char *blk,
                         unsigned int *ret_blk_num);

extern int fs_log_checkpoint(Context *context, int free_previous);
extern int fs_log_format(Context *context, unsigned int log_seed);
extern int fs_log_mount(Context *context);


// Layout
extern int fs_lay_create_inode(Context *context, unsigned int inode_num,
                               unsigned int attr, unsigned int *ret_blk_num,
                               unsigned int atomic_flags);
extern int fs_lay_delete_inode(Context *context, unsigned int inode_num,
                               int is_file, unsigned int atomic_flags);

extern int fs_lay_find_last_entry(Context *context, unsigned int inode_blk_num,
                                  unsigned char *ret_inode_blk,
                                  unsigned char *ret_dirdata_blk,
                                  unsigned int *ret_dirdata_blk_num,
                                  int *ret_dirdata_index,
                                  int *ret_direntry_index);
extern int fs_lay_find_entry(Context *context, unsigned int inode_blk_num,
                             char *name, unsigned char *ret_inode_blk,
                             unsigned char *ret_dirdata_blk,
                             unsigned int *ret_dirdata_blk_num,
                             int *ret_dirdata_index,
                             int *ret_direntry_index);
extern int fs_lay_add_dir_entry(Context *context, unsigned char *inode_blk,
                                unsigned int inode_blk_num,
                                unsigned char *last_dirdata_blk,
                                unsigned int last_dirdata_blk_num,
                                int last_dirdata_index,
                                int last_direntry_index,
                                char *new_name, unsigned int new_inode_num,
                                unsigned int atomic_flags);
extern int fs_lay_delete_dir_entry(Context *context, unsigned char *inode_blk,
                                   unsigned int inode_blk_num,
                                   unsigned char *dirdata_blk,
                                   unsigned int dirdata_blk_num,
                                   int dirdata_index, int direntry_index,
                                   unsigned char *last_dirdata_blk,
                                   unsigned int last_dirdata_blk_num,
                                   int last_dirdata_index, 
                                   int last_direntry_index,
                                   unsigned int atomic_flags);

extern int fs_lay_file_seek(Context *context, unsigned int inode_num,
                            int data_blk_index, unsigned char *ret_inode_blk,
                            unsigned int *ret_inode_blk_num,
                            unsigned int *ret_data_blk_index,
                            unsigned char *ret_data_blk,
                            unsigned int *ret_data_blk_num,
                            unsigned char *ret_indirect_blks,
                            unsigned int *ret_indirect_blk_nums,
                            int *ret_num_indirects);
extern int fs_lay_file_read(Context *context, unsigned int num_bytes, 
                            unsigned char *inode_blk,
                            unsigned int inode_blk_num, int data_blk_index,
                            int data_blk_offset, unsigned char *data_blk,
                            unsigned int data_blk_num, 
                            unsigned char *indirect_blks,
                            unsigned int *indirect_blk_nums, int num_indirects,
                            unsigned char *ret_buf, unsigned int *ret_num_bytes,
                            int *ret_num_indirects);
extern int fs_lay_file_append(Context *context, unsigned char *buf,
                              int num_buf_bytes, unsigned char *inode_blk,
                              unsigned int inode_blk_num, int data_blk_index,
                              unsigned char *data_blk,
                              unsigned int data_blk_num,
                              unsigned char *indirect_blks,
                              unsigned int *indirect_blk_nums,
                              int num_indirects, unsigned int atomic_flags);


// Abstract
extern int fs_abs_create_entity(Context *context, unsigned int parent_inode_num,
                                char *new_name, unsigned int new_inode_num,
                                unsigned int attr);
extern int fs_abs_delete_entity(Context *context, unsigned int parent_inode_num,
                                char *name, int is_file);

extern int fs_abs_file_read(Context *context, unsigned int inode_num,
                            unsigned int num_bytes, unsigned int seek_pos, 
                            unsigned char *ret_buf, 
                            unsigned int *ret_num_bytes);
extern int fs_abs_file_append(Context *context, unsigned int inode_num,
                              unsigned char *buf, int num_buf_bytes,
                              unsigned int *ret_new_size);

extern int fs_abs_search_dir(Context *context, unsigned int dir_inode_num, 
                             char *name, unsigned int *ret_inode_num, 
                             unsigned int *ret_inode_attr,
                             unsigned int *ret_size);


// Top-level
extern int fs_ls_minus_lr(Context *context, unsigned int inode_num, char *name);
extern void fs_test();


// ============================================================================
// VFS routines
// ============================================================================
extern int fs_vfs_open(Context *context, char *pathname, unsigned int flags, 
                       int *ret_fdesc);
extern int fs_vfs_unlink(Context *context, char *pathname);

extern int fs_vfs_seek(Context *context, int fdesc, long int offset,
                       unsigned int whence, long int *ret_offset);
extern int fs_vfs_read(Context *context, int fdesc, unsigned int num_bytes,
                       unsigned char *ret_buf, unsigned int *ret_read_bytes);
extern int fs_vfs_write(Context *context, int fdesc, unsigned char *buf, 
                        unsigned int num_bytes, unsigned int *ret_write_bytes);

extern int fs_vfs_mkdir(Context *context, char *pathname);
extern int fs_vfs_rmdir(Context *context, char *pathname);



#endif
