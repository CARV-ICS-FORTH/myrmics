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
// Abstract      : Virtual file system layer routines
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: vfs.c,v $
// CVS revision  : $Revision: 1.4 $
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
// fs_vfs_open()                Opens a file and assigns a new file descriptor
//                              to it.
//
//                              Depending on the flags given, it will either:
//
//                              (a) open an existing file to read and/or write
//                              (b) do (a) and also seek to the end of file
//                              (c) create a new file
//                              (d) do (a) but truncate it to 0 length
// ===========================================================================
// * INPUTS
//   char *pathname             Absolute pathname of file to be opened
//   unsigned int flags         Flags to be used. See FS_VFS_FLAG_* in
//                              include/filesystem.h
//
// * INOUTS
//   Context *context           The current context, whose file descriptor 
//                              table is updated
//
// * OUTPUTS
//   int *ret_fdesc             File fdesc of the opened file, in the range
//                              3 ... FS_MAX_DESCRIPTORS. Descriptors 0, 1 and 
//                              2 are reserved for stdin, stdout, stderr.
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_open(Context *context, char *pathname, unsigned int flags, 
                int *ret_fdesc) {

  char *len;
  char *beg;
  char *end;
  char buf[FS_MAX_DIRENTRY_NAME];
  unsigned int cur_inode_num;
  unsigned int nxt_inode_num;
  unsigned int nxt_inode_attr;
  unsigned int nxt_inode_size;
  int ret;
  int create;
  int truncate;
  int i;
  int fdesc;

  // Sanity check
  if (!ret_fdesc) {
    ar_panic("VFS: open called with no return fdesc");
  }


  // Pathname should begin from the root of the filesystem
  if (pathname[0] != '/') {
    ar_panic("VFS: open called with relative pathname");
  }

  // Mark the end of the pathname and sanity-check it
  len = pathname + kt_strlen(pathname);
  if (len - pathname < 2) {
    ar_panic("VFS: open called with no file");
  }
  if (*(len - 1) == '/') {
    ar_panic("VFS: open called with trailing slash");
  }

  // Initialize stuff, starting from / (inode number 0)
  cur_inode_num = 0;
  create = 0;
  truncate = 0;

  // For all path components
  beg = pathname;
  while (*beg) {

    // Get next component
    end = kt_strchr(beg + 1, '/');
    if (!end) {
      end = len;
    }
    if (end - beg - 1 >= FS_MAX_DIRENTRY_NAME) {
      return -ENAMETOOLONG;
    }
    kt_strncpy(buf, beg + 1, end - beg - 1);
    buf[end - beg - 1] = 0;

    // Search for it in the current hierarchy
    if ((ret = fs_abs_search_dir(context, cur_inode_num, buf, &nxt_inode_num,
                                 &nxt_inode_attr, &nxt_inode_size))) {
      // We proceed only if search failed because of a non-existent 
      // path component.
      if (ret != -ENOENT) {
        return ret;
      }
    }

    // If not found and it wasn't the last path component, sorry
    if (ret && (end != len)) {
      return -ENOENT;
    }

    // If not found and it's the last component 
    if (ret && (end == len)) {

      // Do we need to create it?
      if (flags & FS_VFS_FLAG_CREAT) {
        // Other flags incompatible with CREAT
        if (flags & (FS_VFS_FLAG_TRUNC)) {
          return -EINVAL;
        }
        create = 1;
        break;
      }
      else {
        return -ENOENT;
      }
    }

    // Ok, it exists. Is it the last component?
    if (end == len) {
      
      // Are they asking us to create it? Fail the call if so, although this
      // is non-POSIX behavior (better safe than sorry...)
      if (flags & FS_VFS_FLAG_CREAT) {
        return -EEXIST;
      }

      // Make sure it's a file
      if (nxt_inode_attr & FS_ATTR_DIRECTORY) {
        return -EISDIR;
      }

      // Do we need to truncate it?
      if (flags & FS_VFS_FLAG_TRUNC) {
        truncate = 1;
        break;
      }

      // Then it simply exists and we must open it for read and/or write.
      break;
    }

    // It exists and it's not the last component. Make sure it's a 
    // directory and move on.
    if (!(nxt_inode_attr & FS_ATTR_DIRECTORY)) {
      return -ENOTDIR;
    }
    cur_inode_num = nxt_inode_num;
    beg = end;
  }

  // Initialize fdesc search to invalid fdesc 0 (reserved for 
  // always opened stdin)
  fdesc = 0; 

  // Run all the physical fdesc numbers. We try to find a free one, but
  // we must also check if the file is already opened.
  for (i = 3; i < FS_MAX_DESCRIPTORS; i++) {
    if (context->fs_fdesc[i].valid) {
      if ((!create) && (context->fs_fdesc[i].inode_num == nxt_inode_num)) {
        return -EALREADY;
      }
    }
    else {
      if (!fdesc) {
        fdesc = i;
      }
    }
  }

  // Maybe all the file descriptors are busy?
  if (!fdesc) {
    return -ENFILE;
  }

  // Open the descriptor
  context->fs_fdesc[fdesc].valid = 1;
  context->fs_fdesc[fdesc].physical = 1;
  context->fs_fdesc[fdesc].pos = 0;
  context->fs_fdesc[fdesc].flags = flags;
  
  // Assign existing inode number and size...
  if (!create) {
    context->fs_fdesc[fdesc].inode_num = nxt_inode_num;
    context->fs_fdesc[fdesc].size = nxt_inode_size;
  }
  // ... or do the creation
  else {
    // Find a new inode number
    if (fs_mem_find_free_inode(context, 
                               &(context->fs_fdesc[fdesc].inode_num))) {
      context->fs_fdesc[fdesc].valid = 0;
      return -EMAXINODES;
    }

    // Create the file
    if ((ret = fs_abs_create_entity(context, cur_inode_num, buf, 
                                    context->fs_fdesc[fdesc].inode_num, 
                                    FS_ATTR_FILE | FS_ATTR_U_R | FS_ATTR_U_W))){
      context->fs_fdesc[fdesc].valid = 0;
      return ret;
    }

    // No size as of yet
    context->fs_fdesc[fdesc].size = 0;
  }

  // Do we need to truncate?
  if (truncate) {
    ar_panic("FIXME");
  }

  // Return file descriptor
  *ret_fdesc = fdesc;

  return 0;
}


// ===========================================================================
// fs_vfs_seek()                Sets the current position of an opened file
// ===========================================================================
// * INPUTS
//   int fdesc                  The file descriptor to be seeked. Must be 
//                              already opened.
//   long int offset            Offset value in bytes
//   unsigned int whence        One of:
//                              FS_VFS_SEEK_SET: new pos set to exactly offset
//                                               bytes
//                              FS_VFS_SEEK_CUR: new pos set to current pos
//                                               plus offset bytes
//                              FS_VFS_SEEK_END: new pos set to file size plus
//                                               offset bytes
//
// * OUTPUTS
//   unsigned char *ret_offset  Actual offset returned
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_seek(Context *context, int fdesc, long int offset,
                unsigned int whence, long int *ret_offset) {
  
  FSDescriptor *f;
  long int new_offset;

  // Sanity checks
  if ((fdesc < 0) || (fdesc >= FS_MAX_DESCRIPTORS)) {
    return -EBADF;
  }
  f = &(context->fs_fdesc[fdesc]);
  if (!f->valid) {
    return -EBADF;
  }
  if (!f->physical) {
    return -EPERM;
  }
  if (f->pos > f->size) {
    ar_panic("VFS: read found inconsistent seek pos and file size");
  }

  // Translate what offset means
  switch (whence) {
    case FS_VFS_SEEK_SET: new_offset = offset; break;
    case FS_VFS_SEEK_CUR: new_offset = f->pos + offset; break;
    case FS_VFS_SEEK_END: new_offset = f->size + offset; break;
    default: return -EINVAL;
  }

  // Check the new offset
  if ((new_offset < 0) || (new_offset > f->size)) {
    return -EINVAL;
  }

  // Set it and return it
  f->pos = new_offset;
  if (ret_offset) {
    *ret_offset = new_offset;
  }

  return 0;
}


// ===========================================================================
// fs_vfs_read()                Reads up to a specified amount of bytes from
//                              an opened file
// ===========================================================================
// * INPUTS
//   int fdesc                  The file descriptor to be read. Must be 
//                              already opened.
//   unsigned int num_bytes     Bytes to read
//
// * OUTPUTS
//   unsigned char *ret_buf     Buffer to hold the returned read bytes
//   unsigned int               Actual number of bytes read. Might be
//      *ret_read_bytes         less than num_bytes if end-of-file occured.
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_read(Context *context, int fdesc, unsigned int num_bytes,
                unsigned char *ret_buf, unsigned int *ret_read_bytes) {
  
  FSDescriptor *f;
  int ret;


  // Sanity checks
  if ((fdesc < 0) || (fdesc >= FS_MAX_DESCRIPTORS)) {
    return -EBADF;
  }
  f = &(context->fs_fdesc[fdesc]);
  if (!f->valid) {
    return -EBADF;
  }
  if (f->flags & FS_VFS_FLAG_WRONLY) {
    return -EPERM;
  }
  if (!ret_buf) {
    ar_panic("VFS: read called with no return buffer");
  }
  if (!ret_read_bytes) {
    ar_panic("VFS: read called with no return read_bytes");
  }
  if (f->pos > f->size) {
    ar_panic("VFS: read found inconsistent seek pos and file size");
  }

  // Differentiate physical files from stdin, stdout and stderr
  if (f->physical) {

    // Read the file
    *ret_read_bytes = 0;
    ret = fs_abs_file_read(context, f->inode_num, num_bytes, f->pos,
                           ret_buf, ret_read_bytes);
    
    // Advance seek position, even if an error occured
    f->pos += *ret_read_bytes;
    if (ret) {
      return ret;
    }
  }
  // Handle stdin reads
  else if (fdesc == 0) {
    ar_panic("FIXME read stdin");
  }
  // stdout, stderr only for writing
  else if ((fdesc == 1) || (fdesc == 2)) {
    return -EPERM;
  }
  else {
    ar_panic("VFS: read called with non-physical, valid descriptor > 3");
  }

  return 0;
}


// ===========================================================================
// fs_vfs_write()               Writes a specified amount of bytes to an
//                              opened file
// ===========================================================================
// * INPUTS
//   int fdesc                  The file descriptor to be written. Must be 
//                              already opened.
//   unsigned char buf          Buffer of data to write
//   unsigned int num_bytes     Bytes in the buffer to write
//
// * OUTPUTS
//   unsigned int               Actual number of bytes written. Might be
//      *ret_write_bytes        less than num_bytes if some error occurs.
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_write(Context *context, int fdesc, unsigned char *buf, 
                 unsigned int num_bytes, unsigned int *ret_write_bytes) {
  
  FSDescriptor *f;
  unsigned int old_size;
  int ret;


  // Sanity checks
  if ((fdesc < 0) || (fdesc >= FS_MAX_DESCRIPTORS)) {
    return -EBADF;
  }
  f = &(context->fs_fdesc[fdesc]);
  if (!f->valid) {
    return -EBADF;
  }
  if (f->flags & FS_VFS_FLAG_RDONLY) {
    return -EPERM;
  }
  if (!buf) {
    ar_panic("VFS: write called with no buffer");
  }
  if (!ret_write_bytes) {
    ar_panic("VFS: write called with no return write_bytes");
  }
  if (f->pos > f->size) {
    ar_panic("VFS: write found inconsistent seek pos and file size");
  }

  // Differentiate physical files from stdin, stdout and stderr
  if (f->physical) {

    // We only do appends for now
    if (!(f->flags & FS_VFS_FLAG_APPEND)) {
      return -ENOTSUP;
    }

    // Remember old size
    old_size = f->size;
    
    // Append to the file
    if ((ret = fs_abs_file_append(context, f->inode_num, buf, num_bytes,
                                  &(f->size)))) {
      return ret;
    }
    f->pos = f->size;

    // Return bytes written
    if (ret_write_bytes) {
      *ret_write_bytes = f->size - old_size;
    }
  }
  // Handle stdout writes
  else if (fdesc == 1) {
    ar_panic("FIXME write stdout");
  }
  else if (fdesc == 2) {
    ar_panic("FIXME write stderr");
  }
  // stdin only for reading
  else if (fdesc == 0) {
    return -EPERM;
  }
  else {
    ar_panic("VFS: write called with non-physical, valid descriptor > 3");
  }

  return 0;
}


// ===========================================================================
// fs_vfs_mkdir()               Creates a new directory
// ===========================================================================
// * INPUTS
//   char *pathname             Absolute pathname of the directory to be 
//                              created
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_mkdir(Context *context, char *pathname) {

  char *len;
  char *beg;
  char *end;
  char buf[FS_MAX_DIRENTRY_NAME];
  unsigned int cur_inode_num;
  unsigned int nxt_inode_num;
  unsigned int nxt_inode_attr;
  unsigned int nxt_inode_size;
  int ret;

  
  // Pathname should begin from the root of the filesystem
  if (pathname[0] != '/') {
    ar_panic("VFS: mkdir called with relative pathname");
  }

  // Mark the end of the pathname and sanity-check it
  len = pathname + kt_strlen(pathname);
  if (len - pathname < 2) {
    ar_panic("VFS: mkdir called with no file");
  }
  if (*(len - 1) == '/') {
    ar_panic("VFS: mkdir called with trailing slash");
  }

  // Initialize stuff, starting from / (inode number 0)
  cur_inode_num = 0;

  // For all path components
  beg = pathname;
  while (*beg) {

    // Get next component
    end = kt_strchr(beg + 1, '/');
    if (!end) {
      end = len;
    }
    if (end - beg - 1 >= FS_MAX_DIRENTRY_NAME) {
      return -ENAMETOOLONG;
    }
    kt_strncpy(buf, beg + 1, end - beg - 1);
    buf[end - beg - 1] = 0;

    // Search for it in the current hierarchy
    if ((ret = fs_abs_search_dir(context, cur_inode_num, buf, &nxt_inode_num,
                                 &nxt_inode_attr, &nxt_inode_size))) {
      // We proceed only if search failed because of a non-existent 
      // path component.
      if (ret != -ENOENT) {
        return ret;
      }
    }

    // If not found and it wasn't the last path component, sorry
    if (ret && (end != len)) {
      return -ENOENT;
    }

    // If not found and it's the last component, it's time to create it.
    if (ret && (end == len)) {
      break;
    }

    // Ok, it exists. Is it the last component?
    if (end == len) {
      return -EEXIST;
    }

    // It exists and it's not the last component. Make sure it's a 
    // directory and move on.
    if (!(nxt_inode_attr & FS_ATTR_DIRECTORY)) {
      return -ENOTDIR;
    }
    cur_inode_num = nxt_inode_num;
    beg = end;
  }

  // Find a new inode number
  if (fs_mem_find_free_inode(context, &nxt_inode_num)) {
    return -EMAXINODES;
  }

  // Create the directory
  if ((ret = fs_abs_create_entity(context, cur_inode_num, buf, 
                                  nxt_inode_num, FS_ATTR_DIRECTORY | 
                                  FS_ATTR_U_R | FS_ATTR_U_W | FS_ATTR_U_X))){
    return ret;
  }


  return 0;
}


// ===========================================================================
// fs_vfs_rmdir()               Deletes an empty directory
// ===========================================================================
// * INPUTS
//   char *pathname             Absolute pathname of the directory to be 
//                              deleted
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_rmdir(Context *context, char *pathname) {

  char *len;
  char *beg;
  char *end;
  char buf[FS_MAX_DIRENTRY_NAME];
  unsigned int cur_inode_num;
  unsigned int nxt_inode_num;
  unsigned int nxt_inode_attr;
  unsigned int nxt_inode_size;
  int ret;

  
  // Pathname should begin from the root of the filesystem
  if (pathname[0] != '/') {
    ar_panic("VFS: rmdir called with relative pathname");
  }

  // Mark the end of the pathname and sanity-check it
  len = pathname + kt_strlen(pathname);
  if (len - pathname < 2) {
    ar_panic("VFS: rmdir called with no file");
  }
  if (*(len - 1) == '/') {
    ar_panic("VFS: rmdir called with trailing slash");
  }

  // Initialize stuff, starting from / (inode number 0)
  cur_inode_num = 0;

  // For all path components
  beg = pathname;
  while (*beg) {

    // Get next component
    end = kt_strchr(beg + 1, '/');
    if (!end) {
      end = len;
    }
    if (end - beg - 1 >= FS_MAX_DIRENTRY_NAME) {
      return -ENAMETOOLONG;
    }
    kt_strncpy(buf, beg + 1, end - beg - 1);
    buf[end - beg - 1] = 0;

    // If we found the last component, break
    if (end == len) {
      break;
    }

    // Search component it in the current hierarchy
    if ((ret = fs_abs_search_dir(context, cur_inode_num, buf, &nxt_inode_num,
                                 &nxt_inode_attr, &nxt_inode_size))) {
      return ret;
    }

    // Make sure it's a directory and move on.
    if (!(nxt_inode_attr & FS_ATTR_DIRECTORY)) {
      return -ENOTDIR;
    }
    cur_inode_num = nxt_inode_num;
    beg = end;
  }

  // Delete the directory. It will be checked for emptiness.
  if ((ret = fs_abs_delete_entity(context, cur_inode_num, buf, 0))) {
    return ret;
  }


  return 0;
}


// ===========================================================================
// fs_vfs_unlink()              Deletes a file from a directory. No hard
//                              links are supported, so this can only be
//                              used on closed files which are going to
//                              be deleted, unlike the POSIX unlink call.
// ===========================================================================
// * INPUTS
//   char *pathname             Absolute pathname of the file to be deleted
//
// * INOUTS
//   Context *context           The current context
//
// * RETURN VALUE
//   int                        0 for success
//                              < 0: see include/errno.h
// ===========================================================================
int fs_vfs_unlink(Context *context, char *pathname) {

  char *len;
  char *beg;
  char *end;
  char buf[FS_MAX_DIRENTRY_NAME];
  unsigned int cur_inode_num;
  unsigned int nxt_inode_num;
  unsigned int nxt_inode_attr;
  unsigned int nxt_inode_size;
  int ret;
  int i;

  
  // Pathname should begin from the root of the filesystem
  if (pathname[0] != '/') {
    ar_panic("VFS: unlink called with relative pathname");
  }

  // Mark the end of the pathname and sanity-check it
  len = pathname + kt_strlen(pathname);
  if (len - pathname < 2) {
    ar_panic("VFS: unlink called with no file");
  }
  if (*(len - 1) == '/') {
    ar_panic("VFS: unlink called with trailing slash");
  }

  // Initialize stuff, starting from / (inode number 0)
  cur_inode_num = 0;

  // For all path components
  beg = pathname;
  while (*beg) {

    // Get next component
    end = kt_strchr(beg + 1, '/');
    if (!end) {
      end = len;
    }
    if (end - beg - 1 >= FS_MAX_DIRENTRY_NAME) {
      return -ENAMETOOLONG;
    }
    kt_strncpy(buf, beg + 1, end - beg - 1);
    buf[end - beg - 1] = 0;

    // Search for it in the current hierarchy
    if ((ret = fs_abs_search_dir(context, cur_inode_num, buf, &nxt_inode_num,
                                 &nxt_inode_attr, &nxt_inode_size))) {
      return ret;
    }

    // Ok, it exists. Is it the last component?
    if (end == len) {
      
      // Make sure it's a file
      if (nxt_inode_attr & FS_ATTR_DIRECTORY) {
        return -EISDIR;
      }

      // Exit the loop
      break;
    }

    // It exists and it's not the last component. Make sure it's a 
    // directory and move on.
    if (!(nxt_inode_attr & FS_ATTR_DIRECTORY)) {
      return -ENOTDIR;
    }
    cur_inode_num = nxt_inode_num;
    beg = end;
  }

  // Run all the physical fdesc numbers to see if it's already opened
  for (i = 3; i < FS_MAX_DESCRIPTORS; i++) {
    if (context->fs_fdesc[i].valid) {
      if (context->fs_fdesc[i].inode_num == nxt_inode_num) {
        return -EALREADY;
      }
    }
  }

  // Delete the file
  // FIXME: This call will rescan the directory to find the file.
  //        This is already done here, because we needed to check if the
  //        file is opened. Need to remove this redundancy, by separating
  //        fs_abs_delete_entity to a file and a directory function.
  if ((ret = fs_abs_delete_entity(context, cur_inode_num, buf, 1))) {
    return ret;
  }

  return 0;
}


