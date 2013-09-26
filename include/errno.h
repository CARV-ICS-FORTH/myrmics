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
// Abstract      : Error codes
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: errno.h,v $
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2013/03/15 18:09:00 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#ifndef _ERRNO_H
#define _ERRNO_H


#define EIO           1  // I/O error
#define ECRC          2  // CRC error
#define ENOMEM        3  // Memory allocation error
#define EINVAL        4  // Invalid argument error
#define ENOSPC        5  // No space left on filesystem
#define ENOENT        6  // No such file or directory
#define EOVERFLOW     7  // Too many directory entries / seek pos too big
#define ENAMETOOLONG  8  // Filename too long
#define EEXIST        9  // File already exists
#define EISDIR       10  // Pathname is a directory
#define ENFILE       11  // Maximum number of files opened
#define ENOTDIR      12  // Path component is not a directory
#define EACCES       13  // Access denied
#define EALREADY     14  // File already open / cannot unlink opened file
#define EMAXINODES   15  // Maximum number of inodes reached
#define EBADF        16  // Bad file descriptor
#define EPERM        17  // Permission denied
#define ENOTSUP      18  // Operation not supported
#define ENOTEMPTY    19  // Directory not empty

#endif
