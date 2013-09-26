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
// Abstract      : Header file for the kernel-side implementation of 
//                 system calls
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: syscall.h,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/11/15 08:49:39 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _SYSCALL_H
#define _SYSCALL_H

#include <stdarg.h>

#include <types.h>

// Flags for _sys_spawn. Keep these in accordance to the app compiler
// (e.g. SCOOP) if one is used. Note they are stored as a character, so
// use values < 0xFF.
#define SYS_TYPE_IN_ARG            (1 << 0)     
#define SYS_TYPE_OUT_ARG           (1 << 1)       
#define SYS_TYPE_NO_TRANSFER_ARG   (1 << 2)     
#define SYS_TYPE_SAFE_ARG          (1 << 3)     
#define SYS_TYPE_REGION_ARG        (1 << 4)     
#define SYS_TYPE_DEP_PENDING       (1 << 7) // Internal use: dependence must
                                            // be still waited for

#define SYS_TYPE_INOUT_ARG         (SYS_TYPE_IN_ARG | SYS_TYPE_OUT_ARG)
#define SYS_TYPE_BYVALUE_ARG       (SYS_TYPE_SAFE_ARG | SYS_TYPE_IN_ARG)


// ===========================================================================
// Memory allocation system calls
// ===========================================================================
extern void *_sys_alloc(char *filename, int line_nr, size_t size, rid_t region);
extern void _sys_balloc(char *filename, int line_nr, size_t size, rid_t region, 
                        int num_objects, void **objects);
extern void _sys_free(char *filename, int line_nr, void *ptr);
extern void *_sys_realloc(char *filename, int line_nr, void *old_ptr, 
                          size_t new_size, rid_t new_region);

extern rid_t _sys_ralloc(char *filename, int line_nr, rid_t parent, 
                         int level_hint);
extern void _sys_rfree(char *filename, int line_nr, rid_t region);


// ===========================================================================
// Task-related system calls and task ID table
// ===========================================================================
extern void (**_sys_task_table)();

extern void _sys_spawn(char *filename, int line_nr, unsigned int idx, 
                       void **args, unsigned int *types, 
                       unsigned int num_args);

extern void _sys_wait_on(char *filename, int line_nr, void **ptrs, 
                         unsigned int num_ptrs);


// ===========================================================================
// Miscellaneous system calls & helper functions
// ===========================================================================
#if 0
extern void _sys_zerostats(char *filename, int line_nr);
extern void _sys_getstats(char *filename, int line_nr,
                          unsigned int *ret_sched_sec, 
                          unsigned int *ret_sched_usec,
                          unsigned int *ret_sched_calls, 
                          unsigned int *ret_worker_sec, 
                          unsigned int *ret_worker_usec, 
                          unsigned long *ret_worker_sent, 
                          unsigned long *ret_worker_recv,
                          unsigned int *ret_barrier_sec, 
                          unsigned int *ret_barrier_usec);
#endif

extern int _sys_warning(char *filename, int line_nr, char *format, ...);
extern void _sys_error(char *filename, int line_nr, char *format, ...);



#endif
