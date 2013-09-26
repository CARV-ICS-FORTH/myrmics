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
// Author        : Spyros Lyberis
// Abstract      : Myrmics header file for application code. Defines the
//                 system calls that can be used and provides some common
//                 utility functions as well.
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: myrmics.h,v $
// CVS revision  : $Revision: 1.11 $
// Last modified : $Date: 2013/02/22 16:34:51 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#ifndef _MYRMICS_H
#define _MYRMICS_H

#include <stdarg.h>

#include <syscall.h>


// ===========================================================================
// System calls are just linked together for the moment (no user/kernel space
// separation yet). They are provided through macros that prepend the source
// code filename and line number
// ===========================================================================

// Memory allocation system calls
#define  sys_alloc(size, region) \
        _sys_alloc(__FILE__, __LINE__, (size), (region))

#define  sys_balloc(size, region, num_objects, objects) \
        _sys_balloc(__FILE__, __LINE__, (size), (region), (num_objects), \
                    (objects))

#define  sys_free(ptr) \
        _sys_free(__FILE__, __LINE__, (ptr))

#define  sys_realloc(old_ptr, new_size, new_region) \
        _sys_realloc(__FILE__, __LINE__, (old_ptr), (new_size), (new_region))

#define  sys_ralloc(parent, level_hint) \
        _sys_ralloc(__FILE__, __LINE__, (parent), (level_hint))

#define  sys_rfree(region) \
        _sys_rfree(__FILE__, __LINE__, (region))


// Miscellaneous system calls
#if 0
#define  sys_zerostats() \
        _sys_zerostats(__FILE__, __LINE__)

#define  sys_getstats(ret_sched_sec, ret_sched_usec, ret_sched_calls, \
                      ret_worker_sec, ret_worker_usec, ret_worker_sent,  \
                      ret_worker_recv, ret_barrier_sec, ret_barrier_usec) \
        _sys_getstats(__FILE__, __LINE__, (ret_sched_sec), (ret_sched_usec), \
                      (ret_sched_calls), (ret_worker_sec), (ret_worker_usec), \
                      (ret_worker_sent), (ret_worker_recv), (ret_barrier_sec), \
                      (ret_barrier_usec))
#endif

// ===========================================================================
// Compiler-generated calls, included here if anyone uses the system without
// the compiler and writes all calls by hand.
// ===========================================================================
#define  sys_spawn(idx, args, types, num_args) \
        _sys_spawn(__FILE__, __LINE__, (idx), (args), (types), (num_args))

#define  sys_wait_on(ptrs, num_ptrs) \
        _sys_wait_on(__FILE__, __LINE__, (ptrs), (num_ptrs))


// ===========================================================================
// Common utility functions are just mapped to the kernel toolset ones
// (no user/kernel separation yet).
// ===========================================================================
#define vsprintf                        kt_vsprintf
#define vprintf                         kt_vprintf
#define printf                          kt_printf
#define sprintf                         kt_sprintf
#define memcpy                          kt_memcpy
#define memset                          kt_memset
#define memmove                         kt_memmove

extern int kt_vsprintf(char *buf, const char *fmt, va_list args);
extern int kt_vprintf(const char *format, va_list ap);
extern int kt_printf(const char *format, ...);
extern int kt_sprintf(char *buf, const char *format, ...);
extern void *kt_memcpy(void *dst, const void *src, int num_bytes);
extern void *kt_memset(void *buf, char val, int num_bytes);
extern void *kt_memmove(void *dst, const void *src, int num_bytes);

static inline unsigned int rand(unsigned int seed) {
  return (1103515245 * seed + 12345);
}

static inline int int_log2(int n) {
  int tmp1 = n;
  int tmp2 = 0;

  while (tmp1 >>= 1) {
    tmp2++;
  }

  return tmp2;
}

static inline int int_log10(int n) {
  int tmp1 = n;
  int tmp2 = 0;

  while (tmp1 /= 10) {
    tmp2++;
  }

  return tmp2;
}



// ===========================================================================
// Other useful functions we want to expose from the architectural layer
// ===========================================================================
#define sys_timer_reset                 ar_timer_reset
#define sys_timer_get_cycles            ar_timer_get_cycles
#define sys_timer_get_msec              ar_timer_get_msec
#define sys_timer_busy_wait_cycles      ar_timer_busy_wait_cycles
#define sys_timer_busy_wait_msec        ar_timer_busy_wait_msec
#define sys_free_timer_get_ticks        ar_free_timer_get_ticks
#define sys_free_timer_get_cycles       ar_free_timer_get_cycles

#define sqrt                            ar_float_sqrt
#define pow                             ar_float_pow

extern void         ar_timer_reset();
extern unsigned int ar_timer_get_cycles();
extern unsigned int ar_timer_get_msec();
extern void         ar_timer_busy_wait_cycles(unsigned int cycles);
extern void         ar_timer_busy_wait_msec(unsigned int msec);
extern unsigned int ar_free_timer_get_ticks();
extern unsigned int ar_free_timer_get_cycles();

extern float        ar_float_sqrt(float num);
extern float        ar_float_pow(float x, float y);

extern void __attribute__ ((noreturn)) ar_panic(char *msg);

#define sys_assert(c) {\
                       if (!(c)) {\
                         char _ar_panic_buf[1024]; \
                         kt_sprintf(_ar_panic_buf, \
                                    "Assertion failed: %s line %d", \
                                    __FILE__, __LINE__); \
                         ar_panic(_ar_panic_buf); \
                       }\
                     }

#define sys_abort() {\
                     char _ar_panic_buf[1024]; \
                     kt_sprintf(_ar_panic_buf, \
                                "Internal error: %s line %d", \
                                __FILE__, __LINE__); \
                         ar_panic(_ar_panic_buf); \
                     }



// ===========================================================================
// Other useful functions we want to expose from the processing layer
// ===========================================================================
#define sys_get_worker_id               pr_get_worker_id
#define sys_get_num_workers             pr_get_num_workers
#define sys_get_num_schedulers          pr_get_num_schedulers

extern int pr_get_worker_id();
extern int pr_get_num_workers();
extern int pr_get_num_schedulers();


#endif
