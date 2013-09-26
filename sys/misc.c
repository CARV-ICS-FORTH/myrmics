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
// Abstract      : Miscellaneous helper routines for system calls
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: misc.c,v $
// CVS revision  : $Revision: 1.3 $
// Last modified : $Date: 2013/01/13 17:03:07 $
// Last author   : $Author: lyberis-home $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>
#include <memory_management.h>
#include <syscall.h>



// ===========================================================================
// _sys_warning()               Issue a warning message
// ===========================================================================
// * INPUTS
//   char *filename             Filename that caused the warning
//   int line_nr                Line number in filename that caused the 
//                              warning
//   char *format               printf() format string
//   ...                        printf() vararg arguments
//
// * RETURN VALUE
//   int                        printf() return status
// ===========================================================================
int _sys_warning(char *filename, int line_nr, char *format, ...) {
  
  Context *context;
  va_list args;
  char    buf[512];    // FIXME: implement kt_snprintf() and use it here
  int     ret;

  context = mm_get_context(ar_get_core_id());

  //kt_sprintf(buf, "Core %*d: Warning: %s [%s:%d]\n",
  //    kt_int_log10(context->pr_num_cores) + 1, context->pr_core_id,
  //    format, filename, line_nr);
  kt_sprintf(buf, "Core %d: Warning: %s [%s:%d]\n",
      context->pr_core_id, format, filename, line_nr);

  va_start(args, format);
  ret = kt_vprintf(buf, args);
  va_end(args);

  return ret;
}


// ===========================================================================
// _sys_error()                 Issue an error message and panic
// ===========================================================================
// * INPUTS
//   char *filename             Filename that caused the warning
//   int line_nr                Line number in filename that caused the 
//                              warning
//   char *format               printf() format string
//   ...                        printf() vararg arguments
// ===========================================================================
void _sys_error(char *filename, int line_nr, char *format, ...) {
  
  Context *context;
  va_list args;
  char    buf[512];    // FIXME: implement kt_snprintf() and use it here
  int     ret;

  context = mm_get_context(ar_get_core_id());

  //kt_sprintf(buf, "Core %*d: Error: %s [%s:%d]\n",
  //    kt_int_log10(context->pr_num_cores) + 1, context->pr_core_id,
  //    format, filename, line_nr);
  kt_sprintf(buf, "Core %d: Error: %s [%s:%d]\n",
      context->pr_core_id, format, filename, line_nr);

  va_start(args, format);
  ret = kt_vprintf(buf, args);
  va_end(args);

  ar_halt();
}

#if 0
// ===========================================================================
// _sys_track_sched_time()      Get system time, compute difference with
//                              given old time and credit it as scheduler
//                              communication time.
// ===========================================================================
// * INPUTS
//   unsigned int old_sec       Old time in seconds
//   unsigned int old_usec      Old time fractional part in microseconds
// ===========================================================================
static inline void _sys_track_sched_time(unsigned int old_sec, 
                                        unsigned int old_usec) {

  Context      *context;
  unsigned int now_sec;
  unsigned int now_usec;
  int          tmp;


  // Get current time
  ar_get_time(&now_sec, &now_usec);

  // Get context
  context = mm_get_context();

  // Microsec difference
  tmp = now_usec - old_usec;
  if (tmp >= 0) {
    // Straight case
    context->sys_sched_sec  += now_sec - old_sec;
    context->sys_sched_usec += tmp;
  }
  else {
    // Microsec rollover case
    context->sys_sched_sec  += now_sec - old_sec - 1;
    context->sys_sched_usec += 1000000 + tmp;
  }

  // Fix microsec overflow in context timer
  if (context->sys_sched_usec > 1000000) {
    context->sys_sched_sec++;
    context->sys_sched_usec -= 1000000;
  }
}


// ===========================================================================
// _sys_track_worker_time()     Get system time, compute difference with
//                              given old time and credit it as worker
//                              communication time.
// ===========================================================================
// * INPUTS
//   unsigned int old_sec       Old time in seconds
//   unsigned int old_usec      Old time fractional part in microseconds
// ===========================================================================
static inline void _sys_track_worker_time(unsigned int old_sec, 
                                         unsigned int old_usec) {

  Context      *context;
  unsigned int now_sec;
  unsigned int now_usec;
  int          tmp;


  // Get current time
  ar_get_time(&now_sec, &now_usec);

  // Get context
  context = mm_get_context();

  // Microsec difference
  tmp = now_usec - old_usec;
  if (tmp >= 0) {
    // Straight case
    context->sys_worker_sec  += now_sec - old_sec;
    context->sys_worker_usec += tmp;
  }
  else {
    // Microsec rollover case
    context->sys_worker_sec  += now_sec - old_sec - 1;
    context->sys_worker_usec += 1000000 + tmp;
  }

  // Fix microsec overflow in context timer
  if (context->sys_worker_usec > 1000000) {
    context->sys_worker_sec++;
    context->sys_worker_usec -= 1000000;
  }
}


// ===========================================================================
// _sys_zerostats()             Initializes all profiling stats to 0
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
// ===========================================================================
void _sys_zerostats(char *filename, int line_nr) {

  Context *context;

  context = mm_get_context();

  context->sys_sched_sec = 0;
  context->sys_sched_usec = 0;
  context->sys_sched_calls = 0;
  context->sys_worker_sec = 0;
  context->sys_worker_usec = 0;
  context->sys_worker_sent = 0;
  context->sys_worker_recv = 0;
  context->sys_barrier_sec = 0;
  context->sys_barrier_usec = 0;
}


// ===========================================================================
// _sys_getstats()              Returns all available statistics
// ===========================================================================
// * INPUTS
//   char *filename             Source code filename where this call is done
//   int line_nr                Line number in filename this call is done
//
// * OUTPUTS
//   unsigned int *ret_sched_sec     If not NULL, time spent in scheduler
//                                   communication, in seconds
//   unsigned int *ret_sched_usec    If not NULL, fractional part of the 
//                                   above, in microseconds
//   unsigned int *ret_sched_calls   If not NULL, how many calls to the
//                                   scheduler were performed
//   unsigned int *ret_worker_sec    If not NULL, time spent in worker
//                                   communication, in seconds
//   unsigned int *ret_worker_usec   If not NULL, fractional part of the
//                                   above, in microseconds
//   unsigned long *ret_worker_sent  If not NULL, total data sent to all other
//                                   workers, in bytes
//   unsigned long *ret_worker_recv  If not NULL, total data received from
//                                   all other workers, in bytes
//   unsigned int *ret_barrier_sec   If not NULL, time spent waiting on
//                                   barriers, in seconds
//   unsigned int *ret_barrier_usec  If not NULL, fractional part of the
//                                   above, in microseconds
// ===========================================================================
void _sys_getstats(char *filename, int line_nr, unsigned int *ret_sched_sec, 
                  unsigned int *ret_sched_usec, unsigned int *ret_sched_calls, 
                  unsigned int *ret_worker_sec, unsigned int *ret_worker_usec, 
                  unsigned long *ret_worker_sent, 
                  unsigned long *ret_worker_recv,
                  unsigned int *ret_barrier_sec, 
                  unsigned int *ret_barrier_usec) {
  
  Context               *context;

  context = mm_get_context();

  if (ret_sched_sec) {
    *ret_sched_sec = context->sys_sched_sec;
  }
  if (ret_sched_usec) {
    *ret_sched_usec = context->sys_sched_usec;
  }
  if (ret_sched_calls) {
    *ret_sched_calls = context->sys_sched_calls;
  }

  if (ret_worker_sec) {
    *ret_worker_sec = context->sys_worker_sec;
  }
  if (ret_worker_usec) {
    *ret_worker_usec = context->sys_worker_usec;
  }
  if (ret_worker_sent) {
    *ret_worker_sent = context->sys_worker_sent;
  }
  if (ret_worker_recv) {
    *ret_worker_recv = context->sys_worker_recv;
  }

  if (ret_barrier_sec) {
    *ret_barrier_sec = context->sys_barrier_sec;
  }
  if (ret_barrier_usec) {
    *ret_barrier_usec = context->sys_barrier_usec;
  }
}


#endif
