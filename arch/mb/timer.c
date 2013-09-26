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
// Abstract      : Timer-related functions
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: timer.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2013/07/05 15:41:14 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <mbs_regs.h>
#include <arch.h>
#include <types.h>


// ###########################################################################
// ###                                                                     ###
// ###                    Separated low-level interface                    ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ar_priv_timer_set()          Sets the Private timer to a value and
//                              enables countdown. Each tick of this counter
//                              is 1 processor clock cycle (i.e. 10 MHz).
// ===========================================================================
// * INPUTS
//   unsigned int value         Value to load
// ===========================================================================
void ar_priv_timer_set(unsigned int value) {

  // Load timer
  *MBS_TMR_PRIVATE = (1 << 31) | (0 << 30) | value;
}


// ===========================================================================
// ar_priv_timer_read()         Reads the Private timer value
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value
// ===========================================================================
unsigned int ar_priv_timer_read() {
  return (*MBS_TMR_PRIVATE & 0x3FFFFFFF);
}


// ===========================================================================
// ar_priv_timer_busy_wait()    Uses the Private timer to busy wait for the
//                              given value of milliseconds
// ===========================================================================
// * INPUTS
//   unsigned int value         Milliseconds to wait
// ===========================================================================
void ar_priv_timer_busy_wait(unsigned int msec) {

  // Start counter with timeout value
  *MBS_TMR_PRIVATE = (1 << 31) | (0 << 30) | (msec * 10000);

  // Wait
  while (*MBS_TMR_PRIVATE & 0x3FFFFFFF) {
    ;
  }
}


// ===========================================================================
// ar_glob_timer_read()         Reads the global timer value. Each tick of this
//                              counter is 1 processor clock cycle (100 ns or
//                              10 MHz)
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value
// ===========================================================================
unsigned int ar_glob_timer_read() {
  return *MBS_TMR_GLOBAL;
}


// ###########################################################################
// ###                                                                     ###
// ###                    Unified high-level interface                     ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ar_timer_reset()             Sets the Private timer to the max value and
//                              enables countdown. Each tick of this counter
//                              is 1 processor clock cycle (i.e. 10 MHz).
// ===========================================================================
void ar_timer_reset() {
  *MBS_TMR_PRIVATE = (1 << 31) | (0 << 30) | 0x3FFFFFFF;
}


// ===========================================================================
// ar_timer_get_ticks()         Reads the Private timer value from the last
//                              reset time and returns it in ticks
//                              (each tick = 1 CPU cycle = 100 ns = 10 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in counter ticks
// ===========================================================================
unsigned int ar_timer_get_ticks() {
  return (0x3FFFFFFF - (*MBS_TMR_PRIVATE & 0x3FFFFFFF));
}


// ===========================================================================
// ar_timer_get_cycles()        Reads the Private timer value from the last
//                              reset time and returns it in clock cycles
//                              (each tick = 1 CPU cycle = 100 ns = 10 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in clock cycles
// ===========================================================================
unsigned int ar_timer_get_cycles() {
  return (0x3FFFFFFF - (*MBS_TMR_PRIVATE & 0x3FFFFFFF));
}


// ===========================================================================
// ar_timer_get_msec()          Reads the Private timer value from the last
//                              reset time and returns it in milliseconds
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in milliseconds
// ===========================================================================
unsigned int ar_timer_get_msec() {
  unsigned int cycles;
  unsigned int msec;

  cycles = 0x3FFFFFFF - (*MBS_TMR_PRIVATE & 0x3FFFFFFF);
  ar_uint_divide(cycles, 10000, &msec, NULL);
  return (msec);
}


// ===========================================================================
// ar_timer_get_usec()          Reads the Private timer value from the last
//                              reset time and returns it in microseconds
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in microseconds
// ===========================================================================
unsigned int ar_timer_get_usec() {
  unsigned int cycles;
  unsigned int usec;

  cycles = 0x3FFFFFFF - (*MBS_TMR_PRIVATE & 0x3FFFFFFF);
  ar_uint_divide(cycles, 10, &usec, NULL);
  return (usec);
}


// ===========================================================================
// ar_timer_busy_wait_ticks()   Uses the Private timer to busy-wait for
//                              a number of counter ticks (1 tick = 1 CPU
//                              cycle = 100 ns = 10 MHz).
// ===========================================================================
// * INPUTS
//   unsigned int ticks         Counter ticks to wait
// ===========================================================================
void ar_timer_busy_wait_ticks(unsigned int ticks) {

  // Load timer
  *MBS_TMR_PRIVATE = (1 << 31) | (0 << 30) | ticks;

  // Wait
  while (*MBS_TMR_PRIVATE & 0x3FFFFFFF) {
    ;
  }
}


// ===========================================================================
// ar_timer_busy_wait_cycles()  Uses the Private timer to busy-wait for
//                              a number of CPU clock cycles (1 tick = 1 CPU
//                              cycle = 100 ns = 10 MHz).
// ===========================================================================
// * INPUTS
//   unsigned int cycles        Cycles to wait
// ===========================================================================
void ar_timer_busy_wait_cycles(unsigned int cycles) {

  // Load timer
  *MBS_TMR_PRIVATE = (1 << 31) | (0 << 30) | cycles;

  // Wait
  while (*MBS_TMR_PRIVATE & 0x3FFFFFFF) {
    ;
  }
}


// ===========================================================================
// ar_timer_busy_wait_msec()    Uses the Private timer to busy-wait for
//                              a number of milliseconds
// ===========================================================================
// * INPUTS
//   unsigned int msec          Milliseconds to wait
// ===========================================================================
void ar_timer_busy_wait_msec(unsigned int msec) {

  // Load timer
  *MBS_TMR_PRIVATE = (1 << 31) | (0 << 30) | (msec * 10000);

  // Wait
  while (*MBS_TMR_PRIVATE & 0x3FFFFFFF) {
    ;
  }
}


// ===========================================================================
// ar_free_timer_get_ticks()    Reads the global timer value. Each tick of this
//                              counter is 1 processor clock cycle (10 ns or
//                              10 MHz)
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value in ticks
// ===========================================================================
unsigned int ar_free_timer_get_ticks() {
  return *MBS_TMR_GLOBAL;
}


// ===========================================================================
// ar_free_timer_get_cycles()   Reads the global timer value. Each tick of this
//                              counter is 1 processor clock cycle (10 ns or
//                              10 MHz)
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value in CPU clock cycles
// ===========================================================================
unsigned int ar_free_timer_get_cycles() {
  return *MBS_TMR_GLOBAL;
}
