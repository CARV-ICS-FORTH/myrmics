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
// Abstract      : Counters, timers, watchdogs and anything else that counts
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: timer.c,v $
// CVS revision  : $Revision: 1.10 $
// Last modified : $Date: 2013/06/05 03:15:50 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <arch.h>
#include <arm_regs.h>
#include <types.h>


// ###########################################################################
// ###                                                                     ###
// ###                    Separated low-level interface                    ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ar_timer0_set()              Sets the Private timer to a value and
//                              enables countdown. Each tick of this counter
//                              is 2 processor clock cycles (i.e. 200 MHz).
// ===========================================================================
// * INPUTS
//   unsigned int value         Value to load
// ===========================================================================
void ar_timer0_set(unsigned int value) {
  // Disable counter
  *ARM_PRIV_TIMER_CONTROL = 0;

  // Load it with value
  *ARM_PRIV_TIMER_LOAD = value;

  // Enable countdown
  *ARM_PRIV_TIMER_CONTROL = 1;
}

// ===========================================================================
// ar_timer0_read()             Reads the Private timer value
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value
// ===========================================================================
unsigned int ar_timer0_read() {
  return *ARM_PRIV_TIMER_COUNTER;
}

// ===========================================================================
// ar_timer0_busy_wait()        Uses the Private timer to busy wait for the
//                              given value of milliseconds
// ===========================================================================
// * INPUTS
//   unsigned int value         Milliseconds to wait
// ===========================================================================
void ar_timer0_busy_wait(unsigned int msec) {

  // Start counter with timeout value
  *ARM_PRIV_TIMER_CONTROL = 0;
  *ARM_PRIV_TIMER_LOAD = msec * 200000;
  *ARM_PRIV_TIMER_CONTROL = 1;

  // Wait
  while (*ARM_PRIV_TIMER_COUNTER) {
    ;
  }
}


// ===========================================================================
// ar_timer1_set()              Sets the Watchdog timer to a value and
//                              enables countdown. Each tick of this counter
//                              is 2 processor clock cycles (i.e. 200 MHz).
// ===========================================================================
// * INPUTS
//   unsigned int value         Value to load
// ===========================================================================
void ar_timer1_set(unsigned int value) {

  // Disable counter
  *ARM_WATCHDOG_TIMER_CONTROL = 0;

  // Ensure watchdog is in timer mode
  *ARM_WATCHDOG_TIMER_DISABLE = 0x12345678;
  *ARM_WATCHDOG_TIMER_DISABLE = 0x87654321;

  // Load it with value
  *ARM_WATCHDOG_TIMER_LOAD = value;

  // Enable countdown
  *ARM_WATCHDOG_TIMER_CONTROL = 1;
}

// ===========================================================================
// ar_timer1_read()             Reads the Watchdog timer value
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value
// ===========================================================================
unsigned int ar_timer1_read() {
  return *ARM_WATCHDOG_TIMER_COUNTER;
}

// ===========================================================================
// ar_timer1_busy_wait()        Uses the Watchdog timer to busy wait for the
//                              given value of milliseconds.
// ===========================================================================
// * INPUTS
//   unsigned int value         Milliseconds to wait
// ===========================================================================
void ar_timer1_busy_wait(unsigned int msec) {

  // Start counter with timeout value
  *ARM_WATCHDOG_TIMER_CONTROL = 0;
  *ARM_WATCHDOG_TIMER_DISABLE = 0x12345678;
  *ARM_WATCHDOG_TIMER_DISABLE = 0x87654321;
  *ARM_WATCHDOG_TIMER_LOAD = msec * 200000;
  *ARM_WATCHDOG_TIMER_CONTROL = 1;

  // Wait
  while (*ARM_WATCHDOG_TIMER_COUNTER) {
    ;
  }
}


// ===========================================================================
// ar_timer2_read()             Reads the global timer value (each timer tick
//                              is 40 clock cycles, i.e. 100 ns or 10 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value
// ===========================================================================
unsigned int ar_timer2_read() {
  return *ARS_TMR_GLOBAL(ar_get_core_id());
}


// ###########################################################################
// ###                                                                     ###
// ###                    Unified high-level interface                     ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ar_timer_reset()             Sets the Private timer to the max value and
//                              enables countdown. Each tick of this counter
//                              is 2 processor clock cycles (i.e. 200 MHz).
// ===========================================================================
void ar_timer_reset() {
  // Disable counter
  *ARM_PRIV_TIMER_CONTROL = 0;

  // Load it with max value
  *ARM_PRIV_TIMER_LOAD = 0xFFFFFFFF;

  // Enable countdown
  *ARM_PRIV_TIMER_CONTROL = 1;
}


// ===========================================================================
// ar_timer_get_ticks()         Reads the Private timer value from the last
//                              reset time and returns it in counter ticks
//                              (each timer tick is 2 clock cycles, i.e. 5 ns
//                              or 200 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in counter ticks
// ===========================================================================
unsigned int ar_timer_get_ticks() {
  return (0xFFFFFFFF - *ARM_PRIV_TIMER_COUNTER);
}


// ===========================================================================
// ar_timer_get_cycles()        Reads the Private timer value from the last
//                              reset time and returns it in clock cycles
//                              (each timer tick is 2 clock cycles, i.e. 5 ns
//                              or 200 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in clock cycles
// ===========================================================================
unsigned int ar_timer_get_cycles() {
  unsigned int ticks;

  ticks = 0xFFFFFFFF - *ARM_PRIV_TIMER_COUNTER;
  if (ticks > 0x7FFFFFFF) {
    return 0xFFFFFFFF;
  }
  else {
    return 2 * ticks;
  }
}


// ===========================================================================
// ar_timer_get_msec()          Reads the Private timer value from the last
//                              reset time and returns it in milliseconds
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Elapsed time in milliseconds
// ===========================================================================
unsigned int ar_timer_get_msec() {
  unsigned int ticks;
  unsigned int msec;

  ticks = 0xFFFFFFFF - *ARM_PRIV_TIMER_COUNTER;
  ar_uint_divide(ticks, 200000, &msec, NULL);
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
  unsigned int ticks;
  unsigned int usec;

  ticks = 0xFFFFFFFF - *ARM_PRIV_TIMER_COUNTER;
  ar_uint_divide(ticks, 200, &usec, NULL);
  return (usec);
}


// ===========================================================================
// ar_timer_busy_wait_ticks()   Uses the Private timer to busy-wait for
//                              a number of counter ticks (each 400-MHz clock
//                              cycle is 2 counter ticks)
// ===========================================================================
// * INPUTS
//   unsigned int ticks        Ticks to wait
// ===========================================================================
void ar_timer_busy_wait_ticks(unsigned int ticks) {

  // Disable counter
  *ARM_PRIV_TIMER_CONTROL = 0;

  // Load it
  *ARM_PRIV_TIMER_LOAD = ticks;

  // Enable countdown
  *ARM_PRIV_TIMER_CONTROL = 1;

  // Wait
  while (*ARM_PRIV_TIMER_COUNTER) {
    ;
  }
}


// ===========================================================================
// ar_timer_busy_wait_cycles()  Uses the Private timer to busy-wait for
//                              a number of clock cycles (each 400-MHz clock
//                              cycle is 2 counter ticks)
// ===========================================================================
// * INPUTS
//   unsigned int cycles        Cycles to wait
// ===========================================================================
void ar_timer_busy_wait_cycles(unsigned int cycles) {

  // Disable counter
  *ARM_PRIV_TIMER_CONTROL = 0;

  // Load it
  *ARM_PRIV_TIMER_LOAD = cycles / 2;

  // Enable countdown
  *ARM_PRIV_TIMER_CONTROL = 1;

  // Wait
  while (*ARM_PRIV_TIMER_COUNTER) {
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

  // Disable counter
  *ARM_PRIV_TIMER_CONTROL = 0;

  // Load it
  *ARM_PRIV_TIMER_LOAD = msec * 200000;

  // Enable countdown
  *ARM_PRIV_TIMER_CONTROL = 1;

  // Wait
  while (*ARM_PRIV_TIMER_COUNTER) {
    ;
  }
}


// ===========================================================================
// ar_free_timer_get_ticks()    Reads the free-running timer value
//                              (each timer tick is 40 clock cycles, i.e. 100 ns
//                              or 10 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value in counter ticks
// ===========================================================================
unsigned int ar_free_timer_get_ticks() {
  return *ARS_TMR_GLOBAL(ar_get_core_id());
}


// ===========================================================================
// ar_free_timer_get_cycles()   Reads the free-running timer value
//                              (each timer tick is 40 clock cycles, i.e. 100 ns
//                              or 10 MHz).
// ===========================================================================
// * RETURN VALUE
//   unsigned int               Current value in CPU clock cycles
// ===========================================================================
unsigned int ar_free_timer_get_cycles() {
  return *ARS_TMR_GLOBAL(ar_get_core_id()) * 40;
}
