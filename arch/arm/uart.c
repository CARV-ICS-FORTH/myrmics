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
// Abstract      : ARM UART device driver. Note that UART initialization has
//                 already happened during boot.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: uart.c,v $
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2013/03/16 12:47:03 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <arm_regs.h>
#include <memory_management.h>


// ============================================================================
// Waits until the UART FIFO is completely emptied
// ============================================================================
void ar_uart_flush() {

  volatile unsigned int flag;

  flag = *ARM_UART0_FR;
  while (!((flag >> 7) & 0x1)) {
    flag = *ARM_UART0_FR;
  }
}


// ============================================================================
// Sends a character to the UART and waits if the FIFO gets full
// ============================================================================
void ar_uart_send_char(int c) {

  volatile unsigned int flag;

  flag = *ARM_UART0_FR;
  while ((flag >> 5) & 0x1) {
    flag = *ARM_UART0_FR;
  }

  *ARM_UART0_DR = c;
}


// ============================================================================
// Sends a string to the UART
// ============================================================================
void ar_uart_send_str(char *s) {

  int core_id = ar_get_core_id();


  // Get UART lock
  while (ar_lock_get((void *) MM_ARM_VA_PRINT_LOCK, core_id)) {
    ;
  }

  // Send string, char-by-char
  while (*s) {
    ar_uart_send_char(*s++);
  }


  // Release UART lock
  while (ar_lock_release((void *) MM_ARM_VA_PRINT_LOCK)) {
    ;
  }
}


// ============================================================================
// Prints out a number in hexadecimal
// ============================================================================
void ar_uart_send_hex(unsigned int hex) {
  int i, j;

  for (i = 28; i >= 0; i -= 4) {
    j = (hex >> i) & 0xf;
    if (j < 10) {
      j += '0';
    }
    else {
      j += 'A' - 10;
    }
    ar_uart_send_char(j);
  }
}


// ============================================================================
// Prints out a character in hexadecimal
// ============================================================================
void ar_uart_send_char_hex(unsigned char hex) {
  int i, j;

  for (i = 4; i >= 0; i -= 4) {
    j = (hex >> i) & 0xf;
    if (j < 10) {
      j += '0';
    }
    else {
      j += 'A' - 10;
    }
    ar_uart_send_char(j);
  }
}


// ============================================================================
// Prints out a number in decimal (FIXME: do it with ar_uint_divide())
// ============================================================================
//void ar_uart_send_dec(unsigned int dec) {
//  int i, j;
//  unsigned char val[10];
//
//  j = 0;
//  for (i = 0; i <= 10; i++) {
//    val[j++] = dec % 10;
//    dec = dec / 10;
//    if (dec == 0) {
//      break;
//    }
//  }
//
//  for (j-- ; j >= 0; j--) {
//    ar_uart_send_char(val[j] + '0');
//  }
//}


// ============================================================================
// Reads a character from the UART, if it has any. Returns the character (if
// available) or -1 (if no character was found)
// ============================================================================
int ar_uart_get_char() {

  if ((*ARM_UART0_FR >> 4) & 0x1) {
    return -1;
  }

  return *ARM_UART0_DR & 0xFF;
}

