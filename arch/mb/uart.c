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
// Abstract      : Formic UART-accessing functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: uart.c,v $
// CVS revision  : $Revision: 1.11 $
// Last modified : $Date: 2013/03/16 12:47:03 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <mbs_regs.h>
#include <formic_regs.h>
#include <kernel_toolset.h>
#include <noc.h>
#include <memory_management.h>
#include <arch.h>


// ============================================================================
// Waits until the UART FIFO is completely emptied
// ============================================================================
void ar_uart_flush() {

  unsigned int core_id;
  unsigned int board_id;
  unsigned int space;

  core_id = *MBS_CPU_STATUS & 0x7;
  board_id = *MBS_CPU_STATUS >> 24;

  do {
    *MBS_MNI_SRC_ADR  = (int) FORMIC_UART_STATUS;
    *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
    *MBS_MNI_DMA_SIZE = 4;
    *MBS_MNI_BRD_NODE = (AR_FORMIC_UART_BID << 20) | (0xF << 16);
    *MBS_MNI_OPCODE   = (board_id << 20) | (core_id << 16) | 0x2;

    space = *MBS_MSL_ACCESS >> 16;
  } while (space < 1024);
}


// ============================================================================
// Sends a string to the UART
// ============================================================================
void ar_uart_send_str(char *s) {
  unsigned int core_id;
  unsigned int board_id;
  volatile char *buf;
  volatile char *c;
  int len;
  int total_len;
  unsigned int space;

  // Get us a properly aligned per-core char buffer
  core_id = *MBS_CPU_STATUS & 0x7;
  board_id = *MBS_CPU_STATUS >> 24;
  buf = (char *) (MM_MB_VA_PRINT_BUF + MM_MB_PRINT_BUF_SIZE * core_id);
  c = buf + 1;

#define HIGH_THRESHOLD 512
#define LOW_THRESHOLD  128

  // Wait until there's enough space in the UART TX FIFO. We leave some
  // additional space (HIGH_THRESHOLD), so that a few processors can print
  // intermixed messages by using this function, without overflowing the FIFO.
  // UART space will be rechecked at every LOW_THRESHOLD bytes sent.
  //
  // Note that in general ar_uart_send_str() does not guarantee overflow
  // protection: if many processors read the FIFO level from FORMIC_UART_STATUS
  // and it is near the end, they may all decide that there's enough space
  // for a single message; however, many simultaneous messages will result
  // in an overflow.
  len = kt_strlen(s) + 64;
  if (len > MM_MB_PRINT_BUF_SIZE - 2) {
    ar_uart_send_str("ERROR: Message too large\r\n");
    return;
  }

  do {
    *MBS_MNI_SRC_ADR  = (int) FORMIC_UART_STATUS;
    *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
    *MBS_MNI_DMA_SIZE = 4;
    *MBS_MNI_BRD_NODE = (AR_FORMIC_UART_BID << 20) | (0xF << 16);
    *MBS_MNI_OPCODE   = (board_id << 20) | (core_id << 16) | 0x2;

    space = *MBS_MSL_ACCESS >> 16;
  } while (space <= HIGH_THRESHOLD);

  len = 0;
  total_len = 0;
  while (*s) {
    *c++ = *s++;
    len++;

    // Send each 63-byte chunk
    if (len == 63) {
      buf[0] = 63;

      ar_cnt_set(core_id, NOC_COUNTER_UART, -64);
      ar_dma_with_ack(core_id,
                      board_id,           core_id, (int) buf,
                      AR_FORMIC_UART_BID, 0xF,     (int) FORMIC_UART_BULK_WRITE,
                      board_id,           core_id, NOC_COUNTER_UART, 
                      64, 0, 0, 0);
      while (ar_cnt_get(core_id, NOC_COUNTER_UART)) {
        ;
      }

      len = 0;
      total_len += 64;
      c = buf + 1;

      if (total_len >= LOW_THRESHOLD) {
        do {
          *MBS_MNI_SRC_ADR  = (int) FORMIC_UART_STATUS;
          *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
          *MBS_MNI_DMA_SIZE = 4;
          *MBS_MNI_BRD_NODE = (AR_FORMIC_UART_BID << 20) | (0xF << 16);
          *MBS_MNI_OPCODE   = (board_id << 20) | (core_id << 16) | 0x2;

          space = *MBS_MSL_ACCESS >> 16;
        } while (space <= HIGH_THRESHOLD);
        total_len = 0;
      }
    }
  }

  // Send the remainder
  if (len) {
    buf[0] = len;

    ar_cnt_set(core_id, NOC_COUNTER_UART, -64);
    ar_dma_with_ack(core_id,
                    board_id,           core_id, (int) buf,
                    AR_FORMIC_UART_BID, 0xF,     (int) FORMIC_UART_BULK_WRITE,
                    board_id,           core_id, NOC_COUNTER_UART, 
                    64, 0, 0, 0);
    while (ar_cnt_get(core_id, NOC_COUNTER_UART)) {
      ;
    }
  }
}


// ============================================================================
// Sends a character to the UART
// ============================================================================
void ar_uart_send_char(int c) {

  unsigned int core_id;
  unsigned int board_id;
  unsigned int space;


  core_id = *MBS_CPU_STATUS & 0x7;
  board_id = *MBS_CPU_STATUS >> 24;

  // Wait until there's enough space in the UART TX FIFO
  do {
    *MBS_MNI_SRC_ADR  = (int) FORMIC_UART_STATUS;
    *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
    *MBS_MNI_DMA_SIZE = 4;
    *MBS_MNI_BRD_NODE = (AR_FORMIC_UART_BID << 20) | (0xF << 16);
    *MBS_MNI_OPCODE   = (board_id << 20) | (core_id << 16) | 0x2;

    space = *MBS_MSL_ACCESS >> 16;
  } while (!space);

  // Send character
  *MBS_MNI_DST_ADR  = (int) FORMIC_UART_BYTE;
  *MBS_MNI_MSG0     = c;
  *MBS_MNI_OPCODE   = (board_id << 20) | (0xF << 16) | 0;

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
// Not implemented for MicroBlaze. If needed, imlpement this by reading the
// UART board's RX status and fetching a byte from there if it exists.
// ============================================================================
int ar_uart_get_char() {
  ar_abort();
}
