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
// Abstract      : Exceptional cases where the kernel must crash
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: panic.c,v $
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2013/05/21 06:14:50 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <mbs_regs.h>
#include <formic_regs.h>
#include <arch.h>


#define AR_MAX_BACKTRACE 10000  // Max number of instructions to search
                                // for frame size during backtrace

// ===========================================================================
// ar_backtrace()               Dumps all the return addresses that lead
//                              to the current stack frame
// ===========================================================================
// * INPUTS
//   int my_bid                 Board ID
//   int my_cid                 Core ID
//   unsigned int pc            Program counter to start searching
//   unsigned int sp            Stack pointer of initial frame
// ===========================================================================
void ar_backtrace(int my_bid, int my_cid, unsigned int pc, unsigned int sp) {

  unsigned int  *cur;
  unsigned int  instr;
  int           found;
  int           frame_size;
  int           depth;
  int           i;

  cur = (unsigned int *) pc;
  depth = 0;
  while (1) {
    //kt_printf("0x%02X/%d: Backtrace depth %d, PC = 0x%08X\r\n",
    //          my_bid, my_cid, depth, cur);
    kt_printf("%02X/%d [%d] %8X\r\n", my_bid, my_cid, depth, cur);

    // Find out the frame size by inspecting previous instructions
    found = 0;
    for (i = 0; i < (AR_MAX_BACKTRACE) && cur; i++, cur--) {
      instr = *cur;
      // Search "addik r1, r1, *" instruction that creates the stack frame.
      // The immediate should be negative (bit 15 == 1).
      if ((instr & 0xFFFF8000) == 0x30218000) {
        found = 1;
        break;
      }
    }

    if (found) {
      // Frame size is the immediate part of the instruction, inversed
      frame_size = instr & 0xFFFF;
      frame_size |= 0xFFFF0000; // sign extend to 32 bits
      frame_size = -frame_size; // inverse
    }
    else {
      // Something went wrong, we're stopping...
      break;
    }

    // Proceed with link register of the current frame
    cur = (unsigned int *) *((int *) sp);

    // Pop stack frame
    sp += frame_size;

    // Did we reach *sp == 0?
    if (!cur) {
      break;
    }

    depth++;
  }
}


// ===========================================================================
// ar_halt()                    Halts the system
// ===========================================================================
// * RETURN VALUE
//   -                          Does not return
// ===========================================================================
void ar_halt() {

  kt_printf("0x%02X/%d: System halted.\r\n",
      ar_get_board_id(), ar_get_core_id());

  // Stall
  while (1) {
    ;
  }
}


// ===========================================================================
// ar_panic()                   Prints a panic message and halts the system
// ===========================================================================
// * INPUTS
//   char *msg                  Message to print
//
// * RETURN VALUE
//   -                          Does not return
// ===========================================================================
void ar_panic(char *msg) {

  int my_bid;
  int my_cid;
  unsigned int pc;
  unsigned int sp;


  // Get board/core ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();

  // Blink board lights
  ar_write_brd_control(my_cid, my_bid, 0xFFFFFF00);

  /*
  // Print dead ant and message
  ar_uart_send_str("\r\n\r\n"
    "    :::::::\r\n"
    "    :::::'\r\n"
    "    :::'                     46    4a    ?? jf\r\n"
    "    :'                    )Wa_P )4aj(    jf ]f\r\n"
    "  _TT_                    ]f  _[Q   j'   ]f <f    ?,\r\n"
    " /____\\                 )??y'  __P   J   ]f )[   ]/     '\r\n"
    " |    |               )WQQQQ4 )6  J )\\  ?P4a   ?QQQ\?\?/\r\n"
    " |    |               mQQQQQQD6   4]a  4QQ/J'4?QQQQQ\r\n"
    " |    |               ]QQQQQQQQQPQ/Q)4jQ/  f4QQga jL'\r\n"
    " |    |               _QQQQQQQQQQQQQQ  aj? QQQQW  Q /\r\n"
    " |    |                jQQQQQQQQQQQQQ'?QQQQQQ5Qa?QQ)\r\n"
    " |    |                 aQQQQQQQQQQQQQQQQQQQQQaayQ//\r\n"
    " |    |                    _aaaa sQQQ6QQQQ6/]     a\"'\r\n"
    " |    |                            aa _aa/  ]'       aa)\?\?)\r\n"
    " |    |                                      J' ?'\r\n"
    " |    |                                       aaa\r\n"
    " |    |  Kernel panic:\r\n"
    " |____|  "
  );
  */
  kt_printf("0x%02X/%d: Kernel panic: %s\r\n", my_bid, my_cid, msg);

  // Attempt a backtrace of calls that led us here
  asm("mfs %0, rpc" : "=r"(pc) );
  asm("add %0, r1, r0" : "=r"(sp) );
  ar_backtrace(my_bid, my_cid, pc, sp);

  // Halt
  ar_halt();
}
