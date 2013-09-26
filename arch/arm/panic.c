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
// CVS revision  : $Revision: 1.7 $
// Last modified : $Date: 2012/07/18 16:50:07 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>

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
  while (1)
    ;
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
  
  // Print dead ant and message
  //ar_uart_send_str("\r\n"
  //  "    :::::::\r\n"
  //  "    :::::'\r\n"
  //  "    :::'                     46    4a    ?? jf\r\n"
  //  "    :'                    )Wa_P )4aj(    jf ]f\r\n"
  //  "  _TT_                    ]f  _[Q   j'   ]f <f    ?,\r\n"
  //  " /____\\                 )??y'  __P   J   ]f )[   ]/     '\r\n"
  //  " |    |               )WQQQQ4 )6  J )\\  ?P4a   ?QQQ\?\?/\r\n"
  //  " |    |               mQQQQQQD6   4]a  4QQ/J'4?QQQQQ\r\n"
  //  " |    |               ]QQQQQQQQQPQ/Q)4jQ/  f4QQga jL'\r\n"
  //  " |    |               _QQQQQQQQQQQQQQ  aj? QQQQW  Q /\r\n"
  //  " |    |                jQQQQQQQQQQQQQ'?QQQQQQ5Qa?QQ)\r\n"
  //  " |    |                 aQQQQQQQQQQQQQQQQQQQQQaayQ//\r\n"
  //  " |    |                    _aaaa sQQQ6QQQQ6/]     a\"'\r\n"
  //  " |    |                            aa _aa/  ]'       aa)\?\?)\r\n"
  //  " |    |                                      J' ?'\r\n"
  //  " |    |                                       aaa\r\n"
  //  " |    |  Kernel panic:\r\n"
  //  " |____|  "
  //);
  //ar_uart_send_str(msg);
  //ar_uart_send_str("\r\n");
  kt_printf("0x%02X/%d: Kernel panic: %s\r\n", 
            ar_get_board_id(), ar_get_core_id(), msg);

  // Halt
  ar_halt();
}













         


