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
// ===========================================================================
// The code in this file is based on the pyascii85 implementation by 
// Alexandre Fiori. The original code is licensed under the GPL v2.
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : ASCII85 encoding
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: ascii85.c,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2012/10/25 12:07:50 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>

#define KT_ASCII85_LINE_WIDTH   78      // Max chars per line


// ===========================================================================
// ===========================================================================
void kt_encode85_tuple(unsigned int tuple, char *line_buf, int *pos) {

  int i;
  char buf[5];
  char *s = buf;

  i = 5;
  do {
    *s++ = tuple % 85;
    tuple /= 85;
  } while (--i > 0);

  i = 5;
  do {
    --s;
    line_buf[(*pos)++] = *s + '!';
    if ((*pos) >= KT_ASCII85_LINE_WIDTH) {
      line_buf[(*pos)++] = '\r';
      line_buf[(*pos)++] = '\n';
      line_buf[(*pos)++] = 0;
      ar_uart_send_str(line_buf);
      *pos = 0;
    }
  } while (--i > 0);
}


// ===========================================================================
// kt_encode85()                Prints a binary array of 32-bit words
//                              using the ASCII85 encoding (Adobe variant). All
//                              output characters are in the range of ASCII
//                              33-126, i.e. the 94 printable characters.
// ===========================================================================
// * INPUTS
//   unsigned int *array        The binary array to encode
//   int length                 Number of 32-bit words in array
// ===========================================================================
void kt_encode85(unsigned int *array, int length) {

  int           pos;
  char          line_buf[KT_ASCII85_LINE_WIDTH + 3];
  int           i;

  pos = 0;

  // Init tag
  line_buf[pos++] = '<';
  line_buf[pos++] = '~';

  // Encode
  for (i = 0; i < length; i++) {
    if (array[i] == 0) {
      line_buf[pos++] = 'z';
      if (pos >= KT_ASCII85_LINE_WIDTH) {
        line_buf[pos++] = '\r';
        line_buf[pos++] = '\n';
        line_buf[pos++] = 0;
        ar_uart_send_str(line_buf);
        pos = 0;
      }
    } 
    else {
      kt_encode85_tuple(array[i], line_buf, &pos);
    }
  }

  // End tag
  if (pos + 2 > KT_ASCII85_LINE_WIDTH) {
    line_buf[pos++] = '\r';
    line_buf[pos++] = '\n';
    line_buf[pos++] = 0;
    ar_uart_send_str(line_buf);
    pos = 0;
  }

  line_buf[pos++] = '~';
  line_buf[pos++] = '>';
  line_buf[pos++] = '\r';
  line_buf[pos++] = '\n';
  line_buf[pos++] = 0;
  ar_uart_send_str(line_buf);
}

