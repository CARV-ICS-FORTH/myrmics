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
// The code in this file is derived from the original Linux (version 0.01)
// code, which was under the following license notice:
//
// This kernel is (C) 1991 Linus Torvalds, but all or part of it may be
// redistributed provided you do the following:
// 
//         - Full source must be available (and free), if not with the
//           distribution then at least on asking for it.
// 
//         - Copyright notices must be intact. (In fact, if you distribute
//           only parts of it you may have to add copyrights, as there aren't
//           (C)'s in all files.) Small partial excerpts may be copied
//           without bothering with copyrights.
// 
//         - You may not distibute this for a fee, not even "handling"
//           costs.
// 
// Mail me at "torvalds@kruuna.helsinki.fi" if you have any questions.
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : Print-related functions
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: print.c,v $
// CVS revision  : $Revision: 1.12 $
// Last modified : $Date: 2013/04/09 16:28:12 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <stdarg.h>     // GCC built-in

#include <arch.h>
#include <kernel_toolset.h>


// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
static inline int kt_skip_atoi(const char **s) {
  int i=0;

  while (**s >= '0' && **s <= '9') {
    i = i*10 + *((*s)++) - '0';
  }
  return i;
}

//
// FIXME: redefine with KT_PRINT_ prefix

#define ZEROPAD 1               /* pad with zero */
#define SIGN    2               /* unsigned/signed long */
#define PLUS    4               /* show plus */
#define SPACE   8               /* space if plus */
#define LEFT    16              /* left justified */
#define SPECIAL 32              /* 0x */
#define SMALL   64              /* use 'abcdef' instead of 'ABCDEF' */

// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
static char *kt_parse_number(char *str, int num, int base, int size,
                             int precision, int type) {
  char c;
  char sign;
  char tmp[36];
  const char *digits="0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ";
  int i;

  if (type & SMALL) {
    digits="0123456789abcdefghijklmnopqrstuvwxyz";
  }
  if (type & LEFT) {
    type &= ~ZEROPAD;
  }
  if (base < 2 || base > 36) {
    return 0;
  }
  c = (type & ZEROPAD) ? '0' : ' ' ;
  if ((type & SIGN) && (num < 0)) {
    sign='-';
    num = -num;
  } else {
    sign = (type & PLUS) ? '+' :
                           ( (type & SPACE) ? ' ' : 0);
  }
  if (sign) {
    size--;
  }
  if (type & SPECIAL) {
    if (base==16) {
      size -= 2;
    }
    else if (base == 8) {
      size--;
    }
  }

  i = 0;
  if (num==0) {
    tmp[i++]='0';
  }
  else {
    while (num!=0) {
      unsigned int quot, rem;
      ar_uint_divide(num, base, &quot, &rem);
      tmp[i++]=digits[rem];
      num = quot;
    }
  }
  if (i > precision) {
    precision = i;
  }
  size -= precision;
  if (!(type & (ZEROPAD + LEFT))) {
    while (size-- > 0) {
      *str++ = ' ';
    }
  }
  if (sign) {
    *str++ = sign;
  }
  if (type & SPECIAL) {
    if (base==8) {
      *str++ = '0';
    }
    else if (base==16) {
      *str++ = '0';
      *str++ = digits[33];
    }
  }
  if (!(type & LEFT)) {
    while (size-- > 0) {
      *str++ = c;
    }
  }
  while (i < precision--) {
    *str++ = '0';
  }
  while(i-- > 0) {
    *str++ = tmp[i];
  }
  while(size-- > 0) {
    *str++ = ' ';
  }

  return str;
}


// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
static char *kt_parse_float(char *str, float num, int size, int precision,
                            int type) {
  int           whole;
  unsigned int  quot;
  unsigned int  rem;
  char          tmp[36];
  int           i;


  if (num < 0) {
    num = -num;
    *str++ = '-';
  }


  // Do the integer part
  whole = num;
  i = 0;
  if (whole == 0) {
    tmp[i++] = '0';
  }
  else {
    while (whole != 0) {
      ar_uint_divide(whole, 10, &quot, &rem);
      tmp[i++] = rem + '0';
      whole = quot;
    }
  }
  while (i-- > 0) {
    *str++ = tmp[i];
  }

  // Do the fractional part
  *str++ = '.';
  if (precision <= 0) {
    precision = 6;
  }
  for (i = 0; i < precision; i++) {
    num -= (int) num;
    num *= 10.0;
    whole = num;
    *str++ = whole + '0';
  }


  return str;
}


// FIXME: make it sNprintf and fix kt_vprintf

// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int kt_vsprintf(char *buf, const char *fmt, va_list args) {
  int len;
  int i;
  char * str;
  char *s;
  int *ip;

  int flags;              // flags to kt_parse_number()

  int field_width;        // width of output field
  int precision;          // min. # of digits for integers; max
                          // number of chars for from string
  int qualifier;          // 'h', 'l', or 'L' for integer fields

  for (str = buf ; *fmt ; ++fmt) {
    if (*fmt != '%') {
      *str++ = *fmt;
      continue;
    }

    // process flags
    flags = 0;
    repeat:
      ++fmt;          // this also skips first '%'
      switch (*fmt) {
        case '-': flags |= LEFT; goto repeat;
        case '+': flags |= PLUS; goto repeat;
        case ' ': flags |= SPACE; goto repeat;
        case '#': flags |= SPECIAL; goto repeat;
        case '0': flags |= ZEROPAD; goto repeat;
      }

    // get field width
    field_width = -1;
    if (*fmt >= '0' && *fmt <= '9') {
      field_width = kt_skip_atoi(&fmt);
    }
    else if (*fmt == '*') {
      // it's the next argument
      field_width = va_arg(args, int);
      if (field_width < 0) {
        field_width = -field_width;
        flags |= LEFT;
      }
    }

    // get the precision
    precision = -1;
    if (*fmt == '.') {
      ++fmt;
      if (*fmt >= '0' && *fmt <= '9') {
        precision = kt_skip_atoi(&fmt);
      }
      else if (*fmt == '*') {
        // it's the next argument
        precision = va_arg(args, int);
      }
      if (precision < 0) {
        precision = 0;
      }
    }

    // get the conversion qualifier
    qualifier = -1;
    if (*fmt == 'h' || *fmt == 'l' || *fmt == 'L') {
      qualifier = *fmt;
      ++fmt;
    }

    switch (*fmt) {
      case 'c':
        if (!(flags & LEFT)) {
          while (--field_width > 0) {
            *str++ = ' ';
          }
        }
        *str++ = (unsigned char) va_arg(args, int);
        while (--field_width > 0) {
          *str++ = ' ';
        }
        break;

      case 's':
        s = va_arg(args, char *);
        len = kt_strlen(s);
        if (precision < 0) {
          precision = len;
        }
        else if (len > precision) {
          len = precision;
        }

        if (!(flags & LEFT)) {
          while (len < field_width--) {
            *str++ = ' ';
          }
        for (i = 0; i < len; ++i) {
          *str++ = *s++;
        }
        while (len < field_width--) {
          *str++ = ' ';
        }
        break;

      case 'o':
        str = kt_parse_number(str, va_arg(args, unsigned long), 8,
                field_width, precision, flags);
        break;

      case 'p':
        if (field_width == -1) {
          field_width = 8;
          flags |= ZEROPAD;
        }
        str = kt_parse_number(str, (unsigned long) va_arg(args, void *), 16,
                field_width, precision, flags);
        break;

      case 'x':
        flags |= SMALL;
      case 'X':
        str = kt_parse_number(str, va_arg(args, unsigned long), 16,
                field_width, precision, flags);
        break;

      case 'd':
      case 'i':
        flags |= SIGN;
      case 'u':
        str = kt_parse_number(str, va_arg(args, unsigned long), 10,
                field_width, precision, flags);
        break;

      case 'n':
        ip = va_arg(args, int *);
        *ip = (str - buf);
        break;

      case 'f':
        str = kt_parse_float(str, va_arg(args, double), field_width,
                precision, flags);
        break;

      default:
        if (*fmt != '%') {
          *str++ = '%';
        }
        if (*fmt) {
          *str++ = *fmt;
        }
        else {
          --fmt;
        }
        break;
      }
    }
  }

  *str = '\0';
  return str-buf;
}





// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int kt_vprintf(const char *format, va_list ap) {
  char buf[8192]; // FIXME
  int ret;

  ret = kt_vsprintf(buf, format, ap);

  ar_uart_send_str(buf);

  return ret;
}


// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int kt_printf(const char *format, ...) {
  va_list ap;
  int     ret;

  // Initialize variable argument list
  va_start (ap, format);

  // Type the actual message
  ret = kt_vprintf (format, ap);

  // Close list
  va_end (ap);

  // Return value
  return ret;
}


// ===========================================================================
// function()                   FIXME comments
// ===========================================================================
// * INPUTS
//   unsigned char *arg1        Describe arg1
//   int arg2                   Describe arg2
//
// * OUTPUTS
//   int *arg3                  Describe arg3
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int kt_sprintf(char *buf, const char *format, ...) {
  va_list ap;
  int     ret;

  // Initialize variable argument list
  va_start (ap, format);

  // Write to buffer
  ret = kt_vsprintf(buf, format, ap);

  // Close list
  va_end (ap);

  // Return value
  return ret;
}
