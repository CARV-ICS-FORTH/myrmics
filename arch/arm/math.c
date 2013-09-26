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
// The software for the integer devision routine (ar_uint_divide()) was 
// written by Ian Kaplan (http://www.bearcave.com/software/divide.htm), with
// the following copyright notice:
//
// Copyright stuff
//
// Use of this program, for any purpose, is granted the author,
// Ian Kaplan, as long as this copyright notice is included in
// the source code or any source code derived from this program.
// The user assumes all responsibility for using this code.
//
// Ian Kaplan, October 1996
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : Math-related functions
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: math.c,v $
// CVS revision  : $Revision: 1.12 $
// Last modified : $Date: 2013/06/07 18:55:17 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <arch.h>
#include <types.h>

// ===========================================================================
// ar_uint_divide()             Perform unsigned integer division to compute
//                              a quotient and/or a remainder
// ===========================================================================
// * INPUTS
//   unsigned int num           The number to be divided (dividend)
//   unsigned int div           The divisor
//
// * OUTPUTS
//   unsigned int *quot         If not NULL, the quotient is returned here
//   unsigned int *rem          If not NULL, the remainder is returned here
// ===========================================================================
/*
 * The code below uses the NEON FPU to do a fast division and then tries
 * to fix inaccuracies. It won't work for some extreme cases (when the ratio of
 * divisor/dividend is too big), so either debug it or forget about it.
 *
void ar_uint_divide(unsigned int num, unsigned int div,
                    unsigned int *quot, unsigned int *rem) {

  float         fl_quot;
  unsigned int  int_quot;
  unsigned int  int_rem;
  unsigned int  cur;


  // Sanity check
  if (!div) {
    ar_panic("Division by zero");
  }

  // ARM does not have integer division capabilities. We'll do this with
  // floating point instructions, but we have to correct for any possible
  // rounding problems or floating point accuracy errors.
  fl_quot = ((float) num) / ((float) div);

  // Start with cutting off decimal digits
  int_quot = fl_quot;

  // Loop until we get this right
  while (1) {

    cur = int_quot * div;

    // Perfect division
    if (cur == num) {
      int_rem = 0;
      break;
    }

    // Accuracy error overshot this. Decrement quotient.
    if (cur > num) {
      int_quot--;
      continue;
    }

    // cur < num. Check remainder.
    int_rem = num - cur;
    if (int_rem < div) {
      break;
    }
    int_quot++;
  }

  // Return results
  if (quot) {
    *quot = int_quot;
  }
  if (rem) {
    *rem = int_rem;
  }
}
*/


// ===========================================================================
// ar_uint_divide()             Perform unsigned integer division to compute
//                              a quotient and/or a remainder
// ===========================================================================
// * INPUTS
//   unsigned int num           The number to be divided (dividend)
//   unsigned int div           The divisor
//
// * OUTPUTS
//   unsigned int *ret_quot     If not NULL, the quotient is returned here
//   unsigned int *ret_rem      If not NULL, the remainder is returned here
// ===========================================================================
void ar_uint_divide(unsigned int num, unsigned int div,
                    unsigned int *ret_quot, unsigned int *ret_rem) {
  unsigned int quot, rem;
  unsigned int t, num_bits;
  unsigned int q, bit, d;
  int i;

  rem = 0;
  quot = 0;
  d = 0;

  if (div == 0) {
    ar_panic("Division by zero");
  }

  if (div > num) {
    rem = num;
    goto end;
  }

  if (div == num) {
    quot = 1;
    goto end;
  }

  num_bits = 32;

  while (rem < div) {
    bit = (num & 0x80000000) >> 31;
    rem = (rem << 1) | bit;
    d = num;
    num = num << 1;
    num_bits--;
  }


  /* The loop, above, always goes one iteration too far.
     To avoid inserting an "if" statement inside the loop
     the last iteration is simply reversed. */

  num = d;
  rem = rem >> 1;
  num_bits++;

  for (i = 0; i < num_bits; i++) {
    bit = (num & 0x80000000) >> 31;
    rem = (rem << 1) | bit;
    t = rem - div;
    q = !((t & 0x80000000) >> 31);
    num = num << 1;
    quot = (quot << 1) | q;
    if (q) {
       rem = t;
     }
  }

end:
  if (ret_quot) {
    *ret_quot = quot;
  }
  if (ret_rem) {
    *ret_rem = rem;
  }
}


// ===========================================================================
// ar_int_divide()              Perform signed integer division to compute
//                              a quotient and/or a remainder
// ===========================================================================
// * INPUTS
//   int num                    The number to be divided (dividend)
//   int div                    The divisor
//
// * OUTPUTS
//   int *ret_quot              If not NULL, the quotient is returned here
//   int *ret_rem               If not NULL, the remainder is returned here
// ===========================================================================
void ar_int_divide(int num, int div, int *ret_quot, int *ret_rem) {
  char inv_quot;
  unsigned int unum, udiv, uquot;
  int final_quot;

  // Track sign changes
  inv_quot = 0;

  // Convert num to positive
  if (num < 0) {
    inv_quot = 1;
    unum = -num;
  }
  else {
    unum = num;
  }

  // Convert div to positive
  if (div < 0) {
    inv_quot = 1 - inv_quot;
    udiv = -div;
  }
  else {
    udiv = div;
  }

  // Do the unsigned division
  ar_uint_divide(unum, udiv, &uquot, NULL);

  // Inverse quotient
  final_quot = (inv_quot) ? -uquot : uquot;

  // Return results
  if (ret_quot) {
    *ret_quot = final_quot;
  }
  if (ret_rem) {
    // Remainder is tricky with negative numbers. It must always satisfy
    // the expression
    //
    //          (quotient * divisor) + remainder = dividend
    //
    // which makes its sign really counter-intuitive. E.g.:
    // (-7) / (-2) -> quotient = +3, remainder = -1
    //
    *ret_rem = num - (final_quot * div);
  }
}


// ===========================================================================
// __aeabi_idiv()               GCC integer division software wrappers
// __aeabi_idivmod()
// __aeabi_uidiv()
// __aeabi_uidivmod()
// ===========================================================================
// Quotient of signed division
int __aeabi_idiv(int a, int b) {
  int q;
  ar_int_divide(a, b, &q, NULL);
  return q;
}
// Remainder of signed division. ARM EABI demands the return of this struct:
//   typedef struct {
//     int quot;
//     int rem;
//   } idiv_return;
// passed directly to the registers (r0 and r1).
int __aeabi_idivmod(int a, int b) {
  int q, r;
  ar_int_divide(a, b, &q, &r);
  asm("mov r1, %0" : : "r"(r));
  return q;
}
// Quotient of unsigned division
unsigned int __aeabi_uidiv(unsigned int a, unsigned int b) {
  unsigned int q;
  ar_uint_divide(a, b, &q, NULL);
  return q;
}
// Remainder of unsigned division
// Remainder of signed division. ARM EABI demands the return of this struct:
//typedef struct {
//  unsigned int quot;
//  unsigned int rem;
//} uidiv_return;
// passed directly to the registers (r0 and r1).
unsigned int __aeabi_uidivmod(unsigned int a, unsigned int b) {
  unsigned int q, r;
  ar_uint_divide(a, b, &q, &r);
  asm("mov r1, %0" : : "r"(r));
  return q;
}


// ===========================================================================
// ar_float_sqrt()              Perform floating-point square root
// ===========================================================================
// * INPUT
//   float num                  The number
//
// * OUTPUT
//   float                      square root of num
// ===========================================================================
float ar_float_sqrt(float num) {
  float ret;
  asm("vmov s15, %0" : : "r"(num) );
  asm("vsqrt.f32 s14, s15");
  asm("vmov %0, s14" : "=r"(ret) );
  return ret;
}


// ===========================================================================
// ar_float_pow()               Perform floating-point raise to power
// ===========================================================================
// * INPUT
//   float x                    The base
//   float y                    The exponent
//
// * OUTPUT
//   float                      x ** y
// ===========================================================================
float ar_float_pow(float x, float y) {

  int   inverse;
  float tmp;
  float low;
  float high;
  float acc;
  float mid;

  if (x == 0.0f) {
    ar_assert(y > 0.0f);
    return 0.0f;
  }

  if (y < 0.0f) {
    inverse = 1;
    y = -y;
  }
  else {
    inverse = 0;
  }

  if (y >= 1.0f) {
    tmp = ar_float_pow(x, y / 2.0f);
    tmp *= tmp;
    if (inverse) {
      return 1.0f / tmp;
    }
    else {
      return tmp;
    }
  }

  low = 0.0f;
  high = 1.0f;
  tmp = ar_float_sqrt(x);
  acc = tmp;
  mid = high / 2.0f;

  while (((mid > y) ? (mid - y) : (y - mid)) > 0.0001){
    tmp = ar_float_sqrt(tmp);
    if (mid <= y) {
      low = mid;
      acc *= tmp;
    }
    else {
      high = mid;
      acc *= (1.0f / tmp);
    }
    mid = (low + high) / 2.0f;
  }

  if (inverse) {
    return 1.0f / acc;
  }
  else {
    return acc;
  }
}


// ===========================================================================
// ar_sin()                     Computes the sine of a number
// ===========================================================================
// * INPUTS
//   float x                    The number (in radians)
//
// * RETURN VALUE
//   float                      sin(x)
// ===========================================================================
float ar_sin(float x) {

  float ret;
  float xsq;
  float prod;
  int   inv;


  // Bring to 0...PI section
  if (x < 0.0f) {
    x += (((int) (-x / 6.283185f)) + 1) * 6.283185f;
  }
  else if (x > 6.283185f) {
    x -= ((int) (x / 6.283185f)) * 6.283185f;
  }

  // Bring to 0...PI/4 section, mirror by PI/4 in even quadrants, remember
  // if we need to invert the result
  if (x > 4.712388f) {
    x -= 4.712388f;
    x = 1.570796f - x;
    inv = 1;
  }
  else if (x > 3.141592f) {
    x -= 3.141592f;
    inv = 1;
  }
  else if (x > 1.570796f) {
    x -= 1.570796f;
    x = 1.570796f - x;
    inv = 0;
  }
  else {
    inv = 0;
  }

  // Do Taylor series
  ret = x;
  xsq = x * x;          // x^2
  prod = x * xsq;       // x^3
  ret -= prod / 6.0f;
  prod *= xsq;          // x^5
  ret += prod / 120.0f;
  prod *= xsq;          // x^7
  ret -= prod / 5040.0f;

  // Inverse result if needed and return
  if (inv) {
    ret = -ret;
  }
  return ret;
}


// ===========================================================================
// ar_cos()                     Computes the cosine of a number
// ===========================================================================
// * INPUTS
//   float x                    The number (in radians)
//
// * RETURN VALUE
//   float                      cos(x)
// ===========================================================================
float ar_cos(float x) {
  return ar_sin(x + 1.570796f);
}


// ===========================================================================
// ar_tan()                     Computes the tangent of a number
// ===========================================================================
// * INPUTS
//   float x                    The number (in radians)
//
// * RETURN VALUE
//   float                      tan(x)
// ===========================================================================
float ar_tan(float x) {
  float sin;
  float cos;
  sin = ar_sin(x);
  cos = ar_cos(x);
  ar_assert(cos != 0.0f);
  return sin / cos;
}
