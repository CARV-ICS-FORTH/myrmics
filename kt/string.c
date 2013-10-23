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
// Abstract      : String-related functions
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: string.c,v $
// CVS revision  : $Revision: 1.12 $
// Last modified : $Date: 2013/04/09 14:41:41 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <kernel_toolset.h>
#include <types.h>

// ===========================================================================
// kt_isspace()                 Checks if the given character is a space
//                              (' ', '\t', '\v', '\f', '\r', '\n')
// ===========================================================================
// * INPUTS
//   const char c               The character
//
// * RETURN VALUE
//   int                        1 if it is a space character, else 0
// ===========================================================================
int kt_isspace(const char c) {
  switch (c) {
  case ' ':
  case '\t':
  case '\v':
  case '\f':
  case '\r':
  case '\n':
    return 1;
  default:
    return 0;
  }
}

// ===========================================================================
// kt_atoi()                    Converts the initial portion of the string
//                              to int. Works for decimal octal and hex.
// ===========================================================================
// * INPUTS
//   const char *s              The string
//
// * RETURN VALUE
//   int                        The converted value
// ===========================================================================
int kt_atoi(const char *s) {
  int res    = 0;
  int sign   = 1;
  int base   = 10;
  char limit = '0';

  /* Ignore white space */
  while ( kt_isspace(*s) )
    s++;

  if ( *s=='-' ) {
    sign = -1;
    s++;
  } else if ( *s=='+' ) {
    s++;
  }

  /* Check sign */
  if ( *s=='0' ) {
    s++;
    base = 8;
    if ( *s=='x' ) {
      base <<= 1; /* base = 16 */
      s++;
    }
  }

  /* Do the conversion */
  if ( base==16 ) {
    while ( (*s>='0' && *s<='9') ||
            (*s>='a' && *s<='f') ||
            (*s>='A' && *s<='F') ) {
      if (*s<='9')
        res = res*base+(*s)-'0';
      else if (*s<='F')
        res = res*base+10+(*s)-'A';
      else if (*s<='f')
        res = res*base+10+(*s)-'a';
      s++;
    }
  } else {
    limit+=base;
    while ( *s>='0' && *s<=limit ) {
      res = res*base+(*s)-'0';
      s++;
    }
  }

  return sign*res;
}

// ===========================================================================
// kt_strlen()                  Returns length of string
// ===========================================================================
// * INPUTS
//   const char *s              The string
//
// * RETURN VALUE
//   int                        Its length
// ===========================================================================
unsigned int kt_strlen(const char *s) {
  const char *c;

  for (c = s; *c; c++)
    ;

  return (c - s);
}

// ===========================================================================
// kt_strstr()                  Finds the first occurrence of a string which
//                              appears in another string
// ===========================================================================
// * INPUTS
//   const char *s1             The first string
//   const char *s2             The second string
//
// * RETURN VALUE
//   char *                     pointer to the first occurence of s2 in s1
//                              NULL if not match
// ===========================================================================
char *kt_strstr(const char *s1, const char *s2)
{
  char *cp = (char *) s1;
  char *str1, *str2;

  if (!*s2) return (char *) s1;

  while (*cp)
  {
    str1 = cp;
    str2 = (char *) s2;

    while (*str1 && *str2 && !(*str1 - *str2)) str1++, str2++;
    if (!*str2) return cp;
    cp++;
  }

  return NULL;
}


// ===========================================================================
// kt_strcmp()                  Compares two NULL-terminated strings for
//                              equality
// ===========================================================================
// * INPUTS
//   const char *s1             The first string
//   const char *s2             The second string
//
// * RETURN VALUE
//   int                         0: strings are equal
//                              >0: s1 lexicographically precedes s2
//                              <0: s1 lexicographically succeeds s2
// ===========================================================================
int kt_strcmp(const char *s1, const char *s2) {

  for ( ; *s1 == *s2; s1++, s2++) {
    if (!(*s1)) {
      return 0;
    }
  }

  return *(const unsigned char *)s1 - *(const unsigned char *)s2;
}


// ===========================================================================
// kt_strncmp()                 Similar to kt_strcmp, except it only compares
//                              the first at most n characters of the two
//                              NULL-terminated strings for
// ===========================================================================
// * INPUTS
//   const char *s1             The first string
//   const char *s2             The second string
//   int num_bytes              The number of elements to be compared
//
//
// * RETURN VALUE
//   int                         0: strings are equal
//                              >0: s1 lexicographically precedes s2
//                              <0: s1 lexicographically succeeds s2
// ===========================================================================
int kt_strncmp(const char *s1, const char *s2, int num_bytes) {
  if (num_bytes == 0)
    return 0;

  for ( ; *s1 == *s2 && num_bytes > 1; s1++, s2++, num_bytes--) {
    if (!(*s1)) {
      return 0;
    }
  }

  return *(const unsigned char *)s1 - *(const unsigned char *)s2;
}


// ===========================================================================
// kt_strcpy()                  Copies a string from a source to a
//                              destination buffer
// ===========================================================================
// * INPUTS
//   const char *src            The source string
//
// * OUTPUTS
//   const char *dst            The destination string
//
// * RETURN VALUE
//   char *                     Returns the destination string
// ===========================================================================
char *kt_strcpy(char *dst, const char *src) {
  char *ret = dst;

  while ((*dst++ = *src++))
    ;

  return ret;
}


// ===========================================================================
// kt_strncpy()                 Copies a string from a source to a
//                              destination buffer. Bytes of transfer are
//                              capped.
// ===========================================================================
// * INPUTS
//   const char *src            The source string
//   int num_bytes              Maximum number of bytes to copy
//
// * OUTPUTS
//   const char *dst            The destination string
//
// * RETURN VALUE
//   char *                     Returns the destination string
// ===========================================================================
char *kt_strncpy(char *dst, const char *src, int num_bytes) {
  char *ret = dst;

  while ((*dst++ = *src++) && (num_bytes--))
    ;

  return ret;
}


// ===========================================================================
// kt_strdup()                  Allocates memory and duplicates string
// ===========================================================================
// * INPUTS
//   const char *src            The source string
//
// * OUTPUTS
//
// * RETURN VALUE
//   char *                     Returns the destination string
// ===========================================================================
char *kt_strdup(const char *s)
{
  char *t;
  int len;

  if (!s) return NULL;
  len = kt_strlen(s);
  t = (char *) kt_malloc(len + 1);
  kt_memcpy(t, s, len + 1);
  return t;
}



// ===========================================================================
// kt_strchr()                  Returns a pointer to the first occurence
//                              of a character in a string
// ===========================================================================
// * INPUTS
//   const char *s              The string
//
// * RETURN VALUE
//   int                        Its length
// ===========================================================================
char *kt_strchr(const char *s, char c) {
  const char *r;

  for (r = s; *r && (*r != c); r++)
    ;

  return (char *) ((*r) ? r : NULL);
}


// ===========================================================================
// kt_bzero()                   Fills a buffer with 0
// ===========================================================================
// * INPUTS
//   void *buf                  The buffer to be filled
//   unsigned int num_bytes     Number of bytes to fill
// ===========================================================================
void kt_bzero(void *buf, int num_bytes) {
  unsigned int ptr;
  int i;
  unsigned char offset;

  ptr = (unsigned int) buf;

  // Do byte-by-byte until we're aligned to 4B pointer
  offset = ptr & 0x3;
  i = 0;
  if (offset) {
    for ( ; i < num_bytes && i < (4 - offset); i++) {
      *((char *) ptr++) = 0;
    }
  }

  // Do as many words as possible
  for ( ; i < num_bytes - 3; i += 4) {
    *((unsigned int *) ptr) = 0;
    ptr += 4;
  }

  // Do the rest in bytes
  for ( ; i < num_bytes; i++) {
    *((char *) ptr++) = 0;
  }
}


// ===========================================================================
// kt_memset()                  Fills a buffer with a constant byte value
// ===========================================================================
// * INPUTS
//   void *buf                  The buffer to be filled
//   char val                   Byte value to be written
//   unsigned int num_bytes     Number of bytes to fill
//
// * RETURN VALUE
//   void *                     Returns the buffer
// ===========================================================================
void *kt_memset(void *buf, char val, int num_bytes) {
  unsigned int ptr;
  unsigned int val4;
  int i;
  unsigned char offset;

  val4 = (val << 24) | (val << 16) | (val << 8) | val;
  ptr = (unsigned int) buf;

  // Do byte-by-byte until we're aligned to 4B pointer
  offset = ptr & 0x3;
  i = 0;
  if (offset) {
    for ( ; i < num_bytes && i < (4 - offset); i++) {
      *((char *) ptr++) = val;
    }
  }

  // Do as many words as possible
  for ( ; i < num_bytes - 3; i += 4) {
    *((unsigned int *) ptr) = val4;
    ptr += 4;
  }

  // Do the rest in bytes
  for ( ; i < num_bytes; i++) {
    *((char *) ptr++) = val;
  }

  return buf;
}


// ===========================================================================
// kt_memcpy()                  Copies data from a source to a destination
//                              buffer
// ===========================================================================
// * INPUTS
//   const void *src            The source buffer
//   unsigned int num_bytes     Number of bytes to copy
//
// * OUTPUTS
//   const void *dst            The destination buffer
//
// * RETURN VALUE
//   void *                     Returns the buffer
// ===========================================================================
void *kt_memcpy(void *dst, const void *src, int num_bytes) {
  unsigned int i;

  for (i = 0; i < num_bytes; i++) {
    ((char *) dst)[i] = ((char *) src)[i];
  }

  return dst;
}


// ===========================================================================
// kt_memmove()                 Moves data from a source to a destination
//                              buffer
// ===========================================================================
// * INPUTS
//   const void *src            The source buffer
//   unsigned int num_bytes     Number of bytes to move
//
// * OUTPUTS
//   const void *dst            The destination buffer
//
// * RETURN VALUE
//   void *                     Returns the buffer
// ===========================================================================
void *kt_memmove(void *dst_void, const void *src_void, int length) {
  char *dst = dst_void;
  const char *src = src_void;

  if (src < dst && dst < src + length) {
    /* Have to copy backwards */
    src += length;
    dst += length;
    while (length--)
      *--dst = *--src;
  } else
    while (length--)
      *dst++ = *src++;

  return dst_void;
}


// ===========================================================================
// kt_memcmp()                  Compares two buffers for equality
// ===========================================================================
// * INPUTS
//   const void *s1             The first buffer
//   const void *s2             The second buffer
//   unsigned int num_bytes     Number of bytes to compare
//
// * RETURN VALUE
//   int                        0 if buffers are equal, 1 otherwise
// ===========================================================================
int kt_memcmp(const void *s1, const void *s2, int num_bytes) {
  unsigned int i;

  for (i = 0; i < num_bytes; i++) {
    if (((char *) s1)[i] != ((char *) s2)[i]) {
      return 1;
    }
  }

  return 0;
}
