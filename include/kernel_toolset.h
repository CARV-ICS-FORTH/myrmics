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
// The code in this file is derived from the MMA version 1.0 project, which 
// was licensed under the following copyright:
//
// Copyright (c) 2011, Lawrence Livermore National Security, LLC.
// Produced at the Lawrence Livermore National Laboratory
// Written by Spyros Lyberis <lymperis1@llnl.gov>
// LLNL-CODE-637218, OCEC-13-184
// All rights reserved.
// 
// This file is part of MMA, version 1.0. 
// For details, please see http://myrmics.com/download.php
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
// - Redistributions of source code must retain the above copyright notice,
//   this list of conditions and the disclaimer below.
// - Redistributions in binary form must reproduce the above copyright notice,
//   this list of conditions and the disclaimer (as noted below) in the
//   documentation and/or other materials provided with the distribution.
// - Neither the name of the LLNS/LLNL nor the names of its contributors may be
//   used to endorse or promote products derived from this software without
//   specific prior written permission.
// 
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
// THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
// INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
// (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
// LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
// THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// Additional BSD Notice
// 
// 1. This notice is required to be provided under our contract with the U.S.
//    Department of Energy (DOE). This work was produced at Lawrence Livermore
//    National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.
// 2. Neither the United States Government nor Lawrence Livermore National
//    Security, LLC nor any of their employees, makes any warranty, express or
//    implied, or assumes any liability or responsibility for the accuracy,
//    completeness, or usefulness of any information, apparatus, product, or
//    process disclosed, or represents that its use would not infringe
//    privately-owned rights.
// 3. Also, reference herein to any specific commercial products, process, or
//    services by trade name, trademark, manufacturer or otherwise does not
//    necessarily constitute or imply its endorsement, recommendation, or
//    favoring by the United States Government or Lawrence Livermore National
//    Security, LLC. The views and opinions of authors expressed herein do not
//    necessarily state or reflect those of the United States Government or
//    Lawrence Livermore National Security, LLC, and shall not be used for
//    advertising or product endorsement purposes.
//
// ==========================[ Static Information ]===========================
//
// Author        : Spyros Lyberis
// Abstract      : Utility functions for in-kernel use
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: kernel_toolset.h,v $
// CVS revision  : $Revision: 1.20 $
// Last modified : $Date: 2013/03/02 22:41:02 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#ifndef _KERNEL_TOOLSET_H
#define _KERNEL_TOOLSET_H

#include <stdarg.h>

#include <types.h>


// ===========================================================================
// String functions
// ===========================================================================
extern unsigned int kt_strlen(const char *s);
extern char *kt_strstr(const char *s1, const char *s2);
extern int kt_strcmp(const char *s1, const char *s2);
extern int kt_strncmp(const char *s1, const char *s2, int num_bytes);
extern char *kt_strcpy(char *dst, const char *src);
extern char *kt_strncpy(char *dst, const char *src, int num_bytes);
extern char *kt_strdup(const char *s);
extern char *kt_strchr(const char *s, char c);
extern void kt_bzero(void *buf, int num_bytes);
extern void *kt_memset(void *buf, char val, int num_bytes);
extern void *kt_memcpy(void *dst, const void *src, int num_bytes);
extern int kt_memcmp(const void *s1, const void *s2, int num_bytes);
extern int kt_atoi(const char *s);


// ===========================================================================
// Print functions
// ===========================================================================
extern int kt_vsprintf(char *buf, const char *fmt, va_list args);
extern int kt_vprintf(const char *format, va_list ap);
extern int kt_printf(const char *format, ...);
extern int kt_sprintf(char *buf, const char *format, ...);


// ===========================================================================
// ASCII85 encoding
// ===========================================================================
void kt_encode85(unsigned int *array, int length);


// ===========================================================================
// Kernel allocation functions
// ===========================================================================
extern void kt_mem_init();

extern void *kt_malloc(size_t size);
extern void *kt_zalloc(size_t size);
extern void kt_free(void *ptr);
extern void *kt_realloc(void *old_ptr, size_t new_size);


// ===========================================================================
// Trie data structure
// ===========================================================================

// Types
typedef struct TrieNodeStruct TrieNode;
struct TrieNodeStruct {

  TrieNode      *children[8];   // 8 children, 0 being the lowest in raw key
                                // and 7 being the highest.
  TrieNode      *parent;        // Parent node
  int           bitmap;         // Bitmap of children that are present. A
                                // value of 0, means the node is a leaf.
  void          *data;          // User data pointer (valid when no children)
  size_t        raw_key;        // Raw key (valid when no children)

};

typedef struct {

  int           msb;            // Most significant bit of key values (0...63)
  int           lsb;            // Least significant bit of key values (0...63)
  TrieNode      *roots;         // Array of roots (one root per possible
                                // number of leading triads of zeroes)
  int           num_roots;      // Number of roots
  TrieNode      *walk;          // Current trie walking position (leaf node)
  TrieNode      *walk_int;      // Same, but for internals function use. We
                                // use this so that user can safely walk
                                // and insert/delete stuff while walking.
  unsigned int  num_keys;       // Number of keys in the Trie

} Trie;


// Exported functions
Trie            *kt_alloc_trie(int msb, int lsb);
void            kt_free_trie(Trie *t, void (*free_data)(void *));
void            kt_reset_trie(Trie *t, void (*free_data)(void *));

int             kt_trie_insert(Trie *t, size_t raw_key, void *data);
int             kt_trie_delete(Trie *t, size_t raw_key,
                               void (*free_data)(void *));

size_t          kt_trie_find(Trie *t, size_t raw_key, void **ret_data);
size_t          kt_trie_find_approx(Trie *t, int greater, size_t raw_key,
                                    void **ret_data);
size_t          kt_trie_find_minmax(Trie *t, int max, void **ret_data);
size_t          kt_trie_find_next(Trie *t, int greater, void **ret_data);
int             kt_trie_update_data(Trie *t, size_t raw_key, void *new_data,
                                    void (*free_old_data)(void *));

unsigned int    kt_trie_size(Trie *t);

// Debugging functions
//void kt_trie_dbg_dotty(Trie *t, char *filename);


// ===========================================================================
// List data structure
// ===========================================================================

// Types
typedef struct ListNodeStruct ListNode;
struct ListNodeStruct {

  ListNode      *prev;          // Previous in list
  ListNode      *next;          // Next in list
  void          *data;          // User data pointer

};

typedef struct {

  ListNode      *head;          // Head of the list
  ListNode      *tail;          // Tail of the list
  unsigned int  num_nodes;      // Number of list nodes

} List;


// Exported functions
extern ListNode *kt_list_insert(List *l, void *data, ListNode *insert_ref,
                                int insert_after);
extern void kt_list_delete(List *l, ListNode *n, void (*free_data)(void *));
extern List *kt_alloc_list();
extern void kt_free_list(List *l, void (*free_data)(void *));

// Inlined functions
static inline ListNode *kt_list_head(List *l) {
  return (l->head);
}

static inline ListNode *kt_list_tail(List *l) {
  return (l->tail);
}

static inline ListNode *kt_list_next(ListNode *n) {
  return (n->next);
}

static inline ListNode *kt_list_prev(ListNode *n) {
  return (n->prev);
}

static inline unsigned int kt_list_size(List *l) {
  return (l->num_nodes);
}


// ===========================================================================
// Math helper macros
// ===========================================================================

// Log base 2
static inline int kt_int_log2(int n) {
  int tmp1 = n;
  int tmp2 = 0;

  while (tmp1 >>= 1) {
    tmp2++;
  }

  return tmp2;
}

// Log base 10
static inline int kt_int_log10(int n) {
  int tmp1 = n;
  int tmp2 = 0;

  while (tmp1 /= 10) {
    tmp2++;
  }

  return tmp2;
}

// Integer power
static inline int kt_int_pow(int base, int exp) {
  int tmp1;
  int tmp2 = 1;
  for (tmp1 = 0; tmp1 < exp; tmp1++) {
    tmp2 *= base;
  }
  return tmp2;
}


// Random number generator: call with returned result as next seed to get
// up to 2^32 different numbers.
static inline unsigned int kt_rand(unsigned int seed) {
  return (1103515245 * seed + 12345);
}

#endif
