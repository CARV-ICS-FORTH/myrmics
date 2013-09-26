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
// Author        : Spyros LYBERIS
// Abstract      : Implementation of 8-degree Trie data structure. Supports
//                 approximate searches, as well as keys that only a part of
//                 their bits is meaningful.
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: trie.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/09/27 14:14:10 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>


// ===========================================================================
// Helper, inlined functions and macros
// ===========================================================================
#define kt_trie_is_leaf(n)      ((n)->bitmap == 0)

// Convert raw key to cooked, by chopping MSBs above what we care and shifting
// until LSB is on 0
static inline size_t kt_trie_cook_key(Trie *t, size_t raw_key) {
  size_t cooked_key;

  if (t->msb == ((sizeof(size_t) * 8) - 1)) {
    cooked_key = raw_key;
  }
  else {
    cooked_key = (1UL << (t->msb + 1)) - 1;
    cooked_key &= raw_key;
  }
  cooked_key >>= t->lsb;

  return cooked_key;
}

// Find appropriate root for cooked key
static inline int kt_trie_find_bin(Trie *t, size_t cooked_key) {
  int bin = 0;

  while (1) {
    cooked_key >>= 3;
    if (!cooked_key) {
      ar_assert(bin < t->num_roots);
      return t->num_roots - 1 - bin;
    }
    bin++;
  }
}

static inline TrieNode *kt_trie_find_leftmost_child(TrieNode *mid) {

  int i;
  int found;


  // Starting node should not be a leaf
  ar_assert(!kt_trie_is_leaf(mid));

  while (1) {

    // Mid-node cases
    if (!kt_trie_is_leaf(mid)) {
      found = 0;
      for (i = 0; i < 8; i++) {
        if (mid->bitmap & (1 << i)) {
          mid = mid->children[i];
          ar_assert(mid);
          found = 1;
          break;
        }
      }
      // A mid node should always have children
      if (!found) {
        ar_abort();
      }
    }
    // We reached a leaf
    else {
      return mid;
    }
  }
}

static inline TrieNode *kt_trie_find_rightmost_child(TrieNode *mid) {

  int i;
  int found;


  // Starting node should not be a leaf
  ar_assert(!kt_trie_is_leaf(mid));

  while (1) {

    // Mid-node cases
    if (!kt_trie_is_leaf(mid)) {
      found = 0;
      for (i = 7; i >= 0; i--) {
        if (mid->bitmap & (1 << i)) {
          mid = mid->children[i];
          ar_assert(mid);
          found = 1;
          break;
        }
      }
      // A mid node should always have children
      if (!found) {
        ar_abort();
      }
    }
    // We reached a leaf
    else {
      return mid;
    }
  }
}


// ===========================================================================
// kt_trie_insert()             Inserts a key into the Trie, creating as
//                              many mid- and leaf-nodes as needed.
// ===========================================================================
// * INPUTS
//   size_t raw_key             The raw key value to insert
//   void *data                 User data to hook to the new node for this key
//                              (can be NULL)
//
// * INOUTS
//   Trie *t                    The Trie
//
// * RETURN VALUE
//   int                        0: success
//                              1: key already exists
// ===========================================================================
int kt_trie_insert(Trie *t, size_t raw_key, void *data) {

  size_t        cooked_key;
  int           bin;
  TrieNode      *mid;
  TrieNode      *new_mid;
  TrieNode      *new_leaf;
  int           idx;
  int           i;


  // Sanity checks
  ar_assert(t);
  ar_assert(raw_key); // disallow raw NULL keys (cooked NULL is ok)

  // Cook key to chop of bits above MSB and below LSB
  cooked_key = kt_trie_cook_key(t, raw_key);

  // Fast track to correct root or leaf
  bin = kt_trie_find_bin(t, cooked_key);
  mid = &t->roots[bin];

  // Traverse tree, creating mid nodes as needed
  for (i = t->num_roots - 1 - bin; i >= 0; i--) {

    idx = (cooked_key >> (3 * i)) & 0x7;

    // Mid-node cases
    if (i > 0) {
      new_mid = mid->children[idx];
      if (!new_mid) {
        ar_assert(new_mid = kt_zalloc(sizeof(TrieNode)));
        mid->children[idx] = new_mid;
        mid->bitmap |= 1 << idx;
        new_mid->parent = mid;
      }
      else {
        ar_assert(!kt_trie_is_leaf(new_mid));
      }
      mid = new_mid;
    }
    // Leaf nodes
    else {
      new_leaf = mid->children[idx];
      if (!new_leaf) {
        ar_assert(new_leaf = kt_malloc(sizeof(TrieNode)));
        mid->children[idx] = new_leaf;
        mid->bitmap |= 1 << idx;
        new_leaf->parent = mid;
        new_leaf->bitmap = 0;
        new_leaf->data = data;
        new_leaf->raw_key = raw_key;
      }
      else {
        ar_assert(kt_trie_is_leaf(new_leaf));
        return 1;
      }
    }
  }

  // Increase counter
  t->num_keys++;

  return 0;
}


// ===========================================================================
// kt_trie_internal_find()      Searches the trie for a given key and returns
//                              the leaf node. This is the internal version
//                              of the function, used also for deletions.
//
//                              Note that a find can return a key which may
//                              not be exactly like the raw key it was told
//                              to lookup: bits above MSB or below LSB may
//                              be different. It is the user's responsibility
//                              to check this out.
//              
//                              This function saves the result internally
//                              for future use with find_next().
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   size_t raw_key             The raw key to be searched for
//
// * OUTPUTS
//   TrieNode **ret_abandon     If not NULL and the exact search fails,
//                              return the mid-node or root where the search
//                              was abandoned.
//   int *ret_idx               If not NULL and the exact search fails,
//                              return the index of the missing child.
//   
// * RETURN VALUE
//   TrieNode *                 The leaf node of the key, or NULL if the
//                              key does not exist
// ===========================================================================
TrieNode *kt_trie_internal_find(Trie *t, size_t raw_key, 
                                TrieNode **ret_abandon, int *ret_idx) {

  size_t        cooked_key;
  int           bin;
  int           idx;
  TrieNode      *mid;
  TrieNode      *new_mid;
  TrieNode      *leaf;
  int           i;


  // Sanity checks
  ar_assert(t);
  ar_assert(raw_key); // disallow raw NULL keys (cooked NULL is ok)

  // Cook key to chop of bits above MSB and below LSB
  cooked_key = kt_trie_cook_key(t, raw_key);

  // Fast track to correct root or leaf
  bin = kt_trie_find_bin(t, cooked_key);
  mid = &t->roots[bin];
  ar_assert(!mid->parent);

  // Traverse tree
  for (i = t->num_roots - 1 - bin; i >= 0; i--) {

    idx = (cooked_key >> (3 * i)) & 0x7;

    // Mid-node cases
    if (i > 0) {
      new_mid = mid->children[idx];
      if (!new_mid) {
        t->walk_int = NULL;
        if (ret_abandon) {
          *ret_abandon = mid;
        }
        if (ret_idx) {
          *ret_idx = idx;
        }
        return NULL;
      }
      ar_assert(!kt_trie_is_leaf(new_mid));
      mid = new_mid;
    }

    // We reached a leaf. Existence of a leaf means a key is present.
    else {
      leaf = mid->children[idx];
      if (!leaf) {
        t->walk_int = NULL;
        if (ret_abandon) {
          *ret_abandon = mid;
        }
        if (ret_idx) {
          *ret_idx = idx;
        }
        return NULL;
      }
      ar_assert(kt_trie_is_leaf(leaf));
      ar_assert(leaf->raw_key);
      t->walk_int = leaf;
      return leaf;
    }
  }

  // We should not be here
  ar_abort();
  return NULL;
}


// ===========================================================================
// kt_trie_find()               Searches the trie for a given key and returns
//                              the key. Just a user wrapper of
//                              kt_trie_internal_find().
//
//                              Note that a find can return a key which may
//                              not be exactly like the raw key it was told
//                              to lookup: bits above MSB or below LSB may
//                              be different. It is the user's responsibility
//                              to check this out.
//              
//                              This function saves the result internally
//                              for future use with kt_trie_find_next().
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   size_t raw_key             The raw key to be searched for
//
// * OUTPUTS
//   void **ret_data            If not NULL, returns the data associated with
//                              the search key (if found; otherwise, NULL is
//                              written to *ret_data)
//
// * RETURN VALUE
//   size_t                     0: key does not exist
//                              otherwise: the raw key that was found
// ===========================================================================
size_t kt_trie_find(Trie *t, size_t raw_key, void **ret_data) {

  TrieNode *n;

  // Lookup
  n = kt_trie_internal_find(t, raw_key, NULL, NULL);

  if (n) {

    // Should be a leaf node
    ar_assert(kt_trie_is_leaf(n));

    // Return data
    if (ret_data) {
      *ret_data = n->data;
    }

    // Return key
    return n->raw_key;
  }

  // Not found; return NULL and 0
  if (ret_data) {
    *ret_data = NULL;
  }
  return 0;
}


// ===========================================================================
// kt_trie_update_data()        Searches the trie for a given key and when
//                              it finds it, updates its data structure with
//                              a new one.
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   size_t raw_key             The raw key to be searched for
//   void *new_data             The new data to be associated with the key
//   void (*free_old_data)(void *)  If not NULL, and the node found has user
//                              data which is also not NULL, this function will
//                              be called to free the old user data
//
// * RETURN VALUE
//   int                        0: success
//                              1: key not found
// ===========================================================================
int kt_trie_update_data(Trie *t, size_t raw_key, void *new_data,
                        void (*free_old_data)(void *)) {

  TrieNode *n;


  // Lookup
  n = kt_trie_internal_find(t, raw_key, NULL, NULL);

  // Found?
  if (!n) {
    return 1;
  }

  // Should be a leaf node
  ar_assert(kt_trie_is_leaf(n));

  // Free old data?
  if (n->data && free_old_data) {
    free_old_data(n->data);
  }

  // Hook new data
  n->data = (void *) new_data;

  return 0;
}


// ===========================================================================
// kt_trie_internal_delete()    Deletes an already located leaf from a trie.
//                              This is an internal function, used also for
//                              freeing the whole trie.
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   TrieNode *leaf             The leaf to be deleted
//   void (*free_data)(void *)  If not NULL, and the node found has user data
//                              which is also not NULL, this function will
//                              be called to free the user data
//
// * RETURN VALUE
//   size_t                     0: success
//                              1: key did not exist
// ===========================================================================
int kt_trie_internal_delete(Trie *t, TrieNode *leaf, 
                            void (*free_data)(void *)) {

  TrieNode      *mid;
  TrieNode      *src;
  int           i;
  int           found;


  // Sanity check
  ar_assert(kt_trie_is_leaf(leaf));

  // User data?
  if (leaf->data) {
    if (free_data) {
      free_data(leaf->data);
    }
    leaf->data = NULL;
  }

  // Traverse upwards
  src = leaf;
  mid = leaf->parent; 
  while (1) {

    // Free the child
    kt_free(src);

    // Update parent
    found = 0;
    for (i = 0; i < 8; i++) {
      if (mid->children[i] == src) {
        mid->children[i] = NULL;
        mid->bitmap &= ~(1 << i);
        found = 1;
        break;
      }
    }
    if (!found) {
      // corruption: parent does not point to child
      ar_abort();
    }

    // We must continue freeing stuff until we either:
    // - reach a parent who has another child
    // - reach a root
    if (mid->bitmap || !mid->parent) {
      break;
    }

    // Continue
    src = mid;
    mid = mid->parent;
    ar_assert(!kt_trie_is_leaf(mid));
  }

  // Decrease counter
  t->num_keys--;

  return 0;
}
      

// ===========================================================================
// kt_trie_delete()             Searches the trie for a given key and deletes
//                              it. Mainly a user wrapper of
//                              kt_trie_internal_delete.
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   size_t raw_key             The raw key to be deleted
//   void (*free_data)(void *)  If not NULL, and the node found has user data
//                              which is also not NULL, this function will
//                              be called to free the user data
//
// * RETURN VALUE
//   size_t                     0: success
//                              1: key did not exist
// ===========================================================================
int kt_trie_delete(Trie *t, size_t raw_key, void (*free_data)(void *)) {

  TrieNode      *leaf;
  int           ret;

  // Sanity checks
  ar_assert(raw_key); // Disallow NULL raw keys

  // Find the leaf node
  leaf = kt_trie_internal_find(t, raw_key, NULL, NULL);
  if (!leaf) {
    return 1;
  }

  // Delete it
  ret = kt_trie_internal_delete(t, leaf, free_data);

  return ret;
}
      

// ===========================================================================
// kt_trie_find_minmax()        Returns the minimum or maximum key from the 
//                              trie and remembers it for future calls of 
//                              kt_trie_find_next().
// ===========================================================================
// * INOUTS
//   Trie *t                    The Trie
//   int max                    0: search for the minimum key
//                              1: search for the maximum key
//
// * OUTPUTS
//   void **ret_data            If not NULL, returns the data associated with
//                              the minimum key (if found; otherwise, NULL is
//                              written to *ret_data)
//
// * RETURN VALUE
//   size_t                     Minimum key of trie, or 0 if the trie is empty
// ===========================================================================
size_t kt_trie_find_minmax(Trie *t, int max, void **ret_data) {
  
  TrieNode   *mid;
  TrieNode   *leaf;
  int        i;
 

  // Sanity checks
  ar_assert(t);
  ar_assert((max == 0) || (max == 1));

  // Is it empty?
  if (!t->num_keys) {
    t->walk = NULL;
    if (ret_data) {
      *ret_data = NULL;
    }
    return 0;
  }

  // Locate first non-empty root, starting from:
  // - max == 0: the smallest value (as many leading zeroes as possible)
  // - max == 1: the biggest  value (as few leading zeroes as possible)
  mid  = NULL;
  i = (!max) ? (t->num_roots - 1) : 0;
  while (((!max) && (i >= 0)) ||
         (( max) && (i < t->num_roots))) {

    // Does root have any children?
    mid = &t->roots[i];
    ar_assert(!mid->parent);
    if (mid->bitmap) {
      break;
    }

    // Go to next one
    if (!max) {
      i--;
    }
    else {
      i++;
    }
  }

  if (!mid) {
    // This should not be reached: if t is empty, it's handled above.
    ar_abort();
  }
    
  // Find its leftmost (min) or rightmost (max) leaf and return it
  if (!max) {
    leaf = kt_trie_find_leftmost_child(mid);
  }
  else {
    leaf = kt_trie_find_rightmost_child(mid);
  }
  ar_assert(leaf);
  ar_assert(kt_trie_is_leaf(leaf));
  t->walk = leaf;
  
  // Return data
  if (ret_data) {
    *ret_data = leaf->data;
  }

  // Return key
  return leaf->raw_key;
}


// ===========================================================================
// kt_trie_internal_find_next() Returns the next key (bigger or smaller) 
//                              from a leaf or non-leaf starting point.
//                              Updates t->walk_int on what it found.
// ===========================================================================
// * INPUTS
//   TrieNode *src              Starting node from where the search begins
//   int greater                0: search for the next smaller key available
//                              1: search for the next greater key available
//   int src_idx                If src is a mid node, index of the child that
//                              is the start of this search (i.e. start from
//                              the next child of src_idx). If src is a leaf,
//                              src_idx is don't care.
// * INOUTS
//   Trie *t                    The Trie
//
// * OUTPUTS
//   void **ret_data            If not NULL, returns the data associated with
//                              the next key (if found; otherwise, NULL is
//                              written to *ret_data)
//
// * RETURN VALUE
//   size_t                     Next key of trie, or 0 if trie walking
//                              is finished (no next value exists)
// ===========================================================================
size_t kt_trie_internal_find_next(Trie *t, TrieNode *src, int greater, 
                                  int src_idx, void **ret_data) {
  TrieNode      *mid;
  int           found;
  int           i;


  // Sanity checks
  ar_assert(src);

  mid = NULL;
  // If we're in a mid node or a root, simply try to continue from the next
  // available child, without going up:
  // - greater == 1: take a righter child
  // - greater == 0: take a lefter  child
  if (!kt_trie_is_leaf(src)) {
    ar_assert((src_idx >= 0) && (src_idx < 8));
    if (greater) {
      for (i = src_idx; i < 8; i++) {
        if (src->children[i]) {
          mid = src->children[i];
          break;
        }
      }
    }
    else {
      for (i = src_idx; i >= 0; i--) {
        if (src->children[i]) {
          mid = src->children[i];
          break;
        }
      }
    }
  }

  // If we're at a leaf, or if there are no other children, go up
  if (!mid) {
    // Try to walk upwards from the current node (leaf or mid) to take:
    // - greater == 1: a right turn
    // - greater == 0: a left  turn
    while (1) {
      mid = src->parent;
      if (!mid) {
        break;
      }
      ar_assert(!kt_trie_is_leaf(mid));
      found = 0;
      if (greater) {
        // Try to find a child more to the right than us
        for (i = 0; i < 8; i++) {
          if (mid->children[i] == src) {
            ar_assert(!found);
            found = 1;
            continue;
          }
          if (found && mid->children[i]) {
            mid = mid->children[i];
            found = 2;
            break;
          }
        }
      }
      else {
        // Try to find a child more to the left than us
        for (i = 7; i >= 0; i--) {
          if (mid->children[i] == src) {
            ar_assert(!found);
            found = 1;
            continue;
          }
          if (found && mid->children[i]) {
            mid = mid->children[i];
            found = 2;
            break;
          }
        }
      }
      // Did we find such a child?
      ar_assert(found >= 1);
      if (found == 2) {
        break;
      }
      else {
        // If not, continue upwards
        src = mid;
      }
    }
  }

  // If we managed to take a correct turn above, return the:
  // - greater == 1: leftmost child of this subtree
  // - greater == 0: rightmost child of this subtree
  if (mid) {

    // If it's a leaf, just return it
    if (kt_trie_is_leaf(mid)) {
      t->walk_int = mid;

      // Return data
      if (ret_data) {
        *ret_data = mid->data;
      }

      // Return key
      return mid->raw_key;
    }

    // Mid-node: traverse down the subtree
    else {
      if (greater) {
        t->walk_int = kt_trie_find_leftmost_child(mid);
      }
      else {
        t->walk_int = kt_trie_find_rightmost_child(mid);
      }
      ar_assert(t->walk_int);
      ar_assert(kt_trie_is_leaf(t->walk_int));

      // Return data
      if (ret_data) {
        *ret_data = t->walk_int->data;
      }

      // Return key
      return t->walk_int->raw_key;
    }
  }


  // Otherwise, we reached a root. We must proceed to the next non-empty root.
  // We search linearly among the roots to find the next one. We could have
  // stored some info to do this quicker, but changing roots is a rare case 
  // and is very bounded anyway.
  //
  // Direction of search:
  // - greater == 1: search for a bigger  root (less leading zeroes)
  // - greater == 0: search for a smaller root (more leading zeroes)
  found = 0;
  i = (greater) ? (t->num_roots - 1) : 0;
  while ((( greater) && (i >= 0)) ||
         ((!greater) && (i < t->num_roots))) {

    ar_assert(!t->roots[i].parent);

    // Locate the current root
    if (&t->roots[i] == src) {
      ar_assert(!found);
      found = 1;
      goto next;
    }

    // Once located, search for the next non-empty one
    if (found) {

      // Is it non-empty?
      if (t->roots[i].bitmap) {

        // Take its leftmost or rightmost child
        if (greater) {
          t->walk_int = kt_trie_find_leftmost_child(&t->roots[i]);
        }
        else {
          t->walk_int = kt_trie_find_rightmost_child(&t->roots[i]);
        }
        ar_assert(t->walk_int);
        ar_assert(kt_trie_is_leaf(t->walk_int));

        // Return data
        if (ret_data) {
          *ret_data = t->walk_int->data;
        }

        // Return key
        return t->walk_int->raw_key;
      }
    }

next:

    if (greater) {
      i--;
    }
    else {
      i++;
    }
  }

  // We should have found src by now...
  ar_assert(found);

  // ... but still, we may have finished walking. If so, return nothing.
  if (ret_data) {
    *ret_data = NULL;
  }
  t->walk_int = NULL;
  return 0;
}


// ===========================================================================
// kt_trie_find_next()          Returns the next key (bigger or smaller) from
//                              the trie. The function must be called for the
//                              first time only after any other find function
//                              has been called and has successfully returned
//                              an element.
//
//                              This is a user wrapper for
//                              kt_trie_internal_find_next().
// ===========================================================================
// * INPUTS
//   int greater                0: search for the next smaller key available
//                              1: search for the next greater key available
//
// * INOUTS
//   Trie *t                    The Trie
//
// * OUTPUTS
//   void **ret_data            If not NULL, returns the data associated with
//                              the minimum key (if found; otherwise, NULL is
//                              written to *ret_data)
//
// * RETURN VALUE
//   size_t                     Minimum key of trie, or 0 if trie walking
//                              is finished (no greater value exists)
// ===========================================================================
size_t kt_trie_find_next(Trie *t, int greater, void **ret_data) {

  size_t key;

  // Sanity checks: this function is always called for a leaf node which
  //                must be the result of another find function
  ar_assert(t);
  ar_assert(t->walk);
  ar_assert(kt_trie_is_leaf(t->walk));

  // Call the internal function. Update user version (t->walk) with what the
  // node it found (or NULL).
  key = kt_trie_internal_find_next(t, t->walk, greater, -1, ret_data);
  t->walk = t->walk_int;

  return key;
}


// ===========================================================================
// kt_trie_find_approx()        Searches the trie for a given key and returns
//                              either the exact key (if found), or the
//                              immediately adjacent (smaller or bigger) one.
//              
//                              This function saves the result internally
//                              for future use with find_next().
// ===========================================================================
// * INPUTS
//   int greater                0: search for the next smaller key available
//                              1: search for the next greater key available
//
//   size_t raw_key             The raw key to be searched for
//
// * INOUTS
//   Trie *t                    The Trie
//
// * OUTPUTS
//   void **ret_data            If not NULL, returns the data associated with
//                              the minimum key (if found; otherwise, NULL is
//                              written to *ret_data)
//
// * RETURN VALUE
//   TrieNode *                 The leaf node of the key, or NULL if the
//                              key does not exist
// ===========================================================================
size_t kt_trie_find_approx(Trie *t, int greater, size_t raw_key, 
                           void **ret_data) {

  TrieNode      *node;
  TrieNode      *abandon_node;
  int           abandon_idx;
  size_t        key;


  // Search for an exact match
  node = kt_trie_internal_find(t, raw_key, &abandon_node, &abandon_idx);

  if (node) {

    // Should be a leaf node
    ar_assert(kt_trie_is_leaf(node));

    // Return data
    if (ret_data) {
      *ret_data = node->data;
    }

    // Remember node (user view) for find_next
    t->walk = t->walk_int;

    // Return key
    return node->raw_key;
  }

  // If not found, get the next available node (if any), from the point
  // where the exact search was abandoned. If we abandoned from a leaf
  // (actually a root leaf), go upwards (next root). If we abandoned from 
  // a mid node, just take the other turn.
  key = kt_trie_internal_find_next(t, abandon_node, greater, 
                                   abandon_idx, ret_data);
  t->walk = t->walk_int;
  return key;
}


// ===========================================================================
// kt_trie_size()               Returns the number of keys in the trie
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//
// * RETURN VALUE
//   unsigned int               Number of keys
// ===========================================================================
unsigned int kt_trie_size(Trie *t) {
  ar_assert(t);
  return (t->num_keys);
}



// ===========================================================================
// kt_alloc_trie()              Allocates a new Trie structure
// ===========================================================================
// * INPUTS
//   int msb                    Most significant bit position of Trie keys
//   int lsb                    Least significant bit position of Trie keys
//
// * RETURN VALUE
//   Trie *                     The newly allocated Trie
// ===========================================================================
Trie *kt_alloc_trie(int msb, int lsb) {

  Trie  *t;


  // Allocate trie struct, everything zeroed out
  t = kt_zalloc(sizeof(Trie));
  ar_assert(t);

  // Remember key MSB/LSB values
  ar_assert(lsb <= msb);
  ar_assert((lsb >= 0) && (lsb < sizeof(size_t) * 8));
  ar_assert((msb >= 0) && (msb < sizeof(size_t) * 8));
  t->msb = msb;
  t->lsb = lsb;

  // Allocate roots. We need (msb - lsb + 1) / 3, rounded up (2 is added to
  // make the rounding up).
  ar_int_divide(msb - lsb + 3, 3, &(t->num_roots), NULL);
  ar_assert(t->roots = kt_zalloc(t->num_roots * sizeof(TrieNode)));


  return t;
}


// ===========================================================================
// kt_reset_trie()              Deallocates all the trie nodes, but keeps the
//                              (now empty) trie intact.
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   void (*free_data)(void *)  If not NULL, and for any nodes found that have 
//                              user data which are also not NULL, this 
//                              function will be called to free these data
// ===========================================================================
void kt_reset_trie(Trie *t, void (*free_data)(void *)) {

  TrieNode      *prv;
  size_t        key;


  // Sanity checks
  ar_assert(t);

  // Do we have nodes to free?
  if (t->num_keys) {

    // Find first key
    ar_assert(kt_trie_find_minmax(t, 0, NULL));
    
    // Delete loop
    do {

      // Sanity check this node
      prv = t->walk;
      ar_assert(prv);
      ar_assert(kt_trie_is_leaf(prv));

      // Walk to next node
      key = kt_trie_find_next(t, 1, NULL);

      // Free the previous one
      kt_trie_internal_delete(t, prv, free_data);

    } while (key);

    // We must have nothing left
    ar_assert(!t->walk);
    ar_assert(!t->num_keys);
  }
}


// ===========================================================================
// kt_free_trie()               Deallocates the trie, along with all its
//                              internal nodes.
// ===========================================================================
// * INPUTS
//   Trie *t                    The Trie
//   void (*free_data)(void *)  If not NULL, and for any nodes found that have 
//                              user data which are also not NULL, this 
//                              function will be called to free these data
// ===========================================================================
void kt_free_trie(Trie *t, void (*free_data)(void *)) {

  TrieNode      *prv;
  size_t        key;


  // Sanity checks
  ar_assert(t);

  // Do we have nodes to free?
  if (t->num_keys) {

    // Find first key
    ar_assert(kt_trie_find_minmax(t, 0, NULL));
    
    // Delete loop
    do {

      // Sanity check this node
      prv = t->walk;
      ar_assert(prv);
      ar_assert(kt_trie_is_leaf(prv));

      // Walk to next node
      key = kt_trie_find_next(t, 1, NULL);

      // Free the previous one
      kt_trie_internal_delete(t, prv, free_data);

    } while (key);

    // We must have nothing left
    ar_assert(!t->walk);
    ar_assert(!t->num_keys);
  }

  // Free the trie structure itself
  kt_free(t->roots);
  kt_free(t);
}


// ===========================================================================
// kt_trie_dbg_dotty()          Debug function to create a DOT graph file for 
//                              a Trie
// ===========================================================================
// * INPUTS
//   Trie *t                    The trie
//   char *filename             Filename for the graph output
// ===========================================================================
#if 0

#define KT_TRIE_DBG_STACK_SIZE 1024 // Debug stack size

void kt_trie_dbg_dotty(Trie *t, char *filename) {

  FILE *f;
  int i;
  TrieNode *stack[KT_TRIE_DBG_STACK_SIZE];
  int stack_size;
  TrieNode *node;
  int found;


  // Open file
  ar_assert(f = fopen(filename, "w"));

  // Print out header
  fprintf(f, 
"digraph \"Trie\" {\n\n"
"  root [shape=record, label=\"{MSB = %d|LSB = %d}\"];\n"
"\n", 
          t->msb, t->lsb);

  // Create a stack with the roots and link them to trie "root" 
  for (i = 0; i < t->num_roots; i++) {
    stack[i] = &t->roots[i];
    fprintf(f,
"  root->node%p [label=%d];\n",
            &t->roots[i], i);
  }
  stack_size = i;

  

  // Iterate on the stack
  while (stack_size) {

    // Pop
    node = stack[--stack_size];
  
    // Do the node
    if (!node->parent) {
      // root
      fprintf(f,
"  node%p [shape=box, style=filled, fillcolor=red, label=\"\"];\n",
              node);
    }

    else {
      if (kt_trie_is_leaf(node)) {
        // leaf
        fprintf(f,
"  node%p [shape=oval, style=filled, fillcolor=palegreen, label=\"0x%lX\"];\n",
                node, kt_trie_cook_key(t, node->raw_key));
      }
      else {
      // mid
        ar_assert(node->bitmap);
        fprintf(f,
"  node%p [shape=box, style=filled, fillcolor=green, label=\"\"];\n",
                node);
      }
    }


    // Add its children to the stack and link to parent
    if (!kt_trie_is_leaf(node)) {
      found = 0;
      for (i = 0; i < 8; i++) {
        if (node->children[i]) {
          stack[stack_size++] = node->children[i];
          ar_assert(stack_size < KT_TRIE_DBG_STACK_SIZE);
          fprintf(f,
"  node%p -> node%p [label=%d, headport=n, tailport=s];\n", 
                  node, node->children[i], i);
          found = 1;
        }
      }
      ar_assert(found);
    }
  }

  fprintf(f, "}\n");
  fflush(f);
  fclose(f);
}

#endif
