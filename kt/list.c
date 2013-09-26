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
// Author        : Spyros LYBERIS
// Abstract      : Generic double-linked list data structure
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: list.c,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2013/03/14 12:56:51 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <kernel_toolset.h>
#include <arch.h>


// ===========================================================================
// kt_list_insert()             Creates a new list node and inserts it after
//                              or before a reference node
// ===========================================================================
// * INPUTS
//   void *data                 User data to hook to the new node (can be NULL)
//   ListNode *insert_ref       List node to be used as reference for the
//                              insertion position (can be NULL only if the
//                              list is empty)
//   int insert_after           0: Create new node before insert_ref
//                              1: Create new node after insert_ref
//
// * INOUTS
//   List *l                    The List
//
// * RETURN VALUE
//   ListNode *                 The newly created list node
// ===========================================================================
ListNode *kt_list_insert(List *l, void *data, ListNode *insert_ref,     
                         int insert_after) {

  ListNode      *n;


  // Sanity checks
  ar_assert(l);
  ar_assert(insert_ref || (!l->num_nodes && !l->head && !l->tail));

  // Create new node
  n = kt_malloc(sizeof(ListNode));

  // Assign user data
  n->data = data;

  // First node?
  if (!l->num_nodes) {
    n->prev = NULL;
    n->next = NULL;

    l->head = n;
    l->tail = n;
  }

  // Insert after insert_ref
  else if (insert_after) {

    n->prev = insert_ref;
    n->next = insert_ref->next;

    insert_ref->next = n;

    if (insert_ref == l->tail) {
      l->tail = n;
      ar_assert(!n->next);
    }
    else {
      ar_assert(n->next);
      ar_assert(n->next->prev == insert_ref);
      n->next->prev = n;
    }
  }
  
  // Insert before insert_ref
  else {

    n->prev = insert_ref->prev;
    n->next = insert_ref;

    insert_ref->prev = n;
    
    if (insert_ref == l->head) {
      l->head = n;
      ar_assert(!n->prev);
    }
    else {
      ar_assert(n->prev);
      ar_assert(n->prev->next == insert_ref);
      n->prev->next = n;
    }
  }

  // Count new node
  l->num_nodes++;

  // Success
  return n;
}


// ===========================================================================
// kt_list_delete()             Deletes a node from a list
// ===========================================================================
// * INPUTS
//   ListNode *n                The node to be deleted
//   void (*free_data)(void *)  If not NULL, and the node found has user data
//                              which is also not NULL, this function will
//                              be called to free the user data
// * INOUTS
//   List *l                    The List
// ===========================================================================
void kt_list_delete(List *l, ListNode *n, void (*free_data)(void *)) {
  
  // Sanity checks
  ar_assert(l);
  ar_assert(n);

  // Update previous node
  if (n->prev) {
    n->prev->next = n->next;
  }
  else {
    ar_assert(l->head == n);
    l->head = n->next;
  }

  // Update next node
  if (n->next) {
    n->next->prev = n->prev;
  }
  else {
    ar_assert(l->tail == n);
    l->tail = n->prev;
  }


  // Free user data
  if (n->data && free_data) {
    free_data(n->data);
  }

  // Free this node
  kt_free(n);


  // Count it
  l->num_nodes--;
}


// ===========================================================================
// kt_alloc_list()              Allocates a new List structure
// ===========================================================================
// * RETURN VALUE
//   List *                     The newly allocated List
// ===========================================================================
List *kt_alloc_list() {

  List *l;

  // Create list with all fields zeroed
  l = kt_zalloc(sizeof(List));

  return l;
}


// ===========================================================================
// kt_free_list()               Deallocates the list along with all its nodes
// ===========================================================================
// * INPUTS
//   List *l                    The List
//   void (*free_data)(void *)  If not NULL, and for any nodes found that have 
//                              user data which are also not NULL, this 
//                              function will be called to free these data
// ===========================================================================
void kt_free_list(List *l, void (*free_data)(void *)) {

  ListNode *n;


  // Sanity checks
  ar_assert(l);

  // For all nodes
  for (n = l->head; n; n = n->next) {

    // Free user data
    if (n->data && free_data) {
      free_data(n->data);
    }

    // Free list node
    kt_free(n);
  }

  // Free list structure
  kt_free(l);
}


