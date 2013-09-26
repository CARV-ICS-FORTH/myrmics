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
// Various Myrmics test routines -- Code by Spyros Lyberis
// ============================================================================



#include <myrmics.h>

// ###########################################################################
// ###########################################################################
// Delayed task spawning speed test
// ###########################################################################
// ###########################################################################
#if 1

// ===========================================================================
// ===========================================================================
void test_myrmics_finish(rid_t r, unsigned int time_start, int delay) {

  unsigned int time_stop;
  unsigned int time;


  // Compute elapsed times
  time_stop = sys_free_timer_get_ticks();

  if (time_stop > time_start) {
    time = time_stop - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_stop);
  }
  printf("%3d workers, %10d delay: time %10u cycles\r\n", 
         sys_get_num_workers(), delay, time);
}


// ===========================================================================
// ===========================================================================
void test_myrmics_task(int *obj, int delay) {
  sys_timer_busy_wait_cycles(delay);
}


// ===========================================================================
// ===========================================================================
void test_myrmics(int num_tasks, int delay) {

  rid_t         r;
  int           **obj;
  int           *scoop_obj;
  unsigned int  time_start;
  unsigned int  time_spawn;
  int           i;


  // Allocate objects in a region
  r = sys_ralloc(0, 99);
  obj = sys_alloc(num_tasks * sizeof(int *), r);
  sys_balloc(sizeof(int), r, num_tasks, obj);

  // Start time
  time_start = sys_free_timer_get_ticks();

  // Spawn tasks
  for (i = 0; i < num_tasks; i++) {
    scoop_obj = obj[i];
    #pragma myrmics task inout(scoop_obj) in(delay)
    test_myrmics_task(scoop_obj, delay);
  }

  // Wait for them
  #pragma myrmics task region in(r) in(time_start, delay)
  test_myrmics_finish(r, time_start, delay);
}
#endif



// ###########################################################################
// ###########################################################################
// Task spawning speed test
// ###########################################################################
// ###########################################################################
#if 0

// ===========================================================================
// ===========================================================================
void test_myrmics_finish(int *obj, unsigned int time_start,
                         unsigned int time_spawn) {

  unsigned int time_stop;
  unsigned int time;


  // Compute elapsed times
  time_stop = sys_free_timer_get_ticks();

  if (time_spawn > time_start) {
    time = time_spawn - time_start;
  }
  else {
    time = 0xFFFFFFFF - (time_start - time_spawn);
  }
  printf("Time to spawn:   %10u cycles (%6u msec)\r\n", time, time / 10000);

  if (time_stop > time_spawn) {
    time = time_stop - time_spawn;
  }
  else {
    time = 0xFFFFFFFF - (time_spawn - time_stop);
  }
  printf("Time to execute: %10u cycles (%6u msec)\r\n", time, time / 10000);

  // See what's in there
  printf("Object value:    %10d\r\n", *obj);
}


// ===========================================================================
// ===========================================================================
void test_myrmics_task(int *obj) {
  (*obj)++;
}


// ===========================================================================
// ===========================================================================
void test_myrmics(int num_tasks) {

  int           *obj;
  unsigned int  time_start;
  unsigned int  time_spawn;
  int           i;


  // Allocate object
  obj = sys_alloc(sizeof(int), 0);
  *obj = 0;

  printf("Spawning %d tasks...\r\n", num_tasks);

  // Start time
  time_start = sys_free_timer_get_ticks();

  // Spawn tasks
  for (i = 0; i < num_tasks; i++) {
    #pragma myrmics task inout(obj)
    test_myrmics_task(obj);
  }

  // Time for spawning
  time_spawn = sys_free_timer_get_ticks();

  // Wait for them
  #pragma myrmics task in(obj, time_start, time_spawn)
  test_myrmics_finish(obj, time_start, time_spawn);

}

#endif



// ###########################################################################
// ###########################################################################
// Test for Jacobi read-only deadlock bug
// ###########################################################################
// ###########################################################################
#if 0

#define DELAY 50

// ===========================================================================
// ===========================================================================
void test_myrmics_finish(rid_t r2) {
  printf("finish: ok\r\n");
  sys_timer_busy_wait_msec(DELAY);
}


// ===========================================================================
// ===========================================================================
void test_myrmics_object_task(int *aux_in, int *obj_in, int *obj_out) {

  //printf("object_task: begin, in = 0x%X 0x%X, out = 0x%X\r\n", 
  //       aux_in, obj_in, obj_out);
  
  sys_timer_busy_wait_msec(DELAY);
  //printf("========\r\n");

  //printf("object_task: end\r\n");
}


// ===========================================================================
// ===========================================================================
void test_myrmics_region_task(int *aux_in, rid_t r_in, int *obj_in, 
                              rid_t r_out, int *obj_out) {

  //printf("region_task: begin, in = (0x%X) %d 0x%X, out = %d 0x%X\r\n", 
  //       aux_in, r_in, obj_in, r_out, obj_out);
  
  sys_timer_busy_wait_msec(DELAY);
  //printf("========\r\n");

  //printf("region_task: spawning in = 0x%X 0x%X, out = 0x%X\r\n", 
  //       aux_in, obj_in, obj_out);

  #pragma myrmics task in(aux_in, obj_in) inout(obj_out)
  test_myrmics_object_task(aux_in, obj_in, obj_out);

  sys_timer_busy_wait_msec(DELAY);
  //printf("========\r\n");

  //printf("region_task: end\r\n");

}


// ===========================================================================
// ===========================================================================
void test_myrmics() {

  rid_t r2;
  rid_t r3;
  rid_t r4;
  rid_t r5;
  rid_t r6;
  int   *o3a;
  int   *o3b;
  int   *o3c;
  int   *o4a;
  int   *o4b;
  int   *o4c;
  int   *o5a;
  int   *o5b;
  int   *o5c;
  int   *o6a;
  int   *o6b;
  int   *o6c;
  int   i;

  // all-holding
  r2 = sys_ralloc(0, 99); 

  // buffers
  r3 = sys_ralloc(r2, 0);
  r4 = sys_ralloc(r2, 0);
  r5 = sys_ralloc(r2, 0);
  r6 = sys_ralloc(r2, 0);

  // objects
  o3a = sys_alloc(sizeof(int), r3);
  o3b = sys_alloc(sizeof(int), r3);
  o3c = sys_alloc(sizeof(int), r3);
  o4a = sys_alloc(sizeof(int), r4);
  o4b = sys_alloc(sizeof(int), r4);
  o4c = sys_alloc(sizeof(int), r4);
  o5a = sys_alloc(sizeof(int), r5);
  o5b = sys_alloc(sizeof(int), r5);
  o5c = sys_alloc(sizeof(int), r5);
  o6a = sys_alloc(sizeof(int), r6);
  o6b = sys_alloc(sizeof(int), r6);
  o6c = sys_alloc(sizeof(int), r6);

  printf("r3 = %d || objects o3* = 0x%X 0x%X 0x%X\r\n", r3, o3a, o3b, o3c);
  printf("r4 = %d || objects o4* = 0x%X 0x%X 0x%X\r\n", r4, o4a, o4b, o4c);
  printf("r5 = %d || objects o5* = 0x%X 0x%X 0x%X\r\n", r5, o5a, o5b, o5c);
  printf("r6 = %d || objects o6* = 0x%X 0x%X 0x%X\r\n", r6, o6a, o6b, o6c);

  for (i = 0; i < 15; i++) {

    #pragma myrmics task in(o5a) region in(r3) in(o3a) safe(o3a) \
                         region inout(r4) in(o4a) safe(o4a)
    test_myrmics_region_task(o5a, r3, o3a, r4, o4a);


    #pragma myrmics task in(o3a) region in(r5) in(o5a) safe(o5a) \
                         region inout(r6) in(o6a) safe(o6a)
    test_myrmics_region_task(o3a, r5, o5a, r6, o6a);

    // ---

    #pragma myrmics task in(o6a) region in(r4) in(o4a) safe(o4a) \
                         region inout(r3) in(o3a) safe(o3a)
    test_myrmics_region_task(o6a, r4, o4a, r3, o3a);

    #pragma myrmics task in(o4a) region in(r6) in(o6a) safe(o6a) \
                         region inout(r5) in(o5a) safe(o5a)
    test_myrmics_region_task(o4a, r6, o6a, r5, o5a);
  }


  #pragma myrmics task region inout(r2)
  test_myrmics_finish(r2);

}

#endif



// ###########################################################################
// ###########################################################################
// Array test for Jacobi read-only corruption bug
// ###########################################################################
// ###########################################################################
#if 0

#define PARALLEL

// Fails in 5_HET_FLAT with printfs (sometimes)
//#define NUM_ARRAYS      65
//#define NUM_ELEMENTS    128
//#define STALL_MSEC      20

// Fails in 9_HET_FLAT with printfs
//#define NUM_ARRAYS      80
//#define NUM_ELEMENTS    512
//#define STALL_MSEC      0

// Fails in 9_HET_FLAT without printfs
#define NUM_ARRAYS      60
#define NUM_ELEMENTS    1024
#define STALL_MSEC      5

// ===========================================================================
// ===========================================================================
void test_myrmics_fill3(int *ar0, int fil0, int *ar1, int fil1, 
                        int *ar2, int fil2) {
  
  int i;
  int ar_pos;
  int mark0;
  int mark1;
  int mark2;

  // Stall a little bit
  sys_timer_busy_wait_msec(STALL_MSEC);

  // Check integrity of arrays
  mark0 = fil0 & 0xFFFF0000;
  fil0 &= 0xFFFF;
  mark1 = fil1 & 0xFFFF0000;
  fil1 &= 0xFFFF;
  mark2 = fil2 & 0xFFFF0000;
  fil2 &= 0xFFFF;

  for (i = 0; i < NUM_ELEMENTS; i++) {
    if (((i <  fil0) && (ar0[i] != (mark0 | i))) ||
        ((i >= fil0) && (ar0[i] != mark0))) {
      printf("%d %X %d Error %X %d %X\r\n", sys_get_worker_id(), ar0, fil0,
                mark0 | fil0, i, ar0[i]);
      sys_abort();
    }
    if (((i <  fil1) && (ar1[i] != (mark1 | i))) ||
        ((i >= fil1) && (ar1[i] != mark1))) {
      printf("%d %X %d Error %X %d %X\r\n", sys_get_worker_id(), ar1, fil1,
                mark1 | fil1, i, ar1[i]);
      sys_abort();
    }
    if (((i <  fil2) && (ar2[i] != (mark2 | i))) ||
        ((i >= fil2) && (ar2[i] != mark2))) {
      printf("%d %X %d Error %X %d %X\r\n", sys_get_worker_id(), ar2, fil2,
                mark2 | fil2, i, ar2[i]);
      sys_abort();
    }
  }

  // Fill new elements of ar2
  ar2[fil2] = mark2 | fil2;

  // Success
  //printf("%d %X %d %X %d %X %d\r\n", sys_get_worker_id(), ar0, fil0, ar1, fil1, ar2, fil2);

  if (fil0 % 50 == 0) {
    printf("%d at %d\r\n", sys_get_worker_id(), fil0);
  }
}


// ===========================================================================
// ===========================================================================
void test_myrmics_do_checksum(int *array, int filled, int *checksum) {
  int i;
  int mark;
 
  mark = filled & 0xFFFF0000;
  filled &= 0xFFFF;
  for (i = 0; i < NUM_ELEMENTS; i++) {
    if (((i <  filled) && (array[i] != (mark | i))) ||
        ((i >= filled) && (array[i] != mark))) {
      printf("* %d %X %d Error %X %d %X\r\n", 
                sys_get_worker_id(), array, filled,
                mark | filled, i, array[i]);
      sys_abort();
    }
    if (i < filled) {
      *checksum += array[i];
    }
  }
}


// ===========================================================================
// ===========================================================================
void test_myrmics_print_checksum(int *checksum, int exp_checksum) {
  if (*checksum == exp_checksum) {
    printf("Checksum = %d [ [1;32mPASS[0m ]\r\n", *checksum);
  }
  else {
    printf("Checksum = %d, expected %d [ [1;31mFAIL[0m ]\r\n", 
           *checksum, exp_checksum);
  }
}


// ===========================================================================
// ===========================================================================
void test_myrmics() {

  int           **arrays;
  int           *filled;
  int           *checksum;
  int           exp_checksum;
  unsigned      seed = 42;
  unsigned int  idx0;
  unsigned int  idx1;
  unsigned int  idx2;
  int           *ar0;
  int           *ar1;
  int           *ar2;
  int           fil0;
  int           fil1;
  int           fil2;
  int           mark;
  int           i;
  int           j;


  // Allocate arrays
  arrays = sys_alloc(NUM_ARRAYS * sizeof(int *), 0);
  filled = sys_alloc(NUM_ARRAYS * sizeof(int), 0);
  sys_balloc(NUM_ELEMENTS * sizeof(int), 0, NUM_ARRAYS, arrays);
  

  // Initialize them all to their fill mark, except pos 1 which is set to 1
  for (i = 0; i < NUM_ARRAYS; i++) {
    mark = (((int) arrays[i]) << 4) & 0xFFFF0000;
    for (j = 0; j < NUM_ELEMENTS; j++) {
      arrays[i][j] = (j == 1) ? (mark | 1) : mark;
    }
    filled[i] = mark | 2;
  }

  idx0 = 0;
  idx1 = 1;
  idx2 = NUM_ARRAYS / 2;
  while (1) {

    // We're done when one of the arrays is completely filled
    if (((filled[idx0] & 0xFFFF) == NUM_ELEMENTS) ||
        ((filled[idx1] & 0xFFFF) == NUM_ELEMENTS) ||
        ((filled[idx2] & 0xFFFF) == NUM_ELEMENTS)) {
      break;
    }

    ar0  = arrays[idx0];
    ar1  = arrays[idx1];
    ar2  = arrays[idx2];

    fil0 = filled[idx0];
    fil1 = filled[idx1];
    fil2 = filled[idx2];


#ifdef PARALLEL
    #pragma myrmics task in(ar0, fil0, ar1, fil1) inout(ar2) in(fil2)
#endif
    test_myrmics_fill3(ar0, fil0, ar1, fil1, ar2, fil2);
    filled[idx2]++;

    idx0 = (idx0 + 1) % NUM_ARRAYS;
    idx1 = (idx1 + 1) % NUM_ARRAYS;
    idx2 = (idx2 + 1) % NUM_ARRAYS;
    
    if (!idx0 && ((filled[idx0] & 0xFFFF) % 50 == 0)) {
      printf("Now at %d\r\n", filled[idx0] & 0xFFFF);
    }

  }

  // Checksum
  exp_checksum = 0;
  for (i = 0; i < NUM_ARRAYS; i++) {
    mark = (((int) arrays[i]) << 4) & 0xFFFF0000;
    for (j = 0; j < (filled[i] & 0xFFFF); j++) {
      exp_checksum += mark | j;
    }
  }

  checksum = sys_alloc(sizeof(int), 0);
  *checksum = 0;
  for (i = 0; i < NUM_ARRAYS; i++) {
    ar0 = arrays[i];
    fil0 = filled[i];

#ifdef PARALLEL
    #pragma myrmics task in(ar0, fil0) inout(checksum)
#endif
    test_myrmics_do_checksum(ar0, fil0, checksum);
  }

#ifdef PARALLEL
    #pragma myrmics task in(checksum, exp_checksum)
#endif
  test_myrmics_print_checksum(checksum, exp_checksum);

  printf("All spawns done.\r\n");
}
#endif




// ###########################################################################
// ###########################################################################
// Array test for Jacobi read-write corruption bug
// ###########################################################################
// ###########################################################################
#if 0

#define PARALLEL

#define NUM_ARRAYS      60
#define NUM_ELEMENTS    1024

#define STALL_MSEC      0

// Bug in cur 5_HET_FLAT (1 printf per fill + 1 per DMA)
//#define NUM_ARRAYS      127
//#define NUM_ELEMENTS    513
//
//#define STALL_MSEC      30

// Bug in cur 9_HET_FLAT (1 printf per fill + 1 per DMA)
//#define NUM_ARRAYS      60
//#define NUM_ELEMENTS    513
//
//#define STALL_MSEC      50

// Bug in cur 9_HET_FLAT (no printfs, with checksum integrity code)
//#define NUM_ARRAYS      60
//#define NUM_ELEMENTS    513
//
//#define STALL_MSEC      10

// Bug in cur 9_HET_FLAT (no printfs, without checksum array integrity code)
//#define NUM_ARRAYS      60
//#define NUM_ELEMENTS    512
//
//#define STALL_MSEC      10

// Bug in old 9_HET_FLAT (1 printf per fill + 1 per DMA)
//
//#define NUM_ARRAYS      60
//#define NUM_ELEMENTS    512
//
//#define STALL_MSEC      100


// ===========================================================================
// ===========================================================================
void test_myrmics_fill1(int *ar0, int fil0) {
  
  int i;
  int ar_pos;
  int mark;

  // Stall a little bit
  sys_timer_busy_wait_msec(STALL_MSEC);

  // Check integrity of arrays
  mark = fil0 & 0xFFFF0000;
  fil0 &= 0xFFFF;
  for (i = 0; i < NUM_ELEMENTS; i++) {
    if (((i <  fil0) && (ar0[i] != (mark | i))) ||
        ((i >= fil0) && (ar0[i] != mark))) {
      printf("%d %X %d Error %X %d %X\r\n", sys_get_worker_id(), ar0, fil0,
                mark | fil0, i, ar0[i]);
      sys_abort();
    }
  }

  // Fill new elements
  ar0[fil0] = mark | fil0;

  // Success
  //printf("%d %X %d\r\n", sys_get_worker_id(), ar0, fil0);

  if (fil0 % 50 == 0) {
    printf("%d at %d\r\n", sys_get_worker_id(), fil0);
  }
}


// ===========================================================================
// ===========================================================================
void test_myrmics_do_checksum(int *array, int filled, int *checksum) {
  int i;
  int mark;
 
  mark = filled & 0xFFFF0000;
  filled &= 0xFFFF;
  for (i = 0; i < NUM_ELEMENTS; i++) {
    if (((i <  filled) && (array[i] != (mark | i))) ||
        ((i >= filled) && (array[i] != mark))) {
      printf("* %d %X %d Error %X %d %X\r\n", 
                sys_get_worker_id(), array, filled,
                mark | filled, i, array[i]);
      sys_abort();
    }
    if (i < filled) {
      *checksum += array[i];
    }
  }
}


// ===========================================================================
// ===========================================================================
void test_myrmics_print_checksum(int *checksum, int exp_checksum) {
  if (*checksum == exp_checksum) {
    printf("Checksum = %d [ [1;32mPASS[0m ]\r\n", *checksum);
  }
  else {
    printf("Checksum = %d, expected %d [ [1;31mFAIL[0m ]\r\n", 
           *checksum, exp_checksum);
  }
}


// ===========================================================================
// ===========================================================================
void test_myrmics() {

  int           **arrays;
  int           *filled;
  int           *checksum;
  int           exp_checksum;
  unsigned      seed = 42;
  unsigned int  idx0;
  int           *ar0;
  int           fil0;
  int           mark;
  int           i;
  int           j;


  // Allocate arrays
  arrays = sys_alloc(NUM_ARRAYS * sizeof(int *), 0);
  filled = sys_alloc(NUM_ARRAYS * sizeof(int), 0);
  sys_balloc(NUM_ELEMENTS * sizeof(int), 0, NUM_ARRAYS, arrays);
  

  // Initialize them all to their fill mark, except pos 1 which is set to 1
  for (i = 0; i < NUM_ARRAYS; i++) {
    mark = (((int) arrays[i]) << 4) & 0xFFFF0000;
    for (j = 0; j < NUM_ELEMENTS; j++) {
      arrays[i][j] = (j == 1) ? (mark | 1) : mark;
    }
    filled[i] = mark | 2;
  }

  idx0 = 0;
  while (1) {

    // We're done when one of the arrays is completely filled
    if ((filled[idx0] & 0xFFFF) == NUM_ELEMENTS) {
      break;
    }

    ar0  = arrays[idx0];

    fil0 = filled[idx0];

#ifdef PARALLEL
    #pragma myrmics task inout(ar0) in(fil0)
#endif
    test_myrmics_fill1(ar0, fil0);
    
    filled[idx0]++;

    idx0 = (idx0 + 1) % NUM_ARRAYS;
    
    if (!idx0 && ((filled[idx0] & 0xFFFF) % 50 == 0)) {
      printf("Now at %d\r\n", filled[idx0] & 0xFFFF);
    }

  }

  // Checksum
  exp_checksum = 0;
  for (i = 0; i < NUM_ARRAYS; i++) {
    mark = (((int) arrays[i]) << 4) & 0xFFFF0000;
    for (j = 0; j < (filled[i] & 0xFFFF); j++) {
      exp_checksum += mark | j;
    }
  }

  checksum = sys_alloc(sizeof(int), 0);
  *checksum = 0;
  for (i = 0; i < NUM_ARRAYS; i++) {
    ar0 = arrays[i];
    fil0 = filled[i];

#ifdef PARALLEL
    #pragma myrmics task inout(ar0) in(fil0) inout(checksum)
#endif
    test_myrmics_do_checksum(ar0, fil0, checksum);
  }

#ifdef PARALLEL
    #pragma myrmics task inout(checksum) in(exp_checksum)
#endif
  test_myrmics_print_checksum(checksum, exp_checksum);

  printf("All spawns done.\r\n");
}

#endif




// ###########################################################################
// ###########################################################################
// Dummy functions test
// ###########################################################################
// ###########################################################################
#if 0

typedef struct NodeStruct Node;

struct NodeStruct {
  int   val;
  Node  *next;
};


void test_task1(rid_t r, Node *head) {

  Node  *new1;
  Node  *new2;

  new1 = sys_alloc(sizeof(Node), r);
  new2 = sys_alloc(sizeof(Node), r);

  head->next = new1;

  new1->val = 666;
  new1->next = new2;

  new2->val = 999;
  new2->next = NULL;
}

void test_task2(int v1, int v2, int v3, int *v4, int v5, int v6, int v7,
                int v8, int v9, int *v10) {
  printf("task2 %d %d %d %d %d %d %d %d\r\n", v1, v2, v3, v5, v6, v7, v8, v9);
  *v4 = v1 + v2 + v3;
  *v10 = v5 + v6 + v7 + v8 + v9;
}


void test_task3(int *foo, int safe) {

  *foo += safe;

  printf("task3: foo = 0x%p, safe = %d, *foo = %d\r\n", foo, safe, *foo);

  //sys_timer_busy_wait_cycles(1000000);

  //printf("*** after wait, *foo = %d\r\n", *foo);
}


void test_task4(rid_t r, int *ch) {
  printf("task4: r = %d\r\n", r);

  printf("%d: Freeing child 0x%08X\r\n", sys_get_worker_id(), ch);
  sys_free(ch);
  printf("%d: Done\r\n", sys_get_worker_id());

  //sys_timer_busy_wait_cycles(10000000); // 1.0 sec
}

void test_task5(int *foo, int safe) {
  int *objects[6];

  //printf("task3: foo = 0x%p, *foo = %d, safe = %d\r\n", foo, *foo, safe);
  *foo = safe + 1;

  sys_balloc(sizeof(int), 0, 6, (void *) objects);
  
  printf("0x%p 0x%p 0x%p 0x%p 0x%p 0x%p\r\n", objects[0], objects[1],
      objects[2], objects[3], objects[4], objects[5]);

  sys_timer_busy_wait_cycles(1000000);

  //printf("*** after wait, *foo = %d\r\n", *foo);

}


void test_task6(rid_t a, rid_t b, int *c) {
  printf("task6: a = %d, b = %d, c = 0x%p\r\n", a, b, c);
  sys_timer_busy_wait_cycles(10000000); // 1.0 sec
}

void test_task7(int *ch) {
  printf("%d: task7: ch = 0x%p, *ch = %d\r\n", sys_get_worker_id(), ch, *ch);
  (*ch)++;
  sys_timer_busy_wait_cycles(20000000); // 2.0 sec
}

void test_rfree(rid_t top, rid_t r) {
  printf("%d: Freeing rid %d\r\n", sys_get_worker_id(), r);
  sys_rfree(r);
  printf("%d: Done\r\n", sys_get_worker_id());
}

//void test_myrmics(int arg0, int arg1, int arg2, int arg3, int arg4, int arg5,
//                  int arg6, int arg7, int arg8, int arg9) {
void test_myrmics() {

  rid_t         top;
  rid_t         r1;
  rid_t         r2;
  rid_t         r3;
  rid_t         r4;
  rid_t         r5;
  rid_t         r6;
  rid_t         r7;
  rid_t         l6;
  rid_t         l7;
  int           *ch4;
  char          *ch6ra;
  char          *ch6rb;
  char          *ch6rc;
  char          *ch6la;
  char          *ch6lb;
  Node          *head;
  Node          *n;
  int           *foo;
  int           *bar;
  int           v1, v2, v3, *v4, v5, v6, v7, v8, v9, *v10;
  int           i;

  //printf("Hallo Myrmics world!\r\n");
  //printf("args: %d %d %d %d %d %d %d %d %d %d\r\n", 
  //       arg0, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9);




  // Create two objects
  //foo = sys_alloc(sizeof(int), 0);
  //*foo = 42;

  //bar = sys_alloc(sizeof(int), 0);
  //*bar = 43;
  //printf("foo = 0x%p, *foo = %d\r\n", foo, *foo);

  // // Spawn a task
  // #pragma myrmics task inout(foo)
  // test_task3(foo, 30);



  // // Create a tree of regions
  // r1 = sys_ralloc(0, 3);     // L3 scheduler
  // printf("r1 = %d\r\n", r1); 

  // r2 = sys_ralloc(r1, 2);    // L2 scheduler
  // printf("r2 = %d\r\n", r2);
  // r3 = sys_ralloc(r2, 2);
  // printf("r3 = %d\r\n", r3);

  // r4 = sys_ralloc(r3, 1);    // L1 scheduler
  // printf("r4 = %d\r\n", r4);
  // ch4 = sys_alloc(sizeof(int), r4);
  // printf("ch4 = 0x%p\r\n", ch4);
  // r5 = sys_ralloc(r4, 1);
  // printf("r5 = %d\r\n", r5);

  // r6 = sys_ralloc(r5, 0);    // L0 scheduler (left)
  // printf("r6 = %d\r\n", r6);
  // ch6ra = sys_alloc(123, r6);
  // ch6rb = sys_alloc(123, r6);
  // ch6rc = sys_alloc(123, r6);

  // r7 = sys_ralloc(r6, 0);
  // printf("r7 = %d\r\n", r7);
  // l6 = sys_ralloc(r5, 0);    // L0 scheduler (right)
  // printf("l6 = %d\r\n", l6);
  // ch6la = sys_alloc(1, l6);
  // ch6lb = sys_alloc(42, l6);
  // l7 = sys_ralloc(l6, 0);
  // printf("l7 = %d\r\n", l7);



  // #pragma myrmics task inout (ch4) region inout(r6, l6)
  // test_task6(r6, l6, ch4); // Task 0x00020001
  // 
  // //sys_timer_busy_wait_cycles(1000000); // 0.1 sec

  // #pragma myrmics task region inout(r2)
  // test_task4(r2); // Task 0x000C0001
  // 
  // //sys_timer_busy_wait_cycles(1000000); // 0.1 sec
  //   
  // #pragma myrmics task region inout(r5)
  // test_task4(r5); // Task 0x00020002


  int *ch1a;
  int *ch1b;

  top = sys_ralloc(0, 999);

  r1 = sys_ralloc(top, 0);
  ch1a = sys_alloc(sizeof(int), r1);
  ch1b = sys_alloc(sizeof(int), r1);
  *ch1a = 333;
  *ch1b = 666;

  #pragma myrmics task inout(ch1a)
  test_task7(ch1a);

  #pragma myrmics task inout(ch1b)
  test_task7(ch1b);

  #pragma myrmics task inout(ch1a)
  test_task7(ch1a);

  #pragma myrmics task inout(ch1b)
  test_task7(ch1b);

  #pragma myrmics task region inout(r1) in(ch1a) safe(ch1a)
  test_task4(r1, ch1a);

  #pragma myrmics task region inout(top) in(r1) safe(r1)
  test_rfree(top, r1);



  //*ch4 = 0;

  //for (i = 0; i < 20; i++) {
  //  #pragma myrmics task inout (ch4)
  //  test_task7(ch4);
  //}



  //*ch4 = 0;

  //for (i = 0; i < 10; i++) {
  //  #pragma myrmics task inout (ch4) in (i)
  //  test_task3(ch4, i);
  //}



  //sys_timer_busy_wait_cycles(30000000); // 3 sec



  //for (i = 0; i < 48; i++) {

  //  // Create a region
  //  //r = sys_ralloc(0, 0);
  //  //printf("%d: r = %d\r\n", i, r);

  //  #pragma myrmics task inout(foo)
  //  test_task3(foo, i);

  //  sys_timer_busy_wait_cycles(100000);
  //}



  //// Create 2nd region
  //r2 = sys_ralloc(0, 0);
  //printf("r2 = %d\r\n", r2);

  //// Spawn a 2nd task
  //#pragma myrmics task region inout(r2)
  //test_task4(r2);
  //

  //sys_timer_busy_wait_cycles(1000000);



  // // Create a region
  // r = sys_ralloc(0, 0);

  // // Create a head node
  // head = sys_alloc(sizeof(Node), r);
  // head->val = 42;
  // head->next = NULL;

  // // Call a func to add some more nodes
  // #pragma myrmics task region inout(r) inout(head)
  // test_task1(r, head);

  // v1 = 1;
  // v2 = 2;
  // v3 = 3;
  // v4 = sys_alloc(sizeof(int), r);
  // v5 = 5;
  // v6 = 6;
  // v7 = 7;
  // v8 = 8;
  // v9 = 9;
  // v10 = sys_alloc(sizeof(int), r);
  // #pragma myrmics task in(v1, v2, v3, v5, v6, v7, v8, v9) out(v4, v10)
  // test_task2(v1, v2, v3, v4, v5, v6, v7, v8, v9, v10);

  // // Wait for it
  // #pragma myrmics wait on(r)
  // for (n = head; n; n = n->next) {
  //   //printf("Node %p: val = %d\r\n", n, n->val);
  // }
}
#endif
