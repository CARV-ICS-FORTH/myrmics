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
// Abstract      : Bare-metal communication measurements
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: test_bare.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/10/16 06:15:17 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <debug.h>

#define TEST_PASS "[32;1mPASS[0m"
#define TEST_FAIL "[31;1mFAIL[0m"

#define TEST_MASTER_BID 0x00
#define TEST_MASTER_CID 0

#define CONTROL_COUNTER 123
#define DONE_COUNTER    124

#define DMA_XPUT_TOTAL_SIZE     (32 * 1024 * 1024)

#define LOAD_LATENCY_MESSAGES       10000
#define LOAD_LATENCY_BUF_SIZE       (1 * 1024)
#define LOAD_LATENCY_MAX_INTERVAL   1000
#define LOAD_LATENCY_INTERVAL_STEP  25


// ===========================================================================
// Cache-to-cache and cache-to-mem xput measurement
// ===========================================================================
int test_bare_cache_dma_xput(int my_bid, int my_cid, int dst_bid, 
                             int dst_cid) {

  unsigned int  src_buf;
  unsigned int  dst_buf;
  int           size;
  int           total_size;
  int           i;
  unsigned int  cycles;
  unsigned int  xput;
  int           warmup;
  int           checksum;


  // Assign buffers
  src_buf = mm_va_kernel_base(my_cid);
  dst_buf = mm_va_kernel_base(dst_cid);

  // Initialize source buffer
  kt_printf("(initializing source buffer)\r\n");
  for (i = 0; i < 1024 * 1024; i += 4) {
    *((unsigned int *) (src_buf + i)) = i;
  }

  checksum = 0;
  //for (warmup = 1; warmup >= 0; warmup--) {
  warmup = 0;
    for (size = 64; size <= 1024 * 1024; size *= 2) {

      // Bring set to caches
      for (i = 0; i < size; i += 4) {
        checksum += *((unsigned int *) (src_buf + i));
      }

      ar_timer_reset();

      for (i = 0, total_size = 0; 
           total_size < DMA_XPUT_TOTAL_SIZE; 
           total_size += size, i++) {

        // DMA engine full?
        while ((ar_ni_status_get(my_cid) & 0xFF) < 2) {
          ;
        }
        ar_dma_no_ack(my_cid, my_bid,  my_cid,  src_buf, 
                              dst_bid, dst_cid, dst_buf,
                              size, 0, 0, 0);

      }

      // Wait until DMA engine empties
      while ((ar_ni_status_get(my_cid) & 0xFF) != 32) {
        ;
      }

      cycles = ar_timer_get_cycles();

      if (!warmup) {
        ar_uint_divide(total_size * 8, cycles, &xput, NULL);

        kt_printf("0x%02X/%X->0x%02X/%X: "
                  "%6d DMAs x %7d B = %8d B, %8u cc, %4u bpcc\r\n",
                  my_bid, my_cid, dst_bid, dst_cid, 
                  i, size, total_size, cycles, xput);
      }
      else if (size == 64) {
        kt_printf("(running once for warm-up)\r\n");
      }
    }
  //}
  return checksum;
}


// ===========================================================================
// Mem-to-mem xput measurement
// ===========================================================================
void test_bare_mem_dma_xput(int my_bid, int my_cid, int dst_bid) {

  unsigned int  src_buf;
  unsigned int  dst_buf;
  int           size;
  int           total_size;
  int           i;
  unsigned int  cycles;
  unsigned int  xput;


  // Assign buffers
  src_buf = mm_va_kernel_base(my_cid);
  dst_buf = src_buf + 1024 * 1024;

  // No initialization of buffers needed since it's a DRAM source (and not
  // a cache). No warmup for the same reason.

  for (size = 64; size <= 1024 * 1024; size *= 2) {

    ar_timer_reset();

    for (i = 0, total_size = 0; 
         total_size < DMA_XPUT_TOTAL_SIZE; 
         total_size += size, i++) {

      // DMA engine full?
      while ((ar_ni_status_get(my_cid) & 0xFF) < 2) {
        ;
      }
      ar_dma_no_ack(my_cid, my_bid,  0xC, src_buf, 
                            dst_bid, 0xC, dst_buf,
                            size, 0, 0, 0);

    }

    // Wait until DMA engine empties
    while ((ar_ni_status_get(my_cid) & 0xFF) != 32) {
      ;
    }

    cycles = ar_timer_get_cycles();

      ar_uint_divide(total_size * 8, cycles, &xput, NULL);

      kt_printf("0x%02X/C->0x%02X/C: "
                "%6d DMAs x %7d B = %8d B, %8u cc, %4u bpcc\r\n",
                my_bid, dst_bid, i, size, total_size, cycles, xput);
  }
}

// ===========================================================================
// DMA load-latency, 1 board, 8 cores
// ===========================================================================
void test_bare_dma_load_latency_1brd(int my_bid, int my_cid) {

  int           interval;
  unsigned int  seed;
  int           rnd;
  unsigned int  my_buf;
  int           time;
  int           i;
  int           j;
  

  // Initialize buffer
  my_buf = mm_va_kernel_base(my_cid);
  for (i = 0; i < 1024 * 1024 / 4; i++) {
    ((int *) my_buf)[i] = 0;
  }

  seed = 42 * (my_cid) + 1;

  for (interval = LOAD_LATENCY_MAX_INTERVAL; 
       interval >= 0; 
       interval -= LOAD_LATENCY_INTERVAL_STEP) {

    // Synchronize
    if (!my_cid) {
      for (j = 1; j < 8; j++) {
        ar_assert(ar_mbox_get(my_cid) == interval + 666);
      }
      for (j = 1; j < 8; j++) {
        ar_mbox_send(my_cid, my_bid, j, interval);
      }
    }
    else {
      ar_mbox_send(my_cid, my_bid, 0, interval + 666);
      ar_assert(ar_mbox_get(my_cid) == interval);
    }

    if (my_cid < 7) {
      time = 0;
      for (i = 0; i < LOAD_LATENCY_MESSAGES; i++) {

        ar_cnt_set(my_cid, 0, -LOAD_LATENCY_BUF_SIZE);

        // Write to core #7 buffer
        ar_timer_reset();
        ar_dma_with_ack(my_cid,
                        my_bid, my_cid,   my_buf,
                        my_bid, 7, mm_va_kernel_base(7),
                        my_bid, my_cid,   0,
                        LOAD_LATENCY_BUF_SIZE, 0, 0, 0);

        while (ar_cnt_get(my_cid, 0)) {
          ;
        }
        time += ar_timer_get_cycles();

        // Wait interval clock cycles
        rnd = ((seed = kt_rand(seed)) % 128) - 64;
        rnd *= (interval / 1024);
        rnd = interval + rnd;
        if (rnd < 0) {
          rnd = 0;
        }
        ar_timer_busy_wait_cycles(rnd);
      }

      // Print latency
      if (my_cid == 0) {
        kt_printf(
         "0x%02X/%d: interval = %7d cc, latency = %12d cc, %.3f Bpcc\r\n", 
                  my_bid, my_cid, interval, time, 
                  (float) (LOAD_LATENCY_MESSAGES * LOAD_LATENCY_BUF_SIZE) / (float) time);
      }
    }
  }
}


// ===========================================================================
// Ping-pong latency between two cores
// ===========================================================================
void test_bare_ping_pong(int src_bid, int src_cid, 
                         int dst_bid, int dst_cid) {

  int           my_cid;
  int           my_bid;
  int           peer_cid = -1;
  int           peer_bid = -1;
  unsigned int  time;
  int           i;


  // Who are we?
  my_cid = ar_get_core_id();
  my_bid = ar_get_board_id();

  // Who is our peer?
  if ((my_bid == src_bid) && (my_cid == src_cid)) {
    peer_bid = dst_bid;
    peer_cid = dst_cid;
  }
  else if ((my_bid == dst_bid) && (my_cid == dst_cid)) {
    peer_bid = src_bid;
    peer_cid = src_cid;
  }


  if ((my_bid == src_bid) && (my_cid == src_cid)) {
    ar_mbox_send(my_cid, peer_bid, peer_cid, 999);
  }

  for (i = 0; i < 10; i++) {
    ar_timer_reset();
    ar_mbox_get(my_cid);
    ar_mbox_send(my_cid, peer_bid, peer_cid, i);
    time = ar_timer_get_cycles();
  }

  if ((my_bid == dst_bid) && (my_cid == dst_cid)) {
    ar_mbox_get(my_cid);
  }
  
  if ((my_bid == src_bid) && (my_cid == src_cid)) {
    kt_printf("0x%02X/%d -> 0x%02X/%d: ping-pong latency = %d cc\r\n", 
        src_bid, src_cid, dst_bid, dst_cid, time);
  }
}


// ===========================================================================
// Board controller read latency
// ===========================================================================
void test_bare_brd_ctl_read(int my_bid, int my_cid) {

  int i;
  unsigned int time;
  
  for (i = 0; i < 10; i++) {
    ar_timer_reset();
#ifdef ARCH_MB
    *MBS_MNI_SRC_ADR  = (int) FORMIC_BRD_LINK_STATUS;
    *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
    *MBS_MNI_DMA_SIZE = 4;
    *MBS_MNI_BRD_NODE = (my_bid << 20) | (0xF << 16);
    *MBS_MNI_OPCODE   = (my_bid << 20) | (my_cid << 16) | 0x2;
    *MBS_MSL_ACCESS;
#else
  kt_printf("ERROR: not supported in ARM mode\r\n");
  return;
#endif
    time = ar_timer_get_cycles();
  }

  kt_printf("0x%02X/%d: board controller read latency = %d cc\r\n", 
      my_bid, my_cid, time);
}


// ===========================================================================
// Centralized mailbox-based barrier
// ===========================================================================
void test_bare_mbox_barrier(int my_bid, int my_cid) {

  unsigned int  time[10];
  int           i;
  int           j;
  

  // Master
  if (!my_bid && !my_cid) {

    for (i = 0; i < 10; i++) {

      ar_timer_reset();
      
      // Collect
      for (j = 0; j < 511; j++) {
        ar_mbox_get(my_cid);
      }

      // Broadcast
      for (j = 1; j < 512; j++) {
        ar_mbox_send(my_cid, j / 8, j % 8, 0);
      }

      time[i] = ar_timer_get_cycles();
    }

    kt_printf("Mailbox barrier = %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d\r\n"
              "                  %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }

  // Slaves
  else {
    for (i = 0; i < 10; i++) {
      ar_mbox_send(my_cid, 0, 0, 0);
      ar_mbox_get(my_cid);
    }
  }

}


// ===========================================================================
// Centralized counter-based barrier
// ===========================================================================
void test_bare_ccounter_barrier(int my_bid, int my_cid) {

  unsigned int  time[10];
  int           cnt;
  int           i;
  int           j;
  

  // Master
  if (!my_bid && !my_cid) {

    ar_cnt_set(my_cid, 0, -511);
    ar_cnt_set(my_cid, 1, -511);

    // Do a mailbox-based barrier to make sure counters are setup
    for (j = 0; j < 511; j++) {
      ar_mbox_get(my_cid);
    }
    for (j = 1; j < 512; j++) {
      ar_mbox_send(my_cid, j / 8, j % 8, 0);
    }

    

    for (i = 0; i < 10; i++) {

      ar_timer_reset();

      cnt = i % 2;

      while (ar_cnt_get(my_cid, cnt)) {
        ;
      }

      for (j = 1; j < 512; j++) {
        ar_cnt_incr(my_cid, j / 8, j % 8, cnt, 1);
      }

      ar_cnt_set(my_cid, cnt, -511);

      time[i] = ar_timer_get_cycles();
    }

    kt_printf("Centralized counter barrier = %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d\r\n"
              "                              %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }

  // Slaves
  else {

    ar_cnt_set(my_cid, 0, -1);
    ar_cnt_set(my_cid, 1, -1);

    // Do a mailbox-based barrier to make sure counters are setup
    ar_mbox_send(my_cid, 0, 0, 0);
    ar_mbox_get(my_cid);


    for (i = 0; i < 10; i++) {
      cnt = i % 2;
      ar_cnt_incr(my_cid, 0, 0, cnt, 1);
      while (ar_cnt_get(my_cid, cnt)) {
        ;
      }
      ar_cnt_set(my_cid, cnt, -1);
    }
  }

}


// ===========================================================================
// Hierarchical counter-based barrier
// ===========================================================================

// tree_direction:   0 = forward tree (barrier enter), 
//                   1 = reverse tree (barrier exit)
// level:            0 = leaf,
//                   9 = log2(512) = root
// id:               rank we're searching for, from 0 to 511
// *ret_cnt:         returned counter number (valid for rank and its sibling)
// *ret_id:          returned rank of *ret_cnt
// *ret_sibling_id:  returned rank of sibling in the tree
void test_bare_find_counter(int tree_direction, int level, int id, 
                            int *ret_cnt, int *ret_id,
                            int *ret_sibling_id) {

  int exp_level;
  int cnt;
  int whose;
  int whose_sibling;
  int i;


  // Forward tree has no counters for leaves
  ar_assert((tree_direction > 0) || (level > 0));

  // Find whose rank's the counter is
  whose = 0;
  for (i = 0; i < 9; i++) {
    if (i >= level) {
      whose |= id & (1 << i);
    }
  }

  // Find its sibling id
  whose_sibling = whose ^ (1 << level);

  // Expand level so that 0 ... 9 = forward tree, 9 ... 18 = reverse tree
  // (9 = root for both)
  if (tree_direction) {
    exp_level = 9 + (9 - level);
  } 
  else {
    exp_level = level;
  }

  // Counter depends on the expanded level (no leaf for forward)
  cnt = exp_level - 1;


  // Return stuff
  if (ret_cnt) {
    *ret_cnt = cnt;
  }
  if (ret_id) {
    *ret_id = whose;
  }
  if (ret_sibling_id) {
    *ret_sibling_id = whose_sibling;
  }
}

void test_bare_hcounter_barrier(int my_bid, int my_cid) {

  unsigned int  time[10];
  int           i;
  int           j;
  int           cnt;
  int           id;
  int           my_id;
  int           tmp_cnt;
  int           tmp_id;
  int           tmp_sibling;
  int           counters_to_init[20];
  int           values_to_init[20];
  int           num_counters;
  int           counter_to_enter;
  int           id_to_enter;
  int           bid_to_enter;
  int           cid_to_enter;
  int           counter_to_exit;
  

  // Compute which counters comprise our part of the trees, set their
  // notifications and store their info so we can init them again.
  // 
  // NOTE: we know it's a maximum of 18 counters per core, so we maintain 2
  // sets of trees by adding 32 to all computed counter numbers.
  num_counters = 0;
  my_id = (my_bid << 3) | my_cid;

  // Forward tree: counters wait 2 arrivals and notify next level
  for (i = 1; i < 9; i++) {

    // Find counter for this tree level
    test_bare_find_counter(0, i, my_id, &cnt, &id, NULL);

    // Are we responsible for it?
    if (id != my_id) {
      continue;
    }

    // Find counter of next tree level
    test_bare_find_counter(0, i + 1, my_id, &tmp_cnt, &tmp_id, NULL);

    // Setup our counter
    ar_cnt_set_notify_cnt(my_cid, cnt, -2, 
                          tmp_id >> 3, tmp_id & 0x7, tmp_cnt);
    ar_cnt_set_notify_cnt(my_cid, cnt + 32, -2, 
                          tmp_id >> 3, tmp_id & 0x7, tmp_cnt + 32);

    // Remember init value
    counters_to_init[num_counters] = cnt;
    values_to_init[num_counters] = -2;
    num_counters++;
  }


  // Roots of both trees: counter waits 2 arrivals and notifies two
  // tree siblings
  test_bare_find_counter(0, 9, my_id, &cnt, &id, NULL);

  // Are we responsible for it?
  if (id == my_id) {

    // Find counters of next tree level siblings
    test_bare_find_counter(1, 8, my_id, &tmp_cnt, &tmp_id, &tmp_sibling);

    // Setup our counter
    ar_cnt_set_dbl_notify_cnt(my_cid, cnt, -2, 
                              tmp_id >> 3,      tmp_id & 0x7,      tmp_cnt,
                              tmp_sibling >> 3, tmp_sibling & 0x7, tmp_cnt);
    ar_cnt_set_dbl_notify_cnt(my_cid, cnt + 32, -2, 
                              tmp_id >> 3,      tmp_id & 0x7,      tmp_cnt + 32,
                              tmp_sibling >> 3, tmp_sibling & 0x7, tmp_cnt + 32
                             );

    // Remember init value
    counters_to_init[num_counters] = cnt;
    values_to_init[num_counters] = -2;
    num_counters++;
  }


  // Reverse tree: counter waits 1 arrival and notifies two tree siblings
  for (i = 8; i >= 1; i--) {

    // Find counter for this tree level
    test_bare_find_counter(1, i, my_id, &cnt, &id, NULL);

    // Are we responsible for it?
    if (id != my_id) {
      continue;
    }

    // Find counters of next tree level siblings
    test_bare_find_counter(1, i - 1, my_id, &tmp_cnt, &tmp_id, &tmp_sibling);

    // Setup our counter
    ar_cnt_set_dbl_notify_cnt(my_cid, cnt, -1, 
                              tmp_id >> 3,      tmp_id & 0x7,      tmp_cnt,
                              tmp_sibling >> 3, tmp_sibling & 0x7, tmp_cnt);
    ar_cnt_set_dbl_notify_cnt(my_cid, cnt + 32, -1, 
                              tmp_id >> 3,      tmp_id & 0x7,      tmp_cnt + 32,
                              tmp_sibling >> 3, tmp_sibling & 0x7, tmp_cnt + 32
                             );

    // Remember init value
    counters_to_init[num_counters] = cnt;
    values_to_init[num_counters] = -1;
    num_counters++;
  }


  // Finally, compute which counter we need to notify in order to enter the 
  // barrier...
  test_bare_find_counter(0, 1, my_id, &counter_to_enter, &id_to_enter, NULL);
  bid_to_enter = id_to_enter >> 3;
  cid_to_enter = id_to_enter & 0x7;

  // ... and which counter to spin onto in order to exit the barrier. We need
  // to initialize this one too.
  test_bare_find_counter(1, 0, my_id, &counter_to_exit, NULL, NULL);

  ar_cnt_set(my_cid, counter_to_exit, -1);
  ar_cnt_set(my_cid, counter_to_exit + 32, -1);

  counters_to_init[num_counters] = counter_to_exit;
  values_to_init[num_counters] = -1;
  num_counters++;


  // Sanity checks
  ar_assert(num_counters < 20);
  for (i = 0; i < num_counters; i++) {
    ar_assert(counters_to_init[i] < 32);
    ar_assert(values_to_init[i] < 0);
  }


  // Do a mailbox barrier once, so that we know everyone has finished setting
  // up its counters.
  if (!my_bid && !my_cid) {
    for (j = 0; j < 511; j++) {
      ar_mbox_get(my_cid);
    }
    for (j = 1; j < 512; j++) {
      ar_mbox_send(my_cid, j / 8, j % 8, 0);
    }
  }
  else {
    ar_mbox_send(my_cid, 0, 0, 0);
    ar_mbox_get(my_cid);
  }



  // Now enter the barriers in rapid succession
  for (i = 0; i < 10; i++) {

    ar_timer_reset();

    // Even or odd phase?
    if (i % 2 == 0) {

      // Even-phased barrier
      ar_cnt_incr(my_cid, bid_to_enter, cid_to_enter, counter_to_enter, 1);
      while (ar_cnt_get(my_cid, counter_to_exit)) {
        ;
      }

      // Reset even-phased barrier counters
      for (j = 0; j < num_counters; j++) {
#ifdef ARCH_ARM
        *ARS_CNT_VAL(my_cid, counters_to_init[j]) = values_to_init[j];
#endif
#ifdef ARCH_MB
        *MBS_CNT_VAL(counters_to_init[j]) = values_to_init[j];
#endif
      }
    }

    else {

      // Odd-phased barrier
      ar_cnt_incr(my_cid, bid_to_enter, cid_to_enter, counter_to_enter + 32, 1);
      while (ar_cnt_get(my_cid, counter_to_exit + 32)) {
        ;
      }

      // Reset odd-phased barrier counters
      for (j = 0; j < num_counters; j++) {
#ifdef ARCH_ARM
        *ARS_CNT_VAL(my_cid, counters_to_init[j] + 32) = values_to_init[j];
#endif
#ifdef ARCH_MB
        *MBS_CNT_VAL(counters_to_init[j] + 32) = values_to_init[j];
#endif
      }
    }

    time[i] = ar_timer_get_cycles();
  }

  if (!my_bid && !my_cid) {
    kt_printf("Hierarchical Counter barrier = %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d\r\n"
              "                               %12d cycles\r\n",
              time[0], time[1], time[2], time[3], time[4],
              time[5], time[6], time[7], time[8], time[9]);
  }


}



// ===========================================================================
// ===========================================================================
void test_bare_phase_enter(int my_bid, int my_cid, int phase, char *title) {
  int i;

  // Broadcast phase number
  if ((my_bid == TEST_MASTER_BID) && (my_cid == TEST_MASTER_CID)) {
    kt_printf("\nPhase %d: %s\r\n", phase, title);
    for (i = 0; i < 512; i++) {
      ar_cnt_incr(my_cid, i / 8, i % 8, CONTROL_COUNTER, 1);
    }
  }

  // Phase barrier
  while (ar_cnt_get(my_cid, CONTROL_COUNTER) < phase) {
    ;
  }
}


// ===========================================================================
// ===========================================================================
void test_bare_phase_exit(int my_bid, int my_cid) {

  // Test master: make sure everybody finished the phase
  if ((my_bid == TEST_MASTER_BID) && (my_cid == TEST_MASTER_CID)) {
    while (ar_cnt_get(my_cid, DONE_COUNTER) < 511) {
      ;
    }

    // Prepare for next phase
    ar_cnt_set(my_cid, DONE_COUNTER, 0);
  }

  // Test slaves: notify master
  else {
    ar_cnt_incr(my_cid, TEST_MASTER_BID, TEST_MASTER_CID, DONE_COUNTER, 1);
  }
}


// ===========================================================================
// test_bare()                  Various bare-metal communication test patterns
// ===========================================================================
void test_bare() {

  int           my_bid;
  int           my_cid;
  int           bid;
  int           cid;
  int           i;


  // Who are we?
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();





  // =========================================================================
  // =========================================================================
#ifdef ARCH_MB

  // Relevant test_bare() code in 0x12500-0x12800

#define BASE            0x00312680              // multiple of 0x40
#define NUM_ARRAYS      60  
#define NUM_ELEMENTS    512                     // multiple of 16
#define PAD_ELEMENTS    (8192 - NUM_ELEMENTS)   // multiple of 16

#define CMD_WORK        0xCAFE0000
#define CMD_FLUSH       0xCAFE1111
#define CMD_FLUSH_DONE  0xCAFE2222

#define IGN_D           1


  int   *arrays[NUM_ARRAYS];
  int   mark;
  int   nxt_cid;
  int   prv_cid;
  int   cmd;
  int   rep;
  int   pos;
  int   j;
  int   k;


  if (my_bid == 0x00) {

    // Install common "user" region (just as in Myrmics) replacing the
    // existing private kernel one
    ar_art_install_region(4,                         // region ID
                          0x00300000,                // base
                          (93 * 1024 * 1024),        // size
                          1,                         // cacheable
                          0,                         // read only
                          0,                         // executable
                          1,                         // user accessible
                          1);                        // privileged accessible


    // Enable performance registers
    *MBS_CPU_CONTROL |= (1 << 24);

    // Allocate arrays
    for (i = 0; i < NUM_ARRAYS; i++) {
      arrays[i] = (int *) (BASE + i * (NUM_ELEMENTS + PAD_ELEMENTS) * 4);
    }

    // Cores in a chain
    nxt_cid = (my_cid + 1) % 8;
    prv_cid = (my_cid + 7) % 8;

    // Core 0: initialize arrays and send to next one
    if (!my_cid) {

      // Init counter for all arrays
      *MBS_CNT_VAL(0) = -NUM_ARRAYS * NUM_ELEMENTS * 4;
      *MBS_CNT_NTF_BRD_NODE(0) = 0;

      // Init
      for (i = 0; i < NUM_ARRAYS; i++) {
        mark = 0xA0000000 | (i << 16);
        for (j = 0; j < NUM_ELEMENTS; j++) {
          arrays[i][j] = mark;
        }

        // Send to 1
        ar_dma_with_ack(my_cid,
                        my_bid, my_cid, (int) arrays[i],
                        my_bid, nxt_cid, (int) arrays[i],
                        my_bid, my_cid, 0, 
                        NUM_ELEMENTS * 4, IGN_D, 0, 0);
      }

      // Wait until DMAs are done
      while (*MBS_CNT_VAL(0)) {
        ;
      }

      // Notify next mailbox to work on the buffers
      ar_mbox_send(my_cid, my_bid, nxt_cid, CMD_WORK);
    }


    // =======================================================================

    rep = (my_cid) ? 0 : 1;
    while (1) {
      
      // Wait for a mailbox command
      cmd = ar_mbox_get(my_cid);

      // Flush?
      if (cmd == CMD_FLUSH) {

        // Flush L2
        *MBS_CHE_MAINT = 0x20000;
        //kt_printf("%d: Flushed L2\r\n", my_cid);

        // Flush/clear all caches
        //*MBS_CHE_MAINT = 0x30101;
        //kt_printf("%d: Cleared caches\r\n", my_cid);

        // Tell next core we did it
        ar_mbox_send(my_cid, my_bid, nxt_cid, CMD_FLUSH_DONE);

        // Wait for next command
        continue;
      }
      
      ar_assert(cmd == CMD_WORK);
      
      // Decide new position
      pos = rep * 8 + my_cid;
      if (pos == NUM_ELEMENTS) {
        kt_printf("%d: All done\r\n", my_cid);
        ar_uart_flush();
        ar_timer_busy_wait_msec(100);
        kt_printf("%s", DBG_KILL_CONNECTION);
        break;
      }

      // Init counter for all arrays
      *MBS_CNT_VAL(0) = -NUM_ARRAYS * NUM_ELEMENTS * 4;
      *MBS_CNT_NTF_BRD_NODE(0) = 0;



      //*MBS_PRF_DL2 = 0;

      for (k = 0; k < NUM_ARRAYS; k++) {
        
        // Pick an array in a different order per core
        i = (my_cid * 7 * rep + 12345 + k) % NUM_ARRAYS;

        mark = 0xA0000000 | (i << 16);

        // Check up to current pos
        for (j = 0; j < NUM_ELEMENTS; j++) {
          if (((j < pos) && (arrays[i][j] != (mark | j))) ||
              ((j >= pos) && (arrays[i][j] != mark))) {
            kt_printf("%d: Error: arrays[%d][%d] = 0x%08X\r\n", 
                my_cid, i, j, arrays[i][j]);
            while (1) {
              ;
            }
          }
        }

        // Update new pos
        arrays[i][pos] = mark | pos;

        // Send to the next one
        ar_dma_with_ack(my_cid,
                        my_bid, my_cid, (int) arrays[i],
                        my_bid, nxt_cid, (int) arrays[i],
                        my_bid, my_cid, 0, 
                        NUM_ELEMENTS * 4, IGN_D, 0, 0);
      }

      //kt_printf("%d: Work mode: Data L2 %d misses, %d hits\r\n", my_cid, 
      //          *MBS_PRF_DL2 >> 16, *MBS_PRF_DL2 & 0xFFFF);

      // Wait until DMAs are done
      while (*MBS_CNT_VAL(0)) {
        ;
      }

      //*MBS_PRF_DL2 = 0;

      //// Rewrite all arrays with dummy values
      //for (i = 0; i < NUM_ARRAYS; i++) {
      //  for (j = 0; j < NUM_ELEMENTS; j++) {
      //    arrays[i][j] = 0x11111111;
      //  }
      //}

      //kt_printf("%d: Zero mode: Data L2 %d misses, %d hits\r\n", my_cid, 
      //          *MBS_PRF_DL2 >> 16, *MBS_PRF_DL2 & 0xFFFF);
      
      // Done
      rep++;
      kt_printf("%d: done %d\r\n", my_cid, pos);

      // Notify previous mailbox to flush its cache and wait till it's done
      ar_mbox_send(my_cid, my_bid, prv_cid, CMD_FLUSH);
      ar_assert(ar_mbox_get(my_cid) == CMD_FLUSH_DONE);
      
      // Notify next mailbox
      ar_mbox_send(my_cid, my_bid, nxt_cid, CMD_WORK);
    }

  }

#endif

  // Finished
  while (1) {
    ;
  }

  // =========================================================================
  // =========================================================================








  // Test master?
  if ((my_bid == TEST_MASTER_BID) && (my_cid == TEST_MASTER_CID)) {

    // Synchronize everybody through mailbox, so we know their control counters
    // are all initialized to 0
    for (i = 0; i < 511; i++) {
      ar_assert(ar_mbox_get(my_cid) == 0x12345678);
    }
    ar_cnt_set(my_cid, CONTROL_COUNTER, 0);
    ar_cnt_set(my_cid, DONE_COUNTER, 0);
    for (i = 0; i < 512; i++) {
      bid = i / 8;
      cid = i % 8;
      if ((bid != my_bid) || (cid != my_cid)) {
        ar_mbox_send(my_cid, bid, cid, 0x87654321);
      }
    }
  }
  else {
    // All slaves initialize their control counter
    ar_cnt_set(my_cid, CONTROL_COUNTER, 0);

    // Initial synchronization through mailbox
    ar_mbox_send(my_cid, TEST_MASTER_BID, TEST_MASTER_CID, 0x12345678);
    ar_assert(ar_mbox_get(my_cid) == 0x87654321);
  }


  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 1, 
                        "Cache-to-cache DMA throughput on-board");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_cache_dma_xput(my_bid, my_cid, 0x00, 0x01);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 2, 
                        "Cache-to-cache DMA throughput neighbor board");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_cache_dma_xput(my_bid, my_cid, 0x10, 0x00);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 3, 
                        "Cache-to-mem DMA throughput on-board");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_cache_dma_xput(my_bid, my_cid, 0x00, 0xC);
  }

  test_bare_phase_exit(my_bid, my_cid);


  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 4, 
                        "Cache-to-mem DMA throughput neighbor board");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_cache_dma_xput(my_bid, my_cid, 0x10, 0xC);
  }

  test_bare_phase_exit(my_bid, my_cid);


  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 5, 
                        "Mem-to-mem DMA throughput on-board");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_mem_dma_xput(my_bid, my_cid, 0x00);
  }

  test_bare_phase_exit(my_bid, my_cid);


  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 6, 
                  "Mem-to-mem DMA throughput neighbor board");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_mem_dma_xput(my_bid, my_cid, 0x10);
  }

  test_bare_phase_exit(my_bid, my_cid);


  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 7, 
                  "DMA load-latency, 8 cores");
  
  if (my_bid == 0x00) {
    test_bare_dma_load_latency_1brd(my_bid, my_cid);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 8, 
                  "Ping-pong intra-board latency");
  
  if ((my_bid == 0x00) && (my_cid < 2)) {
    test_bare_ping_pong(0x00, 0, 0x00, 1);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 9, 
                  "Ping-pong nearest board (1 hop) latency");
  
  if (((my_bid == 0x00) || (my_bid == 0x01)) && (my_cid == 0)) {
    test_bare_ping_pong(0x00, 0, 0x01, 0);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 10, 
                  "Ping-pong furthest board (9 hops) latency");
  
  if (((my_bid == 0x00) || (my_bid == 0x3F)) && (my_cid == 0)) {
    test_bare_ping_pong(0x00, 0, 0x3F, 0);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 11, 
                  "Board controller read latency");
  
  if ((my_bid == 0x00) && (my_cid == 0)) {
    test_bare_brd_ctl_read(0x00, 0);
  }

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 12, 
                  "Centralized Mailbox-based barrier");
  
  test_bare_mbox_barrier(my_bid, my_cid);

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 13, 
                  "Hierarchical Counter-based barrier");
  
  test_bare_hcounter_barrier(my_bid, my_cid);

  test_bare_phase_exit(my_bid, my_cid);

  // =========================================================================
  test_bare_phase_enter(my_bid, my_cid, 14, 
                  "Centralized Counter-based barrier");
  
  test_bare_ccounter_barrier(my_bid, my_cid);

  test_bare_phase_exit(my_bid, my_cid);




  // That's the end
  if ((my_bid == TEST_MASTER_BID) && (my_cid == TEST_MASTER_CID)) {
    kt_printf("\r\nAll done.\r\n");
  }
  while (1) {
    ;
  }
}
