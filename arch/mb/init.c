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
// Abstract      : Arch-specific initialization, after basic boot has finished
//
// =============================[ CVS Variables ]=============================
//
// File name     : $RCSfile: init.c,v $
// CVS revision  : $Revision: 1.11 $
// Last modified : $Date: 2013/04/09 14:49:24 $
// Last author   : $Author: zakkak $
//
// ===========================================================================

#include <arch.h>
#include <mbs_regs.h>
#include <memory_management.h>
#include <fmpi.h>


// ===========================================================================
// ar_check_formic_links()      Reads all Formic board link status and
//                              compares to the expected value of how many
//                              GTP links should be up.
// ===========================================================================
// * INPUTS
//   int my_bid                 Local board ID
//   int my_cid                 Local core ID
//   int reset_link_error       If 1, clear link error bits when accessing
//                              the Formic board
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int ar_check_formic_links(int my_bid, int my_cid, int reset_link_error) {
  int bid;
  int link_status;
  int exp_links;
  int x;
  int y;
  int z;
  int i;
  int ret;
  int done_boards;
  int done_links;
  unsigned int word;

  ret = 0;
  done_boards = 0;
  done_links = 0;

  for (z = AR_FORMIC_MAX_Z; z >= AR_FORMIC_MIN_Z; z--) {
    for (y = AR_FORMIC_MAX_Y; y >= AR_FORMIC_MIN_Y; y--) {
      for (x = AR_FORMIC_MAX_X; x >= AR_FORMIC_MIN_X; x--) {

        bid = (x << 4) | (y << 2) | (z << 0);

        // Print where we're trying to connect
        kt_printf("\r%02d/%02d: Checking board 0x%02X [X=%d, Y=%d, Z=%d]...",
                  ++done_boards,
                  (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1) *
                  (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
                  (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1),
                  bid, x, y, z);

        // Read link status
        link_status = ar_read_link_status(my_bid, my_cid, bid);

        // Are the expected links up?
        exp_links = 0;
        if (x > AR_FORMIC_MIN_X)      exp_links |= 1 << 4; // x-
        if (x < AR_FORMIC_MAX_X)      exp_links |= 1 << 7; // x+
        if (y > AR_FORMIC_MIN_Y)      exp_links |= 1 << 6; // y-
        if (y < AR_FORMIC_MAX_Y)      exp_links |= 1 << 0; // y+
        if (z > AR_FORMIC_MIN_Z)      exp_links |= 1 << 3; // z-
        if (z < AR_FORMIC_MAX_Z)      exp_links |= 1 << 2; // z+

        if ((AR_ARM0_BID > -1) &&
            (x >= AR_FORMIC_ARM0_MIN_X) &&
            (x <= AR_FORMIC_ARM0_MAX_X) &&
            (y == AR_FORMIC_ARM0_Y) &&
            (z == AR_FORMIC_ARM0_Z))  exp_links |= 1 << 1; // w+ (ARM0 board)

        if ((AR_ARM1_BID > -1) &&
            (x >= AR_FORMIC_ARM1_MIN_X) &&
            (x <= AR_FORMIC_ARM1_MAX_X) &&
            (y == AR_FORMIC_ARM1_Y) &&
            (z == AR_FORMIC_ARM1_Z))  exp_links |= 1 << 1; // w+ (ARM1 board)

        if ((AR_XUP_BID > -1) &&
            (x == AR_FORMIC_XUP_X) &&
            (y == AR_FORMIC_XUP_Y) &&
            (z == AR_FORMIC_XUP_Z))   exp_links |= 1 << 1; // w+ (XUP board)

        for (i = 0; i < 8; i++) {
          if (exp_links & (1 << i)) {
            done_links++;
          }
        }

        // Print status
        if (exp_links != (link_status & 0xFF)) {
          kt_printf(" [ [1;31mFAIL[0m ]\r\n");
          ret = 1;
          for (i = 0; i < 8; i++) {
            if ((exp_links & (1 << i)) && !(link_status & (1 << i))) {
              kt_printf("                                              "
                        "> Check connector S%d (%s)\r\n", i,
                        (i == 0) ? "Y+" :
                        (i == 1) ? "W+" :
                        (i == 2) ? "Z+" :
                        (i == 3) ? "Z-" :
                        (i == 4) ? "X-" :
                        (i == 6) ? "Y-" :
                        (i == 7) ? "X+" : "?");
            }
          }
        }

        // Set the board LEDs to blink on errors on the expected links
        if (reset_link_error) {
          ar_write_link_status(my_cid, bid, 0);
          word = 0;
          for (i = 0; i < 8; i++) {
            if (exp_links & (1 << i)) {
              word |= 0x3 << (8 + 2*i);
            }
          }
          word |= (0x01 << 30) | (0x3 << 4);
          ar_write_brd_control(my_cid, bid, word);
        }
      }
    }
  }
  if (!ret) {
    kt_printf(" checked %d endpoints [ [1;32mPASS[0m ]\r\n", done_links);
  }
  else {
    kt_printf("\r\n");
  }

  return ret;
}


// ===========================================================================
// ar_check_formic_version()    Reads all Formic board version status and
//                              compares to the expected value.
// ===========================================================================
// * INPUTS
//   int my_bid                 Local board ID
//   int my_cid                 Local core ID
//
// * RETURN VALUE
//   int                        0 for success
// ===========================================================================
int ar_check_formic_version(int my_bid, int my_cid) {
  int bid;
  int x;
  int y;
  int z;
  int ret;
  int done_boards;
  unsigned int status;

  ret = 0;
  done_boards = 0;

  for (z = AR_FORMIC_MAX_Z; z >= AR_FORMIC_MIN_Z; z--) {
    for (y = AR_FORMIC_MAX_Y; y >= AR_FORMIC_MIN_Y; y--) {
      for (x = AR_FORMIC_MAX_X; x >= AR_FORMIC_MIN_X; x--) {

        bid = (x << 4) | (y << 2) | (z << 0);

        // Print where we're trying to connect
        kt_printf("\r%02d/%02d: Checking board 0x%02X [X=%d, Y=%d, Z=%d]...",
                  ++done_boards,
                  (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1) *
                  (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
                  (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1),
                  bid, x, y, z);

        // Read board status
        status = ar_read_status(my_bid, my_cid, bid);

        // Is it of the expected type and version?
        if ((AR_FORMIC_BOARD_TYPE    != (status >> 28)) ||
            (AR_FORMIC_BOARD_VERSION != ((status >> 12) & 0xFFFF))) {
          kt_printf(" [ [1;31mFAIL[0m ]\r\n");
          ret = 1;
          kt_printf("                                              "
                    "> Type %s, version %d\r\n",
                    (status >> 28 == 0) ? "formic_m1" :
                    (status >> 28 == 1) ? "formic_m4g8" :
                    (status >> 28 == 2) ? "formic_m8" :
                    (status >> 28 == 3) ? "formic_m8g8" :
                    (status >> 28 == 8) ? "versatile_axigw" :
                                          "UNKNOWN",
                    (status >> 12) & 0xFFFF);
        }
      }
    }
  }
  if (!ret) {
    kt_printf(" checked %d boards [ [1;32mPASS[0m ]\r\n", done_boards);
  }
  else {
    kt_printf("\r\n");
  }

  return ret;
}


// ===========================================================================
// ar_wake_up_formic_boards()   Checks link status, checks board version,
//                              transfers the included MicroBlaze ELF code
//                              to all Formic boards and wakes all their
//                              cores up
// ===========================================================================
// * INPUTS
//   int my_bid                 Local board ID
//   int my_cid                 Local core ID
// ===========================================================================
void ar_wake_up_formic_boards(int my_bid, int my_cid) {

  extern int    __linker_start_code;
  extern int    __linker_end_bss;
  int           mb_code_size;
  int           bid;
  int           x;
  int           y;
  int           z;
  int           i;


  // Check all Formic cube links
  kt_printf("\nChecking Formic cube links:\r\n");
  ar_check_formic_links(my_bid, my_cid, 1);


  // Check all Formic cube board versions
  kt_printf("\nChecking Formic cube board versions:\r\n");
  ar_check_formic_version(my_bid, my_cid);


  // Compute MicroBlaze code size to be transferred
  mb_code_size = (int) &__linker_end_bss - (int) &__linker_start_code;

  // Align it to 64-B cache-line boundaries, because the MicroBlaze linker
  // script does not do that
  if (mb_code_size & 63) {
    mb_code_size = (mb_code_size & 0xFFFFFFC0) + 64;
  }

  kt_printf("\r\nMicroblaze code from %p - %p [%d bytes]\r\n",
            &__linker_start_code, &__linker_end_bss, mb_code_size);


  // Wake up all Formic boards
  for (x = AR_FORMIC_MIN_X; x <= AR_FORMIC_MAX_X; x++) {
    for (y = AR_FORMIC_MIN_Y; y <= AR_FORMIC_MAX_Y; y++) {
      for (z = AR_FORMIC_MIN_Z; z <= AR_FORMIC_MAX_Z; z++) {

        // Find board ID
        bid = (x << 4) | (y << 2) | z;
        kt_printf("\rWaking up board x%dy%dz%d (0x%02X)...", x, y, z, bid);

        // If it's our board, don't transfer code
        if (bid != my_bid) {

          // Suspend core #0 (the only active by default)
          ar_suspend_core(my_bid, my_cid, bid, 0);

          // Transfer code
          ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP1, -mb_code_size);
          ar_dma_with_ack(my_cid,
                          my_bid, my_cid, (int) &__linker_start_code,
                          bid,    0xC,    0,
                          my_bid, my_cid, NOC_COUNTER_WAKEUP1,
                          mb_code_size, 0, 0, 0);

          while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP1)) {
            ;
          }

          // Wake up core #0
          ar_wake_up_core(my_bid, my_cid, bid, 0);
        }

        // Wake up the other board cores
        for (i = 1; i < AR_FORMIC_CORES_PER_BOARD; i++) {
          ar_wake_up_core(my_bid, my_cid, bid, i);
        }
      }
    }
  }

  kt_printf("\rSuccessfully woke up %d boards     \r\n",
      (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1) *
      (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
      (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1));
}


// ===========================================================================
// ar_init()                    Basic initialization functions. Initialize
//                              caches and install ART regions
// ===========================================================================
// * INPUTS
//   int my_bid                 Local board ID
//   int my_cid                 Local core ID
// ===========================================================================
void ar_init(int my_bid, int my_cid) {

  extern int    __linker_end_bss;
  int           code_end;


  // Make sure code page is enough for the compiled code
  code_end = (int) &__linker_end_bss;
  ar_assert(MM_MB_VA_CODE_PAGE == 0);
  ar_assert(code_end < MM_MB_VA_STACK_PAGE);

  // Clear & enable I/D caches
#ifdef MYRMICS
  ar_cache_enable(1);
#else
  ar_cache_enable(FMPI_CACHE_EPOCH_MODE);
#endif

  // Re-install ART map for write/execute anywhere memory accesses
  //*MBS_ART_ENTRY1 = (0x95 << 24) |  // C=1, I=0, R=0, X=1, U=0, P=1, V=1
  //                  (0xFFE << 12) | // Bound = 0xFFE
  //                  (0x000);        // Base = 0x000

  // Install dummy last all-inclusive region, so we don't hang execution
  // with an incomplete region setup
  ar_art_install_region(4,
                        0x0,                       // base
                        128 * 1024 * 1024,         // size (all DRAM)
                        1,                         // cacheable
                        0,                         // read only
                        1,                         // executable
                        0,                         // user accessible
                        1);                        // privileged accessible

  // Install ART code region
  ar_art_install_region(1,                         // region ID
                        MM_MB_VA_CODE_PAGE,        // base
                        0x100000,                  // size (1 MB)
                        1,                         // cacheable
                        1,                         // read only
                        1,                         // executable
                        1,                         // user accessible
                        1);                        // privileged accessible

  // Install ART stack/print buf region
  ar_art_install_region(2,                         // region ID
                        MM_MB_VA_STACK_PAGE,       // base
                        0x200000,                  // size (2 MB)
                        1,                         // cacheable
                        0,                         // read only
                        0,                         // executable
                        1,                         // user accessible
                        1);                        // privileged accessible

#ifdef MYRMICS
  // Install ART user region (only for Myrmics)
  ar_art_install_region(3,                         // region ID
                        MM_MB_VA_USER_BASE,        // base
                        MM_MB_USER_SIZE,           // size
                        1,                         // cacheable
                        0,                         // read only
                        0,                         // executable
                        1,                         // user accessible
                        1);                        // privileged accessible
#endif

  // Install ART kernel region (reinstall over dummy region)
  ar_art_install_region(4,                         // region ID
                        mm_va_kernel_base(my_cid), // base
                        MM_MB_KERNEL_SIZE,         // size
                        1,                         // cacheable
                        0,                         // read only
                        0,                         // executable
                        0,                         // user accessible
                        1);                        // privileged accessible

}
