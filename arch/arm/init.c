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
// CVS revision  : $Revision: 1.8 $
// Last modified : $Date: 2012/11/01 17:20:26 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <arm_regs.h>

// Page table section attributes. Reference: ARMARM B3-8
#define SECTION_NON_GLOBAL      (1 << 17)
#define SECTION_SHAREABLE       (1 << 16)
//#define SECTION_SHAREABLE       (0)
#define SECTION_PRIVILEGED_ONLY (1 << 10)
#define SECTION_FULL_ACCESS     ((1 << 11) | (1 << 10))
#define SECTION_EXECUTE_NEVER   (1 << 4)
#define SECTION_CACHEABLE       (1 << 3)
#define SECTION_BUFFERABLE      (1 << 2)
#define SECTION_TEX_ALLOC       (1 << 12)
#define SECTION_IS_SECTION      (2)



// ===========================================================================
// ar_tlb_map_section()         Maps a 1-MB page into the page table
// ===========================================================================
// * INPUTS
//   unsigned int virtual_adr   Virtual address of the page
//   unsigned int physical_adr  Physical address of the page
//   unsigned int attr          Section attributes bitmap of the translation
// ===========================================================================
void ar_tlb_map_section(unsigned int virtual_adr, unsigned int physical_adr, 
                        unsigned int attr) {

  volatile unsigned int *entry;

  // Find page table entry address
  entry = (volatile unsigned int *) 
    ((MM_ARM_VA_PAGE_TABLE & 0xFFFFC000) | ((virtual_adr >> 18) & 0x3FFC));

  // Write entry. Reference: ARMARM B3.3.1
  *entry = (physical_adr & 0xFFF00000) | SECTION_IS_SECTION | attr;

}


// ===========================================================================
// ar_tlb_unmap_all_sections()  Sets all page table entries to invalid
// ===========================================================================
void ar_tlb_unmap_all_sections() {
  unsigned int i;

  for (i = 0; i < 16 * 1024; i += 4) {
    *((volatile unsigned int *) (MM_ARM_VA_PAGE_TABLE + i)) = 0;
  }
}


// ===========================================================================
// ar_tlb_install_page_table()  Translate the memory map into an initial 
//                              page table on the DRAM
// ===========================================================================
void ar_tlb_install_page_table() {
  
  int          cid;
  unsigned int i;


  // Prepare Page Table. All pages are invalid to begin with.
  ar_tlb_unmap_all_sections();

  // Code page
  ar_tlb_map_section(MM_ARM_VA_CODE_PAGE, MM_ARM_PA_CODE_PAGE, 
                    SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                    SECTION_CACHEABLE | SECTION_BUFFERABLE |
                    SECTION_TEX_ALLOC);

  // Stack page
  ar_tlb_map_section(MM_ARM_VA_STACK_PAGE, MM_ARM_PA_STACK_PAGE,
                    SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                    SECTION_CACHEABLE | SECTION_BUFFERABLE | 
                    SECTION_EXECUTE_NEVER | SECTION_TEX_ALLOC);

  // Page table
  ar_tlb_map_section(MM_ARM_VA_PAGE_TABLE, MM_ARM_PA_PAGE_TABLE,
                    SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                    SECTION_CACHEABLE | SECTION_BUFFERABLE | 
                    SECTION_EXECUTE_NEVER | SECTION_TEX_ALLOC);

  // Kernel space pages, private per core
  for (cid = 0; cid < 4; cid++) {
    for (i = 0; i < MM_KERNEL_SIZE; i += 1024*1024) {
      ar_tlb_map_section(mm_va_kernel_base(cid) + i, mm_pa_kernel_base(cid) + i,

      // This works (non-cacheable)
                  SECTION_SHAREABLE | SECTION_FULL_ACCESS |
                  SECTION_EXECUTE_NEVER | SECTION_TEX_ALLOC);

      // This does not work (fully cacheable)
      //        SECTION_SHAREABLE | SECTION_FULL_ACCESS |
      //        SECTION_CACHEABLE | SECTION_BUFFERABLE |
      //        SECTION_TEX_ALLOC);

      // This works (Inner non-cacheable, Outer cacheable)
      //        SECTION_SHAREABLE | SECTION_FULL_ACCESS |
      //        (1 << 14) |             // TEX[2] = 1
      //        (0 << 13) | (1 << 12) | // TEX[1:0] = 2'b01 (outer = wb, wa)
      //        (0 << 3)  | (0 << 2)    // C/B = 2'b00 (inner = uncacheable)
      //        );

      // This does not work (Outer non-cacheable, Inner cacheable)
      //        SECTION_SHAREABLE | SECTION_FULL_ACCESS |
      //        (1 << 14) |             // TEX[2] = 1
      //        (0 << 13) | (0 << 12) | // TEX[1:0] = 2'b00 (outer = uncach.)
      //        (0 << 3)  | (1 << 2)    // C/B = 2'b01 (outer = non cacheable)
      //        );
    }
  }

#ifdef MYRMICS
  // User space pages, common for all cores (Myrmics only)
  for (i = 0; i < MM_USER_SIZE; i += 1024*1024) {
    ar_tlb_map_section(MM_VA_USER_BASE + i, MM_PA_USER_BASE + i,
                       SECTION_SHAREABLE | SECTION_FULL_ACCESS |
                       SECTION_TEX_ALLOC);
  }
#endif

  // Device map 0 (motherboard and peripherals)
  ar_tlb_map_section(MM_ARM_VA_DEV0_PAGE, MM_ARM_PA_DEV0_PAGE,
                     SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                     SECTION_EXECUTE_NEVER);

  // Device map 1 (testchip private memory)
  ar_tlb_map_section(MM_ARM_VA_DEV1_PAGE, MM_ARM_PA_DEV1_PAGE,
                     SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                     SECTION_EXECUTE_NEVER);

  // Test-routines 1-1 maps
  //map_section(0x48000000, 0x48000000,   // SRAM
  //map_section(0x44000000, 0x44000000,   // NOR1
  //map_section(0xE0000000, 0xE0000000,   // FPGA Reg
  //map_section(0xF0000000, 0xF0000000,   // FPGA ZBT

  ar_tlb_map_section(MM_ARM_VA_ARS0_PAGE, MM_ARM_PA_ARS0_PAGE,   // ARS #0
                     SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                     SECTION_EXECUTE_NEVER);

  ar_tlb_map_section(MM_ARM_VA_ARS1_PAGE, MM_ARM_PA_ARS1_PAGE,   // ARS #1
                     SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                     SECTION_EXECUTE_NEVER);

  ar_tlb_map_section(MM_ARM_VA_ARS2_PAGE, MM_ARM_PA_ARS2_PAGE,   // ARS #2
                     SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                     SECTION_EXECUTE_NEVER);

  ar_tlb_map_section(MM_ARM_VA_ARS3_PAGE, MM_ARM_PA_ARS3_PAGE,   // ARS #3
                     SECTION_SHAREABLE | SECTION_PRIVILEGED_ONLY |
                     SECTION_EXECUTE_NEVER);
}

  
// ===========================================================================
// ar_init_fpu_mmu_caches()     Init FPU, caches, branch predictor and MMU
// ===========================================================================
void ar_init_fpu_mmu_caches() {

  volatile unsigned int u32;


  // Notify MMU where the Page Table is
  ARM_TTBR0_WRITE(MM_ARM_PA_PAGE_TABLE | 0xB);    // Full cacheable (?)
  //ARM_TTBR0_WRITE(MM_ARM_PA_PAGE_TABLE | 0x0);  // Non cacheable
  //ARM_TTBR1_WRITE(MM_ARM_PA_PAGE_TABLE | 0x0);  // Non cacheable

  // Enable domain 0
  ARM_DACR_WRITE(1); 


  // Enable SMP bit
  ARM_ACTLR_READ(u32);
  u32 |= (1 << 6);
  ARM_ACTLR_WRITE (u32);

  // Enable FPU
  ar_enable_fpu();

  // Now we turned on MMU and FPU, invalidate caches once again. 
  // Remember that D-cache was never turned on by boot sequence.
  ar_disable_icache();
  ar_invalidate_caches();
  asm (
      "dsb\n"
      "isb\n"
      );

  // Turn on MMU
  ARM_SCTLR_READ(u32);                    
  u32 |= (1 << 0);
  ARM_SCTLR_WRITE (u32);

  // Enable caches
  ar_enable_dcache();
  ar_enable_icache();
  asm (
      "dsb\n"
      "isb\n"
      );

  // Enable branch predictor
  ar_enable_branch_pred();
  asm (
      "dsb\n"
      "isb\n"
      );
}


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

  extern int    __mb_start_include;
  extern int    __mb_end_include;
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
  mb_code_size = (int) &__mb_end_include - (int) &__mb_start_include;
  kt_printf("\r\nMicroblaze code from %p - %p [%d bytes]\r\n",
            &__mb_start_include, &__mb_end_include, mb_code_size);


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
                          my_bid, my_cid, (int) &__mb_start_include,
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

  kt_printf("\rSuccessfully woke up %d boards    \r\n", 
      (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1) * 
      (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) * 
      (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1));
}


// ===========================================================================
// ar_init()                    Basic initialization functions. Initialize
//                              caches, FPU, install page table, etc.
// ===========================================================================
// * INPUTS
//   int my_bid                 Local board ID
//   int my_cid                 Local core ID 
// ===========================================================================
void ar_init(int my_bid, int my_cid) {
  
  extern int    __mb_end_include;
  int           code_end;
  int           my_other_bid;
  int           local_active_cores;
  int           all_active_cores;
  unsigned int  sp;
  int           i;

  // Board-specific inits
  my_other_bid = my_bid + 0x10;  // That's the X+1 link
  if (my_bid == AR_ARM0_BID) {
    local_active_cores = AR_ARM0_CORES_PER_BOARD;
  }
  else if (my_bid == AR_ARM1_BID) {
    local_active_cores = AR_ARM1_CORES_PER_BOARD;
  }
  else {
    ar_panic("Unexpected ARM board ID");
  }
  
  all_active_cores = AR_ARM0_CORES_PER_BOARD;
  if (AR_ARM1_BID > -1) {
    all_active_cores += AR_ARM1_CORES_PER_BOARD;
  }

  // Core 0
  if (!my_cid) {

    // Make sure code page is enough for the compiled code
    code_end   = (int) &__mb_end_include;
    ar_assert(MM_ARM_VA_CODE_PAGE == 0);
    ar_assert(code_end < MM_ARM_VA_STACK_PAGE);

    // Write page table to main memory
    ar_tlb_install_page_table();

    // Enable caches, FPU, MMU
    ar_init_fpu_mmu_caches();

    // Enable SCU
    *ARM_SCU_CONTROL = 1;

    // Initialize UART lock
    ar_lock_release((void *) MM_ARM_VA_PRINT_LOCK);
    
    // Greet
    ARM_SP_READ(sp);
    kt_printf("ARM Core ID = %d, board IDs = 0x%02X / 0x%02X, SP = 0x%08X\r\n", 
              my_cid, my_bid, my_other_bid, sp);

    // Wake up core 1, who still waits on the mailbox from boot.s
    if (local_active_cores > 1) {
      ar_mbox_send(my_cid, my_bid, 1, 0xCAFE0001);
    }
  }

  // Cores 1-3
  else {

    // Enable caches, FPU, MMU
    ar_init_fpu_mmu_caches();

    // Greet
    ARM_SP_READ(sp);
    kt_printf("ARM Core ID = %d, board IDs = 0x%02X / 0x%02X, SP = 0x%08X\r\n", 
              my_cid, my_bid, my_other_bid, sp);
    
    // Wake up next core, who still waits on the mailbox from boot.s
    if (my_cid < local_active_cores - 1) {
      ar_mbox_send(my_cid, my_bid, (my_cid + 1) % 4, 0xCAFE0000 | (my_cid + 1));
    }
    else if ((my_bid == AR_ARM0_BID) && (AR_ARM1_BID > -1)) {
      ar_mbox_send(my_cid, AR_ARM1_BID, 0, 0xCAFE0000);
    }
  }

  // Do a mailbox-based barrier until all cores are woken up, so that the
  // initial greeting printfs are all printed together
  if ((my_bid == AR_ARM0_BID) && (!my_cid)) {

    // Wait for all others to finish
    for (i = 0; i < all_active_cores - 1; i++) {
      ar_assert(ar_mbox_get(my_cid) == 0xDEADBABE);
    }

    // Send them the all-ok
    for (i = 1; i < AR_ARM0_CORES_PER_BOARD; i++) {
      ar_mbox_send(my_cid, AR_ARM0_BID, i, 0xBABE0000 + i);
    }
    if (AR_ARM1_BID > -1) {
      for (i = 0; i < AR_ARM1_CORES_PER_BOARD; i++) {
        ar_mbox_send(my_cid, AR_ARM1_BID, i, 0xBABE0000 + i);
      }
    }
  }
  else {
    ar_mbox_send(my_cid, AR_ARM0_BID, 0, 0xDEADBABE);
    ar_assert(ar_mbox_get(my_cid) == 0xBABE0000 + my_cid);
  }
}
