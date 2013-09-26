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
// Abstract      : Main C language entry point for the MPI version
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: mpi_main.c,v $
// CVS revision  : $Revision: 1.26 $
// Last modified : $Date: 2013/03/15 15:54:52 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <memory_management.h>
#include <fmpi.h>
#include <bench.h>
#include <video.h>
#include <debug.h>



// ==========================================================================
// Main: entry point after boot.s
// ==========================================================================
int main() {

  Context       *context;
  int           total_cores;
  int           my_bid;
  int           my_cid;
  unsigned int  word;
  int           ret;
  int           x;
  int           y;
  int           z;
  int           bid;
  int           i;


  // Get core & board ID
  my_bid = ar_get_board_id();
  my_cid = ar_get_core_id();

  // Arch-specific initialization (FPU/caches/MMU, etc)
  ar_init(my_bid, my_cid);


  // Boot master
  if ((my_bid == AR_BOOT_MASTER_BID) && (!my_cid)) {

    // Greet
    kt_printf("Boot master is 0x%02X/%d\r\n", my_bid, my_cid);

    // Check links, board versions, transfer code and wake up Formic boards
    ar_wake_up_formic_boards(my_bid, my_cid);

    // Compute how many cores should be activated
    total_cores = (AR_FORMIC_MAX_X - AR_FORMIC_MIN_X + 1) *
                  (AR_FORMIC_MAX_Y - AR_FORMIC_MIN_Y + 1) *
                  (AR_FORMIC_MAX_Z - AR_FORMIC_MIN_Z + 1) * 
                                                AR_FORMIC_CORES_PER_BOARD +
                  (AR_ARM0_BID != -1) * AR_ARM0_CORES_PER_BOARD +
                  (AR_ARM1_BID != -1) * AR_ARM1_CORES_PER_BOARD;
    
    // Collect from all cores their "boot ok" message
    for (i = 0; i < total_cores - 1; i++) {
      word = ar_mbox_get(my_cid);
      if ((word >> 16) != 0x1111) {
        kt_printf("Boot master: invalid response 0x%08X\r\n", word);
        ar_abort();
      }
    }
    kt_printf("Collected %d boot answers\r\n\n", i);

    // Flush the UART, so multiple MicroBlaze cores can print debug stuff with
    // lower probability that the UART will get full
    ar_uart_flush();

    // Make Formic slaves reach their boot barrier
    for (x = AR_FORMIC_MIN_X; x <= AR_FORMIC_MAX_X; x++) {
      for (y = AR_FORMIC_MIN_Y; y <= AR_FORMIC_MAX_Y; y++) {
        for (z = AR_FORMIC_MIN_Z; z <= AR_FORMIC_MAX_Z; z++) {
          for (i = 0; i < AR_FORMIC_CORES_PER_BOARD; i++) {

            bid = (x << 4) | (y << 2) | z;

            if ((bid != my_bid) || (i != my_cid)) {
              ar_cnt_incr(my_cid, bid, i, NOC_COUNTER_WAKEUP1, 1);
            }
          }
        }
      }
    }

    // Make ARM slaves reach their boot barrier
    if (AR_ARM0_BID != -1) {
      for (i = 0; i < AR_ARM0_CORES_PER_BOARD; i++) {
        if ((AR_ARM0_BID != my_bid) || (i != my_cid)) {
          ar_cnt_incr(my_cid, AR_ARM0_BID, i, NOC_COUNTER_WAKEUP1, 1);
        }
      }
    }
    if (AR_ARM1_BID != -1) {
      for (i = 0; i < AR_ARM1_CORES_PER_BOARD; i++) {
        if ((AR_ARM1_BID != my_bid) || (i != my_cid)) {
          ar_cnt_incr(my_cid, AR_ARM1_BID, i, NOC_COUNTER_WAKEUP1, 1);
        }
      }
    }
  }

  // Boot slave
  else {

    // Initialize a boot barrier counter
    ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP1, -1);

    // Signal master we've booted
    ar_mbox_send(my_cid, AR_BOOT_MASTER_BID, 0, 
                 0x11110000 | (my_bid << 8) | my_cid);

    // Block on barrier, until master indicates all cores have booted
    while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP1)) {
      ;
    }
  }



  // =========================================================================
  // Bare-metal benchmarks. If enabled, they must run here because they
  // assume any kernel memory is available (and thus cannot run after 
  // kt_mem_init() or fmpi_rank_init() etc.
  // =========================================================================
  //test_bare();

  //if (!my_bid && !my_cid) {
  //  for (i = 100; i < 2000; i += 100) {
  //    stream_bare(i);
  //  }
  //  for (i = 2000; i <= 44000; i += 1000) {
  //    stream_bare(i);
  //  }
  //}

  // =========================================================================
  // Jacobi (bare-metal)
  // =========================================================================
  //while (1) {
  //  jacobi_bare(512, 1024, 256, 50);
  //}

  //jacobi_bare(512, 1024, 1024, 50);
  //jacobi_bare(256, 1024, 1024, 50);
  //jacobi_bare(128, 1024, 1024, 50);
  //jacobi_bare(64,  1024, 1024, 50);
  //jacobi_bare(32,  1024, 1024, 50);
  //jacobi_bare(16,  1024, 1024, 50);
  //jacobi_bare(8,   1024, 1024, 50);
  //jacobi_bare(4,   1024, 1024, 50);
  //jacobi_bare(2,   1024, 1024, 50);
  //jacobi_bare(1,   1024, 1024, 50);


  // =========================================================================
  // ================                                         ================
  // ========                                                         ========
  // ====                End of bare metal -- Begin of MPI                ====
  // ========                                                         ========
  // ================                                         ================
  // =========================================================================
  //while (1) {
  //  ;
  //}

  
  // =========================================================================
  // Memory management and FMPI initialization
  // =========================================================================

  // Initialize kernel memory allocator
  kt_mem_init();

  // We now have a context
  context = mm_get_context(my_cid);

  // Assign FMPI ranks
  //ret = fmpi_select_rank_init(FMPI_CFG_1_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_2_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_4_ARM);
  //ret = fmpi_select_rank_init(FMPI_CFG_4_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_8_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_8_HET);
  //ret = fmpi_select_rank_init(FMPI_CFG_12_HET);
  //ret = fmpi_select_rank_init(FMPI_CFG_16_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_20_HET);
  //ret = fmpi_select_rank_init(FMPI_CFG_64_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_128_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_256_MB);
  ret = fmpi_select_rank_init(FMPI_CFG_512_MB);
  //ret = fmpi_select_rank_init(FMPI_CFG_516_HET);

  // For video demo, enable ARM0 and uncomment this
  if (my_bid == AR_ARM0_BID) {
    demo_mpi();
  }

  // If we're not part of the FMPI setup, block here.
  if (ret) {
    while (1) {
      ;
    }
  }


  // Print a separator 
  if (!context->fmpi->rank) {
    kt_printf("=========================================================\r\n");
  }

  // Initialize MPI
  MPI_Init(NULL, NULL);
 

  // =========================================================================
  // Video demo (MPI). Don't forget to enable ARM0 cores and use the video
  // clause above, to make them enter demo_mpi() on their own.
  // =========================================================================
  demo_mpi();


  // =========================================================================
  // Matrix multiplication (MPI)
  // =========================================================================

  // 128x128 with verification
  //
  // matrix_mult_mpi(1,    1,  128, 1);
  // matrix_mult_mpi(4,    2,   64, 1);
  // matrix_mult_mpi(16,   4,   32, 1);
  // matrix_mult_mpi(64,   8,   16, 1);
  // matrix_mult_mpi(256, 16,    8, 0); // FIXME: verif=1 "out of counters"

  // 256x256, no verification
  //
  // matrix_mult_mpi(256, 16,  16, 0);
  // matrix_mult_mpi(64,  8,   32, 0);
  // matrix_mult_mpi(16,  4,   64, 0);
  // matrix_mult_mpi(4,   2,  128, 0);
  // matrix_mult_mpi(1,   1,  256, 0);

  // 384x384, no verification
  //
  // matrix_mult_mpi(256, 16,  24, 0);
  // matrix_mult_mpi(64,  8,   48, 0);
  // matrix_mult_mpi(16,  4,   96, 0);
  // matrix_mult_mpi(4,   2,  192, 0);
  // matrix_mult_mpi(1,   1,  384, 0);

  // Strong scaling
  //
  // matrix_mult_mpi(1,   1,  768, 0); // 4161.475 sec (296.011 + 9*ovflow)
  // matrix_mult_mpi(4,   2,  384, 0); //  439.735 sec (10.239 + ovflow)
  // matrix_mult_mpi(16,  4,  192, 0); //   86.629 sec
  // matrix_mult_mpi(64,  8,   96, 0); //   21.628 sec
  // matrix_mult_mpi(256, 16,  48, 0); //    4.421 sec

  // Weak scaling
  //
  // matrix_mult_mpi(1,   1,  48, 0); //  0.272 sec
  // matrix_mult_mpi(4,   2,  48, 0); //  0.549 sec
  // matrix_mult_mpi(16,  4,  48, 0); //  1.096 sec
  // matrix_mult_mpi(64,  8,  48, 0); //  2.196 sec
  // matrix_mult_mpi(256, 16, 48, 0); //  4.424 sec

  // =========================================================================
  // Jacobi (MPI)
  // =========================================================================
  // jacobi_mpi(512, 1024, 1024, 50);
  // jacobi_mpi(256, 1024, 1024, 50);
  // jacobi_mpi(128, 1024, 1024, 50);
  // jacobi_mpi(64,  1024, 1024, 50);
  // jacobi_mpi(32,  1024, 1024, 50);
  // jacobi_mpi(16,  1024, 1024, 50);
  // jacobi_mpi(8,   1024, 1024, 50);
  // jacobi_mpi(4,   1024, 1024, 50);
  // jacobi_mpi(2,   1024, 1024, 50);
  // jacobi_mpi(1,   1024, 1024, 50);
 
  // Strong scaling
  //
  // jacobi_mpi(4,   4096, 1024, 16); // 117.558 sec
  // jacobi_mpi(8,   4096, 1024, 16); //  58.817 sec
  // jacobi_mpi(16,  4096, 1024, 16); //  29.439 sec
  // jacobi_mpi(32,  4096, 1024, 16); //  14.724 sec
  // jacobi_mpi(64,  4096, 1024, 16); //   7.245 sec
  // jacobi_mpi(128, 4096, 1024, 16); //   3.341 sec
  // jacobi_mpi(256, 4096, 1024, 16); //   1.628 sec
  // jacobi_mpi(512, 4096, 1024, 16); //   0.823 sec
 
  // Weak scaling
  //
  // jacobi_mpi(4,     32, 1024, 16); // 0.820 sec
  // jacobi_mpi(8,     64, 1024, 16); // 0.821 sec
  // jacobi_mpi(16,   128, 1024, 16); // 0.822 sec
  // jacobi_mpi(32,   256, 1024, 16); // 0.823 sec
  // jacobi_mpi(64,   512, 1024, 16); // 0.823 sec
  // jacobi_mpi(128, 1024, 1024, 16); // 0.823 sec
  // jacobi_mpi(256, 2048, 1024, 16); // 0.823 sec
  // jacobi_mpi(512, 4096, 1024, 16); // 0.823 sec
 
  // =========================================================================
  // Smith-Waterman (MPI)
  // =========================================================================
  // swat_mpi(512, 4096, 8192, 16, 4);
  // swat_mpi(256, 4096, 8192, 16, 4);
  // swat_mpi(128, 4096, 8192, 16, 4);
  // swat_mpi(64,  4096, 8192, 16, 8);
  // swat_mpi(32,  4096, 8192, 16, 16);
  // swat_mpi(16,  4096, 8192, 16, 32);
  // swat_mpi(8,   4096, 8192, 16, 64);
  // swat_mpi(4,   4096, 8192, 16, 128);
  // swat_mpi(2,   4096, 8192, 16, 128);
  // swat_mpi(1,   4096, 8192, 16, 128);
  
  // =========================================================================
  // Bitonic sort (MPI)
  // =========================================================================
  // bitonic_mpi(512, 524288, 1);
  // bitonic_mpi(256, 524288, 1);
  // bitonic_mpi(128, 524288, 1);
  // bitonic_mpi(64,  524288, 1);
  // bitonic_mpi(32,  524288, 1);
  // bitonic_mpi(16,  524288, 1);
  // bitonic_mpi(8,   524288, 1);
  // bitonic_mpi(4,   524288, 1);
  // bitonic_mpi(2,   524288, 1);
  // bitonic_mpi(1,   524288, 1);
  
  // Strong scaling
  //
  // bitonic_mpi(4,   2097152, 0); // 234.978 sec
  // bitonic_mpi(8,   2097152, 0); // 111.449 sec
  // bitonic_mpi(16,  2097152, 0); //  52.418 sec
  // bitonic_mpi(32,  2097152, 0); //  22.651 sec
  // bitonic_mpi(64,  2097152, 0); //  11.152 sec
  // bitonic_mpi(128, 2097152, 0); //   5.074 sec
  // bitonic_mpi(256, 2097152, 0); //   2.585 sec
  // bitonic_mpi(512, 2097152, 0); //   1.344 sec
  
  // Weak scaling
  //
  // bitonic_mpi(4,     16384, 0); // 0.718 sec
  // bitonic_mpi(8,     32768, 0); // 0.761 sec
  // bitonic_mpi(16,    65536, 0); // 0.820 sec
  // bitonic_mpi(32,   131072, 0); // 0.893 sec
  // bitonic_mpi(64,   262144, 0); // 0.980 sec
  // bitonic_mpi(128,  524288, 0); // 1.089 sec
  // bitonic_mpi(256, 1048576, 0); // 1.205 sec
  // bitonic_mpi(512, 2097152, 0); // 1.343 sec
  
  // =========================================================================
  // 2D-FFT (MPI)
  // =========================================================================
  // Strong scaling
  //
  // fft_mpi(512, 2048);
  // fft_mpi(256, 2048);
  // fft_mpi(128, 2048);
  // fft_mpi(64,  2048);
  // fft_mpi(32,  2048);
  // fft_mpi(16,  2048);
  // fft_mpi(8,   2048);
  // fft_mpi(4,   2048);

  // Weak scaling
  //
  // fft_mpi(512, 8192);
  // fft_mpi(128, 4096);
  // fft_mpi(32,  2048);
  // fft_mpi(8,   1024);
  // fft_mpi(2,   512);


  // =========================================================================
  // c-ray (MPI)
  // =========================================================================
  
  // Strong scaling
  //
  // cray_mpi(1,   800, 512); // 502.107 sec (72.611 + ovflow)
  // cray_mpi(2,   800, 512); // 284.682 sec
  // cray_mpi(4,   800, 512); // 158.388 sec
  // cray_mpi(8,   800, 512); //  80.068 sec
  // cray_mpi(16,  800, 512); //  44.095 sec
  // cray_mpi(32,  800, 512); //  22.583 sec
  // cray_mpi(64,  800, 512); //  11.363 sec
  // cray_mpi(128, 800, 512); //   5.755 sec
  // cray_mpi(256, 800, 512); //   2.889 sec
  // cray_mpi(512, 800, 512); //   1.462 sec

  // Weak scaling
  //
  // cray_mpi(4,   800,   4); // 1.260 sec
  // cray_mpi(8,   800,   8); // 1.290 sec
  // cray_mpi(16,  800,  16); // 1.342 sec
  // cray_mpi(32,  800,  32); // 1.414 sec
  // cray_mpi(64,  800,  64); // 1.439 sec
  // cray_mpi(128, 800, 128); // 1.440 sec
  // cray_mpi(256, 800, 256); // 1.444 sec
  // cray_mpi(512, 800, 512); // 1.451 sec


  // =========================================================================
  // k-means (MPI)
  // =========================================================================

  // Strong scaling
  //
  // kmeans_mpi(4,   16, 2097152, 8); // 630.375 sec (200.879 + ovflow)
  // kmeans_mpi(8,   16, 2097152, 8); // 315.303 sec
  // kmeans_mpi(16,  16, 2097152, 8); // 157.759 sec
  // kmeans_mpi(32,  16, 2097152, 8); //  78.980 sec
  // kmeans_mpi(64,  16, 2097152, 8); //  39.530 sec
  // kmeans_mpi(128, 16, 2097152, 8); //  19.746 sec
  // kmeans_mpi(256, 16, 2097152, 8); //   9.952 sec
  // kmeans_mpi(512, 16, 2097152, 8); //   5.085 sec

  // Weak scaling
  //
  // kmeans_mpi(4,   16,   16384, 8); // 5.075 sec
  // kmeans_mpi(8,   16,   32768, 8); // 5.074 sec
  // kmeans_mpi(16,  16,   65536, 8); // 5.074 sec
  // kmeans_mpi(32,  16,  131072, 8); // 5.075 sec
  // kmeans_mpi(64,  16,  262144, 8); // 5.074 sec
  // kmeans_mpi(128, 16,  524288, 8); // 5.078 sec
  // kmeans_mpi(256, 16, 1048576, 8); // 5.079 sec
  // kmeans_mpi(512, 16, 2097152, 8); // 5.081 sec


  // =========================================================================
  // Sandia MPI Benchmark suite
  // =========================================================================
 
  // Host overhead test (needs CFG_2_MB)
  //for (i = 1; i <= 1024 * 1024; i = i * 2) { 
  //  smb_overhead_mpi(0, 20, i, 1.5, 1.02, !(i == 1), 0);
  //}

  // Message rate test (needs CFG_8_MB)
  //if (!my_bid && !my_cid) kt_printf("\n>> Varying number of messages\r\n");
  //for (i = 1; i <= 16; i *= 2) {
  //  smb_msgrate_mpi(4, 8, i,  64, 256 * 1024 / 4, 1, 0);
  //}

  //if (!my_bid && !my_cid) kt_printf("\n>> Varying number of bytes\r\n");
  //for (i = 64; i <= 128 * 1024; i *= 2) {
  //  smb_msgrate_mpi(4, 8, 8, i, 256 * 1024 / 4, 1, 0);
  //}

  // =========================================================================
  // MPI self-test routine
  // =========================================================================
  //test_mpi();




  // Stall
  if (!context->fmpi->rank) {
    kt_printf("\nAll done.\r\n");
  
    // Instruct server to kill the power
    ar_uart_flush();
    ar_timer_busy_wait_msec(200);
    kt_printf("%s", DBG_KILL_CONNECTION);
  }
  while (1) {
    ;
  }

  
}

