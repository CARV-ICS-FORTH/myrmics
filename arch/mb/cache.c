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
// Abstract      : Microblaze CPU/Cache/ART arch-specific functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: cache.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/12/12 16:37:19 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <mbs_regs.h>
#include <formic_regs.h>



// ===========================================================================
// ar_art_install_region()      Installs an ART region
// ===========================================================================
// * INPUTS
//   int region                 ART region (1 to 4), 0 is reserved for I/O
//   unsigned int base          Base address (must be 1-MB aligned)
//   unsigned int size          Region size in bytes (must be 1-MB aligned)
//   int cacheable              When 1, region is cacheable
//   int read_only              When 1, region is read-only
//   int executable             When 1, region can be used for code
//   int user_accessible        When 1, region is accessible from user mode
//   int privileged_accessible  When 1, region is accessible from privileged
//                              mode
// ===========================================================================
void ar_art_install_region(int region, unsigned int base, unsigned int size,
                           int cacheable, int read_only, int executable,
                           int user_accessible, int privileged_accessible) {
  unsigned int value;

  
  // Compute entry
  value = (cacheable                      << 31) | 
          (read_only                      << 29) |
          (executable                     << 28) |
          (user_accessible                << 27) |
          (privileged_accessible          << 26) |
          (1                              << 24) | // valid bit
          (((base + size) / 0x100000 - 1) << 12) | // bound
          ((base / 0x100000)              << 0);   // base
  
  // Install entry
  switch (region) {
    case 1: *MBS_ART_ENTRY1 = value; break;
    case 2: *MBS_ART_ENTRY2 = value; break;
    case 3: *MBS_ART_ENTRY3 = value; break;
    case 4: *MBS_ART_ENTRY4 = value; break;
    default: ar_abort();
  }
}


// ===========================================================================
// ar_art_uninstall_region()    Disable an ART region (set to invalid)
// ===========================================================================
// * INPUTS
//   int region                 ART region (1 to 4), 0 is reserved for I/O
// ===========================================================================
void ar_art_uninstall_region(int region) {
  switch (region) {
    case 1: *MBS_ART_ENTRY1 = 0; break;
    case 2: *MBS_ART_ENTRY2 = 0; break;
    case 3: *MBS_ART_ENTRY3 = 0; break;
    case 4: *MBS_ART_ENTRY4 = 0; break;
    default: ar_abort();
  }
}


// ===========================================================================
// ar_cache_enable()            Enable the cache hierarchy
// ===========================================================================
// * INPUTS
//   int epoch_mode             0: Set L2 cache in LRU mode
//                              1: Set L2 cache in Epoch mode
// ===========================================================================
void ar_cache_enable(int epoch_mode) {

  // Clear caches
  *MBS_CHE_MAINT = 0x10101;

  // Enable caches
  if (epoch_mode) {
    *MBS_CHE_CONTROL = (AR_CACHE_EPOCH_CPU_WAYS << 24) | // CPU ways
                       (0 << 16) |                       // Epoch number
                       (0 << 4)  |                       // LRU
                       (1 << 0);                         // Enable
  }
  else {
    *MBS_CHE_CONTROL = (0 << 24) |  // CPU ways
                       (0 << 16) |  // Epoch number
                       (1 << 4)  |  // LRU
                       (1 << 0);    // Enable
  }
}


// ==========================================================================
// ar_cache_disable()           Disable the cache hierarchy
// ==========================================================================
void ar_cache_disable() {

  // Disable caches
  *MBS_CHE_CONTROL = 0x0;

  // Flush L2
  *MBS_CHE_MAINT = 0x20000;
}


// ==========================================================================
// ar_cache_flush()            Flush the cache hierarchy
// ==========================================================================
void ar_cache_flush() {

  // Flush L2
  *MBS_CHE_MAINT = 0x20000;
}

