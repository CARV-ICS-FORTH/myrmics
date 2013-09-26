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
// Abstract      : ARM CompactFlash device driver
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: flash.c,v $
// CVS revision  : $Revision: 1.6 $
// Last modified : $Date: 2012/03/15 15:31:47 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <arm_regs.h>
#include <types.h>
#include <kernel_toolset.h>
#include <memory_management.h>




// ===========================================================================
// ar_flash_wait_data_req()     Busy-wait for DRQ to be set
// ===========================================================================
// * INPUTS
//   unsigned int msec_timeout  Timeout period to wait in milliseconds
//
// * RETURN VALUE
//   int                        0 for success, 1 otherwise
// ===========================================================================
int ar_flash_wait_data_req(unsigned int msec_timeout) {
  
  // See if it's already ready without setting the timer
  if (*ARM_FLASH_ACCESS_STATUS & 0x8) {
    return 0;
  }

  // Start the timer with the timeout
  ar_timer0_set(msec_timeout * 200000);

  // Wait for ready or timeout
  while (ar_timer0_read() != 0) {
    if (*ARM_FLASH_ACCESS_STATUS & 0x8) {
      return 0;
    }
  }

  // Failure
  return 1;
}


// ===========================================================================
// ar_flash_wait_ready()        Busy-wait for Ready to be set
// ===========================================================================
// * INPUTS
//   unsigned int msec_timeout  Timeout period to wait in milliseconds
//
// * RETURN VALUE
//   int                        0 for success, 1 otherwise
// ===========================================================================
int ar_flash_wait_ready(unsigned int msec_timeout) {
  
  // See if it's already ready without setting the timer
  if (*ARM_FLASH_ACCESS_STATUS & 0x40) {
    return 0;
  }

  // Start the timer with the timeout
  ar_timer0_set(msec_timeout * 200000);

  // Wait for ready or timeout
  while (ar_timer0_read() != 0) {
    if (*ARM_FLASH_ACCESS_STATUS & 0x40) {
      return 0;
    }
  }

  // Failure
  return 1;
}


// ===========================================================================
// ar_flash_wait_busy()         Busy-wait for Busy to be set
// ===========================================================================
// * INPUTS
//   unsigned int msec_timeout  Timeout period to wait in milliseconds
//
// * RETURN VALUE
//   int                        0 for success, 1 otherwise
// ===========================================================================
int ar_flash_wait_busy(unsigned int msec_timeout) {
  
  // See if it's already busy without setting the timer
  if (*ARM_FLASH_ACCESS_STATUS & 0x80) {
    return 0;
  }

  // Start the timer with the timeout
  ar_timer0_set(msec_timeout * 200000);

  // Wait for busy or timeout
  while (ar_timer0_read() != 0) {
    if (*ARM_FLASH_ACCESS_STATUS & 0x80) {
      return 0;
    }
  }

  // Failure
  return 1;
}


// ===========================================================================
// ar_flash_wait_not_busy()     Busy-wait for Busy to be cleared
// ===========================================================================
// * INPUTS
//   unsigned int msec_timeout  Timeout period to wait in milliseconds
//
// * RETURN VALUE
//   int                        0 for success, 1 otherwise
// ===========================================================================
int ar_flash_wait_not_busy(unsigned int msec_timeout) {
  
  // See if it's already not busy without setting the timer
  if (!(*ARM_FLASH_ACCESS_STATUS & 0x80)) {
    return 0;
  }

  // Start the timer with the timeout
  ar_timer0_set(msec_timeout * 200000);

  // Wait for busy or timeout
  while (ar_timer0_read() != 0) {
    if (!(*ARM_FLASH_ACCESS_STATUS & 0x80)) {
      return 0;
    }
  }

  // Failure
  return 1;
}


// ===========================================================================
// ar_flash_init_device()       Initializes the CompactFlash card and 
//                              queries it.
// ===========================================================================
// * RETURN VALUE
//   int                        Number of 512B-sectors on the card, or 
//                              -1 on failure
// ===========================================================================
int ar_flash_init_device() {
  unsigned int buf[256];
  int i;
  int sectors;
  Context *context;

  kt_printf("Resetting CompactFlash peripheral...\r\n");

  // Reset on, power off  
  *ARM_FLASH_CONTROL = 0x40000;
  ar_timer0_busy_wait(1000);
  // Power on  
  *ARM_FLASH_CONTROL = 0x40001;
  ar_timer0_busy_wait(100);
  // Reset off
  *ARM_FLASH_CONTROL = 0x40003;
  ar_timer0_busy_wait(100);


  kt_printf("Identifying flash card...\r\n");

  // Wait for ready
  if (ar_flash_wait_ready(1000)) {
    kt_printf("Flash ready timed out.\r\n");
    return -1;
  }

  // Issue Identify Device command
  *ARM_FLASH_ACCESS_DEVICE  = 0xE0;
  *ARM_FLASH_ACCESS_NSECT   = 0x01;
  *ARM_FLASH_ACCESS_LBAL    = 0x01;
  *ARM_FLASH_ACCESS_LBAM    = 0x00;
  *ARM_FLASH_ACCESS_LBAH    = 0x00;
  *ARM_FLASH_ACCESS_COMMAND = 0xEC;

  // Wait for DRQ
  if (ar_flash_wait_data_req(1000)) {
    kt_printf("Flash Identify timed out.\r\n");
    return -1;
  }

  // Check for error
  if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
    kt_printf("Flash Identify gave an error.\r\n");
    return -1;
  }

  // Read identification sector
  for (i = 0; i < 256; i++) {
    buf[i] = *ARM_FLASH_ACCESS_DATA & 0xFFFF;
  }

  // Collect & print out useful stuff from the ID process
  sectors = (buf[61] << 16) | buf[60];

  kt_printf("  Model nr:    ");
  for (i = 27; i < 47; i++) {
    kt_printf("%c", buf[i] >> 8);
    kt_printf("%c", buf[i] & 0xFF);
  }
  kt_printf("\r\n");

  kt_printf("  Serial nr:   ");
  for (i = 10; i < 20; i++) {
    kt_printf("%c", buf[i] >> 8);
    kt_printf("%c", buf[i] & 0xFF);
  }
  kt_printf("\r\n");

  kt_printf("  Firmware nr: ");
  for (i = 23; i < 27; i++) {
    kt_printf("%c", buf[i] >> 8);
    kt_printf("%c", buf[i] & 0xFF);
  }
  kt_printf("\r\n");

  kt_printf("  LBA sectors: %d\r\n", sectors);
  kt_printf("  Capacity:    %d MiB\r\n", sectors / 2048);

  // Store capacity for filesystem usage
  context = mm_get_context(ar_get_core_id());
  context->fs_flash_num_sectors = sectors;

  // Success
  return sectors;
}


// ===========================================================================
// ar_flash_read_sector()       Reads a 512-B sector from the Compact Flash
// ===========================================================================
// * INPUT ARGUMENTS
//   int sector                 LBA address of the sector to be read
//
// * OUTPUT ARGUMENTS
//   u8 *buf                    512 byte buffer to be filled with read data
//
// * RETURN VALUE
//   int                        0 on success
// ===========================================================================
int ar_flash_read_sector(int sector, unsigned char *buf) {
  int i;
  unsigned int tmp;
  
  // Wait for ready
  if (ar_flash_wait_ready(1000)) {
    kt_printf("Flash ready timed out.\r\n");
    return -1;
  }

  // Issue Read Sector command
  *ARM_FLASH_ACCESS_DEVICE  = 0xE0 | ((sector & 0xF000000) >> 24);
  *ARM_FLASH_ACCESS_LBAH    = (sector & 0xFF0000) >> 16;
  *ARM_FLASH_ACCESS_LBAM    = (sector & 0xFF00) >> 8;
  *ARM_FLASH_ACCESS_LBAL    = sector & 0xFF;
  *ARM_FLASH_ACCESS_NSECT   = 0x01;
  *ARM_FLASH_ACCESS_COMMAND = 0x20;

  // Wait for DRQ
  if (ar_flash_wait_data_req(1000)) {
    kt_printf("Flash DRQ timed out.\r\n");
    return -1;
  }

  // Check for error
  if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
    kt_printf("Flash read sector gave an error.\r\n");
    return -1;
  }

  // Read the sector
  for (i = 0; i < 256; i++) {
    tmp = *ARM_FLASH_ACCESS_DATA & 0xFFFF;
    buf[2 * i]     = tmp >> 8;
    buf[2 * i + 1] = tmp & 0xFF;
  }

  return 0;
}



// ===========================================================================
// ar_flash_write_sector()      Writes (after erasing) a 512-byte sector to 
//                              the CompactFlash card and verifies it
// ===========================================================================
// * INPUT ARGUMENTS
//   int sector                 LBA address of the sector to be written
//   u8 *buf                    512 byte buffer to be written
//
// * RETURN VALUE
//   int                        0 on success
// ===========================================================================
int ar_flash_write_sector(int sector, unsigned char *buf) {
  int i;
  
  // Wait for ready
  if (ar_flash_wait_ready(1000)) {
    kt_printf("Flash ready timed out.\r\n");
    return -1;
  }

  // Issue Write Verify Sector command
  *ARM_FLASH_ACCESS_DEVICE  = 0xE0 | ((sector & 0xF000000) >> 24);
  *ARM_FLASH_ACCESS_LBAH    = (sector & 0xFF0000) >> 16;
  *ARM_FLASH_ACCESS_LBAM    = (sector & 0xFF00) >> 8;
  *ARM_FLASH_ACCESS_LBAL    = sector & 0xFF;
  *ARM_FLASH_ACCESS_NSECT   = 0x01;
  *ARM_FLASH_ACCESS_COMMAND = 0x3C;

  // Wait for DRQ
  if (ar_flash_wait_data_req(1000)) {
    kt_printf("Flash DRQ timed out.\r\n");
    return -1;
  }

  // Check for error
  if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
    kt_printf("Flash write sector gave an early error.\r\n");
    return -1;
  }

  // Write the sector
  for (i = 0; i < 256; i++) {
    *ARM_FLASH_ACCESS_DATA = (buf[2*i] << 8) | buf[2*i + 1];
  }

  // Wait for not busy
  if (ar_flash_wait_not_busy(1000)) {
    kt_printf("Flash not_busy timed out.\r\n");
    return -1;
  }

  // Check for error
  if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
    kt_printf("Flash write sector gave a late error.\r\n");
    return -1;
  }

  return 0;
}


// ===========================================================================
// ar_flash_read_sectors()      Reads a number of 512-byte sectors from the 
//                              Compact Flash
// ===========================================================================
// * INPUT ARGUMENTS
//   int start_sector           LBA address of the first sector to be read
//   int num_sectors            Number of sectors to read
//
// * OUTPUT ARGUMENTS
//   unsigned char *buf         (num_sectors * 512) byte buffer to be 
//                              filled with read data
//
// * RETURN VALUE
//   int                        0 on success
// ===========================================================================
int ar_flash_read_sectors(int start_sector, int num_sectors,
                          unsigned char *buf) {
  int i, j;
  unsigned int tmp;
  
  // Sanity check
  if ((num_sectors < 1) || (num_sectors > 256)) {
    return -1;
  }
  // Wait for ready
  if (ar_flash_wait_ready(1000)) {
    kt_printf("Flash ready timed out.\r\n");
    return -1;
  }

  // Issue Read Sector command
  *ARM_FLASH_ACCESS_DEVICE  = 0xE0 | ((start_sector & 0xF000000) >> 24);
  *ARM_FLASH_ACCESS_LBAH    = (start_sector & 0xFF0000) >> 16;
  *ARM_FLASH_ACCESS_LBAM    = (start_sector & 0xFF00) >> 8;
  *ARM_FLASH_ACCESS_LBAL    = start_sector & 0xFF;
  *ARM_FLASH_ACCESS_NSECT   = (num_sectors == 256) ? 0 : num_sectors;
  *ARM_FLASH_ACCESS_COMMAND = 0x20;

  // For each sector
  for (j = 0; j < num_sectors; j++) {

    // Wait for DRQ
    if (ar_flash_wait_data_req(1000)) {
      kt_printf("Flash DRQ timed out.\r\n");
      return -1;
    }

    // Check for error
    if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
      kt_printf("Flash read sector gave an error.\r\n");
      return -1;
    }

    // Read the sector
    for (i = 0; i < 256; i++) {
      tmp = *ARM_FLASH_ACCESS_DATA & 0xFFFF;
      buf[512*j + 2*i]     = tmp >> 8;
      buf[512*j + 2*i + 1] = tmp & 0xFF;
    }
  }

  return 0;
}


// ===========================================================================
// ar_flash_write_sectors()     Writes (after erasing) a number of 512-byte 
//                              sectors to the CompactFlash card and verifies
//                              them.
// ===========================================================================
// * INPUT ARGUMENTS
//   int start_sector           LBA address of the first sector to be written
//   int num_sectors            Number of sectors to write
//   unsigned char *buf         (num_sectors * 512) bytes buffer to be written
//
// * RETURN VALUE
//   int                        0 on success
// ===========================================================================
int ar_flash_write_sectors(int start_sector, int num_sectors, 
                           unsigned char *buf) {
  int i, j;
  
  // Sanity check
  if ((num_sectors < 1) || (num_sectors > 256)) {
    return -1;
  }

  // Wait for ready
  if (ar_flash_wait_ready(1000)) {
    kt_printf("Flash ready timed out.\r\n");
    return -1;
  }

  // Issue Write Verify Sector command
  *ARM_FLASH_ACCESS_DEVICE  = 0xE0 | ((start_sector & 0xF000000) >> 24);
  *ARM_FLASH_ACCESS_LBAH    = (start_sector & 0xFF0000) >> 16;
  *ARM_FLASH_ACCESS_LBAM    = (start_sector & 0xFF00) >> 8;
  *ARM_FLASH_ACCESS_LBAL    = start_sector & 0xFF;
  *ARM_FLASH_ACCESS_NSECT   = (num_sectors == 256) ? 0 : num_sectors;
  *ARM_FLASH_ACCESS_COMMAND = 0x3C;

  // For each sector
  for (j = 0; j < num_sectors; j++) {

    // Wait for DRQ
    if (ar_flash_wait_data_req(1000)) {
      kt_printf("Flash DRQ timed out.\r\n");
      return -1;
    }

    // Check for error
    if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
      kt_printf("Flash write sector gave an early error.\r\n");
      return -1;
    }

    // Write the sector
    for (i = 0; i < 256; i++) {
      *ARM_FLASH_ACCESS_DATA = (buf[512*j + 2*i] << 8) | buf[512*j + 2*i + 1];
    }
  }

  // Wait for not busy
  if (ar_flash_wait_not_busy(1000)) {
    kt_printf("Flash not_busy timed out.\r\n");
    return -1;
  }

  // Check for error
  if (*ARM_FLASH_ACCESS_STATUS & 0x1) {
    kt_printf("Flash write sector gave a late error.\r\n");
    return -1;
  }

  // Success
  return 0;
}

