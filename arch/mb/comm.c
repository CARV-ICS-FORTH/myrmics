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
// Abstract      : Microblaze communication-specific functions
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: comm.c,v $
// CVS revision  : $Revision: 1.10 $
// Last modified : $Date: 2013/03/19 14:52:54 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <mbs_regs.h>
#include <formic_regs.h>
#include <arch.h>
#include <noc.h>
#include <memory_management.h>


// ==========================================================================
// Do a remote read of a Formic board link status register using the mailslot
// ==========================================================================
unsigned int ar_read_link_status(int my_bid, int my_cid, int src_bid) {
  *MBS_MNI_SRC_ADR  = (int) FORMIC_BRD_LINK_STATUS;
  *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
  *MBS_MNI_DMA_SIZE = 4;
  *MBS_MNI_BRD_NODE = (src_bid << 20) | (0xF << 16);
  *MBS_MNI_OPCODE   = (my_bid << 20) | (my_cid << 16) | 0x2;

  return (*MBS_MSL_ACCESS);
}


// ==========================================================================
// Write a Formic board link status register
// ==========================================================================
void ar_write_link_status(int my_cid, int dst_bid, unsigned int val) {
  *MBS_MNI_MSG0    = val;
  *MBS_MNI_DST_ADR = (int) FORMIC_BRD_LINK_STATUS;
  *MBS_MNI_OPCODE  = (dst_bid << 20) | (0xF << 16);
}


// ==========================================================================
// Do a remote read of a Formic board status register using the mailslot
// ==========================================================================
unsigned int ar_read_status(int my_bid, int my_cid, int src_bid) {
  *MBS_MNI_SRC_ADR  = (int) FORMIC_BRD_STATUS;
  *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
  *MBS_MNI_DMA_SIZE = 4;
  *MBS_MNI_BRD_NODE = (src_bid << 20) | (0xF << 16);
  *MBS_MNI_OPCODE   = (my_bid << 20) | (my_cid << 16) | 0x2;

  return (*MBS_MSL_ACCESS);
}


// ==========================================================================
// Write a Formic board control register
// ==========================================================================
void ar_write_brd_control(int my_cid, int dst_bid, unsigned int val) {
  *MBS_MNI_MSG0    = val;
  *MBS_MNI_DST_ADR = (int) FORMIC_BRD_CONTROL;
  *MBS_MNI_OPCODE  = (dst_bid << 20) | (0xF << 16);
}


// ==========================================================================
// Wake up a remote Formic MicroBlaze core
// ==========================================================================
void ar_wake_up_core(int my_bid, int my_cid, int dst_bid, int dst_cid) {

  // Set counter
  ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP1, -4);

  // Send message with ack
  *MBS_MNI_MSG0     = 0x101;
  *MBS_MNI_DST_ADR  = (int) MBS_CPU_CONTROL;
  *MBS_MNI_ACK_ADR  = (int) MBS_CNT_VAL(NOC_COUNTER_WAKEUP1);
  *MBS_MNI_BRD_NODE = (my_bid  << 4) | (my_cid << 0);
  *MBS_MNI_OPCODE   = (dst_bid << 20) | (dst_cid << 16) | 0x4;

  // Wait until message is received
  while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP1)) {
    ;
  }
}


// ==========================================================================
// Disable a remote Formic MicroBlaze core
// ==========================================================================
void ar_suspend_core(int my_bid, int my_cid, int dst_bid, int dst_cid) {

  // Set counter
  ar_cnt_set(my_cid, NOC_COUNTER_WAKEUP1, -4);

  // Send message with ack
  *MBS_MNI_MSG0     = 0x100;
  *MBS_MNI_DST_ADR  = (int) MBS_CPU_CONTROL;
  *MBS_MNI_ACK_ADR  = (int) MBS_CNT_VAL(NOC_COUNTER_WAKEUP1);
  *MBS_MNI_BRD_NODE = (my_bid  << 4) | (my_cid << 0);
  *MBS_MNI_OPCODE   = (dst_bid << 20) | (dst_cid << 16) | 0x4;

  // Wait until message is received
  while (ar_cnt_get(my_cid, NOC_COUNTER_WAKEUP1)) {
    ;
  }
}


// ==========================================================================
// Do a remote read of a XUP video peripheral register using the mailslot
// ==========================================================================
unsigned int ar_read_xup_register(int my_bid, int my_cid, int xup_reg_adr) {
  ar_assert(AR_XUP_BID > -1);

  *MBS_MNI_SRC_ADR  = xup_reg_adr;
  *MBS_MNI_DST_ADR  = (int) MBS_MSL_ACCESS;
  *MBS_MNI_DMA_SIZE = 4;
  *MBS_MNI_BRD_NODE = (AR_XUP_BID << 20) | (0x0 << 16);
  *MBS_MNI_OPCODE   = (my_bid << 20) | (my_cid << 16) | 0x2;

  return (*MBS_MSL_ACCESS);
}


// ==========================================================================
// Write a XUP video peripheral register
// ==========================================================================
void ar_write_xup_register(int my_cid, int xup_reg_adr, unsigned int val) {
  *MBS_MNI_MSG0    = val;
  *MBS_MNI_DST_ADR = xup_reg_adr;
  *MBS_MNI_OPCODE  = (AR_XUP_BID << 20) | (0x0 << 16);
}


// ==========================================================================
// Send a 800x600 video buffer to a XUP DRAM location.  We use a special
// routine for this transfer, to avoid expanding the ar_virt_to_phys() routine
// to exclude translation for the XUP boards at runtime.
// ==========================================================================
void ar_send_xup_frame(int my_bid, int my_cid, unsigned int my_adr, 
                       unsigned int xup_adr, int my_cnt) {

  // Set counter
  ar_cnt_set(my_cid, my_cnt, -800 * 600 * 4);

  // Wait until DMA engine can support 1 more DMA
  while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
    ;
  }

  // DMA
  *MBS_MNI_SRC_ADR  = ar_virt_to_phys(my_bid, my_adr);
  *MBS_MNI_DST_ADR  = xup_adr;
  *MBS_MNI_ACK_ADR  = ar_addr_of_cnt(my_bid, my_cid, my_cnt);
  *MBS_MNI_DMA_SIZE = 800 * 600 * 4;
  *MBS_MNI_BRD_NODE = (my_bid << 20) | (my_cid << 16) |
                      (my_bid  << 4) | (my_cid << 0);
  *MBS_MNI_OPCODE   = (AR_XUP_BID << 20) | (0x0 << 16) | (1 << 2) | 0x2;
}


// ==========================================================================
// Receive a 800x600 video buffer from a XUP DRAM location. See 
// ar_send_xup_frame() above for comments.
// ==========================================================================
void ar_receive_xup_frame(int my_bid, int my_cid, unsigned int my_adr, 
                          unsigned int xup_adr, int my_cnt) {

  // Set counter
  ar_cnt_set(my_cid, my_cnt, -800 * 600 * 4);

  // DMA
  *MBS_MNI_SRC_ADR  = xup_adr;
  *MBS_MNI_DST_ADR  = ar_virt_to_phys(my_bid, my_adr);
  *MBS_MNI_ACK_ADR  = ar_addr_of_cnt(my_bid, my_cid, my_cnt);
  *MBS_MNI_DMA_SIZE = 800 * 600 * 4;
  *MBS_MNI_BRD_NODE = (AR_XUP_BID << 20) | (0x0 << 16) |
                      (my_bid  << 4) | (my_cid << 0);
  *MBS_MNI_OPCODE   = (my_bid << 20) | (my_cid << 16) | (1 << 2) | 0x2;
}


// ==========================================================================
// Send an arbitrary amount to a XUP frame (warning: counter not initialized)
// ==========================================================================
void ar_send_xup_arbitrary(int my_bid, int my_cid, unsigned int my_adr, 
                           unsigned int xup_adr, int num_pixels, int my_cnt) {

  // Sanity checks
  ar_assert(num_pixels > 0);
  ar_assert(num_pixels % 16 == 0);
  
  // Wait until DMA engine can support 1 more DMA
  while ((ar_ni_status_get(my_cid) & 0xFF) < 1) {
    ;
  }

  // DMA
  *MBS_MNI_SRC_ADR  = ar_virt_to_phys(my_bid, my_adr);
  *MBS_MNI_DST_ADR  = xup_adr;
  *MBS_MNI_ACK_ADR  = ar_addr_of_cnt(my_bid, my_cid, my_cnt);
  *MBS_MNI_DMA_SIZE = num_pixels * 4;
  *MBS_MNI_BRD_NODE = (my_bid << 20) | (my_cid << 16) |
                      (my_bid  << 4) | (my_cid << 0);
  *MBS_MNI_OPCODE   = (AR_XUP_BID << 20) | (0x0 << 16) | (1 << 2) | 0x2;
}


