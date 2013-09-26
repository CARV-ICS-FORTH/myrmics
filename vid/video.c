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
// Abstract      : XUP board video-related functionality
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: video.c,v $
// CVS revision  : $Revision: 1.11 $
// Last modified : $Date: 2013/03/21 17:36:37 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#include <arch.h>
#include <kernel_toolset.h>
#include <video.h>


// ###########################################################################
// ###                                                                     ###
// ###                            XUP board I2C master                     ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// ===========================================================================
void vid_i2c_init(int my_bid, int my_cid) {

  // Prescale registers
  ar_write_xup_register(my_cid, VID_I2C_SL, 0x31);
  ar_write_xup_register(my_cid, VID_I2C_SH, 0x00);

  // Enable
  ar_write_xup_register(my_cid, VID_I2C_CT, 0x80);
}


// ===========================================================================
// ===========================================================================
void vid_i2c_write_byte(int my_bid, int my_cid, int i2c_slave_addr, 
                        int i2c_slave_offset, int val) {

  // Slave address
  ar_write_xup_register(my_cid, VID_I2C_TX_RX, i2c_slave_addr << 1);

  // Start condition & write to slave
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x90);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }

  // Wait RX ack
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 0x80)) {
    ;
  }

  
  // Slave register offset
  ar_write_xup_register(my_cid, VID_I2C_TX_RX, i2c_slave_offset);

  // ??
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x10);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }


  // Value and EOP
  ar_write_xup_register(my_cid, VID_I2C_TX_RX, val);

  // ??
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x50);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }
}


// ===========================================================================
// ===========================================================================
unsigned int vid_i2c_read_bytes(int my_bid, int my_cid, int i2c_slave_addr, 
                                int i2c_slave_offset, int num_bytes) {

  unsigned int value;
  int i;


  // Slave address
  ar_write_xup_register(my_cid, VID_I2C_TX_RX, i2c_slave_addr << 1);
  
  // Start condition & write to slave
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x90);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }

  // Wait RX ack
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 0x80)) {
    ;
  }

  // Slave register offset
  ar_write_xup_register(my_cid, VID_I2C_TX_RX, i2c_slave_offset);

  // ??
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x10);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }


  // Slave address (read)
  ar_write_xup_register(my_cid, VID_I2C_TX_RX, (i2c_slave_addr << 1) | 0x1);
  
  // Start condition & write to slave
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x90);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }


  // For all bytes - 1
  for (i = 1; i < num_bytes; i++) {

    // End (RD, ACK)
    ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x28);

    // Transaction in progress
    while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
      ;
    }

    // Read byte
    value = ar_read_xup_register(my_bid, my_cid, VID_I2C_TX_RX);
  }

  // End (RD, NACK, STOP)
  ar_write_xup_register(my_cid, VID_I2C_CR_SR, 0x68);

  // Transaction in progress
  while ((ar_read_xup_register(my_bid, my_cid, VID_I2C_CR_SR) & 2)) {
    ;
  }

  // Read byte
  value = ar_read_xup_register(my_bid, my_cid, VID_I2C_TX_RX);

  return value;
}



// ###########################################################################
// ###                                                                     ###
// ###                        XUP board Video Out module                   ###
// ###                                                                     ###
// ###########################################################################

// I2C configuration (address, value) for the DVI out chip (Chrontel CH7301C).
// This configuration is for 1024x768 resolution.
unsigned char vid_out_config_regs[] = {
   0x21, 0x09,
   0x33, 0x08,
   0x34, 0x16,
   0x36, 0x60,
   0x49, 0xC0
};

// All I2C register addresses for dumping purposes
unsigned char vid_out_regs_addrs[] = {
  0x1C, 0x1D, 0x1E, 0x1F, 0x20, 0x21, 0x23, 0x31, 0x33, 0x34,
  0x35, 0x36, 0x37, 0x48, 0x49, 0x4A, 0x4B, 0x4D, 0x56
};


// ===========================================================================
// ===========================================================================
void vid_out_configure_i2c(int my_bid, int my_cid) {

  int i;

  for (i = 0; i < sizeof(vid_out_config_regs); i += 2) {
    vid_i2c_write_byte(my_bid, my_cid, VID_I2C_VID_OUT_SLV_ADR, 
                       vid_out_config_regs[i], vid_out_config_regs[i+1]);
  }
}


// ===========================================================================
// ===========================================================================
void vid_out_dump_i2c(int my_bid, int my_cid) {

  int          i;
  unsigned int value;

  for (i = 0; i < sizeof(vid_out_regs_addrs); i++) {
    value = vid_i2c_read_bytes(my_bid, my_cid, VID_I2C_VID_OUT_SLV_ADR,
                               vid_out_regs_addrs[i], 1);
    kt_printf("Video Out %02x = %02x\n\r", vid_out_regs_addrs[i], value);
  }
}

// ===========================================================================
// ===========================================================================
void vid_out_set_1024x768(int my_bid, int my_cid, unsigned int mask) {

  // Disable controller
  ar_write_xup_register(my_cid, VID_OUT_CONTROL, 0x0);

  // Video buffer
  ar_write_xup_register(my_cid, VID_OUT_FB_ADR, VID_DRAM_BASE);

  // Buffer mask
  ar_write_xup_register(my_cid, VID_OUT_FB_INDICES, mask);

  // HSync len
  ar_write_xup_register(my_cid, VID_OUT_HOR_SYNC, 136);

  // HSync fporch
  ar_write_xup_register(my_cid, VID_OUT_HOR_FRONT_PORCH, 24);

  // HSync bporch
  ar_write_xup_register(my_cid, VID_OUT_HOR_BACK_PORCH, 160);

  // HSync active
  ar_write_xup_register(my_cid, VID_OUT_HOR_PIXEL_TIME, 1024);

  // VSync len
  ar_write_xup_register(my_cid, VID_OUT_VER_SYNC, 6);

  // VSync fporch
  ar_write_xup_register(my_cid, VID_OUT_VER_FRONT_PORCH, 3);

  // VSync bporch
  ar_write_xup_register(my_cid, VID_OUT_VER_BACK_PORCH, 29);

  // VSync active
  ar_write_xup_register(my_cid, VID_OUT_VER_FRAME, 768);

  // Enable controller: 4 = ycbcr , 3:2 = clk_sel, 1 = soft rst, 0 = enable
  ar_write_xup_register(my_cid, VID_OUT_CONTROL, 0x09);
}


// ===========================================================================
// ===========================================================================
void vid_out_set_800x600(int my_bid, int my_cid, unsigned int mask) {

  // Disable controller
  ar_write_xup_register(my_cid, VID_OUT_CONTROL, 0x0);

  // Video buffer
  ar_write_xup_register(my_cid, VID_OUT_FB_ADR, VID_DRAM_BASE);

  // Buffer mask
  ar_write_xup_register(my_cid, VID_OUT_FB_INDICES, mask);

  // HSync len
  ar_write_xup_register(my_cid, VID_OUT_HOR_SYNC, 128);

  // HSync fporch
  ar_write_xup_register(my_cid, VID_OUT_HOR_FRONT_PORCH, 40);

  // HSync bporch
  ar_write_xup_register(my_cid, VID_OUT_HOR_BACK_PORCH, 88);

  // HSync active
  ar_write_xup_register(my_cid, VID_OUT_HOR_PIXEL_TIME, 800);

  // VSync len
  ar_write_xup_register(my_cid, VID_OUT_VER_SYNC, 4);

  // VSync fporch
  ar_write_xup_register(my_cid, VID_OUT_VER_FRONT_PORCH, 1);

  // VSync bporch
  ar_write_xup_register(my_cid, VID_OUT_VER_BACK_PORCH, 23);

  // VSync active
  ar_write_xup_register(my_cid, VID_OUT_VER_FRAME, 600);

  // Enable controller: 4 = ycbcr , 3:2 = clk_sel, 1 = soft rst, 0 = enable
  ar_write_xup_register(my_cid, VID_OUT_CONTROL, 0x05);
}


// ===========================================================================
// ===========================================================================
void vid_out_reset(int my_bid, int my_cid) {

  // Disable and reset controller
  ar_write_xup_register(my_cid, VID_OUT_CONTROL, 0x2);

  // Enable controller: 4 = ycbcr , 3:2 = clk_sel, 1 = soft rst, 0 = enable
  ar_write_xup_register(my_cid, VID_OUT_CONTROL, 0x09);
}


// ===========================================================================
// ===========================================================================
void vid_out_set_frame_buffer(int my_bid, int my_cid, void *framebuf){

  // Video buffer
  ar_write_xup_register(my_cid, VID_OUT_FB_ADR, (unsigned int) framebuf);
}


// ===========================================================================
// ===========================================================================
void vid_out_set_frame_mask(int my_bid, int my_cid, unsigned int mask) {

  // Buffer mask
  ar_write_xup_register(my_cid, VID_OUT_FB_INDICES, mask);
}


// ###########################################################################
// ###                                                                     ###
// ###                        XUP board Video In module                    ###
// ###                                                                     ###
// ###########################################################################

// I2C configuration (address, value) for the VGA in chip (Analog Devices
// AD9980). This configuration is for 1024x768 resolution.
unsigned char vid_in_config_regs_1024x768[] = {

  0x1F, 0x16, // Output_Control {[7:5]=output_mode, [4]=primary_out_en, 
              //                 [3]=secondary_out_en, [2:1]=drive_strength,
              //                 [0]=clk_inv}   
  0x20, 0x05, // Output Select 2
             
  0x05, 0x40, // Red Gain
  0x06, 0x00,
  0x07, 0x40, // Green Gain
  0x08, 0x00,
  0x09, 0x40, // Blue Gain
  0x0A, 0x00,

  // RGB
  0x0B, 0x00, // Red Offset
  0x0C, 0x80,
  0x0D, 0x00, // Green Offset
  0x0E, 0x80,
  0x0F, 0x00, // Blue Offset
  0x10, 0x80,
             
  0x1B, 0x23, // Clamp and Offset (update on every clamp)
             
  0x28, 0xBF, // Required Test Write
  0x29, 0x02, // Required Test Write
  0x2D, 0xE8, // Required Test Write
  0x2E, 0xE0, // Required Test Write
             
  // 1024x768 P60
  0x01, 0x53, // PLL_div_msb {[7:0]=pll_div[11:4]}    
  0x02, 0x00, // PLL_div_lsb {[7:4]=pll_div[3:0]}  
  0x03, 0xB0, // Clk_Gen_Ctrl {[7:6]=vco_range, [5:3]=charge_pump_current, 
              //               [2]=ext_clk_en}   
  0x04, 0xA0, // Phase_Adjust {[7:3]=phase_adjust} 
  0x12, 0x18, // Hsync Control {...,[3]= hsync_out polarity}
  0x13, 0x88, // Hsync_Duration {[7:0]= hsync_pulsewidth_duration}
  0x14, 0x18, // Vsync Control {...,[3]= vsync_out polarity}
  0x19, 0x06, // Clamp_Placement {[7:0]= clamp_placement}   
  0x1A, 0x20  // Clamp_Duration {[7:0]= clamp_duration}
};


// I2C configuration (address, value) for the VGA in chip (Analog Devices
// AD9980). This configuration is for 800x600 resolution.
unsigned char vid_in_config_regs_800x600[] = {

  0x1F, 0x16, // Output_Control {[7:5]=output_mode, [4]=primary_out_en, 
              //                 [3]=secondary_out_en, [2:1]=drive_strength,
              //                 [0]=clk_inv}   
  0x20, 0x05, // Output Select 2
             
  0x05, 0x40, // Red Gain
  0x06, 0x00,
  0x07, 0x40, // Green Gain
  0x08, 0x00,
  0x09, 0x40, // Blue Gain
  0x0A, 0x00,

  // RGB
  0x0B, 0x00, // Red Offset
  0x0C, 0x80,
  0x0D, 0x00, // Green Offset
  0x0E, 0x80,
  0x0F, 0x00, // Blue Offset
  0x10, 0x80,
             
  0x1B, 0x23, // Clamp and Offset (update on every clamp)
             
  0x28, 0xBF, // Required Test Write
  0x29, 0x02, // Required Test Write
  0x2D, 0xE8, // Required Test Write
  0x2E, 0xE0, // Required Test Write
             
  // 800x600
  0x01, 0x42, // PLL_div_msb {[7:0]=pll_div[11:4]}     
  0x02, 0x00, // PLL_div_lsb {[7:4]=pll_div[3:0]}   
  0x03, 0x68, // Clk_Gen_Ctrl {[7:6]=vco_range, [5:3]=charge_pump_current, [2]=ext_clk_en}  
  0x04, 0xA0, // Phase_Adjust {[7:3]=phase_adjust} 
  0x12, 0x18, // Hsync Control {...,[3]= hsync_out polarity}
  0x13, 0x80, // Hsync_Duration {[7:0]= hsync_pulsewidth_duration}
  0x14, 0x18, // Vsync Control {...,[3]= vsync_out polarity}
  0x19, 0x04, // Clamp_Placement {[7:0]= clamp_placement}   
  0x1A, 0x3C  // Clamp_Duration {[7:0]= clamp_duration}
};


// ===========================================================================
// mode: 0 = 800x600, 1 = 1024x768
// ===========================================================================
void vid_in_configure_i2c(int my_bid, int my_cid, int mode) {
  int i;

  if (mode) {
    for (i = 0; i < sizeof(vid_in_config_regs_1024x768); i += 2) {
      vid_i2c_write_byte(my_bid, my_cid, VID_I2C_VID_IN_SLV_ADR, 
                         vid_in_config_regs_1024x768[i], 
                         vid_in_config_regs_1024x768[i+1]);
    }
  }
  else {
    for (i = 0; i < sizeof(vid_in_config_regs_800x600); i += 2) {
      vid_i2c_write_byte(my_bid, my_cid, VID_I2C_VID_IN_SLV_ADR, 
                         vid_in_config_regs_800x600[i], 
                         vid_in_config_regs_800x600[i+1]);
    }
  }
}


// ===========================================================================
// ===========================================================================
void vid_in_dump_i2c(int my_bid, int my_cid) {
  int i;
  unsigned int value;

  for (i = 0 ; i < 47; i++) {
    value = vid_i2c_read_bytes(my_bid, my_cid, VID_I2C_VID_IN_SLV_ADR, i, 1);
    kt_printf("Video In %02x = %02x\n\r", i, value);
  }

}


// ===========================================================================
// ===========================================================================
void vid_in_set_1024x768(int my_bid, int my_cid, unsigned int mask) {

  // Disable controller
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x0);

  // Video buffer
  ar_write_xup_register(my_cid, VID_IN_FB_ADR, VID_DRAM_BASE);

  // Buffer mask
  ar_write_xup_register(my_cid, VID_IN_FB_INDICES, mask);

  // Pipe latency
  ar_write_xup_register(my_cid, VID_IN_HOR_PIPE, 6);

  // HSync len
  ar_write_xup_register(my_cid, VID_IN_HOR_SYNC, 136);

  // HSync fporch
  ar_write_xup_register(my_cid, VID_IN_HOR_FRONT_PORCH, 24);

  // HSync bporch 70Hz
  ar_write_xup_register(my_cid, VID_IN_HOR_BACK_PORCH, 144);

  // HSync active
  ar_write_xup_register(my_cid, VID_IN_HOR_PIXEL_TIME, 1024);


  // VSync len
  ar_write_xup_register(my_cid, VID_IN_VER_SYNC, 6);

  // VSync fporch
  ar_write_xup_register(my_cid, VID_IN_VER_FRONT_PORCH, 3);

  // VSync bporch
  ar_write_xup_register(my_cid, VID_IN_VER_BACK_PORCH, 29);

  // VSync active
  ar_write_xup_register(my_cid, VID_IN_VER_FRAME, 768);


  // Controller reset & enable: 
  // 4 = ycbcr , 2 = dcm rst, 1 = soft rst, 0 = enable
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x04);
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x01);
}


// ===========================================================================
// ===========================================================================
void vid_in_set_800x600(int my_bid, int my_cid, unsigned int mask) {

  // Disable controller
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x0);

  // Video buffer
  ar_write_xup_register(my_cid, VID_IN_FB_ADR, VID_DRAM_BASE);

  // Buffer mask
  ar_write_xup_register(my_cid, VID_IN_FB_INDICES, mask);

  // Pipe latency
  ar_write_xup_register(my_cid, VID_IN_HOR_PIPE, 6);

  // HSync len
  ar_write_xup_register(my_cid, VID_IN_HOR_SYNC, 128);

  // HSync fporch
  ar_write_xup_register(my_cid, VID_IN_HOR_FRONT_PORCH, 40);

  // HSync bporch 70Hz
  ar_write_xup_register(my_cid, VID_IN_HOR_BACK_PORCH, 88);

  // HSync active
  ar_write_xup_register(my_cid, VID_IN_HOR_PIXEL_TIME, 800);


  // VSync len
  ar_write_xup_register(my_cid, VID_IN_VER_SYNC, 4);

  // VSync fporch
  ar_write_xup_register(my_cid, VID_IN_VER_FRONT_PORCH, 1);

  // VSync bporch
  ar_write_xup_register(my_cid, VID_IN_VER_BACK_PORCH, 23);

  // VSync active
  ar_write_xup_register(my_cid, VID_IN_VER_FRAME, 600);


  // Controller reset & enable: 
  // 4 = ycbcr , 2 = dcm rst, 1 = soft rst, 0 = enable
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x04);
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x01);
}


// ===========================================================================
// ===========================================================================
void vid_in_reset(int my_bid, int my_cid) {

  // Controller reset & enable: 
  // 4 = ycbcr , 2 = dcm rst, 1 = soft rst, 0 = enable
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x06);
  ar_write_xup_register(my_cid, VID_IN_CONTROL, 0x01);
}


// ===========================================================================
// ===========================================================================
void vid_in_set_frame_buffer(int my_bid, int my_cid, void *framebuf) {

  // Video buffer
  ar_write_xup_register(my_cid, VID_IN_FB_ADR, (unsigned int) framebuf);
}


// ===========================================================================
// ===========================================================================
void vid_in_set_frame_mask(int my_bid, int my_cid, unsigned int mask) {

  // Buffer mask
  ar_write_xup_register(my_cid, VID_IN_FB_INDICES, mask);
}


// ===========================================================================
// ===========================================================================
int vid_in_hold_last_frame(int my_bid, int my_cid, 
                           unsigned int full_vid_in_mask) {

  unsigned int mask;
  unsigned int avail;
  unsigned int candidates;
  unsigned int test;
  unsigned int previous;
  

  // Buffer mask available
  avail = ar_read_xup_register(my_bid, my_cid, VID_IN_FB_INDICES) & 0xFF;

  // Buffer mask status
  mask = ar_read_xup_register(my_bid, my_cid, VID_IN_STATUS) & 0xFF;

  // Candidates
  candidates = avail ^ mask;

  // Find the last fully written frame
  previous = 0;
  for (test = mask >> 1; test > 0; ) {
    if (test & candidates) {
      previous = test;
      break;
    }
    else {
      test = test >> 1;
    }
  }
  if (!previous) { 
    for (test = 0x80; test > 0; ) {
      if (test & candidates) {
	previous = test;
	break;
      }
      else {
	test = test >> 1;
      }
    }
  }

  // Found anything?
  if (previous == 0) {
    return -1;
  }

  // Disable it from the available buffers for video in
  if (previous ^ mask) {
    ar_write_xup_register(my_cid, VID_IN_FB_INDICES, 
                          (full_vid_in_mask ^ previous) & 0xFF);
  }

  // Return result
  switch (previous) {
    case 0x01:	return 0;
    case 0x02:	return 1;
    case 0x04:	return 2;
    case 0x08:	return 3;
    case 0x10:	return 4;
    case 0x20:	return 5;
    case 0x40:	return 6;
    case 0x80:	return 7;
    // error
    default:	return -1;
  }
}


// ###########################################################################
// ###                                                                     ###
// ###                           Font routines                             ###
// ###                                                                     ###
// ###########################################################################

// A 8-pixel wide by 12-pixel tall font (96 glyphs, ASCII slots 32-127)
//
// vid_font_8x12_data[] size: 128 x 72 = 9216 bits = 1152 chars = 288 words
unsigned int vid_font_8x12_data[] = {
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFBFAFFF,0xEFDBCFBF,0xEFBFEFFF,
  0xFFFFFFF7,0xFFBFAFFF,0xC3ABB7BF,0xDFDFABFF,0xFFFFFFF7,0xFFBFAFD7,0xAFD7B7BF,
  0xDFDFC7FF,0xFFFFFFEF,0xFFBFFF83,0xC7F7CFFF,0xBFEFEFEF,0xFFFFFFEF,0xFFBFFFD7,
  0xEBEFABFF,0xBFEFC7EF,0xFFFFFFDF,0xFFBFFF83,0xEBEBB7FF,0xBFEFAB83,0xFF83FFDF,
  0xFFFFFFD7,0x87D5B7FF,0xDFDFEFEF,0xFFFFFFBF,0xFFBFFFFF,0xEFDBCBFF,0xDFDFFFEF,
  0xDFFFBFBF,0xFFFFFFFF,0xFFFFFFFF,0xEFBFFFFF,0xBFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xCFDFCFCF,0xB787CF87,0xCFCFFFFF,0xFFFFFFC7,
  0xB79FB7B7,0xB7BFBFF7,0xB7B7FFFF,0xF7FFBFBB,0xB7DFF7F7,0xB7BFBFF7,0xB7B7FFFF,
  0xEFFFDFFB,0xB7DFEFEF,0x878F8FEF,0xCFC7BFDF,0xDF87EFF7,0xB7DFDFF7,0xF7F7B7EF,
  0xB7F7FFFF,0xBFFFF7EF,0xB7DFBFF7,0xF7F7B7DF,0xB7F7FFFF,0xDF87EFEF,0xB7DFBFB7,
  0xF7F7B7DF,0xB7EFBFDF,0xEFFFDFFF,0xCF8F87CF,0xF78FCFDF,0xCFDFFFBF,0xF7FFBFEF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xC7EF8FC7,0x8F8787C7,0xB78FC7B7,0xBFBBBBCF,0xBBD7B7BF,
  0xB7BFBFBF,0xB7DFEFB7,0xBF939BB7,0xB3D7B7BF,0xB7BFBFBF,0xB7DFEFAF,0xBFAB9BB7,
  0xABC78FBF,0xB78F8FA7,0x87DFEF9F,0xBFBBABB7,0xABBBB7BF,0xB7BFBFB7,0xB7DFEFAF,
  0xBFBBABB7,0xB3BBB7BF,0xB7BFBFB7,0xB7DFEFAF,0xBFBBB3B7,0xBFBBB7BF,0xB7BFBFB7,
  0xB7DFEFB7,0xBFBBB3B7,0xC7BB8FC7,0x8F87BFC7,0xB78F9FB7,0x87BBBBCF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0x8FCF8FC7,0x83B7BBBB,0xBBBB878F,0xBF8FEFFF,0xB7B7B7BF,0xEFB7BBBB,
  0xD7D7F7BF,0xBFEFD7FF,0xB7B7B7BF,0xEFB7BBBB,0xD7D7EFBF,0xDFEFBBFF,0x8FB78FCF,
  0xEFB7D7AB,0xEFEFEFBF,0xDFEFFFFF,0xBFB7AFF7,0xEFB7D7AB,0xEFEFDFBF,0xEFEFFFFF,
  0xBFB7B7F7,0xEFB7D7AB,0xD7EFDFBF,0xEFEFFFFF,0xBFA7B7F7,0xEFB7EFD7,0xD7EFBFBF,
  0xF7EFFFFF,0xBFC3B78F,0xEFCFEFD7,0xBBEF87BF,0xF7EFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFF8F,0xFF8FFF87,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xBFFFBFFF,0xF7FFEFFF,0xBFFFFFBF,0xBFFFFFFF,0xDFFFBFFF,0xF7FFDFFF,0xBFBFDFBF,
  0xBFFFFFFF,0xEFFFBFFF,0xF7FFDFFF,0xBFFFFFBF,0xBFFFFFFF,0xFFCF8FC7,0xC7CF8FC7,
  0x8FBFDFB7,0xBF978FCF,0xFFB7B7BF,0xB7B7DFB7,0xB7BFDFAF,0xBFABB7B7,0xFFB7B7BF,
  0xB787DFB7,0xB7BFDF9F,0xBFABB7B7,0xFFB7B7BF,0xB7BFDFB7,0xB7BFDFAF,0xBFABB7B7,
  0xFFC78FC7,0xC7C7DFC7,0xB7BFDFB7,0xBFABB7CF,0xFFFFFFFF,0xFFFFFFF7,0xFFFFDFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFCF,0xFFFFBFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFEF,0xBFBFCBFF,0xFFFFFFFF,0xDFFFFFFF,0xFFFFFFDF,0xBFDFB7FF,
  0xFFFFFFFF,0xDFFFFFFF,0xFFFFFFDF,0xBFDFFFFF,0x8FC78FC7,0x87B7BBBB,0xBBBB87DF,
  0xBFDFFFFF,0xB7B7B7BF,0xDFB7BBAB,0xD7BBF7BF,0xBFEFFFFF,0xB7B7BFCF,0xDFB7D7AB,
  0xEFD7CFDF,0xBFDFFFFF,0xB7B7BFF7,0xDFB7D7AB,0xD7D7BFDF,0xBFDFFFFF,0x8FC7BF8F,
  0xEFC7EFD7,0xBBEF87DF,0xBFDFFFFF,0xBFF7FFFF,0xFFFFFFFF,0xFFEFFFEF,0xBFBFFFFF,
  0xBFF3FFFF,0xFFFFFFFF,0xFF9FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF
};

// A 19-pixel wide by 32-pixel tall font (96 glyphs, ASCII slots 32-127)
//
// vid_font_19x32_data[] size: 304 x 192 = 58368 bits = 7296 chars = 1824 words
unsigned int vid_font_19x32_data[] = {
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFEFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xF3FFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFC7,0xFFFFFFFF,0xFFFFFFFF,0x33FFF3FF,0xFFFFFFFF,0xFFFFFFFF,0xFBFF7FFF,
  0xFEFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFCFFFFF,0xFFE7FFC1,0x07FF33FF,0xF1FFF07F,
  0xFFFFFFF0,0x7FFFF1FE,0x3FFFFC7F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF8F,0xFFFFFFE7,
  0xFFC10FFF,0x33FFC01F,0xE03FFFFF,0xFFF07FFF,0xF3FF3FFF,0xFC7FFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFF9FFFFF,0xFFE7FFC1,0x0FFF33FF,0x841FE73F,0xFFC7FFF0,0x7FFFE3FF,
  0x1FFFFC7F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF1F,0xFFFFFFE7,0xFFE10FFF,0x33FF1F1F,
  0xCF9FFF01,0xFFF07FFF,0xE7FF9FFF,0xCC67FF9F,0xFFFFFFFF,0xFFFFFFFF,0xFF3FFFFF,
  0xFFE7FFE3,0x0FFE33FF,0x3F9FE73F,0xFE33FFF8,0xFFFFC7FF,0x8FFFC007,0xFF9FFFFF,
  0xFFFFFFFF,0xFFFFFF3F,0xFFFFFFE7,0xFFE38FFE,0x73FF3FFF,0xE23FFE7F,0xFFF8FFFF,
  0xCFFFCFFF,0xE00FFF9F,0xFFFFFFFF,0xFFFFFFFF,0xFE7FFFFF,0xFFE7FFE3,0x8FFE63FF,
  0x3FFFF079,0xFE7FFFF8,0xFFFF8FFF,0xC7FFFC7F,0xFF9FFFFF,0xFFFFFFFF,0xFFFFFE7F,
  0xFFFFFFE7,0xFFE38FF0,0x00FF0FFF,0xFDE1FE7F,0xFFF8FFFF,0x8FFFC7FF,0xF87FFF9F,
  0xFFFFFFFF,0xFFFFFFFF,0xFC7FFFFF,0xFFE7FFE3,0x9FF000FF,0x80FFFF03,0xFE7FFFF8,
  0xFFFF8FFF,0xC7FFF93F,0xFF9FFFFF,0xFFFFFFFF,0xFFFFFCFF,0xFFFFFFE7,0xFFF7DFFE,
  0x67FFE03F,0xF81FFF3F,0xFFFDFFFF,0x8FFFC7FF,0xF11FF000,0x7FFFFF80,0x01FFFFFF,
  0xF8FFFFFF,0xFFE7FFFF,0xFFFE67FF,0xFF1FE07F,0xFC1FFFFF,0xFFFF9FFF,0xE7FFF39F,
  0xF0007FFF,0xFF8001FF,0xFFFFF9FF,0xFFFFFFE7,0xFFFFFFF0,0x00FFFF9F,0xC3DFF81C,
  0x7FFFFFFF,0x9FFFE7FF,0xFFFFFF9F,0xFFFFFFFF,0xFFFFFFFF,0xF1FFFFFF,0xFFE7FFFF,
  0xFFF000FF,0xFF8FDF07,0xF98CFFFF,0xFFFF8FFF,0xC7FFFFFF,0xFF9FFFFF,0xFFFFFFFF,
  0xFFFFF3FF,0xFFFFFFFF,0xFFFFFFFE,0x67FE7F9F,0xFE23F9C1,0xFFFFFFFF,0x8FFFC7FF,
  0xFFFFFF9F,0xFFFFFFFF,0xFFFFFFFF,0xE3FFFFFF,0xFFFFFFFF,0xFFFE67FE,0x3F9FFE73,
  0xF9E1FFFF,0xFFFF8FFF,0xC7FFFFFF,0xFF9FFFC1,0xFFFFFFFF,0xFFFFE7FF,0xFFFFFFE7,
  0xFFFFFFFE,0x67FE1E1F,0xFCF9F9E1,0xFFFFFFFF,0x8FFFC7FF,0xFFFFFF9F,0xFFC1FFFF,
  0xFFFF83FF,0xE7FFFFFF,0xFFC3FFFF,0xFFFC67FE,0x003FFE73,0xF8E3FFFF,0xFFFFCFFF,
  0xCFFFFFFF,0xFF9FFF83,0xFFFFFFFF,0x83FFCFFF,0xFFFFFFC3,0xFFFFFFFC,0xE7FFE0FF,
  0xFE03FC00,0x7FFFFFFF,0xC7FF8FFF,0xFFFFFF9F,0xFF87FFFF,0xFFFF83FF,0xCFFFFFFF,
  0xFFC3FFFF,0xFFFCC7FF,0xF3FFFF07,0xFE087FFF,0xFFFFC7FF,0x8FFFFFFF,0xFFFFFF87,
  0xFFFFFFFF,0x83FF9FFF,0xFFFFFFFF,0xFFFFFFFC,0xCFFFF3FF,0xFFFFFFFF,0xFFFFFFFF,
  0xE7FF9FFF,0xFFFFFFFF,0xFF8FFFFF,0xFFFFFFFF,0x9FFFFFFF,0xFFFFFFFF,0xFFFCCFFF,
  0xF3FFFFFF,0xFFFFFFFF,0xFFFFE3FF,0x1FFFFFFF,0xFFFFFF0F,0xFFFFFFFF,0xFFFF1FFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFF3FF,0xFFFFFFFF,0xFFFFFFFF,0xF3FF3FFF,0xFFFFFFFF,
  0xFF1FFFFF,0xFFFFFFFF,0x3FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFF1FE,0x3FFFFFFF,0xFFFFFF1F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFF3FFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFE7,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFE0F,
  0xFF07FFF0,0x3FFE03FF,0xFC7FF003,0xFFE03F00,0x0FFC07FF,0xC0FFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xF803FC07,0xFFE00FF8,0x01FFF87F,0xF003FFC0,0x3F000FF8,
  0x03FF007F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xC03FF8E3,0xFC67FFC7,0xCFF8F8FF,
  0xF87FF7FF,0xFF8FFF3F,0xCFF1F1FF,0x3E3FFFFF,0xFFFFFFFF,0xE7FFFFF3,0xFFFF801F,
  0xF1F1FFE7,0xFF8FE7FB,0xFCFFF27F,0xF7FFFF1F,0xFF7FCFF3,0xF9FE7F3F,0xFFFFFFFF,
  0xFFFF87FF,0xFFF0FFFF,0x9F8FF3F9,0xFFE7FF9F,0xE7FFFCFF,0xE27FF7FF,0xFE3FFFFF,
  0xCFF3F9FE,0x7F9FFFFF,0xFFFFFFFE,0x1FFFFFF8,0x3FFF9FCF,0xF3F9FFE7,0xFFFFE7FF,
  0xFCFFE67F,0xF7FFFC7F,0xFFFF9FE3,0xF9FE7F9F,0xFFFFFFFF,0xFFF83FFF,0xFFFE1FFF,
  0xBFCFF3F9,0xFFE7FFFF,0xE7FFF8FF,0xC67FF77F,0xFCFFFFFF,0x9FF3F9FE,0x7F9FF83F,
  0xFE0FFFF0,0xFFFFFFFF,0x87FFFFCF,0xE3F8FFE7,0xFFFFC7FF,0xC1FFCE7F,0xF00FFCF7,
  0xFFFF9FF1,0xF1FE7F1F,0xF83FFE0F,0xFFC3FF80,0x007FE1FF,0xFFCFE7FC,0xFFE7FFFF,
  0x8FFF83FF,0x8E7FF087,0xFC80FFFF,0x3FF803FE,0x3E1FF83F,0xFE0FFF0F,0xFF80007F,
  0xF07FFF8F,0xE7FCFFE7,0xFFFF1FFF,0xC1FF9E7F,0xF7E3FC08,0x7FFF3FFC,0x07FF081F,
  0xF83FFE0F,0xFC1FFFFF,0xFFFFFC1F,0xFE1FE7FC,0xFFE7FFFE,0x3FFFFCFF,0x3E7FFFF3,
  0xFC3E3FFF,0x3FF843FF,0x809FFFFF,0xFFFFFC1F,0xFFFFFFFF,0xFC1FF87F,0xE3F8FFE7,
  0xFFFC7FFF,0xFE7E3E7F,0xFFF3FC7F,0x3FFE7FF1,0xF1FFF79F,0xFFFFFFFF,0xFF0FFF80,
  0x007FF07F,0xF9FFF3F9,0xFFE7FFF8,0xFFFFFE7E,0x003FFFF9,0xFCFF3FFE,0x7FF3F9FF,
  0xFF9FFFFF,0xFFFFFFC3,0xFF80007F,0xE1FFF9FF,0xF3F9FFE7,0xFFF1FFFF,0xFE7E001F,
  0xFFF3FCFF,0x3FFE7FE3,0xF9FFFF1F,0xFFFFFFFF,0xFFF0FFFF,0xFFFF87FF,0xFFFFF3F9,
  0xFFE7FFE3,0xFFFFFE7F,0xFE7FFFF3,0xFCFF3FFC,0xFFF3F9FF,0xFE3FFFFF,0xFE0FFFF8,
  0x3FFFFFFE,0x0FFFFFFF,0xF1F1FFE7,0xFFC7F7F7,0xFC7FFE7F,0xCFF3FE7F,0x3FFCFFF3,
  0xF9FFFC7F,0xF83FFE0F,0xFFFE1FFF,0xFFF83FFF,0xFFFFF8E3,0xFFE7FF8F,0xE7F1F8FF,
  0xFE7FC7C7,0xFE3E7FFC,0xFFF1F1FF,0xF8FFF83F,0xFC1FFFFF,0x87FFFFF0,0xFFFFF07F,
  0xF803FE00,0x3F0007F0,0x01FFF81F,0xE007FF00,0x7FF9FFF8,0x03FE01FF,0xF83FFC1F,
  0xFFFFE7FF,0xFFF3FFFF,0xF07FFE0F,0xFC001F00,0x07FC03FF,0xF01FF81F,0xFF81FFF9,
  0xFFFC07FE,0x03FFF83F,0xFC3FFFFF,0xFFFFFFFF,0xFFFFF07F,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFC3F,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xF87FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFF8FF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xF8FFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFF9FF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFC0F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xF807F803,0xFE001FFE,
  0x037800FF,0x0001E000,0x1FE017C0,0xE07C001F,0xF800C038,0x3C01FE0F,0xF041F81F,
  0xE03FF1F3,0xFC03FE00,0x0FF80038,0x003F0001,0xE0001F80,0x03C0E07E,0x003FF800,
  0xC0783C01,0xFE0FF060,0xF81FC01F,0xF3F3FF81,0xFF9FC7F0,0xFC3E7E1F,0xCFF9F9FF,
  0x9F8FC3F3,0xF9FFE7FF,0xFFE7F3F8,0xFF9FFF87,0xE0F8FF3F,0x8F8FE7F9,0xFF99FF9F,
  0xE7F3FE3E,0x7F9FCFF9,0xF9FF9F1F,0xE3F3F9FF,0xE7FFFFE7,0xF3F1FF9F,0xFF87E4F8,
  0x7F3F1FC7,0xE7F9FF99,0xFF9FF3E3,0xFE3E7FCF,0xCFF9F9FF,0x9E3FF3F3,0xF9FFE7FF,
  0xFFE7F3E3,0xFF9FFF83,0xE4F83F3E,0x3FE3E7C1,0xFF38FF9F,0xE7E7FF7E,0x7FCFCFBB,
  0xF9F7BE7F,0xFFF3F9FF,0xE7FFFFE7,0xF3C7FF9F,0xFF93CCF9,0x3F3E7FF3,0xE781FF3C,
  0xFF9FE7E7,0xFFFE7FCF,0xCF3FF9E7,0xFE7FFFF3,0xF9FFE7FF,0xFFE7F30F,0xFF9FFF91,
  0xCCF91F3E,0x7FF3E739,0xFF3C7F9F,0x0FE7FFFE,0x7FE7CF3F,0xF9E7FE7F,0xFFF3F9FF,
  0xE7FFFFE7,0xF21FFF9F,0xFF998CF9,0x9F3E7FF3,0xE739FE7E,0x7F800FE7,0xFFFE7FE7,
  0xC03FF807,0xFE7FFFF0,0x01FFE7FF,0xFFE7F00F,0xFF9FFF99,0x9CF9CF3C,0xFFF9E679,
  0xFE7E7F80,0x03E7FFFE,0x7FE7C03F,0xF807FE7F,0xFFF001FF,0xE7FFBFE7,0xF0C7FF9F,
  0xFF9C1CF9,0xC73CFFF9,0xE739FC00,0x3F9FF1E7,0xFFFE7FE7,0xCF3FF9E7,0xFE7C01F3,
  0xF9FFE7FF,0x3FE7F1E3,0xFF9FF79C,0x3CF9E73E,0x7FF3E719,0xFC003F9F,0xF9E7FFFE,
  0x7FCFCFBF,0xF9F7FE7E,0x03F3F9FF,0xE7FF3FE7,0xF3F1FF9F,0xE79C3CF9,0xE33E7FF3,
  0xE781FCFF,0x1F9FF9E7,0xFFFE7FCF,0xCFF9F9FF,0xFE7FF3F3,0xF9FFE7FF,0x3FE7F3F9,
  0xFF9FE79E,0x7CF9F33E,0x7FF3E7E3,0xF9FF9F9F,0xF9E3FF7E,0x7F8FCFF9,0xF9FFFE7F,
  0xF3F3F9FF,0xE7FF3FCF,0xF3F8FF9F,0xE79FFCF9,0xF13E3FE3,0xE7FFF9FF,0x9F9FF9F1,
  0xFE7E7F9F,0xCFF9F9FF,0xFF3FF3F3,0xF9FFE7FF,0x1FCFF3FC,0xFF9FE79F,0xFCF9F83F,
  0x1FC7E3FF,0xF9FF8F9F,0xF1F8FC7E,0x7E1FCFF9,0xF9FFFF0F,0xE3F3F9FF,0xE7FF8F8F,
  0xF3FCFF9F,0xE79FFCF9,0xFC3F8F8F,0xF3FFE07E,0x020003FC,0x00F8003F,0x0001E00F,
  0xFF8007C0,0xE07E003F,0xC01FC07C,0x1C000607,0xE0607C3F,0xC01FF1F3,0xE07E0200,
  0x07FE01F8,0x00FF0001,0xE007FFE0,0x0F80E03C,0x001FE03F,0xC03E1C00,0x0603E060,
  0x3E3FE03F,0xF803FFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFC07,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFF7FFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFE7F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFF3FFFFF,0xFFFFEFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFC1FF3F,0xFFF07FFF,
  0xC7FFFFFF,0x8007FF80,0xFE001FFE,0x027C0007,0x01C0C0FC,0x0C0F8181,0xF0383E0F,
  0x8007FFC1,0xFF1FFFF0,0x7FFF83FF,0xFFFF8003,0xFF007E00,0x0FFC007C,0x000703C0,
  0xC0FC0C0F,0x8181F038,0x3E0F8007,0xFFCFFF9F,0xFFFE7FFF,0x13FFFFFF,0xE7F1FE3E,
  0x3F9FC7F8,0xF87CF3E7,0xCFF9F3FF,0x3E7FF3E3,0xF8FE7F3F,0x9FCFFFCF,0xFF8FFFFE,
  0x7FFF39FF,0xFFFFE7F9,0xFC7F1F9F,0xE7F9FC7C,0xF3E7CFF9,0xF3FF3E7F,0xF3F1F1FE,
  0x3E3F9F8F,0xFFCFFFCF,0xFFFE7FFE,0x38FFFFFF,0xE7F9F8FF,0x8F9FE7F1,0xFE7DF3EF,
  0xCFF9F1FE,0x3E78F3F8,0xE3FF3E7F,0x9F1FFFCF,0xFFC7FFFE,0x7FFC7C7F,0xFFFFE7FC,
  0xF9FFCF9F,0xE3F9FEFF,0xF3FFCFF9,0xF9FE7E78,0xF3F8E3FF,0x9CFF9F3F,0xFFCFFFE7,
  0xFFFE7FF8,0xFE7FFFFF,0xE7F9F9FF,0xCF9FE7F9,0xFFFFF3FF,0xCFF9F9FE,0x7E7873FC,
  0x47FF88FF,0xFE3FFFCF,0xFFE3FFFE,0x7FF9FF3F,0xFFFFE7F9,0xF9FFCF9F,0xC7F83FFF,
  0xF3FFCFF9,0xFCFCFE30,0x73FE0FFF,0xC9FFFC7F,0xFFCFFFF3,0xFFFE7FFF,0xFFFFFFFF,
  0xE7F1F3FF,0xE7800FFE,0x03FFF3FF,0xCFF9FCFC,0xFF3267FF,0x1FFFC3FF,0xFCFFFFCF,
  0xFFF1FFFE,0x7FFFFFFF,0xFFFFE003,0xF3FFE780,0x3FFFC0FF,0xF3FFCFF9,0xFC78FF32,
  0x67FE0FFF,0xE3FFF9FF,0xFFCFFFF9,0xFFFE7FFF,0xFFFFFFFF,0xE00FF9FF,0xCF9E1FFF,
  0xFC7FF3FF,0xCFF9FE79,0xFF2227FE,0x0FFFE7FF,0xF1FFFFCF,0xFFF9FFFE,0x7FFFFFFF,
  0xFFFFE7FF,0xF9FFCF9F,0x8FFFFE7F,0xF3FFCFF9,0xFE79FF27,0x27FC47FF,0xE7FFF3F7,
  0xFFCFFFFC,0xFFFE7FFF,0xFFFFFFFF,0xE7FFF9FF,0xCF9F8FF7,0xFE7FF3FF,0xCFF1FE33,
  0xFF2727F8,0xE3FFE7FF,0xE7F3FFCF,0xFFFCFFFE,0x7FFFFFFF,0xFFFFE7FF,0xF8FF8F9F,
  0xC7F3FE7F,0xF3FFCFF3,0xFF33FF27,0x07F1F1FF,0xE7FFC7F3,0xFFCFFFFC,0x7FFE7FFF,
  0xFFFFFFFF,0xE7FFFCFF,0x9F9FE7F3,0xFE7FF3FF,0xC7F3FF23,0xFF0F87E3,0xF8FFE7FF,
  0xCFF3FFCF,0xFFFE7FFE,0x7FFFFFFF,0xFFFFE7FF,0xFE3E3F9F,0xF3F1FC7F,0xF3FFE3C7,
  0xFF07FF0F,0x87E3F8FF,0xE7FF9FF3,0xFFCFFFFE,0x3FFE7FFF,0xFFFFFFFF,0x803FFE00,
  0x7E03F0F0,0x00FF003F,0xE007FF87,0xFF0F8781,0xF03F007F,0x8003FFCF,0xFFFF3FFE,
  0x7FFFFFFF,0xFFFF803F,0xFF80FE03,0xF8F603FF,0x003FF81F,0xFF8FFF9F,0x8F81F03F,
  0x007F8003,0xFFCFFFFF,0x1FFE7FFF,0xFFFFFFFF,0xFFFFFFC7,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFCF,0xFFFF9FFE,0x7FFFFFFF,0xFFFFFFFF,
  0xFF83FFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFCFFFFF,
  0x8FFE7FFF,0xFFFFFFFF,0xFFFFFF80,0x0FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFC1,0xFFFFDFF0,0x7FFFFFFE,0x0001FFFF,0xFFC40FFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFC1FFFF,0xFFF07FFF,0xFFFE0001,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFBFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xF9FFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFF87F,0xFFFFFC1F,0xFFFFFFFF,0xFF0FFFFF,0xFFE01FFF,0xFF87FFFF,
  0xE7FFFF9F,0xE1FFFF03,0xFFFFFFFF,0xFFFFFFFF,0xFC3FFFFF,0xFE1FFFFF,0xFFFFFF0F,
  0xFFFFFF80,0x1FFFFF87,0xFFFFE7FF,0xFF9FE1FF,0xFF03FFFF,0xFFFFFFFF,0xFFFFFF3F,
  0xFFFFFF9F,0xFFFFFFFF,0xFFCFFFFF,0xFF9FFFFF,0xFFE7FFFF,0xE7FFFF9F,0xF9FFFFF3,
  0xFFFFFFFF,0xFFFFFFFF,0xFFBFFFFF,0xFF9FFFFF,0xFFFFFFCF,0xFFFFFF9F,0xFFFFFFE7,
  0xFFFFFFFF,0xFFFFF9FF,0xFFF3FFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF9F,0xFFFFFFFF,
  0xFFCFFFFF,0xFF9FFFFF,0xFFE7FFFF,0xFFFFFFFF,0xF9FFFFF3,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFF9FFFFF,0xFFFFFFCF,0xFFFFFF9F,0xFFFFFFE7,0xFFFFFFFF,0xFFFFF9FF,
  0xFFF3FFFF,0xFFFFFFFF,0xFFFFFFFF,0xFE00FF90,0x0FFE017F,0x804FF00F,0xF8007FC0,
  0x41E607FE,0x07FFC007,0xF9E07FF3,0xFE00C3F0,0x81FFC03F,0xFFFFFC00,0x7F8007FC,
  0x007F000F,0xE007F800,0x7F8001E0,0x03FE07FF,0xE007F9E0,0x7FF3FE00,0x01F0007F,
  0x801FFFFF,0xFEFE7F87,0xE3F8FC7E,0x3F0FC7E3,0xFF9FFF1F,0x0FE1F1FF,0xE7FFFFE7,
  0xF9C3FFF3,0xFF8E19FC,0x7C7F1F8F,0xFFFFFFFF,0x3F8FF1F1,0xFE7C7F8F,0x8FF1FF9F,
  0xFE3F8FE3,0xF9FFE7FF,0xFFE7F98F,0xFFF3FF9E,0x3CFCFE7E,0x3FC7FFFF,0xFFC73F9F,
  0xF9F3FF7C,0xFFCF9FF9,0xFF9FFE7F,0xCFE7F9FF,0xE7FFFFE7,0xF91FFFF3,0xFF9E7CFC,
  0xFE7E7FE7,0xFFFFFE00,0x3F9FF9F3,0xFFFCFFCF,0x8001FF9F,0xFE7FCFE7,0xF9FFE7FF,
  0xFFE7F83F,0xFFF3FF9E,0x7CFCFE7E,0x7FE7FFFF,0xFC383F9F,0xF9E7FFFC,0xFFCF8001,
  0xFF9FFE7F,0xCFE7F9FF,0xE7FFFFE7,0xF81FFFF3,0xFF9E7CFC,0xFE7E7FE7,0xFFFFF8FF,
  0x3F9FF9E3,0xFFFCFFCF,0x9FFFFF9F,0xFE7FCFE7,0xF9FFE7FF,0xFFE7F88F,0xFFF3FF9E,
  0x7CFCFE7E,0x7FE7FFFF,0xF9FF3F9F,0xF9F3FFFC,0xFFCF9FFF,0xFF9FFE7F,0xCFE7F9FF,
  0xE7FFFFE7,0xF9C7FFF3,0xFF9E7CFC,0xFE7E7FE7,0xFFFFF9FE,0x3F8FF1F1,0xFF3C7F8F,
  0x8FFFFF9F,0xFE3F8FE7,0xF9FFE7FF,0xFFE7F9E3,0xFFF3FF9E,0x7CFCFE7E,0x3FC7FFFF,
  0xF9F83F87,0xE3F8FC3E,0x3F0FC7F1,0xFF9FFF1F,0x0FE7F9FF,0xE7FFFFE7,0xF9F1FFF3,
  0xFF9E7CFC,0xFE7F1F8F,0xFFFFF800,0x0C0007FC,0x007F0001,0xE001F800,0x7F800FC1,
  0xF07C001F,0xFFE7E1F0,0x3E000E0E,0x1C707C1F,0x801FFFFF,0xFC030C10,0x0FFE01FF,
  0x8041F807,0xF8007FC0,0x4F80E038,0x001FFFE7,0xE1F03E00,0x0E061C60,0x3C1FC03F,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFCFFF,0xFFFFFFFF,0xFFE7FFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xCFFFFFFF,0xFFFFFFC7,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFCFFF,0xFFFFFFFF,0xFFCFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x8FFFFFFF,0xFFFFFFCF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFE01FFF,0xFFFFFFFF,0xE01FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFC0,0x3FFFFFFF,0xFFFFC03F,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFF1,0xFFFFFFF1,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFE3FFF3,0xFFF8FFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFBFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFC7,
  0xFFF3FFFC,0x7FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3FFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFCFFFF3,0xFFFE7FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFF3FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF9F,0xFFF3FFFE,0x7FFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x3FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFF9FFFF3,0xFFFE7FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFF3FFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF9F,0xFFF3FFFE,0x7FFFFFFF,0xFFFF0403,0xFF008383,
  0x81FC02F8,0x003F0F83,0xE0781C1F,0xC1E0E0F0,0x3E078003,0xFF9FFFF3,0xFFFE7FFF,
  0xFFFFFFFF,0x0001FE00,0x03C300F8,0x00FC003F,0x0FC3E078,0x1C1FC1E0,0xE0F83E0F,
  0x8007FF9F,0xFFF3FFFE,0x7FFFFFFF,0xFFFFE1F8,0xFC7E1FF0,0x3DF9F8FF,0x3FFFCFF3,
  0xF9FE7E7F,0xF3F9F3FC,0xFF9F9F8F,0xFF9FFFF3,0xFFFE7FFC,0x1F3FFFFF,0xE3FC78FF,
  0x1FF0FFF1,0xFCFF3FFF,0xCFF3F8FC,0x7E78E3FC,0xE7FE7F3F,0xDF0FFF3F,0xFFF3FFFF,
  0x3FF80F3F,0xFFFFE7FE,0x79FF9FF1,0xFFF8FFFF,0x3FFFCFF3,0xFCFCFF38,0xE7FE4FFE,
  0x7F3FFF1F,0xFE3FFFF3,0xFFFF1FF9,0xC03FFFFF,0xE7FE79FF,0x9FF3FFFC,0x07FF3FFF,
  0xCFF3FC7C,0xFF3067FE,0x1FFF3E7F,0xFE3FFC3F,0xFFF3FFFF,0x0FFBF07F,0xFFFFE7FE,
  0x79FF9FF3,0xFFFE01FF,0x3FFFCFF3,0xFE79FF32,0x67FF1FFF,0x3E7FFC7F,0xFE3FFFF3,
  0xFFFF1FFF,0xFFFFFFFF,0xE7FE79FF,0x9FF3FFFF,0xF8FF3FFF,0xCFF3FE79,0xFF3267FE,
  0x0FFF9CFF,0xF8FFFF3F,0xFFF3FFFF,0x3FFFFFFF,0xFFFFE7FE,0x79FF9FF3,0xFFF7FCFF,
  0x3FFFCFF3,0xFF33FF82,0x0FFC47FF,0x9CFFF1FF,0xFF9FFFF3,0xFFFE7FFF,0xFFFFFFFF,
  0xE3FC78FF,0x1FF3FFF3,0xFCFF1FFF,0xC7F3FF33,0xFF870FF8,0xE3FFC9FF,0xE3F7FF9F,
  0xFFF3FFFE,0x7FFFFFFF,0xFFFFE1F8,0xFC7E1FF3,0xFFF1FCFF,0x9F8FE7C3,0xFF03FF87,
  0x0FF1F1FF,0xC9FFC7F3,0xFF9FFFF3,0xFFFE7FFF,0xFFFFFFFF,0xE001FE00,0x1F800FF0,
  0x00FF800F,0xE000FF87,0xFF878FC0,0xE07FC3FF,0xC003FF9F,0xFFF3FFFE,0x7FFFFFFF,
  0xFFFFE403,0xFF009F80,0x07F803FF,0xC03FF010,0xFF87FF8F,0x8FC0E07F,0xE3FF8003,
  0xFF9FFFF3,0xFFFE7FFF,0xFFFFFFFF,0xE7FFFFFF,0x9FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFE7FF,0xFFFFFF9F,0xFFF3FFFE,0x7FFFFFFF,0xFFFFE7FF,0xFFFF9FFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xE7FFFFFF,0xFFCFFFF3,0xFFFE7FFF,
  0xFFFFFFFF,0xE7FFFFFF,0x9FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFC7FF,
  0xFFFFFFC7,0xFFF3FFF8,0xFFFFFFFF,0xFFFFE7FF,0xFFFF9FFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xCFFFFFFF,0xFFE3FFF3,0xFFF1FFFF,0xFFFFFFFF,0x80FFFFFC,
  0x07FFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFC07FF,0xFFFFFFFF,0xFFFBFFFF,
  0xFFFFFFFF,0xFFFF007F,0xFFF803FF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFF8,
  0x03FFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF
};


// ===========================================================================
// vid_font_putchar()           Draw a font glyph onto a frame in memory
// ===========================================================================
// * INPUTS
//   int buf_width              Buffer width (i.e. number of pixels/line)
//   int buf_x                  Left origin of frame to draw the font
//   int buf_y                  Top origin of frame to draw the font
//   int font_large             1: Use the large font (19x32)
//                              0: Use the small font (8x12)
//   unsigned int fg_color      Foreground font color in BBGGRRxx format
//   unsigned int bg_color      Background color in BBGGRRxx format
//   int transparency           0: no transparency (use fg/bg colors normally)
//                              1-256: don't use bg color, blend the font on
//                                     the existing frame using the 
//                                     transparency value: 1 is almost fully 
//                                     transparent, 256 is fully opaque.
//   int ch                     The character to draw
//
// * INOUTS   
//   unsigned int *buf          Frame in memory, in BBGGRRxx format (1 pixel
//                              per integer)
// ===========================================================================
void vid_font_putchar(unsigned int *buf, int buf_width, int buf_x, int buf_y, 
                      int font_large, unsigned int fg_color, 
                      unsigned int bg_color, int transparency, int ch) {

  int   font_x;
  int   font_y;
  int   font_index;
  int   buf_index;
  int   x;
  int   y;
  unsigned int b;
  unsigned int g;
  unsigned int r;
  

  // Get font origin
  if (font_large) {
    font_x = (ch - 32) % 16 * 19;
    font_y = (ch - 32) / 16 * 32;
  }
  else {
    font_x = (ch - 32) % 16 * 8;
    font_y = (ch - 32) / 16 * 12;
  }

  // Blit
  for (x = 0; x < (font_large ? 19 : 8); x++) {
    for (y = 0; y < (font_large ? 32 : 12); y++) {

      // Compute indices
      font_index = (font_y + y) * (font_large? 304 : 128) + (font_x + x);
      buf_index = (buf_y + y) * buf_width + (buf_x + x);

      // Font bit inactive?
      if (((font_large ? vid_font_19x32_data[font_index / 32] :
                         vid_font_8x12_data[font_index / 32]) & 
           (1 << (31 - (font_index % 32))))) {
        
        // Draw background pixel
        if (!transparency) {
          buf[buf_index] = bg_color;
        }
      }

      // Font bit active?
      else {
        // Draw fully opaque...
        if (!transparency) {
          buf[buf_index] = fg_color;
        }
        // ... or blit semi-transparent value
        else {

          // Compute new pixel values
          b = (((fg_color >> 24) * transparency) + 
               ((buf[buf_index] >> 24) * (256 - transparency))) / 256;
          g = ((((fg_color >> 16) & 0xFF) * transparency) + 
               (((buf[buf_index] >> 16) & 0xFF) * (256 - transparency))) / 256;
          r = ((((fg_color >> 8) & 0xFF) * transparency) + 
               (((buf[buf_index] >> 8) & 0xFF) * (256 - transparency))) / 256;
          
          // Clamp
          if (b > 255) {
            b = 255;
          }
          if (g > 255) {
            g = 255;
          }
          if (r > 255) {
            r = 255;
          }

          // Replace
          buf[buf_index] = (b << 24) | (g << 16) | (r << 8);
        }
      }
    }
  }
}


// ===========================================================================
// condense: 0 for normal font char distance, >0 to shorten char distance
//           or <0 to expand char distance
//
// See vid_font_putchar() above for other parameters explanation
// ===========================================================================
void vid_font_putstr(unsigned int *buf, int buf_width, int buf_x, int buf_y, 
                     int font_large, int condense, unsigned int fg_color, 
                     unsigned int bg_color, int transparency, char *str) {

  while (*str) {
    if (buf_x > buf_width - (font_large ? 19 : 8)) {
      return;
    }
    vid_font_putchar(buf, buf_width, buf_x, buf_y, font_large, 
                     fg_color, bg_color, transparency, *str++);
    if (font_large) {
      buf_x += 19 - condense;
    }
    else {
      buf_x += 8 - condense;
    }
  }
}


// ###########################################################################
// ###                                                                     ###
// ###                           Blit routines                             ###
// ###                                                                     ###
// ###########################################################################

// The Myrmics logo, 108x76 pixels. The bitmap contains 2-bit color, where
// color 0 = fg, color3 = bg and others are in between.
//
// vid_blit_myrmics_data[] size: 2x108x76 = 16416 bits = 2052 chars = 513 words
unsigned int vid_blit_myrmics_data[] = {
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFE56FFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFE,0x0007FFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFF814,0x01FFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFE07D00,0xBFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x82FF90BF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF47,0xFFFBFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFE0FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFF,0xFE2FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFFFD,0x3FFFFFFF,0xFFFFFFFF,0x9556FFFF,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFFFFD7F,0xFFFFFFFF,0xFFFFA401,0x502FFFFF,0xFFFFFFFF,0xFFFFFFFB,
  0xFFFFAAFF,0xFFFCBFFF,0xFFFFFFFF,0xE501AFF4,0x0BFFFFFF,0xFFFFFFFF,0xFFFF91BF,
  0xFD006FFF,0xF8BFFFFF,0xFFFFFD00,0x6BFFF807,0xFFFFFFFF,0xFFFFFFFF,0xFD001FF4,
  0x0007FFF8,0xBFFFFFFF,0xFF900AFF,0xFFF80BFF,0xFFFFFFFF,0xFFFFFFF4,0x000BE000,
  0x01FFF8BF,0xFFFFFFFD,0x01BFFFFF,0xFD1FFFFF,0xFFFFFFFF,0xFFFFD000,0x03C00000,
  0x7FF8BFFF,0xFFFFE01F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFF800001,0x4000002F,
  0xF8BFFFFF,0xFF41FFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFEAFFFF,0x40000000,0x00000BF8,
  0xBFFFFFFD,0x1FFFFFFF,0xFFFFFFFF,0xFFFFFA50,0x0006A900,0x00000000,0x0007F8BF,
  0xFFFFF0BF,0xFFFFFFFF,0xFFFFFFFF,0xFE400000,0x00000000,0x00000000,0x01B8BFFF,
  0xFF82FFFF,0xFFFFFFFF,0xFFFFFFE0,0x00000000,0x00000000,0x00000000,0x142BFFFE,
  0x07FFFFFF,0xFFFFFFFF,0xFFFF4000,0x00000000,0x00000000,0x00000000,0x01A5500F,
  0xFFFFFFFF,0xFFFFFFFF,0xFD000000,0x00000000,0x00000000,0x00000000,0x00001FFF,
  0xFFFFFFFF,0xFFFFFFF4,0x00000000,0x00000000,0x00000000,0x00000000,0x000BAFFF,
  0xFFFFFFFF,0xFFFFE000,0x00000000,0x00000000,0x00000000,0x00005900,0x015FFFFF,
  0xFFFFFFFF,0xFF800000,0x00000000,0x00000000,0x00000000,0x05BE0000,0x6BFFFFFF,
  0xFFFFFFFF,0x40000000,0x00000000,0x001A4000,0x00000015,0x5F800067,0xFFFFFFFF,
  0xFFFFFE00,0x00000000,0x00000000,0x1FE00000,0x00005007,0xE0007AFF,0xFFFFFFFF,
  0xFFFE0000,0x00000000,0x0000000F,0xF9000000,0x004002F8,0x007AFFFF,0xFFFFFFFF,
  0xFC000000,0x00000000,0x00001FFF,0x90000001,0x4001FE00,0x7EFFFFFF,0xFFFFFFF8,
  0x00000000,0x00000000,0x002FFFFD,0x00000140,0x01FF807D,0xFFFFFFFF,0xFFFFF400,
  0x00000000,0x00000000,0x3FFFF901,0xB8014001,0xFFE07EFF,0xFFFFFFFF,0xFFF00000,
  0x00000000,0x0000002F,0xFF9007F8,0x014001FF,0xF07EFFFF,0xFFFFFFFF,0xE0000000,
  0x00000000,0x00041FFE,0x002FF801,0x4001FFF4,0x7EBFFFFF,0xFFFFFFE0,0x00000000,
  0x00000000,0x090FF800,0xBFF98140,0x01FFF47E,0xBFFFFFFF,0xFFFFD000,0x00000000,
  0x0000024D,0x0BF802FF,0xF6C04002,0xFFF83EBF,0xFFFFFFFF,0xFFD00000,0x00000000,
  0x001F4E02,0xA40BFFF6,0xC05002FF,0xF82EBFFF,0xFFFFFFFF,0xC0000000,0x00000000,
  0x2F4F4000,0x2FFFF5A4,0x10007FF8,0x1EBFFFFF,0xFFFFFFC0,0x00000000,0x0000101F,
  0x4FD0001F,0xFFF0B800,0x001FF41A,0xBFFFFFFF,0xFFFF8000,0x00000000,0x01B90B8B,
  0xF90002FF,0xF5F90000,0x07F001FF,0xFFFFFFFF,0xFF800000,0x0000002B,0xFF478BFF,
  0xA4002FEB,0xAF400001,0x9002FFFF,0xFFFFFFFF,0x80000000,0x0001FFFF,0xD2C3FFFD,
  0x0007EB5F,0x50000000,0x02FFFFFF,0xFFFFFF80,0x00000000,0x02FFFFF5,0xD1FFFD0A,
  0x41E6EAB4,0x00000002,0xFFFFFFFF,0xFFFF4000,0x00000001,0xFFFFF9E0,0xBFFE0690,
  0xF2F5FA40,0x000002FF,0xFFFFFFFF,0xFF400000,0x000050BF,0xFFF9F02F,0xFF400075,
  0x9A9B8000,0x0002FFFF,0xFFFFFFFF,0x40000000,0x00B42FFF,0xF9F81FFF,0xD0002A45,
  0x06400000,0x02FBFFFF,0xFFFFFF80,0x00000001,0xF81FFFF5,0xFD0BFFF9,0x0006BEAA,
  0xA9000000,0x007FFFFF,0xFFFFC000,0x000005FD,0x0BFFE2FF,0x43FFFFE4,0x006FFFFF,
  0x80000040,0x0BFFFFFF,0xFFE00000,0x0026FF42,0xFFD7FFD1,0xFFFFFC14,0x02FFFFE0,
  0x0002F902,0xFFFFFFFF,0xF4000000,0xA7FFD1FF,0x8BFFF5FF,0xFFFC7E01,0xFFFFF800,
  0x02FFD07F,0xFFFFFFFD,0x00000287,0xFFF4FF0F,0xFFFDFFFF,0xFC7FD0BF,0xFFFD1007,
  0xFFFD1FFF,0xFFFFFF90,0x000B4FFF,0xFCBE2FFF,0xFDFFFFF8,0x7FF4BFFF,0xFE295BFF,
  0xFF87FFFF,0xFFFFF941,0x7E1FFFFD,0xB83FFFF9,0xFFFFF87F,0xFC7FFFFF,0x1FFFFFFF,
  0xE6FFFFFF,0xFFFFEBFC,0x3FFFF8F4,0xBFFFF9FF,0xFFF87FFD,0x7FFFFF4F,0xFFFFFFFE,
  0xFFFFFFFF,0xFFFFF47F,0xFFF9E0FF,0xFFF5FFFF,0xF87FFD7F,0xFFFF8BFF,0xFFFFFFFF,
  0xFFFFFFFF,0xFFE0FFFF,0xF6D1FFFF,0xE2FFFFF4,0x7FFD7FFF,0xFFD3FFFF,0xFFFFFFFF,
  0xFFFFFFFF,0xD1FFFFE2,0x82FFFFC7,0xFFFFF47F,0xFD7FFFFF,0xE2FFFFFF,0xFFFFFFFF,
  0xFFFFFF82,0xFFFFD747,0xFFFF4BFF,0xFFF47FFD,0x7FFFFFF4,0xFFFFFFFF,0xFFFFFFFF,
  0xFFFF43FF,0xFF8B07FF,0xFE1FFFFF,0xF07FFC3F,0xFFFFF8BF,0xFFFFFFFF,0xFFFFFFFF,
  0xFF07FFFE,0x1E07FFFC,0x2FFFFFF0,0x7FFC7FFF,0xFFFE7FFF,0xFFFFFFFF,0xFFFFFFFE,
  0x07FFFD2E,0x03FFF47F,0xFFFFE07F,0xF87FFFFF,0xFFAFFFFF,0xFFFFFFFF,0xFFFFFE02,
  0xFFF87F02,0xFFE0BFFF,0xFFE07FF8,0x7FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFF41FF,
  0xE0BF80BF,0xD1FFFFFF,0xE07FF87F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFF807FD1,
  0xFFD01F82,0xFFFFFFE0,0xBFF47FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xD01F82FF,
  0xF90242FF,0xFFFFD0BF,0xF47FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFF9,0x0242FFFF,
  0x9003FFFF,0xFFD0BFF4,0x7FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFF90,0x02FFFFF9,
  0x03FFFFFF,0xD07FF47F,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFF902,0xFFFFFF02,
  0xFFFFFFD0,0x1FF47FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFF02FF,0xFFFF41FF,
  0xFFFFD007,0xE07FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFF41FFFF,0xFF807FFF,
  0xFFE401E0,0x7FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0x807FFFFF,0xD01FFFFF,
  0xFEAAE0BF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFE0,0x1BFFFFF8,0x02FFFFFF,
  0xFFE0BFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFF902,0xFFFFFF90,0xBFFFFFFF,
  0xE0BFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFF907F,0xFFFFFA6F,0xFFFFFFE0,
  0x7FFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFA6FFF,0xFFFFFFFF,0xFFFFD01F,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFD007FF,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xE001BFFF,
  0xFFFFFFFF,0xFFFFFFFF
};

// The FORTH logo, 72x76 pixels. The bitmap contains 2-bit color, where
// color 0 = fg, color3 = bg and others are in between.
//
// vid_blit_forth_data[] size: 2x72x76 = 10944 bits = 1368 chars = 342 words
unsigned int vid_blit_forth_data[] = {
  0xFFFFFFFF,0xFFFFFFEA,0xAAAAAFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xE5FFFFFE,
  0xFD5FFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFD7FFF,0xFBFFFFF9,0x7FFFFFFF,0xFFFFFFFF,
  0xFFFFFFD7,0xF9F1FBF6,0xD4FFD7FF,0xFFFFFFFF,0xFFFFFFFF,0xFD7FF9C5,0x9FFBAEFF,
  0xF5BFFFFF,0xFFFFFFFF,0xFFFFE7FF,0xFCC1AF93,0xFBD0FBDF,0xFFFFFFFF,0xFFFFFFFF,
  0x6FFE2CF6,0x7F971CBF,0xCFF7FFFF,0xFFFFFFFF,0xFFFDFFFC,0xDCFFF7AB,0x807B4DFE,
  0xBFFFFFFF,0xFFFFFFEB,0xFD1EDFBF,0xFB47FFEA,0x2CFFEFFF,0xFFFFFFFF,0xFFAFFAFB,
  0xEBFF6F4F,0xFFFE3E13,0xEBFFFFFF,0xFFFFFEBF,0xEEAB9BFF,0x0F5E9BFD,0x6AAFFEFF,
  0xFFFFFFFF,0xF9FF9EAF,0xFFFFEF6E,0xAE6FF657,0xFEBFFFFF,0xFFFFFABE,0xEEFFFAA9,
  0xFE3EEFB6,0xF6FA5F6F,0xFFFFFFFF,0xEFEF7C7F,0xBFFFAE77,0xBEF7AAFA,0xFAEFFFFF,
  0xFFFFBFEB,0xC0FAFBFF,0xEBE06DEF,0xFBF7EFFB,0xFFFFFFFE,0xFEAAD7DF,0x3796FF79,
  0x1AEF64BB,0xAFBEFFFF,0xFFFAFEFF,0xBF7F3B7E,0xEEFFFECE,0xFDEFB92F,0xBFFFFFF7,
  0xCAFBEBB7,0x776BDEFF,0xFE9CFFE9,0xF997BFFF,0xFFEFDAAF,0xDBE801EF,0xCF7BFFFB,
  0xFE6DF593,0xEFFFFFDF,0xFFBFBF5C,0xFBFFCBB7,0xBFE3B0F9,0xE25BEFFF,0xFFBF6ABE,
  0xBEAABFBF,0xFBDBEAF0,0x83E7E5AF,0xFBFFFFBA,0x9AFEE7EA,0x6AFAFF6E,0xADBE8F5B,
  0xFBFBABFF,0xFDF9EBFA,0xEECDBFF7,0x9EFD3F6F,0xDD2EFFEA,0x6DFFFDFE,0xA7FBBCBF,
  0xFF7BAEAF,0x1BB7FBEA,0x7EEAFEFF,0xFBF89BBA,0xFBFFFF3F,0x7FEB85FB,0xFEEFEEFF,
  0xFF7FF2B7,0xDBBDFEEE,0x7F1F3FFF,0xBEB7EEEE,0xFFBEAF7F,0xFBFBEFFE,0xE77F7F0E,
  0x3FBFEFF3,0x676EEFBE,0xDFEFDFBB,0xEF6EAEFF,0x9A5C3AAB,0xF7F6FBB9,0xBBBBDFCF,
  0xDFD03AD6,0xFEBF9EFC,0x72F9BCEB,0xA6EDFDFF,0x99DFDFF1,0xFFFEAB5B,0xEFE9FB94,
  0x7EF5FFEE,0x7AFE7FFB,0xFEFAF83F,0xF3F49AB6,0xEAF5FF7F,0xF8FB74AB,0xFAF7B9AA,
  0xFFE01E7D,0xB7FBABF7,0x43BC63F6,0x9EBFEEF7,0xBF5BF2FD,0x1FFAEFF5,0xFBFD1EDD,
  0xFB7BFBF7,0xD4FBBFF9,0xF3FFF28F,0x0FFC2E6F,0x3CDEBFF6,0xAABBFFFF,0xBDBAF3FE,
  0x7D66F803,0xFFF5FEDF,0xE9FBFFFB,0xEBFEBD4A,0xEAAF3F25,0x2FD7D79F,0xE6DFBCBF,
  0x49BF688A,0xBFFFEEE6,0x3CB3FC7F,0xD3AFBFED,0x91FF2FFD,0x6E5EBD6B,0xEFEFBDFD,
  0xEFA98265,0x2F2FAFFB,0x0E3BFABE,0xBDF7EFF7,0xE800D19F,0xA3FFBA8F,0x7BFBCFFB,
  0xDF7EBDF7,0xEC5A6ABD,0xFCF2FFBF,0xBFCCB5F7,0xE1FBDFBE,0xBFABFFFF,0x7E5EFCF1,
  0xBEFC2DAE,0xFCF7EEBB,0xF01ABFF0,0xFEA7FAFB,0x7D4FEBF7,0xFFEEBFFB,0xDEF9BFFE,
  0xBCACBFF5,0xFBF7BD29,0x0FE7E379,0xD9EFCAFB,0xF93EBFBC,0xBBF0EE60,0xE6EF97FE,
  0xBB7AFEDD,0xAFEFBBBE,0xBF08FE1D,0xAFF837BA,0x784396FD,0x55EEFBEF,0x4EFEBFFF,
  0xEBFCAF5F,0x3ED6FDD3,0xBEFFFEBF,0x83EEDFFE,0xEFAABF8A,0xF7B7BBA0,0xF9F3FB6B,
  0xFFFE8BBF,0x6CFECAF5,0xBFFFFBE7,0xEFEEF9F7,0xDCB1FCE7,0xFFBFCCFE,0xDF80BDAF,
  0xDDFF95FB,0xFDFFBC33,0x2FEDFEBF,0xF9A3DFFE,0xFEF43EFB,0xBFBB79F9,0xFF9BCF49,
  0x2EAB82F7,0xFBFBAF9F,0xFFAFAA9F,0xF9AF5AB5,0x5C2B9BFF,0x09F7F75E,0xBFBF9FDF,
  0xE96FFFFA,0x8EBE3CB9,0xEBFD3DFB,0xFBFEFFBF,0x2DFBFEBA,0x5EFBDDFE,0xFCAEEFFA,
  0xAA9FFDFF,0xD7EF2EFE,0xFCFFEE67,0xDDF7FFF7,0xEEAFBADF,0xFDFFFBEB,0x78AFBEF7,
  0xB7E3EFEF,0xEBAFBE9B,0xBBBFFFB9,0x9FEBEEF7,0xFFF92BFB,0xFEBFABBE,0xFFA67F7F,
  0xFF7E2FBE,0xFDF87EBF,0xFEF7EBFA,0xABBBFF29,0x7EFFFFDF,0xFEAFBEAC,0x0F3FFFF6,
  0xBF73FAEB,0xF6BDBDFF,0xFFEFEAAF,0xAFFEC7BE,0xAEFFF9D7,0xBAEFFBBF,0xFBFFFFF6,
  0xFEFFEBEB,0xDB6DEEFF,0xBDDFEE7F,0xAF4BFBFF,0xFFFDF5A9,0xF9FB3F79,0xF8230EDF,
  0xDBEF8FEF,0xEFFFFFFE,0xFBFEFFBF,0xEDB6E9B3,0x4DD7DFFF,0xDBFFBFFF,0xFFFFBFFA,
  0xC3EBBDFB,0xFEB7CDBF,0x7FF7BD2E,0xBFFFFFFF,0xAF9B44FE,0x3EF2EFB2,0xCFFAFF27,
  0xEF7AFFFF,0xFFFFEBFE,0xBDFFE7F1,0xE036FF9F,0x8B7EFBEB,0xFFFFFFFF,0xFAF9D25B,
  0xFEAFFEBE,0xEAFFCB2F,0xAE6FFFFF,0xFFFFFEFA,0xA9F7FFEE,0x5AA9AFFE,0xF2FBBEBF,
  0xFFFFFFFF,0xFF7CCEF7,0x07EFFFFF,0xEFFCA9F9,0xFAFFFFFF,0xFFFFFFEE,0x6E5B10EE,
  0xFFFDEFAE,0xF2FFDFFF,0xFFFFFFFF,0xFFFAFFFF,0xBDC669E7,0x4EBADD1E,0xBFFFFFFF,
  0xFFFFFFFF,0x7FFFBAEA,0xFFFFECF7,0x7FF7FFFF,0xFFFFFFFF,0xFFFFE7FF,0x16F6BF7F,
  0xFDF7FF6F,0xFFFFFFFF,0xFFFFFFFF,0xFD7FFFE0,0xF7FFBB9F,0xF6FFFFFF,0xFFFFFFFF,
  0xFFFFFFD7,0xFFFFF6E1,0x37FD6FFF,0xFFFFFFFF,0xFFFFFFFF,0xFFFE5FFF,0xFFFFFD5B,
  0xFFFFFFFF,0xFFFFFFFF,0xFFFFFFFF,0xFEABFFFE,0xAFFFFFFF,0xFFFFFFFF
};


// ===========================================================================
// vid_blit_logo()              Draw the Myrmics or FORTH logo onto a frame 
//                              in memory
// ===========================================================================
// * INPUTS
//   int buf_width              Buffer width (i.e. number of pixels/line)
//   int buf_x                  Left origin of frame to draw the logo
//   int buf_y                  Top origin of frame to draw the logo
//   int which_logo             0: Myrmics logo
//                              1: FORTH logo
//   unsigned int colors        Array of 4 colors, where colors[0] is the
//                              foreground, colors[3] the background and
//                              colors[1] and colors[2] the in-between shades
// * INOUTS   
//   unsigned int *buf          Frame in memory, in BBGGRRxx format (1 pixel
//                              per integer)
// ===========================================================================
void vid_blit_logo(unsigned int *buf, int buf_width, int buf_x, int buf_y, 
                   int which_logo, unsigned int *colors) {

  int           logo_idx;
  int           x;
  int           y;
  unsigned int  palette;
  int           max_x;
  int           max_y;
  unsigned int  *logo_ptr;
  
  switch (which_logo) {
    case 0: max_x = 108; max_y = 76; logo_ptr = vid_blit_myrmics_data; break;
    case 1: max_x =  72; max_y = 76; logo_ptr = vid_blit_forth_data; break;
    default: ar_abort();
  }


  for (x = 0; x < max_x; x++) {
    for (y = 0; y < max_y; y++) {

      logo_idx = y * max_x + x;

      palette = (logo_ptr[logo_idx / 16] >> 
                 (30 - 2 * (logo_idx % 16))) & 3;

      buf[(buf_y + y) * buf_width + (buf_x + x)] = colors[palette];
    }
  }
}



// ###########################################################################
// ###                                                                     ###
// ###                         Drawing routines                            ###
// ###                                                                     ###
// ###########################################################################

// ===========================================================================
// vid_draw_line()              Draw a line onto a frame in memory using
//                              Bresenham's algorithm.
//
//                              The main functionality is taken from the
//                              Wikipedia article:
//
//                       http://en.wikipedia.org/wiki/Bresenham%27s_algorithm
// ===========================================================================
// * INPUTS
//   int buf_width              Buffer width (i.e. number of pixels/line)
//   int buf_height             Buffer width (i.e. number of lines)
//   int x1                     First line point, X coordinate
//   int y1                     First line point, Y coordinate
//   int x2                     Second line point, X coordinate
//   int y2                     Second line point, Y coordinate
//   unsigned int color         Color in BBGGRRxx format
//   int transparency           0: no transparency (use color normally)
//                              1-256: blend the line on the existing frame
//                                     buffer using the transparency value: 1
//                                     is almost fully transparent, 255 is
//                                     almost fully opaque.
//
// * INOUTS   
//   unsigned int *buf          Frame in memory, in BBGGRRxx format (1 pixel
//                              per integer)
// ===========================================================================
void vid_draw_line(unsigned int *buf, int buf_width, int buf_height, 
                   int x1, int y1, int x2, int y2, unsigned int color,
                   int transparency) {

  int dx;
  int dy;
  int x;
  int y;
  int xstep;
  int ystep;
  int err;
  int err2;
  unsigned int b;
  unsigned int g;
  unsigned int r;
  int buf_index;


  // Sanity check
  ar_assert((x1 >= 0) && (x1 < buf_width) && (y1 >= 0) && (y1 < buf_height) &&
            (x2 >= 0) && (x2 < buf_width) && (y2 >= 0) && (y2 < buf_height));

  // Init
  if (x1 < x2) {
    xstep = 1;
    dx = x2 - x1;
  }
  else {
    xstep = -1;
    dx = x1 - x2;
  }
  if (y1 < y2) {
    ystep = 1;
    dy = y2 - y1;
  }
  else {
    ystep = -1;
    dy = y1 - y2;
  }
  err = dx - dy;
  x = x1;
  y = y1;

  // Loop
  while (1) {
    buf_index = y * buf_width + x;

    // Draw fully opaque...
    if (!transparency) {
      buf[buf_index] = color;
    }
    // ... or blit semi-transparent value
    else {

      // Compute new pixel values
      b = (((color >> 24) * transparency) + 
           ((buf[buf_index] >> 24) * (256 - transparency))) / 256;
      g = ((((color >> 16) & 0xFF) * transparency) + 
           (((buf[buf_index] >> 16) & 0xFF) * (256 - transparency))) / 256;
      r = ((((color >> 8) & 0xFF) * transparency) + 
           (((buf[buf_index] >> 8) & 0xFF) * (256 - transparency))) / 256;
      
      // Clamp
      if (b > 255) {
        b = 255;
      }
      if (g > 255) {
        g = 255;
      }
      if (r > 255) {
        r = 255;
      }

      // Replace
      buf[buf_index] = (b << 24) | (g << 16) | (r << 8);
    }

    // End?
    if ((x == x2) && (y == y2)) {
      break;
    }

    // Advance
    err2 = 2 * err;
    if (err2 > -dy) {
      err -= dy;
      x += xstep;
    }
    if (err2 < dx) {
      err += dx;
      y += ystep;
    }
  }
}


// ===========================================================================
// vid_draw_circle()            Draw a circle onto a frame in memory using
//                              the midpoint circle algorithm.
//
//                              The main functionality is taken from the
//                              Wikipedia article:
//
//                     http://en.wikipedia.org/wiki/Midpoint_circle_algorithm
// ===========================================================================
// * INPUTS
//   int buf_width              Buffer width (i.e. number of pixels/line)
//   int buf_height             Buffer width (i.e. number of lines)
//   int x0                     Center, X coordinate
//   int y0                     Center, Y coordinate
//   int radius                 Radius
//   unsigned int color         Color in BBGGRRxx format
//
// * INOUTS   
//   unsigned int *buf          Frame in memory, in BBGGRRxx format (1 pixel
//                              per integer)
// ===========================================================================
void vid_draw_circle(unsigned int *buf, int buf_width, int buf_height, 
                     int x0, int y0, int radius, unsigned int color) {

  int x;
  int y;
  int xChange;
  int yChange;
  int radiusError;

  // Sanity check
  ar_assert((x0 - radius >= 0) && (x0 + radius < buf_width) && 
            (y0 - radius >= 0) && (y0 + radius < buf_height));

  // Init
  x = radius;
  y = 0;
  xChange = 1 - (radius << 1);
  yChange = 0;
  radiusError = 0;
 
  // Loop
  while (x >= y) {

    buf[(y + y0 ) * buf_width + (x + x0)]  = color;
    buf[(x + y0 ) * buf_width + (y + x0)]  = color;
    buf[(y + y0 ) * buf_width + (-x + x0)] = color;
    buf[(x + y0 ) * buf_width + (-y + x0)] = color;
    buf[(-y + y0) * buf_width + (-x + x0)] = color;
    buf[(-x + y0) * buf_width + (-y + x0)] = color;
    buf[(-y + y0) * buf_width + (x + x0)]  = color;
    buf[(-x + y0) * buf_width + (y + x0)]  = color;
 
    y++;
    radiusError += yChange;
    yChange += 2;
    if (((radiusError << 1) + xChange) > 0) {
      x--;
      radiusError += xChange;
      xChange += 2;
    }
  }
}


// ===========================================================================
// vid_draw_floodfill()         Fills an area, starting from a given point
//                              and expanding until it reaches a boundary
//                              color
// ===========================================================================
// * INPUTS
//   int buf_width              Buffer width (i.e. number of pixels/line)
//   int buf_height             Buffer width (i.e. number of lines)
//   int x                      Starting point, X coordinate
//   int y                      Starting point, Y coordinate
//   unsigned int color         Color to paint, in BBGGRRxx format
//   unsigned int boundary      Color to recognize as boundary, (also BBGGRRxx)
//
// * INOUTS   
//   unsigned int *buf          Frame in memory, in BBGGRRxx format (1 pixel
//                              per integer)
// ===========================================================================
#define STACK_SIZE 50176
// WARNING: 196 KB stack, don't use this from a MicroBlaze with the current
//          per-core stack. We call this from the ARM cores, which have
//          a per-core stack of 256 KB.
void vid_draw_floodfill(unsigned int *buf, int buf_width, int buf_height,
                        int x, int y, unsigned int color, 
                        unsigned int boundary) {

  int           stack[STACK_SIZE];
  int           sp;
  int           new_x;
  int           new_y;
  unsigned int  *ptr;
  unsigned int  word;


  // Validate the starting point
  ar_assert((x >= 0) && (y >= 0) && (x < buf_width) && (y < buf_height));
  ptr = buf + y * buf_width + x;
  if ((*ptr == color) || (*ptr == boundary)) {
    return;
  }

  // Initialize the stack
  sp = 0;

  // Push the starting point
  stack[sp++] = (x << 16) | y;

  // Loop
  while (sp) {

    // Pop next point
    word = stack[--sp];
    y = word & 0xFFFF;
    x = (word >> 16) & 0xFFFF;

    // Paint it
    buf[y * buf_width + x] = color;

    // Examine its north, east, south and west neighbors for coloring
    ar_assert(sp + 8 < STACK_SIZE);
    if (y > 0) {
      new_x = x;
      new_y = y - 1;
      ptr = buf + new_y * buf_width + new_x;
      if ((*ptr != color) && (*ptr != boundary)) {
        stack[sp++] = (new_x << 16) | new_y;
      }
    }
    if (x < buf_width - 1) {
      new_x = x + 1;
      new_y = y;
      ptr = buf + new_y * buf_width + new_x;
      if ((*ptr != color) && (*ptr != boundary)) {
        stack[sp++] = (new_x << 16) | new_y;
      }
    }
    if (x > 0) {
      new_x = x - 1;
      new_y = y;
      ptr = buf + new_y * buf_width + new_x;
      if ((*ptr != color) && (*ptr != boundary)) {
        stack[sp++] = (new_x << 16) | new_y;
      }
    }
    if (y < buf_height - 1) {
      new_x = x;
      new_y = y + 1;
      ptr = buf + new_y * buf_width + new_x;
      if ((*ptr != color) && (*ptr != boundary)) {
        stack[sp++] = (new_x << 16) | new_y;
      }
    }
  }
}



// ###########################################################################
// ###                                                                     ###
// ###                     Top-level video management                      ###
// ###                                                                     ###
// ###########################################################################


// ===========================================================================
// vid_init()                   Initialize XUP video peripherals
// ===========================================================================
// * INPUTS
//   int my_bid                 Caller board ID
//   int my_cid                 Caller core ID
// ===========================================================================
void vid_init(int my_bid, int my_cid) {

  // Initialize XUP board I2C master
  vid_i2c_init(my_bid, my_cid);

  // Configure video in/out I2C regs
  vid_out_configure_i2c(my_bid, my_cid);
  vid_in_configure_i2c(my_bid, my_cid, 0);

  // Configure video in/out FPGA peripherals
  vid_out_set_800x600(my_bid, my_cid, 0x01);
  vid_in_set_800x600(my_bid, my_cid, 0x01);
  //vid_out_set_1024x768(my_bid, my_cid, 0x01);
  //vid_in_set_1024x768(my_bid, my_cid, 0x01);

  kt_printf("XUP Video initialized\r\n");
}
