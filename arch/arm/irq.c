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
// Abstract      : ARM interrupt handler
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: irq.c,v $
// CVS revision  : $Revision: 1.3 $
// Last modified : $Date: 2012/03/15 15:31:47 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================


#include <arch.h>
#include <kernel_toolset.h>
#include <filesystem.h>
#include <memory_management.h>
#include <ars_regs.h>



// ============================================================================
// Interrupts
// ============================================================================
#define INT_SGI0                 0      // Software-generated, 0 - 15
#define INT_SGI1                 1 
#define INT_SGI2                 2 
#define INT_SGI3                 3 
#define INT_SGI4                 4 
#define INT_SGI5                 5 
#define INT_SGI6                 6 
#define INT_SGI7                 7 
#define INT_SGI8                 8 
#define INT_SGI9                 9 
#define INT_SGI10               10
#define INT_SGI11               11
#define INT_SGI12               12
#define INT_SGI13               13
#define INT_SGI14               14
#define INT_SGI15               15

#define INT_PRIV_LEGACY_NFIQ    28      // PPI 0 - 3
#define INT_PRIV_TIMER          29
#define INT_PRIV_WATCHDOG       30
#define INT_PRIV_LEGACY_NIRQ    31

#define INT_MB_WATCHDOG         32      // Motherboard SB_IRQ 0 - 42
#define INT_MB_SOFTWARE         33
#define INT_MB_TIMER01          34
#define INT_MB_TIMER23          35
#define INT_MB_REAL_TIME_CLOCK  36
#define INT_MB_UART0            37
#define INT_MB_UART1            38
#define INT_MB_UART2            39
#define INT_MB_UART3            40
#define INT_MB_MEDIA_CARD0      41
#define INT_MB_MEDIA_CARD1      42
#define INT_MB_AUDIO_CODEC      43
#define INT_MB_KEYB_MOUSE0      44
#define INT_MB_KEYB_MOUSE1      45
#define INT_MB_DISPLAY          46
#define INT_MB_ETHERNET         47
#define INT_MB_USB              48
#define INT_MB_PCI_EXPRESS      49
#define INT_MB_SITE1_0          64
#define INT_MB_SITE1_1          65
#define INT_MB_SITE1_2          66
#define INT_MB_SITE1_3          67
#define INT_MB_SITE2_0          68
#define INT_MB_SITE2_1          69
#define INT_MB_SITE2_2          70
#define INT_MB_SITE2_3          71

#define INT_TC_L2_CACHE         75      // Test-chip peripherals
#define INT_TC_DISPLAY          76
#define INT_TC_SMC0             77
#define INT_TC_SMC1             78
#define INT_TC_NMC              79
#define INT_TC_TIMER0           80
#define INT_TC_TIMER1           81
#define INT_TC_GPIO             82
#define INT_TC_WATCHDOG         83
#define INT_TC_UART             84
#define INT_TC_PERF_CPU0        92
#define INT_TC_PERF_CPU1        93
#define INT_TC_PERF_CPU2        94
#define INT_TC_PERF_CPU3        95

#define INT_SPURIOUS          1023      // Special ID: spurious interrupt



  // ==========================================================================
  // Interrupts
  // ==========================================================================
  // asm (
  //       "mrs     r0, CPSR\n"
  //       "bic     r0, r0, #(0xC0)\n"  // Enable IRQs
  //       "msr     CPSR_c, r0\n"
  //       : : : "r0"
  //     );

  // send_banner("= Interrupts enabled\r\n");

  // // Enable UART0 receive interrupt, at 1/8 FIFO level = 2 chars
  // u32 = *ARM_UART0_IMSC;
  // *ARM_UART0_IMSC = (u32 | (1 << 4));
  // *ARM_UART0_IFLS = 0;

  // // Activate MB SWINT
  // //u32 = *((vu32_t *) 0x10000060);
  // //u32 |= (1 << 19); // SWINT
  // //*((vu32_t *) 0x10000060) = u32;
  // u32 = *((vu32_t *) 0x10000060);
  // ar_uart_send_str("MB SYS_MISC = [");
  // ar_uart_send_hex(u32);
  // ar_uart_send_str("]\r\n");
  // 
  // // Disable distributor
  // *ARM_ICDDCR = 0;

  // // Disable processor interface
  // *ARM_ICPICR = 0;

  // // Processor at lowest priority
  // *ARM_ICCIPMR = 0xFF;
  // 
  // // Read controller type
  // ar_uart_send_str("ICDICTR is [");
  // ar_uart_send_hex(*ARM_ICDICTR);
  // ar_uart_send_str("]\r\n");

  // // Set all interrupt targets to core0 only
  // for (i = 0; i < 96/4; i++) {
  //   *(ARM_ICDIPTR + i) = 0x01010101;
  // }

  // // Set all set-enable bits to 1
  // //for (i = 0; i < 96/32; i++) {
  // //  *((vu32_t *) (0x1E001100 + i * 4)) = 0xFFFFFFFF;
  // //  ar_uart_send_str("ICDISER");
  // //  ar_uart_send_dec(i);
  // //  ar_uart_send_str(" = [");
  // //  ar_uart_send_hex(*((vu32_t *) (0x1E001100 + i * 4)));
  // //  ar_uart_send_str("]\r\n");
  // //}

  // // Disable all interrupts
  // for (i = 0; i < 96/32; i++) {
  //   *(ARM_ICDICER + i) = 0xFFFFFFFF;
  //   ar_uart_send_str("ICDICER");
  //   ar_uart_send_dec(i);
  //   ar_uart_send_str(" = [");
  //   ar_uart_send_hex(*(ARM_ICDICER + i));
  //   ar_uart_send_str("]\r\n");
  // }

  // // Clear all pending interrupts
  // ar_uart_send_str("\r\n");
  // for (i = 0; i < 96/32; i++) {
  //   *(ARM_ICDICPR + i) = 0xFFFFFFFF;
  //   ar_uart_send_str("ICDICPR");
  //   ar_uart_send_dec(i);
  //   ar_uart_send_str(" = [");
  //   ar_uart_send_hex(*(ARM_ICDICPR + i));
  //   ar_uart_send_str("]\r\n");
  // }


  // // Activate MB SWINT
  // ar_enable_interrupt(INT_MB_SOFTWARE);


  // // Activate MB UART0
  // ar_enable_interrupt(INT_MB_UART0);



  // // Enable processor interface
  // *ARM_ICPICR = 1;

  // // Enable distributor
  // *ARM_ICDDCR = 1;



  // //*((vu32_t *) 0x1E000000) = 1;
  // ar_uart_send_str("SCU Control is 0x");
  // ar_uart_send_hex(*((vu32_t *) 0x1E000000));
  // ar_uart_send_str("\r\n");
  // ar_uart_send_str("SCU Config is 0x");
  // ar_uart_send_hex(*((vu32_t *) 0x1E000004));
  // ar_uart_send_str("\r\n");


  //asm("mov r0, #0x60000000\n"
  //    "ldr r1, [r0, #1]\n"
  //    : : : "r0", "r1"
  //   );


  // Shutdown = 8, Reboot = 9
  // ar_uart_send_str("Rebooting...\r\n");
  // *((vu32_t *) 0x100000A4) = 
  //                               (1 << 31) | // Start bit
  //                               (1 << 30) | // Write bit
  //                               (8 << 20) | // Function bits
  //                               (0 << 16);  // 0 = Motherboard
  // ar_uart_send_str("Should not be here!\r\n");
  // 
  // while (1);
    
    
    
    
    

/*
void ar_enable_interrupt(u32_t interrupt) {
  int    offset;
  int    bit;
  vu32_t reg;

  // Find register offset in distributor registers (32 interrupts/reg)
  offset = interrupt / 32;
  bit    = interrupt % 32;

  // Read set-enable reg
  reg = *(ARM_ICDISER + offset);

  // Activate bit
  reg |= 1 << bit;

  // Write set-enable reg
  *(ARM_ICDISER + offset) = reg;
}



void ar_irq_handler() {
  u32_t ack;

  ar_uart_send_str("Enter IRQ handler\r\n");

  // Acknowledge interrupt
  ack = *ARM_ICCIAR;

  ar_uart_send_str("Interrupt ID = ");
  ar_uart_send_dec(ack & 0x1FF);
  ar_uart_send_str(", CPU ID = ");
  ar_uart_send_dec((ack >> 10) & 0xF);
  ar_uart_send_str("\r\n");

  switch (ack & 0x1FF) {
    case 37: // UART0
      ar_uart_send_str("RX FIFO: [");
      ar_uart_send_dec(*((volatile unsigned int *) 0x10009000) & 0xFF);
      ar_uart_send_str("] and [");
      ar_uart_send_dec(*((volatile unsigned int *) 0x10009000) & 0xFF);
      ar_uart_send_str("\r\n");
      break;
    default:
      ar_uart_send_str(">> Interrupt unsupported by IRQ handler\r\n");
  }


  // End of interrupt
  *ARM_ICCEOIR = ack;

  ar_uart_send_str("Exit IRQ handler\r\n");
}
*/
