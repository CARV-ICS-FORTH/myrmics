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
// Abstract      : ARM on-chip and motherboard registers
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: arm_regs.h,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/03/26 13:14:25 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _ARM_REGS_H
#define _ARM_REGS_H

// Coprocessor assembly macros
#define ARM_SCTLR_READ(arg)   asm("mrc p15, 0, %0, c1, c0, 0" : "=r"(arg) );
#define ARM_SCTLR_WRITE(arg)  asm("mcr p15, 0, %0, c1, c0, 0" : : "r"(arg) );

#define ARM_ACTLR_READ(arg)   asm("mrc p15, 0, %0, c1, c0, 1" : "=r"(arg) );
#define ARM_ACTLR_WRITE(arg)  asm("mcr p15, 0, %0, c1, c0, 1" : : "r"(arg) );

#define ARM_CPACR_READ(arg)   asm("mrc p15, 0, %0, c1, c0, 2" : "=r"(arg) );
#define ARM_CPACR_WRITE(arg)  asm("mcr p15, 0, %0, c1, c0, 2" : : "r"(arg) );

#define ARM_TTBR0_READ(arg)   asm("mrc p15, 0, %0, c2, c0, 0" : "=r"(arg) );
#define ARM_TTBR0_WRITE(arg)  asm("mcr p15, 0, %0, c2, c0, 0" : : "r"(arg) );
#define ARM_TTBR1_READ(arg)   asm("mrc p15, 0, %0, c2, c0, 1" : "=r"(arg) );
#define ARM_TTBR1_WRITE(arg)  asm("mcr p15, 0, %0, c2, c0, 1" : : "r"(arg) );
#define ARM_TTBCR_READ(arg)   asm("mrc p15, 0, %0, c2, c0, 2" : "=r"(arg) );
#define ARM_TTBCR_WRITE(arg)  asm("mcr p15, 0, %0, c2, c0, 2" : : "r"(arg) );

#define ARM_DACR_READ(arg)    asm("mrc p15, 0, %0, c3, c0, 0" : "=r"(arg) );
#define ARM_DACR_WRITE(arg)   asm("mcr p15, 0, %0, c3, c0, 0" : : "r"(arg) );

#define ARM_CLIDR_READ(arg)   asm("mrc p15, 1, %0, c0, c0, 1" : "=r"(arg) );

#define ARM_SP_READ(arg)      asm("mov %0, sp" : "=r"(arg) );


// Snoop Control Unit
#define ARM_SCU_CONTROL               ((volatile int *) 0x1E000000)
#define ARM_SCU_CONFIG                ((volatile int *) 0x1E000004)
#define ARM_SCU_CPU_POWER_STATUS      ((volatile int *) 0x1E000008)
#define ARM_SCU_INVALIDATE_SECURE     ((volatile int *) 0x1E00000C)
#define ARM_SCU_FILTERING_START       ((volatile int *) 0x1E000040)
#define ARM_SCU_FILTERING_END         ((volatile int *) 0x1E000044)
#define ARM_SCU_ACCESS_CONTROL        ((volatile int *) 0x1E000050)
#define ARM_SCU_ACCESS_CONTROL_NS     ((volatile int *) 0x1E000054)

// Timers
#define ARM_PRIV_TIMER_LOAD           ((volatile unsigned int *) 0x1E000600)
#define ARM_PRIV_TIMER_COUNTER        ((volatile unsigned int *) 0x1E000604)
#define ARM_PRIV_TIMER_CONTROL        ((volatile unsigned int *) 0x1E000608)
#define ARM_PRIV_TIMER_INT_STATUS     ((volatile unsigned int *) 0x1E00060C)

#define ARM_WATCHDOG_TIMER_LOAD       ((volatile unsigned int *) 0x1E000620)
#define ARM_WATCHDOG_TIMER_COUNTER    ((volatile unsigned int *) 0x1E000624)
#define ARM_WATCHDOG_TIMER_CONTROL    ((volatile unsigned int *) 0x1E000628)
#define ARM_WATCHDOG_TIMER_INT_STATUS ((volatile unsigned int *) 0x1E00062C)
#define ARM_WATCHDOG_TIMER_RST_STATUS ((volatile unsigned int *) 0x1E000630)
#define ARM_WATCHDOG_TIMER_DISABLE    ((volatile unsigned int *) 0x1E000634)


// Interrupt Controller
#define ARM_ICPICR      ((volatile int *) 0x1E000100) // Processor interface ctl
#define ARM_ICCIPMR     ((volatile int *) 0x1E000104) // Priority mask
#define ARM_ICCBPR      ((volatile int *) 0x1E000108) // Binary point
#define ARM_ICCIAR      ((volatile int *) 0x1E00010C) // Interrupt acknowledge
#define ARM_ICCEOIR     ((volatile int *) 0x1E000110) // End-of-interrupt
#define ARM_ICCRPR      ((volatile int *) 0x1E000114) // Running priority
#define ARM_ICCHPIR     ((volatile int *) 0x1E000118) // Highest pending intrpt
#define ARM_ICCABPR     ((volatile int *) 0x1E00011C) // Aliased insecure bin pt
#define ARM_ICPIIR      ((volatile int *) 0x1E0001FC) // Pr interface ident ID

#define ARM_ICDDCR      ((volatile int *) 0x1E001000) // Distributor control
#define ARM_ICDICTR     ((volatile int *) 0x1E001004) // Interrupt ctl type
#define ARM_ICDDIR      ((volatile int *) 0x1E001008) // Distributor impl ID
#define ARM_ICDISR      ((volatile int *) 0x1E001080) // Interrupt security
#define ARM_ICDISER     ((volatile int *) 0x1E001100) // Enable set
#define ARM_ICDICER     ((volatile int *) 0x1E001180) // Enable clear
#define ARM_ICDISPR     ((volatile int *) 0x1E001200) // Pending set
#define ARM_ICDICPR     ((volatile int *) 0x1E001280) // Pending clear
#define ARM_ICDABR      ((volatile int *) 0x1E001300) // Active status
#define ARM_ICDIPR      ((volatile int *) 0x1E001400) // Priority level
#define ARM_ICDIPTR     ((volatile int *) 0x1E001800) // SPI target
#define ARM_ICDICR      ((volatile int *) 0x1E001C00) // Interrupt configuration
#define ARM_ICDPPISTAT  ((volatile int *) 0x1E001D00) // PPI status
#define ARM_ICDSPISTAT  ((volatile int *) 0x1E001D04) // SPI status
#define ARM_ICDSGIR     ((volatile int *) 0x1E001F00) // Software generated int
#define ARM_PERIPH_ID   ((volatile int *) 0x1E001000) // Peripheral ID
#define ARM_COMP_ID     ((volatile int *) 0x1E001000) // Component ID


// Motherboard UART
#define ARM_UART0_DR    ((volatile int *) 0x10009000) // Data
#define ARM_UART0_RSR   ((volatile int *) 0x10009004) // Rcv status/error clear
#define ARM_UART0_FR    ((volatile int *) 0x10009018) // Flag
#define ARM_UART0_ILPR  ((volatile int *) 0x10009020) // IrDA low-power counter
#define ARM_UART0_IBRD  ((volatile int *) 0x10009024) // Integer baud rate
#define ARM_UART0_FBRD  ((volatile int *) 0x10009028) // Fractional baud rate
#define ARM_UART0_LCR_H ((volatile int *) 0x1000902C) // Line control
#define ARM_UART0_CR    ((volatile int *) 0x10009030) // Control
#define ARM_UART0_IFLS  ((volatile int *) 0x10009034) // Interrupt FIFO lvl sel
#define ARM_UART0_IMSC  ((volatile int *) 0x10009038) // Interrupt mask
#define ARM_UART0_RIS   ((volatile int *) 0x1000903C) // Raw interrupt status
#define ARM_UART0_MIS   ((volatile int *) 0x10009040) // Masked interrupt status
#define ARM_UART0_ICR   ((volatile int *) 0x10009044) // Interrupt clear
#define ARM_UART0_DMACR ((volatile int *) 0x10009048) // DMA control
#define ARM_UART0_PID0  ((volatile int *) 0x10009FE0) // Periperal ID0
#define ARM_UART0_PID1  ((volatile int *) 0x10009FE4) // Periperal ID1
#define ARM_UART0_PID2  ((volatile int *) 0x10009FE8) // Periperal ID2
#define ARM_UART0_PID3  ((volatile int *) 0x10009FEC) // Periperal ID3
#define ARM_UART0_CID0  ((volatile int *) 0x10009FF0) // PrimeCell ID0
#define ARM_UART0_CID1  ((volatile int *) 0x10009FF4) // PrimeCell ID1
#define ARM_UART0_CID2  ((volatile int *) 0x10009FF8) // PrimeCell ID2
#define ARM_UART0_CID3  ((volatile int *) 0x10009FFC) // PrimeCell ID3


// Motherboard CompactFlash
#define ARM_FLASH_CONTROL           ((volatile unsigned int *) 0x1001A300)

#define ARM_FLASH_ACCESS_DATA       ((volatile unsigned int *) 0x1001A000)
#define ARM_FLASH_ACCESS_ERROR      ((volatile unsigned int *) 0x1001A004)
#define ARM_FLASH_ACCESS_NSECT      ((volatile unsigned int *) 0x1001A008)
#define ARM_FLASH_ACCESS_LBAL       ((volatile unsigned int *) 0x1001A00C)
#define ARM_FLASH_ACCESS_LBAM       ((volatile unsigned int *) 0x1001A010)
#define ARM_FLASH_ACCESS_LBAH       ((volatile unsigned int *) 0x1001A014)
#define ARM_FLASH_ACCESS_DEVICE     ((volatile unsigned int *) 0x1001A018)
#define ARM_FLASH_ACCESS_STATUS     ((volatile unsigned int *) 0x1001A01C)
#define ARM_FLASH_ACCESS_CTL        ((volatile unsigned int *) 0x1001A038)
#define ARM_FLASH_ACCESS_FEATURE    ARM_FLASH_ACCESS_ERROR
#define ARM_FLASH_ACCESS_COMMAND    ARM_FLASH_ACCESS_STATUS

#endif
