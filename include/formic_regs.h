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
// Abstract      : Formic FPGA board controller registers
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: formic_regs.h,v $
// CVS revision  : $Revision: 1.1 $
// Last modified : $Date: 2012/01/27 12:30:15 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _FORMIC_REGS_H
#define _FORMIC_REGS_H

#define FORMIC_BRD_CONTROL      ((volatile unsigned int *) 0xFFF00000)
#define FORMIC_BRD_STATUS       ((volatile unsigned int *) 0xFFF00004)
#define FORMIC_BRD_LINK_STATUS  ((volatile unsigned int *) 0xFFF00008)
#define FORMIC_BRD_CORES_STATUS ((volatile unsigned int *) 0xFFF0000C)
#define FORMIC_BRD_DRAM_STATUS  ((volatile unsigned int *) 0xFFF00010)

#define FORMIC_TLB_CONTROL      ((volatile unsigned int *) 0xFFF20000)
#define FORMIC_TLB_ENTRY(x)     ((volatile unsigned int *) \
                                 (0xFFF21000 + ((x) << 2)))

#define FORMIC_TMR_BRD_GLOBAL   ((volatile unsigned int *) 0xFFF40000)

#define FORMIC_UART_CONTROL     ((volatile unsigned int *) 0xFFF60000)
#define FORMIC_UART_STATUS      ((volatile unsigned int *) 0xFFF60004)
#define FORMIC_UART_BYTE        ((volatile unsigned int *) 0xFFF60008)
#define FORMIC_UART_BULK_WRITE  ((volatile unsigned int *) 0xFFF60100)

#define FORMIC_TRC_BASE         ((volatile unsigned int *) 0xFFF80000)
#define FORMIC_TRC_BOUND        ((volatile unsigned int *) 0xFFF80004)
#define FORMIC_TRC_CONTROL      ((volatile unsigned int *) 0xFFF80008)

#endif
