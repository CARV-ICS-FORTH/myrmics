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
// Abstract      : Versatile FPGA ARS block registers
//
// =============================[ CVS Variables ]=============================
// 
// File name     : $RCSfile: ars_regs.h,v $
// CVS revision  : $Revision: 1.2 $
// Last modified : $Date: 2012/01/27 12:30:15 $
// Last author   : $Author: lyberis-spree $
// 
// ===========================================================================

#ifndef _ARS_REGS_H
#define _ARS_REGS_H

#define ARS_CPU_STATUS(x)          ((volatile unsigned int *) \
                                    (0xFFC00008 + ((x) << 20)))

#define ARS_PRF_REGS(x)            ((volatile unsigned int *) \
                                    (0xFFC60010 + ((x) << 20)))
#define ARS_PRF_NET_OPS(x)         ((volatile unsigned int *) \
                                    (0xFFC60014 + ((x) << 20)))
#define ARS_PRF_PACKETS0(x)        ((volatile unsigned int *) \
                                    (0xFFC60020 + ((x) << 20)))
#define ARS_PRF_PACKETS1(x)        ((volatile unsigned int *) \
                                    (0xFFC60024 + ((x) << 20)))
#define ARS_PRF_PACKETS2(x)        ((volatile unsigned int *) \
                                    (0xFFC60028 + ((x) << 20)))

#define ARS_TMR_GLOBAL(x)          ((volatile unsigned int *) \
                                    (0xFFC80000 + ((x) << 20)))

#define ARS_MNI_OPCODE(x)          ((volatile unsigned int *) \
                                    (0xFFCA0000 + ((x) << 20)))
#define ARS_MNI_DMA_SIZE(x)        ((volatile unsigned int *) \
                                    (0xFFCA0004 + ((x) << 20)))
#define ARS_MNI_DST_ADR(x)         ((volatile unsigned int *) \
                                    (0xFFCA0008 + ((x) << 20)))
#define ARS_MNI_BRD_NODE(x)        ((volatile unsigned int *) \
                                    (0xFFCA000C + ((x) << 20)))
#define ARS_MNI_ACK_ADR(x)         ((volatile unsigned int *) \
                                    (0xFFCA0010 + ((x) << 20)))
#define ARS_MNI_SRC_ADR(x)         ((volatile unsigned int *) \
                                    (0xFFCA0014 + ((x) << 20)))
#define ARS_MNI_MSG0(x)            ((volatile unsigned int *) \
                                    (0xFFCA0018 + ((x) << 20)))
#define ARS_MNI_MSG1(x)            ((volatile unsigned int *) \
                                    (0xFFCA001C + ((x) << 20)))
#define ARS_MNI_STATUS(x)          ((volatile unsigned int *) \
                                    (0xFFCA0020 + ((x) << 20)))

#define ARS_CNT_VAL(x, y)          ((volatile unsigned int *) \
                                    (0xFFCC1000 + ((x) << 20) + ((y) << 4)))
#define ARS_CNT_NTF_ADR0(x, y)     ((volatile unsigned int *) \
                                    (0xFFCC1004 + ((x) << 20) + ((y) << 4)))
#define ARS_CNT_NTF_ADR1(x, y)     ((volatile unsigned int *) \
                                    (0xFFCC1008 + ((x) << 20) + ((y) << 4)))
#define ARS_CNT_NTF_BRD_NODE(x, y) ((volatile unsigned int *) \
                                    (0xFFCC100C + ((x) << 20) + ((y) << 4)))

#define ARS_MBX_ACCESS(x)          ((volatile unsigned int *) \
                                    (0xFFCE0000 + ((x) << 20)))
#define ARS_MBX_STATUS(x)          ((volatile unsigned int *) \
                                    (0xFFCE0004 + ((x) << 20)))
#define ARS_MSL_ACCESS(x)          ((volatile unsigned int *) \
                                    (0xFFCF0000 + ((x) << 20)))
#define ARS_MSL_STATUS(x)          ((volatile unsigned int *) \
                                    (0xFFCF0004 + ((x) << 20)))

#endif
