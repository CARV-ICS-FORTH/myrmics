# ===========================================================================
# Copyright (c) 2013, FORTH-ICS / CARV 
#                     (Foundation for Research & Technology -- Hellas,
#                      Institute of Computer Science,
#                      Computer Architecture & VLSI Systems Laboratory)
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# 
# ==========================[ Static Information ]===========================
#
# Author        : Spyros Lyberis
# Abstract      : MicroBlaze boot code, starting at address 0x50
#
# =============================[ CVS Variables ]=============================
#
# File name     : $RCSfile: boot.s,v $
# CVS revision  : $Revision: 1.4 $
# Last modified : $Date: 2013/05/21 06:14:50 $
# Last author   : $Author: zakkak $
#
# ===========================================================================

.text

.global ar_boot

ar_boot:

    # Enable microblaze interrupts
    msrset  r0, 2

    # Setup the stack
    addi    r1, r0, 0x120000

    # Increment stack pointer by 0x20000 * core_id
    lwi     r3, r0, 0xFFF00008
    andi    r3, r3, 7
    bslli   r4, r3, 17
    add     r1, r1, r4

    # Clear BSS (core_id 0 only, visible to all because caches are still off)
    bnei    r3, proceed
    addi    r3, r0, __linker_start_bss
    addi    r4, r0, __linker_end_bss
clear_bss:
    cmp     r5, r4, r3
    bgei    r5, proceed
    sw      r0, r3, r0
    addi    r3, r3, 4
    bri     clear_bss

proceed:
    # Hic sunt dracones...
    brai    main
    nop
