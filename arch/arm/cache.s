@ ===========================================================================
@ Copyright (c) 2013, FORTH-ICS / CARV 
@                     (Foundation for Research & Technology -- Hellas,
@                      Institute of Computer Science,
@                      Computer Architecture & VLSI Systems Laboratory)
@ 
@ Licensed under the Apache License, Version 2.0 (the "License");
@ you may not use this file except in compliance with the License.
@ You may obtain a copy of the License at
@ 
@     http://www.apache.org/licenses/LICENSE-2.0
@ 
@ Unless required by applicable law or agreed to in writing, software
@ distributed under the License is distributed on an "AS IS" BASIS,
@ WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
@ See the License for the specific language governing permissions and
@ limitations under the License.
@ 
@ ==========================[ Static Information ]===========================
@
@ Author        : Spyros Lyberis
@ Abstract      : Cache, Branch predictor, FPU initialization, control and
@                 maintenance functions
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: cache.s,v $
@ CVS revision  : $Revision: 1.5 $
@ Last modified : $Date: 2012/01/27 16:32:59 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================

.text
.align 4

.global ar_invalidate_caches

.global ar_disable_icache
.global ar_enable_icache

.global ar_disable_dcache
.global ar_enable_dcache

.global ar_disable_branch_pred
.global ar_enable_branch_pred

.global ar_enable_fpu

.global ar_get_core_id
.global ar_get_board_id


@ =============================================================================
@ Get core ID
@
@ Reference: Multiprocessor Affinity Register (MPIDR), ARMARM B3.12.11
@ =============================================================================
ar_get_core_id:
        mrc     p15, 0, r0, c0, c0, 5           @ r0 = MPIDR
        and     r0, r0, #0xf                    @ r0 = core id
        bx      lr

@ =============================================================================
@ Get board ID
@
@ Reference: ARS register CPU_CONTROL (0xFFC00008 + 0x00100000 * core_id)
@ =============================================================================
ar_get_board_id:
        push    {LR}
        bl      ar_get_core_id                  @ r0 = core_id
        lsl     r0, r0, #20                     @ r0 = core_id << 20
        ldr     r1, =0xFFC00008
        add     r0, r0, r1                      @ r0 += 0xFFC00000
        ldr     r0, [r0]                        @ r0 = *r0
        lsr     r0, #24                         @ r0 >>= 24
        pop     {PC}
        bx      lr


@ =============================================================================
@ Invalidate all caches
@
@ Reference: Cache enabling and disabling, ARMARM B2-8
@ =============================================================================
ar_invalidate_caches:
        push    {r4-r8,r10-r11}            @ Save registers on stack

        mov     r0, #0
        mcr     p15, 0, r0, c8, c7, 0      @ Invalidate Inst TLB and Data TLB
        mcr     p15, 0, r0, c7, c5, 0      @ Invalidate I Cache
        mcr     p15, 0, r0, c7, c5, 6      @ Invalidate branch predictor array
    
        @ Must iterate over the caches in order to synthesise a complete clean
        @ of data/unified cache
        mrc     p15, 1, r0, c0, c0, 1      @ read Cache Level ID register (clidr)
        ands    r3, r0, #0x7000000         @ extract level of coherency from clidr
        mov     r3, r3, lsr #23            @ left align level of coherency bit field
        beq     ar_inv_finished            @ if loc is 0, then no need to clean
        
        mov     r10, #0                    @ start clean at cache level 0 (in r10)
ar_inv_loop1:
        add     r2, r10, r10, lsr #1       @ work out 3x current cache level
        mov     r1, r0, lsr r2             @ extract cache type bits from clidr
        and     r1, r1, #7                 @ mask of the bits for current cache only
        cmp     r1, #2                     @ see what cache we have at this level
        blt     ar_inv_skip                @ skip if no cache, or just i-cache
        mcr     p15, 2, r10, c0, c0, 0     @ select current cache level in cssr
        mov     r1, #0
        mcr     p15, 0, r1, c7, c5, 4      @ prefetchflush to synch the new cssr&csidr
        mrc     p15, 1, r1, c0, c0, 0      @ read the new csidr
        and     r2, r1, #7                 @ extract the length of the cache lines
        add     r2, r2, #4                 @ add 4 (line length offset)
        ldr     r6, =0x3ff
        ands    r6, r6, r1, lsr #3         @ find maximum number on the way size
        clz     r5, r6                     @ find bit position of way size increment
        ldr     r7, =0x7fff
        ands    r7, r7, r1, lsr #13        @ extract max number of the index size
ar_inv_loop2:
        mov     r8, r6                     @ create working copy of max way size
ar_inv_loop3:
        orr     r11, r10, r8, lsl r5       @ factor way and cache number into r11
        orr     r11, r11, r7, lsl r2       @ factor index number into r11
        mcr     p15, 0, r11, c7, c6, 2     @ invalidate by set/way
        subs    r8, r8, #1                 @ decrement the way
        bge     ar_inv_loop3
        subs    r7, r7, #1                 @ decrement the index
        bge     ar_inv_loop2
ar_inv_skip:
        add     r10, r10, #2               @ increment cache number
        cmp     r3, r10
        bgt     ar_inv_loop1
ar_inv_finished:
        pop     {r4-r8,r10-r11}            @ Load registers from stack
        bx      lr                         @ Return


@ =============================================================================
@ Turn off instruction cache
@ =============================================================================
ar_disable_icache:
        mrc     p15, 0, r0, c1, c0, 0      @ Read SCTLR
        bic     r0, r0, #(1 << 12)         @ Turn off I bit
        mcr     p15, 0, r0, c1, c0, 0      @ Write SCTLR
        bx      lr                         @ Return

@ =============================================================================
@ Turn on instruction cache
@ =============================================================================
ar_enable_icache:
        mrc     p15, 0, r0, c1, c0, 0      @ Read SCTLR
        orr     r0, r0, #(1 << 12)         @ Turn on I bit
        mcr     p15, 0, r0, c1, c0, 0      @ Write SCTLR
        bx      lr                         @ Return

@ =============================================================================
@ Turn off data cache
@ =============================================================================
ar_disable_dcache:
        mrc     p15, 0, r0, c1, c0, 0      @ Read SCTLR
        bic     r0, r0, #(1 << 2)          @ Turn off C bit
        mcr     p15, 0, r0, c1, c0, 0      @ Write SCTLR
        bx      lr                         @ Return

@ =============================================================================
@ Turn on data cache
@ =============================================================================
ar_enable_dcache:
        mrc     p15, 0, r0, c1, c0, 0      @ Read SCTLR
        orr     r0, r0, #(1 << 2)          @ Turn on C bit
        mcr     p15, 0, r0, c1, c0, 0      @ Write SCTLR
        bx      lr                         @ Return

@ =============================================================================
@ Turn off branch prediction
@ =============================================================================
ar_disable_branch_pred:
        mrc     p15, 0, r0, c1, c0, 0      @ Read SCTLR
        bic     r0, r0, #(1 << 11)         @ Turn off Z bit
        mcr     p15, 0, r0, c1, c0, 0      @ Write SCTLR
        bx      lr                         @ Return

@ =============================================================================
@ Turn on branch prediction
@ =============================================================================
ar_enable_branch_pred:
        mrc     p15, 0, r0, c1, c0, 0      @ Read SCTLR
        orr     r0, r0, #(1 << 11)         @ Turn on Z bit
        mcr     p15, 0, r0, c1, c0, 0      @ Write SCTLR
        bx      lr                         @ Return

@ =====================================================================
@ Enable the FPU
@ =====================================================================
ar_enable_fpu:
        @ Reference: Coprocessor Access Control Register, ARMARM B3-104
        @ Enable CP10 and CP11 access
        mov     r0, #0xF00000
        mcr     p15, 0, r0, c1, c0, 2

        @ Enable the FPEXC EN bit
        @ Reference: Cortex-A9 FPU manual
        mov     r0, #0x40000000
        msr     fpexc, r0

        bx lr

