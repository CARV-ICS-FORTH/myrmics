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
@ Abstract      : Exception handlers, to be relocated to the Exception Base
@                 Address (0x0 currently). Booting starts here, by taking the 
@                 exception_reset branch, which simply points to the boot.s 
@                 entry point.
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: exceptions.s,v $
@ CVS revision  : $Revision: 1.6 $
@ Last modified : $Date: 2013/03/16 15:32:35 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================


@ =============================================================================
@ Messages
@ =============================================================================
.data
msg_undefined:  .asciz "\r\n\nUndefined instruction exception taken\r\n"
msg_svc:        .asciz "\r\n\nSVC exception taken\r\n"
msg_prefetch:   .asciz "\r\n\nPrefetch abort exception taken\r\n"
msg_data:       .asciz ": Data abort exception taken\r\n"
msg_unused:     .asciz "\r\n\nUnused exception taken\r\n"
msg_irq:        .asciz "\r\n\nIRQ exception taken\r\n"
msg_fiq:        .asciz "\r\n\nFIQ exception taken\r\n"

msg_fin:        .asciz "\r\n\nIRQ exception finished\r\n"

@ =============================================================================
@ The following 8 instructions are the entry points for the 8 defined
@ exceptions. This code is loaded exactly at the Exception Base Address.
@ =============================================================================

.text
.align 4
.global ar_vectors

ar_vectors:
        b       exception_reset
        b       exception_undefined
        b       exception_svc
        b       exception_prefetch
        b       exception_data
        b       exception_unused
        b       exception_irq

exception_fiq:
        @ =====================================================================
        @ Exception 0x1C: Fast interrupt.
        @ This is the last exception and is handled directly on entry, so that
        @ a branch instruction is saved.
        @ =====================================================================
        ldr     r0, =msg_fiq
        bl      ar_asm_print_str
        b       busy_wait

exception_reset:
        @ =====================================================================
        @ Exception 0x00: Reset. This is also the entry point for normal
        @ boot, so we assume that everything is at reset state, all 
        @ processor state (caches, interrupts, etc.) is at default 'off'
        @ positions and so on. Hence, we just jump to normal program entry.
        @ =====================================================================
        b       ar_boot

exception_undefined:
        @ =====================================================================
        @ Exception 0x04: Undefined instruction
        @ =====================================================================
        ldr     r0, =msg_undefined
        bl      ar_asm_print_str
        b       busy_wait

exception_svc:
        @ =====================================================================
        @ Exception 0x08: Supervisor Call (SVC)
        @ =====================================================================
        ldr     r0, =msg_svc
        bl      ar_asm_print_str
        b       busy_wait

exception_prefetch:
        @ =====================================================================
        @ Exception 0x0C: Prefetch abort
        @ =====================================================================
        ldr     r0, =msg_prefetch
        bl      ar_asm_print_str
        b       busy_wait

exception_data:
        @ =====================================================================
        @ Exception 0x10: Data abort
        @ =====================================================================
        mrc     p15, 0, r4, c0, c0, 5           @ r4 = cpu id
        and     r4, r4, #0xf

        lsl     r0, r4, #20                     @ r0 = core_id << 20
        ldr     r1, =0xFFC00008
        add     r0, r0, r1                      @ r0 = ARS_CPU_STATUS(core_id)
        ldr     r5, [r0]
        lsr     r5, r5, #24                     @ r5 = board id

        @ Print board and core ID
        mov     r0, #'\r'
        bl      ar_asm_print_char
        mov     r0, #'\n'
        bl      ar_asm_print_char
        mov     r0, #'\n'
        bl      ar_asm_print_char
        mov     r0, #'0'
        bl      ar_asm_print_char
        mov     r0, #'x'
        bl      ar_asm_print_char
        mov     r0, r5
        bl      ar_asm_print_hex
        mov     r0, #'/'
        bl      ar_asm_print_char
        add     r0, r4, #48
        bl      ar_asm_print_char

        @ Print abort message and stall
        ldr     r0, =msg_data
        bl      ar_asm_print_str
        b       busy_wait

exception_unused:
        @ =====================================================================
        @ Exception 0x14: Unused (this exception is undefined)
        @ =====================================================================
        ldr     r0, =msg_unused
        bl      ar_asm_print_str
        b       busy_wait

exception_irq:
        @ =====================================================================
        @ Exception 0x18: IRQ (normal interrupt)
        @ =====================================================================
        sub     lr, #4                  @ Decrement return address, so that
                                        @ interrupted instruction executes 
                                        @ again
        srsdb   sp!, #0x12              @ Push LR on IRQ stack
        cpsid   i, #0x12                @ Disable further interrupts

        ldr     r0, =msg_irq
        bl      ar_asm_print_str

        b       busy_wait

        @mov     r0, #'S'
        @bl      ar_asm_print_char
        @mov     r0, #'P'
        @bl      ar_asm_print_char
        @mov     r0, #'='
        @bl      ar_asm_print_char
        @mov     r0, #'['
        @bl      ar_asm_print_char
        @mov     r0, sp
        @bl      ar_asm_print_hex
        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #'\r'
        @bl      ar_asm_print_char
        @mov     r0, #'\n'
        @bl      ar_asm_print_char

        @bl      ar_irq_handler

        @ldr     r0, =msg_fin
        @bl      ar_asm_print_str

        @rfeia   sp!                     @ Pop LR from IRQ stack and return


        @ Do nothing
busy_wait:
        b       busy_wait

