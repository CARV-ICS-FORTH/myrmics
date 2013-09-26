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
@ Abstract      : First instructions to be executed upon boot
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: boot.s,v $
@ CVS revision  : $Revision: 1.10 $
@ Last modified : $Date: 2012/11/01 17:20:26 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================

.data
msg_smc:        .asciz     "SMC init ok.\r\n"
msg_dmc:        .asciz     "DMC init ok.\r\n"
msg_bus:        .asciz     "AXI Bus init ok.\r\n"
msg_scu:        .asciz     "SCU Config = ["
msg_remap:      .asciz     "Remap = ["
msg_tst:        .asciz     "DDR 0x1000 = ["
msg_tst2:       .asciz     "UART_FR = ["
msg_reloc1:     .asciz     "Relocated 0x"
msg_reloc2:     .asciz     " bytes.\r\n"

.text
.align 4
.global ar_boot
.global ar_asm_print_str
.global ar_asm_print_char
.global ar_asm_print_hex


@ =============================================================================
@ Definitions
@ =============================================================================

.equ    UART0_BASE,          0x10009000
.equ    UART_FR,             0x18

.equ    BRD_BASE,            0x10000000
.equ    SYS_PROCID0,         0x84
.equ    SYS_PROCID1,         0x88


.equ    CPSR_MODE_USR,       0x10               @ Processor modes
.equ    CPSR_MODE_FIQ,       0x11
.equ    CPSR_MODE_IRQ,       0x12
.equ    CPSR_MODE_SVC,       0x13
.equ    CPSR_MODE_MON,       0x16
.equ    CPSR_MODE_ABT,       0x17
.equ    CPSR_MODE_UND,       0x1B
.equ    CPSR_MODE_SYS,       0x1F


.equ    PA_CODE_BASE,        0x80000000         @ Physical address 64-MB window
                                                @ remapped to address 0x0

.equ    VA_STACK0_SVC_BASE,  0x0013FFFC         @ Core0 236 KB Supervisor stack
.equ    VA_STACK0_MON_BASE,  0x00104FFC         @ Core0   4 KB Monitor stack
.equ    VA_STACK0_ABT_BASE,  0x00103FFC         @ Core0   4 KB Abort stack
.equ    VA_STACK0_UND_BASE,  0x00102FFC         @ Core0   4 KB Undefined stack
.equ    VA_STACK0_IRQ_BASE,  0x00101FFC         @ Core0   4 KB IRQ stack
.equ    VA_STACK0_FIQ_BASE,  0x00100FFC         @ Core0   4 KB FIQ stack

                                                @ Stacks for cores 1-3 are at
                                                @ 0x40000 * core_id above the
                                                @ core 0 stack.


@ =============================================================================
@ Entry point (after reset exception has been taken, it branches here)
@ =============================================================================

ar_boot:
        @ =====================================================================
        @ Disable interrupts, turn off MMU and d-cache. Enable i-cache,
        @ otherwise reads from 0x1000_0000 motherboard peripherals hang
        @ for some reason.
        @ =====================================================================
        mrs     r0, CPSR
        orr     r0, r0, #(0xC0)                 @ Disable IRQ & FIQ
        msr     CPSR_c, r0

        @ Reference: System Control Register, ARMARM B3-96
        mrc     p15, 0, r0, c1, c0, 0           @ Get control register
        bic     r0, r0, #(1 << 0)               @ Disable MMU
        bic     r0, r0, #(1 << 2)               @ Disable D Cache
        @orr     r0, r0, #(1 << 2)               @ Enable D Cache
        @bic     r0, r0, #(1 << 12)              @ Disable I Cache
        orr     r0, r0, #(1 << 12)              @ Enable I Cache
        mcr     p15, 0, r0, c1, c0, 0           @ Write control register

        @ =====================================================================
        @ All cores are booting. Make cores 1-3 of ARM0 (board ID 0x6B) and
        @ cores 0-3 of ARM1 (board ID 0x4B) block on their mailbox, so that
        @ only a single ARM core is booting even if both ARM boards are on.
        @ =====================================================================

        @ Reference: Multiprocessor Affinity Register, ARMARM B3.12.11
        mrc     p15, 0, r3, c0, c0, 5           @ r3 = cpu id
        and     r3, r3, #0xf

        lsl     r0, r3, #20                     @ r0 = core_id << 20
        ldr     r1, =0xFFC00008
        add     r0, r0, r1                      @ r0 = ARS_CPU_STATUS(core_id)
        ldr     r6, [r0]
        lsr     r6, r6, #24                     @ r6 = board id

        cmp     r3, #0                          @ Branch accordingly
        bne     cores123567

        b       cores04


        @ =====================================================================
        @ Initialize UART0, so we can print debug info through the serial port
        @ =====================================================================
cores04:
        ldr     r0, =UART0_BASE                 @ UART0 base adr

        mov     r1, #0
        str     r1, [r0, #0x30]                 @ Disable UART

        mov     r1, #0x27                       @ IBRD = d_39
        str     r1, [r0, #0x24]

        mov     r1, #0x5                        @ FBRD = d_5
        str     r1, [r0, #0x28]

        mov     r1, #0x60                       @ LCRH: disable FEN
        str     r1, [r0, #0x2C]
        mov     r1, #0x70                       @ LCRH: enable FEN, 8-N-1
        str     r1, [r0, #0x2C]

        mov     r1, #0x0300                     @ CR: enable TX/RX
        str     r1, [r0, #0x30]

        orr     r1, #0x01                       @ Enable UART
        str     r1, [r0, #0x30]


        @ =====================================================================
        @ Make core 0 of ARM1 to also wait on mailbox (as cores 1-3 of both
        @ boards) before continuing with the ARM1 initialization
        @ =====================================================================
        lsl     r0, r3, #20                     @ r0 = core_id << 20
        ldr     r1, =0xFFCE0000
        add     r0, r0, r1                      @ r0 = ARS_MBX_ACCESS(core_id)
        cmp     r6, #0x4B                       @ block on mailbox read, if 
        ldreq   r0, [r0]                        @ we're ARM1


        @ =====================================================================
        @ Greet, print PC, SP, FP
        @ =====================================================================
        ldr     r0, =msg_greet
        cmp     r6, #0x6B                       @ ARM0 prints the Myrmics ant
        bleq    ar_asm_print_str

        @mov     r0, #'P'
        @bl      ar_asm_print_char
        @mov     r0, #'C'
        @bl      ar_asm_print_char
        @mov     r0, #'='
        @bl      ar_asm_print_char
        @mov     r0, #'['
        @bl      ar_asm_print_char
        @mov     r0, pc
        @bl      ar_asm_print_hex
        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #','
        @bl      ar_asm_print_char
        @mov     r0, #' '
        @bl      ar_asm_print_char

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
        @mov     r0, #','
        @bl      ar_asm_print_char
        @mov     r0, #' '
        @bl      ar_asm_print_char

        @mov     r0, #'S'
        @bl      ar_asm_print_char
        @mov     r0, #'C'
        @bl      ar_asm_print_char
        @mov     r0, #'T'
        @bl      ar_asm_print_char
        @mov     r0, #'L'
        @bl      ar_asm_print_char
        @mov     r0, #'R'
        @bl      ar_asm_print_char
        @mov     r0, #'='
        @bl      ar_asm_print_char
        @mov     r0, #'['
        @bl      ar_asm_print_char
        @mrc     p15, 0, r0, c1, c0, 0
        @bl      ar_asm_print_hex
        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #'\r'
        @bl      ar_asm_print_char
        @mov     r0, #'\n'
        @bl      ar_asm_print_char


        @ =====================================================================
        @ Initialize Static Memory Bus
        @ =====================================================================
        bl      ar_bus_smc

        @ldr     r0, =msg_smc
        @bl      ar_asm_print_str


        @ =====================================================================
        @ Initialize Dynamic Memory Controller
        @ =====================================================================
        bl      ar_bus_dmc

        @ldr     r0, =msg_dmc
        @bl      ar_asm_print_str


        @ =====================================================================
        @ Initialize On-Chip PL301 AXI Bus
        @ =====================================================================
        bl      ar_bus_axi

        @ldr     r0, =msg_bus
        @bl      ar_asm_print_str


        @ =====================================================================
        @ Check board IDs for expected values
        @ =====================================================================

        @@ Check SCU config
        @ldr     r0, =msg_scu
        @bl      ar_asm_print_str

        @ldr     r0, =0x1E000000
        @ldr     r0, [r0, #4]

        @bl      ar_asm_print_hex

        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #','
        @bl      ar_asm_print_char
        @mov     r0, #' '
        @bl      ar_asm_print_char

        @@ Check remap bits
        @ldr     r0, =msg_remap
        @bl      ar_asm_print_str

        @ldr     r0, =0x100E2000
        @ldr     r0, [r0, #4]

        @bl      ar_asm_print_hex

        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #'\r'
        @bl      ar_asm_print_char
        @mov     r0, #'\n'
        @bl      ar_asm_print_char

        @@ ddr tst
        @ldr     r0, =msg_tst
        @bl      ar_asm_print_str

        @ldr     r0, =0x80001000
        @ldr     r0, [r0, #0]

        @bl      ar_asm_print_hex

        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #','
        @bl      ar_asm_print_char
        @mov     r0, #' '
        @bl      ar_asm_print_char


        @ =====================================================================
        @ Relocate kernel code to DDR. 
        @
        @ Note that all the following relocation functions may copy 4 bytes
        @ more than necessary, since they use doublewords.
        @ =====================================================================
        ldr     r0, =__linker_start_code        @ Start of code from remapped
                                                @ FLASH space (PA 0x0000_0000)
        ldr     r1, =__mb_end_include           @ End of included MB ELF
        ldr     r2, =PA_CODE_BASE               @ Relocation start point...
        sub     r2, r2, r0                      @ ... minus start of code.

code_reloc_loop:
        ldrd    r4, r5, [r0]                    @ Load 2 words from src adr...
        strd    r4, r5, [r0, r2]                @ ... and store them to dst.
        add     r0, #8                          @ Add 8 bytes
        cmp     r0, r1                          @ Loop
        blt     code_reloc_loop

        ldr     r1, =__linker_start_code        @ Say how many bytes we copied
        sub     r3, r0, r1
        ldr     r0, =msg_reloc1
        bl      ar_asm_print_str
        mov     r0, r3
        bl      ar_asm_print_hex
        ldr     r0, =msg_reloc2
        bl      ar_asm_print_str


        @ =====================================================================
        @ Remap first 16 MB to DDR and prepare core0 stacks for all modes.
        @
        @ We are still using physical addresses, but we are now going to
        @ remap the first MBs to point to the DDR space, and not FLASH.
        @ Therefore, we will be ready to switch to virtual addresses later,
        @ as long as they match the physical ones.
        @ =====================================================================
        ldr     r0, =0x100E2000                 @ Remap first MBs to DDR
        ldr     r1, [r0, #0x04]
        orr     r1, r1, #0x40000000  
        str     r1, [r0, #0x04]

        msr     cpsr_c, #(CPSR_MODE_MON)        @ Set up stack pointers for
        ldr     sp, =VA_STACK0_MON_BASE         @ all processor modes of 
        msr     cpsr_c, #(CPSR_MODE_ABT)        @ core 0
        ldr     sp, =VA_STACK0_ABT_BASE
        msr     cpsr_c, #(CPSR_MODE_UND)
        ldr     sp, =VA_STACK0_UND_BASE
        msr     cpsr_c, #(CPSR_MODE_IRQ)
        ldr     sp, =VA_STACK0_IRQ_BASE
        msr     cpsr_c, #(CPSR_MODE_FIQ)
        ldr     sp, =VA_STACK0_FIQ_BASE
        msr     cpsr_c, #(CPSR_MODE_SVC)
        ldr     sp, =VA_STACK0_SVC_BASE

        @mov     r0, #'P'                        @ PC should now be reading
        @bl      ar_asm_print_char               @ from DDR, not FLASH
        @mov     r0, #'C'
        @bl      ar_asm_print_char
        @mov     r0, #'='
        @bl      ar_asm_print_char
        @mov     r0, #'['
        @bl      ar_asm_print_char
        @mov     r0, pc
        @bl      ar_asm_print_hex
        @mov     r0, #']'
        @bl      ar_asm_print_char
        @mov     r0, #','
        @bl      ar_asm_print_char
        @mov     r0, #' '
        @bl      ar_asm_print_char

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

        
        @ =====================================================================
        @ Hic sunt dracones...
        @ =====================================================================
        b       main




        @ =====================================================================
        @ Cores 1-3 of both ARM boards block here, until core 0 of ARM0 sets up
        @ everything, jumps in main() and starts the waking up of all other
        @ cores, by sending something in our mailbox.
        @ =====================================================================
cores123567:
        lsl     r0, r3, #20                     @ r0 = core_id << 20
        ldr     r1, =0xFFCE0000
        add     r0, r0, r1                      @ r0 = ARS_MBX_ACCESS(core_id)
        ldr     r0, [r0]                        @ block on mailbox read

        lsl     r3, r3, #18                     @ r3 = core_id * 0x40000

        msr     cpsr_c, #(CPSR_MODE_MON)        @ Set up stack pointers for
        ldr     sp, =VA_STACK0_MON_BASE         @ all processor modes of
        add     sp, sp, r3                      @ this core (spaced at
        msr     cpsr_c, #(CPSR_MODE_ABT)        @ 0x40000 per core)
        ldr     sp, =VA_STACK0_ABT_BASE
        add     sp, sp, r3
        msr     cpsr_c, #(CPSR_MODE_UND)
        ldr     sp, =VA_STACK0_UND_BASE
        add     sp, sp, r3
        msr     cpsr_c, #(CPSR_MODE_IRQ)
        ldr     sp, =VA_STACK0_IRQ_BASE
        add     sp, sp, r3
        msr     cpsr_c, #(CPSR_MODE_FIQ)
        ldr     sp, =VA_STACK0_FIQ_BASE
        add     sp, sp, r3
        msr     cpsr_c, #(CPSR_MODE_SVC)
        ldr     sp, =VA_STACK0_SVC_BASE
        add     sp, sp, r3

        b       main                            @ Etiam hic sunt dracones...


@ =============================================================================
@ ar_asm_print_char()  Prints a character to the debug serial port
@ =============================================================================
@ Args:     r0: character (input)
@ Clobbers: r0, r1
@ =============================================================================
ar_asm_print_char:
        ldr     r1, =UART0_BASE                 @ UART0 base adr (also DR)

        str     r0, [r1, #0]                    @ write byte to serial

print_char_loop:
        ldr     r0, [r1, #0x18]                 @ r0 = flag register
        tst     r0, #0x20                       @ get TXFF bit
        bne     print_char_loop

        bx      lr                              @ return


@ =============================================================================
@ ar_asm_print_str()   Prints a string to the debug serial port. Cannot use a 
@                      stack, we don't have one :)
@ =============================================================================
@ Args:     r0: null-terminated string to be printed (input)
@ Clobbers: r0, r1, r2
@ =============================================================================
ar_asm_print_str:
        ldr     r1, =UART0_BASE                 @ UART0 base adr (also DR)

print_str_outer:
        ldrb    r2, [r0]                        @ get byte from r0
        cmp     r2, #0                          @ if null, exit loop
        beq     print_str_return

        str     r2, [r1, #0]                    @ write byte to serial

print_str_inner:
        ldr     r2, [r1, #0x18]                 @ r0 = flag register
        tst     r2, #0x20                       @ get TXFF bit
        bne     print_str_inner

        add     r0, #1                          @ go to next byte
        b       print_str_outer                 @ loop until end of string

print_str_return:
        bx      lr                              @ return


@ =============================================================================
@ ar_asm_print_hex()   Prints a 32-bit value in hex to the debug serial port
@ =============================================================================
@ Args:     r0: 32-bit value (input)
@ Clobbers: r1, r2, r3
@ =============================================================================
ar_asm_print_hex:
        ldr     r1, =UART0_BASE                 @ UART0 base adr (also DR)
        mov     r2, #28                         @ shift position init
print_hex_outer:
        lsr     r3, r0, r2                      @ r3 = (r0 >> r2) & 0xf
        and     r3, #0xf

        add     r3, #0x30                       @ r3 += '0'
        cmp     r3, #0x39                       @ if in A-F range...
        ble     print_hex_val_ok
        add     r3, #0x7                        @ ... add 7 more.

print_hex_val_ok:
        str     r3, [r1, #0]                    @ write byte to serial

print_hex_inner:
        ldr     r3, [r1, #0x18]                 @ r0 = flag register
        tst     r3, #0x20                       @ get TXFF bit
        bne     print_hex_inner

        sub     r2, #4                          @ go to next shift amount
        cmp     r2, #0                          @ if finished, exit
        blt     print_hex_return
        b       print_hex_outer                 @ loop until end of shift

print_hex_return:
        bx      lr                              @ return



