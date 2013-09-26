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
@ Abstract      : Arbitrary code execution functionality
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: exec.s,v $
@ CVS revision  : $Revision: 1.2 $
@ Last modified : $Date: 2012/07/13 12:18:44 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================

.text
.align 4

.global ar_exec



@ ===========================================================================
@ ar_exec()                     Executes a function call, putting its
@                               arguments into the proper places (combination
@                               of registers and stack (if needed)).
@
@                               LIMITATIONS:
@                               - Does not support a return value
@                               - Does not support doublewords (64 bit) or
@                                 bigger co-processor arguments
@                               - Has not been tested with weird combinations
@                                 like passing structs by-value etc.
@ ===========================================================================
@ * INPUTS
@   void (*func_adr)()  [r0]    Address of function to be executed
@   int num_args        [r1]    Number of arguments to pass to the function
@   void **args         [r2]    Array of arguments to give to the function
@   int *types          [r3]    Array of argument types. FIXME -- define or
@                               remove
@ ===========================================================================
ar_exec:
        push    {r4-r8, r10-r11, LR}            @ save old regs on stack

        
        mov     r4, r0                          @ Move r0-r3 to r4-r7,
        mov     r5, r1                          @ because we want to overwrite
        mov     r6, r2                          @ r0-r3
        mov     r7, r3

        @ ===================================================================
        @ Pass 4 first args directly to r0-r3
        @ ===================================================================
        mov     r8, #0                          @ r8 = arg loop counter
        mov     r10, #0                         @ stack grow size = 0

        cmp     r8, r5                          @ No args? call function
        beq     call

        ldr     r0, [r6, #0]                    @ first arg -> r0
        add     r8, r8, #1
        cmp     r8, r5                          @ 1 arg? call function
        beq     call

        ldr     r1, [r6, #4]                    @ second arg -> r1
        add     r8, r8, #1
        cmp     r8, r5                          @ 2 args? call function
        beq     call

        ldr     r2, [r6, #8]                    @ third arg -> r2
        add     r8, r8, #1
        cmp     r8, r5                          @ 3 args? call function
        beq     call

        ldr     r3, [r6, #12]                   @ fourth arg -> r3
        add     r8, r8, #1
        cmp     r8, r5                          @ 4 args? call function
        beq     call


        @ ===================================================================
        @ Enlarge stack for the remaining arguments and put them there
        @ ===================================================================
        sub     r10, r5, r8                     @ r10 = num of remaining args
        lsl     r10, r10, #2                    @ r10 *= 4
        sub     sp, sp, r10                     @ increase stack size by r10

loop:
        @ Arg copy: sp[(arg nr - 4) * 4] = r6[arg nr * 4]

        lsl     ip, r8, #2                      @ ip = arg_nr * 4
        ldr     ip, [r6, ip]                    @ ip = r6[arg_nr * 4]

        sub     r11, r8, #4                     @ r11 = arg_nr - 4
        lsl     r11, r11, #2                    @ r11 = (arg_nr -4) * 4
        str     ip, [sp, r11]                   @ sp[(arg_nr -4) * 4] = ip

        add     r8, r8, #1                      @ we did one more arg
        cmp     r8, r5                          @ finished args? call function
        beq     call
        b       loop

        @ ===================================================================
        @ Call the function
        @ ===================================================================
call:
        blx     r4                              @ call func

        add     sp, sp, r10                     @ restore stack pointer, if
                                                @ we messed with it

        pop     {r4-r8, r10-r11, PC}            @ restore registers and 
                                                @ return

