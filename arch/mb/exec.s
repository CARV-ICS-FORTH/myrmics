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
# Abstract      : Arbitrary code execution functionality
#
# =============================[ CVS Variables ]=============================
# 
# File name     : $RCSfile: exec.s,v $
# CVS revision  : $Revision: 1.2 $
# Last modified : $Date: 2012/07/13 12:18:44 $
# Last author   : $Author: lyberis-spree $
# 
# ===========================================================================

.text

.global ar_exec



# ===========================================================================
# ar_exec()                     Executes a function call, putting its
#                               arguments into the proper places (combination
#                               of registers and stack (if needed)).
#
#                               LIMITATIONS:
#                               - Does not support a return value
#                               - Has not been tested with weird combinations
#                                 like passing structs by-value etc.
# ===========================================================================
# * INPUTS
#   void (*func_adr)()  [r5]    Address of function to be executed
#   int num_args        [r6]    Number of arguments to pass to the function
#   void **args         [r7]    Array of arguments to give to the function
#   int *types          [r8]    Array of argument types. FIXME -- define or
#                               remove
# ===========================================================================
ar_exec:
        addi    r3, r6, 8                       # r3 = num_args + 8 regs
        bslli   r3, r3, 2                       # r3 *= 4

        swi     r22, r1, -4                     # save old regs on the top of
        swi     r23, r1, -8                     # our expanded stack frame
        swi     r24, r1, -12
        swi     r25, r1, -16
        swi     r26, r1, -20
        swi     r27, r1, -24
        swi     r28, r1, -28
        rsub    r1, r3, r1                      # increase stack size by r3
        swi     r15, r1, 0                      # save link register on bottom

        add     r22, r5, r0                     # Move r5-r8 to r22-r25,
        add     r23, r6, r0                     # because we want to overwrite
        add     r24, r7, r0                     # r5-r10
        add     r25, r8, r0

        add     r27, r3, r0                     # r27 = stack frame size


        # ===================================================================
        # Pass 6 first args directly to r5-r10
        # ===================================================================
        add     r28, r0, r0                     # r28 = arg loop counter

        cmp     r12, r28, r23                   # No args? call function
        beqi    r12, call

        lwi     r5, r24, 0                      # first arg -> r5
        addi    r28, r28, 1
        cmp     r12, r28, r23                   # 1 arg? call function
        beqi    r12, call

        lwi     r6, r24, 4                      # second arg -> r6
        addi    r28, r28, 1
        cmp     r12, r28, r23                   # 2 args? call function
        beqi    r12, call

        lwi     r7, r24, 8                      # third arg -> r7
        addi    r28, r28, 1
        cmp     r12, r28, r23                   # 3 args? call function
        beqi    r12, call

        lwi     r8, r24, 12                     # fourth arg -> r8
        addi    r28, r28, 1
        cmp     r12, r28, r23                   # 4 args? call function
        beqi    r12, call

        lwi     r9, r24, 16                     # fifth arg -> r9
        addi    r28, r28, 1
        cmp     r12, r28, r23                   # 5 args? call function
        beqi    r12, call

        lwi     r10, r24, 20                    # sixth arg -> r10
        addi    r28, r28, 1
        cmp     r12, r28, r23                   # 6 args? call function
        beqi    r12, call


        # ===================================================================
        # Enlarge stack for the remaining arguments and put them there.
        # Note that Microblaze GCC places arg #7 in r1 + 28,
        # arg #8 in r1 + 32, etc. although the first six are NOT placed
        # in (r1 + 0) ... (r1 + 24) but instead in regs r5-r10.
        # ===================================================================
loop:
        # Arg copy: r1[arg_nr * 4] = r24[arg_nr * 4]

        bslli   r11, r28, 2                     # r11 = arg_nr * 4
        lw      r26, r24, r11                   # r26 = r24[arg_nr * 4]

        addi    r11, r11, 4                     # r11 += 4
        sw      r26, r1, r11                    # r1[(arg_nr + 1) * 4] = r26

        addi    r28, r28, 1                     # we did one more arg
        cmp     r12, r28, r23                   # finished args? call function
        beqi    r12, call
        bri     loop

        # ===================================================================
        # Call the function
        # ===================================================================
call:
        brald   r15, r22                        # call func
        nop

        lwi     r15, r1, 0                      # restore link register

        add     r1, r1, r27                     # restore stack pointer

        lwi     r22, r1, -4                     # load back old regs
        lwi     r23, r1, -8
        lwi     r24, r1, -12
        lwi     r25, r1, -16
        lwi     r26, r1, -20
        lwi     r27, r1, -24
        lwi     r28, r1, -28

        rtsd    r15, 8                          # return
        nop
