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
@ Abstract      : Cache-coherent lock functions, based on 
@                 load-link/store-conditional primitives
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: lock.s,v $
@ CVS revision  : $Revision: 1.1 $
@ Last modified : $Date: 2012/03/26 13:14:24 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================

.text
.align 4

.global ar_lock_get
.global ar_lock_release



@ ===========================================================================
@ ar_lock_get()                 Tries to get a lock. Does not block, it
@                               simply returns if the lock is not acquired.
@
@                               The core id is written in the lock variable,
@                               if the lock acquisition is successful.
@ ===========================================================================
@ * INPUTS
@   void *lock_adr      [r0]    Address of lock variable
@   int core_id         [r1]    Core ID of the caller
@
@ * RETURN VALUE
@   int                 [r0]    0: lock acquired
@                               > 0: lock not acquired
@ ===========================================================================
ar_lock_get:
        ldrex   r2, [r0]                        @ Load-link the lock

        cmp     r2, #-1                         @ If lock != -1, we failed
        bne     get_fail1

        strex   r2, r1, [r0]                    @ Store-conditional our core
                                                @ ID in the lock

        cmp     r2, #0                          @ Check if store-conditional
        bne     get_fail2                       @ failed

        dmb                                     @ Memory barrier

        mov     r0, #0                          @ Success
        bx      lr

get_fail1:
        mov     r0, #1
        bx      lr

get_fail2:
        mov     r0, #2
        bx      lr


@ ===========================================================================
@ ar_lock_release()             Releases a lock, by writing -1 to the
@                               variable. 
@ ===========================================================================
@ * INPUTS
@   void *lock_adr      [r0]    Address of lock variable
@
@ * RETURN VALUE
@   int                 [r0]    0: release successful
@                               1: release failed
@ ===========================================================================
ar_lock_release:
        ldrex   r2, [r0]                        @ Load-link the lock

        mov     r1, #-1                         @ Store-conditional -1 
        strex   r2, r1, [r0]                    @ in the lock

        cmp     r2, #0                          @ Check if store-conditional
        bne     rel_fail                        @ failed

        dmb                                     @ Memory barrier

        mov     r0, #0                          @ Success
        bx      lr

rel_fail:
        mov     r0, #1
        bx      lr



