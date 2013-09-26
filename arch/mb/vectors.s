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
# Abstract      : MicroBlaze vectors, starting at binary position 0x0
#
# =============================[ CVS Variables ]=============================
# 
# File name     : $RCSfile: vectors.s,v $
# CVS revision  : $Revision: 1.2 $
# Last modified : $Date: 2012/01/27 16:32:59 $
# Last author   : $Author: lyberis-spree $
# 
# ===========================================================================

.text

.global ar_vectors

ar_vectors:
    brai    ar_boot             # 0x00: reset exception entry
    brai    ar_abort_handler    # 0x08: user exception entry
    brai    ar_irq_handler      # 0x10: interrupt entry
    brai    ar_abort_handler    # 0x18: break entry
    brai    ar_abort_handler    # 0x20: hardware exception entry


