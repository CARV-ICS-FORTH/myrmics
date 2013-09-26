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
@ Author        : Spyros Lyberis / ARM reference DVD
@ Abstract      : Routines to perform initialization of basic ARM hardware
@                 busses. Code based on ARM Reference DVD for the Versatile
@                 Express boards, as noted per each related routine below.
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: bus.s,v $
@ CVS revision  : $Revision: 1.1 $
@ Last modified : $Date: 2012/03/15 15:31:47 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================

.text
.align 4
.global ar_bus_smc
.global ar_bus_dmc
.global ar_bus_axi

@ =============================================================================
@ Definitions
@ =============================================================================

@ SMC
.equ SMC_BASE,          0x100E1000

.equ SMC_DIRECT_CMD,    0x10
.equ SMC_SET_CYCLES,    0x14
.equ SMC_SET_OPMODE,    0x18


@ DMC
.equ DMC_BASE,                 0x100E0000  @ DMC base address

.equ DMC_STATUS,               0x00        @ Status register
.equ DMC_COMMAND,              0x04        @ Command register
.equ DMC_DIRECT_CMD,           0x08        @ Direct Command register
.equ DMC_MEMORY_CONFIG,        0x0C        @ Memory Configuration register

.equ DMC_REFRESH_PRD,          0x10        @ Refresh Period register
.equ DMC_CAS_LATENCY,          0x14        @ CAS Latency Register
.equ DMC_WRITE_LATENCY,        0x18        @ t_dqss register
.equ DMC_T_MRD,                0x1C        @ t_mrd register

.equ DMC_T_RAS,                0x20        @ t_ras register
.equ DMC_T_RC,                 0x24        @ t_rc register
.equ DMC_T_RCD,                0x28        @ t_rcd register
.equ DMC_T_RFC,                0x2C        @ t_rfc register

.equ DMC_T_RP,                 0x30        @ t_rp register
.equ DMC_T_RRD,                0x34        @ t_rrd register
.equ DMC_T_WR,                 0x38        @ t_wr register
.equ DMC_T_WTR,                0x3C        @ t_wtr register

.equ DMC_T_XP,                 0x40        @ t_xp register
.equ DMC_T_XSR,                0x44        @ t_xsr register
.equ DMC_T_ESR,                0x48        @ t_esr register

.equ DMC_MEMORY_CFG2,          0x4C        @ memory_cfg2 register
.equ DMC_MEMORY_CFG3,          0x50        @ memory_cfg3 register
.equ DMC_T_FAW,                0x54        @ t_faw register

.equ DMC_ID_0_CFG,             0x100       @ id_cfg registers
.equ DMC_ID_1_CFG,             0x104       @ id_cfg registers
.equ DMC_ID_2_CFG,             0x108       @ id_cfg registers
.equ DMC_ID_3_CFG,             0x10c       @ id_cfg registers
.equ DMC_ID_4_CFG,             0x110       @ id_cfg registers
.equ DMC_ID_5_CFG,             0x114       @ id_cfg registers
.equ DMC_ID_6_CFG,             0x118       @ id_cfg registers
.equ DMC_ID_7_CFG,             0x11c       @ id_cfg registers
.equ DMC_ID_8_CFG,             0x120       @ id_cfg registers
.equ DMC_ID_9_CFG,             0x124       @ id_cfg registers
.equ DMC_ID_10_CFG,            0x128       @ id_cfg registers
.equ DMC_ID_11_CFG,            0x12c       @ id_cfg registers
.equ DMC_ID_12_CFG,            0x130       @ id_cfg registers
.equ DMC_ID_13_CFG,            0x134       @ id_cfg registers
.equ DMC_ID_14_CFG,            0x138       @ id_cfg registers
.equ DMC_ID_15_CFG,            0x13c       @ id_cfg registers

.equ DMC_CHIP_0_CFG,           0x200       @ chip_cfg 0 registers
.equ DMC_CHIP_1_CFG,           0x204       @ chip_cfg 1 registers

.equ DMC_USER_STATUS,          0x300       @ user_status register (RO)

.equ DMC_USER_0_CFG,           0x304       @ user_cfg0 registers (WO)
.equ DMC_USER_1_CFG,           0x308       @ user_cfg1 registers (WO)
.equ DMC_FEATURE_CRTL,         0x30C       @ Feature Control (RW)
.equ DMC_USER_2_CFG,           0x310       @ user_cfg registers (RO)

.equ DMC_USR2_BIGMEM_MASK,     0x40000000
.equ DMC_USR2_BIGMEM_SHIFT,    30
.equ DMC_USR2_UODTSW_MASK,     0x01800000
.equ DMC_USR2_UODTSW_SHIFT,    23

.equ TC_UIOLHNC_MASK,          0x000003C0
.equ TC_UIOLHNC_SHIFT,         0x6
.equ TC_UIOLHPC_MASK,          0x0000003F
.equ TC_UIOLHPC_SHIFT,         0x2
.equ TC_UIOHOCT_MASK,          0x2
.equ TC_UIOHOCT_SHIFT,         0x1

.equ TC_UIOHSTOP_SHIFT,        0x0
.equ TC_UIOLHXC_VALUE,         0x4

.equ DDR_EMR_OCD_MASK,         0x0000380
.equ DDR_EMR_OCD_SHIFT,        0x7
.equ DDR_EMR_RTT_MASK,         0x00000044  @ DDR2 Device RTT (ODT) settings
.equ DDR_EMR_RTT_SHIFT,        0x2     
.equ DDR_EMR_ODS_MASK,         0x00000002  @ DDR2 Output Drive Strength
.equ DDR_EMR_ODS_SHIFT,        0x0001

.equ DDR_EMR_RTT_50,           0x00000044  @ DDR2 50 Ohm termination
.equ DDR_EMR_RTT_75R,          0x00000004  @ DDR2 75 Ohm termination
.equ DDR_EMR_RTT_150,          0x00000040  @ DDR2 150 Ohm termination
.equ DDR_EMR_ODS_FULL,         0x0         @ DDR2 Full Drive Strength
.equ DDR_EMR_ODS_HALF,         0x1         @ DDR2 Half Drive Strength
.equ DDR_EMR_OCD_DEFAULT,      0x7
.equ DDR_EMR_OCD_NS,           0x0

.equ DDR_EMR_RTT_VAL,          DDR_EMR_RTT_75R         
.equ DDR_EMR_ODS_VAL,          DDR_EMR_ODS_FULL


@ AXI Bus
.equ PL301_FAXI_BASE,               0x100E9000


.equ PL301_QOS_TIDEMARK_MI_0,       0x400    @ QoS Tidemark for MI 0
.equ PL301_QOS_ACCESSCONTROL_MI_0,  0x404    @ QoS Tidemark for MI 0

.equ PL301_QOS_TIDEMARK_MI_1,       0x420    @ QoS Tidemark for MI 0
.equ PL301_QOS_ACCESSCONTROL_MI_1,  0x424    @ QoS Tidemark for MI 0

.equ PL301_QOS_TIDEMARK_MI_2,       0x440    @ QoS Tidemark for MI 0
.equ PL301_QOS_ACCESSCONTROL_MI_2,  0x444    @ QoS Tidemark for MI 0

.equ PL301_AR_ARB_MI_0,             0x408    @ QoS Tidemark for MI 0
.equ PL301_AW_ARB_MI_0,             0x40C    @ QoS Tidemark for MI 0

.equ PL301_AR_ARB_MI_1,             0x428    @ QoS Tidemark for MI 0
.equ PL301_AW_ARB_MI_1,             0x42C    @ QoS Tidemark for MI 0

.equ PL301_AR_ARB_MI_2,             0x448    @ QoS Tidemark for MI 0
.equ PL301_AW_ARB_MI_2,             0x44C    @ QoS Tidemark for MI 0

.equ PL301_MI_1_OFFSET,             0x20     @ Offset for MI 1 from MI 1
.equ PL301_MI_2_OFFSET,             0x40     @ Offset for MI 1 from MI 2
.equ PL301_MI_3_OFFSET,             0x60     @ Offset for MI 1 from MI 3
.equ PL301_MI_4_OFFSET,             0x80     @ Offset for MI 1 from MI 4
.equ PL301_MI_5_OFFSET,             0xa0     @ Offset for MI 1 from MI 5


.equ V2P_CA9_FAXI_MI0_TIDEMARK_VAL, 0x6      @ QoS takes control after 2 outstanding transactions
.equ V2P_CA9_FAXI_MI0_ACCESSCNTRL_VAL, 0x1   @ Give spare slots to Slave 0 = CLCD

.equ V2P_CA9_FAXI_MI1_TIDEMARK_VAL, 0x6      @ QoS takes control after 2 outstanding transactions
.equ V2P_CA9_FAXI_MI1_ACCESSCNTRL_VAL, 0x1   @ Give spare slots to Slave 0 = CLCD

.equ V2P_CA9_FAXI_MI2_TIDEMARK_VAL,    0x6   @ QoS takes control after 2 outstanding transactions
.equ V2P_CA9_FAXI_MI2_ACCESSCNTRL_VAL, 0x1   @ Give spare slots to Slave 0 = CLCD


@ =============================================================================
@ Initialize Static Memory Bus. Code adapted from ARM Boot Monitor code
@ located in Boot_Monitor/Platform/Source/sys_smc.s
@ =============================================================================

ar_bus_smc:
        ldr     r0, = SMC_BASE                 @ r0 = SMC base adr

@ cs0-0 setup (nor0)
        ldr     r1, = 0x0002393A
        str     r1, [r0, #SMC_SET_CYCLES]

        ldr     r1, = 0x00000AAA
        str     r1, [r0, #SMC_SET_OPMODE]

        ldr     r1, = 0x00400000
        str     r1, [r0, #SMC_DIRECT_CMD]

@ cs0-2 cycles (psram)
        ldr     r1, = 0x00027158
        str     r1, [r0, #SMC_SET_CYCLES]

        ldr     r1, = 0x00000802
        str     r1, [r0, #SMC_SET_OPMODE]

        ldr     r1, = 0x01400000
        str     r1, [r0, #SMC_DIRECT_CMD]

@ cs0-3 cycles (usb/eth/vram)
        ldr     r1, = 0x000CD2AA
        str     r1, [r0, #SMC_SET_CYCLES]

        ldr     r1, = 0x00000046
        str     r1, [r0, #SMC_SET_OPMODE]

        ldr     r1, = 0x01C00000
        str     r1, [r0, #SMC_DIRECT_CMD]

@ cs1-0  cycles (nor1)
        ldr     r1, = 0x0002393A
        str     r1, [r0, #SMC_SET_CYCLES]

        ldr     r1, = 0x00000AAA
        str     r1, [r0, #SMC_SET_OPMODE]

        ldr     r1, = 0x02400000
        str     r1, [r0, #SMC_DIRECT_CMD]

@ cs1-3 cycles (peripherals)
        ldr     r1, = 0x00025156
        str     r1, [r0, #SMC_SET_CYCLES]

        ldr     r1, = 0x00000046
        str     r1, [r0, #SMC_SET_OPMODE]

        ldr     r1, = 0x03C00000
        str     r1, [r0, #SMC_DIRECT_CMD]

@ cs0-1 VRAM

        ldr     r1, =  0x00049249
        str     r1, [r0, #SMC_SET_CYCLES]

        ldr     r1, = 0x00000046
        str     r1, [r0, #SMC_SET_OPMODE]

        ldr     r1, = 0x00C00000
        str     r1, [r0, #SMC_DIRECT_CMD]

@ page mode setup for VRAM
        ldr     r0, = 0x4CFFFFFC

@ read current state 
        ldr     r1, [r0, #0]
        ldr     r1, [r0, #0]
        ldr     r1, = 0x00000000
        str     r1, [r0, #0]
        ldr     r1, [r0, #0]

@ enable page mode 
        ldr     r1, [r0, #0]
        ldr     r1, [r0, #0]
        ldr     r1, = 0x00000000
        str     r1, [r0, #0]
        ldr     r1, = 0x00900090
        str     r1, [r0, #0]

@ confirm page mode enabled
        ldr     r1, [r0, #0]
        ldr     r1, [r0, #0]
        ldr     r1, = 0x00000000
        str     r1, [r0, #0]
        ldr     r1, [r0, #0]

@ return
        bx      lr



@ =============================================================================
@ Initialize Dynamic Memory Controller. Code adapted from ARM Boot Monitor code
@ located in Boot_Monitor/Platform/Source/sys_dmc_v2p_ca9.s
@ =============================================================================
ar_bus_dmc:

        LDR     r0, =DMC_BASE
    
@ set config mode
        MOV     r1, #0x4
        STR     r1, [r0, #DMC_COMMAND]
    
@ Setup the QoS AXI ID bits    
        
        @ CLCD AXIID = 000
        LDR     r1, =0x3
        STR     r1, [r0, #DMC_ID_0_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_1_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_2_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_3_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_4_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_5_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_6_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_7_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_8_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_9_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_10_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_11_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_12_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_13_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_14_CFG]

        @ Default disable QoS
        LDR     r1, =0x0
        STR     r1, [r0, #DMC_ID_15_CFG]


@ initialise memory controlller

        @ refresh period
        LDR     r1, =0x3D0
        STR     r1, [r0, #DMC_REFRESH_PRD]

        @ cas latency
        MOV     r1, #0x8 
        STR     r1, [r0, #DMC_CAS_LATENCY]

        @ write latency
        MOV     r1, #0x3
        STR     r1, [r0, #DMC_WRITE_LATENCY]

        @ t_mrd
        MOV     r1, #0x2
        STR     r1, [r0, #DMC_T_MRD]

        @ t_ras
        MOV     r1, #0xa 
        STR     r1, [r0, #DMC_T_RAS]

        @ t_rc
        MOV     r1, #0xE 
        STR     r1, [r0, #DMC_T_RC]

        @ t_rcd
        LDR     r1, =0x104
        STR     r1, [r0,#DMC_T_RCD]

        @ t_rfc
        LDR     r1, =0x2f32     
        STR     r1, [r0, #DMC_T_RFC]
    
        @ t_rp
        LDR     r1, =0x14
        STR     r1, [r0, #DMC_T_RP]

        @ t_rrd
        MOV     r1, #0x2
        STR     r1, [r0, #DMC_T_RRD]

        @ t_wr
        MOV     r1, #0x4
        STR     r1, [r0, #DMC_T_WR]

        @ t_wtr
        MOV     r1, #0x2
        STR     r1, [r0, #DMC_T_WTR]

        @ t_xp
        MOV     r1, #0x2
        STR     r1, [r0, #DMC_T_XP]

        @ t_xsr
        MOV     r1, #0xC8
        STR     r1, [r0, #DMC_T_XSR]

        @ t_esr
        MOV     r1, #0x14
        STR     r1, [r0, #DMC_T_ESR]

        @ t_faw
        LDR     r1, =0x00000a0d 
        STR     r1, [r0, #DMC_T_FAW]

@ =======================================================================
@ Initialise PL341 Mem Config Registers
@ =======================================================================

        @ |======================================
        @ | Set PL341 Memory Config
        @ |======================================
        @ | ACLK async to MCLK
        @ | dqm & cke = 1
        @ |======================================
        LDR     r1, =0x00010022         
        STR     r1, [r0, #DMC_MEMORY_CONFIG]

        @ |======================================
        @ | Set PL341 Memory Config 2
        @ |======================================
        @ | ACLK async to MCLK
        @ | dqm & cke = 1
        @ |======================================
        LDR     r1, =0x0000007C   
        STR     r1, [r0, #DMC_MEMORY_CFG2]

        @ |======================================
        @ | Set PL341 Chip Select 0 
        @ |======================================
        LDR     r1, =0x00010000                         
        STR     r1, [r0, #DMC_CHIP_0_CFG]

        @ |======================================
        @ | Set PL341 Memory Config 3 
        @ |======================================
        @ | Default value used
        LDR     r1, =0x00000001
        STR     r1, [r0, #DMC_MEMORY_CFG3]

        @ |========================================================
        @ |Set Test Chip PHY Registers via PL341 User Config Reg
        @ |Note that user_cfgX registers are Write Only
        @ |
        @ |DLL Freq set = 250MHz - 266MHz
        @ |======================================================== 
        LDR     r1, =0x7C924924         
        STR     r1, [r0, #DMC_USER_0_CFG]
 
        @ Default values of user_cfg1 are most suitable
        ;LDR     r1, =0x008A28A2         
        ;STR     r1, [r0, #DMC_USER_1_CFG]
 
        @ user_config2
        @ ------------
        @ Set defaults before calibrating the DDR2 buffer impendence
        @ -Disable ODT
        @ -Default drive strengths

        LDR     r1, =0x40000198          
        STR     r1, [r0, #DMC_USER_2_CFG]
 
        @ |=======================================================
        @ |Auto calibrate the DDR2 buffers impendence 
        @ |=======================================================
label_99:             
        LDR     r3, [r0, #DMC_USER_STATUS]    @ Read PL341 USER STATUS register.
        TST     r3, #0x100                    @ Has task completed?
        BEQ     label_99
        LDR     r1, =0x40800199               @ Now stop the auto adjustment
        STR     r1, [r0, #DMC_USER_2_CFG]

        @ Set the output driven strength
        @LDR     r1, =0x40800000:OR:((TC_UIOLHXC_VALUE:SHL:TC_UIOLHNC_SHIFT):OR:(TC_UIOLHXC_VALUE:SHL:TC_UIOLHPC_SHIFT):OR:(0x1:SHL:TC_UIOHOCT_SHIFT):OR:(0x1:SHL:TC_UIOHSTOP_SHIFT)) 
        LDR     r1, =(0x40800000|((TC_UIOLHXC_VALUE << TC_UIOLHNC_SHIFT)|(TC_UIOLHXC_VALUE << TC_UIOLHPC_SHIFT)|(0x1 << TC_UIOHOCT_SHIFT)|(0x1 << TC_UIOHSTOP_SHIFT)))
        STR     r1, [r0, #DMC_USER_2_CFG]

        @ |======================================
        @ | Set PL341 Feature Control Register 
        @ |======================================
        @ | Disable early BRESP - use to optimise CLCD performance
        LDR     r1, =0x00000001
        STR     r1, [r0, #DMC_FEATURE_CRTL]
 
@=================
@ Config memories
@=================
        
        @ send nop
        LDR     r1, =0x000C0000
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ pre-charge all
        MOV     r1, #0x0
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ delay
        MOV     r1, #0             
label_1:
        LDR     r3, [r0, #DMC_STATUS]    @ read status register
        ADD     r1, r1, #1
        CMP     r1, #10
        BLT     label_1

        @ set (EMR2) extended mode register 2
        LDR     r1, =0x000A0000
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ set (EMR3) extended mode register 3
        LDR     r1, =0x000B0000
        STR     r1, [r0, #DMC_DIRECT_CMD]   

        @ =================================
        @  set (EMR) Extended Mode Register
        @ ==================================
        @ Put into OCD default state
        LDR     r1, =0x00090000
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ ===========================================================        
        @ set (MR) mode register - With DLL reset
        @ ===========================================================
        @ Burst Length = 4 (010)
        @ Burst Type   = Seq (0)
        @ Latency      = 4 (100)
        @ Test mode    = Off (0)
        @ DLL reset    = Yes (1)
        @ Wr Recovery  = 4  (011)      
        @ PD           = Normal (0)
        LDR     r1, =0x00080742                                 
        STR     r1, [r0, #DMC_DIRECT_CMD]
        
        @ pre-charge all 
        MOV     r1, #0x0
        STR     r1, [r0, #DMC_DIRECT_CMD]
        
        @ auto-refresh 
        LDR     r1, =0x00040000
        STR     r1, [r0, #DMC_DIRECT_CMD]
        
        @ delay
        MOV     r1, #0                
label_2:
        LDR     r3, [r0, #DMC_STATUS]    @ read status register
        ADD     r1, r1, #1
        CMP     r1, #10
        BLT     label_2

        @ auto-refresh
        LDR     r1, =0x00040000
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ ===========================================================        
        @ set (MR) mode register - Without DLL reset
        @ ===========================================================
        LDR     r1, =0x00080642
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ delay
        MOV     r1, #0             
label_3:
        LDR     r3, [r0, #DMC_STATUS]    @ read status register
        ADD     r1, r1, #1
        CMP     r1, #10
        BLT     label_3

        @ ======================================================        
        @ set (EMR) extended mode register - Enable OCD defaults
        @ ====================================================== 
        NOP
        
        @LDR     r1, =0x00090000 :OR:((DDR_EMR_OCD_DEFAULT:SHL:DDR_EMR_OCD_SHIFT):OR:(DDR_EMR_RTT_VAL)(DDR_EMR_ODS_VAL:SHL:DDR_EMR_ODS_MASK))         
        LDR     r1, =(0x00090000 | ((DDR_EMR_OCD_DEFAULT << DDR_EMR_OCD_SHIFT ) | (DDR_EMR_RTT_VAL) | (DDR_EMR_ODS_VAL << DDR_EMR_ODS_MASK)))
        STR     r1, [r0, #DMC_DIRECT_CMD]

        @ delay
        MOV     r1, #0             
label_4:
        LDR     r3, [r0, #DMC_STATUS]    @ read status register
        ADD     r1, r1, #1
        CMP     r1, #10
        BLT     label_4
    
        @ Set (EMR) extended mode register - OCD Exit
        NOP
        @LDR     r1, =0x00090000 :OR:((DDR_EMR_OCD_NS:SHL:DDR_EMR_OCD_SHIFT) \ :OR:(DDR_EMR_RTT_VAL) \ :OR:(DDR_EMR_ODS_VAL:SHL:DDR_EMR_ODS_MASK))         
        LDR     r1, =(0x00090000 | ((DDR_EMR_OCD_NS << DDR_EMR_OCD_SHIFT) | (DDR_EMR_RTT_VAL) | (DDR_EMR_ODS_VAL << DDR_EMR_ODS_MASK)) )
        STR     r1, [r0, #DMC_DIRECT_CMD]
        
        @ delay
        MOV     r1, #0             
label_5:
        LDR     r3, [r0, #DMC_STATUS]    @ read status register
        ADD     r1, r1, #1
        CMP     r1, #10
        BLT     label_5
    
@----------------------------------------    
@ go command
        MOV     r1, #0x0
        STR     r1, [r0, #DMC_COMMAND]

@ wait for ready
label_31:
        LDR     r1, [r0,#DMC_STATUS]
        TST     r1,#1
        BEQ     label_31
        
@ return
        bx      lr


@ =============================================================================
@ Initialize PL301 on-chip AXI Bus. Code adapted from ARM Boot Monitor code
@ located in Boot_Monitor/Platform/Source/sys_bus.s
@ =============================================================================

ar_bus_axi:
        @ ======================================
        @ Enable the PL301 Fast AXI QoS settings
        @ ======================================
                        
        @ Configure Tidemark Register for Master Port 0 (MI 0)
        LDR             r0, =PL301_FAXI_BASE
        LDR     r1, =V2P_CA9_FAXI_MI0_TIDEMARK_VAL
        STR     r1, [r0, #PL301_QOS_TIDEMARK_MI_0]

        @ Configure the Access Control Register (MI 0)
        LDR     r1, =V2P_CA9_FAXI_MI0_ACCESSCNTRL_VAL
        STR     r1, [r0, #PL301_QOS_ACCESSCONTROL_MI_0]

        @ Configure Tidemark Register for Master Port 0 (MI 0)
        ;LDR            r0, =PL301_FAXI_BASE
        ;LDR     r1, =V2P_CA9_FAXI_MI1_TIDEMARK_VAL
        ;STR     r1, [r0, #PL301_QOS_TIDEMARK_MI_1]

        @ Configure the Access Control Register (MI 0)
        ;LDR     r1, =V2P_CA9_FAXI_MI1_ACCESSCNTRL_VAL
        ;STR     r1, [r0, #PL301_QOS_ACCESSCONTROL_MI_1]

        @ Configure Tidemark Register for Master Port 0 (MI 0)
        ;LDR            r0, =PL301_FAXI_BASE
        ;LDR     r1, =V2P_CA9_FAXI_MI2_TIDEMARK_VAL
        ;STR     r1, [r0, #PL301_QOS_TIDEMARK_MI_2]

        @ Configure the Access Control Register (MI 0)
        ;LDR     r1, =V2P_CA9_FAXI_MI2_ACCESSCNTRL_VAL
        ;STR     r1, [r0, #PL301_QOS_ACCESSCONTROL_MI_2]

        @ Configure Arbitration Priority for Master Port 0 
        @ We need to set the CLCD master port (SP 0) to 0 and all the others to 1
        @ CLCD will then have highest priority (i.e. lowest priority value)
        LDR     r0, =PL301_FAXI_BASE
        
        @ MP0 
        @ Set priority for Read
        LDR     r1, = 0x00000100
        STR     r1, [r0,#PL301_AR_ARB_MI_0]             @ Set SP0 to priority 0   
        
        LDR     r1, = 0x01000200
        STR     r1, [r0,#PL301_AR_ARB_MI_0]             @ Set SP1 to priority 1        
        
        LDR     r1, = 0x02000200
        STR     r1, [r0,#PL301_AR_ARB_MI_0]             @ Set SP2 to priority 1        
        
        LDR     r1, = 0x03000200
        STR     r1, [r0,#PL301_AR_ARB_MI_0]             @ Set SP3 to priority 1        
        
        LDR     r1, = 0x04000200
        STR     r1, [r0,#PL301_AR_ARB_MI_0]             @ Set SP4 to priority 1 
        
        @ Set priority for Write
        LDR     r1, = 0x00000100
        STR     r1, [r0,#PL301_AW_ARB_MI_0]             @ Set SP0 to priority 0   
        
        LDR     r1, = 0x01000200
        STR     r1, [r0,#PL301_AW_ARB_MI_0]             @ Set SP1 to priority 1        
        
        LDR     r1, = 0x02000200
        STR     r1, [r0,#PL301_AW_ARB_MI_0]             @ Set SP2 to priority 1        
        
        LDR     r1, = 0x03000200
        STR     r1, [r0,#PL301_AW_ARB_MI_0]             @ Set SP3 to priority 1        
        
        LDR     r1, = 0x04000200
        STR     r1, [r0,#PL301_AW_ARB_MI_0]             @ Set SP4 to priority 1 

        @ MP1
        @ Set priority for Read
        LDR     r1, = 0x00000100
        STR     r1, [r0,#PL301_AR_ARB_MI_1]             @ Set SP0 to priority 0   
        
        LDR     r1, = 0x01000200
        STR     r1, [r0,#PL301_AR_ARB_MI_1]             @ Set SP1 to priority 1        
        
        LDR     r1, = 0x02000200
        STR     r1, [r0,#PL301_AR_ARB_MI_1]             @ Set SP2 to priority 1        
        
        LDR     r1, = 0x03000200
        STR     r1, [r0,#PL301_AR_ARB_MI_1]             @ Set SP3 to priority 1        
        
        LDR     r1, = 0x04000200
        STR     r1, [r0,#PL301_AR_ARB_MI_1]             @ Set SP4 to priority 1 

        @ Set priority for Write
        LDR     r1, = 0x00000100
        STR     r1, [r0,#PL301_AW_ARB_MI_1]             @ Set SP0 to priority 0   
        
        LDR     r1, = 0x01000200
        STR     r1, [r0,#PL301_AW_ARB_MI_1]             @ Set SP1 to priority 1        
        
        LDR     r1, = 0x02000200
        STR     r1, [r0,#PL301_AW_ARB_MI_1]             @ Set SP2 to priority 1        
        
        LDR     r1, = 0x03000200
        STR     r1, [r0,#PL301_AW_ARB_MI_1]             @ Set SP3 to priority 1        
        
        LDR     r1, = 0x04000200
        STR     r1, [r0,#PL301_AW_ARB_MI_1]             @ Set SP4 to priority 1 

        @ MP2
        @ Set priority for Read
        LDR     r1, = 0x00000100
        STR     r1, [r0,#PL301_AR_ARB_MI_2]             @ Set SP0 to priority 0   
        
        LDR     r1, = 0x01000100
        STR     r1, [r0,#PL301_AR_ARB_MI_2]             @ Set SP1 to priority 1        
        
        LDR     r1, = 0x02000100
        STR     r1, [r0,#PL301_AR_ARB_MI_2]             @ Set SP2 to priority 1        
        
        LDR     r1, = 0x03000100
        STR     r1, [r0,#PL301_AR_ARB_MI_2]             @ Set SP3 to priority 1        
        
        LDR     r1, = 0x04000100
        STR     r1, [r0,#PL301_AR_ARB_MI_2]             @ Set SP4 to priority 1 
        
        @ Set priority for Write
        LDR     r1, = 0x00000100
        STR     r1, [r0,#PL301_AW_ARB_MI_2]             @ Set SP0 to priority 0   
        
        LDR     r1, = 0x01000200
        STR     r1, [r0,#PL301_AW_ARB_MI_2]             @ Set SP1 to priority 1        
        
        LDR     r1, = 0x02000200
        STR     r1, [r0,#PL301_AW_ARB_MI_2]             @ Set SP2 to priority 1        
        
        LDR     r1, = 0x03000200
        STR     r1, [r0,#PL301_AW_ARB_MI_2]             @ Set SP3 to priority 1        
        
        LDR     r1, = 0x04000200
        STR     r1, [r0,#PL301_AW_ARB_MI_2]             @ Set SP4 to priority 1 
   
        BX      lr
                
