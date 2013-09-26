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
@ Abstract      : Myrmics logo with auto-generated build date
@
@ =============================[ CVS Variables ]=============================
@ 
@ File name     : $RCSfile: version_myrmics.s,v $
@ CVS revision  : $Revision: 1.1 $
@ Last modified : $Date: 2012/03/28 16:25:55 $
@ Last author   : $Author: lyberis-spree $
@ 
@ ===========================================================================

.data
.global msg_greet

msg_greet: .asciz "                                       _aa/\r\n                                     _m?)?Qf\r\n                                     Q'\r\n                      aa    aaa.    ]f            aaaP??\"\"QQ\r\n                    aQQQQa_QQQWQ6.  ]         _aP?'       )?\r\n        _aaaaaaaaaajQQQQQQQQQQQQQQc ]       aJ'\r\n     .wQQQQWQWQWQWQQQQQQQQQQQQQQQQQQQQmaaamQP\r\n    _QQQQQQQQQQQQQQQQQQQQPQQQQQQQQQQQQ4Q?4QQQQ)/\r\n   _QQQQQQQQQQQQQQQQQQQQQ  )?4QQQQQQQ$QQQ/)4QQ[\\\r\n   jQQQQQQQQQQQQQQQQQQQQ6   _aQP' jQQQQQQf  )W ]\r\n  _QQQQQQQQQQQQQQQQQP'Q)Wa/]QP'   f4QJQQQ    Qf)\r\n  ]QQQQQQQQQQQQQP??46 Q )?4QQaa   P_?QQQQQ6,_QQ7\r\n  ]QQQQQQQQQQQ6     )[]6.   QQP46]/?a'PQQQQQQQQf\r\n  ]QQQQQQQQQQ 4Q/   _f 4Q/   ?4QQ6jP??')?4QQQQQgaaa\r\n   ?QQQQQQP\\'  )4/ _P   \"4      Q-?4Q,    )WQQQf  ??w,\r\n     )???'jP    _f_D     j      Q   ]f     )6        )\r\n         jP     J_Q'    w'     ]Q   ]f      )6\r\n        ]Q,   _m']Q    j[      ]Q   ]f       )`\r\n        -4Qa _Q' )4Qa/jP       jQ   mf\r\n           ??QE      9Qf       QQa  Qf                 Myrmics rev BUILD_VER\r\n             9Qa      ?Q6/     )???_Qf                 build BUILD_DATE\r\n               ??\\       ??        ]Q6                 (c) FORTH-ICS/CARV\r\n\nThis is Myrmics.\r\n\n"


