!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: dlaed6_init_mod
!
! Copyright(C) 2012-2018 RIKEN.
! Copyright(C) 2011-2012 Toshiyuki Imamura
!                        Graduate School of Informatics and Engineering,
!                        The University of Electro-Communications.
! Copyright (C) 2011- 2015 Japan Atomic Energy Agency.
! 
! Redistribution  and  use  in  source and binary forms, with or without
! modification,  are  permitted  provided  that the following conditions
! are met:
! 
! * Redistributions  of  source  code  must  retain  the above copyright
!   notice,  this  list  of  conditions  and  the  following  disclaimer.
! * Redistributions  in  binary  form must reproduce the above copyright
!   notice,  this list of conditions and the following disclaimer in the
!   documentation  and/or other materials provided with the distribution.
! 
! THIS  SOFTWARE  IS  PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
! ``AS IS''  AND  ANY  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
! LIMITED  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
! A  PARTICULAR  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
! HOLDERS OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
! SPECIAL,  EXEMPLARY,  OR  CONSEQUENTIAL  DAMAGES  (INCLUDING,  BUT NOT
! LIMITED  TO,  PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
! DATA,  OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
! THEORY  OF  LIABILITY,  WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
! (INCLUDING  NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF  THIS  SOFTWARE,  EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      MODULE DLAED6_INIT_MOD

      USE DLAED6_COMMON_MOD
      IMPLICIT NONE
      PUBLIC

      CONTAINS

      SUBROUTINE DLAED6_INIT()
!     
!     
!     =====================================================================
!     
!     .. Parameters ..
      REAL(8), PARAMETER  :: ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     &     THREE = 3.0D0, FOUR = 4.0D0, EIGHT = 8.0D0
!     ..
!     .. External Functions ..
      REAL(8), EXTERNAL   :: DLAMCH
!     ..
!     .. Local Scalars ..
      REAL(8)             :: BASE
!     ..
!     .. Intrinsic Functions ..
      INTRINSIC           :: INT, LOG
!     ..
!     .. Executable Statements ..
!     
!$OMP CRITICAL
      IF( FIRST ) THEN
         EPS = DLAMCH( 'Epsilon' )
         BASE = DLAMCH( 'Base' )
         SMALL1 = BASE**( INT( LOG( DLAMCH( 'SafMin' ) ) / LOG( BASE ) /
     &        THREE ) )
         SMINV1 = ONE / SMALL1
         SMALL2 = SMALL1*SMALL1
         SMINV2 = SMINV1*SMINV1
         FIRST = .FALSE.
      END IF
!$OMP END CRITICAL
!     
      RETURN
!     
!     End of DLAED6_INIT
!     
      END SUBROUTINE DLAED6_INIT

      END MODULE DLAED6_INIT_MOD
