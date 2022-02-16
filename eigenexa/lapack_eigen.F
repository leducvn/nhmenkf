!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: lapack_eigen_mod
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

      MODULE LAPACK_EIGEN_MOD

      IMPLICIT NONE

      CONTAINS

      SUBROUTINE LAPACK_EIGEN2( N_GLOBAL, N, HBW, ID, D, E, LDE, Q,
     &     LDQ )
      INTEGER, intent(IN)    :: N_GLOBAL, N, HBW, ID, LDE, LDQ
      REAL(8), intent(INOUT) :: D(*), E(LDE,*), Q(LDQ,*)

      INTEGER            :: i, j, info, lwork, liwork
      REAL(8)            :: temp
      REAL(8), PARAMETER :: ZERO = 0D0
      REAL(8), POINTER   :: work(:)
      INTEGER, POINTER   :: iwork(:)
      CHARACTER*1        :: JOBU, UPLO


      JOBU = 'V'; UPLO= 'U'

      lwork  = -1
      liwork = -1

      call dsyevd(JOBU, UPLO, N, Q, LDQ, D(1),
     &     temp, lwork, i, liwork, info)

      lwork  = int(temp)
      liwork = i

      DO J=1,N
         DO I=1,N
            Q(I,J) = ZERO
         END DO
      END DO

      DO I=1,N
         Q(I,I) = D(I)
      END DO

      DO J=1,HBW
         DO I=1,N-J
            Q(I,I+J) = E(I,J)
            Q(I+J,I) = E(I,J)
         END DO
      END DO

      allocate ( work(lwork), iwork(liwork) )

      call dsyevd(JOBU, UPLO, N, Q, LDQ, D(1),
     &     work, lwork, iwork, liwork, info)

      deallocate ( work, iwork )


      END SUBROUTINE LAPACK_EIGEN2

      END MODULE LAPACK_EIGEN_MOD
