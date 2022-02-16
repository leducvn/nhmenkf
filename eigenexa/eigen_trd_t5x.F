!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t5x_mod
!
! Purpose
! =======
!
! eigen_trd_t5_mod manages the modules for the special update to the
! next target vector for creating a Householder reflector.
!
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

      module eigen_trd_t5x_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_local_2update0

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_local_2update0
!
! Purpose
! =======
!
! The special updater to the next target vector for creating a
! Householder reflector.
!
!
! Arguments
! =========
!
! uxx    (input/output) real(8) array, dimension(*)
!        uxx returns the next target reflector vector updated.
!
! ux     (input) real(8) array, dimension(*)
!        ux contains the reflector vector.
!
! vx     (input) real(8) array, dimension(*)
!        vx contains the counterpart vector to ux.
!
! m      (input) integer
!        Block factor
!
! i      (input) integer
!        The global loop counter
!
! u_n    (input) real(8)
!        u_n contains replicated elements of ux(i,i:i+m-1)
!
! v_n    (input) real(8)
!        v_n contains replicated elements of vx(i,i:i+m-1)
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_local_2update0(
     &     uxx, ux, vx,
     &     m, i, u_n, v_n)

      real(8), intent(inout) ::  uxx(*)
      real(8), intent(in)    ::  ux(*)
      real(8), intent(in)    ::  vx(*)
      integer, intent(in)    ::  m
      integer, intent(in)    ::  i
      real(8), intent(in)    ::  u_n
      real(8), intent(in)    ::  v_n

      integer                ::  j, k, L
      integer                ::  j_1, j_2, j_3

      include 'CSTAB.h'


      if (m < 1) return
     
      L = i - 1

      j_2 = eigen_loop_start(1, x_nnod,x_inod)
      j_3 = eigen_loop_end  (L, x_nnod,x_inod)

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
      do j_1=j_2,j_3
        uxx(j_1) = uxx(j_1)
     &           - (ux(j_1)*v_n + vx(j_1)*u_n)
      end do                    ! j_1

      return

      end subroutine eigen_trd_local_2update0

      end module eigen_trd_t5x_mod
