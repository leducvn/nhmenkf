!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t5_mod
!
! Purpose
! =======
!
! eigen_trd_t5_mod manages the modules for the local update.
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

      module eigen_trd_t5_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_local_2update

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_local_2update
!
! Purpose
! =======
!
! The local updater.
!
!
! Arguments
! =========
!
! u_x    (input/output) real(8) array, dimension(ldv,*)
!        u_x returns the reflector vectors updated.
!
! ldv    (input) integer
!        The leading dimension of the working arrays, u_x, u_y, ...
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
! u_nn   (input) real(8) array, dimension(*)
!        u_nn contains replicated elements of ux(i,i:i+m-1)
!
! v_nn   (input) real(8) array, dimension(*)
!        v_nn contains replicated elements of vx(i,i:i+m-1)
!
! local_size (input) integer
!        Thread ID
!
! local_rank (input) integer
!        The number of threads in the main team
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_local_2update(
     &     u_x, ldv,
     &     ux, vx,
     &     m, i, u_nn, v_nn,
     &     local_size, local_rank)

      real(8), intent(inout) ::  u_x(1:ldv,*)
      integer, intent(in)    ::  ldv
      real(8), intent(in)    ::  ux(*)
      real(8), intent(in)    ::  vx(*)
      integer, intent(in)    ::  m
      integer, intent(in)    ::  i
      real(8), intent(in)    ::  u_nn(*)
      real(8), intent(in)    ::  v_nn(*)
      integer, intent(in)    ::  local_size
      integer, intent(in)    ::  local_rank

      integer                ::  j, k, L
      integer                ::  j_1, j_2, j_3
      integer                ::  l_1
      integer                ::  jj_1, jj_2, jj_3
      integer                ::  LX

      include 'CSTAB.h'

      real(8)                :: u_n, v_n

      integer                ::  ll_size
      integer                ::  ll_rank


      if (m < 1) return

      ll_size = local_size - min(local_size-1,1)
      ll_rank = local_rank - min(local_size-1,1)
      if (ll_rank < 0) return
     
      LX  = 64
     
      L = i - 1
      jj_2 = eigen_loop_start(1, x_nnod,x_inod)
      jj_3 = eigen_loop_end  (L, x_nnod,x_inod)

      do jj_1=jj_2+LX*ll_rank,jj_3,LX*ll_size
        j_2 = jj_1; j_3 = min(jj_1+LX-1, jj_3)

        do l_1=1,m

          u_n = u_nn(l_1)
          v_n = v_nn(l_1)

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
!IBM* INDEPENDENT
!IBM* ASSERT(NODEPS)
          do j_1=j_2,j_3
            u_x(j_1, l_1) = u_x(j_1, l_1)
     &                    - (ux(j_1)*v_n + vx(j_1)*u_n)
          end do                ! j_1

        end do                  ! l_1

      end do                    ! jj_1

      return

      end subroutine eigen_trd_local_2update

      end module eigen_trd_t5_mod
