!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_prd_t4x_mod
!
! Purpose
! =======
!
! eigen_prd_t4x_mod manages the modules for the generation of seed
! Householder reflectors.
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

      module eigen_prd_t4x_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use eigen_house_mod, only : tsave
      use mpi
!$    use omp_lib

      implicit none
      private


      public :: eigen_prd_compute_u

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_compute_u
!
! Purpose
! =======
!
! Generator of two contiguous Householder reflectors
!
!
! Arguments
! =========
!
! a       (input/output) real(8) array, dimension(nm,*)
!         a contains the elements of the target panel.
!         On exit, it returns reflector vectors broadcasted.
!
! nm      (input) integer 
!         The leading dimension of the array a.
!
! u_x     (output) real(8) array, dimension(nv,*)
!         u_x returns the reflector vectors broadcasted.
!
! u_y     (output) real(8) array, dimension(nv,*)
!         u_y returns the reflector vectors broadcasted.
!
! nv      (input) integer 
!         The leading dimension of the working arrays, u_x, u_y, ...
!
! u_t     (output) real(8) array, dimension(*)
!         working buffer
!
! v_t     (output) real(8) array, dimension(*)
!         working buffer
!
! i       (input) integer 
!         The global loop counter
!
! c       (output) real(8) array, dimension(MBAND,MBAND)
!         c is an interface array containing the reciprocal of
!         beta(1:2).
!
! e       (output) real(8) array, dimension(ne,*)
!         e contains off-diagonal elements of the output penta-diagonal
!         matrix.
!
! ne      (input) integer
!         The leading dimension of the working array e.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_prd_compute_u(
     &     a, nm,
     &     u_x, u_y, nv,
     &     u_t, v_t, i, c, e, ne)

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(out)   :: u_x(1:nv,*)
      real(8), intent(out)   :: u_y(1:nv,*)
      integer, intent(in)    :: nv
      real(8), intent(out)   :: u_t(*)
      real(8), intent(out)   :: v_t(*)
      integer, intent(in)    :: i
      real(8), intent(out)   :: c(mband,mband)
      real(8), intent(out)   :: e(1:ne,*)
      integer, intent(in)    :: ne

!     integer, parameter     :: VLEN = 2048
      integer, parameter     :: VLEN = 512

      real(8)                :: sgm(MBAND), bet(MBAND)

      integer                :: x_pos, y_owner_nod
      integer                :: loop_sta, loop_end
      integer                :: j_1, j_2, j_3
      integer                :: L, k, itr
      logical                :: mask(MBAND)

      real(8)                :: t, s
      real(8)                :: t11, t12, t22
      real(8)                :: s11, s12, s22
      real(8)                :: r12, rr(MBAND)


      L  = i-MBAND

!$OMP MASTER

      y_owner_nod = eigen_owner_node   (L+1, y_nnod, y_inod)
      x_pos       = eigen_translate_g2l(L+1, x_nnod, x_inod)
      call bcastw_dbl(a(1, 1), x_pos, y_owner_nod,
     &     nm, 2, u_t, 1, y_COMM_WORLD)


!     eps = get_constant_eps()

      loop_sta = eigen_loop_start(1, x_nnod, x_inod)
      loop_end = eigen_loop_end  (L, x_nnod, x_inod)
      do j_1=loop_sta,loop_end
        u_x(j_1, 2) = a(j_1, 2)
        u_x(j_1, 1) = a(j_1, 1)
      end do                    ! j_1

      r12 = ZERO
      do itr=1,2

                                ! compute a Gram matrix
        t11 = ZERO
        t12 = ZERO
        t22 = ZERO
        loop_sta = eigen_loop_start(1, x_nnod, x_inod)
        loop_end = eigen_loop_end  (L, x_nnod, x_inod)
        do j_1=loop_sta,loop_end
          t = u_x(j_1, 2)
          s = u_x(j_1, 1)
          t11 = t11 + t * t
          t12 = t12 + s * t
          t22 = t22 + s * s
        end do                  ! j_1

        u_t(1) = t11
        u_t(2) = t12
        u_t(3) = t22
        j_1 = eigen_owner_index(L-1, x_nnod, x_inod)
        j_2 = eigen_owner_index(L,   x_nnod, x_inod)
        j_3 = eigen_owner_index(L+1, x_nnod, x_inod)
        if (j_1 > 0) then
          u_t(4) = u_x(j_1, 1)
        else
          u_t(4) = ZERO
        end if
        if (j_2 > 0) then
          u_t(5) = u_x(j_2, 1)
        else
          u_t(5) = ZERO
        end if
        if (itr == 1) then
          if (j_1 > 0) then
            u_t(6) = u_x(j_1, 2)
          else
            u_t(6) = ZERO
          end if
          if (j_2 > 0) then
            u_t(7) = u_x(j_2, 2)
          else
            u_t(7) = ZERO
          end if
          if (j_3 > 0) then
            u_t(8) = a(j_3, 2)
          else
            u_t(8) = ZERO
          end if
          k = 8
        else
          k = 5
        end if
        call reduce_dbl(u_t, v_t, k, 3, x_COMM_WORLD)
        t11 = u_t(1)
        t12 = u_t(2)
        t22 = u_t(3)

                                !
                                ! do modified Cholesky factorization
                                ! [1    ][c1  ][1 s12] [c1 0     ][1 s12
                                !]
                                ! [s12 1][  c2][    1]=[s12*c1 c2][    1
                                !]
                                !                      [c1              
                                !    ]
                                !                     =[c1*s12 c1*s12*s1
                                !2+c2]
                                ! if g11 /= 0
                                ! -> c1=g11, s12=g12/c1, c2=g22-c1*s12**
                                !2
                                ! if g11 == 0
                                ! -> c1=0, s12=0, c2=g22
        if (t11 == ZERO) then
          mask(2) = .false.
          s11 = ZERO
          s12 = ZERO
          s22 = t22
        else
          mask(2) = .true.
          s11 = t11
          s12 = t12 / t11
          s22 = t22 - s11 * s12**2
        end if
        if (s22 == ZERO) then
          mask(1) = .false.
        else
          mask(1) = .true.
        end if

                                !
                                ! do multiply R^{-1} from right and comp
                                !ute Q
                                !
        loop_sta = eigen_loop_start(1, x_nnod, x_inod)
        loop_end = eigen_loop_end  (L, x_nnod, x_inod)
        if (mask(2)) then
          if (mask(1)) then
            u_x(loop_sta:loop_end,1) = u_x(loop_sta:loop_end,1)
     &           - s12 * u_x(loop_sta:loop_end,2)
            u_t(4) = u_t(4) - s12 * u_t(6)
            u_t(5) = u_t(5) - s12 * u_t(7)
          end if
        else
          u_x(loop_sta:loop_end,2) = ZERO
          u_t(6) = ZERO
          u_t(7) = ZERO
        end if
        if (mask(1)) then
          continue
        else
          u_x(loop_sta:loop_end,1) = ZERO
          u_t(4) = ZERO
          u_t(5) = ZERO
        end if

                                !
                                ! R_updated <= R_new * R_old
                                !  [1   s12][1 r12] = [1 r12+s12]
                                !  [    1  ][  1  ]   [  1      ]
                                !
        r12 = r12 + s12

        if (itr == 1) then

                                !
                                ! if the condition number is good, skip 
                                !second QR step.
                                !
          s = max(s22, s11)
          t = min(s22, s11)
          if (t * (2**20) > s) then
            exit
!     !     else
!     !        if (TRD_inod == 1) then
!     !           print*,"Do the second modified CholeskyQR", i
!     !        end if
          end if

        end if

      end do

      rr(1) = sqrt(s22)
      rr(2) = sqrt(s11)

      bet(1:2) = ONE
      sgm(1:2) = ZERO

                                !
                                ! Compute Householder reflectors
                                !
      if (mask(2)) then
        sgm(2) = -sign(rr(2), u_t(7))
        j_1 = eigen_owner_index(L, x_nnod, x_inod)
        if (j_1 > 0) then
          u_x(j_1, 2) = u_x(j_1, 2) - sgm(2)
        end if
        u_t(7) =  u_t(7) - sgm(2)
        bet(2) = -u_t(7) * sgm(2)
        if (mask(1)) then
          s = sgm(2) * u_t(5) / bet(2)
          loop_sta = eigen_loop_start(1  , x_nnod, x_inod)
          loop_end = eigen_loop_end  (L-1, x_nnod, x_inod)
          if (loop_sta <= loop_end) then
            u_x(loop_sta:loop_end, 1) = u_x(loop_sta:loop_end, 1)
     &           + s * u_x(loop_sta:loop_end, 2)
          end if
          u_t(4) = u_t(4) + s * u_t(6)
        end if
      end if
      if (mask(1)) then
        sgm(1) = -sign(rr(1), u_t(4))
        j_1 = eigen_owner_index(L-1, x_nnod, x_inod)
        if (j_1 > 0) then
          u_x(j_1, 1) = u_x(j_1, 1) - sgm(1)
        end if
        u_t(4) =  u_t(4) - sgm(1)
        bet(1) = -u_t(4) * sgm(1)
      end if
      if (mask(2)) then
        j_1 = eigen_owner_index(L, x_nnod, x_inod)
        if (j_1 > 0) then
          u_x(j_1, 1) = ZERO
        end if
        j_1 = eigen_owner_index(L+1, x_nnod, x_inod)
        if (j_1 > 0) then
          u_x(j_1, 2) = ZERO
          a  (j_1, 2) = ZERO
        end if
      end if


                                ! [sgm(1)      0][1   0]
                                ! [ 0     sgm(2)][r12 1]
      e(i-1, 1) = sgm(2) * r12
      e(i  , 1) = u_t(8)
      e(i-1, 2) = sgm(1)
      e(i  , 2) = sgm(2)

      c(1,1) = ONE / bet(1)
      c(1,2) = ZERO
      c(2,1) = ZERO
      c(2,2) = ONE / bet(2)


      loop_sta = eigen_loop_start(1, x_nnod, x_inod)
      loop_end = eigen_loop_end  (L, x_nnod, x_inod)
      if (loop_sta <= loop_end) then
        a(loop_sta:loop_end, 2) = u_x(loop_sta:loop_end, 2)
        a(loop_sta:loop_end, 1) = u_x(loop_sta:loop_end, 1)
      end if


      x_pos = eigen_translate_g2l(L+1, x_nnod, x_inod)
      call datacast_dbl2(u_y(1, 1), u_y(1, 2),
     &     u_x(1, 1), u_x(1, 2), v_t, u_t, L+1, 1)

      call eigen_vector_zeropad_x(u_x(1, 1), L+1)
      call eigen_vector_zeropad_x(u_x(1, 2), L+1)
      call eigen_vector_zeropad_y(u_y(1, 1), L+1)
      call eigen_vector_zeropad_y(u_y(1, 2), L+1)

!$OMP END MASTER

      return

      end subroutine eigen_prd_compute_u

      end module eigen_prd_t4x_mod
