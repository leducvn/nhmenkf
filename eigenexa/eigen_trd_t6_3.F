!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t6_3_mod
!
! Purpose
! =======
!
! eigen_trd_t6_3_mod manages the modules for the vx generator.
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

      module eigen_trd_t6_3_mod

      use eigen_house_mod
      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_compute_v

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_compute_v
!
! Purpose
! =======
!
! the vx generator.
!
!
! Arguments
! =========
!
! u_x    (input) real(8) array, dimension(*)
!        u_x contains the reflector vector.
!
! v_x    (output) real(8) array, dimension(*)
!        v_x returns the counterpart vector to the reflector u_x.
!
! v_y    (output) real(8) array, dimension(*)
!        v_y returns the counterpart vector to the reflector u_x.
!
! ux     (input) real(8) arrya, dimension(ldv,*)
!        ux contains U, in which previously calculated u's are stored.
!
! vx     (input) real(8) arrya, dimension(ldv,*)
!        vx contains V, in which previously calculated v's are stored.
!
! ldv    (input) integer
!        The leading dimension of the working arrays.
!
! u_t    (input/output) real(8) array, dimension(*)
!        working buffer
!
! v_t    (input/output) real(8) array, dimension(*)
!        working buffer
!
! beta   (input) real(8)
!        beta contains the sum of square of the reflector vector u_x.
!
! u_n    (input) real(8)
!        n-th element of u_x
!
! v_n    (input/output) real(8)
!        n-th element of v_x
!
! loop_info (input) integer array, dimension(*)
!        Information regarding loop and so on.
!
! flag_putoff (input) logical
!        Flag for putting off broadcast
!
! v_nn   (output) real(8) array, dimension(*)
!        u_nn contains a part of u_x(:) which must be commonly shared.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_compute_v(
     &     u_x, v_x, v_y,
     &     ux, vx, ldv,
     &     u_t, v_t, beta, u_n, v_n,
     &     loop_info, flag_putoff, v_nn)
     
      real(8), intent(in)    ::  u_x(*)
      real(8), intent(out)   ::  v_x(*)
      real(8), intent(out)   ::  v_y(*)
      real(8), intent(in)    ::  ux(1:ldv,*)
      real(8), intent(in)    ::  vx(1:ldv,*)
      integer, intent(in)    ::  ldv
      real(8), intent(inout) ::  u_t(*)
      real(8), intent(inout) ::  v_t(*)
      real(8), intent(in)    ::  beta
      real(8), intent(in)    ::  u_n
      real(8), intent(inout) ::  v_n
      integer, intent(in)    ::  loop_info(*)
      logical, intent(in)    ::  flag_putoff
      real(8), intent(out)   ::  v_nn(*)

      integer                ::  i, i_base, m0, num_v_ptr, num_v_done
      integer                ::  j, k, L, n, LL
      integer                ::  k_1, k_2, k_3
      integer                ::  j_1, j_2, j_3
      integer                ::  l_1, l_2, l_3, l_4
      integer                ::  jj_1
      integer                ::  LX, blk_0

      include 'CSTAB.h'

      real(8)                ::  w0
      real(8)                ::  u0, v0
      real(8)                ::  u1, v1
      real(8)                ::  u2, v2
      real(8)                ::  ux0, vx0
      real(8)                ::  ux1, vx1
      real(8)                ::  ux2, vx2

      real(8)                ::  alpha, prod_uv
      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod


!     v=v-(UV+VU)u
!     
      m0         = loop_info(LOOP_INFO_V_WIDTH)
      i_base     = loop_info(LOOP_INFO_I_BASE)
      num_v_ptr  = loop_info(LOOP_INFO_V_CURNT)
      i          = i_base + num_v_ptr
      num_v_done = m0 - num_v_ptr

      L = i-1
      n = eigen_translate_g2l(L, x_nnod,x_inod)

!$OMP MASTER
      call eigen_vector_zeropad_x(v_x(1), L)
!$OMP END MASTER

!$OMP BARRIER

      prod_uv = u_t(1+2*num_v_done)

      if (num_v_done > 0) then

!  do l_2,l_3
!   =>
!  do l_2,l_4-1
!  do l_4,l_3,3
!
        l_2 = num_v_ptr+1
        l_3 = m0
        l_4 = mod(num_v_done, 3)+l_2

        alpha = ddot(num_v_done, u_t(1), 1, u_t(1+num_v_done), 1)
        prod_uv = prod_uv - 2*alpha

        LX  = 64

!$OMP DO
        do jj_1=1,n,LX

          j_2 = jj_1; j_3 = min(jj_1+LX-1, n)

          if (l_4-1==l_2) then
            l_1 = l_2

            j = l_1-num_v_ptr

            u0 = u_t(0+j+num_v_done*0)
            v0 = u_t(0+j+num_v_done*1)
            do j_1=j_2,j_3
              w0 = v_x(j_1)
              ux0 = ux(j_1, l_1+0)
              vx0 = vx(j_1, l_1+0)
              w0 = w0
     &             -ux0*u0
     &             -vx0*v0
              v_x(j_1) = w0
            end do              ! j_1
          end if

          if (l_4-1==l_2+1) then
            l_1 = l_2

            j = l_1-num_v_ptr

            u0 = u_t(0+j+num_v_done*0)
            v0 = u_t(0+j+num_v_done*1)
            u1 = u_t(1+j+num_v_done*0)
            v1 = u_t(1+j+num_v_done*1)
            do j_1=j_2,j_3
              w0 = v_x(j_1)
              ux0 = ux(j_1, l_1+0)
              vx0 = vx(j_1, l_1+0)
              ux1 = ux(j_1, l_1+1)
              vx1 = vx(j_1, l_1+1)
              w0 = w0
     &             -ux0*u0
     &             -vx0*v0
     &             -ux1*u1
     &             -vx1*v1
              v_x(j_1) = w0
            end do              ! j_1
          end if

          do l_1=l_4,m0,3       ! 2

            j = l_1-num_v_ptr

            u0 = u_t(0+j+num_v_done*0)
            v0 = u_t(0+j+num_v_done*1)
            u1 = u_t(1+j+num_v_done*0)
            v1 = u_t(1+j+num_v_done*1)
            u2 = u_t(2+j+num_v_done*0)
            v2 = u_t(2+j+num_v_done*1)
            do j_1=j_2,j_3
              w0 = v_x(j_1)
              ux0 = ux(j_1, l_1+0)
              vx0 = vx(j_1, l_1+0)
              ux1 = ux(j_1, l_1+1)
              vx1 = vx(j_1, l_1+1)
              ux2 = ux(j_1, l_1+2)
              vx2 = vx(j_1, l_1+2)
              w0 = w0
     &             -ux0*u0
     &             -vx0*v0
     &             -ux1*u1
     &             -vx1*v1
     &             -ux2*u2
     &             -vx2*v2
              v_x(j_1) = w0
            end do              ! j_1
          end do                ! l_1

        end do                  ! jj_1
!$OMP ENDDO

      end if

      if (beta /= ZERO) then
        alpha = prod_uv/(2*beta)
!$OMP DO
        do j_1=1,n
          v_x(j_1) = (v_x(j_1)-alpha*u_x(j_1))/beta
        end do                  ! j_1
!$OMP ENDDO NOWAIT
      end if

!$OMP MASTER

!
! u_n and v_n must be necessary as global variables
! to update u_x(next index) in eigen_trd_local_2update0.
!
      v_n = (v_n-alpha*u_n)/beta

!$OMP END MASTER

!$OMP BARRIER

!$OMP MASTER

!
! v_nn is partly stored from v_x, but not consistent as a global array.
! allreduce called in eigen_trd_au (or other) will make consistent.
!
      v_nn(1:m0) = ZERO
      if (num_v_ptr > 0) then
!        v_nn(1:num_v_ptr) = ZERO
        k_2 = eigen_loop_start(i_base+1, x_nnod,x_inod)
        k_3 = eigen_loop_end  (i,        x_nnod,x_inod)
        do k_1=k_2,k_3
          j = eigen_translate_l2g(k_1, x_nnod,x_inod)
          v_nn(j-i_base) = v_x(k_1)
        end do
      end if

!
! Re-distribution of v_x (column-vector) to v_y (row-vector)
! if this is not the last re-distribution, it will be
! put off and combined with u_y (next index) in eigen_trd_u,
! otherwise done here.
!
      if (.not. flag_putoff) then
        if (flag_oncache) then
          l_2 = num_v_done + 1
          x_pos = eigen_translate_g2l(L+l_2-1, x_nnod,x_inod)
          call datacast_dblx(l_2, v_y, v_x, ldv,
     &                       v_t(1), u_t(1), L+l_2-1, 2)
          do k_1=1,l_2
          call eigen_vector_zeropad_y(v_y(1+ldv*(k_1-1)), L+k_1-1)
          end do
        else
          x_pos = eigen_translate_g2l(L, x_nnod,x_inod)
          call datacast_dbl(v_y(1), v_x(1), v_t(1), u_t(1), L, 2)
          call eigen_vector_zeropad_y(v_y(1), L)
        end if
      end if

      call eigen_vector_zeropad_x(v_x(1), L)

!$OMP END MASTER

      return

      end subroutine eigen_trd_compute_v

      end module eigen_trd_t6_3_mod
