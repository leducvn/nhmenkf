!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_prd_t8_mod
!
! Purpose
! =======
!
! eigen_prd_t8_mod manages the modules for initializer and finalizer
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

      module eigen_prd_t8_mod

      use eigen_libs_mod
      use comm_mod
      use eigen_devel_mod
      use eigen_house_mod
      use CSTAB_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_prd_init
      public :: eigen_prd_final

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_init
!
! Purpose
! =======
!
! initializer
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(nm,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! nm     (input) integer
!        The leading dimension of the array a.
!
! n      (input) integer
!        The dimension of matrix A
!
! d_out  (output) real(8) array, dimension(*)
!        d_out contains the diagonal elements of the pena-diagonal
!        matrix.
!
! e_out  (output) real(8) array, dimension(ne,*)
!        e_out contains the off-diagonal elements of the penta-
!        diagonal matrix.
!
! ne     (input) integer
!        The leading dimension of the array e
!
! u_t    (input) real(8) array, dimension(*)
!        working buffer
!
! v_t    (input) real(8) array, dimension(*)
!        working buffer
!
! nv     (input) integer
!        The leading diemnsion of the working arrays.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_prd_init(a, nm, n,
     &     d_out, e_out, ne,
     &     u_t, v_t, nv)

      real(8), intent(inout) :: a(1:nm, *)
      integer, intent(in)    :: nm
      integer, intent(in)    :: n
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne,*)
      integer, intent(in)    :: ne
      real(8), intent(in)    :: u_t(*)
      real(8), intent(in)    :: v_t(*)
      integer, intent(in)    :: nv

      include 'CSTAB.h'

      integer                ::  iloop_sta, iloop_end
      integer                ::  jloop_sta, jloop_end
      integer                ::  i_1, i_4, j_1, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  j

      integer                :: local_size, local_rank


      d_out(1:n)   = ZERO
      e_out(1:n,1) = ZERO
      e_out(1:n,2) = ZERO

      if (diag_0 > 0) then
        iloop_sta = diag_0
        iloop_end = eigen_loop_end(n, y_nnod, y_inod)
        if (iloop_sta <= iloop_end) then
          i_4 = x_nnod/n_common
          jloop_sta = diag_1
          jloop_end = eigen_loop_end(n, x_nnod, x_inod)
          j_4 = y_nnod/n_common
          k_2 = 0
          k_3 = (iloop_end-iloop_sta)/i_4
          do k_1 = k_2, k_3
            i_1 = iloop_sta + k_1 * i_4
            j_1 = jloop_sta + k_1 * j_4
            j   = (i_1-1)*y_nnod+y_inod
            d_out(j) = a(j_1, i_1)
          end do                ! k_1
        end if
      end if

      iloop_sta = eigen_loop_start(1, y_nnod,y_inod)
      iloop_end = eigen_loop_end  (n, y_nnod,y_inod)
      do i_1=iloop_sta,iloop_end
        j = eigen_translate_l2g(i_1, y_nnod,y_inod)
        jloop_sta = eigen_loop_start(j+1, x_nnod,x_inod)
        if (j <= n) then
          a(jloop_sta:nm, i_1) = ZERO
        else
          a(1:nm, i_1) = ZERO
        end if
      end do                    ! i_1

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( local_size > 0 ) then

        allocate(u0_z(nv*local_size+n_columns))
        allocate(v0_z(nv*local_size+n_columns))
        allocate(u1_z(nv*local_size+n_columns))
        allocate(v1_z(nv*local_size+n_columns))
        call CSTAB_adjust_base(u0_z(1), u_t(1), offset1)
        call CSTAB_adjust_base(v0_z(1), v_t(1), offset2)
        call CSTAB_adjust_base(u1_z(1), u_t(1), offset3)
        call CSTAB_adjust_base(v1_z(1), v_t(1), offset4)
        k_1 = (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &       +(L1_LSIZE)
     &       +(L2_LSIZE/8)
        offset1 = offset1 + k_1 * 1
        offset2 = offset2 + k_1 * 2
        offset3 = offset3 + k_1 * 3
        offset4 = offset4 + k_1 * 4
        call CSTAB_round_offset(offset1)
        call CSTAB_round_offset(offset2)
        call CSTAB_round_offset(offset3)
        call CSTAB_round_offset(offset4)

      end if

      allocate (tsave(0:local_size-1, 0:MBAND))

      return

      end subroutine eigen_prd_init

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_final
!
! Purpose
! =======
!
! finalizer
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(nm,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! nm     (input) integer
!        The leading dimension of the array a.
!
! n      (input) integer
!        The dimension of matrix A
!
! d_out  (output) real(8) array, dimension(*)
!        d_out contains the diagonal elements of the pena-diagonal
!        matrix.
!
! e_out  (output) real(8) array, dimension(ne,*)
!        e_out contains the off-diagonal elements of the penta-
!        diagonal matrix.
!
! ne     (input) integer
!        The leading dimension of the array e
!
! u_t    (output) real(8) array, dimension(*)
!        working buffer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_prd_final(a, nm, n, d_out, e_out, ne, u_t)

      real(8), intent(inout) :: a(1:nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: n
      real(8), intent(out)   :: d_out(*)
      real(8), intent(out)   :: e_out(1:ne,*)
      integer, intent(in)    :: ne
      real(8), intent(out)   :: u_t(*)

      integer                ::  x_owner_nod, y_owner_nod
      integer                ::  i_1, iloop_sta, iloop_end, i_4
      integer                ::  j_1
      integer                ::  i, j, L
      real(8)                ::  t

      integer                :: local_size, local_rank


      if (n >= 2) then
        do i = 2+mod(n,2), 2, -1
          L = i-1
          i_1 = eigen_owner_index(i, y_nnod,y_inod)
          j_1 = eigen_owner_index(L, x_nnod,x_inod)
          if (i_1 > 0 .and. j_1 > 0) then
            e_out(i,1) = a(j_1,i_1)
            a(j_1,i_1) = ZERO
          else
            e_out(i,1) = ZERO
          end if
          L = i-2; if (L < 1) cycle
          i_1 = eigen_owner_index(i, y_nnod,y_inod)
          j_1 = eigen_owner_index(L, x_nnod,x_inod)
          if (i_1 > 0 .and. j_1 > 0) then
            e_out(i,2) = a(j_1,i_1)
            a(j_1,i_1) = ZERO
          else
            e_out(i,2) = ZERO
          end if
        end do
      end if

      if (diag_0 > 0) then
        iloop_sta = diag_0
        iloop_end = eigen_loop_end(n, y_nnod,y_inod)
        i_4 = x_nnod/n_common
        do i_1=iloop_sta,iloop_end,i_4
          j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
          j_1 = eigen_translate_g2l(j, x_nnod,x_inod)
          t          = d_out(j)
          d_out(j)   = a(j_1,i_1)
          a(j_1,i_1) = t
        end do                  ! i_1
      end if


      call reduce_dbl(d_out(1),   u_t(1), n, 6, TRD_COMM_WORLD)
      if (mod(n,2) == 1) then
!     x * *
!     x x * o
!     x x x o o
!     o o o o
!     o o o
        u_t(1) = e_out(2,1)
        u_t(2) = e_out(3,1)
        u_t(3) = e_out(3,2)
        call reduce_dbl(u_t(1), u_t(4), 3, 6, TRD_COMM_WORLD)
        e_out(2,1) = u_t(1)
        e_out(3,1) = u_t(2)
        e_out(3,2) = u_t(3)
      else
!     x * o
!     x x o o
!     o o o o
!     o o o
        y_owner_nod = eigen_owner_node(2, y_nnod, y_inod)
        x_owner_nod = eigen_owner_node(1, x_nnod, x_inod)
        if (y_owner_nod == y_inod) then
          call bcast_dbl(e_out(2,1), 1, x_owner_nod,
     &         2, x_COMM_WORLD)
        end if
        call bcast_dbl(e_out(2,1), 1, y_owner_nod,
     &       2, y_COMM_WORLD)
      end if

      e_out(1,1) = ZERO
      e_out(1,2) = ZERO
      e_out(2,2) = ZERO

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if ( local_size > 0 ) then

        deallocate(u0_z, v0_z)
        deallocate(u1_z, v1_z)

      end if

      deallocate(tsave)

      return

      end subroutine eigen_prd_final

      end module eigen_prd_t8_mod
