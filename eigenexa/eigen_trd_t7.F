!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_trd_t7_mod
!
! Purpose
! =======
!
! eigen_trd_t7_mod manages the modules for data loader and restorer
! of the panel
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

      module eigen_trd_t7_mod

      use eigen_libs_mod
      use eigen_devel_mod
      use comm_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_trd_panel_load
      public :: eigen_trd_panel_restore

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_panel_load
!
! Purpose
! =======
!
! data loader of the panel
!
!
! Arguments
! =========
!
! a      (input/output) real(8) array, dimension(lda,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! lda    (input) integer
!        The leading dimension of the array a.
!
! u_x    (output) real(8) array, dimension(ldv,*)
!        u_x returns the panel data.
!
! ldv    (input) integer
!        The leading dimension of the array u_x.
!
! d_t    (output) real(8) array, dimension(*)
!        d_t contains the diagonal elements of A
!
! buffer (output) real(8) array, dimension(*)
!        working buffer
!
! m      (input) integer
!        Block factor
!
! i_base (input) integer
!        The offset pointer
!
! i_block(input) integer
!        not used
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_panel_load(a, lda, u_x, ldv,
     &     d_t, buffer,
     &     m, i_base, i_block)

      real(8), intent(inout) ::  a(1:lda,*)
      integer, intent(in)    ::  lda
      real(8), intent(out)   ::  u_x(1:ldv,*)
      integer, intent(in)    ::  ldv
      real(8), intent(out)   ::  d_t(*)
      real(8), intent(out)   ::  buffer(*)
      integer, intent(in)    ::  m
      integer, intent(in)    ::  i_base
      integer, intent(in)    ::  i_block

      integer                ::  x_pos, x_owner_nod
      integer                ::  y_pos, y_owner_nod
      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j, L, k0, i_root


      x_pos = eigen_translate_g2l(i_base+m, x_nnod, x_inod)

      i_root = eigen_owner_node  (i_base+1, y_nnod, y_inod)
      i_2    = eigen_loop_start  (i_base+1, y_nnod, y_inod)
      i_3    = eigen_loop_end    (i_base+m, y_nnod, y_inod)

!$OMP DO
      do i_1=i_2,i_3
        j = eigen_translate_l2g(i_1, y_nnod,y_inod)
        u_x(1:x_pos, j-i_base) = a(1:x_pos, i_1)
      end do                    ! i_1
!$OMP END DO

      do j=1,m,y_nnod
        k0 = min(m+1-j,y_nnod)
        if (k0 == y_nnod) then
          i_1 = j + mod(y_inod-i_root+y_nnod,y_nnod)
!$OMP MASTER
          call allgather_dbl(u_x(1, i_1), buffer, x_pos,
     &         1, y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
!$OMP DO
          do j_2=1,y_nnod
            j_3 = j + mod(j_2-i_root+y_nnod, y_nnod)
!            if ( i_1 /= j_3 ) then
              k_2 = (j_2-1)*x_pos
              u_x(1:x_pos, j_3) = buffer(k_2+1:k_2+x_pos)
!            end if
          end do
!$OMP END DO
        else
!$OMP MASTER
          y_owner_nod = eigen_owner_node (i_base+j, y_nnod, y_inod)
          call bcastw_dbl(u_x(1, j), x_pos, y_owner_nod,
     &         ldv, k0, buffer, 1, y_COMM_WORLD)
!$OMP END MASTER
        end if
      end do

!$OMP BARRIER

      if (diag_0 > 0) then
        i_2 = diag_0
        i_3 = eigen_loop_end(i_base+m-1, y_nnod,y_inod)
        if (i_2 <= i_3) then
          i_4 = x_nnod/n_common
          j_2 = diag_1
          j_3 = eigen_loop_end(i_base+m-1, x_nnod, x_inod)
          j_4 = y_nnod/n_common
          k_2 = 0
          k_3 = (i_3-i_2)/i_4
!$OMP DO
          do  k_1 = k_2, k_3
            i_1 = i_2 + k_1 * i_4
            j_1 = j_2 + k_1 * j_4
            d_t(i_1)    = a(j_1, i_1)
            a(j_1, i_1) = ZERO
          end do                ! k_1
!$OMP END DO
        end if
      end if


      return

      end subroutine  eigen_trd_panel_load

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_trd_panel_restore
!
! Purpose
! =======
!
! data restorer of the panel
!
!
! Arguments
! =========
!
! a      (output) real(8) array, dimension(lda,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! lda    (input) integer
!        The leading dimension of the array a.
!
! u_x    (input) real(8) array, dimension(ldv,*)
!        u_x is filled in zero.
!
! ldv    (input) integer
!        The leading dimension of the array u_x.
!
! d_t    (output) real(8) array, dimension(*)
!        d_t contains the diagonal elements of A
!
! m      (input) integer
!        Block factor
!
! i_base (input) integer
!        The offset pointer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_trd_panel_restore(a, lda, u_x, ldv,
     &     d_t,
     &     m, i_base)

      real(8), intent(out)   ::  a(1:lda,*)
      integer, intent(in)    ::  lda
      real(8), intent(in)    ::  u_x(1:ldv,*)
      integer, intent(in)    ::  ldv
      real(8), intent(in)    ::  d_t(*)
      integer, intent(in)    ::  m
      integer, intent(in)    ::  i_base

      integer                ::  i_1, i_2, i_3, i_4
      integer                ::  j_1, j_2, j_3, j_4
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j, L


      i_2 = eigen_loop_start(i_base+1, y_nnod,y_inod)
      i_3 = eigen_loop_end  (i_base+m, y_nnod,y_inod)
!$OMP DO
      do i_1=i_2,i_3
        j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
        j_3 = eigen_loop_end  (j, x_nnod,x_inod)
        a(1:j_3, i_1) = u_x(1:j_3, j-i_base)
      end do                    ! i_1
!$OMP END DO

      if (diag_0 > 0) then
        i_2 = diag_0
        i_3 = eigen_loop_end(i_base, y_nnod,y_inod)
        if (i_2 <= i_3) then
          i_4 = x_nnod/n_common
          j_2 = diag_1
          j_3 = eigen_loop_end(i_base, x_nnod, x_inod)
          j_4 = y_nnod/n_common
          k_2 = 0
          k_3 = (i_3-i_2)/i_4
!$OMP DO
          do  k_1 = k_2, k_3
            i_1 = i_2 + k_1 * i_4
            j_1 = j_2 + k_1 * j_4
            a(j_1, i_1) = d_t(i_1)
          end do                ! k_1
!$OMP END DO
        end if
      end if


      return

      end subroutine eigen_trd_panel_restore

      end module eigen_trd_t7_mod
