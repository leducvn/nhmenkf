!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_prd_t7_mod
!
! Purpose
! =======
!
! eigen_prd_t7_mod manages the modules for data loader and restorer
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

      module eigen_prd_t7_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private


      public :: eigen_prd_panel_load
      public :: eigen_prd_panel_store

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_panel_load
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
! a      (input/output) real(8) array, dimension(nm,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! w      (output) real(8) array, dimension(nm,*)
!        W returns the target panel data.
!
! nm     (input) integer
!        The leading dimension of the array w.
!
! d_t    (output) real(8) array, dimension(*)
!        d_t contains the diagonal elements of A.
!
! u_x    (output) real(8) array, dimension(nv,*)
! u_y    (output) real(8) array, dimension(nv,*)
! v_x    (output) real(8) array, dimension(nv,*)
! v_y    (output) real(8) array, dimension(nv,*)
!        filled in zeros.
!
! nv     (input) integer
!        The leading dimension of the arrays.
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

      subroutine eigen_prd_panel_load(a, w, nm,
     &     d_t,
     &     u_x, u_y, v_x, v_y, nv,
     &     m, i_base, i_block)

      real(8), intent(inout) :: a(1:nm,*)
      real(8), intent(out)   :: w(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(out)   :: d_t(*)
      real(8), intent(out)   :: u_x(1:nv,*)
      real(8), intent(out)   :: u_y(1:nv,*)
      real(8), intent(out)   :: v_x(1:nv,*)
      real(8), intent(out)   :: v_y(1:nv,*)
      integer, intent(in)    :: nv
      integer, intent(in)    :: m
      integer, intent(in)    :: i_base
      integer, intent(in)    :: i_block

      integer                ::  x_pos, y_pos
      integer                ::  iloop_sta, iloop_end
      integer                ::  jloop_sta, jloop_end
      integer                ::  i_1, i_4, j_1
      integer                ::  k_1, k_2, k_3
      integer                ::  i, j, L


      iloop_sta = eigen_loop_start(i_base+1, y_nnod,y_inod)
      iloop_end = eigen_loop_end  (i_base+m, y_nnod,y_inod)
      jloop_sta = eigen_loop_start(1,        x_nnod,x_inod)
      jloop_end = eigen_loop_end  (i_base+m, x_nnod,x_inod)

      do i_1=iloop_sta,iloop_end
        j = eigen_translate_l2g(i_1, y_nnod,y_inod)
        do j_1=jloop_sta,jloop_end
          w(j_1,j-i_base) = a(j_1,i_1)
        end do                  ! j_1
      end do                    ! i_1

      if (diag_0 > 0) then
        iloop_sta = diag_0
        iloop_end = eigen_loop_end(i_base+m-1, y_nnod,y_inod)
        i_4 = x_nnod/n_common
        if (y_nnod == x_nnod) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
          do i_1=iloop_sta,iloop_end,i_4
            d_t(i_1)   = a(i_1,i_1)
            a(i_1,i_1) = ZERO
          end do                ! i_1
        else
          do i_1=iloop_sta,iloop_end,i_4
            j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
            j_1 = eigen_translate_g2l(j, x_nnod,x_inod)
            d_t(i_1)   = a(j_1,i_1)
            a(j_1,i_1) = ZERO
          end do                ! i_1
        end if
      end if

      i = i_base+m
      L = i - 2
      y_pos  = eigen_translate_g2l(i, y_nnod,y_inod)
      x_pos  = eigen_translate_g2l(i, x_nnod,x_inod)

      k_2 = m
      k_3 = max(1,3*(2-i_block))
      do k_1=k_3,k_2
        do j_1=1,x_pos
          u_x(j_1,k_1) = ZERO
          v_x(j_1,k_1) = ZERO
        end do                  ! j_1
        do j_1=1,y_pos
          u_y(j_1,k_1) = ZERO
          v_y(j_1,k_1) = ZERO
        end do                  ! j_1
      end do                    ! k_1

      return

      end subroutine eigen_prd_panel_load

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_prd_panel_store
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
! a      (output) real(8) array, dimension(nm,*)
!        a contains the element of upper triangluar part of the
!        symmetric matrix.
!
! w      (input) real(8) array, dimension(nm,*)
!        W returns the target panel data.
!
! nm     (input) integer
!        The leading dimension of the array w.
!
! d_t    (input) real(8) array, dimension(*)
!        d_t contains the diagonal elements of A.
!
! m      (input) integer
!        Block factor
!
! i_base (input) integer
!        The offset pointer
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_prd_panel_store(a, w, nm,
     &     d_t,
     &     m, i_base)

      real(8), intent(out)   :: a(1:nm,*)
      real(8), intent(in)    :: w(1:nm,*)
      integer, intent(in)    :: nm
      real(8), intent(in)    :: d_t(*)
      integer, intent(in)    :: m
      integer, intent(in)    :: i_base

      integer                ::  iloop_sta, iloop_end
      integer                ::  jloop_sta, jloop_end
      integer                ::  i_1, i_4, j, j_1


      iloop_sta = eigen_loop_start(i_base+1, y_nnod,y_inod)
      iloop_end = eigen_loop_end  (i_base+m, y_nnod,y_inod)
      do i_1=iloop_sta,iloop_end
        j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
        jloop_end = eigen_loop_end(j, x_nnod,x_inod)
        do j_1=1,jloop_end
          a(j_1,i_1) = w(j_1,j-i_base)
        end do                  ! j_1
      end do                    ! i_1

      if (diag_0 > 0) then
        iloop_sta = diag_0
        iloop_end = eigen_loop_end(i_base, y_nnod,y_inod)
        i_4 = x_nnod/n_common
        if (iloop_sta <= iloop_end) then
          if (y_nnod == x_nnod) then
            do i_1=iloop_sta,iloop_end,i_4
              j_1 = i_1
              a(j_1,i_1) = d_t(i_1)
            end do              ! i_1
          else
            do i_1=iloop_sta,iloop_end,i_4
!     j   = eigen_translate_l2g(i_1, y_nnod,y_inod)
!     j_1 = eigen_translate_g2l(j, x_nnod,x_inod)
              j   = (i_1-1)*y_nnod+y_inod
              j_1 = (j-1)/x_nnod+1
              a(j_1,i_1) = d_t(i_1)
            end do              ! i_1
          end if
        end if
      end if

      return

      end subroutine  eigen_prd_panel_store

      end module eigen_prd_t7_mod
