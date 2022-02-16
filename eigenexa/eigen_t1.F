!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_t1_mod
!
! Purpose
! =======
!
! eigen_t1_mod manages the modules for the Rank-2k update routine.
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

      module eigen_t1_mod

      use eigen_libs_mod
      use eigen_devel_mod
!$    use omp_lib

      implicit none
      private

      public :: eigen_common_2update

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_common_2update
!
! Purpose
! =======
!
! parallel Rank-2k update, A := A - Vy*U - Uy*V
!
!
! Arguments
! =========
!
! ar      (input/output) real(8) array, dimension(nma,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! nma     (input) integer
!         The leading dimenstion of the array a.
!
! ur      (input) real(8) array, dimension(nmv,*)
!         ur contains m_size vectors of U.
!
! uyr     (input) real(8) array, dimension(nmv,*)
!         uyr contains m_size vectors of U^t.
!
! vr      (input) real(8) array, dimension(nmv,*)
!         vr contains m_size vectors of V.
!
! vyr     (input) real(8) array, dimension(nmv,*)
!         vyr contains m_size vectors of V^t.
!
! nmv     (input) integer
!         The leading dimenstion of the array U and V.
!
! m_size  (input) integer
!         The number of vectors bunched.
!
! i_base  (input) integer
!         The dimenion of the sub-matrix of A
!     
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_common_2update(
     &     ar, nma,
     &     ur, uyr, vr, vyr, nmv,
     &     m_size, i_base)

      real(8), intent(inout) :: ar (nma, *)
      integer, intent(in)    :: nma
      real(8), intent(in)    :: ur (nmv, *)
      real(8), intent(in)    :: uyr(nmv, *)
      real(8), intent(in)    :: vr (nmv, *)
      real(8), intent(in)    :: vyr(nmv, *)
      integer, intent(in)    :: nmv
      integer, intent(in)    :: m_size
      integer, intent(in)    :: i_base

#if __FUJITSU
      integer                :: BLAS_CHUNK_N
      integer                :: BLAS_CHUNK_M
#else
      integer, parameter     :: BLAS_CHUNK_N = 8*4
      integer, parameter     :: BLAS_CHUNK_M = 96
#endif

      integer                :: k1, k2, m, n
      integer                :: i_1,i_2,i_3,i_4
      integer                :: j_1,j_2,j_3,j_4
      integer                :: i_loop_sta, i_loop_end
      integer                :: j_loop_sta, j_loop_end
      integer                :: blk_size1, blk_size2
      integer                :: ii_step
      integer                :: local_rank, local_size

      intrinsic              :: min, max
      external               :: dgemm


      if (i_base <= 0) return

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      m  = m_size
      n  = eigen_translate_g2l(i_base, y_nnod, y_inod)

!     12 + 8 + 4 = 24 = 8*3
#if __FUJITSU
!     BLAS_CHUNK_M = max((n-1)/ 3+1, 32)
!     BLAS_CHUNK_N = max((n-1)/12+1,  8)
!     BLAS_CHUNK_M = ((BLAS_CHUNK_M-1)/4+1)*4
!     BLAS_CHUNK_N = ((BLAS_CHUNK_N-1)/2+1)*2
      BLAS_CHUNK_M = 256*3
      BLAS_CHUNK_N = 7*8
#endif

      j_loop_sta = 1
      j_loop_end = eigen_loop_end(i_base, x_nnod, x_inod)

      ii_step = 0
      do j_1 = j_loop_sta, j_loop_end, BLAS_CHUNK_M

        j_4 = min(j_1+BLAS_CHUNK_M-1, j_loop_end)

        k1 = eigen_translate_l2g(j_1, x_nnod, x_inod)

        i_loop_sta = eigen_loop_start(k1, y_nnod, y_inod)
        i_loop_sta = max(1, i_loop_sta)

        i_loop_end = eigen_loop_end(i_base, y_nnod, y_inod)

        do i_1 = i_loop_sta, i_loop_end, BLAS_CHUNK_N

          i_4 = min(i_1+BLAS_CHUNK_N-1, i_loop_end)
          i_2 = i_1; i_3 = i_4

          k2  = eigen_translate_l2g(i_4, y_nnod, y_inod)

          j_2 = j_1
          j_2 = max(1, j_2)

          j_3 = eigen_loop_end(k2, x_nnod, x_inod)
          j_3 = min(j_4, j_3)

          blk_size1 = j_3-j_2+1
          blk_size2 = i_3-i_2+1

          if (blk_size1 > 0 .and. blk_size2 > 0) then

            if (mod(ii_step, local_size) == local_rank) then

              call dgemm('N','T',
     &             blk_size1, blk_size2, m,
     &             MONE, ur (j_1, 1),   nmv,
     &             vyr(i_1, 1),   nmv,
     &             ONE,  ar (j_1, i_1), nma)

              call dgemm('N','T',
     &             blk_size1, blk_size2, m,
     &             MONE, vr (j_1, 1),   nmv,
     &             uyr(i_1, 1),   nmv,
     &             ONE,  ar (j_1, i_1), nma)

            end if

            ii_step = ii_step + 1

          end if

        end do                  ! i_1

      end do                    ! j_1

      return

      end subroutine  eigen_common_2update

      end module eigen_t1_mod
