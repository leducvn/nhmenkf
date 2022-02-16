!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_scaling_mod
!
! Purpose
! =======
!
! eigen_scaling_mod manages the modules for scaling of the input
! matrix.
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

      module eigen_scaling_mod

      use eigen_libs_mod
      use eigen_devel_mod
      use mpi
!$    use omp_lib

      implicit none
      private


      public :: eigen_scaling


      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_scaling
!
! Purpose
! =======
!
! eigen_scaling scales the input matrix in order to avoid overflow and
! underflow numerically.
!
!
! Arguments
! =========
!
! n       (input) integer
!         The dimension of the matrix. N >= 0.
!
! a       (input/output) real(8) array, dimension(lda,*)
!         a contains the element of upper triangluar part of the
!         symmetric matrix.
!
! lda     (input) integer
!         The leading dimension of the array a.
!
! x       (output) real(8)
!         The scaling factor
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_scaling(n, a, lda, x)

      integer, intent(in)    :: n
      real(8), intent(inout) :: a(lda, *)
      integer, intent(in)    :: lda
      real(8), intent(out)   :: x

      real(8)                :: SAFMIN, EPS, SMLNUM, BIGNUM
      real(8)                :: RMIN, RMAX, ANRM, SIGMA
      real(8)                :: t
      real(8), external      :: DLAMCH

      integer                :: ERR_MPI
      integer                :: i, i_1, iloop_sta, iloop_end
      integer                :: j, j_1, jloop_sta, jloop_end


      SAFMIN = DLAMCH('Safe minimum')
      EPS    = DLAMCH('Precision')
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN   = sqrt(SMLNUM)
      RMAX   = min(sqrt(BIGNUM), ONE / sqrt(sqrt(SAFMIN)))

      iloop_sta = eigen_loop_start(1, y_nnod, y_inod)
      iloop_end = eigen_loop_end  (n, y_nnod, y_inod)

      ANRM = ZERO
!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1,jloop_sta,jloop_end)
!$OMP+         REDUCTION(max: ANRM)
      do i_1=iloop_sta,iloop_end

        i   = eigen_translate_l2g(i_1, y_nnod, y_inod)
        j   = i

        jloop_sta = eigen_loop_start(1, x_nnod, x_inod)
        jloop_end = eigen_loop_end  (j, x_nnod, x_inod)

        do j_1=jloop_sta,jloop_end
          ANRM = max(ANRM, abs(a(j_1, i_1)))
        end do

      end do
!$OMP END PARALLEL DO

      t = ANRM
      call MPI_Allreduce(t, ANRM,
     &     1, MPI_DOUBLE_PRECISION,
     &     MPI_MAX, TRD_COMM_WORLD, ERR_MPI)

!-
!     
!     SIGMA is calcurated in order to avoid Under/over-flow
!     
      SIGMA = ONE
      if (ANRM /= ZERO .and. ANRM < RMIN) then
         SIGMA = RMIN / ANRM
      else if (ANRM > RMAX) then
         SIGMA = RMAX / ANRM
      end if

      x = SIGMA
      if (x == ONE) return

      iloop_sta = eigen_loop_start(1, y_nnod, y_inod)
      iloop_end = eigen_loop_end  (n, y_nnod, y_inod)

!$OMP PARALLEL DO
!$OMP+         PRIVATE(i,j,i_1,j_1,jloop_sta,jloop_end)
      do i_1=iloop_sta,iloop_end

        i   = eigen_translate_l2g(i_1, y_nnod, y_inod)
        j   = i

        jloop_sta = eigen_loop_start(1, x_nnod, x_inod)
        jloop_end = eigen_loop_end  (j, x_nnod, x_inod)

        call dscal(jloop_end-jloop_sta+1, x, a(jloop_sta, i_1), 1)
!     do j_1=jloop_sta,jloop_end
!       a(j_1, i_1) = a(j_1, i_1) * x
!     end do

      end do
!$OMP END PARALLEL DO

      return

      end subroutine eigen_scaling

      end module eigen_scaling_mod
