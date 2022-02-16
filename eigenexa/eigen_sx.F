!--------*---------*---------*---------*---------*---------*---------*-*
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

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_sx
!
! Purpose
! =======
!
! implementation of interface eigen_sx (see eigen_lib.F module)
!
!
! Arguments
! =========
!
! n, nvec, lda, ldz    (input) integer
!  lda and ldz are the leading dimensions corresponding to a and z,
!  respectively. These values must be equal or larger than one obtained 
!  by calling eigen_get_matdims.
!
! m_forward, m_backward  (input) integer, optional
!  default values are m_forward = 48 and m_backward = 128.
!     
! mode  (input) character, optional
!  = 'A' : eigenvalues and corresponding eigenvectors (default)
!  = 'N' : eigenvalues only
!  = 'X' : mode 'A' + improve accuracy of eigenvalues
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_sx(n, nvec, a, lda, w, z, ldz,
     &     m_forward, m_backward, mode)

      use eigen_scaling_mod
      use trbakwy4_mod
      use eigen_prd_mod
      use bisect2_mod
      use dcx_mod
      use eigen_libs_mod, only : eigen_get_matdims
      use eigen_devel_mod
      use comm_mod, only : barrier
      use mpi
!$    use omp_lib

      implicit none

      integer,   intent(in)           :: n
      integer,   intent(in), optional :: nvec
      real(8),   intent(inout)        :: a(lda,*)
      integer,   intent(in)           :: lda
      real(8),   intent(out)          :: w(*)
      real(8),   intent(out)          :: z(ldz,*)
      integer,   intent(in)           :: ldz
      integer,   intent(in), optional :: m_forward
      integer,   intent(in), optional :: m_backward
      character, intent(in), optional :: mode

      real(8), allocatable   :: d(:), e(:), e2(:)

      integer                :: nvec_, mf, mb
      character              :: mode_

      real(8)                :: SIGMA
      integer                :: nm, ny, nb, nme
      integer                :: my_rank, world_size
      integer                :: INFO
      logical                :: flag

      real(8)                :: hs0, hs1, s0, s1
      real(8)                :: ret, ret_1, ret_2, ret_3
!-
!-------------------------------------------------------------------
!-
      call eigen_get_initialized(flag)
      if (.not. flag) then
        return
      end if
      if (TRD_COMM_WORLD == MPI_COMM_NULL) then
        return
      end if
!-
!-------------------------------------------------------------------
!-
      if (n <= 0) then
        print*,"Warining: Negative dimesion is invalid!"
        return
      end if

      call eigen_get_matdims(n, nm, ny)
      if (nm <= 0 .or. ny <= 0) then
        print*,"Warining: Problem size is too large for"
     &       //  " 32bit fortarn integer binary."
        return
      end if
!-
!-------------------------------------------------------------------
!-
      if (present(nvec)) then
        nvec_ = nvec
      else
        nvec_ = n
      end if

      if (present(mode)) then
        mode_ = mode
      else
        mode_ = 'A'
      end if

      if (present(m_forward)) then
        mf = m_forward
      else
        mf = 48
      end if
      if (mode_ /= 'N') then
        if (present(m_backward)) then
          mb = m_backward
        else
          mb = 128
        end if
      end if
!-
!-------------------------------------------------------------------
!-
      hs0 = eigen_get_wtime()
      ret_1 = ZERO; ret_2 = ZERO; ret_3 = ZERO
!-
!-------------------------------------------------------------------
!-
      nme = ((n-1)/4+1)*4
      allocate(d(1:n), e(1:nme*2), e2(1:nme*2))

      world_size = TRD_nnod
      my_rank    = TRD_inod-1
!-
!-------------------------------------------------------------------
!-
      SIGMA = ONE
!     if (present(mode)) then
!     if (mode == 'X' ) then
      call eigen_scaling(n, a(1,1), lda, SIGMA)
!     end if
!     end if
!-
!-------------------------------------------------------------------
!-
#if TIMER_PRINT
      s0 = eigen_get_wtime()
#endif

      call eigen_prd(n, a(1,1), lda, d(1), e(1), nme, mf)
      ret_1 = dble(n)**3*4/3

#if TIMER_PRINT
      s1 = eigen_get_wtime()
      if (my_rank == 0) then
        print 10000,"TRD-BLK ", n, s1-s0, 1D-9*ret_1/(s1-s0), "GFLOPS"
      end if
      call flush(6)
      call barrier(TRD_COMM_WORLD)
#endif
!-
!-------------------------------------------------------------------
!-
#if TIMER_PRINT
      s0 = eigen_get_wtime()
#endif

      if ((mode_ == 'A' .or. mode_ == 'X') .and. nvec_ > 0) then
        e2(0*nme+1:0*nme+N-1) = e(0*nme+2:0*nme+N)
        e2(0*nme+N) = 0
        e2(1*nme+1:1*nme+N-2) = e(1*nme+3:1*nme+N)
        e2(1*nme+N-1) = 0
        e2(1*nme+N) = 0
        w(1:n)=d(1:n)
        call eigen_dcx(n, nvec_, w, e2, nme, z, ldz, INFO, ret_2)
        if (mode_ == 'X') then
          call eigen_bisect2(d, e, e(1+nme), w, n, 1)
        end if
      end if
      if (mode_ == 'N' .OR. nvec_ <= 0) then
        call eigen_bisect2(d, e, e(1+nme), w, n, 0)
      end if

#if TIMER_PRINT
      s1 = eigen_get_wtime()
      if (my_rank == 0) then
        print 10000,"D&C     ", n, s1-s0, 1D-9*ret_2/(s1-s0), "GFLOPS"
      end if
      call flush(6)
      call barrier(TRD_COMM_WORLD)
#endif
      if (mode_ == 'N' .OR. nvec_ <= 0) then
        goto 99999
      end if
!-
!-------------------------------------------------------------------
!-
#if TIMER_PRINT
      s0 = eigen_get_wtime()
#endif

      nb = MBAND
      call eigen_common_trbakwy(n, nvec_, a(1,1), lda, z(1,1), ldz,
     &     e(1+nme), mb, nb)
      ret_3 = 2*dble(nvec_)*dble(n)**2

#if TIMER_PRINT
      s1 = eigen_get_wtime()
      if (my_rank == 0) then
        print 10000,"TRDBAK  ", n, s1-s0, 1D-9*ret_3/(s1-s0), "GFLOPS"
      end if
      call flush(6)
#endif
!-
!-------------------------------------------------------------------
!-
      if (SIGMA /= ONE) then
        if (my_rank == 0) then
          print*,"Scaling done ",SIGMA
        end if
        SIGMA = ONE / SIGMA
        call dscal(n, w, SIGMA, 1)
      end if
!-
!-------------------------------------------------------------------
!-
#if TIMER_PRINT
      hs1 = eigen_get_wtime()
      ret = ret_1 + ret_2 + ret_3
      if (my_rank == 0) then
        print*,"Total(eigen_sx)",
     &       hs1-hs0, 1D-9*ret/(hs1-hs0), "GFLOPS"
      end if
      call flush(6)
#endif
!-
!-------------------------------------------------------------------
!-
99999 continue
!-
!-------------------------------------------------------------------
!-
      hs1 = eigen_get_wtime()
      ret = ret_1 + ret_2 + ret_3
      if (ret_2 == 0) ret = -ret
      a(1, 1) = ret
      a(2, 1) = hs1-hs0
#if TIMER_PRINT
      a(3, 1) = comm_time_reduction
     &     + comm_time_dc
     &     + comm_time_backtrafo
#else
      a(3, 1) = -1d0
#endif
!-
!-------------------------------------------------------------------
!-
      deallocate(d, e, e2)
!-
!-------------------------------------------------------------------
!-
10000 format (X, A8, I8, E25.16e2, E25.16e2, X, A)
!-
!-------------------------------------------------------------------
!-
      end subroutine eigen_sx
