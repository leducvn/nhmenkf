!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: trbakwy4_mod
!
! Purpose
! =======
!
! trbakwy4_mod manages the modules for back-transformation routines.
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

      module trbakwy4_mod

      use trbakwy4_body_mod
      use comm_mod
      use eigen_libs_mod
      use eigen_devel_mod
      use eigen_trbak_mod
      use CSTAB_mod
      use mpi
!$    use omp_lib

      implicit none
      private
      

      public  :: eigen_common_trbakwy
      private :: eigen_trbakwy_body

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_common_trbakwy
!
! Purpose
! =======
!
! eigen_common_trbakwy is the main body of the Householder back
! transformation routine by taking use of the WY representation
! algorithm.
!
!
! Arguments
! =========
!
! n       (input) integer
!         The dimension of the input matrix. N >= 0.
!
! nvec    (input) integer
!         The number of vectors to be transformed. N >= NVEC >= 0.
!
! a       (input) real(8) array, dimension(*)
!         a contains the data of the householder reflectors.
!
! nma0    (input) integer
!         The leading dimension of the array a.
!
! z       (output) real(8) array, dimension(*)
!         z contains the elements of the eigenvectors.
!
! nmz0    (input) integer
!         The leading dimension of the array z.
!
! e       (input) real(8) array, dimension(*)
!         e contains the data of the householder transformation.
!
! m0      (input) integer
!         Block factor
!
! iblk    (input) integer
!         Auxiliary data for back-transformation (tri-/penta-/...).
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_common_trbakwy(n, nvec, a, nma0, z, nmz0,
     &     e, m0, iblk)

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      real(8), intent(in)    :: a(*)
      integer, intent(in)    :: nma0
      real(8), intent(out)   :: z(*)
      integer, intent(in)    :: nmz0
      real(8), intent(in)    :: e(*)
      integer, intent(in)    :: m0
      integer, intent(in)    :: iblk

      real(8),  allocatable  :: d(:)
      real(8),  allocatable  :: v1(:), ss1(:)
      real(8),  allocatable  :: v2(:), ss2(:)
      real(8),  allocatable  :: v3(:), ss3(:)
      real(8),  allocatable  :: tt(:)

      include 'CSTAB.h'

      integer                :: nma, nmz, m
      integer                :: na, nm, kx
      integer                :: i_v1, i_v2, i_v3
      integer                :: i_s1, i_s2, i_s3
      integer                :: i_t
      integer                :: i_stat, ierr


      nma = nma0
      nmz = nmz0
      m   = min(nsm, m0)
      if (m < 1) m = 1

      call eigen_timer_reset(2, 4, 0, 1)

      na   = (n-1)/y_nnod+1
      na   = na + mod(na - 1,2)
      call CSTAB_get_optdim(nma,9,16*4,16*6,nm)

      allocate(
     &     d(1:n),
     &     v1(1:max(nm*m,n)+n_columns),
     &     ss1(1:na*m+ns0+n_columns),
     &     v2(1:max(nm*m,n)+n_columns),
     &     ss2(1:na*m+ns0+n_columns),
     &     v3(1:max(nm*m,n)+n_columns),
     &     ss3(1:na*m+ns0+n_columns),
     &     tt(1:na*m+ns0+n_columns),
     &     stat=i_stat)
      if (i_stat /= 0) then
        if (TRD_inod == 1) print*, "Memory allocation error."
        call MPI_Abort(TRD_COMM_WORLD, MPI_ERR_OTHER, ierr)
      end if

      d = ZERO
      v1 = ZERO
      ss1 = ZERO
      v2 = ZERO
      ss2 = ZERO
      v3 = ZERO
      ss3 = ZERO
      tt = ZERO

      call CSTAB_adjust_base(v1(1), z(1), i_v1)
      call CSTAB_adjust_base(ss1(1), z(1), i_s1)
      call CSTAB_adjust_base(v2(1), z(1), i_v2)
      call CSTAB_adjust_base(ss2(1), z(1), i_s2)
      call CSTAB_adjust_base(v3(1), z(1), i_v3)
      call CSTAB_adjust_base(ss3(1), z(1), i_s3)
      call CSTAB_adjust_base(tt(1), z(1), i_t)
      kx = (L1_WINDOW/8)
!     &           +(L1_WINDOW)
!     &           +(L1_LSIZE/8)
     &     +(L1_LSIZE)
     &     +(L2_LSIZE/8)
      i_v1 = i_v1 + kx*5
      i_s1 = i_s1 + kx*1
      i_v2 = i_v2 + kx*5
      i_s2 = i_s2 + kx*1
      i_v3 = i_v3 + kx*5
      i_s3 = i_s3 + kx*1
      i_t = i_t + kx*1
      call CSTAB_round_offset(i_v1)
      call CSTAB_round_offset(i_s1)
      call CSTAB_round_offset(i_v2)
      call CSTAB_round_offset(i_s2)
      call CSTAB_round_offset(i_v3)
      call CSTAB_round_offset(i_s3)
      call CSTAB_round_offset(i_t)

      allocate(TRBK_buf(1:m*nm))
     
      call MPI_Barrier(TRD_COMM_WORLD,ierr)

      TRBK_TIME_REDUC_OVERHEAD_x = Reduce_cont_Overhead_x

      call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP PARALLEL
      call eigen_trbakwy_body(n, nvec,
     &     a(1), nma,
     &     z(1), nmz,
     &     e(1), d(1),
     &     v1(1+i_v1), v2(1+i_v2), v3(1+i_v3), nm, m,
     &     ss1(1+i_s1), ss2(1+i_s2), ss3(1+i_s3),
     &     tt(1+i_t), iblk, na)
!$OMP END PARALLEL
      call MPI_Barrier(TRD_COMM_WORLD,ierr)

      deallocate(TRBK_buf)

      deallocate(d)
      deallocate(v1)
      deallocate(ss1)
      deallocate(v2)
      deallocate(ss2)
      deallocate(v3)
      deallocate(ss3)
      deallocate(tt)

      comm_time_backtrafo =
     &     eigen_timer_print('EigenExa (Back-transformation)')

      return

      end subroutine eigen_common_trbakwy

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine eigen_trbakwy_body(n, nvec, a, nma, z, nmz, e, d,
     &     v1, v2, v3, nm, m, ss1, ss2, ss3, tt, iblk, nss)

      integer, intent(in)    :: n
      integer, intent(in)    :: nvec
      real(8), intent(in)    :: a(1:nma,*)
      integer, intent(in)    :: nma
      real(8), intent(inout) :: z(1:nmz,*)
      integer, intent(in)    :: nmz
      real(8), intent(in)    :: e(1:n)
      real(8), intent(inout) :: d(1:n)
      real(8), intent(out)   :: v1(1:nm,*)
      real(8), intent(out)   :: v2(1:nm,*)
      real(8), intent(out)   :: v3(1:nm,*)
      integer, intent(in)    :: nm
      integer, intent(in)    :: m
      real(8), intent(out)   :: ss1(*)
      real(8), intent(out)   :: ss2(*)
      real(8), intent(out)   :: ss3(*)
      real(8), intent(out)   :: tt(*)
      integer, intent(in)    :: iblk
      integer, intent(in)    :: nss

      real(8), pointer       :: wk(:)
      integer                :: nodes(0:nsm-1)
      integer                :: nx, lwk, ierr
      integer                :: i, L, mode
      integer                :: x_root
      integer                :: iloop_sta, iloop_end
      integer                :: jloop_sta, jloop_end
      integer                :: i_1, i_4, j_1
      real(8)                :: s0, s1, s2, s3
      real(8)                :: d0, d1, d2, ds, de, dcom, dx, dy, dz


#if TIMER_PRINT
      dx = ZERO; dy = ZERO; dz = ZERO
#endif

!$OMP MASTER
      call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
#if TIMER_PRINT
      d1 = eigen_get_wtime()
#endif
      dcom = ZERO
!-
#if AT_BCAST_OVERLAP
!$    if(.TRUE.) DO_OVERLAP_BCAST_LEVEL  = 2
!$    TRBK_TIME_COUNTER  = 0
!$    TRBK_TIME_INTERVAL = 1
!$    TRBK_TIME_NEXT     = 2
!$    TRBK_SWITCHED      = 0
!$    TRBK_TIME_TRBK1    = ZERO
!$    TRBK_TIME_TRBK1_   = ZERO
!$    TRBK_TIME_TRBK1x   = ZERO
!$    TRBK_TIME_TRBK1x_  = ZERO
!$    TRBK_TIME_TRBK2    = ZERO
!$    TRBK_TIME_BCAST    = ZERO
!$    TRBK_TIME_REDUC    = ONE
!$    TRBK_TIME_FR       = ZERO
!===  
!$    call sync_other_than_master_init(TRBK_LOCK, TRBK_MASK)
!===  
#else
!$    if(.true.) DO_OVERLAP_BCAST_LEVEL  = 0
#endif
!-
      lwk=max(((m-1)/y_nnod+1)*((n-1)/x_nnod+1),n)
      allocate(wk(lwk))
!$OMP END MASTER
!$OMP BARRIER

      iloop_sta = eigen_loop_start(1+iblk, y_nnod,y_inod)
      iloop_end = eigen_loop_end  (n,      y_nnod,y_inod)

!$OMP MASTER
      d(1:n) = ZERO
      do i_1=iloop_sta,iloop_end
        i = (i_1-1)*y_nnod+y_inod
        l = i-iblk
        x_root = mod(l-1,x_nnod)+1
        if (x_root == x_inod) then
          j_1  = (l-1)/x_nnod+1
          d(i) = a(j_1,i_1)
        end if
      end do                    ! i_1
      call reduce_dbl(d, wk, n, 1, TRD_COMM_WORLD)
      do i=1+iblk,n
        if (e(i)*d(i) == ZERO) then
          s0 = ZERO
        else
          s0 = (ONE/d(i))/e(i)
        end if
        d(i) = s0
      end do                    ! i
!$OMP END MASTER

      iloop_sta = eigen_loop_start(1,    y_nnod,y_inod)
      iloop_end = eigen_loop_end  (nvec, y_nnod,y_inod)

!----------------------------------------------------------------

      nx = min(mod(n-(1+iblk)+1,m)+(1+iblk)-1,n)

      do i=(1+iblk),nx

!$OMP BARRIER

        if (e(i) == ZERO) cycle

        jloop_sta = eigen_loop_start(1,      x_nnod,x_inod)
        jloop_end = eigen_loop_end  (i-iblk, x_nnod,x_inod)
        i_4=mod(iloop_end-iloop_sta+1,4)+iloop_sta

!$OMP MASTER
        ds = eigen_get_wtime()
        nodes(0) = eigen_owner_node(i, y_nnod, y_inod)
        if (nodes(0) == y_inod) then
          i_1 = eigen_translate_g2l(i, y_nnod, y_inod)
          do j_1=jloop_sta,jloop_end
            v1(j_1,1) = a(j_1,i_1)
          end do                ! j_1
        end if
        call bcast_dbl(v1(jloop_sta,1), jloop_end-jloop_sta+1,
     &       nodes(0), 1, y_COMM_WORLD)

        de = eigen_get_wtime()
        dcom = dcom + (de - ds)
!$OMP END MASTER

!$OMP BARRIER


#define	USE_BLAS	1

!$OMP MASTER
        call barrier(y_COMM_WORLD)
        call barrier(x_COMM_WORLD)
!$OMP END MASTER
        if (jloop_end >= jloop_sta) then
          i_4=mod(iloop_end-iloop_sta+1,8)+iloop_sta
          if (i_4 > iloop_sta) then
!$OMP DO SCHEDULE(STATIC)
            do i_1=iloop_sta,i_4-1
#if USE_BLAS
              ss1(i_1+0)=ddot(jloop_end-jloop_sta+1, z(jloop_sta,i_1),1,
     &             v1(jloop_sta,1),1)
#else
              s0 = ZERO
              do j_1=jloop_sta,jloop_end
                s0=s0+v1(j_1,1)*z(j_1,i_1+0)
              end do            ! j_1
              ss1(i_1+0)=s0
#endif
            end do
!$OMP END DO
          end if
!$OMP DO SCHEDULE(STATIC)
#if USE_BLAS
          do i_1=i_4,iloop_end,8
            call dgemv('T', jloop_end-jloop_sta+1, 8,
     &           ONE, z(jloop_sta,i_1), nmz, v1(jloop_sta,1), 1,
     &           ZERO, ss1(i_1), 1)
          end do                ! i_1
#else
          do i_1=i_4,iloop_end,4
            s0 = ZERO
            s1 = ZERO
            s2 = ZERO
            s3 = ZERO
            do j_1=jloop_sta,jloop_end
              s0=s0+v1(j_1,1)*z(j_1,i_1+0)
              s1=s1+v1(j_1,1)*z(j_1,i_1+1)
              s2=s2+v1(j_1,1)*z(j_1,i_1+2)
              s3=s3+v1(j_1,1)*z(j_1,i_1+3)
            end do              ! j_1
            ss1(i_1+0)=s0
            ss1(i_1+1)=s1
            ss1(i_1+2)=s2
            ss1(i_1+3)=s3
          end do                ! i_1
#endif
!$OMP END DO
        else
!$OMP DO
          do i_1=iloop_sta,iloop_end
            ss1(i_1) = ZERO
          end do
!$OMP END DO
        end if

!$OMP BARRIER

!$OMP MASTER
        ds = eigen_get_wtime()
        call reduce_dbl(ss1(iloop_sta),tt, iloop_end-iloop_sta+1, 2, 
     &       x_COMM_WORLD)
        de = eigen_get_wtime()
        dcom = dcom + (de - ds)

#if USE_BLAS
        call dscal(iloop_end-iloop_sta+1, d(i), ss1(iloop_sta), 1)
#else
        s0 = d(i)
        do i_1=iloop_sta,iloop_end
          ss1(i_1) = ss1(i_1) * s0
        end do                  ! i_1
#endif
!$OMP END MASTER

!$OMP BARRIER

        jloop_sta = eigen_loop_start(1,      x_nnod,x_inod)
        jloop_end = eigen_loop_end  (i-iblk, x_nnod,x_inod)
        i_4=mod(iloop_end-iloop_sta+1,8)+iloop_sta

        if (jloop_end >= jloop_sta) then
          if (i_4 > iloop_sta) then
!$OMP DO SCHEDULE(STATIC)
            do i_1=iloop_sta, i_4-1
#if USE_BLAS
              call daxpy(jloop_end-jloop_sta+1,
     &             ss1(i_1), v1(jloop_sta,1), 1, z(jloop_sta,i_1), 1)
#else
              s0 = ss1(i_1+0)
              do j_1=jloop_sta,jloop_end
                z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
              end do            ! j_1
#endif
            end do
!$OMP END DO
          end if
!$OMP DO SCHEDULE(STATIC)
#if USE_BLAS
          do i_1=i_4,iloop_end,8
            call dger(jloop_end-jloop_sta+1, 8,
     &           ONE, v1(jloop_sta,1), 1, 
     &           ss1(i_1), 1, z(jloop_sta,i_1), nmz)
          end do                ! i_1
#else
          do i_1=i_4,iloop_end,4
            s0 = ss1(i_1+0)
            s1 = ss1(i_1+1)
            s2 = ss1(i_1+2)
            s3 = ss1(i_1+3)
            do j_1=jloop_sta,jloop_end
              z(j_1,i_1+0) = z(j_1,i_1+0) + s0 * v1(j_1,1)
              z(j_1,i_1+1) = z(j_1,i_1+1) + s1 * v1(j_1,1)
              z(j_1,i_1+2) = z(j_1,i_1+2) + s2 * v1(j_1,1)
              z(j_1,i_1+3) = z(j_1,i_1+3) + s3 * v1(j_1,1)
            end do              ! j_1
          end do                ! i_1
#endif
!$OMP ENDDO
        end if

!$OMP BARRIER

      end do

!----------------------------------------------------------------

      d2 = eigen_get_wtime()
      d0 = d2 -d1

#define	USE_BCASTW	1

      mode = 0; i = nx+1

!$OMP MASTER
      if (nx+1 <= n) then
        ds = eigen_get_wtime()
        call trbakwy_datacast(iloop_end, m, i, a(1,1), nma,
     &       v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     &       wk(1),
#endif
     &       iblk)
        de = eigen_get_wtime()
        dcom = dcom + (de - ds)
      end if
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
      if (nx+1+m <= n) then
        ds = eigen_get_wtime()
        call trbakwy_datacast(iloop_end, m, i+m, a(1,1), nma,
     &       v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     &       wk(1),
#endif
     &       iblk)
        de = eigen_get_wtime()
        dcom = dcom + (de - ds)
      end if
!$OMP END MASTER

      do i=nx+1, n, m

        if (mode == 0) then
          call eigen_trbakwy_block_body(iloop_end, z(1,1), nmz,
     &         d(1), e(1), v1(1,1), nm, m,
     &         i, ss1(1), tt, nss, iblk,
     &         dcom, dx, dy, dz)
        end if
        if (mode == 1) then
          call eigen_trbakwy_block_body(iloop_end, z(1,1), nmz,
     &         d(1), e(1), v2(1,1), nm, m,
     &         i, ss2(1), tt, nss, iblk,
     &         dcom, dx, dy, dz)
        end if
        if (mode == 2) then
          call eigen_trbakwy_block_body(iloop_end, z(1,1), nmz,
     &         d(1), e(1), v3(1,1), nm, m,
     &         i, ss3(1), tt, nss, iblk,
     &         dcom, dx, dy, dz)
        end if

!$OMP MASTER
        if (i+2*m <= n) then
          ds = eigen_get_wtime()
          if (mode == 0) then
            call trbakwy_datacast(iloop_end, m, i+2*m, a(1,1), nma,
     &           v3(1,1), nm, ss3(1),
#if !USE_BCASTW
     &           wk(1),
#endif
     &           iblk)
          end if
          if (mode == 1) then
            call trbakwy_datacast(iloop_end, m, i+2*m, a(1,1), nma,
     &           v1(1,1), nm, ss1(1),
#if !USE_BCASTW
     &           wk(1),
#endif
     &           iblk)
          end if
          if (mode == 2) then
            call trbakwy_datacast(iloop_end, m, i+2*m, a(1,1), nma,
     &           v2(1,1), nm, ss2(1),
#if !USE_BCASTW
     &           wk(1),
#endif
     &           iblk)
          end if
          de = eigen_get_wtime()
          dcom = dcom + (de - ds)
        end if
!$OMP END MASTER

        mode = mod(mode + 1, 3)

      end do

!----------------------------------------------------------------

!$    if (omp_get_num_threads() > 1) then
!$OMP BARRIER
!$       if (omp_get_thread_num() == 1) then
!$          ss1(1) = dx; ss1(2) = dy
!$       end if
!$OMP BARRIER
!$       if (omp_get_thread_num() == 0) then
!$          dx = ss1(1); dy = ss1(2)
!$       end if
!$    end if
!$OMP BARRIER
!     
#if AT_BCAST_OVERLAP
!===  
!$OMP MASTER
!$    call sync_other_than_master_finalize(TRBK_LOCK)
!$OMP END MASTER
!===  
#endif
!$OMP MASTER
      call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER
!$OMP MASTER
      deallocate(wk)

#if TIMER_PRINT>1
      d2 = eigen_get_wtime()

      if (1 == TRD_inod) then
        print*,"TRBAK=",(d2-d1)
        print*,"COMM=",dcom
        print*,"   ",d2-d1,(2d0*nvec*n*n)/(d2-d1)*1d-9,"GFLOPS"
        if (dx > 0)
     &       print*,"   ",dx,(1d0*nvec*n*n)/(dx)*1d-9,"GFLOPS"
        if (dy > 0)
     &       print*,"   ",dy,(1d0*nvec*n*n)/(dy)*1d-9,"GFLOPS"
        print*,"   ",d0,dz
      end if
#endif
!$OMP END MASTER

      return

      end subroutine eigen_trbakwy_body

!------------------------------------------------------------------
!------------------------------------------------------------------

      subroutine trbakwy_datacast(iloop_end, m, i, a, nma, v, nm,
     &     ss,
#if !USE_BCASTW
     &     wk,
#endif
     &     iblk)

      integer, intent(in)   :: iloop_end
      integer, intent(in)   :: m
      integer, intent(in)   :: i
      real(8), intent(in)   :: a(nma,*)
      integer, intent(in)   :: nma
      real(8), intent(out)  :: v(nm,*)
      integer, intent(in)   :: nm
      real(8), intent(out)  :: ss(*)
#if !USE_BCASTW
      real(8), intent(out)  :: wk(*)
#endif
      integer, intent(in)   :: iblk

      integer               :: j, iy, k0
      integer               :: i_1, j_1, j_4, j_5
      integer               :: jloop_sta, jloop_end
      integer               :: nodes(0:nsm-1)

!$    real(8)               :: ds, de


!$    ds = omp_get_wtime()

      ss(1:iloop_end*m+ns0) = ZERO

#if USE_BCASTW
      jloop_end = min(eigen_loop_end(i+m-1-iblk, x_nnod, x_inod), nm)

      if (mod(m, y_nnod) == 0) then

        do j=0,m-1
          if (y_inod == eigen_owner_node(i+j, y_nnod, y_inod)) then
            iy  = eigen_owner_node(i+j, y_nnod, y_inod)
            i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
            k0  = (j/y_nnod)*jloop_end
            v(1+k0:jloop_end+k0, 1) = a(1:jloop_end, i_1)
          end if
        end do

        k0  = (m/y_nnod)*jloop_end
        call allgather_dbl(v, TRBK_buf, k0, 1, y_COMM_WORLD)

        do j=0,m-1
          iy = eigen_owner_node(i+j, y_nnod, y_inod)
          k0 = ((iy-1)*(m/y_nnod)+(j/y_nnod))*jloop_end
          v(1:jloop_end, j+1) = TRBK_buf(1+k0:jloop_end+k0)
          jloop_sta = eigen_loop_start(i+j,  x_nnod,x_inod)
          if (jloop_sta <= jloop_end)
     &         v(jloop_sta:jloop_end, j+1) = ZERO
        end do

      else

        do j=0,m-1
          if (y_inod == eigen_owner_node(i+j, y_nnod, y_inod)) then
            iy  = eigen_owner_node(i+j, y_nnod, y_inod)
            i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
            v(1:jloop_end, j+1) = a(1:jloop_end, i_1)
            jloop_sta = eigen_loop_start(i+j,  x_nnod,x_inod)
            if (jloop_sta <= jloop_end)
     &           v(jloop_sta:jloop_end, j+1) = ZERO
          end if
        end do
        do j=0,m-1,y_nnod
          iy  = eigen_owner_node(i+j, y_nnod, y_inod)
          k0  = min(m-j, y_nnod)
          call bcastw_dbl(v(1,j+1), jloop_end, iy,
     &         nm, k0, TRBK_buf, 2, y_COMM_WORLD)
        end do

      end if
#else

      if (m > y_nnod .and. y_nnod > 1) then

        do j=0,m-1
          nodes(j) = eigen_owner_node(i+j, y_nnod, y_inod)
        end do

        do iy=1,y_nnod

          jloop_sta = eigen_loop_start(1,         x_nnod,x_inod)
          jloop_end = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)

          k0=0
          do j=0,m-1
            if (nodes(j) == iy) then
              i_1 = eigen_translate_g2l(i+j, y_nnod, iy)
              if (y_inod == iy) then
                do j_1=jloop_sta,jloop_end
                  wk(k0+j_1) = a(j_1, i_1)
                end do          ! k
              end if
              k0=k0+(jloop_end-jloop_sta+1)
            end if
          end do

          call bcast_dbl(wk, k0, iy, 2, y_COMM_WORLD)

          k0=0
          do j=0,m-1
            if (nodes(j) == iy) then
              do j_1=jloop_sta,jloop_end
                v(j_1, j+1) = wk(k0+j_1)
              end do            ! k
              k0=k0+(jloop_end-jloop_sta+1)
            end if
          end do

          do j=0,m-1
            if (nodes(j) == iy) then
              j_4 = eigen_loop_start(i+j,       x_nnod,x_inod)
              j_5 = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
              do j_1=j_4,j_5
                v(j_1, j+1) = ZERO
              end do
            end if
          end do

        end do

      else

        do j=0,m-1
          nodes(j) = eigen_owner_node(i+j, y_nnod, y_inod)
          if (nodes(j) == y_inod) then
            i_1 = eigen_translate_g2l(i+j, y_nnod, y_inod)
            jloop_sta = eigen_loop_start(1,         x_nnod,x_inod)
            jloop_end = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
            do j_1=jloop_sta,jloop_end
              v(j_1, j+1) = a(j_1, i_1)
            end do              ! k
            jloop_sta = eigen_loop_start(i+j,       x_nnod,x_inod)
            jloop_end = eigen_loop_end  (i+m-1-iblk,x_nnod,x_inod)
!     CDIR NOVECTOR
            do j_1=jloop_sta,jloop_end
              v(j_1, j+1) = ZERO
            end do
          end if
        end do
        jloop_sta = eigen_loop_start(1,          x_nnod,x_inod)
        jloop_end = eigen_loop_end  (i+m-1-iblk, x_nnod,x_inod)
        do j=0,m-1
          call bcast_dbl(v(1,j+1), jloop_end-jloop_sta+1,
     &         nodes(j), 2, y_COMM_WORLD)
        end do

      end if

#endif

!$    de = omp_get_wtime()
!$    TRBK_TIME_BCAST = de - ds

      return

      end subroutine

      end module trbakwy4_mod
