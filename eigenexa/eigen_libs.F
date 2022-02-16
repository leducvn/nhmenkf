!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Module: eigen_libs_mod
!
! Purpose
! =======
!
! eigen_libs_mod manages the public modules for the EigenExa library
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

      module eigen_libs_mod

      use comm_mod
      use eigen_blacs_mod
      use eigen_devel_mod
      use CSTAB_mod
      use mpi
!$    use omp_lib

      implicit none   
      private

      character(32)   ::  CodeNAME = 'EigenExa'
      type, public :: version_t
        integer       :: Major_Version ! Major version
        integer       :: Minor_Version ! Minor version
        integer       :: Patch_Level   ! Patchlevel 0=none,1=a,2=b,...
        character(32) :: date          ! Release date
        character(32) :: vcode         ! Version code name
      end type  version_t

      type(version_t) :: Eigen_Version ! 2.4b
     &     = version_t (
     &     2, 4, 2,                    ! Major, Minor, Patchlevel
     &     'August 20, 2018',          ! Release date
     &     'Akashi'                    ! Version code
     &     )

      integer, public, parameter ::  eigen_NB    = 64

      interface

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_sx
!
! Purpose
! =======
!
! Main driver for solving eigenproblem of a symmetric matrix via
! transformation of a penta-diagonal format.
!
!
! Arguments
! =========
!
! n, lda, ldz    (input) integer
!  lda and ldz are the leading dimensions corresponding to a and z,
!  respectively. These values must be equal or larger than one obtained 
!  by calling eigen_get_matdims.
!
! nvec  (input) integer, optional
!  the number of eigenvectors computed
!  default case, it is set as n
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
     &        m_forward, m_backward, mode)
      integer,   intent(in)           :: n
      integer,   intent(in), optional :: nvec
      real(8),   intent(inout)        :: a(1:lda, *)
      integer,   intent(in)           :: lda
      real(8),   intent(out)          :: w(*)
      real(8),   intent(out)          :: z(1:ldz, *)
      integer,   intent(in)           :: ldz
      integer,   intent(in), optional :: m_forward
      integer,   intent(in), optional :: m_backward
      character, intent(in), optional :: mode
      end subroutine eigen_sx

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_s
!
! Purpose
! =======
!
! Main driver for solving eigenproblem of a symmetric matrix via
! transformation of a tridiagonal format.
!
!
! Arguments
! =========
!
! n, lda, ldz    (input) integer
!  lda and ldz are the leading dimensions corresponding to a and z,
!  respectively. These values must be equal or larger than one obtained
!  by calling eigen_get_matdims.
!
! nvec  (input) integer, optional
!  the number of eigenvectors computed
!  default case, it is set as n
!
! m_forward, m_backward  (input) integer, optional
!  default values are m_forward = 48, and m_backward = 128.
!
! mode  (input) character, optional
!  = 'A' : eigenvalues and corresponding eigenvectors (default)
!  = 'N' : eigenvalues only
!  = 'X' : mode 'A' + improve accuracy of eigenvalues
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_s(n, nvec, a, lda, w, z, ldz,
     &        m_forward, m_backward, mode)
      integer,   intent(in)           :: n
      integer,   intent(in), optional :: nvec
      real(8),   intent(inout)        :: a(1:lda, *)
      integer,   intent(in)           :: lda
      real(8),   intent(out)          :: w(*)
      real(8),   intent(out)          :: z(1:ldz, *)
      integer,   intent(in)           :: ldz
      integer,   intent(in), optional :: m_forward
      integer,   intent(in), optional :: m_backward
      character, intent(in), optional :: mode
      end subroutine  eigen_s

      end interface

      public :: eigen_sx
      public :: eigen_s
      public :: eigen_get_version
      public :: eigen_show_version
      public :: eigen_initialized
      public :: eigen_init
      public :: eigen_free
      public :: eigen_get_matdims
      public :: eigen_memory_internal
      public :: eigen_get_procs
      public :: eigen_get_id
      public :: eigen_get_comm
      public :: eigen_loop_start
      public :: eigen_loop_end
      public :: eigen_translate_l2g
      public :: eigen_translate_g2l
      public :: eigen_owner_node
      public :: eigen_owner_index
      public :: eigen_vector_zeropad_x
      public :: eigen_vector_zeropad_y
      public :: get_constant_eps
      public :: get_constant_pai
      public :: get_constant_2pai
      public :: get_constant_pai_2

      contains

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_version
!
! Purpose
! =======
!
! Get version information
!
!
! Arguments
! =========
!
! version    (output) integer
!            current version number
!
! date       (output) character*(*), optional
!            date string
!
! vcode      (output) character*(*), optional
!            version code string
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_version(version, date, vcode)

      integer,       intent(out)           :: version
      character*(*), intent(out), optional :: date
      character*(*), intent(out), optional :: vcode


      version = Eigen_Version%Major_Version * 100
     &     + Eigen_Version%Minor_Version * 10
     &     + Eigen_Version%Patch_Level

      if (present(date)) then
        date = Eigen_Version%date
      end if

      if (present(vcode)) then
        vcode = Eigen_Version%vcode
      end if

      return

      end subroutine eigen_get_version

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_show_version
!
! Purpose
! =======
!
! Show version information
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_show_version()

      character*(256) :: version
      character*(1  ) :: patchlevel
      integer         :: i


      i = min(26, Eigen_Version%Patch_Level) + 1
      patchlevel = (" abcdefghijklmnopqrstuvwxyz*" (i:i))

      write(version, '(I1,A,I1,A)')
     &     Eigen_Version%Major_Version,
     &     '.',Eigen_Version%Minor_Version, trim(patchlevel)

      if (TRD_inod == 1) then
        print*, "## EigenExa version (", trim(version),
     &       ") / (", trim(Eigen_Version%date),
     &       ") / (", trim(Eigen_Version%vcode), ")"
      end if

      return

      end subroutine eigen_show_version

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_initialized
!
! Purpose
! =======
!
! Get the flag that means initialized.
!
!
! Arguments
! =========
!
! flag   (output) logical
!        flag that means whether or not initialized.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_initialized(flag)

      logical, intent(out) :: flag


      call eigen_get_initialized(flag)

      return

      end subroutine eigen_initialized

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_init
!
! Purpose
! =======
!
! Initialize the utility of the EigenExa library.
!
!
! Arguments
! =========
!
! comm   (input) integer, optional
!        communicator
!
! order  (input) character*(*), optional
!        Grid major
!        'R': row-major
!        'C': column-major
!
! scalapack_context  (input) integer, optional
!         This option will be supported in future version.
!
! gridmap (input) integer, optional
!         This option will be supported in a future version.
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_init(comm, order,
     &     scalapack_context, gridmap)

      integer,       intent(in), optional :: comm
      character*(*), intent(in), optional :: order
      integer,       intent(in), optional :: scalapack_context
      integer,       intent(in), optional :: gridmap(:,:)

      logical          :: flag
      integer          :: n1, n2, n3, i, j, k, ierr
      integer          :: topo_type, cart_dim, dims(2), coords(2)
      logical          :: periods(2)
      integer          :: old_grp, new_grp

      integer          :: comm0
      integer          :: group0, group1, kk0(1), kk1(1)

      integer          :: local_size
!$    integer          :: th0(2), th1(2)

      real(8)          :: s1, s2, d1, d2, buff(4)
!      real(8), pointer :: buff1(:), buff2(:)
      real(8), allocatable :: buff1(:), buff2(:)

      character*1      :: GRID_major = 'C'

!--------------------

      call eigen_timer_reset(0, 0, 0, 0)

!--------------------

      call eigen_initialized(flag)

      if (flag) then
        if (TRD_inod == 1) then
          print*,"*************"
          print*,"** CAUTION **"
          print*,"*************"
          print*,"You are going to initialize EigenExa,"
          print*,"while EigenExa was not freed at last call."
          print*,"EigenExa restarts again by itself."
        end if
        call eigen_free()
      end if

!---- Comm setup ---

      call eigen_init_comm_setup()

!---- OpenMP thread check ---

      call eigen_init_omp_setup()

!---- Cartesian check ---

      call eigen_init_cartesian_check()

!---- BLACS setup ---

      call eigen_init_blacs_setup()

!---- Preliminary Sampling Collective comms. ---

      call eigen_init_collective_comms()

!---- Final setup ---

      call eigen_set_initialized()

!--------------------

      return

      contains

!--------------------

      subroutine eigen_init_comm_setup()

      if (present(comm)) then
        comm0 = comm
      else
        comm0 = MPI_COMM_WORLD
      end if
      call MPI_Comm_Test_Inter(comm0, flag, ierr)
      if (flag .or. ierr /= MPI_SUCCESS) then
        print*,"*************"
        print*,"** CAUTION **"
        print*,"*************"
        print*,"You are going to initialize EigenExa with"
        print*,"an invalid communicator."
        print*,"EigenExa terminates this run."
        call sleep(1)
        call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
      end if
      if (comm0 == MPI_COMM_NULL) then
        TRD_COMM_WORLD = MPI_COMM_NULL
        x_COMM_WORLD   = MPI_COMM_NULL
        y_COMM_WORLD   = MPI_COMM_NULL
        z_COMM_WORLD   = MPI_COMM_NULL
        w_COMM_WORLD   = MPI_COMM_NULL
        TRD_nnod = 0; TRD_inod = 0
        x_nnod   = 0; x_inod   = 0
        y_nnod   = 0; y_inod   = 0
        z_inod   = 0; z_nnod   = 0
        w_inod   = 0; w_nnod   = 0
      else
        call MPI_Comm_dup(comm0, TRD_COMM_WORLD, ierr)
      end if

      return

      end subroutine eigen_init_comm_setup

!--------------------

      subroutine eigen_init_omp_setup()

!$    if (TRD_COMM_WORLD /= MPI_COMM_NULL) then
!$      call MPI_Query_thread(MPI_THREAD_MODE, ierr)
!$      local_size = 1
!$OMP PARALLEL
!$OMP MASTER
!$      local_size = omp_get_num_threads()
!$OMP END MASTER
!$OMP END PARALLEL
!$      th0(1) = local_size
!$      th0(2) = -th0(1)
!$      call MPI_Allreduce(th0, th1, 2, MPI_INTEGER, MPI_MAX,
!$   &        TRD_COMM_WORLD, ierr)
!$      j = th1(1) + th1(2)
!$      if (j /= 0) then
!$        call MPI_Barrier(TRD_COMM_WORLD, ierr)
!$        call flush(6)
!$        if (TRD_inod == 1) then
!$          print*,"*************"
!$          print*,"** CAUTION **"
!$          print*,"*************"
!$          print*,"EigenExa supports only homogeneous thread"
!$   &           // " setting !"
!$          print*,"EigenExa terminates this run."
!$          print*,"*************"
!$        end if
!$        do j=1,2
!$          call MPI_Barrier(TRD_COMM_WORLD, ierr)
!$          call flush(6)
!$        end do
!$        call MPI_Barrier(TRD_COMM_WORLD, ierr)
!$        call sleep(1)
!$        call MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER, ierr)
!$      end if
!$    end if

      return

      end subroutine eigen_init_omp_setup

!--------------------

      subroutine eigen_init_cartesian_check()

      if (TRD_COMM_WORLD /= MPI_COMM_NULL) then

!---- Cartesian check ---

        call MPI_Comm_size(TRD_COMM_WORLD, TRD_nnod, ierr)
        call MPI_Comm_rank(TRD_COMM_WORLD, TRD_inod, ierr)
        TRD_inod = TRD_inod+1

        call MPI_Topo_test(TRD_COMM_WORLD, topo_type, ierr)
        if (topo_type == MPI_CART) then
          call MPI_Cartdim_get(TRD_COMM_WORLD, cart_dim, ierr)
        else
          cart_dim = 1
        end if

!---- 2D Cartesian setup ---

        if (cart_dim == 2) then

          call MPI_Cart_get(TRD_COMM_WORLD, cart_dim,
     &         dims, periods, coords, ierr)
          x_nnod = dims(1)
          y_nnod = dims(2)
          x_inod = coords(1) +1 ! fixed
          y_inod = coords(2) +1 ! fixed

          GRID_major = 'R'
          if (present(order)) then
            if (x_inod == 1 .and. y_inod == 1) then
              if (order(1:1) == 'C' .or. order(1:1) == 'c') then
                print*,"*************"
                print*,"** CAUTION **"
                print*,"*************"
                print*,"The MPI_CART you specified is based on"
     &               // " R-major,"
                print*,"but you also specified C-major option."
                print*,"EigenExa solve this conflict by taking"
     &               // " R-major,"
              end if
            end if
          end if

        else                    ! if (cart_dim == 2

!---- Setup 2D process map ---

          x_nnod = int(sqrt(dble(TRD_nnod)))
          i = 1                 ! minimum factor, x_nnod must be
!     multiple of k
          if (mod(TRD_nnod, i) == 0) then
            k = i
          else
            k = 1
          end if
          do
            if (x_nnod <= k) exit
            if (mod(x_nnod, k) == 0 .and.
     &          mod(TRD_nnod, x_nnod) == 0) exit
            x_nnod = x_nnod-1
          end do                !!
          y_nnod = TRD_nnod/x_nnod

          if (present(order)) then
            GRID_major = order(1:1)
          else
            GRID_major = 'C'
          end if

          if (GRID_major == 'R' .or. GRID_major == 'r') then
            GRID_major = 'R'
          else
            GRID_major = 'C'
          end if

          call eigen_set_grid_major(GRID_major)

          if (GRID_major == 'R') then
!     row-major
            x_inod =    (TRD_inod-1)/y_nnod +1
            y_inod = mod(TRD_inod-1, y_nnod)+1
          else
!     column-major
!     ** EigenK adopts column-major in default
!     The process ordering on the Earth Simulator was done in the way
!     such as incrmental from internal processors to external nodes.
!     We want to make the processes, which possess a distributed vector,
!     close physically on the network connection.
!     These are the historical reason to adopt column-major.
            x_inod = mod(TRD_inod-1, x_nnod)+1
            y_inod =    (TRD_inod-1)/x_nnod +1
          end if

        end if                  ! if (cart_dim == 2) then

!---- x-Comm & y-Comm setup ---

        call MPI_Comm_split(TRD_COMM_WORLD, y_inod, x_inod,
     &       x_COMM_WORLD, ierr)
        call MPI_Comm_split(TRD_COMM_WORLD, x_inod, y_inod,
     &       y_COMM_WORLD, ierr)

!---- Parameter setup for datacast_dbl (row-column vector redistribution) ---

        n1 = max(x_nnod, y_nnod)
        n2 = min(x_nnod, y_nnod)
        do
          if (n1 == n2) then
            n_common = n1
            exit
          end if
          n3 = n1-n2
          n1 = max(n2, n3)
          n2 = min(n2, n3)
        end do                  !!

        if (associated(p0_)) then
          deallocate(p0_)
        end if
        allocate(p0_(1:max(x_nnod, y_nnod)))

        if (associated(q0_)) then
          deallocate(q0_)
        end if
        allocate(q0_(1:max(x_nnod, y_nnod)))

        p0_(:) = -1
        q0_(:) = -1
        do i = 1, x_nnod
          if (mod(i-1, n_common) == mod(y_inod-1, n_common)) then
            n1 = y_inod-i
            if (n1 >= 0) then
              do j = 1, x_nnod
                k = +n1+(j-1)*y_nnod
                if (mod(k, x_nnod) == 0) then
                  p0_(i) = k/x_nnod
                  q0_(i) = (j-1)
                  exit
                end if
              end do            ! j
            else
              do j = 1, y_nnod
                k = -n1+(j-1)*x_nnod
                if (mod(k, y_nnod) == 0) then
                  q0_(i) = k/y_nnod
                  p0_(i) = (j-1)
                  exit
                end if
              end do            ! j
            end if
          end if
        end do                  ! i
        p0_(:) = p0_(:)+1
        q0_(:) = q0_(:)+1


        diag_0 = 0
        diag_1 = 0
        do i = 1, y_nnod/n_common
          j = (i-1)*y_nnod+y_inod
          k = mod(j-1, x_nnod)+1
          if (k == x_inod) then
            diag_0 = i
            diag_1 = (j-1)/x_nnod+1
            exit
          end if
        end do                  ! i_1

!---- extra-Comm setup for datacast_dbl ---

        if (n_common > 1) then
          w_inod =     (x_inod-1)/ n_common +1
          z_inod = mod((x_inod-1), n_common)+1
          w_nnod = x_nnod/n_common
          z_nnod = n_common
          call MPI_Comm_split(x_COMM_WORLD, w_inod-1, z_inod-1,
     &         z_COMM_WORLD, ierr)
          call MPI_Comm_split(x_COMM_WORLD, z_inod-1, w_inod-1,
     &         w_COMM_WORLD, ierr)
        else
          z_COMM_WORLD = MPI_COMM_SELF
          w_COMM_WORLD = x_COMM_WORLD
          z_inod = 1; z_nnod = 1
          w_inod = x_inod; w_nnod = x_nnod
        end if

!---- Exceptional setup ---

      else                      !  if (TRD_COMM_WORLD /= MPI_COMM_NULL

        x_COMM_WORLD = MPI_COMM_NULL
        y_COMM_WORLD = MPI_COMM_NULL
        z_COMM_WORLD = MPI_COMM_NULL
        w_COMM_WORLD = MPI_COMM_NULL
        x_nnod = 0; x_inod = 0
        y_nnod = 0; y_inod = 0
        z_inod = 0; z_nnod = 0
        w_inod = 0; w_nnod = 0

        GRID_major = 'C'
        call eigen_set_grid_major(GRID_major)

      end if                    !  if (TRD_COMM_WORLD /= MPI_COMM_NULL

      return

      end subroutine eigen_init_cartesian_check

!--------------------

      subroutine eigen_init_blacs_setup()

#if TIMER_PRINT>1
      if (TRD_inod <= 1) then
        print*,"GRID major ",GRID_major," is specified."
      end if
#endif

      if (TRD_COMM_WORLD == MPI_COMM_NULL) then
        comm0  = MPI_COMM_SELF
        x_nnod = 1; x_inod = 1
        y_nnod = 1; y_inod = 1
      else
        comm0  = TRD_COMM_WORLD
      end if

      call Eigen_BLACS_Init(comm0, x_nnod, y_nnod, GRID_Major)

      if (TRD_COMM_WORLD == MPI_COMM_NULL) then
        x_nnod = 0; x_inod = 0
        y_nnod = 0; y_inod = 0
      end if

      return

      end subroutine eigen_init_blacs_setup

!--------------------

      subroutine eigen_init_collective_comms()

      if (TRD_COMM_WORLD /= MPI_COMM_NULL) then

!---- Measure the overhead of collective communications ---

        call MPI_Barrier(TRD_COMM_WORLD,ierr)
        call MPI_Barrier(TRD_COMM_WORLD,ierr)

!$OMP PARALLEL
!$OMP+          PRIVATE(i,j,d1,d2,buff1,buff2)

        allocate(buff1(2048), buff2(2048))

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do  i = 1, 10
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        call MPI_Barrier(TRD_COMM_WORLD,ierr)
        call MPI_Barrier(TRD_COMM_WORLD,ierr)
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        call MPI_Barrier(TRD_COMM_WORLD,ierr)
        call MPI_Barrier(TRD_COMM_WORLD,ierr)

        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Barrier_Overhead_x = buff(3) / TRD_nnod
        Barrier_Overhead_y = buff(4) / TRD_nnod
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),2048,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),2048,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Reduce_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
        Reduce_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call reduce_dbl(buff2(1),buff1(1),1024,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do i = 1, 10
!$OMP MASTER
          call reduce_dbl(buff1(1),buff2(1),1024,0,x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call reduce_dbl(buff1(1),buff2(1),1024,0,y_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call reduce_dbl(buff2(1),buff1(1),1024,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Reduce_cont_Overhead_x = buff(3) / TRD_nnod
     &       - Barrier_Overhead_x - Reduce_Overhead_x
        Reduce_cont_Overhead_y = buff(4) / TRD_nnod
     &       - Barrier_Overhead_y - Reduce_Overhead_y
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),2048,1,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),2048,1,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &        MPI_SUM, TRD_COMM_WORLD, ierr)
        Bcast_Overhead_x = buff(3) / TRD_nnod - Barrier_Overhead_x
        Bcast_Overhead_y = buff(4) / TRD_nnod - Barrier_Overhead_y
!$OMP END MASTER

!--------------------

!$OMP MASTER
        s1 = ZERO
        s2 = ZERO
!$OMP END MASTER
!$OMP BARRIER

        do i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call bcast_dbl(buff2(1),1024,1,0,x_COMM_WORLD)
          call MPI_Barrier(x_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s1 = s1 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

        do  i = 1, 10
!$OMP MASTER
          call bcast_dbl(buff1(1),1024,1,0,x_COMM_WORLD)
!$OMP END MASTER
          continue
!$OMP MASTER
          call MPI_Barrier(TRD_COMM_WORLD,ierr)
          d1 = eigen_get_wtime()
          call bcast_dbl(buff1(1),1024,1,0,y_COMM_WORLD)
!$OMP END MASTER
!$OMP BARRIER
          continue
!$OMP MASTER
          call bcast_dbl(buff2(1),1024,1,0,y_COMM_WORLD)
          call MPI_Barrier(y_COMM_WORLD,ierr)
          d2 = eigen_get_wtime()
          if (i > 5) s2 = s2 + (d2-d1)
!$OMP END MASTER
!$OMP BARRIER
        end do

!$OMP MASTER
        buff(1) = s1 / 5
        buff(2) = s2 / 5
        call MPI_Allreduce(buff(1), buff(3), 2, MPI_DOUBLE_PRECISION,
     &       MPI_SUM, TRD_COMM_WORLD, ierr)
        Bcast_cont_Overhead_x = buff(3) / TRD_nnod
     &       - Barrier_Overhead_x - Bcast_Overhead_x
        Bcast_cont_Overhead_y = buff(4) / TRD_nnod
     &       - Barrier_Overhead_y - Bcast_Overhead_y
!$OMP END MASTER

        deallocate(buff1, buff2)

!$OMP END PARALLEL

!--------------------

        if (Reduce_cont_Overhead_x < ZERO) then
          Reduce_cont_Overhead_x = ZERO
        end if
        if (Reduce_cont_Overhead_y < ZERO) then
          Reduce_cont_Overhead_y = ZERO
        end if
        if (Bcast_cont_Overhead_x  < ZERO) then
          Bcast_cont_Overhead_x  = ZERO
        end if
        if (Bcast_cont_Overhead_y  < ZERO) then
          Bcast_cont_Overhead_y  = ZERO
        end if

#if TIMER_PRINT>0
        if (TRD_inod == 1) then
          print*,";----- MPI in-(x/y)-group overhead -----"
          print*,"Overhead Barrier(x) :: ", Barrier_Overhead_x
          print*,"Overhead Barrier(y) :: ", Barrier_Overhead_y
          print*,"Overhead Bcast(x)   :: ", Bcast_Overhead_x
          print*,"Overhead Bcast(y)   :: ", Bcast_Overhead_y
          print*,"Overhead Reduce(x)  :: ", Reduce_Overhead_x
          print*,"Overhead Reduce(y)  :: ", Reduce_Overhead_y
          print*,"Overhead Bcast_cont(x)   :: ",
     &         Bcast_cont_Overhead_x
          print*,"Overhead Bcast_cont(y)   :: ",
     &         Bcast_cont_Overhead_y
          print*,"Overhead Reduce_cont(x)  :: ",
     &         Reduce_cont_Overhead_x
          print*,"Overhead Reduce_cont(y)  :: ",
     &         Reduce_cont_Overhead_y
        end if
#endif

!--------------------

      end if

      end subroutine  eigen_init_collective_comms

      end subroutine  eigen_init

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_free
!
! Purpose
! =======
!
! Free the utility of the EigenExa library
!
!
! Arguments
! =========
!
! flag   (input) integer, optional
!        TRUE:  prit out the stats of timer table.
!        FALSE: 
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_free(flag)

      integer, intent(in), optional ::  flag

      integer                ::  ierr
      logical                ::  local_flag
      real(8)                ::  tmp


      call eigen_get_initialized(local_flag)

      if (.not. local_flag) then
        return
      end if

      if (TRD_COMM_WORLD /= MPI_COMM_NULL) then

        call MPI_Comm_free(x_COMM_WORLD, ierr)
        call MPI_Comm_free(y_COMM_WORLD, ierr)

        if (n_common > 1) then
          call MPI_Comm_free(z_COMM_WORLD, ierr)
          call MPI_Comm_free(w_COMM_WORLD, ierr)
        end if

        call MPI_Comm_free(TRD_COMM_WORLD, ierr)

        deallocate(p0_)
        deallocate(q0_)

        if (present(flag)) then
          if (flag == 1) then
            tmp = eigen_timer_print('EigenExa(finalized)')
          end if
        end if
        call eigen_timer_reset(0, 0, 0, 0)

      end if

!     call Eigen_BLACS_Exit()

      TRD_COMM_WORLD = MPI_COMM_WORLD

      call eigen_unset_initialized()

      return

      end subroutine eigen_free

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_matdims
!
! Purpose
! =======
!
! Returns the recommended array size for the distributed matrix,
! such as (1:nx, 1:ny).
!
!
! Arguments
! =========
!
! n      (input) integer
!        Global matrix dimension
!
! nx     (output) integer
!        Array size for the first index
!
! ny     (output) integer
!        Array size for the second index
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_matdims(n, nx, ny)

      integer, intent(in)    :: n
      integer, intent(out)   :: nx, ny

      integer :: NPROW, NPCOL, NB
      integer :: n1, nm, nmz, nmw, nn, larray

      integer :: lddw, lddz, nxx

      integer(8) :: LX1, LX2


      if (n <= 0) then
        nx = -1; ny = -1
        return
      end if

      NPROW = x_nnod
      NPCOL = y_nnod

      n1 = ((n-1)/NPROW+1)
      call CSTAB_get_optdim(n1, 6, 16*4, 16*4*2, nm)

      NB  = eigen_NB

      nmz = ((n-1)/NPROW+1)
      nmz = ((nmz-1)/NB+1)*NB+1
      nn  = nmz
      nmz = (n-1)/NB+1
      nmz = ((nmz-1)/NPROW+1)*NB
!     Fix on version 2.2b
!     to avoid unexpected SIGSEGV,
!     use the maximum of nn and nmz.
      nmz = max(nn, nmz)

      nmw = ((n-1)/NPCOL+1)
      nmw = ((nmw-1)/NB+1)*NB+1
      nn  = nmw
      nmw = (n-1)/NB+1
      nmw = ((nmw-1)/NPCOL+1)*NB
!     Fix on version 2.2b
!     to avoid unexpected SIGSEGV,
!     use the maximum of nn and nmw.
      nmw = max(nn, nmw)


      larray = max(nmz, nm)*nmw

      nx = nm
      ny = (larray-1)/nm+1

#if 0
!     
!     check whether intermediate datasize is valid with in integer(4).
!     
      NB = eigen_NB
      NB = min(NB, N)

      lddz = (n-1)/NPROW+1
      lddz = ((lddz-1)/NB+1)*NB
      lddw = (n-1)/NPCOL+1
      lddw = ((lddw-1)/NB+1)*NB
      nxx  = (N-1)/NPCOL+1

      LX1 = lddw; LX1 = LX1 * lddw
      LX2 = nxx;  LX2 = LX2 * n

      if (LX1 >= (2_8)**31 .OR. LX2 >= (2_8)**31) then
        print*, "Warning :: oversized problem !!"
        larray = -1
        nx     = -1
        ny     = -1
      end if
#endif

      return

      end subroutine  eigen_get_matdims

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function : eigen_memory_internal
!
! Purpose
! =======
!
! Returns size of the required memory internally while proceeding the
! EigenExa library.
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*
!
!     ******* ! important ! API has been changed ! *********************
!     Since the size of intermediate data may overfolow integer(4),
!     version 2.3b has changed the datatype of the function  result
!     from integer to integer(8).
!     ******************************************************************
!
      integer(8) function eigen_memory_internal(n, lda, ldz,
     &     m1_opt, m0_opt)
     &     result( byte )

      integer, intent(in)           :: n, lda, ldz
      integer, intent(in), optional :: m1_opt, m0_opt

      integer, parameter :: nm_max_L1 = 16*4
      integer, parameter :: nm_max_L2 = 16*6
      integer, parameter :: REAL_SIZE = 8_8 ! must be in integer(8)
      integer, parameter :: INT_SIZE  = 4_8 ! must be in integer(8)

      integer    :: m1, m0
      integer    :: MYROW, MYCOL, NPROW, NPCOL, NP, NQ, NB
      integer    :: local_rank, local_size
      integer    :: m , nv, nm, ierr

      integer    :: lddz, lddw, LWORK, LIWORK
      integer    :: na, nx, ny, nz, kx
      integer    :: l_array_1, l_array_2

      integer(8) :: nax(2)
      integer(8) :: byte_tridi
      integer(8) :: byte_dcx
      integer(8) :: byte_trbak

      integer, external  :: NUMROC

      include 'CSTAB.h'


      if (n <= 0) then
        byte = -1
        return
      end if

      local_rank = 0
      local_size = 1
!$    local_rank = omp_get_thread_num()
!$    local_size = omp_get_num_threads()

      if (present(m1_opt)) then
        m1 = m1_opt
      else
        m1 = 48
      end if

      if (present(m0_opt)) then
        m0 = m0_opt
      else
        m0 = 128
      end if

!     ===== for TRIDI ====

      nx = (n-1)/x_nnod+1 +2
      call CSTAB_get_optdim(nx, 6, nm_max_L1, nm_max_L2, nv)

      nx = (n-1)/x_nnod+1
      kx = y_nnod / n_common
      nz = ((nx-1)/kx+1)
      ny = nz*MAX(x_nnod,y_nnod)
!
      l_array_1 = MAX(MAX(m,2)*ny, nx+4*m+6)
      l_array_2 = MAX(MAX(m,2)*nz, nx+4*m+6)

      LWORK =
     &       (l_array_1+n_columns)          ! u_t
     &     + (l_array_2+n_columns)          ! v_t
     &     + MAX(3*m1, (nm*m1+n_columns))   ! w
     &     + (nv+n_columns)                 ! d_t
     &     + 2*(nv*2*m1+2*n_columns)        ! u_x, v_x, u_y, v_y
     &     + 4*(nv*local_size+n_columns)    ! u0_z, v0_z

      byte_tridi = LWORK * REAL_SIZE

!     ===== for DCX ====

      NB = eigen_NB
      NB = min(NB, N)

      NPROW = x_nnod
      NPCOL = y_nnod
      MYROW = x_inod-1
      MYCOL = y_inod-1
      NP = NUMROC(n, NB, MYROW, 0, NPROW)
      NQ = NUMROC(n, NB, MYCOL, 0, NPCOL)

      lddz = (n-1)/NPROW+1
      lddz = ((lddz-1)/NB+1)*NB
      lddw = (n-1)/NPCOL+1
      lddw = ((lddw-1)/NB+1)*NB

      nx     = (N-1)/NPCOL+1
      LWORK  = max(1+6*N+2*NP*(NQ+max(NQ,NB)), lddz*lddw, ldz*nx)
     &     + n
      LIWORK = 2+7*n+8*NPCOL

      byte_dcx = LWORK * REAL_SIZE + LIWORK * INT_SIZE

!     ===== for TRBAK ====

      m    = min(nsm, m0)
      if (m < 1) m = 1

      na   = (n-1)/y_nnod+1
      na   = na  + mod(na  -1,2)
      call CSTAB_get_optdim(lda, 9, 16*4, 16*6, nm)

      LWORK = n
     &     + 3*(max(nm*m,n)+n_columns)
     &     + 4*(na*m+ns0+n_columns)
     &     + m*nm
      LIWORK = max(((m-1)/y_nnod+1)*((n-1)/x_nnod+1),n)

      byte_trbak = LWORK * REAL_SIZE + LIWORK * INT_SIZE


      nax(1) = max(byte_tridi, byte_dcx, byte_trbak)
      call MPI_Allreduce(nax(1), nax(2), 1, MPI_INTEGER8,
     &     MPI_MAX, TRD_COMM_WORLD, ierr)
      byte = nax(2)

      return

      end function eigen_memory_internal

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_procs
!
! Purpose
! =======
!
! Returns the number of processes corresponding to three comms.
!
!
! Arguments
! =========
!
! procs  (output) integer
!        Num of Global communicator
!
! x_procs (output) integer
!        Num of x- communicator
!
! y_procs (output) integer
!        Num of y- communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_procs(procs, x_procs, y_procs)

      integer, intent(out)   ::  procs
      integer, intent(out)   ::  x_procs
      integer, intent(out)   ::  y_procs


      procs   = TRD_nnod
      x_procs = x_nnod
      y_procs = y_nnod

      return

      end subroutine  eigen_get_procs

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_id
!
! Purpose
! =======
!
! Returns process ID corresponding to three comms.
!
!
! Arguments
! =========
!
! id     (output) integer
!        Global ID
!
! x_id   (output) integer
!        x-ID
!
! y_id   (output) integer
!        y-ID
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine  eigen_get_id(id, x_id, y_id)

      integer, intent(out)   ::  id
      integer, intent(out)   ::  x_id
      integer, intent(out)   ::  y_id


      id   = TRD_inod
      x_id = x_inod
      y_id = y_inod

      return

      end subroutine  eigen_get_id

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_get_comm
!
! Purpose
! =======
!
! Returns communicators.
!
!
! Arguments
! =========
!
! comm   (output) integer
!        global communicator
!
! x_comm (output) integer
!        x-dir communicator, column-group
!
! y_comm (output) integer
!        y-dir communicator, row-group
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_get_comm(comm, x_comm, y_comm)

      integer, intent(out)   ::  comm
      integer, intent(out)   ::  x_comm
      integer, intent(out)   ::  y_comm


      comm   = TRD_COMM_WORLD
      x_comm = x_COMM_WORLD
      y_comm = y_COMM_WORLD

      return

      end subroutine eigen_get_comm

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_loop_start
!
! Purpose
! =======
!
! Return the local loop initial value of the specific global loop
! interval
!
!
! Arguments
! =========
!
! istart (input) integer
!        global loop initial value
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_loop_start(istart, nnod, inod)
     &     result(ret)

      integer, intent(in)    ::  istart
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod


      ret = (istart + nnod - 1 - inod) / nnod + 1

      return

      end function eigen_loop_start

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_loop_end
!
! Purpose
! =======
!
! Return the local loop temination value of the specific global loop
! interval
!
!
! Arguments
! =========
!
! iend   (input) integer
!        global loop termination value
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_loop_end(iend, nnod, inod)
     &     result(ret)

      integer, intent(in)    ::  iend
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = (iend + nnod - 0 - inod) / nnod + 0

      return

      end function eigen_loop_end

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_translate_l2g
!
! Purpose
! =======
!
! Translate the local index to the global index
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_translate_l2g(ictr, nnod, inod)
     &     result(ret)

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = (ictr-1)*nnod+inod

      return

      end function eigen_translate_l2g

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_translate_g2l
!
! Purpose
! =======
!
! Translate the global index to the local index
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_translate_g2l(ictr, nnod, inod)
     &     result(ret)

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = (ictr-1)/nnod+1

      return

      end function eigen_translate_g2l

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_owner_node
!
! Purpose
! =======
!
! Returns the owner ID of the specific global index
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_owner_node(ictr, nnod, inod)
     &     result(ret)

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      ret = mod(ictr-1, nnod)+1

      return

      end function eigen_owner_node

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: eigen_owner_index
!
! Purpose
! =======
!
! Returns the index if owner but -1 if not owner.
!
!
! Arguments
! =========
!
! ictr   (input) integer
!        global index
!
! nnod   (input) integer
!        The number of nodes in a communicator
!
! inod   (input) integer
!        ID of the node in a communicator
!
!--------*---------*---------*---------*---------*---------*---------*-*

      integer function eigen_owner_index(ictr, nnod, inod)
     &     result(ret)

      integer, intent(in)    ::  ictr
      integer, intent(in)    ::  nnod
      integer, intent(in)    ::  inod

      integer  ::  j_2, j_3

      j_2 = eigen_loop_start(ictr, nnod, inod)
      j_3 = eigen_loop_end  (ictr, nnod, inod)
      if (j_2 == j_3) then
        ret = j_2
      else
        ret = -1
      end if

      return

      end function eigen_owner_index

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_vector_zeropad_x
!
! Purpose
! =======
!
! Fill zeros in padding area of the working buffer
!
!
! Arguments
! =========
!
! v_x    (input/output) real(8) array, dimension(*)
!        working buffer
!
! L      (input) integer
!        Length of global vector
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_vector_zeropad_x(v_x, L)

      real(8), intent(inout) :: v_x(*)
      integer, intent(in)    :: L

      integer                :: j_3, n

      j_3 = eigen_loop_end     (L, x_nnod, x_inod)
      n   = eigen_translate_g2l(L, x_nnod, x_inod)
      if (j_3 < n) then
        v_x(j_3+1:n) = ZERO     ! in case
      end if

      end subroutine eigen_vector_zeropad_x

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Subroutine: eigen_vector_zeropad_y
!
! Purpose
! =======
!
! Fill zeros in padding area of the working buffer
!
!
! Arguments
! =========
!
! v_y    (input/output) real(8) array, dimension(*)
!        working buffer
!
! L      (input) integer
!        Length of global vector
!
!--------*---------*---------*---------*---------*---------*---------*-*

      subroutine eigen_vector_zeropad_y(v_y, L)

      real(8), intent(inout) :: v_y(*)
      integer, intent(in)    :: L
      integer                :: j_3, n

      j_3 = eigen_loop_end     (L, y_nnod, y_inod)
      n   = eigen_translate_g2l(L, y_nnod, y_inod)
      if (j_3 < n) then
        v_y(j_3+1:n) = ZERO     ! in case
      end if

      end subroutine eigen_vector_zeropad_y

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_eps
!
! Purpose
! =======
!
! Constant number machine_epsilon
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_eps()
     &     result(r)

      real(8) :: eps
      data eps /z'3CB0000000000000'/

      r = eps

      return

      end function get_constant_eps
   
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_pai
!
! Purpose
! =======
!
! Constant number PAI
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_pai()
     &     result(r)

      real(8) :: const_pai
      data const_pai /z'400921FB54442D18'/

      r = const_pai

      return

      end function get_constant_pai

!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_2pai
!
! Purpose
! =======
!
! Constant number 2*PAI
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_2pai()
     &     result(r)

      real(8) :: const_2pai
      data const_2pai /z'401921FB54442D18'/

      r = const_2pai

      return

      end function  get_constant_2pai
   
!--------*---------*---------*---------*---------*---------*---------*-*
!     
! Function: get_constant_pai_2
!
! Purpose
! =======
!
! Constant number PAI/2
!
!
! Arguments
! =========
!
!--------*---------*---------*---------*---------*---------*---------*-*

      real(8) function get_constant_pai_2()
     &     result(r)

      real(8) :: const_pai2
      data const_pai2 /z'3ff921FB54442D18'/

      r = const_pai2

      return

      end function get_constant_pai_2

      end module eigen_libs_mod
