module NodeMPI
! Author: Le Duc
! Created date: 27 Jul 2015
   use variable, only : r_size, r_sngl
   use NodeInfo_class
   use mpi
   implicit none
   !
   integer :: r_type, r_2type
   integer :: MPI_COMM_XWORLD, MPI_COMM_YWORLD, MPI_COMM_EWORLD, MPI_COMM_XYWORLD
   !
contains
   !
   subroutine initialize_mpi(info)
      implicit none
      type(NodeInfo), intent(in) :: info
      integer :: myid, myidx, myidy, myide, nproc, nxpe, nype, nepe, id, ierror
      !
      if (r_size == kind(0.0e0)) then
         r_type = MPI_REAL
         r_2type = MPI_2REAL
      else if (r_size == kind(0.0d0)) then
         r_type = MPI_DOUBLE_PRECISION
         r_2type = MPI_2DOUBLE_PRECISION
      end if
      !
      call get_myid(info, myid, nproc)
      call get_myidx(info, myidx)
      call get_myidy(info, myidy)
      call get_myide(info, myide)
      call get_nxpe(info, nxpe)
      call get_nype(info, nype)
      call get_nepe(info, nepe)
      id = myide*nype + myidy
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, id, myid, MPI_COMM_XWORLD, ierror)
      id = myide*nxpe + myidx
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, id, myid, MPI_COMM_YWORLD, ierror)
      id = myidy*nxpe + myidx
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, id, myid, MPI_COMM_EWORLD, ierror)
      id = myide
      call MPI_COMM_SPLIT(MPI_COMM_WORLD, id, myid, MPI_COMM_XYWORLD, ierror)
      !
      return
   end subroutine initialize_mpi
   !
   !
   !
   subroutine int_broadcast0D(header, datum, source)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source
      integer, intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, 1, MPI_INTEGER, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, 1, MPI_INTEGER, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, 1, MPI_INTEGER, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, 1, MPI_INTEGER, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, 1, MPI_INTEGER, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine int_broadcast0D
   !
   !
   !
   subroutine real_broadcast0D(header, datum, source)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source
      real, intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, 1, MPI_REAL, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, 1, MPI_REAL, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, 1, MPI_REAL, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, 1, MPI_REAL, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, 1, MPI_REAL, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine real_broadcast0D
   !
   !
   !
   subroutine broadcast0D(header, datum, source)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source
      real(kind=r_size), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, 1, r_type, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, 1, r_type, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, 1, r_type, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, 1, r_type, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, 1, r_type, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine broadcast0D
   !
   !
   !
   subroutine int_broadcast1D(header, datum, source, nx0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0
      integer, dimension(nx0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0, MPI_INTEGER, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0, MPI_INTEGER, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0, MPI_INTEGER, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0, MPI_INTEGER, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0, MPI_INTEGER, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine int_broadcast1D
   !
   !
   !
   subroutine real_broadcast1D(header, datum, source, nx0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0
      real, dimension(nx0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0, MPI_REAL, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0, MPI_REAL, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0, MPI_REAL, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0, MPI_REAL, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0, MPI_REAL, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine real_broadcast1D
   !
   !
   !
   subroutine broadcast1D(header, datum, source, nx0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0
      real(kind=r_size), dimension(nx0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0, r_type, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0, r_type, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0, r_type, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0, r_type, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0, r_type, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine broadcast1D
   !
   !
   !
   subroutine int_broadcast2D(header, datum, source, nx0, ny0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0, ny0
      integer, dimension(nx0,ny0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0*ny0, MPI_INTEGER, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0*ny0, MPI_INTEGER, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0*ny0, MPI_INTEGER, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0*ny0, MPI_INTEGER, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0*ny0, MPI_INTEGER, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine int_broadcast2D
   !
   !
   !
   subroutine broadcast2D(header, datum, source, nx0, ny0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0, ny0
      real(kind=r_size), dimension(nx0,ny0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0*ny0, r_type, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0*ny0, r_type, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0*ny0, r_type, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0*ny0, r_type, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0*ny0, r_type, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine broadcast2D
   !
   !
   !
   subroutine int_broadcast3D(header, datum, source, nx0, ny0, nz0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0, ny0, nz0
      integer, dimension(nx0,ny0,nz0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0*ny0*nz0, MPI_INTEGER, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0*ny0*nz0, MPI_INTEGER, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0*ny0*nz0, MPI_INTEGER, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0*ny0*nz0, MPI_INTEGER, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0*ny0*nz0, MPI_INTEGER, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine int_broadcast3D
   !
   !
   !
   subroutine broadcast3D(header, datum, source, nx0, ny0, nz0)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0, ny0, nz0
      real(kind=r_size), dimension(nx0,ny0,nz0), intent(inout) :: datum
      integer :: ierror
      !
      if (trim(header) == 'x') then
         call MPI_BCAST(datum, nx0*ny0*nz0, r_type, source, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
	 call MPI_BCAST(datum, nx0*ny0*nz0, r_type, source, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
	 call MPI_BCAST(datum, nx0*ny0*nz0, r_type, source, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
	 call MPI_BCAST(datum, nx0*ny0*nz0, r_type, source, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BCAST(datum, nx0*ny0*nz0, r_type, source, MPI_COMM_WORLD, ierror)
      end if
      !
      return
   end subroutine broadcast3D
   !
   !
   !
   subroutine int_scatter(info, header, datum, source, nx0, ny0, nz0)
      implicit none
      type(NodeInfo), intent(in) :: info
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0, ny0, nz0
      integer, dimension(nx0,ny0,nz0), intent(inout) :: datum
      integer :: id, idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: myid, myidx, myidy, nproc, nxpe, nype, count, nx, ny, nz, ierror
      integer, dimension(:,:,:), allocatable :: datum_r
      integer, dimension(:,:,:,:), allocatable :: datum_s
      !
      if (trim(header) == 'xy') then
	 call get_nxpe(info, nxpe); call get_nype(info, nype)
         call get_myidx(info, myidx); call get_myidy(info, myidy)
	 myid = myidy*nxpe+myidx
	 nproc = nxpe*nype
      else
         call get_myid(info, myid, nproc)
      end if
      nx = 0
      ny = 0
      nz = nz0
      do id = 1, nproc
	 call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	 nx = max(nx, ie-is+1)
	 ny = max(ny, je-js+1)
      end do
      if (myid == source) allocate(datum_s(nx,ny,nz,nproc))
      allocate(datum_r(nx,ny,nz))
      count = nx*ny*nz
      !
      if (myid == source) then
	 datum_s = 0
	 do id = 1, nproc
	    call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    nx = ie-is+1
	    ny = je-js+1
	    datum_s(1:nx,1:ny,:,id) = datum(is:ie,js:je,:)
	 end do
      end if
      if (trim(header) == 'xy') then
	 call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
	 call MPI_SCATTER(datum_s, count, MPI_INTEGER, datum_r, count, MPI_INTEGER, source, MPI_COMM_XYWORLD, ierror)
	 !call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BARRIER(MPI_COMM_WORLD, ierror)
	 call MPI_SCATTER(datum_s, count, MPI_INTEGER, datum_r, count, MPI_INTEGER, source, MPI_COMM_WORLD, ierror)
	 !call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      end if
      call get_nx(info, nx)
      call get_ny(info, ny)
      datum(1:nx,1:ny,:) = datum_r(1:nx,1:ny,:)
      deallocate(datum_r)
      if (myid == source) deallocate(datum_s)
      !
      return
   end subroutine int_scatter
   !
   !
   !
   subroutine scatter(info, header, datum, source, nx0, ny0, nz0)
      implicit none
      type(NodeInfo), intent(in) :: info
      character(len=*), intent(in) :: header
      integer, intent(in) :: source, nx0, ny0, nz0
      real(kind=r_size), dimension(nx0,ny0,nz0), intent(inout) :: datum
      integer :: id, idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: myid, myidx, myidy, nproc, nxpe, nype, count, nx, ny, nz, ierror
      real(kind=r_size), dimension(:,:,:), allocatable :: datum_r
      real(kind=r_size), dimension(:,:,:,:), allocatable :: datum_s
      !
      if (trim(header) == 'xy') then
	 call get_nxpe(info, nxpe); call get_nype(info, nype)
         call get_myidx(info, myidx); call get_myidy(info, myidy)
	 myid = myidy*nxpe+myidx
	 nproc = nxpe*nype
      else
         call get_myid(info, myid, nproc)
      end if
      nx = 0
      ny = 0
      nz = nz0
      do id = 1, nproc
	 call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	 nx = max(nx, ie-is+1)
	 ny = max(ny, je-js+1)
      end do
      if (myid == source) allocate(datum_s(nx,ny,nz,nproc))
      allocate(datum_r(nx,ny,nz))
      count = nx*ny*nz
      !
      if (myid == source) then
	 datum_s = 0.
	 do id = 1, nproc
	    call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    nx = ie-is+1
	    ny = je-js+1
	    datum_s(1:nx,1:ny,:,id) = datum(is:ie,js:je,:)
	 end do
      end if
      if (trim(header) == 'xy') then
	 call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
	 call MPI_SCATTER(datum_s, count, r_type, datum_r, count, r_type, source, MPI_COMM_XYWORLD, ierror)
	 !call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_BARRIER(MPI_COMM_WORLD, ierror)
	 call MPI_SCATTER(datum_s, count, r_type, datum_r, count, r_type, source, MPI_COMM_WORLD, ierror)
	 !call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      end if
      call get_nx(info, nx)
      call get_ny(info, ny)
      datum(1:nx,1:ny,:) = datum_r(1:nx,1:ny,:)
      deallocate(datum_r)
      if (myid == source) deallocate(datum_s)
      !
      return
   end subroutine scatter
   !
   !
   !
   subroutine int_gather(info, header, datum, destination, nx0, ny0, nz0)
      implicit none
      type(NodeInfo), intent(in) :: info
      character(len=*), intent(in) :: header
      integer, intent(in) :: destination, nx0, ny0, nz0
      integer, dimension(nx0,ny0,nz0), intent(inout) :: datum
      integer :: id, idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: myid, myidx, myidy, nproc, nxpe, nype, count, nx, ny, nz, ierror
      integer, dimension(:,:,:), allocatable :: datum_s
      integer, dimension(:,:,:,:), allocatable :: datum_r
      !
      if (trim(header) == 'xy') then
	 call get_nxpe(info, nxpe); call get_nype(info, nype)
         call get_myidx(info, myidx); call get_myidy(info, myidy)
	 myid = myidy*nxpe+myidx
	 nproc = nxpe*nype
      else
         call get_myid(info, myid, nproc)
      end if
      nx = 0; ny = 0; nz = nz0
      do id = 1, nproc
	 call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	 nx = max(nx, ie-die-is-dis+1)
	 ny = max(ny, je-dje-js-djs+1)
      end do
      if (myid == destination) allocate(datum_r(nx,ny,nz,nproc))
      allocate(datum_s(nx,ny,nz))
      count = nx*ny*nz
      !
      datum_s = 0
      call get_nx(info, nx); call get_ny(info, ny)
      call get_di(info, dis, die); call get_dj(info, djs, dje)
      datum_s(1:nx-dis-die,1:ny-djs-dje,:) = datum(1+dis:nx-die,1+djs:ny-dje,:)
      if (trim(header) == 'xy') then
	 call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
	 call MPI_GATHER(datum_s, count, MPI_INTEGER, datum_r, count, MPI_INTEGER, destination, MPI_COMM_XYWORLD, ierror)
	 !call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
      else
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         call MPI_GATHER(datum_s, count, MPI_INTEGER, datum_r, count, MPI_INTEGER, destination, MPI_COMM_WORLD, ierror)
         !call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      end if
      if (myid == destination) then
	 do id = 1, nproc
	    call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    nx = ie-die-is-dis+1
	    ny = je-dje-js-djs+1
	    datum(is+dis:ie-die,js+djs:je-dje,:) = datum_r(1:nx,1:ny,:,id)
	 end do
	 deallocate(datum_r)
      end if
      deallocate(datum_s)
      !
      return
   end subroutine int_gather
   !
   !
   !
   subroutine gather(info, header, datum, destination, nx0, ny0, nz0)
      implicit none
      type(NodeInfo), intent(in) :: info
      character(len=*), intent(in) :: header
      integer, intent(in) :: destination, nx0, ny0, nz0
      real(kind=r_size), dimension(nx0,ny0,nz0), intent(inout) :: datum
      integer :: id, idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: myid, myidx, myidy, nproc, nxpe, nype, count, nx, ny, nz, ierror
      real(kind=r_size), dimension(:,:,:), allocatable :: datum_s
      real(kind=r_size), dimension(:,:,:,:), allocatable :: datum_r
      !
      if (trim(header) == 'xy') then
	 call get_nxpe(info, nxpe); call get_nype(info, nype)
         call get_myidx(info, myidx); call get_myidy(info, myidy)
	 myid = myidy*nxpe+myidx
	 nproc = nxpe*nype
      else
         call get_myid(info, myid, nproc)
      end if
      nx = 0; ny = 0; nz = nz0
      do id = 1, nproc
	 call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	 nx = max(nx, ie-die-is-dis+1)
	 ny = max(ny, je-dje-js-djs+1)
      end do
      if (myid == destination) allocate(datum_r(nx,ny,nz,nproc))
      allocate(datum_s(nx,ny,nz))
      count = nx*ny*nz
      !
      datum_s = 0.
      call get_nx(info, nx); call get_ny(info, ny)
      call get_di(info, dis, die); call get_dj(info, djs, dje)
      datum_s(1:nx-dis-die,1:ny-djs-dje,:) = datum(1+dis:nx-die,1+djs:ny-dje,:)
      if (trim(header) == 'xy') then
	 call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
	 call MPI_GATHER(datum_s, count, r_type, datum_r, count, r_type, destination, MPI_COMM_XYWORLD, ierror)
	 !call MPI_BARRIER(MPI_COMM_XYWORLD, ierror)
      else
         call MPI_BARRIER(MPI_COMM_WORLD, ierror)
         call MPI_GATHER(datum_s, count, r_type, datum_r, count, r_type, destination, MPI_COMM_WORLD, ierror)
         !call MPI_BARRIER(MPI_COMM_WORLD, ierror)
      end if
      if (myid == destination) then
	 do id = 1, nproc
	    call get_info(info, id-1, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    nx = ie-die-is-dis+1
	    ny = je-dje-js-djs+1
	    datum(is+dis:ie-die,js+djs:je-dje,:) = datum_r(1:nx,1:ny,:,id)
	 end do
	 deallocate(datum_r)
      end if
      deallocate(datum_s)
      !
      return
   end subroutine gather
   !
   !
   !
   subroutine int_allreduce0Dmax(header, datum)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(inout) :: datum
      integer :: ierror, maxvalue
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, maxvalue, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, maxvalue, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, maxvalue, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, maxvalue, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, maxvalue, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
      end if
      datum = maxvalue
      !
      return
   end subroutine int_allreduce0Dmax
   !
   !
   !
   subroutine allreduce0Dminloc(header, datum)
      implicit none
      character(len=*), intent(in) :: header
      real(kind=r_size), dimension(2), intent(inout) :: datum
      integer :: ierror
      real(kind=r_size), dimension(2) :: minvalue
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, minvalue, 1, r_2type, MPI_MINLOC, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, minvalue, 1, r_2type, MPI_MINLOC, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, minvalue, 1, r_2type, MPI_MINLOC, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, minvalue, 1, r_2type, MPI_MINLOC, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, minvalue, 1, r_2type, MPI_MINLOC, MPI_COMM_WORLD, ierror)
      end if
      datum(:) = minvalue(:)
      !
      return
   end subroutine allreduce0Dminloc
   !
   !
   !
   subroutine int_allreduce0D(header, datum)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(inout) :: datum
      integer :: ierror, sum
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, sum, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
      end if
      datum = sum
      !
      return
   end subroutine int_allreduce0D
   !
   !
   !
   subroutine allreduce0D(header, datum)
      implicit none
      character(len=*), intent(in) :: header
      real(kind=r_size), intent(inout) :: datum
      integer :: ierror
      real(kind=r_size) :: sum
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, sum, 1, r_type, MPI_SUM, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, sum, 1, r_type, MPI_SUM, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, sum, 1, r_type, MPI_SUM, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, sum, 1, r_type, MPI_SUM, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, sum, 1, r_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      end if
      datum = sum
      !
      return
   end subroutine allreduce0D
   !
   !
   !
   subroutine allreduce1D(header, datum, nx)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: nx
      real(kind=r_size), dimension(nx), intent(inout) :: datum
      integer :: ierror
      real(kind=r_size), dimension(nx) :: sum
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, sum, nx, r_type, MPI_SUM, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, sum, nx, r_type, MPI_SUM, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, sum, nx, r_type, MPI_SUM, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, sum, nx, r_type, MPI_SUM, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, sum, nx, r_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      end if
      datum(:) = sum(:)
      !
      return
   end subroutine allreduce1D
   !
   !
   !
   subroutine int_allreduce2D(header, datum, nx, ny)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: nx, ny
      integer, dimension(nx,ny), intent(inout) :: datum
      integer :: ierror
      integer, dimension(nx,ny) :: sum
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, MPI_INTEGER, MPI_SUM, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, MPI_INTEGER, MPI_SUM, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, MPI_INTEGER, MPI_SUM, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, MPI_INTEGER, MPI_SUM, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, sum, nx*ny, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
      end if
      datum(:,:) = sum(:,:)
      !
      return
   end subroutine int_allreduce2D
   !
   !
   !
   subroutine allreduce2D(header, datum, nx, ny)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: nx, ny
      real(kind=r_size), dimension(nx,ny), intent(inout) :: datum
      integer :: ierror
      real(kind=r_size), dimension(nx,ny) :: sum
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, r_type, MPI_SUM, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, r_type, MPI_SUM, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, r_type, MPI_SUM, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, sum, nx*ny, r_type, MPI_SUM, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, sum, nx*ny, r_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      end if
      datum(:,:) = sum(:,:)
      !
      return
   end subroutine allreduce2D
   !
   !
   !
   subroutine allreduce3D(header, datum, nx, ny, nz)
      implicit none
      character(len=*), intent(in) :: header
      integer, intent(in) :: nx, ny, nz
      real(kind=r_size), dimension(nx,ny,nz), intent(inout) :: datum
      integer :: ierror
      real(kind=r_size), dimension(nx,ny,nz) :: sum
      !
      if (trim(header) == 'x') then
         call MPI_ALLREDUCE(datum, sum, nx*ny*nz, r_type, MPI_SUM, MPI_COMM_XWORLD, ierror)
      else if (trim(header) == 'y') then
         call MPI_ALLREDUCE(datum, sum, nx*ny*nz, r_type, MPI_SUM, MPI_COMM_YWORLD, ierror)
      else if (trim(header) == 'e') then
         call MPI_ALLREDUCE(datum, sum, nx*ny*nz, r_type, MPI_SUM, MPI_COMM_EWORLD, ierror)
      else if (trim(header) == 'xy') then
         call MPI_ALLREDUCE(datum, sum, nx*ny*nz, r_type, MPI_SUM, MPI_COMM_XYWORLD, ierror)
      else
	 call MPI_ALLREDUCE(datum, sum, nx*ny*nz, r_type, MPI_SUM, MPI_COMM_WORLD, ierror)
      end if
      datum(:,:,:) = sum(:,:,:)
      !
      return
   end subroutine allreduce3D
   !
   !
   !
   subroutine add_halo(info, datum, nx, ny, nz)
      implicit none
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nz
      real(kind=r_size), dimension(nx,ny,nz), intent(inout) :: datum
      integer :: myid, nproc, myidx, myidy, nxpe, nype, halo, dis, die, djs, dje
      integer :: is, ie, js, je, n, count, ierror
      integer, dimension(6) :: irequest
      integer, dimension(MPI_STATUS_SIZE,6) :: istatus
      real(kind=r_size), dimension(:), allocatable :: left, right, bottom, top, bottomleft, topright
      !
      call get_myid(info, myid, nproc)
      call get_myidx(info, myidx)
      call get_myidy(info, myidy)
      call get_nxpe(info, nxpe)
      call get_nype(info, nype)
      call get_halo(info, halo)
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      ! Send
      n = 0
      if (myidx < nxpe-1) then
	 n = n + 1
	 count = halo*(ny-djs-dje)*nz
         allocate(right(count))
	 right = reshape(datum(nx-halo+1:nx,js:je,:),(/count/))
         call MPI_ISEND(right, count, r_type, myid+1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidy < nype-1) then
	 n = n + 1
	 count = (nx-dis-die)*halo*nz
         allocate(top(count))
	 top = reshape(datum(is:ie,ny-halo+1:ny,:),(/count/))
         call MPI_ISEND(top, count, r_type, myid+nxpe, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidx < nxpe-1 .and. myidy < nype-1) then
	 n = n + 1
	 count = halo*halo*nz
         allocate(topright(count))
	 topright = reshape(datum(nx-halo+1:nx,ny-halo+1:ny,:),(/count/))
         call MPI_ISEND(topright, count, r_type, myid+nxpe+1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      !
      ! Receive
      if (myidx > 0) then
	 n = n + 1
	 count = halo*(ny-djs-dje)*nz
         allocate(left(count))
	 call MPI_IRECV(left, count, r_type, myid-1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidy > 0) then
	 n = n + 1
	 count = (nx-dis-die)*halo*nz
         allocate(bottom(count))
	 call MPI_IRECV(bottom, count, r_type, myid-nxpe, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidx > 0 .and. myidy > 0) then
	 n = n + 1
	 count = halo*halo*nz
         allocate(bottomleft(count))
	 call MPI_IRECV(bottomleft, count, r_type, myid-nxpe-1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      !
      call MPI_WAITALL(n, irequest, istatus, ierror)
      if (myidx > 0) then
	 datum(1+dis:halo+dis,js:je,:) = datum(1+dis:halo+dis,js:je,:) + reshape(left,(/halo,ny-djs-dje,nz/))
      end if
      if (myidy > 0) then
	 datum(is:ie,1+djs:halo+djs,:) = datum(is:ie,1+djs:halo+djs,:) + reshape(bottom,(/nx-dis-die,halo,nz/))
      end if
      if (myidx > 0 .and. myidy > 0) then
	 datum(1+dis:halo+dis,1+djs:halo+djs,:) = datum(1+dis:halo+dis,1+djs:halo+djs,:) + reshape(bottomleft,(/halo,halo,nz/))
      end if
      if (allocated(left)) deallocate(left)
      if (allocated(right)) deallocate(right)
      if (allocated(bottom)) deallocate(bottom)
      if (allocated(top)) deallocate(top)
      if (allocated(bottomleft)) deallocate(bottomleft)
      if (allocated(topright)) deallocate(topright)
      !
      return
   end subroutine add_halo
   !
   !
   !
   subroutine update_halo(info, datum, nx, ny, nz)
      implicit none
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nz
      real(kind=r_size), dimension(nx,ny,nz), intent(inout) :: datum
      integer :: myid, nproc, myidx, myidy, nxpe, nype, halo, dis, die, djs, dje
      integer :: is, ie, js, je, n, count, ierror
      integer, dimension(6) :: irequest
      integer, dimension(MPI_STATUS_SIZE,6) :: istatus
      real(kind=r_size), dimension(:), allocatable :: left, right, bottom, top, bottomleft, topright
      !
      call get_myid(info, myid, nproc)
      call get_myidx(info, myidx)
      call get_myidy(info, myidy)
      call get_nxpe(info, nxpe)
      call get_nype(info, nype)
      call get_halo(info, halo)
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      ! Send
      n = 0
      if (myidx > 0) then
	 n = n + 1
	 count = halo*(ny-djs-dje)*nz
         allocate(left(count))
	 left = reshape(datum(1+dis:halo+dis,js:je,:),(/count/))
	 call MPI_ISEND(left, count, r_type, myid-1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidy > 0) then
	 n = n + 1
	 count = (nx-dis-die)*halo*nz
         allocate(bottom(count))
	 bottom = reshape(datum(is:ie,1+djs:halo+djs,:),(/count/))
	 call MPI_ISEND(bottom, count, r_type, myid-nxpe, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidx > 0 .and. myidy > 0) then
	 n = n + 1
	 count = halo*halo*nz
         allocate(bottomleft(count))
	 bottomleft = reshape(datum(1+dis:halo+dis,1+djs:halo+djs,:),(/count/))
	 call MPI_ISEND(bottomleft, count, r_type, myid-nxpe-1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      !
      ! Receive
      if (myidx < nxpe-1) then
	 n = n + 1
	 count = halo*(ny-djs-dje)*nz
         allocate(right(count))
         call MPI_IRECV(right, count, r_type, myid+1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidy < nype-1) then
	 n = n + 1
	 count = (nx-dis-die)*halo*nz
         allocate(top(count))
         call MPI_IRECV(top, count, r_type, myid+nxpe, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      if (myidx < nxpe-1 .and. myidy < nype-1) then
	 n = n + 1
	 count = halo*halo*nz
         allocate(topright(count))
         call MPI_IRECV(topright, count, r_type, myid+nxpe+1, 0, MPI_COMM_WORLD, irequest(n), ierror)
      end if
      !
      call MPI_WAITALL(n, irequest, istatus, ierror)
      if (myidx < nxpe-1) then
	 datum(nx-halo+1:nx,js:je,:) = reshape(right,(/halo,ny-djs-dje,nz/))
      end if
      if (myidy < nype-1) then
	 datum(is:ie,ny-halo+1:ny,:) = reshape(top,(/nx-dis-die,halo,nz/))
      end if
      if (myidx < nxpe-1 .and. myidy < nype-1) then
	 datum(nx-halo+1:nx,ny-halo+1:ny,:) = reshape(topright,(/halo,halo,nz/))
      end if
      if (allocated(left)) deallocate(left)
      if (allocated(right)) deallocate(right)
      if (allocated(bottom)) deallocate(bottom)
      if (allocated(top)) deallocate(top)
      if (allocated(bottomleft)) deallocate(bottomleft)
      if (allocated(topright)) deallocate(topright)
      !
      return
   end subroutine update_halo
   !
   !
   !
end module NodeMPI
