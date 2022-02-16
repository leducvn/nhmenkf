program MVSensitivity1000
   use variable
   use nhmlib
   use NodeInfo_class
   use NodeControl_class
   use NodeScorr_class
   use NodeMPI
   implicit none
   !
   character(2) :: header
   integer :: nproc, myid, np, nx, ny, nz, nt, ne, myide
   integer :: ierror, ip, jp, is, ie
   integer, dimension(:), allocatable :: i0, j0
   real(r_size), dimension(:), allocatable :: ymean, ystd
   real(r_sngl), dimension(:,:), allocatable :: ytmp
   real(r_size), dimension(:,:), allocatable :: y0, y
   !
   type(NodeInfo) :: info
   type(NodeControl), dimension(:), allocatable :: xnhm, xpert
   type(NodeControl) :: xmean, xstd, xtmp, xy, xymean, xyrmse
   !
   !
   !
   ! MPI Initialization
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, nproc, ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD, myid,  ierror)
   !
   ! Read namelist
   call initialize_namelist()
   open(10, file='namelist.txt')
   read(10, control_nl)
   read(10, model_nl)
   close(10)
   !
   ! MPI Decomposition
   call new(info, 3, nxpe, nype, nepe, myid, nx0, ny0, ne0, 0)
   call initialize_mpi(info)
   call get_myide(info, myide)
   call get_nx(info, nx)
   call get_ny(info, ny)
   call get_ne(info, ne)
   nz = nz0
   nt = nt0
   !
   !
   !
   ! Y
   if (myid == 0) then
      open(90, file='y', form='binary')
      read(90) np
      allocate(i0(np), j0(np), ymean(np), ystd(np), ytmp(np,ne0), y0(np,ne0))
      read(90) i0
      read(90) j0
      do ip = 1, np
         read(90) ytmp(ip,:)
      end do
      close(90)
      y0(:,:) = ytmp(:,:)
      print*, np, i0(1), j0(1), y0(1,:)
      do ip = 1, np
         ymean(ip) = sum(y0(ip,:))/ne0
	 y0(ip,:) = y0(ip,:) - ymean(ip)
	 ystd(ip) = sqrt(sum(y0(ip,:)**2)/(ne0-1))
	 if (ystd(ip) < 1.e-12) then
	    y0(ip,:) = 0.d0
	 else
	    y0(ip,:) = y0(ip,:)/ystd(ip)
	 end if
	 y0(ip,:) = y0(ip,:)/sqrt(1.d0*(ne0-1))
      end do
      deallocate(ytmp)
   end if
   call int_broadcast0D('all', np, 0)
   if (myid > 0) allocate(i0(np), j0(np), ymean(np), ystd(np), y0(np,ne0))
   allocate(y(np,ne))
   call int_broadcast1D('all', i0, 0, np)
   call int_broadcast1D('all', j0, 0, np)
   call broadcast1D('all', ymean, 0, np)
   call broadcast1D('all', ystd, 0, np)
   call broadcast2D('all', y0, 0, np, ne0)
   is = myide*ne + 1; ie = is + ne - 1
   y(:,1:ne) = y0(:,is:ie)
   deallocate(y0)
   !
   !
   !
   ! X
   allocate(xnhm(ne))
   do ie = 1, ne
      call new_nhm_NodeControl(xnhm(ie), info, .True., nz0, nt0)
   end do
   call read_ensemble(xnhm, info, 'x', ne)
   call new_sens1_NodeControl(xmean, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xstd, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xtmp, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xy, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xymean, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xyrmse, info, nz0, nt0, xnhm(1))
   allocate(xpert(ne))
   do ie = 1, ne
      call new_sens1_NodeControl(xpert(ie), info, nz0, nt0, xnhm(ie))
      call destroy(xnhm(ie))
   end do
   deallocate(xnhm)
   !
   ! Forecast perturbations:
   call set_const(xmean, 0.d0)
   do ie = 1, ne
      call add(xmean, xpert(ie))
   end do
   call allreduce_ens(xmean)
   call divide(xmean, 1.d0*ne0)
   do ie = 1, ne
      call subtract(xpert(ie), xmean)
   end do
   call set_const(xstd, 0.d0)
   do ie = 1, ne
      xtmp = xpert(ie)
      call power(xtmp, 2.d0)
      call add(xstd, xtmp)
   end do
   call allreduce_ens(xstd)
   call divide(xstd, 1.d0*(ne0-1))
   call power(xstd, 0.5d0)
   do ie = 1, ne
      call divide(xpert(ie), xstd)
      call divide(xpert(ie), sqrt(1.d0*(ne0-1)))
   end do
   !
   !
   !
   !
   ! Changes of y
   call set_const(xymean, 0.d0)
   call set_const(xyrmse, 0.d0)
   do ip = 1, np
      call set_const(xy, 0.d0)
      do ie = 1, ne
         xtmp = xpert(ie)
         call multiply(xtmp, y(ip,ie))
         call add(xy, xtmp)
      end do
      call allreduce_ens(xy)
      call multiply(xy, ystd(ip))
      call add(xymean, xy)
      call power(xy, 2.d0)
      call add(xyrmse, xy)
   end do
   call divide(xymean, 1.d0*np)
   call divide(xyrmse, 1.d0*np)
   call power(xyrmse, 0.5d0)
   !
   !
   !
   ! Output
   call write_control(xmean, 1, info, 'xmean')
   call write_control(xstd, 1, info, 'xstd')
   call write_control(xymean, 1, info, 'xymean')
   call write_control(xyrmse, 1, info, 'xyrmse')
   !
   ! deallocation
   call destroy(xmean); call destroy(xstd); call destroy(xtmp)
   call destroy(xy); call destroy(xymean); call destroy(xyrmse)
   do ie = 1, ne
      call destroy(xpert(ie))
   end do
   deallocate(xpert)
   deallocate(i0, j0, ymean, ystd, y)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program MVSensitivity1000
