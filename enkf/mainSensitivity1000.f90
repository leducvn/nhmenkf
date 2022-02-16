program sensitivity1000
   use variable
   use nhmlib
   use NodeInfo_class
   use NodeControl_class
   use NodeMPI
   implicit none
   !
   integer :: nproc, myid, nx, ny, nz, nt, ne, nlev, myide
   integer :: i0, j0, ierror, is, ie 
   real(r_sngl) :: Rtmp
   real(r_size) :: R0, ymean, ystd, norm, ustd
   real(r_sngl), dimension(:), allocatable :: ytmp
   real(r_size), dimension(:), allocatable :: y0, y
   !
   type(NodeInfo) :: info
   type(NodeControl), dimension(:), allocatable :: xnhm, xpert
   type(NodeControl) :: x, u, xmean, xstd
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
   allocate(y0(ne0), y(ne))
   if (myid == 0) then
      allocate(ytmp(ne0))
      open(90, file='y', form='binary')
      read(90) i0, j0, Rtmp
      read(90) ytmp
      close(90)
      print*, i0, j0, Rtmp, ytmp
      R0 = 1000.d0*Rtmp/resolution
      y0(:) = ytmp(:)
      ymean = sum(y0)/ne0
      y0 = y0 - ymean
      ystd = sqrt(sum(y0**2)/(ne0-1))
      y0 = y0/ystd
      y0 = y0/sqrt(1.d0*(ne0-1))
      deallocate(ytmp)
   end if
   call int_broadcast0D('all', i0, 0)
   call int_broadcast0D('all', j0, 0)
   call broadcast0D('all', R0, 0)
   call broadcast1D('all', y0, 0, ne0)
   is = myide*ne + 1; ie = is + ne - 1
   y(1:ne) = y0(is:ie)
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
   call new_sens1_NodeControl(x, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(u, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xmean, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xstd, info, nz0, nt0, xnhm(1))
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
      x = xpert(ie)
      call power(x, 2.d0)
      call add(xstd, x)
   end do
   call allreduce_ens(xstd)
   call divide(xstd, 1.d0*(ne0-1))
   call power(xstd, 0.5d0)
   do ie = 1, ne
      call divide(xpert(ie), xstd)
      call divide(xpert(ie), sqrt(1.d0*(ne0-1)))
   end do
   call write_control(xmean, 1, info, 'xmean')
   call write_control(xstd, 1, info, 'xstd')
   !
   !
   !
   ! Sensitivities
   call set_const(u, 0.d0)
   do ie = 1, ne
      x = xpert(ie)
      call multiply(x, y(ie))
      call add(u, x)
   end do
   call allreduce_ens(u)
   call write_control(u, 1, info, 's')
   call innerproduct(u, u, info, norm)
   call allreduce0D('xy', norm)
   norm = sqrt(norm)
   call divide(u, norm)
   if (myid == 0) print*, norm
   !
   ! Project perturbations on u
   ustd = 0.d0
   do ie = 1, ne
      call innerproduct(xpert(ie), u, info, norm)
      call allreduce0D('xy', norm)
      ustd = ustd + norm**2
   end do
   call allreduce0D('e', ustd)
   ustd = sqrt(ustd/(ne0-1))
   if (myid == 0) then
      open(90, file="ustd.txt")
      write(90,'(e14.7)') ustd
      close(90)
   end if
   !
   ! deallocation
   call destroy(x); call destroy(u); call destroy(xmean); call destroy(xstd)
   do ie = 1, ne
      call destroy(xpert(ie))
   end do
   deallocate(xpert)
   deallocate(y)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program sensitivity1000
