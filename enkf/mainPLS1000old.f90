program PLS1000
   use variable
   use enkflib
   use nhmlib
   use NodeInfo_class
   use NodeControl_class
   use NodeScorr_class
   use NodeMPI
   implicit none
   !
   character(2) :: header
   integer :: nproc, myid, np, nx, ny, nz, nt, ne, nrank, myide
   integer :: ierror, ip, jp, is, ie, irank
   integer, dimension(:), allocatable :: i0, j0
   real(r_size), parameter :: tolerance = 0.99
   real(r_size) :: total, explained
   real(r_size), dimension(:), allocatable :: ymean, ystd, sv, umean, vmean
   real(r_sngl), dimension(:,:), allocatable :: ytmp
   real(r_size), dimension(:,:), allocatable :: y0, y, YYT, v
   !
   type(NodeInfo) :: info
   type(NodeControl), dimension(:), allocatable :: xnhm, xpert
   type(NodeControl) :: xmean, xstd, u, xy1, xy2
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
   hscale = 1000.d0*hscale/resolution
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
   call new_sens1_NodeControl(u, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xy1, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(xy2, info, nz0, nt0, xnhm(1))
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
      u = xpert(ie)
      call power(u, 2.d0)
      call add(xstd, u)
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
   ! YYT
   allocate(sv(np), v(np,np), YYT(np,np))
   do ip = 1, np
      call set_const(xy1, 0.d0)
      do ie = 1, ne
         u = xpert(ie)
         call multiply(u, y(ip,ie))
         call add(xy1, u)
      end do
      call allreduce_ens(xy1)
      !call localize(xy1, info, i0(ip), j0(ip))
      call innerproduct(xy1, xy1, info, YYT(ip,ip))
      call allreduce0D('xy', YYT(ip,ip))
      !if (myid == 0) print*, ip, ip, YYT(ip,ip)
      !
      do jp = ip+1, np
         call set_const(xy2, 0.d0)
         do ie = 1, ne
            u = xpert(ie)
            call multiply(u, y(jp,ie))
            call add(xy2, u)
         end do
         call allreduce_ens(xy2)
         !call localize(xy2, info, i0(jp), j0(jp))
	 call innerproduct(xy1, xy2, info, YYT(ip,jp))
	 call allreduce0D('xy', YYT(ip,jp))
         YYT(jp,ip) = YYT(ip,jp)
	 !if (myid == 0) print*, ip, jp, YYT(ip,jp)
      end do
   end do
   !
   ! Eigen-decomposition for sv and V
   call mtx_eigen(1, np, YYT, sv, v, nrank)
   where (sv < 0.) sv = 0.
   sv = sqrt(sv)
   total = sum(sv)
   explained = 0.d0
   do ip = 1, np
      explained = explained + sv(ip)
      if (explained/total > tolerance) then
         nrank = ip
	 exit
      end if
   end do
   if (myid == 0) print*, nrank, sv(1:10)
   nrank = min(10,nrank)
   allocate(umean(nrank), vmean(nrank))
   !
   !
   !
   ! Find U
   do irank = 1, nrank
      call set_const(u, 0.d0)
      do ip = 1, np
         call set_const(xy1, 0.d0)
         do ie = 1, ne
            xy2 = xpert(ie)
            call multiply(xy2, y(ip,ie))
            call add(xy1, xy2)
         end do
         call allreduce_ens(xy1)
         !call localize(xy1, info, i0(ip), j0(ip))
	 call multiply(xy1, v(ip,irank))
	 call add(u, xy1)
      end do
      call divide(u, sv(irank))
      !
      umean(irank) = 0.d0; vmean(irank) = 0.d0
      do ie = 1, ne
         call innerproduct(xpert(ie), u, info, total)
         call allreduce0D('xy', total)
	 umean(irank) = umean(irank) + abs(total)
	 vmean(irank) = vmean(irank) + abs(sum(y(:,ie)*v(:,irank))) 
      end do
      call allreduce0D('e', umean(irank)); call allreduce0D('e', vmean(irank))
      umean(irank) = umean(irank)/ne0; vmean(irank) = vmean(irank)/ne0
      call multiply(u, umean(irank)*sqrt(1.d0*(ne0-1)))
      call multiply(u, xstd)
      call add(u, xmean)
      write(header(1:2),'(I2.2)') irank
      call write_control(u, 1, info, 'u'//header)
      v(:,irank) = vmean(irank)*sqrt(1.d0*(ne0-1))*v(:,irank)
      v(:,irank) = v(:,irank)*ystd + ymean
      if (myid == 0) print*, irank, sv(irank), umean(irank), vmean(irank) 
   end do
   !
   !
   !
   ! Output
   !call write_control(xmean, info, 'xmean')
   !call write_control(xstd, info, 'xstd')
   if (myid == 0) then
      open(90, file='v', form='binary')
      write(90) np,nrank
      !write(90) ymean(1:np), ystd(1:np)
      write(90) sv(1:nrank), umean, vmean
      write(90) v(1:np,1:nrank)
      close(90)
   end if
   !
   ! deallocation
   call destroy(xmean); call destroy(xstd); call destroy(u)
   call destroy(xy1); call destroy(xy2)
   do ie = 1, ne
      call destroy(xpert(ie))
   end do
   deallocate(xpert)
   deallocate(i0, j0, ymean, ystd, y)
   deallocate(YYT, sv, v, umean, vmean)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program PLS1000
