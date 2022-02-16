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
   real(r_size) :: norm, total, explained
   real(r_size), dimension(:), allocatable :: ymean, ystd, sv, ustd, vstd, coef, w
   real(r_sngl), dimension(:,:), allocatable :: ytmp
   real(r_size), dimension(:,:), allocatable :: y0, y, v, Vlcz, Tlcz
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
   ! Lanczos
   allocate(w(np), sv(niteration), v(np,niteration), Vlcz(np,niteration), Tlcz(niteration,niteration))
   Tlcz(:,:) = 0.d0
   if (myid == 0) then
      call random_number(Vlcz(:,1))
      norm = sqrt(sum(Vlcz(:,1)**2))
      Vlcz(:,1) = Vlcz(:,1)/norm
   endif
   call broadcast1D('all', Vlcz(:,1), 0, np)
   ! Av
   call set_const(xy2, 0.d0)
   do ip = 1, np
      call set_const(xy1, 0.d0)
      do ie = 1, ne
         u = xpert(ie)
         call multiply(u, y(ip,ie))
         call add(xy1, u)
      end do
      call allreduce_ens(xy1)
      call multiply(xy1, Vlcz(ip,1))
      call add(xy2, xy1)
   end do
   do ip = 1, np
      call set_const(xy1, 0.d0)
      do ie = 1, ne
         u = xpert(ie)
         call multiply(u, y(ip,ie))
         call add(xy1, u)
      end do
      call allreduce_ens(xy1)
      call innerproduct(xy1, xy2, info, w(ip))
      call allreduce0D('xy', w(ip))
      !if (myid == 0) print*, 'iter = 1', ip, w(ip)
   end do
   Tlcz(1,1) = sum(Vlcz(:,1)*w)
   w = w - Tlcz(1,1)*Vlcz(:,1)
   if (myid == 0) print*, 'iter = 1', Tlcz(1,1)
   !
   do jp = 2, niteration
      Tlcz(jp,jp-1) = sqrt(sum(w**2))
      if (Tlcz(jp,jp-1) >= 1.e-12) then
         Vlcz(:,jp) = w/Tlcz(jp,jp-1)
      else
         total = Tlcz(jp,jp-1)
	 do while (total < 1.e-12)
	    if (myid == 0) then
               call random_number(Vlcz(:,jp))
               norm = sqrt(sum(Vlcz(:,jp)**2))
               Vlcz(:,jp) = Vlcz(:,jp)/norm
            endif
	    call broadcast1D('all', Vlcz(:,jp), 0, np)
	    do ip = 1, jp-1
	       norm = sum(Vlcz(:,jp)*Vlcz(:,ip))
	       Vlcz(:,jp) = Vlcz(:,jp) - norm*Vlcz(:,ip)
	    end do
	    total = sqrt(sum(Vlcz(:,jp)**2))
	 end do
      end if
      !
      ! Av
      call set_const(xy2, 0.d0)
      do ip = 1, np
         call set_const(xy1, 0.d0)
         do ie = 1, ne
            u = xpert(ie)
            call multiply(u, y(ip,ie))
            call add(xy1, u)
         end do
         call allreduce_ens(xy1)
	 call multiply(xy1, Vlcz(ip,jp))
	 call add(xy2, xy1)
      end do
      do ip = 1, np
         call set_const(xy1, 0.d0)
         do ie = 1, ne
            u = xpert(ie)
            call multiply(u, y(ip,ie))
            call add(xy1, u)
         end do
         call allreduce_ens(xy1)
         call innerproduct(xy1, xy2, info, w(ip))
         call allreduce0D('xy', w(ip))
	 !if (myid == 0) print*, 'iter = ', jp, ip, w(ip)
      end do
      Tlcz(jp,jp) = sum(Vlcz(:,jp)*w)
      w = w - Tlcz(jp,jp)*Vlcz(:,jp) - Tlcz(jp,jp-1)*Vlcz(:,jp-1)
      !
      ! Reorthogonalization
      do ip = 1, jp
	 norm = sum(w*Vlcz(:,ip))
	 w = w - norm*Vlcz(:,ip)
      end do
      Tlcz(jp-1,jp) = Tlcz(jp,jp-1)
      if (myid == 0) print*, 'iter = ', jp, Tlcz(jp,jp)
   end do
   !
   !
   !
   ! Eigen-decomposition for sv and V
   call mtx_eigen(1, niteration, Tlcz, sv, v(1:niteration,:), nrank)
   where (sv < 0.) sv = 0.
   sv = sqrt(sv)
   total = sum(sv)
   explained = 0.d0
   do ip = 1, niteration
      explained = explained + sv(ip)
      if (explained/total > tolerance) then
         nrank = ip
	 exit
      end if
   end do
   nrank = min(10,nrank)
   if (myid == 0) print*, nrank, sv(1:nrank)
   allocate(ustd(nrank), vstd(nrank), coef(nrank))
   !
   !
   !
   ! Find U
   do irank = 1, nrank
      write(header(1:2),'(I2.2)') irank
      w(:) = 0.d0
      do ip = 1, niteration
         w(:) = w(:) + v(ip,irank)*Vlcz(:,ip)
      end do
      !
      call set_const(u, 0.d0)
      do ip = 1, np
         call set_const(xy1, 0.d0)
         do ie = 1, ne
            xy2 = xpert(ie)
            call multiply(xy2, y(ip,ie))
            call add(xy1, xy2)
         end do
         call allreduce_ens(xy1)
	 call multiply(xy1, w(ip))
	 call add(u, xy1)
      end do
      call divide(u, sv(irank))
      !
      ustd(irank) = 0.d0; vstd(irank) = 0.d0; coef(irank) = 0.d0
      do ie = 1, ne
         call innerproduct(xpert(ie), u, info, norm)
         call allreduce0D('xy', norm)
	 total = sum(y(:,ie)*w(:))
	 ustd(irank) = ustd(irank) + norm**2
	 vstd(irank) = vstd(irank) + total**2
	 coef(irank) = coef(irank) + norm*total
      end do
      call allreduce0D('e', ustd(irank)); call allreduce0D('e', vstd(irank)); call allreduce0D('e', coef(irank))
      ustd(irank) = sqrt(ustd(irank)/(ne0-1)); vstd(irank) = sqrt(vstd(irank)/(ne0-1))
      coef(irank) = coef(irank)/(ne0-1); coef(irank) = coef(irank)/ustd(irank)**2
      call write_control(u, 1, info, 'u'//header)
      xy1 = u
      call multiply(xy1, ustd(irank)*sqrt(1.d0*(ne0-1)))
      call multiply(xy1, xstd)
      !call add(xy1, xmean)
      call write_control(xy1, 1, info, 'x'//header)
      v(:,irank) = w(:)
      !v(:,irank) = vstd(irank)*sqrt(1.d0*(ne0-1))*w(:)
      !v(:,irank) = v(:,irank)*ystd + ymean
      if (myid == 0) print*, irank, sv(irank), ustd(irank), vstd(irank), coef(irank) 
   end do
   !
   !
   !
   ! Output
   call write_control(xmean, 1, info, 'xmean')
   call write_control(xstd, 1, info, 'xstd')
   if (myid == 0) then
      open(90, file='v', form='binary')
      write(90) np,nrank
      write(90) ymean, ystd
      write(90) sv(1:nrank), ustd(1:nrank), vstd(1:nrank), coef(1:nrank)
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
   deallocate(w, sv, v, Vlcz, Tlcz, ustd, vstd, coef)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program PLS1000
