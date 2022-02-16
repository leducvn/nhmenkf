program sensitivity1000
   use variable
   use nhmlib
   use NodeInfo_class
   use NodeControl_class
   use NodeScorr_class
   use NodeMPI
   implicit none
   !
   integer :: nproc, myid, nx, ny, nz, nt, ne, nlev0, nlev, myide
   integer :: i0, j0, ierror, is, ie, iflag, iter, jter, kter 
   real(r_sngl) :: Rtmp
   real(r_sngl), dimension(:), allocatable :: ytmp
   real(r_size), dimension(:), allocatable :: y0, y
   !
   type(NodeScorr) :: Slocx, Slocy
   !
   type(NodeInfo) :: info
   type(NodeControl), dimension(:), allocatable :: xnhm, xpert, w
   type(NodeControl) :: x, r, b
   type(NodeControl), dimension(:), allocatable :: q
   !
   real(r_size), parameter :: tolerance = 1.e-6, lambda = 0.001
   real(r_size) :: R0, rnorm0, ratio, norm, rmat1, rmat2
   real(r_size), dimension(:), allocatable :: rnorm, p, s
   real(r_size), dimension(:,:), allocatable :: hmat, rmat, g
   real(r_size), dimension(:,:,:,:), allocatable :: ftmp
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
      rmat1 = sum(y0)/ne0
      y0 = y0 - rmat1
      rmat1 = sqrt(sum(y0**2)/(ne0-1))
      y0 = y0/rmat1
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
   call new_sens1_NodeControl(r, info, nz0, nt0, xnhm(1))
   call new_sens1_NodeControl(b, info, nz0, nt0, xnhm(1))
   allocate(q(niteration))
   do iter = 1, niteration
      call new_sens1_NodeControl(q(iter), info, nz0, nt0, xnhm(1))
   end do
   allocate(xpert(ne))
   do ie = 1, ne
      call new_sens1_NodeControl(xpert(ie), info, nz0, nt0, xnhm(ie))
      call destroy(xnhm(ie))
   end do
   deallocate(xnhm)
   !
   ! Forecast perturbations: use x and b temporarily
   call set_const(x, 0.d0)
   do ie = 1, ne
      call add(x, xpert(ie))
   end do
   call allreduce_ens(x)
   call divide(x, 1.d0*ne0)
   do ie = 1, ne
      call subtract(xpert(ie), x)
   end do
   call set_const(x, 0.d0)
   do ie = 1, ne
      b = xpert(ie)
      call power(b, 2.d0)
      call add(x, b)
   end do
   call allreduce_ens(x)
   call divide(x, 1.d0*(ne0-1))
   call power(x, 0.5d0)
   do ie = 1, ne
      call divide(xpert(ie), x)
      call divide(xpert(ie), sqrt(1.d0*(ne0-1)))
   end do
   !
   ! Control
   allocate(w(ne))
   do ie = 1, ne
      call new_enkf0_NodeControl(w(ie), info, .True., nz0, 1, 1)
   end do
   allocate(rnorm(0:niteration), p(niteration+1), s(niteration))
   allocate(hmat(niteration+1,niteration), rmat(niteration,niteration), g(niteration,2))
   !
   !
   !
   ! Localization
   call get_nlev(w(1), nlev)
   allocate(ftmp(nx,ny,nlev,1))
   call new(Slocx, nx0, nlev, 'x')
   call set_Scorr(Slocx, info, 'loc')
   call new(Slocy, ny0, nlev, 'y')
   call set_Scorr(Slocy, info, 'loc')
   !
   !
   !
   ! Preprocessing: calculate b
   call set_const(b, 0.d0)
   do ie = 1, ne
      x = xpert(ie)
      call multiply(x, y(ie))
      call add(b, x)
   end do
   call allreduce_ens(b)
   !call localize(b, info, i0, j0, R0)
   call write_control(b, 1, info, 'us')
   !
   ! GMRES Ax = b
   !x = b ! initial guess
   !call apply_UT1000(w, ne, x, xpert)
   !do ie = 1, ne
      !call get_field(w(ie), ftmp, nlev)
      !call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
      !call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
      !call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
      !call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
      !call set_field(w(ie), ftmp, nlev)
   !end do
   !call apply_U1000(r, ne, w, xpert)
   !call allreduce_ens(r)
   !call multiply(b, lambda)
   !call add(r, b)
   !call multiply(r, -1.d0)
   !call add(r, x)
   call set_const(x, 0.d0)
   r = b
   call innerproduct(r, r, info, rnorm(0))
   call allreduce0D('xy', rnorm(0))
   rnorm(0) = sqrt(rnorm(0))
   rnorm0 = rnorm(0)
   p(:) = 0.d0
   p(1) = rnorm(0)
   if (myid == 0) print*, rnorm0
   !
   ! GMRES niteration steps
   do iter = 1, niteration
	 call divide(r, rnorm(iter-1))
	 q(iter) = r
	 !
	 ! Calculate new r
         call apply_UT1000(w, ne, r, xpert)
	 do ie = 1, ne
	    call get_field(w(ie), ftmp, nlev)
	    call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call set_field(w(ie), ftmp, nlev)
	 end do
	 call apply_U1000(r, ne, w, xpert)
	 call allreduce_ens(r)
	 b = q(iter)
	 call multiply(b, lambda)
         call add(r, b)
	 !
	 ! Arnoldi (Gram-Schmidt): use b temporarily
	 do jter = 1, iter
	    call innerproduct(q(jter), r, info, hmat(jter,iter))
	    call allreduce0D('xy', hmat(jter,iter))
	    b = q(jter)
	    call multiply(b, hmat(jter,iter))
	    call subtract(r, b)
	 end do
	 !
	 ! Calculate new norm
	 call innerproduct(r, r, info, rnorm(iter))
	 call allreduce0D('xy', rnorm(iter))
	 rnorm(iter) = sqrt(rnorm(iter))
	 hmat(iter+1,iter) = rnorm(iter)
	 !
	 if (myid == 0) then
	    ! calculate new G and R
	    rmat(1:iter,iter) = hmat(1:iter,iter)
	    do jter = 1, iter-1
	       rmat1 = g(jter,1)*rmat(jter,iter) + g(jter,2)*rmat(jter+1,iter)
	       rmat2 = -g(jter,2)*rmat(jter,iter) + g(jter,1)*rmat(jter+1,iter)
	       rmat(jter,iter) = rmat1
	       rmat(jter+1,iter) = rmat2
	    end do
	    rmat1 = sqrt(rmat(iter,iter)**2+hmat(iter+1,iter)**2)
	    g(iter,1) = rmat(iter,iter)/rmat1
	    g(iter,2) = hmat(iter+1,iter)/rmat1
	    rmat(iter,iter) = rmat1
	    !
	    ! calculate new p
	    rmat1 = g(iter,1)*p(iter) + g(iter,2)*p(iter+1)
	    rmat2 = -g(iter,2)*p(iter) + g(iter,1)*p(iter+1)
	    p(iter) = rmat1
	    p(iter+1) = rmat2
	 end if
	 !
	 ! Check convergence
	 if (myid == 0) then
	    ratio = abs(p(iter+1))/rnorm0
	    if (ratio < tolerance .or. abs(p(iter+1)) < tolerance) then
	       iflag = 0
	    else
	       iflag = -1
	    end if
	    ! Monitoring
	    open(90, file="iteration.log", position="append")
	    if (iter == 1) write(90,'(a6,2a14)') "iter", "residual", "ratio"
	    write(90,'(i4,2e14.7)') iter, abs(p(iter+1)), ratio
	    close(90)
	 end if
	 call int_broadcast0D('all', iflag, 0)
	 if (iflag == 0) then
	    if (myid == 0) print*, 'Linear equation system has reached its tolerance.', abs(p(iter+1)), rnorm0
	    exit
	 end if
	 if (iter == niteration) then
	    if (myid == 0) print*, 'The iteration number exceeds.', abs(p(iter+1)), rnorm0
	    exit
	 end if
   end do
   !
   ! solve Rs = p
   if (myid == 0) then
      do jter = iter, 1, -1
	 s(jter) = p(jter)
	 do kter = iter, jter+1, -1
	    s(jter) = s(jter) - s(kter)*rmat(jter,kter)
	 end do
	 s(jter) = s(jter)/rmat(jter,jter)
      end do
   end if
   call broadcast1D('all', s(1:iter), 0, iter)
   !
   ! Calculate x
   do jter = 1, iter
      call multiply(q(jter), s(jter))
      call add(x, q(jter))
   end do
   !
   !
   !
   ! Output
   call write_control(x, 1, info, 'ms')
   !
   ! deallocation
   call destroy(x); call destroy(r); call destroy(b)
   do ie = 1, ne
      call destroy(xpert(ie))
      call destroy(w(ie))
   end do
   do iter = 1, niteration
      call destroy(q(iter))
   end do
   deallocate(xpert, w, q)
   call destroy(Slocx); call destroy(Slocy)
   deallocate(rnorm, p, s, hmat, rmat, g, ftmp)
   deallocate(y)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program sensitivity1000
