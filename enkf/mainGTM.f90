program GTM
   use variable
   use enkflib
   use nhmlib
   use NodeInfo_class
   use NodeControl_class
   use NodeMPI
   implicit none
   !
   character(2) :: header
   integer, parameter :: nlow = 2
   integer :: nproc, myid, np, nx, ny, nz, nt, ne, nrank, myide, nhigh, nlatent, nrbf, nlcz
   integer :: ierror, ip, jp, is, ie, irank, i, j, k, ilatent, irbf, iter, monitoring
   integer, dimension(:), allocatable :: i0, j0
   real(r_size) :: norm, drbf, sigmarbf, dx, dy, xs, ys, dist, dmin, alpha, beta, logp, logp0, logp1, logp2, logp3, logp4
   real(r_size), dimension(:), allocatable :: ymean, ystd, xv, Elcz, ygtm, G, E
   real(r_sngl), dimension(:,:), allocatable :: ytmp
   real(r_size), dimension(:,:), allocatable :: y0, y, u, urbf, phi, wy, Vylcz, Tlcz, Vlcz, D, R, FGF, V, FTR, upost, ypost
   !
   type(NodeInfo) :: info
   type(NodeControl), dimension(:), allocatable :: xnhm, xpert, wx, Vxlcz, xpost
   type(NodeControl) :: xmean, xstd, xtmp, xgtm
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
   ! GTM parameters
   nlatent = ngroup
   alpha = qcthreshold
   drbf = hscale
   sigmarbf = vscale
   monitoring = inflmode
   nlcz = 3*(nlow+1)
   ! Latent grid
   allocate(u(nlow,nlatent**nlow))
   dx = 2.d0/(nlatent-1); dy = 2.d0/(nlatent-1)
   xs = -1.d0; ys = -1.d0
   do j = 1, nlatent
      do i = 1, nlatent
	 ilatent = (j-1)*nlatent + i-1 + 1
	 u(1,ilatent) = xs + (i-1)*dx
	 u(2,ilatent) = ys + (j-1)*dy
      end do
   end do
   !
   ! RBF grid
   drbf = drbf*2.d0/(nlatent-1)
   nrbf = 2*floor(1.d0/drbf) + 1
   allocate(urbf(nlow,nrbf**nlow))
   dx = drbf; dy = drbf
   xs = -(nrbf-1)/2*dx; ys = -(nrbf-1)/2*dy
   do j = 1, nrbf
      do i = 1, nrbf
	 irbf = (j-1)*nrbf + i-1 + 1
	 urbf(1,irbf) = xs + (i-1)*dx
	 urbf(2,irbf) = ys + (j-1)*dy
      end do
   end do
   nlatent = nlatent**nlow
   nrbf = nrbf**nlow + nlow + 1 ! nlow: linear, 1: bias
   !
   ! Phi
   sigmarbf = 1.d0/(sigmarbf*drbf)**2
   allocate(phi(nrbf,nlatent))
   do k = 1, nlatent
      do i = 1, nrbf-nlow-1
	 dist = sum((u(:,k)-urbf(:,i))**2)
	 phi(i,k) = exp(-0.5d0*sigmarbf*dist)
      end do
      do i = 1, nlow
	 phi(nrbf-nlow-1+i,k) = u(i,k)
      end do
      phi(nrbf,k) = 1.d0
   end do
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
      !print*, np, i0(1), j0(1), y0(1,:)
      do ip = 1, np
         ymean(ip) = sum(y0(ip,:))/ne0
	 y0(ip,:) = y0(ip,:) - ymean(ip)
	 ystd(ip) = sqrt(sum(y0(ip,:)**2)/(ne0-1))
	 if (ystd(ip) < 1.e-12) then
	    y0(ip,:) = 0.d0
	 else
	    y0(ip,:) = y0(ip,:)/ystd(ip)
	 end if
      end do
      deallocate(ytmp)
   end if
   call int_broadcast0D('all', np, 0)
   if (myid > 0) allocate(i0(np), j0(np), ymean(np), ystd(np), y0(np,ne0))
   allocate(y(np,ne), wy(np,nrbf), ygtm(np))
   call int_broadcast1D('all', i0, 0, np)
   call int_broadcast1D('all', j0, 0, np)
   call broadcast1D('all', ymean, 0, np)
   call broadcast1D('all', ystd, 0, np)
   call broadcast2D('all', y0, 0, np, ne0)
   is = myide*ne + 1; ie = is + ne - 1
   y(:,1:ne) = y0(:,is:ie)
   wy(:,:) = 0.d0
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
   call new_sens1_NodeControl(xgtm, info, nz0, nt0, xnhm(1))
   allocate(Vxlcz(nlcz))
   do ie = 1, nlcz
      call new_sens1_NodeControl(Vxlcz(ie), info, nz0, nt0, xnhm(1))
   end do
   allocate(wx(nrbf))
   do irbf = 1, nrbf
      call new_sens1_NodeControl(wx(irbf), info, nz0, nt0, xnhm(1))
      call set_const(wx(irbf), 0.d0)
   end do
   allocate(xpost(nlatent))
   do ilatent = 1, nlatent
      call new_sens1_NodeControl(xpost(ilatent), info, nz0, nt0, xnhm(1))
   end do
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
   end do
   !
   !
   !
   !
   ! Lanczos
   allocate(xv(ne), Vylcz(np,nlcz), Tlcz(nlcz,nlcz))
   Tlcz(:,:) = 0.d0
   call randomize(Vxlcz(1), info)
   if (myid == 0) call random_number(Vylcz(:,1))      
   call broadcast1D('all', Vylcz(:,1), 0, np)
   call innerproduct(Vxlcz(1), Vxlcz(1), info, norm)
   call allreduce0D('xy', norm)
   norm = norm + sum(Vylcz(:,1)**2)
   norm = sqrt(norm)
   call divide(Vxlcz(1), norm)
   Vylcz(:,1) = Vylcz(:,1)/norm
   !
   ! XXTv
   call set_const(wx(1), 0.d0) ! use wx(1) and wy(1) temporarily as w
   wy(:,1) = 0.d0
   do ie = 1, ne
      call innerproduct(xpert(ie), Vxlcz(1), info, xv(ie))
      call allreduce0D('xy', xv(ie))
      xv(ie) = xv(ie) + sum(y(:,ie)*Vylcz(:,1))
   end do
   do ie = 1, ne  
      xtmp = xpert(ie)
      call multiply(xtmp, xv(ie))
      call add(wx(1), xtmp)
      wy(:,1) = wy(:,1) + xv(ie)*y(:,ie)
   end do
   call divide(wx(1), 1.d0*(ne0-1))
   wy(:,1) = wy(:,1)/(ne0-1)
   call allreduce_ens(wx(1))
   call allreduce1D('e', wy(:,1), np)
   call innerproduct(wx(1), Vxlcz(1), info, Tlcz(1,1))
   call allreduce0D('xy', Tlcz(1,1))
   Tlcz(1,1) = Tlcz(1,1) + sum(wy(:,1)*Vylcz(:,1))
   xtmp = Vxlcz(1)
   call multiply(xtmp, Tlcz(1,1))
   call subtract(wx(1), xtmp)
   wy(:,1) = wy(:,1) - Tlcz(1,1)*Vylcz(:,1)
   if (myid == 0) print*, 'iter = 1', Tlcz(1,1)
   !
   do jp = 2, nlcz
      call innerproduct(wx(1), wx(1), info, norm)
      call allreduce0D('xy', norm)
      norm = norm + sum(wy(:,1)**2)
      norm = sqrt(norm)
      Tlcz(jp,jp-1) = norm
      Vxlcz(jp) = wx(1)
      call divide(Vxlcz(jp), norm)
      Vylcz(:,jp) = wy(:,1)/norm
      !
      ! XXTv
      call set_const(wx(1), 0.d0)
      wy(:,1) = 0.d0
      do ie = 1, ne
         call innerproduct(xpert(ie), Vxlcz(jp), info, xv(ie))
         call allreduce0D('xy', xv(ie))
         xv(ie) = xv(ie) + sum(y(:,ie)*Vylcz(:,jp))
      end do
      do ie = 1, ne  
         xtmp = xpert(ie)
         call multiply(xtmp, xv(ie))
         call add(wx(1), xtmp)
         wy(:,1) = wy(:,1) + xv(ie)*y(:,ie)
      end do
      call divide(wx(1), 1.d0*(ne0-1))
      wy(:,1) = wy(:,1)/(ne0-1)
      call allreduce_ens(wx(1))
      call allreduce1D('e', wy(:,1), np)
      call innerproduct(wx(1), Vxlcz(jp), info, Tlcz(jp,jp))
      call allreduce0D('xy', Tlcz(jp,jp))
      Tlcz(jp,jp) = Tlcz(jp,jp) + sum(wy(:,1)*Vylcz(:,jp))
      xtmp = Vxlcz(jp)
      call multiply(xtmp, Tlcz(jp,jp))
      call subtract(wx(1), xtmp)
      wy(:,1) = wy(:,1) - Tlcz(jp,jp)*Vylcz(:,jp)
      xtmp = Vxlcz(jp-1)
      call multiply(xtmp, Tlcz(jp,jp-1))
      call subtract(wx(1), xtmp)
      wy(:,1) = wy(:,1) - Tlcz(jp,jp-1)*Vylcz(:,jp-1)
      !
      ! Reorthogonalization
      do ip = 1, jp
         call innerproduct(wx(1), Vxlcz(ip), info, norm)
         call allreduce0D('xy', norm)
         norm = norm + sum(wy(:,1)*Vylcz(:,ip))
	 xtmp = Vxlcz(ip)
         call multiply(xtmp, norm)
         call subtract(wx(1), xtmp)
         wy(:,1) = wy(:,1) - norm*Vylcz(:,ip)
      end do
      Tlcz(jp-1,jp) = Tlcz(jp,jp-1)
      if (myid == 0) print*, 'iter = ', jp, Tlcz(jp,jp)
   end do
   !
   !
   !
   ! Eigen-decomposition for Tlcz
   allocate(Elcz(nlcz), Vlcz(nlcz,nlcz))
   call mtx_eigen(1, nlcz, Tlcz, Elcz, Vlcz, nrank)
   where (Elcz < 0.) Elcz = 0.
   Elcz = sqrt(Elcz)
   nrank = nlow
   if (myid == 0) print*, nrank, Elcz(1:nrank)
   !
   !
   !
   ! Firstguesses from the eigenvalues and eigenvectors of Lanczos
   call set_const(wx(1), 0.d0)
   wy(:,1) = 0.d0
   do irank = 1, nrank
      irbf = nrbf-nlow-1+irank
      call set_const(wx(irbf), 0.d0)
      wy(:,irbf) = 0.d0
      do ip = 1, nlcz
         xtmp = Vxlcz(ip)
         call multiply(xtmp, Vlcz(ip,irank))
         call add(wx(irbf), xtmp)
         wy(:,irbf) = wy(:,irbf) + Vlcz(ip,irank)*Vylcz(:,ip)
      end do
      call multiply(wx(irbf), Elcz(irank))
      wy(:,irbf) = Elcz(irank)*wy(:,irbf)
   end do
   beta = 1.d0/Elcz(nlow+1)**2
   !
   !
   !
   ! GTM
   call get_ndim(xgtm, nhigh)
   allocate(D(nlatent,ne), R(nlatent,ne), G(nlatent))
   allocate(FGF(nrbf,nrbf), E(nrbf), V(nrbf,nrbf), FTR(nrbf,ne))
   do k = 1, nlatent
      call set_const(xgtm, 0.d0)
      ygtm(:) = 0.d0
      do i = 1, nrbf
         xtmp = wx(i)
         call multiply(xtmp, phi(i,k))
         call add(xgtm, xtmp)
         ygtm(:) = ygtm(:) + phi(i,k)*wy(:,i)
      end do
      do ie = 1, ne
         xtmp = xgtm
         call subtract(xtmp, xpert(ie))
	 call innerproduct(xtmp, xtmp, info, norm)
         call allreduce0D('xy', norm)
	 norm = norm + sum((ygtm(:)-y(:,ie))**2)
         D(k,ie) = norm
      end do
   end do
   !
   do iter = 1, niteration
      ! R
      do ie = 1, ne
	 dmin = minval(D(:,ie))
	 R(:,ie) = exp(-0.5d0*beta*(D(:,ie)-dmin))
	 norm = sum(R(:,ie))
	 R(:,ie) = R(:,ie)/norm
      end do
      do k = 1, nlatent
	 G(k) = sum(R(k,:))
      end do
      call allreduce1D('e', G, nlatent)
      !
      ! Log-likelihood
      if (monitoring == 1) then
         logp0 = log(1.d0/nlatent)
	 logp1 = 0.5d0*nhigh*log(0.5d0*beta/pi)
	 logp2 = -0.5d0*beta*sum(R*D)
	 call allreduce0D('e', logp2)
	 logp3 = -sum(log(R**R))
	 call allreduce0D('e', logp3)
	 logp4 = 0.d0
	 do i = 1, nrbf
	    call innerproduct(wx(i), wx(i), info, norm)
            call allreduce0D('xy', norm)
	    norm = norm + sum(wy(:,i)**2)
	    logp4 = logp4 + norm
	 end do
	 logp4 = -0.5d0*alpha*logp4
	 logp = ne0*(logp0+logp1) + logp2 + logp3 + logp4
	 if (myid == 0) print*, iter, logp, beta
      end if
      !
      ! Solve W
      ! FGF
      do i = 1, nrbf
	 do j = i, nrbf
	    FGF(i,j) = sum(phi(i,:)*phi(j,:)*G(:))
	    if (j == i) FGF(i,j) = FGF(i,j) + alpha/beta
	    if (j > i) FGF(j,i) = FGF(i,j)
	 end do
      end do
      call mtx_eigen(1, nrbf, FGF, E, V, nrank)
      !if (myid == 0) print*, nrank, nrbf, E
      !where (E < 0.d0) E = 0.d0
      ! FTR
      do i = 1, nrbf
	 do ie = 1, ne
	    FTR(i,ie) = sum(phi(i,:)*R(:,ie))
	 end do
      end do
      do ie = 1, ne
         do irank = 1, nrank
	    FGF(irank,1) = sum(V(:,irank)*FTR(:,ie))/E(irank) ! use FGF(:,1) temporarily
	 end do
	 FTR(:,ie) = 0.d0
	 do irank = 1, nrank
	    FTR(:,ie) = FTR(:,ie) + FGF(irank,1)*V(:,irank)
	 end do
      end do
      ! Compute W
      do i = 1, nrbf
         call set_const(wx(i), 0.d0)
         wy(:,i) = 0.d0
	 do ie = 1, ne
            xtmp = xpert(ie)
            call multiply(xtmp, FTR(i,ie))
            call add(wx(i), xtmp)
            wy(:,i) = wy(:,i) + FTR(i,ie)*y(:,ie)
	 end do
	 call allreduce_ens(wx(i))
         call allreduce1D('e', wy(:,i), np)
      end do
      !
      ! Solve beta
      do k = 1, nlatent
         call set_const(xgtm, 0.d0)
         ygtm(:) = 0.d0
         do i = 1, nrbf
            xtmp = wx(i)
            call multiply(xtmp, phi(i,k))
            call add(xgtm, xtmp)
            ygtm(:) = ygtm(:) + phi(i,k)*wy(:,i)
         end do
         do ie = 1, ne
	    xtmp = xgtm
            call subtract(xtmp, xpert(ie))
	    call innerproduct(xtmp, xtmp, info, norm)
            call allreduce0D('xy', norm)
	    norm = norm + sum((ygtm(:)-y(:,ie))**2)
            D(k,ie) = norm
         end do
      end do
      norm = sum(R*D)
      call allreduce0D('e', norm)
      beta = ne0*nhigh/norm
   end do
   !
   !
   !
   ! Post GTM
   allocate(upost(nlow,ne0), ypost(np,nlatent))
   upost(:,:) = 0.d0
   ! D
   do k = 1, nlatent
      call set_const(xpost(k), 0.d0)
      ypost(:,k) = 0.d0
      do i = 1, nrbf
         xtmp = wx(i)
         call multiply(xtmp, phi(i,k))
         call add(xpost(k), xtmp)
         ypost(:,k) = ypost(:,k) + phi(i,k)*wy(:,i)
      end do
      do ie = 1, ne
         xtmp = xpost(k)
         call subtract(xtmp, xpert(ie))
	 call innerproduct(xtmp, xtmp, info, norm)
         call allreduce0D('xy', norm)
	 norm = norm + sum((ygtm(:)-y(:,ie))**2)
         D(k,ie) = norm
      end do
   end do
   ! R
   is = myide*ne + 1
   do ie = 1, ne
      dmin = minval(D(:,ie))
      R(:,ie) = exp(-0.5d0*beta*(D(:,ie)-dmin))
      norm = sum(R(:,ie))
      R(:,ie) = R(:,ie)/norm
      do i = 1, nlow
	 upost(i,ie+is-1) = sum(R(:,ie)*u(i,:))
      end do
   end do
   call allreduce2D('e', upost, nlow, ne0)
   !
   !
   !
   ! Output
   do i = 1, nrbf
      write(header(1:2),'(I2.2)') i
      !call multiply(wx(i), xstd)
      !call add(wx(i), xmean)
      call write_control(wx(i), 1, info, 'wx'//header)
   end do
   do i = 1, nlatent
      write(header(1:2),'(I2.2)') i
      call multiply(xpost(i), xstd)
      !call add(xpost(i), xmean)
      call write_control(xpost(i), 1, info, 'xpost'//header)
      ypost(:,i) = ypost(:,i)*ystd
      !ypost(:,i) = ypost(:,i) + ymean
   end do
   call write_control(xmean, 1, info, 'xmean')
   call write_control(xstd, 1, info, 'xstd')
   if (myid == 0) then
      open(90, file='ypost', form='binary')
      write(90) np, nrbf, nlatent
      write(90) ymean, ystd
      write(90) wy
      write(90) ypost
      write(90) upost
      close(90)
   end if
   !
   ! deallocation
   call destroy(xmean); call destroy(xstd); call destroy(xtmp); call destroy(xgtm)
   do i = 1, nlcz
      call destroy(Vxlcz(i))
   end do
   do irbf = 1, nrbf
      call destroy(wx(irbf))
   end do
   do ilatent = 1, nlatent
      call destroy(xpost(ilatent))
   end do
   do ie = 1, ne
      call destroy(xpert(ie))
   end do
   deallocate(Vxlcz, wx, xpost, xpert)
   deallocate(i0, j0, ymean, ystd, xv, Elcz, ygtm, G, E)
   deallocate(y, u, urbf, phi, wy, Vylcz, Vlcz, Tlcz, D, R, FGF, V, FTR, upost, ypost)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program GTM
