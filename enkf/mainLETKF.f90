program letkf
   use variable
   use nhmlib
   use enkflib
   use interpolate, only : setGridInf
   use NodeInfo_class
   use NodeObsSpaceControl_class
   use NodeObsValidControl_class
   use NodeObsVLocControl_class
   use NodeObsControl_class
   use NodeProfileControl_class
   use NodeControl_class
   use NodeBCControl_class
   use NodeHControl_class
   use NodePreH_class
   use Localization_class
   use NodeMPI
   implicit none
   !
   integer :: nproc, myid, nx, ny, nz, nt, ne, nlev, hrange
   integer :: ierror, i, j, k, it, ip, jp, iobs, jobs, ii, jj
   integer :: i0, j0, k0, it0
   real(r_size), dimension(:), allocatable :: vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw
   real(r_size), dimension(:,:), allocatable :: zs
   real(r_sngl), dimension(:,:), allocatable :: zstmp, landsea, lon, lat
   real(r_size), dimension(:,:,:), allocatable :: zs0
   !
   type(NodeBCControl) :: bc
   type(NodePreH) :: preH
   type(NodeHControl) :: h
   type(Localization) :: loc
   !
   type(NodeInfo) :: info
   type(NodeControl) :: xana, xinc, rho
   type(NodeControl), dimension(:), allocatable :: xapert
   type(NodeProfileControl) :: xbck, xobs, xpert
   type(NodeObsSpaceControl) :: obsspace, global_obsspace
   type(NodeObsValidControl) :: valid, validtmp
   type(NodeObsControl) :: obs, error, ybck, xyloc
   type(NodeObsControl), dimension(:), allocatable :: ypert
   type(NodeObsVLocControl) :: logp, corr, vproftmp, vloc
   type(NodeObsVLocControl), dimension(:), allocatable :: vprof
   !
   integer :: mmax, mobs, nrank, inflana, adapmode, inflpoly, nlevnhm, nxyt, ncol0, ncol, myidx, myidy, myide
   integer :: id, idx, idy, is, ie, js, je, dis, die, djs, dje, ixyt, ivar, ivarnhm, iter, icol, ia, ja
   real(r_size) :: rhoML, rhoMLold, dh
   logical, dimension(:,:,:), allocatable :: processed, updated
   integer, dimension(:), allocatable :: ix0, jy0, ix, jy, istart, iend, jstart, jend
   integer, dimension(:,:), allocatable :: nobs, k2ijt
   real(r_size), dimension(:), allocatable :: X, VTX, Y0, gammab, d, lambda, flambda, gammaa, gammainfl
   real(r_size), dimension(:,:), allocatable :: U, V, Y, YYT
   real(r_size), dimension(:,:,:), allocatable :: Ye
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
   !if (myid == 0) write(6, control_nl)
   close(10)
   hscale = 1000.d0*hscale/resolution
   hrange = int(floor(hscale))
   !
   ! Optimal decomposition
   call setGridInf(nx0, ny0, nz0, proj, resolution, slon, xi, xj, xlat, xlon, slat, slat2)
   call new(info, 1, nxpe, nype, nepe, myid, nx0, ny0, ne0, hrange)
   call initialize_mpi(info)
   call destroy(info)
   call new_conv(global_obsspace, nx0, ny0, nt0)
   call read_obs(global_obsspace, myid)
   allocate(istart(nxpe), iend(nxpe), jstart(nype), jend(nype))
   if (myid == 0) then
      allocate(nobs(nx0,ny0))
      do i = 1, nx0
      do j = 1, ny0
         call get_mobs(global_obsspace, i, j, mobs)
         nobs(i,j) = mobs
      end do
      end do
      call load_balance(nxpe, nype, hrange, nobs, istart, iend, jstart, jend, nx0, ny0)
      deallocate(nobs)
   end if
   call int_broadcast1D('all', istart, 0, nxpe)
   call int_broadcast1D('all', iend, 0, nxpe)
   call int_broadcast1D('all', jstart, 0, nype)
   call int_broadcast1D('all', jend, 0, nype)
   !
   ! MPI Decomposition
   call new(info, 0, nxpe, nype, nepe, myid, nx0, ny0, ne0, hrange, istart, iend, jstart, jend)
   call get_nx(info, nx)
   call get_ny(info, ny)
   call get_ne(info, ne)
   call get_myidx(info, myidx)
   call get_myidy(info, myidy)
   call get_myide(info, myide)
   call get_di(info, dis, die)
   call get_dj(info, djs, dje)
   nz = nz0
   nt = nt0
   !
   ! Grid metric
   allocate(vdz(nz), vrdz(nz), vrdz2(nz))
   allocate(zrp(nz), vctransp(nz), dvtransp(nz))
   allocate(zrw(nz), vctransw(nz), dvtransw(nz))
   if (myid == 0) then
      call vrgdis(nz, nz, 1, nz, nz, 1, dz2, dz1, dz2, vdz, vrdz, vrdz2)
      call setzrp(vrdz2, nz, zrp)
      call calc_zcoordinate(vctrans_type, n_vctrans, ztop, zl_vctrans, zh_vctrans, zrp, nz, vctransp, dvtransp)
      call setzrw(vrdz, nz, zrw)
      call calc_zcoordinate(vctrans_type, n_vctrans, real(zrw(nz-1),r_sngl), zl_vctrans, zh_vctrans, zrw, nz, vctransw, dvtransw)
   end if
   call broadcast1D('all', vdz, 0, nz)
   call broadcast1D('all', vrdz, 0, nz)
   call broadcast1D('all', vrdz2, 0, nz)
   call broadcast1D('all', zrp, 0, nz)
   call broadcast1D('all', vctransp, 0, nz)
   call broadcast1D('all', dvtransp, 0, nz)
   call broadcast1D('all', zrw, 0, nz)
   call broadcast1D('all', vctransw, 0, nz)
   call broadcast1D('all', dvtransw, 0, nz)
   !
   ! Const
   allocate(zs0(nx0,ny0,1), zs(nx,ny))
   if (myid == 0) then
      allocate(zstmp(nx0,ny0), landsea(nx0,ny0), lon(nx0,ny0), lat(nx0,ny0))
      call read_nhmcst('mfhm', nx0, ny0, lon, lat, zstmp, landsea)
      zs0(:,:,1) = zstmp
      deallocate(zstmp, landsea, lon, lat)
   end if
   if (myide == 0) call scatter(info, 'xy', zs0(:,:,1), 0, nx0, ny0, 1)
   call broadcast2D('e', zs0(1:nx,1:ny,1), 0, nx, ny)
   zs(1:nx,1:ny) = zs0(1:nx,1:ny,1)
   deallocate(zs0)
   !
   !
   !
   ! Observation space
   call new_conv(obsspace, nx, ny, nt0)
   call scatter_obs(global_obsspace, info, obsspace)
   ! Observation and observational errors
   call new(obs, obsspace)
   call new(error, obsspace)
   call set_obs(obs, obsspace)
   call set_error(error, obsspace)
   ! Observation validity
   call new(valid, obsspace)
   call new(validtmp, obsspace)
   ! Pre-observation operator
   call new_NodePreH(preH, .False., vdz, zrp, vctransp, dvtransp, dvtransw, zs, nx, ny, nz0, nt0)
   ! Observation operator
   call new(h, info, obsspace)
   call new(xyloc, obsspace, 2)
   call get_xyloc(h, info, xyloc)
   ! Observation perturbation
   allocate(ypert(ne))
   do ie = 1, ne
      call new(ypert(ie), obsspace)
   end do
   ! Temporary variables
   call new(ybck, obsspace)
   call new(ytmp, obsspace)
   !
   !
   !
   ! Model space
   ! Index for profiles
   allocate(updated(nx,ny,nt))
   nxyt = 0
   updated(:,:,:) = .False.
   do i = 1, nx-1
      do j = 1, ny-1
	 do it = 1, nt
	    call get_mobs(obsspace, i, j, it, mobs)
	    if (mobs == 0) cycle
	    do ii = i, i+1
	       do jj = j, j+1
		  if (updated(ii,jj,it)) cycle
		  nxyt = nxyt + 1
		  updated(ii,jj,it) = .True.
	       end do
	    end do
	 end do
      end do
   end do
   allocate(k2ijt(nxyt,3))
   ixyt = 0
   updated(:,:,:) = .False.
   do i = 1, nx-1
      do j = 1, ny-1
	 do it = 1, nt
	    call get_mobs(obsspace, i, j, it, mobs)
	    if (mobs == 0) cycle
	    do ii = i, i+1
	       do jj = j, j+1
		  if (updated(ii,jj,it)) cycle
	          ixyt = ixyt + 1
	          k2ijt(ixyt,1) = ii; k2ijt(ixyt,2) = jj; k2ijt(ixyt,3) = it
		  updated(ii,jj,it) = .True.
	       end do
	    end do
	 end do
      end do
   end do
   deallocate(updated)
   ! Allocate
   call new_nhm_NodeControl(xana, info, .True., nz0, 2)
   call new_nhm_NodeControl(xinc, info, .True., nz0, 2)
   allocate(xapert(ne))
   do ie = 1, ne
      call new_nhm_NodeControl(xapert(ie), info, .True., nz0, 2)
   end do
   call new_nhm_NodeProfileControl(xbck, nxyt, nz0, nx, ny)
   call new_obs_NodeProfileControl(xobs, nxyt, nz0, nx, ny)
   call new_nhm_NodeProfileControl(xpert, nxyt, nz0, nx, ny)
   ! Boundary conditions
   call new(bc, xana)
   ! Control
   call new(rho, control_mode, info, .True., nz0, nz0, 1)
   call new(loc, hscale, hscaleqv, rho, ybck)
   !
   !
   !
   ! Forecast perturbations Y
   ! Numbers of obs per column
   ncol0 = 0; mmax = 0
   do i = 1, nx-1
      do j = 1, ny-1
         call get_mobs(obsspace, i, j, mobs)
         if (mobs > 0) then
	    ncol0 = ncol0 + 1
	    mmax = mmax + mobs
	 end if
      end do
   end do
   allocate(processed(loc%nvar,nx,ny))
   processed(:,:,:) = .True.
   if (ncol0 > 0) then
      allocate(ix0(ncol0), jy0(ncol0), ix(ncol0), jy(ncol0))
      icol = 0
      do i = 1, nx-1
	 do j = 1, ny-1
	    call get_mobs(obsspace, i, j, mobs)
	    if (mobs > 0) then
	       icol = icol + 1
	       ix0(icol) = i; jy0(icol) = j
	    end if
	 end do
      end do
      do i = 1, ncol0
         processed(:,ix0(i),jy0(i)) = .False.
      end do
   end if
   !
   ! First vertical localization: set vertical perturbations
   if (myid == 0) print*, 'VLOC0'
   call new(vloc, obsspace, nz0, 1)
   call new(logp, obsspace, nz0, 1)
   call new(vproftmp, obsspace, nz0)
   call new(corr, obsspace, nz0)
   allocate(vprof(ne))
   call set_field(logp, 0.d0)
   call set_field(corr, 0.d0)
   call set_field(xbck, 0.d0)
   call set_field(valid, 1)
   call set_field(ybck, 0.d0)
   do ie = 1, ne
      call read_ensemble(xpert, info, ie, k2ijt, xapert(ie), 'fg', nt, ne, nxyt)
      call apply_preH(preH, control_mode, processed(1,:,:), k2ijt, xpert, xbck, xobs, nxyt)
      ! mean logp: use vloc as temporary variable
      call set_logp(vloc, info, processed(1,:,:), k2ijt, obsspace, xobs, nx, ny, nxyt)
      call add(logp, vloc)
      ! mean vprof: use corr as ensemble mean
      call new(vprof(ie), obsspace, nz0)
      call set_profile(vprof(ie), info, processed(1,:,:), k2ijt, obsspace, xobs, nx, ny, nxyt)
      call add(corr, vprof(ie))
      ! mean y: use ybck as ensemble mean
      call apply_H(h, info, processed(1,:,:), k2ijt, xobs, obsspace, validtmp, ypert(ie), nx, ny, nxyt)
      call and(valid, validtmp)
      call add_obs(ybck, ypert(ie))
   end do
   call allreduce_ensvloc(logp)
   call allreduce_ensvloc(corr)
   call divide_vloc(logp, 1.d0*ne0)
   call divide_vloc(corr, 1.d0*ne0)
   call allreduce_ensvalid(valid)
   call allreduce_ensobs(ybck)
   call divide_obs(ybck, 1.d0*ne0)
   do ie = 1, ne
      call subtract_vloc(vprof(ie), corr)
      call subtract_obs(ypert(ie), ybck)
   end do
   deallocate(processed)
   ! Quality control
   if (iflg_incremental == 1) call add_obs(obs, ybck)
   call subtract_obs(obs, ybck)
   call divide_obs(obs, error)
   call qccheck(obs, qcthreshold, validtmp)
   call and(valid, validtmp)
   !
   ! Second vertical localization: normalize perturbations
   if (myid == 0) print*, 'VLOC1'
   call set_field(corr, 0.d0)
   call set_field(ybck, 0.d0)
   do ie = 1, ne
      vproftmp = vprof(ie)
      call power_vloc(vproftmp, 2.d0)
      call add_vloc(corr, vproftmp)
      ytmp = ypert(ie)
      call power_obs(ytmp, 2.d0)
      call add_obs(ybck, ytmp)
   end do
   call allreduce_ensvloc(corr)
   call divide_vloc(corr, 1.d0*(ne0-1))
   call power_vloc(corr, 0.5d0)
   call allreduce_ensobs(ybck)
   call divide_obs(ybck, 1.d0*(ne0-1))
   call power_obs(ybck, 0.5d0)
   do ie = 1, ne
      call divide_vloc(vprof(ie), corr)
      call divide_obs(ypert(ie), ybck)
   end do
   !
   ! Third vertical localization: correlation
   if (myid == 0) print*, 'VLOC2'
   call set_field(corr, 0.d0)
   do ie = 1, ne
      call multiply_vloc(vprof(ie), ypert(ie))
      call add_vloc(corr, vprof(ie))
   end do
   call allreduce_ensvloc(corr)
   call divide_vloc(corr, 1.d0*(ne0-1))
   call ecorap(corr, vscale, logp, vloc)
   call destroy(vproftmp)
   call destroy(corr)
   call destroy(logp)
   do ie = 1, ne
      call destroy(vprof(ie))
   end do
   deallocate(vprof)
   ! restore perturbations in observation space
   do ie = 1, ne
      call multiply_obs(ypert(ie), ybck)
      call divide_obs(ypert(ie), sqrt(1.d0*(ne0-1)))
      call divide_obs(ypert(ie), error)
   end do
   if (myid == 0) print*, 'VLOC3'
   !
   !
   !
   ! Forecast perturbations X
   ! X
   call set_const(xana, 0.d0)
   do ie = 1, ne
      call add(xana, xapert(ie))
   end do
   call allreduce_ens(xana)
   call divide(xana, 1.d0*ne0)
   do ie = 1, ne
      call subtract(xapert(ie), xana)
      call divide(xapert(ie), sqrt(1.d0*(ne0-1)))
      call apply_bc(bc, xapert(ie))
   end do
   !
   !
   !
   ! Pre EnSRF
   inflana = inflmode/100
   adapmode = mod(inflmode,100)/10
   inflpoly = mod(inflmode,10)
   ! Allocate
   if (ncol0 > 0) then
      allocate(X(ne), Y(mmax,ne), Y0(mmax))
      allocate(Ye(mmax,ne,nepe))
      allocate(YYT(ne0,ne0), gammab(ne0), U(ne0,ne0), V(ne,ne0))
      allocate(d(ne0), lambda(ne0), flambda(ne0), gammaa(ne0), gammainfl(ne0), VTX(ne0))
   end if
   call set_const(rho, 1.d0)
   call set_const(xinc, 0.d0)
   !
   !
   !
   ! EnSRF
   do i = 1+dis, nx-die
   do j = 1+djs, ny-dje
   !do i = 1+dis, 1+dis
   !do j = 1+djs, 1+djs
      ! Find relevant obs
      ncol = 0
      do icol = 1, ncol0
	 dh = sqrt(1.d0*(ix0(icol)-i)**2+1.d0*(jy0(icol)-j)**2)
	 if (dh > hscale) cycle
	 ncol = ncol + 1
	 ix(ncol) = ix0(icol); jy(ncol) = jy0(icol)
      end do
      if (ncol == 0) cycle
      !
      ia = i-dis; ja = j-djs
      do ivar = 1, loc%nvar
         call get_nlev(rho, ivar, nlev)
	 if (nlev > 1) nlev = nlev-1
	 do k = 1, nlev
	    ! Extract Y
	    call get_field(vloc, k, ytmp)
	    do ie = 1, ne
	       call extract_y(ypert(ie), i, j, ix(1:ncol), jy(1:ncol), loc.correlated(ivar,:), obsspace, valid, xyloc, ytmp, mobs, Y(:,ie), ncol, mmax)
	    end do
	    call extract_y(obs, i, j, ix(1:ncol), jy(1:ncol), loc.correlated(ivar,:), obsspace, valid, xyloc, ytmp, mobs, Y0, ncol, mmax)
	    if (mobs == 0) cycle
	    !if (myid == 41) print*, myid, minval(Y0(1:mobs)), maxval(Y0(1:mobs)), minval(Y(1:mobs,:)), maxval(Y(1:mobs,:))
	    !
	    !
	    !
	    ! Eigen-decomposition
	    if (mobs <= ne0) then
	       ! Covariance
	       do iobs = 1, mobs
	          do jobs = iobs, mobs
		     YYT(iobs,jobs) = sum(Y(iobs,:)*Y(jobs,:))
		     if (jobs > iobs) YYT(jobs,iobs) = YYT(iobs,jobs)
	          end do
	       end do
	       call allreduce2D('e', YYT(1:mobs,1:mobs), mobs, mobs)
	       do iobs = 1, mobs
	          YYT(iobs,iobs) = YYT(iobs,iobs) + 1.d0 
	       end do
	       call mtx_eigen(1, mobs, YYT(1:mobs,1:mobs), gammab(1:mobs), U(1:mobs,1:mobs), nrank)
	       gammab(1:mobs) = gammab(1:mobs) - 1.d0
	       ! Find rank
	       do iobs = 1, mobs
	          if (gammab(iobs) < 1.e-12) exit
	       end do
	       nrank = iobs - 1
	       if (nrank == 0) cycle
	       gammab(1:nrank) = sqrt(gammab(1:nrank))
	       ! V
	       do ie = 1, nrank
	          do je = 1, ne
		     V(je,ie) = sum(Y(1:mobs,je)*U(1:mobs,ie))/gammab(ie)
	          end do
	       end do
	       ! Project Y0 on U
	       do iobs = 1, nrank
	          d(iobs) = sum(U(1:mobs,iobs)*Y0(1:mobs))
	       end do
	    else
	       call MPI_BARRIER(MPI_COMM_EWORLD, ierror)
               call MPI_GATHER(Y(1:mobs,:), mobs*ne, r_type, Ye(1:mobs,:,:), mobs*ne, r_type, 0, MPI_COMM_EWORLD, ierror)
	       ! Covariance
	       if (myide == 0) then
	          do ie = 1, ne0
	             idx = (ie-1)/ne+1; is = mod(ie-1,ne)+1
		     do je = ie, ne0
		        idy = (je-1)/ne+1; js = mod(je-1,ne)+1
		        YYT(ie,je) = sum(Ye(1:mobs,is,idx)*Ye(1:mobs,js,idy))
		        if (je > ie) YYT(je,ie) = YYT(ie,je)
		     end do
	          end do
	          do ie = 1, ne0
		     YYT(ie,ie) = YYT(ie,ie) + 1.d0
	          end do
	          call mtx_eigen(1, ne0, YYT, gammab, U, nrank) ! use U as V temporarily
	          gammab(:) = gammab(:) - 1.d0
	          ! Find rank
	          do ie = 1, ne0
		     if (gammab(ie) < 1.e-12) exit
	          end do
	          nrank = ie - 1
	          gammab(1:nrank) = sqrt(gammab(1:nrank))
	          !if (myide == 0) print*, myid, minval(YYT), maxval(YYT), gammab(1)
	       end if
	       call int_broadcast0D('e', nrank, 0)
	       if (nrank == 0) cycle
	       call broadcast1D('e', gammab(1:nrank), 0, nrank)
	       ! V
	       call broadcast2D('e', U(:,1:nrank), 0, ne0, nrank)
	       is = myide*ne + 1; ie = is + ne - 1
	       V(:,1:nrank) = U(is:ie,1:nrank)
	       ! Project Y0 on U
	       do ie = 1, ne
	          X(ie) = sum(Y(1:mobs,ie)*Y0(1:mobs))
	       end do
	       do ie = 1, nrank
	          d(ie) = sum(V(:,ie)*X(:))/gammab(ie)
	       end do
	       call allreduce1D('e', d(1:nrank), nrank)
	    end if
	    !
	    !
	    !
	    ! Inflation functions
	    lambda(1:nrank) = 1.d0/sqrt(1.d0+gammab(1:nrank)**2)
	    !if (myid == 41) print*, i, j, mobs, nrank, lambda(1), d(1)
	    ! Non-parametric prior
	    if (adapmode == 1) then 
	       do iobs = 1, nrank
	          if (abs(d(iobs)) > 1.d0) then
		     gammainfl(iobs) = sqrt(d(iobs)**2-1.d0)
	          else
		     gammainfl(iobs) = 0.d0
	          end if
	          gammainfl(iobs) = max(gammab(iobs),gammainfl(iobs))
	          if (iobs > 1) then
		     if (gammainfl(iobs) > gammainfl(iobs-1)) gammainfl(iobs) = gammainfl(iobs-1)
	          end if
	          flambda(iobs) = gammainfl(iobs)/(gammab(iobs)*sqrt(1.d0+gammainfl(iobs)**2))
	       end do
	       rho%control(ivar)%p%field(ia,ja,1,1) = gammainfl(1)/gammab(1)
	    ! Non-parametric posterior
	    else if (adapmode == 2) then
	       do iobs = 1, nrank
	          gammainfl(iobs) = gammab(iobs)*abs(d(iobs))*lambda(iobs)**2
	          gammaa(iobs) = gammab(iobs)*lambda(iobs)
	          gammainfl(iobs) = max(gammaa(iobs),gammainfl(iobs))
	          if (iobs > 1) then
		     if (gammainfl(iobs) > gammainfl(iobs-1)) gammainfl(iobs) = gammainfl(iobs-1)
	          end if
	          flambda(iobs) = gammainfl(iobs)/gammab(iobs)
	       end do
	       rho%control(ivar)%p%field(ia,ja,1,1) = gammainfl(1)/gammaa(1)
	    ! EM algorithm: single rho
	    else if (adapmode == 3) then
	       inflana = 1
	       rhoML = 1.d0
	       do iter = 1, niteration
	          rhoMLold = rhoML
	          gammainfl(1:nrank) = sqrt(rhoML)*gammab(1:nrank)
	          flambda(1:nrank) = 1.d0/sqrt(1.d0+gammainfl(1:nrank)**2)
	          rhoML = rhoML*sum(flambda(1:nrank)**2+gammainfl(1:nrank)**2*flambda(1:nrank)**4*d(1:nrank)**2)/nrank
	          if (rhoML <= 1.d0) then
		     rhoML = rhoMLold
		     exit
	          end if
	       end do
	       gammainfl(1:nrank) = sqrt(rhoML)*gammab(1:nrank)
	       flambda(1:nrank) = 1.d0/sqrt(1.d0+gammainfl(1:nrank)**2) 
	       flambda(1:nrank) = gammainfl(1:nrank)*flambda(1:nrank)/gammab(1:nrank)
	       rho%control(ivar)%p%field(ia,ja,1,1) = sqrt(rhoML)
	    ! EM algorithm: multiple rho
	    else if (adapmode == 4) then
	       VTX(1:nrank) = 1.d0 ! usr VTX, gammaa as rhoML, rhoMLold
	       do iter = 1, niteration
	          gammaa(1:nrank) = VTX(1:nrank)
	          gammainfl(1:nrank) = sqrt(VTX(1:nrank))*gammab(1:nrank)
	          flambda(1:nrank) = 1.d0/sqrt(1.d0+gammainfl(:)**2)
	          VTX(1:nrank) = VTX(1:nrank)*(flambda(1:nrank)**2+gammainfl(1:nrank)**2*flambda(1:nrank)**4*d(1:nrank)**2)
	          if (maxval(VTX(1:nrank)) <= 1.d0) then
		     VTX(1:nrank) = gammaa(1:nrank)
		     exit
	          end if
	       end do
	       do iobs = 1, nrank
	          gammainfl(iobs) = sqrt(VTX(iobs))*gammab(iobs)
	          gammainfl(iobs) = max(gammab(iobs),gammainfl(iobs))
	          if (iobs > 1) then
		     if (gammainfl(iobs) > gammainfl(iobs-1)) gammainfl(iobs) = gammainfl(iobs-1)
	          end if
	          flambda(iobs) = gammainfl(iobs)/(gammab(iobs)*sqrt(1.d0+gammainfl(iobs)**2))
	       end do
	       rho%control(ivar)%p%field(ia,ja,1,1) = sqrt(VTX(1))
	    end if
	    !
	    !
	    !
	    ! Update X
	    if (inflana == 1) then
	       lambda(1:nrank) = gammab(1:nrank)*flambda(1:nrank)**2*d(1:nrank)
	    else
	       lambda(1:nrank) = gammab(1:nrank)*lambda(1:nrank)**2*d(1:nrank)
	    end if
	    flambda(1:nrank) = (1.d0-flambda(1:nrank))
	    do iter = 1, rho%nvarnhm(ivar)
	       ivarnhm = rho%ivarnhm(ivar,iter)
	       ! Project X on V
	       do ie = 1, ne
		  X(ie) = xapert(ie)%control(ivarnhm)%p%field(ia,ja,k,2)
	       end do
	       do ie = 1, nrank
		  VTX(ie) = sum(V(:,ie)*X(:))
	       end do
	       call allreduce1D('e', VTX(1:nrank), nrank)
	       ! Analysis
	       xinc%control(ivarnhm)%p%field(ia,ja,k,2) = sum(VTX(1:nrank)*lambda(1:nrank))
	       ! Analysis perturbations
	       VTX(1:nrank) = flambda(1:nrank)*VTX(1:nrank)
	       X(1:ne) = 0.d0
	       do ie = 1, nrank
		  X(1:ne) = X(1:ne) + VTX(ie)*V(:,ie)
	       end do
	       do ie = 1, ne
		  xapert(ie)%control(ivarnhm)%p%field(ia,ja,k,2) = xapert(ie)%control(ivarnhm)%p%field(ia,ja,k,2) - X(ie)
	       end do
	       !
	       if (itout == nt) cycle
	       ! Project X on V
	       do ie = 1, ne
		  X(ie) = xapert(ie)%control(ivarnhm)%p%field(ia,ja,k,1)
	       end do
	       do ie = 1, nrank
		  VTX(ie) = sum(V(:,ie)*X(:))
	       end do
	       call allreduce1D('e', VTX(1:nrank), nrank)
	       ! Analysis
	       xinc%control(ivarnhm)%p%field(ia,ja,k,1) = sum(VTX(1:nrank)*lambda(1:nrank))
	    end do
	 end do
      end do
   end do
   end do
   !
   !
   !
   ! Output
   call apply_bc(bc, xinc)
   if (iflg_outana == 0) then
      if (itout == nt) then
         call write_control(xinc, 2, info, 'inc')
      else
         call write_control(xinc, 1, info, 'inc')
      end if
   else
      call add(xana, xinc)
      if (itout == nt) then
         call write_control(xana, 2, info, 'ana')
      else
         call write_control(xana, 1, info, 'ana')
      end if
   end if
   call write_control(rho, 1, info, 'rho')
   do ie = 1, ne
      call multiply(xapert(ie), sqrt(1.d0*(ne0-1)))
      call apply_bc(bc, xapert(ie))
   end do
   call write_ensemble(xapert, 2, imember1, imember2, info, 'pert', ne)
   !
   ! deallocation
   call destroy(xana); call destroy(xinc); call destroy(rho)
   call destroy(xbck); call destroy(xobs); call destroy(xpert)
   do ie = 1, ne
      call destroy(xapert(ie))
      call destroy(ypert(ie))
   end do
   deallocate(xapert, ypert)
   call destroy(obsspace); call destroy(global_obsspace)
   call destroy(valid); call destroy(validtmp)
   call destroy(obs); call destroy(error); call destroy(ybck)
   call destroy(h); call destroy(preH)
   if (ncol0 > 0) then
      deallocate(Ye)
      deallocate(k2ijt, ix0, jy0, ix, jy)
      deallocate(X, Y, YYT, VTX, gammab, U, V)
      deallocate(Y0, d, lambda, flambda, gammaa, gammainfl)
   end if
   deallocate(vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw, zs)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program letkf
