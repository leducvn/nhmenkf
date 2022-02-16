program etkf
   use variable
   use nhmlib
   use enkflib
   use rttovlib
   use interpolate, only : setGridInf
   use NodeInfo_class
   use NodeObsSpaceControl_class
   use NodeObsValidControl_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeBCControl_class
   use NodeHControl_class
   use NodePreH_class
   use NodeMPI
   use eigen_libs_mod
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
   !
   type(NodeInfo) :: info
   type(NodeControl) :: xtmp, xbck, xobs, rho
   type(NodeControl), dimension(:), allocatable :: xpert
   type(NodeObsSpaceControl) :: obsspace, global_obsspace
   type(NodeObsValidControl) :: valid, validtmp
   type(NodeObsControl) :: obs, error, ytmp
   type(NodeObsControl), dimension(:), allocatable :: ypert
   !
   integer :: mmax, mobs, nrank, inflana, adapmode, inflpoly, nlevnhm, myidx, myidy, myide
   integer :: id, idx, idy, is, ie, js, je, dis, die, djs, dje, ivar, ivarnhm, iter
   real(r_size) :: rhoML, rhoMLold
   real(r_size), dimension(:), allocatable :: X, VTX, Y0, gammab, d, lambda, flambda, gammaa, gammainfl
   real(r_size), dimension(:,:), allocatable :: U, V, Y, YYT
   real(r_size), dimension(:,:,:), allocatable :: Ye
   ! EigenExa
   integer :: eigennx, eigenny
   integer :: eigen_nnod, eigenx_nnod, eigeny_nnod, eigen_inod, eigenx_inod, eigeny_inod 
   integer :: eigen_istart, eigen_iend, eigen_il, eigen_ig, eigen_jstart, eigen_jend, eigen_jl, eigen_jg
   real(8), dimension(:), allocatable :: eigenW
   real(8), dimension(:,:), allocatable :: eigenA, eigenZ
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
   !
   ! MPI Decomposition
   call new(info, 3, nxpe, nype, nepe, myid, nx0, ny0, ne0, 1)
   call initialize_mpi(info)
   call get_myidx(info, myidx)
   call get_myidy(info, myidy)
   call get_myide(info, myide)
   call get_nx(info, nx)
   call get_ny(info, ny)
   call get_ne(info, ne)
   nz = nz0
   nt = nt0
   call get_di(info, dis, die)
   call get_dj(info, djs, dje)
   !
   ! EigenExa
   call eigen_init(MPI_COMM_WORLD)
   call eigen_get_matdims(ne0, eigennx, eigenny)
   allocate(eigenW(ne0), eigenA(eigennx,eigenny), eigenZ(eigennx,eigenny))
   call eigen_get_procs(eigen_nnod, eigenx_nnod, eigeny_nnod)
   call eigen_get_id(eigen_inod, eigenx_inod, eigeny_inod)
   eigen_istart = eigen_loop_start(1, eigenx_nnod, eigenx_inod)
   eigen_iend = eigen_loop_end(ne0, eigenx_nnod, eigenx_inod)
   eigen_jstart = eigen_loop_start(1, eigeny_nnod, eigeny_inod)
   eigen_jend = eigen_loop_end(ne0, eigeny_nnod, eigeny_inod)
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
   call setGridInf(nx0, ny0, nz0, proj, resolution, slon, xi, xj, xlat, xlon, slat, slat2)
   call initialize_rttov(myid)
   call new_convrad(global_obsspace, nx0, ny0, nt0)
   call read_obs(global_obsspace, myid)
   call new_convrad(obsspace, nx, ny, nt0)
   call scatter_obs(global_obsspace, info, obsspace)
   call MPI_BARRIER(MPI_COMM_WORLD, ierror)
   call destroy(global_obsspace)
   ! Observation and observational errors
   call new(obs, obsspace)
   call new(error, obsspace)
   call set_obs(obs, obsspace)
   call set_error(error, obsspace)
   ! Observation validity
   call new(valid, obsspace)
   call new(validtmp, obsspace)
   ! Pre-observation operator
   call new_obs_NodeControl(xobs, info, nz0, nt0)
   call new_NodePreH(preH, .False., vdz, zrp, vctransp, dvtransp, dvtransw, zs, nx, ny, nz0, nt0)
   ! Observation operator
   call new(h, info, obsspace)
   ! Observation perturbation
   allocate(ypert(ne))
   do ie = 1, ne
      call new(ypert(ie), obsspace)
   end do
   ! Temporary variables
   call new(ytmp, obsspace)
   !
   !
   !
   ! Background
   call new_nhm_NodeControl(xtmp, info, .False., nz0, nt0)
   if (iflg_outana == 1) call new_nhm_NodeControl(xbck, info, .False., nz0, nt0)
   allocate(xpert(ne))
   do ie = 1, ne
      call new_nhm_NodeControl(xpert(ie), info, .False., nz0, nt0)
   end do
   call read_ensemble(xpert, info, 'fg', ne)
   ! Boundary conditions
   call new(bc, xtmp)
   ! Control
   call new(rho, control_mode, info, .False., nz0, 1, 1)
   !
   ! Y
   call set_const(xtmp, 0.d0)
   call set_field(valid, 1)
   call set_field(ytmp, 0.d0)
   do ie = 1, ne
      call apply_PreH(preH, control_mode, xpert(ie), xtmp, xobs)
      call apply_H(h, info, xobs, obsspace, validtmp, ypert(ie))
      call and(valid, validtmp)
      call add_obs(ytmp, ypert(ie))
   end do
   call allreduce_ensvalid(valid)
   call allreduce_ensobs(ytmp)
   call divide_obs(ytmp, 1.d0*ne0)
   do ie = 1, ne
      call subtract_obs(ypert(ie), ytmp)
      call divide_obs(ypert(ie), sqrt(1.d0*(ne0-1)))
      call divide_obs(ypert(ie), error)
   end do
   ! Quality control
   if (iflg_incremental == 1) call add_obs(obs, ytmp)
   call subtract_obs(obs, ytmp)
   call divide_obs(obs, error)
   call qccheck(obs, qcthreshold, validtmp)
   call and(valid, validtmp)
   !
   ! X
   call set_const(xtmp, 0.d0)
   do ie = 1, ne
      call add(xtmp, xpert(ie))
   end do
   call allreduce_ens(xtmp)
   call divide(xtmp, 1.d0*ne0)
   do ie = 1, ne
      call subtract(xpert(ie), xtmp)
      call divide(xpert(ie), sqrt(1.d0*(ne0-1)))
      call apply_bc(bc, xpert(ie))
   end do
   !
   !
   !
   ! Pre ETKF
   inflana = inflmode/100
   adapmode = mod(inflmode,100)/10
   !inflpoly = mod(inflmode,10)
   ! Allocate
   ! Numbers of obs per column
   mmax = 0
   do i = 1, nx-1
      do j = 1, ny-1
         call get_mobs(obsspace, i, j, mobs)
         if (mobs > 0) mmax = mmax + mobs
      end do
   end do
   allocate(X(ne), Y(mmax,ne), Y0(mmax))
   allocate(Ye(mmax,ne,nepe))
   allocate(YYT(ne0,ne0), gammab(ne0), U(ne0,ne0), V(ne,ne0))
   allocate(d(ne0), lambda(ne0), flambda(ne0), gammaa(ne0), gammainfl(ne0), VTX(ne0))
   !
   !
   !
   ! ETKF
   do ie = 1, ne
      call extract_y(ypert(ie), obsspace, valid, mobs, Y(:,ie), mmax)
   end do
   call extract_y(obs, obsspace, valid, mobs, Y0, mmax)
   !if (mobs == 0) cycle
   !
   ! Covariance
   call MPI_BARRIER(MPI_COMM_EWORLD, ierror)
   call MPI_GATHER(Y(1:mobs,:), mobs*ne, r_type, Ye(1:mobs,:,:), mobs*ne, r_type, 0, MPI_COMM_EWORLD, ierror)
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
      call allreduce2D('xy', YYT, ne0, ne0)
   end if
   call broadcast2D('e', YYT(:,:), 0, ne0, ne0)
   !
   ! EigenExa
   do eigen_jl = eigen_jstart, eigen_jend
      eigen_jg = eigen_translate_l2g(eigen_jl, eigeny_nnod, eigeny_inod)
      do eigen_il = eigen_istart, eigen_iend
         eigen_ig = eigen_translate_l2g(eigen_il, eigenx_nnod, eigenx_inod)
         eigenA(eigen_il,eigen_jl) = YYT(eigen_ig,eigen_jg)
      end do
   end do
   call eigen_s(ne0, ne0, eigenA, eigennx, eigenW, eigenZ, eigennx, m_forward=32, m_backward=128)
   YYT(:,:) = 0.
   do eigen_jl = eigen_jstart, eigen_jend
      eigen_jg = eigen_translate_l2g(eigen_jl, eigeny_nnod, eigeny_inod)
      do eigen_il = eigen_istart, eigen_iend
         eigen_ig = eigen_translate_l2g(eigen_il, eigenx_nnod, eigenx_inod)
         YYT(eigen_ig,eigen_jg) = eigenZ(eigen_il,eigen_jl)
      end do
   end do
   call allreduce2D('all', YYT, ne0, ne0) ! use U as V temporarily
   !
   ! Find rank
   do ie = 1, ne0
      gammab(ie) = eigenW(ne0-ie+1) - 1.d0
      U(:,ie) = YYT(:,ne0-ie+1)
   end do
   do ie = 1, ne0
      if (gammab(ie) < 1.e-12) exit
   end do
   nrank = ie - 1
   !if (nrank == 0) cycle
   gammab(1:nrank) = sqrt(gammab(1:nrank))
   ! V
   is = myide*ne + 1; ie = is + ne - 1
   V(:,1:nrank) = U(is:ie,1:nrank)
   ! Project Y0 on U
   do ie = 1, ne
      X(ie) = sum(Y(1:mobs,ie)*Y0(1:mobs))
   end do
   call allreduce1D('xy', X, ne)
   do ie = 1, nrank
      d(ie) = sum(V(:,ie)*X(:))/gammab(ie)
   end do
   call allreduce1D('e', d(1:nrank), nrank)
   !
   ! Inflation functions
   lambda(1:nrank) = 1.d0/sqrt(1.d0+gammab(1:nrank)**2)
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
   ! EM algorithm: single rho
   else if (adapmode == 3) then
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
   end if
   if (inflana == 1) then
      lambda(1:nrank) = gammab(1:nrank)*flambda(1:nrank)**2*d(1:nrank)
   else
      lambda(1:nrank) = gammab(1:nrank)*lambda(1:nrank)**2*d(1:nrank)
   end if
   flambda(1:nrank) = (1.d0-flambda(1:nrank))
   !
   ! Update X
   ivar = 1
   if (iflg_outana == 1) xbck = xtmp
   call set_const(xtmp, 0.d0)
   do i = 1+dis, nx-die
   do j = 1+djs, ny-dje
      do iter = 1, rho%nvarnhm(ivar)
         ivarnhm = rho%ivarnhm(ivar,iter)
         nlevnhm = rho%nlevnhm(ivar,iter)
         if (nlevnhm > 1) nlevnhm = nlevnhm - 1
         do k = 1, nlevnhm
            ! Project X on V
            do ie = 1, ne
               X(ie) = xpert(ie)%control(ivarnhm)%p%field(i,j,k,nt)
            end do
            do ie = 1, nrank
               VTX(ie) = sum(V(:,ie)*X(:))
            end do
            call allreduce1D('e', VTX(1:nrank), nrank)
            ! Analysis
            xtmp%control(ivarnhm)%p%field(i,j,k,nt) = sum(VTX(1:nrank)*lambda(1:nrank))
            ! Analysis perturbations
            VTX(1:nrank) = flambda(1:nrank)*VTX(1:nrank)
            X(1:ne) = 0.d0
            do ie = 1, nrank
               X(1:ne) = X(1:ne) + VTX(ie)*V(:,ie)
            end do
            do ie = 1, ne
               xpert(ie)%control(ivarnhm)%p%field(i,j,k,nt) = xpert(ie)%control(ivarnhm)%p%field(i,j,k,nt) - X(ie)
            end do
            !
            if (itout == nt) cycle
            ! Project X on V
            do ie = 1, ne
               X(ie) = xpert(ie)%control(ivarnhm)%p%field(i,j,k,1)
            end do
            do ie = 1, nrank
               VTX(ie) = sum(V(:,ie)*X(:))
            end do
            call allreduce1D('e', VTX(1:nrank), nrank)
            ! Analysis
            xtmp%control(ivarnhm)%p%field(i,j,k,1) = sum(VTX(1:nrank)*lambda(1:nrank))
         end do
      end do
   end do
   end do
   !
   !
   !
   ! Output
   call apply_bc(bc, xtmp)
   if (iflg_outana == 0) then
      if (itout == nt) then
         call write_control(xtmp, nt, info, 'inc')
      else
         call write_control(xtmp, 1, info, 'inc')
      end if
   else
      call add(xtmp, xbck)
      if (itout == nt) then
         call write_control(xtmp, nt, info, 'ana')
      else
         call write_control(xtmp, 1, info, 'ana')
      end if
   end if
   do ie = 1, ne
      call multiply(xpert(ie), sqrt(1.d0*(ne0-1)))
      call apply_bc(bc, xpert(ie))
   end do
   call write_ensemble(xpert, nt, imember1, imember2, info, 'pert', ne)
   !
   ! deallocation
   call destroy(xtmp); call destroy(xobs); call destroy(rho)
   if (iflg_outana == 1) call destroy(xbck)
   do ie = 1, ne
      call destroy(xpert(ie))
      call destroy(ypert(ie))
   end do
   deallocate(xpert, ypert)
   call destroy(obsspace) 
   call destroy(valid); call destroy(validtmp)
   call destroy(obs); call destroy(error); call destroy(ytmp)
   call destroy(h); call destroy(preH)
   deallocate(Ye)
   deallocate(X, Y, YYT, VTX, gammab, U, V)
   deallocate(Y0, d, lambda, flambda, gammaa, gammainfl)
   deallocate(vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw, zs)
   call destroy(info)
   !
   deallocate(eigenW, eigenA, eigenZ)
   call eigen_free()
   call finalize_rttov()
   call MPI_FINALIZE(ierror)
   !
   !stop
end program etkf
