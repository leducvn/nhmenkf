program envardpert
   use variable
   use nhmlib
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
   use NodeScorr_class
   use NodeSigma_class
   use NodeMPI
   implicit none
   !
   integer :: nproc, myid, nx, ny, nz, nt, ne, nlev0, nlev, nrestart, ns, myide
   integer :: ierror, ie, je, ke, iflag, iter, jter, kter, irestart, istep, ide, k1, k2
   integer :: i0, j0, k0, it0
   real(r_size), dimension(:), allocatable :: vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw
   real(r_size), dimension(:,:), allocatable :: zs
   real(r_sngl), dimension(:,:), allocatable :: zstmp, landsea, lon, lat
   real(r_size), dimension(:,:,:), allocatable :: zs0
   !
   type(NodeBCControl) :: bc
   type(NodePreH) :: preH
   type(NodeHControl) :: h
   type(NodeScorr) :: Scorx, Scory, Scorz, Slocx, Slocy, Slocz
   type(NodeSigma) :: Sigma
   !
   type(NodeInfo) :: info
   type(NodeControl) :: x, xbck, xobs, xobsbck, w0, std
   type(NodeControl), dimension(:), allocatable :: xpert, xapert, w
   type(NodeObsControl) :: Az
   type(NodeObsControl), dimension(:), allocatable :: r, z, b
   type(NodeObsControl), dimension(:,:), allocatable :: q
   type(NodeObsSpaceControl) :: obsspace, global_obsspace
   type(NodeObsValidControl) :: valid, validtmp
   type(NodeObsControl) :: obs, error, y, ybck
   !
   character(8) :: pertfile = 'pert0000'
   real(r_size), parameter :: tolerance = 1.e-2
   real(r_size) :: rnorm0, ratio, rnorm, rmat1, rmat2, rmat3
   real(r_size), dimension(:), allocatable :: rmata, rmatb, rmatc, rmatd
   real(r_size), dimension(:,:), allocatable :: beta
   real(r_size), dimension(:,:,:), allocatable :: p, s
   real(r_size), dimension(:,:,:,:), allocatable :: hmat, rmat, g
   real(r_size), dimension(:,:,:,:), allocatable :: ftmp0, ftmp
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
   nrestart = 1
   ns = imember2 - imember1 + 1
   weight_clim = sqrt(weight_clim)
   weight_ens = sqrt(weight_ens)
   !
   ! MPI Decomposition
   call new(info, 3, nxpe, nype, nepe, myid, nx0, ny0, ne0, 1)
   call initialize_mpi(info)
   call get_myide(info, myide)
   call get_nx(info, nx)
   call get_ny(info, ny)
   call get_ne(info, ne)
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
   call new_obs_NodeControl(xobsbck, info, nz0, nt0)
   call new_NodePreH(preH, .True., vdz, zrp, vctransp, dvtransp, dvtransw, zs, nx, ny, nz0, nt0)
   ! Observation operator
   call new(h, info, obsspace)
   ! Temporary variables
   call new(y, obsspace)
   call new(ybck, obsspace)
   !
   !
   !
   ! Background
   call new_nhm_NodeControl(x, info, .False., nz0, nt0)
   call new_nhm_NodeControl(xbck, info, .False., nz0, nt0)
   allocate(xpert(ne), xapert(ne))
   do ie = 1, ne
      call new_nhm_NodeControl(xpert(ie), info, .False., nz0, nt0)
      call new_nhm_NodeControl(xapert(ie), info, .False., nz0, 1)
   end do
   call read_ensemble(xpert, info, 'fg', ne)
   ! Boundary conditions
   call new(bc, x)
   ! Control
   allocate(w(ne))
   do ie = 1, ne
      call new(w(ie), control_mode, info, .False., nz0, nz0, 1)
   end do
   call new(Az, obsspace)
   allocate(z(imember1:imember2), r(imember1:imember2), b(imember1:imember2), q(niteration+1,imember1:imember2))
   do ke = imember1, imember2
      call new(z(ke), obsspace)
      call new(r(ke), obsspace)
      call new(b(ke), obsspace)
   end do
   do iter = 1, niteration+1
      do ke = imember1, imember2
         call new(q(iter,ke), obsspace)
      end do
   end do
   allocate(rmata(imember1:imember2), rmatb(imember1:imember2), rmatc(imember1:imember2), rmatd(imember1:imember2))
   allocate(beta(imember1:imember2,imember1:imember2))
   allocate(p(niteration+1,imember1:imember2,imember1:imember2), s(niteration,imember1:imember2,imember1:imember2))
   allocate(hmat(niteration+1,niteration,imember1:imember2,imember1:imember2))
   allocate(rmat(niteration,niteration,imember1:imember2,imember1:imember2))
   allocate(g(niteration,imember1:2*imember2-imember1+1,imember1:imember2,2))
   !
   !
   !
   ! Climatological
   if (weight_clim > 0.) then
      call new(w0, control_mode, info, .False., nz0, nz0, 1)
      call new(std, control_mode, info, .False., nz0, nz0, 1)
      call get_nlev(w0, nlev0)
      allocate(ftmp0(nx,ny,nlev0,1))
      call new(Scorx, nx0, nlev0, 'x')
      call set_Scorr(Scorx, info, 'cor')
      call new(Scory, ny0, nlev0, 'y')
      call set_Scorr(Scory, info, 'cor')
      call new(Scorz, nlev0, 1, 'z')
      call set_Scorr(Scorz, info, 'cor')
      call new(Sigma, nx, ny, nlev0)
      call set_Sigma1D(Sigma, info)
      !call set_Sigma3D(Sigma, info)
      call get_sigma(Sigma, ftmp0(:,:,:,1))
      call set_field(std, ftmp0, nlev0)
   end if
   !
   ! Localization
   call get_nlev(w(1), nlev)
   allocate(ftmp(nx,ny,nlev,1))
   call new(Slocx, nx0, nlev, 'x')
   call set_Scorr(Slocx, info, 'loc')
   call new(Slocy, ny0, nlev, 'y')
   call set_Scorr(Slocy, info, 'loc')
   call new(Slocz, nlev, 1, 'z')
   call set_Scorr(Slocz, info, 'loc')
   !
   !
   !
   ! Forecast perturbations
   call set_const(x, 0.d0)
   do ie = 1, ne
      call add(x, xpert(ie))
   end do
   call allreduce_ens(x)
   call divide(x, 1.d0*ne0)
   do ie = 1, ne
      call subtract(xpert(ie), x)
      call divide(xpert(ie), sqrt(1.d0*(ne0-1)))
      call apply_bc(bc, xpert(ie))
   end do
   ! Background
   if (iflg_usectl == 1) then
      call read_control(xbck, info, 'ctl')
   else
      xbck = x
   end if
   call set_const(x, 0.d0)
   call set_field(valid, 1)
   call apply_PreH(preH, control_mode, xbck, x, xobsbck)
   call apply_H(h, info, xobsbck, obsspace, validtmp, ybck)
   call and(valid, validtmp)
   if (iflg_incremental == 1) call add_obs(obs, ybck)
   ! Quality control
   y = obs
   call subtract_obs(y, ybck)
   call divide_obs(y, error)
   call qccheck(y, qcthreshold, validtmp)
   call and(valid, validtmp)
   !
   !
   !
   ! Pre EnVAR
   do ie = 1, ne
      call apply_DpreH(preH, control_mode, xpert(ie), xobs)
      call update(xobs, info)
      call set_index_nohalo(ne0, nepe, myide, k1, k2)
      ke = ie+k1-1
      call apply_DH(h, info, xobsbck, xobs, obsspace, validtmp, b(ke))
      call divide_obs(b(ke), error)
      call multiply_obs(b(ke), -inflfactor1) ! fixed quadratic inflation
   end do
   do ke = imember1, imember2
      call set_field(z(ke), 0.d0)
      do ide = 0, nepe-1
         call set_index_nohalo(ne0, nepe, ide, k1, k2)
         if (k1 <= ke .and. ke <= k2) exit
      end do
      call broadcast_ensobs(b(ke), ide)
   end do
   !
   ! EnVAR: Block GMRES Az = b
   istep = 0
   do irestart = 1, nrestart
      ! Calculate r = b - Az
      do ke = imember1, imember2
	 r(ke) = b(ke)
	 if (irestart > 1) then
	    ! Calculate Az
	    Az = z(ke)
	    call divide_obs(Az, error)
	    call set_const(xobs, 0.d0)
	    call apply_DHT(h, info, xobsbck, obsspace, valid, Az, xobs)
	    call apply_DPreHT(preH, control_mode, xobs, x)
	    call apply_bcT(bc, x)
	    !
	    ! C*w
	    if (weight_clim > 0.) then
	       call apply_USITclim(w0, std, x)
	       call get_field(w0, ftmp0, nlev0)
	       call apply_Scorr(Scorz, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scory, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scorx, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scorx, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scory, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scorz, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call set_field(w0, ftmp0, nlev0)
	    end if
	    call apply_USITens(w, ne, x, xpert)
	    do ie = 1, ne
	       call get_field(w(ie), ftmp, nlev)
	       call apply_Scorr(Slocz, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocz, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call set_field(w(ie), ftmp, nlev)
	    end do
	    !
	    call apply_USIens(x, ne, w, xpert)
	    call allreduce_ens(x)
	    if (weight_clim > 0.) call apply_USIclim(x, std, w0)
	    call apply_bc(bc, x)
	    call apply_DpreH(preH, control_mode, x, xobs)
	    call update(xobs, info)
	    call apply_DH(h, info, xobsbck, xobs, obsspace, validtmp, Az)
	    !call and(validtmp, valid)
	    call divide_obs(Az, error)
	    call add_obs(Az, z(ke))
	    call subtract_obs(r(ke), Az)
	 end if
      end do
      !
      ! First norm: QR decomposition
      iter = 1
      beta(:,:) = 0.d0
      do ke = imember1, imember2
         q(iter,ke) = r(ke)
      end do
      do ke = imember1, imember2
         call innerproduct(q(iter,ke), q(iter,ke), obsspace, valid, beta(ke,ke))
         call allreduce0D('xy', beta(ke,ke))
         beta(ke,ke) = sqrt(beta(ke,ke))
         call divide_obs(q(iter,ke), beta(ke,ke))
         do je = ke+1,imember2
            call innerproduct(q(iter,ke), q(iter,je), obsspace, valid, beta(ke,je))
            call allreduce0D('xy', beta(ke,je))
            Az = q(iter,ke)
            call multiply_obs(Az, beta(ke,je))
            call subtract_obs(q(iter,je), Az)
         end do
      end do
      p(:,:,:) = 0.d0
      p(1,:,:) = beta
      rnorm0 = sqrt(sum(beta**2))
      if (myid == 0) print*, 'iter=0', rnorm0
      !
      !
      !
      ! Block GMRES niteration steps
      do iter = 1, niteration
         istep = istep + 1
	 ! Calculate new r
	 do ke = imember1, imember2
	    r(ke) = q(iter,ke)
	    call divide_obs(r(ke), error)
	    call set_const(xobs, 0.d0)
	    call apply_DHT(h, info, xobsbck, obsspace, valid, r(ke), xobs)
	    call apply_DPreHT(preH, control_mode, xobs, x)
	    call apply_bcT(bc, x)
	    !
	    ! C*w
	    if (weight_clim > 0.) then
	       call apply_USITclim(w0, std, x)
	       call get_field(w0, ftmp0, nlev0)
	       call apply_Scorr(Scorz, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scory, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scorx, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scorx, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scory, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call apply_Scorr(Scorz, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	       call set_field(w0, ftmp0, nlev0)
	    end if
	    call apply_USITens(w, ne, x, xpert)
	    do ie = 1, ne
	       call get_field(w(ie), ftmp, nlev)
	       call apply_Scorr(Slocz, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call apply_Scorr(Slocz, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	       call set_field(w(ie), ftmp, nlev)
	    end do
	    !
	    call apply_USIens(x, ne, w, xpert)
	    call allreduce_ens(x)
	    if (weight_clim > 0.) call apply_USIclim(x, std, w0)
	    call apply_bc(bc, x)
	    call apply_DpreH(preH, control_mode, x, xobs)
	    call update(xobs, info)
	    call apply_DH(h, info, xobsbck, xobs, obsspace, validtmp, r(ke))
	    !call and(validtmp, valid)
	    call divide_obs(r(ke), error)
	    call add_obs(r(ke), q(iter,ke))
	 end do
	 !
	 ! Block Arnoldi: Gram-Schmidt
	 do jter = 1, iter
	    do ke = imember1, imember2
	       do je = imember1, imember2
	          call innerproduct(q(jter,je), r(ke), obsspace, valid, hmat(jter,iter,je,ke))
	          call allreduce0D('xy', hmat(jter,iter,je,ke))
	       end do
	       do je = imember1, imember2
	          Az = q(jter,je)
	          call multiply_obs(Az, hmat(jter,iter,je,ke))
	          call subtract_obs(r(ke), Az)
	       end do
	    end do
	 end do
	 !
	 ! Calculate new norm: QR decomposition
         do ke = imember1, imember2
            q(iter+1,ke) = r(ke)
         end do
         do ke = imember1, imember2
            call innerproduct(q(iter+1,ke), q(iter+1,ke), obsspace, valid, hmat(iter+1,iter,ke,ke))
            call allreduce0D('xy', hmat(iter+1,iter,ke,ke))
            hmat(iter+1,iter,ke,ke) = sqrt(hmat(iter+1,iter,ke,ke))
            call divide_obs(q(iter+1,ke), hmat(iter+1,iter,ke,ke))
            do je = ke+1,imember2
               call innerproduct(q(iter+1,ke), q(iter+1,je), obsspace, valid, hmat(iter+1,iter,ke,je))
               call allreduce0D('xy', hmat(iter+1,iter,ke,je))
               Az = q(iter+1,ke)
               call multiply_obs(Az, hmat(iter+1,iter,ke,je))
               call subtract_obs(q(iter+1,je), Az)
            end do
         end do
	 !
	 if (myid == 0) then
	    ! apply G to the new block at column iter
	    rmat(1:iter,iter,:,:) = hmat(1:iter,iter,:,:)
	    do jter = 1, iter-1
	       ! g(jter) now consists of s**2 rotation
	       do ke = imember1, imember2
	          do je = ke+ns, ke+1, -1
	             if (je-1 > imember2) then
	                rmata = rmat(jter+1,iter,je-1-imember2,:)
	             else
	                rmata = rmat(jter,iter,je-1,:)
	             end if
	             if (je > imember2) then
	                rmatb = rmat(jter+1,iter,je-imember2,:)
	             else
	                rmatb = rmat(jter,iter,je,:)
	             end if
	             rmatc = g(jter,je,ke,1)*rmata + g(jter,je,ke,2)*rmatb
	             rmatd = -g(jter,je,ke,2)*rmata + g(jter,je,ke,1)*rmatb
	             if (je-1 > imember2) then
	                rmat(jter+1,iter,je-1-imember2,:) = rmatc
	             else
	                rmat(jter,iter,je-1,:) = rmatc
	             end if
	             if (je > imember2) then
	                rmat(jter+1,iter,je-imember2,:) = rmatd
	             else
	                rmat(jter,iter,je,:) = rmatd
	             end if
	          end do
	       end do
	    end do
	    !
	    ! calculate new G and R
	    do ke = imember1, imember2
	       do je = ke+ns, ke+1, -1
	          if (je-1 > imember2) then
	             rmat1 = hmat(iter+1,iter,je-1-imember2,ke)
	          else
	             rmat1 = rmat(iter,iter,je-1,ke)
	          end if
	          if (je > imember2) then
	             rmat2 = hmat(iter+1,iter,je-imember2,ke)
	          else
	             rmat2 = rmat(iter,iter,je,ke)
	          end if
	          rmat3 = sqrt(rmat1**2+rmat2**2)
	          g(iter,je,ke,1) = rmat1/rmat3
	          g(iter,je,ke,2) = rmat2/rmat3
	          !
	          ! Now apply this new rotation
	          if (je-1 > imember2) then
	             rmata = hmat(iter+1,iter,je-1-imember2,:)
	          else
	             rmata = rmat(iter,iter,je-1,:)
	          end if
	          if (je > imember2) then
	             rmatb = hmat(iter+1,iter,je-imember2,:)
	          else
	             rmatb = rmat(iter,iter,je,:)
	          end if
	          rmatc = g(iter,je,ke,1)*rmata + g(iter,je,ke,2)*rmatb
	          rmatd = -g(iter,je,ke,2)*rmata + g(iter,je,ke,1)*rmatb
	          if (je-1 > imember2) then
	             hmat(iter+1,iter,je-1-imember2,:) = rmatc
	             hmat(iter+1,iter,je-1-imember2,ke) = rmat3
	          else
	             rmat(iter,iter,je-1,:) = rmatc
	             rmat(iter,iter,je-1,ke) = rmat3
	          end if
	          if (je > imember2) then
	             hmat(iter+1,iter,je-imember2,:) = rmatd
	             hmat(iter+1,iter,je-imember2,ke) = 0.d0
	          else
	             rmat(iter,iter,je,:) = rmatd
	             rmat(iter,iter,je,ke) = 0.d0
	          end if
	       end do
	    end do
	    !
	    ! Now apply the new rotation to p
	    do ke = imember1, imember2
	       do je = ke+ns, ke+1, -1
	          if (je-1 > imember2) then
	             rmata = p(iter+1,je-1-imember2,:)
	          else
	             rmata = p(iter,je-1,:)
	          end if
	          if (je > imember2) then
	             rmatb = p(iter+1,je-imember2,:)
	          else
	             rmatb = p(iter,je,:)
	          end if
	          rmatc = g(iter,je,ke,1)*rmata + g(iter,je,ke,2)*rmatb
	          rmatd = -g(iter,je,ke,2)*rmata + g(iter,je,ke,1)*rmatb
	          if (je-1 > imember2) then
	             p(iter+1,je-1-imember2,:) = rmatc
	          else
	             p(iter,je-1,:) = rmatc
	          end if
	          if (je > imember2) then
	             p(iter+1,je-imember2,:) = rmatd
	          else
	             p(iter,je,:) = rmatd
	          end if
	       end do
	    end do
	 end if
	 !
	 ! Check convergence
	 if (myid == 0) then
	    rnorm = sqrt(sum(p(iter+1,:,:)**2))
	    ratio = rnorm/rnorm0
	    print*, 'iter=', iter, rnorm, ratio
	    if (ratio < tolerance .or. rnorm < tolerance) then
	       iflag = 0
	    else
	       iflag = -1
	    end if
	    ! Monitoring
	    open(90, file="iteration.log", position="append")
	    if (irestart == 1 .and. iter == 1) then
	       write(90,'(a6,2a14)') "istep", "residual", "ratio"
	    end if
	    write(90,'(i4,2e14.7)') istep, rnorm, ratio
	    close(90)
	 end if
	 call int_broadcast0D('all', iflag, 0)
	 if (iflag == 0) then
	    if (myid == 0) print*, 'Linear equation system has reached its tolerance.', rnorm, rnorm0
	    exit
	 end if
	 if (iter == niteration) exit
      end do
      !
      !
      !
      ! solve Rs = p
      if (myid == 0) then
	 do jter = iter, 1, -1
	    s(jter,:,:) = p(jter,:,:)
	    do kter = iter, jter+1, -1
	       do ke = imember1, imember2
	          do je = imember1, imember2
	             beta(je,ke) = sum(rmat(jter,kter,je,:)*s(kter,:,ke))
	          end do
	       end do
	       s(jter,:,:) = s(jter,:,:) - beta
	    end do
	    beta = s(jter,:,:)
	    do je = imember2, imember1, -1
	       s(jter,je,:) = beta(je,:)
	       do ke = imember2, je+1, -1
	          s(jter,je,:) = s(jter,je,:) - s(jter,ke,:)*rmat(jter,jter,je,ke)
	       end do
	       s(jter,je,:) = s(jter,je,:)/rmat(jter,jter,je,je)
	    end do
	 end do
      end if
      call broadcast3D('all', s(1:iter,:,:), 0, iter, ns, ns)
      !
      ! Calculate new z
      do ke = imember1, imember2
         do jter = 1, iter
            do je = imember1, imember2
	       Az = q(jter,je)
	       call multiply_obs(Az, s(jter,je,ke))
	       call add_obs(z(ke), Az)
	    end do
         end do
      end do
      !
      if (iflag == 0) exit
      if (irestart == nrestart) then
	 iflag = 1
	 if (myid == 0) print*, 'The iteration number exceeds.', rnorm, rnorm0
	 exit
      end if
   end do
   !
   !
   !
   ! Transform back to the weight space and output
   do ke = imember1, imember2
      call divide_obs(z(ke), error)
      call set_const(xobs, 0.d0)
      call apply_DHT(h, info, xobsbck, obsspace, valid, z(ke), xobs)
      call apply_DPreHT(preH, control_mode, xobs, x)
      call apply_bcT(bc, x)
      !
      ! C*w
      if (weight_clim > 0.) then
	 call apply_USITclim(w0, std, x)
	 call get_field(w0, ftmp0, nlev0)
	 call apply_Scorr(Scorz, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	 call apply_Scorr(Scory, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	 call apply_Scorr(Scorx, .True., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	 call apply_Scorr(Scorx, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	 call apply_Scorr(Scory, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	 call apply_Scorr(Scorz, .False., info, ftmp0(:,:,:,1), nx, ny, nlev0)
	 call set_field(w0, ftmp0, nlev0)
      end if
      call apply_USITens(w, ne, x, xpert)
      do ie = 1, ne
	 call get_field(w(ie), ftmp, nlev)
	 call apply_Scorr(Slocz, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	 call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	 call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	 call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	 call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	 call apply_Scorr(Slocz, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	 call set_field(w(ie), ftmp, nlev)
      end do
      !
      call apply_USIens(x, ne, w, xpert)
      call allreduce_ens(x)
      if (weight_clim > 0.) call apply_USIclim(x, std, w0)
      do ide = 0, nepe-1
         call set_index_nohalo(ne0, nepe, ide, k1, k2)
         if (k1 <= ke .and. ke <= k2) exit
      end do
      if (myide == ide) then
	 call add(x, xpert(ke-k1+1))
         call get_field(x, nt, xapert(ke-k1+1))
         !call multiply(xapert(ke-k1+1), sqrt(1.d0*(ne0-1)))
         !call apply_bc(bc, xapert(ke-k1+1))
      end if
   end do
   !
   do ie = 1, ne
      call multiply(xapert(ie), sqrt(1.d0*(ne0-1)))
      call apply_bc(bc, xapert(ie))
   end do
   call write_ensemble(xapert, 1, imember1, imember2, info, 'pert', ne)
   !do ie = 1, ne
      !call multiply(xpert(ie), sqrt(1.d0*(ne0-1)))
      !call apply_bc(bc, xpert(ie))
   !end do
   !call write_ensemble(xpert, nt, imember1, imember2, info, 'pertf', ne)
   !
   !
   !
   ! deallocation
   call destroy(x); call destroy(xbck); call destroy(xobs); call destroy(xobsbck)
   do ie = 1, ne
      call destroy(xpert(ie))
      call destroy(xapert(ie))
      call destroy(w(ie))
   end do
   call destroy(Az)
   do ke = imember1, imember2
      call destroy(z(ke)); call destroy(r(ke)); call destroy(b(ke))
   end do
   do iter = 1, niteration
      do ke = imember1, imember2
         call destroy(q(iter,ke))
      end do
   end do
   deallocate(xpert, xapert, w, z, r, b, q)
   call destroy(obsspace)
   call destroy(valid); call destroy(validtmp)
   call destroy(obs); call destroy(error); call destroy(y); call destroy(ybck)
   call destroy(Slocx); call destroy(Slocy); call destroy(Slocz)
   call destroy(h); call destroy(preH)
   if (weight_clim > 0.) then
      call destroy(w0); call destroy(std); deallocate(ftmp0)
      call destroy(Scorx); call destroy(Scory); call destroy(Scorz); call destroy(Sigma)
   end if
   deallocate(rmata, rmatb, rmatc, rmatd, beta, p, s, hmat, rmat, g)
   deallocate(vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw, zs)
   call destroy(info)
   !
   call finalize_rttov()
   call MPI_FINALIZE(ierror)
end program envardpert

