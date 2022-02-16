program envard1000
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
   integer :: nproc, myid, nx, ny, nz, nt, ne, nlev0, nlev, nrestart, myide
   integer :: ierror, ie, iflag, iter, jter, kter, irestart, istep
   integer :: i0, j0, k0, it0
   real(r_size), dimension(:), allocatable :: vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw
   real(r_size), dimension(:,:), allocatable :: zs
   real(r_sngl), dimension(:,:), allocatable :: zstmp, landsea, lon, lat
   real(r_size), dimension(:,:,:), allocatable :: zs0
   !
   type(NodeBCControl) :: bc
   type(NodePreH) :: preH
   type(NodeHControl) :: h
   type(NodeScorr) :: Scorx, Scory, Scorz, Slocx, Slocy
   type(NodeSigma) :: Sigma
   !
   type(NodeInfo) :: info
   type(NodeControl) :: x, xbck, xobs, xobsbck, w0, std
   type(NodeControl), dimension(:), allocatable :: xpert, w
   type(NodeObsControl) :: r, z, b, Az
   type(NodeObsControl), dimension(:), allocatable :: q
   type(NodeObsSpaceControl) :: obsspace, global_obsspace
   type(NodeObsValidControl) :: valid, validtmp
   type(NodeObsControl) :: obs, error, y, ybck
   !
   real(r_size), parameter :: tolerance = 1.e-5
   real(r_size) :: rnorm0, ratio, norm, rmat1, rmat2
   real(r_size), dimension(:), allocatable :: rnorm, p, s
   real(r_size), dimension(:,:), allocatable :: hmat, rmat, g
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
   call new_convgnss(obsspace, nx, ny, nt0)
   call new_convgnss(global_obsspace, nx0, ny0, nt0)
   call read_obs(global_obsspace, myid)
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
   allocate(xpert(ne))
   do ie = 1, ne
      call new_nhm_NodeControl(xpert(ie), info, .False., nz0, nt0)
   end do
   call read_ensemble(xpert, info, 'fg', ne)
   ! Boundary conditions
   call new(bc, x)
   ! Control
   allocate(w(ne))
   do ie = 1, ne
      call new(w(ie), control_mode, info, .False., nz0, 1, 1)
   end do
   call new(z, obsspace)
   call new(r, obsspace)
   call new(b, obsspace)
   call new(Az, obsspace)
   allocate(q(niteration))
   do iter = 1, niteration
      call new(q(iter), obsspace)
   end do
   allocate(rnorm(0:niteration), p(niteration+1), s(niteration))
   allocate(hmat(niteration+1,niteration), rmat(niteration,niteration), g(niteration,2))
   !
   !
   !
   ! Climatological
   if (weight_clim > 0.) then
      call new(w0, enkf5, info, .False., nz0, nz0, 1)
      call new(std, enkf5, info, .False., nz0, nz0, 1)
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
   !include "checkSbT.h"
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
   !include "checkPreH.h"
   call apply_H(h, info, xobsbck, obsspace, validtmp, ybck)
   call and(valid, validtmp)
   !include "checkDH.h"
   if (iflg_incremental == 1) call add_obs(obs, ybck)
   ! Quality control
   y = obs
   call subtract_obs(y, ybck)
   call divide_obs(y, error)
   call qccheck(y, qcthreshold, validtmp)
   call and(valid, validtmp)
   !include "checkUSI.h"
   !
   !
   !
   ! Pre EnVAR
   call set_field(z, 0.d0)
   ! Calculate b
   b = y
   !
   ! EnVAR: GMRES Az = b
   istep = 0
   do irestart = 1, nrestart
      ! Calculate r = b - Az
      r = b
      if (irestart > 1) then
         ! Calculate Az
         Az = z
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
	    call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
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
         call add_obs(Az, z)
         call subtract_obs(r, Az)
      end if
      ! First rnorm
      call innerproduct(r, r, obsspace, valid, rnorm(0))
      call allreduce0D('xy', rnorm(0))
      rnorm(0) = sqrt(rnorm(0))
      if (irestart == 1) rnorm0 = rnorm(0)
      p(:) = 0.d0
      p(1) = rnorm(0)
      !
      !
      !
      ! GMRES niteration steps
      do iter = 1, niteration
         istep = istep + 1
	 call divide_obs(r, rnorm(iter-1))
	 q(iter) = r
	 !
	 ! Calculate new r
	 r = q(iter)
	 call divide_obs(r, error)
	 call set_const(xobs, 0.d0)
         call apply_DHT(h, info, xobsbck, obsspace, valid, r, xobs)
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
	    call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
	    call set_field(w(ie), ftmp, nlev)
	 end do
	 !
	 call apply_USIens(x, ne, w, xpert)
	 call allreduce_ens(x)
	 if (weight_clim > 0.) call apply_USIclim(x, std, w0)
         call apply_bc(bc, x)
         call apply_DpreH(preH, control_mode, x, xobs)
	 call update(xobs, info)
         call apply_DH(h, info, xobsbck, xobs, obsspace, validtmp, r)
	 !call and(validtmp, valid)
         call divide_obs(r, error)
         call add_obs(r, q(iter))
	 !
	 ! Arnoldi: Gram-Schmidt
	 do jter = 1, iter
	    call innerproduct(q(jter), r, obsspace, valid, hmat(jter,iter))
	    call allreduce0D('xy', hmat(jter,iter))
	    Az = q(jter)
	    call multiply_obs(Az, hmat(jter,iter))
	    call subtract_obs(r, Az)
	 end do
	 !
	 ! Calculate new norm
	 call innerproduct(r, r, obsspace, valid, rnorm(iter))
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
	    if (irestart == 1 .and. iter == 1) then
	       write(90,'(a6,2a14)') "istep", "residual", "ratio"
	    end if
	    write(90,'(i4,2e14.7)') istep, abs(p(iter+1)), ratio
	    close(90)
	 end if
	 call int_broadcast0D('all', iflag, 0)
	 if (iflag == 0) then
	    if (myid == 0) print*, 'Linear equation system has reached its tolerance.', abs(p(iter+1)), rnorm0
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
	    s(jter) = p(jter)
	    do kter = iter, jter+1, -1
	       s(jter) = s(jter) - s(kter)*rmat(jter,kter)
	    end do
	    s(jter) = s(jter)/rmat(jter,jter)
	 end do
      end if
      call broadcast1D('all', s(1:iter), 0, iter)
      !
      ! Calculate new z
      do jter = 1, iter
	 Az = q(jter)
	 call multiply_obs(Az, s(jter))
	 call add_obs(z, Az)
      end do
      !
      if (iflag == 0) exit
      if (irestart == nrestart) then
	 iflag = 1
	 if (myid == 0) print*, 'The iteration number exceeds.', abs(p(iter+1)), rnorm0
	 exit
      end if
   end do
   !include "checkA.h"
   !
   !
   !
   ! Transform back to the weight space
   call divide_obs(z, error)
   call set_const(xobs, 0.d0)
   call apply_DHT(h, info, xobsbck, obsspace, valid, z, xobs)
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
      call apply_Scorr(Slocy, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
      call apply_Scorr(Slocx, .True., info, ftmp(:,:,:,1), nx, ny, nlev)
      call apply_Scorr(Slocx, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
      call apply_Scorr(Slocy, .False., info, ftmp(:,:,:,1), nx, ny, nlev)
      call set_field(w(ie), ftmp, nlev)
   end do
   !
   !
   !
   ! Output
   call apply_USIens(x, ne, w, xpert)
   call allreduce_ens(x)
   if (weight_clim > 0.) call apply_USIclim(x, std, w0)
   call apply_bc(bc, x)
   if (iflg_outana == 1) then
      call add(x, xbck)
      call write_control(x, itout, info, 'ana')
   else
      call write_control(x, itout, info, 'inc')
   end if
   !
   ! deallocation
   call destroy(x); call destroy(xbck); call destroy(xobs); call destroy(xobsbck)
   do ie = 1, ne
      call destroy(xpert(ie))
      call destroy(w(ie))
   end do
   call destroy(z); call destroy(r); call destroy(b); call destroy(Az)
   do iter = 1, niteration
      call destroy(q(iter))
   end do
   deallocate(xpert, w, q)
   call destroy(obsspace)
   call destroy(valid); call destroy(validtmp)
   call destroy(obs); call destroy(error); call destroy(y); call destroy(ybck)
   call destroy(Slocx); call destroy(Slocy)
   call destroy(h); call destroy(preH)
   if (weight_clim > 0.) then
      call destroy(w0); call destroy(std); deallocate(ftmp0)
      call destroy(Scorx); call destroy(Scory); call destroy(Scorz); call destroy(Sigma)
   end if
   deallocate(rnorm, p, s, hmat, rmat, g, ftmp)
   deallocate(vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw, zs)
   call destroy(info)
   !
   call finalize_rttov()
   call MPI_FINALIZE(ierror)
   !
   !stop
end program envard1000
