program letkfinfl
   use variable
   use nhmlib
   use interpolate, only : setGridInf
   use NodeInfo_class
   use NodeObsSpaceControl_class
   use NodeObsValidControl_class
   use NodeObsVLocControl_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeBCControl_class
   use NodeHControl_class
   use NodePreH_class
   use Localization_class
   use enkf
   use NodeMPI
   implicit none
   !
   integer :: nproc, myid, nx, ny, nz, nt, ne
   integer :: ierror, ie, ke, iflag, ivar
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
   type(NodeControl) :: x, xbck, xobs, rho, muy, xa, xb
   type(NodeControl), dimension(:), allocatable :: xpert, xapert
   type(NodeObsSpaceControl) :: obsspace, global_obsspace
   type(NodeObsValidControl) :: valid, global_valid, validtmp
   type(NodeObsControl) :: obs, error, y, ybck, ytmp, global_y
   type(NodeObsControl), dimension(:), allocatable :: ypert, global_ypert
   type(NodeObsControl) :: xyloc, global_xyloc
   type(NodeObsVLocControl) :: logp, corr, vproftmp, vloc, global_vloc
   type(NodeObsVLocControl), dimension(:), allocatable :: vprof
   !
   character(8) :: pertfile = 'pert0000'
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
   nepe = 1
   hscale = 1000.*hscale/resolution
   hscaleqv = 1000.*hscaleqv/resolution
   !
   ! MPI Decomposition
   call new_NodeInfo(info, nxpe, nype, nepe, myid, nx0, ny0, ne0, 1)
   call initialize_mpi(info)
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
   call broadcast1D(vdz, 0, nz)
   call broadcast1D(vrdz, 0, nz)
   call broadcast1D(vrdz2, 0, nz)
   call broadcast1D(zrp, 0, nz)
   call broadcast1D(vctransp, 0, nz)
   call broadcast1D(dvtransp, 0, nz)
   call broadcast1D(zrw, 0, nz)
   call broadcast1D(vctransw, 0, nz)
   call broadcast1D(dvtransw, 0, nz)
   !
   ! Const
   allocate(zs0(nx0,ny0,1), zs(nx,ny))
   if (myid == 0) then
      allocate(zstmp(nx0,ny0), landsea(nx0,ny0), lon(nx0,ny0), lat(nx0,ny0))
      call read_nhmcst('mfhm', nx0, ny0, lon, lat, zstmp, landsea)
      zs0(:,:,1) = zstmp
      deallocate(zstmp, landsea, lon, lat)
   end if
   call scatter(info, zs0(:,:,1), 0, nx0, ny0, 1)
   zs(1:nx,1:ny) = zs0(1:nx,1:ny,1)
   deallocate(zs0)
   !
   !
   !
   ! Observation space
   call setGridInf(nx0, ny0, nz0, proj, resolution, slon, xi, xj, xlat, xlon, slat, slat2)
   call new_conv(obsspace, nx, ny, nt0)
   call new_conv(global_obsspace, nx0, ny0, nt0)
   call read_obs(global_obsspace, info)
   call scatter_obs(global_obsspace, info, obsspace)
   ! Observation and observational errors
   call new(obs, obsspace)
   call new(error, obsspace)
   call set_obs(obs, obsspace)
   call set_error(error, obsspace)
   ! Observation validity
   call new(valid, obsspace)
   call new(validtmp, obsspace)
   call new(global_valid, global_obsspace)
   ! Pre-observation operator
   call new_obs_NodeControl(xobs, info, nz0, nt0)
   call new_NodePreH(preH, vdz, zrp, vctransp, dvtransp, dvtransw, zs, nx, ny, nz0, nt0)
   ! Observation operator
   call new(h, info, obsspace)
   call new(xyloc, obsspace, 2)
   call new(global_xyloc, global_obsspace, 2)
   call get_xyloc(h, xyloc)
   call gather_obs(xyloc, info, obsspace, global_obsspace, global_xyloc)
   call destroy(xyloc)
   ! Observation perturbation
   allocate(ypert(ne))
   allocate(global_ypert(ne))
   do ie = 1, ne
      call new(ypert(ie), obsspace)
      call new(global_ypert(ie), global_obsspace)
   end do
   ! Temporary variables
   call new(y, obsspace)
   call new(ytmp, obsspace)
   call new(ybck, obsspace)
   call new(global_y, global_obsspace)
   !
   !
   !
   ! Background
   call new_nhm_NodeControl(x, info, nz0, nt0)
   call new_nhm_NodeControl(xbck, info, nz0, nt0)
   call new_nhm_NodeControl(xa, info, nz0, 1)
   allocate(xpert(ne))
   allocate(xapert(ne))
   do ie = 1, ne
      call new_nhm_NodeControl(xpert(ie), info, nz0, nt0)
      call new_nhm_NodeControl(xapert(ie), info, nz0, 1)
   end do
   call read_ensemble(xpert, info, 'fg', ne)
   ! Boundary conditions
   call new(bc, x)
   ! Control
   call new(rho, control_mode, info, nz0, nz0, 1)
   call new(muy, control_mode, info, nz0, nz0, 1)
   call new(loc, hscale, hscaleqv, rho, y)
   !
   !
   !
   ! First vertical localization: set vertical perturbations
   call new(vloc, obsspace, nz0, 1)
   call new(global_vloc, global_obsspace, nz0, 1)
   call new(logp, obsspace, nz0, 1)
   call new(vproftmp, obsspace, nz0)
   call new(corr, obsspace, nz0)
   allocate(vprof(ne))
   call set_field(logp, 0.d0)
   call set_field(corr, 0.d0)
   call set_field(ybck, 0.d0)
   call set_field(valid, 1)
   call set_const(x, 0.d0)
   do ie = 1, ne
      call apply_preH(preH, control_mode, xpert(ie), x, xobs)
      ! mean logp: use vloc as temporary variable
      call set_logp(vloc, info, obsspace, xobs)
      call add(logp, vloc)
      ! mean vprof: use corr as ensemble mean
      call new(vprof(ie), obsspace, nz0)
      call set_profile(vprof(ie), info, obsspace, xobs)
      call add(corr, vprof(ie))
      ! mean y: use y as ensemble mean
      call apply_H(h, info, xobs, obsspace, validtmp, ypert(ie))
      call and(valid, validtmp)
      call add(ybck, ypert(ie))
   end do
   call divide(logp, 1.d0*ne)
   call divide(corr, 1.d0*ne)
   call divide(ybck, 1.d0*ne)
   do ie = 1, ne
      call subtract(vprof(ie), corr)
      call subtract(ypert(ie), ybck)
   end do
   !
   ! Second vertical localization: normalize perturbations
   call set_field(corr, 0.d0)
   call set_field(y, 0.d0)
   do ie = 1, ne
      vproftmp = vprof(ie)
      call power(vproftmp, 2.d0)
      call add(corr, vproftmp)
      ytmp = ypert(ie)
      call power(ytmp, 2.d0)
      call add(y, ytmp)
   end do
   call divide(corr, 1.d0*(ne-1))
   call power(corr, 0.5d0)
   call divide(y, 1.d0*(ne-1))
   call power(y, 0.5d0)
   do ie = 1, ne
      call divide(vprof(ie), corr)
      call divide(ypert(ie), y)
   end do
   !
   ! Third vertical localization: correlation
   call set_field(corr, 0.d0)
   do ie = 1, ne
      call multiply(vprof(ie), ypert(ie))
      call add(corr, vprof(ie))
   end do
   call divide(corr, 1.d0*(ne-1))
   call ecorap(corr, vscale, logp, vloc)
   call gather_obs(vloc, info, obsspace, global_obsspace, global_vloc)
   do ie = 1, ne
      call destroy(vprof(ie))
   end do
   call destroy(vproftmp)
   call destroy(corr)
   call destroy(logp)
   call destroy(vloc)
   deallocate(vprof)
   ! restore perturbations in observation space
   do ie = 1, ne
      call multiply(ypert(ie), y)
      call divide(ypert(ie), sqrt(1.d0*(ne-1)))
   end do
   !
   !
   !
   ! Forecast perturbations: X and Y
   call set_const(x, 0.d0)
   do ie = 1, ne
      call add(x, xpert(ie))
   end do
   call divide(x, 1.d0*ne)
   do ie = 1, ne
      call subtract(xpert(ie), x)
      call divide(xpert(ie), sqrt(1.d0*(ne-1)))
      call apply_bc(bc, xpert(ie))
   end do
   ! Background
   if (iflg_usectl == 1) then
      call read_control(xbck, info, 'ctl')
      call set_const(x, 0.d0)
      call apply_preH(preH, control_mode, xbck, x, xobs)
      call apply_H(h, info, xobs, obsspace, validtmp, ybck)
      call and(valid, validtmp)
   else
      xbck = x
   end if
   if (iflg_incremental == 1) call add(obs, ybck)
   !
   ! Normalization
   y = obs
   call subtract(y, ybck)
   call divide(y, error)
   call qccheck(y, qcthreshold, validtmp)
   call and(valid, validtmp)
   call gather_obs(valid, info, obsspace, global_obsspace, global_valid)
   call gather_obs(y, info, obsspace, global_obsspace, global_y)
   do ie = 1, ne
      call divide(ypert(ie), error)
      call gather_obs(ypert(ie), info, obsspace, global_obsspace, global_ypert(ie))
   end do
   !
   !
   !
   ! LETKF
   call set_const(xa, 0.d0)
   do ie = 1, ne
      call get_field(xpert(ie), nt, xapert(ie))
   end do
   call solve_LETKFinfl(loc, info, itout, nt, ne, inflmode, inflfactor1, inflfactor2, global_obsspace, global_valid, global_xyloc, global_vloc, global_y, global_ypert, xpert, rho, muy, xa, xapert)
   call apply_bc(bc, xa)
   if (iflg_outana == 1) then
      call new_nhm_NodeControl(xb, info, nz0, 1)
      call get_field(xbck, itout, xb)
      call add(xa, xb)
      call write_control(xa, 1, info, 'ana')
   else
      call write_control(xa, 1, info, 'inc')
   end if
   call write_control(rho, 1, info, 'rho')
   call write_control(muy, 1, info, 'muy')
   do ie = 1, ne
      call multiply(xapert(ie), sqrt(1.d0*(ne-1)))
      call apply_bc(bc, xapert(ie))
   end do
   call write_ensemble(xapert, 1, imember1, imember2, info, 'pert', ne)
   !
   !
   !
   ! deallocation
   call destroy(x); call destroy(xbck); call destroy(xobs); call destroy(rho); call destroy(muy); call destroy(xa)
   if (iflg_outana == 1) call destroy(xb)
   do ie = 1, ne
      call destroy(xpert(ie)); call destroy(xapert(ie))
   end do
   deallocate(xpert, xapert)
   call destroy(obsspace); call destroy(global_obsspace)
   call destroy(valid); call destroy(validtmp); call destroy(global_valid)
   call destroy(obs); call destroy(error); call destroy(y); call destroy(ybck); call destroy(ytmp); call destroy(global_y)
   do ie = 1, ne
      call destroy(ypert(ie)); call destroy(global_ypert(ie))
   end do
   deallocate(ypert, global_ypert)
   call destroy(global_xyloc); call destroy(global_vloc)
   call destroy(loc); call destroy(h); call destroy(preH)
   deallocate(vdz, vrdz, vrdz2, zrp, vctransp, dvtransp, zrw, vctransw, dvtransw, zs)
   call destroy(info)
   !
   call MPI_FINALIZE(ierror)
   !
   !stop
end program letkfinfl
