module variable
   implicit none
   !
   integer, parameter :: r_size=kind(0.0d0), r_dble=kind(0.0d0), r_sngl=kind(0.0e0)
   integer, parameter :: enkf0 = 0, enkf1 = 1, enkf2 = 2, enkf3 = 3, enkf4 = 4, enkf5 = 5, &
                       & enkf11 = 11, enkf12 = 12, enkf13 = 13, enkf14 = 14, nsoil = 4
   real(kind=r_dble), parameter :: rd = 287.05d0, cp = 7.0d0/2.d0*rd, lqv = 2.51d6, rdvcp = rd/cp, cpdvrd = cp/rd
   real(kind=r_dble), parameter :: tkelvn = 273.15d0, g0 = 9.80665d0, gvrd = g0/rd, gvcp = g0/cp
   real(kind=r_dble), parameter :: gamma = 0.0065, gmrgh = gamma*rd/g0, gmrgi = 1/gmrgh
   real(kind=r_dble), parameter :: tcw = 0.0d0, e0cw = 6.11d2, tetn1w = 17.27d0, tetn2w = 273.15d0, tetn3w = 35.85d0
   real(kind=r_dble), parameter :: tci = -36.0d0, e0ci = 611.d0, tetn1i = 21.875d0, tetn2i = 273.15d0, tetn3i = 7.65d0
   real(kind=r_dble), parameter :: ra = 6371.E+3, ptrf = 300., presrf = 100000., pi = 4.d0*atan(1.d0)
   real(kind=r_sngl), parameter :: slat2 = 60.0
   real(kind=r_dble), parameter :: factr = 1.d+10, epsmch = epsilon(1.d0), pgtol = 1.d-5, eps = 1.0d-10
   !
   character(len=3) :: cnvobs_format, tcobs_format, gnssobs_format, radobs_format
   integer :: control_mode, iflg_perturbed_obs, iflg_incremental, iflg_surface_obs, iflg_skip_psobs, &
            & iflg_usectl, iflg_outana, imember1, imember2, itout
   integer :: nxpe, nype, nepe, niteration, ne0, ngroup, inflmode, nzmode
   character(len=4) :: proj
   integer :: nx0, ny0, nz0, nt0, vctrans_type, n_vctrans
   real(kind=r_sngl) :: resolution, xi, xj, xlon, xlat, slon, slat, dz1, dz2, ztop, zl_vctrans, zh_vctrans
   real(kind=r_dble) :: hscale, hscaleqv, vscale, qcthreshold, inflfactor1, inflfactor2, weight_clim, weight_ens
   real(kind=r_dble) :: ramin = 0.08
   !
   namelist /control_nl/ control_mode, iflg_perturbed_obs, iflg_incremental, iflg_surface_obs, iflg_skip_psobs, iflg_usectl, iflg_outana, &
                       & cnvobs_format, tcobs_format, gnssobs_format, radobs_format, &
                       & imember1, imember2, itout, nxpe, nype, nepe, niteration, ne0, ngroup, inflmode, nzmode, &
                       & hscale, hscaleqv, vscale, qcthreshold, inflfactor1, inflfactor2, weight_clim, weight_ens
   namelist /model_nl/ proj, resolution, nx0, ny0, nz0, nt0, xi, xj, xlon, xlat, slon, slat, &
                     & dz1, dz2, ztop, vctrans_type, n_vctrans, zl_vctrans, zh_vctrans
   !
   !
   !
contains
   !
   !
   !
   subroutine initialize_namelist()
      implicit none
      !
      iflg_perturbed_obs = 0
      iflg_incremental = 0
      iflg_surface_obs = 1
      iflg_skip_psobs = 1
      iflg_usectl = 0
      iflg_outana = 0
      cnvobs_format = 'cda'
      tcobs_format = 'txt'
      gnssobs_format = 'cda'
      radobs_format = 'cda'
      imember1 = 0
      imember2 = 0
      niteration = 50
      ngroup = 0
      itout = 1
      inflmode = 60
      nzmode = 1
      hscale = 550.d0
      hscaleqv = 400.d0
      vscale = 1.45d0
      qcthreshold = 8.d0
      inflfactor1 = 0.9d0
      inflfactor2 = 0.0d0
      weight_clim = 0.d0
      weight_ens = 1.d0
      !
      return
   end subroutine initialize_namelist
end module variable
