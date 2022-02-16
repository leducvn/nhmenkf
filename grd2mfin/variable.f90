module variable
   implicit none
   !
   integer, parameter :: ngm = 4
   real, parameter :: rd = 287.05d0, cp = 7.0d0/2.d0*rd, rdvcp = rd/cp, cpdvrd = cp/rd
   real, parameter :: tkelvn = 273.15d0, g0 = 9.80665d0, gvcp = g0/cp
   real, parameter :: gamma = 0.0065, gmrgh = gamma*rd/g0, gmrgi = 1/gmrgh
   real, parameter :: tcw = 0.0d0, e0cw = 6.11d2, tetn1w = 17.27d0, tetn2w = 273.15d0, tetn3w = 35.85d0
   real, parameter :: tci = -36.0d0, e0ci = 611.d0, tetn1i = 21.875d0, tetn2i = 273.15d0, tetn3i = 7.65d0
   real, parameter :: ra = 6371.E+3, ptrf = 300., presrf = 100000.
   !
   integer :: hydrostatic = 0, iyear, imonth, iday, ihour, iminute
   character(len=4) :: proj
   integer :: nx, ny, nz, vctrans_type, n_vctrans
   real :: resolution, xi, xj, xlon, xlat, slon, slat, dz1, dz2, ztop, zl_vctrans, zh_vctrans
   !
   real, dimension(:), allocatable :: vdz, vrdz, vrdz2, zrp, zrw, vctransp, vctransw, dvtransp, dvtransw
   real, dimension(:,:), allocatable :: zs, landsea, lon, lat
   real, dimension(:,:,:), allocatable :: fmap
   real, dimension(:,:,:), allocatable :: tsoil, dtsoil
   real, dimension(:,:,:), allocatable :: u, v, w, dua, dva, dwa, du, dv, dw
   real, dimension(:,:,:), allocatable :: t, pnh, dt, dpnh
   real, dimension(:,:,:), allocatable :: qv, qc, qi, qr, qs, qg, dqv, dqc, dqi, dqr, dqs, dqg
   real, dimension(:,:), allocatable :: pb, pts, ps, dptsdt, sst, ptop, pmsl, dps
   real, dimension(:,:,:), allocatable :: p, ph, pt, rh, h
   !
   namelist /control/ hydrostatic, iyear, imonth, iday, ihour, iminute
   namelist /model/ proj, resolution, nx, ny, xi, xj, xlon, xlat, slon, slat, &
                  & nz, dz1, dz2, ztop, vctrans_type, n_vctrans, zl_vctrans, zh_vctrans
end module variable
