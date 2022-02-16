module variable
   implicit none
   !
   integer, parameter :: mfin = 0, nhm = 10, ngm = 4
   real, parameter :: rd = 287.05d0, cp = 7.0d0/2.d0*rd, rdvcp = rd/cp, cpdvrd = cp/rd
   real, parameter :: tkelvn = 273.15d0, g0 = 9.80665d0, gvcp = g0/cp
   real, parameter :: gamma = 0.0065, gmrgh = gamma*rd/g0, gmrgi = 1/gmrgh
   real, parameter :: tcw = 0.0d0, e0cw = 6.11d2, tetn1w = 17.27d0, tetn2w = 273.15d0, tetn3w = 35.85d0
   real, parameter :: tci = -36.0d0, e0ci = 611.d0, tetn1i = 21.875d0, tetn2i = 273.15d0, tetn3i = 7.65d0
   real, parameter :: ra = 6371.E+3, ptrf = 300., presrf = 100000.
   !
   integer :: input_mode, iyear, imonth, iday, ihour, iminute
   character(len=4) :: proj
   integer :: nx, ny, nz, vctrans_type, n_vctrans
   real :: resolution, xi, xj, xlon, xlat, slon, slat, dz1, dz2, ztop, zl_vctrans, zh_vctrans
   character(len=4) :: proj0
   integer :: nx0, ny0, nz0, vctrans_type0, n_vctrans0
   real :: resolution0, xi0, xj0, xlon0, xlat0, slon0, slat0, dz10, dz20, ztop0, zl_vctrans0, zh_vctrans0
   !
   ! inner
   real, dimension(:), allocatable :: vdz, vrdz, vrdz2, zrp, zrw, vctransp, vctransw, dvtransp, dvtransw
   real, dimension(:,:), allocatable :: zs, landsea, lon, lat
   real, dimension(:,:,:), allocatable :: fmap
   ! outer
   real, dimension(:), allocatable :: vdz0, vrdz0, vrdz20, zrp0, zrw0, vctransp0, vctransw0, dvtransp0, dvtransw0
   real, dimension(:,:), allocatable :: zs0, landsea0, lon0, lat0, lonu0, latu0, lonv0, latv0
   real, dimension(:,:), allocatable :: fi, fj, fiu, fju, fiv, fjv
   ! outer fields
   real, dimension(:,:), allocatable :: ps0, pb0, pmsl0, ptop0, ts0, pts0, us0, vs0, qvs0
   real, dimension(:,:,:), allocatable :: u0, v0, w0, pt0, p0
   real, dimension(:,:,:), allocatable :: qv0, qc0, qi0, qr0, qs0, qg0
   real, dimension(:,:,:), allocatable :: tsoil0, h0, t0, pnh0
   ! inner fields
   real, dimension(:,:), allocatable :: deltazs, sst, dptsdt, ps, ts, ptop, pmsl, pts, us, vs, qvs
   real, dimension(:,:,:), allocatable :: u, v, w, pt, pnh
   real, dimension(:,:,:), allocatable :: qv, qc, qi, qr, qs, qg
   real, dimension(:,:,:), allocatable :: tsoil, p, t, rh
   !
   namelist /control/ input_mode, iyear, imonth, iday, ihour, iminute
   namelist /model/ proj, resolution, nx, ny, xi, xj, xlon, xlat, slon, slat, &
                  & nz, dz1, dz2, ztop, vctrans_type, n_vctrans, zl_vctrans, zh_vctrans
   namelist /driving/ proj0, resolution0, nx0, ny0, xi0, xj0, xlon0, xlat0, slon0, slat0, &
                    & nz0, dz10, dz20, ztop0, vctrans_type0, n_vctrans0, zl_vctrans0, zh_vctrans0
end module variable
