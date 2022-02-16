module variable
   implicit none
   !
   character(len=2) :: dtype = 'R4'
   character(len=4) :: member = '    '
   integer, parameter :: mfin = 0, enkf1 = 1, enkf2 = 2, enkf3 = 3, enkf4 = 4, &
                       & nhm = 10, enkf11 = 11, enkf12 = 12, enkf13 = 13, enkf14 = 14
   real, parameter :: rd = 287.05d0, cp = 7.0d0/2.d0*rd, rdvcp = rd/cp, cpdvrd = cp/rd
   real, parameter :: tkelvn = 273.15d0, g = 9.80665d0, gvcp = g/cp
   real, parameter :: tcw = 0.0d0, e0cw = 6.11d2, tetn1w = 17.27d0, tetn2w = 273.15d0, tetn3w = 35.85d0
   real, parameter :: tci = -36.0d0, e0ci = 611.d0, tetn1i = 21.875d0, tetn2i = 273.15d0, tetn3i = 7.65d0
   !
   character(len=8) :: type1, type1s
   character(len=4) :: type2, type2s, type3, type3s
   integer :: iyear, imonth, iday, ihour, iminute, interval, nslot
   !
   character(len=4) :: proj
   integer :: ibase
   integer :: nx, ny, nz
   real :: slat, ptrf, presrf
   !
   real, dimension(:), allocatable :: zrp, zrw, vctrans_p, vctrans_w, dvtrans_p, dvtrans_w, vrdz2
   real, dimension(:,:), allocatable :: zs, lat, pb, us, vs, ts, ps, tds, rhs, qvs
   real, dimension(:,:,:), allocatable :: fmap, g2, tsoil, pairf
   real, dimension(:,:,:), allocatable :: dnsg2, rhou, rhov, rhow, uc, vc, wc, prs
   real, dimension(:,:,:), allocatable :: u, v, w, pt, t, p, pnh, qv, qc, qi, qr, qs, qg
   !
   namelist /control/ type1, type2, type3, type1s, type2s, type3s, iyear, imonth, iday, ihour, iminute, interval, nslot
end module variable
