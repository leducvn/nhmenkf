module nhmlib
   use variable, only : r_size, r_dble, ra, pi
contains
   !
   !
   !
   subroutine read_nhmcst(filename, nx, ny, lon, lat, terrain, landsea_mask)
      implicit none
      !
      character(*), intent(in) :: filename
      integer, intent(in) :: nx, ny
      real, dimension(nx,ny), intent(out) :: lon, lat, terrain, landsea_mask
      integer :: istart, iend, jstart, jend
      real, dimension(nx,ny) :: tmp
      !
      open(40, file=trim(filename), form='unformatted')
      read(40) istart, iend, jstart, jend, terrain, landsea_mask, tmp, tmp, lat, lon
      close(40)
      !
      return
   end subroutine read_nhmcst
   !
   !
   !
   subroutine vrgdis(ii1, ii2, imin, imax, nd, msw, dx, dxl, dxr, vdx, vrdx, vrdx2)
      implicit none
      integer, intent(in) :: ii1, ii2, imin, imax, nd, msw
      real, intent(in) :: dx, dxl, dxr
      real(kind=r_size), dimension(nd), intent(out) :: vdx, vrdx, vrdx2
      integer :: i, i1, i2, ixi1, ix
      real :: c1, c2
      !
      if (imax <= 4) then
         do i = 1, nd
            vdx(i) = dx
            vrdx(i) = 1 / dx
            vrdx2(i) = vrdx(i)
         end do
         return
      end if
      !
      i1 = min(nd, max(ii1, 1))
      i2 = max(2, min(ii2, nd))
      c1 = 0.0
      if (i1 > 2) c1 = (dxl - dx) / abs(float(i1 - 2))
      c2 = 0.0
      if (imax - 1 - i2 > 0) c2 = (dxr - dx) / abs(float(imax - 1 - i2))
      if (msw == 1) then
         ! W-POSITION
	 ixi1 = max(2, min(i1, nd))
	 do ix = 2, ixi1
	    vrdx2(ix) = 1.0 / (c1 * abs(float(i1 - ix)) + dx)
	    vrdx(ix) = 1.0 / (c1 * (abs(float(i1 - ix)) + abs(float(ix - 1 - i1))) / 2.0 + dx)
	 end do
	 do ix = i1, i2
	    vrdx2(ix) = 1 / dx
	 end do
	 ixi1 = min(i2, i1 + 1)
	 do ix = ixi1, i2
	    vrdx(ix) = 1 / dx
	 end do
	 do ix = i2 + 1, imax - 1
	    vrdx2(ix) = 1.0 / (c2 * abs(float(ix - i2)) + dx)
	    vrdx(ix) = 1.0 / (c2 * (abs(float(ix - i2)) + abs(float(ix - 1 - i2))) / 2.0 + dx)
	 end do
	 vrdx2(imax) = vrdx2(imax - 1)
	 vrdx(imax) = vrdx2(imax - 1)
	 vrdx2(1) = vrdx2(2)
	 vrdx(2) = vrdx2(2)
	 vrdx(1) = vrdx2(2)
      else
         ! U,V -- POSITION
	 c1 = 0.0
	 if (i1 > 2) c1 = (dxl - dx) / abs(float(i1 - 2))
	 c2 = 0.0
	 if (imax - 1 - i2 > 0) c2 = (dxr - dx) / abs(float(imax - 1 - i2))
	 ixi1 = max(2, i1 - 1)
	 do ix = 2, ixi1
	    vrdx2(ix) = 1.0 / (c1 * abs(float(ix - i1)) + dx)
	    vrdx(ix) = 1.0 / (c1 * (abs(float(ix - i1)) + abs(float(ix + 1 - i1))) / 2.0 + dx)
	 end do
	 do ix = i1, i2 - 1
	    vrdx2(ix) = 1 / dx
	 end do
	 do ix = i1, i2 - 1
	    vrdx(ix) = 1 / dx
	 end do
	 do ix = i2, imax - 1
	    vrdx2(ix) = 1.0 / (c2 * abs(float(ix - i2)) + dx)
	    vrdx(ix) = 1.0 / (c2 * (abs(float(ix - i2)) + abs(float(ix + 1 - i2))) / 2.0 + dx)
	 end do
	 vrdx2(imax) = vrdx2(imax - 1)
	 vrdx(imax - 1) = vrdx2(imax - 1)
	 vrdx(imax) = vrdx(imax - 1)
	 vrdx2(1) = vrdx2(2)
	 vrdx(1) = vrdx2(2)
      end if
      do ix = 1, nd
         vdx(ix) = 1.0 / vrdx(ix)
      end do
      !
      return
   end subroutine vrgdis
   !
   !
   !
   subroutine setzrp(vrdz2, nz, zrp)
      implicit none
      integer, intent(in) :: nz
      real(kind=r_size), dimension(nz), intent(in) :: vrdz2
      real(kind=r_size), dimension(nz), intent(out) :: zrp
      integer :: kz
      !
      zrp(1) = -0.5 * 1. / vrdz2(1)
      do kz = 2, nz
         zrp(kz) = zrp(kz - 1) + 1. / vrdz2(kz - 1)
      end do
      !write(6, '(1x, "*I SETZRP", /, 5(i3, e11 .4))') (kz, zrp(kz), kz = 1, nz)
      !
      return
   end subroutine setzrp
   !
   !
   !
   subroutine setzrw(vrdz, nz, zrw)
      implicit none
      integer, intent(in) :: nz
      real(kind=r_size), dimension(nz), intent(in) :: vrdz
      real(kind=r_size), dimension(nz), intent(out) :: zrw
      integer :: kz
      !
      zrw(1) = 0.
      do kz = 2, nz
         zrw(kz) = zrw(kz - 1) + 1. / vrdz(kz)
      end do
      !write(6, '(1x, "*I SETZRP", /, 5(i3, e11 .4))') (kz, zrp(kz), kz = 1, nz)
      !
      return
   end subroutine setzrw
   !
   !
   !
   subroutine calc_zcoordinate(mode_in, n_in, zt, zl, zh, zeta, nz, fp, dfdz)
      implicit none
      integer, intent(in) :: mode_in, n_in, nz
      real, intent(in) :: zt, zl, zh
      real(kind=r_size), dimension(nz), intent(in) :: zeta
      real(kind=r_size), dimension(nz), intent(out) :: fp, dfdz
      integer :: kz, mode, n
      real :: p, zmzt, c
      !
      if (mode_in == 0) then
         mode = 1
         n = 1
      else
         mode = mode_in
         n = n_in
      end if
      !
      ! f(p) =  (1-p)^n
      if (mode == 1) then
         do kz = 1, nz
            p = zeta(kz)/zt
            fp(kz) = (1.0-p)**n
            dfdz(kz) = -n*(1.0-p)**(n-1)
            dfdz(kz) = dfdz(kz)/zt
         end do
      ! f(p) = 20p^7 - 70p^6 + 84p^5 -35p^4 + 1
      elseif (mode == 2) then
         do kz = 1, nz
            p = (zeta(kz)-zl)/(zh-zl)
            if (p < 0.0) then
               fp(kz) = 1.0
               dfdz(kz) = 0.0
            elseif (p < 1.0) then
               fp(kz) = 1.0 - 35.*p**4 + 84.*p**5 - 70.*p**6 + 20.*p**7
               dfdz(kz) = -140.*p**3 + 420.*p**4 - 420.*p**5 + 140.*p**6
            else
               fp(kz) = 0.0
               dfdz(kz) = 0.0
            endif
            dfdz(kz) = dfdz(kz)/(zh-zl)
         end do
      ! f(p) = k(1 - p^n) / (k + p^n)
      elseif (mode == 3) then
         zmzt = (zl+zh)*0.5/zt
         c = zmzt**n/(1.0-2*zmzt**n)
         do kz = 1, nz
            p = zeta(kz)/zt
            fp(kz) = c*(1-p**n)/(c+p**n)
            dfdz(kz) = -n*c*p**(n-1)*(1+c)/(c+p**n)**2
            dfdz(kz) = dfdz(kz)/zt
         end do
      endif
      !
      return
   end subroutine calc_zcoordinate
   !
   !
   !
   subroutine gauss(a, w, k)
    !*
    !  A; COSINE OF COLATITUDE
    !  W; GAUSSIAN WEIGHT
    !  K; ORDER OF LEGENDRE FUNCTIONS
    !*
    implicit none
    integer(4), intent(in):: k
    real(r_dble), intent(out):: a(k), w(k)
    real(r_dble):: esp = 1.d-14
    real(r_dble):: c, fk, fn, pk, pkm1, pkm2, pkmrk, xz, sp
    integer(4):: kk, is, iter, n
    integer(4):: itermax = 10
    !*
    c = (1.d0 - (2.d0 / pi) ** 2) * 0.25d0
    fk = real(k)
    kk = k / 2
    call bsslz1(a, kk)
    do is = 1, kk
      xz = cos(a(is) / sqrt((fk + 0.5D0) ** 2 + c))
      do iter = 1, itermax
        pkm2 = 1.d0
        pkm1 = xz
        do n = 2, k
          fn = real(n)
          pk = ((2.d0 * fn - 1.d0) * xz * pkm1 - (fn - 1.d0) * pkm2) / fn
          pkm2 = pkm1
          pkm1 = pk
        end do
        pkm1 = pkm2
        pkmrk = (fk * (pkm1 - xz * pk)) / (1.d0 - xz ** 2)
        sp = pk / pkmrk
        xz = xz - sp
      end do
      if(abs(sp) > esp) then
        write(6,*)'error in gauss'
        stop 1
      end if
      a(is) = xz
      w(is) = (2.d0 * (1.d0 - xz ** 2)) / (fk * pkm1) ** 2
    end do
    if(k /= kk * 2) then
      a(kk + 1) = 0.D0
      pk = 2.d0 / fk ** 2
      do n = 2, k, 2
        fn = real(n)
        pk = pk * fn ** 2 / (fn - 1.d0) ** 2
      end do
      w(kk + 1) = pk
    end if
    do n = 1, kk
      a(k + 1 - n) = -a(n)
      w(k + 1 - n) =  w(n)
    end do
    return
   end subroutine gauss
   !
   !
   !
   subroutine bsslz1(bes, n)
    implicit none
    integer(4), intent(in):: n
    real(r_dble), intent(out):: bes(n)
    real(r_dble):: bz(50) = (/             2.4048255577D0,  5.5200781103D0, &
      &  8.6537279129D0, 11.7915344391D0, 14.9309177086D0, 18.0710639679D0, &
      & 21.2116366299D0, 24.3524715308D0, 27.4934791320D0, 30.6346064684D0, &
      & 33.7758202136D0, 36.9170983537D0, 40.0584257646D0, 43.1997917132D0, &
      & 46.3411883717D0, 49.4826098974D0, 52.6240518411D0, 55.7655107550D0, &
      & 58.9069839261D0, 62.0484691902D0, 65.1899648002D0, 68.3314693299D0, &
      & 71.4729816036D0, 74.6145006437D0, 77.7560256304D0, 80.8975558711D0, &
      & 84.0390907769D0, 87.1806298436D0, 90.3221726372D0, 93.4637187819D0, &
      & 96.6052679510D0, 99.7468198587D0, 102.888374254D0, 106.029930916D0, &
      & 109.171489649D0, 112.313050280D0, 115.454612653D0, 118.596176630D0, &
      & 121.737742088D0, 124.879308913D0, 128.020877005D0, 131.162446275D0, &
      & 134.304016638D0, 137.445588020D0, 140.587160352D0, 143.728733573D0, &
      & 146.870307625D0, 150.011882457D0, 153.153458019D0, 156.295034268D0 /)
    integer(4):: i
    !*
    if(n <= 50) then
      do i = 1, n
        bes(i) = bz(i)
      end do
    else
      do i = 1, 50
        bes(i) = bz(i)
      end do
      do i = 51, n
        bes(i) = bes(i - 1) + pi
      end do 
    endif
   end subroutine bsslz1
   !
   !
   !
   subroutine mapconv_xy(fi, fj, rlat, rlon, sw_conv, &
    & im, jm, nprojc, delsx, delsy, slat, slon, xi, xj, xlat, xlon)
    !---------------------------------------------------------------------!
    !    GRID POINT OF LATITUDE AND LONGITUDE                             !
    !                                                                     !
    !IN/OUT   (FI,FJ)      GRID POINT                                     !
    !         RLAT         LATITUDE   (DEGREE)   -90.<   < 90.            !
    !         RLON         LONGITUDE  (DEGREE)     0.<   <360.            !
    !                                                                     !
    !IN-PUT   SW_CONV      CONVERTER SWITCH                               !
    !            'LLTOIJ'                                                 !
    !                IN-PUT  RLAT, RLON                                   !
    !                OUTPUT  FI, FJ                                       !
    !            'IJTOLL'                                                 !
    !                IN-PUT  FI, FJ                                       !
    !                OUTPUT  RLAT, RLON                                   !
    !         NPRO       MAP PROJECTION                                   !
    !            'PSN ' POLAR STEREO     SLAT= 60 N                       !
    !            'PSS ' POLAR STEREO     SLAT=-60 N                       !
    !            'MER ' MERCATOR         SLAT= 22.5 N                     !
    !            'LMN ' LAMBERT          SLAT= 30 N , 60 N                !
    !            'LMS ' LAMBERT          SLAT=-30 N ,-60 N                !
    !            'LL  ' LATITUDE LONGITUDE                                !
    !            'CE  ' Cylinderical Equidistant projection
    !         DELS         GRID PARAM (     M)                            !
    !                      IF NPRO='LL  ','GAUS' (DEG)                    !
    !         SLAT         STANDARD LATITUDE   ; X-COORDINATE -90.<  < 90.!
    !         SLON         STANDARD LONGITUDE  ; Y-COORDINATE   0.<  <360.!
    !                                                                     !
    !         (XI,XJ) <----------> (XLAT,XLON)                            !
    !                      STANDARD POINT                                 !
    !---------------------------------------------------------------------!
    !*
    implicit none
    integer(4), intent(in):: im, jm
    character(4), intent(in):: nprojc
    character(6), intent(in):: sw_conv
    real(r_size), intent(in):: delsx, delsy, slat, slon, xi, xj, xlat, xlon
    real(r_size), intent(inout):: fi(im, jm), fj(im, jm)
    real(r_size), intent(inout):: rlat(im, jm), rlon(im, jm)
    !*
    real(r_dble):: slat1, slon1, xlat1, xlon1, x0, y0, xp, yp, poi, rpoi, pi4
    real(r_dble):: rl0, ac, acn, ck, rck, t1, fxlon, rrlon, xslon
    real(r_dble):: slata, slatb, slata1, slatb1, slat2
    real(r_dble):: ali, di, dj, rlat1, rlon1, rdelsx, rdelsy, r0, cyclei
    real(r_dble):: dcoslt(1400), dgw(1400), glat(1400), eps, xxxj
    integer(4):: i, j, k, jmg, jx0, jx
    integer(4), parameter ::is_debug = 0
    !*
    if (sw_conv(1:6) /= 'LLTOIJ' .and. sw_conv(1:6) /= 'IJTOLL' ) then
      write(6,*)'warning ! : sw_conv', sw_conv,'is wrong'
      return
    end if
    poi = pi/180.d0
    rpoi = 1.d0/poi
    pi4 = pi/4.d0
    !*
    !-----------------------------------POLAR STEREO NORTH-----------------
    !*
    if (nprojc(1:3) == 'PSN') then
      slat1 = +slat * poi
      slon1 = +slon * poi
      xlat1 = +xlat * poi
      xlon1 = +xlon * poi
      rl0 = ra * (1.d0 + sin(slat1)) * tan(0.5d0 * (0.5d0 * pi - xlat1))
      x0 = rl0 * sin(xlon1 - slon1)
      y0 = rl0 * cos(xlon1 - slon1)
      !*
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
            rlat1 = +rlat(i, j) * poi
            rlon1 = +rlon(i, j) * poi
            ali = ra * (1.d0 + sin(slat1)) * tan(0.5d0 * (0.5d0 * pi - rlat1))
            di = ali * sin(rlon1 - slon1)
            dj = ali * cos(rlon1 - slon1)
            fi(i, j) = (di - x0) / delsx + xi
            fj(i, j) = (dj - y0) / delsy + xj
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            xp = x0 + (fi(i, j) - xi) * delsx
            yp = y0 + (fj(i, j) - xj) * delsy
            if ((xp == 0._r_dble) .and. (yp == 0._r_dble)) then
              rlat(i, j) = + 90._r_dble
              rlon(i, j) = slon
            else
              rlat(i, j) = + 90. - rpoi * 2. * atan( &
                & sqrt(xp * xp + yp * yp) / (ra * (1. + sin(slat1))))
              rlon(i, j) = slon + rpoi * atan2(xp, yp)
            end if
          end do
        end do
      end if
    end if
    !*
    !-----------------------------------POLAR STEREO SOUTH-----------------
    !*
    if (nprojc(1:3) == 'PSS') then
      slat1 = -slat * poi
      slon1 = +slon * poi
      xlat1 = -xlat * poi
      xlon1 = +xlon * poi
      rl0 = ra * (1.d0 + sin(slat1)) * tan(0.5d0 * (0.5d0 * pi - xlat1))
      x0 = rl0 * sin(xlon1 - slon1)
      y0 = rl0 * cos(xlon1 - slon1)
      !*
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
            rlat1 = -rlat(i, j) * poi
            rlon1 = +rlon(i, j) * poi
            ali = ra * (1. + sin(slat1)) * tan(0.5d0 * (0.5d0 * pi - rlat1))
            di = ali * sin(rlon1 - slon1)
            dj = ali * cos(rlon1 - slon1)
            fi(i, j) = (di - x0) / delsx + xi
            fj(i, j) = (y0 - dj) / delsy + xj
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            xp = x0 + (fi(i, j) - xi) * delsx
            yp = y0 + (xj - fj(i, j)) * delsy
            if ((xp == 0._r_dble) .and. (yp == 0._r_dble)) then
              rlat(i, j) = - 90._r_dble
              rlon(i, j) = slon
            else
              rlat(i, j) = - 90. + rpoi * 2.* atan( &
                & sqrt(xp * xp + yp * yp) / (ra * (1. + sin(slat1))))
              rlon(i, j) = slon + rpoi * atan2(xp, yp)
            end if
          end do
        end do
      end if
    end if
    !*
    !-----------------------------------MERCATOR---------------------------
    !*
    if (nprojc(1:3) == 'MER') then
      slat1 = slat * poi
      slon1 = slon * poi
      xlat1 = xlat * poi
      xlon1 = xlon * poi
      ac = ra * cos(slat1)
      x0 = 0.d0
      y0 = ac * log((1.d0 + sin(xlat1)) / cos(xlat1))
      !*
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
          rlat1 = rlat(i, j) * poi
          rlon1 = rlon(i, j) * poi
          di = rlon(i, j) - xlon
          di = mod(di + 900.d0, 360.d0) - 180.d0
          di = ac * (di * poi)
          dj = ac * log((1.d0 + sin(rlat1)) / cos(rlat1))
          fi(i, j) = (di - x0) / delsx + xi
          fj(i, j) = (y0 - dj) / delsy + xj
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            rlon(i, j) = xlon + rpoi * (fi(i, j) - xi) * delsx / ac
            t1 = exp((xj - fj(i, j)) * delsx / ac) * (1.d0 + sin(xlat1)) &
              & / cos(xlat1)
            rlat(i, j) = rpoi * asin((t1 * t1 - 1.d0) / (t1 * t1 + 1.d0))
          end do
        end do
      end if
    end if
    !*
    !-----------------------------------LAMBERT NORTH----------------------
    !*
    if (nprojc(1:3) == 'LMN') then
      slata = +30.
      slatb = +60.
      slata1 = +slata * poi
      slatb1 = +slatb * poi
      slat1 = pi4 - slata1 * 0.5d0
      slat2 = pi4 - slatb1 * 0.5d0
      slon1 = slon * poi
      xlat1 = +xlat * poi
      xlon1 = +xlon * poi
      ck = log(cos(slata1) / cos(slatb1)) / log(tan(slat1) / tan(slat2))
      rck = 1.d0 / ck
      acn = ra * cos(slata1) * rck
      r0 = acn / (tan(slat1)) ** ck
      rl0 = r0 * (tan(pi4 - xlat1 * 0.5d0)) ** ck
      xslon = xlon - slon
      xslon = mod(xslon + 900.d0, 360.d0) - 180.d0
      xslon = xslon * poi * ck
      x0 = rl0 * sin(xslon)
      y0 = rl0 * cos(xslon)
      !*
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
            rlat1 = +rlat(i, j) * poi
            rlon1 = +rlon(i, j) * poi
            ali = r0 * (tan(pi4 - rlat1 * 0.5d0)) ** ck
            fxlon = rlon(i, j) - xlon
            fxlon = mod(fxlon + 900.d0, 360.d0) - 180.d0
            fxlon = fxlon * poi * ck
            di = ali * sin(fxlon + xslon)
            dj = ali * cos(fxlon + xslon)
            fi(i, j) = (di - x0) / delsx + xi
            fj(i, j) = (dj - y0) / delsy + xj
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            xp = x0 + (fi(i, j) - xi) * delsx
            yp = y0 + (fj(i, j) - xj) * delsy
            if ((xp == 0._r_dble) .and. (yp == 0._r_dble)) then
              rlat(i, j) = + 90._r_dble
              rlon(i, j) = slon
            else
              rlat(i, j) = + 90. - 2. * rpoi * atan( &
                & (sqrt(xp * xp + yp * yp) / acn) ** rck * tan(slat1))
              rlon(i, j) = slon + rck * rpoi * atan2(xp, yp)
              rrlon = rlon(i, j)
              rrlon = rrlon - xlon
              rrlon = mod(rrlon + 900., 360._r_dble) - 180.
              rlon(i, j) = rrlon + xlon
            end if
          end do
        end do
      end if
    end if
    !*
    !-----------------------------------LAMBERT SOUTH----------------------
    !*
    if (nprojc(1:3) == 'LMS') then
      slata = -30.
      slatb = -60.
      slata1 = -slata * poi
      slatb1 = -slatb * poi
      slat1 = pi4 - slata1 * 0.5d0
      slat2 = pi4 - slatb1 * 0.5d0
      slon1 = slon * poi
      xlat1 = -xlat * poi
      xlon1 = +xlon * poi
      ck = log(cos(slata1) / cos(slatb1)) / log(tan(slat1) / tan(slat2))
      rck = 1.d0 / ck
      acn = ra * cos(slata1) * rck
      r0 = acn / (tan(slat1)) ** ck
      rl0 = r0 * (tan(pi4 - xlat1 * 0.5d0)) ** ck
      xslon = xlon - slon
      xslon = mod(xslon + 900.d0, 360.d0) - 180.d0
      xslon = xslon * poi * ck
      x0 = rl0 * sin(xslon)
      y0 = rl0 * cos(xslon)
      !*
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
            rlat1 = -rlat(i, j) * poi
            rlon1 = +rlon(i, j) * poi
            ali = r0 * (tan(pi4 - rlat1 * 0.5d0)) ** ck
            fxlon = rlon(i, j) - xlon
            fxlon = mod(fxlon + 900.d0, 360.d0) - 180.d0
            fxlon = fxlon * poi * ck
            di = ali * sin(fxlon + xslon)
            dj = ali * cos(fxlon + xslon)
            fi(i, j) = (di - x0) / delsx + xi
            fj(i, j) = (y0 - dj) / delsy + xj
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            xp = x0 + (fi(i, j) - xi) * delsx
            yp = y0 + (xj - fj(i, j)) * delsy
            if ((xp == 0._r_dble) .and. (yp == 0._r_dble)) then
              rlat(i, j) = - 90._r_dble
              rlon(i, j) = slon
            else
              rlat(i, j) = - 90. + 2. * rpoi * atan( &
                & (sqrt(xp * xp + yp * yp) / acn) ** rck * tan(slat1))
              rlon(i, j) = slon + rck * rpoi * atan2(xp, yp)
              rrlon = rlon(i, j)
              rrlon = rrlon - xlon
              rrlon = mod(rrlon + 900., 360._r_dble) - 180.
              rlon(i, j) = rrlon + xlon
            end if
          end do
        end do
      end if
    end if
    !*
    !-----------------------------------LATITUDE-LONGITUDE-----------------
    !*
    if (nprojc(1:2) == 'LL') then
      rdelsx = 1.d0 / delsx
      rdelsy = 1.d0 / delsy
      cyclei = 360.d0 / delsx
      !*
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
            fi(i, j) = xi + rdelsx * (rlon(i, j) - xlon)
            fj(i, j) = xj - rdelsy * (rlat(i, j) - xlat)
            if (fi(i, j) <  0.0d0) fi(i, j) = fi(i, j) + cyclei
            if (fi(i, j) > cyclei) fi(i, j) = fi(i, j) - cyclei
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            rlon(i, j) = xlon + (fi(i, j) - xi) * delsx
            rlat(i, j) = xlat - (fj(i, j) - xj) * delsy
          end do
        end do
      end if
    end if
    !*
    ! --- Cylinderical Equidistant projection ---
    !*
    if (nprojc == 'CE  ') then
      ac = ra * cos(slat * poi)
      if (sw_conv == 'LLTOIJ' ) then
        do j = 1, jm
          do i = 1, im
            di = rlon(i, j) - xlon
            di = mod(di + 900., 360._r_dble) - 180.
            di = ac * (di * poi)
            dj = rlat(i, j) - xlat
            dj = ra * (dj * poi)
            fi(i, j) = xi + di / delsx
            fj(i, j) = xj - dj / delsy
          end do
        end do
      else if (sw_conv == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            rlon(i, j) = xlon + rpoi * (fi(i, j) - xi) * delsx / ac
            rlat(i, j) = xlat - rpoi * (fj(i, j) - xj) * delsy / ra
          end do
        end do
      end if
    end if
    !*
    ! ----------------------------------GAUSSIAN GRID----------------------
    !*
    if (nprojc(1:1) == 'G') then
      jmg = int(slon + 0.01d0)
      if (mod(jmg, 2) == 1 .or. jmg > 1400) then
        write(6, *) 'DIMENSION OF GAUSSIAN LATITUDE IS INVALID'
        write(6, *) 'JMG=', jmg
        stop 90
      end if
      call gauss(dcoslt, dgw, jmg)
      if (is_debug == 1) then
        do j = 1, jmg
          write(6,*)dcoslt(j), dgw(j)
        end do
      end if
      do j = 1, jmg
        glat(j) = 90.d0 - rpoi * acos(dcoslt(j))
      end do
      if (abs(xlat) < 0.01d0) then
        jx0 = jmg / 2
      else
        write(6, *) ' ----- SORRY ! ------------------------- '
        write(6, *) '  THIS PROGRAM ONLY SUPORT FOR XLAT=0.0  '
        write(6, *) '        AT GASSIAN LATITUDE              '
        write(6, *) ' --------------------------------------- '
        stop 91
        ! --------------------- FOLLOWING COMMENT PROGRAM HAS NOT CHECKED,YET.
        !         ---- SEARCH JX0
        !         DO 62 J=1,JMG
        !             IF( ABS(GLAT(J)-XLAT).LT.0.1 ) THEN
        !                  JX0=J
        !                  GOTO 63
        !             ENDIF
        !  62     CONTINUE
        !             WRITE (6,*) ' I CAN NOT FIND JX0 FOR  GAUSIAN GRID '
        !             STOP 900
        !  63     CONTINUE
        ! ---------------------------------------------------------------------
      end if
      eps = 0.01d0
      if (xj < 0.5d0) eps = -0.01d0
      rdelsx = 1.d0 / delsx
      rdelsy = 1.d0 / delsy
      cyclei = 360.d0 / delsx
      if (sw_conv(1:6) == 'LLTOIJ' ) then
        ! ---- FI
        do j = 1, jm
          do i = 1, im
          fi(i, j) = xi + rdelsx * (rlon(i, j) - xlon)
          if (fi(i, j) <  0.0d0) fi(i, j) = fi(i, j) + cyclei
          if (fi(i, j) > cyclei) fi(i, j) = fi(i, j) - cyclei
          end do
        end do
        ! ---- FJ
        do j = 1, jm
          do i = 1, im
            xxxj = (glat(1) - rlat(i, j)) / (glat(1) - glat(jmg)) * &
              & (jmg - 1.d0) + 1
            jx = int(xxxj)
            if (xxxj <= 1.0d0) then
              jx0 = 1
            elseif(xxxj >= float(jmg) - 1.d0) then
              jx0 = jmg - 1
            else
              if(rlat(i, j) >= glat(jx + 1) .and. &
               & rlat(i, j) <  glat(jx    )) then
                jx0 = jx
              end if
              if(rlat(i, j) >= glat(jx    ) .and. &
               & rlat(i, j) <  glat(jx - 1)) then
                jx0 = jx - 1
              end if
              if(rlat(i, j) >= glat(jx + 2) .and. &
               & rlat(i, j) <  glat(jx + 1)) then
                jx0 = jx + 1
              end if
            end if
            fj(i, j) = jx0 + (rlat(i, j) - glat(jx0)) / &
             & (glat(jx0 + 1) - glat(jx0))
          end do
        end do
      end if
      if (sw_conv(1:6) == 'IJTOLL' ) then
        do j = 1, jm
          do i = 1, im
            rlon(i, j) = xlon + (fi(i, j) - xi) * delsx
            rlat(i, j) = glat(nint(fj(i, j)) - int(xj + eps) + jx0)
          end do
        end do
      end if
    end if
    !*
    !----------------------------------------------------------------------
    !*
    if (sw_conv(1:6) == 'IJTOLL' ) then
      do j = 1, jm
        do i = 1, im
          if (rlat(i, j) > +90.d0) rlat(i, j) = +90.d0
          if (rlat(i, j) < -90.d0) rlat(i, j) = -90.d0
        end do
      end do
      do k = 1, 2
        do j = 1, jm
          do i = 1, im
            if (rlon(i, j) <   0.d0) rlon(i, j) = rlon(i, j) + 360.d0
            if (rlon(i, j) > 360.d0) rlon(i, j) = rlon(i, j) - 360.d0
          end do
        end do
      end do
    end if
    !*
    return
   end subroutine mapconv_xy
   !
   !
   !
   subroutine ij2lonlat(flon, flat, fi, fj, nprojc, dels, xi, xj, xlon, xlat, slon, slat)
      implicit none
      real(r_size), intent(in):: fi, fj
      character(4), intent(in):: nprojc
      real(r_size), intent(in):: dels, xi, xj, xlon, xlat, slon, slat
      real(r_size), intent(out):: flon, flat
      !*
      real(r_size) :: gi(1,1), gj(1,1), rlon(1,1), rlat(1,1)
      gi(1,1) = fi
      gj(1,1) = fj
      !*
      call mapconv_xy(gi, gj, rlat, rlon, 'IJTOLL', 1, 1, nprojc, dels, dels, slat, slon, xi, xj, xlat, xlon)
      flon = rlon(1,1)
      flat = rlat(1,1)
      !*
      return
   end subroutine ij2lonlat
   !
   !
   !
   subroutine compute_qvs(p, t, qvs, nx, ny ,nz)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: p, t
      real(r_size), dimension(nx,ny,nz), intent(out) :: qvs
      !
      integer :: i, j, k
      real(r_size) :: e0c, al, bl, e0i, ai, bi, tc, qvs0, qvs1
      !
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tc = t(i,j,k) - tkelvn
               if (tc >= tcw ) then
                  qvs(i,j,k) = e0c*exp(al*tc/(bl+tc))
	       else if(tc <= tci) then
                  qvs(i,j,k) = e0i*exp(ai*tc/(bi+tc))
	       else
		  qvs0 = e0i*exp(ai*tc/(bi+tc))
	          qvs1 = e0c*exp(al*tc/(bl+tc))
                  qvs(i,j,k) = qvs1*(tc-tci)/(tcw-tci) + qvs0*(tcw-tc)/(tcw-tci)
	       end if
               qvs(i,j,k) = 0.622d0*qvs(i,j,k)/p(i,j,k)
	    end do
         end do
      end do
      !
      return
   end subroutine compute_qvs
   !
   !
   !
   subroutine pst2p(zrp, dvtrans, zs, ps, t, qv, p, nx, ny, nz)
   ! modified so that t, qv at 1st level are t, qv at 2nd level
      use variable, only : g0, rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nz), intent(in) :: zrp, dvtrans
      real(r_size), dimension(nx,ny), intent(in) :: zs, ps
      real(r_size), dimension(nx,ny,nz), intent(in) :: t, qv
      real(r_size), dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real(r_size), dimension(nx,ny) :: factor
      real(r_size), dimension(nx,ny,nz) :: g2w, deltaz
      !
      do k = 1, nz-1
         g2w(:,:,k) = 1. + dvtrans(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      !factor = 0.5*(t(:,:,1)*(1.+0.608*qv(:,:,1))+t(:,:,2)*(1.+0.608*qv(:,:,2)))
      factor = t(:,:,2)*(1.+0.608*qv(:,:,2))
      factor = g0/rd*deltaz(:,:,1)/factor(:,:)
      p(:,:,1) = ps(:,:)*exp(-factor(:,:))
      do k = 2, nz
         if (k == 2) then
	    factor = t(:,:,k)*(1.+0.608*qv(:,:,k))
	 else
	    factor = 0.5*(t(:,:,k-1)*(1.+0.608*qv(:,:,k-1))+t(:,:,k)*(1.+0.608*qv(:,:,k)))
	 end if
	 factor = g0/rd*deltaz(:,:,k)/factor(:,:)
	 p(:,:,k) = p(:,:,k-1)*exp(-factor(:,:))
      end do
      p(:,:,1) = ps(:,:)
      !
      return
   end subroutine pst2p
   !
   !
   !
   subroutine pbpt2p(zrp, dvtransw, zs, pb, pt, qv, p, nx, ny, nz)
      use variable, only : presrf, rdvcp, gvcp, cpdvrd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nz), intent(in) :: zrp, dvtransw
      real(r_size), dimension(nx,ny), intent(in) :: zs, pb
      real(r_size), dimension(nx,ny,nz), intent(in) :: pt, qv
      real(r_size), dimension(nx,ny,nz), intent(out) :: p
      !
      integer :: k
      real(r_size), dimension(nx,ny) :: paib, ptv
      real(r_size), dimension(nx,ny,nz) :: g2w, deltaz, pai, ptm, ptratio
      !
      do k = 1, nz-1
	 g2w(:,:,k) = 1. + dvtransw(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      ! ptm
      ptratio = 1.d0+0.608d0*qv
      ptm = pt*ptratio
      ! paib
      paib = (pb/presrf)**rdvcp
      ! hydrostatic pai
      pai(:,:,1) = paib
      do k = 2, nz
	 ptv = 0.5d0*(ptm(:,:,k)+ptm(:,:,k-1))
	 pai(:,:,k) = pai(:,:,k-1) - gvcp*deltaz(:,:,k)/ptv
      end do
      p = presrf*pai**cpdvrd
      !
      return
   end subroutine pbpt2p
   !
   !
   !
   subroutine compute_t(p, pt, t, nx, ny, nz)
      use variable, only : presrf, rdvcp
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: p, pt
      real(r_size), dimension(nx,ny,nz), intent(out) :: t
      !
      t = pt*(p/presrf)**rdvcp
      !
      return
   end subroutine compute_t
   !
   !
   !
   subroutine compute_pt(p, t, pt, nx, ny, nz)
      use variable, only : presrf, rdvcp
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: p, t
      real(r_size), dimension(nx,ny,nz), intent(out) :: pt
      !
      pt = t*(presrf/p)**rdvcp
      !
      return
   end subroutine compute_pt
   !
   !
   !
   subroutine compute_rh(p, t, qv, rh, nx, ny ,nz)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: p, t, qv
      real(r_size), dimension(nx,ny,nz), intent(out) :: rh
      !
      real(r_size), parameter :: alfa = 0.378d0, beta = 0.622d0
      integer :: i, j, k
      real(r_size) :: e0c, al, bl, e0i, ai, bi, tc, es0, es1
      real(r_size), dimension(nx,ny,nz) :: es, e
      !
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tc = t(i,j,k) - tkelvn
               if (tc >= tcw) then
                  es(i,j,k) = e0c*exp(al*tc/(bl+tc))
	       else if(tc <= tci) then
                  es(i,j,k) = e0i*exp(ai*tc/(bi+tc))
	       else
		  es0 = e0i*exp(ai*tc/(bi+tc))
	          es1 = e0c*exp(al*tc/(bl+tc))
                  es(i,j,k) = es1*(tc-tci)/(tcw-tci) + es0*(tcw-tc)/(tcw-tci)
	       end if
            end do
         end do
      end do
      !
      e = qv*p/(alfa*qv+beta)
      rh = e/es
      !
      return
   end subroutine compute_rh
   !
   !
   !
   subroutine compute_rhs(p, t, qv, rh, nx, ny)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny
      real(r_size), dimension(nx,ny), intent(in) :: p, t, qv
      real(r_size), dimension(nx,ny), intent(out) :: rh
      !
      real(r_size), parameter :: alfa = 0.378d0, beta = 0.622d0
      integer :: i, j
      real(r_size) :: e0c, al, bl, e0i, ai, bi, tc, es0, es1
      real(r_size), dimension(nx,ny) :: es, e
      !
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      do j = 1, ny
	 do i = 1, nx
	    tc = t(i,j) - tkelvn
	    if (tc >= tcw) then
	       es(i,j) = e0c*exp(al*tc/(bl+tc))
	    else if(tc <= tci) then
	       es(i,j) = e0i*exp(ai*tc/(bi+tc))
	    else
	       es0 = e0i*exp(ai*tc/(bi+tc))
	       es1 = e0c*exp(al*tc/(bl+tc))
	       es(i,j) = es1*(tc-tci)/(tcw-tci) + es0*(tcw-tc)/(tcw-tci)
	    end if
	 end do
      end do
      !
      e = qv*p/(alfa*qv+beta)
      rh = e/es
      !
      return
   end subroutine compute_rhs
   !
   !
   !
   subroutine compute_pwv(vdz, dvtransp, zs, p, t, qv, pwv, nx, ny, nz)
      use variable, only : rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nz), intent(in) :: vdz, dvtransp
      real(r_size), dimension(nx,ny), intent(in) :: zs
      real(r_size), dimension(nx,ny,nz), intent(in) :: p, t, qv
      real(r_size), dimension(nx,ny), intent(out) :: pwv
      !
      integer :: k
      real(r_size), dimension(nx,ny) :: tmp
      real(r_size), dimension(nx,ny,nz) :: g2, rhog2
      !
      do k = 1, nz
         g2(:,:,k) = 1. + zs*dvtransp(k)
      end do
      rhog2 = p/(rd*t)*g2
      pwv = 0.d0
      do k = 2, nz-1
         pwv = pwv + qv(:,:,k)*rhog2(:,:,k)*vdz(k)
      end do
      !
      return
   end subroutine compute_pwv
   !
   !
   !
   subroutine compute_ps(zrp, dvtransw, zs, pb, pt, qv, ps, nx, ny, nz)
      use variable, only : presrf, rdvcp, gvcp, cpdvrd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nz), intent(in) :: zrp, dvtransw
      real(r_size), dimension(nx,ny), intent(in) :: zs, pb
      real(r_size), dimension(nx,ny,nz), intent(in) :: pt, qv
      real(r_size), dimension(nx,ny), intent(out) :: ps
      !
      real(r_size), dimension(nx,ny) :: deltaz, ptv, paib, pais
      real(r_size), dimension(nx,ny,2) :: ptm, ptratio
      !
      deltaz = zrp(1)*(1.+dvtransw(1)*zs)
      ptratio = 1.d0+0.608d0*qv(:,:,1:2)
      ptm = pt(:,:,1:2)*ptratio
      ptv = 0.5d0*(ptm(:,:,1)+ptm(:,:,2))
      paib = (pb/presrf)**rdvcp
      pais = paib + gvcp*deltaz/ptv
      ps = presrf*pais**cpdvrd
      !
      return
   end subroutine compute_ps
   !
   !
   !
   subroutine compute_ts(ps, pts, ts, nx, ny)
      use variable, only : presrf, rdvcp
      implicit none
      !
      integer, intent(in) :: nx, ny
      real(r_size), dimension(nx,ny), intent(in) :: ps, pts
      real(r_size), dimension(nx,ny), intent(out) :: ts
      !
      ts = pts*(ps/presrf)**rdvcp
      !
      return
   end subroutine compute_ts
   !
   !
   !
   subroutine compute_pts(ps, ts, pts, nx, ny)
      use variable, only : presrf, rdvcp
      implicit none
      !
      integer, intent(in) :: nx, ny
      real(r_size), dimension(nx,ny), intent(in) :: ps, ts
      real(r_size), dimension(nx,ny), intent(out) :: pts
      !
      pts = ts*(presrf/ps)**rdvcp
      !
      return
   end subroutine compute_pts
   !
   !
   !
end module nhmlib