module outer2innerlib
contains
   !
   !
   !
   subroutine read_nhmcst(fileunit, nx, ny, lon, lat, terrain, landsea_mask)
      implicit none
      !
      integer, intent(in) :: fileunit, nx, ny
      real, dimension(nx,ny), intent(out) :: lon, lat, terrain, landsea_mask
      integer :: istart, iend, jstart, jend
      real, dimension(nx,ny) :: tmp
      !
      read(fileunit) istart, iend, jstart, jend, terrain, landsea_mask, tmp, tmp, lat, lon
      close(fileunit)
      !
      return
   end subroutine read_nhmcst
   !
   !
   !
   subroutine bsslz1(bes, n)
      implicit none
      integer, intent(in):: n
      real, intent(out):: bes(n)
      real :: bz(50) = (/             2.4048255577,  5.5200781103, &
      &  8.6537279129, 11.7915344391, 14.9309177086, 18.0710639679, &
      & 21.2116366299, 24.3524715308, 27.4934791320, 30.6346064684, &
      & 33.7758202136, 36.9170983537, 40.0584257646, 43.1997917132, &
      & 46.3411883717, 49.4826098974, 52.6240518411, 55.7655107550, &
      & 58.9069839261, 62.0484691902, 65.1899648002, 68.3314693299, &
      & 71.4729816036, 74.6145006437, 77.7560256304, 80.8975558711, &
      & 84.0390907769, 87.1806298436, 90.3221726372, 93.4637187819, &
      & 96.6052679510, 99.7468198587, 102.888374254, 106.029930916, &
      & 109.171489649, 112.313050280, 115.454612653, 118.596176630, &
      & 121.737742088, 124.879308913, 128.020877005, 131.162446275, &
      & 134.304016638, 137.445588020, 140.587160352, 143.728733573, &
      & 146.870307625, 150.011882457, 153.153458019, 156.295034268 /)
      integer :: i
      real :: pi = 3.141592653589793
      !
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
      !
      return
   end subroutine bsslz1
   !
   !
   !
   subroutine gauss(a, w, k)
   !  A; COSINE OF COLATITUDE
   !  W; GAUSSIAN WEIGHT
   !  K; ORDER OF LEGENDRE FUNCTIONS
      implicit none
      integer, intent(in) :: k
      real, intent(out) :: a(k), w(k)
      real :: esp = 1.e-14, pi = 3.141592653589793
      real :: c, fk, fn, pk, pkm1, pkm2, pkmrk, xz, sp
      integer :: kk, is, iter, n
      integer :: itermax = 10
      !
      c = (1. - (2. / pi) ** 2) * 0.25
      fk = real(k)
      kk = k / 2
      call bsslz1(a, kk)
      do is = 1, kk
	 xz = cos(a(is) / sqrt((fk + 0.5) ** 2 + c))
	 do iter = 1, itermax
	   pkm2 = 1.
	   pkm1 = xz
	   do n = 2, k
	     fn = real(n)
	     pk = ((2. * fn - 1.) * xz * pkm1 - (fn - 1.) * pkm2) / fn
	     pkm2 = pkm1
	     pkm1 = pk
	   end do
	   pkm1 = pkm2
	   pkmrk = (fk * (pkm1 - xz * pk)) / (1. - xz ** 2)
	   sp = pk / pkmrk
	   xz = xz - sp
	 end do
	 if(abs(sp) > esp) then
	    write(6,*)'error in gauss'
	    stop 1
	 end if
	 a(is) = xz
	 w(is) = (2. * (1. - xz ** 2)) / (fk * pkm1) ** 2
      end do
      if(k /= kk * 2) then
	 a(kk + 1) = 0.
	 pk = 2. / fk ** 2
	 do n = 2, k, 2
	   fn = real(n)
	   pk = pk * fn ** 2 / (fn - 1.) ** 2
	 end do
	 w(kk + 1) = pk
      end if
      do n = 1, kk
         a(k + 1 - n) = -a(n)
         w(k + 1 - n) =  w(n)
      end do
      !
      return
   end subroutine gauss
   !
   !
   !
   subroutine ij2ll(fi, fj, rlat, rlon, im, jm, nprojc, delsx, delsy, slat, slon, xi, xj, xlat, xlon)
      use variable, only : ra
      implicit none
      integer, intent(in):: im, jm
      character(4), intent(in):: nprojc
      real, intent(in):: delsx, delsy, slat, slon, xi, xj, xlat, xlon
      real, dimension(im,jm), intent(in):: fi, fj
      real, dimension(im,jm), intent(out):: rlat, rlon
      !
      real :: slat1, slon1, xlat1, xlon1, x0, y0, xp, yp
      real :: rl0, ac, acn, ck, rck, t1, fxlon, rrlon, xslon
      real :: slata, slatb, slata1, slatb1, slat2
      real :: ali, di, dj, rlat1, rlon1, rdelsx, rdelsy, r0, cyclei
      real :: dcoslt(1400), dgw(1400), glat(1400), eps, xxxj
      integer :: i, j, k, jmg, jx0, jx
      integer, parameter :: is_debug = 0
      real, parameter :: pi = 3.141592653589793, pi4 = pi*0.25, poi = pi/180., rpoi = 180./pi
      !
      ! MERCATOR
      if (nprojc(1:3) == 'MER') then
	 slat1 = slat * poi
	 slon1 = slon * poi
	 xlat1 = xlat * poi
	 xlon1 = xlon * poi
	 ac = ra * cos(slat1)
	 x0 = 0.
	 y0 = ac * log((1. + sin(xlat1)) / cos(xlat1))
	 do j = 1, jm
	    do i = 1, im
	       rlon(i, j) = xlon + rpoi * (fi(i, j) - xi) * delsx / ac
	       t1 = exp((xj - fj(i, j)) * delsx / ac) * (1. + sin(xlat1))/cos(xlat1)
	       rlat(i, j) = rpoi * asin((t1 * t1 - 1.) / (t1 * t1 + 1.))
	    end do
	 end do
      end if
      !
      ! LAMBERT NORTH
      if (nprojc(1:3) == 'LMN') then
	 slata = +30.
	 slatb = +60.
	 slata1 = +slata * poi
	 slatb1 = +slatb * poi
	 slat1 = pi4 - slata1 * 0.5
	 slat2 = pi4 - slatb1 * 0.5
	 slon1 = slon * poi
	 xlat1 = +xlat * poi
	 xlon1 = +xlon * poi
	 ck = log(cos(slata1) / cos(slatb1)) / log(tan(slat1) / tan(slat2))
	 rck = 1. / ck
	 acn = ra * cos(slata1) * rck
	 r0 = acn / (tan(slat1)) ** ck
	 rl0 = r0 * (tan(pi4 - xlat1 * 0.5)) ** ck
	 xslon = xlon - slon
	 xslon = mod(xslon + 900., 360.) - 180.
	 xslon = xslon * poi * ck
	 x0 = rl0 * sin(xslon)
	 y0 = rl0 * cos(xslon)
	 !
	 do j = 1, jm
            do i = 1, im
               xp = x0 + (fi(i, j) - xi) * delsx
               yp = y0 + (fj(i, j) - xj) * delsy
               if ((xp == 0.) .and. (yp == 0.)) then
                  rlat(i, j) = + 90.
                  rlon(i, j) = slon
               else
                  rlat(i, j) = + 90. - 2. * rpoi * atan((sqrt(xp * xp + yp * yp) / acn) ** rck * tan(slat1))
		  rlon(i, j) = slon + rck * rpoi * atan2(xp, yp)
		  rrlon = rlon(i, j)
		  rrlon = rrlon - xlon
		  rrlon = mod(rrlon + 900., 360.) - 180.
		  rlon(i, j) = rrlon + xlon
	       end if
	    end do
         end do
      endif
      !
      ! LATITUDE-LONGITUDE
      if (nprojc(1:2) == 'LL') then
	 rdelsx = 1. / delsx
	 rdelsy = 1. / delsy
	 cyclei = 360. / delsx
	 do j = 1, jm
	    do i = 1, im
	       rlon(i, j) = xlon + (fi(i, j) - xi) * delsx
	       rlat(i, j) = xlat - (fj(i, j) - xj) * delsy
	    end do
	 end do
      end if
      !
      ! Cylinderical Equidistant projection ---
      if (nprojc == 'CE  ') then
	 ac = ra * cos(slat * poi)
	 do j = 1, jm
	    do i = 1, im
	       rlon(i, j) = xlon + rpoi * (fi(i, j) - xi) * delsx / ac
	       rlat(i, j) = xlat - rpoi * (fj(i, j) - xj) * delsy / ra
	    end do
	 end do
      end if
      !
      ! GAUSSIAN GRID
      if (nprojc(1:1) == 'G') then
	 jmg = int(slon + 0.01)
	 if (mod(jmg, 2) == 1 .or. jmg > 1400) then
	    write(6, *) 'DIMENSION OF GAUSSIAN LATITUDE IS INVALID'
	    write(6, *) 'JMG=', jmg
	    stop 90
	 end if
	 call gauss(dcoslt, dgw, jmg)
	 do j = 1, jmg
	    glat(j) = 90. - rpoi * acos(dcoslt(j))
	 end do
	 if (abs(xlat) < 0.01) then
	    jx0 = jmg / 2
	 else
	    write(6, *) ' ----- SORRY ! ------------------------- '
	    write(6, *) '  THIS PROGRAM ONLY SUPORT FOR XLAT=0.0  '
	    write(6, *) '        AT GASSIAN LATITUDE              '
	    write(6, *) ' --------------------------------------- '
	 end if
	 eps = 0.01
	 if (xj < 0.5) eps = -0.01
	 rdelsx = 1. / delsx
	 rdelsy = 1. / delsy
	 cyclei = 360. / delsx
	 !
	 do j = 1, jm
            do i = 1, im
               rlon(i, j) = xlon + (fi(i, j) - xi) * delsx
               rlat(i, j) = glat(nint(fj(i, j)) - int(xj + eps) + jx0)
            end do
         end do
      end if
      !
      do j = 1, jm
	 do i = 1, im
	    if (rlat(i, j) > +90.) rlat(i, j) = +90.
	    if (rlat(i, j) < -90.) rlat(i, j) = -90.
	 end do
      end do
      do k = 1, 2
	 do j = 1, jm
	    do i = 1, im
	       if (rlon(i, j) <   0.) rlon(i, j) = rlon(i, j) + 360.
	       if (rlon(i, j) > 360.) rlon(i, j) = rlon(i, j) - 360.
	    end do
	 end do
      end do
      !
      return
   end subroutine ij2ll
   !
   !
   !
   subroutine ll2ij(fi, fj, rlat, rlon, im, jm, nprojc, delsx, delsy, slat, slon, xi, xj, xlat, xlon)
      use variable, only : ra
      implicit none
      integer, intent(in):: im, jm
      character(4), intent(in):: nprojc
      real, intent(in):: delsx, delsy, slat, slon, xi, xj, xlat, xlon
      real, dimension(im,jm), intent(in):: rlat, rlon
      real, dimension(im,jm), intent(out):: fi, fj
      !
      real :: slat1, slon1, xlat1, xlon1, x0, y0, xp, yp
      real :: rl0, ac, acn, ck, rck, t1, fxlon, rrlon, xslon
      real :: slata, slatb, slata1, slatb1, slat2
      real :: ali, di, dj, rlat1, rlon1, rdelsx, rdelsy, r0, cyclei
      real :: dcoslt(1400), dgw(1400), glat(1400), eps, xxxj
      integer :: i, j, k, jmg, jx0, jx
      integer, parameter :: is_debug = 0
      real, parameter :: pi = 3.141592653589793, pi4 = pi*0.25, poi = pi/180., rpoi = 180./pi
      !
      ! MERCATOR
      if (nprojc(1:3) == 'MER') then
	 slat1 = slat * poi
	 slon1 = slon * poi
	 xlat1 = xlat * poi
	 xlon1 = xlon * poi
	 ac = ra * cos(slat1)
	 x0 = 0.
	 y0 = ac * log((1. + sin(xlat1)) / cos(xlat1))
	 do j = 1, jm
	    do i = 1, im
	       rlat1 = rlat(i, j) * poi
	       rlon1 = rlon(i, j) * poi
	       di = rlon(i, j) - xlon
	       di = mod(di + 900., 360.) - 180.
	       di = ac * (di * poi)
	       dj = ac * log((1. + sin(rlat1)) / cos(rlat1))
	       fi(i, j) = (di - x0) / delsx + xi
	       fj(i, j) = (y0 - dj) / delsy + xj
	    end do
	 end do
      end if
      !
      ! LAMBERT NORTH
      if (nprojc(1:3) == 'LMN') then
	 slata = +30.
	 slatb = +60.
	 slata1 = +slata * poi
	 slatb1 = +slatb * poi
	 slat1 = pi4 - slata1 * 0.5
	 slat2 = pi4 - slatb1 * 0.5
	 slon1 = slon * poi
	 xlat1 = +xlat * poi
	 xlon1 = +xlon * poi
	 ck = log(cos(slata1) / cos(slatb1)) / log(tan(slat1) / tan(slat2))
	 rck = 1. / ck
	 acn = ra * cos(slata1) * rck
	 r0 = acn / (tan(slat1)) ** ck
	 rl0 = r0 * (tan(pi4 - xlat1 * 0.5)) ** ck
	 xslon = xlon - slon
	 xslon = mod(xslon + 900., 360.) - 180.
	 xslon = xslon * poi * ck
	 x0 = rl0 * sin(xslon)
	 y0 = rl0 * cos(xslon)
	 !
	 do j = 1, jm
            do i = 1, im
	       rlat1 = +rlat(i, j) * poi
	       rlon1 = +rlon(i, j) * poi
	       ali = r0 * (tan(pi4 - rlat1 * 0.5)) ** ck
	       fxlon = rlon(i, j) - xlon
	       fxlon = mod(fxlon + 900., 360.) - 180.
	       fxlon = fxlon * poi * ck
	       di = ali * sin(fxlon + xslon)
	       dj = ali * cos(fxlon + xslon)
	       fi(i, j) = (di - x0) / delsx + xi
	       fj(i, j) = (dj - y0) / delsy + xj
            end do
         end do
      endif
      !
      ! LATITUDE-LONGITUDE
      if (nprojc(1:2) == 'LL') then
	 rdelsx = 1. / delsx
	 rdelsy = 1. / delsy
	 cyclei = 360. / delsx
	 do j = 1, jm
	    do i = 1, im
	       fi(i, j) = xi + rdelsx * (rlon(i, j) - xlon)
	       fj(i, j) = xj - rdelsy * (rlat(i, j) - xlat)
	       if (fi(i, j) <  0.0) fi(i, j) = fi(i, j) + cyclei
	       if (fi(i, j) > cyclei) fi(i, j) = fi(i, j) - cyclei
	    end do
	 end do
      end if
      !
      ! Cylinderical Equidistant projection ---
      if (nprojc == 'CE  ') then
	 ac = ra * cos(slat * poi)
	 do j = 1, jm
	    do i = 1, im
	       di = rlon(i, j) - xlon
	       di = mod(di + 900., 360.) - 180.
	       di = ac * (di * poi)
	       dj = rlat(i, j) - xlat
	       dj = ra * (dj * poi)
	       fi(i, j) = xi + di / delsx
	       fj(i, j) = xj - dj / delsy
	    end do
	 end do
      end if
      !
      ! GAUSSIAN GRID
      if (nprojc(1:1) == 'G') then
	 jmg = int(slon + 0.01)
	 if (mod(jmg, 2) == 1 .or. jmg > 1400) then
	    write(6, *) 'DIMENSION OF GAUSSIAN LATITUDE IS INVALID'
	    write(6, *) 'JMG=', jmg
	    stop 90
	 end if
	 call gauss(dcoslt, dgw, jmg)
	 do j = 1, jmg
	    glat(j) = 90. - rpoi * acos(dcoslt(j))
	 end do
	 if (abs(xlat) < 0.01) then
	    jx0 = jmg / 2
	 else
	    write(6, *) ' ----- SORRY ! ------------------------- '
	    write(6, *) '  THIS PROGRAM ONLY SUPORT FOR XLAT=0.0  '
	    write(6, *) '        AT GASSIAN LATITUDE              '
	    write(6, *) ' --------------------------------------- '
	 end if
	 eps = 0.01
	 if (xj < 0.5) eps = -0.01
	 rdelsx = 1. / delsx
	 rdelsy = 1. / delsy
	 cyclei = 360. / delsx
	 !
	 ! FI
         do j = 1, jm
            do i = 1, im
               fi(i, j) = xi + rdelsx * (rlon(i, j) - xlon)
	       if (fi(i, j) <  0.0d0) fi(i, j) = fi(i, j) + cyclei
	       if (fi(i, j) > cyclei) fi(i, j) = fi(i, j) - cyclei
	    end do
         end do
         ! FJ
         do j = 1, jm
            do i = 1, im
	       xxxj = (glat(1) - rlat(i, j)) / (glat(1) - glat(jmg)) * (jmg - 1.) + 1
	       jx = int(xxxj)
	       if (xxxj <= 1.0) then
		  jx0 = 1
	       elseif(xxxj >= float(jmg) - 1.) then
		  jx0 = jmg - 1
	       else
		  if(rlat(i, j) >= glat(jx + 1) .and. rlat(i, j) <  glat(jx    )) then
		     jx0 = jx
		  end if
		  if(rlat(i, j) >= glat(jx    ) .and. rlat(i, j) <  glat(jx - 1)) then
		     jx0 = jx - 1
		  end if
		  if(rlat(i, j) >= glat(jx + 2) .and. rlat(i, j) <  glat(jx + 1)) then
		     jx0 = jx + 1
		  end if
	       end if
	       fj(i, j) = jx0 + (rlat(i, j) - glat(jx0)) / (glat(jx0 + 1) - glat(jx0))
	    end do
         end do
      endif
      !
      return
   end subroutine ll2ij
   !
   !
   !
   subroutine setorg(npro, dels, slat, slon, xi, xj, xlat, xlon, flat, flon, nx, ny, flatu, flonu, flatv, flonv)
      implicit none
      character(len = 4), intent(in) :: npro
      integer(4), intent(in):: nx, ny
      real, intent(in) :: dels, slat, slon, xi, xj, xlat, xlon
      real, dimension(nx,ny), intent(in ) :: flat, flon
      real, dimension(nx,ny), intent(out) :: flatu, flonu, flatv, flonv
      real, dimension(nx,ny) :: fiu, fju, fiv, fjv
      integer :: ix, jy
      !
      do jy = 1, ny
	 do ix = 1, nx
	    fiu(ix, jy) = real(2 * ix - 1) / 2.
	    fju(ix, jy) = real(jy)
	    fiv(ix, jy) = real(ix)
	    fjv(ix, jy) = real(2 * jy - 1) / 2.
	 end do
      end do
      fju(:, :) = float(ny + 1) - fju(:, :)
      fjv(:, :) = float(ny + 1) - fjv(:, :)
      call ij2ll(fiu, fju, flatu, flonu, nx, ny, npro, dels, dels, slat, slon, xi, xj, xlat, xlon)
      flonu(1:nx, 1:ny) = mod(flonu(1:nx, 1:ny) + 720., 360.)
      call ij2ll(fiv, fjv, flatv, flonv, nx, ny, npro, dels, dels, slat, slon, xi, xj, xlat, xlon)
      flonv(1:nx, 1:ny) = mod(flonv(1:nx, 1:ny) + 720., 360.)
      !
      return
   end subroutine setorg
   !
   !
   !
   subroutine cfmap(nprojc, stdlat, flati, lx, ly, fmap)
      implicit none
      character(len=4), intent(in) :: nprojc
      integer, intent(in):: lx, ly
      real, intent(in) :: stdlat
      real, dimension(lx,ly), intent(in) :: flati
      real, intent(out) :: fmap(lx,ly,6)
      real, parameter :: pi = 3.14159265358979, rearth = 6371000.
      integer :: imap
      real :: rad, fhms, pi4, slata, slatb, slata1, slatb1, slat1, slat2, ck, cns
      !
      rad = pi/180.d0
      if (stdlat >= 0.) then
         fhms = +1.d0
      else
         fhms = -1.d0
      end if
      !
      if ((nprojc == 'PSN ') .or. (nprojc == 'PSS ')) then
         do imap = 1, 3
            fmap(1:lx,1:ly,imap) = (1.+fhms*sin(stdlat*rad))/(1.d0+fhms*sin(flati(1:lx,1:ly)*rad))
         end do
      else if ((nprojc == 'LMN ') .or. (nprojc == 'LMS ')) then
	 pi4 = 0.25d0*pi
	 slata = 30.d0
	 slatb = 60.d0
	 slata1 = slata*rad
	 slatb1 = slatb*rad
	 slat1 = pi4-slata1*0.5d0
	 slat2 = pi4-slatb1*0.5d0
	 ck = log(cos(slata1)/cos(slatb1))/log(tan(slat1)/tan(slat2))
	 cns = cos(slata1)/(tan(slat1))**ck
	 do imap = 1, 3
	   fmap(1:lx,1:ly,imap) = cns*(tan(pi4-fhms*flati(1:lx,1:ly)*rad*0.5))**ck/cos(flati(1:lx,1:ly)*rad)
	 end do
      else if (nprojc == 'MER ') then
	 do imap = 1, 3
	    fmap(1:lx,1:ly,imap) = cos(stdlat*rad)/cos(flati(1:lx,1:ly)*rad)
	 end do
      else if (nprojc == 'DES ') then
         fmap(1:lx, 1:ly, 1:3) = 1.
      else if (nprojc == 'CE  ') then 
         fmap(1:lx, 1:ly,   1) = cos(stdlat * rad) / cos(flati(1:lx, 1:ly) * rad)
         fmap(1:lx, 1:ly, 2:3) = 1.
      else if (nprojc == 'LL  ') then 
         fmap(1:lx, 1:ly, 1) = 1. / (rearth * cos(flati(1:lx, 1:ly) * rad))
         fmap(1:lx, 1:ly, 2) = 1. / rearth
         fmap(1:lx, 1:ly, 3) = 1.
      end if
      fmap(1:lx, 1:ly, 4:6) = 1. / fmap(1:lx, 1:ly, 1:3)
      !
      return
   end subroutine cfmap
   !
   !
   !
   subroutine vrgdis(ii1, ii2, imin, imax, nd, msw, dx, dxl, dxr, vdx, vrdx, vrdx2)
      implicit none
      integer, intent(in) :: ii1, ii2, imin, imax, nd, msw
      real, intent(in) :: dx, dxl, dxr
      real, dimension(nd), intent(out) :: vdx, vrdx, vrdx2
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
      real, dimension(nz), intent(in) :: vrdz2
      real, dimension(nz), intent(out) :: zrp
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
      real, dimension(nz), intent(in) :: vrdz
      real, dimension(nz), intent(out) :: zrw
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
      real, dimension(nz), intent(in) :: zeta
      real, dimension(nz), intent(out) :: fp, dfdz
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
   subroutine setfifj(npro, nxo, nyo, delsx, delsy, slat, slon, xi, xj, xlat, xlon, &
                      flatu, flonu, flatv, flonv, flat, flon, nx, ny, fi, fj, fiu, fju, fiv, fjv)
      implicit none
      character(len=4), intent(in):: npro
      integer, intent(in):: nx, ny, nxo, nyo
      real, intent(in):: delsx, delsy, slat, slon, xi, xj, xlat, xlon
      real, dimension(nx,ny), intent(in):: flat, flon, flatu, flonu, flatv, flonv
      real, dimension(nx,ny), intent(out):: fi, fj, fiu, fju, fiv, fjv
      !
      call ll2ij(fi , fj , flat , flon , nx, ny, npro, delsx, delsy, slat, slon, xi, xj, xlat, xlon)
      call ll2ij(fiu, fju, flatu, flonu, nx, ny, npro, delsx, delsy, slat, slon, xi, xj, xlat, xlon)
      call ll2ij(fiv, fjv, flatv, flonv, nx, ny, npro, delsx, delsy, slat, slon, xi, xj, xlat, xlon)
      fj (:, :) = float(nyo + 1) - fj (:, :)
      fju(:, :) = float(nyo + 1) - fju(:, :)
      fjv(:, :) = float(nyo + 1) - fjv(:, :)
      if (.false.) then
	 where (fi <= 1)
	    fi = 1 + 1e-4
	 end where
	 where (fi >= nxo)
	    fi = nxo - 1e-4
	 end where
	 where (fiu <= 1)
	    fiu = 1 + 1e-4
	 end where
	 where (fiu >= nxo)
	    fiu = nxo - 1e-4
	 end where
	 where (fiv <= 1)
	    fiv = 1 + 1e-4
	 end where
	 where (fiv >= nxo)
	    fiv = nxo - 1e-4
	 end where
	 ! j
	 where (fj <= 1)
	    fj = 1 + 1e-4
	 end where
	 where (fj >= nyo)
	    fj = nyo - 1e-4
	 end where
	 where (fju <= 1)
	    fju = 1 + 1e-4
	 end where
	 where (fju >= nyo)
	    fju = nyo - 1e-4
	 end where
	 where (fjv <= 1)
	    fjv = 1 + 1e-4
	 end where
	 where (fjv >= nyo)
	    fjv = nyo - 1e-4
	 end where
      end if
      !
      return
   end subroutine setfifj
   !
   !
   !
   subroutine average2D(dx0, dy0, fi, fj, u0, nx0, ny0, nx, ny, u)
      implicit none
      integer, intent(in) :: nx0, ny0, nx, ny
      real, intent(in) :: dx0, dy0
      real, dimension(nx0,ny0), intent(in) :: fi, fj
      real, dimension(nx0,ny0), intent(in) :: u0
      real, dimension(nx,ny), intent(out) :: u
      !
      integer :: i, j, i0, j0
      real :: xb, yb
      real, dimension(0:nx+1) :: xe
      real, dimension(0:ny+1) :: ye
      real, dimension(0:nx+1,0:ny+1) :: weight, total, umean
      real, dimension(nx0,ny0) :: fil, fiu, fjl, fju
      !
      do i = 0, nx+1
         xe(i) = i*1.d0
      end do
      do j = 0, ny+1
         ye(j) = j*1.d0
      end do
      umean(:,:) = 0.d0
      total(:,:) = 0.d0
      fil(:,:) = fi(:,:) - 0.5*dx0
      fiu(:,:) = fi(:,:) + 0.5*dx0
      fjl(:,:) = fj(:,:) - 0.5*dy0
      fju(:,:) = fj(:,:) + 0.5*dy0
      !
      do i0 = 1, nx0
	 do j0 = 1, ny0
	    if (fi(i0,j0) <= xe(0) .or. fi(i0,j0) >= xe(nx+1) .or. &
	      & fj(i0,j0) <= ye(0) .or. fj(i0,j0) >= ye(ny+1)) cycle
	    i = ifix(fi(i0,j0))
	    if (i == nx+1) i = i - 1
	    j = ifix(fj(i0,j0))
	    if (j == ny+1) j = j - 1
	    weight(i,j) = 1.d0
	    weight(i+1,j) = 1.d0
	    weight(i+1,j+1) = 1.d0
	    weight(i,j+1) = 1.d0
	    !
	    xb = i + 0.5
	    if (fiu(i0,j0) <= xb) then
	       weight(i,j) = weight(i,j)*dx0
	       weight(i+1,j) = weight(i+1,j)*0.d0
	       weight(i+1,j+1) = weight(i+1,j+1)*0.d0
	       weight(i,j+1) = weight(i,j+1)*dx0
	    else if (fil(i0,j0) >= xb) then
	       weight(i,j) = weight(i,j)*0.d0
	       weight(i+1,j) = weight(i+1,j)*dx0
	       weight(i+1,j+1) = weight(i+1,j+1)*dx0
	       weight(i,j+1) = weight(i,j+1)*0.d0
	    else
	       weight(i,j) = weight(i,j)*(xb-fil(i0,j0))
	       weight(i+1,j) = weight(i+1,j)*(fiu(i0,j0)-xb)
	       weight(i+1,j+1) = weight(i+1,j+1)*(fiu(i0,j0)-xb)
	       weight(i,j+1) = weight(i,j+1)*(xb-fil(i0,j0))
	    end if
	    !
	    yb = j + 0.5
	    if (fju(i0,j0) <= yb) then
	       weight(i,j) = weight(i,j)*dy0
	       weight(i+1,j) = weight(i+1,j)*dy0
	       weight(i+1,j+1) = weight(i+1,j+1)*0.d0
	       weight(i,j+1) = weight(i,j+1)*0.d0
	    else if (fjl(i0,j0) >= yb) then
	       weight(i,j) = weight(i,j)*0.d0
	       weight(i+1,j) = weight(i+1,j)*0.d0
	       weight(i+1,j+1) = weight(i+1,j+1)*dy0
	       weight(i,j+1) = weight(i,j+1)*dy0
	    else
	       weight(i,j) = weight(i,j)*(yb-fjl(i0,j0))
	       weight(i+1,j) = weight(i+1,j)*(yb-fjl(i0,j0))
	       weight(i+1,j+1) = weight(i+1,j+1)*(fju(i0,j0)-yb)
	       weight(i,j+1) = weight(i,j+1)*(fju(i0,j0)-yb)
	    end if
	    !
	    umean(i,j) = umean(i,j) + weight(i,j)*u0(i0,j0)
	    umean(i+1,j) = umean(i+1,j) + weight(i+1,j)*u0(i0,j0)
	    umean(i+1,j+1) = umean(i+1,j+1) + weight(i+1,j+1)*u0(i0,j0)
	    umean(i,j+1) = umean(i,j+1) + weight(i,j+1)*u0(i0,j0)
	    total(i,j) = total(i,j) + weight(i,j)
	    total(i+1,j) = total(i+1,j) + weight(i+1,j)
	    total(i+1,j+1) = total(i+1,j+1) + weight(i+1,j+1)
	    total(i,j+1) = total(i,j+1) + weight(i,j+1)
	    !if (i == 100 .and. j == 100) print*, fi(i0,j0), fj(i0,j0), weight(i,j)
	    !if (i+1 == 100 .and. j == 100) print*, fi(i0,j0), fj(i0,j0), weight(i+1,j)
	    !if (i+1 == 100 .and. j+1 == 100) print*, fi(i0,j0), fj(i0,j0), weight(i+1,j+1)
	    !if (i == 100 .and. j+1 == 100) print*, fi(i0,j0), fj(i0,j0), weight(i,j+1)
	 end do
      end do
      u(:,:) = umean(1:nx,1:ny)/total(1:nx,1:ny)
      !
      return
   end subroutine average2D
   !
   !
   !
   subroutine average3D0(hsdiff, dx0, dy0, zrp0, vctrans0, zrp, vctrans, hs0, hs, fi, fj, u0, nx0, ny0, nz0, nx, ny, nz, u)
      implicit none
      integer, intent(in) :: hsdiff, nx0, ny0, nz0, nx, ny, nz
      real, intent(in) :: dx0, dy0
      real, dimension(nz), intent(in) :: zrp, vctrans
      real, dimension(nx,ny), intent(in) :: hs
      real, dimension(nz0), intent(in) :: zrp0, vctrans0
      real, dimension(nx0,ny0), intent(in) :: hs0, fi, fj
      real, dimension(nx0,ny0,nz0), intent(in) :: u0
      real, dimension(nx,ny,nz), intent(out) :: u
      !
      integer :: i, j, k, k0, ks, ke
      real :: weight
      real, dimension(0:nz) :: hb
      real, dimension(0:nz0) :: hintb
      real, dimension(nz) :: total, umean
      real, dimension(nx,ny) :: hsint
      real, dimension(nx,ny,nz0) :: uint, hint
      real, dimension(nx,ny,nz) :: h
      !
      call average2D(dx0, dy0, fi, fj, hs0, nx0, ny0, nx, ny, hsint)
      do k = 1, nz0
         call average2D(dx0, dy0, fi, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 hint(:,:,k) = zrp0(k) + hsint(:,:)*vctrans0(k)
	 if (hsdiff == 0) hint(:,:,k) = hint(:,:,k) - hsint(:,:)
      end do
      do k = 1, nz
         h(:,:,k) = zrp(k) + hs(:,:)*vctrans(k)
	 if (hsdiff == 0) h(:,:,k) = h(:,:,k) - hs(:,:)
      end do
      !
      do i = 1, nx
         do j = 1, ny
	    umean(:) = 0.d0
            total(:) = 0.d0
	    hb(0) = 1.5*h(i,j,1) - 0.5*h(i,j,2)
	    hb(1:nz-1) = 0.5*(h(i,j,1:nz-1)+h(i,j,2:nz))
	    hb(nz) = 1.5*h(i,j,nz) - 0.5*h(i,j,nz-1)
	    hintb(0) = 1.5*hint(i,j,1) - 0.5*hint(i,j,2)
	    hintb(1:nz0-1) = 0.5*(hint(i,j,1:nz0-1)+hint(i,j,2:nz0))
	    hintb(nz0) = 1.5*hint(i,j,nz0) - 0.5*hint(i,j,nz0-1)
	    do k0 = 1, nz0
	       if (hintb(k0) <= hb(0) .or. hintb(k0-1) >= hb(nz)) cycle
	       do k = 1, nz
	          if (hb(k) > hintb(k0-1)) exit
	       end do
	       ks = k
	       do k = 1, nz
	          if (hb(k) >= hintb(k0)) exit
	       end do
	       ke = k
	       do k = ks, ke
	          weight = min(hb(k),hintb(k0)) - max(hb(k-1),hintb(k0-1))
		  umean(k) = umean(k) + weight*uint(i,j,k0)
		  total(k) = total(k) + weight
		  !if (i == 100 .and. j == 100 .and. k0 == 10) print*, k, hb(k-1), hb(k), weight, hintb(k0-1), hintb(k0)
	       end do
	    end do
	    !
	    do k = 1, nz
	       if (total(k) > 0.) exit
	    end do
	    ks = k
	    do k = nz, 1, -1
	       if (total(k) > 0.) exit
	    end do
	    ke = k
	    u(i,j,ks:ke) = umean(ks:ke)/total(ks:ke)
	    do k = 1, ks-1
	       u(i,j,k) = u(i,j,ks)
	    end do
	    do k = ke+1, nz
	       u(i,j,k) = u(i,j,ke)
	    end do
         end do
      end do
      !
      return
   end subroutine average3D0
   !
   !
   !
   subroutine average3D(hsdiff, dx0, dy0, zrp0, vctrans0, zrp, vctrans, hs0, hs, fi, fj, u0, nx0, ny0, nz0, nx, ny, nz, u)
      implicit none
      integer, intent(in) :: hsdiff, nx0, ny0, nz0, nx, ny, nz
      real, intent(in) :: dx0, dy0
      real, dimension(nz), intent(in) :: zrp, vctrans
      real, dimension(nx,ny), intent(in) :: hs
      real, dimension(nz0), intent(in) :: zrp0, vctrans0
      real, dimension(nx0,ny0), intent(in) :: hs0, fi, fj
      real, dimension(nx0,ny0,nz0), intent(in) :: u0
      real, dimension(nx,ny,nz), intent(out) :: u
      !
      integer :: i, j, k, kk
      real :: weight
      real, dimension(nx,ny) :: hsint
      real, dimension(nx,ny,nz0) :: uint, hint
      real, dimension(nx,ny,nz) :: h
      !
      call average2D(dx0, dy0, fi, fj, hs0, nx0, ny0, nx, ny, hsint)
      do k = 1, nz0
         call average2D(dx0, dy0, fi, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 hint(:,:,k) = zrp0(k) + hsint(:,:)*vctrans0(k)
	 if (hsdiff == 0) hint(:,:,k) = hint(:,:,k) - hsint(:,:)
      end do
      do k = 1, nz
         h(:,:,k) = zrp(k) + hs(:,:)*vctrans(k)
	 if (hsdiff == 0) h(:,:,k) = h(:,:,k) - hs(:,:)
      end do
      !
      do i = 1, nx
         do j = 1, ny
	    do k = 1, nz
	       if (h(i,j,k) < hint(i,j,1)) then
		  u(i,j,k) = uint(i,j,1)
	       else if (h(i,j,k) > hint(i,j,nz0)) then
		  u(i,j,k) = uint(i,j,nz0)
	       else
		  do kk = 2, nz0
		     if (h(i,j,k) < hint(i,j,kk)) exit
		  end do
		  if (kk > nz0) kk = nz0
		  kk = kk - 1
		  weight = (h(i,j,k)-hint(i,j,kk))/(hint(i,j,kk+1)-hint(i,j,kk))
		  u(i,j,k) = (1.-weight)*uint(i,j,kk) + weight*uint(i,j,kk+1)
	       end if
	    end do
         end do
      end do
      !
      return
   end subroutine average3D
   !
   !
   !
   subroutine average3Duv0(hsdiff, dx0, dy0, stag0, stag, zrp0, vctrans0, dvtrans0, zrp, vctrans, dvtrans, hs0, hs, fi, fj, u0, nx0, ny0, nz0, nx, ny, nz, u)
      implicit none
      integer, intent(in) :: hsdiff, stag0, stag, nx0, ny0, nz0, nx, ny, nz
      real, intent(in) :: dx0, dy0
      real, dimension(nz), intent(in) :: zrp, vctrans, dvtrans
      real, dimension(nx,ny), intent(in) :: hs
      real, dimension(nz0), intent(in) :: zrp0, vctrans0, dvtrans0
      real, dimension(nx0,ny0), intent(in) :: hs0, fi, fj ! Note: fi, fj w.r.t mass points
      real, dimension(nx0,ny0,nz0), intent(in) :: u0
      real, dimension(nx,ny,nz), intent(out) :: u
      !
      integer :: i, j, k, k0, ks, ke
      real :: weight
      real, dimension(0:nz) :: hb
      real, dimension(0:nz0) :: hintb
      real, dimension(nz) :: total, umean
      real, dimension(nx,ny) :: hsint, hstmp
      real, dimension(nx,ny,nz0) :: uint, hint, g2int
      real, dimension(nx,ny,nz) :: h, g2
      !
      call average2D(dx0, dy0, fi, fj, hs0, nx0, ny0, nx, ny, hsint)
      do k = 1, nz0
         if (stag0 == 1) then
	    call average2D(dx0, dy0, fi+0.5, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 else if (stag0 == 2) then
	    call average2D(dx0, dy0, fi, fj+0.5, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 else
            call average2D(dx0, dy0, fi, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 end if
	 hint(:,:,k) = zrp0(k) + hsint(:,:)*vctrans0(k)
	 if (hsdiff == 0) hint(:,:,k) = hint(:,:,k) - hsint(:,:)
	 g2int(:,:,k) = 1. + hsint(:,:)*dvtrans0(k)
      end do
      !
      if (stag == 1) then
	 do i = 2, nx
	    hstmp(i,:) = 0.5*(hs(i-1,:)+hs(i,:))
	 end do
	 hstmp(1,:) = 2*hstmp(2,:) - hstmp(3,:)
      else if (stag == 2) then
	 do j = 2, ny
	    hstmp(:,j) = 0.5*(hs(:,j-1)+hs(:,j))
	 end do
	 hstmp(:,1) = 2*hstmp(:,2) - hstmp(:,3)
      else
	 hstmp(:,:) = hs
      end if
      do k = 1, nz
         h(:,:,k) = zrp(k) + hstmp(:,:)*vctrans(k)
	 if (hsdiff == 0) h(:,:,k) = h(:,:,k) - hstmp(:,:)
	 g2(:,:,k) = 1. + hstmp(:,:)*dvtrans(k)
      end do
      !
      do i = 1, nx
         do j = 1, ny
	    umean(:) = 0.d0
            total(:) = 0.d0
	    hb(0) = 1.5*h(i,j,1) - 0.5*h(i,j,2)
	    hb(1:nz-1) = 0.5*(h(i,j,1:nz-1)+h(i,j,2:nz))
	    hb(nz) = 1.5*h(i,j,nz) - 0.5*h(i,j,nz-1)
	    hintb(0) = 1.5*hint(i,j,1) - 0.5*hint(i,j,2)
	    hintb(1:nz0-1) = 0.5*(hint(i,j,1:nz0-1)+hint(i,j,2:nz0))
	    hintb(nz0) = 1.5*hint(i,j,nz0) - 0.5*hint(i,j,nz0-1)
	    do k0 = 1, nz0
	       if (hintb(k0) <= hb(0) .or. hintb(k0-1) >= hb(nz)) cycle
	       do k = 1, nz
	          if (hb(k) > hintb(k0-1)) exit
	       end do
	       ks = k
	       do k = 1, nz
	          if (hb(k) >= hintb(k0)) exit
	       end do
	       ke = k
	       do k = ks, ke
	          weight = min(hb(k),hintb(k0)) - max(hb(k-1),hintb(k0-1))
		  umean(k) = umean(k) + weight*uint(i,j,k0)/g2int(i,j,k0)
		  total(k) = total(k) + weight
	       end do
	    end do
	    u(i,j,:) = umean/total*g2(i,j,:)
         end do
      end do
      !if (stag == 1) then
	 !u(1,:,:) = 2*u(2,:,:) - u(3,:,:)
      !else if (stag == 2) then
	 !u(:,1,:) = 2*u(:,2,:) - u(:,3,:)
      !end if
      !
      return
   end subroutine average3Duv0
   !
   !
   !
   subroutine average3Duv(hsdiff, dx0, dy0, stag0, stag, zrp0, vctrans0, dvtrans0, zrp, vctrans, dvtrans, hs0, hs, fi, fj, u0, nx0, ny0, nz0, nx, ny, nz, u)
      implicit none
      integer, intent(in) :: hsdiff, stag0, stag, nx0, ny0, nz0, nx, ny, nz
      real, intent(in) :: dx0, dy0
      real, dimension(nz), intent(in) :: zrp, vctrans, dvtrans
      real, dimension(nx,ny), intent(in) :: hs
      real, dimension(nz0), intent(in) :: zrp0, vctrans0, dvtrans0
      real, dimension(nx0,ny0), intent(in) :: hs0, fi, fj ! Note: fi, fj w.r.t mass points
      real, dimension(nx0,ny0,nz0), intent(in) :: u0
      real, dimension(nx,ny,nz), intent(out) :: u
      !
      integer :: i, j, k, kk
      real :: weight
      real, dimension(nx,ny) :: hsint, hstmp
      real, dimension(nx,ny,nz0) :: uint, hint, g2int
      real, dimension(nx,ny,nz) :: h, g2
      !
      call average2D(dx0, dy0, fi, fj, hs0, nx0, ny0, nx, ny, hsint)
      do k = 1, nz0
         if (stag0 == 1) then
	    call average2D(dx0, dy0, fi+0.5, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 else if (stag0 == 2) then
	    call average2D(dx0, dy0, fi, fj+0.5, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 else
            call average2D(dx0, dy0, fi, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 end if
	 hint(:,:,k) = zrp0(k) + hsint(:,:)*vctrans0(k)
	 if (hsdiff == 0) hint(:,:,k) = hint(:,:,k) - hsint(:,:)
	 g2int(:,:,k) = 1. + hsint(:,:)*dvtrans0(k)
      end do
      !
      if (stag == 1) then
	 do i = 2, nx
	    hstmp(i,:) = 0.5*(hs(i-1,:)+hs(i,:))
	 end do
	 hstmp(1,:) = 2*hstmp(2,:) - hstmp(3,:)
      else if (stag == 2) then
	 do j = 2, ny
	    hstmp(:,j) = 0.5*(hs(:,j-1)+hs(:,j))
	 end do
	 hstmp(:,1) = 2*hstmp(:,2) - hstmp(:,3)
      else
	 hstmp(:,:) = hs
      end if
      do k = 1, nz
         h(:,:,k) = zrp(k) + hstmp(:,:)*vctrans(k)
	 if (hsdiff == 0) h(:,:,k) = h(:,:,k) - hstmp(:,:)
	 g2(:,:,k) = 1. + hstmp(:,:)*dvtrans(k)
      end do
      !
      do i = 1, nx
         do j = 1, ny
	    do k = 1, nz
	       if (h(i,j,k) < hint(i,j,1)) then
		  u(i,j,k) = uint(i,j,1)
	       else if (h(i,j,k) > hint(i,j,nz0)) then
		  u(i,j,k) = uint(i,j,nz0)
	       else
		  do kk = 2, nz0
		     if (h(i,j,k) < hint(i,j,kk)) exit
		  end do
		  if (kk > nz0) kk = nz0
		  kk = kk - 1
		  weight = (h(i,j,k)-hint(i,j,kk))/(hint(i,j,kk+1)-hint(i,j,kk))
		  u(i,j,k) = (1.-weight)*uint(i,j,kk)/g2int(i,j,kk) + weight*uint(i,j,kk+1)/g2int(i,j,kk+1)
		  u(i,j,k) = u(i,j,k)*g2(i,j,k)
	       end if
	    end do
         end do
      end do
      !
      return
   end subroutine average3Duv
   !
   !
   !
   subroutine pspt2p(zrp, dvtrans, zs, ps, pt, qv, nx, ny, nz, p)
      use variable, only : rdvcp, cpdvrd, gvcp, presrf
      implicit none
      !
      integer, intent(in)  :: nx, ny, nz
      real, dimension(nz), intent(in)  :: zrp, dvtrans
      real, dimension(nx,ny), intent(in)  :: zs, ps
      real, dimension(nx,ny,nz), intent(in)  :: pt, qv
      real, dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real, dimension(nx,ny) :: virtpt, pai
      real, dimension(nx,ny,nz) :: g2w, deltaz
      !
      do k = 1, nz
	 g2w(:,:,k) = 1. + dvtrans(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      pai = (ps/presrf)**rdvcp
      virtpt = 0.5*(pt(:,:,1)*(1.+0.608*qv(:,:,1)) + pt(:,:,2)*(1.+0.608*qv(:,:,2)))
      pai = pai - gvcp/virtpt*deltaz(:,:,1)
      p(:,:,1) = presrf*pai(:,:)**cpdvrd
      do k = 2, nz
	 virtpt = 0.5*(pt(:,:,k)*(1.+0.608*qv(:,:,k)) + pt(:,:,k-1)*(1.+0.608*qv(:,:,k-1)))
	 pai = pai - gvcp/virtpt*deltaz(:,:,k)
	 p(:,:,k) = presrf*pai(:,:)**cpdvrd
      end do
      !
      return
   end subroutine pspt2p
   !
   !
   !
   subroutine pbpt2p(zrp, dvtrans, zs, pb, pt, qv, nx, ny, nz, p)
      use variable, only : rdvcp, cpdvrd, gvcp, presrf
      implicit none
      !
      integer, intent(in)  :: nx, ny, nz
      real, dimension(nz), intent(in)  :: zrp, dvtrans
      real, dimension(nx,ny), intent(in)  :: zs, pb
      real, dimension(nx,ny,nz), intent(in)  :: pt, qv
      real, dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real, dimension(nx,ny) :: virtpt, pai
      real, dimension(nx,ny,nz) :: g2w, deltaz
      !
      do k = 1, nz
	 g2w(:,:,k) = 1. + dvtrans(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      p(:,:,1) = pb(:,:)
      pai = (pb/presrf)**rdvcp
      do k = 2, nz
	 virtpt = 0.5*(pt(:,:,k)*(1.+0.608*qv(:,:,k)) + pt(:,:,k-1)*(1.+0.608*qv(:,:,k-1)))
	 pai = pai - gvcp/virtpt*deltaz(:,:,k)
	 p(:,:,k) = presrf*pai(:,:)**cpdvrd
      end do
      !
      return
   end subroutine pbpt2p
   !
   !
   !
   subroutine t2es(t, nx, ny, nz, es)
      use variable, only : tkelvn, &
                         & tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         & tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,nz), intent(in) :: t
      real, dimension(nx,ny,nz), intent(out) :: es
      integer :: i, j, k
      real :: e0c, al, bl, e0i, ai, bi, tc
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
                  es(i,j,k) = e0c*exp(al*tc/(bl+tc))
               else if(tc <= tci) then
                  es(i,j,k) = e0i*exp(ai*tc/(bi+tc))
               else
                  es(i,j,k) = e0c*exp(al*tc/(bl+tc))*(tc-tci)/(tcw-tci) + &
                            & e0i*exp(ai*tc/(bi+tc))*(tcw-tc)/(tcw-tci)
               end if
            end do
         end do
      end do
      !
      return
   end subroutine t2es
   !
   !
   !
   subroutine qv2rh(p, t, qv, nx, ny, nz, rh)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,nz), intent(in) :: p, t, qv
      real, dimension(nx,ny,nz), intent(out) :: rh
      real, dimension(nx,ny,nz) :: es
      !
      call t2es(t, nx, ny, nz, es)
      where (es >= p)
         rh = 100*(qv*p/(0.622d0+0.378d0*qv))/p
      elsewhere
         rh = 100*(qv*p/(0.622d0+0.378d0*qv))/es
      end where
      where (rh < 0.)
         rh = 0.
      elsewhere (rh > 100.)
         rh = 100.
      end where
      !
      return
   end subroutine qv2rh
   !
   !
   !
   subroutine rh2qv(p, t, rh, nx, ny, nz, qv)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,nz), intent(in) :: p, t, rh
      real, dimension(nx,ny,nz), intent(out) :: qv
      real, dimension(nx,ny,nz) :: es, e
      !
      call t2es(t, nx, ny, nz, es)
      where (es >= p)
         e = 0.01*rh*p
      elsewhere
         e = 0.01*rh*es
      end where
      qv = 0.622d0*e/(p-0.378d0*e)
      !
      return
   end subroutine rh2qv
   !
   !
   !
   subroutine read_sst(fileunit, nx, ny, sst)
      implicit none
      !
      integer, intent(in) :: fileunit, nx, ny
      real, dimension(nx,ny), intent(out) :: sst
      integer :: istart, iend, jstart, jend
      !
      read(fileunit) istart, iend, jstart, jend, sst
      close(fileunit)
      !
      return
   end subroutine read_sst
   !
   !
   !
   subroutine calptop(horg, phi, t, plev, nx, ny, nz, ptop)
      use variable, only : g0, rd
      implicit none
      integer, intent(in)  :: nx, ny, nz
      real, intent(in)  :: horg
      real, dimension(nx,ny,nz), intent(in)  :: phi, t, plev
      real, dimension(nx,ny), intent(out) :: ptop
      integer :: ix, jy, kz
      integer, dimension(nx,ny) :: kz_top
      real :: zeta1, zeta2, gmt, plev2, ttop, tmean
      !
      kz_top(1:nx, 1:ny) = nz
      do kz = 2, nz
         do jy = 1, ny
            do ix = 1, nx
               zeta1 = phi(ix, jy, kz - 1)
               zeta2 = phi(ix, jy, kz)
               if (zeta1 < horg .and. zeta2 >= horg) then
                  kz_top(ix, jy) = kz
               end if
            end do
         end do
      end do
      !
      do jy = 1, ny
         do ix = 1, nx
            kz = kz_top(ix,jy)
            zeta1 = phi(ix,jy,kz-1)
            zeta2 = phi(ix,jy,kz)
            gmt = -(t(ix,jy,kz)-t(ix,jy,kz-1))/(zeta2-zeta1)
            plev2 = plev(ix,jy,kz)
            if (abs(gmt) > 0.001) then
               ttop = t(ix,jy,kz) + (zeta2-horg)*gmt
               ptop(ix,jy) = plev2*(ttop/t(ix,jy,kz))**(g0/rd/gmt)
            else
               tmean = t(ix,jy,kz) + 0.5*(zeta2-horg)*gmt
               ptop(ix, jy) = plev2*exp(g0*(zeta2-horg)/rd/tmean)
            end if
         end do
      end do
      !
      return
   end subroutine calptop
   !
   !
   !
   subroutine calpmean(dxi, dyi, pseai, prsgrd, ptopi, fmap, nx, ny, pseam, pgrdm, ptopm)
      implicit none
      integer, intent(in) :: nx, ny
      real, intent(in) :: dxi, dyi
      real, dimension(nx,ny), intent(in) :: pseai, prsgrd, ptopi
      real, dimension(nx,ny,6), intent(in) :: fmap
      real, intent(out) :: pseam, pgrdm, ptopm
      integer :: ix, jy
      real :: celarea, sumarea
      !
      pseam = 0.
      pgrdm = 0.
      ptopm = 0.
      sumarea = 0.
      do jy = 2, ny - 1
         do ix = 2, nx - 1
	    celarea = dxi * fmap(ix, jy, 4) * dyi * fmap(ix, jy, 5)
            pseam  = pseam + pseai (ix, jy) * celarea
            pgrdm  = pgrdm + prsgrd(ix, jy) * celarea
            ptopm  = ptopm + ptopi (ix, jy) * celarea
            sumarea = sumarea + celarea
         end do
      end do
      pseam = pseam / sumarea
      pgrdm = pgrdm / sumarea
      ptopm = ptopm / sumarea
      !
      return
   end subroutine calpmean
   !
   !
   !
   subroutine read_nhm(input_file, nx, ny, nz, u, v, w, t, pnh, qv, qc, qi, qr, qs, qg, tsoil)
      implicit none
      !
      character(*), intent(in) :: input_file
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,4), intent(out) :: tsoil
      real, dimension(nx,ny,nz), intent(out) :: u, v, w, t, pnh, qv, qc, qi, qr, qs, qg
      !
      open(80, file=trim(input_file), form='unformatted')
      read(80) u
      read(80) v
      read(80) w
      read(80) t
      read(80) pnh
      read(80) qv
      read(80) qc
      read(80) qi
      read(80) qr
      read(80) qs
      read(80) qg
      read(80) tsoil
      close(80)
      !
      return
   end subroutine read_nhm
   !
   !
   !
   subroutine read_nonhydro_icbc(input_file, readw, readq, readp, nx, ny, nz, pseam, ptopm, pgrdm, ptgrd, ptgrdt, sst, psurf, tin, u, v, w, pt, p, qv, qc, qi, qr, qs, qg)
      implicit none
      !
      character(*), intent(in) :: input_file
      integer, intent(in) :: readw, readq, readp, nx, ny, nz
      real, intent(out) :: pseam, ptopm, pgrdm
      real, dimension(nx,ny), intent(out) :: ptgrd, ptgrdt, sst, psurf
      real, dimension(nx,ny,4), intent(out) :: tin
      real, dimension(nx,ny,nz), intent(out) :: u, v, w, pt, p, qv, qc, qi, qr, qs, qg
      integer :: kt, mtuv, idate(5), flgq, flgqn, flgv, nx0, ny0, nz0, ngm, iz1, iz2, vctrans_type, n_vctrans, ix1, ix2, iy1, iy2
      real :: dtratio, dz, dzl, dzr, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr
      !
      open(80, file=trim(input_file), form='unformatted')
      read(80) kt, mtuv, pseam, ptopm, idate, pgrdm, &
      & flgq, flgqn, flgv, &
      & dtratio, nx0, ny0, nz0, ngm, dz, dzl, dzr, iz1, iz2, &
      & vctrans_type, zl_vctrans, zh_vctrans, n_vctrans, &
      & dx, dxl, dxr, ix1, ix2, dy, dyl, dyr, iy1, iy2
      read(80) u
      if (readw == 1) then
	 read(80) v,w
      else
         read(80) v
      endif
      read(80) pt
      if (readq == 1) then
	 read(80) qv, qc, qi, qr, qs, qg
      else
         read(80) qv
      endif
      !
      if (readp == 1) then
	 read(80) ptgrd,ptgrdt,sst,psurf,tin,p
      else
         read(80) ptgrd,ptgrdt,sst,psurf,tin
      endif
      close(80)
      !
      return
   end subroutine read_nonhydro_icbc
   !
   !
   !
   subroutine write_nhm(output_file, u, v, w, t, pnh, qv, qc, qi, qr, qs, qg, tsoil, nx, ny, nz)
      implicit none
      !
      character(*), intent(in) :: output_file
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,4), intent(in) :: tsoil
      real, dimension(nx,ny,nz), intent(in) :: u, v, w, t, pnh, qv, qc, qi, qr, qs, qg
      !
      open(90, file=trim(output_file), form='unformatted')
      write(90) u
      write(90) v
      write(90) w
      write(90) t
      write(90) pnh
      write(90) qv
      write(90) qc
      write(90) qi
      write(90) qr
      write(90) qs
      write(90) qg
      write(90) tsoil
      close(90)
      !
      return
   end subroutine write_nhm
   !
   !
   !
   subroutine write_ic(output_file, kt, idate, vctrans_type, n_vctrans, ix1, ix2, iy1, iy2, iz1, iz2, &
                     & pseam, ptopm, pgrdm, dtratio, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr, dz, dzl, dzr, &
		     & ptgrd, ptgrdt, sst, psurf, tin, u, v, w, pt, p, qv, qc, qi, qr, qs, qg, nx, ny, nz)
      implicit none
      !
      character(*), intent(in) :: output_file
      integer, intent(in) :: kt, nx, ny, nz, idate(5), vctrans_type, n_vctrans, ix1, ix2, iy1, iy2, iz1, iz2
      real, intent(in) :: pseam, ptopm, pgrdm, dtratio, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr, dz, dzl, dzr
      real, dimension(nx,ny), intent(in) :: ptgrd, ptgrdt, sst, psurf
      real, dimension(nx,ny,4), intent(in) :: tin
      real, dimension(nx,ny,nz), intent(in) :: u, v, w, pt, p, qv, qc, qi, qr, qs, qg
      integer :: mtuv, flgq, flgqn, flgv
      !
      mtuv = 22023
      flgq = 11111
      flgqn = 0
      flgv = 100001010
      !
      open(90, file=trim(output_file), form='unformatted')
      write(90) kt, mtuv, pseam, ptopm, idate, pgrdm, &
      & flgq, flgqn, flgv, &
      & dtratio, nx, ny, nz, 4, dz, dzl, dzr, iz1, iz2, &
      & vctrans_type, zl_vctrans, zh_vctrans, n_vctrans, &
      & dx, dxl, dxr, ix1, ix2, dy, dyl, dyr, iy1, iy2
      write(90) u
      write(90) v,w
      write(90) pt
      write(90) qv,qc,qi,qr,qs,qg
      write(90) ptgrd,ptgrdt,sst,psurf,tin,p
      close(90)
      !
      return
   end subroutine write_ic
   !
   !
   !
end module outer2innerlib