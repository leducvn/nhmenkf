module inner2outerlib
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
   subroutine nintr2d(fi, fj, zso, nxo, nyo, nx, ny, zsw)
      implicit none
      integer, intent(in) :: nxo, nyo, nx, ny
      real, dimension(nx,ny), intent(in) :: fi, fj
      real, dimension(nxo,nyo), intent(in) :: zso
      real, dimension(nx,ny), intent(out) :: zsw
      integer :: ix, jy
      real :: xx, yy, xo1, xo2, yo1, yo2
      integer, dimension(nx,ny) :: ixo, jyo
      real, dimension(nx,ny) :: t11, t21, t12, t22
      !
      do jy = 1, ny
	 do ix = 1, nx
	    ixo(ix, jy) = ifix(fi(ix, jy))
	    jyo(ix, jy) = ifix(fj(ix, jy))
	    ! Duc modified
	    if ((ixo(ix, jy) <= 0) .or. (ixo(ix, jy) >= nxo) .or. (jyo(ix, jy) <= 0) .or. (jyo(ix, jy) >= nyo)) then
	       write(6, '(1x, "THE REGION OF OUTER MODEL IS TOO SMALL")')
	       write(6, '(2 i4, 2 f9 .3)') ixo(ix,jy), jyo(ix,jy), fi(ix,jy), fj(ix,jy)
	       !stop 99
	    end if
	    if (ixo(ix,jy) <= 0) then
	       ixo(ix,jy) = 1
	    else if (ixo(ix,jy) >= nxo) then
	       ixo(ix,jy) = nxo - 1
	    end if
	    if (jyo(ix,jy) <= 0) then
	       jyo(ix,jy) = 1
	    else if (jyo(ix,jy) >= nyo) then
	       jyo(ix,jy) = nyo - 1
	    end if
	 end do
      end do
      !
      do jy = 1, ny
	 do ix = 1, nx
	    xx = fi(ix, jy)
	    yy = fj(ix, jy)
	    xo1 = float(ixo(ix, jy))
	    xo2 = float(ixo(ix, jy) + 1)
	    yo1 = float(jyo(ix, jy))
	    yo2 = float(jyo(ix, jy) + 1)
	    t11(ix, jy) = (xo2 - xx) * (yo2 - yy)
	    t21(ix, jy) = (xx - xo1) * (yo2 - yy)
	    t12(ix, jy) = (xo2 - xx) * (yy - yo1)
	    t22(ix, jy) = (xx - xo1) * (yy - yo1)
	 end do
      end do
      do jy = 1, ny
	 do ix = 1, nx
	    zsw(ix,jy) = zso(ixo(ix,jy),jyo(ix,jy))*t11(ix,jy) + zso(ixo(ix,jy)+1,jyo(ix,jy))*t21(ix,jy) &
		       + zso(ixo(ix,jy),jyo(ix,jy)+1)*t12(ix,jy) + zso(ixo(ix,jy)+1,jyo(ix,jy)+1)*t22(ix,jy)
	 end do
      end do
      !
      return
   end subroutine nintr2d
   !
   !
   !
   subroutine nintrp(iuv_swit, zrpo, zrp, fi, fj, pto, nxo, nyo, nzo, nx, ny, nz, pt)
      implicit none
      integer, intent(in) :: nxo, nyo, nzo, nx, ny, nz, iuv_swit
      real, dimension(nzo), intent(in) :: zrpo
      real, dimension(nz), intent(in) :: zrp
      real, dimension(nx, ny), intent(in) :: fi, fj
      real, dimension(nxo, nyo, nzo), intent(in) :: pto
      real, dimension(nx, ny, nz), intent(out) :: pt
      integer :: ix, jy, kz, ixo, jyo, kzo
      real :: xx, yy, zz, xo1, xo2, yo1, yo2, zo1, zo2, xo05, yo05
      real :: t111, t112, t121, t122, t211, t212, t221, t222
      !
      xo05 = 0.0
      yo05 = 0.0
      if (iuv_swit == 1) then
         xo05 = 0.5
         yo05 = 0.0
      elseif(iuv_swit == 2) then
         xo05 = 0.0
         yo05 = 0.5
      end if
      !
      kzo = 1
      do kz = 1, nz
         zz = zrp(kz)
         kzo = kzo - 1
1100     kzo = max0(kzo + 1, 1)
         if (zrpo(kzo) <= zz) then
            if(kzo < nzo) goto 1100
            if (zrpo(kzo) < zz) then
	       ! Duc modified
               write(6, '(1x, "WARNING: THE TOP OF OUTER MODEL IS TOO SHORT")')
               write(6, '(i3, f9.2, i3, f9.2)') kzo, zrpo(kzo), kz, zrp(kz)
               !stop 99
	       kzo = nzo + 1
            end if
         end if
         kzo = kzo - 1
         if (kzo > 0 .and. kzo < nzo) then
            zo1 = zrpo(kzo)
            zo2 = zrpo(kzo + 1)
            do jy = 1, ny
               do ix = 1, nx
		  ixo = ifix(fi(ix, jy) + xo05)
		  jyo = ifix(fj(ix, jy) + yo05)
		  if ((ixo <= 0) .or. (ixo >= nxo) .or. (jyo <= 0) .or. (jyo >= nyo)) then
		     write(6, '(1x, "WARNING: THE REGION OF OUTER MODEL IS TOO SMALL")')
		     write(6, '(4i4, 2f9.3)') ixo, jyo, ix, jy, fi(ix, jy), fj(ix, jy)
		     !stop 99
		  end if
		  if (ixo <= 0) then
		     ixo = 1
		     xx = 1. - xo05
		  else if (ixo >= nxo) then
		     ixo = nxo - 1
		     xx = nxo - xo05
		  else
		     xx = fi(ix, jy)
		  end if
		  if (jyo <= 0) then
		     jyo = 1
		     yy = 1. - yo05
		  else if (jyo >= nyo) then
		     jyo = nyo - 1
		     yy = nyo - yo05
		  else
		     yy = fj(ix, jy)
		  end if
		  xo1 = float(ixo) - xo05
		  xo2 = float(ixo + 1) - xo05
		  yo1 = float(jyo) - yo05
		  yo2 = float(jyo + 1) - yo05
		  t111 = (xo2 - xx) * (yo2 - yy) * (zo2 - zz)
		  t211 = (xx - xo1) * (yo2 - yy) * (zo2 - zz)
		  t121 = (xo2 - xx) * (yy - yo1) * (zo2 - zz)
		  t221 = (xx - xo1) * (yy - yo1) * (zo2 - zz)
		  t112 = (xo2 - xx) * (yo2 - yy) * (zz - zo1)
		  t212 = (xx - xo1) * (yo2 - yy) * (zz - zo1)
		  t122 = (xo2 - xx) * (yy - yo1) * (zz - zo1)
		  t222 = (xx - xo1) * (yy - yo1) * (zz - zo1)
		  pt(ix,jy,kz) = (pto(ixo,jyo,kzo)*t111 + pto(ixo+1,jyo,kzo)*t211 + pto(ixo,jyo+1,kzo)*t121 &
		                + pto(ixo+1,jyo+1,kzo)*t221 + pto(ixo,jyo,kzo+1)*t112 + pto(ixo+1,jyo,kzo+1)*t212 &
		                + pto(ixo,jyo+1,kzo+1)*t122 + pto(ixo+1,jyo+1,kzo+1)*t222) &
		               / (t111 + t211 + t121 + t221 + t112 + t212 + t122 + t222)
	       end do
	    end do
         else
	    if (kzo == 0) kzo = 1
            do jy = 1, ny
               do ix = 1, nx
		  ixo = ifix(fi(ix, jy) + xo05)
		  jyo = ifix(fj(ix, jy) + yo05)
		  if ((ixo <= 0) .or. (ixo >= nxo) .or. (jyo <= 0) .or. (jyo >= nyo)) then
		     write(6, '(1x, "WARNING: THE REGION OF OUTER MODEL IS TOO SMALL")')
		     write(6, '(4i4, 2f9.3)') ixo, jyo, ix, jy, fi(ix, jy), fj(ix, jy)
		     !stop 99
		  end if
		  if (ixo <= 0) then
		     ixo = 1
		     xx = 1. - xo05
		  else if (ixo >= nxo) then
		     ixo = nxo - 1
		     xx = nxo - xo05
		  else
		     xx = fi(ix, jy)
		  end if
		  if (jyo <= 0) then
		     jyo = 1
		     yy = 1. - yo05
		  else if (jyo >= nyo) then
		     jyo = nyo - 1
		     yy = nyo - yo05
		  else
		     yy = fj(ix, jy)
		  end if
		  xo1 = float(ixo) - xo05
		  xo2 = float(ixo + 1) - xo05
		  yo1 = float(jyo) - yo05
		  yo2 = float(jyo + 1) - yo05
		  t111 = (xo2 - xx) * (yo2 - yy)
		  t211 = (xx - xo1) * (yo2 - yy)
		  t121 = (xo2 - xx) * (yy - yo1)
		  t221 = (xx - xo1) * (yy - yo1)
		  pt(ix,jy,kz) = (pto(ixo,jyo,kzo)*t111 + pto(ixo+1,jyo,kzo)*t211 + pto(ixo,jyo+1,kzo)*t121 &
		                + pto(ixo+1,jyo+1,kzo)*t221)/(t111+t211+t121+t221)
               end do
            end do
         end if
      end do
      !
      return
   end subroutine nintrp
   !
   !
   !
   subroutine nintrpt(zrpo, vctrans_o, zrp, vctrans_i, zso, zsi, fi, fj, pto, nxo, nyo, nzo, nx, ny, nz, pt)
      implicit none
      integer, intent(in) :: nxo, nyo, nzo, nx, ny, nz
      real, dimension(nzo), intent(in) :: zrpo, vctrans_o
      real, dimension(nz), intent(in) :: zrp, vctrans_i
      real, dimension(nxo, nyo), intent(in) :: zso
      real, dimension(nx, ny), intent(in) :: zsi, fi, fj
      real, dimension(nxo, nyo, nzo), intent(in) :: pto
      real, dimension(nx,ny,nz), intent(out) :: pt
      !
      integer :: ix, jy, kz, ixo, jyo, ixq, jyq, kzq, kzq1
      real :: xx, yy, xo1, xo2, yo1, yo2, zzri, zzro, zzro1
      real :: ptint11, ptint21, ptint12, ptint22, t11, t21, t12, t22
      real, dimension(nx, ny) :: t11sav, t12sav, t21sav, t22sav, rtsumsav
      integer, dimension(nx, ny, nz, 4) :: kzqall
      real, dimension(nx, ny, nz, 4) :: dif1zrp, dif2zrp, rdifzrp
      !
      do jy = 1, ny
         do ix = 1, nx
            ixo = ifix(fi(ix, jy))
            jyo = ifix(fj(ix, jy))
            do kz = 1, nz
               zzri = zrp(kz) + zsi(ix, jy) * vctrans_i(kz)
	       ixq = ixo
	       jyq = jyo
	       do kzq = 1, nzo
		  zzro = zrpo(kzq) + zso(ixq, jyq) * vctrans_o(kzq)
		  if(zzro > zzri) exit
	       end do
	       kzqall(ix, jy, kz, 1) = kzq
	       !
	       ixq = ixo + 1
	       jyq = jyo
	       do kzq = 1, nzo
		  zzro = zrpo(kzq) + zso(ixq, jyq) * vctrans_o(kzq)
		  if(zzro > zzri) exit
	       end do
	       kzqall(ix, jy, kz, 2) = kzq
	       !
	       ixq = ixo
	       jyq = jyo + 1
	       do kzq = 1, nzo
		  zzro = zrpo(kzq) + zso(ixq, jyq) * vctrans_o(kzq)
		  if(zzro > zzri) exit
	       end do
	       kzqall(ix, jy, kz, 3) = kzq
	       !
	       ixq = ixo + 1
	       jyq = jyo + 1
	       do kzq = 1, nzo
		  zzro = zrpo(kzq) + zso(ixq, jyq) * vctrans_o(kzq)
		  if(zzro > zzri) exit
	       end do
	       kzqall(ix, jy, kz, 4) = kzq
            end do
         end do
      end do
      !
      do kz = 1, nz
         do jy = 1, ny
            do ix = 1, nx
	       ixo = ifix(fi(ix, jy))
	       jyo = ifix(fj(ix, jy))
	       zzri = zrp(kz) + zsi(ix, jy) * vctrans_i(kz)
	       ixq = ixo
	       jyq = jyo
	       kzq = kzqall(ix, jy, kz, 1)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  zzro  = zrpo(kzq ) + zso(ixq, jyq) * vctrans_o(kzq )
		  zzro1 = zrpo(kzq1) + zso(ixq, jyq) * vctrans_o(kzq1)
		  dif1zrp(ix, jy, kz, 1) = zzri - zzro1
		  dif2zrp(ix, jy, kz, 1) = zzro - zzri
		  rdifzrp(ix, jy, kz, 1) = 1.0 / (zzro - zzro1)
               else
                  dif1zrp(ix, jy, kz, 1) = 0.0
                  dif2zrp(ix, jy, kz, 1) = 0.0
                  rdifzrp(ix, jy, kz, 1) = 0.0
               end if
               !
	       ixq = ixo + 1
	       jyq = jyo
	       kzq = kzqall(ix, jy, kz, 2)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  zzro  = zrpo(kzq ) + zso(ixq, jyq) * vctrans_o(kzq )
		  zzro1 = zrpo(kzq1) + zso(ixq, jyq) * vctrans_o(kzq1)
		  dif1zrp(ix, jy, kz, 2) = zzri - zzro1
		  dif2zrp(ix, jy, kz, 2) = zzro - zzri
		  rdifzrp(ix, jy, kz, 2) = 1.0 / (zzro - zzro1)
	       else
		  dif1zrp(ix, jy, kz, 2) = 0.0
		  dif2zrp(ix, jy, kz, 2) = 0.0
		  rdifzrp(ix, jy, kz, 2) = 0.0
	       end if
	       !
	       ixq = ixo
	       jyq = jyo + 1
	       kzq = kzqall(ix, jy, kz, 3)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  zzro  = zrpo(kzq ) + zso(ixq, jyq) * vctrans_o(kzq )
		  zzro1 = zrpo(kzq1) + zso(ixq, jyq) * vctrans_o(kzq1)
		  dif1zrp(ix, jy, kz, 3) = zzri - zzro1
		  dif2zrp(ix, jy, kz, 3) = zzro - zzri
		  rdifzrp(ix, jy, kz, 3) = 1.0 / (zzro - zzro1)
	       else
		  dif1zrp(ix, jy, kz, 3) = 0.0
		  dif2zrp(ix, jy, kz, 3) = 0.0
		  rdifzrp(ix, jy, kz, 3) = 0.0
	       end if
	       !
	       ixq = ixo + 1
	       jyq = jyo + 1
	       kzq = kzqall(ix, jy, kz, 4)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  zzro  = zrpo(kzq ) + zso(ixq, jyq) * vctrans_o(kzq )
		  zzro1 = zrpo(kzq1) + zso(ixq, jyq) * vctrans_o(kzq1)
		  dif1zrp(ix, jy, kz, 4) = zzri - zzro1
		  dif2zrp(ix, jy, kz, 4) = zzro - zzri
		  rdifzrp(ix, jy, kz, 4) = 1.0 / (zzro - zzro1)
	       else
		  dif1zrp(ix, jy, kz, 4) = 0.0
		  dif2zrp(ix, jy, kz, 4) = 0.0
		  rdifzrp(ix, jy, kz, 4) = 0.0
	       end if
            end do
         end do
      end do
      !
      do jy = 1, ny
         do ix = 1, nx
            ixo = ifix(fi(ix, jy))
            jyo = ifix(fj(ix, jy))
            xx = fi(ix, jy)
            yy = fj(ix, jy)
            xo1 = float(ixo)
            xo2 = float(ixo + 1)
            yo1 = float(jyo)
            yo2 = float(jyo + 1)
            t11sav(ix, jy) = (xo2 - xx) * (yo2 - yy)
            t21sav(ix, jy) = (xx - xo1) * (yo2 - yy)
            t12sav(ix, jy) = (xo2 - xx) * (yy - yo1)
            t22sav(ix, jy) = (xx - xo1) * (yy - yo1)
            rtsumsav(ix, jy) = 1.0/(t11sav(ix,jy)+t21sav(ix,jy)+t12sav(ix,jy)+t22sav(ix,jy))
         end do
      end do
      !
      do kz = 1, nz
         do jy = 1, ny
            do ix = 1, nx
               ixo = ifix(fi(ix, jy))
               jyo = ifix(fj(ix, jy))
               ixq = ixo
               jyq = jyo
               kzq = kzqall(ix, jy, kz, 1)
               if (kzq > 1 .and. kzq < nzo+1) then
                  kzq1 = kzq - 1
                  ptint11 = (pto(ixq,jyq,kzq)*dif1zrp(ix,jy,kz,1)+pto(ixq,jyq,kzq1)*dif2zrp(ix,jy,kz,1))*rdifzrp(ix,jy,kz,1)
               else if (kzq == 1) then
                  ptint11 = pto(ixq, jyq, kzq)
	       else
		  ptint11 = pto(ixq, jyq, kzq-1)
               end if
	       !
	       ixq = ixo + 1
	       jyq = jyo
	       kzq = kzqall(ix, jy, kz, 2)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  ptint21 = (pto(ixq,jyq,kzq)*dif1zrp(ix,jy,kz,2)+pto(ixq,jyq,kzq1)*dif2zrp(ix,jy,kz,2))*rdifzrp(ix,jy,kz,2)
	       else if (kzq == 1) then
	          ptint21 = pto(ixq, jyq, kzq)
	       else
		  ptint21 = pto(ixq, jyq, kzq-1)
	       end if
	       !
	       ixq = ixo
	       jyq = jyo + 1
	       kzq = kzqall(ix, jy, kz, 3)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  ptint12 = (pto(ixq,jyq,kzq)*dif1zrp(ix,jy,kz,3)+pto(ixq,jyq,kzq1)*dif2zrp(ix,jy,kz,3))*rdifzrp(ix,jy,kz,3)
	       else if (kzq == 1) then
	          ptint12 = pto(ixq, jyq, kzq)
	       else
		  ptint12 = pto(ixq, jyq, kzq-1)
	       end if
	       !
	       ixq = ixo + 1
	       jyq = jyo + 1
	       kzq = kzqall(ix, jy, kz, 4)
	       if (kzq > 1 .and. kzq < nzo+1) then
		  kzq1 = kzq - 1
		  ptint22 = (pto(ixq,jyq,kzq)*dif1zrp(ix,jy,kz,4)+pto(ixq,jyq,kzq1)*dif2zrp(ix,jy,kz,4))*rdifzrp(ix,jy,kz,4)
	       else if (kzq == 1) then
	          ptint22 = pto(ixq, jyq, kzq)
	       else
		  ptint22 = pto(ixq, jyq, kzq-1)
	       end if
	       t11 = t11sav(ix, jy)
	       t21 = t21sav(ix, jy)
	       t12 = t12sav(ix, jy)
	       t22 = t22sav(ix, jy)
	       pt(ix, jy, kz) = (ptint11 * t11 + ptint21 * t21 + ptint12 * t12 + ptint22 * t22) * rtsumsav(ix, jy)
            end do
         end do
      end do
      !
      return
   end subroutine nintrpt
   !
   !
   !
   subroutine interpolate2D(fi, fj, u0, nx0, ny0, nx, ny, u)
      implicit none
      integer, intent(in) :: nx0, ny0, nx, ny
      real, dimension(nx,ny), intent(in) :: fi, fj
      real, dimension(nx0,ny0), intent(in) :: u0
      real, dimension(nx,ny), intent(out) :: u
      integer :: i, j, ii, jj
      real :: weight, wx, wy
      !
      do i = 1, nx
	 do j = 1, ny
	    ii = ifix(fi(i,j))
	    jj = ifix(fj(i,j))
	    if (fi(i,j) < 1. .and. fj(i,j) < 1.) then
	       u(i,j) = u0(1,1)
	    else if (fi(i,j) > 1.*nx0 .and. fj(i,j) < 1.) then
	       u(i,j) = u0(nx0,1)
	    else if (fi(i,j) > 1.*nx0 .and. fj(i,j) > 1.*ny0) then
	       u(i,j) = u0(nx0,ny0)
	    else if (fi(i,j) < 1. .and. fj(i,j) > 1.*ny0) then
	       u(i,j) = u0(1,ny0)
	    else if (fj(i,j) < 1.) then
	       if (ii >= nx0) ii = nx0 - 1
	       weight = fi(i,j) - ii
	       u(i,j) = (1.-weight)*u0(ii,1) + weight*u0(ii+1,1)
	    else if (fj(i,j) > 1.*ny0) then
	       if (ii >= nx0) ii = nx0 - 1
	       weight = fi(i,j) - ii
	       u(i,j) = (1.-weight)*u0(ii,ny0) + weight*u0(ii+1,ny0)
	    else if (fi(i,j) < 1.) then
	       if (jj >= ny0) jj = ny0 - 1
	       weight = fj(i,j) - jj
	       u(i,j) = (1.-weight)*u0(1,jj) + weight*u0(1,jj+1)
	    else if (fi(i,j) > 1.*nx0) then
	       if (jj >= ny0) jj = ny0 - 1
	       weight = fj(i,j) - jj
	       u(i,j) = (1.-weight)*u0(nx0,jj) + weight*u0(nx0,jj+1)
	    else
	       if (ii >= nx0) ii = nx0 - 1
	       if (jj >= ny0) jj = ny0 - 1
	       wx = fi(i,j) - ii
	       wy = fj(i,j) - jj
	       u(i,j) = (1.-wx)*(1.-wy)*u0(ii,jj) + wx*(1.-wy)*u0(ii+1,jj) + wx*wy*u0(ii+1,jj+1) + (1.-wx)*wy*u0(ii,jj+1)
	    end if
	 end do
      end do
      !
      return
   end subroutine interpolate2D
   !
   !
   !
   subroutine interpolate3D(hsdiff, zrp0, vctrans0, zrp, vctrans, hs0, hs, fi, fj, u0, nx0, ny0, nz0, nx, ny, nz, u)
      implicit none
      integer, intent(in) :: hsdiff, nx0, ny0, nz0, nx, ny, nz
      real, dimension(nz0), intent(in) :: zrp0, vctrans0
      real, dimension(nz), intent(in) :: zrp, vctrans
      real, dimension(nx0,ny0), intent(in) :: hs0
      real, dimension(nx,ny), intent(in) :: hs, fi, fj
      real, dimension(nx0,ny0,nz0), intent(in) :: u0
      real, dimension(nx,ny,nz), intent(out) :: u
      !
      integer :: i, j, k, kk
      real :: weight
      real, dimension(nx,ny) :: hsint
      real, dimension(nx,ny,nz0) :: uint, hint
      real, dimension(nx,ny,nz) :: h
      !
      call interpolate2D(fi, fj, hs0, nx0, ny0, nx, ny, hsint)
      do k = 1, nz0
         call interpolate2D(fi, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
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
   end subroutine interpolate3D
   !
   !
   !
   subroutine interpolate3Duv(hsdiff, stag0, stag, zrp0, vctrans0, dvtrans0, zrp, vctrans, dvtrans, hs0, hs, fi, fj, u0, nx0, ny0, nz0, nx, ny, nz, u)
      implicit none
      integer, intent(in) :: hsdiff, stag0, stag, nx0, ny0, nz0, nx, ny, nz
      real, dimension(nz0), intent(in) :: zrp0, vctrans0, dvtrans0
      real, dimension(nz), intent(in) :: zrp, vctrans, dvtrans
      real, dimension(nx0,ny0), intent(in) :: hs0
      real, dimension(nx,ny), intent(in) :: hs, fi, fj ! Note: fi, fj w.r.t mass points
      real, dimension(nx0,ny0,nz0), intent(in) :: u0
      real, dimension(nx,ny,nz), intent(out) :: u
      !
      integer :: i, j, k, kk
      real :: weight
      real, dimension(nx,ny) :: hsint, hstmp
      real, dimension(nx,ny,nz0) :: uint, hint, g2int
      real, dimension(nx,ny,nz) :: h, g2
      !
      call interpolate2D(fi, fj, hs0, nx0, ny0, nx, ny, hsint)
      do k = 1, nz0
         if (stag0 == 1) then
	    call interpolate2D(fi+0.5, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 else if (stag0 == 2) then
	    call interpolate2D(fi, fj+0.5, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
	 else
            call interpolate2D(fi, fj, u0(:,:,k), nx0, ny0, nx, ny, uint(:,:,k))
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
	 !hstmp(1,:) = hstmp(2,:)
	 hstmp(1,:) = 2*hstmp(2,:) - hstmp(3,:)
      else if (stag == 2) then
	 do j = 2, ny
	    hstmp(:,j) = 0.5*(hs(:,j-1)+hs(:,j))
	 end do
	 !hstmp(:,1) = hstmp(:,2)
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
      !if (stag == 1) then
	 !u(1,:,:) = 2*u(2,:,:) - u(3,:,:)
      !else if (stag == 2) then
	 !u(:,1,:) = 2*u(:,2,:) - u(:,3,:)
      !end if
      !
      return
   end subroutine interpolate3Duv
   !
   !
   !
   subroutine read_nhm(input_file, nx, ny, nz, u, v, w, t, pnh, qv, qc, qi, qr, qs, qg, tsoil)
      use variable, only : ngm
      implicit none
      !
      character(*), intent(in) :: input_file
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,ngm), intent(out) :: tsoil
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
   subroutine write_nhm(output_file, u, v, w, t, pnh, qv, qc, qi ,qr, qs, qg, tsoil, nx, ny, nz)
      use variable, only : ngm
      implicit none
      !
      character(*), intent(in) :: output_file
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,ngm), intent(in) :: tsoil
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
end module inner2outerlib