module nus2mfinlib
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
   subroutine read_field(type1, type2, type3, ibase, member, ivalid, level, elem, dtype, nx, ny, found, f)
      character(len=8), intent(in) :: type1
      character(len=4), intent(in) :: type2, type3, member
      character(len=6), intent(in) :: level, elem
      character(len=2), intent(in) :: dtype
      integer, intent(in):: ibase, ivalid
      integer, intent(in) :: nx, ny
      integer, intent(out) :: found
      real, dimension(nx,ny), intent(out) :: f
      integer :: ir, j
      real, dimension(nx,ny) :: work
      !
      call nusdas_read(type1, type2, type3, ibase, member, ivalid, level, elem, work, dtype, nx*ny, ir)
      if (ir <= 0) then
         found = 0
         return
      end if
      found = 1
      do j = 1, ny
          f(:,j) = work(:,ny+1-j)
      end do
      !print*, minval(f), maxval(f)
      !
      return
   end subroutine read_field
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
   subroutine uvw_ru2u(fmap, dnsg2, ru, rv, rw, nx, ny, nz, u, v, w)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,6), intent(in) :: fmap
      real, dimension(nx,ny,nz), intent(in) :: dnsg2, ru, rv, rw
      real, dimension(nx,ny,nz), intent(out) :: u, v, w
      integer :: i, j, k
      !
      do k = 1, nz
         do j = 1, ny
            do i = 2, nx
               u(i,j,k) = ru(i,j,k)*(fmap(i,j,2)+fmap(i-1,j,2))/(dnsg2(i,j,k)+dnsg2(i-1,j,k))
            end do
            u(1,j,k) = ru(1,j,k)*fmap(1,j,2)/dnsg2(1,j,k)
         end do
	 do i = 1, nx
            do j = 2, ny
               v(i,j,k) = rv(i,j,k)*(fmap(i,j,1)+fmap(i,j-1,1))/(dnsg2(i,j,k)+dnsg2(i,j-1,k))
            end do
	    v(i,1,k) = rv(i,1,k)*fmap(i,1,1)/dnsg2(i,1,k)
         end do
      end do
      !
      do j = 1, ny
         do i = 1, nx
            do k = 1, nz-1
               w(i,j,k) = rw(i,j,k)*fmap(i,j,3)/(0.5*(dnsg2(i,j,k)+dnsg2(i,j,k+1)))
            end do
	    w(i,j,nz) = rw(i,j,nz)*fmap(i,j,3)/dnsg2(i,j,nz)
         end do
      end do
      !
      return
   end subroutine uvw_ru2u
   !
   !
   !
   subroutine uvw_c2a(vrdz2, uu, vv, ww, nx, ny, nz, u, v, w)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nz), intent(in) :: vrdz2
      real, dimension(nx,ny,nz), intent(in) :: uu, vv, ww
      real, dimension(nx,ny,nz), intent(out) :: u, v, w
      integer :: i, j, k
      !
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx-1
               u(i,j,k) = (uu(i,j,k)+uu(i+1,j,k))*0.5
            end do
            u(nx,j,k) = uu(nx,j,k)
         end do
	 do i = 1, nx
            do j = 1, ny-1
               v(i,j,k) = (vv(i,j,k)+vv(i,j+1,k))*0.5
            end do
            v(i,ny,k) = vv(i,ny,k)
         end do
      end do
      !
      do j = 1, ny
         do i = 1, nx
	    do k = 2, nz
               w(i,j,k) = (vrdz2(k-1)*ww(i,j,k-1)+vrdz2(k)*ww(i,j,k))/(vrdz2(k-1)+vrdz2(k))
            end do
	    w(i,j,1) = ww(i,j,1)
         end do
      end do
      !
      return
   end subroutine uvw_c2a
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
   subroutine write_ic(fileout, idate, vctrans_type, n_vctrans, ix1, ix2, iy1, iy2, iz1, iz2, &
                     & pseam, ptopm, pgrdm, dtratio, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr, dz, dzl, dzr, &
		     & ptgrd, ptgrdt, sst, psurf, tin, u, v, w, pt, p, qv, qc, qi, qr, qs, qg, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: fileout, nx, ny, nz, idate(5), vctrans_type, n_vctrans, ix1, ix2, iy1, iy2, iz1, iz2
      real, intent(in) :: pseam, ptopm, pgrdm, dtratio, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr, dz, dzl, dzr
      real, dimension(nx,ny), intent(in) :: ptgrd, ptgrdt, sst, psurf
      real, dimension(nx,ny,4), intent(in) :: tin
      real, dimension(nx,ny,nz), intent(in) :: u, v, w, pt, p, qv, qc, qi, qr, qs, qg
      integer :: kt, mtuv, flgq, flgqn, flgv
      !
      kt = 0
      mtuv = 22023
      flgq = 11111
      flgqn = 0
      flgv = 100001010
      !
      write(fileout) kt, mtuv, pseam, ptopm, idate, pgrdm, &
      & flgq, flgqn, flgv, &
      & dtratio, nx, ny, nz, 4, dz, dzl, dzr, iz1, iz2, &
      & vctrans_type, zl_vctrans, zh_vctrans, n_vctrans, &
      & dx, dxl, dxr, ix1, ix2, dy, dyl, dyr, iy1, iy2
      write(fileout) u
      write(fileout) v,w
      write(fileout) pt
      write(fileout) qv,qc,qi,qr,qs,qg
      write(fileout) ptgrd,ptgrdt,sst,psurf,tin,p
      close(fileout)
      !
      return
   end subroutine write_ic
   !
   !
   !
end module nus2mfinlib