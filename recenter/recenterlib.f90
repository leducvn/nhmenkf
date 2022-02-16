module recenterlib
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
   subroutine read_hydro_icbc(filein, nx, ny, nz, pseam, ptopm, pgrdm, ptgrd, ptgrdt, sst, psurf, tin, u, v, pt, qv)
      implicit none
      !
      integer, intent(in) :: filein, nx, ny, nz
      real, intent(out) :: pseam, ptopm, pgrdm
      real, dimension(nx,ny), intent(out) :: ptgrd, ptgrdt, sst, psurf
      real, dimension(nx,ny,4), intent(out) :: tin
      real, dimension(nx,ny,nz), intent(out) :: u, v, pt, qv
      integer :: kt, mtuv, idate(5), flgq, flgqn, flgv, nx0, ny0, nz0, ngm, iz1, iz2, vctrans_type, n_vctrans, ix1, ix2, iy1, iy2
      real :: dtratio, dz, dzl, dzr, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr
      real, dimension(nx,ny,nz) :: qc, qci, qr, qs, qg
      !
      read(filein) kt, mtuv, pseam, ptopm, idate, pgrdm, &
      & flgq, flgqn, flgv, &
      & dtratio, nx0, ny0, nz0, ngm, dz, dzl, dzr, iz1, iz2, &
      & vctrans_type, zl_vctrans, zh_vctrans, n_vctrans, &
      & dx, dxl, dxr, ix1, ix2, dy, dyl, dyr, iy1, iy2
      read(filein) u
      read(filein) v
      read(filein) pt
      read(filein) qv
      !read(filein) qv,qc,qci,qr,qs,qg
      read(filein) ptgrd,ptgrdt,sst,psurf,tin
      !
      return
   end subroutine read_hydro_icbc
   !
   !
   !
   subroutine read_nonhydro_icbc(filein, readw, readq, readp, nx, ny, nz, pseam, ptopm, pgrdm, ptgrd, ptgrdt, sst, psurf, tin, u, v, w, pt, p, qv, qc, qi, qr, qs, qg)
      implicit none
      !
      integer, intent(in) :: filein, readw, readq, readp, nx, ny, nz
      real, intent(out) :: pseam, ptopm, pgrdm
      real, dimension(nx,ny), intent(out) :: ptgrd, ptgrdt, sst, psurf
      real, dimension(nx,ny,4), intent(out) :: tin
      real, dimension(nx,ny,nz), intent(out) :: u, v, w, pt, p, qv, qc, qi, qr, qs, qg
      integer :: kt, mtuv, idate(5), flgq, flgqn, flgv, nx0, ny0, nz0, ngm, iz1, iz2, vctrans_type, n_vctrans, ix1, ix2, iy1, iy2
      real :: dtratio, dz, dzl, dzr, zl_vctrans, zh_vctrans, dx, dxl, dxr, dy, dyl, dyr
      !
      read(filein) kt, mtuv, pseam, ptopm, idate, pgrdm, &
      & flgq, flgqn, flgv, &
      & dtratio, nx0, ny0, nz0, ngm, dz, dzl, dzr, iz1, iz2, &
      & vctrans_type, zl_vctrans, zh_vctrans, n_vctrans, &
      & dx, dxl, dxr, ix1, ix2, dy, dyl, dyr, iy1, iy2
      read(filein) u
      if (readw == 1) then
	 read(filein) v,w
      else
         read(filein) v
      endif
      read(filein) pt
      if (readq == 1) then
	 read(filein) qv, qc, qi, qr, qs, qg
      else
         read(filein) qv
      endif
      !
      if (readp == 1) then
	 read(filein) ptgrd,ptgrdt,sst,psurf,tin,p
      else
         read(filein) ptgrd,ptgrdt,sst,psurf,tin
      endif
      close(filein)
      !
      return
   end subroutine read_nonhydro_icbc
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
      open(90, file=trim(input_file), form='unformatted')
      read(90) u
      read(90) v
      read(90) w
      read(90) t
      read(90) pnh
      read(90) qv
      read(90) qc
      read(90) qi
      read(90) qr
      read(90) qs
      read(90) qg
      read(90) tsoil
      close(90)
      !
      return
   end subroutine read_nhm
   !
   !
   !
   subroutine compute_paihydro(zrp, dvtrans, zs, pb, pt, qv, nx, ny, nz, pai)
      use variable, only : rdvcp, cpdvrd, gvcp, presrf
      implicit none
      !
      integer, intent(in)  :: nx, ny, nz
      real, dimension(nz), intent(in)  :: zrp, dvtrans
      real, dimension(nx,ny), intent(in)  :: zs, pb
      real, dimension(nx,ny,nz), intent(in)  :: pt, qv
      real, dimension(nx,ny,nz), intent(out) :: pai
      integer :: k
      real, dimension(nx,ny) :: virtpt
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
      pai(:,:,1) = (pb/presrf)**rdvcp
      do k = 2, nz
	 virtpt(:,:) = 0.5*(pt(:,:,k)*(1.+0.608*qv(:,:,k)) + pt(:,:,k-1)*(1.+0.608*qv(:,:,k-1)))
	 pai(:,:,k) = pai(:,:,k-1) - gvcp/virtpt*deltaz(:,:,k)
      end do
      !
      return
   end subroutine compute_paihydro
   !
   !
   !
   subroutine compute_qvs(p, t, nx, ny, nz, sqv)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nx,ny,nz), intent(in) :: p, t
      real, dimension(nx,ny,nz), intent(out) :: sqv
      !
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
                  sqv(i,j,k) = e0c*exp(al*tc/(bl+tc))
               else if(tc <= tci) then
                  sqv(i,j,k) = e0i*exp(ai*tc/(bi+tc))
               else
                  sqv(i,j,k) = e0c*exp(al*tc/(bl+tc))*(tc-tci)/(tcw-tci) + &
                             & e0i*exp(ai*tc/(bi+tc))*(tcw-tc)/(tcw-tci)
               end if
               sqv(i,j,k) = 0.622d0*sqv(i,j,k)/p(i,j,k)
            end do
         end do
      end do
      !
      return
   end subroutine compute_qvs
   !
   !
   !
   subroutine compute_qv(zrp, dvtrans, zs, pb, pt, pnh, twtr, nx, ny, nz, qv, qc)
      use variable, only : rdvcp, cpdvrd, gvcp, presrf, tkelvn, &
                         & tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         & tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in)  :: nx, ny, nz
      real, dimension(nz), intent(in)  :: zrp, dvtrans
      real, dimension(nx,ny), intent(in)  :: zs, pb
      real, dimension(nx,ny,nz), intent(in)  :: pt, pnh, twtr
      real, dimension(nx,ny,nz), intent(out) :: qv, qc
      !
      integer :: i, j, k, iter
      real :: e0c, al, bl, e0i, ai, bi, tc, virtpt, error, epsilon = 1.e-12
      real, dimension(nz) :: pai, p, t, qvs, qvold
      real, dimension(nx,ny,nz) :: g2w, deltaz
      !
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      do k = 1, nz
	 g2w(:,:,k) = 1. + dvtrans(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      do j = 1, ny
	 do i = 1, nx
	    qv(i,j,:) = twtr(i,j,:)
	    error = 1.
	    iter = 0
	    do while (error > epsilon)
	       iter = iter + 1
	       ! pai, p
	       pai(1) = (pb(i,j)/presrf)**rdvcp
	       do k = 2, nz
		  virtpt = 0.5*(pt(i,j,k)*(1.+0.608*qv(i,j,k)) + pt(i,j,k-1)*(1.+0.608*qv(i,j,k-1)))
		  pai(k) = pai(k-1) - gvcp/virtpt*deltaz(i,j,k)
	       end do
	       p(:) = presrf*pai(:)**cpdvrd
	       p(:) = p(:) + pnh(i,j,:)
	       t(:) = pt(i,j,:)*(p/presrf)**rdvcp
	       !
	       ! qvs
	       qvold(:) = qv(i,j,:)
	       do k= 1, nz
		  tc = t(k) - tkelvn
		  if (tc >= tcw) then
		     qvs(k) = e0c*exp(al*tc/(bl+tc))
		  else if(tc <= tci) then
		     qvs(k) = e0i*exp(ai*tc/(bi+tc))
		  else
		     qvs(k) = e0c*exp(al*tc/(bl+tc))*(tc-tci)/(tcw-tci) + &
			    & e0i*exp(ai*tc/(bi+tc))*(tcw-tc)/(tcw-tci)
		  end if
		  qvs(k) = 0.622d0*qvs(k)/p(k)
		  qv(i,j,k) = min(twtr(i,j,k), qvs(k))
	       end do
	       error = maxval(abs(qv(i,j,:)-qvold(:))/qvold(:))
	       if (iter > 10000) then
		  !print*, i, j, error
		  !print*, 1000*twtr(i,j,:)
		  !print*, 1000*qvs(:)
		  !print*, 1000*qv(i,j,:)
	          !stop
		  exit
	       end if
	    end do
	 end do
      end do
      qc = twtr - qv
      !
      return
   end subroutine compute_qv
   !
   !
   !
   subroutine compute_p(zrp, dvtrans_p, zs, psfc, t, qv, p, nx, ny, nz)
      use variable, only : g0, rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nz), intent(in) :: zrp, dvtrans_p
      real, dimension(nx,ny), intent(in) :: zs, psfc
      real, dimension(nx,ny,nz), intent(in) :: t, qv
      real, dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real, dimension(nx,ny) :: factor
      real, dimension(nx,ny,nz) :: g2w, deltaz
      !
      do k = 1, nz-1
         g2w(:,:,k) = 1. + zs(:,:)*(dvtrans_p(k+1)+dvtrans_p(k))*0.5
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      do k = 1, nz
	 if (k == 1) then
	    factor(:,:) = 0.5*(t(:,:,1)*(1.+0.608*qv(:,:,1))+t(:,:,2)*(1.+0.608*qv(:,:,2)))
	    factor(:,:) = g0*deltaz(:,:,1)/rd/factor(:,:)
	    p(:,:,1) = psfc(:,:)*exp(-factor(:,:))
	 else
	    factor(:,:) = 0.5*(t(:,:,k-1)*(1.+0.608*qv(:,:,k-1))+t(:,:,k)*(1.+0.608*qv(:,:,k)))
	    factor(:,:) = g0*deltaz(:,:,k)/rd/factor(:,:)
	    p(:,:,k) = p(:,:,k-1)*exp(-factor(:,:))
	 end if
      end do
      !
      return
   end subroutine compute_p
   !
   !
   !
   subroutine pst2p(zrp, dvtrans, zs, ps, t, qv, nx, ny, nz, p)
      use variable, only : g0, rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nz), intent(in) :: zrp, dvtrans
      real, dimension(nx,ny), intent(in) :: zs, ps
      real, dimension(nx,ny,nz), intent(in) :: t, qv
      real, dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real, dimension(nx,ny) :: factor
      real, dimension(nx,ny,nz) :: g2w, deltaz
      !
      do k = 1, nz-1
         g2w(:,:,k) = 1. + dvtrans(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      factor = 0.5*(t(:,:,1)*(1.+0.608*qv(:,:,1))+t(:,:,2)*(1.+0.608*qv(:,:,2)))
      factor = g0/rd*deltaz(:,:,1)/factor(:,:)
      p(:,:,1) = ps(:,:)*exp(-factor(:,:))
      do k = 2, nz
	 factor = 0.5*(t(:,:,k-1)*(1.+0.608*qv(:,:,k-1))+t(:,:,k)*(1.+0.608*qv(:,:,k)))
	 factor = g0/rd*deltaz(:,:,k)/factor(:,:)
	 p(:,:,k) = p(:,:,k-1)*exp(-factor(:,:))
      end do
      !
      return
   end subroutine pst2p
   !
   !
   !
   subroutine pbt2p(zrp, dvtrans, zs, pb, t, qv, nx, ny, nz, p)
      use variable, only : g0, rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nz), intent(in) :: zrp, dvtrans
      real, dimension(nx,ny), intent(in) :: zs, pb
      real, dimension(nx,ny,nz), intent(in) :: t, qv
      real, dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real, dimension(nx,ny) :: factor
      real, dimension(nx,ny,nz) :: g2w, deltaz
      !
      do k = 1, nz-1
         g2w(:,:,k) = 1. + dvtrans(k)*zs(:,:)
      end do
      deltaz(:,:,1) = zrp(1)*g2w(:,:,1)
      do k = 2, nz
	 deltaz(:,:,k) = (zrp(k)-zrp(k-1))*g2w(:,:,k-1)
      end do
      !
      p(:,:,1) = pb(:,:)
      do k = 2, nz
	 factor = 0.5*(t(:,:,k-1)*(1.+0.608*qv(:,:,k-1))+t(:,:,k)*(1.+0.608*qv(:,:,k)))
	 factor = g0/rd*deltaz(:,:,k)/factor(:,:)
	 p(:,:,k) = p(:,:,k-1)*exp(-factor(:,:))
      end do
      !
      return
   end subroutine pbt2p
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
   subroutine uvw_a2c(vrdz, uu, vv, ww, nx, ny, nz, u, v, w)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real, dimension(nz), intent(in) :: vrdz
      real, dimension(nx,ny,nz), intent(in) :: uu, vv, ww
      real, dimension(nx,ny,nz), intent(out) :: u, v, w
      integer :: i, j, k
      !
      do k = 1, nz
         do j = 1, ny
            do i = 2, nx
               u(i,j,k) = (uu(i-1,j,k)+uu(i,j,k))*0.5
            end do
            u(1,j,k) = uu(1,j,k)
         end do
	 do i = 1, nx
            do j = 2, ny
               v(i,j,k) = (vv(i,j-1,k)+vv(i,j,k))*0.5
            end do
            v(i,1,k) = vv(i,1,k)
         end do
      end do
      !
      do j = 1, ny
         do i = 1, nx
	    do k = 1, nz-1
               w(i,j,k) = (vrdz(k)*ww(i,j,k)+vrdz(k+1)*ww(i,j,k+1))/(vrdz(k)+vrdz(k+1))
            end do
	    w(i,j,nz) = ww(i,j,nz)
         end do
      end do
      w(1:nx,1:ny,1) = 0.0 ! ground
      !
      return
   end subroutine uvw_a2c
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
end module recenterlib