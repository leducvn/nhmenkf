module nus2grdlib
contains
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
   subroutine compute_phydro(zrp, dvtrans, zs, pb, pt, qv, nx, ny, nz, p)
      use variable, only : rdvcp, cpdvrd, gvcp, presrf
      implicit none
      !
      integer, intent(in)  :: nx, ny, nz
      real, dimension(nz), intent(in)  :: zrp, dvtrans
      real, dimension(nx,ny), intent(in)  :: zs, pb
      real, dimension(nx,ny,nz), intent(in)  :: pt, qv
      real, dimension(nx,ny,nz), intent(out) :: p
      integer :: k
      real, dimension(nx,ny) :: virtpt
      real, dimension(nx,ny,nz) :: g2w, deltaz, pai
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
      p(:,:,1) = pb
      do k = 2, nz
	 virtpt(:,:) = 0.5*(pt(:,:,k)*(1.+0.608*qv(:,:,k)) + pt(:,:,k-1)*(1.+0.608*qv(:,:,k-1)))
	 pai(:,:,k) = pai(:,:,k-1) - gvcp/virtpt*deltaz(:,:,k)
	 p(:,:,k) = presrf*pai(:,:,k)**cpdvrd
      end do
      !
      return
   end subroutine compute_phydro
   !
   !
   !
   subroutine compute_rhs(p, t, td, nx, ny, rh)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny
      real, dimension(nx,ny), intent(in) :: p, t, td
      real, dimension(nx,ny), intent(out) :: rh
      !
      integer :: i, j
      real :: e0c, al, bl, e0i, ai, bi, tc, es, e
      !
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      do j = 1, ny
	 do i = 1, nx
	    tc = td(i,j) - tkelvn
	    if (tc >= tcw ) then
	       e = e0c*exp(al*tc/(bl+tc))
	    else if(tc <= tci) then
	       e = e0i*exp(ai*tc/(bi+tc))
	    else
	       e = e0c*exp(al*tc/(bl+tc))*(tc-tci)/(tcw-tci) + &
		  & e0i*exp(ai*tc/(bi+tc))*(tcw-tc)/(tcw-tci)
	    end if
	    tc = t(i,j) - tkelvn
	    if (tc >= tcw ) then
	       es = e0c*exp(al*tc/(bl+tc))
	    else if(tc <= tci) then
	       es = e0i*exp(ai*tc/(bi+tc))
	    else
	       es = e0c*exp(al*tc/(bl+tc))*(tc-tci)/(tcw-tci) + &
		  & e0i*exp(ai*tc/(bi+tc))*(tcw-tc)/(tcw-tci)
	    end if
	    rh(i,j) = min(1.0d0, max(0.0d0,e/es))
	 end do
      end do
      !
      return
   end subroutine compute_rhs
   !
   !
   !
   subroutine t2es(t, nx, ny, es)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      integer, intent(in) :: nx, ny
      real, dimension(nx,ny), intent(in) :: t
      real, dimension(nx,ny), intent(out) :: es
      integer :: i, j, k
      real :: e0c, al, bl, e0i, ai, bi, tc
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
	    if (tc >= tcw ) then
	       es(i,j) = e0c*exp(al*tc/(bl+tc))
	    else if(tc <= tci) then
	       es(i,j) = e0i*exp(ai*tc/(bi+tc))
	    else
	       es(i,j) = e0c*exp(al*tc/(bl+tc))*(tc-tci)/(tcw-tci) + &
		       & e0i*exp(ai*tc/(bi+tc))*(tcw-tc)/(tcw-tci)
	    end if
	 end do
      end do
      !
      return
   end subroutine t2es
   !
   !
   !
   subroutine rh2qv(p, t, rh, nx, ny, qv)
      ! p(Pa), t(K), rh(fraction), qv(kg/kg)
      implicit none
      !
      integer, intent(in) :: nx, ny
      real, dimension(nx,ny), intent(in) :: p, t, rh
      real, dimension(nx,ny), intent(out) :: qv
      real, dimension(nx,ny) :: es, e
      !
      call t2es(t, nx, ny, es)
      where (es >= p)
         e = rh*p
      elsewhere
         e = rh*es
      end where
      qv = 0.622d0*e/(p-0.378d0*e)
      !
      return
   end subroutine rh2qv
   !
   !
   !
   subroutine write_nhm(output_file, u, v, w, t, p, qv, qc, qi, qr, qs, qg, tsoil, nx, ny, nz, ngm)
      implicit none
      !
      character(*), intent(in) :: output_file
      integer, intent(in) :: nx, ny, nz, ngm
      real, dimension(nx,ny,ngm), intent(in) :: tsoil
      real, dimension(nx,ny,nz), intent(in) :: u, v, w, t, p, qv, qc, qi, qr, qs, qg
      !
      open(90, file=trim(output_file), form='unformatted')
      write(90) u
      write(90) v
      write(90) w
      write(90) t
      write(90) p
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
end module nus2grdlib