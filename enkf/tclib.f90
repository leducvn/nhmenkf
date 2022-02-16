module tclib
! Author: Le Duc
! Created date: 30 Nov 2016
   use variable, only : r_size, pi
   implicit none
   integer, parameter :: tcloc = 1 ! localization flag
   real(r_size), parameter :: dbuffer = 200., dloc = 800. ! km
   !
   !
   !
contains
   !
   !
   !
   subroutine search_max(lon, lat, field, nlon, nlat, lonmax, latmax, fmax)
      implicit none
      integer, intent(in) :: nlon, nlat
      real(r_size), dimension(nlon), intent(in) :: lon
      real(r_size), dimension(nlat), intent(in) :: lat
      real(r_size), dimension(nlon,nlat), intent(in) :: field
      real(r_size), intent(out) :: lonmax, latmax, fmax
      !
      integer, parameter :: norder = 3, nout = 3000
      ! hdown: the first radius (in degree) for finding maximum
      ! hdhmin: the last radius (in degree) for finding maximum
      real(r_size) :: std_tunning = 0.1, hdown = 0.5, hdhmin = 0.05
      integer :: i, j, nex, iex_tmp, jex_tmp, niter, kk, ii, nnr, invalid
      real(r_size) :: dlon, sum, mean, std, fex_tmp, ww(3), wo, hin
      complex(r_size) :: zs, aa, vv, uu(7), z(3)
      integer, dimension(:), allocatable :: iex, jex
      real(r_size), dimension(:), allocatable :: fex
      real(r_size), dimension(16,nlon-1,nlat-1) :: a
      !
      ! Preprocessing
      dlon = (lon(nlon)-lon(1))/(nlon-1)
      hdown = hdown*dlon
      hdhmin = hdhmin*dlon
      call bicubic_itpl_coefs(lat, lon, field, nlat, nlon, a)
      !
      ! Compute mean and standard deviation
      sum = 0.
      do i = 1, nlon
         do j = 1, nlat
            sum = sum + field(i,j)
         enddo
      enddo
      mean = sum/real(nlon*nlat)
      sum = 0.
      do i = 1, nlon
         do j = 1, nlat
            sum = sum + (field(i,j)-mean)*(field(i,j)-mean)
         enddo
      enddo
      std = sqrt(sum/real(nlon*nlat-1))
      !write(6,*) 'Staticstics:'
      !write(6,*) 'Mean: ', mean
      !write(6,*) 'Standard deviation: ', std
      !
      ! Check for extrema existence
      nex=0
      do j = 2, nlat-1
         do i = 2, nlon-1
            fmax = max(field(i+1,j-1), field(i,j-1), field(i-1,j-1), field(i+1,j), &
                       field(i-1,j), field(i+1,j+1), field(i,j+1), field(i-1,j+1))
            if (fmax.le.field(i,j)) nex=nex+1
      enddo
      enddo
      lonmax = -9999.
      latmax = -9999.
      fmax = -9999.
      if (nex.eq.0) then
         !write(6,*) 'No extreme found'
         return
      endif
      !
      allocate(fex(nex))
      allocate(iex(nex))
      allocate(jex(nex))
      nex = 0
      do j = 2, nlat-1
         do i = 2, nlon-1
            fmax = max(field(i+1,j-1), field(i,j-1), field(i-1,j-1), field(i+1,j), &
                       field(i-1,j), field(i+1,j+1), field(i,j+1), field(i-1,j+1))
            if (fmax.le.field(i,j)) then
	       nex=nex+1
	       fex(nex) = field(i,j)
	       iex(nex) = i
	       jex(nex) = j
            endif
         enddo
      enddo
      !write(6,*) 'Find local extrema:'
      !write(6,*) 'Number of local extrema found: ', nex
      !
      ! Find the extreme
      if ((nex.eq.1) .and. (fex(1).lt.(mean+std_tunning*std))) then
         !write(6,*) 'No distinct extreme found'
         return
      endif
      if (nex.gt.1) then
	 ! sorting in ascending order
	 do 30 j = 2, nex
            fex_tmp = fex(j)
            iex_tmp = iex(j)
            jex_tmp = jex(j)
            do 20 i = j-1, 1, -1
               if (fex(i).le.fex_tmp) goto 10
               fex(i+1) = fex(i)
               iex(i+1) = iex(i)
               jex(i+1) = jex(i)
   20       continue
            i = 0
   10       fex(i+1) = fex_tmp
            iex(i+1) = iex_tmp
            jex(i+1) = jex_tmp
   30    continue
         if (fex(nex).lt.(mean+std_tunning*std)) then
            !write(6,*) 'No distinct extreme found'
            return
         endif
         do while (fex(1).lt.(mean+std_tunning*std))
            do j=2, nex
               fex(j-1)=fex(j)
               iex(j-1)=iex(j)
               jex(j-1)=jex(j)
            enddo
            nex=nex-1
         enddo
         !print*, 'Remaining extrema: ', nex
      endif
      !write(6,*) 'From the grid:'
      !write(6,*) 'lonmax: ', lon(iex(nex))
      !write(6,*) 'latmax: ', lat(jex(nex))
      !write(6,*) 'fmax: ', fex(nex)
      !
      ! Initial conditions for DOWNHILL method
      zs = cmplx(lon(iex(nex)), lat(jex(nex)))
      uu(1) = (1.,0.)
      uu(2) = (0.8660254,.5)
      uu(3) = (0.,1.)
      uu(4) = (0.9659258,.2588190)
      uu(5) = (.7071068,.7071068)
      uu(6) = (.2588190,.9659258)
      uu(7) = (-2.588190,.9659258)
      hin   = hdown
      niter = 0
      call bicubic_itpl(lat, lon, a, aimag(zs), real(zs), nlat, nlon, invalid, wo)
      !
      ! Search for maximum
  201 kk = 1
      ii = 0
  202 vv = (-1.,0.)
  203 aa = (-.5,.8660254)
  204 z(1) = zs + hin*vv*aa
      z(2) = zs + hin*vv
      z(3) = zs + hin*conjg(aa)*vv
      call bicubic_itpl(lat, lon, a, aimag(z(1)), real(z(1)), nlat, nlon, invalid, ww(1))
      if (invalid == 1) goto 118
      call bicubic_itpl(lat, lon, a, aimag(z(2)), real(z(2)), nlat, nlon, invalid, ww(2))
      if (invalid == 1) goto 118
      call bicubic_itpl(lat, lon, a, aimag(z(3)), real(z(3)), nlat, nlon, invalid, ww(3))
      if (invalid == 1) goto 118
      niter = niter+1
      if (niter.gt.nout) then
         !write(6,*) 'Iteration number exceeded'
         if (wo-ww(1)) 224,223,223
  223    if (wo-ww(2)) 227,225,225
  225    if (wo-ww(3)) 222,219,219
  224    if (ww(1)-ww(2)) 227,227,226
  226    if (ww(1)-ww(3)) 222,220,220
  227    if (ww(2)-ww(3)) 222,221,221
  219    goto 118
  220    zs = z(1)
         goto 118
  221    zs = z(2)
         goto 118
  222    zs = z(3)
         goto 118
      endif
      if (ww(1)-ww(3)) 206,205,205
  205 if (ww(1)-ww(2)) 208,208,207
  206 if (ww(2)-ww(3)) 209,208,208
  207 nnr = 1
      goto 210
  208 nnr = 2
      goto 210
  209 nnr = 3
  210 if (wo-ww(nnr)) 212,212,211
  211 goto (213,214,215), kk
  212 kk = 1
      ii = 0
      aa = (.7071068,.7071068)
      vv = (z(nnr)-zs)/hin
      wo = ww(nnr)
      zs = z(nnr)
      goto 204
  213 kk = 2
      if (hin.lt.hdhmin) goto 118
      hin = hin*.25
      goto 203
  214 kk = 3
      hin = hin*4.
      goto 202
  215 ii = ii+1
      if (ii-7) 216,216,217
  216 vv = uu(ii)
      goto 203
  217 if (hin.lt.hdhmin) goto 118
      hin = hin*.25
      ii = 0
      goto 202
  118 lonmax = real(zs)
      latmax = aimag(zs)
      fmax = wo
      !
      !write(6,*) 'After Downhill method:'
      !write(6,*) 'Number of iteration: ', niter
      !write(6,*) 'lonmax: ', lonmax
      !write(6,*) 'latmax: ', latmax
      !write(6,*) 'fmax: ', fmax
      deallocate(fex)
      deallocate(iex)
      deallocate(jex)
      !
      return
   end subroutine search_max
   !
   !
   !
   subroutine bicubic_itpl_coefs(x, y, f, nx, ny, a)
      integer, intent(in) :: nx, ny
      real(r_size), dimension(nx), intent(in) :: x
      real(r_size), dimension(ny), intent(in) :: y
      real(r_size), dimension(ny,nx), intent(in) :: f
      real(r_size), dimension(16,ny-1,nx-1), intent(out) :: a
      integer :: i, j
      real(r_size) :: dx, dy
      real(r_size), dimension(ny,nx) :: dfdx, dfdy, d2fdxdy
      !
      dx = (x(nx)-x(1))/(nx-1)
      dy = (y(ny)-y(1))/(ny-1)
      do j = 2, ny-1
         do i = 2, nx-1
	    dfdx(j,i) = 0.5*(f(j,i+1)-f(j,i-1))/dx
	    dfdy(j,i) = 0.5*(f(j+1,i)-f(j-1,i))/dy
            d2fdxdy(j,i) = 0.25*((f(j+1,i+1)-f(j-1,i+1))-(f(j+1,i-1)-f(j-1,i-1)))/(dx*dy)
         enddo
      enddo
      !
      do j = 2, ny-1
         dfdy(j,1) = 0.5*(f(j+1,1)-f(j-1,1))/dy
	 dfdy(j,nx) = 0.5*(f(j+1,nx)-f(j-1,nx))/dy
	 dfdx(j,1) = 0.5*(-3*f(j,1)+4*f(j,2)-f(j,3))/dx
	 dfdx(j,nx) = 0.5*(f(j,nx-2)-4*f(j,nx-1)+3*f(j,nx))/dx
	 d2fdxdy(j,1) = 0.25*(-3*(f(j+1,1)-f(j-1,1))+4*(f(j+1,2)-f(j-1,2))-(f(j+1,3)-f(j-1,3)))/(dx*dy)
	 d2fdxdy(j,nx) = 0.25*((f(j+1,nx-2)-f(j-1,nx-2))-4*(f(j+1,nx-1)-f(j-1,nx-1))+3*(f(j+1,nx)-f(j-1,nx)))/(dx*dy)
      enddo
      do i = 2, nx-1
         dfdx(1,i) = 0.5*(f(1,i+1)-f(1,i-1))/dx
	 dfdx(ny,i) = 0.5*(f(ny,i+1)-f(ny,i-1))/dx
	 dfdy(1,i) = 0.5*(-3*f(1,i)+4*f(2,i)-f(3,i))/dy
	 dfdy(ny,i) = 0.5*(f(ny-2,i)-4*f(ny-1,i)+3*f(ny,i))/dy
	 d2fdxdy(1,i) = 0.25*(-3*(f(1,i+1)-f(1,i-1))+4*(f(2,i+1)-f(2,i-1))-(f(3,i+1)-f(3,i-1)))/(dx*dy)
	 d2fdxdy(ny,i) = 0.25*((f(ny-2,i+1)-f(ny-2,i-1))-4*(f(nx-1,i+1)-f(nx-1,i-1))+3*(f(nx,i+1)-f(nx,i-1)))/(dx*dy)
      enddo
      !
      dfdx(1,1) = 0.5*(-3*f(1,1)+4*f(1,2)-f(1,3))/dx
      dfdx(ny,1) = 0.5*(-3*f(ny,1)+4*f(ny,2)-f(ny,3))/dx
      dfdx(1,nx) = 0.5*(f(1,nx-2)-4*f(1,nx-1)+3*f(1,nx))/dx
      dfdx(ny,nx) = 0.5*(f(ny,nx-2)-4*f(ny,nx-1)+3*f(ny,nx))/dx
      dfdy(1,1) = 0.5*(-3*f(1,1)+4*f(2,1)-f(3,1))/dy
      dfdy(1,nx) = 0.5*(-3*f(1,nx)+4*f(2,nx)-f(3,nx))/dy
      dfdy(ny,1) = 0.5*(f(ny-2,1)-4*f(ny-1,1)+3*f(ny,1))/dy
      dfdy(ny,nx) = 0.5*(f(ny-2,nx)-4*f(ny-1,nx)+3*f(ny,nx))/dy
      d2fdxdy(1,1) = 0.25*(-3*(-3*f(1,1)+4*f(2,1)-f(3,1))+4*(-3*f(1,2)+4*f(2,2)-f(3,2))-(-3*f(1,3)+4*f(2,3)-f(3,3)))/(dx*dy)
      d2fdxdy(ny,1) = 0.25*(-3*(f(ny-2,1)-4*f(ny-1,1)+3*f(ny,1))+4*(f(ny-2,2)-4*f(ny-1,2)+3*f(ny,2))-(f(ny-2,3)-4*f(ny-1,3)+3*f(ny,3)))/(dx*dy)
      d2fdxdy(1,nx) = 0.25*((-3*f(1,nx-2)+4*f(2,nx-2)-f(3,nx-2))-4*(-3*f(1,nx-1)+4*f(2,nx-1)-f(3,nx-1))+3*(-3*f(1,nx)+4*f(2,nx)-f(3,nx)))/(dx*dy)
      d2fdxdy(ny,nx) = 0.25*((f(ny-2,nx-2)-4*f(ny-1,nx-2)+3*f(ny,nx-2))-4*(f(ny-2,nx-1)-4*f(ny-1,nx-1)+3*f(ny,nx-1))+3*(f(ny-2,nx)-4*f(ny-1,nx)+3*f(ny,nx)))/(dx*dy)
      !
      dfdx = dx*dfdx
      dfdy = dy*dfdy
      d2fdxdy = dx*dy*d2fdxdy
      do j = 1, ny-1
         do i = 1, nx-1
	    a(1,j,i) = f(j,i)
	    a(2,j,i) = dfdx(j,i)
	    a(3,j,i) = -3*f(j,i)+3*f(j,i+1)-2*dfdx(j,i)-dfdx(j,i+1)
	    a(4,j,i) = 2*f(j,i)-2*f(j,i+1)+dfdx(j,i)+dfdx(j,i+1)
	    a(5,j,i) = dfdy(j,i)
	    a(6,j,i) = d2fdxdy(j,i)
	    a(7,j,i) = -3*dfdy(j,i)+3*dfdy(j,i+1)-2*d2fdxdy(j,i)-d2fdxdy(j,i+1)
	    a(8,j,i) = 2*dfdy(j,i)-2*dfdy(j,i+1)+d2fdxdy(j,i)+d2fdxdy(j,i+1)
	    a(9,j,i) = -3*f(j,i)+3*f(j+1,i)-2*dfdy(j,i)-dfdy(j+1,i)
	    a(10,j,i) = -3*dfdx(j,i)+3*dfdx(j+1,i)-2*d2fdxdy(j,i)-d2fdxdy(j+1,i)
	    a(11,j,i) = 9*f(j,i)-9*f(j,i+1)-9*f(j+1,i)+9*f(j+1,i+1) + 6*dfdx(j,i)+3*dfdx(j,i+1)-6*dfdx(j+1,i)-3*dfdx(j+1,i+1) + &
	                6*dfdy(j,i)-6*dfdy(j,i+1)+3*dfdy(j+1,i)-3*dfdy(j+1,i+1) + 4*d2fdxdy(j,i)+2*d2fdxdy(j,i+1)+2*d2fdxdy(j+1,i)+d2fdxdy(j+1,i+1)
	    a(12,j,i) = -6*f(j,i)+6*f(j,i+1)+6*f(j+1,i)-6*f(j+1,i+1) - 3*dfdx(j,i)-3*dfdx(j,i+1)+3*dfdx(j+1,i)+3*dfdx(j+1,i+1) - &
	                4*dfdy(j,i)+4*dfdy(j,i+1)-2*dfdy(j+1,i)+2*dfdy(j+1,i+1) - 2*d2fdxdy(j,i)-2*d2fdxdy(j,i+1)-d2fdxdy(j+1,i)-d2fdxdy(j+1,i+1)
	    a(13,j,i) = 2*f(j,i)-2*f(j+1,i)+dfdy(j,i)+dfdy(j+1,i)
	    a(14,j,i) = 2*dfdx(j,i)-2*dfdx(j+1,i)+d2fdxdy(j,i)+d2fdxdy(j+1,i)
	    a(15,j,i) = -6*f(j,i)+6*f(j,i+1)+6*f(j+1,i)-6*f(j+1,i+1) - 4*dfdx(j,i)-2*dfdx(j,i+1)+4*dfdx(j+1,i)+2*dfdx(j+1,i+1) - &
	                3*dfdy(j,i)+3*dfdy(j,i+1)-3*dfdy(j+1,i)+3*dfdy(j+1,i+1) - 2*d2fdxdy(j,i)-d2fdxdy(j,i+1)-2*d2fdxdy(j+1,i)-d2fdxdy(j+1,i+1)
	    a(16,j,i) = 4*f(j,i)-4*f(j,i+1)-4*f(j+1,i)+4*f(j+1,i+1) + 2*dfdx(j,i)+2*dfdx(j,i+1)-2*dfdx(j+1,i)-2*dfdx(j+1,i+1) + &
	                2*dfdy(j,i)-2*dfdy(j,i+1)+2*dfdy(j+1,i)-2*dfdy(j+1,i+1) + d2fdxdy(j,i)+d2fdxdy(j,i+1)+d2fdxdy(j+1,i)+d2fdxdy(j+1,i+1)
	 enddo
      enddo
      !
      return
   end subroutine bicubic_itpl_coefs
   !
   !
   !
   subroutine bicubic_itpl(x, y, a, xstar, ystar, nx, ny, invalid, fstar)
      integer, intent(in) :: nx, ny
      real(r_size), dimension(nx), intent(in) :: x
      real(r_size), dimension(ny), intent(in) :: y
      real(r_size), dimension(16,ny-1,nx-1), intent(in) :: a
      real(r_size), intent(in) :: xstar, ystar
      integer, intent(out) :: invalid
      real(r_size), intent(out) :: fstar
      real(r_size), parameter :: undefined = -9999.
      integer :: i, j, k, ii, jj
      real(r_size) :: dx, dy, t, u
      !
      dx = (x(nx)-x(1))/(nx-1)
      dy = (y(ny)-y(1))/(ny-1)
      !
      if (xstar < x(1) .or. xstar > x(nx) .or. ystar < y(1) .or. ystar > y(ny)) then
	 invalid = 1
	 fstar = undefined
      else
	 invalid = 0
	 j = int((ystar-y(1))/dy) + 1
	 if (j == ny) j = j-1
	 u = (ystar-y(j))/dy
	 i = int((xstar-x(1))/dx) + 1
	 if (i == nx) i = i-1
	 t = (xstar-x(i))/dx
	 fstar = 0.
	 do k = 1, 16
	    ii = mod(k-1,4)
	    jj = (k-1)/4
	    fstar = fstar + a(k,j,i)*t**ii*u**jj
	 enddo
      endif
      !
      return
   end subroutine bicubic_itpl
   !
   !
   !
   subroutine distance(rlon1, rlat1, rlon2, rlat2, d)
      implicit none
      real(r_size), intent(in) :: rlon1, rlat1
      real(r_size), intent(in) :: rlon2, rlat2
      real(r_size), intent(out) :: d
      real(r_size) :: lon1, lat1, lon2, lat2, dlon, w1, w2, w3, w4
      !
      lon1 = rlon1*pi/180.
      lat1 = rlat1*pi/180.
      lon2 = rlon2*pi/180.
      lat2 = rlat2*pi/180.
      dlon = lon1-lon2
      w1 = cos(lat2)*sin(dlon)
      w2 = cos(lat1)*sin(lat2)-sin(lat1)*cos(lat2)*cos(dlon)
      w3 = sin(lat1)*sin(lat2)
      w4 = cos(lat1)*cos(lat2)*cos(dlon)
      d = 180./pi*atan2(sqrt(w1*w1+w2*w2),(w3+w4))
      !
      return
   end subroutine distance
   !
   !
   !
end module tclib
