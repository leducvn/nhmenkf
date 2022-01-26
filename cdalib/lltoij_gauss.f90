subroutine lltoij_gauss(fi, fj, x, y, lat, lon, glat, nx0, ny0)
  implicit none
  real(4), intent(out)         :: fi, fj, x, y
  real(4), intent(in)          :: lat, lon, glat(nx0, ny0)
  integer(4), intent(in)       :: nx0, ny0
!
  integer(4)                   :: j
!
  fi = lon * nx0 / 360. + 1
  x = fi - int(fi)
  if (lat <= glat(1,1) .and. lat > glat(1,ny0)) then
    do j = 2, ny0
      if (lat > glat(1,j)) then
        fj = (j - 1) + (glat(1,j-1) - lat) / (glat(1,j-1) - glat(1,j))
        exit
      end if
    end do
    y = fj - int(fj)
  else if (lat > glat(1,1)) then  ! near the north pole
    fj = (90. - lat) / (90. - glat(1,1))
    if (fj < 0.) then
      fj = 0.
    end if
    y = fj - int(fj)
  else if (lat <= glat(1,ny0)) then ! near the south pole
    fj = ny0 &
        & + (glat(1,ny0) - lat) / (glat(1,ny0) - (-90.))
    fj = ny0 &
        & + (glat(1,ny0) - lat) / (glat(1,ny0) - (-90.))
    if (fj >= ny0 + 1) then ! at the south pole
      fj = ny0
      y = 1.
    else
      y = fj - int(fj)
    end if
  end if
end subroutine lltoij_gauss
