! $Id: occ_point.f90 2019 2009-01-14 10:20:26Z frhl $

!****f* Coordinates/occ_point
!
! NAME
!    occ_point - Determine the occultation point
!
! SYNOPSIS
!    call occ_point(r_leo, r_gns, lat, lon, r_coc, roc, azimuth, undulation, 
!                   cfile, efile)
! 
! DESCRIPTION
!    This subroutine calculates the lowest occultation perigree point projected
!    to the Earth's surface
!
! INPUTS
!    r_leo         cartesian LEO position vector (relative to ECF frame)
!    r_gns         cartesian GPS position vector (relative to ECF frame)
!    cfile         path to geoid potential coefficients file (optional)
!    efile         path to geoid potential corrections file (optional)
!
! OUTPUT
!    lat           Occultation point latitude
!    lon           Occultation point longitude
!    r_coc         Cartesian centre of curvature vector for occ point (ECF)
!    roc           Radius of curvature value for occultation point
!    azimuth       GPS to LEO azimuth direction wrt true North (deg)
!    undulation    Difference between ellipsoid (WGS-84) and EGM-96 geoid (m)
!
! SEE ALSO
!
! REFERENCES
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

subroutine occ_point(r_leo, r_gns, lat, lon, r_coc, roc, azimuth, &
                     undulation, cfile,efile)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
! use coordinates, not_this => occ_point
  use coordinates
  use EarthMod, only: datum_hmsl

  implicit none

  real(wp), dimension(:,:), intent(in) :: r_leo   ! LEO position vector (ECF)
  real(wp), dimension(:,:), intent(in) :: r_gns   ! GPS position vector (ECF)
  real(wp), intent(out)                :: lat     ! Occultation point latitude
  real(wp), intent(out)                :: lon     ! Occultation point longitude
  real(wp), dimension(size(r_leo,2)), intent(out) :: r_coc  ! Centre curvature
  real(wp), intent(out)                :: roc               ! Radius curvature
  real(wp), intent(out)                :: azimuth           ! Azimuth (deg)
  real(wp), intent(out)                :: undulation        ! Undulation
  character(len=*), optional, intent(in) :: cfile   ! Coefficient file path
  character(len=*), optional, intent(in) :: efile   ! Corrections file path

  real(wp), dimension(size(r_leo,1),size(r_leo,2)) :: perigee
  real(wp)                             :: slta    ! Straight line tangent ht
  real(wp)                             :: ro      ! Length of r_leo
  real(wp)                             :: alpha   ! Angle r_leo and perigee
  real(wp)                             :: theta   ! Cross section azimuth
  real(wp), dimension(size(r_leo,1))   :: height_per
  real(wp), dimension(size(r_leo,1))   :: lat_per
  real(wp), dimension(size(r_leo,1))   :: lon_per
  real(wp), dimension(3), parameter    :: pa = (/0,0,1/)   ! Polar axis

  integer:: i, iocc

  do i=1,size(r_leo,1)

! 1.2 Determine ray tangent points
! --------------------------------
     
     slta = impact_parameter(r_leo(i,:), r_gns(i,:))
     ro = Sqrt(Dot_Product(r_leo(i,:), r_leo(i,:)))   
     alpha = acos(slta/ro)
     
     perigee(i,:) = rotate(r_leo(i,:), vector_product(r_leo(i,:), r_gns(i,:)), alpha) * (slta/ro)

  enddo
  
! 1.3 Convert cartesian to geodetic points
! ---------------------------------------- 

  call cart2geod(perigee, lat_per, lon_per, height_per)

! 1.4 Find the lowest perigee
! ---------------------------
 
  iocc = Sum(MinLoc(Abs(height_per)))

! 1.5 Define occultation point latitude and longitude
! ---------------------------------------------------

  lat = lat_per(iocc)
  lon = lon_per(iocc)

! 1.6 Determine occultation point centre of curvature and radius
! --------------------------------------------------------------

  ! 1.6.1 Cross-section azimuth

  theta = vector_angle( vector_product(perigee(iocc,:), pa),          &
                        vector_product(r_gns(iocc,:), r_leo(iocc,:)), &
                        -perigee(iocc,:) )

  azimuth = theta * rad2deg

  ! 1.6.2 Compute curvature

  call curvature(lat, lon, theta, r_coc, roc)

! 1.7 Find undulation - height difference between local ellipsoid and geoid
! -------------------------------------------------------------------------

  if (present(cfile) .and. present(efile)) then
    call datum_hmsl("WGS84", (/lat, lon, 0.0_wp/), undulation, cfile, efile)
  else
    call datum_hmsl("WGS84", (/lat, lon, 0.0_wp/), undulation)
  endif
  if (undulation > -999999.000) undulation = -1.0_wp * undulation

end subroutine occ_point



