! $Id: occ_point.f90 2019 2009-01-14 10:20:26Z frhl $

!****f* Coordinates/tangent_point
!
! NAME
!    tangent_point - Determine tangent point coordinates
!
! SYNOPSIS
!    call occ_point(r_leo, r_gns, lat_tp, lon_tp, azimuth_tp)
! 
! DESCRIPTION
!    This subroutine calculates the latitude, longitude and azimuth of
!    each tangent point for an occultation
!
! INPUTS
!    r_leo         cartesian LEO position vector (relative to ECF frame)
!    r_gns         cartesian GPS position vector (relative to ECF frame)
!
! OUTPUT
!    lat_tp        tangent point latitude
!    lon_tp        tangent point longitude
!    azimuth_tp    GPS to LEO azimuth direction wrt true North at tangent (deg)
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

subroutine tangent_point(r_leo, r_gns, lat_tp, lon_tp, azimuth_tp)

! 1.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
! use coordinates, not_this => tangent_point
  use coordinates, only: rad2deg, &
                       & rotate, &
                       & impact_parameter, &
                       & vector_angle, &
                       & vector_product

  implicit none

  real(wp), dimension(:,:), intent(in) :: r_leo   ! LEO position vector (ECF)
  real(wp), dimension(:,:), intent(in) :: r_gns   ! GPS position vector (ECF)
  real(wp), dimension(:), intent(out)  :: lat_tp  ! Tangent point latitude
  real(wp), dimension(:), intent(out)  :: lon_tp  ! Tangent point longitude
  real(wp), dimension(:), intent(out)  :: azimuth_tp  ! Tangent point azimuth

  real(wp), dimension(size(r_leo,1),size(r_leo,2)) :: perigee
  real(wp)                             :: slta    ! Straight line tangent ht
  real(wp)                             :: ro      ! Length of r_leo
  real(wp)                             :: alpha   ! Angle r_leo and perigee
  real(wp)                             :: theta   ! Cross section azimuth
  real(wp), dimension(size(r_leo,1))   :: height_per
  real(wp), dimension(3), parameter    :: pa = (/0,0,1/)   ! Polar axis
  

  integer:: i

  do i=1,size(r_leo,1)

! 1.2 Determine ray tangent points
! --------------------------------
     
    slta = impact_parameter(r_leo(i,:), r_gns(i,:))
    ro = Sqrt(Dot_Product(r_leo(i,:), r_leo(i,:)))   
    alpha = acos(slta/ro)
     
    perigee(i,:) = rotate(r_leo(i,:), vector_product(r_leo(i,:), r_gns(i,:)), alpha) * (slta/ro)

  enddo
  
! 1.3 Convert cartesian to geodetic points (TP latitude and longitude)
! ---------------------------------------- 

  call cart2geod(perigee, lat_tp, lon_tp, height_per)

! 1.4 Cross-section azimuth at tangent point
! ------------------------------------------ 

  do i=1,size(r_leo,1)
  
    theta = vector_angle( vector_product(perigee(i,:), pa),       &
                          vector_product(r_gns(i,:), r_leo(i,:)), &
                          -perigee(i,:) )

    azimuth_tp(i) = theta * rad2deg

    if (azimuth_tp(i) < 0.0_wp ) azimuth_tp(i) = azimuth_tp(i) + 360.0_wp

  enddo

end subroutine tangent_point



