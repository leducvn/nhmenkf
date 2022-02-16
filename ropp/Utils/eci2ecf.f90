! $Id: eci2ecf.f90 2019 2009-01-14 10:20:26Z frhl $

!****f* Coordinates/eci2ecf
!
! NAME
!    eci2ecf - Transform occultation geometry relative to ECI
!              reference frame to geometry relative to ECF frame.
!
! SYNOPSIS
!    r_ecf = eci2ecf(year, month, day, hour, minute, sec, dsec, r_eci)
! 
! DESCRIPTION
!    This subroutine calculates occultation coordinates relative to the ECF 
!    frame.
!
! INPUTS
!    year           occultation start year
!    month          occultation start month
!    day            occultation start day
!    hour           occultation start hour
!    minute         occultation start minute
!    sec            occultation start second
!    dsec           time since occultation start
!    r_eci          cartesian position vector (relative to ECI frame)
!
! OUTPUT
!    r_ecf          cartesian position vector (relative to ECF frame)  
!
! NOTES
!
! SEE ALSO
!    ecf2eci
!
! REFERENCES
!     Astronomical Alamanus, 1993
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

!-------------------------------------------------------------------------------
! 1. Double, scalar arguments 
!-------------------------------------------------------------------------------

function eci2ecf_0d(year, month, day, hour, minute, sec, dsec, r_eci) result(r_ecf)

! 1.1 Declarations
! ----------------

  use typesizes,      only: wp => EightByteReal
  use datetimeprogs,  only: CalTojul
  use coordinates,    only: rotate, Pi

  implicit none

  integer,  intent(inout)            :: year     ! occultation start year
  integer,  intent(inout)            :: month    ! occultation start month
  integer,  intent(inout)            :: day      ! occultation start day
  integer,  intent(in)               :: hour     ! occultation start hour
  integer,  intent(in)               :: minute   ! occultation start minute
  integer , intent(in)               :: sec      ! occultation start second
  real(wp), intent(in)               :: dsec     ! time since occ start
  real(wp), dimension(3), intent(in) :: r_eci    ! position vector (ECI)
  real(wp), dimension(3)             :: r_ecf    ! position vector (ECF)

  integer,  dimension(8)             :: cdt      ! date/time array
  real(wp)                           :: jdf      ! Julian Day & fraction
  real(wp)                           :: phi             ! GPS frame rotation
  real(wp), dimension(3), parameter  :: pa = (/0,0,1/)  ! Polar axis
  real(wp)                           :: tu, gmst, utco, utc

! 1.2 Calculate Greenwich Mean Sidereal Time
! ------------------------------------------

  cdt = (/year,month,day,0,0,0,0,0/)
  call CalToJul ( cdt, jdf, 1 )
  tu   = (jdf - 2451545.0_wp)/36525.0_wp
  gmst = 24110.54841_wp + 8640184.812866_wp*tu + &
         0.093104_wp*tu**2 - 6.2e-6_wp*tu**3
  utco = hour*3600.0_wp + minute*60.0_wp + sec
  utc  = (utco + dsec)*1.0027379093_wp
  gmst = Modulo(gmst + utc, 86400.0_wp)

! 1.3 Calculate Greenwich Mean Sidereal Time Angle
! ------------------------------------------------
 
  phi = gmst*2.0_wp*Pi/86400.0_wp        

! 1.2 Frame rotation
! ------------------

  r_ecf = rotate(r_eci, pa, -phi)     !! N.B. ECI -> ECF, negative phi

end function eci2ecf_0d


!-------------------------------------------------------------------------------
! 2. Double, array argument for position 
!-------------------------------------------------------------------------------

function eci2ecf_1d(year, month, day, hour, minute, sec, dsec, r_eci) result(r_ecf)

! 2.1 Declarations
! ----------------

  use typesizes, only: wp => EightByteReal
  use coordinates, only: eci2ecf

  implicit none

  integer,  intent(inout)            :: year     ! occultation start year
  integer,  intent(inout)            :: month    ! occultation start month
  integer,  intent(inout)            :: day      ! occultation start day
  integer,  intent(in)               :: hour     ! occultation start hour
  integer,  intent(in)               :: minute   ! occultation start minute
  integer,  intent(in)               :: sec      ! occultation start second
  real(wp), dimension(:), intent(in) :: dsec     ! time since occ start
  real(wp), dimension(:,:), intent(in)             :: r_eci   ! position (ECI)
  real(wp), dimension(size(r_eci,1),size(r_eci,2)) :: r_ecf   ! position (ECF)

  integer :: i
  
! 2.2 Frame transformation
! ------------------------

  do i=1,size(dsec)
    r_ecf(i,:) =  eci2ecf(year, month, day, hour, minute, &
                          sec, dsec(i), r_eci(i,:))
  enddo
  
end function eci2ecf_1d



