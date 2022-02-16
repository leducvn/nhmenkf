! $Id: ropp_pp_set_coordinates.f90 1882 2008-10-27 15:45:52Z frhl $

!****s* Coordinates/ropp_pp_set_coordinates *
!
! NAME
!    ropp_pp_set_coordinatess - Set internally used position and velocity 
!                               reference frames used in ropp_pp processing
!
! SYNOPSIS
!    use ropp_io
!    use ropp_pp
!      ...
!    type(ROprof) :: rodata
!      ...
!    call ropp_pp_set_coordinates(rodata)
!
! DESCRIPTION
!    This subroutine sets the reference coordinate frame within an ROprof
!    data structure to the reference frame used internally in the ropp_pp 
!    package. For each variable to be defined a call to ropp_utils function
!    coordinates_eci2ecef is made to transform the data if required. 
!
! INPUTS
!    rodata  Radio occultation profile data structure
!
! OUTPUT
!    rodata  As above, but with Level1a data coordinate reference frame 
!             modified to reflect standard reference frames as used within 
!             ropp_pp.
!
! NOTES
!    Default coordinate frames assumed by ropp_pp processing are:
!         r_leo : ECF (Earth Centred Fixed)
!         r_gns : ECF (Earth Centred Fixed)
!         v_leo : ECI (Earth Centred Inertial)
!         v_gns : ECI (Earth Centred Inertial)
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. ROprof, scalar arguments
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_set_coordinates_sca(rodata)

! 1.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
! USE ropp_pp_preproc,  not_this => ropp_pp_set_coordinates_sca

  IMPLICIT NONE

  TYPE(ROprof), INTENT(inout) :: rodata

  IF(rodata%lev1a%reference_frame%r_leo == "ECI")THEN

     rodata%Lev1a%r_leo = eci2ecf(rodata%dtocc%year,  rodata%dtocc%month,  & 
                                  rodata%dtocc%day,   rodata%dtocc%hour,   &
                                  rodata%dtocc%minute,rodata%dtocc%second, &
                                  rodata%Lev1a%dtime +                     &
                                  rodata%dtocc%msec/1000.0_wp,             &
                                  rodata%Lev1a%r_leo)
     rodata%lev1a%reference_frame%r_leo = "ECF"
     
  ENDIF
  
  IF(rodata%lev1a%reference_frame%r_gns == "ECI")THEN
     
     rodata%Lev1a%r_gns = eci2ecf(rodata%dtocc%year,  rodata%dtocc%month,  &
                                  rodata%dtocc%day,   rodata%dtocc%hour,   &
                                  rodata%dtocc%minute,rodata%dtocc%second, &
                                  rodata%Lev1a%dtime +                     &
                                  rodata%dtocc%msec/1000.0_wp,             &
                                  rodata%Lev1a%r_gns)
     rodata%lev1a%reference_frame%r_gns = "ECF"

   ENDIF
      
   IF(rodata%lev1a%reference_frame%v_gns .NE. "ECI")THEN
     CALL message(msg_info, "Expecting v_gns data with reference to ECI frame")
   ENDIF
   
   IF(rodata%lev1a%reference_frame%v_leo .NE. "ECI")THEN
     CALL message(msg_info, "Expecting v_leo data with reference to ECI frame")
   ENDIF
   
 END SUBROUTINE ropp_pp_set_coordinates_sca
                             

!-------------------------------------------------------------------------------
! 2. ROprof, array arguments
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_set_coordinates_arr(rodata_set)

! 2.1 Declarations
! ----------------

  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
! USE ropp_pp_preproc,  not_this => ropp_pp_set_coordinates_arr

  IMPLICIT NONE

  TYPE(ROprof), DIMENSION(:), INTENT(inout) :: rodata_set

  INTEGER                                   :: i

! 2.2 Loop over all elements
! --------------------------

  DO i = 1, SIZE(rodata_set)
     CALL ropp_pp_set_coordinates(rodata_set(i))
  ENDDO

END SUBROUTINE ropp_pp_set_coordinates_arr


