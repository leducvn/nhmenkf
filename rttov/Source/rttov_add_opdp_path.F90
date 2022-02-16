!
SUBROUTINE rttov_add_opdp_path(opdp_path, opdp_path1, opdp_path2)
! Description:
! allocation/deallocation of a opdp_path structure
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2007, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code
!  1.1       11/10/2007  Add addclouds, addaerosl, init logicals
!                        nullify unused pointers P.Marguinaud
!  1.2       03/11/2009  Transmittances / optical depths on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY : opdp_path_Type
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(opdp_path_Type), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type), INTENT(IN)    :: opdp_path1
  TYPE(opdp_path_Type), INTENT(IN)    :: opdp_path2
!INTF_END
!- End of header --------------------------------------------------------
  opdp_path%atm_level = opdp_path1%atm_level + opdp_path2%atm_level
  opdp_path%sun_level = opdp_path1%sun_level + opdp_path2%sun_level
END SUBROUTINE 
