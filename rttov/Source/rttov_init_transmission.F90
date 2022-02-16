!
SUBROUTINE rttov_init_transmission( transmission )
! Description:
! allocation/deallocation of a profile structure
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
  USE rttov_types, ONLY : transmission_Type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(transmission_type), INTENT(INOUT)           :: transmission
!INTF_END
  transmission%tau_total  = 0._jprb
  transmission%tau_levels = 0._jprb
END SUBROUTINE 
