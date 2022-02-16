SUBROUTINE rttov_init_rad(rad)
! Description:
!   Initialise radiance structure
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
  USE rttov_types, ONLY : radiance_type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(radiance_type), INTENT(INOUT) :: rad
!INTF_END
  rad%clear     = 0._jprb
  rad%cloudy    = 0._jprb
  rad%bt_clear  = 0._jprb
  rad%upclear   = 0._jprb
  rad%dnclear   = 0._jprb
  rad%reflclear = 0._jprb
  rad%overcast  = 0._jprb
  rad%up        = 0._jprb
  rad%down      = 0._jprb
  rad%surf      = 0._jprb
  rad%bt        = 0._jprb
  rad%total     = 0._jprb
END SUBROUTINE 
