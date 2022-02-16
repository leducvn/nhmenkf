SUBROUTINE rttov_init_pccomp( pccomp )
! Description:
!   Initialise PC components structure
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
  USE rttov_types, ONLY : rttov_pccomp
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(rttov_pccomp), INTENT(INOUT) :: pccomp
!INTF_END
    IF (Associated(pccomp%pcscores)) THEN
      pccomp%pcscores = 0._jprb
    ENDIF
    IF (Associated(pccomp%bt_pccomp)) THEN
      pccomp%bt_pccomp = 0._jprb
    ENDIF
    IF (Associated(pccomp%clear_pccomp)) THEN
      pccomp%clear_pccomp = 0._jprb
    ENDIF
END SUBROUTINE 
