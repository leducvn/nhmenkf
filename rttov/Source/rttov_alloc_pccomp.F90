SUBROUTINE rttov_alloc_pccomp( &
            & ERR,          &
            & pccomp,       &
            & npcscores,    &
            & asw,          &
            & init,         &
            & nchannels_rec)
! Description:
!   Allocates/deallocates the PC component structure
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

!INTF_OFF
#include "throw.h"
!INTF_ON
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_pccomp
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  INTEGER(KIND=jpim), INTENT(OUT)             :: ERR
  TYPE(rttov_pccomp), INTENT(INOUT)           :: pccomp
  INTEGER(KIND=jpim), INTENT(IN)              :: npcscores
  INTEGER(KIND=jpim), INTENT(IN)              :: asw
  LOGICAL(KIND=jplm), OPTIONAL   , INTENT(IN) :: init
  INTEGER(KIND=jpim), OPTIONAL   , INTENT(IN) :: nchannels_rec
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_pccomp.h"
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY


  init1 = .FALSE.
  IF (Present(init)) init1 = init
  IF (asw .EQ. 1) THEN
    NULLIFY (pccomp%pcscores, pccomp%bt_pccomp, pccomp%clear_pccomp)
    ALLOCATE (pccomp%pcscores(npcscores), STAT = err)
    THROWM( ERR .NE. 0 , "allocation of pcscores")
    IF (Present(nchannels_rec)) THEN
      IF (nchannels_rec .GT. 0) THEN
        ALLOCATE (pccomp%bt_pccomp(nchannels_rec), pccomp%clear_pccomp(nchannels_rec), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of bt_pccomp, clear_pccomp")
      ENDIF
    ENDIF
    IF (init1) THEN
      CALL rttov_init_pccomp(pccomp)
    ENDIF
  ENDIF
  IF (asw .EQ. 0) THEN
    IF (Associated(pccomp%pcscores)) THEN
      DEALLOCATE (pccomp%pcscores, STAT = err)
      THROWM( ERR .NE. 0 , "deallocation of pcscores")
    ENDIF
    IF (Associated(pccomp%bt_pccomp)) THEN
      DEALLOCATE (pccomp%bt_pccomp, STAT = err)
      THROWM( ERR .NE. 0 , "deallocation of bt_pccomp")
    ENDIF
    IF (Associated(pccomp%clear_pccomp)) THEN
      DEALLOCATE (pccomp%clear_pccomp, STAT = err)
      THROWM( ERR .NE. 0 , "deallocation of clear_pccomp")
    ENDIF
    NULLIFY (pccomp%pcscores, pccomp%bt_pccomp, pccomp%clear_pccomp)
  ENDIF
  CATCH
END SUBROUTINE 
