!
SUBROUTINE rttov_alloc_transmission( &
            & ERR,          &
            & transmission, &
            & nlayers,      &
            & nchannels,    &
            & asw,          &
            & init)
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
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : transmission_Type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)     , INTENT(OUT)             :: ERR         ! return code
  TYPE(transmission_type), INTENT(INOUT)           :: transmission
  INTEGER(KIND=jpim)     , INTENT(IN)              :: nlayers     ! number of levels
  INTEGER(KIND=jpim)     , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)     , INTENT(IN)              :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_transmission.h"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlevels
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  nlevels = nlayers + 1
  init1   = .FALSE.
  IF (Present(init)) init1 = init
  IF (asw .EQ. 1) THEN
    ALLOCATE (transmission%tau_levels(nlevels, nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission%tau_levels")
    ALLOCATE (transmission%tau_total(nchannels), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of transmission%tau_total")
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (transmission%tau_levels, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission%tau_levels")
    DEALLOCATE (transmission%tau_total, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of transmission%tau_total")
    NULLIFY (transmission%tau_levels)
    NULLIFY (transmission%tau_total)
  ENDIF
  IF (init1 .AND. (asw .EQ. 1)) THEN
    CALL rttov_init_transmission(transmission)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_transmission
