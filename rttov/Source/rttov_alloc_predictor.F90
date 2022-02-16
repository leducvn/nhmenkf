!
SUBROUTINE rttov_alloc_predictor( &
            & ERR,         &
            & npredictors, &
            & predictors,  &
            & coef,        &
            & asw,         &
            & addsolar,    &
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
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Declarations:
! Modules used:
! Imported Parameters:
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_coef, predictors_Type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  INTEGER(KIND=jpim)   , INTENT(OUT)             :: ERR                  ! return code
  INTEGER(KIND=jpim)   , INTENT(IN)              :: npredictors
  TYPE(predictors_Type), INTENT(INOUT)           :: predictors
  TYPE(rttov_coef     ), INTENT(IN)              :: coef
  INTEGER(KIND=jpim)   , INTENT(IN)              :: asw                  ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)   , INTENT(IN)              :: addsolar
  LOGICAL(KIND=jplm)   , INTENT(IN)   , OPTIONAL :: init
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_init_predictor.h"
! Local Arrays and Scalars:
  INTEGER(KIND=jpim) :: nlayers
  LOGICAL(KIND=jplm) :: init1
!- End of header --------------------------------------------------------
  TRY
  init1 = .FALSE.
  IF (Present(init)) init1 = init
  nlayers   = coef%nlayers
!  Allocate section
  IF (asw .EQ. 1) THEN
    NULLIFY (predictors%mixedgas)
    NULLIFY (predictors%watervapour)
    NULLIFY (predictors%ozone)
    NULLIFY (predictors%wvcont)
    NULLIFY (predictors%co2)
    NULLIFY (predictors%n2o)
    NULLIFY (predictors%co)
    NULLIFY (predictors%ch4)
    NULLIFY (predictors%clw)
    NULLIFY (predictors%mixedgas_sun)
    NULLIFY (predictors%watervapour_sun)
    NULLIFY (predictors%ozone_sun)
    NULLIFY (predictors%wvcont_sun)
    NULLIFY (predictors%co2_sun)
    NULLIFY (predictors%n2o_sun)
    NULLIFY (predictors%co_sun)
    NULLIFY (predictors%ch4_sun)
    predictors%nlevels = coef%nlevels
    predictors%nmixed  = coef%nmixed
    predictors%nwater  = coef%nwater
    predictors%nozone  = coef%nozone
    predictors%nwvcont = coef%nwvcont
    predictors%nco2    = coef%nco2
    predictors%nn2o    = coef%nn2o
    predictors%nco     = coef%nco
    predictors%nch4    = coef%nch4
    predictors%ncloud  = 0
  ENDIF
  IF (asw .EQ. 1) THEN
    ALLOCATE (predictors%mixedgas(coef%nmixed, nlayers, npredictors), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of predictors%mixedgas")
    ALLOCATE (predictors%watervapour(coef%nwater, nlayers, npredictors), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of predictors%watervapour")
    ALLOCATE (predictors%clw(nlayers, npredictors), STAT = ERR)
    THROWM( ERR .NE. 0 , "allocation of predictors%clw")
    IF (coef%nozone > 0) THEN
      ALLOCATE (predictors%ozone(coef%nozone, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%ozone")
    ENDIF
    IF (coef%nwvcont > 0) THEN
      ALLOCATE (predictors%wvcont(coef%nwvcont, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%wvcont")
    ENDIF
    IF (coef%nco2 > 0) THEN
      ALLOCATE (predictors%co2(coef%nco2, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%co2")
    ENDIF
    IF (coef%nn2o > 0) THEN
      ALLOCATE (predictors%n2o(coef%nn2o, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%n2o")
    ENDIF
    IF (coef%nco > 0) THEN
      ALLOCATE (predictors%co(coef%nco, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%co")
    ENDIF
    IF (coef%nch4 > 0) THEN
      ALLOCATE (predictors%ch4(coef%nch4, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%ch4")
    ENDIF
    IF (addsolar) THEN
      ALLOCATE (predictors%mixedgas_sun(coef%nmixed, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%mixedgas_sun")
      ALLOCATE (predictors%watervapour_sun(coef%nwater, nlayers, npredictors), STAT = ERR)
      THROWM( ERR .NE. 0 , "allocation of predictors%watervapour_sun")
      IF (coef%nozone > 0) THEN
        ALLOCATE (predictors%ozone_sun(coef%nozone, nlayers, npredictors), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of predictors%ozone_sun")
      ENDIF
      IF (coef%nwvcont > 0) THEN
        ALLOCATE (predictors%wvcont_sun(coef%nwvcont, nlayers, npredictors), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of predictors%wvcont_sun")
      ENDIF
      IF (coef%nco2 > 0) THEN
        ALLOCATE (predictors%co2_sun(coef%nco2, nlayers, npredictors), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of predictors%co2_sun")
      ENDIF
      IF (coef%nn2o > 0) THEN
        ALLOCATE (predictors%n2o_sun(coef%nn2o, nlayers, npredictors), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of predictors%n2o_sun")
      ENDIF
      IF (coef%nco > 0) THEN
        ALLOCATE (predictors%co_sun(coef%nco, nlayers, npredictors), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of predictors%co_sun")
      ENDIF
      IF (coef%nch4 > 0) THEN
        ALLOCATE (predictors%ch4_sun(coef%nch4, nlayers, npredictors), STAT = ERR)
        THROWM( ERR .NE. 0 , "allocation of predictors%ch4_sun")
      ENDIF
    ENDIF
  ENDIF
  IF (asw .EQ. 0) THEN
    DEALLOCATE (predictors%mixedgas, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of predictors%mixedgas")
    DEALLOCATE (predictors%watervapour, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of predictors%watervapour")
    DEALLOCATE (predictors%clw, STAT = ERR)
    THROWM( ERR .NE. 0 , "deallocation of predictors%clw")
    IF (coef%nozone > 0) THEN
      DEALLOCATE (predictors%ozone, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%ozone")
    ENDIF
    IF (coef%nwvcont > 0) THEN
      DEALLOCATE (predictors%wvcont, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%wvcont")
    ENDIF
    IF (coef%nco2 > 0) THEN
      DEALLOCATE (predictors%co2, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%co2")
    ENDIF
    IF (coef%nn2o > 0) THEN
      DEALLOCATE (predictors%n2o, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%n2o")
    ENDIF
    IF (coef%nco > 0) THEN
      DEALLOCATE (predictors%co, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%co")
    ENDIF
    IF (coef%nch4 > 0) THEN
      DEALLOCATE (predictors%ch4, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%ch4")
    ENDIF
    IF (addsolar) THEN
      DEALLOCATE (predictors%mixedgas_sun, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%mixedgas_sun")
      DEALLOCATE (predictors%watervapour_sun, STAT = ERR)
      THROWM( ERR .NE. 0 , "deallocation of predictors%watervapour_sun")
      IF (coef%nozone > 0) THEN
        DEALLOCATE (predictors%ozone_sun, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of predictors%ozone_sun")
      ENDIF
      IF (coef%nwvcont > 0) THEN
        DEALLOCATE (predictors%wvcont_sun, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of predictors%wvcont_sun")
      ENDIF
      IF (coef%nco2 > 0) THEN
        DEALLOCATE (predictors%co2_sun, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of predictors%co2_sun")
      ENDIF
      IF (coef%nn2o > 0) THEN
        DEALLOCATE (predictors%n2o_sun, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of predictors%n2o_sun")
      ENDIF
      IF (coef%nco > 0) THEN
        DEALLOCATE (predictors%co_sun, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of predictors%co_sun")
      ENDIF
      IF (coef%nch4 > 0) THEN
        DEALLOCATE (predictors%ch4_sun, STAT = ERR)
        THROWM( ERR .NE. 0 , "deallocation of predictors%ch4_sun")
      ENDIF
    ENDIF
  ENDIF
  IF (init1 .AND. (asw .EQ. 1)) THEN
    CALL rttov_init_predictor(predictors)
  ENDIF
  IF (asw .EQ. 0) THEN
    NULLIFY (predictors%mixedgas)
    NULLIFY (predictors%watervapour)
    NULLIFY (predictors%ozone)
    NULLIFY (predictors%wvcont)
    NULLIFY (predictors%co2)
    NULLIFY (predictors%n2o)
    NULLIFY (predictors%co)
    NULLIFY (predictors%ch4)
    NULLIFY (predictors%clw)
    NULLIFY (predictors%mixedgas_sun)
    NULLIFY (predictors%watervapour_sun)
    NULLIFY (predictors%ozone_sun)
    NULLIFY (predictors%wvcont_sun)
    NULLIFY (predictors%co2_sun)
    NULLIFY (predictors%n2o_sun)
    NULLIFY (predictors%co_sun)
    NULLIFY (predictors%ch4_sun)
  ENDIF
  CATCH
END SUBROUTINE rttov_alloc_predictor
