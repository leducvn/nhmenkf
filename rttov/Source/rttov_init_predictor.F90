SUBROUTINE rttov_init_predictor(predictors)
! Description:
!   Initialise predictors structure
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
  USE rttov_types, ONLY : predictors_type
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
  TYPE(predictors_type), INTENT(INOUT) :: predictors
!INTF_END
  IF (Associated(predictors%mixedgas)       ) predictors%mixedgas        = 0._jprb
  IF (Associated(predictors%watervapour)    ) predictors%watervapour     = 0._jprb
  IF (Associated(predictors%ozone)          ) predictors%ozone           = 0._jprb
  IF (Associated(predictors%wvcont)         ) predictors%wvcont          = 0._jprb
  IF (Associated(predictors%co2)            ) predictors%co2             = 0._jprb
  IF (Associated(predictors%n2o)            ) predictors%n2o             = 0._jprb
  IF (Associated(predictors%co)             ) predictors%co              = 0._jprb
  IF (Associated(predictors%ch4)            ) predictors%ch4             = 0._jprb
  IF (Associated(predictors%clw)            ) predictors%clw             = 0._jprb
  IF (Associated(predictors%mixedgas_sun)   ) predictors%mixedgas_sun    = 0._jprb
  IF (Associated(predictors%watervapour_sun)) predictors%watervapour_sun = 0._jprb
  IF (Associated(predictors%ozone_sun)      ) predictors%ozone_sun       = 0._jprb
  IF (Associated(predictors%wvcont_sun)     ) predictors%wvcont_sun      = 0._jprb
  IF (Associated(predictors%co2_sun)        ) predictors%co2_sun         = 0._jprb
  IF (Associated(predictors%n2o_sun)        ) predictors%n2o_sun         = 0._jprb
  IF (Associated(predictors%co_sun)         ) predictors%co_sun          = 0._jprb
  IF (Associated(predictors%ch4_sun)        ) predictors%ch4_sun         = 0._jprb
END SUBROUTINE 
