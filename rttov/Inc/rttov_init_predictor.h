Interface
SUBROUTINE rttov_init_predictor(predictors)
  USE rttov_types, ONLY : predictors_type
  IMPLICIT NONE
  TYPE(predictors_type), INTENT(INOUT) :: predictors
End Subroutine
End Interface
