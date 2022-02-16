Interface
SUBROUTINE rttov_alloc_predictor( &
            & ERR,         &
            & npredictors, &
            & predictors,  &
            & coef,        &
            & asw,         &
            & addsolar,    &
            & init)
  USE rttov_types, ONLY : rttov_coef, predictors_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)   , INTENT(OUT)             :: ERR                  ! return code
  INTEGER(KIND=jpim)   , INTENT(IN)              :: npredictors
  TYPE(predictors_Type), INTENT(INOUT)           :: predictors
  TYPE(rttov_coef     ), INTENT(IN)              :: coef
  INTEGER(KIND=jpim)   , INTENT(IN)              :: asw                  ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)   , INTENT(IN)              :: addsolar
  LOGICAL(KIND=jplm)   , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
