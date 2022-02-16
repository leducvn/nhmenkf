Interface
SUBROUTINE rttov_setpredictors_9_solar( &
            & opts,       &
            & prof,       &
            & coef,       &
            & predictors, &
            & raytracing)
  USE rttov_types, ONLY :  &
       & rttov_options,   &
       & rttov_coef,      &
       & profile_Type,    &
       & predictors_Type, &
       & raytracing_type
  IMPLICIT NONE
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)         ! profile
  TYPE(rttov_coef     ), INTENT(IN)    :: coef            ! coefficients
  TYPE(predictors_Type), INTENT(INOUT) :: predictors      ! predictors
  TYPE(raytracing_Type), INTENT(INOUT) :: raytracing      !
End Subroutine
End Interface
