Interface
SUBROUTINE rttov_setpredictors_9_solar_k( &
            & opts,         &
            & chanprof,     &
            & prof,         &
            & prof_k,       &
            & coef,         &
            & predictors,   &
            & predictors_k, &
            & raytracing,   &
            & raytracing_k)
  USE rttov_types, ONLY :  &
       & rttov_options,   &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & profile_Type,    &
       & predictors_Type, &
       & raytracing_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options)  , INTENT(IN)    :: opts
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type   ), INTENT(IN)    :: prof    (:)
  TYPE(profile_Type   ), INTENT(INOUT) :: prof_k  (size(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
End Subroutine
End Interface
