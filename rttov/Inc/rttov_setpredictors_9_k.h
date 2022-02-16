Interface
SUBROUTINE rttov_setpredictors_9_k( &
            & opts,         &
            & nlayers,      &
            & angles,       &
            & coef_pccomp,  &
            & chanprof,     &
            & prof,         &
            & prof_k,       &
            & coef,         &
            & predictors,   &
            & predictors_k, &
            & raytracing,   &
            & raytracing_k)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_options,     &
       & rttov_coef_pccomp, &
       & profile_Type,      &
       & geometry_Type,     &
       & predictors_Type,   &
       & raytracing_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options)    , INTENT(IN)    :: opts
  INTEGER(KIND=jpim)     , INTENT(IN)    :: nlayers               ! Number of layers
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(profile_Type     ), INTENT(IN)    :: prof  (:)
  TYPE(profile_Type     ), INTENT(INOUT) :: prof_k(size(chanprof))
  TYPE(geometry_Type    ), INTENT(IN)    :: angles(size(prof)    )
  TYPE(predictors_Type  ), INTENT(IN)    :: predictors
  TYPE(predictors_Type  ), INTENT(INOUT) :: predictors_k
  TYPE(raytracing_type  ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type  ), INTENT(INOUT) :: raytracing_k
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
End Subroutine
End Interface
