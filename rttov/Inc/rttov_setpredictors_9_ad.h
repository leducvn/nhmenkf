Interface
SUBROUTINE rttov_setpredictors_9_ad( &
            & opts,          &
            & prof,          &
            & prof_ad,       &
            & geom,          &
            & coef_pccomp,   &
            & coef,          &
            & predictors,    &
            & predictors_ad, &
            & raytracing,    &
            & raytracing_ad)
  USE rttov_types, ONLY :  &
       & rttov_coef,        &
       & rttov_options,     &
       & rttov_coef_pccomp, &
       & profile_Type,      &
       & geometry_Type,     &
       & predictors_Type,   &
       & raytracing_type
  IMPLICIT NONE
  TYPE(rttov_options    ), INTENT(IN)    :: opts
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  TYPE(profile_Type     ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type     ), INTENT(INOUT) :: prof_ad(size(prof))
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(geometry_Type    ), INTENT(IN)    :: geom(size(prof))
  TYPE(predictors_Type  ), INTENT(IN)    :: predictors
  TYPE(predictors_Type  ), INTENT(INOUT) :: predictors_ad
  TYPE(raytracing_type  ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type  ), INTENT(INOUT) :: raytracing_ad
End Subroutine
End Interface
