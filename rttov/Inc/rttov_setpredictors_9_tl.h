Interface
SUBROUTINE rttov_setpredictors_9_tl( &
            & opts,          &
            & prof,          &
            & prof_tl,       &
            & geom,          &
            & coef_pccomp,   &
            & coef,          &
            & predictors,    &
            & predictors_tl, &
            & raytracing,    &
            & raytracing_tl)
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
  TYPE(profile_Type     ), INTENT(IN)    :: prof_tl(size(prof))
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(geometry_Type    ), INTENT(IN)    :: geom(size(prof))
  TYPE(predictors_Type  ), INTENT(IN)    :: predictors
  TYPE(predictors_Type  ), INTENT(INOUT) :: predictors_tl      ! in because of mem allocation
  TYPE(raytracing_type  ), INTENT(IN)    :: raytracing
  TYPE(raytracing_type  ), INTENT(IN)    :: raytracing_tl
End Subroutine
End Interface
