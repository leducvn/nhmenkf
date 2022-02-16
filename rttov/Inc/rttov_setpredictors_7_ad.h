Interface
SUBROUTINE rttov_setpredictors_7_ad( &
            & opts,          &
            & prof,          &
            & prof_ad,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & predictors,    &
            & predictors_ad, &
            & raytracing,    &
            & raytracing_ad)
  USE rttov_types, ONLY :  &
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & predictors_Type, &
       & raytracing_type
  IMPLICIT NONE
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type   ), INTENT(INOUT) :: prof_ad(size(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_ad
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
  TYPE(profile_aux    ), INTENT(IN)    :: aux                ! auxillary profiles info.
End Subroutine
End Interface
