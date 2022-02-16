Interface
SUBROUTINE rttov_setpredictors_7( &
            & opts,       &
            & prof,       &
            & geom,       &
            & coef,       &
            & aux,        &
            & predictors, &
            & raytracing)
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
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)         ! profile
  TYPE(rttov_coef     ), INTENT(IN)    :: coef            ! coefficients
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))! geometry
  TYPE(predictors_Type), INTENT(INOUT) :: predictors      ! predictors
  TYPE(profile_aux    ), INTENT(IN)    :: aux             ! auxillary profiles info.
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing
End Subroutine
End Interface
