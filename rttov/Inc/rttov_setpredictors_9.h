Interface
SUBROUTINE rttov_setpredictors_9( &
            & opts,        &
            & prof,        &
            & geom,        &
            & coef_pccomp, &
            & coef,        &
            & predictors,  &
            & raytracing)
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
  TYPE(profile_Type     ), INTENT(IN)    :: prof(:)         ! profile
  TYPE(rttov_coef       ), INTENT(IN)    :: coef            ! coefficients
  TYPE(geometry_Type    ), INTENT(IN)    :: geom(size(prof))! geometry
  TYPE(predictors_Type  ), INTENT(INOUT) :: predictors      ! predictors
  TYPE(raytracing_Type  ), INTENT(INOUT) :: raytracing      !
End Subroutine
End Interface
