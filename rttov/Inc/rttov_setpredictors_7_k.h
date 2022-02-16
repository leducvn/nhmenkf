Interface
SUBROUTINE rttov_setpredictors_7_k( &
            & opts,         &
            & nlayers,      &
            & angles,       &
            & chanprof,     &
            & profiles,     &
            & profiles_k,   &
            & coef,         &
            & aux_prof,     &
            & predictors,   &
            & predictors_k, &
            & raytracing,   &
            & raytracing_k)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & predictors_Type, &
       & raytracing_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options)  , INTENT(IN)    :: opts
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                   ! Number of layers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type   ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type   ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(geometry_Type  ), INTENT(IN)    :: angles    (size(profiles))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
  TYPE(profile_aux    ), INTENT(IN)    :: aux_prof                  ! auxillary profiles info.
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
End Subroutine
End Interface
