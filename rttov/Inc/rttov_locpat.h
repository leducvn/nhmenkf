Interface
SUBROUTINE rttov_locpat( &
            & OPTS,       &
            & profiles,   &
            & aux,        &
            & coef,       &
            & angles,     &
            & raytracing)
  USE rttov_types, ONLY :  &
       & rttov_options,   &
       & rttov_coef,      &
       & profile_aux,     &
       & profile_Type,    &
       & geometry_Type,   &
       & raytracing_type
  IMPLICIT NONE
  TYPE(RTTOV_OPTIONS  ), INTENT(IN)    :: OPTS
  TYPE(PROFILE_TYPE   ), INTENT(IN)    :: PROFILES(:)
  TYPE(PROFILE_AUX    ), INTENT(IN)    :: AUX
  TYPE(GEOMETRY_TYPE  ), INTENT(IN)    :: ANGLES(size(profiles))
  TYPE(RTTOV_COEF     ), INTENT(IN)    :: COEF
  TYPE(RAYTRACING_TYPE), INTENT(INOUT) :: RAYTRACING
End Subroutine
End Interface
