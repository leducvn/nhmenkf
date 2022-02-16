Interface
SUBROUTINE RTTOV_LOCPAT_TL( &
            & OPTS,          &
            & PROFILES,      &
            & PROFILES_TL,   &
            & AUX,           &
            & COEF,          &
            & ANGLES,        &
            & RAYTRACING,    &
            & RAYTRACING_TL)
  USE rttov_types, ONLY :  &
       & rttov_options,   &
       & rttov_coef,      &
       & profile_aux,     &
       & profile_Type,    &
       & geometry_Type,   &
       & raytracing_type
  IMPLICIT NONE
  TYPE(RTTOV_OPTIONS  ), INTENT(IN)    :: OPTS
  TYPE(PROFILE_TYPE   ), INTENT(IN)    :: PROFILES   (:)             ! Atmospheric profiles
  TYPE(PROFILE_TYPE   ), INTENT(INOUT) :: PROFILES_TL(size(profiles))! Atmospheric profiles
  TYPE(PROFILE_AUX    ), INTENT(IN)    :: AUX
  TYPE(GEOMETRY_TYPE  ), INTENT(IN)    :: ANGLES(size(profiles))     ! angles
  TYPE(RTTOV_COEF     ), INTENT(IN)    :: COEF
  TYPE(RAYTRACING_TYPE), INTENT(IN)    :: RAYTRACING
  TYPE(RAYTRACING_TYPE), INTENT(INOUT) :: RAYTRACING_TL
End Subroutine
End Interface
