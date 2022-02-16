Interface
SUBROUTINE rttov_setgeometry_ad( &
            & opts,          &
            & profiles,      &
            & profiles_ad,   &
            & aux,           &
            & coef,          &
            & angles,        &
            & raytracing,    &
            & raytracing_ad)
  USE rttov_types, ONLY :  &
       & rttov_coef,      &
       & profile_Type,    &
       & profile_aux,     &
       & geometry_Type,   &
       & raytracing_type, &
       & rttov_options
  IMPLICIT NONE
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: profiles   (:)             ! profile
  TYPE(profile_Type   ), INTENT(INOUT) :: profiles_ad(size(profiles))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                       ! coefficient
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(geometry_Type  ), INTENT(IN)    :: angles(size(profiles))     ! angles
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_ad
End Subroutine
End Interface
