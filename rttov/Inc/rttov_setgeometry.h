Interface
SUBROUTINE rttov_setgeometry( &
            & opts,       &
            & profiles,   &
            & aux,        &
            & coef,       &
            & angles,     &
            & raytracing)! inout, optional
  USE rttov_types, ONLY :  &
       & rttov_coef,      &
       & profile_Type,    &
       & profile_aux,     &
       & geometry_Type,   &
       & raytracing_type, &
       & rttov_options
  IMPLICIT NONE
  TYPE(rttov_options  ), INTENT(IN)              :: opts
  TYPE(profile_Type   ), INTENT(IN)              :: profiles(:)           ! profile
  TYPE(rttov_coef     ), INTENT(IN)              :: coef                  ! coefficient
  TYPE(profile_aux    ), INTENT(IN)   , OPTIONAL :: aux
  TYPE(geometry_Type  ), INTENT(OUT)             :: angles(size(profiles))! angles
  TYPE(raytracing_type), INTENT(INOUT), OPTIONAL :: raytracing
End Subroutine
End Interface
