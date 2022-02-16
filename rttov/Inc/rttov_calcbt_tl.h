Interface
SUBROUTINE rttov_calcbt_tl( &
            & chanprof, &
            & coeffs,   &
            & rad,      &
            & rad_tl)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, radiance_type
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs
  TYPE(radiance_Type ), INTENT(IN)    :: rad        ! rad%total rad%clear
  TYPE(radiance_Type ), INTENT(INOUT) :: rad_tl
End Subroutine
End Interface
