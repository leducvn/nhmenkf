Interface
SUBROUTINE rttov_calcbt(chanprof, coeffs, rad)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, radiance_Type
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)! Array of channel indices.
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs     ! Coefficients
  TYPE(radiance_Type ), INTENT(INOUT) :: rad        ! input radiances and output BT
End Subroutine
End Interface
