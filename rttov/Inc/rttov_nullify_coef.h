Interface
SUBROUTINE rttov_nullify_coef(coef)
  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_coef
  IMPLICIT NONE
  TYPE(rttov_coef), INTENT(INOUT) :: coef
End Subroutine
End Interface
