Interface
SUBROUTINE rttov_nullify_coef_pccomp(coef_pccomp)
  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_coef_pccomp
  IMPLICIT NONE
  TYPE(rttov_coef_pccomp), INTENT(INOUT) :: coef_pccomp
End Subroutine
End Interface
