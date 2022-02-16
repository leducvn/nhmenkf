Interface
SUBROUTINE rttov_nullify_coef_scatt_ir(coef_scatt_ir)
  USE parkind1, ONLY : jpim, jprb
  USE rttov_types, ONLY : rttov_coef_scatt_ir
  IMPLICIT NONE
  TYPE(rttov_coef_scatt_ir), INTENT(INOUT) :: coef_scatt_ir
End Subroutine
End Interface
