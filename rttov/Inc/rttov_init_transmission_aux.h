Interface
SUBROUTINE rttov_init_transmission_aux(transmission_aux)
  USE rttov_types, ONLY : transmission_Type_aux
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(transmission_type_aux), INTENT(INOUT) :: transmission_aux
End Subroutine
End Interface
