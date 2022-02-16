Interface
SUBROUTINE rttov_copy_transmission(transmission1, transmission2)
  USE rttov_types, ONLY : transmission_Type
  IMPLICIT NONE
  TYPE(transmission_type), INTENT(INOUT) :: transmission1
  TYPE(transmission_type), INTENT(IN)    :: transmission2
End Subroutine
End Interface
