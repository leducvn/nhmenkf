Interface
SUBROUTINE rttov_copy_raytracing(raytracing1, raytracing2)
  USE rttov_types, ONLY : raytracing_type
  IMPLICIT NONE
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing1
  TYPE(raytracing_type), INTENT(IN)    :: raytracing2
End Subroutine
End Interface
