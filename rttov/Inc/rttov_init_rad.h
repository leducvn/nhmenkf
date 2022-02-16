Interface
SUBROUTINE rttov_init_rad(rad)
  USE rttov_types, ONLY : radiance_type
  IMPLICIT NONE
  TYPE(radiance_type), INTENT(INOUT) :: rad
End Subroutine
End Interface
