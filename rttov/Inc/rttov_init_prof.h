Interface
SUBROUTINE rttov_init_prof(profiles, p)
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : profile_Type
  IMPLICIT NONE
  TYPE(profile_Type), INTENT(INOUT)           :: profiles(:)
  REAL(KIND=jprb)   , INTENT(IN)   , OPTIONAL :: p(:)
End Subroutine
End Interface
