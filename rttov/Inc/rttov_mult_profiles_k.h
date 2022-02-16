Interface
SUBROUTINE rttov_mult_profiles_k(profiles_k_rec, profiles_k, clear_k_pc)
  USE parkind1, ONLY : jprb
  USE rttov_types, ONLY : profile_type
  IMPLICIT NONE
  TYPE(profile_type), INTENT(INOUT) :: profiles_k_rec(:)
  TYPE(profile_type), INTENT(IN)    :: profiles_k   (:)
  REAL(KIND=jprb)   , INTENT(IN)    :: clear_k_pc   (:, :, :)! dBT(PC)/dRadiance(RTTOV)
End Subroutine
End Interface
