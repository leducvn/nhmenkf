Interface
SUBROUTINE rttov_intavg_prof( &
            & opts,    &
            & kni,     &
            & kno,     &
            & ProfIn,  &
            & ProfOut)
  USE rttov_types, ONLY : rttov_options, profile_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options), INTENT(IN)    :: opts
  INTEGER(KIND=jpim) , INTENT(IN)    :: kni , kno ! number of levels
  TYPE(profile_type) , INTENT(IN)    :: ProfIn (:)! atmospheric profiles
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut(:)! atmospheric profiles
End Subroutine
End Interface
