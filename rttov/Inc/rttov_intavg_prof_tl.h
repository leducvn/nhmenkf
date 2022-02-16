Interface
SUBROUTINE rttov_intavg_prof_tl( &
            & opts,       &
            & kni,        &
            & kno,        &
            & ProfIn,     &
            & ProfIn_tl,  &
            & ProfOut,    &
            & ProfOut_tl)
  USE rttov_types, ONLY : rttov_options, profile_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options), INTENT(IN)    :: opts
  INTEGER(KIND=jpim) , INTENT(IN)    :: kni , kno               ! number of levels
  TYPE(profile_type) , INTENT(IN)    :: ProfIn    (:)           ! atmospheric profiles
  TYPE(profile_type) , INTENT(IN)    :: ProfIn_tl (size(ProfIn))
  TYPE(profile_type) , INTENT(IN)    :: ProfOut   (size(ProfIn))! atmospheric profiles
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut_tl(size(ProfIn))
End Subroutine
End Interface
