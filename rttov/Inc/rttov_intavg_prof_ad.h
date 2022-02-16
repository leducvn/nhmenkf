Interface
SUBROUTINE rttov_intavg_prof_ad( &
            & opts,       &
            & kni,        &
            & kno,        &
            & ProfIn,     &
            & ProfIn_ad,  &
            & ProfOut,    &
            & ProfOut_ad)
  USE rttov_types, ONLY : rttov_options, profile_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options), INTENT(IN)    :: opts
  INTEGER(KIND=jpim) , INTENT(IN)    :: kni , kno               ! number of levels
  TYPE(profile_type) , INTENT(IN)    :: ProfIn    (:)           ! atmospheric profiles
  TYPE(profile_type) , INTENT(INOUT) :: ProfIn_ad (size(ProfIn))
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut   (size(ProfIn))! atmospheric profiles
  TYPE(profile_type) , INTENT(INOUT) :: ProfOut_ad(size(ProfIn))
End Subroutine
End Interface
