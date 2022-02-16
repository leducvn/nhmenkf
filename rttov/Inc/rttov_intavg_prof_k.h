Interface
SUBROUTINE rttov_intavg_prof_k( &
            & opts,      &
            & chanprof,  &
            & kni,       &
            & kno,       &
            & ProfIn,    &
            & ProfIn_k,  &
            & ProfOut,   &
            & ProfOut_k)
  USE rttov_types, ONLY : rttov_options, rttov_chanprof, profile_type
  USE parkind1, ONLY : jpim
  IMPLICIT NONE
  TYPE(rttov_options ), INTENT(IN)    :: opts
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)  , INTENT(IN)    :: kni, kno                 ! number of levels
  TYPE(profile_type)  , INTENT(IN)    :: ProfIn   (:)             ! atmospheric profiles
  TYPE(profile_type)  , INTENT(INOUT) :: ProfIn_k (size(chanprof))
  TYPE(profile_type)  , INTENT(INOUT) :: ProfOut  (size(ProfIn)  )! atmospheric profiles
  TYPE(profile_type)  , INTENT(INOUT) :: ProfOut_k(size(chanprof))
End Subroutine
End Interface
