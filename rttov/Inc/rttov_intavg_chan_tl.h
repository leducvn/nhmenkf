Interface
SUBROUTINE rttov_intavg_chan_tl( &
            & lgradp,      &
            & sun,         &
            & kni,         &
            & kno,         &
            & chanprof,    &
            & ProfIn,      &
            & ProfOut,     &
            & ProfOut_tl,  &
            & OpdepIn,     &
            & OpdepIn_tl,  &
            & OpdepOut,    &
            & OpdepOut_tl)
  USE rttov_types, ONLY : rttov_chanprof, profile_type, opdp_path_type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof(:)
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: lgradp
  LOGICAL(KIND=jplm)  , INTENT(IN)    :: sun(size(chanprof))     ! switch for solar beam
  INTEGER(KIND=jpim)  , INTENT(IN)    :: kni, kno                ! number of levels
  TYPE(profile_type  ), INTENT(IN)    :: ProfIn    (:)           ! atmospheric profiles
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut   (size(ProfIn))! atmospheric profiles
  TYPE(profile_type  ), INTENT(IN)    :: ProfOut_tl(size(ProfIn))! TL of atmospheric profiles
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn                 ! optical depths
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepIn_tl
  TYPE(opdp_path_type), INTENT(IN)    :: OpdepOut                ! optical depths
  TYPE(opdp_path_type), INTENT(INOUT) :: OpdepOut_tl
End Subroutine
End Interface
