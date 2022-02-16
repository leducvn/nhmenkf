Interface
SUBROUTINE rttov_add_prof( &
            & profiles,  &
            & profiles1, &
            & profiles2, &
            & lair,      &
            & lground)
  USE rttov_types, ONLY : profile_type
  USE parkind1, ONLY : jplm
  TYPE(profile_type), INTENT(INOUT)           :: profiles (:)
  TYPE(profile_type), INTENT(IN)              :: profiles1(size(profiles))
  TYPE(profile_type), INTENT(IN)              :: profiles2(size(profiles))
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: lair
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: lground
End Subroutine
End Interface
