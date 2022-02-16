Interface
SUBROUTINE rttov_cldstr_ad( &
            & profiles,    &
            & profiles_ad, &
            & ircld,       &
            & ircld_ad)
  USE rttov_types, ONLY : profile_Type, ircld_type
  IMPLICIT NONE
  TYPE(profile_type), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type), INTENT(INOUT) :: profiles_ad(size(profiles))
  TYPE(ircld_type  ), INTENT(INOUT) :: ircld
  TYPE(ircld_type  ), INTENT(INOUT) :: ircld_ad
End Subroutine
End Interface
