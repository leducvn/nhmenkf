Interface
SUBROUTINE rttov_cldstr_tl( &
            & profiles,    &
            & profiles_tl, &
            & ircld,       &
            & ircld_tl)
  USE rttov_types, ONLY : profile_Type, ircld_type
  IMPLICIT NONE
  TYPE(profile_type), INTENT(IN)    :: profiles   (:)
  TYPE(profile_type), INTENT(IN)    :: profiles_tl(size(profiles))
  TYPE(ircld_type  ), INTENT(INOUT) :: ircld
  TYPE(ircld_type  ), INTENT(INOUT) :: ircld_tl
End Subroutine
End Interface
