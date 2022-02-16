Interface
SUBROUTINE rttov_cldstr_k( &
            & chanprof,   &
            & profiles,   &
            & profiles_k, &
            & ircld,      &
            & ircld_k)
  USE rttov_types, ONLY :  &
       & rttov_chanprof, &
       & profile_Type,   &
       & ircld_type
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_type  ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_type  ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(ircld_type    ), INTENT(INOUT) :: ircld
  TYPE(ircld_type    ), INTENT(INOUT) :: ircld_k
End Subroutine
End Interface
