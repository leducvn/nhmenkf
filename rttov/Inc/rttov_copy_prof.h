Interface
SUBROUTINE rttov_copy_prof( &
            & profiles1, &
            & profiles2, &
            & larray,    &
            & lscalar)
  USE rttov_types, ONLY : profile_type
  USE parkind1, ONLY : jplm
  TYPE(profile_type), INTENT(INOUT)           :: profiles1(:)
  TYPE(profile_type), INTENT(IN)              :: profiles2(size(profiles1))
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: larray
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: lscalar
End Subroutine
End Interface
