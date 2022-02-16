Interface
SUBROUTINE rttov_cldstr( &
            & profiles,         &
            & cldstr_threshold, &
            & ircld,            &
            & nstreams)
  USE rttov_types, ONLY : profile_Type, ircld_type
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(KIND=jpim), INTENT(OUT)   :: nstreams
  TYPE(PROFILE_TYPE), INTENT(IN)    :: PROFILES(:)
  REAL(KIND=jprb),    INTENT(IN)    :: cldstr_threshold
  TYPE(IRCLD_TYPE  ), INTENT(INOUT) :: IRCLD
End Subroutine
End Interface
