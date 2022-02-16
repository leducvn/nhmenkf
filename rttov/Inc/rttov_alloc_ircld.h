Interface
SUBROUTINE rttov_alloc_ircld( &
            & ERR,       &
            & nircld,    &
            & irclds,    &
            & nlayers,   &
            & asw,       &
            & addaerosl, &
            & addclouds, &
            & init)
  USE rttov_types, ONLY : ircld_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim), INTENT(OUT)             :: ERR      ! return code
  INTEGER(KIND=jpim), INTENT(IN)              :: nircld   !
  INTEGER(KIND=jpim), INTENT(IN)              :: nlayers  ! number of layers
  TYPE(ircld_Type)  , INTENT(INOUT)           :: irclds
  INTEGER(KIND=jpim), INTENT(IN)              :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: addclouds
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
