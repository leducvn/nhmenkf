Interface
SUBROUTINE rttov_alloc_transmission( &
            & ERR,          &
            & transmission, &
            & nlayers,      &
            & nchannels,    &
            & asw,          &
            & init)
  USE rttov_types, ONLY : transmission_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)     , INTENT(OUT)             :: ERR         ! return code
  TYPE(transmission_type), INTENT(INOUT)           :: transmission
  INTEGER(KIND=jpim)     , INTENT(IN)              :: nlayers     ! number of levels
  INTEGER(KIND=jpim)     , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)     , INTENT(IN)              :: asw         ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)     , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
