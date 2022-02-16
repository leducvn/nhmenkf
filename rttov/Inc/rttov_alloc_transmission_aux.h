Interface
SUBROUTINE rttov_alloc_transmission_aux( &
            & ERR,              &
            & transmission_aux, &
            & nlayers,          &
            & nchannels,        &
            & asw,              &
            & nstreams,         &
            & init)
  USE rttov_types, ONLY : transmission_Type_aux
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)         , INTENT(OUT)             :: ERR             ! return code
  TYPE(transmission_type_aux), INTENT(INOUT)           :: transmission_aux
  INTEGER(KIND=jpim)         , INTENT(IN)              :: nlayers         ! number of levels
  INTEGER(KIND=jpim)         , INTENT(IN)              :: nchannels
  INTEGER(KIND=jpim)         , INTENT(IN)              :: asw             ! 1=allocate, 0=deallocate
  INTEGER(KIND=jpim)         , INTENT(IN)              :: nstreams
  LOGICAL(KIND=jplm)         , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
