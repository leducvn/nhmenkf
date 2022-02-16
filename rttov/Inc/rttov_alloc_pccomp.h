Interface
SUBROUTINE rttov_alloc_pccomp( &
            & ERR,          &
            & pccomp,       &
            & npcscores,    &
            & asw,          &
            & init,         &
            & nchannels_rec)
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : rttov_pccomp
  INTEGER(KIND=jpim), INTENT(OUT)             :: ERR
  TYPE(rttov_pccomp), INTENT(INOUT)           :: pccomp
  INTEGER(KIND=jpim), INTENT(IN)              :: npcscores
  INTEGER(KIND=jpim), INTENT(IN)              :: asw
  LOGICAL(KIND=jplm), OPTIONAL   , INTENT(IN) :: init
  INTEGER(KIND=jpim), OPTIONAL   , INTENT(IN) :: nchannels_rec
End Subroutine
End Interface
