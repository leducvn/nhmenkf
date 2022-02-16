Interface
SUBROUTINE rttov_alloc_opdp_path( &
            & ERR,       &
            & opdp_path, &
            & nlevels,   &
            & nchannels, &
            & asw,       &
            & init)
  USE rttov_types, ONLY : opdp_path_Type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)  , INTENT(OUT)             :: ERR      ! return code
  INTEGER(KIND=jpim)  , INTENT(IN)              :: nlevels
  INTEGER(KIND=jpim)  , INTENT(IN)              :: nchannels
  TYPE(opdp_path_Type), INTENT(INOUT)           :: opdp_path
  INTEGER(KIND=jpim)  , INTENT(IN)              :: asw      ! 1=allocate, 0=deallocate
  LOGICAL(KIND=jplm)  , INTENT(IN)   , OPTIONAL :: init
End Subroutine
End Interface
