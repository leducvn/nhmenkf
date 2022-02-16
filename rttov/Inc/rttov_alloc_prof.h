Interface
SUBROUTINE rttov_alloc_prof( &
            & ERR,      &
            & nprof,    &
            & profiles, &
            & nlevels,  &
            & opts,     &
            & asw,      &
            & coefs,    &
            & init,     &
            & blob)
  USE rttov_types, ONLY : rttov_options, profile_Type, rttov_coefs, blob_type
  USE parkind1, ONLY : jpim, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim) , INTENT(OUT)             :: ERR            ! return code
  INTEGER(KIND=jpim) , INTENT(IN)              :: nprof          ! number of profiles
  INTEGER(KIND=jpim) , INTENT(IN)              :: nlevels        ! number of levels
  TYPE(profile_Type ), INTENT(INOUT)           :: profiles(nprof)! profiles
  TYPE(rttov_options), INTENT(IN)              :: opts
  INTEGER(KIND=jpim) , INTENT(IN)              :: asw            ! 1=allocate, 0=deallocate
  TYPE(rttov_coefs)  , INTENT(IN)   , OPTIONAL :: coefs
  LOGICAL(KIND=jplm) , INTENT(IN)   , OPTIONAL :: init
  TYPE(blob_type)    , INTENT(INOUT), OPTIONAL :: blob
End Subroutine
End Interface
