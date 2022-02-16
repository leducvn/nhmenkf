Interface
SUBROUTINE rttov_alloc_sunglint( &
            & err,       &
            & sunglint,  &
            & nprofiles, &
            & nomega,    &
            & asw,       &
            & init)
  USE parkind1, ONLY : jpim, jplm
  USE rttov_types, ONLY : sunglint_type
  INTEGER(KIND=jpim) , INTENT(INOUT)             :: err
  TYPE(sunglint_type), INTENT(INOUT)             :: sunglint
  INTEGER(KIND=jpim) , INTENT(IN)                :: nprofiles
  INTEGER(KIND=jpim) , INTENT(IN)                :: nomega
  INTEGER(KIND=jpim) , INTENT(IN)                :: asw
  LOGICAL(KIND=jplm) , OPTIONAL     , INTENT(IN) :: init
End Subroutine
End Interface
