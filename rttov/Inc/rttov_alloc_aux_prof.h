Interface
SUBROUTINE rttov_alloc_aux_prof( &
            & err,       &
            & nprofiles, &
            & nlevels,   &
            & id_sensor, &
            & aux_prof,  &
            & opts,      &
            & asw,       &
            & init)
  USE rttov_types, ONLY : rttov_options, profile_aux
  USE parkind1, ONLY : jpim, jplm
  INTEGER(KIND=jpim) , INTENT(OUT)               :: err
  INTEGER(KIND=jpim) , INTENT(IN)                :: nprofiles
  INTEGER(KIND=jpim) , INTENT(IN)                :: nlevels
  INTEGER(KIND=jpim) , INTENT(IN)                :: id_sensor
  TYPE(profile_aux  ), INTENT(INOUT)             :: aux_prof
  TYPE(rttov_options), INTENT(IN)                :: opts
  INTEGER(KIND=jpim) , INTENT(IN)                :: asw
  LOGICAL(KIND=jplm) , OPTIONAL     , INTENT(IN) :: init
End Subroutine
End Interface
