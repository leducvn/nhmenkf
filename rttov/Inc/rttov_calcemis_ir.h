Interface
SUBROUTINE rttov_calcemis_ir( &
            & profiles,    &
            & geometry,    &
            & coef,        &
            & addpc,       &
            & coef_pccomp, &
            & chanprof,    &
            & calcemis,    &
            & emissivity,  &
            & err)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_coef_pccomp, &
       & profile_Type,      &
       & geometry_Type
  USE parkind1, ONLY : jprb, jpim, jplm
  IMPLICIT NONE
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: addpc
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type     ), INTENT(IN)    :: profiles(:)
  TYPE(geometry_Type    ), INTENT(IN)    :: geometry(:)
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: calcemis  (size(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT) :: emissivity(size(chanprof))
  INTEGER(KIND=jpim)     , INTENT(OUT)   :: err
End Subroutine
End Interface
