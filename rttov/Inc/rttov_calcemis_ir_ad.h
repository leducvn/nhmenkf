Interface
SUBROUTINE rttov_calcemis_ir_ad( &
            & profiles,      &
            & profiles_ad,   &
            & coef,          &
            & addpc,         &
            & coef_pccomp,   &
            & chanprof,      &
            & calcemis,      &
            & emissivity_ad)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,    &
       & rttov_coef,        &
       & rttov_coef_pccomp, &
       & profile_Type
  USE parkind1, ONLY : jprb, jplm
  IMPLICIT NONE
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: addpc
  TYPE(rttov_chanprof   ), INTENT(IN)    :: chanprof   (:)
  TYPE(profile_Type     ), INTENT(INOUT) :: profiles_ad(:)
  TYPE(profile_Type     ), INTENT(IN)    :: profiles   (:)
  TYPE(rttov_coef       ), INTENT(IN)    :: coef
  TYPE(rttov_coef_pccomp), INTENT(IN)    :: coef_pccomp
  LOGICAL(KIND=jplm)     , INTENT(IN)    :: calcemis     (size(chanprof))
  REAL(KIND=jprb)        , INTENT(INOUT) :: emissivity_ad(size(chanprof))
End Subroutine
End Interface
