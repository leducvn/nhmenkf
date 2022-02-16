Interface
SUBROUTINE rttov_calcemis_mw( &
            & opts,             &
            & profiles,         &
            & geometry,         &
            & coef,             &
            & chanprof,         &
            & transmission_aux, &
            & calcemis,         &
            & emissivity,       &
            & reflectivity,     &
            & err)
  USE rttov_types, ONLY :  &
       & rttov_options,         &
       & rttov_chanprof,        &
       & rttov_coef,            &
       & profile_Type,          &
       & transmission_Type_aux, &
       & geometry_Type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  TYPE(rttov_options        ), INTENT(IN)            :: opts
  TYPE(rttov_chanprof       ), INTENT(IN)            :: chanprof(:)
  TYPE(profile_Type         ), INTENT(IN)   , TARGET :: profiles(:)
  TYPE(geometry_Type        ), INTENT(IN)   , TARGET :: geometry(size(profiles))
  TYPE(rttov_coef           ), INTENT(IN)            :: coef
  TYPE(transmission_Type_aux), INTENT(IN)            :: transmission_aux
  LOGICAL(KIND=jplm)         , INTENT(IN)            :: calcemis    (size(chanprof))
  REAL   (KIND=jprb)         , INTENT(INOUT)         :: emissivity  (size(chanprof))
  REAL   (KIND=jprb)         , INTENT(OUT)           :: reflectivity(size(chanprof))
  INTEGER(KIND=jpim)         , INTENT(OUT)           :: err
End Subroutine
End Interface
