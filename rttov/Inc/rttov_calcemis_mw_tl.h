Interface
SUBROUTINE rttov_calcemis_mw_tl( &
            & opts,                &
            & profiles,            &
            & profiles_tl,         &
            & geometry,            &
            & coef,                &
            & chanprof,            &
            & transmission_aux,    &
            & transmission_aux_tl, &
            & calcemis,            &
            & emissivity_tl,       &
            & reflectivity_tl)
  USE rttov_types, ONLY :  &
       & rttov_options,         &
       & rttov_chanprof,        &
       & rttov_coef,            &
       & profile_Type,          &
       & transmission_Type_aux, &
       & geometry_Type
  USE parkind1, ONLY : jprb, jplm
  IMPLICIT NONE
  TYPE(rttov_options        ), INTENT(IN)            :: opts
  TYPE(rttov_chanprof       ), INTENT(IN)            :: chanprof(:)
  TYPE(profile_Type         ), INTENT(IN)   , TARGET :: profiles(:)
  TYPE(geometry_Type        ), INTENT(IN)   , TARGET :: geometry(size(profiles))
  TYPE(rttov_coef           ), INTENT(IN)            :: coef
  TYPE(transmission_Type_aux), INTENT(IN)            :: transmission_aux
  LOGICAL(KIND=jplm)         , INTENT(IN)            :: calcemis   (size(chanprof))
  TYPE(profile_Type         ), INTENT(IN)   , TARGET :: profiles_tl(size(profiles))
  TYPE(transmission_Type_aux), INTENT(IN)            :: transmission_aux_tl
  REAL(KIND=jprb)            , INTENT(INOUT)         :: emissivity_tl  (size(chanprof))
  REAL(KIND=jprb)            , INTENT(OUT)           :: reflectivity_tl(size(chanprof))
End Subroutine
End Interface
