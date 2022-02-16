Interface
SUBROUTINE rttov_tl( &
            & errorstatus,       &
            & chanprof,          &
            & opts,              &
            & profiles,          &
            & profiles_tl,       &
            & coefs,             &
            & calcemis,          &
            & emissivity,        &
            & emissivity_tl,     &
            & emissivity_out,    &
            & emissivity_out_tl, &
            & transmission,      &
            & transmission_tl,   &
            & radiancedata,      &
            & radiancedata_tl,   &
            & traj,              &
            & traj_tl,           &
            & pccomp,            &
            & pccomp_tl,         &
            & channels_rec)
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & radiance_Type,     &
       & transmission_Type, &
       & profile_Type,      &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & rttov_chanprof,    &
       & rttov_traj
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  TYPE(profile_Type)     , INTENT(IN)                           :: profiles   (:)
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof   (:)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(profile_Type  )   , INTENT(INOUT)                        :: profiles_tl(size(profiles))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis         (size(chanprof))
  REAL(KIND=jprb)        , INTENT(IN)                           :: emissivity       (size(chanprof))
  REAL(KIND=jprb)        , INTENT(IN)                           :: emissivity_tl    (size(chanprof))
  REAL(KIND=jprb)        , INTENT(OUT)                          :: emissivity_out   (size(chanprof))
  REAL(KIND=jprb)        , INTENT(OUT)                          :: emissivity_out_tl(size(chanprof))
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission                     ! in because of meme allocation
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission_tl                  ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                     ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata_tl                  ! in because of meme allocation
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL     , TARGET :: traj           , traj_tl         ! target is needed here
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp_tl
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
End Subroutine
End Interface
