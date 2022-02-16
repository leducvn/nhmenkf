Interface
SUBROUTINE rttov_ad( &
            & errorstatus,       &
            & chanprof,          &
            & opts,              &
            & profiles,          &
            & profiles_ad,       &
            & coefs,             &
            & calcemis,          &
            & emissivity,        &
            & emissivity_ad,     &
            & emissivity_out,    &
            & emissivity_out_ad, &
            & transmission,      &
            & transmission_ad,   &
            & radiancedata,      &
            & radiancedata_ad,   &
            & traj,              &
            & traj_ad,           &
            & pccomp,            &
            & pccomp_ad,         &
            & channels_rec)
  USE rttov_types, ONLY :  &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & profile_Type,      &
       & transmission_Type, &
       & radiance_Type,     &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_traj
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  TYPE(profile_Type  )   , INTENT(IN)                           :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(profile_Type  )   , INTENT(INOUT)                        :: profiles_ad(size(profiles))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis         (size(chanprof))
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus
  REAL   (KIND=jprb)     , INTENT(IN)                           :: emissivity       (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_ad    (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(OUT)                          :: emissivity_out   (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_out_ad(size(chanprof))
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission                     ! in because of meme allocation
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission_ad                  ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                     ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata_ad                  ! in because of meme allocation
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL     , TARGET :: traj           , traj_ad         ! target is needed here, because
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp_ad
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
End Subroutine
End Interface
