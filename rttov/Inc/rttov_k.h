Interface
SUBROUTINE rttov_k( &
            & errorstatus,      &
            & chanprof,         &
            & opts,             &
            & profiles,         &
            & profiles_k,       &
            & coefs,            &
            & calcemis,         &
            & emissivity,       &
            & emissivity_k,     &
            & emissivity_out,   &
            & emissivity_out_k, &
            & transmission,     &
            & transmission_k,   &
            & radiancedata,     &
            & radiancedata_k,   &
            & traj,             &
            & traj_k,           &
            & pccomp,           &
            & pccomp_k,         &
            & profiles_k_pc,    &
            & profiles_k_rec,   &
            & channels_rec)
  USE rttov_types, ONLY :  &
       & rttov_options,     &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & transmission_Type, &
       & profile_Type,      &
       & radiance_Type,     &
       & rttov_chanprof,    &
       & rttov_traj
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  TYPE(profile_Type  )   , INTENT(IN)                           :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(profile_Type  )   , INTENT(INOUT)                        :: profiles_k(size(chanprof))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis        (size(chanprof))
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus
  REAL   (KIND=jprb)     , INTENT(IN)                           :: emissivity      (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_k    (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(OUT)                          :: emissivity_out  (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_out_k(size(chanprof))
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission                    ! in because of meme allocation
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission_k                  ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                    ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata_k
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL     , TARGET :: traj          , traj_k          ! Target is needed here (see rttov_check_temp)
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp_k
  TYPE(profile_Type     ), OPTIONAL     , INTENT(INOUT)         :: profiles_k_pc(:)
  TYPE(profile_Type     ), OPTIONAL     , INTENT(INOUT)         :: profiles_k_rec(:)
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
End Subroutine
End Interface
