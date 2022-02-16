Interface
SUBROUTINE rttov_direct( &
            & errorstatus,    &
            & chanprof,       &
            & opts,           &
            & profiles,       &
            & coefs,          &
            & calcemis,       &
            & emissivity,     &
            & emissivity_out, &
            & transmission,   &
            & radiancedata,   &
            & traj,           &
            & traj_dyn,       &
            & traj_sta,       &
            & pccomp,         &
            & channels_rec,   &
            & lbl_check)
  USE rttov_types, ONLY :  &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & profile_Type,      &
       & transmission_type, &
       & radiance_Type,     &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_traj,        &
       & rttov_traj_dyn,    &
       & rttov_traj_sta,    &
       & rttov_lbl_check
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  TYPE(profile_Type  )   , INTENT(IN)                           :: profiles(:)                   ! Atmospheric profiles
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof(:)                   ! Channel indices (nchannels)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs                         ! It is necessary to have "Target" attribute here
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis      (size(chanprof))! switch for emmissivity calc.
  REAL(KIND=jprb)        , INTENT(IN)                           :: emissivity    (size(chanprof))! surface emmissivity
  REAL(KIND=jprb)        , INTENT(INOUT)                        :: emissivity_out(size(chanprof))! surface emmissivity
  TYPE(transmission_type), INTENT(INOUT)                        :: transmission                  ! transmittances and singlelayer optical depths (on User levels)
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                  ! radiances (mw/cm-1/ster/sq.m) and degK
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus                   ! return flag
  TYPE(rttov_traj  )     , INTENT(INOUT), OPTIONAL     , TARGET :: traj                          ! Target is *NEEDED* here (see rttov_check_temp)
  TYPE(rttov_traj_dyn)   , INTENT(INOUT), OPTIONAL     , TARGET :: traj_dyn                      
  TYPE(rttov_traj_sta)   , INTENT(INOUT), OPTIONAL     , TARGET :: traj_sta
  TYPE(rttov_pccomp)     , OPTIONAL     , INTENT(INOUT)         :: pccomp
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
  TYPE(rttov_lbl_check)  , OPTIONAL     , INTENT(IN)            :: lbl_check
End Subroutine
End Interface
