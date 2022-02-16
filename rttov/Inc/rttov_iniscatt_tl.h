Interface
Subroutine rttov_iniscatt_tl (&
      & errorstatus,       &! out
      & nlevels,           &! in
      & nchannels,         &! in
      & nprofiles,         &! in
      & chanprof,          &! in
      & frequencies,       &! in
      & profiles,          &! in  
      & profiles_tl,       &! in  
      & cld_profiles,      &! in 
      & cld_profiles_tl,   &! in 
      & coef_rttov,        &! in
      & coef_scatt,        &! in
      & transmission,      &! in
      & transmission_tl,   &! in
      & calcemiss,         &! in
      & usenewcld,         &! in
      & usercfrac,         &! in
      & angles,            &! out
      & scatt_aux,         &! inout
      & scatt_aux_tl)       ! inout 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :     &
       & rttov_coef          ,&
       & rttov_scatt_coef    ,&
       & transmission_type   ,&
       & transmission_type_aux  ,&
       & geometry_Type       ,&
       & profile_scatt_aux   ,&
       & profile_Type        ,&
       & profile_cloud_Type  ,&
       & rttov_chanprof      ,&
       & rttov_options 
  Use rttov_const, Only:      &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & gravity,             &
       & pressure_top,        &
       & rgp,                 &
       & rm,                  &
       & rho_rain,            &
       & rho_snow,            &
       & ccthres 
  Use parkind1, Only : jpim, jprb, jplm
  Implicit None
  Integer (Kind=jpim), Intent (in) :: nlevels    ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles  ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent(out) :: errorstatus             ! Error return code
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Logical (Kind=jplm), Intent (in) :: calcemiss   (nchannels) ! Emissivity flags
  Logical (Kind=jplm), Intent (in) :: usenewcld               ! New or old cloud partition
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction
  Type (profile_Type),        Intent (in)    :: profiles        (nprofiles)   ! Atmospheric profiles
  Type (profile_Type),        Intent (in)    :: profiles_tl     (nprofiles)   ! Atmospheric profiles
  Type (rttov_coef),          Intent (in)    :: coef_rttov                    ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                    ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (nprofiles)   ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles_tl (nprofiles)   ! Cloud profiles 
  Type (transmission_Type),   Intent (in)    :: transmission                  ! Transmittances and optical depths
  Type (transmission_Type),   Intent (in)    :: transmission_tl               ! Transmittances and optical depths
  Type (geometry_Type),       Intent (out)   :: angles          (nprofiles)   ! Zenith angles
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux                     ! Auxiliary profile variables
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux_tl                  ! Auxiliary profile variables
End Subroutine
End Interface
