Interface
Subroutine rttov_scatt_tl( &
     & errorstatus,        &! out
     & nlevels,            &! in
     & chanprof,           &! in
     & frequencies,        &! in
     & profiles,           &! in  
     & cld_profiles,       &! in
     & coef_rttov,         &! in
     & coef_scatt,         &! in
     & calcemiss,          &! in
     & emissivity_in,      &! in
     & profiles_tl,        &! inout
     & cld_profiles_tl,    &! in
     & emissivity_in_tl,   &! in
     & radiance,           &! inout
     & radiance_tl,        &! inout 
     & lnewcld,            &! in, optional
     & lusercfrac)          ! in, optional
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_const, Only :   &
       & errorstatus_success ,&
       & errorstatus_fatal, &
       & sensor_id_mw     , &
       & ccthres          ,&  
       & sensor_id_po
  Use rttov_types, Only :    &
       & rttov_coefs          ,&
       & rttov_scatt_coef     ,&
       & geometry_Type        ,&
       & profile_Type         ,&
       & profile_cloud_Type   ,&
       & profile_scatt_aux    ,&
       & transmission_Type    ,&
       & radiance_Type        ,&
       & rttov_chanprof       ,&
       & rttov_options          
  Use parkind1, Only : jpim, jprb, jplm
  Implicit None
  Integer (Kind=jpim), Intent (in)  :: nlevels ! Number of levels
  Type(rttov_chanprof),Intent (in)  :: chanprof(:)             ! Indices
  Type (profile_Type), Intent (in)  :: profiles(:)             ! Atmospheric profiles
  Integer (Kind=jpim), Intent (in)  :: frequencies (size(chanprof)) ! Frequency indices
  Integer (Kind=jpim), Intent (out) :: errorstatus                  ! Error return flag
  Logical (Kind=jplm), Intent (in)  :: calcemiss        (size(chanprof))         ! Switch for emmissivity calculation
  Real    (Kind=jprb), Intent (in)  :: emissivity_in    (size(chanprof))         ! Surface emmissivity 
  Real    (Kind=jprb), Intent (in)  :: emissivity_in_tl (size(chanprof))         ! Surface emmissivity 
  Type (profile_Type),        Intent (inout) :: profiles_tl     (size(profiles))    
  Type (rttov_coefs),         Intent (in)    :: coef_rttov                  ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                  ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (size(profiles)) ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles_tl (size(profiles))   
  Type (radiance_Type),       Intent (inout) :: radiance                    ! Radiances
  Type (radiance_Type),       Intent (inout) :: radiance_tl          
  Logical (Kind=jplm), optional, Intent(in)  :: lnewcld ! T = revised cloud/rain treatment; F = old treatment; Default T
  Logical (Kind=jplm), optional, Intent(in)  :: lusercfrac ! T = take av. cloud fraction from that supplied; Default F
End Subroutine
End Interface
