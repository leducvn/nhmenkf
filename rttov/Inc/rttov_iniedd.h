Interface
Subroutine rttov_iniedd (&        
     & nlevels,       &! in
     & nchannels ,    &! in
     & nprofiles ,    &! in
     & lprofiles ,    &! in
     & angles    ,    &! in
     & scatt_aux)      ! inout 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :    &
       & geometry_Type        ,&
       & profile_cloud_Type   ,&
       & profile_scatt_aux 
  Use rttov_const, Only : min_ssa, ccthres
  Use parkind1, Only : jpim     ,jprb
  Implicit None
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of channels*profiles
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (geometry_Type),     Intent (in)    :: angles (nprofiles) ! Zenith angles
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux     ! Auxiliary profile variables
End Subroutine
End Interface
