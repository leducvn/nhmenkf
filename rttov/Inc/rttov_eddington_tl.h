Interface
Subroutine rttov_eddington_tl ( &
     & nlevels,        &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & lprofiles,         &! in
     & angles,            &! in
     & profiles,          &! in
     & profiles_tl,       &! in
     & scatt_aux,         &! in
     & scatt_aux_tl,      &! in
     & cld_bt,            &! out
     & cld_bt_tl)          ! out 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :      &
       & geometry_Type        ,&
       & profile_Type         ,&
       & profile_cloud_Type   ,&
       & profile_scatt_aux    
  Use rttov_const, Only:       &
       & tcosmic, ccthres
  Use parkind1, Only : jpim     ,jprb
  Implicit None
  Integer (Kind=jpim), Intent (in) :: nlevels            ! Number of NWP-levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (geometry_Type),       Intent (in)    :: angles       (nprofiles) ! Zenith angles
  Type (profile_Type),        Intent (in)    :: profiles     (nprofiles) ! Profiles 
  Type (profile_Type),        Intent (in)    :: profiles_tl  (nprofiles) ! Profiles 
  Type (profile_scatt_aux),   Intent (in)    :: scatt_aux                ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_scatt_aux),   Intent (in)    :: scatt_aux_tl             ! Auxiliary profile variables for RTTOV_SCATT
  Real (Kind=jprb),           Intent (out)   :: cld_bt       (nchannels) ! Radiances
  Real (Kind=jprb),           Intent (out)   :: cld_bt_tl    (nchannels) ! Radiances
End Subroutine
End Interface
