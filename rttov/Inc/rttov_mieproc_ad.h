Interface
Subroutine rttov_mieproc_ad (&        
     & nlevels,           &! in
     & nchannels,         &! in
     & nprofiles,         &! in
     & nprofilesad,       &! in
     & frequencies,       &! in
     & lprofiles,         &! in
     & profiles,          &! in
     & coef_scatt,        &! in
     & scatt_aux,         &! inout
     & scatt_aux_ad)       ! inout 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :    &
       & profile_scatt_aux    ,&
       & profile_Type   ,&
       & rttov_scatt_coef     
  Use rttov_const, Only: adk_adjoint, adk_k, ccthres
  Use parkind1, Only : jpim     ,jprb
  Implicit None
  Integer (Kind=jpim), Intent (in) :: nlevels                 ! Number oflevels
  Integer (Kind=jpim), Intent (in) :: nchannels               ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofiles               ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad             ! Number of profiles in adjoint variables
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Integer (Kind=jpim), Intent (in) :: lprofiles   (nchannels) ! Profile indices
  Type (profile_Type),      Intent (in)    :: profiles    (nprofiles)  
  Type (rttov_scatt_coef),  Intent (in)    :: coef_scatt               ! RTTOV_SCATT Coefficients
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux                ! Auxiliary profile variables
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux_ad             ! Auxiliary profile variables
End Subroutine
End Interface
