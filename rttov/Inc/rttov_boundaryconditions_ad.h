Interface
Subroutine rttov_boundaryconditions_ad (& 
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & nprofilesad,   &! in
     & lprofiles,     &! in
     & scatt_aux,     &! in
     & scatt_aux_ad,  &! inout
     & profiles ,     &! in
     & profiles_ad ,  &! inout
     & ftop,          &! in
     & ftop_ad,       &! inout
     & dp,            &! out
     & dp_ad,         &! inout
     & dm,            &! out
     & dm_ad)          ! inout 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :    &
       & profile_Type         ,&
       & profile_scatt_aux 
  Use rttov_const, Only: adk_adjoint, adk_k, ccthres
  Use parkind1, Only : jpim     ,jprb
  Implicit none
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad           ! Number of profiles in adjoint variables
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (profile_scatt_aux), Intent    (in) :: scatt_aux               ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux_ad            ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_Type),      Intent    (in) :: profiles    (nprofiles) ! Profiles on RTTOV levels
  Type (profile_Type),      Intent (inout) :: profiles_ad (nprofilesad) ! Profiles on RTTOV levels
  Real (Kind=jprb), Intent    (in), dimension (nchannels)            :: ftop
  Real (Kind=jprb), Intent (inout), dimension (nchannels)            :: ftop_ad
  Real (Kind=jprb), Intent   (out), dimension (nchannels,nlevels) :: dp   , dm
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: dp_ad, dm_ad
End Subroutine
End Interface
