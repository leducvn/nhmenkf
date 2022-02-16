Interface
Subroutine rttov_integratesource_ad (&        
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & nprofilesad,   &! in
     & lprofiles,     &! in
     & angles,        &! in
     & scatt_aux,     &! in
     & scatt_aux_ad,  &! inout
     & dp,            &! in
     & dp_ad,         &! inout
     & dm,            &! in
     & dm_ad,         &! inout
     & j_do,          &! inout
     & j_do_ad,       &! inout
     & j_up,          &! inout
     & j_up_ad)        ! inout 
  Use rttov_types, Only :    &
       & profile_scatt_aux    ,&
       & geometry_Type 
  Use rttov_const, Only: adk_adjoint, adk_k, ccthres
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use parkind1, Only : jpim     ,jprb
  Implicit None
  Integer (Kind=jpim), Intent (in) :: nlevels      ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles    ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels    ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (profile_scatt_aux), Intent (in)    :: scatt_aux          ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_scatt_aux), Intent (inout) :: scatt_aux_ad       ! Auxiliary profile variables for RTTOV_SCATT
  Type (geometry_Type),     Intent (in)    :: angles (nprofiles) ! Zenith angles  
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dp       ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dm       ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do     ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up     ! Upward source terms
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: dp_ad    ! D+ for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: dm_ad    ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do_ad  ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up_ad  ! Upward source terms
End Subroutine
End Interface
