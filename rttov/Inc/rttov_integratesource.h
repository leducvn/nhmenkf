Interface
Subroutine rttov_integratesource (&        
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & lprofiles,     &! in
     & angles,        &! in
     & scatt_aux,     &! in
     & dp,            &! in
     & dm,            &! in
     & j_do,          &! inout
     & j_up)           ! inout 
  Use rttov_types, Only :    &
       & profile_scatt_aux    ,&
       & geometry_Type 
  Use rttov_const, Only: ccthres
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use parkind1, Only : jpim     ,jprb
  Implicit None
  Integer (Kind=jpim), Intent (in) :: nlevels   ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (profile_scatt_aux), Intent (in) :: scatt_aux
  Type (geometry_Type),     Intent (in) :: angles (nprofiles) ! Zenith angles  
  Real (Kind=jprb), Intent (in)    :: dp  (nchannels,nlevels) ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)    :: dm  (nchannels,nlevels) ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout) :: j_do(nchannels,nlevels) ! Downward source terms 
  Real (Kind=jprb), Intent (inout) :: j_up(nchannels,nlevels) ! Upward source terms
End Subroutine
End Interface
