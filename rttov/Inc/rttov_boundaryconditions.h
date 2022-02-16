Interface
Subroutine rttov_boundaryconditions (&
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & lprofiles,     &! in
     & scatt_aux,     &! in
     & profiles ,     &! in
     & ftop,          &! in
     & dp,            &! out
     & dm)             ! out 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :    &
       & profile_Type       ,&
       & profile_scatt_aux 
  Use rttov_const, Only : ccthres
  Use parkind1, Only : jpim     ,jprb
  Implicit none
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (profile_scatt_aux), Intent (in) :: scatt_aux 
  Type (profile_Type), Intent (in) :: profiles (nprofiles)
  Real (Kind=jprb), Intent (in), dimension (nchannels) :: ftop
  Real (Kind=jprb), Intent (out), dimension (nchannels,nlevels) :: dp, dm
End Subroutine
End Interface
