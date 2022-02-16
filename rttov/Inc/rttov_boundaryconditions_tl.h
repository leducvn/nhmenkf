Interface
Subroutine rttov_boundaryconditions_tl (&         
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & lprofiles,     &! in
     & scatt_aux,     &! in
     & scatt_aux_tl,  &! in
     & profiles ,     &! in
     & profiles_tl ,  &! in
     & ftop,          &! in
     & ftop_tl,       &! in
     & dp,            &! out
     & dp_tl,         &! out
     & dm,            &! out
     & dm_tl)          ! out 
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  Use rttov_types, Only :    &
       & profile_Type         ,&
       & profile_scatt_aux 
  Use rttov_const, Only : ccthres 
  Use parkind1, Only : jpim     ,jprb
  Implicit none
  Integer (Kind=jpim), Intent (in) :: nlevels    ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles  ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices
  Type (profile_scatt_aux), Intent (in)    :: scatt_aux      ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_scatt_aux), Intent (in)    :: scatt_aux_tl   ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_Type),      Intent (in)    :: profiles    (nprofiles) ! Profiles on RTTOV levels
  Type (profile_Type),      Intent (in)    :: profiles_tl (nprofiles) ! Profiles on RTTOV levels
  Real (Kind=jprb), Intent  (in), dimension (nchannels)            :: ftop
  Real (Kind=jprb), Intent  (in), dimension (nchannels)            :: ftop_tl
  Real (Kind=jprb), Intent (out), dimension (nchannels,nlevels) :: dp   , dm
  Real (Kind=jprb), Intent (out), dimension (nchannels,nlevels) :: dp_tl, dm_tl
End Subroutine
End Interface
