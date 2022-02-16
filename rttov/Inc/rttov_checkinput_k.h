Interface
Subroutine rttov_checkinput_k( &
       & opts,        &
       & chanprof,    &
       & prof,        &
       & prof_k,      &
       & coef) ! out
  Use rttov_types, Only : &
         & rttov_chanprof, &
         & rttov_coef     ,&
         & rttov_options,  &
         & profile_Type
  Implicit None
  Type(rttov_options),Intent(in)  :: opts
  Type(rttov_chanprof),Intent(in) :: chanprof(:)
  Type(profile_Type), Intent (in) :: prof(:)    ! input profiles
  Type(profile_Type), Intent (inout) :: prof_k(:)    ! input profiles
  Type( rttov_coef ), Intent (in) :: coef    ! coefficients
End Subroutine
End Interface
