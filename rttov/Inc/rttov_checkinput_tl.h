Interface
Subroutine rttov_checkinput_tl( &
       & opts,        &
       & prof,        &
       & prof_tl,     &
       & coef,        &
       & coef_pccomp) ! out
  Use rttov_types, Only : &
         & rttov_coef     ,&
         & rttov_options,  &
         & rttov_coef_pccomp   ,&
         & profile_Type
  Implicit None
  Type(rttov_options),Intent(in)  :: opts
  Type(profile_Type), Intent (in) :: prof(:)    ! input profiles
  Type(profile_Type), Intent (inout) :: prof_tl(:)    ! input profiles
  Type( rttov_coef ), Intent (in) :: coef    ! coefficients
  Type( rttov_coef_pccomp ), Intent (in) :: coef_pccomp  ! coefficients
End Subroutine
End Interface
