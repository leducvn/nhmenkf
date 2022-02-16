Interface
Subroutine rttov_checkinput( &
       & opts,        &
       & prof,        &
       & coef,        &! in
       & coef_pccomp, &
       & ERR          ) ! out
  Use rttov_types, Only : &
         & rttov_coef     ,&
         & rttov_options,  &
         & rttov_coef_pccomp,   &
         & profile_Type
  Use parkind1, Only : jpim
  Implicit None
  Type(rttov_options),Intent(in)  :: opts
  Type(profile_Type), Intent (inout) :: prof(:)    ! input profiles
  Type( rttov_coef ), Intent (in)    :: coef    ! coefficients
  Type( rttov_coef_pccomp ), Intent (in)  :: coef_pccomp
  Integer(Kind=jpim), Intent (out) :: ERR       ! return code
End Subroutine
End Interface
