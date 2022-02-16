Interface
Subroutine rttov_user_profile_checkinput( &
       & opts,               &
       & prof,               &
       & coefs,              &
       & ERR                 )
  Use rttov_types, Only :      &
         & rttov_coefs       , &
         & rttov_options     , &
         & profile_Type
  Use parkind1, Only : jpim
  Implicit None
  Type(rttov_options),Intent (in) :: opts    ! rttov options
  Type(profile_Type), Intent (in) :: prof    ! input profiles
  Type( rttov_coefs), Intent (in) :: coefs   ! coefficients
  Integer(Kind=jpim), Intent (out) :: ERR    ! return code
End Subroutine
End Interface
