Interface
Subroutine rttov_errorreport (ErrStatus, ErrMessage, NameOfRoutine)
  Use parkind1, Only : jpim
  Implicit None
  Integer(Kind=jpim) , Intent (in) ::      ErrStatus     ! +ve => fatal error, -ve => warning 
  Character (len=*) , Intent (in) :: ErrMessage    ! ..to output 
  Character (len=*) , Intent (in) :: NameOfRoutine ! ..calling this one
End Subroutine
End Interface
