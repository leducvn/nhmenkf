Interface
Subroutine rttov_coeffname (ERR, instrument, coeffname, type, lbinary)
#include "throw.h"
  Use parkind1, Only : jpim, jplm
  Implicit None
  Character(len=*), Intent(in) :: type
  Integer(Kind=jpim), Intent (in) :: instrument(3)  ! (platform, sat_id, inst) numbers
  Logical(Kind=jplm), Optional, Intent (in) :: lbinary  ! if binary file wanted
  Integer(Kind=jpim), Intent (out)       :: ERR! return code
  Character (*), Intent (out) :: coeffname  ! filename of the coefficient file
End Subroutine
End Interface
