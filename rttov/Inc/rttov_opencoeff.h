Interface
Subroutine rttov_opencoeff (&
       & ERR,&
       & coeffname,  &
       & file_id,    &
       & for_output, &
       & lbinary     ) 
#include "throw.h"
  Use parkind1, Only : jpim, jplm
  Implicit None
  Character (*), Intent (in) :: coeffname  ! filename of the coefficient file
  Logical(Kind=jplm), Optional, Intent (in) :: for_output     ! file access mode
  Logical(Kind=jplm), Optional, Intent (in) :: lbinary        ! if binary file wanted
  Integer(Kind=jpim), Intent(inout) :: file_id
  Integer(Kind=jpim), Intent(out) :: ERR
End Subroutine
End Interface
