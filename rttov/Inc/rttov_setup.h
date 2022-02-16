Interface
Subroutine rttov_setup (&
       & ERR,      &! out
       & Err_unit,         &! in
       & verbosity_level,  &! in
       & opts,             &
       & coefs,            &
       & instrument,       &! in
       & channels  ,       &! in optional
       & channels_rec)      ! in optional
#include "throw.h"
  Use rttov_types, Only : &
        & rttov_coefs,           &
        & rttov_options 
  Use parkind1, Only : jpim
  Implicit None
  Type(rttov_options), Intent(in) :: opts
  Integer(Kind=jpim), Intent (in) :: Err_Unit        ! Logical error unit (<0 for default) 
  Integer(Kind=jpim), Intent (in) :: verbosity_level ! (<0 for default)
  Integer(Kind=jpim), Intent (in) :: instrument(3) ! Instrument triplet
  Integer(Kind=jpim), Optional, Intent (in) :: channels    (:)   ! list of channels to extract (channels,msat)
  Integer(Kind=jpim), Optional, Intent (in) :: channels_rec(:)
  Integer(Kind=jpim),  Intent (out) :: ERR
  Type( rttov_coefs ), Intent (out) :: coefs
End Subroutine
End Interface
