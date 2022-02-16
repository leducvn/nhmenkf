Interface
subroutine rttov_read_coefs( err, coefs, opts, &
                             channels,      &
                             channels_rec,  &
                             form_coef,     &
                             form_scaer,    &
                             form_sccld,    &
                             form_pccoef,   &
                             file_coef,     &
                             file_scaer,    &
                             file_sccld,    &
                             file_pccoef,   &
                             file_id_coef,  &
                             file_id_scaer, &
                             file_id_sccld, &
                             file_id_pccoef,&
                             instrument)
#include "throw.h"
Use rttov_types, Only : rttov_coefs, rttov_options
Use parkind1, Only : jpim
Implicit None
Integer(Kind=jpim), Intent(out) :: err
Type(rttov_coefs),  Intent(out) :: coefs
Type(rttov_options),Intent(in)  :: opts
Integer(Kind=jpim), Intent(in), Optional :: channels(:)
Integer(Kind=jpim), Intent(in), Optional :: channels_rec(:)
Character(len=*),   Intent(in), Optional :: form_coef
Character(len=*),   Intent(in), Optional :: form_scaer
Character(len=*),   Intent(in), Optional :: form_sccld
Character(len=*),   Intent(in), Optional :: form_pccoef
Character(len=*),   Intent(in), Optional :: file_coef
Character(len=*),   Intent(in), Optional :: file_scaer
Character(len=*),   Intent(in), Optional :: file_sccld
Character(len=*),   Intent(in), Optional :: file_pccoef
Integer(Kind=jpim), Intent(in), Optional :: file_id_coef
Integer(Kind=jpim), Intent(in), Optional :: file_id_scaer
Integer(Kind=jpim), Intent(in), Optional :: file_id_sccld
Integer(Kind=jpim), Intent(in), Optional :: file_id_pccoef
Integer(Kind=jpim), Intent(in), Optional :: instrument(3)
End Subroutine
End Interface
