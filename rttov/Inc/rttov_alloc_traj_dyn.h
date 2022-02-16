Interface
SUBROUTINE rttov_alloc_traj_dyn (err, traj_dyn, opts, nchannels, nlayers, nstreams, ncldtyp, asw)
USE rttov_types, ONLY : &
    & rttov_traj_dyn,   &
    & rttov_options
USE parkind1, ONLY : &
    & jpim
#include "throw.h"
IMPLICIT NONE
INTEGER(KIND=jpim),   INTENT(OUT)   :: err
TYPE(rttov_traj_dyn), INTENT(INOUT) :: traj_dyn
TYPE(rttov_options),  INTENT(IN)    :: opts
INTEGER(KIND=jpim),   INTENT(IN)    :: nchannels
INTEGER(KIND=jpim),   INTENT(IN)    :: nlayers
INTEGER(KIND=jpim),   INTENT(IN)    :: nstreams
INTEGER(KIND=jpim),   INTENT(IN)    :: ncldtyp
INTEGER(KIND=jpim),   INTENT(IN)    :: asw
End Subroutine
End Interface
