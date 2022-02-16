Interface
SUBROUTINE rttov_alloc_traj_sta (err, traj_sta, opts, coef, nlevels, nchannels, nprofiles, asw, npcscores, channels_rec)
#include "throw.h"
USE rttov_types, ONLY : &
 & rttov_traj_sta, &
 & rttov_coef,     &
 & rttov_options
USE parkind1, ONLY : &
  & jpim, jprb
IMPLICIT NONE
INTEGER(KIND=jpim),   INTENT(OUT)   :: err
TYPE(rttov_traj_sta), INTENT(INOUT) :: traj_sta
TYPE(rttov_options),  INTENT(IN)    :: opts
TYPE(rttov_coef),     INTENT(IN)    :: coef
INTEGER(KIND=jpim),   INTENT(IN)    :: nlevels
INTEGER(KIND=jpim),   INTENT(IN)    :: nchannels
INTEGER(KIND=jpim),   INTENT(IN)    :: nprofiles
INTEGER(KIND=jpim),   INTENT(IN)    :: asw
INTEGER(KIND=jpim),   INTENT(IN), OPTIONAL :: npcscores
INTEGER(KIND=jpim),   INTENT(IN), OPTIONAL :: channels_rec(:)
End Subroutine
End Interface
