Interface
SUBROUTINE rttov_alloc_traj( &
            & ERR,       &
            & nprofiles, &
            & nchannels, &
            & opts,      &
            & nlevels,   &
            & coefs,     &
            & asw,       &
            & traj,      &
            & traj_tl,   &
            & traj_ad,   &
            & traj_k)
  USE parkind1, ONLY : jpim
  USE rttov_types, ONLY :  &
       & rttov_options,  &
       & rttov_coefs,    &
       & rttov_traj
  IMPLICIT NONE
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: nprofiles
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: nchannels
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: nlevels
  TYPE(rttov_options) , INTENT(IN)                 :: opts
  TYPE(rttov_coefs  ) , INTENT(IN) , TARGET        :: coefs                             ! Target attribute needed
  INTEGER(KIND=jpim)  , INTENT(OUT)                :: ERR
  TYPE(rttov_traj)    , OPTIONAL   , INTENT(INOUT) :: traj    , traj_tl, traj_ad, traj_k
  INTEGER(KIND=jpim)  , INTENT(IN)                 :: asw
End Subroutine
End Interface
