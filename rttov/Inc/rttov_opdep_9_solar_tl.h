Interface
SUBROUTINE rttov_opdep_9_solar_tl( &
            & nlayers,       &
            & chanprof,      &
            & profiles,      &
            & sun,           &
            & predictors,    &
            & predictors_tl, &
            & coef,          &
            & opdp_path,     &
            & opdp_path_tl,  &
            & opdpsun_ref)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_Type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  TYPE(rttov_chanprof) , INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type  ) , INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: sun(size(chanprof))
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(IN)    :: predictors_tl
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_tl
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  REAL(KIND=jprb)      , INTENT(IN)    :: opdpsun_ref(nlayers, size(chanprof))
End Subroutine
End Interface
