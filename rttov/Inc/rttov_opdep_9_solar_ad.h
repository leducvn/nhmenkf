Interface
SUBROUTINE rttov_opdep_9_solar_ad( &
            & nlayers,       &
            & chanprof,      &
            & profiles,      &
            & sun,           &
            & predictors,    &
            & predictors_ad, &
            & coef,          &
            & opdp_path,     &
            & opdp_path_ad,  &
            & opdpsun_ref)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type   ), INTENT(IN)    :: profiles(:)
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: sun(size(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_ad
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_ad
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  REAL(KIND=jprb)      , INTENT(IN)    :: opdpsun_ref(nlayers, size(chanprof))
End Subroutine
End Interface
