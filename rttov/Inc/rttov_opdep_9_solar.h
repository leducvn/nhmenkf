Interface
SUBROUTINE rttov_opdep_9_solar( &
            & nlayers,     &
            & chanprof,    &
            & profiles,    &
            & sun,         &
            & predictors,  &
            & coef,        &
            & opdp_path,   &
            & opdpsun_ref)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_Type
  USE parkind1, ONLY : jpim, jprb, jplm
  IMPLICIT NONE
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                             ! Number of pressure levels
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)                         ! Channel indices
  TYPE(profile_Type   ), INTENT(IN)    :: profiles(:)                         ! Atmospheric profiles
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: sun(size(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors                          ! Predictors
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                                ! Coefficients
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path
  REAL(KIND=jprb)      , INTENT(OUT)   :: opdpsun_ref(nlayers, size(chanprof))! layer optical depth
End Subroutine
End Interface
