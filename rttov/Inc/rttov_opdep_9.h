Interface
SUBROUTINE rttov_opdep_9( &
            & nlayers,    &
            & chanprof,   &
            & predictors, &
            & aux,        &
            & coef,       &
            & opdp_path,  &
            & opdp_ref)
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                          ! Number of pressure levels
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)                      ! Channel indices
  TYPE(predictors_Type), INTENT(IN)    :: predictors                       ! Predictors
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                             ! Coefficients
  TYPE(profile_aux    ), INTENT(IN)    :: aux                              ! auxillary profiles informations
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path                        ! optical depths
  REAL(KIND=jprb)      , INTENT(OUT)   :: opdp_ref(nlayers, size(chanprof))! layer optical depth
End Subroutine
End Interface
