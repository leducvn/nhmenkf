Interface
SUBROUTINE rttov_calcrad_k( &
            & chanprof,      &
            & profiles,      &
            & profiles_k,    &
            & coeffs,        &
            & rad_skin,      &
            & rad_surfair,   &
            & rad_air,       &
            & rad_skin_k,    &
            & rad_surfair_k, &
            & rad_air_k)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_type
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type  ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type  ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin     (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair  (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_air      (profiles(1)%nlevels, size(chanprof))
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin_k   (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair_k(size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_air_k    (profiles(1)%nlevels, size(chanprof))
End Subroutine
End Interface
