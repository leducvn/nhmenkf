Interface
SUBROUTINE rttov_calcrad_ad( &
            & chanprof,       &
            & profiles,       &
            & profiles_ad,    &
            & coeffs,         &
            & rad_skin,       &
            & rad_surfair,    &
            & rad_air,        &
            & rad_skin_ad,    &
            & rad_surfair_ad, &
            & rad_air_ad)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_type
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(profile_Type  ), INTENT(IN)    :: profiles   (:)
  TYPE(profile_Type  ), INTENT(INOUT) :: profiles_ad(size(profiles))
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof      (:)
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin      (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair   (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_air       (profiles(1)%nlevels, size(chanprof))
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin_ad   (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair_ad(size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_air_ad    (profiles(1)%nlevels, size(chanprof))
End Subroutine
End Interface
