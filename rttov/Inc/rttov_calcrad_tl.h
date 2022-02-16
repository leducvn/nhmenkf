Interface
SUBROUTINE rttov_calcrad_tl( &
            & chanprof,       &
            & profiles,       &
            & profiles_tl,    &
            & coeffs,         &
            & rad_skin,       &
            & rad_surfair,    &
            & rad_air,        &
            & rad_skin_tl,    &
            & rad_surfair_tl, &
            & rad_air_tl)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_type
  USE parkind1, ONLY : jprb
  IMPLICIT NONE
  TYPE(profile_Type  ), INTENT(IN)  :: profiles   (:)
  TYPE(profile_Type  ), INTENT(IN)  :: profiles_tl(size(profiles))
  TYPE(rttov_coef    ), INTENT(IN)  :: coeffs
  TYPE(rttov_chanprof), INTENT(IN)  :: chanprof      (:)
  REAL(KIND=jprb)     , INTENT(IN)  :: rad_skin      (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)  :: rad_surfair   (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)  :: rad_air       (profiles(1)%nlevels, size(chanprof))
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_skin_tl   (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_surfair_tl(size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_air_tl    (profiles(1)%nlevels, size(chanprof))
End Subroutine
End Interface
