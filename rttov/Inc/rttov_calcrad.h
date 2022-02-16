Interface
SUBROUTINE rttov_calcrad( &
            & addcosmic,   &
            & chanprof,    &
            & profiles,    &
            & coeffs,      &
            & rad_cosmic,  &
            & rad_skin,    &
            & rad_surfair, &
            & rad_air)
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_Type
  USE parkind1, ONLY : jprb, jplm
  IMPLICIT NONE
  TYPE(profile_Type  ), INTENT(IN)  :: profiles(:)                                     ! profiles
  TYPE(rttov_coef    ), INTENT(IN)  :: coeffs                                          ! coefficients (Planck)
  TYPE(rttov_chanprof), INTENT(IN)  :: chanprof(:)                                     ! Array of channel indices.
  LOGICAL(KIND=jplm)  , INTENT(IN)  :: addcosmic                                       ! switch for adding cosmic background
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_cosmic (size(chanprof)                     )! cosmic background radiance
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_skin   (size(chanprof)                     )! surface skin radiance
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_surfair(size(chanprof)                     )! 2m air radiance
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_air    (profiles(1)%nlevels, size(chanprof))! level radiances
End Subroutine
End Interface
