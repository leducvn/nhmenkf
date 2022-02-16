Interface
SUBROUTINE RTTOV_LAYERAVG_AD( &
            & PX1,    &
            & PX1_AD, &
            & PX2,    &
            & PX2_AD, &
            & KN1,    &
            & KN2,    &
            & PZ,     &
            & PZ_AD,  &
            & kstart, &
            & kend)
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(KIND=jpim), INTENT(IN)    :: KN1, KN2
  REAL   (KIND=jprb), INTENT(IN)    :: PX1   (KN1     ), PX2   (KN2     )
  REAL   (KIND=jprb), INTENT(INOUT) :: PX1_AD(KN1     ), PX2_AD(KN2     )
  REAL   (KIND=jprb), INTENT(INOUT) :: PZ    (KN2, KN1), PZ_AD (KN2, KN1)
  INTEGER(KIND=jpim), INTENT(OUT)   :: kstart(KN1     ), kend  (KN1     )
End Subroutine
End Interface
