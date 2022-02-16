Interface
SUBROUTINE RTTOV_LAYERAVG_TL( &
            & PX1,    &
            & PX1_TL, &
            & PX2,    &
            & PX2_TL, &
            & KN1,    &
            & KN2,    &
            & PZ,     &
            & PZ_TL,  &
            & kstart, &
            & kend)
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(KIND=jpim), INTENT(IN)  :: KN1, KN2
  REAL   (KIND=jprb), INTENT(IN)  :: PX1   (KN1     ), PX1_TL(KN1     )
  REAL   (KIND=jprb), INTENT(IN)  :: PX2   (KN2     ), PX2_TL(KN2     )
  REAL   (KIND=jprb), INTENT(OUT) :: PZ    (KN2, KN1), PZ_TL (KN2, KN1)
  INTEGER(KIND=jpim), INTENT(OUT) :: kstart(KN1     ), kend  (KN1     )
End Subroutine
End Interface
