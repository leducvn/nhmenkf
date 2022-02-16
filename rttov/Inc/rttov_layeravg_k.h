Interface
SUBROUTINE RTTOV_LAYERAVG_K( &
            & PX1,    &
            & PX1_K,  &
            & PX2,    &
            & PX2_K,  &
            & KN1,    &
            & KN2,    &
            & PZ,     &
            & PZ_K,   &
            & kstart, &
            & kend)
  USE parkind1, ONLY : jpim, jprb
  IMPLICIT NONE
  INTEGER(KIND=jpim), INTENT(IN)    :: KN1, KN2
  REAL   (KIND=jprb), INTENT(IN)    :: PX1   (KN1     ), PX2  (KN2     )
  REAL   (KIND=jprb), INTENT(INOUT) :: PX1_K (KN1     ), PX2_K(KN2     )
  REAL   (KIND=jprb), INTENT(INOUT) :: PZ    (KN2, KN1), PZ_K (KN2, KN1)
  INTEGER(KIND=jpim), INTENT(OUT)   :: kstart(KN1     ), kend (KN1     )
End Subroutine
End Interface
