Interface
  SUBROUTINE rttov_fastem5( fastem_version, &  ! Input
                              Frequency   , &  ! Input
                              Zenith_Angle, &  ! Input
                              Temperature , &  ! Input
                              Salinity    , &  ! Input
                              Wind_Speed  , &  ! Input
                              Emissivity  , &  ! Output
                              Reflectivity, &  ! Output
                              Transmittance,&  ! Input, may not be used
                              Rel_Azimuth  )   ! Input, may not be used
    USE mod_rttov_fastem5_coef, ONLY : FresnelVariables_type, PermittivityVariables_type,&
        fp, ZERO, POINT_5, ONE, TWO, THREE,PI,DEGREES_TO_RADIANS,transmittance_limit_lower,&
        transmittance_limit_upper, e0_4, e0_5, min_f, max_f, min_wind, max_wind, A_COEF, Lcoef4, Lcoef5,&
        Scoef, t_c4, t_c5, b_coef, FR_COEFF,x,y
    USE parkind1, ONLY : jpim
    IMPLICIT NONE
    INTEGER(jpim),   INTENT(IN)     :: fastem_version
    REAL(fp),        INTENT(IN)     :: Frequency
    REAL(fp),        INTENT(IN)     :: Zenith_Angle
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Wind_Speed
    REAL(fp),        INTENT(OUT)    :: Emissivity(4), Reflectivity(4)
    REAL(fp), INTENT(IN)  :: Transmittance
    REAL(fp), INTENT(IN)  :: Rel_Azimuth
End Subroutine
End Interface
