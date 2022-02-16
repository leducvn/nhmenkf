Interface
  SUBROUTINE rttov_fastem5_tl(fastem_version,  &  ! Input
                              Frequency   ,    &  ! Input
                              Zenith_Angle,    &  ! Input
                              Temperature ,    &  ! Input
                              Salinity    ,    &  ! Input
                              Wind_Speed  ,    &  ! Input
                              Temperature_tl,  &  ! Input
                              Salinity_tl ,    &  ! Input
                              Wind_Speed_tl,   &  ! Input
                              Emissivity  ,    &  ! Output
                              Reflectivity,    &  ! Output
                              Emissivity_tl,   &  ! Output
                              Reflectivity_tl, &  ! Output
                              Transmittance,   &  ! Input, may not be used
                              Rel_Azimuth,     &  ! Input, may not be used
                              Transmittance_tl,&  ! Input, may not be used
                              Rel_Azimuth_tl )    ! Input, may not be used
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
    REAL(fp),                    INTENT(IN)     :: Temperature_tl
    REAL(fp),                    INTENT(IN)     :: Salinity_tl
    REAL(fp),                    INTENT(IN)     :: Wind_Speed_tl
    REAL(fp),        INTENT(INOUT)    :: Emissivity_tl(4), Reflectivity_tl(4)
    REAL(fp), OPTIONAL, INTENT(IN)  :: Transmittance_tl
    REAL(fp), OPTIONAL, INTENT(IN)  :: Rel_Azimuth_tl
End Subroutine
End Interface
