Interface
  SUBROUTINE rttov_fastem5_k( fastem_version,  &  ! Input
                              Frequency   ,    &  ! Input
                              Zenith_Angle,    &  ! Input
                              Temperature ,    &  ! Input
                              Salinity    ,    &  ! Input
                              Wind_Speed  ,    &  ! Input
                              Emissivity_k,   &  ! Input
                              Reflectivity_k, &  ! Input
                              Temperature_k,  &  ! Output
                              Salinity_k ,    &  ! Output
                              Wind_Speed_k,   &  ! Output
                              Emissivity  ,    &  ! Output
                              Reflectivity,    &  ! Output
                              Transmittance,   &  ! Input, may not be used
                              Rel_Azimuth  ,   &  ! Input, may not be used
                              Transmittance_k,&  ! Output
                              Rel_Azimuth_k )    ! Output
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
    REAL(fp),        INTENT(IN)  :: Transmittance
    REAL(fp),        INTENT(IN)  :: Rel_Azimuth
    REAL(fp),        INTENT(INOUT)  :: Emissivity_k(4), Reflectivity_k(4)
    REAL(fp),        INTENT(OUT)    :: Emissivity(4), Reflectivity(4)
    REAL(fp),        INTENT(INOUT)  :: Temperature_k,Salinity_k,Wind_Speed_k
    REAL(fp),        INTENT(INOUT)  :: Transmittance_k, Rel_Azimuth_k
End Subroutine
End Interface
