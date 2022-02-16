!
  SUBROUTINE rttov_fastem5_ad(fastem_version,  &  ! Input
                              Frequency   ,    &  ! Input
                              Zenith_Angle,    &  ! Input
                              Temperature ,    &  ! Input
                              Salinity    ,    &  ! Input
                              Wind_Speed  ,    &  ! Input
                              Emissivity_ad,   &  ! Input
                              Reflectivity_ad, &  ! Input
                              Temperature_ad,  &  ! Output
                              Salinity_ad ,    &  ! Output
                              Wind_Speed_ad,   &  ! Output
                              Emissivity  ,    &  ! Output
                              Reflectivity,    &  ! Output
                              Transmittance,   &  ! Input, may not be used
                              Rel_Azimuth  ,   &  ! Input, may not be used
                              Transmittance_ad,&  ! Output
                              Rel_Azimuth_ad )    ! Output
  ! Description:
  ! Adjoint of rttov_fastem5
  ! To compute FASTEM-4 emissivities and reflectivities.
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2009, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  ! An improved fast microwave sea surface emissivity model, FASTEM4
  ! Liu, Q., S. English, F. Weng, 2009: Report in prepare
  !
  ! It is an extension of the FASTEM-3 English 2003.
  ! http://www.metoffice.com/research/interproj/nwpsaf/rtm/evalfastems.pdf
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       27/08/2009  New F90 code (Q. Liu)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
  ! Imported Parameters:
    USE mod_rttov_fastem5_coef, ONLY : FresnelVariables_type, PermittivityVariables_type,&
        fp, ZERO, POINT_5, ONE, TWO, THREE,PI,DEGREES_TO_RADIANS,transmittance_limit_lower,&
        transmittance_limit_upper, e0_4, e0_5, min_f, max_f, min_wind, max_wind, A_COEF, Lcoef4, Lcoef5,&
        Scoef, t_c4, t_c5, b_coef, FR_COEFF,x,y
    USE parkind1, ONLY : jpim
    IMPLICIT NONE
    ! Arguments
    INTEGER(jpim),   INTENT(IN)     :: fastem_version
    REAL(fp),        INTENT(IN)     :: Frequency
    REAL(fp),        INTENT(IN)     :: Zenith_Angle
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Wind_Speed
    REAL(fp),        INTENT(IN)  :: Transmittance
    REAL(fp),        INTENT(IN)  :: Rel_Azimuth
    REAL(fp),        INTENT(INOUT)  :: Emissivity_ad(4), Reflectivity_ad(4)
    REAL(fp),        INTENT(OUT)    :: Emissivity(4), Reflectivity(4)
    REAL(fp),        INTENT(INOUT)  :: Temperature_ad,Salinity_ad,Wind_Speed_ad
    REAL(fp),        INTENT(INOUT)  :: Transmittance_ad, Rel_Azimuth_ad

!INTF_END

  !local variable

    REAL(fp) :: e0
    REAL(fp) :: cos_z, Foam_Cover
    REAL(fp) :: scor, small_corr, Azimuth_Emi(4),RV_Fresnel,RH_Fresnel
    REAL(fp) :: Ev,Eh,RvL,RhL,RvS,RhS
    REAL(fp) :: zreflmod_v,zreflmod_h,zrough_v,zrough_h
    INTEGER :: i, j, L, m

    REAL(fp) :: Foam_Rv,Foam_Rh
    ! Local variables for the permittivity model
    REAL(fp) :: einf, sigma25
    REAL(fp) :: tau1, tau2, es, e1
    REAL(fp) :: perm_Real, perm_imag
!    REAL(fp) ::
    TYPE(PermittivityVariables_type) :: iVar
    COMPLEX( fp ) :: Permittivity

    ! Local variables for Fresnel reflectance
    COMPLEX(fp) :: zRv ! Vertical
    COMPLEX(fp) :: zRh ! Horizontal
    TYPE(FresnelVariables_type) :: frVar

    ! Local variables for small-scale
    REAL(fp) :: windspeed, freq_S
    ! Local variables for large-scale
    REAL(fp) :: seczen, zc(12)

    ! Local variables for including transmittance
    REAL(fp) :: variance,varm,opdpsfc,zx(9)

    ! Foam reflectance
    REAL(fp) :: Fv, Fh, Foam_ref

    ! Local variables for azimuth angle
    REAL(fp) :: asc(4,3), fre_c, phi
    ! Local arrays to hold FASTEM-4/5 coefs
    REAL(fp) :: Lcoef(size(Lcoef5)), t_c(size(t_c5))

  ! ad part
    REAL(fp) :: scor_ad, small_corr_ad, Azimuth_Emi_ad(4)
    REAL(fp) :: Foam_Cover_ad,RV_Fresnel_ad,RH_Fresnel_ad
    REAL(fp) :: Ev_ad,Eh_ad,RvL_ad,RhL_ad,RvS_ad,RhS_ad
    COMPLEX( fp ) ::  Permittivity_AD

    REAL( fp ) :: ac_ad, sc_ad, phi_ad,wind10_ad,zreflmod_v_ad,zreflmod_h_ad

    REAL(fp) :: sigma_AD
    REAL(fp) :: einf_AD, e1_AD, es_AD, ce1_AD, ces_AD
    REAL(fp) :: t_AD, t_sq_AD, t_cu_AD, S_AD, beta_AD, delta_AD, sigma25_AD
    REAL(fp) :: tau1_AD, tau2_AD, ctau1_AD, ctau2_AD, f1_AD, f2_AD, del1_AD, del2_AD
    REAL(fp) :: perm_Real_AD, perm_imag_AD
    ! Local variables
    COMPLEX(fp) :: z1_AD, z2_AD
    COMPLEX(fp) :: zRv_AD           ! Vertical
    COMPLEX(fp) :: zRh_AD           ! Horizontal
    REAL(fp)    :: rzRv_AD,izRv_AD  ! Vertical
    REAL(fp)    :: rzRh_AD,izRh_AD  ! Horizontal
    COMPLEX(fp) :: denom
    REAL(fp) :: variance_ad,varm_ad,opdpsfc_ad,zx_ad(9),zrough_v_ad,zrough_h_ad

    IF (fastem_version == 4) THEN
      e0 = e0_4
      Lcoef = Lcoef4
      t_c = t_c4
    ELSE
      e0 = e0_5
      Lcoef = Lcoef5
      t_c = t_c5
    ENDIF
    cos_z = cos( Zenith_Angle*DEGREES_TO_RADIANS )

  ! Permittivity Calculation
  ! ------------------------
    !1.2 calculate permittivity using double-debye formula
    !-----------------------------------------------------
    !Set values for temperature polynomials (convert from kelvin to celsius)
    iVar%t = Temperature - 273.15_fp
    iVar%t_sq = iVar%t * iVar%t     !quadratic
    iVar%t_cu = iVar%t_sq * iVar%t  !cubic
    iVar%S = Salinity
    !-----------------------------------------------------
    !1.2 Pure or fresh water
    !-----------------------------------------------------
    einf = A_COEF(0) + A_COEF(1)*iVar%t
    es   = A_COEF(2) + A_COEF(3)*iVar%t  + A_COEF(4)*iVar%t_sq + A_COEF(5)*iVar%t_cu
    e1   = A_COEF(9) + A_COEF(10)*iVar%t + A_COEF(11)*iVar%t_sq
    tau1 = A_COEF(15) + A_COEF(16)*iVar%t + A_COEF(17)*iVar%t_sq + A_COEF(18)*iVar%t_cu
    tau2 = A_COEF(22) + A_COEF(23)*iVar%t + A_COEF(24)*iVar%t_sq + A_COEF(25)*iVar%t_cu

    iVar%es_k = es
    iVar%e1_k = e1
    iVar%tau1_k = tau1
    iVar%tau2_k = tau2
    perm_imag = ZERO

    IF( iVar%S > ZERO ) THEN
      iVar%delta = 25.0_fp - iVar%t
      iVar%beta  = A_COEF(29) +A_COEF(30)*iVar%delta +A_COEF(31)*iVar%delta**2  &
            + iVar%S*(A_COEF(32) +A_COEF(33)*iVar%delta +A_COEF(34)*iVar%delta**2)
      sigma25 = iVar%S*(A_COEF(35) +A_COEF(36)*iVar%S +A_COEF(37)*iVar%S**2  &
              +A_COEF(38)*iVar%S**3)
      iVar%sigma = sigma25*exp(-iVar%delta*iVar%beta)

      iVar%ces = ONE + iVar%S*(A_COEF(6) + A_COEF(7)*iVar%S + A_COEF(8)*iVar%t )
      iVar%ce1 = ONE + iVar%S*(A_COEF(12) + A_COEF(13)*iVar%S +A_COEF(14)*iVar%t )
      iVar%ctau1 = ONE + iVar%S*(A_COEF(19) +A_COEF(20)*iVar%t + A_COEF(21)*iVar%t_sq)
      iVar%ctau2 = ONE + iVar%S*(A_COEF(26) + A_COEF(27)*iVar%t + A_COEF(28)*iVar%S**2 )
      es = iVar%es_k * iVar%ces
      e1 = iVar%e1_k * iVar%ce1
      tau1 = iVar%tau1_k * iVar%ctau1
      tau2 = iVar%tau2_k * iVar%ctau2
      perm_imag = -iVar%sigma/(TWO*PI*e0*Frequency)
    END IF
    !Define two relaxation frequencies, f1 and f2
    iVar%f1 = Frequency*tau1
    iVar%f2 = Frequency*tau2
    iVar%del1 = es - e1
    iVar%del2 = e1 - einf

    perm_Real = einf + iVar%del1/(ONE + iVar%f1**2) + iVar%del2/(ONE + iVar%f2**2)
    perm_imag = -perm_imag + iVar%del1*iVar%f1/(ONE + iVar%f1**2)  &
              + iVar%del2*iVar%f2/(ONE + iVar%f2**2)
    Permittivity = Cmplx(perm_Real,-perm_imag,fp)

  ! Compute Fresnel reflectance code, adopted from Masahiro Kazumori, JMA
  !
    ! Compute the complex reflectivity components
    frVar%z1 = SQRT(permittivity - ONE + (cos_z*cos_z))
    frVar%z2 = permittivity * cos_z
    zRh = (cos_z  -frVar%z1) / (cos_z  +frVar%z1)
    zRv = (frVar%z2-frVar%z1) / (frVar%z2+frVar%z1)

    ! The square of the vertical abs value
    frVar%rzRv = REAL(zRv,fp)
    frVar%izRv = AIMAG(zRv)
    Rv_Fresnel = frVar%rzRv**2 + frVar%izRv**2

    ! The square of the horizontal abs value
    frVar%rzRh = REAL(zRh,fp)
    frVar%izRh = AIMAG(zRh)
    Rh_Fresnel = frVar%rzRh**2 + frVar%izRh**2


  ! Apply small-scale correction
  ! --------------------------------
    windspeed = Wind_Speed
    IF( windspeed < min_wind ) windspeed = min_wind
    IF( windspeed > max_wind ) windspeed = max_wind

    freq_S = Frequency
    IF( freq_S < min_f ) freq_S = min_f
    IF( freq_S > max_f ) freq_S = max_f

    scor = Scoef(1) *windspeed*freq_S +Scoef(2) *windspeed*freq_S**2 &
           + Scoef(3) *windspeed**2* freq_S +Scoef(4) *windspeed**2* freq_S**2 &
           + Scoef(5) *windspeed**2 /freq_S +Scoef(6) *windspeed**2 /freq_S**2 &
           + Scoef(7) *windspeed + Scoef(8) *windspeed**2

    small_corr = exp(-scor*cos_z*cos_z )
    RvS = Rv_Fresnel * small_corr

    RhS = Rh_Fresnel * small_corr

  ! Large Scale Correction Calculation
  ! ----------------------------------
    seczen = ONE/cos_z
    ! compute fitting coefficients for a given frequency
    DO j = 1, 12
      zc(j) = Lcoef(j*3-2) + Lcoef(j*3-1)*frequency + Lcoef(j*3)*frequency**2
    END DO


    RvL = zc(1) + zc(2)*seczen + zc(3)*seczen**2 + zc(4)*Wind_Speed &
      + zc(5)*Wind_Speed**2 + zc(6)*Wind_Speed*seczen
    RhL = zc(7) + zc(8)*seczen + zc(9)*seczen**2 + zc(10)*Wind_Speed &
      + zc(11)*Wind_Speed**2 + zc(12)*Wind_Speed*seczen

    ! change foam coverage back to FASTEM1,2,3 for FASTEM5
    IF (fastem_version == 4) THEN
      ! Compute foam coverage after Tang, 1974
      Foam_Cover = 7.75E-06_fp * Wind_Speed ** 3.231_fp
    ELSE
      ! Monahan et al., 1986 without surface stability term
      Foam_Cover = 1.95E-05_fp * Wind_Speed ** 2.55_fp
    END IF
    
  ! The foam vertical and horizontal reflectanc codes, adopted from Masahiro Kazumori, JMA
  ! ----------------------------------
    Fv = ONE + Zenith_Angle*(FR_COEFF(1)+ Zenith_Angle*(FR_COEFF(2)  &
       + Zenith_Angle*FR_COEFF(3))) + FR_COEFF(4)*Zenith_Angle**10
    Foam_Rv = FR_COEFF(5)
    Fh = ONE + Zenith_Angle*(FR_COEFF(6) +  Zenith_Angle*(FR_COEFF(7)  &
       + Zenith_Angle*FR_COEFF(8)))
    Foam_Rh = ONE + FR_COEFF(9)*Fh

  ! Added frequency dependence derived from Stogry model, 1971
  ! ----------------------------------
    Foam_ref = 0.4_fp * exp(-0.05_fp*Frequency )
    Foam_Rv = Foam_Rv * Foam_ref
    Foam_Rh = Foam_Rh * Foam_ref

    Ev = (ONE-Foam_Cover)*(ONE - RvS + RvL) + Foam_Cover*(ONE-Foam_Rv)
    Eh = (ONE-Foam_Cover)*(ONE - RhS + RhL) + Foam_Cover*(ONE-Foam_Rh)
    Emissivity(1) = Ev
    Emissivity(2) = Eh

    zreflmod_v = ONE
    zreflmod_h = ONE

  ! correction for anisotropic downward radiation, adopted from the FASTEM3
  ! ----------------------------------
    IF( Transmittance > transmittance_limit_lower .and. Transmittance < transmittance_limit_upper) THEN
        !Using the Cox and Munk model to compute slope variance
        variance = 0.00512_fp * Wind_Speed + 0.0030_fp
        varm     = variance * t_c(43)
        variance = varm * ( t_c(44) * Frequency + t_c(45) )

        IF ( variance >= varm ) THEN
          variance = varm
        END IF
        IF ( variance <= ZERO  ) THEN
          variance = ZERO
        END IF
        !Compute surface to space optical depth
        opdpsfc = -log(Transmittance ) * cos_z
        !Define nine predictors for the effective angle calculation
        zx(1) = ONE
        zx(2) = variance
        zx(4) = ONE / cos_z
        zx(3) = zx(2) * zx(4)
        zx(5) = zx(3) * zx(3)
        zx(6) = zx(4) * zx(4)
        zx(7) = zx(2) * zx(2)
        zx(8) = log(opdpsfc)
        zx(9) = zx(8) * zx(8)

        zrough_h = ONE
        zrough_v = ONE
        DO i = 1, 7
           j = i-1
          !Switched h to v Deblonde SSMIS june 7, 2001
          zrough_h = zrough_h + zx(i) *(t_c(1+j*3) + zx(8)*t_c(2+j*3) + zx(9)*t_c(3+j*3))
          zrough_v = zrough_v + zx(i) *(t_c(22+j*3)+ zx(8)*t_c(23+j*3)+ zx(9)*t_c(24+j*3))

        END DO
        zreflmod_v = (ONE-Transmittance ** zrough_v)/(ONE-Transmittance )
        zreflmod_h = (ONE-Transmittance ** zrough_h)/(ONE-Transmittance )

    END IF

  ! azimuthal component, emissivity for full Stokes parameters
  ! --------------------------------
    Azimuth_Emi = ZERO
    IF( abs(Rel_Azimuth) <= 360.0_fp ) THEN
      Fre_C = ZERO
      IF( Frequency >= min_f .or. Frequency <= max_f ) THEN
        DO i = 1, 8
          IF( Frequency >= x(i) .and. Frequency < x(i+1) ) THEN
            Fre_C = y(i) + (y(i+1)-y(i))/(x(i+1)-x(i))*(Frequency-x(i))
          END IF
        END DO
      END IF

      phi = Rel_Azimuth * DEGREES_TO_RADIANS

      DO m = 1, 3
        L = 10*(m-1)
        asc(1,m) = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen &
          +b_coef(L+4)*seczen*Frequency +b_coef(L+5)*wind_speed  &
          +b_coef(L+6)*wind_speed*Frequency +b_coef(L+7)*wind_speed**2  &
          +b_coef(L+8)*Frequency*wind_speed**2 +b_coef(L+9)*wind_speed*seczen &
          +b_coef(L+10)*wind_speed*seczen*Frequency
        Azimuth_Emi(1) = Azimuth_Emi(1) + asc(1,m)*cos(m*phi)

        L = 10*(m-1) + 30
        asc(2,m) = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen &
          +b_coef(L+4)*seczen*Frequency +b_coef(L+5)*wind_speed  &
          +b_coef(L+6)*wind_speed*Frequency +b_coef(L+7)*wind_speed**2  &
          +b_coef(L+8)*Frequency*wind_speed**2 +b_coef(L+9)*wind_speed*seczen &
          +b_coef(L+10)*wind_speed*seczen*Frequency
        Azimuth_Emi(2) = Azimuth_Emi(2) + asc(2,m)*cos(m*phi)

        L = 10*(m-1) + 60
        asc(3,m) = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen &
          +b_coef(L+4)*seczen*Frequency +b_coef(L+5)*wind_speed &
          +b_coef(L+6)*wind_speed*Frequency +b_coef(L+7)*wind_speed**2  &
          +b_coef(L+8)*Frequency*wind_speed**2 +b_coef(L+9)*wind_speed*seczen &
          +b_coef(L+10)*wind_speed*seczen*Frequency
        Azimuth_Emi(3) = Azimuth_Emi(3) + asc(3,m)*sin(m*phi)

        L = 10*(m-1) + 90
        asc(4,m) = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen &
          +b_coef(L+4)*seczen*Frequency +b_coef(L+5)*wind_speed  &
          +b_coef(L+6)*wind_speed*Frequency +b_coef(L+7)*wind_speed**2  &
          +b_coef(L+8)*Frequency*wind_speed**2 +b_coef(L+9)*wind_speed*seczen &
          +b_coef(L+10)*wind_speed*seczen*Frequency
        Azimuth_Emi(4) = Azimuth_Emi(4) + asc(4,m)*sin(m*phi)

      END DO
      Azimuth_Emi = Azimuth_Emi * Fre_C

    END IF

    Emissivity(1) = Emissivity(1) + Azimuth_Emi(1)
    Emissivity(2) = Emissivity(2) + Azimuth_Emi(2)
    Emissivity(3) = Azimuth_Emi(3)
    Emissivity(4) = Azimuth_Emi(4)
    Reflectivity(1)  = zreflmod_v * (ONE-Emissivity(1))
    Reflectivity(2)  = zreflmod_h * (ONE-Emissivity(2))
    ! Reflectivities not computed for 3rd or 4th elements of Stokes vector, 
    ! as never used subsequently, as atmospheric source term = zero.
    Reflectivity(3:4)  = ZERO

  !  ****** End of Forward Part ******
  !  ****** start of adjoint code ******

    wind10_ad = ZERO
    phi_ad = 0.0_fp
    zreflmod_v_ad = ZERO
    zreflmod_h_ad = ZERO

    Emissivity_ad(2) = Emissivity_ad(2) -zreflmod_h*Reflectivity_ad(2)
    zreflmod_h_ad = zreflmod_h_ad + Reflectivity_ad(2)* (ONE-Emissivity(2))

    Emissivity_ad(1) = Emissivity_ad(1) -zreflmod_v*Reflectivity_ad(1)
    zreflmod_v_ad = zreflmod_v_ad + Reflectivity_ad(1)* (ONE-Emissivity(1))
    Azimuth_Emi_ad(4) = Emissivity_ad(4)
    Azimuth_Emi_ad(3) = Emissivity_ad(3)
    Azimuth_Emi_ad(2) = Emissivity_ad(2)
    Azimuth_Emi_ad(1) = Emissivity_ad(1)

    ! azimuthal component
    ! --------------------------------
  IF( abs(Rel_Azimuth) <= 360.0_fp ) THEN

    seczen = ONE/cos_z
    phi = Rel_Azimuth * DEGREES_TO_RADIANS

    Azimuth_Emi_ad = Azimuth_Emi_ad * Fre_C

    DO m = 1, 3
      phi_ad = phi_ad +m*Azimuth_Emi_ad(4)*asc(4,m)*cos(m*phi)
      sc_ad = Azimuth_Emi_ad(4)*sin(m*phi)
      L = 10*(m-1) + 90
      wind10_ad = wind10_ad + ( b_coef(L+5) +b_coef(L+6)*Frequency )*sc_ad  &
          +( TWO*(b_coef(L+7) +b_coef(L+8)*Frequency)*wind_speed &
          +b_coef(L+9)*seczen +b_coef(L+10)*seczen*Frequency )*sc_ad

      phi_ad = phi_ad +m*Azimuth_Emi_ad(3)*asc(3,m)*cos(m*phi)
      sc_ad = Azimuth_Emi_ad(3)*sin(m*phi)
      L = 10*(m-1) + 60
      wind10_ad = wind10_ad + ( b_coef(L+5) +b_coef(L+6)*Frequency )*sc_ad  &
           +( TWO*(b_coef(L+7) +b_coef(L+8)*Frequency)*wind_speed  &
           +b_coef(L+9)*seczen +b_coef(L+10)*seczen*Frequency )*sc_ad

      phi_ad = phi_ad -m*Azimuth_Emi_ad(2)*asc(2,m)*sin(m*phi)
      ac_ad = Azimuth_Emi_ad(2)*cos(m*phi)
      L = 10*(m-1) + 30
      wind10_ad = wind10_ad + ( b_coef(L+5) +b_coef(L+6)*Frequency )*ac_ad  &
           +( TWO*(b_coef(L+7) +b_coef(L+8)*Frequency)*wind_speed  &
           +b_coef(L+9)*seczen +b_coef(L+10)*seczen*Frequency )*ac_ad

      phi_ad = phi_ad -m*Azimuth_Emi_ad(1)*asc(1,m)*sin(m*phi)
      ac_ad = Azimuth_Emi_ad(1)*cos(m*phi)
      L = 10*(m-1)
      wind10_ad = wind10_ad + ( b_coef(L+5) +b_coef(L+6)*Frequency )*ac_ad  &
           +( TWO*(b_coef(L+7) +b_coef(L+8)*Frequency)*wind_speed  &
           +b_coef(L+9)*seczen +b_coef(L+10)*seczen*Frequency )*ac_ad

    END DO

    Rel_Azimuth_ad = Rel_Azimuth_ad + phi_ad* DEGREES_TO_RADIANS

  END IF


    Azimuth_Emi = ZERO
    Azimuth_Emi_ad = ZERO
  !  print *, zreflmod_v_ad,zreflmod_h_ad
  ! correction for anisotropic downward radiation
  IF( Transmittance > transmittance_limit_lower .and. Transmittance < transmittance_limit_upper) THEN

      zrough_v_ad = ZERO
      zrough_h_ad = ZERO

      zrough_h_ad = - zreflmod_h_ad * &
            & ( Transmittance**zrough_h *log(Transmittance) ) / &
            & (ONE-Transmittance)

      Transmittance_AD = Transmittance_AD + zreflmod_h_ad/(ONE-Transmittance)* &
        & (-zrough_h * Transmittance**(zrough_h-ONE) + zreflmod_h)

      zrough_v_ad = -zreflmod_v_ad * &
        & ( Transmittance**zrough_v * log(Transmittance) ) / &
        & (ONE-Transmittance)

      Transmittance_AD = Transmittance_AD + zreflmod_v_ad/(ONE-Transmittance)* &
        & (-zrough_v * Transmittance**(zrough_v-ONE) + zreflmod_v)

      zx_ad(:) = ZERO

      DO i = 1,7
         j = i-1
         !Switched h to v Deblonde SSMIS june 7, 2001

         zx_ad(9) = zx_ad(9) + zrough_v_ad * zx(i) * t_c(24+j*3)
         zx_ad(8) = zx_ad(8) + zrough_v_ad * zx(i) * t_c(23+j*3)
         zx_ad(i) = zx_ad(i) + zrough_v_ad *(t_c(22+j*3)+zx(8)*t_c(23+j*3)+zx(9) &
                  *t_c(24+j*3) )

         zx_ad(9) = zx_ad(9) + zrough_h_ad * zx(i) * t_c(3+j*3)
         zx_ad(8) = zx_ad(8) + zrough_h_ad * zx(i) * t_c(2+j*3)
         zx_ad(i) = zx_ad(i) + zrough_h_ad*(t_c(1+j*3) +zx(8)*t_c(2+j*3)   &
                  + zx(9)  *  t_c(3+j*3) )
      END DO
        zrough_v_ad = ZERO
        zrough_h_ad = ZERO

        !Define nine predictors for the effective angle calculation
        zx_ad(8) = zx_ad(8) + zx_ad(9) * 2 * zx(8)
        opdpsfc_ad = zx_ad(8) / opdpsfc
        zx_ad(2) = zx_ad(2) + zx_ad(7) * 2 * zx(2)
        zx_ad(3) = zx_ad(3) + zx_ad(5) * 2 * zx(3)
        zx_ad(2) = zx_ad(2) + zx_ad(3) * zx(4)
        zx_ad(4) = ZERO
        variance_ad = zx_ad(2)
        zx_ad(1) = ZERO

        !Compute surface to space optical depth
        Transmittance_AD = Transmittance_AD - opdpsfc_ad*cos_z /Transmittance

        IF ( variance < varm ) THEN
           varm_ad = variance_ad * ( t_c(44) *Frequency + t_c(45) )
        ELSE
           varm_ad = variance_ad
        END IF

        variance_ad = varm_ad * t_c(43)
        wind10_ad = wind10_ad + variance_ad * 0.00512_fp

  END IF

    zreflmod_v_ad = ZERO
    zreflmod_h_ad = ZERO
    Eh_ad = Emissivity_ad(2)
    Ev_ad = Emissivity_ad(1)
    RhL_ad = (ONE-Foam_Cover)*Eh_ad
    RhS_ad = -(ONE-Foam_Cover)*Eh_ad
    Foam_Cover_ad = Eh_ad*(RhS - RhL-Foam_Rh)

    RvL_ad = (ONE-Foam_Cover)*Ev_ad
    RvS_ad = -(ONE-Foam_Cover)*Ev_ad

    Foam_Cover_ad = Foam_Cover_ad  + Ev_ad*(RvS - RvL-Foam_Rv)
    IF (fastem_version == 4) THEN
      wind10_ad = wind10_ad + 7.75E-06_fp*3.231_fp*Wind_Speed**2.231_fp*Foam_Cover_ad
    ELSE
      wind10_ad = wind10_ad + 2.55_fp*1.95E-05_fp*Wind_Speed**1.55_fp*Foam_Cover_ad
    END IF


    ! Large Scale Correction Calculation
    ! ----------------------------------
    wind10_ad = wind10_ad + zc(10)*RhL_ad &
      + TWO*zc(11)*Wind_Speed*RhL_ad + zc(12)*RhL_ad*seczen

    wind10_ad = wind10_ad + zc(4)*RvL_ad &
      + TWO*zc(5)*Wind_Speed*RvL_ad + zc(6)*RvL_ad*seczen


    small_corr_ad = Rh_Fresnel *RhS_ad
    Rh_Fresnel_ad = RhS_ad * small_corr

    small_corr_ad = small_corr_ad + Rv_Fresnel *RvS_ad
    Rv_Fresnel_ad = RvS_ad * small_corr

    ! Small scale Correction Calculation
    scor_ad = -small_corr_ad*cos_z*cos_z *small_corr

    wind10_ad = wind10_ad + ( Scoef(1)*frequency +Scoef(2)*frequency**2 &
            +TWO*Scoef(3) *windspeed*frequency +TWO*Scoef(4)*windspeed* frequency**2 &
            +TWO*Scoef(5) *windspeed /frequency +TWO*Scoef(6) *windspeed/frequency**2 &
            +Scoef(7) + TWO*Scoef(8) *windspeed  )*scor_ad

    ! Fresnel reflectivity calculation
    ! --------------------------------
    permittivity_ad = ZERO

     ! The adjoint of the horizontal reflectivity
    izRh_AD = TWO*frVar%izRh*Rh_Fresnel_ad
    rzRh_AD = TWO*frVar%rzRh*Rh_Fresnel_ad
    Rh_Fresnel_ad = ZERO
    zRh_AD = CMPLX(rzRh_AD, -izRh_AD, fp)  ! complex conjugate

    ! The adjoint of the vertical reflectivity
    izRv_AD = TWO*frVar%izRv*Rv_Fresnel_ad
    rzRv_AD = TWO*frVar%rzRv*Rv_Fresnel_ad
    Rv_Fresnel_ad  = ZERO
    zRv_AD = CMPLX(rzRv_AD, -izRv_AD, fp)  ! complex conjugate

    ! The adjoint of the complex vertical polarised component
    denom = (frVar%z2+frVar%z1)**2
    z1_AD = -TWO*frVar%z2*zRv_AD / denom
    z2_AD =  TWO*frVar%z1*zRv_AD / denom

    ! The adjoint of the complex horizontal polarised component
    z1_AD = z1_AD - ( TWO*cos_z*zRh_AD / (cos_z+frVar%z1)**2 )

    ! The adjoint of the preserved variables
    permittivity_AD = permittivity_AD + CONJG(cos_z*z2_AD)
    permittivity_AD = permittivity_AD + CONJG(POINT_5*z1_AD/frVar%z1)


    ! Permittivity Calculation
    ! ------------------------
    t_AD = ZERO
    S_AD = ZERO
    t_sq_AD = ZERO
    einf_AD = ZERO

!    perm_imag_AD = -AIMAG(Permittivity_AD)

    perm_imag_AD = -AIMAG(Permittivity_AD)
    perm_Real_AD = REAL(Permittivity_AD,fp)
    Permittivity_AD = ZERO

    f2_AD = (iVar%del2-TWO*iVar%del2*iVar%f2**2/(ONE + iVar%f2**2)) &
          *perm_imag_AD/(ONE + iVar%f2**2)

    del2_AD = iVar%f2*perm_imag_AD/(ONE + iVar%f2**2)

    f1_AD = (iVar%del1-TWO*iVar%del1*iVar%f1**2/(ONE + iVar%f1**2))  &
          *perm_imag_AD/(ONE + iVar%f1**2)
    del1_AD = iVar%f1*perm_imag_AD/(ONE + iVar%f1**2)

    perm_imag_AD = -perm_imag_AD

    f2_AD = f2_AD -TWO*iVar%del2*iVar%f2*perm_Real_AD/(ONE + iVar%f2**2)**2
    del2_AD = del2_AD + perm_Real_AD/(ONE + iVar%f2**2)
    f1_AD = f1_AD -TWO*iVar%del1*iVar%f1*perm_Real_AD/(ONE + iVar%f1**2)**2
    del1_AD = del1_AD + perm_Real_AD/(ONE + iVar%f1**2)
    einf_AD = perm_Real_AD


    einf_AD = einf_AD -del2_AD
    e1_AD = del2_AD
    e1_AD = e1_AD - del1_AD
    es_AD = del1_AD

    tau2_AD = Frequency*f2_AD
    tau1_AD = Frequency*f1_AD
    IF( iVar%S > ZERO ) THEN
      sigma_AD = -perm_imag_AD/(TWO*PI*e0*Frequency)
      ctau2_AD = iVar%tau2_k * tau2_AD
      tau2_AD = tau2_AD * iVar%ctau2
      ctau1_AD = iVar%tau1_k * tau1_AD
      tau1_AD = tau1_AD * iVar%ctau1
      ce1_AD = iVar%e1_k * e1_AD
      e1_AD = e1_AD * iVar%ce1
      ces_AD = iVar%es_k * es_AD
      es_AD = es_AD * iVar%ces
      S_AD = S_AD + (A_COEF(26) + A_COEF(27)*iVar%t + THREE*A_COEF(28)*iVar%S**2)*ctau2_AD
      t_AD = t_AD + iVar%S*A_COEF(27)*ctau2_AD
      t_sq_AD = t_sq_AD + iVar%S*A_COEF(21)*ctau1_AD
      t_AD = t_AD + iVar%S*A_COEF(20)*ctau1_AD
      S_AD = S_AD + ctau1_AD*(A_COEF(19) +A_COEF(20)*iVar%t + A_COEF(21)*iVar%t_sq)
      t_AD = t_AD + A_COEF(14)*iVar%S*ce1_AD
      S_AD = S_AD + (A_COEF(12)+TWO*A_COEF(13)*iVar%S+A_COEF(14)*iVar%t)*ce1_AD
      t_AD = t_AD + A_COEF(8)*iVar%S*ces_AD
      S_AD = S_AD + (A_COEF(6) + TWO*A_COEF(7)*iVar%S + A_COEF(8)*iVar%t)*ces_AD
      beta_AD = -iVar%delta*iVar%sigma*sigma_AD
      delta_AD = -sigma_AD*iVar%beta*iVar%sigma
      sigma25_AD = sigma_AD*exp(-iVar%delta*iVar%beta)
      S_AD = S_AD + sigma25_AD*(A_COEF(35) +A_COEF(36)*iVar%S   &
           +A_COEF(37)*iVar%S**2 +A_COEF(38)*iVar%S**3)
      S_AD = S_AD + iVar%S*(A_COEF(36) +TWO*A_COEF(37)*iVar%S   &
           +THREE*A_COEF(38)*iVar%S**2)*sigma25_AD
      delta_AD = delta_AD + (A_COEF(30)+TWO*A_COEF(31)*iVar%delta+iVar%S*A_COEF(33)  &
               +iVar%S*TWO*A_COEF(34)*iVar%delta)*beta_AD
      S_AD = S_AD + beta_AD*(A_COEF(32) +A_COEF(33)*iVar%delta +A_COEF(34)*iVar%delta**2)
      t_AD = t_AD -delta_AD
    END IF

      t_cu_AD = A_COEF(25)*tau2_AD
      t_sq_AD = t_sq_AD + A_COEF(24)*tau2_AD
      t_AD = t_AD + A_COEF(23)*tau2_AD
      t_cu_AD = t_cu_AD + A_COEF(18)*tau1_AD
      t_sq_AD = t_sq_AD + A_COEF(17)*tau1_AD
      t_AD = t_AD + A_COEF(16)*tau1_AD

      t_sq_AD = t_sq_AD + A_COEF(11)*e1_AD
      t_AD = t_AD + A_COEF(10)*e1_AD
      t_cu_AD = t_cu_AD + A_COEF(5)*es_AD
      t_sq_AD = t_sq_AD + A_COEF(4)*es_AD
      t_AD = t_AD + A_COEF(3)*es_AD
      t_AD = t_AD + A_COEF(1)*einf_AD

      t_AD = t_AD + THREE * iVar%t_sq * t_cu_AD
      t_AD = t_AD + TWO * iVar%t * t_sq_AD
      Temperature_AD = Temperature_AD + t_AD
      Salinity_AD = Salinity_ad + S_AD

    wind_speed_ad = Wind_Speed_ad + wind10_ad

   ! set to zero
   Emissivity_ad = ZERO
   Reflectivity_ad = ZERO
   RETURN

  END SUBROUTINE rttov_fastem5_ad
!
