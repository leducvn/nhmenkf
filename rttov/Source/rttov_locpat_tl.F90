!     Calculate the secant of the zenith angle at the lower boundary
!     of each atmospheric layer.
SUBROUTINE RTTOV_LOCPAT_TL( &
            & OPTS,          &
            & PROFILES,      &
            & PROFILES_TL,   &
            & AUX,           &
            & COEF,          &
            & ANGLES,        &
            & RAYTRACING,    &
            & RAYTRACING_TL)
!     Description:
!     The calculation of the secant of the zenith angle at the lower
!     boundary of each atmospheric layer is performed by evaluating
!     the height of the pressure levels.
!
!     Method.
!     The calculation of the height of the pressure levels is based
!     on the hydrostatic equation. To account for the presence of
!     water vapour, virtual temperature are used. The variation of
!     gravity with latitude is introduced using the international
!     gravity formula. The bending of rays as they traverse the atmosphere
!     is fully accounted for applying the Snell's law. For the computation
!     of the refractive index of air, an updated version of Edlen's formula
!     is used.
!
!     K.P. Birch and M.J.Downs:'An Updated Edlen Equation for the
!     raefractive index of air'. Metrologia, 1993, 30, 155-162.
!
!     K.P. Birch and M.J.Downs:'Correction to the Updated Edlen Equation
!     for the refractive index of air'. Metrologia, 1994, 31, 315-316.
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1            14/02/2003  Original code. RTIASI.4. Marco Matricardi. ECMWF.
!     2            30/07/2004  RTIASI.5. Marco Matricardi. ECMWF.
!                              Cloud variables removed
!     3            30/11/2005  Changed intent of raytracing to inout (R Saunders)
!     4            04/02/2008  lgradp option for TL/AD of pressure levels, Niels Bormann
!     5            15/07/2009  Bug corrected in defining DISPCO2 (P. Rayer)
!     6            15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!     7            02/12/2009  The computation of the local zenith angle has been
!                              changed to take into account the extra top level.
!                              PATHSAT and PATHSUN are now layer arrays (Marco Matricardi).
!     8            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!     9            02/09/2010  Bug fixes (J Hocking)
!    10            15/10/2010  Ensure refraction assumed when addpc true (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!     Module used:
  USE rttov_types, ONLY :  &
       & rttov_options,   &
       & rttov_coef,      &
       & profile_aux,     &
       & profile_Type,    &
       & geometry_Type,   &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY :  &
       & deg2rad, &
       & D1,      &
       & D2,      &
       & D3,      &
       & D4,      &
       & D5,      &
       & DCO2,    &
       & ED1,     &
       & ED2,     &
       & ED3,     &
       & ED4,     &
       & EW1,     &
       & EW2,     &
       & HTOP,    &
       & CTOM,    &
       & WAVER,   &
       & RGC,     &
       & MAIR,    &
       & MH2O,    &
       & FLATT,   &
       & EQRAD,   &
       & OMEGA,   &
       & GRAVE,   &
       & T0,      &
       & MAX_SOL_ZEN
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!       Scalar arguments with intent in:
  TYPE(RTTOV_OPTIONS  ), INTENT(IN)    :: OPTS
!       Arrays arguments with intent in:
  TYPE(PROFILE_TYPE   ), INTENT(IN)    :: PROFILES   (:)             ! Atmospheric profiles
  TYPE(PROFILE_TYPE   ), INTENT(INOUT) :: PROFILES_TL(size(profiles))! Atmospheric profiles
  TYPE(PROFILE_AUX    ), INTENT(IN)    :: AUX
  TYPE(GEOMETRY_TYPE  ), INTENT(IN)    :: ANGLES(size(profiles))     ! angles
  TYPE(RTTOV_COEF     ), INTENT(IN)    :: COEF
  TYPE(RAYTRACING_TYPE), INTENT(IN)    :: RAYTRACING
  TYPE(RAYTRACING_TYPE), INTENT(INOUT) :: RAYTRACING_TL
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: J, IR   , I, K
  INTEGER(KIND=jpim) :: JPLEV , JPLAY
  REAL   (KIND=jprb) :: REARTH
  REAL   (KIND=jprb) :: GRAVL                             ! Gravity at latitude LAT              [m/s^2]
  REAL   (KIND=jprb) :: GRAVH
  REAL   (KIND=jprb) :: ETA   , BETA                      ! Coefficients of the international gravity
! formula
  REAL   (KIND=jprb) :: DISP                              ! The value of the refractive index given by
! the dispersion equation.
  REAL   (KIND=jprb) :: RLH
  REAL   (KIND=jprb) :: LATR
  REAL   (KIND=jprb) :: DFLAT
  REAL   (KIND=jprb) :: FAC
!       Local arrays:
  REAL   (KIND=jprb) :: PRES   ((profiles(1)%NLEVELS)    )! Profile level pressures      [hPa]
  REAL   (KIND=jprb) :: TEMPP  ((profiles(1)%NLEVELS) + 1)! Temperature profile          [K]
  REAL   (KIND=jprb) :: WATER  ((profiles(1)%NLEVELS) + 1)! Profile volume mixing ratio  [ppmv]
  REAL   (KIND=jprb) :: PRES_D ((profiles(1)%NLEVELS)    )! Profile level pressures      [hPa]
  REAL   (KIND=jprb) :: TEMPP_D((profiles(1)%NLEVELS) + 1)! Temperature profile          [K]
  REAL   (KIND=jprb) :: WATER_D((profiles(1)%NLEVELS) + 1)! Profile volume mixing ratio  [ppmv]
  INTEGER(KIND=jpim) :: nprofiles                         ! Number of profiles
!-----End of header-------------------------------------------------------------
  nprofiles = size(profiles)
  JPLEV     = profiles(1)%nlevels
  JPLAY     = JPLEV - 1
  DO J = 1, NPROFILES
    RAYTRACING_TL%PATHSUN(:, J) = 0_jprb
!-------1.  Caculate the earth's radius at latitude LAT assuming the-----------
!       Earth is an ellipsoid of revolution                                    |
!------------------------------------------------------------------------------
    DFLAT = (1.0_jprb - FLATT) ** 2_jpim
    LATR = DEG2RAD * ABS(PROFILES(J)%LATITUDE)
    REARTH = SQRT(EQRAD ** 2_jpim * DFLAT / (SIN(LATR) ** 2_jpim + DFLAT * COS(LATR) ** 2_jpim))
!-------2.  The valute of earth's gravity at surface at latitude LAT is--------
!       computed using the international gravity formula.                      |
!------------------------------------------------------------------------------
    FAC = (omega ** 2._jprb * (EQRAD * 1000._jprb)) / (GRAVE)
    BETA = 5._jprb * FAC / 2._jprb - FLATT - 17._jprb * FAC * FLATT / 14._jprb
    ETA = FLATT * (5._jprb * FAC - FLATT) / 8._jprb
    GRAVL = GRAVE * (1 + BETA * (SIN(DEG2RAD * PROFILES(J)%LATITUDE)) ** 2_jpim +      &
      & ETA * (SIN(2 * DEG2RAD * PROFILES(J)%LATITUDE)) ** 2_jpim)
!-------3.  The value of the gravity as a function of altitude H can be--------
!       expressed using the inverse-square law of gravitation:                 |
!                                                                              |
!       GRAVL(H)=GRAVL*(REARTH/(H+REARTH))**2._jprb=                                 |
!                GRAVL*(1-2*H/REARTH+3*H**2._jprb/REARTH**2._jprb+terms of higher order)   |
!                                                                              |
!       If we eliminate the second and higher order terms we can write:        |
!                                                                              |
!       GRAVL(H)=GRAVL-2*GRAVL*H/REARTH=GRAVL-GRAVH*H                          |
!------------------------------------------------------------------------------
    GRAVH = 2.0E-3_jprb * GRAVL / REARTH
!-------4.  Unpack input profile------------------------------------------------
    DO I = 1, JPLEV
      PRES_D(I)  = PROFILES(J)%P(profiles(1)%NLEVELS - I + 1)
      WATER_D(I) = PROFILES(J)%Q(profiles(1)%NLEVELS - I + 1)
      TEMPP_D(I) = PROFILES(J)%T(profiles(1)%NLEVELS - I + 1)
    ENDDO
    DO I = 1, JPLEV
      WATER(I) = PROFILES_TL(J)%Q(profiles(1)%NLEVELS - I + 1)
      TEMPP(I) = PROFILES_TL(J)%T(profiles(1)%NLEVELS - I + 1)
    ENDDO
! to agree with RT9
    WATER(jplev) = PROFILES_TL(J)%Q(2)
    TEMPP(jplev) = PROFILES_TL(J)%T(2)
    IF (OPTS%LGRADP) THEN
      DO I = 1, JPLEV
        PRES(I) = PROFILES_TL(J)%P(profiles(1)%NLEVELS - I + 1)
      ENDDO
    ELSE
      PRES = 0.0_JPRB
    ENDIF
!-------5.  Calculate density for dry air (ideal gas)---------------------------
    DO IR = 1, JPLEV
      RAYTRACING_TL%DAIR(IR, J) =                                                                           &
        & ( - 1.E+02_jprb * PRES_D(IR) * MAIR / (1000._jprb * RGC * TEMPP_D(IR) ** 2._jprb)) * TEMPP(IR) +  &
        & 1.E+02_jprb * PRES(IR) * MAIR / (1000._jprb * RGC * TEMPP_D(IR))
    ENDDO
!-------6.  The density of dry air is adjusted to account for the presence ----
!       of water vapour by introducing a scaling factor (i.e. the              |
!       temperature is replaced by the virtual temperature).The density        |
!       for moist air is thus obtained. It is assumed that the partial         |
!       pressure for water vapour can be written as PRES*WATER*1.E-6           |
!------------------------------------------------------------------------------
    DO IR = 1, JPLEV
      K = JPLEV - IR + 1
      RAYTRACING_TL%DMAIR(K, J) =                                                                                       &
        & RAYTRACING_TL%DAIR(IR, J) * (MAIR * (1.0E6_jprb - WATER_D(IR)) + MH2O * WATER_D(IR)) / (1.0E6_jprb * MAIR) +  &
        & RAYTRACING%DAIR(IR, J) * (MAIR * ( - WATER(IR)) + MH2O * WATER(IR)) / (1.0E6_jprb * MAIR)
    ENDDO
!-------7.  Set up the height of the surface pressure level---------------------
    RAYTRACING_TL%HGPL(JPLEV - AUX%s(J)%NEARESTLEV_SURF + 1, J) = 0._jprb
!-------8  Compute height of pressure levels:levels above surface. ------------
!       The height of pressure levels H is obtained by integrating             |
!       the hydrostatic equation dPRES=-GRAVL(H)*DMAIR*dH between              |
!       two adjacent pressure levels                                           |
!                                                                              |
!           -P2         - H2                                                   |
!          |           |                                                       |
!          | dP/DMAIR= | GRAVL(H)*dH                                           |
!          |           |                                                       |
!         -  P1       -   H1                                                   |
!                                                                              |
!       The integration id DMAIR is carried out assuming this                  |
!       quantity varies linearly with pressure. The integral is                |
!       then computed applying the trapezoidal formula.                        |
!------------------------------------------------------------------------------
    DO IR = JPLEV - AUX%s(J)%NEARESTLEV_SURF + 2, JPLEV
      RAYTRACING_TL%INT(IR, J)  = 0.5_jprb * (PRES_D(IR) - PRES_D(IR - 1)) * (                                    &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 1, J) ** 2._jprb) * RAYTRACING_TL%DMAIR(JPLEV - IR + 1, J) +   &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 2, J) ** 2._jprb) * RAYTRACING_TL%DMAIR(JPLEV - IR + 2, J)) -  &
        & 0.5_jprb * (PRES(IR) - PRES(IR - 1)) *                                                                  &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 1, J) + 1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 2, J))
      RAYTRACING_TL%INT(IR, J)  = RAYTRACING_TL%INT(IR, J) * 1.0E2_jprb! Integral in units [m2.s-1]
      RLH = GRAVL / GRAVH
      RAYTRACING_TL%HL(IR, J)   = 1.0E3_jprb * RAYTRACING_TL%HGPL(IR - 1, J)
      RAYTRACING_TL%HGPL(IR, J) =                                                                                           &
        & 1.0E-3_jprb * (RAYTRACING_TL%HL(IR, J) * (RLH - RAYTRACING%HL(IR, J)) + RAYTRACING_TL%INT(IR, J) / GRAVH) / SQRT( &
        & RLH ** 2._jprb - RAYTRACING%HL(IR, J) * (2.0_jprb * RLH - RAYTRACING%HL(IR, J)) -                                 &
        & 2.0_jprb * RAYTRACING%INT(IR, J) / GRAVH)
    ENDDO
!-------8.1  Compute height of pressure levels:levels below surface-------------
    DO IR = JPLEV - AUX%s(J)%NEARESTLEV_SURF, 1,  - 1
      RAYTRACING_TL%INT(IR, J)  =  - 0.5_jprb * (PRES_D(IR + 1) - PRES_D(IR - 1 + 1)) * (                         &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR, J) ** 2._jprb) * RAYTRACING_TL%DMAIR(JPLEV - IR, J) +           &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 1, J) ** 2._jprb) * RAYTRACING_TL%DMAIR(JPLEV - IR + 1, J)) +  &
        & 0.5_jprb * (PRES(IR + 1) - PRES(IR - 1 + 1)) *                                                          &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR, J) + 1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 1, J))
      RAYTRACING_TL%INT(IR, J)  = RAYTRACING_TL%INT(IR, J) * 1.0E2_jprb! Integral in units [m2.s-1]
      RLH = GRAVL / GRAVH
      RAYTRACING_TL%HL(IR, J)   = 1.0E3_jprb * RAYTRACING_TL%HGPL(IR + 1, J)
      RAYTRACING_TL%HGPL(IR, J) = 1.0E-3_jprb * (                                                                            &
        & RLH * RAYTRACING_TL%HL(IR, J) - RAYTRACING%HL(IR, J) * RAYTRACING_TL%HL(IR, J) + RAYTRACING_TL%INT(IR, J) / GRAVH) &
        &  / SQRT(RLH ** 2._jprb - RAYTRACING%HL(IR, J) * (2.0_jprb * RLH - RAYTRACING%HL(IR, J)) -                          &
        & 2.0_jprb * RAYTRACING%INT(IR, J) / GRAVH)
    ENDDO
!-------9.  Compute refractive index of air using a updated  version of Edlen -
!       equation.                                                              |
!------------------------------------------------------------------------------
!-------9.1  Compute the refractivity index as a function of wave-number -------
!       using the dispersion equation
    DISP = 1.0E-8_jprb * (                                                                                                 &
      & D1 + D2 * (D3 - (waver * CTOM) ** 2._jprb) ** ( - 1._jprb) + D4 * (D5 - (waver * CTOM) ** 2._jprb) ** ( - 1._jprb) &
      & )
    DO IR = 1, JPLEV
      RAYTRACING_TL%PPW(IR, J) = PRES_D(IR) * WATER(IR) * 1.0E-6_jprb + PRES(IR) * WATER_D(IR) * 1.0E-6_jprb
!-------9.2  A correction factor is applied to the dispersion equation to-------
!       account for the presence of carbon dioxide in the air
      IF (opts%CO2_DATA) THEN
        RAYTRACING_TL%DISPCO2(IR, J) = DISP * (DCO2 * PROFILES_TL(J)%CO2(jplev - IR + 1) * 1E-6)
      ELSE
        RAYTRACING_TL%DISPCO2(IR, J) = 0._jprb
      ENDIF
!-------9.3  Compute the refractivity index for dry air at temperature TEMPP  --
!       and pressure (PRES-PPW).
      RAYTRACING_TL%R(IR, J) = RAYTRACING_TL%DISPCO2(IR, J) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP / ED1 *            &
        & (1._jprb + 1.0E-8_jprb * (ED2 - ED3 * (TEMPP_D(IR) - T0)) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP) /         &
        & (1._jprb + ED4 * (TEMPP_D(IR) - T0)) +                                                                             &
        & (PRES(IR) - RAYTRACING_TL%PPW(IR, J)) * RAYTRACING%DISPCO2(IR, J) * HTOP / ED1 *                                   &
        & (1._jprb + 1.0E-8_jprb * (ED2 - ED3 * (TEMPP_D(IR) - T0)) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP) /         &
        & (1._jprb + ED4 * (TEMPP_D(IR) - T0)) -                                                                             &
        & TEMPP(IR) * ED4 * RAYTRACING%DISPCO2(IR, J) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP / ED1 *                  &
        & (1._jprb + 1.0E-8_jprb * (ED2 - ED3 * (TEMPP_D(IR) - T0)) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP) /         &
        & (1._jprb + ED4 * (TEMPP_D(IR) - T0)) ** 2._jprb +                                                                  &
        & (PRES(IR) - RAYTRACING_TL%PPW(IR, J)) * HTOP * RAYTRACING%DISPCO2(IR, J) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) *  &
        & HTOP / ED1 * (1.0E-8_jprb * (ED2 - ED3 * (TEMPP_D(IR) - T0))) / (1._jprb + ED4 * (TEMPP_D(IR) - T0)) -   &
        & TEMPP(IR) * 1.0E-8_jprb * ED3 * RAYTRACING%DISPCO2(IR, J) * (PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP / ED1 *    &
        & ((PRES_D(IR) - RAYTRACING%PPW(IR, J)) * HTOP) / (1._jprb + ED4 * (TEMPP_D(IR) - T0))
!-------9.4  Compute the refractivity index for air containig a partial --------
!       pressure PPW of water vapour.
      RAYTRACING_TL%R(IR, J) =      &
        & RAYTRACING_TL%R(IR, J) - RAYTRACING_TL%PPW(IR, J) * HTOP * (EW1 - EW2 * (waver * CTOM) ** 2._jprb) * 1.0E-10_jprb
    ENDDO
!------------------------------------------------------------------------------
!-------10.  Compute secant of the zenith angle at the lower boundary of each -
!       layer for the satellite-to-surface line of view.                       |
!                                                                              |
!       The atmospheric layers are considered as concentric rings. If we trace |
!       a ray across these rings at any angle other than nadir, the local      |
!       angle relative to the outward radial direction at the point of         |
!       intersection will be different at each ring because due to the         |
!       curvature of the Earth and to atmospheric refraction. The secant of    |
!       the local zenith angle PATHSAT is thus computed taking into account    |
!       the geometry of the situation and the bending of rays as they traverse |
!       the atmosphere (by application of Snell's law).                        |
!------------------------------------------------------------------------------
    IF (OPTS%ADDREFRAC .OR. OPTS%ADDPC) THEN
      DO IR = 1, JPLAY
!            RAYTRACING_TL%RATOESAT(IR,J)     =(REARTH+COEF%FC_SAT_HEIGHT)   / &
!                   (REARTH+RAYTRACING%HGPL(JPLEV-IR,J))                   * &
!                   (RAYTRACING_TL%R(JPLEV,J)/RAYTRACING%R(JPLEV-IR,J))   - &
!                   (REARTH+COEF%FC_SAT_HEIGHT)                               / &
!                   (REARTH+RAYTRACING%HGPL(JPLEV-IR,J))                   * &
!                   (RAYTRACING%R(JPLEV,J)*RAYTRACING_TL%R(JPLEV-IR,J)    / &
!                   RAYTRACING%R(JPLEV-IR,J)**2._jprb)                           - &
!                   (REARTH+COEF%FC_SAT_HEIGHT)                               / &
!                   (REARTH+RAYTRACING%HGPL(JPLEV-IR,J))**2._jprb                * &
!                   (RAYTRACING%R(JPLEV,J)/RAYTRACING%R(JPLEV-IR,J))      * &
!                   RAYTRACING_TL%HGPL(JPLEV-IR,J)
! to agree with RT9
! NB here RAYTRACING%R(JPLEV-1,J) agrees with RTTOV-9 value RAYTRACING%R(JPLEV,J)
        RAYTRACING_TL%RATOESAT(IR, J) = (REARTH + COEF%FC_SAT_HEIGHT) / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) *      &
          & (RAYTRACING_TL%R(JPLAY, J) / RAYTRACING%R(JPLAY + 1 - IR, J)) -                                                  &
          & (REARTH + COEF%FC_SAT_HEIGHT) / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) *                                  &
          & (RAYTRACING%R(JPLAY, J) * RAYTRACING_TL%R(JPLAY + 1 - IR, J) / RAYTRACING%R(JPLAY + 1 - IR, J) ** 2._jprb) -     &
          & (REARTH + COEF%FC_SAT_HEIGHT) / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) ** 2._jprb *                       &
          & (RAYTRACING%R(JPLAY, J) / RAYTRACING%R(JPLAY + 1 - IR, J)) * RAYTRACING_TL%HGPL(JPLAY + 1 - IR, J)
        RAYTRACING_TL%ZASAT(IR, J)    = SIN(ANGLES(J)%VIEWANG * DEG2RAD) * RAYTRACING_TL%RATOESAT(IR, J)
        RAYTRACING_TL%PATHSAT(IR, J)  = RAYTRACING_TL%ZASAT(IR, J) * RAYTRACING%ZASAT(IR, J) /      &
          & (1._jprb - RAYTRACING%ZASAT(IR, J) * RAYTRACING%ZASAT(IR, J)) ** 1.5_jprb
      ENDDO
    ELSE
      DO IR = 1, JPLAY
        RAYTRACING_TL%RATOESAT(IR, J) =  -                                                              &
          & (REARTH + COEF%FC_SAT_HEIGHT) / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) ** 2._jprb *  &
          & RAYTRACING_TL%HGPL(JPLAY + 1 - IR, J)
        RAYTRACING_TL%ZASAT(IR, J)    = SIN(ANGLES(J)%VIEWANG * DEG2RAD) * RAYTRACING_TL%RATOESAT(IR, J)
        RAYTRACING_TL%PATHSAT(IR, J)  = RAYTRACING_TL%ZASAT(IR, J) * RAYTRACING%ZASAT(IR, J) /      &
          & (1._jprb - RAYTRACING%ZASAT(IR, J) * RAYTRACING%ZASAT(IR, J)) ** 1.5_jprb
      ENDDO
    ENDIF
!-------10.1  Compute secant of the zenith angle at the lower boundary of each-
!       layer for the sun-to-surface line of view.                             |
!------------------------------------------------------------------------------
    IF (OPTS%ADDSOLAR .AND. PROFILES(J)%SUNZENANGLE >= 0.0 .AND. &
                            PROFILES(J)%SUNZENANGLE < MAX_SOL_ZEN) THEN
      IF (OPTS%ADDREFRAC .OR. OPTS%ADDPC) THEN
        DO IR = 1, JPLAY
          IF (IR .NE. JPLAY) THEN
            RAYTRACING_TL%RATOESUN(IR, J) = (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J))) *                    &
              & (RAYTRACING_TL%R(1, J) / RAYTRACING%R(JPLAY + 1 - IR, J)) -                                               &
              & (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J))) *                                                &
              & (RAYTRACING%R(1, J) * RAYTRACING_TL%R(JPLAY + 1 - IR, J) / RAYTRACING%R(JPLAY + 1 - IR, J) ** 2._jprb) -  &
              & (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) ** 2._jprb) *                                     &
              & (RAYTRACING%R(1, J) / RAYTRACING%R(JPLAY + 1 - IR, J)) * RAYTRACING_TL%HGPL(JPLAY + 1 - IR, J)
            RAYTRACING_TL%ZASUN(IR, J)    = SIN(PROFILES(J)%SUNZENANGLE * DEG2RAD) * RAYTRACING_TL%RATOESUN(IR, J)
            RAYTRACING_TL%PATHSUN(IR, J)  = RAYTRACING_TL%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J) /      &
              & (1._jprb - RAYTRACING%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J)) ** 1.5_jprb
          ELSE
            RAYTRACING_TL%RATOESUN(IR, J) =  -      &
              & (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) ** 2._jprb) * RAYTRACING_TL%HGPL(JPLAY + 1 - IR, J)
            RAYTRACING_TL%ZASUN(IR, J)    = SIN(PROFILES(J)%SUNZENANGLE * DEG2RAD) * RAYTRACING_TL%RATOESUN(IR, J)
            RAYTRACING_TL%PATHSUN(IR, J)  = RAYTRACING_TL%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J) /      &
              & (1._jprb - RAYTRACING%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J)) ** 1.5_jprb
          ENDIF
        ENDDO
      ELSE
        DO IR = 1, JPLAY
          RAYTRACING_TL%RATOESUN(IR, J) =      &
            &  - (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)) ** 2._jprb) * RAYTRACING_TL%HGPL(JPLAY + 1 - IR, J)
          RAYTRACING_TL%ZASUN(IR, J)    = SIN(PROFILES(J)%SUNZENANGLE * DEG2RAD) * RAYTRACING_TL%RATOESUN(IR, J)
          RAYTRACING_TL%PATHSUN(IR, J)  = RAYTRACING_TL%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J) /      &
            & (1._jprb - RAYTRACING%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J)) ** 1.5_jprb
        ENDDO
      ENDIF
    ENDIF
!-------------------------------------------------------------------------------
! to agree with RT9
    RAYTRACING_TL%HGPL(JPLEV, J) = RAYTRACING_TL%HGPL(JPLEV - 1, J)
    DO IR = 1, JPLAY
      RAYTRACING_TL%LTICK(IR, J) = RAYTRACING_TL%HGPL(IR + 1, J) - RAYTRACING_TL%HGPL(IR, J)
    ENDDO
  ENDDO
END SUBROUTINE RTTOV_LOCPAT_TL
