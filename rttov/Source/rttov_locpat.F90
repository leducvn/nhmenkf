!     Calculate the secant of the zenith angle at the lower boundary
!     of each atmospheric layer.
SUBROUTINE rttov_locpat( &
            & OPTS,       &
            & profiles,   &
            & aux,        &
            & coef,       &
            & angles,     &
            & raytracing)
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
!     1            14/02/2003  Oroginal code. RTIASI.4. Marco Matricardi. ECMWF.
!     2            30/07/2004  RTIASI.5. Marco Matricardi. ECMWF.
!                              Cloud variables removed
!     3            27/02/2009  Profile levels to include ToA. All arrays in
!                              RAYTRACING are on levels except the layer arrays
!                              INT, DMAIR, LTICK. Change RATOESAT and LTICK to
!                              to allow agreement with RTTOV-9.
!     4            27/02/2009  Bug corrected in defining DISPCO2 (P. Rayer)
!     5            02/12/2009  The computation of the local zenith angle has been
!                              changed to take into account the extra top level.
!                              PATHSAT and PATHSUN are now layer arrays (Marco Matricardi).
!     6            05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!     7            15/10/2010  Ensure refraction assumed when addpc true (J Hocking)
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
  TYPE(PROFILE_TYPE   ), INTENT(IN)    :: PROFILES(:)
  TYPE(PROFILE_AUX    ), INTENT(IN)    :: AUX
  TYPE(GEOMETRY_TYPE  ), INTENT(IN)    :: ANGLES(size(profiles))
  TYPE(RTTOV_COEF     ), INTENT(IN)    :: COEF
  TYPE(RAYTRACING_TYPE), INTENT(INOUT) :: RAYTRACING
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: J, IR   , I, K
  INTEGER(KIND=jpim) :: JPLEV , JPLAY
  REAL   (KIND=jprb) :: REARTH
  REAL   (KIND=jprb) :: GRAVL                           ! Gravity at latitude LAT           [m/s^2]
  REAL   (KIND=jprb) :: GRAVH
  REAL   (KIND=jprb) :: ETA   , BETA                    ! Coefficients of the international gravity
! formula
  REAL   (KIND=jprb) :: DISP                            ! The value of the refractive index given
! by the dispersion equation.
  REAL   (KIND=jprb) :: RLH
  REAL   (KIND=jprb) :: LATR
  REAL   (KIND=jprb) :: DFLAT
  REAL   (KIND=jprb) :: FAC
!       Local arrays:
  REAL   (KIND=jprb) :: TEMPP((profiles(1)%NLEVELS) + 1)
  REAL   (KIND=jprb) :: WATER((profiles(1)%NLEVELS) + 1)
  REAL   (KIND=jprb) :: PRES ((profiles(1)%NLEVELS)    )
  INTEGER(KIND=jpim) :: nprofiles                       ! Number of profiles
!-----End of header-------------------------------------------------------------
  nprofiles = size(profiles)
  JPLEV     = profiles(1)%nlevels
  JPLAY     = JPLEV - 1
  DO J = 1, NPROFILES
    RAYTRACING%PATHSUN(:, J) = 1.E+38_jprb
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
    GRAVL = GRAVE * (1._jprb + BETA * (SIN(DEG2RAD * PROFILES(J)%LATITUDE)) ** 2_jpim +      &
      & ETA * (SIN(2._jprb * DEG2RAD * PROFILES(J)%LATITUDE)) ** 2_jpim)
!-------3.  The value of the gravity as a function of altitude H can be--------
!       expressed using the inverse-square law of gravitation:                 |
!                                                                              |
!       GRAVL(H)=GRAVL*(REARTH/(H+REARTH))**2=                                 |
!                GRAVL*(1-2*H/REARTH+3*H**2/REARTH**2+terms of higher order)   |
!                                                                              |
!       If we eliminate the second and higher order terms we can write:        |
!                                                                              |
!       GRAVL(H)=GRAVL-2*GRAVL*H/REARTH=GRAVL-GRAVH*H                          |
!------------------------------------------------------------------------------
    GRAVH = 2.0E-3_jprb * GRAVL / REARTH
!-------4.  Unpack input profile------------------------------------------------
    DO I = 1, JPLEV
! profile arrays - levels, bottom to top
      PRES(I)  = PROFILES(J)%P(profiles(1)%nlevels - I + 1)
      WATER(I) = PROFILES(J)%Q(profiles(1)%nlevels - I + 1)
      TEMPP(I) = PROFILES(J)%T(profiles(1)%nlevels - I + 1)
    ENDDO
!-------5.  Calculate density for dry air (ideal gas)---------------------------
    DO IR = 1, JPLEV
! dry air density profile
! on levels, bottom to top (follows PRES)
      RAYTRACING%DAIR(IR, J) = 1.E+02_jprb * PRES(IR) * MAIR / (1000._jprb * RGC * TEMPP(IR))
    ENDDO
!-------6.  The density of dry air is adjusted to account for the presence ----
!       of water vapour by introducing a scaling factor (i.e. the              |
!       temperature is replaced by the virtual temperature).The density        |
!       for moist air is thus obtained. It is assumed that the partial         |
!       pressure for water vapour can be written as PRES*WATER*1.E-6           |
!------------------------------------------------------------------------------
    DO IR = 1, JPLEV
! reverse direction
      K = JPLEV - IR + 1
! moist air density profile
! on levels, top to bottom (reverses DAIR & WATER)
      RAYTRACING%DMAIR(K, J) =      &
        & RAYTRACING%DAIR(IR, J) * (MAIR * (1.0E6_jprb - WATER(IR)) + MH2O * WATER(IR)) / (1.0E6_jprb * MAIR)
    ENDDO
!-------7.  Set up the height of the surface pressure level---------------------
! height of nearest level below surf
! counted from bottom of profile grid
    RAYTRACING%HGPL(JPLEV - AUX%s(J)%NEARESTLEV_SURF + 1, J) = PROFILES(J)%ELEVATION
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
! integrated layer values for DMAIR,
! lower boundary of first layer is nearest level below surface
! upper boundary of last layer is ToA
      RAYTRACING%INT(IR, J)  =  - 0.5_jprb * (PRES(IR) - PRES(IR - 1)) *      &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 1, J) + 1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 2, J))
      RAYTRACING%INT(IR, J)  = RAYTRACING%INT(IR, J) * 1.0E2_jprb
      RLH = GRAVL / GRAVH
      RAYTRACING%HL(IR, J)   = 1.0E3_jprb * RAYTRACING%HGPL(IR - 1, J)
      RAYTRACING%HGPL(IR, J) = (RLH - SQRT(                                                  &
        & RLH ** 2._jprb - RAYTRACING%HL(IR, J) * (2.0_jprb * RLH - RAYTRACING%HL(IR, J)) -  &
        & 2.0_jprb * RAYTRACING%INT(IR, J) / GRAVH)) * 1.0E-3_jprb
    ENDDO
!-------8.1  Compute height of pressure levels:levels below surface-------------
! integrated layer values for DMAIR,
! upper boundary of first layer is nearest level below surface
! lower boundary of last layer is bottom of profile grid
    DO IR = JPLEV - AUX%s(J)%NEARESTLEV_SURF, 1,  - 1
      RAYTRACING%INT(IR, J)  = 0.5_jprb * (PRES(IR + 1) - PRES(IR - 1 + 1)) *      &
        & (1._jprb / RAYTRACING%DMAIR(JPLEV - IR, J) + 1._jprb / RAYTRACING%DMAIR(JPLEV - IR + 1, J))
      RAYTRACING%INT(IR, J)  = RAYTRACING%INT(IR, J) * 1.0E2_jprb
      RLH = GRAVL / GRAVH
      RAYTRACING%HL(IR, J)   = 1.0E3_jprb * RAYTRACING%HGPL(IR + 1, J)
      RAYTRACING%HGPL(IR, J) = (RLH - SQRT(                                                  &
        & RLH ** 2._jprb - RAYTRACING%HL(IR, J) * (2.0_jprb * RLH - RAYTRACING%HL(IR, J)) -  &
        & 2.0_jprb * RAYTRACING%INT(IR, J) / GRAVH)) * 1.0E-3_jprb
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
! partial pressure of water vapour
! on levels, bottom to top (follows PRES & WATER)
      RAYTRACING%PPW(IR, J) = PRES(IR) * WATER(IR) * 1.0E-6_jprb
!                       *(1._jprb/(1._jprb+WATER(IR)*1.0E-6_jprb*0.378_jprb))
!-------9.2  A correction factor is applied to the dispersion equation to-------
!       account for the presence of carbon dioxide in the air
      IF (opts%CO2_DATA) THEN
! correction factor profile
! on levels, bottom to top (reverses (PROFILES(J)%CO2)
        RAYTRACING%DISPCO2(IR, J) = DISP * (1 + DCO2 * (PROFILES(J)%CO2(jplev - IR + 1) * 1E-6 - 0.0003))
      ELSE
! correction factor profile
! on levels, const
        RAYTRACING%DISPCO2(IR, J) = DISP * (1._jprb + DCO2 * (376._jprb * 1E-6_jprb - 0.0003_jprb))
      ENDIF
!-------9.3  Compute the refractivity index for dry air at temperature TEMPP  --
!       and pressure (PRES-PPW).
! refractive-index-of-dry-air profile
! on levels, bottom to top (follows PRES & TEMPP & PPW)
      RAYTRACING%R(IR, J) = RAYTRACING%DISPCO2(IR, J) * (PRES(IR) - RAYTRACING%PPW(IR, J)) * HTOP / ED1 *         &
        & (1._jprb + 1.0E-8_jprb * (ED2 - ED3 * (TEMPP(IR) - T0)) * (PRES(IR) - RAYTRACING%PPW(IR, J)) * HTOP) /  &
        & (1._jprb + ED4 * (TEMPP(IR) - T0))
!-------9.4  Compute the refractivity index for air containig a partial --------
!       pressure PPW of water vapour.
! refractive-index-of-moist-air profile
! on levels, bottom to top (follows PRES & TEMPP & PPW)
      RAYTRACING%R(IR, J) =      &
        & RAYTRACING%R(IR, J) - RAYTRACING%PPW(IR, J) * HTOP * (EW1 - EW2 * (waver * CTOM) ** 2._jprb) * 1.0E-10_jprb
      RAYTRACING%R(IR, J) = RAYTRACING%R(IR, J) + 1._jprb
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
! intermediate trig variable
! on levels, top to bottom (reverses HGPL & R)
!---------------RT10--------------------------------------------------------------------
!             RAYTRACING%RATOESAT(IR,J) = (REARTH+COEF%FC_SAT_HEIGHT)                / &
!                             (REARTH+(RAYTRACING%HGPL(JPLEV-IR,J)))                 * &
!                             (RAYTRACING%R(JPLEV,J)/RAYTRACING%R(JPLEV-IR,J))
!---------------------------------------------------------------------------------------
! to agree with RT9---------------------------------------------------------------------
! intermediate trig variable
! on levels, top to bottom (reverses HGPL & R)
        RAYTRACING%RATOESAT(IR, J) = (REARTH + COEF%FC_SAT_HEIGHT) / (REARTH + (RAYTRACING%HGPL(JPLAY + 1 - IR, J))) *      &
          & (RAYTRACING%R(JPLAY, J) / RAYTRACING%R(JPLAY + 1 - IR, J))
!---------------------------------------------------------------------------------------
! intermediate trig variable
! on levels, top to bottom (uses RATOESAT)
        RAYTRACING%ZASAT(IR, J)    = SIN(ANGLES(J)%VIEWANG * DEG2RAD) * RAYTRACING%RATOESAT(IR, J)
! secant of satellite zenith angle
! on levels, top to bottom (uses ZASAT)
        RAYTRACING%PATHSAT(IR, J)  = 1._jprb / SQRT(1._jprb - RAYTRACING%ZASAT(IR, J) * RAYTRACING%ZASAT(IR, J))
!print*,acos(1/RAYTRACING%PATHSAT(IR,J))/deg2rad,IR,RAYTRACING%R(IR,J)
      ENDDO
    ELSE
      DO IR = 1, JPLAY
!----------RT10-------------------------------------------------------------------------
!            RAYTRACING%RATOESAT(IR,J) = (REARTH+COEF%FC_SAT_HEIGHT)                  / &
!                             (REARTH+(RAYTRACING%HGPL(JPLEV-IR,J)))
!---------------------------------------------------------------------------------------
! to agree with RT9---------------------------------------------------------------------
        RAYTRACING%RATOESAT(IR, J) = (REARTH + COEF%FC_SAT_HEIGHT) / (REARTH + (RAYTRACING%HGPL(JPLAY + 1 - IR, J)))
!---------------------------------------------------------------------------------------
        RAYTRACING%ZASAT(IR, J)    = SIN(ANGLES(J)%VIEWANG * DEG2RAD) * RAYTRACING%RATOESAT(IR, J)
! secant of satellite zenith angle - levels, top to bottom (reverses HGPL & R)
        RAYTRACING%PATHSAT(IR, J)  = 1._jprb / SQRT(1._jprb - RAYTRACING%ZASAT(IR, J) * RAYTRACING%ZASAT(IR, J))
!print*,acos(1/RAYTRACING%PATHSAT(IR,J))/deg2rad,IR
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
!--------RT10---------------------------------------------------------------------------
!                RAYTRACING%RATOESUN(IR,J)     = (REARTH/(REARTH                      + &
!                                      RAYTRACING%HGPL(JPLEV-IR,J)))                  * &
!                                (RAYTRACING%R(1,J)/RAYTRACING%R(JPLEV-IR,J))
!---------------------------------------------------------------------------------------
! to agree with RT9---------------------------------------------------------------------
            RAYTRACING%RATOESUN(IR, J) = (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J))) *      &
              & (RAYTRACING%R(1, J) / RAYTRACING%R(JPLAY + 1 - IR, J))
!---------------------------------------------------------------------------------------
            RAYTRACING%ZASUN(IR, J)    = SIN(PROFILES(J)%SUNZENANGLE * DEG2RAD) * RAYTRACING%RATOESUN(IR, J)
! secant of sun zenith angle
! on levels, top to bottom (reverses HGPL & R)
! excluding ToA
            RAYTRACING%PATHSUN(IR, J)  = 1._jprb / SQRT(1. - RAYTRACING%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J))
!            print*,acos(1/RAYTRACING%PATHSUN(IR,J))/deg2rad,IR,'pathsun'
          ELSE
!--------------RT10---------------------------------------------------------------------
!                RAYTRACING%RATOESUN(IR,J)     = REARTH                      / &
!                                         (REARTH+RAYTRACING%HGPL(JPLEV-IR,J))
!---------------------------------------------------------------------------------------
! to agree with RT9---------------------------------------------------------------------
            RAYTRACING%RATOESUN(IR, J) = REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J))
!---------------------------------------------------------------------------------------
            RAYTRACING%ZASUN(IR, J)    = SIN(PROFILES(J)%SUNZENANGLE * DEG2RAD) * RAYTRACING%RATOESUN(IR, J)
! secant of sun zenith angle
! ToA
            RAYTRACING%PATHSUN(IR, J)  = 1._jprb / SQRT(1. - RAYTRACING%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J))
!            print*,acos(1/RAYTRACING%PATHSUN(IR,J))/deg2rad,IR,'jplay'
          ENDIF
        ENDDO
      ELSE
        DO IR = 1, JPLAY
!------------RT10-----------------------------------------------------------------------
!                RAYTRACING%RATOESUN(IR,J)     = (REARTH                              / &
!                                        (REARTH+RAYTRACING%HGPL(JPLEV-IR,J)))
!---------------------------------------------------------------------------------------
! to agree with RT9---------------------------------------------------------------------
          RAYTRACING%RATOESUN(IR, J) = (REARTH / (REARTH + RAYTRACING%HGPL(JPLAY + 1 - IR, J)))
!---------------------------------------------------------------------------------------
          RAYTRACING%ZASUN(IR, J)    = SIN(PROFILES(J)%SUNZENANGLE * DEG2RAD) * RAYTRACING%RATOESUN(IR, J)
! secant of sun zenith angle
! on levels, top to bottom (reverses HGPL & R)
          RAYTRACING%PATHSUN(IR, J)  = 1._jprb / SQRT(1._jprb - RAYTRACING%ZASUN(IR, J) * RAYTRACING%ZASUN(IR, J))
!            print*,acos(1/RAYTRACING%PATHSUN(IR,J))/deg2rad,IR
        ENDDO
      ENDIF
    ENDIF
!-------------------------------------------------------------------------------
! to agree with RT9
    RAYTRACING%HGPL(JPLEV, J) = RAYTRACING%HGPL(JPLEV - 1, J)
    DO IR = 1, JPLAY
      RAYTRACING%LTICK(IR, J) = RAYTRACING%HGPL(IR + 1, J) - RAYTRACING%HGPL(IR, J)
!          print*,RAYTRACING%LTICK(IR,J),IR,'thick'
    ENDDO
  ENDDO
END SUBROUTINE RTTOV_LOCPAT
