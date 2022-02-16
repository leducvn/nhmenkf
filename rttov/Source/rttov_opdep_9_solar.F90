
SUBROUTINE rttov_opdep_9_solar( &
            & nlayers,     &
            & chanprof,    &
            & profiles,    &
            & sun,         &
            & predictors,  &
            & coef,        &
            & opdp_path,   &
            & opdpsun_ref)
!
! Description:
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEP and RTTAU from previous versions of RTTOV
! Only one profile per call
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
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    29/01/2007  Modified to remove polarisation (R Saunders)
!  1.2    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - size, index
!                     labels, looping (P. Rayer)
!  1.3    30/06/2009  Use rttov9 constant intervals from rttov_const
!                     and change test on upper bounds from <= to <
!                     Philippe Marguinaud
!  1.4    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.5    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.6    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  1.7    14/12/2010  Use traj0_sta%sun array to flag channels for which solar calculations
!                     should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_Type
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & rttov9_wv2000_00, &
       & rttov9_wv2295_25, &
       & rttov9_wv2360_00, &
       & max_sol_zen
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers                             ! Number of pressure levels
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)                         ! Channel indices
  TYPE(profile_Type   ), INTENT(IN)    :: profiles(:)                         ! Atmospheric profiles
  LOGICAL(KIND=jplm)   , INTENT(IN)    :: sun(size(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors                          ! Predictors
  TYPE(rttov_coef     ), INTENT(IN)    :: coef                                ! Coefficients
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path
  REAL(KIND=jprb)      , INTENT(OUT)   :: opdpsun_ref(nlayers, size(chanprof))! layer optical depth
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: opticaldepth(nlayers, size(chanprof))    ! raw layer optical depth
  REAL   (KIND=jprb) :: chanx
  INTEGER(KIND=jpim) :: lev         , lay, chan, j, prof, nlevels! loop variables
! cloud liquid water local variables
  INTEGER(KIND=jpim) :: ii
  INTEGER(KIND=jpim) :: nchannels                                ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_SOLAR', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
!-----------------------------------------
!1. calculate layer gaseous optical depths
!-----------------------------------------
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    prof  = chanprof(j)%prof
    IF (sun(j)) THEN
      DO lay = 1, nlayers
        opticaldepth(lay, j) = coef%mixedgas(lay, chan, 1) * predictors%mixedgas_sun(1, lay, prof)
        DO ii = 2, 8
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%mixedgas(lay, chan, ii) * predictors%mixedgas_sun(ii, lay, prof)
        ENDDO
        opticaldepth(lay, j) =      &
          & opticaldepth(lay, j) + coef%mixedgas(lay, chan, 9) * predictors%mixedgas_sun(9, lay, prof)
        opticaldepth(lay, j) =      &
          & opticaldepth(lay, j) + coef%mixedgas(lay, chan, 10) * predictors%mixedgas_sun(10, lay, prof)
      ENDDO
    ENDIF
  ENDDO
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    chan  = chanprof(j)%chan
    chanx = coef%ff_cwn(chan)
    prof  = chanprof(j)%prof
    IF (sun(j)) THEN
      IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
        DO lay = 1, nlayers
          DO ii = 1, 13
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%watervapour(lay, chan, ii) * predictors%watervapour_sun(ii, lay, prof)
          ENDDO
          IF (coef%nco2 > 0) THEN
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%watervapour(lay, chan, 14) * predictors%co2_sun(1, lay, prof)
          ENDIF
          IF (coef%nco > 0) THEN
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%watervapour(lay, chan, 15) * predictors%co_sun(1, lay, prof)
          ENDIF
        ENDDO
      ELSE IF (chanx >= rttov9_wv2295_25 .AND. chanx < rttov9_wv2360_00) THEN
        DO lay = 1, nlayers
          DO ii = 1, 9
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%watervapour(lay, chan, ii) * predictors%watervapour_sun(ii, lay, prof)
          ENDDO
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 10) * predictors%watervapour_sun(10, lay, prof)
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 11) * predictors%watervapour_sun(11, lay, prof)
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 12) * predictors%watervapour_sun(12, lay, prof)
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 13) * predictors%watervapour_sun(13, lay, prof)
        ENDDO
      ELSE
        DO lay = 1, nlayers
          DO ii = 1, 9
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%watervapour(lay, chan, ii) * predictors%watervapour_sun(ii, lay, prof)
          ENDDO
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 10) * predictors%watervapour_sun(14, lay, prof)
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 11) * predictors%watervapour_sun(11, lay, prof)
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 12) * predictors%watervapour_sun(12, lay, prof)
          opticaldepth(lay, j) =      &
            & opticaldepth(lay, j) + coef%watervapour(lay, chan, 13) * predictors%watervapour_sun(15, lay, prof)
          IF (coef%nch4 > 0) THEN
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%watervapour(lay, chan, 14) * predictors%ch4_sun(1, lay, prof) ** 1.25
            opticaldepth(lay, j) = opticaldepth(lay, j) +                                     &
              & coef%watervapour(lay, chan, 15) * predictors%ch4_sun(1, lay, prof) ** 0.25 *  &
              & predictors%ch4_sun(3, lay, prof)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!-------------
!1.3 add ozone
!-------------
  IF (coef%nozone > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = 1, nlayers
          DO ii = 1, 13
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%ozone(lay, chan, ii) * predictors%ozone_sun(ii, lay, prof)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!------------------------------
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = 1, nlayers
          DO ii = 1, coef%nwvcont
            opticaldepth(lay, j) =      &
              & opticaldepth(lay, j) + coef%wvcont(lay, chan, ii) * predictors%wvcont_sun(ii, lay, prof)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            DO ii = 1, 14
              opticaldepth(lay, j) =      &
                & opticaldepth(lay, j) + coef%co2(lay, chan, ii) * predictors%co2_sun(ii, lay, prof)
            ENDDO
          ENDDO
        ELSE
          DO lay = 1, nlayers
            DO ii = 1, 7
              opticaldepth(lay, j) =      &
                & opticaldepth(lay, j) + coef%co2(lay, chan, ii) * predictors%co2_sun(ii, lay, prof)
            ENDDO
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co2(lay, chan, 8) * predictors%co2_sun(9, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co2(lay, chan, 9) * predictors%co2_sun(10, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co2(lay, chan, 10) * predictors%co2_sun(11, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co2(lay, chan, 11) * predictors%co2_sun(12, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co2(lay, chan, 12) * predictors%co2_sun(13, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co2(lay, chan, 13) * predictors%co2_sun(14, lay, prof)
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.6 add N2O
!-----------
  IF (coef%nn2o > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      chanx = coef%ff_cwn(chan)
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        IF (chanx >= rttov9_wv2000_00 .AND. chanx < rttov9_wv2295_25) THEN
          DO lay = 1, nlayers
            DO ii = 1, 10
              opticaldepth(lay, j) =      &
                & opticaldepth(lay, j) + coef%n2o(lay, chan, ii) * predictors%n2o_sun(ii, lay, prof)
            ENDDO
            IF (coef%nco > 0) THEN
              opticaldepth(lay, j) = opticaldepth(lay, j) + coef%n2o(lay, chan, 11) * predictors%co_sun(1, lay, prof)
              opticaldepth(lay, j) = opticaldepth(lay, j) +      &
                & coef%n2o(lay, chan, 12) * predictors%co_sun(1, lay, prof) ** 3 / predictors%co_sun(8, lay, prof)
            ENDIF
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%n2o(lay, chan, 13) * predictors%n2o_sun(11, lay, prof)
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%n2o(lay, chan, 14) * predictors%n2o_sun(12, lay, prof)
          ENDDO
        ELSE
          DO lay = 1, nlayers
            DO ii = 1, 12
              opticaldepth(lay, j) =      &
                & opticaldepth(lay, j) + coef%n2o(lay, chan, ii) * predictors%n2o_sun(ii, lay, prof)
            ENDDO
          ENDDO
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.7 add CO
!-----------
  IF (coef%nco > 0) THEN
    DO j = 1, nchannels
      chan  = chanprof(j)%chan
      prof  = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = 1, nlayers
          DO ii = 1, 12
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%co(lay, chan, ii) * predictors%co_sun(ii, lay, prof)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!-----------
!1.8 add CH4
!-----------
  IF (coef%nch4 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      IF (sun(j)) THEN
        DO lay = 1, nlayers
          DO ii = 1, coef%nch4
            opticaldepth(lay, j) = opticaldepth(lay, j) + coef%ch4(lay, chan, ii) * predictors%ch4_sun(ii, lay, prof)
          ENDDO
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!----------------------------------------
!2. Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
! gamma correction done in rttov_transmit_9_solar
! single layer
! put raw opticaldepth in reference array for TL, AD and K calculations
! then check if opticaldepth is sensible and constrain if necessary.
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lay = 1, nlayers
        opdpsun_ref(lay, j) = opticaldepth(lay, j)
        IF (opticaldepth(lay, j) > 0.0_JPRB) THEN
          opticaldepth(lay, j) = 0.0_JPRB
        ENDIF
      ENDDO
    ENDIF
  ENDDO
! path to space
  opdp_path%sun_level = 0.0_jprb
  DO j = 1, nchannels
    IF (sun(j)) THEN
      opdp_path%sun_level(1, j) = 0.0_jprb
      DO lev = 2, nlevels
        lay = lev - 1
        opdp_path%sun_level(lev, j) = opdp_path%sun_level(lev - 1, j) + opticaldepth(lay, j)
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_9_SOLAR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_9_solar
