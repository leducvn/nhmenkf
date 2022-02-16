!
SUBROUTINE rttov_setpredictors_9_solar_k( &
            & opts,         &
            & chanprof,     &
            & prof,         &
            & prof_k,       &
            & coef,         &
            & predictors,   &
            & predictors_k, &
            & raytracing,   &
            & raytracing_k)
! Description
! RTTOV-8 Model
! K of rttov_setpredictors_88_solar
! To calculate and store the profile variables (predictors) required
! in subsequent transmittance calculations.
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
! see RTTOV7 science and validation report pages 18/19
! variable names are close to the documentation
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0   17/08/2006  Roger Saunders
!           --       New routine based on rttov_setpredictors_k.F90.
!           --       Variable trace gases CO2, N2O,CO and CH4
!           --       introduced for IASI and AIRS.
!           --       Altitude dependent local zenith angle also
!           --       introduced.
!  1.1   12/02/2007  Removed polarisation index (R Saunders)
!  1.2   15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.3   02/12/2009  Pathsat, Pathsun and related quantities are now
!                    layer arrays (Marco Matricardi).
!  1.4   05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.5   04/08/2010  Move addsolar check to calling routine (J Hocking)
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
       & rttov_options,   &
       & rttov_chanprof,  &
       & rttov_coef,      &
       & profile_Type,    &
       & predictors_Type, &
       & raytracing_type
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY : max_sol_zen
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options)  , INTENT(IN)    :: opts
  TYPE(rttov_chanprof ), INTENT(IN)    :: chanprof(:)
  TYPE(profile_Type   ), INTENT(IN)    :: prof    (:)
  TYPE(profile_Type   ), INTENT(INOUT) :: prof_k  (size(chanprof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_k
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer
  INTEGER(KIND=jpim) :: i, j
! user profile
  REAL   (KIND=Jprb) :: t(prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: w(prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: o(prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: co2  (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: co   (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: n2o  (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4  (prof(1)%nlayers, size(prof))
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: tro  (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: or   (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: co2r (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: n2or (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: cor  (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4r (prof(1)%nlayers, size(prof))
! user - reference
  REAL   (KIND=Jprb) :: dt   (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: dto  (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: dtabs(prof(1)%nlayers, size(prof))
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: twr  (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: tuw  (prof(1)%nlayers, size(prof))
  REAL   (KIND=Jprb) :: tuwr (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: ww   (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: ow   (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: co2w (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: n2ow (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: cow  (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4w (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: n2owr(prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: cowr (prof(1)%nlayers, size(prof))
  REAL   (KIND=jprb) :: ch4wr(prof(1)%nlayers, size(prof))
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: sum2_ww   (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_wwr  (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_ow   (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_twr  (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_tuw  (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_co2w (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_ch4w (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_ch4wr(prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_n2ow (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_n2owr(prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_cow  (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: sum2_cowr (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: tr_sq     (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=Jprb) :: tr_4      (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sec_wrwr  (prof(1)%nlayers, size(prof)    )
  REAL   (KIND=jprb) :: sec_wr    (prof(1)%nlayers, size(prof)    )
! K variables
  REAL   (KIND=Jprb) :: t_k       (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: w_k       (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: o_k       (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co2_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: n2o_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ch4_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tr_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tro_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: wr_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: or_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: wwr_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co2r_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: n2or_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: cor_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ch4r_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: twr_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tuw_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tuwr_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: dt_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: dto_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: tw_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ww_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ow_k      (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: co2w_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: n2ow_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: cow_k     (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ch4w_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: n2owr_k   (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: cowr_k    (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: ch4wr_k   (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: sec_or_k  (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: sec_wr_k  (prof(1)%nlayers, size(chanprof))
  REAL   (KIND=Jprb) :: sec_wrwr_k(prof(1)%nlayers, size(chanprof))
  INTEGER(KIND=jpim) :: nprofiles                          ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                          ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_SOLAR_K', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  nchannels = size(chanprof)
!-------------------------------------------------------------------------------
! Recompute Direct variables
!-------------------------------------------------------------------------------
!1) Profile layer quantities
!-------------------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N

  DO j = 1, nprofiles
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    DO layer = 1, prof(j)%nlayers
      level       = layer + 1
!-Temperature--------------------------------------------------------------------
      t(layer, j) = (prof(j)%t(level - 1) + prof(j)%t(level)) / 2._JPRB
!-H2O----------------------------------------------------------------------------
      w(layer, j) = (prof(j)%q(level - 1) + prof(j)%q(level)) / 2._JPRB
!-O3-----------------------------------------------------------------------------
      IF (opts%ozone_Data .AND. coef%nozone > 0) o(layer, j)   = (prof(j)%o3(level - 1) + prof(j)%o3(level)) / 2._JPRB
!-CO2----------------------------------------------------------------------------
      IF (opts%co2_Data .AND. coef%nco2 > 0    ) co2(layer, j) = (prof(j)%co2(level - 1) + prof(j)%co2(level)) / 2._JPRB
!-N2O----------------------------------------------------------------------------
      IF (opts%n2o_Data .AND. coef%nn2o > 0    ) n2o(layer, j) = (prof(j)%n2o(level - 1) + prof(j)%n2o(level)) / 2._JPRB
!-CO-----------------------------------------------------------------------------
      IF (opts%co_Data .AND. coef%nco > 0      ) co(layer, j)  = (prof(j)%co(level - 1) + prof(j)%co(level)) / 2._JPRB
!-CH4----------------------------------------------------------------------------
      IF (opts%ch4_Data .AND. coef%nch4 > 0    ) ch4(layer, j) = (prof(j)%ch4(level - 1) + prof(j)%ch4(level)) / 2._JPRB
    ENDDO

!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:, j)      = t(:, j) - coef%tstar(:)
    dtabs(:, j)   = Abs(dt(:, j))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr(:, j)      = t(:, j) / coef%tstar(:)
    tr_sq(:, j)   = tr(:, j) * tr(:, j)
    tr_4(:, j)    = tr_sq(:, j) * tr_sq(:, j)
    wr(:, j)      = w(:, j) / coef%wstar(:)

    IF (coef%nozone > 0) THEN
      dto(:, j) = t(:, j) - coef%to3star(:)
      tro(:, j) = t(:, j) / coef%to3star(:)
      IF (opts%ozone_Data) THEN
        or(:, j)  = o(:, j) / coef%ostar(:)
      ELSE
        or(:, j)  = 1._JPRB
      ENDIF
    ENDIF
    
!-CO2----------------------------------------------------------------------------

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      co2r(:, j) = co2(:, j) / coef%co2star(:)
    ELSE
      co2r(:, j) = 1._JPRB
    ENDIF

!-N2O----------------------------------------------------------------------------

    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or(:, j) = n2o(:, j) / coef%n2ostar(:)
    ELSE
      n2or(:, j) = 1._JPRB
    ENDIF

!-CO----------------------------------------------------------------------------

    IF (opts%co_Data .AND. coef%nco > 0) THEN
      cor(:, j) = co(:, j) / coef%costar(:)
    ELSE
      cor(:, j) = 1._JPRB
    ENDIF

!-CH4---------------------------------------------------------------------------

    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r(:, j) = ch4(:, j) / coef%ch4star(:)
    ELSE
      ch4r(:, j) = 1._JPRB
    ENDIF

!-------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw wwr
!--------------------------------------------------------------------
    tw(1, j) = 0._jprb

    DO layer = 2, prof(j)%nlayers
      tw(layer, j) = tw(layer - 1, j) + coef%dpp(layer - 1) * tr(layer - 1, j)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + t(layer, j)
      sum2 = sum2 + coef%tstar(layer)
      sum2_tuw(layer, j) = sum2
      tuw(layer, j)      = sum1 / sum2
    ENDDO

    tuwr(1, j)                 = coef%dpp(0) * t(1, j) / (coef%dpp(0) * coef%tstar(1))
    tuwr(2:prof(j)%nlayers, j) = tuw(2:prof(j)%nlayers, j)

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * t(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        sum2_twr(layer, j) = sum2
        twr(layer, j)      = sum1 / sum2
      ENDDO
    ENDIF

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + coef%dpp(layer - 1) * w(layer, j) * t(layer, j)
      sum2 = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      sum2_wwr(layer, j) = sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(j)%nlayers
      sum1 = sum1 + coef%dpp(layer - 1) * w(layer, j)
      sum2 = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      sum2_ww(layer, j) = sum2
      ww(layer, j)      = sum1 / sum2
    ENDDO


    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * o(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        sum2_ow(layer, j) = sum2
        ow(layer, j)      = sum1 / sum2
      ENDDO

    ELSE
      ow(:,:) = 1._JPRB
    ENDIF


    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(j)%nlayers
        sum1 = sum1 + coef%dpp(layer - 1) * co2(layer, j)
        sum2 = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        sum2_co2w(layer, j) = sum2
        co2w(layer, j)      = sum1 / sum2
      ENDDO

    ELSE
      co2w(:,:) = 1._JPRB
    ENDIF

!-N2O---------------------------------------------------------------------------

    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * n2o(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          sum2_n2ow(layer, j) = sum2
          n2ow(layer, j)      = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * n2o(layer, j) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          sum2_n2owr(layer, j) = sum2
          n2owr(layer, j)      = sum1 / sum2
        ENDDO
  
      ELSE
        n2ow(:,:)  = 1._JPRB

        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          sum2_n2owr(layer, j) = sum2
          n2owr(layer, j)      = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF

!-CO---------------------------------------------------------------------------

    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * co(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          sum2_cow(layer, j) = sum2
          cow(layer, j)      = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * co(layer, j) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          sum2_cowr(layer, j) = sum2
          cowr(layer, j)      = sum1 / sum2
        ENDDO
  
      ELSE
        cow(:,:)  = 1._JPRB
        
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * coef%costar(layer) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          sum2_cowr(layer, j) = sum2
          cowr(layer, j)      = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF
    
!-CH4---------------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * ch4(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          sum2_ch4w(layer, j) = sum2
          ch4w(layer, j)      = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * ch4(layer, j) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          sum2_ch4wr(layer, j) = sum2
          ch4wr(layer, j)      = sum1 / sum2
        ENDDO
  
      ELSE
        ch4w(:,:)  = 1._JPRB
        
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(j)%nlayers
          sum1 = sum1 + coef%dpp(layer - 1) * coef%ch4star(layer) * t(layer, j)
          sum2 = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          sum2_ch4wr(layer, j) = sum2
          ch4wr(layer, j)      = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF
    
  ENDDO

!-------------------------------------------------------------------------
! K code
!-------------------------------------------------------------------------
  w_k(:,:)        = 0._JPRB
  wr_k(:,:)       = 0._JPRB
  ww_k(:,:)       = 0._JPRB
  wwr_k(:,:)      = 0._JPRB
  sec_wr_k(:,:)   = 0._JPRB
  sec_wrwr_k(:,:) = 0._JPRB
  dt_k(:,:)       = 0._JPRB
  dto_k(:,:)      = 0._JPRB
  t_k(:,:)        = 0._JPRB
  tr_k(:,:)       = 0._JPRB
  tro_k(:,:)      = 0._JPRB
  tw_k(:,:)       = 0._JPRB
  twr_k(:,:)      = 0._JPRB
  tuw_k(:,:)      = 0._JPRB
  tuwr_k(:,:)     = 0._JPRB
  o_k(:,:)        = 0._JPRB
  or_k(:,:)       = 0._JPRB
  ow_k(:,:)       = 0._JPRB
  sec_or_k(:,:)   = 0._JPRB
  ch4_k(:,:)      = 0._JPRB
  ch4r_k(:,:)     = 0._JPRB
  ch4w_k(:,:)     = 0._JPRB
  ch4wr_k(:,:)    = 0._JPRB
  co_k(:,:)       = 0._JPRB
  cor_k(:,:)      = 0._JPRB
  cow_k(:,:)      = 0._JPRB
  cowr_k(:,:)     = 0._JPRB
  n2o_k(:,:)      = 0._JPRB
  n2or_k(:,:)     = 0._JPRB
  n2ow_k(:,:)     = 0._JPRB
  n2owr_k(:,:)    = 0._JPRB
  co2_k(:,:)      = 0._JPRB
  co2r_k(:,:)     = 0._JPRB
  co2w_k(:,:)     = 0._JPRB
!-------------------------------------------------
!5.9 ch4            transmittance based on RTIASI
!-------------------------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    IF (coef%nch4 > 0) THEN

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        ch4r_k(layer, i)               =      &
          & ch4r_k(layer, i) + predictors_k%ch4_sun(11, layer, i) * 1.5 * predictors%ch4_sun(2, layer, j) / ch4w(layer, j)
        ch4w_k(layer, i)               = ch4w_k(layer, i) -                              &
          & predictors_k%ch4_sun(11, layer, i) * predictors%ch4_sun(2, layer, j) ** 3 /  &
          & (raytracing%patheff(layer, j) * ch4w(layer, j) ** 2)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +       &
          & predictors_k%ch4_sun(11, layer, i) * 0.5 * ch4r(layer, j) ** 1.5 /  &
          & (raytracing%patheff(layer, j) ** 0.5 * ch4w(layer, j))
        ch4w_k(layer, i)               =      &
          & ch4w_k(layer, i) + raytracing%patheff(layer, j) * predictors_k%ch4_sun(10, layer, i)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ch4_sun(10, layer, i) * ch4w(layer, j)
        ch4w_k(layer, i)               = ch4w_k(layer, i) +      &
          & predictors_k%ch4_sun(9, layer, i) * predictors%ch4_sun(10, layer, j) * 2 * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ch4_sun(9, layer, i) * 2 * predictors%ch4_sun(10, layer, j) * ch4w(layer, j)
        ch4wr_k(layer, i)              = ch4wr_k(layer, i) + predictors_k%ch4_sun(8, layer, i)
        ch4wr_k(layer, i)              =      &
          & ch4wr_k(layer, i) + predictors_k%ch4_sun(7, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ch4_sun(7, layer, i) * ch4wr(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) +      &
          & predictors_k%ch4_sun(6, layer, i) * raytracing%patheff(layer, j) * 0.25 / predictors%ch4_sun(6, layer, j) ** 3
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ch4_sun(6, layer, i) * ch4r(layer, j) * 0.25 / predictors%ch4_sun(6, layer, j) ** 3
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%ch4_sun(5, layer, i) * ch4r(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) + predictors_k%ch4_sun(5, layer, i) * dt(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) +      &
          & predictors_k%ch4_sun(4, layer, i) * predictors%ch4_sun(1, layer, j) * raytracing%patheff(layer, j) * 2
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ch4_sun(4, layer, i) * 2 * predictors%ch4_sun(1, layer, j) * ch4r(layer, j)
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%ch4_sun(3, layer, i) * predictors%ch4_sun(1, layer, j)
        ch4r_k(layer, i)               =      &
          & ch4r_k(layer, i) + predictors_k%ch4_sun(3, layer, i) * raytracing%patheff(layer, j) * dt(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ch4_sun(3, layer, i) * ch4r(layer, j) * dt(layer, j)
        ch4r_k(layer, i)               = ch4r_k(layer, i) +      &
          & predictors_k%ch4_sun(2, layer, i) * 0.5 * raytracing%patheff(layer, j) / predictors%ch4_sun(2, layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ch4_sun(2, layer, i) * 0.5 * ch4r(layer, j) / predictors%ch4_sun(2, layer, j)
        ch4r_k(layer, i)               =      &
          & ch4r_k(layer, i) + predictors_k%ch4_sun(1, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ch4_sun(1, layer, i) * ch4r(layer, j)
      ENDDO

    ENDIF

  ENDDO

!------------------------------------------------
!5.8 CO             transmittance based on RTIASI
!-------------------------------------------------
!

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    IF (coef%nco > 0) THEN

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        cowr_k(layer, i)               = cowr_k(layer, i) +                            &
          & predictors_k%co_sun(12, layer, i) * 0.25 * raytracing%patheff(layer, j) *  &
          & (raytracing%patheff(layer, j) * cowr(layer, j)) ** ( - 0.75)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(12, layer, i) * 0.25 * cowr(layer, j) *        &
          & (raytracing%patheff(layer, j) * cowr(layer, j)) ** ( - 0.75)
        cowr_k(layer, i)               = cowr_k(layer, i) +                           &
          & predictors_k%co_sun(11, layer, i) * raytracing%patheff(layer, j) * 0.4 *  &
          & (raytracing%patheff(layer, j) * cowr(layer, j)) ** ( - 0.6)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(11, layer, i) * cowr(layer, j) * 0.4 *         &
          & (raytracing%patheff(layer, j) * cowr(layer, j)) ** ( - 0.6)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co_sun(10, layer, i) * 2 * predictors%co_sun(1, layer, j) / cow(layer, j) ** 0.5
        cow_k(layer, i)                =      &
          & cow_k(layer, i) - predictors_k%co_sun(10, layer, i) * 0.5 * predictors%co_sun(10, layer, j) / cow(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co_sun(10, layer, i) * cor(layer, j) ** 2 / SQRT(cow(layer, j))
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co_sun(9, layer, i) * 1.5 * predictors%co_sun(2, layer, j) / cow(layer, j)
        cow_k(layer, i)                = cow_k(layer, i) -                            &
          & predictors_k%co_sun(9, layer, i) * predictors%co_sun(2, layer, j) ** 3 /  &
          & (raytracing%patheff(layer, j) * cow(layer, j) ** 2)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(9, layer, i) * 0.5 * cor(layer, j) ** 1.5 /    &
          & (raytracing%patheff(layer, j) ** 0.5 * cow(layer, j))
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co_sun(8, layer, i) * 2 * predictors%co_sun(1, layer, j) / cow(layer, j)
        cow_k(layer, i)                =      &
          & cow_k(layer, i) - predictors_k%co_sun(8, layer, i) * predictors%co_sun(8, layer, j) / cow(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co_sun(8, layer, i) * cor(layer, j) ** 2 / cow(layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +      &
          & predictors_k%co_sun(7, layer, i) * dtabs(layer, j) * raytracing%patheff(layer, j) * dt(layer, j)
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%co_sun(7, layer, i) * dtabs(layer, j) * predictors%co_sun(1, layer, j) * 2
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(7, layer, i) * dtabs(layer, j) * cor(layer, j) * dt(layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +      &
          & predictors_k%co_sun(6, layer, i) * raytracing%patheff(layer, j) * 0.25 / predictors%co_sun(6, layer, j) ** 3
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(6, layer, i) * cor(layer, j) * 0.25 / predictors%co_sun(6, layer, j) ** 3
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%co_sun(5, layer, i) * predictors%co_sun(2, layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +                                          &
          & predictors_k%co_sun(5, layer, i) * 0.5 * dt(layer, j) * raytracing%patheff(layer, j) /  &
          & predictors%co_sun(2, layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(5, layer, i) * 0.5 * dt(layer, j) * cor(layer, j) / predictors%co_sun(2, layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +      &
          & predictors_k%co_sun(4, layer, i) * 2 * predictors%co_sun(1, layer, j) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(4, layer, i) * 2 * predictors%co_sun(1, layer, j) * cor(layer, j)
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%co_sun(3, layer, i) * predictors%co_sun(1, layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co_sun(3, layer, i) * dt(layer, j) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co_sun(3, layer, i) * cor(layer, j) * dt(layer, j)
        cor_k(layer, i)                = cor_k(layer, i) +      &
          & predictors_k%co_sun(2, layer, i) * 0.5 * raytracing%patheff(layer, j) / predictors%co_sun(2, layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co_sun(2, layer, i) * 0.5 * cor(layer, j) / predictors%co_sun(2, layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co_sun(1, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co_sun(1, layer, i) * cor(layer, j)
      ENDDO

    ENDIF

  ENDDO

!-------------------------------------------------
!5.7 N2O            transmittance based on RTIASI
!-------------------------------------------------
!

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    IF (coef%nn2o > 0) THEN

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%n2o_sun(12, layer, i) * raytracing%patheff(layer, j) ** 2 * n2owr(layer, j)
        n2owr_k(layer, i)              =      &
          & n2owr_k(layer, i) + predictors_k%n2o_sun(12, layer, i) * raytracing%patheff(layer, j) ** 2 * dt(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%n2o_sun(12, layer, i) * 2 * raytracing%patheff(layer, j) * n2owr(layer, j) * dt(layer, j)
        n2owr_k(layer, i)              = n2owr_k(layer, i) +      &
          & predictors_k%n2o_sun(11, layer, i) * 3 * predictors%n2o_sun(8, layer, j) ** 2 * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%n2o_sun(11, layer, i) * 3 * predictors%n2o_sun(8, layer, j) ** 2 * n2owr(layer, j)
        n2or_k(layer, i)               =      &
          & n2or_k(layer, i) + predictors_k%n2o_sun(10, layer, i) * 1.5 * predictors%n2o_sun(2, layer, j) / n2ow(layer, j)
        n2ow_k(layer, i)               = n2ow_k(layer, i) -                              &
          & predictors_k%n2o_sun(10, layer, i) * predictors%n2o_sun(2, layer, j) ** 3 /  &
          & (raytracing%patheff(layer, j) * n2ow(layer, j) ** 2)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +       &
          & predictors_k%n2o_sun(10, layer, i) * 0.5 * n2or(layer, j) ** 1.5 /  &
          & (raytracing%patheff(layer, j) ** 0.5 * n2ow(layer, j))
        n2owr_k(layer, i)              = n2owr_k(layer, i) +      &
          & predictors_k%n2o_sun(9, layer, i) * 2 * predictors%n2o_sun(8, layer, j) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%n2o_sun(9, layer, i) * predictors%n2o_sun(8, layer, j) * 2 * n2owr(layer, j)
        n2owr_k(layer, i)              =      &
          & n2owr_k(layer, i) + predictors_k%n2o_sun(8, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%n2o_sun(8, layer, i) * n2owr(layer, j)
        n2ow_k(layer, i)               =      &
          & n2ow_k(layer, i) + predictors_k%n2o_sun(7, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%n2o_sun(7, layer, i) * n2ow(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) +      &
          & predictors_k%n2o_sun(6, layer, i) * raytracing%patheff(layer, j) * 0.25 / predictors%n2o_sun(6, layer, j) ** 3
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%n2o_sun(6, layer, i) * n2or(layer, j) * 0.25 / predictors%n2o_sun(6, layer, j) ** 3
        dt_k(layer, i)                 = dt_k(layer, i) + predictors_k%n2o_sun(5, layer, i) * n2or(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) + predictors_k%n2o_sun(5, layer, i) * dt(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) +      &
          & predictors_k%n2o_sun(4, layer, i) * 2 * predictors%n2o_sun(1, layer, j) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%n2o_sun(4, layer, i) * 2 * predictors%n2o_sun(1, layer, j) * n2or(layer, j)
        dt_k(layer, i)                 =      &
          & dt_k(layer, i) + predictors_k%n2o_sun(3, layer, i) * predictors%n2o_sun(1, layer, j)
        n2or_k(layer, i)               =      &
          & n2or_k(layer, i) + predictors_k%n2o_sun(3, layer, i) * raytracing%patheff(layer, j) * dt(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%n2o_sun(3, layer, i) * n2or(layer, j) * dt(layer, j)
        n2or_k(layer, i)               = n2or_k(layer, i) +      &
          & predictors_k%n2o_sun(2, layer, i) * 0.5 * raytracing%patheff(layer, j) / predictors%n2o_sun(2, layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%n2o_sun(2, layer, i) * 0.5 * n2or(layer, j) / predictors%n2o_sun(2, layer, j)
        n2or_k(layer, i)               =      &
          & n2or_k(layer, i) + predictors_k%n2o_sun(1, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%n2o_sun(1, layer, i) * n2or(layer, j)
      ENDDO

    ENDIF

  ENDDO

!------------------------------------------------------------------------------------
!5.6 CO2
!------------------------------------------------------------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    IF (coef%nco2 > 0) THEN

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2_sun(14, layer, i) * 2 * SQRT(predictors%co2_sun(14, layer, j)) * twr(layer, j)
        twr_k(layer, i)                =      &
          & twr_k(layer, i) + predictors_k%co2_sun(14, layer, i) * 2 * SQRT(predictors%co2_sun(14, layer, j)) * tr(layer, j)
        twr_k(layer, i)                = twr_k(layer, i) +                                     &
          & predictors_k%co2_sun(13, layer, i) * 3 * twr(layer, j) ** 2 * tr(layer, j) ** 2 *  &
          & raytracing%patheff(layer, j) ** 0.5
        tr_k(layer, i)                 = tr_k(layer, i) +      &
          & predictors_k%co2_sun(13, layer, i) * 2 * tr(layer, j) * twr(layer, j) ** 3 * SQRT(raytracing%patheff(layer, j))
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +                        &
          & predictors_k%co2_sun(13, layer, i) * tr(layer, j) ** 2 * twr(layer, j) ** 3 * 0.5 /  &
          & SQRT(raytracing%patheff(layer, j))
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2_sun(12, layer, i) * 3 * tr(layer, j) ** 2 * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(12, layer, i) * tr(layer, j) ** 3
        tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%co2_sun(11, layer, i) * 3 * tr(layer, j) ** 2
        co2r_k(layer, i)               = co2r_k(layer, i) +      &
          & predictors_k%co2_sun(10, layer, i) * 0.5 * raytracing%patheff(layer, j) / SQRT(predictors%co2_sun(1, layer, j))
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co2_sun(10, layer, i) * 0.5 * co2r(layer, j) / SQRT(predictors%co2_sun(1, layer, j))
        twr_k(layer, i)                = twr_k(layer, i) +      &
          & predictors_k%co2_sun(9, layer, i) * raytracing%patheff(layer, j) * SQRT(predictors%co2_sun(5, layer, j))
        tr_k(layer, i)                 = tr_k(layer, i) +                                &
          & predictors_k%co2_sun(9, layer, i) * 0.5 * predictors%co2_sun(6, layer, j) /  &
          & SQRT(predictors%co2_sun(5, layer, j))
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(9, layer, i) * twr(layer, j) * tr(layer, j) ** 0.5
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(8, layer, i) * cor(layer, j)
        cor_k(layer, i)                =      &
          & cor_k(layer, i) + predictors_k%co2_sun(8, layer, i) * raytracing%patheff(layer, j)
        co2w_k(layer, i)               = co2w_k(layer, i) +      &
          & predictors_k%co2_sun(7, layer, i) * 2 * raytracing%patheff(layer, j) * SQRT(predictors%co2_sun(7, layer, j))
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%co2_sun(7, layer, i) * 2 * co2w(layer, j) * SQRT(predictors%co2_sun(7, layer, j))
        twr_k(layer, i)                =      &
          & twr_k(layer, i) + predictors_k%co2_sun(6, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(6, layer, i) * twr(layer, j)
        tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%co2_sun(5, layer, i)
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2_sun(4, layer, i) * 2 * predictors%co2_sun(3, layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(4, layer, i) * tr(layer, j) ** 2
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2_sun(3, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(3, layer, i) * tr(layer, j)
        tr_k(layer, i)                 =      &
          & tr_k(layer, i) + predictors_k%co2_sun(2, layer, i) * 2 * predictors%co2_sun(5, layer, j)
        co2r_k(layer, i)               =      &
          & co2r_k(layer, i) + predictors_k%co2_sun(1, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%co2_sun(1, layer, i) * co2r(layer, j)
      ENDDO

    ENDIF

  ENDDO

!----------------------------------------------------------------------------------
!5.4 ozone
!----------------------------------------------------------------------------------
! One can pack all ow_k lines in one longer statement
! same for sec_or_k and dto_k

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    IF (coef%nozone > 0) THEN

      DO layer = 1, prof_k(i)%nlayers
        level = layer + 1
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ozone_sun(13, layer, i) * tro(layer, j) ** 3
        tro_k(layer, i)                =      &
          & tro_k(layer, i) + predictors_k%ozone_sun(13, layer, i) * 3 * tro(layer, j) ** 2 * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +                                                &
          & predictors_k%ozone_sun(12, layer, i) * 0.5 * raytracing%patheff(layer, j) ** ( - 0.5) * ow(layer, j) ** 2 *  &
          & dto(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) +      &
          & predictors_k%ozone_sun(12, layer, i) * 2 * ow(layer, j) * raytracing%patheff(layer, j) ** 0.5 * dto(layer, j)
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone_sun(12, layer, i) * raytracing%patheff(layer, j) ** 0.5 * ow(layer, j) ** 2
        ow_k(layer, i)                 = ow_k(layer, i) +      &
          & predictors_k%ozone_sun(11, layer, i) * 2 * raytracing%patheff(layer, j) * predictors%ozone_sun(10, layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ozone_sun(11, layer, i) * 2 * ow(layer, j) * predictors%ozone_sun(10, layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) + predictors_k%ozone_sun(10, layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ozone_sun(10, layer, i) * ow(layer, j)
        or_k(layer, i)                 = or_k(layer, i) +                                              &
          & predictors_k%ozone_sun(9, layer, i) * SQRT(raytracing%patheff(layer, j) * ow(layer, j)) *  &
          & raytracing%patheff(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) +                                                                   &
          & predictors_k%ozone_sun(9, layer, i) * predictors%ozone_sun(1, layer, j) * 0.5 * raytracing%patheff(layer, j) /  &
          & SQRT(raytracing%patheff(layer, j) * ow(layer, j))
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ozone_sun(9, layer, i) * 1.5 * or(layer, j) * SQRT(predictors%ozone_sun(10, layer, j))
        ow_k(layer, i)                 = ow_k(layer, i) +                                &
          & predictors_k%ozone_sun(8, layer, i) * 1.75 * raytracing%patheff(layer, j) *  &
          & (raytracing%patheff(layer, j) * ow(layer, j)) ** 0.75
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ozone_sun(8, layer, i) * 1.75 * ow(layer, j) *        &
          & (raytracing%patheff(layer, j) * ow(layer, j)) ** 0.75
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone_sun(7, layer, i) * 1.5 * predictors%ozone_sun(2, layer, j) / ow(layer, j)
        ow_k(layer, i)                 =      &
          & ow_k(layer, i) - predictors_k%ozone_sun(7, layer, i) * predictors%ozone_sun(7, layer, j) / ow(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
          & predictors_k%ozone_sun(7, layer, i) * 0.5 * or(layer, j) ** 1.5 /  &
          & (raytracing%patheff(layer, j) ** 0.5 * ow(layer, j))
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone_sun(6, layer, i) * 2 * predictors%ozone_sun(1, layer, j) * ow(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) +      &
          & predictors_k%ozone_sun(6, layer, i) * predictors%ozone_sun(4, layer, j) / raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ozone_sun(6, layer, i) * or(layer, j) ** 2 * ow(layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) +      &
          & predictors_k%ozone_sun(5, layer, i) * 0.5_JPRB * dto(layer, j) / (predictors%ozone_sun(2, layer, j))
        dto_k(layer, i)                =      &
          & dto_k(layer, i) + predictors_k%ozone_sun(5, layer, i) * predictors%ozone_sun(2, layer, j)
        sec_or_k(layer, i)             =      &
          & sec_or_k(layer, i) + predictors_k%ozone_sun(4, layer, i) * 2 * predictors%ozone_sun(1, layer, j)
        raytracing_k%patheff(layer, i) =      &
          & raytracing_k%patheff(layer, i) + predictors_k%ozone_sun(3, layer, i) * or(layer, j) / ow(layer, j)
        or_k(layer, i)                 =      &
          & or_k(layer, i) + predictors_k%ozone_sun(3, layer, i) * raytracing%patheff(layer, j) / ow(layer, j)
        ow_k(layer, i)                 = ow_k(layer, i) -      &
          & predictors_k%ozone_sun(3, layer, i) * or(layer, j) * raytracing%patheff(layer, j) / ow(layer, j) ** 2
        sec_or_k(layer, i)             =      &
          & sec_or_k(layer, i) + predictors_k%ozone_sun(2, layer, i) * 0.5_JPRB / predictors%ozone_sun(2, layer, j)
        sec_or_k(layer, i)             = sec_or_k(layer, i) + predictors_k%ozone_sun(1, layer, i)
        or_k(layer, i)                 = or_k(layer, i) + sec_or_k(layer, i) * raytracing%patheff(layer, j)
        raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) + sec_or_k(layer, i) * or(layer, j)
      ENDDO

    ENDIF

  ENDDO

!-------------------------------------------------
!5.3 Water Vapour Continuum based on RTIASI
!-------------------------------------------------

  IF (coef%nwvcont > 0) THEN

    DO i = 1, nchannels
      j = chanprof(i)%prof
      IF (prof(j)%sunzenangle < 0.0 .OR. &
          prof(j)%sunzenangle >= max_sol_zen) CYCLE

      DO layer = 1, prof_k(i)%nlayers
        level                = layer + 1
        sec_wr_k(layer, i)   = sec_wr_k(layer, i) + predictors_k%wvcont_sun(4, layer, i) / tr_sq(layer, j)
        tr_k(layer, i)       = tr_k(layer, i) -                                                   &
          & 2 * predictors_k%wvcont_sun(4, layer, i) * predictors%watervapour_sun(7, layer, j) /  &
          & (tr_sq(layer, j) * tr(layer, j))
        sec_wr_k(layer, i)   = sec_wr_k(layer, i) + predictors_k%wvcont_sun(3, layer, i) / tr(layer, j)
        tr_k(layer, i)       = tr_k(layer, i) -      &
          & predictors_k%wvcont_sun(3, layer, i) * predictors%watervapour_sun(7, layer, j) / tr_sq(layer, j)
        sec_wrwr_k(layer, i) = sec_wrwr_k(layer, i) + predictors_k%wvcont_sun(2, layer, i) / tr_4(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - 4 * predictors_k%wvcont_sun(2, layer, i) * predictors%wvcont_sun(1, layer, j) / tr_4(layer, j)
        sec_wrwr_k(layer, i) = sec_wrwr_k(layer, i) + predictors_k%wvcont_sun(1, layer, i) / tr(layer, j)
        tr_k(layer, i)       =      &
          & tr_k(layer, i) - predictors_k%wvcont_sun(1, layer, i) * predictors%wvcont_sun(1, layer, j) / tr(layer, j)
      ENDDO

    ENDDO

  ENDIF

!------------------------------------------------------------------
!5.2 water vapour based on RTIASI
!------------------------------------------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    DO layer = 1, prof_k(i)%nlayers
      level = layer + 1
      sec_wr(layer, j)               = raytracing%patheff(layer, j) * wr(layer, j)
      sec_wrwr(layer, j)             = sec_wr(layer, j) * wr(layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) +      &
        & predictors_k%watervapour_sun(15, layer, i) * 2 * predictors%watervapour_sun(7, layer, j) / ww(layer, j)
      ww_k(layer, i)                 = ww_k(layer, i) -                                           &
        & predictors_k%watervapour_sun(15, layer, i) * predictors%watervapour_sun(1, layer, j) /  &
        & (raytracing%patheff(layer, j) * ww(layer, j) ** 2)
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%watervapour_sun(15, layer, i) * wr(layer, j) ** 2 / ww(layer, j)
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour_sun(14, layer, i) * predictors%watervapour_sun(7, layer, j) ** 1.5
      wr_k(layer, i)                 = wr_k(layer, i) +                                                 &
        & predictors_k%watervapour_sun(14, layer, i) * 1.5 * predictors%watervapour_sun(5, layer, j) *  &
        & raytracing%patheff(layer, j) * dt(layer, j)
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +                                                &
        & predictors_k%watervapour_sun(14, layer, i) * 1.5 * predictors%watervapour_sun(5, layer, j) * wr(layer, j) *  &
        & dt(layer, j)
      wr_k(layer, i)                 =      &
        & wr_k(layer, i) + predictors_k%watervapour_sun(13, layer, i) * 1.5 * predictors%watervapour_sun(5, layer, j)
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
        & predictors_k%watervapour_sun(13, layer, i) * 0.5 * wr(layer, j) ** 1.5 / raytracing%patheff(layer, j) ** 0.5
      ww_k(layer, i)                 = ww_k(layer, i) +                                       &
        & predictors_k%watervapour_sun(12, layer, i) * 1.25 * raytracing%patheff(layer, j) *  &
        & (predictors%watervapour_sun(2, layer, j)) ** 0.25
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +       &
        & predictors_k%watervapour_sun(12, layer, i) * 1.25 * ww(layer, j) *  &
        & (predictors%watervapour_sun(2, layer, j)) ** 0.25
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour_sun(11, layer, i) * predictors%watervapour_sun(5, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) +      &
        & 0.5_JPRB * predictors_k%watervapour_sun(11, layer, i) * dt(layer, j) / predictors%watervapour_sun(5, layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + predictors_k%watervapour_sun(10, layer, i) * dtabs(layer, j) * dt(layer, j)
      dt_k(layer, i)                 = dt_k(layer, i) +      &
        & 2 * predictors_k%watervapour_sun(10, layer, i) * predictors%watervapour_sun(7, layer, j) * dtabs(layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) +                                     &
        & predictors_k%watervapour_sun(9, layer, i) * 1.5 * raytracing%patheff(layer, j) *  &
        & SQRT(predictors%watervapour_sun(7, layer, j))
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
        & predictors_k%watervapour_sun(9, layer, i) * 1.5 * wr(layer, j) * SQRT(predictors%watervapour_sun(7, layer, j))
      ww_k(layer, i)                 = ww_k(layer, i) +                                     &
        & predictors_k%watervapour_sun(8, layer, i) * 1.5 * raytracing%patheff(layer, j) *  &
        & SQRT(predictors%watervapour_sun(2, layer, j))
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
        & predictors_k%watervapour_sun(8, layer, i) * 1.5 * ww(layer, j) * SQRT(predictors%watervapour_sun(2, layer, j))
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + predictors_k%watervapour_sun(7, layer, i)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) +      &
        & 0.25_JPRB * predictors_k%watervapour_sun(6, layer, i) / predictors%watervapour_sun(6, layer, j) ** 3
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) +      &
        & 0.5_JPRB * predictors_k%watervapour_sun(5, layer, i) / predictors%watervapour_sun(5, layer, j)
      dt_k(layer, i)                 =      &
        & dt_k(layer, i) + predictors_k%watervapour_sun(4, layer, i) * predictors%watervapour_sun(7, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + predictors_k%watervapour_sun(4, layer, i) * dt(layer, j)
      ww_k(layer, i)                 = ww_k(layer, i) +                                              &
        & predictors_k%watervapour_sun(3, layer, i) * 2 * predictors%watervapour_sun(2, layer, j) *  &
        & raytracing%patheff(layer, j)
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
        & predictors_k%watervapour_sun(3, layer, i) * 2 * predictors%watervapour_sun(2, layer, j) * ww(layer, j)
      ww_k(layer, i)                 =      &
        & ww_k(layer, i) + predictors_k%watervapour_sun(2, layer, i) * raytracing%patheff(layer, j)
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%watervapour_sun(2, layer, i) * ww(layer, j)
      sec_wr_k(layer, i)             =      &
        & sec_wr_k(layer, i) + 2 * predictors_k%watervapour_sun(1, layer, i) * predictors%watervapour_sun(7, layer, j)
      sec_wr_k(layer, i)             = sec_wr_k(layer, i) + sec_wrwr_k(layer, i) * wr(layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) + sec_wrwr_k(layer, i) * predictors%watervapour_sun(7, layer, j)
      wr_k(layer, i)                 = wr_k(layer, i) + sec_wr_k(layer, i) * raytracing%patheff(layer, j)
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) + sec_wr_k(layer, i) * wr(layer, j)
    ENDDO

  ENDDO

!-------------------------------------------------------------------------------------------
!5.1 mixed gases
!-------------------------------------------------------------------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    DO layer = 1, prof_k(i)%nlayers
      level = layer + 1
! X10
      tr_k(layer, i)                 = tr_k(layer, i) +      &
        & predictors_k%mixedgas_sun(10, layer, i) * 0.5 * raytracing%patheff(layer, j) ** 1.5 / SQRT(tr(layer, j))
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
        & predictors_k%mixedgas_sun(10, layer, i) * 1.5 * SQRT(predictors%mixedgas_sun(3, layer, j))
! X9
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas_sun(9, layer, i) * 3 * raytracing%patheff(layer, j) * tr(layer, j) ** 2
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%mixedgas_sun(9, layer, i) * tr(layer, j) ** 3
! X8
      tuwr_k(layer, i)               =      &
        & tuwr_k(layer, i) + predictors_k%mixedgas_sun(8, layer, i) * raytracing%patheff(layer, j)
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%mixedgas_sun(8, layer, i) * tuwr(layer, j)
! X7
      tuw_k(layer, i)                =      &
        & tuw_k(layer, i) + predictors_k%mixedgas_sun(7, layer, i) * raytracing%patheff(layer, j)
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%mixedgas_sun(7, layer, i) * tuw(layer, j)
! X6
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas_sun(6, layer, i) * 2 * predictors%mixedgas_sun(5, layer, j)
! X5
      tr_k(layer, i)                 = tr_k(layer, i) + predictors_k%mixedgas_sun(5, layer, i)
! X4
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas_sun(4, layer, i) * 2 * predictors%mixedgas_sun(3, layer, j)
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) +      &
        & predictors_k%mixedgas_sun(4, layer, i) * predictors%mixedgas_sun(5, layer, j) ** 2
! X3
      tr_k(layer, i)                 =      &
        & tr_k(layer, i) + predictors_k%mixedgas_sun(3, layer, i) * raytracing%patheff(layer, j)
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%mixedgas_sun(3, layer, i) * tr(layer, j)
! X2
      raytracing_k%patheff(layer, i) =      &
        & raytracing_k%patheff(layer, i) + predictors_k%mixedgas_sun(2, layer, i) * 2 * raytracing%patheff(layer, j)
! X1
      raytracing_k%patheff(layer, i) = raytracing_k%patheff(layer, i) + predictors_k%mixedgas_sun(1, layer, i)
      raytracing_k%pathsat(layer, i) = raytracing_k%pathsat(layer, i) + raytracing_k%patheff(layer, i)
      raytracing_k%pathsun(layer, i) = raytracing_k%pathsun(layer, i) + raytracing_k%patheff(layer, i)
    ENDDO

  ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile/reference sums
!-------------------------------------------------------------------

  DO i = 1, nchannels
    j = chanprof(i)%prof
    IF (prof(j)%sunzenangle < 0.0 .OR. &
        prof(j)%sunzenangle >= max_sol_zen) CYCLE

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + ch4w_k(layer, i) / sum2_ch4w(layer, j)
          ch4_k(layer, i) = ch4_k(layer, i) + sum1 * coef%dpp(layer - 1)
        ENDDO
  
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + ch4wr_k(layer, i) / sum2_ch4wr(layer, j)
          ch4_k(layer, i) = ch4_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * ch4(layer, j)
        ENDDO
  
      ELSE
        ch4_k(:, i) = 0._JPRB
        
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + ch4wr_k(layer, i) / sum2_ch4wr(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * coef%ch4star(layer)
        ENDDO

      ENDIF
    ENDIF
    

    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1           = sum1 + cow_k(layer, i) / sum2_cow(layer, j)
          co_k(layer, i) = co_k(layer, i) + sum1 * coef%dpp(layer - 1)
        ENDDO
  
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1           = sum1 + cowr_k(layer, i) / sum2_cowr(layer, j)
          co_k(layer, i) = co_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
          t_k(layer, i)  = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * co(layer, j)
        ENDDO
  
      ELSE
        co_k(:, i) = 0._JPRB
        
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1           = sum1 + cowr_k(layer, i) / sum2_cowr(layer, j)
          t_k(layer, i)  = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * coef%costar(layer)
        ENDDO

      ENDIF
    ENDIF
    

    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + n2ow_k(layer, i) / sum2_n2ow(layer, j)
          n2o_k(layer, i) = n2o_k(layer, i) + sum1 * coef%dpp(layer - 1)
        ENDDO
  
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + n2owr_k(layer, i) / sum2_n2owr(layer, j)
          n2o_k(layer, i) = n2o_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * n2o(layer, j)
        ENDDO
  
      ELSE
        n2o_k(:, i) = 0._JPRB
        
        sum1 = 0._JPRB
  
        DO layer = prof_k(i)%nlayers, 1,  - 1
          sum1            = sum1 + n2owr_k(layer, i) / sum2_n2owr(layer, j)
          t_k(layer, i)   = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * coef%n2ostar(layer)
        ENDDO

      ENDIF
    ENDIF
    

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1            = sum1 + co2w_k(layer, i) / sum2_co2w(layer, j)
        co2_k(layer, i) = co2_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      co2_k(:, i) = 0._JPRB
    ENDIF

!
    sum1 = 0._JPRB

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1          = sum1 + ow_k(layer, i) / sum2_ow(layer, j)
        o_k(layer, i) = o_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO

    ELSE
      o_k(:, i) = 0._JPRB
    ENDIF

!
    sum1 = 0._JPRB

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + wwr_k(layer, i) / sum2_wwr(layer, j)
      w_k(layer, i) = w_k(layer, i) + sum1 * coef%dpp(layer - 1) * t(layer, j)
      t_k(layer, i) = t_k(layer, i) + sum1 * coef%dpp(layer - 1) * w(layer, j)
    ENDDO

!
    sum1 = 0._JPRB

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + ww_k(layer, i) / sum2_ww(layer, j)
      w_k(layer, i) = w_k(layer, i) + sum1 * coef%dpp(layer - 1)
    ENDDO

    sum1 = 0._JPRB

    IF (coef%nco2 > 0) THEN

      DO layer = prof_k(i)%nlayers, 1,  - 1
        sum1          = sum1 + twr_k(layer, i) / sum2_twr(layer, j)
        t_k(layer, i) = t_k(layer, i) + sum1 * coef%dpp(layer - 1)
      ENDDO
    
    ENDIF

    tuw_k(2:prof(j)%nlayers, i) = tuw_k(2:prof(j)%nlayers, i) + tuwr_k(2:prof(j)%nlayers, i)
    t_k(1, i)                   = t_k(1, i) + tuwr_k(1, i) * coef%dpp(0) / (coef%dpp(0) * coef%tstar(1))
    sum1 = 0._JPRB

    DO layer = prof_k(i)%nlayers, 1,  - 1
      sum1          = sum1 + tuw_k(layer, i) / sum2_tuw(layer, j)
      t_k(layer, i) = t_k(layer, i) + sum1
    ENDDO


    DO layer = prof_k(i)%nlayers, 2,  - 1
      tw_k(layer - 1, i) = tw_k(layer - 1, i) + tw_k(layer, i)
      tr_k(layer - 1, i) = tr_k(layer - 1, i) + tw_k(layer, i) * coef%dpp(layer - 1)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile deviations
!-------------------------------------------------------------------

    DO layer = 1, prof_k(i)%nlayers

      IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
        ch4_k(layer, i) = ch4_k(layer, i) + ch4r_k(layer, i) / coef%ch4star(layer)
      ELSE
        ch4_k(layer, i) = 0._JPRB
      ENDIF


      IF (opts%co_Data .AND. coef%nco > 0) THEN
        co_k(layer, i) = co_k(layer, i) + cor_k(layer, i) / coef%costar(layer)
      ELSE
        co_k(layer, i) = 0._JPRB
      ENDIF


      IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
        n2o_k(layer, i) = n2o_k(layer, i) + n2or_k(layer, i) / coef%n2ostar(layer)
      ELSE
        n2o_k(layer, i) = 0._JPRB
      ENDIF


      IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
        co2_k(layer, i) = co2_k(layer, i) + co2r_k(layer, i) / coef%co2star(layer)
      ELSE
        co2_k(layer, i) = 0._JPRB
      ENDIF


      IF (coef%nozone > 0) THEN
        t_k(layer, i) = t_k(layer, i) + tro_k(layer, i) / coef%to3star(layer)
        IF (opts%ozone_Data) THEN
          o_k(layer, i) = o_k(layer, i) + or_k(layer, i) / coef%ostar(layer)
        ELSE
          o_k(layer, i) = 0._JPRB
        ENDIF
      ENDIF

      w_k(layer, i) = w_k(layer, i) + wr_k(layer, i) / coef%wstar(layer)
      t_k(layer, i) = t_k(layer, i) + tr_k(layer, i) / coef%tstar(layer)

      IF (coef%nozone > 0) t_k(layer, i) = t_k(layer, i) + dto_k(layer, i)

      t_k(layer, i) = t_k(layer, i) + dt_k(layer, i)
    ENDDO

!-------------------------------------------------------------------
!   calc adjoint of profile layer means
!-------------------------------------------------------------------

    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN

      DO level = 2, prof_k(i)%nlevels
        layer = level - 1
        prof_k(i)%ch4(level - 1) = prof_k(i)%ch4(level - 1) + 0.5_JPRB * ch4_k(layer, i)
        prof_k(i)%ch4(level)     = prof_k(i)%ch4(level) + 0.5_JPRB * ch4_k(layer, i)
      ENDDO

    ENDIF


    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN

      DO level = 2, prof_k(i)%nlevels
        layer = level - 1
        prof_k(i)%n2o(level - 1) = prof_k(i)%n2o(level - 1) + 0.5_JPRB * n2o_k(layer, i)
        prof_k(i)%n2o(level)     = prof_k(i)%n2o(level) + 0.5_JPRB * n2o_k(layer, i)
      ENDDO

    ENDIF


    IF (opts%co_Data .AND. coef%nco > 0) THEN

      DO level = 2, prof_k(i)%nlevels
        layer = level - 1
        prof_k(i)%co(level - 1) = prof_k(i)%co(level - 1) + 0.5_JPRB * co_k(layer, i)
        prof_k(i)%co(level)     = prof_k(i)%co(level) + 0.5_JPRB * co_k(layer, i)
      ENDDO

    ENDIF


    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN

      DO level = 2, prof_k(i)%nlevels
        layer = level - 1
        prof_k(i)%co2(level - 1) = prof_k(i)%co2(level - 1) + 0.5_JPRB * co2_k(layer, i)
        prof_k(i)%co2(level)     = prof_k(i)%co2(level) + 0.5_JPRB * co2_k(layer, i)
      ENDDO

    ENDIF

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN

      DO level = 2, prof_k(i)%nlevels
        layer = level - 1
        prof_k(i)%o3(level - 1) = prof_k(i)%o3(level - 1) + 0.5_JPRB * o_k(layer, i)
        prof_k(i)%o3(level)     = prof_k(i)%o3(level) + 0.5_JPRB * o_k(layer, i)
      ENDDO

    ENDIF

    DO level = 2, prof_k(i)%nlevels
      layer = level - 1
      prof_k(i)%q(level - 1) = prof_k(i)%q(level - 1) + 0.5_JPRB * w_k(layer, i)
      prof_k(i)%q(level)     = prof_k(i)%q(level) + 0.5_JPRB * w_k(layer, i)
      prof_k(i)%t(level - 1) = prof_k(i)%t(level - 1) + 0.5_JPRB * t_k(layer, i)
      prof_k(i)%t(level)     = prof_k(i)%t(level) + 0.5_JPRB * t_k(layer, i)
    ENDDO

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_SOLAR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_9_solar_k
