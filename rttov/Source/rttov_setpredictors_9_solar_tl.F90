!
SUBROUTINE rttov_setpredictors_9_solar_tl( &
            & opts,          &
            & prof,          &
            & prof_tl,       &
            & coef,          &
            & predictors,    &
            & predictors_tl, &
            & raytracing,    &
            & raytracing_tl)
! Description
! RTTOV-8 Model
! TL of rttov_setpredictors_8
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
!  1.0   01/06/2005  Marco Matricardi (ECMWF):
!           --       New routine based on rttov_setpredictors_tl.F90.
!           --       Variable trace gases CO2, N2O,CO and CH4
!           --       introduced for IASI and AIRS.
!           --       Altitude dependent local zenith angle also
!           --       introduced.
!  1.1   15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.2   02/12/2009  Pathsat, Pathsun and related quantities are now
!                    layer arrays (Marco Matricardi).
!  1.3   05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.4   04/08/2010  Move addsolar check to calling routine (J Hocking)
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
       & rttov_coef,      &
       & profile_Type,    &
       & predictors_Type, &
       & raytracing_type
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
  USE rttov_const, ONLY : max_sol_zen
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type   ), INTENT(IN)    :: prof_tl(size(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_tl      ! in because of mem allocation
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_tl
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level, layer, iprof
! user profile
  REAL   (KIND=Jprb) :: t(prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w(prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2o  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4  (prof(1)%nlayers)
! reference profile
  REAL   (KIND=Jprb) :: tr   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tro  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: or   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2r (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2or (prof(1)%nlayers)
  REAL   (KIND=jprb) :: cor  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4r (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuw  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuwr (prof(1)%nlayers)
! user - reference
  REAL   (KIND=Jprb) :: dt   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dtabs(prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=Jprb) :: tw   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: co2w (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2ow (prof(1)%nlayers)
  REAL   (KIND=jprb) :: cow  (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4w (prof(1)%nlayers)
  REAL   (KIND=jprb) :: n2owr(prof(1)%nlayers)
  REAL   (KIND=jprb) :: cowr (prof(1)%nlayers)
  REAL   (KIND=jprb) :: ch4wr(prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=Jprb) :: sum1, sum2
  REAL   (KIND=Jprb) :: tr_sq      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_4       (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wr     (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_wrwr   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_or     (prof(1)%nlayers)
! TL variables
  REAL   (KIND=Jprb) :: t_tl       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: w_tl       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: o_tl       (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2o_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tr_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tro_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wr_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: wwr_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: or_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2r_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2or_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cor_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4r_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: twr_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuw_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tuwr_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dt_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: dto_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: tw_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ww_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ow_tl      (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: co2w_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2ow_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cow_tl     (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4w_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: n2owr_tl   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: cowr_tl    (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: ch4wr_tl   (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_or_tl  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wr_tl  (prof(1)%nlayers)
  REAL   (KIND=Jprb) :: sec_wrwr_tl(prof(1)%nlayers)
  INTEGER(KIND=jpim) :: nprofiles                   ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-------------------------------------------------------------------------------
! Recompute direct variables
!-------------------------------------------------------------------------------
! 1) profile layer mean quantities
!------------------------------------------------------------------------------
! layer N-1 lies between levels N-1 and N
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_SOLAR_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  DO iprof = 1, nprofiles
    IF (prof(iprof)%sunzenangle < 0.0 .OR. &
        prof(iprof)%sunzenangle >= max_sol_zen) CYCLE

    DO layer = 1, prof(1)%nlayers
      level    = layer + 1
!-Temperature
      t(layer) = (prof(iprof)%t(level - 1) + prof(iprof)%t(level)) / 2._JPRB
!-H2O
      w(layer) = (prof(iprof)%q(level - 1) + prof(iprof)%q(level)) * 0.5_JPRB
!-O3
      IF (opts%ozone_Data .AND. coef%nozone > 0)     &
        &  o(layer)   = (prof(iprof)%o3(level - 1) + prof(iprof)%o3(level)) * 0.5_JPRB
!-CO2
      IF (opts%co2_Data .AND. coef%nco2 > 0    )     &
        &  co2(layer) = (prof(iprof)%co2(level - 1) + prof(iprof)%co2(level)) * 0.5_JPRB
!-N2O
      IF (opts%n2o_Data .AND. coef%nn2o > 0    )     &
        &  n2o(layer) = (prof(iprof)%n2o(level - 1) + prof(iprof)%n2o(level)) * 0.5_JPRB
!-CO
      IF (opts%co_Data .AND. coef%nco > 0      )     &
        &  co(layer)  = (prof(iprof)%co(level - 1) + prof(iprof)%co(level)) * 0.5_JPRB
!-CH4
      IF (opts%ch4_Data .AND. coef%nch4 > 0    )     &
        &  ch4(layer) = (prof(iprof)%ch4(level - 1) + prof(iprof)%ch4(level)) * 0.5_JPRB
    ENDDO

! Direct value of o and co2 NOT needed for TL
!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt(:)    = t(:) - coef%tstar(:)
    dtabs(:) = Abs(dt(:))
!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
!-Temperature
    tr(:)    = t(:) / coef%tstar(:)
!-H2O
    wr(:)    = w(:) / coef%wstar(:)
!-CO2

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      co2r(:) = co2(:) / coef%co2star(:)
    ELSE
      co2r(:) = 1._JPRB
    ENDIF

!-Ozone

    IF (coef%nozone > 0) THEN
      dto(:) = t(:) - coef%to3star(:)
      tro(:) = t(:) / coef%to3star(:)
      IF (opts%ozone_Data) THEN
        or(:)  = o(:) / coef%ostar(:)
      ELSE
        or(:)  = 1._JPRB
      ENDIF
    ENDIF
    
!-N2O

    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or(:) = n2o(:) / coef%n2ostar(:)
    ELSE
      n2or(:) = 1._JPRB
    ENDIF

!-CO

    IF (opts%co_Data .AND. coef%nco > 0) THEN
      cor(:) = co(:) / coef%costar(:)
    ELSE
      cor(:) = 1._JPRB
    ENDIF

!-CH4

    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r(:) = ch4(:) / coef%ch4star(:)
    ELSE
      ch4r(:) = 1._JPRB
    ENDIF

!-------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw wwr twr
!--------------------------------------------------------------------
!-Temperature-------------------------------------------------------------------
    tw(1) = 0._jprb

    DO layer = 2, prof(1)%nlayers
      tw(layer) = tw(layer - 1) + coef%dpp(layer - 1) * tr(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1       = sum1 + t(layer)
      sum2       = sum2 + coef%tstar(layer)
      tuw(layer) = sum1 / sum2
    ENDDO

    tuwr(1)                 = coef%dpp(0) * t(1) / (coef%dpp(0) * coef%tstar(1))
    tuwr(2:prof(1)%nlayers) = tuw(2:prof(1)%nlayers)

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1       = sum1 + coef%dpp(layer - 1) * t(layer)
        sum2       = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        twr(layer) = sum1 / sum2
      ENDDO
    ENDIF

!-H2O----------------------------------------------------------------------------
    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1       = sum1 + coef%dpp(layer - 1) * w(layer) * t(layer)
      sum2       = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1      = sum1 + coef%dpp(layer - 1) * w(layer)
      sum2      = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      ww(layer) = sum1 / sum2
    ENDDO

!-O3-----------------------------------------------------------------------------

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1      = sum1 + coef%dpp(layer - 1) * o(layer)
        sum2      = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow(:) = 1._JPRB
    ENDIF

!-CO2---------------------------------------------------------------------------

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1        = sum1 + coef%dpp(layer - 1) * co2(layer)
        sum2        = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w(layer) = sum1 / sum2
      ENDDO

    ELSE
      co2w(:) = 1._JPRB
    ENDIF

!-N2O---------------------------------------------------------------------------

    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * n2o(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          n2ow(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * n2o(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        n2ow(:)  = 1._JPRB

        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr(layer) = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF
    
!-CO---------------------------------------------------------------------------

    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1       = sum1 + coef%dpp(layer - 1) * co(layer)
          sum2       = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          cow(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * co(layer) * t(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        cow(:)  = 1._JPRB
        
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * coef%costar(layer) * t(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr(layer) = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF
    
!-CH4---------------------------------------------------------------------------

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1        = sum1 + coef%dpp(layer - 1) * ch4(layer)
          sum2        = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          ch4w(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * ch4(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        ch4w(:)  = 1._JPRB
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1         = sum1 + coef%dpp(layer - 1) * coef%ch4star(layer) * t(layer)
          sum2         = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr(layer) = sum1 / sum2
        ENDDO
      
      ENDIF
    ENDIF
    
!-------------------------------------------------------------------------------
! Now compute TL variables
!-------------------------------------------------------------------------------
! 1) profile layer mean quantities
!------------------------------------------------------------------------------

    DO layer = 1, prof(1)%nlayers
      level       = layer + 1
      t_tl(layer) = (prof_tl(iprof)%t(level - 1) + prof_tl(iprof)%t(level)) * 0.5_JPRB
      w_tl(layer) = (prof_tl(iprof)%q(level - 1) + prof_tl(iprof)%q(level)) * 0.5_JPRB
      IF (opts%ozone_Data .AND. coef%nozone > 0)     &
        &  o_tl(layer)   = (prof_tl(iprof)%o3(level - 1) + prof_tl(iprof)%o3(level)) * 0.5_JPRB
      IF (opts%co2_Data .AND. coef%nco2 > 0    )     &
        &  co2_tl(layer) = (prof_tl(iprof)%co2(level - 1) + prof_tl(iprof)%co2(level)) * 0.5_JPRB
      IF (opts%n2o_Data .AND. coef%nn2o > 0    )     &
        &  n2o_tl(layer) = (prof_tl(iprof)%n2o(level - 1) + prof_tl(iprof)%n2o(level)) * 0.5_JPRB
      IF (opts%co_Data .AND. coef%nco > 0      )     &
        &  co_tl(layer)  = (prof_tl(iprof)%co(level - 1) + prof_tl(iprof)%co(level)) * 0.5_JPRB
      IF (opts%ch4_Data .AND. coef%nch4 > 0    )     &
        &  ch4_tl(layer) = (prof_tl(iprof)%ch4(level - 1) + prof_tl(iprof)%ch4(level)) * 0.5_JPRB
    ENDDO

!------------------------------------------------------------------------------
!2) calculate deviations from reference profile (layers)
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    dt_tl(:) = t_tl(:)

    IF (coef%nozone > 0) dto_tl(:) = t_tl(:)

!------------------------------------------------------------------------------
!3) calculate (profile / reference profile) ratios; tr_tl wr_tl or_tl
!------------------------------------------------------------------------------
! All assignments 1:prof(1) % nlayers
    tr_tl(:) = t_tl(:) / coef%tstar(:)
    wr_tl(:) = w_tl(:) / coef%wstar(:)

    IF (coef%nozone > 0) THEN
      tro_tl(:) = t_tl(:) / coef%to3star(:)
      IF (opts%ozone_Data) THEN
        or_tl(:)  = o_tl(:) / coef%ostar(:)
      ELSE
        or_tl(:)  = 0._JPRB
      ENDIF
    ENDIF
    

    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      co2r_tl(:) = co2_tl(:) / coef%co2star(:)
    ELSE
      co2r_tl(:) = 0._JPRB
    ENDIF


    IF (opts%n2o_Data .AND. coef%nn2o > 0) THEN
      n2or_tl(:) = n2o_tl(:) / coef%n2ostar(:)
    ELSE
      n2or_tl(:) = 0._JPRB
    ENDIF


    IF (opts%co_Data .AND. coef%nco > 0) THEN
      cor_tl(:) = co_tl(:) / coef%costar(:)
    ELSE
      cor_tl(:) = 0._JPRB
    ENDIF


    IF (opts%ch4_Data .AND. coef%nch4 > 0) THEN
      ch4r_tl(:) = ch4_tl(:) / coef%ch4star(:)
    ELSE
      ch4r_tl(:) = 0._JPRB
    ENDIF

!-------------------------------------------------------------------------
! 4. calculate profile / reference profile sums: tw_tl ww_tl ow_tl co2w_tl
!-------------------------------------------------------------------------
    tw_tl(1) = 0._JPRB

    DO layer = 2, prof(1)%nlayers
      tw_tl(layer) = tw_tl(layer - 1) + coef%dpp(layer - 1) * tr_tl(layer - 1)
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1          = sum1 + t_tl(layer)
      sum2          = sum2 + coef%tstar(layer)
      tuw_tl(layer) = sum1 / sum2
    ENDDO

    tuwr_tl(1)                 = coef%dpp(0) * t_tl(1) / (coef%dpp(0) * coef%tstar(1))
    tuwr_tl(2:prof(1)%nlayers) = tuw_tl(2:prof(1)%nlayers)

    IF (coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1          = sum1 + coef%dpp(layer - 1) * t_tl(layer)
        sum2          = sum2 + coef%dpp(layer - 1) * coef%tstar(layer)
        twr_tl(layer) = sum1 / sum2
      ENDDO
    ENDIF

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1         = sum1 + coef%dpp(layer - 1) * w_tl(layer)
      sum2         = sum2 + coef%dpp(layer - 1) * coef%wstar(layer)
      ww_tl(layer) = sum1 / sum2
    ENDDO

    sum1 = 0._JPRB
    sum2 = 0._JPRB

    DO layer = 1, prof(1)%nlayers
      sum1          = sum1 + coef%dpp(layer - 1) * (w_tl(layer) * t(layer) + w(layer) * t_tl(layer))
      sum2          = sum2 + coef%dpp(layer - 1) * coef%wstar(layer) * coef%tstar(layer)
      wwr_tl(layer) = sum1 / sum2
    ENDDO


    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1         = sum1 + coef%dpp(layer - 1) * o_tl(layer)
        sum2         = sum2 + coef%dpp(layer - 1) * coef%ostar(layer)
        ow_tl(layer) = sum1 / sum2
      ENDDO

    ELSE
      ow_tl(:) = 0._JPRB
    ENDIF


    IF (opts%co2_Data .AND. coef%nco2 > 0) THEN
      sum1 = 0._JPRB
      sum2 = 0._JPRB

      DO layer = 1, prof(1)%nlayers
        sum1           = sum1 + coef%dpp(layer - 1) * co2_tl(layer)
        sum2           = sum2 + coef%dpp(layer - 1) * coef%co2star(layer)
        co2w_tl(layer) = sum1 / sum2
      ENDDO

    ELSE
      co2w_tl(:) = 0._JPRB
    ENDIF


    IF (coef%nn2o > 0) THEN
      IF (opts%n2o_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1           = sum1 + coef%dpp(layer - 1) * n2o_tl(layer)
          sum2           = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer)
          n2ow_tl(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1            = sum1 + coef%dpp(layer - 1) * (n2o_tl(layer) * t(layer) + n2o(layer) * t_tl(layer))
          sum2            = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr_tl(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        n2ow_tl(:)  = 0._JPRB

        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1            = sum1 + coef%dpp(layer - 1) * coef%n2ostar(layer) * t_tl(layer)
          sum2            = sum2 + coef%dpp(layer - 1) * coef%n2ostar(layer) * coef%tstar(layer)
          n2owr_tl(layer) = sum1 / sum2
        ENDDO

      ENDIF
    ENDIF
    

    IF (coef%nco > 0) THEN
      IF (opts%co_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1          = sum1 + coef%dpp(layer - 1) * co_tl(layer)
          sum2          = sum2 + coef%dpp(layer - 1) * coef%costar(layer)
          cow_tl(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1           = sum1 + coef%dpp(layer - 1) * (co_tl(layer) * t(layer) + co(layer) * t_tl(layer))
          sum2           = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr_tl(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        cow_tl(:)  = 0._JPRB
        
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1           = sum1 + coef%dpp(layer - 1) * coef%costar(layer) * t_tl(layer)
          sum2           = sum2 + coef%dpp(layer - 1) * coef%costar(layer) * coef%tstar(layer)
          cowr_tl(layer) = sum1 / sum2
        ENDDO
        
      ENDIF
    ENDIF
    

    IF (coef%nch4 > 0) THEN
      IF (opts%ch4_Data) THEN
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1           = sum1 + coef%dpp(layer - 1) * ch4_tl(layer)
          sum2           = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer)
          ch4w_tl(layer) = sum1 / sum2
        ENDDO
  
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1            = sum1 + coef%dpp(layer - 1) * (ch4_tl(layer) * t(layer) + ch4(layer) * t_tl(layer))
          sum2            = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr_tl(layer) = sum1 / sum2
        ENDDO
  
      ELSE
        ch4w_tl(:)  = 0._JPRB
        
        sum1 = 0._JPRB
        sum2 = 0._JPRB
  
        DO layer = 1, prof(1)%nlayers
          sum1            = sum1 + coef%dpp(layer - 1) * coef%ch4star(layer) * t_tl(layer)
          sum2            = sum2 + coef%dpp(layer - 1) * coef%ch4star(layer) * coef%tstar(layer)
          ch4wr_tl(layer) = sum1 / sum2
        ENDDO
        
      ENDIF
    ENDIF

! End of TL profile calcs
! ATTENTION
!  w_tl(:) = prof_tl(iprof) % q(:)
!5) set predictors for RTTOV-8 options
!--
    raytracing_tl%patheff(:, iprof) = raytracing_tl%pathsat(:, iprof) + raytracing_tl%pathsun(:, iprof)
!5.1 mixed gases
!---

    DO layer = 1, prof(1)%nlayers
      level = layer + 1
      predictors_tl%mixedgas_sun(1, layer, iprof)     = raytracing_tl%patheff(layer, iprof)
      predictors_tl%mixedgas_sun(2, layer, iprof)     =      &
        & raytracing_tl%patheff(layer, iprof) * 2 * raytracing%patheff(layer, iprof)
      predictors_tl%mixedgas_sun(3, layer, iprof)     =      &
        & tr_tl(layer) * raytracing%patheff(layer, iprof) + raytracing_tl%patheff(layer, iprof) * tr(layer)
      predictors_tl%mixedgas_sun(4, layer, iprof)     =                        &
        & 2._JPRB * tr_tl(layer) * predictors%mixedgas_sun(3, layer, iprof) +  &
        & raytracing_tl%patheff(layer, iprof) * predictors%mixedgas_sun(5, layer, iprof) ** 2
      predictors_tl%mixedgas_sun(5, layer, iprof)     = tr_tl(layer)
      predictors_tl%mixedgas_sun(6, layer, iprof)     = 2._JPRB * tr_tl(layer) * tr(layer)
      predictors_tl%mixedgas_sun(7, layer, iprof)     =      &
        & raytracing%patheff(layer, iprof) * tuw_tl(layer) + raytracing_tl%patheff(layer, iprof) * tuw(layer)
      predictors_tl%mixedgas_sun(8, layer, iprof)     =      &
        & raytracing%patheff(layer, iprof) * tuwr_tl(layer) + raytracing_tl%patheff(layer, iprof) * tuwr(layer)
      predictors_tl%mixedgas_sun(9, layer, iprof)     = raytracing_tl%patheff(layer, iprof) * tr(layer) ** 3 +      &
        & tr_tl(layer) * 3 * raytracing%patheff(layer, iprof) * tr(layer) ** 2
      predictors_tl%mixedgas_sun(10, layer, iprof)    =                                                 &
        & raytracing_tl%patheff(layer, iprof) * 1.5 * sqrt(predictors%mixedgas_sun(3, layer, iprof)) +  &
        & tr_tl(layer) * 0.5 * raytracing%patheff(layer, iprof) ** 1.5 / sqrt(tr(layer))
!5.2 water vapour lines based on RTIASI
!--------------------------------------
      sec_wr = raytracing%patheff(layer, iprof) * wr(layer)
      sec_wrwr = sec_wr(layer) * wr(layer)
      sec_wr_tl(layer)                                =      &
        & raytracing%patheff(layer, iprof) * wr_tl(layer) + raytracing_tl%patheff(layer, iprof) * wr(layer)
      sec_wrwr_tl(layer)                              =      &
        & sec_wr_tl(layer) * wr(layer) + predictors%watervapour_sun(7, layer, iprof) * wr_tl(layer)
      predictors_tl%watervapour_sun(:, layer, iprof)  = 0._JPRB
      predictors_tl%watervapour_sun(1, layer, iprof)  =      &
        & 2 * predictors%watervapour_sun(7, layer, iprof) * sec_wr_tl(layer)
      predictors_tl%watervapour_sun(2, layer, iprof)  =      &
        & raytracing%patheff(layer, iprof) * ww_tl(layer) + raytracing_tl%patheff(layer, iprof) * ww(layer)
      predictors_tl%watervapour_sun(3, layer, iprof)  =                                                        &
        & 2 * predictors%watervapour_sun(2, layer, iprof) * raytracing%patheff(layer, iprof) * ww_tl(layer) +  &
        & 2 * predictors%watervapour_sun(2, layer, iprof) * raytracing_tl%patheff(layer, iprof) * ww(layer)
      predictors_tl%watervapour_sun(4, layer, iprof)  =      &
        & predictors%watervapour_sun(7, layer, iprof) * dt_tl(layer) + sec_wr_tl(layer) * dt(layer)
      predictors_tl%watervapour_sun(5, layer, iprof)  =      &
        & 0.5_JPRB * sec_wr_tl(layer) / predictors%watervapour_sun(5, layer, iprof)
      predictors_tl%watervapour_sun(6, layer, iprof)  =      &
        & 0.25_JPRB * sec_wr_tl(layer) / predictors%watervapour_sun(6, layer, iprof) ** 3
      predictors_tl%watervapour_sun(7, layer, iprof)  = sec_wr_tl(layer)
      predictors_tl%watervapour_sun(8, layer, iprof)  =                                                                &
        & raytracing_tl%patheff(layer, iprof) * 1.5 * SQRT(predictors%watervapour_sun(2, layer, iprof)) * ww(layer) +  &
        & ww_tl(layer) * 1.5 * SQRT(predictors%watervapour_sun(2, layer, iprof)) * raytracing%patheff(layer, iprof)
      predictors_tl%watervapour_sun(9, layer, iprof)  =                                                                &
        & raytracing_tl%patheff(layer, iprof) * 1.5 * SQRT(predictors%watervapour_sun(7, layer, iprof)) * wr(layer) +  &
        & wr_tl(layer) * 1.5 * SQRT(predictors%watervapour_sun(7, layer, iprof)) * raytracing%patheff(layer, iprof)
      predictors_tl%watervapour_sun(10, layer, iprof) =      &
        & dtabs(layer) * (sec_wr_tl(layer) * dt(layer) + 2 * predictors%watervapour_sun(7, layer, iprof) * dt_tl(layer))
      predictors_tl%watervapour_sun(11, layer, iprof) = predictors%watervapour_sun(5, layer, iprof) * dt_tl(layer) +      &
        & 0.5_JPRB * dt(layer) * sec_wr_tl(layer) / predictors%watervapour_sun(5, layer, iprof)
      predictors_tl%watervapour_sun(12, layer, iprof) =                                                                     &
        & raytracing_tl%patheff(layer, iprof) * 1.25 * (predictors%watervapour_sun(2, layer, iprof)) ** 0.25 * ww(layer) +  &
        & ww_tl(layer) * 1.25 * (predictors%watervapour_sun(2, layer, iprof)) ** 0.25 * raytracing%patheff(layer, iprof)
      predictors_tl%watervapour_sun(13, layer, iprof) = 1.5 * predictors%watervapour_sun(5, layer, iprof) * wr_tl(layer)     &
        &  + raytracing_tl%patheff(layer, iprof) * 0.5 * wr(layer) ** 1.5 / raytracing%patheff(layer, iprof) ** 0.5
      predictors_tl%watervapour_sun(14, layer, iprof) =                                                                      &
        & dt_tl(layer) * predictors%watervapour_sun(7, layer, iprof) ** 1.5 +                                                &
        & raytracing_tl%patheff(layer, iprof) * 1.5 * predictors%watervapour_sun(5, layer, iprof) * wr(layer) * dt(layer) +  &
        & wr_tl(layer) * 1.5 * predictors%watervapour_sun(5, layer, iprof) * raytracing%patheff(layer, iprof) * dt(layer)
      predictors_tl%watervapour_sun(15, layer, iprof) =                                                                    &
        & 2 * predictors%watervapour_sun(7, layer, iprof) * wr_tl(layer) / ww(layer) -                                     &
        & predictors%watervapour_sun(1, layer, iprof) * ww_tl(layer) / (raytracing%patheff(layer, iprof) * ww(layer) ** 2) &
        &  + raytracing_tl%patheff(layer, iprof) * wr(layer) ** 2 / ww(layer)
!
!5.3 water vapour continuum transmittance based on RTIASI
!--------------------------------------------------------
!

      IF (coef%nwvcont > 0) THEN
        tr_sq(layer)                              = tr(layer) * tr(layer)
        tr_4(layer)                               = tr_sq(layer) * tr_sq(layer)
        predictors_tl%wvcont_sun(1, layer, iprof) =      &
          & sec_wrwr_tl(layer) / tr(layer) - predictors%wvcont_sun(1, layer, iprof) * tr_tl(layer) / tr(layer)
        predictors_tl%wvcont_sun(2, layer, iprof) =      &
          & sec_wrwr_tl(layer) / tr_4(layer) - 4 * predictors%wvcont_sun(1, layer, iprof) * tr_tl(layer) / tr_4(layer)
        predictors_tl%wvcont_sun(3, layer, iprof) =      &
          & sec_wr_tl(layer) / tr(layer) - predictors%watervapour_sun(7, layer, iprof) * tr_tl(layer) / tr_sq(layer)
        predictors_tl%wvcont_sun(4, layer, iprof) = sec_wr_tl(layer) / tr_sq(layer) -      &
          & 2 * predictors%watervapour_sun(7, layer, iprof) * tr_tl(layer) / (tr_sq(layer) * tr(layer))
      ENDIF

!
!5.4 ozone
!---------

      IF (coef%nozone > 0) THEN
        sec_or(layer)                             = raytracing%patheff(layer, iprof) * or(layer)
        sec_or_tl(layer)                          =      &
          & raytracing_tl%patheff(layer, iprof) * or(layer) + or_tl(layer) * raytracing%patheff(layer, iprof)
        predictors_tl%ozone_sun(1, layer, iprof)  = sec_or_tl(layer)
        predictors_tl%ozone_sun(2, layer, iprof)  = 0.5_JPRB * sec_or_tl(layer) / predictors%ozone_sun(2, layer, iprof)
        predictors_tl%ozone_sun(3, layer, iprof)  = raytracing_tl%patheff(layer, iprof) * or(layer) / ow(layer) +      &
          & or_tl(layer) * raytracing%patheff(layer, iprof) / ow(layer) -                                              &
          & or(layer) * raytracing%patheff(layer, iprof) * ow_tl(layer) / ow(layer) ** 2
        predictors_tl%ozone_sun(4, layer, iprof)  = 2 * sec_or_tl(layer) * predictors%ozone_sun(1, layer, iprof)
        predictors_tl%ozone_sun(5, layer, iprof)  =                                               &
          & 0.5_JPRB * sec_or_tl(layer) * dto(layer) / (predictors%ozone_sun(2, layer, iprof)) +  &
          & predictors%ozone_sun(2, layer, iprof) * dto_tl(layer)
        predictors_tl%ozone_sun(6, layer, iprof)  =                                                &
          & sec_or_tl(layer) * or(layer) * ow(layer) + or_tl(layer) * sec_or(layer) * ow(layer) +  &
          & ow_tl(layer) * sec_or(layer) * or(layer)
        predictors_tl%ozone_sun(7, layer, iprof)  = predictors_tl%ozone_sun(2, layer, iprof) * or(layer) / ow(layer) +      &
          & or_tl(layer) * predictors%ozone_sun(2, layer, iprof) / ow(layer) -                                              &
          & ow_tl(layer) * predictors%ozone_sun(2, layer, iprof) * or(layer) / ow(layer) ** 2
        predictors_tl%ozone_sun(8, layer, iprof)  =                                                                         &
          & ow_tl(layer) * 1.75 * raytracing%patheff(layer, iprof) * (raytracing%patheff(layer, iprof) * ow(layer)) ** 0.75 &
          &  +                                                                                                              &
          & raytracing_tl%patheff(layer, iprof) * 1.75 * ow(layer) * (raytracing%patheff(layer, iprof) * ow(layer)) ** 0.75
        predictors_tl%ozone_sun(9, layer, iprof)  =                                                           &
          & raytracing%patheff(layer, iprof) * or_tl(layer) * Sqrt(predictors%ozone_sun(10, layer, iprof)) +  &
          & predictors%ozone_sun(1, layer, iprof) * 0.5 * raytracing%patheff(layer, iprof) * ow_tl(layer) /   &
          & Sqrt(predictors%ozone_sun(10, layer, iprof)) +                                                    &
          & 1.5 * raytracing_tl%patheff(layer, iprof) * or(layer) * Sqrt(predictors%ozone_sun(10, layer, iprof))
        predictors_tl%ozone_sun(10, layer, iprof) =      &
          & raytracing_tl%patheff(layer, iprof) * ow(layer) + raytracing%patheff(layer, iprof) * ow_tl(layer)
        predictors_tl%ozone_sun(11, layer, iprof) =                                                         &
          & 2 * ow_tl(layer) * raytracing%patheff(layer, iprof) * predictors%ozone_sun(10, layer, iprof) +  &
          & 2 * raytracing_tl%patheff(layer, iprof) * ow(layer) * predictors%ozone_sun(10, layer, iprof)
        predictors_tl%ozone_sun(12, layer, iprof) =                                                                      &
          & raytracing_tl%patheff(layer, iprof) * 0.5 * raytracing%patheff(layer, iprof) ** ( - 0.5) * ow(layer) ** 2 *  &
          & dto(layer) + ow_tl(layer) * 2 * ow(layer) * raytracing%patheff(layer, iprof) ** 0.5 * dto(layer) +           &
          & dto_tl(layer) * raytracing%patheff(layer, iprof) ** 0.5 * ow(layer) ** 2
        predictors_tl%ozone_sun(13, layer, iprof) = raytracing_tl%patheff(layer, iprof) * tro(layer) ** 3 +      &
          & tro_tl(layer) * 3 * tro(layer) ** 2 * raytracing%patheff(layer, iprof)
      ENDIF

!
!5.6 carbon dioxide transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nco2 > 0) THEN
        predictors_tl%co2_sun(1, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * co2r_tl(layer) + raytracing_tl%patheff(layer, iprof) * co2r(layer)
        predictors_tl%co2_sun(2, layer, iprof)  = 2 * tr_tl(layer) * predictors%co2_sun(5, layer, iprof)
        predictors_tl%co2_sun(3, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * tr_tl(layer) + raytracing_tl%patheff(layer, iprof) * tr(layer)
        predictors_tl%co2_sun(4, layer, iprof)  =      &
          & 2 * tr_tl(layer) * predictors%co2_sun(3, layer, iprof) + raytracing_tl%patheff(layer, iprof) * tr(layer) ** 2
        predictors_tl%co2_sun(5, layer, iprof)  = tr_tl(layer)
        predictors_tl%co2_sun(6, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * twr_tl(layer) + raytracing_tl%patheff(layer, iprof) * twr(layer)
        predictors_tl%co2_sun(7, layer, iprof)  =                                                                &
          & 2 * raytracing%patheff(layer, iprof) * SQRT(predictors%co2_sun(7, layer, iprof)) * co2w_tl(layer) +  &
          & 2 * co2w(layer) * SQRT(predictors%co2_sun(7, layer, iprof)) * raytracing_tl%patheff(layer, iprof)
        predictors_tl%co2_sun(8, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * cor_tl(layer) + cor(layer) * raytracing_tl%patheff(layer, iprof)
        predictors_tl%co2_sun(9, layer, iprof)  =                                                                   &
          & raytracing%patheff(layer, iprof) * SQRT(predictors%co2_sun(5, layer, iprof)) * twr_tl(layer) +          &
          & 0.5 * predictors%co2_sun(6, layer, iprof) * tr_tl(layer) / SQRT(predictors%co2_sun(5, layer, iprof)) +  &
          & raytracing_tl%patheff(layer, iprof) * twr(layer) * tr(layer) ** 0.5
        predictors_tl%co2_sun(10, layer, iprof) =                                                                  &
          & raytracing_tl%patheff(layer, iprof) * 0.5 * co2r(layer) / SQRT(predictors%co2_sun(1, layer, iprof)) +  &
          & co2r_tl(layer) * 0.5 * raytracing%patheff(layer, iprof) / SQRT(predictors%co2_sun(1, layer, iprof))
        predictors_tl%co2_sun(11, layer, iprof) = tr_tl(layer) * 3 * tr(layer) ** 2
        predictors_tl%co2_sun(12, layer, iprof) = raytracing_tl%patheff(layer, iprof) * tr(layer) ** 3 +      &
          & tr_tl(layer) * 3 * tr(layer) ** 2 * raytracing%patheff(layer, iprof)
        predictors_tl%co2_sun(13, layer, iprof) =                                                      &
          & raytracing_tl%patheff(layer, iprof) * tr(layer) ** 2 * twr(layer) ** 3 * 0.5 /             &
          & SQRT(raytracing%patheff(layer, iprof)) +                                                   &
          & tr_tl(layer) * 2 * tr(layer) * twr(layer) ** 3 * SQRT(raytracing%patheff(layer, iprof)) +  &
          & twr_tl(layer) * 3 * twr(layer) ** 2 * tr(layer) ** 2 * raytracing%patheff(layer, iprof) ** 0.5
        predictors_tl%co2_sun(14, layer, iprof) =                                         &
          & tr_tl(layer) * 2 * SQRT(predictors%co2_sun(14, layer, iprof)) * twr(layer) +  &
          & twr_tl(layer) * 2 * SQRT(predictors%co2_sun(14, layer, iprof)) * tr(layer)
      ENDIF

!5.7 n2o            transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nn2o > 0) THEN
        predictors_tl%n2o_sun(1, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * n2or_tl(layer) + raytracing_tl%patheff(layer, iprof) * n2or(layer)
        predictors_tl%n2o_sun(2, layer, iprof)  =                                                            &
          & 0.5 * raytracing%patheff(layer, iprof) * n2or_tl(layer) / predictors%n2o_sun(2, layer, iprof) +  &
          & 0.5 * raytracing_tl%patheff(layer, iprof) * n2or(layer) / predictors%n2o_sun(2, layer, iprof)
        predictors_tl%n2o_sun(3, layer, iprof)  = predictors%n2o_sun(1, layer, iprof) * dt_tl(layer) +      &
          & raytracing%patheff(layer, iprof) * n2or_tl(layer) * dt(layer) +                                 &
          & raytracing_tl%patheff(layer, iprof) * n2or(layer) * dt(layer)
        predictors_tl%n2o_sun(4, layer, iprof)  =                                                          &
          & 2 * predictors%n2o_sun(1, layer, iprof) * n2or_tl(layer) * raytracing%patheff(layer, iprof) +  &
          & 2 * predictors%n2o_sun(1, layer, iprof) * n2or(layer) * raytracing_tl%patheff(layer, iprof)
        predictors_tl%n2o_sun(5, layer, iprof)  = n2or(layer) * dt_tl(layer) + n2or_tl(layer) * dt(layer)
        predictors_tl%n2o_sun(6, layer, iprof)  =                                                                  &
          & raytracing%patheff(layer, iprof) * n2or_tl(layer) * 0.25 / predictors%n2o_sun(6, layer, iprof) ** 3 +  &
          & raytracing_tl%patheff(layer, iprof) * n2or(layer) * 0.25 / predictors%n2o_sun(6, layer, iprof) ** 3
        predictors_tl%n2o_sun(7, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * n2ow_tl(layer) + raytracing_tl%patheff(layer, iprof) * n2ow(layer)
        predictors_tl%n2o_sun(8, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * n2owr_tl(layer) + raytracing_tl%patheff(layer, iprof) * n2owr(layer)
        predictors_tl%n2o_sun(9, layer, iprof)  =                                                           &
          & raytracing_tl%patheff(layer, iprof) * 2 * predictors%n2o_sun(8, layer, iprof) * n2owr(layer) +  &
          & n2owr_tl(layer) * 2 * predictors%n2o_sun(8, layer, iprof) * raytracing%patheff(layer, iprof)
        predictors_tl%n2o_sun(10, layer, iprof) =                                            &
          & 1.5_JPRB * predictors%n2o_sun(2, layer, iprof) * n2or_tl(layer) / n2ow(layer) -  &
          & predictors%n2o_sun(2, layer, iprof) ** 3 * n2ow_tl(layer) /                      &
          & (raytracing%patheff(layer, iprof) * n2ow(layer) ** 2) +                          &
          & 0.5_JPRB * n2or(layer) ** 1.5_JPRB * raytracing_tl%patheff(layer, iprof) /       &
          & (raytracing%patheff(layer, iprof) ** 0.5_JPRB * n2ow(layer))
        predictors_tl%n2o_sun(11, layer, iprof) =                                                                &
          & raytracing_tl%patheff(layer, iprof) * 3 * predictors%n2o_sun(8, layer, iprof) ** 2 * n2owr(layer) +  &
          & n2owr_tl(layer) * 3 * predictors%n2o_sun(8, layer, iprof) ** 2 * raytracing%patheff(layer, iprof)
        predictors_tl%n2o_sun(12, layer, iprof) =                                                                    &
          & raytracing_tl%patheff(layer, iprof) * 2 * raytracing%patheff(layer, iprof) * n2owr(layer) * dt(layer) +  &
          & n2owr_tl(layer) * raytracing%patheff(layer, iprof) ** 2 * dt(layer) +                                    &
          & dt_tl(layer) * raytracing%patheff(layer, iprof) ** 2 * n2owr(layer)
      ENDIF

!5.8 co             transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nco > 0) THEN
        predictors_tl%co_sun(1, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * cor_tl(layer) + raytracing_tl%patheff(layer, iprof) * cor(layer)
        predictors_tl%co_sun(2, layer, iprof)  =                                                           &
          & 0.5 * raytracing%patheff(layer, iprof) * cor_tl(layer) / predictors%co_sun(2, layer, iprof) +  &
          & 0.5 * raytracing_tl%patheff(layer, iprof) * cor(layer) / predictors%co_sun(2, layer, iprof)
        predictors_tl%co_sun(3, layer, iprof)  =                                                                             &
          & predictors%co_sun(1, layer, iprof) * dt_tl(layer) + dt(layer) * raytracing%patheff(layer, iprof) * cor_tl(layer) &
          &  + raytracing_tl%patheff(layer, iprof) * cor(layer) * dt(layer)
        predictors_tl%co_sun(4, layer, iprof)  =                                                         &
          & 2 * predictors%co_sun(1, layer, iprof) * cor_tl(layer) * raytracing%patheff(layer, iprof) +  &
          & 2 * predictors%co_sun(1, layer, iprof) * cor(layer) * raytracing_tl%patheff(layer, iprof)
        predictors_tl%co_sun(5, layer, iprof)  = predictors%co_sun(2, layer, iprof) * dt_tl(layer) +                   &
          & 0.5 * dt(layer) * raytracing%patheff(layer, iprof) * cor_tl(layer) / predictors%co_sun(2, layer, iprof) +  &
          & 0.5 * dt(layer) * raytracing_tl%patheff(layer, iprof) * cor(layer) / predictors%co_sun(2, layer, iprof)
        predictors_tl%co_sun(6, layer, iprof)  =                                                                 &
          & raytracing%patheff(layer, iprof) * cor_tl(layer) * 0.25 / predictors%co_sun(6, layer, iprof) ** 3 +  &
          & raytracing_tl%patheff(layer, iprof) * cor(layer) * 0.25 / predictors%co_sun(6, layer, iprof) ** 3
        predictors_tl%co_sun(7, layer, iprof)  = ABS(dt(layer)) * (         &
          & raytracing%patheff(layer, iprof) * cor_tl(layer) * dt(layer) +  &
          & predictors%co_sun(1, layer, iprof) * 2 * dt_tl(layer) +         &
          & raytracing_tl%patheff(layer, iprof) * cor(layer) * dt(layer))
        predictors_tl%co_sun(8, layer, iprof)  = (2 * predictors%co_sun(1, layer, iprof) * cor_tl(layer) / cow(layer))     &
          &  - predictors%co_sun(8, layer, iprof) * cow_tl(layer) / cow(layer) +                                           &
          & raytracing_tl%patheff(layer, iprof) * cor(layer) ** 2 / cow(layer)
        predictors_tl%co_sun(9, layer, iprof)  = 1.5 * predictors%co_sun(2, layer, iprof) * cor_tl(layer) / cow(layer)     &
          &  -                                                                                                             &
          & predictors%co_sun(2, layer, iprof) ** 3 * cow_tl(layer) / (raytracing%patheff(layer, iprof) * cow(layer) ** 2) &
          &  + 0.5 * cor(layer) ** 1.5 * raytracing_tl%patheff(layer, iprof) /                                             &
          & (raytracing%patheff(layer, iprof) ** 0.5 * cow(layer))
        predictors_tl%co_sun(10, layer, iprof) =                                            &
          & (2 * predictors%co_sun(1, layer, iprof) * cor_tl(layer) / cow(layer) ** 0.5) -  &
          & 0.5 * predictors%co_sun(10, layer, iprof) * cow_tl(layer) / cow(layer) +        &
          & raytracing_tl%patheff(layer, iprof) * cor(layer) ** 2 / SQRT(cow(layer))
        predictors_tl%co_sun(11, layer, iprof) = raytracing%patheff(layer, iprof) * 0.4 * cowr_tl(layer) *      &
          & (raytracing%patheff(layer, iprof) * cowr(layer)) ** ( - 0.6) +                                      &
          & cowr(layer) * 0.4 * raytracing_tl%patheff(layer, iprof) *                                           &
          & (raytracing%patheff(layer, iprof) * cowr(layer)) ** ( - 0.6)
        predictors_tl%co_sun(12, layer, iprof) = raytracing%patheff(layer, iprof) * 0.25 * cowr_tl(layer) *      &
          & (raytracing%patheff(layer, iprof) * cowr(layer)) ** ( - 0.75) +                                      &
          & cowr(layer) * 0.25 * raytracing_tl%patheff(layer, iprof) *                                           &
          & (raytracing%patheff(layer, iprof) * cowr(layer)) ** ( - 0.75)
      ENDIF

!5.9 ch4            transmittance based on RTIASI
!-------------------------------------------------
!

      IF (coef%nch4 > 0) THEN
        predictors_tl%ch4_sun(1, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * ch4r_tl(layer) + raytracing_tl%patheff(layer, iprof) * ch4r(layer)
        predictors_tl%ch4_sun(2, layer, iprof)  =                                                            &
          & 0.5 * raytracing%patheff(layer, iprof) * ch4r_tl(layer) / predictors%ch4_sun(2, layer, iprof) +  &
          & 0.5 * raytracing_tl%patheff(layer, iprof) * ch4r(layer) / predictors%ch4_sun(2, layer, iprof)
        predictors_tl%ch4_sun(3, layer, iprof)  = predictors%ch4_sun(1, layer, iprof) * dt_tl(layer) +      &
          & raytracing%patheff(layer, iprof) * ch4r_tl(layer) * dt(layer) +                                 &
          & raytracing_tl%patheff(layer, iprof) * ch4r(layer) * dt(layer)
        predictors_tl%ch4_sun(4, layer, iprof)  =                                                          &
          & 2 * predictors%ch4_sun(1, layer, iprof) * ch4r_tl(layer) * raytracing%patheff(layer, iprof) +  &
          & 2 * predictors%ch4_sun(1, layer, iprof) * ch4r(layer) * raytracing_tl%patheff(layer, iprof)
        predictors_tl%ch4_sun(5, layer, iprof)  = ch4r(layer) * dt_tl(layer) + ch4r_tl(layer) * dt(layer)
        predictors_tl%ch4_sun(6, layer, iprof)  =                                                                  &
          & raytracing%patheff(layer, iprof) * ch4r_tl(layer) * 0.25 / predictors%ch4_sun(6, layer, iprof) ** 3 +  &
          & raytracing_tl%patheff(layer, iprof) * ch4r(layer) * 0.25 / predictors%ch4_sun(6, layer, iprof) ** 3
        predictors_tl%ch4_sun(7, layer, iprof)  =      &
          & raytracing%patheff(layer, iprof) * ch4wr_tl(layer) + raytracing_tl%patheff(layer, iprof) * ch4wr(layer)
        predictors_tl%ch4_sun(8, layer, iprof)  = ch4wr_tl(layer)
        predictors_tl%ch4_sun(9, layer, iprof)  =                                                           &
          & 2 * predictors%ch4_sun(10, layer, iprof) * ch4w_tl(layer) * raytracing%patheff(layer, iprof) +  &
          & 2 * predictors%ch4_sun(10, layer, iprof) * ch4w(layer) * raytracing_tl%patheff(layer, iprof)
        predictors_tl%ch4_sun(10, layer, iprof) =      &
          & raytracing%patheff(layer, iprof) * ch4w_tl(layer) + raytracing_tl%patheff(layer, iprof) * ch4w(layer)
        predictors_tl%ch4_sun(11, layer, iprof) =                                       &
          & 1.5 * predictors%ch4_sun(2, layer, iprof) * ch4r_tl(layer) / ch4w(layer) -  &
          & predictors%ch4_sun(2, layer, iprof) ** 3 * ch4w_tl(layer) /                 &
          & (raytracing%patheff(layer, iprof) * ch4w(layer) ** 2) +                     &
          & 0.5 * ch4r(layer) ** 1.5 * raytracing_tl%patheff(layer, iprof) /            &
          & (raytracing%patheff(layer, iprof) ** 0.5 * ch4w(layer))
      ENDIF

    ENDDO
! layers
  ENDDO
! profiles
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_9_SOLAR_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_9_solar_tl
