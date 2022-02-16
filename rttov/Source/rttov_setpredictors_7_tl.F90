!
SUBROUTINE rttov_setpredictors_7_tl( &
            & opts,          &
            & prof,          &
            & prof_tl,       &
            & geom,          &
            & coef,          &
            & aux,           &
            & predictors,    &
            & predictors_tl, &
            & raytracing,    &
            & raytracing_tl)
!
! Description
! RTTOV-7 Model
! TL of rttov_setpredictors
! To calculate and store the profile variables (predictors) required
! in subsequent transmittance calculations.
! Code based on PRFTAU from previous versions of RTTOV
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
! see RTTOV7 science and validation report pages 18/19
! variable names are close to the documentation
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.2       29/03/2005  Add end of header comment (J. Cameron)
!  1.3       08/12/2005  Add surface humidity (R Saunders)
!  1.4       03/03/2006  Marco Matricardi (ECMWF):
!               --       Altitude dependent local zenith angle introduced.
!  1.5       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.6       15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.7       02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                        (Marco Matricardi).
!  1.8       17/06/2010  Combined non-Zeeman and Zeeman predictors for SSMIS for
!                        use with single coefficient file  (P Rayer)                       (Marco Matricardi).
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
       & rttov_coef,      &
       & rttov_options,   &
       & profile_Type,    &
       & geometry_Type,   &
       & profile_aux,     &
       & predictors_Type, &
       & raytracing_type
!INTF_OFF
  USE rttov_const, ONLY : gravity
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type   ), INTENT(IN)    :: prof_tl(size(prof))
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(INOUT) :: predictors_tl      ! in because of mem allocation
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(IN)    :: raytracing_tl
  TYPE(profile_aux    ), INTENT(IN)    :: aux                ! auxillary profiles info.
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer
  INTEGER(KIND=jpim) :: iv2lev, iv3lev, iv2lay, i
! user profile
  REAL   (KIND=jprb) :: t(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: w(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(size(prof)     , prof(1)%nlayers)
! reference profile
  REAL   (KIND=jprb) :: tr(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: wr(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: or(size(prof)     , prof(1)%nlayers)
! user - reference
  REAL   (KIND=jprb) :: dt(size(prof)     , prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=jprb) :: tw(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow(size(prof)     , prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=jprb) :: sum1     (size(prof)                      ), sum2(size(prof))
  REAL   (KIND=jprb) :: deltac   (prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_wr   (prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_or   (prof(1)%nlayers                 )
! TL variables
  REAL   (KIND=jprb) :: t_tl     (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: w_tl     (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: o_tl     (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: tr_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: wr_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: or_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: dt_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: dto_tl   (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: tw_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow_tl    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: sec_or_tl(prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_wr_tl(prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: zsqrt       , zrecip
  INTEGER(KIND=jpim) :: nprofiles                                                    ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
! profile layer quantities
! Direct variables
!CDIR NOLOOPCHG
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles
      t(i, layer) = (prof(i)%t(level - 1) + prof(i)%t(level)) * 0.5_JPRB
      w(i, layer) = (prof(i)%q(level - 1) + prof(i)%q(level)) * 0.5_JPRB
    ENDDO

  ENDDO


  IF (opts%use_q2m) THEN

    DO i = 1, nprofiles
! include surface humidity
      iv3lev = aux%s(i)%nearestlev_surf - 1
      iv2lev = aux%s(i)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        w(i, iv2lay) = (prof(i)%s2m%q + prof(i)%q(iv3lev)) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDIF

!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        o(i, layer) = (prof(i)%o3(level - 1) + prof(i)%o3(level)) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDDO


! TL variables
!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles
      t_tl(i, layer) = (prof_tl(i)%t(level - 1) + prof_tl(i)%t(level)) * 0.5_JPRB
      w_tl(i, layer) = (prof_tl(i)%q(level - 1) + prof_tl(i)%q(level)) * 0.5_JPRB
    ENDDO

  ENDDO


!
! include tl surface humidity

  IF (opts%use_q2m) THEN

    DO i = 1, nprofiles
      iv3lev = aux%s(i)%nearestlev_surf - 1
      iv2lev = aux%s(i)%nearestlev_surf

      IF (iv2lev <= coef%nlevels) THEN
        w_tl(i, iv2lay) = (prof_tl(i)%s2m%q + prof_tl(i)%q(iv3lev)) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDIF

!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles
      IF (opts%ozone_Data .AND. coef%nozone > 0)     &
        &  o_tl(i, layer) = (prof_tl(i)%o3(level - 1) + prof_tl(i)%o3(level)) * 0.5_JPRB
    ENDDO

  ENDDO


! calculate deviations from reference profile (layers)
! if no input O3 profile,  we still use the input temperature profile for dto
! Direct variables & TL variables
!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers

    DO i = 1, nprofiles
      dt(i, layer)    = t(i, layer) - coef%tstar(layer)
      dt_tl(i, layer) = t_tl(i, layer)
! Direct value of dto NOT needed for TL

      IF (coef%nozone > 0) dto_tl(i, layer) = t_tl(i, layer)

! calculate (profile / reference profile) ratios; tr wr or
! if no input O3 profile, set to reference value (or =1)
! Direct variables
      tr(i, layer)    = t(i, layer) / coef%tstar(layer)
      wr(i, layer)    = w(i, layer) / coef%wstar(layer)
! TL variables
      tr_tl(i, layer) = t_tl(i, layer) / coef%tstar(layer)
      wr_tl(i, layer) = w_tl(i, layer) / coef%wstar(layer)
! if no input O3 profile, set to reference value (or =1)

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        or(i, layer)    = o(i, layer) / coef%ostar(layer)
        or_tl(i, layer) = o_tl(i, layer) / coef%ostar(layer)
      ELSE
        or(i, layer)    = 1._JPRB
        or_tl(i, layer) = 0._JPRB
      ENDIF

    ENDDO

  ENDDO

! calculate profile / reference profile sums: tw ww ow
! if no input O3 profile, set to reference value (ow =1)
! Direct variables

  DO i = 1, nprofiles
    tw(i, 1)    = 0._JPRB
    tw_tl(i, 1) = 0._JPRB
  ENDDO


  DO layer = 2, prof(1)%nlayers

    DO i = 1, nprofiles
      tw(i, layer)    = tw(i, layer - 1) + coef%dpp(layer - 1) * tr(i, layer - 1)
      tw_tl(i, layer) = tw_tl(i, layer - 1) + coef%dpp(layer - 1) * tr_tl(i, layer - 1)
    ENDDO

  ENDDO

  sum1 = 0._JPRB
  sum2 = 0._JPRB

  DO layer = 1, prof_tl(1)%nlayers

    DO i = 1, nprofiles
      sum1(i)      = sum1(i) + coef%dpp(layer - 1) * w(i, layer)
      sum2(i)      = sum2(i) + coef%dpp(layer - 1) * coef%wstar(layer)
      ww(i, layer) = sum1(i) / sum2(i)
    ENDDO

  ENDDO

  sum1 = 0._JPRB
  sum2 = 0._JPRB

  DO layer = 1, prof_tl(1)%nlayers

    DO i = 1, nprofiles
      sum1(i)         = sum1(i) + coef%dpp(layer - 1) * w_tl(i, layer)
      sum2(i)         = sum2(i) + coef%dpp(layer - 1) * coef%wstar(layer)
      ww_tl(i, layer) = sum1(i) / sum2(i)
    ENDDO

  ENDDO

! if no input O3 profile, set to reference value (ow =1)
  IF (coef%nozone > 0) THEN
    sum1 = 0._JPRB
    sum2 = 0._JPRB
  
    DO layer = 1, prof_tl(1)%nlayers
  
      DO i = 1, nprofiles
  
        IF (opts%ozone_Data) THEN
          sum1(i)      = sum1(i) + coef%dpp(layer - 1) * o(i, layer)
          sum2(i)      = sum2(i) + coef%dpp(layer - 1) * coef%ostar(layer)
          ow(i, layer) = sum1(i) / sum2(i)
        ELSE
          ow(i, layer) = 1._JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  
    sum1 = 0._JPRB
    sum2 = 0._JPRB
  
    DO layer = 1, prof_tl(1)%nlayers
  
      DO i = 1, nprofiles
  
        IF (opts%ozone_Data) THEN
          sum1(i)         = sum1(i) + coef%dpp(layer - 1) * o_tl(i, layer)
          sum2(i)         = sum2(i) + coef%dpp(layer - 1) * coef%ostar(layer)
          ow_tl(i, layer) = sum1(i) / sum2(i)
        ELSE
          ow_tl(i, layer) = 0._JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  ENDIF

! ATTENTION
!  w_tl(:) = prof_tl % q(:)
!5) set predictors
!--
!5.1 mixed gases
!---

  IF (coef%id_inst == 10 .and. coef%IncZeeman) THEN
  ! SSMIS with Zeeman coefficient file
  ! geomagnetic field variables (Be, cosbk) are part of the user input

! X1 -> X9
    DO layer = 1, prof_tl(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles
! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        predictors_tl%mixedgas(1, layer, i) = raytracing_tl%pathsat(layer, i)
        predictors_tl%mixedgas(2, layer, i) = raytracing_tl%pathsat(layer, i) * 2 * raytracing%pathsat(layer, i)
        predictors_tl%mixedgas(3, layer, i) =      &
          & tr_tl(i, layer) * raytracing%pathsat(layer, i) + raytracing_tl%pathsat(layer, i) * tr(i, layer)
        predictors_tl%mixedgas(4, layer, i) = 2._JPRB * tr_tl(i, layer) * predictors%mixedgas(3, layer, i) +      &
          & raytracing_tl%pathsat(layer, i) * predictors%mixedgas(5, layer, i) ** 2
        predictors_tl%mixedgas(5, layer, i) = tr_tl(i, layer)
        predictors_tl%mixedgas(6, layer, i) = 2._JPRB * tr_tl(i, layer) * tr(i, layer)
        predictors_tl%mixedgas(7, layer, i) =      &
          & raytracing_tl%pathsat(layer, i) * tw(i, layer) + tw_tl(i, layer) * raytracing%pathsat(layer, i)
        predictors_tl%mixedgas(8, layer, i) = raytracing_tl%pathsat(layer, i) * tw(i, layer) / tr(i, layer) +      &
          & tw_tl(i, layer) * raytracing%pathsat(layer, i) / tr(i, layer) -                                        &
          & tr_tl(i, layer) * tw(i, layer) * raytracing%pathsat(layer, i) / tr(i, layer) ** 2._JPIM
        predictors_tl%mixedgas(9, layer, i) =      &
          & raytracing_tl%pathsat(layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, i)

! X10 - NB tw(i,1) set to zero in direct
        IF (layer == 1) THEN
          predictors_tl%mixedgas(10, layer, i) = 0._jprb
        ELSE
          predictors_tl%mixedgas(10, layer, i) =                                                                           &
            & raytracing_tl%pathsat(layer, i) * 0.5_JPRB * sqrt(sqrt(tw(i, layer))) / predictors%mixedgas(9, layer, i) +  &
            & tw_tl(i, layer) * 0.25_JPRB * predictors%mixedgas(9, layer, i) * sqrt(sqrt(tw(i, layer))) ** (-3_jpim)
        ENDIF

! X11 -> X21
! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
        Predictors_tl%mixedgas(11, layer, i)  = raytracing_tl%pathsat(layer, i)
        Predictors_tl%mixedgas(12, layer, i)  =  - (Predictors%mixedgas(12, layer, i) / t(i, layer)) * t_tl(i, layer) +      &
          & (300.0_JPRB / t(i, layer)) * raytracing_tl%pathsat(layer, i)
        Predictors_tl%mixedgas(13, layer, i)  = prof(i)%cosbk ** 2 * raytracing_tl%pathsat(layer, i)
        Predictors_tl%mixedgas(14, layer, i)  = Predictors_tl%mixedgas(12, layer, i) / prof(i)%Be
        Predictors_tl%mixedgas(15, layer, i)  = Predictors_tl%mixedgas(12, layer, i) * prof(i)%cosbk ** 2
        Predictors_tl%mixedgas(16, layer, i)  = raytracing_tl%pathsat(layer, i) / prof(i)%Be
        Predictors_tl%mixedgas(17, layer, i)  = Predictors_tl%mixedgas(16, layer, i) / prof(i)%Be
        Predictors_tl%mixedgas(18, layer, i)  = prof(i)%Be * raytracing_tl%pathsat(layer, i)
        Predictors_tl%mixedgas(19, layer, i)  = prof(i)%Be ** 3 * raytracing_tl%pathsat(layer, i)
        Predictors_tl%mixedgas(20, layer, i) = Predictors_tl%mixedgas(13, layer, i) * prof(i)%Be
        Predictors_tl%mixedgas(21, layer, i) = Predictors_tl%mixedgas(20, layer, i) * prof(i)%Be
      ENDDO

    ENDDO


  ELSE ! all sensors, except when running SSMIS with Zeeman coefficient file

! X1 -> X9
    DO layer = 1, prof_tl(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles
! all sensors, except when running SSMIS or AMSU-A with Zeeman coefficient file
        predictors_tl%mixedgas(1, layer, i) = raytracing_tl%pathsat(layer, i)
        predictors_tl%mixedgas(2, layer, i) = raytracing_tl%pathsat(layer, i) * 2 * raytracing%pathsat(layer, i)
        predictors_tl%mixedgas(3, layer, i) =      &
          & tr_tl(i, layer) * raytracing%pathsat(layer, i) + raytracing_tl%pathsat(layer, i) * tr(i, layer)
        predictors_tl%mixedgas(4, layer, i) = 2._JPRB * tr_tl(i, layer) * predictors%mixedgas(3, layer, i) +      &
          & raytracing_tl%pathsat(layer, i) * predictors%mixedgas(5, layer, i) ** 2
        predictors_tl%mixedgas(5, layer, i) = tr_tl(i, layer)
        predictors_tl%mixedgas(6, layer, i) = 2._JPRB * tr_tl(i, layer) * tr(i, layer)
        predictors_tl%mixedgas(7, layer, i) =      &
          & raytracing_tl%pathsat(layer, i) * tw(i, layer) + tw_tl(i, layer) * raytracing%pathsat(layer, i)
        predictors_tl%mixedgas(8, layer, i) = raytracing_tl%pathsat(layer, i) * tw(i, layer) / tr(i, layer) +      &
          & tw_tl(i, layer) * raytracing%pathsat(layer, i) / tr(i, layer) -                                        &
          & tr_tl(i, layer) * tw(i, layer) * raytracing%pathsat(layer, i) / tr(i, layer) ** 2._JPIM
        predictors_tl%mixedgas(9, layer, i) =      &
          & raytracing_tl%pathsat(layer, i) * 0.5_JPRB / predictors%mixedgas(9, layer, i)

! X10 -  NB tw(i,1) set to zero in direct
        IF (layer == 1) THEN
          predictors_tl%mixedgas(10, layer, i) = 0._jprb
        ELSE
          predictors_tl%mixedgas(10, layer, i) =                                                                           &
            & raytracing_tl%pathsat(layer, i) * 0.5_JPRB * sqrt(sqrt(tw(i, layer))) / predictors%mixedgas(9, layer, i) +  &
            & tw_tl(i, layer) * 0.25_JPRB * predictors%mixedgas(9, layer, i) * sqrt(sqrt(tw(i, layer))) **(-3_jpim)
        ENDIF


! X11 -> X14
        IF (coef%id_inst == 3 .AND. coef%IncZeeman) THEN
	! AMSU-A with Zeeman coefficient file
! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
          predictors_tl%mixedgas(11, layer, i) = prof(i)%cosbk ** 2 * raytracing_tl%pathsat(layer, i)
          predictors_tl%mixedgas(12, layer, i) =      &
            & 2.0_JPRB * prof(i)%Be * raytracing%pathsat(layer, i) * raytracing_tl%pathsat(layer, i)
          predictors_tl%mixedgas(13, layer, i) = prof(i)%Be ** 3 * raytracing_tl%pathsat(layer, i)
          predictors_tl%mixedgas(14, layer, i) =      &
            & 2.0_JPRB * (prof(i)%cosbk * prof(i)%Be) ** 2 * raytracing%pathsat(layer, i) * raytracing_tl%pathsat(layer, i)
! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
        ENDIF

      ENDDO

    ENDDO

  ENDIF

!5.2 water vapour  ( numbers in right hand are predictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------

  DO layer = 1, prof_tl(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles
      sec_wr(layer)                           = raytracing%pathsat(layer, i) * wr(i, layer)
      sec_wr_tl(layer)                        =      &
        & raytracing_tl%pathsat(layer, i) * wr(i, layer) + raytracing%pathsat(layer, i) * wr_tl(i, layer)
      predictors_tl%watervapour(1, layer, i)  = sec_wr_tl(layer)!  7
      predictors_tl%watervapour(2, layer, i)  = 0.5_JPRB * sec_wr_tl(layer) / predictors%watervapour(2, layer, i) !  5
      zrecip = 1.0_JPRB / predictors%watervapour(1, layer, i)
      predictors_tl%watervapour(3, layer, i)  =                                                              &
        & sec_wr_tl(layer) * wr(i, layer) / ww(i, layer) + wr_tl(i, layer) * sec_wr(layer) / ww(i, layer) -  &
        & ww_tl(i, layer) * wr(i, layer) * sec_wr(layer) / ww(i, layer) ** 2
      predictors_tl%watervapour(4, layer, i)  = sec_wr_tl(layer) * predictors%watervapour(4, layer, i) * zrecip +      &
        & predictors%watervapour(1, layer, i) * dt_tl(i, layer)
      predictors_tl%watervapour(5, layer, i)  = 2 * predictors%watervapour(1, layer, i) * sec_wr_tl(layer)
      predictors_tl%watervapour(6, layer, i)  = predictors%watervapour(2, layer, i) * dt_tl(i, layer) +      &
        & 0.5_JPRB * sec_wr_tl(layer) * predictors%watervapour(6, layer, i) * zrecip
      predictors_tl%watervapour(7, layer, i)  = 0.25_JPRB * sec_wr_tl(layer) / predictors%watervapour(7, layer, i) ** 3
      predictors_tl%watervapour(8, layer, i)  = predictors_tl%watervapour(2, layer, i) * wr(i, layer) / ww(i, layer) +      &
        & wr_tl(i, layer) * predictors%watervapour(2, layer, i) / ww(i, layer) -                                            &
        & ww_tl(i, layer) * predictors%watervapour(2, layer, i) * wr(i, layer) / ww(i, layer) ** 2
      predictors_tl%watervapour(9, layer, i)  = 3 * sec_wr_tl(layer) * predictors%watervapour(5, layer, i)
      predictors_tl%watervapour(10, layer, i) = 4 * sec_wr_tl(layer) * predictors%watervapour(9, layer, i)
      predictors_tl%watervapour(11, layer, i) =      &
        & Abs(dt(i, layer)) * (sec_wr_tl(layer) * dt(i, layer) + 2 * sec_wr(layer) * dt_tl(i, layer))
      zsqrt = Sqrt(predictors%watervapour(13, layer, i))
      predictors_tl%watervapour(12, layer, i) =                                                                &
        & 4 * raytracing%pathsat(layer, i) * ww_tl(i, layer) * predictors%watervapour(12, layer, i) / zsqrt +  &
        & 4 * raytracing_tl%pathsat(layer, i) * ww(i, layer) * predictors%watervapour(12, layer, i) / zsqrt
      predictors_tl%watervapour(13, layer, i) = 2 * raytracing%pathsat(layer, i) * ww_tl(i, layer) * zsqrt +      &
        & 2 * raytracing_tl%pathsat(layer, i) * ww(i, layer) * zsqrt
      predictors_tl%watervapour(14, layer, i) =                                                              &
        & sec_wr_tl(layer) * wr(i, layer) / tr(i, layer) + wr_tl(i, layer) * sec_wr(layer) / tr(i, layer) -  &
        & tr_tl(i, layer) * sec_wr(layer) * wr(i, layer) / tr(i, layer) ** 2
      predictors_tl%watervapour(15, layer, i) =                                                                        &
        & sec_wr_tl(layer) * wr(i, layer) / tr(i, layer) ** 4 + wr_tl(i, layer) * sec_wr(layer) / tr(i, layer) ** 4 -  &
        & tr_tl(i, layer) * 4 * sec_wr(layer) * wr(i, layer) / tr(i, layer) ** 5
    ENDDO

  ENDDO

!5.3 ozone
!---------

  IF (coef%nozone > 0) THEN

    DO layer = 1, prof_tl(1)%nlayers
      level = layer + 1

      DO i = 1, nprofiles
        sec_or(layer)                     = raytracing%pathsat(layer, i) * or(i, layer)
        sec_or_tl(layer)                  =      &
          & raytracing_tl%pathsat(layer, i) * or(i, layer) + raytracing%pathsat(layer, i) * or_tl(i, layer)
        predictors_tl%ozone(1, layer, i)  = sec_or_tl(layer)
        predictors_tl%ozone(2, layer, i)  = 0.5_JPRB * sec_or_tl(layer) / predictors%ozone(2, layer, i)
        predictors_tl%ozone(3, layer, i)  =                                                     &
          & sec_or_tl(layer) * predictors%ozone(3, layer, i) / predictors%ozone(1, layer, i) +  &
          & predictors%ozone(1, layer, i) * dto_tl(i, layer)
        predictors_tl%ozone(4, layer, i)  = 2 * sec_or_tl(layer) * predictors%ozone(1, layer, i)
        predictors_tl%ozone(5, layer, i)  = 0.5_JPRB * sec_or_tl(layer) * predictors%ozone(3, layer, i) /      &
          & (predictors%ozone(1, layer, i) * predictors%ozone(2, layer, i)) +                                  &
          & predictors%ozone(2, layer, i) * dto_tl(i, layer)
        predictors_tl%ozone(6, layer, i)  =                                                                    &
          & sec_or_tl(layer) * or(i, layer) * ow(i, layer) + or_tl(i, layer) * sec_or(layer) * ow(i, layer) +  &
          & ow_tl(i, layer) * sec_or(layer) * or(i, layer)
        predictors_tl%ozone(7, layer, i)  = predictors_tl%ozone(2, layer, i) * or(i, layer) / ow(i, layer) +      &
          & or_tl(i, layer) * predictors%ozone(2, layer, i) / ow(i, layer) -                                      &
          & ow_tl(i, layer) * predictors%ozone(2, layer, i) * or(i, layer) / ow(i, layer) ** 2
        predictors_tl%ozone(8, layer, i)  = sec_or_tl(layer) * ow(i, layer) + sec_or(layer) * ow_tl(i, layer)
        zsqrt = Sqrt(predictors%ozone(10, layer, i))
        predictors_tl%ozone(9, layer, i)  = raytracing%pathsat(layer, i) * or_tl(i, layer) * zsqrt +        &
          & predictors%ozone(1, layer, i) * 0.5 * raytracing%pathsat(layer, i) * ow_tl(i, layer) / zsqrt +  &
          & 1.5 * raytracing_tl%pathsat(layer, i) * or(i, layer) * zsqrt
        predictors_tl%ozone(10, layer, i) =      &
          & raytracing_tl%pathsat(layer, i) * ow(i, layer) + raytracing%pathsat(layer, i) * ow_tl(i, layer)
        predictors_tl%ozone(11, layer, i) =                                                        &
          & 2 * ow_tl(i, layer) * raytracing%pathsat(layer, i) * predictors%ozone(10, layer, i) +  &
          & 2 * raytracing_tl%pathsat(layer, i) * ow(i, layer) * predictors%ozone(10, layer, i)
      ENDDO

    ENDDO

  ENDIF

!5.4 cloud
!---------
!CDIR NOLOOPCHG

  DO layer = 1, prof(1)%nlayers
    level         = layer + 1
    deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) / (4.3429_JPRB * gravity)

    DO i = 1, nprofiles

      IF (opts%clw_Data) THEN
        predictors_tl%clw(layer, i) = deltac(layer) * prof_tl(i)%clw(level) * geom(i)%seczen
      ELSE
        predictors_tl%clw(layer, i) = 0._jprb
      ENDIF

    ENDDO

  ENDDO

!CDIR NOLOOPCHG

  DO layer = 2, prof_tl(1)%nlayers
    level = layer + 1

    DO i = 1, nprofiles

      IF (opts%clw_Data) THEN
        predictors_tl%clw(layer, i) =      &
          & 0.5_JPRB * (predictors_tl%clw(layer, i) + deltac(layer) * prof_tl(i)%clw(level - 1) * geom(i)%seczen)
      ENDIF

    ENDDO

  ENDDO


  DO i = 1, nprofiles

    IF (.NOT. opts%clw_Data) THEN
      predictors_tl%ncloud = 0
    ENDIF

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7_tl
