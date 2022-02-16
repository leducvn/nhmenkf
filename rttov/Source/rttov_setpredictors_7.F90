!
SUBROUTINE rttov_setpredictors_7( &
            & opts,       &
            & prof,       &
            & geom,       &
            & coef,       &
            & aux,        &
            & predictors, &
            & raytracing)
!
! Description
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
!  1.3       07/12/2005  Add surface humidity (R. Saunders)
!  1.4       16/01/2006  Marco Matricardi (ECMWF):
!               --       Altitude dependent local zenith angle introduced.
!  1.5       22/08/2007  Optimised (D Salmond)
!  1.6       12/08/2008  Zeeman effect based on Yong Han (P. Rayer)
!  1.7       27/02/2009  Profile levels to include ToA. Distinguish arrays
!                        in raytracing (on levels) from all others (on
!                        layers). Predictors prepared to maintain agreement
!                        with the RTTOV-7 scheme in RTTOV-9 (P. Rayer)
!  1.8       02/12/2009  Pathsat, Pathsun and related quantities are now layer arrays
!                        (Marco Matricardi).
!  1.9       17/06/2010  Combined non-Zeeman and Zeeman predictors for SSMIS for
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
  USE rttov_const, ONLY : gravity, sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options  ), INTENT(IN)    :: opts
  TYPE(profile_Type   ), INTENT(IN)    :: prof(:)         ! profile
  TYPE(rttov_coef     ), INTENT(IN)    :: coef            ! coefficients
  TYPE(geometry_Type  ), INTENT(IN)    :: geom(size(prof))! geometry
  TYPE(predictors_Type), INTENT(INOUT) :: predictors      ! predictors
  TYPE(profile_aux    ), INTENT(IN)    :: aux             ! auxillary profiles info.
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing
!INTF_END
!local variables:
  INTEGER(KIND=jpim) :: level , layer
  INTEGER(KIND=jpim) :: iv2lev, iv3lev, iv2lay, i
! user profile
  REAL   (KIND=jprb) :: t(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: w(size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: o(size(prof)     , prof(1)%nlayers)
! reference profile
  REAL   (KIND=jprb) :: tr    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: wr    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: or    (size(prof)     , prof(1)%nlayers)
! user - reference
  REAL   (KIND=jprb) :: dt    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: dto   (size(prof)     , prof(1)%nlayers)
! pressure weighted
  REAL   (KIND=jprb) :: tw    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ww    (size(prof)     , prof(1)%nlayers)
  REAL   (KIND=jprb) :: ow    (size(prof)     , prof(1)%nlayers)
! intermediate variables
  REAL   (KIND=jprb) :: sum1  (size(prof)                      ), sum2(size(prof)), ztemp
  REAL   (KIND=jprb) :: deltac(prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_or(prof(1)%nlayers                 )
  REAL   (KIND=jprb) :: sec_wr(prof(1)%nlayers                 )
  INTEGER(KIND=jpim) :: nprofiles                                                         ! Number of profiles
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
! 1 profile layer quantities
!   the layer number agrees with the level number of its upper boundary
! layer N-1 lies between levels N-1 and N
  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)

  DO i = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      level       = layer + 1
      t(i, layer) = (prof(i)%t(level - 1) + prof(i)%t(level)) * 0.5_JPRB
      w(i, layer) = (prof(i)%q(level - 1) + prof(i)%q(level)) * 0.5_JPRB
    ENDDO

  ENDDO


  IF (opts%use_q2m) THEN

    DO i = 1, nprofiles
! include surface humidity
      iv3lev = aux%s(i)%nearestlev_surf - 1! nearest level above surface
      iv2lev = aux%s(i)%nearestlev_surf    ! nearest level above surface

      IF (iv2lev <= coef%nlevels) THEN
        iv2lay       = iv2lev - 1
        w(i, iv2lay) = (prof(i)%s2m%q + prof(i)%q(iv3lev)) * 0.5_JPRB
      ENDIF

    ENDDO

  ENDIF


  DO i = 1, nprofiles

    IF (opts%ozone_Data .AND. coef%nozone > 0) THEN

      DO layer = 1, prof(1)%nlayers
        level       = layer + 1
        o(i, layer) = (prof(i)%o3(level - 1) + prof(i)%o3(level)) * 0.5_JPRB
      ENDDO

    ENDIF

  ENDDO


! 2 calculate, for layers, deviations from reference profile
! if no input O3 profile we still use the input temperature profile for dto

  DO i = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      dt(i, layer) = t(i, layer) - coef%tstar(layer)
      IF (coef%nozone > 0) dto(i, layer) = t(i, layer) - coef%to3star(layer)
    ENDDO

! 3 calculate (profile / reference profile) ratios; tr wr or

    DO layer = 1, prof(1)%nlayers
      tr(i, layer) = t(i, layer) / coef%tstar(layer)
      wr(i, layer) = w(i, layer) / coef%wstar(layer)
! if no input O3 profile, set to reference value (or =1)

      IF (opts%ozone_Data .AND. coef%nozone > 0) THEN
        or(i, layer) = o(i, layer) / coef%ostar(layer)
      ELSE
        or(i, layer) = 1._JPRB
      ENDIF

    ENDDO

  ENDDO

! 4 calculate profile / reference profile sums: tw ww ow

  DO i = 1, nprofiles
    tw(i, 1) = 0._JPRB
  ENDDO


  DO i = 1, nprofiles
    tw(i, 1) = 0._JPRB

    DO layer = 2, prof(1)%nlayers
! cumulate overlying layers: weighting tr relates to same layer as dpp
! do not need dpp(0) to start
      tw(i, layer) = tw(i, layer - 1) + coef%dpp(layer - 1) * tr(i, layer - 1)
    ENDDO

  ENDDO

  sum1 = 0._JPRB
  sum2 = 0._JPRB

  DO i = 1, nprofiles
! cumulating column overlying layer and layer itself

    DO layer = 1, prof(1)%nlayers
! cumulate overlying layers: weighting w or wstar relates to layer below dpp
! need dpp(0) to start
      sum1(i)      = sum1(i) + coef%dpp(layer - 1) * w(i, layer)
      sum2(i)      = sum2(i) + coef%dpp(layer - 1) * coef%wstar(layer)
      ww(i, layer) = sum1(i) / sum2(i)
    ENDDO

  ENDDO

! if no input O3 profile, set to reference value (ow =1)
  IF (coef%nozone > 0) THEN
    sum1 = 0._JPRB
    sum2 = 0._JPRB
  
    DO i = 1, nprofiles
  
      DO layer = 1, prof(1)%nlayers
  
        IF (opts%ozone_Data) THEN
! cumulate overlying layers: weighting o or ostar relates to layer below dpp
! need dpp(0) to start
          sum1(i)      = sum1(i) + coef%dpp(layer - 1) * o(i, layer)
          sum2(i)      = sum2(i) + coef%dpp(layer - 1) * coef%ostar(layer)
          ow(i, layer) = sum1(i) / sum2(i)
        ELSE
          ow(i, layer) = 1._JPRB
        ENDIF
  
      ENDDO
  
    ENDDO
  ENDIF
  
! for other minor gases do as for O3 for testing presence of
! coefficients and profile values
!
!5) set predictors
!--
!5.1 mixed gases
!---

  IF (coef%id_inst == 10 .AND. coef%IncZeeman) THEN
  ! SSMIS with Zeeman coefficient file
  ! geomagnetic field variables (Be, cosbk) are part of the user input

    DO i = 1, nprofiles

      DO layer = 1, prof(1)%nlayers
        level = layer + 1! raytracing % pathsat(layer,i) is  angle at lower boundary of layer

 ! only effective for non-Zeeman chans 1-18,23-24 - coefficient file will have zeros for chan 19-22
        predictors%mixedgas(1, layer, i)  = raytracing%pathsat(layer, i)
        predictors%mixedgas(2, layer, i)  = raytracing%pathsat(layer, i) * raytracing%pathsat(layer, i)
        predictors%mixedgas(3, layer, i)  = raytracing%pathsat(layer, i) * tr(i, layer)
        predictors%mixedgas(4, layer, i)  = raytracing%pathsat(layer, i) * tr(i, layer) * tr(i, layer)
        predictors%mixedgas(5, layer, i)  = tr(i, layer)
        predictors%mixedgas(6, layer, i)  = tr(i, layer) * tr(i, layer)
        predictors%mixedgas(7, layer, i)  = raytracing%pathsat(layer, i) * tw(i, layer)
        predictors%mixedgas(8, layer, i)  = raytracing%pathsat(layer, i) * tw(i, layer) / tr(i, layer)
        predictors%mixedgas(9, layer, i)  = Sqrt(raytracing%pathsat(layer, i))
        predictors%mixedgas(10, layer, i) = predictors%mixedgas(9, layer, i) * sqrt(sqrt(tw(i, layer)))

! only effective for Zeeman chans 19-22 - coefficient file will have zeros for chan 1-18,23-24
! NB require prof(i) % Be >0. (divisor)
        Predictors % mixedgas(11, layer,i) = raytracing % pathsat(layer,i)
        Predictors % mixedgas(12, layer,i) = (300.0_JPRB/t(i, layer)) * raytracing % pathsat(layer,i)
        Predictors % mixedgas(13, layer,i) = prof(i) % cosbk**2 * raytracing % pathsat(layer,i)    
        Predictors % mixedgas(14, layer,i) = Predictors % mixedgas(12,layer,i) / prof(i) % Be
        Predictors % mixedgas(15, layer,i) = Predictors % mixedgas(12,layer,i) * prof(i) % cosbk**2 
        Predictors % mixedgas(16, layer,i) = raytracing % pathsat(layer,i) / prof(i) % Be
        Predictors % mixedgas(17, layer,i) = Predictors % mixedgas(16,layer,i) / prof(i) % Be
        Predictors % mixedgas(18, layer,i) = prof(i) % Be * raytracing % pathsat(layer,i)
        Predictors % mixedgas(19, layer,i) = prof(i) % Be**3 * raytracing % pathsat(layer,i)
        Predictors % mixedgas(20,layer,i) = Predictors % mixedgas(13,layer,i) * prof(i) % Be
        Predictors % mixedgas(21,layer,i) = Predictors % mixedgas(20,layer,i) * prof(i) % Be

      ENDDO

    ENDDO
 
  ELSE ! all sensors, except when running SSMIS with Zeeman coefficient file

   DO i = 1, nprofiles

      DO layer = 1, prof(1)%nlayers
        level = layer + 1! raytracing % pathsat(layer,i) is angle at lower boundary of layer
! all sensors, except when running SSMIS or AMSU-A with Zeeman coefficient file
        predictors%mixedgas(1, layer, i)  = raytracing%pathsat(layer, i)
        predictors%mixedgas(2, layer, i)  = raytracing%pathsat(layer, i) * raytracing%pathsat(layer, i)
        predictors%mixedgas(3, layer, i)  = raytracing%pathsat(layer, i) * tr(i, layer)
        predictors%mixedgas(4, layer, i)  = raytracing%pathsat(layer, i) * tr(i, layer) * tr(i, layer)
        predictors%mixedgas(5, layer, i)  = tr(i, layer)
        predictors%mixedgas(6, layer, i)  = tr(i, layer) * tr(i, layer)
        predictors%mixedgas(7, layer, i)  = raytracing%pathsat(layer, i) * tw(i, layer)
        predictors%mixedgas(8, layer, i)  = raytracing%pathsat(layer, i) * tw(i, layer) / tr(i, layer)
        predictors%mixedgas(9, layer, i)  = Sqrt(raytracing%pathsat(layer, i))
        predictors%mixedgas(10, layer, i) = predictors%mixedgas(9, layer, i) * sqrt(sqrt(tw(i, layer)))

        IF (coef%id_inst == 3 .AND. coef%IncZeeman) THEN
	! AMSU-A with Zeeman coefficient file
! only effective for Zeeman chan 14 - coefficient file will have zeros for chan 1-13
          predictors%mixedgas(11, layer, i) = prof(i)%cosbk ** 2 * raytracing%pathsat(layer, i)
          predictors%mixedgas(12, layer, i) = prof(i)%Be * raytracing%pathsat(layer, i) ** 2
          predictors%mixedgas(13, layer, i) = prof(i)%Be ** 3 * raytracing%pathsat(layer, i)
          predictors%mixedgas(14, layer, i) = (prof(i)%cosbk * prof(i)%Be * raytracing%pathsat(layer, i)) ** 2
! NB some of YH's original predictors omitted - effectively duplicated by predictors 1-4 above
        ENDIF

      ENDDO

    ENDDO

  ENDIF

!5.2 water vapour  ( numbers in right hand are predictor numbers
! in the reference document for RTTOV7 (science and validation report)
!----------------

  DO i = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
      level = layer + 1! raytracing % pathsat(layer,i) is angle at lower boundary of layer
      ztemp = 1.0_JPRB / ww(i, layer)
      sec_wr(layer)                        = raytracing%pathsat(layer, i) * wr(i, layer)
      predictors%watervapour(1, layer, i)  = sec_wr(layer)                                             !  7
      predictors%watervapour(2, layer, i)  = Sqrt(sec_wr(layer))                                       !  5
      predictors%watervapour(3, layer, i)  = sec_wr(layer) * wr(i, layer) * ztemp                      ! 12
      predictors%watervapour(4, layer, i)  = sec_wr(layer) * dt(i, layer)                              !  4
      predictors%watervapour(5, layer, i)  = sec_wr(layer) * sec_wr(layer)                             !  1
      predictors%watervapour(6, layer, i)  = predictors%watervapour(2, layer, i) * dt(i, layer)        ! 11
      predictors%watervapour(7, layer, i)  = Sqrt(predictors%watervapour(2, layer, i))                 ! 6
      predictors%watervapour(8, layer, i)  = predictors%watervapour(2, layer, i) * wr(i, layer) * ztemp! 13
      predictors%watervapour(9, layer, i)  = predictors%watervapour(5, layer, i) * sec_wr(layer)       ! 8
      predictors%watervapour(10, layer, i) = predictors%watervapour(9, layer, i) * sec_wr(layer)       ! 9
      predictors%watervapour(11, layer, i) = sec_wr(layer) * dt(i, layer) * Abs(dt(i, layer))          ! 10
      predictors%watervapour(12, layer, i) = (raytracing%pathsat(layer, i) * ww(i, layer)) ** 4        ! 3
      predictors%watervapour(13, layer, i) = (raytracing%pathsat(layer, i) * ww(i, layer)) ** 2        ! 2
      ztemp = 1.0_JPRB / tr(i, layer)
      predictors%watervapour(14, layer, i) = sec_wr(layer) * wr(i, layer) * ztemp                      ! 14
      predictors%watervapour(15, layer, i) = sec_wr(layer) * wr(i, layer) * ztemp ** 4
    ENDDO

  ENDDO

!
!5.3 ozone
!---------
! if no input O3 profile, variables or, ow and dto have been set
! to the reference profile values (1, 1, 0)

  IF (coef%nozone > 0) THEN

    DO i = 1, nprofiles

      DO layer = 1, prof(1)%nlayers
        level = layer + 1! raytracing % pathsat(layer,i) is angle at lower boundary of layer
        sec_or(layer)                  = raytracing%pathsat(layer, i) * or(i, layer)
        predictors%ozone(1, layer, i)  = sec_or(layer)
        predictors%ozone(2, layer, i)  = Sqrt(sec_or(layer))
        predictors%ozone(3, layer, i)  = sec_or(layer) * dto(i, layer)
        predictors%ozone(4, layer, i)  = sec_or(layer) * sec_or(layer)
        predictors%ozone(5, layer, i)  = predictors%ozone(2, layer, i) * dto(i, layer)
        predictors%ozone(6, layer, i)  = sec_or(layer) * or(i, layer) * ow(i, layer)
        predictors%ozone(7, layer, i)  = predictors%ozone(2, layer, i) * or(i, layer) / ow(i, layer)
        predictors%ozone(8, layer, i)  = sec_or(layer) * ow(i, layer)
        predictors%ozone(9, layer, i)  = sec_or(layer) * Sqrt(raytracing%pathsat(layer, i) * ow(i, layer))
        predictors%ozone(10, layer, i) = raytracing%pathsat(layer, i) * ow(i, layer)
        predictors%ozone(11, layer, i) = (raytracing%pathsat(layer, i) * ow(i, layer)) ** 2
      ENDDO

    ENDDO

  ENDIF

!
!5.4 cloud
!---------
  ztemp = 1.0_JPRB / (4.3429_JPRB * gravity)

  DO i = 1, nprofiles

    DO layer = 1, prof(1)%nlayers
! NB in RTTOV-10, layer 1 is the layer above level 2
      level         = layer + 1
      deltac(layer) = 0.1820_JPRB * 100.0_JPRB * coef%dp(layer) * ztemp

      IF (opts%clw_Data .AND. coef%id_sensor == sensor_id_mw) THEN
        predictors%clw(layer, i) = deltac(layer) * prof(i)%clw(level) * geom(i)%seczen
      ELSE
        predictors%clw(layer, i) = 0._jprb
      ENDIF

    ENDDO

  ENDDO


  DO i = 1, nprofiles

    DO layer = 2, prof(1)%nlayers
      level = layer + 1

      IF (opts%clw_Data .AND. coef%id_sensor == sensor_id_mw) THEN
        predictors%clw(layer, i) =      &
          & 0.5_JPRB * (predictors%clw(layer, i) + deltac(layer) * prof(i)%clw(level - 1) * geom(i)%seczen)
      ENDIF

    ENDDO

  ENDDO


  DO i = 1, nprofiles

    IF (opts%clw_Data .AND. coef%id_sensor == sensor_id_mw) THEN
      predictors%ncloud = 1
    ENDIF

  ENDDO

  IF (LHOOK) CALL DR_HOOK('RTTOV_SETPREDICTORS_7', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_setpredictors_7
