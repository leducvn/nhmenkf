SUBROUTINE rttov_transmit_k( &
            & addaerosl,                      &
            & addclouds,                      &
            & nlayers,                        &
            & chanprof,                       &
            & aux,                            &
            & aux_k,                          &
            & coef,                           &
            & ircld,                          &
            & opdp_path,                      &
            & opdp_path_k,                    &
            & od_level,                       &
            & transmission,                   &
            & transmission_k,                 &
            & transmission_aux,               &
            & transmission_aux_k,             &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k, &
            & tau_ref,                        &
            & tau_ref_surf,                   &
            & tau_surf,                       &
            & tau_level)
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
!    Copyright 2010, EUMETSAT, All Rights Reserved.
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit_k.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
! 1.1     06/02/2007  Removed polaristion index R Saunders
! 1.2     11/12/2007  Input also adjoint of tau_total and tau_levels since
!                     these are the external outputs of RTTOV (A Geer)
! 1.3     04/06/2008  Fix od_frac and od_frac_ac calculation near surface level (PB PM)
! 1.4     15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
! 1.5     03/11/2009  Transmittances on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Adjoint variables
! input transmission_aux_k% tau_surf and transmission_aux_k% tau_level set inside integrate_k
!
! input/output aux_k
!
! output predictors_k initialised inside rttov_k (need input
!    intent for memory allocation in calling routine)
!
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type,          &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_hi
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path                                   ! path  optical depths
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_k
  TYPE(transmission_Type         ), INTENT(IN)    :: transmission
  TYPE(transmission_Type         ), INTENT(INOUT) :: transmission_k
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_k
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref_surf(size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_level    (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_surf    (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_level   (nlayers + 1   , size(chanprof))
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_level_k      (nlayers + 1   , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_k       (size(chanprof)                )
  REAL   (KIND=jprb) :: od_surf_ac_k    (size(chanprof)                )
  REAL   (KIND=jprb) :: od_frac_k       (size(chanprof)                )
  REAL   (KIND=jprb) :: tau_surf_k      (size(chanprof)                )! sat to surface transmission_aux at each frequency
  REAL   (KIND=jprb) :: tau_level_k     (nlayers + 1   , size(chanprof))! sat to level transmission_aux at each frequency
  REAL   (KIND=jprb) :: od_singlelayer_k(nlayers       , size(chanprof))
  REAL   (KIND=jprb) :: ztemp
  INTEGER(KIND=jpim) :: lev         , lay    , chan   , j, ist
  INTEGER(KIND=jpim) :: prof        , nlevels, levsurf
  INTEGER(KIND=jpim) :: nchannels                                       ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!---------------------------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_K', 0_jpim, ZHOOK_HANDLE)
  nchannels             = size(chanprof)
  nlevels               = nlayers + 1
!-----------------------------------------------------
!K of store transmittances for other streams
!-----------------------------------------------------
  od_level_k(:,:)       = 0._JPRB
  tau_level_k(:,:)      = 0.0_JPRB
  od_singlelayer_k(:,:) = 0.0_JPRB
  tau_surf_k(:)         = 0.0_JPRB
  od_frac_k(:)          = 0.0_JPRB
  od_surf_k(:)          = 0.0_JPRB
  od_surf_ac_k(:)       = 0.0_JPRB
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      prof = chanprof(j)%prof! Profile index
      DO ist = ircld%nstream(prof), 0,  - 1
        transmission_aux_k%tau_surf_t(ist, j) =      &
          & transmission_aux_k%tau_surf_t(ist, j) + transmission_aux_k%tau_surf(ist, j)
        transmission_aux_k%od_frac_t(ist, j)  =      &
          & transmission_aux_k%od_frac_t(ist, j) + transmission_aux_k%od_sfrac(ist, j)
        DO lay = 1, nlayers
          od_singlelayer_k(lay, j)                            =      &
            & od_singlelayer_k(lay, j) + transmission_aux_k%od_singlelayer(lay, ist, j)
          transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) =      &
            & transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) + transmission_aux_k%od_singlelayer(lay, ist, j)
        ENDDO
        DO lev = nlevels, 1,  - 1
          IF (tau_level(lev, j) >= 0._jprb) THEN
            ztemp = exp( - transmission_scatt_ir_stream%opdpac(ist, j, lev))
            tau_level_k(lev, j)                                =      &
              & tau_level_k(lev, j) + transmission_aux_k%tau_level(lev, ist, j) * ztemp
            transmission_scatt_ir_stream_k%opdpac(ist, j, lev) = transmission_scatt_ir_stream_k%opdpac(ist, j, lev) -      &
              & transmission_aux_k%tau_level(lev, ist, j) * tau_level(lev, j) * ztemp
          ELSE
            tau_level_k(lev, j) = tau_level_k(lev, j) + transmission_aux_k%tau_level(lev, ist, j)
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ELSE
    DO j = 1, nchannels
      prof = chanprof(j)%prof! Profile index
      DO ist = ircld%nstream(prof), 0,  - 1
        DO lev = nlevels, 1,  - 1
          tau_level_k(lev, j) = tau_level_k(lev, j) + transmission_aux_k%tau_level(lev, ist, j)
        ENDDO
        DO lay = nlayers, 1,  - 1
          od_singlelayer_k(lay, j) = od_singlelayer_k(lay, j) + transmission_aux_k%od_singlelayer(lay, ist, j)
        ENDDO
        od_frac_k(j)  = od_frac_k(j) - transmission_aux_k%od_sfrac(ist, j)
        tau_surf_k(j) = tau_surf_k(j) + transmission_aux_k%tau_surf(ist, j)
      ENDDO
      tau_level_k(:, j)               = tau_level_k(:, j) + transmission_k%tau_levels(:, j)
      transmission_k%tau_levels(:, j) = 0
      tau_surf_k(j)                   = tau_surf_k(j) + transmission_k%tau_total(j)
      transmission_k%tau_total(j)     = 0
    ENDDO
  ENDIF
  DO j = 1, nchannels
    DO ist = 0, ircld%nstream(prof)
      transmission_aux_k%tau_level(:, ist, j)      = 0.0_JPRB
      transmission_aux_k%od_singlelayer(:, ist, j) = 0.0_JPRB
      transmission_aux_k%tau_surf(ist, j)          = 0.0_JPRB
      transmission_aux_k%od_sfrac(ist, j)          = 0.0_JPRB
    ENDDO
  ENDDO
!---------------------------------------------------------------------------------------
!K of compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      prof    = chanprof(j)%prof           ! Profile index
      chan    = chanprof(j)%chan
      levsurf = aux%s(prof)%nearestlev_surf
!---Loop over the streams-----------------------------------------------------------------
      DO ist = ircld%nstream(prof), 0,  - 1
        IF (tau_surf(j) >= 0) THEN
          transmission_aux_k%tau_surf_ac(ist, j) =      &
            & transmission_aux_k%tau_surf_ac(ist, j) + transmission_aux_k%tau_surf_t(ist, j) * tau_surf(j)
          tau_surf_k(j)                          =      &
            & tau_surf_k(j) + transmission_aux_k%tau_surf_t(ist, j) * transmission_aux%tau_ref_surf_ac(ist, j)
        ELSE
          tau_surf_k(j) = tau_surf_k(j) + transmission_aux_k%tau_surf_t(ist, j)
        ENDIF
        od_frac_k(j)                          = od_frac_k(j) - transmission_aux_k%od_frac_t(ist, j)
        transmission_aux_k%od_frac_ac(ist, j) =      &
          & transmission_aux_k%od_frac_ac(ist, j) + transmission_aux_k%od_frac_t(ist, j)
        od_surf_ac_k(j)                       =      &
          & od_surf_ac_k(j) - transmission_aux_k%tau_surf_ac(ist, j) * transmission_aux%tau_ref_surf_ac(ist, j)
        IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
          od_surf_ac_k(j)                                            =      &
            & od_surf_ac_k(j) + transmission_aux_k%od_frac_ac(ist, j)
          transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf - 1) =      &
            & transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf - 1) - transmission_aux_k%od_frac_ac(ist, j)
        ELSE
          od_surf_ac_k(j)                                        =      &
            & od_surf_ac_k(j) + transmission_aux_k%od_frac_ac(ist, j)
          transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf) =      &
            & transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf) - transmission_aux_k%od_frac_ac(ist, j)
        ENDIF
        transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf)     =      &
          & transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf) +      &
          & od_surf_ac_k(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
        transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf - 1) =      &
          & transmission_scatt_ir_stream_k%opdpac(ist, j, levsurf - 1) + od_surf_ac_k(j) * aux%s(prof)%pfraction_surf
        aux_k%s(j)%pfraction_surf                                  = aux_k%s(j)%pfraction_surf + od_surf_ac_k(j) * (     &
          & transmission_scatt_ir_stream%opdpac(ist, j, levsurf - 1) - transmission_scatt_ir_stream%opdpac(ist, j, levsurf))
        od_surf_ac_k(j)                                            = 0._JPRB
      ENDDO
    ENDDO
  ENDIF
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tau_ref_surf(j) < coef%tt_a0(chan)) THEN
          tau_surf_k(j) = 0._jprb
        ENDIF
      ENDIF
    ENDDO
  ENDIF
  DO j = 1, nchannels
    prof         = chanprof(j)%prof
    levsurf      = aux%s(prof)%nearestlev_surf
    od_surf_k(j) = od_surf_k(j) + tau_surf_k(j) * tau_ref_surf(j)
    IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
      od_surf_k(j)               = od_surf_k(j) + od_frac_k(j)
      od_level_k(levsurf - 1, j) = od_level_k(levsurf - 1, j) - od_frac_k(j)
    ELSE
      od_surf_k(j)           = od_surf_k(j) + od_frac_k(j)
      od_level_k(levsurf, j) = od_level_k(levsurf, j) - od_frac_k(j)
    ENDIF
    od_level_k(levsurf, j)     = od_level_k(levsurf, j) + od_surf_k(j) * (1._JPRB - aux%s(prof)%pfraction_surf)
    od_level_k(levsurf - 1, j) = od_level_k(levsurf - 1, j) + od_surf_k(j) * aux%s(prof)%pfraction_surf
    aux_k%s(j)%pfraction_surf  =      &
      & aux_k%s(j)%pfraction_surf + od_surf_k(j) * (od_level(levsurf - 1, j) - od_level(levsurf, j))
    od_surf_k(j)               = 0._JPRB
  ENDDO
!-------------------------------------------
!K of assemble layer optical depths
!-------------------------------------------
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tau_ref(lev, j) < coef%tt_a0(chan)) THEN
            tau_level_k(lev, j) = 0._jprb
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    DO lev = 1, nlevels
      od_level_k(lev, j)            = od_level_k(lev, j) + tau_level_k(lev, j) * tau_ref(lev, j)
      od_level_k(lev, j)            = coef%ff_gam(chan) * od_level_k(lev, j)
      opdp_path_k%atm_level(lev, j) = opdp_path_k%atm_level(lev, j) + od_level_k(lev, j)
    ENDDO
    DO lay = 1, nlayers
      od_singlelayer_k(lay, j) = coef%ff_gam(chan) * od_singlelayer_k(lay, j)
    ENDDO
  ENDDO
! ad of level to space optical depths
! opdp_path_k % atm_level(:,:) = opdp_path_k % atm_level(:,:) + od_level_k(:,:)
! ad of single layer optical depths
  DO j = 1, nchannels
    DO lay = nlayers, 1,  - 1
      opdp_path_k%atm_level(lay, j)     = opdp_path_k%atm_level(lay, j) + od_singlelayer_k(lay, j)
      opdp_path_k%atm_level(lay + 1, j) = opdp_path_k%atm_level(lay + 1, j) - od_singlelayer_k(lay, j)
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_k
