SUBROUTINE rttov_transmit_tl( &
            & addaerosl,                       &
            & addclouds,                       &
            & nlayers,                         &
            & chanprof,                        &
            & aux,                             &
            & aux_tl,                          &
            & coef,                            &
            & ircld,                           &
            & opdp_path,                       &
            & opdp_path_tl,                    &
            & od_level,                        &
            & transmission,                    &
            & transmission_tl,                 &
            & transmission_aux,                &
            & transmission_aux_tl,             &
            & transmission_scatt_ir_stream,    &
            & transmission_scatt_ir_stream_tl, &
            & tau_ref,                         &
            & tau_ref_surf,                    &
            & tau_surf,                        &
            & tau_level)
!
! Description:
! Tangent linear of rttov_transmit
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEPTL and RTTAUTL from previous versions of RTTOV
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
!            --       New routine based on rttov_transmit_tl.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    29/01/2007  Removed polarisation R Saunders
!  1.2    11/12/2007  Return also TL of tau_total and tau_levels since
!                     these are the external outputs of RTTOV (A Geer)
!  1.3    04/06/2008  Fix od_frac and od_frac_ac calculation near surface level (PB PM)
!  1.4    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.5    03/11/2009  Transmittances / optical depths on levels (A Geer)
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
  TYPE(rttov_chanprof)            , INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(IN)    :: aux_tl
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path                                   ! path  optical depths
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path_tl                                ! path  optical depths
  TYPE(transmission_Type         ), INTENT(IN)    :: transmission
  TYPE(transmission_Type         ), INTENT(INOUT) :: transmission_tl
  TYPE(transmission_Type_aux     ), INTENT(IN)    :: transmission_aux
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux_tl
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_tl
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_ref_surf(size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(IN)    :: od_level    (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_level   (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(IN)    :: tau_surf    (size(chanprof)                )
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_level_tl      (nlayers + 1   , size(chanprof))
  REAL   (KIND=jprb) :: od_surf_tl       (size(chanprof)                )
  REAL   (KIND=jprb) :: od_surf_ac_tl    (size(chanprof)                )
  REAL   (KIND=jprb) :: od_frac_tl       (size(chanprof)                )
  REAL   (KIND=jprb) :: tau_surf_tl      (size(chanprof)                )! sat to surface transmission_aux at each frequency
  REAL   (KIND=jprb) :: tau_level_tl     (nlayers + 1   , size(chanprof))! sat to level transmission_aux at each frequency
  REAL   (KIND=jprb) :: od_singlelayer_tl(nlayers       , size(chanprof))
  REAL   (KIND=jprb) :: ztemp
  INTEGER(KIND=jpim) :: lev         , lay    , chan, j, levsurf
  INTEGER(KIND=jpim) :: prof        , nlevels
  INTEGER(KIND=jpim) :: ist
  INTEGER(KIND=jpim) :: nchannels                                        ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!--------------------------------------------------------------
!1. Assemble layer optical depths
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that values are sensible
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
! single layer optical depths - local variable
  DO j = 1, nchannels
    DO lay = 1, nlayers
      od_singlelayer_tl(lay, j) =  - (opdp_path_tl%atm_level(lay + 1, j) - opdp_path_tl%atm_level(lay, j))
    ENDDO
  ENDDO
! level to space optical depths - local variable
! od_level_tl(:,:) = opdp_path_tl % atm_level(:,:)
! gamma correction of local variables
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    DO lev = 1, nlevels
      od_level_tl(lev, j)  = coef%ff_gam(chan) * opdp_path_tl%atm_level(lev, j)
! associated transmittance
      tau_level_tl(lev, j) = od_level_tl(lev, j) * tau_ref(lev, j)
    ENDDO
    DO lay = 1, nlayers
      od_singlelayer_tl(lay, j) = coef%ff_gam(chan) * od_singlelayer_tl(lay, j)
    ENDDO
  ENDDO
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tau_ref(lev, j) < coef%tt_a0(chan)) THEN
            tau_level_tl(lev, j) = 0._jprb
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!----------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!----------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof          = chanprof(j)%prof
! as defined in rttov_profaux
    levsurf       = aux%s(prof)%nearestlev_surf
! layer above this
! arrays here based on layers, not levels
    od_surf_tl(j) =                                                                                                    &
      & od_level_tl(levsurf, j) + aux_tl%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j)) +  &
      & aux%s(prof)%pfraction_surf * (od_level_tl(levsurf - 1, j) - od_level_tl(levsurf, j))
    IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
      od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf - 1, j)
    ELSE
      od_frac_tl(j) = od_surf_tl(j) - od_level_tl(levsurf, j)
    ENDIF
! associated transmittance
    tau_surf_tl(j) = od_surf_tl(j) * tau_ref_surf(j)
  ENDDO
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tau_ref_surf(j) < coef%tt_a0(chan)) THEN
          tau_surf_tl(j) = 0._jprb
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!---Loop over the streams-----------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      prof    = chanprof(j)%prof
      levsurf = aux%s(prof)%nearestlev_surf
! layer above this
      DO ist = 0, ircld%nstream(prof)
        od_surf_ac_tl(j) = transmission_scatt_ir_stream_tl%opdpac(ist, j, levsurf) + aux%s(prof)%pfraction_surf * (     &
          & transmission_scatt_ir_stream_tl%opdpac(ist, j, levsurf - 1) -                                               &
          & transmission_scatt_ir_stream_tl%opdpac(ist, j, levsurf)) + aux_tl%s(prof)%pfraction_surf * (                &
          & transmission_scatt_ir_stream%opdpac(ist, j, levsurf - 1) - transmission_scatt_ir_stream%opdpac(ist, j, levsurf))
        IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
          transmission_aux_tl%od_frac_ac(ist, j) =      &
            & od_surf_ac_tl(j) - transmission_scatt_ir_stream_tl%opdpac(ist, j, levsurf - 1)
        ELSE
          transmission_aux_tl%od_frac_ac(ist, j) =      &
            & od_surf_ac_tl(j) - transmission_scatt_ir_stream_tl%opdpac(ist, j, levsurf)
        ENDIF
        transmission_aux_tl%tau_surf_ac(ist, j) =  - od_surf_ac_tl(j) * transmission_aux%tau_ref_surf_ac(ist, j)
        transmission_aux_tl%od_frac_t(ist, j)   =  - od_frac_tl(j) + transmission_aux_tl%od_frac_ac(ist, j)
        IF (tau_surf(j) >= 0) THEN
          transmission_aux_tl%tau_surf_t(ist, j) = tau_surf_tl(j) * transmission_aux%tau_ref_surf_ac(ist, j) +      &
            & tau_surf(j) * transmission_aux_tl%tau_surf_ac(ist, j)
        ELSE
          transmission_aux_tl%tau_surf_t(ist, j) = tau_surf_tl(j)
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!-----------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      prof = chanprof(j)%prof
      DO ist = 0, ircld%nstream(prof)
        DO lev = 1, nlevels
          IF (tau_level(lev, j) >= 0._jprb) THEN
            ztemp = exp( - transmission_scatt_ir_stream%opdpac(ist, j, lev))
            transmission_aux_tl%tau_level(lev, ist, j) = tau_level_tl(lev, j) * ztemp -      &
              & tau_level(lev, j) * transmission_scatt_ir_stream_tl%opdpac(ist, j, lev) * ztemp
          ELSE
            transmission_aux_tl%tau_level(lev, ist, j) = tau_level_tl(lev, j)
          ENDIF
        ENDDO
        DO lay = 1, nlayers
          transmission_aux_tl%od_singlelayer(lay, ist, j) =      &
            & od_singlelayer_tl(lay, j) + transmission_scatt_ir_stream_tl%opdpacl(ist, j, lay)
        ENDDO
      ENDDO
      DO ist = 0, ircld%nstream(prof)
        transmission_aux_tl%od_sfrac(ist, j) = transmission_aux_tl%od_frac_t(ist, j)
        transmission_aux_tl%tau_surf(ist, j) = transmission_aux_tl%tau_surf_t(ist, j)
      ENDDO
    ENDDO
  ELSE
    DO j = 1, nchannels
      prof = chanprof(j)%prof
      transmission_tl%tau_total(j)     = tau_surf_tl(j)
      transmission_tl%tau_levels(:, j) = tau_level_tl(:, j)
      DO ist = 0, ircld%nstream(prof)
        DO lev = 1, nlevels
          transmission_aux_tl%tau_level(lev, ist, j) = tau_level_tl(lev, j)
        ENDDO
        DO lay = 1, nlayers
          transmission_aux_tl%od_singlelayer(lay, ist, j) = od_singlelayer_tl(lay, j)
        ENDDO
      ENDDO
      DO ist = 0, ircld%nstream(prof)
        transmission_aux_tl%tau_surf(ist, j) = tau_surf_tl(j)
        transmission_aux_tl%od_sfrac(ist, j) =  - od_frac_tl(j)
      ENDDO
    ENDDO
  ENDIF
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_tl
