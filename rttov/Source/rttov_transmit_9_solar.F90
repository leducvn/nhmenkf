
SUBROUTINE rttov_transmit_9_solar( &
            & addaerosl,                    &
            & addclouds,                    &
            & nlayers,                      &
            & chanprof,                     &
            & profiles,                     &
            & sun,                          &
            & aux,                          &
            & coef,                         &
            & raytracing,                   &
            & ircld,                        &
            & opdp_path,                    &
            & odsun_level,                  &
            & odsun_singlelayer,            &
            & od_frac,                      &
            & transmission_aux,             &
            & transmission_scatt_ir_stream, &
            & tausun_ref,                   &
            & tausun_ref_surf,              &
            & tausun_surf,                  &
            & tausun_level)
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
!  1.2    28/08/2007  Optimised (D. Salmond)
!  1.3    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - size, index
!                     labels, looping (P. Rayer)
!  1.4    03/11/2009  Transmittances / optical depths on levels (A Geer)
!  1.5    02/12/2009  Pathsat, Pathsun and related quantities are now
!                     layer arrays (Marco Matricardi).
!  1.6    05/07/2010  Remove addsolar flag from profiles structure (J Hocking)
!  1.7    04/08/2010  Move addsolar check to calling routine (J Hocking)
!  1.8    14/12/2010  Use traj0_sta%sun array to flag channels for which solar calculations
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
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type,                 &
       & raytracing_type,            &
       & profile_Type
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : max_optical_depth, max_sol_zen
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers                                          ! Number of pressure levels
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)                                      ! Channel indices
  TYPE(profile_Type              ), INTENT(IN)    :: profiles(:)                                      ! Atmospheric profiles
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(rttov_coef                ), INTENT(IN)    :: coef                                             ! Coefficients
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux                                 ! Transmittances and single-layer od
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(profile_aux               ), INTENT(IN)    :: aux                                              ! auxillary profiles informations
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  TYPE(raytracing_type           ), INTENT(IN)    :: raytracing
  REAL(KIND=jprb)                 , INTENT(OUT)   :: odsun_level      (nlayers + 1   , size(chanprof))! sat to level optical depth
  REAL(KIND=jprb)                 , INTENT(OUT)   :: odsun_singlelayer(nlayers       , size(chanprof))! single layer optical depth
  REAL(KIND=jprb)                 , INTENT(OUT)   :: od_frac          (size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_ref_surf  (size(chanprof)                )! sat to surface transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_ref       (nlayers + 1   , size(chanprof))! sat to level transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_level     (nlayers + 1   , size(chanprof))! sat to level transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tausun_surf      (size(chanprof)                )! sat to surface transmittance
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_surf_ac    (size(chanprof)                )
  REAL   (KIND=jprb) :: od_surf       (size(chanprof)                )   ! sat to surface optical depth
  REAL   (KIND=jprb) :: od_level      (nlayers + 1   , size(chanprof))   ! sat to level optical depth
  REAL   (KIND=jprb) :: od_singlelayer(nlayers       , size(chanprof))   ! single layer optical depth
  INTEGER(KIND=jpim) :: lev         , lay    , chan   , j, prof, ist     ! loop variables
  INTEGER(KIND=jpim) :: nlevels
! cloud liquid water local variables
  INTEGER(KIND=jpim) :: levsurf
  INTEGER(KIND=jpim) :: nchannels                                        ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!--------------------------------------------------------------
!1. Assemble layer optical depths and convert to transmittances
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that values are sensible
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lay = 1, nlayers
        od_singlelayer(lay, j) =  - (opdp_path%sun_level(lay + 1, j) - opdp_path%sun_level(lay, j))
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level(lev, j) = opdp_path%sun_level(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level(lev, j) = coef%ff_gam(chan) * od_level(lev, j)
      ENDDO
      DO lay = 1, nlayers
        od_singlelayer(lay, j) = coef%ff_gam(chan) * od_singlelayer(lay, j)
      ENDDO
    ENDIF
  ENDDO
! On some computers when optical depth is too thick
! there is an underflow during the conversion in
! transmittances. In that case uncomment following line
! and the declaration statement of max_optical_depth
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        od_level(lev, j) = Max(od_level(lev, j),  - max_optical_depth)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    IF (sun(j)) THEN
      DO lev = 1, nlevels
        tausun_level(lev, j) = Exp(od_level(lev, j))
        tausun_ref(lev, j)   = tausun_level(lev, j)
      ENDDO
    ENDIF
  ENDDO
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    IF (sun(j)) THEN
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tausun_level(lev, j) < coef%tt_a0(chan)) THEN
            tausun_level(lev, j) = coef%tt_a1(chan)
          ENDIF
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    chan    = chanprof(j)%chan
    prof    = chanprof(j)%prof
! as defined in rttov_profaux
    levsurf = aux%s(prof)%nearestlev_surf
! layer above this
    IF (sun(j)) THEN
      od_surf(j) =      &
        & od_level(levsurf, j) + aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j))
      IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
        od_frac(J) = od_surf(j) - od_level(levsurf - 1, j)
      ELSE
        od_frac(J) = od_surf(j) - od_level(levsurf, j)
      ENDIF
      tausun_surf(j)     = Exp(od_surf(j))
      tausun_ref_surf(j) = tausun_surf(j)
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tausun_surf(j) < coef%tt_a0(chan)) THEN
          tausun_surf(j) = coef%tt_a1(chan)
        ENDIF
      ENDIF
!---Loop over the streams-----------------------------------------------------------------
      IF (addaerosl .OR. addclouds) THEN
        DO ist = 0, ircld%nstream(prof)
          od_surf_ac(j) = transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf) + aux%s(prof)%pfraction_surf * (     &
            & transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf - 1) -                                            &
            & transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf))
          IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
            transmission_aux%odsun_frac_ac(ist, J) =      &
              & od_surf_ac(j) - transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf - 1)
          ELSE
            transmission_aux%odsun_frac_ac(ist, J) =      &
              & od_surf_ac(j) - transmission_scatt_ir_stream%opdpacsun(ist, j, levsurf)
          ENDIF
          transmission_aux%tau_surf_acsun(ist, j)     = Exp( - od_surf_ac(j))
          transmission_aux%tau_ref_surf_acsun(ist, j) = transmission_aux%tau_surf_acsun(ist, j)
          transmission_aux%odsun_frac_t(ist, J)       =                                                           &
            & ( - od_frac(j) + transmission_aux%odsun_frac_ac(ist, J)) * raytracing%pathsun(levsurf - 1, prof) /  &
            & (raytracing%pathsun(levsurf - 1, prof) + raytracing%pathsat(levsurf - 1, prof))
          IF (tausun_surf(j) >= 0) THEN
            transmission_aux%tausun_surf_t(ist, j) = tausun_surf(j) * transmission_aux%tau_surf_acsun(ist, J)
          ELSE
            transmission_aux%tausun_surf_t(ist, j) = tausun_surf(j)
          ENDIF
          transmission_aux%tausun_ref_surf_t(ist, j) = transmission_aux%tausun_surf_t(ist, j)
        ENDDO
      ENDIF
!-----------------------------------------------------------------------------------------
    ENDIF
  ENDDO
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof = chanprof(j)%prof! Profile index
    IF (sun(j)) THEN
      DO ist = 0, ircld%nstream(prof)
        transmission_aux%tausun_level(:, ist, j)      = 1.0_JPRB
        transmission_aux%tausun_surf(ist, j)          = 1.0_JPRB
        transmission_aux%odsun_singlelayer(:, ist, j) = 1.0_JPRB
        transmission_aux%odsun_sfrac(ist, j)          = 1.0_JPRB
      ENDDO
      odsun_level(:, j)       = 0.0_JPRB
      odsun_singlelayer(:, j) = 0.0_JPRB
      DO ist = 0, ircld%nstream(prof)
        IF (addaerosl .OR. addclouds) THEN
          DO lev = 1, nlevels
            IF (tausun_level(lev, j) >= 0) THEN
              transmission_aux%tausun_level(lev, ist, j) =      &
                & tausun_level(lev, j) * exp( - transmission_scatt_ir_stream%opdpacsun(ist, j, lev))
            ELSE
              transmission_aux%tausun_level(lev, ist, j) = tausun_level(lev, j)
            ENDIF
          ENDDO
          transmission_aux%tausun_surf(ist, j) = transmission_aux%tausun_surf_t(ist, j)
          transmission_aux%odsun_sfrac(ist, j) = transmission_aux%odsun_frac_t(ist, J)
          DO lay = 1, nlayers
            lev = lay + 1
            transmission_aux%odsun_singlelayer(lay, ist, j)   =                                    &
              & (od_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpaclsun(ist, j, lay)) *  &
              & raytracing%pathsun(lay, prof) / (raytracing%pathsun(lay, prof) + raytracing%pathsat(lay, prof))
            transmission_scatt_ir_stream%opdpext(ist, j, lay) =                               &
              & od_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpabs(ist, j, lay) +  &
              & transmission_scatt_ir_stream%opdpsca(ist, j, lay)
            IF (transmission_scatt_ir_stream%opdpext(ist, j, lay) /= 0._jprb) THEN
              transmission_scatt_ir_stream%ssa(ist, j, lay) =      &
                & transmission_scatt_ir_stream%opdpsca(ist, j, lay) / transmission_scatt_ir_stream%opdpext(ist, j, lay)
            ENDIF
          ENDDO
          odsun_level(:, j)       = od_level(:, j)
          odsun_singlelayer(:, j) = od_singlelayer(:, j)
        ELSE
          transmission_aux%tausun_level(:, ist, j) = tausun_level(:, j)
          transmission_aux%tausun_surf(ist, j)     = tausun_surf(j)
          odsun_level(:, j)                        = od_level(:, j)
          odsun_singlelayer(:, j)                  = od_singlelayer(:, j)
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT_9_SOLAR', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit_9_solar
