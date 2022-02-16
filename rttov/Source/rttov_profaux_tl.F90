SUBROUTINE rttov_profaux_tl( &
            & opts,    &
            & prof,    &
            & prof_tl, &
            & coef,    &
            & aux,     &
            & aux_tl)
!
! Description:
! TL of rttov_profaux
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
!  1.0    07/10/2004  Added history
!  1.1    29/03/2005  Add end of header comment (J. Cameron)
!  1.2    04/02/2008  opts%lgradp option for TL/AD of pressure levels (N. Bormann)
!  1.3    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.4    02/12/2009  Fixed a number of bugs due to the wrong assumption that cloud
!                     related quantities are defined on levels (thay are layer
!                     average quantities). Marco Matricardi
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
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
       & rttov_options, &
       & rttov_coef,    &
       & profile_Type,  &
       & profile_aux
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw, sensor_id_po, dcoeff
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jpim, jprb
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_options), INTENT(IN)    :: opts
  TYPE(profile_Type ), INTENT(IN)    :: prof   (:)
  TYPE(profile_Type ), INTENT(IN)    :: prof_tl(:)
  TYPE(rttov_coef   ), INTENT(IN)    :: coef
  TYPE(profile_aux  ), INTENT(IN)    :: aux
  TYPE(profile_aux  ), INTENT(INOUT) :: aux_tl
!INTF_END
! local
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp   (prof(1)%nlayers)
  REAL   (KIND=jprb) :: dp_tl(prof(1)%nlayers)
  REAL   (KIND=jprb) :: v                                                       ! temperature ratio
  REAL   (KIND=jprb) :: v_tl                                                    ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq, zmcfarq_tl
  REAL   (KIND=jprb) :: ztempc      , ztempc_tl
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , bwyser_tl   , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!-----------------------------------------
! TL for cloud top and surface levels
!-----------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  DO iprof = 1, nprofiles
    dp(1:prof(iprof)%nlayers) = prof(iprof)%p(2:prof(iprof)%nlevels) - prof(iprof)%p(1:prof(iprof)%nlevels - 1)
    IF (opts%lgradp) THEN
      dp_tl(1:prof(iprof)%nlayers) = prof_tl(iprof)%p(2:prof(iprof)%nlevels) - &
                                   & prof_tl(iprof)%p(1:prof(iprof)%nlevels - 1)
    ENDIF
!nearest level above surface
    IF (opts%lgradp) THEN
      aux_tl%s(iprof)%pfraction_surf =                                                                                      &
        & (prof_tl(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof_tl(iprof)%s2m%p) / dp(aux%s(iprof)%nearestlev_surf - 1) -  &
        & (prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%s2m%p) / dp(aux%s(iprof)%nearestlev_surf - 1) ** 2 *   &
        & dp_tl(aux%s(iprof)%nearestlev_surf - 1)
    ELSE
      aux_tl%s(iprof)%pfraction_surf =  - prof_tl(iprof)%s2m%p / dp(aux%s(iprof)%nearestlev_surf - 1)
    ENDIF
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
!nearest level above cloud top
      IF (opts%lgradp) THEN
        aux_tl%s(iprof)%pfraction_ctp =                                                                                   &
          & (prof_tl(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof_tl(iprof)%ctp) / dp(aux%s(iprof)%nearestlev_ctp - 1) -  &
          & (prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%ctp) / dp(aux%s(iprof)%nearestlev_ctp - 1) ** 2 *   &
          & dp_tl(aux%s(iprof)%nearestlev_ctp - 1)
      ELSE
        aux_tl%s(iprof)%pfraction_ctp =  - prof_tl(iprof)%ctp / dp(aux%s(iprof)%nearestlev_ctp - 1)
      ENDIF
      aux_tl%s(iprof)%cfraction = prof_tl(iprof)%cfraction
    ELSE
! for micro waves do not consider clouds in the RTTOV basis routines
      aux_tl%s(iprof)%pfraction_ctp = 0._JPRB
      aux_tl%s(iprof)%cfraction     = 0._JPRB
    ENDIF
!--------------------------------------------------
! TL for debye terms for pewrmittivity calculations
!--------------------------------------------------
! Description:
!   To calculate individual debye terms for temperature
!   at each level. There are five debye terms. These
!   will be used in fastem and opdep to calculate
!   permittivity which is required for surface emissivity
!   and cloud modelling
!
! Method:
!   The model is a hybrid of LIEBE MPM 1993 and the PIOM laboratory
!   measurements reported by ELLISON et al. 1999.
  IF (coef%id_sensor == sensor_id_mw .OR. coef%id_sensor == sensor_id_po) THEN
    IF (opts%clw_Data) THEN
      DO lev = coef%mwcldtop, prof(iprof)%nlevels
        v = 300.0_JPRB / prof(iprof)%t(lev) - 1.0_JPRB
        v_tl =  - 300.0_JPRB * prof_tl(iprof)%t(lev) / prof(iprof)%t(lev) ** 2
        aux_tl%debye_prof(1, lev, iprof) =  - dcoeff(2) * v_tl + 2 * dcoeff(3) * v_tl * v
        aux_tl%debye_prof(2, lev, iprof) = dcoeff(4) * aux_tl%debye_prof(1, lev, iprof)
        aux_tl%debye_prof(3, lev, iprof) = dcoeff(5) * v_tl
        aux_tl%debye_prof(4, lev, iprof) = dcoeff(7) * aux_tl%debye_prof(3, lev, iprof)
        aux_tl%debye_prof(5, lev, iprof) = 0._JPRB
      ENDDO
    ENDIF
  ENDIF
!--------------------------------------------------------------
! TL for ice  water content into effective generalized diameter
!--------------------------------------------------------------
    IF (opts%addclouds) THEN
! Calculate upper and lower limits for Ou-Liou effective size
!
      zradipou_upp = 326.3_JPRB + rtou_upp * (12.42_JPRB + rtou_upp * (0.197_JPRB + rtou_upp * 0.0012_JPRB))
      zradipou_low = 326.3_JPRB + rtou_low * (12.42_JPRB + rtou_low * (0.197_JPRB + rtou_low * 0.0012_JPRB))
!
! and convert these to the "generalized" effective size used here (using McFarquhar et al 2003 equation),
! not forgetting the factor of 2 to convert from McFarquhar's radius to a diameter
!
      zradipou_upp =  - 1.56_JPRB + zradipou_upp * (0.388_JPRB + zradipou_upp * 0.00051_JPRB)
      zradipou_upp = 2.0_JPRB * zradipou_upp
      zradipou_low =  - 1.56_JPRB + zradipou_low * (0.388_JPRB + zradipou_low * 0.00051_JPRB)
      zradipou_low = 2.0_JPRB * zradipou_low
      DO lay = 1, prof(iprof)%nlayers
        lev = lay + 1
        IF (prof(iprof)%cloud(6, lay) /= 0._jprb) THEN
          ! Use effective diameter from input profile if specified
          IF (prof(iprof)%icede(lay) /= 0._jprb) THEN
            aux_tl%dg(lay, iprof) = prof_tl(iprof)%icede(lay)
          ELSE
            IF (prof(iprof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc = prof(iprof)%t(lev) - rtt
              ztempc_tl                  = prof_tl(iprof)%t(lev)
              aux_tl%fac1_dg(lay, iprof) =      &
                & ztempc_tl * (12.42_JPRB + 2._JPRB * 0.197_JPRB * ztempc + 3._JPRB * 0.0012_JPRB * ztempc * ztempc)
              aux_tl%fac2_dg(lay, iprof) = 0.388_JPRB * aux_tl%fac1_dg(lay, iprof) +      &
                & 0.00051_JPRB * 2 * aux%fac1_dg(lay, iprof) * aux_tl%fac1_dg(lay, iprof)
              aux_tl%fac3_dg(lay, iprof) = 2.0_JPRB * aux_tl%fac2_dg(lay, iprof)
  !
  ! Take Ou-Liou scheme as being valid only between 20C and 60C
  !
              IF (aux%fac3_dg(lay, iprof) < zradipou_low .OR. aux%fac3_dg(lay, iprof) > zradipou_upp) THEN
                aux_tl%dg(lay, iprof) = 0._JPRB
              ELSE
                aux_tl%dg(lay, iprof) = aux_tl%fac3_dg(lay, iprof)
              ENDIF
            ELSE IF (prof(iprof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser    =  - 2.0_JPRB
              bwyser_tl = 0._JPRB
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                bwyser    = bwyser + (0.001_JPRB * &
                  ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * Log10(prof(iprof)%cloud(6, lay) / 50._JPRB))
                bwyser_tl =  -                                                                                        &
                  & 0.001_JPRB * 1.5_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 0.5_JPRB) * prof_tl(iprof)%t(lev) *  &
                  & Log10(prof(iprof)%cloud(6, lay) / 50._JPRB) +                                                     &
                  & 0.001_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * prof_tl(iprof)%cloud(6, lay) /      &
                  & prof(iprof)%cloud(6, lay) / log(10._JPRB)
              ENDIF
              aux_tl%fac1_dg(lay, iprof) =      &
                & (203.3_JPRB + 2 * bwyser * 37.91_JPRB + 3 * bwyser * bwyser * 2.3696_JPRB) * bwyser_tl
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux_tl%fac2_dg(lay, iprof) = aux_tl%fac1_dg(lay, iprof) / nft
              aux_tl%dg(lay, iprof)      = 2._JPRB * 4._JPRB * aux_tl%fac2_dg(lay, iprof) * sqrt(3._JPRB) / 9._JPRB
            ELSE IF (prof(iprof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc                = prof(iprof)%t(lev) - rtt
              ztempc_tl             = prof_tl(iprof)%t(lev)
              aux_tl%dg(lay, iprof) = 53.005_JPRB * 0.06_JPRB *                                 &
                (prof(iprof)%cloud(6, lay) ** (0.06_JPRB - 1)) * prof_tl(iprof)%cloud(6, lay) * &
                exp(0.013_JPRB * ztempc) + 53.005_JPRB *                                        &
                ((prof(iprof)%cloud(6, lay)) ** 0.06_JPRB) * 0.013_JPRB * ztempc_tl * exp(0.013_JPRB * ztempc)
            ELSE IF (prof(iprof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                    = 1.78449_JPRB
              bmcfarq                    = 0.281301_JPRB
              cmcfarq                    = 0.0177166_JPRB
              zmcfarq                    = prof(iprof)%cloud(6, lay)
              zmcfarq_tl                 = prof_tl(iprof)%cloud(6, lay)
              aux_tl%fac1_dg(lay, iprof) =                                                                         &
                & 10.0_JPRB ** (amcfarq + bmcfarq * Log10(zmcfarq) + cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)) *  &
                & (bmcfarq + 2._JPRB * cmcfarq * Log10(zmcfarq)) / zmcfarq * zmcfarq_tl
              aux_tl%dg(lay, iprof)      = 2.0_JPRB * aux_tl%fac1_dg(lay, iprof)
            ELSE
              aux_tl%dg(lay, iprof) = 0._jprb
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux_tl
