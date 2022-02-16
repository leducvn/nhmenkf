!
SUBROUTINE rttov_profaux( &
            & opts, &
            & prof, &
            & coef, &
            & aux)
!
! Description:
! Calculates some variables related to the input profile.
! variables are nearest surface level, nearest cloud top level
! and Debye terms for MW
! The reason of having a separate structure for these
! variables is that the input profiles should be "read only"
! in RTTOV context.
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
!  1.0       01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1       27/02/2009  Profile levels to include ToA. Distinguish between
!                        layer arrays and level arrays - size, index
!                        labels, looping (P. Rayer)
!  1.2       02/12/2009  Fixed a number of bugs due to the wrong assumption that cloud
!                        related quantities are defined on levels (thay are layer
!                        average quantities). Marco Matricardi
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
! debye coefficients
! Imported Type Definitions:
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
  TYPE(profile_Type ), INTENT(IN)    :: prof(:)! profile
  TYPE(rttov_coef   ), INTENT(IN)    :: coef   ! coefficients
  TYPE(profile_aux  ), INTENT(INOUT) :: aux    ! auxilary profile info
!INTF_END
!local
  INTEGER(KIND=jpim) :: nprofiles
  INTEGER(KIND=jpim) :: iprof
  INTEGER(KIND=jpim) :: lev      , lay
  REAL   (KIND=jprb) :: dp(prof(1)%nlayers)
  REAL   (KIND=jprb) :: v                                           ! temperature ratio
  REAL   (KIND=jprb) :: amcfarq     , bmcfarq     , cmcfarq, zmcfarq
  REAL   (KIND=jprb) :: ztempc
  REAL   (KIND=jprb) :: zradipou_upp, zradipou_low
  REAL   (KIND=jprb) :: bwyser      , nft
  REAL(KIND=jprb), PARAMETER :: rtt      = 273.15_JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_upp =  - 20._JPRB
  REAL(KIND=jprb), PARAMETER :: rtou_low =  - 60._JPRB
  REAL(KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
!-----------------------------------------
! determine cloud top and surface levels
!-----------------------------------------
! in line with coef % dp in rttov_initcoeffs
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX', 0_jpim, ZHOOK_HANDLE)
  nprofiles = size(prof)
  DO iprof = 1, nprofiles
    dp(1:prof(iprof)%nlayers) = prof(iprof)%p(2:prof(iprof)%nlevels) - prof(iprof)%p(1:prof(iprof)%nlevels - 1)
! nearest level above surface
    DO lev = prof(iprof)%nlevels - 1, 1,  - 1
      IF (prof(iprof)%s2m%p > prof(iprof)%p(lev)) EXIT
    ENDDO
! case-1: surf lies above lev=nlevels
!         at exit, lev is first level above surface
! case-2: surf lies below lev=nlevels
!         at exit, lev+1=nlevels, there is no level below surface
! case-1: first level below surface
! case-2: first level above surface
    aux%s(iprof)%nearestlev_surf = lev + 1
    aux%s(iprof)%pfraction_surf  =      &
      & (prof(iprof)%p(aux%s(iprof)%nearestlev_surf) - prof(iprof)%s2m%p) / dp(aux%s(iprof)%nearestlev_surf - 1)
! NB for case-2, aux % s(iprof) % pfraction_surf -ve
!nearest level above cloud top
    IF (coef%id_sensor /= sensor_id_mw .AND. coef%id_sensor /= sensor_id_po) THEN
      DO lev = prof(iprof)%nlevels - 1, 1,  - 1
        IF (prof(iprof)%ctp > prof(iprof)%p(lev)) EXIT
      ENDDO
      IF( lev > 1 ) THEN
        aux%s(iprof)%nearestlev_ctp = lev + 1
        aux%s(iprof)%pfraction_ctp  =      &
          & (prof(iprof)%p(aux%s(iprof)%nearestlev_ctp) - prof(iprof)%ctp) / dp(aux%s(iprof)%nearestlev_ctp - 1)
        aux%s(iprof)%cfraction      = prof(iprof)%cfraction
      ELSE
        aux%s(iprof)%nearestlev_ctp = prof(iprof)%nlevels -1
        aux%s(iprof)%pfraction_ctp  = 0._JPRB
        aux%s(iprof)%cfraction      = 0._JPRB
      ENDIF
    ELSE
! for micro waves do not consider clouds in the RTTOV basis routines
      aux%s(iprof)%nearestlev_ctp = prof(iprof)%nlevels - 1
      aux%s(iprof)%pfraction_ctp  = 0._JPRB
      aux%s(iprof)%cfraction      = 0._JPRB
    ENDIF
!---------------------------------------------
! debye terms for pewrmittivity calculations
!---------------------------------------------
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
      aux%debye_prof(:, :, iprof) = 0._jprb
      IF (opts%clw_Data) THEN
! cycle from from level corresponding to the parameter mwcldtp
        DO lev = coef%mwcldtop, prof(iprof)%nlevels
          v = 300.0_JPRB / prof(iprof)%t(lev) - 1.0_JPRB
          aux%debye_prof(1, lev, iprof) = dcoeff(1) - dcoeff(2) * v + dcoeff(3) * v * v
          aux%debye_prof(2, lev, iprof) = dcoeff(4) * aux%debye_prof(1, lev, iprof)
          aux%debye_prof(3, lev, iprof) = dcoeff(5) * v + dcoeff(6)
          aux%debye_prof(4, lev, iprof) = dcoeff(7) * aux%debye_prof(3, lev, iprof)
          aux%debye_prof(5, lev, iprof) = dcoeff(8)
        ENDDO
      ENDIF
    ENDIF
!-----------------------------------------------------------------------------------------
! If ice clouds, convert ice water content into effective generalized diameter
!-----------------------------------------------------------------------------------------
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
            aux%dg(lay, iprof) = prof(iprof)%icede(lay)
          ELSE
            IF (prof(iprof)%idg == 1) THEN
  !Scheme by Ou and Liou, 1995, Atmos. Res., 35, 127-138.
              ztempc                  = prof(iprof)%t(lev) - rtt
  ! intermediate factors in calculating the generalized effective diameter
              aux%fac1_dg(lay, iprof) = 326.3_JPRB + ztempc * (12.42_JPRB + ztempc * (0.197_JPRB + ztempc * 0.0012_JPRB))
              aux%fac2_dg(lay, iprof) =      &
                &  - 1.56_JPRB + aux%fac1_dg(lay, iprof) * (0.388_JPRB + aux%fac1_dg(lay, iprof) * 0.00051_JPRB)
              aux%fac3_dg(lay, iprof) = 2.0_JPRB * aux%fac2_dg(lay, iprof)
  !
  ! Take Ou-Liou scheme as being valid only between -20C and -60C
  !
              aux%dg(lay, iprof)      = max(aux%fac3_dg(lay, iprof), zradipou_low)
              aux%dg(lay, iprof)      = min(aux%dg(lay, iprof), zradipou_upp)
            ELSE IF (prof(iprof)%idg == 2) THEN
  !Scheme by Wyser et al. (see McFarquhar et al. (2003))
              bwyser =  - 2.0_JPRB
              IF (prof(iprof)%t(lev) < 273._JPRB) THEN
                bwyser = bwyser + (0.001_JPRB * ((273._JPRB - prof(iprof)%t(lev)) ** 1.5_JPRB) * &
                  Log10(prof(iprof)%cloud(6, lay) / 50._JPRB))
              ENDIF
              aux%fac1_dg(lay, iprof) = 377.4_JPRB + bwyser * (203.3_JPRB + bwyser * (37.91_JPRB + bwyser * 2.3696_JPRB))
              nft = (sqrt(3._JPRB) + 4._JPRB) / (3._JPRB * sqrt(3._JPRB))
              aux%fac2_dg(lay, iprof) = aux%fac1_dg(lay, iprof) / nft
              aux%dg(lay, iprof)      = 2._JPRB * 4._JPRB * aux%fac2_dg(lay, iprof) * sqrt(3._JPRB) / 9._JPRB
            ELSE IF (prof(iprof)%idg == 3) THEN
  !Scheme by Boudala et al., 2002, Int. J. Climatol., 22, 1267-1284.
              ztempc             = prof(iprof)%t(lev) - rtt
              aux%dg(lay, iprof) = 53.005_JPRB * ((prof(iprof)%cloud(6, lay)) ** 0.06_JPRB) * exp(0.013_JPRB * ztempc)
            ELSE IF (prof(iprof)%idg == 4) THEN
  ! Scheme by McFarquhar et al. (2003)
              amcfarq                 = 1.78449_JPRB
              bmcfarq                 = 0.281301_JPRB
              cmcfarq                 = 0.0177166_JPRB
              zmcfarq                 = prof(iprof)%cloud(6, lay)
              aux%fac1_dg(lay, iprof) =      &
                & 10.0_JPRB ** (amcfarq + (bmcfarq * Log10(zmcfarq)) + (cmcfarq * Log10(zmcfarq) * Log10(zmcfarq)))
              aux%dg(lay, iprof)      = 2.0_JPRB * aux%fac1_dg(lay, iprof)
            ELSE
              aux%dg(lay, iprof) = 0._jprb
            ENDIF
          ENDIF
        ENDIF
      ENDDO
    ENDIF
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_PROFAUX', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_profaux
