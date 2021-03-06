SUBROUTINE rttov_opdep_tl( &
            & nlayers,       &
            & chanprof,      &
            & predictors,    &
            & predictors_tl, &
            & aux,           &
            & aux_tl,        &
            & coef,          &
            & opdp_path,     &
            & opdp_path_tl,  &
            & opdp_ref)
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
!  1.0    01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.1    29/01/2003  Add WV Continuum and CO2 capability (P Brunel)
!  1.2    04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.3    26/09/2003  Modified to allow for multiple polarisations (S English)
!  1.4    28/02/2005  Improved vectorisation (D Dent)
!  1.5    01/06/2005  Marco Matricardi (ECMWF):
!            --       Computation ov variable OD_SFRAC added.
!  1.6    03/01/2007  Corrected bug in tau_surf_freq_tl assignment (R Saunders)
!  1.7    09/02/2007  Removed polarisation index
!  1.8    15/07/2009  User defined ToA. Layers distinct from levels (P.Rayer)
!  1.9    03/11/2009  Transmittances / optical depths on levels (A Geer)
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
       & rttov_chanprof,  &
       & rttov_coef,      &
       & predictors_Type, &
       & opdp_path_Type,  &
       & profile_aux
  USE parkind1, ONLY : jpim, jprb
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_mw
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof) , INTENT(IN)    :: chanprof(:)
  INTEGER(KIND=jpim)   , INTENT(IN)    :: nlayers
  TYPE(predictors_Type), INTENT(IN)    :: predictors
  TYPE(predictors_Type), INTENT(IN)    :: predictors_tl
  TYPE(rttov_coef     ), INTENT(IN)    :: coef
  TYPE(profile_aux    ), INTENT(IN)    :: aux
  TYPE(profile_aux    ), INTENT(IN)    :: aux_tl
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path                        ! optical depths
  TYPE(opdp_path_Type ), INTENT(INOUT) :: opdp_path_tl
  REAL(KIND=jprb)      , INTENT(IN)    :: opdp_ref(nlayers, size(chanprof))
!INTF_END
!local variables:
  REAL(KIND=jprb) :: opticaldepth_tl(nlayers, size(chanprof))
  REAL(KIND=jprb), POINTER :: debye_prof   (:, :)
  REAL(KIND=jprb), POINTER :: debye_prof_tl(:, :)
  INTEGER(KIND=jpim) :: lev, lay, chan      , j, nlevels
  INTEGER(KIND=jpim) :: prof
! cloud liquid water local variables
  REAL   (KIND=jprb) :: zf, zf_sq          , z34_dif   , z45_dif      , z1_sq     , z2_sq     , z1_div    , z2_div
  REAL   (KIND=jprb) :: z1_den         , z2_den         , zastar    , z1_prod      , z2_prod   , z3_prod   , z4_prod
  REAL   (KIND=jprb) :: zbstar         , zbstar_sq      , za2star   , za2star_sq   , zdiv      , zgstar
  REAL   (KIND=jprb) :: z1f_sq_z1_sq   , z2f_sq_z2_sq
  REAL   (KIND=jprb) :: z34_dif_tl     , z45_dif_tl     , z1_sq_tl  , z2_sq_tl     , z1_div_tl , z2_div_tl
  REAL   (KIND=jprb) :: z1_den_tl      , z2_den_tl      , zastar_tl , z1_prod_tl   , z2_prod_tl, z3_prod_tl, z4_prod_tl
  REAL   (KIND=jprb) :: zbstar_tl      , zbstar_sq_tl   , za2star_tl, za2star_sq_tl, zdiv_tl   , zgstar_tl
  REAL   (KIND=jprb) :: z1f_sq_z1_sq_tl, z2f_sq_z2_sq_tl, t1, t2, ztemp
  INTEGER(KIND=jpim) :: ii
  INTEGER(KIND=jpim) :: nchannels                                                                                      ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!- End of header --------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_TL', 0_jpim, ZHOOK_HANDLE)
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
!-----------------------------------------
!1. calculate layer gaseous optical depths
!-----------------------------------------
!--------------------------
!1.1 start with mixed gases
!--------------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    IF (coef%nmixed == 10) THEN
      DO lay = 1, nlayers
        t1 = coef%mixedgas(lay, chan, 1) * predictors_tl%mixedgas(1, lay, prof) +      &
          & coef%mixedgas(lay, chan, 2) * predictors_tl%mixedgas(2, lay, prof) +       &
          & coef%mixedgas(lay, chan, 3) * predictors_tl%mixedgas(3, lay, prof) +       &
          & coef%mixedgas(lay, chan, 4) * predictors_tl%mixedgas(4, lay, prof) +       &
          & coef%mixedgas(lay, chan, 5) * predictors_tl%mixedgas(5, lay, prof)
        t2 = coef%mixedgas(lay, chan, 6) * predictors_tl%mixedgas(6, lay, prof) +      &
          & coef%mixedgas(lay, chan, 7) * predictors_tl%mixedgas(7, lay, prof) +       &
          & coef%mixedgas(lay, chan, 8) * predictors_tl%mixedgas(8, lay, prof) +       &
          & coef%mixedgas(lay, chan, 9) * predictors_tl%mixedgas(9, lay, prof) +       &
          & coef%mixedgas(lay, chan, 10) * predictors_tl%mixedgas(10, lay, prof)
        opticaldepth_tl(lay, j) = t1 + t2
      ENDDO
    ELSE
      opticaldepth_tl(1:nlayers, j) = 0._JPRB
      DO ii = 1, coef%nmixed
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) =      &
            & opticaldepth_tl(lay, j) + coef%mixedgas(lay, chan, ii) * predictors_tl%mixedgas(ii, lay, prof)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
!--------------------
!1.2 add water vapour
!--------------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    IF (coef%nwater == 15) THEN
      DO lay = 1, nlayers
        t1 = coef%watervapour(lay, chan, 1) * predictors_tl%watervapour(1, lay, prof) +      &
          & coef%watervapour(lay, chan, 2) * predictors_tl%watervapour(2, lay, prof) +       &
          & coef%watervapour(lay, chan, 3) * predictors_tl%watervapour(3, lay, prof) +       &
          & coef%watervapour(lay, chan, 4) * predictors_tl%watervapour(4, lay, prof)
        t2 = coef%watervapour(lay, chan, 5) * predictors_tl%watervapour(5, lay, prof) +      &
          & coef%watervapour(lay, chan, 6) * predictors_tl%watervapour(6, lay, prof) +       &
          & coef%watervapour(lay, chan, 7) * predictors_tl%watervapour(7, lay, prof) +       &
          & coef%watervapour(lay, chan, 8) * predictors_tl%watervapour(8, lay, prof)
        opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + t1 + t2
      ENDDO
      DO lay = 1, nlayers
        t1 = coef%watervapour(lay, chan, 9) * predictors_tl%watervapour(9, lay, prof) +      &
          & coef%watervapour(lay, chan, 10) * predictors_tl%watervapour(10, lay, prof) +     &
          & coef%watervapour(lay, chan, 11) * predictors_tl%watervapour(11, lay, prof) +     &
          & coef%watervapour(lay, chan, 12) * predictors_tl%watervapour(12, lay, prof)
        t2 = coef%watervapour(lay, chan, 13) * predictors_tl%watervapour(13, lay, prof) +      &
          & coef%watervapour(lay, chan, 14) * predictors_tl%watervapour(14, lay, prof) +       &
          & coef%watervapour(lay, chan, 15) * predictors_tl%watervapour(15, lay, prof)
        opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + t1 + t2
      ENDDO
    ELSE
      DO ii = 1, coef%nwater
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) =      &
            & opticaldepth_tl(lay, j) + coef%watervapour(lay, chan, ii) * predictors_tl%watervapour(ii, lay, prof)
        ENDDO
      ENDDO
    ENDIF
  ENDDO
!-------------
!1.3 add ozone
!-------------
  DO j = 1, nchannels
    chan = chanprof(j)%chan
    prof = chanprof(j)%prof
    IF (coef%nozone > 0) THEN
      IF (coef%nozone == 11) THEN
        DO lay = 1, nlayers
          t1 = coef%ozone(lay, chan, 1) * predictors_tl%ozone(1, lay, prof) +      &
            & coef%ozone(lay, chan, 2) * predictors_tl%ozone(2, lay, prof) +       &
            & coef%ozone(lay, chan, 3) * predictors_tl%ozone(3, lay, prof) +       &
            & coef%ozone(lay, chan, 4) * predictors_tl%ozone(4, lay, prof) +       &
            & coef%ozone(lay, chan, 5) * predictors_tl%ozone(5, lay, prof)
          t2 = coef%ozone(lay, chan, 6) * predictors_tl%ozone(6, lay, prof) +      &
            & coef%ozone(lay, chan, 7) * predictors_tl%ozone(7, lay, prof) +       &
            & coef%ozone(lay, chan, 8) * predictors_tl%ozone(8, lay, prof) +       &
            & coef%ozone(lay, chan, 9) * predictors_tl%ozone(9, lay, prof) +       &
            & coef%ozone(lay, chan, 10) * predictors_tl%ozone(10, lay, prof) +     &
            & coef%ozone(lay, chan, 11) * predictors_tl%ozone(11, lay, prof)
          opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + t1 + t2
        ENDDO
      ELSE
        DO ii = 1, coef%nozone
          DO lay = 1, nlayers
            opticaldepth_tl(lay, j) =      &
              & opticaldepth_tl(lay, j) + coef%ozone(lay, chan, ii) * predictors_tl%ozone(ii, lay, prof)
          ENDDO
        ENDDO
      ENDIF
    ENDIF
  ENDDO
!1.4 add Water Vapour Continuum
!------------------------------
  IF (coef%nwvcont > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nwvcont
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) =      &
            & opticaldepth_tl(lay, j) + coef%wvcont(lay, chan, ii) * predictors_tl%wvcont(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!-----------
!1.5 add CO2
!-----------
  IF (coef%nco2 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nco2
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + coef%co2(lay, chan, ii) * predictors_tl%co2(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!-----------
!1.6 add N2O
!-----------
  IF (coef%nn2o > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nn2o
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + coef%n2o(lay, chan, ii) * predictors_tl%n2o(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!-----------
!1.7 add CO
!-----------
  IF (coef%nco > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nco
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + coef%co(lay, chan, ii) * predictors_tl%co(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!-----------
!1.8 add CH4
!-----------
  IF (coef%nch4 > 0) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      DO ii = 1, coef%nch4
        DO lay = 1, nlayers
          opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) + coef%ch4(lay, chan, ii) * predictors_tl%ch4(ii, lay, prof)
        ENDDO
      ENDDO
    ENDDO
  ENDIF
!--------------------
!1.9 add liquid water (MW only)
!--------------------
  IF (coef%id_sensor == sensor_id_mw) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      prof = chanprof(j)%prof
      debye_prof => aux%debye_prof(:, :, prof)
      debye_prof_tl => aux_tl%debye_prof(:, :, prof)
      IF (predictors%ncloud >= 1) THEN
        DO lay = 1, nlayers
          lev = lay + 1
          IF (lev >= coef%mwcldtop) THEN
! Repeat direct code
            zf = coef%frequency_ghz(chan)
            zf_sq = zf * zf
            z1_sq = debye_prof(1, lev) * debye_prof(1, lev)
            z2_sq = debye_prof(2, lev) * debye_prof(2, lev)
            z34_dif                 = debye_prof(3, lev) - debye_prof(4, lev)
            z45_dif                 = debye_prof(4, lev) - debye_prof(5, lev)
            z1f_sq_z1_sq            = zf_sq + z1_sq
            z2f_sq_z2_sq            = zf_sq + z2_sq
            z1_div                  = 1.0_JPRB / z1f_sq_z1_sq
            z2_div                  = 1.0_JPRB / z2f_sq_z2_sq
            z1_den                  = z34_dif * z1_div
            z2_den                  = z45_dif * z2_div
            zastar                  = debye_prof(3, lev) - zf_sq * (z1_den + z2_den)
            z1_prod                 = z34_dif * debye_prof(1, lev)
            z2_prod                 = z1_prod * z1_div
            z3_prod                 = z45_dif * debye_prof(2, lev)
            z4_prod                 = z3_prod * z2_div
            zbstar                  =  - zf * (z2_prod + z4_prod)
            zbstar_sq               = zbstar * zbstar
            za2star                 = zastar + 2.0_JPRB
            za2star_sq              = za2star * za2star
            zdiv = za2star_sq + zbstar_sq
            zgstar                  =  - 3.0_JPRB * zbstar / zdiv
! Now compute tangent-linear code
!zf_tl    = 0
!zf_sq_tl = 0
            z1_sq_tl                = 2.0_JPRB * debye_prof(1, lev) * debye_prof_tl(1, lev)
            z2_sq_tl                = 2.0_JPRB * debye_prof(2, lev) * debye_prof_tl(2, lev)
            z34_dif_tl              = debye_prof_tl(3, lev) - debye_prof_tl(4, lev)
            z45_dif_tl              = debye_prof_tl(4, lev) - debye_prof_tl(5, lev)
            z1f_sq_z1_sq_tl         = z1_sq_tl
            z2f_sq_z2_sq_tl         = z2_sq_tl
            z1_div_tl               =  - z1f_sq_z1_sq_tl / (z1f_sq_z1_sq * z1f_sq_z1_sq)
            z2_div_tl               =  - z2f_sq_z2_sq_tl / (z2f_sq_z2_sq * z2f_sq_z2_sq)
            z1_den_tl               = z34_dif * z1_div_tl + z34_dif_tl * z1_div
            z2_den_tl               = z45_dif * z2_div_tl + z45_dif_tl * z2_div
            zastar_tl               = debye_prof_tl(3, lev) - zf_sq * (z1_den_tl + z2_den_tl)
            z1_prod_tl              = z34_dif_tl * debye_prof(1, lev) + z34_dif * debye_prof_tl(1, lev)
            z2_prod_tl              = z1_prod_tl * z1_div + z1_prod * z1_div_tl
            z3_prod_tl              = z45_dif_tl * debye_prof(2, lev) + z45_dif * debye_prof_tl(2, lev)
            z4_prod_tl              = z3_prod_tl * z2_div + z3_prod * z2_div_tl
            zbstar_tl               =  - zf * (z2_prod_tl + z4_prod_tl)
            zbstar_sq_tl            = 2.0_JPRB * zbstar * zbstar_tl
            za2star_tl              = zastar_tl
            za2star_sq_tl           = 2.0_JPRB * za2star * za2star_tl
            zdiv_tl                 = za2star_sq_tl + zbstar_sq_tl
            zgstar_tl               =  - 3.0_JPRB * (zbstar_tl * zdiv - zbstar * zdiv_tl) / (zdiv * zdiv)
            opticaldepth_tl(lay, j) = opticaldepth_tl(lay, j) -      &
              & 1.5_JPRB * zf * (zgstar_tl * predictors%clw(lay, prof) + zgstar * predictors_tl%clw(lay, prof))
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!----------------------------------------
!2. Assemble layer optical depths
!----------------------------------------
! note that optical depth in the calculations is negative
! single layer
  WHERE (opdp_ref(:,:) > 0.0_JPRB)
    opticaldepth_tl(:,:) = 0.0_JPRB
  ENDWHERE
! level to space optical depths
  DO j = 1, nchannels
    opdp_path_tl%atm_level(1, j) = 0.0_jprb
!  Introduce ztemp to stop store followed by load in this recursive loop (DJS)
    ztemp = 0.0_jprb
    DO lev = 2, nlevels
      lay = lev - 1
      ztemp = ztemp + opticaldepth_tl(lay, j)
      opdp_path_tl%atm_level(lev, j) = ztemp
    ENDDO
  ENDDO
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDEP_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdep_tl
