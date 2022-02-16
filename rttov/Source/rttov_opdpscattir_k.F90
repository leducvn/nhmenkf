!    Set up aerosols optical parameters for a climatological profile
SUBROUTINE rttov_opdpscattir_k( &
            & nlayers,                        &
            & chanprof,                       &
            & opts,                           &
            & aux,                            &
            & aux_k,                          &
            & profiles,                       &
            & profiles_k,                     &
            & sun,                            &
            & coef,                           &
            & coef_scatt_ir,                  &
            & raytracing,                     &
            & raytracing_k,                   &
            & transmission_scatt_ir,          &
            & transmission_scatt_ir_k,        &
            & transmission_scatt_ir_stream,   &
            & transmission_scatt_ir_stream_k, &
            & optp,                           &
            & ircld,                          &
            & ircld_k)
!     Description:
!     To set up profile-dependent variables for subsequent
!     rt calculations by other subroutines of RTIASI.
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
! Current Code Owner: SAF NWP
!
!     Method:
!
!     Owner:
!     EUMETSAT
!
!     History:
!     Version      Date        Comment
!     1           01/6/2004    Marco Matricardi. ECMWF.
!     1.1         05/02/2007   Removed polarisation R Saunders
!     1.2         15/09/2009   User defined ToA. Layers distinct from levels (P.Rayer)
!     1.3         03/11/2009   Transmittances / optical depths on levels (A Geer)
!     1.4         02/12/2009   Fixed a number of bugs due to the wrong assumption that aerosol/cloud
!                              related quantities are defined on levels (thay are layer
!                              average quantities). Marco Matricardi
!     1.5         02/12/2009   Introduced multiple cloud types in a single layer. Pathsun, Pathsat and
!                              related quantities are now layer arrays (Marco Matricardi)
!     1.6         05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!     1.7         02/08/2010   Allow variable pressure levels (lgradp) (N Bormann)
!     1.8         14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                              should be performed (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
!     Code description:
!       Language:              Fortran 90.
!       Software Standards:    "European Standards for Writing and Documenting
!                              Exchangeable Fortran 90 code".
!
!     Module used:
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_options,              &
       & rttov_coef,                 &
       & profile_Type,               &
       & raytracing_type,            &
       & transmission_scatt_ir_type, &
       & rttov_coef_scatt_ir,        &
       & rttov_optpar_ir,            &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jplm
!INTF_OFF
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE parkind1, ONLY : jprb
  USE rttov_const, ONLY :  &
       & pi,      &
       & deg2rad, &
       & jpazn,   &
       & t00,     &
       & ti,      &
       & ncldtyp, &
       & max_sol_zen
!INTF_ON
  IMPLICIT NONE
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof  (:)
  TYPE(rttov_options )            , INTENT(IN)    :: opts
  TYPE(profile_type              ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_type              ), INTENT(INOUT) :: profiles_k(size(chanprof))
  LOGICAL(KIND=jplm)              , INTENT(IN)    :: sun(size(chanprof))
  TYPE(profile_aux               ), INTENT(IN)    :: aux
  TYPE(profile_aux               ), INTENT(INOUT) :: aux_k
  TYPE(rttov_coef                ), INTENT(IN)    :: coef
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_k
  TYPE(transmission_scatt_ir_type), INTENT(IN)    :: transmission_scatt_ir_stream
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream_k
  TYPE(rttov_coef_scatt_ir       ), INTENT(IN)    :: coef_scatt_ir
  TYPE(rttov_optpar_ir           ), INTENT(IN)    :: optp
  TYPE(ircld_type),                 INTENT(IN)    :: ircld
  TYPE(ircld_type),                 INTENT(INOUT) :: ircld_k
  TYPE(raytracing_type), INTENT(IN)    :: raytracing
  TYPE(raytracing_type), INTENT(INOUT) :: raytracing_k
!INTF_END
!     End of subroutine arguments
!       Local scalars:
  INTEGER(KIND=jpim) :: j, ich    , ich1 , i, ipf, ish, k, kk, ityp, iconf, iae
  INTEGER(KIND=jpim) :: lev         , lay   , lctyp
  INTEGER(KIND=jpim) :: iang        , iend  , ist
  REAL   (KIND=jprb) :: opd         , opdsun
  REAL   (KIND=jprb) :: absch
  REAL   (KIND=jprb) :: scach
  REAL   (KIND=jprb) :: bparh
  REAL   (KIND=jprb) :: afac
  REAL   (KIND=jprb) :: sfac
  REAL   (KIND=jprb) :: gfac
  REAL   (KIND=jprb) :: scattangle
  REAL   (KIND=jprb) :: deltapice_d
  REAL   (KIND=jprb) :: deltapice_k
  REAL   (KIND=jprb) :: deltapaer_k
  REAL   (KIND=jprb) :: deltapaer_d
  REAL   (KIND=jprb) :: deltapwcl_d
  REAL   (KIND=jprb) :: deltaaer
  REAL   (KIND=jprb) :: deltaice
  REAL   (KIND=jprb) :: deltawcl
  REAL   (KIND=jprb) :: delth
  REAL   (KIND=jprb) :: frach
  REAL   (KIND=jprb) :: phasint
  REAL   (KIND=jprb) :: musat
  REAL   (KIND=jprb) :: musun
  REAL   (KIND=jprb) :: phup
  REAL   (KIND=jprb) :: phdo
  REAL   (KIND=jprb) :: scattangle_d
  REAL   (KIND=jprb) :: musat_d
  REAL   (KIND=jprb) :: musun_d
  REAL   (KIND=jprb) :: frach_d
  REAL   (KIND=jprb) :: absch_d
  REAL   (KIND=jprb) :: scach_d
  REAL   (KIND=jprb) :: bparh_d
!       Local arrays:
!      REAL OPDPE      (KNCHPF,JPLAY)
!      REAL OPDPEAER   (KNCHPF,JPLAY)
  REAL   (KIND=jprb) :: opdpaerl   (size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpcldl   (size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpcldlsun(size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: opdpaerlsun(size(chanprof)          , nlayers)
  REAL   (KIND=jprb) :: cosan_aer  (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: cosan_wcl  (coef_scatt_ir%fmv_wcl_ph         )
  REAL   (KIND=jprb) :: cosan_icl  (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: pfac       (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: phash      (coef_scatt_ir%fmv_aer_ph         )
  REAL   (KIND=jprb) :: phasice    (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: phasice_d  (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: deltap     (coef_scatt_ir%fmv_icl_ph         )
  REAL   (KIND=jprb) :: phash_d    (coef_scatt_ir%fmv_aer_ph         )
  INTEGER(KIND=jpim) :: nchannels                                                      ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE
!-----End of header-----------------------------------------------------------------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_K', 0_jpim, ZHOOK_HANDLE)
  nchannels        = size(chanprof)
  IANG             = 360_jpim / (JPAZN + 1)
  IEND             = IANG * JPAZN
  OPD = 0._jprb
  OPDSUN           = 0._jprb
  FRACH            = 0._jprb
  ABSCH            = 0._jprb
  SCACH            = 0._jprb
  BPARH            = 0._jprb
  PHASINT          = 0._jprb
  SCATTANGLE       = 0._jprb
  MUSAT            = 0._jprb
  MUSUN            = 0._jprb
  PHDO             = 0._jprb
  PHUP             = 0._jprb
  PHASH            = 0._jprb
  phasice          = 0._jprb
  DELTAPaer_K      = 0._jprb
  DELTAPice_K      = 0._jprb
  opdpaerl(:,:)    = 0._jprb
  opdpaerlsun(:,:) = 0._jprb
  opdpcldl(:,:)    = 0._jprb
  opdpcldlsun(:,:) = 0._jprb
!-----Compute optical parameters for each stream------------------------------------------
  DO j = nchannels, 1,  - 1
    ich = chanprof(j)%chan
    ipf = chanprof(j)%prof
    DO ist = ircld%nstream(ipf), 0,  - 1
      IF (ist == 0) THEN
        IF (opts%addaerosl) THEN
          DO lay = nlayers, 1,  - 1
            lev = lay + 1
            opd = opd + transmission_scatt_ir_stream_k%opdpac(IST, J, lev)
            transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) =      &
              & transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) + opd
            opdpaerl(j, lay)                                    =      &
              & opdpaerl(j, lay) + transmission_scatt_ir_stream_k%opdpacl(ist, j, lay)
            IF (sun(j)) THEN
              opdsun = opdsun + transmission_scatt_ir_stream_k%opdpacsun(ist, j, lev)
              transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) =      &
                & transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) + opdsun
              opdpaerlsun(j, lay)                                    =      &
                & opdpaerlsun(j, lay) + transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay)
              transmission_scatt_ir_k%opdpsaer(j, lay)               = transmission_scatt_ir_k%opdpsaer(j, lay) +      &
                & transmission_scatt_ir_stream_k%opdpsca(ist, j, lay) *                                                &
                & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ich)
              raytracing_k%pathsat(lay, j)                           = raytracing_k%pathsat(lay, j) +             &
                & transmission_scatt_ir_stream_k%opdpsca(ist, j, lay) * transmission_scatt_ir%opdpsaer(j, lay) *  &
                & coef%ff_gam(ich)
              raytracing_k%pathsun(lay, j)                           = raytracing_k%pathsun(lay, j) +             &
                & transmission_scatt_ir_stream_k%opdpsca(ist, j, lay) * transmission_scatt_ir%opdpsaer(j, lay) *  &
                & coef%ff_gam(ich)
              transmission_scatt_ir_k%opdpaaer(j, lay)               = transmission_scatt_ir_k%opdpaaer(j, lay) +      &
                & transmission_scatt_ir_stream_k%opdpabs(ist, j, lay) *                                                &
                & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ich)
              raytracing_k%pathsat(lay, j)                           = raytracing_k%pathsat(lay, j) +             &
                & transmission_scatt_ir_stream_k%opdpabs(ist, j, lay) * transmission_scatt_ir%opdpaaer(j, lay) *  &
                & coef%ff_gam(ich)
              raytracing_k%pathsun(lay, j)                           = raytracing_k%pathsun(lay, j) +             &
                & transmission_scatt_ir_stream_k%opdpabs(ist, j, lay) * transmission_scatt_ir%opdpaaer(j, lay) *  &
                & coef%ff_gam(ich)
              transmission_scatt_ir_k%azphaerdoa(j, lay)             =      &
                & transmission_scatt_ir_k%azphaerdoa(j, lay) + transmission_scatt_ir_stream_k%azphacdo(ist, j, lay)
              transmission_scatt_ir_k%azphaerupa(j, lay)             =      &
                & transmission_scatt_ir_k%azphaerupa(j, lay) + transmission_scatt_ir_stream_k%azphacup(ist, j, lay)
              transmission_scatt_ir_k%gparaer(j, lay)                =      &
                & transmission_scatt_ir_k%gparaer(j, lay) + transmission_scatt_ir_stream_k%bcksp(ist, j, lay)
            ENDIF
          ENDDO
          opd = 0._jprb
          opdsun = 0._jprb
          transmission_scatt_ir_stream_k%opdpac(IST, J, 1)    = 0._jprb
          transmission_scatt_ir_stream_k%opdpacsun(IST, J, 1) = 0._jprb
        ELSE IF ((.NOT. opts%addaerosl) .AND. (opts%addclouds)) THEN
          transmission_scatt_ir_stream_k%opdpaclsun(ist, j, :) = 0._jprb
          transmission_scatt_ir_stream_k%OPDPACSUN(IST, J, :)  = 0._jprb
          transmission_scatt_ir_stream_k%OPDPACL(IST, J, :)    = 0._jprb
          transmission_scatt_ir_stream_k%opdpac(ist, j, :)     = 0._jprb
          transmission_scatt_ir_stream_k%azphacup(ist, j, :)   = 0._jprb
          transmission_scatt_ir_stream_k%azphacdo(ist, j, :)   = 0._jprb
        ENDIF
      ELSE
        IF (opts%addclouds) THEN
          DO lay = nlayers, 1,  - 1
            lev = lay + 1
!                opdsun              = opdsun                                             + &
!                              transmission_scatt_ir_stream_k%opdpacsun(ist,j,lev)
            opd = opd + transmission_scatt_ir_stream_k%opdpac(ist, j, lev)
!                transmission_scatt_ir_stream_k%opdpaclsun(ist,j,lay)=                 &
!                              transmission_scatt_ir_stream_k%opdpaclsun(ist,j,lay)  + &
!                              opdsun
            transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) =      &
              & transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) + opd
!                opdpcldlsun(j,lay)   = opdpcldlsun(j,lay)                                  + &
!                              transmission_scatt_ir_stream_k%opdpaclsun(ist,j,lay)  * &
!                              ircld%icldarr(ist,lay,ipf)
!                opdpaerlsun(j,lay)   = opdpaerlsun(j,lay)                                  + &
!                              transmission_scatt_ir_stream_k%opdpaclsun(ist,j,lay)
            opdpcldl(j, lay)                                    =      &
              & opdpcldl(j, lay) + transmission_scatt_ir_stream_k%opdpacl(ist, j, lay) * ircld%icldarr(ist, lay, ipf)
            opdpaerl(j, lay)                                    =      &
              & opdpaerl(j, lay) + transmission_scatt_ir_stream_k%opdpacl(ist, j, lay)
            IF (sun(j)) THEN
              opdsun = opdsun + transmission_scatt_ir_stream_k%opdpacsun(ist, j, lev)
              transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) =      &
                & transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) + opdsun
              opdpcldlsun(j, lay)                                    = opdpcldlsun(j, lay) +      &
                & transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay) * ircld%icldarr(ist, lay, ipf)
              opdpaerlsun(j, lay)                                    =      &
                & opdpaerlsun(j, lay) + transmission_scatt_ir_stream_k%opdpaclsun(ist, j, lay)
              IF ((ircld%icldarr(ist, lay, ipf) * transmission_scatt_ir%opdps(j, lay) +      &
                & transmission_scatt_ir%opdpsaer(j, lay)) /= 0._jprb) THEN
                raytracing_k%pathsun(lay, j)               = raytracing_k%pathsun(lay, j) +                               &
                  & transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) &
                  &  * transmission_scatt_ir%OPDPS(J, lay)
                raytracing_k%pathsun(lay, j)               = raytracing_k%pathsun(lay, j) +      &
                  & transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay) * coef%ff_gam(ICH) *     &
                  & transmission_scatt_ir%OPDPSAER(J, lay)
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) +      &
                  & transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH)
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) +      &
                  & transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay) *                                    &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH)
                raytracing_k%pathsat(lay, j)               = raytracing_k%pathsat(lay, j) +                               &
                  & transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) &
                  &  * transmission_scatt_ir%OPDPS(J, lay)
                raytracing_k%pathsat(lay, j)               = raytracing_k%pathsat(lay, j) +      &
                  & transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay) * coef%ff_gam(ICH) *     &
                  & transmission_scatt_ir%OPDPSAER(J, lay)
                raytracing_k%pathsun(lay, j)               = raytracing_k%pathsun(lay, j) +                               &
                  & transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) &
                  &  * transmission_scatt_ir%OPDPA(J, lay)
                raytracing_k%pathsun(lay, j)               = raytracing_k%pathsun(lay, j) +      &
                  & transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay) * coef%ff_gam(ICH) *     &
                  & transmission_scatt_ir%OPDPAAER(J, lay)
                transmission_scatt_ir_k%OPDPA(J, lay)      = transmission_scatt_ir_k%OPDPA(J, lay) +      &
                  & transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH)
                transmission_scatt_ir_k%OPDPAAER(J, lay)   = transmission_scatt_ir_k%OPDPAAER(J, lay) +      &
                  & transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay) *                                    &
                  & (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH)
                raytracing_k%pathsat(lay, j)               = raytracing_k%pathsat(lay, j) +                               &
                  & transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay) * coef%ff_gam(ICH) * ircld%ICLDARR(IST, lay, ipf) &
                  &  * transmission_scatt_ir%OPDPA(J, lay)
                raytracing_k%pathsat(lay, j)               = raytracing_k%pathsat(lay, j) +      &
                  & transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay) * coef%ff_gam(ICH) *     &
                  & transmission_scatt_ir%OPDPAAER(J, lay)
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) -       &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * transmission_scatt_ir%GPAR(J, lay) *  &
                  & transmission_scatt_ir%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) / (                    &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                      &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) -      &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *    &
                  & transmission_scatt_ir%GPAR(J, lay) * transmission_scatt_ir%OPDPS(J, lay) / (          &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) +      &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *    &
                  & transmission_scatt_ir%GPAR(J, lay) / (                                                &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%GPAR(J, lay)       = transmission_scatt_ir_k%GPAR(J, lay) +                       &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * transmission_scatt_ir%OPDPS(J, lay) *             &
                  & ircld%ICLDARR(IST, lay, ipf) / (ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +  &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) -          &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * transmission_scatt_ir%GPARAER(J, lay) *  &
                  & transmission_scatt_ir%OPDPSAER(J, lay) / (                                                   &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                         &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) -      &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *    &
                  & transmission_scatt_ir%GPARAER(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (    &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                  &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) +           &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * transmission_scatt_ir%GPARAER(J, lay) / ( &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                          &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%GPARAER(J, lay)    = transmission_scatt_ir_k%GPARAER(J, lay) +             &
                  & transmission_scatt_ir_stream_k%BCKSP(IST, J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / ( &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                           &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) -            &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * transmission_scatt_ir%AZPHDO(J, lay) *  &
                  & transmission_scatt_ir%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) / (                         &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                           &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) -       &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & transmission_scatt_ir%AZPHDO(J, lay) * transmission_scatt_ir%OPDPS(J, lay) / (         &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) +       &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & transmission_scatt_ir%AZPHDO(J, lay) / (                                               &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%AZPHDO(J, lay)     = transmission_scatt_ir_k%AZPHDO(J, lay) +                     &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * transmission_scatt_ir%OPDPS(J, lay) *          &
                  & ircld%ICLDARR(IST, lay, ipf) / (ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +  &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) -                &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * transmission_scatt_ir%AZPHAERDOA(J, lay) *  &
                  & transmission_scatt_ir%OPDPSAER(J, lay) / (                                                         &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                               &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) -       &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & transmission_scatt_ir%AZPHAERDOA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (  &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) +                 &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * transmission_scatt_ir%AZPHAERDOA(J, lay) / ( &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%AZPHAERDOA(J, lay) = transmission_scatt_ir_k%AZPHAERDOA(J, lay) +             &
                  & transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / ( &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                              &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) -            &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * transmission_scatt_ir%AZPHUP(J, lay) *  &
                  & transmission_scatt_ir%OPDPS(J, lay) * ircld%ICLDARR(IST, lay, ipf) / (                         &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                           &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) -       &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & transmission_scatt_ir%AZPHUP(J, lay) * transmission_scatt_ir%OPDPS(J, lay) / (         &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) +       &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & transmission_scatt_ir%AZPHUP(J, lay) / (                                               &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%AZPHUP(J, lay)     = transmission_scatt_ir_k%AZPHUP(J, lay) +                     &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * transmission_scatt_ir%OPDPS(J, lay) *          &
                  & ircld%ICLDARR(IST, lay, ipf) / (ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +  &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) -                &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * transmission_scatt_ir%AZPHAERUPA(J, lay) *  &
                  & transmission_scatt_ir%OPDPSAER(J, lay) / (                                                         &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                               &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPS(J, lay)      = transmission_scatt_ir_k%OPDPS(J, lay) -       &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * ircld%ICLDARR(IST, lay, ipf) *  &
                  & transmission_scatt_ir%AZPHAERUPA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / (  &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                   &
                  & transmission_scatt_ir%OPDPSAER(J, lay)) ** 2
                transmission_scatt_ir_k%OPDPSAER(J, lay)   = transmission_scatt_ir_k%OPDPSAER(J, lay) +                 &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * transmission_scatt_ir%AZPHAERUPA(J, lay) / ( &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                                &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
                transmission_scatt_ir_k%AZPHAERUPA(J, lay) = transmission_scatt_ir_k%AZPHAERUPA(J, lay) +             &
                  & transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) * transmission_scatt_ir%OPDPSAER(J, lay) / ( &
                  & ircld%ICLDARR(IST, lay, ipf) * transmission_scatt_ir%OPDPS(J, lay) +                              &
                  & transmission_scatt_ir%OPDPSAER(J, lay))
              ELSE
                transmission_scatt_ir_stream_k%AZPHACUP(IST, J, lay) = 0._jprb
                transmission_scatt_ir_stream_k%AZPHACDO(IST, J, lay) = 0._jprb
                transmission_scatt_ir_stream_k%BCKSP(IST, J, lay)    = 0._jprb
                transmission_scatt_ir_stream_k%OPDPSCA(IST, J, lay)  = 0._jprb
                transmission_scatt_ir_stream_k%OPDPABS(IST, J, lay)  = 0._jprb
              ENDIF
            ENDIF
          ENDDO
          opd = 0._jprb
          opdsun = 0._jprb
          transmission_scatt_ir_stream_k%opdpac(ist, j, 1)    = 0._jprb
          transmission_scatt_ir_stream_k%opdpacsun(ist, j, 1) = 0._jprb
        ENDIF
      ENDIF
    ENDDO
  ENDDO
!-------------------------------------------------------------------------------
!         2.   CALCULATE OPTICAL DEPTHS OF CLOUDS
!-------------------------------------------------------------------------------
  IF (opts%addclouds) THEN
    DO j = nchannels, 1,  - 1
      ich = chanprof(j)%chan
      ipf = chanprof(j)%prof
      ish = profiles(ipf)%ish
!        IOFF=NUMAE
      DO lay = nlayers, 1,  - 1
        lev = lay + 1
!             ITMP = ITEMP  (lay,IPF)
!             ICONF= CLDTYP (lay,IPF)+ITEMP(lay,IPF)-1
        IF (sun(j)) THEN
          ich1 = ich - coef_scatt_ir%fmv_wcl_pha_ioff + 1
!---------------Compute cloud  optical parameters ----------------------------------------
          transmission_scatt_ir_k%OPDPCLDLA(J, lay) = transmission_scatt_ir_k%OPDPCLDLA(J, lay) +      &
            & OPDPCLDLSUN(J, lay) * (raytracing%pathsat(lay, ipf) + raytracing%pathsun(lay, ipf)) * coef%ff_gam(ICH)
          raytracing_k%pathsat(lay, j)              = raytracing_k%pathsat(lay, j) +      &
            & OPDPCLDLSUN(J, lay) * transmission_scatt_ir%OPDPCLDLA(J, lay) * coef%ff_gam(ICH)
          raytracing_k%pathsun(lay, j)              = raytracing_k%pathsun(lay, j) +      &
            & OPDPCLDLSUN(J, lay) * transmission_scatt_ir%OPDPCLDLA(J, lay) * coef%ff_gam(ICH)
          IF (transmission_scatt_ir%OPDPS(J, lay) /= 0._JPRB) THEN
            transmission_scatt_ir_k%AZPHUPTOT(J, lay) = transmission_scatt_ir_k%AZPHUPTOT(J, lay) +      &
              & transmission_scatt_ir_k%AZPHUP(J, lay) / transmission_scatt_ir%OPDPS(J, lay)
            transmission_scatt_ir_k%OPDPS(J, lay)     = transmission_scatt_ir_k%OPDPS(J, lay) -      &
              & transmission_scatt_ir_k%AZPHUP(J, lay) * transmission_scatt_ir%AZPHUPTOT(J, lay) /   &
              & (transmission_scatt_ir%OPDPS(J, lay) ** 2)
            transmission_scatt_ir_k%AZPHDOTOT(J, lay) = transmission_scatt_ir_k%AZPHDOTOT(J, lay) +      &
              & transmission_scatt_ir_k%AZPHDO(J, lay) / transmission_scatt_ir%OPDPS(J, lay)
            transmission_scatt_ir_k%OPDPS(J, lay)     = transmission_scatt_ir_k%OPDPS(J, lay) -      &
              & transmission_scatt_ir_k%AZPHDO(J, lay) * transmission_scatt_ir%AZPHDOTOT(J, lay) /   &
              & (transmission_scatt_ir%OPDPS(J, lay) ** 2)
          ENDIF
        ENDIF
        transmission_scatt_ir_k%OPDPCLDLA(J, lay) =      &
          & transmission_scatt_ir_k%OPDPCLDLA(J, lay) + OPDPCLDL(J, lay) * raytracing%pathsat(lay, ipf) * coef%ff_gam(ICH)
        raytracing_k%pathsat(lay, j)              =      &
          & raytracing_k%pathsat(lay, j) + OPDPCLDL(J, lay) * transmission_scatt_ir%OPDPCLDLA(J, lay) * coef%ff_gam(ICH)
        transmission_scatt_ir_k%OPDPA(J, lay)     =      &
          & transmission_scatt_ir_k%OPDPA(J, lay) + transmission_scatt_ir_k%OPDPCLDLA(J, lay)
        transmission_scatt_ir_k%OPDPS(J, lay)     = transmission_scatt_ir_k%OPDPS(J, lay) +      &
          & transmission_scatt_ir_k%OPDPCLDLA(J, lay) * transmission_scatt_ir%GPAR(J, lay)
        transmission_scatt_ir_k%GPAR(J, lay)      = transmission_scatt_ir_k%GPAR(J, lay) +      &
          & transmission_scatt_ir_k%OPDPCLDLA(J, lay) * transmission_scatt_ir%OPDPS(J, lay)
        IF (transmission_scatt_ir%OPDPS(J, lay) /= 0._JPRB) THEN
          transmission_scatt_ir_k%GPARTOT(J, lay) = transmission_scatt_ir_k%GPARTOT(J, lay) +      &
            & transmission_scatt_ir_k%GPAR(J, lay) / transmission_scatt_ir%OPDPS(J, lay)
          transmission_scatt_ir_k%OPDPS(J, lay)   = transmission_scatt_ir_k%OPDPS(J, lay) -      &
            & transmission_scatt_ir_k%GPAR(J, lay) * transmission_scatt_ir%GPARTOT(J, lay) /     &
            & (transmission_scatt_ir%OPDPS(J, lay) ** 2)
        ENDIF
!-----------------------------------------------------------------------------------------
!                 For water clouds use stored optical parameters
!-----------------------------------------------------------------------------------------
        DO lctyp = 1, ncldtyp
          IF (ircld%cldtyp(lctyp, lay, ipf) /= 0) THEN
            ityp  = ircld%cldtyp(lctyp, lay, ipf)
            iconf = ityp
!-----------------Average phase function for the downward scattered solar beam------------
            IF (sun(j)) THEN
              IF (ITYP == 6_jpim) THEN
!-----------------If ice clouds are present the phase function for for the current value
!                 of the effective generalized diameter is obtained by linear
!                 interpolation and then the azimuthally averaged value is computed
!-----------------------------------------------------------------------------------------
                ich1 = ich - coef_scatt_ir%fmv_wcl_pha_ioff + 1
                DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                  IF (aux%dg(lay, ipf) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                    & aux%dg(lay, ipf) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                    deltap(:)    = (optp%optpicl(ish)%PHA(ICH1, k + 1, :) - optp%optpicl(ish)%PHA(ICH1, k, :))
                    DELTAice     = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                    phasice_d(:) = optp%optpicl(ish)%PHA(ICH1, k, :) +      &
                      & deltap(:) * (aux%dg(lay, ipf) - coef_scatt_ir%fmv_icl_dg(k, ish)) / DELTAice
                    EXIT
                  ENDIF
                ENDDO
!-----------------Average phase function for the downward scattered solar beam------------
                cosan_icl(1:coef_scatt_ir%fmv_icl_ph)           =      &
                  & cos(coef_scatt_ir%fmv_icl_ph_val(1:coef_scatt_ir%fmv_icl_ph) * deg2rad)
                musat_d =  - 1 / raytracing%pathsat(lay, ipf)
                musun_d =  - 1 / raytracing%pathsun(lay, ipf)
                transmission_scatt_ir_k%AZPHDOCLS(J, lay, ityp) = transmission_scatt_ir_k%AZPHDOCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHDOTOT(J, lay) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
                transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp)  = transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHDOTOT(J, lay) * transmission_scatt_ir%AZPHDOCLS(J, lay, ityp)
                loop4 : DO k = iend, 0,  - iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt(1 - musat_d ** 2) * sqrt(1 - musun_d ** 2) * cos(k * deg2rad)
                  phasint      = phasint + transmission_scatt_ir_k%azphdocls(j, lay, ityp) * iang * deg2rad / (2 * pi)
                  DO kk = coef_scatt_ir%fmv_icl_ph - 1, 1,  - 1
                    IF (scattangle_d >= cosan_iCL(kk + 1) .AND. scattangle_d <= cosan_iCL(kk)) THEN
                      deltapice_d     = (PHAsice_d(KK + 1) - PHAsice_d(KK))
                      deltaice        = (cosan_icl(kk) - cosan_icl(kk + 1))
                      PHAsice(KK)     = PHAsice(KK) + phasint
                      DELTAPice_k     = DELTAPice_k + phasint * (COSAN_icl(KK) - SCATTANGLE_d) / DELTAice
                      SCATTANGLE      = SCATTANGLE - phasint * DELTAPice_d / DELTAice
                      phasint         = 0._jprb
                      PHAsice(KK + 1) = PHAsice(KK + 1) + deltapice_k
                      PHAsice(KK)     = PHAsice(KK) - deltapice_k
                      deltapice_k     = 0._jprb
                      EXIT
                    ENDIF
                  ENDDO
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._JPRB) THEN
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    scattangle = 0._jprb
                  ELSE
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    musat      =      &
                      & musat - scattangle * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1 - musat_d ** 2)
                    musun      =      &
                      & musun - scattangle * musun_d * sqrt(1 - musat_d ** 2) * cos(k * deg2rad) / sqrt(1 - musun_d ** 2)
                    scattangle = 0._jprb
                  ENDIF
                ENDDO loop4
                raytracing_k%pathsat(lay, j)                    =      &
                  & raytracing_k%pathsat(lay, j) + MUSAT / raytracing%pathsat(lay, ipf) ** 2
                raytracing_k%pathsun(lay, j)                    =      &
                  & raytracing_k%pathsun(lay, j) + MUSUN / raytracing%pathsun(lay, ipf) ** 2
                MUSAT = 0._jprb
                MUSUN = 0._jprb
                PHASINT = 0._jprb
!-----------------Average phase function for the upward scattered solar beam--------------
                cosan_icl(1:coef_scatt_ir%fmv_icl_ph)           =      &
                  & cos(coef_scatt_ir%fmv_icl_ph_val(1:coef_scatt_ir%fmv_icl_ph) * deg2rad)
                musat_d = 1._jprb / raytracing%pathsat(lay, ipf)
                musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
                transmission_scatt_ir_k%AZPHUPCLS(J, lay, ityp) = transmission_scatt_ir_k%AZPHUPCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHUPTOT(J, lay) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
                transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp)  = transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHUPTOT(J, lay) * transmission_scatt_ir%AZPHUPCLS(J, lay, ityp)
                loop3 : DO k = iend, 0,  - iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt(1 - musat_d ** 2) * sqrt(1 - musun_d ** 2) * cos(k * deg2rad)
                  phasint      = phasint + transmission_scatt_ir_k%AZPHUPCLS(J, lay, ityp) * iang * deg2rad / (2 * pi)
                  DO kk = coef_scatt_ir%fmv_icl_ph - 1, 1,  - 1
                    IF (scattangle_d >= cosan_icl(kk + 1) .AND. scattangle_d <= cosan_icl(kk)) THEN
                      deltapice_d     = (PHAsice_d(KK + 1) - PHAsice_d(KK))
                      deltaice        = (COSAN_icl(KK) - COSAN_icl(KK + 1))
                      PHAsice(KK)     = PHAsice(KK) + phasint
                      DELTAPice_k     = DELTAPice_k + phasint * (COSAN_icl(KK) - SCATTANGLE_d) / DELTAice
                      SCATTANGLE      = SCATTANGLE - phasint * DELTAPice_d / DELTAice
                      phasint         = 0._jprb
                      PHAsice(KK + 1) = PHAsice(KK + 1) + deltapice_k
                      PHAsice(KK)     = PHAsice(KK) - deltapice_k
                      deltapice_k     = 0._jprb
                      EXIT
                    ENDIF
                  ENDDO
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    scattangle = 0._jprb
                  ELSE
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    musat      =      &
                      & musat - scattangle * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1 - musat_d ** 2)
                    musun      =      &
                      & musun - scattangle * musun_d * sqrt(1 - musat_d ** 2) * cos(k * deg2rad) / sqrt(1 - musun_d ** 2)
                    scattangle = 0._jprb
                  ENDIF
                ENDDO loop3
                raytracing_k%pathsat(lay, j) =      &
                  & raytracing_k%pathsat(lay, j) - MUSAT / raytracing%pathsat(lay, ipf) ** 2
                raytracing_k%pathsun(lay, j) =      &
                  & raytracing_k%pathsun(lay, j) + MUSUN / raytracing%pathsun(lay, ipf) ** 2
                musat = 0._jprb
                musun = 0._jprb
                phasint                      = 0._jprb
!-----------------If ice clouds are present the phase function for for the current value
!                 of the effective generalized diameter is obtained by linear
!                 interpolation and then the azimuthally averaged value is computed
!----------------------------------------------------------------------------------------
                DO k = 1, coef_scatt_ir%fmv_icl_comp - 1
                  IF (aux%dg(lay, ipf) >= coef_scatt_ir%fmv_icl_dg(k, ish) .AND.      &
                    & aux%dg(lay, ipf) <= coef_scatt_ir%fmv_icl_dg(k + 1, ish)) THEN
                    deltap(:) = (optp%optpicl(ish)%PHA(ICH1, k + 1, :) - optp%optpicl(ish)%PHA(ICH1, k, :))
                    DELTAice  = (coef_scatt_ir%fmv_icl_dg(k + 1, ish) - coef_scatt_ir%fmv_icl_dg(k, ish))
                    DO kk = 1, coef_scatt_ir%fmv_icl_ph
                      aux_k%dg(lay, j) = aux_k%dg(lay, j) + phasice(kk) * deltap(kk) / DELTAice
                      phasice(kk)      = 0._jprb
                    ENDDO
                    EXIT
                  ENDIF
                ENDDO
              ELSE
!-----------------------------------------------------------------------------------------
!                     Water clouds
!-----------------------------------------------------------------------------------------
                ich1 = ich - coef_scatt_ir%fmv_wcl_pha_ioff + 1
                cosan_wcl(1:coef_scatt_ir%fmv_wcl_ph)           =      &
                  & cos(coef_scatt_ir%fmv_wcl_ph_val(1:coef_scatt_ir%fmv_wcl_ph) * deg2rad)
                musat_d =  - 1 / raytracing%pathsat(lay, ipf)
                musun_d =  - 1 / raytracing%pathsun(lay, ipf)
                transmission_scatt_ir_k%AZPHDOCLS(J, lay, ityp) = transmission_scatt_ir_k%AZPHDOCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHDOTOT(J, lay) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
                transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp)  = transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHDOTOT(J, lay) * transmission_scatt_ir%AZPHDOCLS(J, lay, ityp)
                loop6 : DO k = iend, 0,  - iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt(1 - musat_d ** 2) * sqrt(1 - musun_d ** 2) * cos(k * deg2rad)
                  phasint      = phasint + transmission_scatt_ir_k%azphdocls(j, lay, ityp) * iang * deg2rad / (2 * pi)
                  DO kk = coef_scatt_ir%fmv_wcl_ph - 1, 1,  - 1
                    IF (scattangle_d >= cosan_WCL(kk + 1) .AND. scattangle_d <= cosan_WCL(kk)) THEN
                      deltapwcl_d = (optp%optpwcl(ityp)%pha(ich1, 1, kk + 1) - optp%optpwcl(ityp)%pha(ich1, 1, kk))
                      deltawcl    = (cosan_wcl(kk) - cosan_wcl(kk + 1))
                      scattangle  = scattangle - phasint * deltapwcl_d / deltawcl
                      phasint     = 0._jprb
                      EXIT
                    ENDIF
                  ENDDO
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._JPRB) THEN
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    scattangle = 0._jprb
                  ELSE
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    musat      =      &
                      & musat - scattangle * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1 - musat_d ** 2)
                    musun      =      &
                      & musun - scattangle * musun_d * sqrt(1 - musat_d ** 2) * cos(k * deg2rad) / sqrt(1 - musun_d ** 2)
                    scattangle = 0._jprb
                  ENDIF
                ENDDO loop6
                raytracing_k%pathsat(lay, j)                    =      &
                  & raytracing_k%pathsat(lay, j) + MUSAT / raytracing%pathsat(lay, ipf) ** 2
                raytracing_k%pathsun(lay, j)                    =      &
                  & raytracing_k%pathsun(lay, j) + MUSUN / raytracing%pathsun(lay, ipf) ** 2
                MUSAT = 0._jprb
                MUSUN = 0._jprb
                PHASINT = 0._jprb
!-----------------Average phase function for the upward scattered solar beam--------------
                cosan_wcl(1:coef_scatt_ir%fmv_wcl_ph)           =      &
                  & cos(coef_scatt_ir%fmv_wcl_ph_val(1:coef_scatt_ir%fmv_wcl_ph) * deg2rad)
                musat_d = 1._jprb / raytracing%pathsat(lay, ipf)
                musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
                transmission_scatt_ir_k%AZPHUPCLS(J, lay, ityp) = transmission_scatt_ir_k%AZPHUPCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHUPTOT(J, lay) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
                transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp)  = transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) +      &
                  & transmission_scatt_ir_k%AZPHUPTOT(J, lay) * transmission_scatt_ir%AZPHUPCLS(J, lay, ityp)
                loop5 : DO k = iend, 0,  - iang
                  scattangle_d =      &
                    & musat_d * musun_d + sqrt(1 - musat_d ** 2) * sqrt(1 - musun_d ** 2) * cos(k * deg2rad)
                  phasint      = phasint + transmission_scatt_ir_k%AZPHUPCLS(J, lay, ityp) * iang * deg2rad / (2 * pi)
                  DO kk = coef_scatt_ir%fmv_wcl_ph - 1, 1,  - 1
                    IF (scattangle_d >= cosan_wcl(kk + 1) .AND. scattangle_d <= cosan_wcl(kk)) THEN
                      deltapwcl_d = (optp%optpwcl(ityp)%PHA(ICH1, 1, KK + 1) - optp%optpwcl(ityp)%PHA(ICH1, 1, KK))
                      deltawcl    = (COSAN_wcl(KK) - COSAN_wcl(KK + 1))
                      scattangle  = scattangle - phasint * deltapwcl_d / deltawcl
                      phasint     = 0._jprb
                      EXIT
                    ENDIF
                  ENDDO
                  IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    scattangle = 0._jprb
                  ELSE
                    musat      = musat + scattangle * musun_d
                    musun      = musun + scattangle * musat_d
                    musat      =      &
                      & musat - scattangle * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1 - musat_d ** 2)
                    musun      =      &
                      & musun - scattangle * musun_d * sqrt(1 - musat_d ** 2) * cos(k * deg2rad) / sqrt(1 - musun_d ** 2)
                    scattangle = 0._jprb
                  ENDIF
                ENDDO loop5
                raytracing_k%pathsat(lay, j) =      &
                  & raytracing_k%pathsat(lay, j) - MUSAT / raytracing%pathsat(lay, ipf) ** 2
                raytracing_k%pathsun(lay, j) =      &
                  & raytracing_k%pathsun(lay, j) + MUSUN / raytracing%pathsun(lay, ipf) ** 2
                musat = 0._jprb
                musun = 0._jprb
                phasint                      = 0._jprb
              ENDIF
            ENDIF
            IF (ityp <= 5_jpim) THEN
              transmission_scatt_ir_k%GPAR(J, lay)           = 0._jprb
!              profiles_k(j)%cloud(Ityp,lay)= profiles_k(j)%cloud(Ityp,lay)       + &
!                          OPDPE(J,lay)         * &
!                          CONFAC(ICONF)                            * &
!                          EXTC(ICH,IOFF+ITYP,ITMP)*LTICK_D(90-lay,IPF)
!              LTICK(90-lay,IPF)    = LTICK(90-lay,IPF)+OPDPE(J,lay)             * &
!                                             PROFCLD_D(ITYP,lay,IPF)          * &
!                                    CONFAC(ICONF)*EXTC(ICH,IOFF+ITYP,ITMP)
              transmission_scatt_ir_k%GPARCLS(J, lay, ityp)  = transmission_scatt_ir_k%GPARCLS(J, lay, ityp) +      &
                & transmission_scatt_ir_k%GPARTOT(J, lay) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
              transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) = transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) +      &
                & transmission_scatt_ir_k%GPARTOT(J, lay) * transmission_scatt_ir%GPARCLS(J, lay, ityp)
              transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) =      &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) + transmission_scatt_ir_k%OPDPS(J, lay)
              profiles_k(j)%cloud(Ityp, lay)                 = profiles_k(j)%cloud(Ityp, lay) +      &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * coef_scatt_ir%CONFAC(ICONF) *     &
                & optp%optpwcl(ityp)%sca(ich, 1) * raytracing%LTICK(nlayers + 1 - lay, IPF)
              raytracing_k%LTICK(nlayers + 1 - lay, j)       = raytracing_k%LTICK(nlayers + 1 - lay, j) +      &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * profiles(IPF)%cloud(Ityp, lay) *            &
                & coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%sca(ich, 1)
              transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) =      &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) + transmission_scatt_ir_k%OPDPA(J, lay)
              profiles_k(j)%cloud(Ityp, lay)                 = profiles_k(j)%cloud(Ityp, lay) +      &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * coef_scatt_ir%CONFAC(ICONF) *     &
                & optp%optpwcl(ityp)%abs(ich, 1) * raytracing%LTICK(nlayers + 1 - lay, IPF)
              raytracing_k%LTICK(nlayers + 1 - lay, j)       = raytracing_k%LTICK(nlayers + 1 - lay, j) +      &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * profiles(IPF)%cloud(Ityp, lay) *            &
                & coef_scatt_ir%CONFAC(ICONF) * optp%optpwcl(ityp)%abs(ich, 1)
            ELSE
!-----------------------------------------------------------------------------------------
!                 For ice clouds optical parameters are computed using regression
!                 coefficients
!-----------------------------------------------------------------------------------------
              transmission_scatt_ir_k%GPARCLS(J, lay, ityp)  = transmission_scatt_ir_k%GPARCLS(J, lay, ityp) +      &
                & transmission_scatt_ir_k%GPARTOT(J, lay) * transmission_scatt_ir%OPDPSCLS(J, lay, ityp)
              transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) = transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) +      &
                & transmission_scatt_ir_k%GPARTOT(J, lay) * transmission_scatt_ir%GPARCLS(J, lay, ityp)
              aux_k%dg(lay, j)                               =      &
                & aux_k%dg(lay, j) + transmission_scatt_ir_k%GPARCLS(J, lay, ityp) * optp%optpicl(ish)%bpr(ich, 2)
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) +      &
                & transmission_scatt_ir_k%GPARCLS(J, lay, ityp) * 2 * aux%dg(lay, ipf) * optp%optpicl(ish)%bpr(ich, 3)
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) +      &
                & transmission_scatt_ir_k%GPARCLS(J, lay, ityp) * 3 * aux%dg(lay, ipf) ** 2 * optp%optpicl(ish)%bpr(ich, 4)
              transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) =      &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) + transmission_scatt_ir_k%OPDPS(J, lay)
              raytracing_k%LTICK(nlayers + 1 - lay, j)       = raytracing_k%LTICK(nlayers + 1 - lay, j) +      &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * profiles(IPF)%cloud(Ityp, lay) * (          &
                & optp%optpicl(ish)%sca(ich, 1) + optp%optpicl(ish)%sca(ich, 2) * aux%dg(lay, ipf) +           &
                & optp%optpicl(ish)%sca(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%sca(ich, 4) / aux%dg(lay, ipf) ** 2)
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) -                        &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * optp%optpicl(ish)%sca(ich, 4) *       &
                & profiles(IPF)%cloud(Ityp, lay) * raytracing%LTICK(nlayers + 1 - lay, IPF) * 2._jprb /  &
                & aux%dg(lay, ipf) ** 3
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) -                   &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * optp%optpicl(ish)%sca(ich, 3) *  &
                & profiles(IPF)%cloud(Ityp, lay) * raytracing%LTICK(nlayers + 1 - lay, IPF) / aux%dg(lay, ipf) ** 2
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) +                   &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * optp%optpicl(ish)%sca(ich, 2) *  &
                & profiles(IPF)%cloud(Ityp, lay) * raytracing%LTICK(nlayers + 1 - lay, IPF)
              profiles_k(j)%cloud(Ityp, lay)                 = profiles_k(j)%cloud(Ityp, lay) +                             &
                & transmission_scatt_ir_k%OPDPSCLS(J, lay, ityp) * (                                                        &
                & optp%optpicl(ish)%sca(ich, 1) + optp%optpicl(ish)%sca(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%sca(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%sca(ich, 4) / aux%dg(lay, ipf) ** 2) &
                &  * raytracing%LTICK(nlayers + 1 - lay, IPF)
              transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) =      &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) + transmission_scatt_ir_k%OPDPA(J, lay)
              raytracing_k%LTICK(nlayers + 1 - lay, j)       = raytracing_k%LTICK(nlayers + 1 - lay, j) +      &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * profiles(IPF)%cloud(Ityp, lay) * (          &
                & optp%optpicl(ish)%abs(ich, 1) + optp%optpicl(ish)%abs(ich, 2) * aux%dg(lay, ipf) +           &
                & optp%optpicl(ish)%abs(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%abs(ich, 4) / aux%dg(lay, ipf) ** 2)
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) -                        &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * optp%optpicl(ish)%abs(ich, 4) *       &
                & profiles(IPF)%cloud(Ityp, lay) * raytracing%LTICK(nlayers + 1 - lay, IPF) * 2._jprb /  &
                & aux%dg(lay, ipf) ** 3
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) -                   &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * optp%optpicl(ish)%abs(ich, 3) *  &
                & profiles(IPF)%cloud(Ityp, lay) * raytracing%LTICK(nlayers + 1 - lay, IPF) / aux%dg(lay, ipf) ** 2
              aux_k%dg(lay, j)                               = aux_k%dg(lay, j) +                   &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * optp%optpicl(ish)%abs(ich, 2) *  &
                & profiles(IPF)%cloud(Ityp, lay) * raytracing%LTICK(nlayers + 1 - lay, IPF)
              profiles_k(j)%cloud(Ityp, lay)                 = profiles_k(j)%cloud(Ityp, lay) +                             &
                & transmission_scatt_ir_k%OPDPACLS(J, lay, ityp) * (                                                        &
                & optp%optpicl(ish)%abs(ich, 1) + optp%optpicl(ish)%abs(ich, 2) * aux%dg(lay, ipf) +                        &
                & optp%optpicl(ish)%abs(ich, 3) / aux%dg(lay, ipf) + optp%optpicl(ish)%abs(ich, 4) / aux%dg(lay, ipf) ** 2) &
                &  * raytracing%LTICK(nlayers + 1 - lay, IPF)
            ENDIF
          ENDIF
        ENDDO
        transmission_scatt_ir_k%azphup(j, lay) = 0._jprb
        transmission_scatt_ir_k%azphdo(j, lay) = 0._jprb
        transmission_scatt_ir_k%OPDPA(j, lay)  = 0._jprb
        transmission_scatt_ir_k%OPDPS(j, lay)  = 0._jprb
        transmission_scatt_ir_k%GPAR(j, lay)   = 0._jprb
!            transmission_scatt_ir_k%OPDPAAER(j,lay)=0._jprb
!            transmission_scatt_ir_k%OPDPSAER(j,lay)=0._jprb
!            transmission_scatt_ir_k%GPARAER(j,lay) =0._jprb
      ENDDO
    ENDDO
  ENDIF
!-------------------------------------------------------------------------------
!         1.   CALCULATE OPTICAL DEPTHS OF AEROSOLS
!-------------------------------------------------------------------------------
  IF (opts%addaerosl) THEN
    DO J = nchannels, 1,  - 1
      ich = chanprof(J)%chan
      ipf = chanprof(J)%prof
!---------Compute final values for optical parameters-------------------------------------
      DO lay = nlayers, 1,  - 1
        lev = lay + 1
        IF (sun(j)) THEN
          transmission_scatt_ir_k%OPDPAERLA(J, lay) = transmission_scatt_ir_k%OPDPAERLA(J, lay) +      &
            & OPDPAERLSUN(J, lay) * (raytracing%pathsun(lay, ipf) + raytracing%pathsat(lay, ipf)) * coef%ff_gam(ich)
          raytracing_k%pathsat(lay, j)              = raytracing_k%pathsat(lay, j) +      &
            & OPDPAERLSUN(J, lay) * transmission_scatt_ir%OPDPAERLA(J, lay) * coef%ff_gam(ich)
          raytracing_k%pathsun(lay, j)              = raytracing_k%pathsun(lay, j) +      &
            & OPDPAERLSUN(J, lay) * transmission_scatt_ir%OPDPAERLA(J, lay) * coef%ff_gam(ich)
        ENDIF
        transmission_scatt_ir_k%OPDPAERLA(J, lay) =      &
          & transmission_scatt_ir_k%OPDPAERLA(J, lay) + OPDPAERL(J, lay) * raytracing%pathsat(lay, ipf) * coef%ff_gam(ich)
        raytracing_k%pathsat(lay, j)              =      &
          & raytracing_k%pathsat(lay, j) + OPDPAERL(J, lay) * transmission_scatt_ir%OPDPAERLA(J, lay) * coef%ff_gam(ich)
        transmission_scatt_ir_k%OPDPAAER(J, lay)  =      &
          & transmission_scatt_ir_k%OPDPAAER(J, lay) + transmission_scatt_ir_k%OPDPAERLA(J, lay)
        transmission_scatt_ir_k%OPDPSAER(J, lay)  = transmission_scatt_ir_k%OPDPSAER(J, lay) +      &
          & transmission_scatt_ir_k%OPDPAERLA(J, lay) * transmission_scatt_ir%GPARAER(J, lay)
        transmission_scatt_ir_k%GPARAER(J, lay)   = transmission_scatt_ir_k%GPARAER(J, lay) +      &
          & transmission_scatt_ir_k%OPDPAERLA(J, lay) * transmission_scatt_ir%OPDPSAER(J, lay)
        IF (transmission_scatt_ir%opdpsaer(j, lay) /= 0._jprb) THEN
          IF (sun(j)) THEN
            transmission_scatt_ir_k%AZPHAERDO(J, lay) = transmission_scatt_ir_k%AZPHAERDO(J, lay) +      &
              & transmission_scatt_ir_k%AZPHAERDOA(J, lay) / transmission_scatt_ir%OPDPSAER(J, lay)
            transmission_scatt_ir_k%OPDPSAER(J, lay)  = transmission_scatt_ir_k%OPDPSAER(J, lay) -      &
              & transmission_scatt_ir_k%AZPHAERDOA(J, lay) * transmission_scatt_ir%AZPHAERDO(J, lay) /  &
              & transmission_scatt_ir%OPDPSAER(J, lay) ** 2
            transmission_scatt_ir_k%AZPHAERUP(J, lay) = transmission_scatt_ir_k%AZPHAERUP(J, lay) +      &
              & transmission_scatt_ir_k%AZPHAERUPA(J, lay) / transmission_scatt_ir%OPDPSAER(J, lay)
            transmission_scatt_ir_k%OPDPSAER(J, lay)  = transmission_scatt_ir_k%OPDPSAER(J, lay) -      &
              & transmission_scatt_ir_k%AZPHAERUPA(J, lay) * transmission_scatt_ir%AZPHAERUP(J, lay) /  &
              & transmission_scatt_ir%OPDPSAER(J, lay) ** 2
          ENDIF
          transmission_scatt_ir_k%GPARAERA(J, lay) = transmission_scatt_ir_k%GPARAERA(J, lay) +      &
            & transmission_scatt_ir_k%GPARAER(J, lay) / transmission_scatt_ir%OPDPSAER(J, lay)
          transmission_scatt_ir_k%OPDPSAER(J, lay) = transmission_scatt_ir_k%OPDPSAER(J, lay) -      &
            & transmission_scatt_ir_k%GPARAER(J, lay) * transmission_scatt_ir%GPARAERA(J, lay) /     &
            & transmission_scatt_ir%OPDPSAER(J, lay) ** 2
        ENDIF
      ENDDO
      DO lay = nlayers, 1,  - 1
        lev = lay + 1
        DO I = aux%iaernum(lay, ipf), 1,  - 1
          iae = aux%iaertyp(i, lay, ipf)
!-------------Repeat direct calculations--------------------------------------------------
          IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1) THEN
            DO k = 1, coef_scatt_ir%fmv_aer_rh(iae) - 1
              IF (aux%relhum(lay, ipf) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                & aux%relhum(lay, ipf) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                delth   = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                afac    = (optp%optpaer(iae)%abs(ich, k + 1) - optp%optpaer(iae)%abs(ich, k)) / delth
                sfac    = (optp%optpaer(iae)%sca(ich, k + 1) - optp%optpaer(iae)%sca(ich, k)) / delth
!                    EFAC=(EXTC(ICH,IAE,K+1)-EXTC(ICH,IAE,K))/DELTH
                gfac    = (optp%optpaer(iae)%bpr(ich, k + 1) - optp%optpaer(iae)%bpr(ich, k)) / delth
                frach_d = (aux%relhum(lay, ipf) - optp%optpaer(iae)%fmv_aer_rh_val(k))
                absch_d = optp%optpaer(iae)%abs(ich, k) + afac * frach_d
                scach_d = optp%optpaer(iae)%sca(ich, k) + sfac * frach_d
                bparh_d = optp%optpaer(iae)%bpr(ich, k) + gfac * frach_d
                IF (sun(j)) THEN
                  ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
                  pfac(1:coef_scatt_ir%fmv_aer_ph)    = (                               &
                    & optp%optpaer(iae)%pha(ich1, K + 1, 1:coef_scatt_ir%fmv_aer_ph) -  &
                    & optp%optpaer(iae)%pha(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph)) / delth
                  phash_d(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph) +      &
                    & pfac(1:coef_scatt_ir%fmv_aer_ph) * frach_d
                ENDIF
                EXIT
              ENDIF
            ENDDO
          ELSE
            absch_d = optp%optpaer(iae)%abs(ich, 1)
            scach_d = optp%optpaer(iae)%sca(ich, 1)
            bparh_d = optp%optpaer(iae)%bpr(ich, 1)
            IF (sun(j)) THEN
              ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
              phash_d(1:coef_scatt_ir%fmv_aer_ph) = optp%optpaer(iae)%pha(ich1, 1, 1:coef_scatt_ir%fmv_aer_ph)
            ENDIF
          ENDIF
!-----------------------------------------------------------------------------------------
          IF (sun(j)) THEN
            ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
!-----------------Average phase function for the downward scattered solar beam------------
            cosan_aer(1:coef_scatt_ir%fmv_aer_ph)    =      &
              & cos(coef_scatt_ir%fmv_aer_ph_val(1:coef_scatt_ir%fmv_aer_ph) * deg2rad)
            musat_d =  - 1 / raytracing%pathsat(lay, ipf)
            musun_d =  - 1 / raytracing%pathsun(lay, ipf)
            profiles_k(j)%aerosols(IAE, lay)         = profiles_k(j)%aerosols(IAE, lay) +                              &
              & transmission_scatt_ir_k%AZPHAERDO(J, lay) * transmission_scatt_ir%PHASINTDOREF(J, lay, I) * SCACH_D *  &
              & raytracing%LTICK(nlayers + 1 - lay, ipf)
            PHASINT = PHASINT + transmission_scatt_ir_k%AZPHAERDO(J, lay) * profiles(ipf)%aerosols(IAE, lay) * SCACH_D     &
              &  * raytracing%LTICK(nlayers + 1 - lay, ipf)
            SCACH = SCACH + transmission_scatt_ir_k%AZPHAERDO(J, lay) * profiles(ipf)%aerosols(IAE, lay) *      &
              & transmission_scatt_ir%PHASINTDOREF(J, lay, I) * raytracing%LTICK(nlayers + 1 - lay, ipf)
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +         &
              & transmission_scatt_ir_k%AZPHAERDO(J, lay) * SCACH_D * profiles(ipf)%aerosols(IAE, lay) *  &
              & transmission_scatt_ir%PHASINTDOREF(J, lay, I)
            loop2 : DO K = iend, 0,  - iang
              scattangle_d = musat_d * musun_d + sqrt(1 - musat_d ** 2) * sqrt(1 - musun_d ** 2) * cos(k * deg2rad)
              phdo         = phdo + phasint * iang * deg2rad / (2 * pi)
              DO kk = coef_scatt_ir%fmv_aer_ph - 1, 1,  - 1
                IF (scattangle_d >= cosan_aer(kk + 1) .AND. scattangle_d <= cosan_aer(kk)) THEN
                  deltaaer      = (cosan_aer(kk) - cosan_aer(kk + 1))
                  deltapaer_d   = (phash_d(kk + 1) - phash_d(kk))
                  phash(kk)     = phash(kk) + phdo
                  deltapaer_k   = deltapaer_k + phdo * (cosan_aer(kk) - scattangle_d) / deltaaer
                  scattangle    = scattangle - phdo * deltapaer_d / deltaaer
                  phdo          = 0._jprb
                  phash(kk + 1) = phash(kk + 1) + deltapaer_k
                  phash(kk)     = phash(kk) - deltapaer_k
                  deltapaer_k   = 0._jprb
                  EXIT
                ENDIF
              ENDDO
              IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                musat      = musat + scattangle * musun_d
                musun      = musun + scattangle * musat_d
                scattangle = 0._jprb
              ELSE
                musat      = musat + scattangle * musun_d
                musun      = musun + scattangle * musat_d
                musat      =      &
                  & musat - scattangle * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1 - musat_d ** 2)
                musun      =      &
                  & musun - scattangle * musun_d * sqrt(1 - musat_d ** 2) * cos(k * deg2rad) / sqrt(1 - musun_d ** 2)
                scattangle = 0
              ENDIF
            ENDDO loop2
            raytracing_k%pathsun(lay, j)             =      &
              & raytracing_k%pathsun(lay, j) + MUSUN / raytracing%pathsun(lay, ipf) ** 2
            raytracing_k%pathsat(lay, j)             =      &
              & raytracing_k%pathsat(lay, j) + MUSAT / raytracing%pathsat(lay, ipf) ** 2
            MUSAT = 0._jprb
            MUSUN = 0._jprb
            PHASINT = 0._jprb
!-------------Average phase function for the upward scattered solar beam------------------
            cosan_aer(1:coef_scatt_ir%fmv_aer_ph)    =      &
              & cos(coef_scatt_ir%fmv_aer_ph_val(1:coef_scatt_ir%fmv_aer_ph) * deg2rad)
            musat_d = 1._jprb / raytracing%pathsat(lay, ipf)
            musun_d =  - 1._jprb / raytracing%pathsun(lay, ipf)
            profiles_k(j)%aerosols(IAE, lay)         = profiles_k(j)%aerosols(IAE, lay) +                              &
              & transmission_scatt_ir_k%AZPHAERUP(J, lay) * transmission_scatt_ir%PHASINTUPREF(J, lay, I) * SCACH_D *  &
              & raytracing%LTICK(nlayers + 1 - lay, ipf)
            PHASINT = PHASINT + transmission_scatt_ir_k%AZPHAERUP(J, lay) * profiles(ipf)%aerosols(IAE, lay) * SCACH_D     &
              &  * raytracing%LTICK(nlayers + 1 - lay, ipf)
            SCACH = SCACH + transmission_scatt_ir_k%AZPHAERUP(J, lay) * profiles(ipf)%aerosols(IAE, lay) *      &
              & transmission_scatt_ir%PHASINTUPREF(J, lay, I) * raytracing%LTICK(nlayers + 1 - lay, ipf)
            raytracing_k%LTICK(nlayers + 1 - lay, j) = raytracing_k%LTICK(nlayers + 1 - lay, j) +         &
              & transmission_scatt_ir_k%AZPHAERUP(J, lay) * SCACH_D * profiles(ipf)%aerosols(IAE, lay) *  &
              & transmission_scatt_ir%PHASINTUPREF(J, lay, I)
            loop1 : DO k = iend, 0,  - iang
              scattangle_d = musat_d * musun_d + sqrt(1 - musat_d ** 2) * sqrt(1 - musun_d ** 2) * cos(k * deg2rad)
              phup         = phup + phasint * iang * deg2rad / (2 * pi)
              DO kk = coef_scatt_ir%fmv_aer_ph - 1, 1,  - 1
                IF (scattangle_d >= cosan_aer(kk + 1) .AND. scattangle_d <= cosan_aer(kk)) THEN
                  deltaaer      = (cosan_aer(kk) - cosan_aer(kk + 1))
                  deltapaer_d   = (phash_d(kk + 1) - phash_d(kk))
                  phash(kk)     = phash(kk) + phup
                  deltapaer_k   = deltapaer_k + phup * (cosan_aer(kk) - scattangle_d) / deltaaer
                  scattangle    = scattangle - phup * deltapaer_d / deltaaer
                  phup          = 0._jprb
                  phash(kk + 1) = phash(kk + 1) + deltapaer_k
                  phash(kk)     = phash(kk) - deltapaer_k
                  deltapaer_k   = 0._jprb
                  EXIT
                ENDIF
              ENDDO
              IF (abs(musat_d) == 1._jprb .OR. abs(musun_d) == 1._jprb) THEN
                musat      = musat + scattangle * musun_d
                musun      = musun + scattangle * musat_d
                scattangle = 0._jprb
              ELSE
                musat      = musat + scattangle * musun_d
                musun      = musun + scattangle * musat_d
                musat      =      &
                  & musat - scattangle * musat_d * sqrt(1 - musun_d ** 2) * cos(k * deg2rad) / sqrt(1 - musat_d ** 2)
                musun      =      &
                  & musun - scattangle * musun_d * sqrt(1 - musat_d ** 2) * cos(k * deg2rad) / sqrt(1 - musun_d ** 2)
                scattangle = 0._jprb
              ENDIF
            ENDDO loop1
            raytracing_k%pathsat(lay, j) = raytracing_k%pathsat(lay, j) - MUSAT / raytracing%pathsat(lay, ipf) ** 2
            raytracing_k%pathsun(lay, j) = raytracing_k%pathsun(lay, j) + MUSUN / raytracing%pathsun(lay, ipf) ** 2
            musat = 0._jprb
            musun = 0._jprb
            phasint                      = 0._jprb
          ENDIF
!-----------------------------------------------------------------------------------------
          IF (i == 1) THEN
            profiles_k(j)%aerosols(iae, lay)         = profiles_k(j)%aerosols(iae, lay) +      &
              & transmission_scatt_ir_k%GPARAERA(J, lay) * SCACH_D * BPARH_D * raytracing%ltick(nlayers + 1 - lay, ipf)
            SCACH = SCACH + transmission_scatt_ir_k%GPARAERA(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf) * BPARH_D
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +      &
              & transmission_scatt_ir_k%GPARAERA(J, lay) * SCACH_D * profiles(ipf)%aerosols(iae, lay) * BPARH_D
            BPARH = BPARH + transmission_scatt_ir_k%GPARAERA(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf) * SCACH_D
            profiles_k(j)%aerosols(iae, lay)         = profiles_k(j)%aerosols(iae, lay) +      &
              & transmission_scatt_ir_k%OPDPSAER(J, lay) * SCACH_D * raytracing%ltick(nlayers + 1 - lay, ipf)
            SCACH = SCACH + transmission_scatt_ir_k%OPDPSAER(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf)
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +      &
              & transmission_scatt_ir_k%OPDPSAER(J, lay) * SCACH_D * profiles(ipf)%aerosols(iae, lay)
            profiles_k(j)%aerosols(iae, lay)         = profiles_k(j)%aerosols(iae, lay) +      &
              & transmission_scatt_ir_k%OPDPAAER(J, lay) * ABSCH_D * raytracing%ltick(nlayers + 1 - lay, ipf)
            ABSCH = ABSCH + transmission_scatt_ir_k%OPDPAAER(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf)
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +      &
              & transmission_scatt_ir_k%OPDPAAER(J, lay) * ABSCH_D * profiles(ipf)%aerosols(iae, lay)
          ELSE
            profiles_k(j)%aerosols(iae, lay)         = profiles_k(j)%aerosols(iae, lay) +      &
              & transmission_scatt_ir_k%GPARAERA(J, lay) * SCACH_D * BPARH_D * raytracing%ltick(nlayers + 1 - lay, ipf)
            SCACH = SCACH + transmission_scatt_ir_k%GPARAERA(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf) * BPARH_D
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +      &
              & transmission_scatt_ir_k%GPARAERA(J, lay) * SCACH_D * profiles(ipf)%aerosols(iae, lay) * BPARH_D
            BPARH = BPARH + transmission_scatt_ir_k%GPARAERA(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf) * SCACH_D
            profiles_k(j)%aerosols(iae, lay)         = profiles_k(j)%aerosols(iae, lay) +      &
              & transmission_scatt_ir_k%OPDPSAER(J, lay) * SCACH_D * raytracing%ltick(nlayers + 1 - lay, ipf)
            SCACH = SCACH + transmission_scatt_ir_k%OPDPSAER(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf)
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +      &
              & transmission_scatt_ir_k%OPDPSAER(J, lay) * SCACH_D * profiles(ipf)%aerosols(iae, lay)
            profiles_k(j)%aerosols(iae, lay)         = profiles_k(j)%aerosols(iae, lay) +      &
              & transmission_scatt_ir_k%OPDPAAER(J, lay) * ABSCH_D * raytracing%ltick(nlayers + 1 - lay, ipf)
            ABSCH = ABSCH + transmission_scatt_ir_k%OPDPAAER(J, lay) * profiles(ipf)%aerosols(iae, lay) *      &
              & raytracing%ltick(nlayers + 1 - lay, ipf)
            raytracing_k%ltick(nlayers + 1 - lay, j) = raytracing_k%ltick(nlayers + 1 - lay, j) +      &
              & transmission_scatt_ir_k%OPDPAAER(J, lay) * ABSCH_D * profiles(ipf)%aerosols(iae, lay)
          ENDIF
          IF (coef_scatt_ir%fmv_aer_rh(iae) /= 1_jpim) THEN
!---------------Interpolate scattering parameters to actual value of relative--
!               humidity.
            DO k = coef_scatt_ir%fmv_aer_rh(iae) - 1, 1,  - 1
              IF (aux%relhum(lay, ipf) >= optp%optpaer(iae)%fmv_aer_rh_val(k) .AND.      &
                & aux%relhum(lay, ipf) <= optp%optpaer(iae)%fmv_aer_rh_val(k + 1)) THEN
                delth = (optp%optpaer(iae)%fmv_aer_rh_val(K + 1) - optp%optpaer(iae)%fmv_aer_rh_val(K))
                afac  = (optp%optpaer(iae)%abs(ich, k + 1) - optp%optpaer(iae)%abs(ich, k)) / delth
                sfac  = (optp%optpaer(iae)%sca(ich, k + 1) - optp%optpaer(iae)%sca(ich, k)) / delth
                gfac  = (optp%optpaer(iae)%bpr(ich, k + 1) - optp%optpaer(iae)%bpr(ich, k)) / delth
                IF (sun(j)) THEN
                  ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
                  pfac(1:coef_scatt_ir%fmv_aer_ph) = (optp%optpaer(iae)%pha(ich1, K + 1, 1:coef_scatt_ir%fmv_aer_ph)     &
                    &  - optp%optpaer(iae)%pha(ICH1, K, 1:coef_scatt_ir%fmv_aer_ph)) / delth
                  DO kk = coef_scatt_ir%fmv_aer_ph, 1,  - 1
                    frach     = frach + phash(kk) * pfac(kk)
                    phash(kk) = 0._jprb
                  ENDDO
                ENDIF
                frach                = frach + bparh * gfac
                frach                = frach + scach * sfac
                frach                = frach + absch * afac
                absch                = 0._jprb
                scach                = 0._jprb
                bparh                = 0._jprb
                aux_k%relhum(lay, j) = aux_k%relhum(lay, j) + frach
                frach                = 0._jprb
                EXIT
              ENDIF
            ENDDO
          ELSE
            IF (sun(j)) THEN
              ich1 = ich - coef_scatt_ir%fmv_aer_pha_ioff + 1
              phash(1:coef_scatt_ir%fmv_aer_ph) = 0._jprb
            ENDIF
            absch = 0._jprb
            scach = 0._jprb
            bparh = 0._jprb
          ENDIF
        ENDDO
      ENDDO
    ENDDO
  ENDIF
  IF (opts%addaerosl .OR. opts%addclouds) THEN
    transmission_scatt_ir_k%OPDPAAER(:,:)    = 0._jprb
    transmission_scatt_ir_k%OPDPSAER(:,:)    = 0._jprb
    transmission_scatt_ir_k%OPDPA(:,:)       = 0._jprb
    transmission_scatt_ir_k%OPDPS(:,:)       = 0._jprb
    transmission_scatt_ir_k%GPARAER(:,:)     = 0._jprb
    transmission_scatt_ir_k%GPAR(:,:)        = 0._jprb
    transmission_scatt_ir_k%OPDPAERLA(:,:)   = 0._jprb
    transmission_scatt_ir_k%GPARAERA(:,:)    = 0._jprb
    transmission_scatt_ir_k%azphaerup(:,:)   = 0._jprb
    transmission_scatt_ir_k%azphaerupa(:,:)  = 0._jprb
    transmission_scatt_ir_k%azphaerdo(:,:)   = 0._jprb
    transmission_scatt_ir_k%azphaerdoa(:,:)  = 0._jprb
    transmission_scatt_ir_k%azphupcls(:,:,:) = 0._jprb
    transmission_scatt_ir_k%azphdocls(:,:,:) = 0._jprb
    transmission_scatt_ir_k%azphdotot(:,:)   = 0._jprb
    transmission_scatt_ir_k%azphuptot(:,:)   = 0._jprb
  ENDIF
!-----Compute relative humidity-----------------------------------------------------------
  IF (opts%addaerosl) THEN
    DO j = 1, nchannels
      ipf = chanprof(j)%prof
      DO lay = nlayers, 1,  - 1
        lev = lay + 1
        IF (aux%RELHUMREF(lay, ipf) > 99._jprb) THEN
          aux_k%RELHUM(lay, j) = 0._jprb
        ENDIF
        ircld_k%wmixave(lay, j) = ircld_k%wmixave(lay, j) +                                         &
          & aux_k%RELHUM(lay, j) * 100._jprb * 1e-6_jprb * 0.622_jprb * ircld%xpresave(lay, ipf) /  &
          & (ircld%PPV(lay, ipf) * (0.622_jprb + ircld%wmixave(lay, ipf) * 1e-6_jprb * 0.622_jprb))
        ircld_k%wmixave(lay, j) = ircld_k%wmixave(lay, j) -                                                &
          & aux_k%RELHUM(lay, j) * 100._jprb * (1e-6_jprb * 0.622_jprb) ** 2 * ircld%xpresave(lay, ipf) *  &
          & ircld%wmixave(lay, ipf) * ircld%PPV(lay, ipf) /                                                &
          & (ircld%PPV(lay, ipf) * (0.622_jprb + ircld%wmixave(lay, ipf) * 1e-6_jprb * 0.622_jprb)) ** 2
        ircld_k%PPV(lay, j)     = ircld_k%PPV(lay, j) -                                                                    &
          & aux_k%RELHUM(lay, j) * 100._jprb * ircld%wmixave(lay, ipf) * 1e-6_jprb * 0.622_jprb * ircld%xpresave(lay, ipf) &
          &  * (0.622_jprb + ircld%wmixave(lay, ipf) * 1e-6_jprb * 0.622_jprb) /                                           &
          & (ircld%PPV(lay, ipf) * (0.622_jprb + ircld%wmixave(lay, ipf) * 1e-6_jprb * 0.622_jprb)) ** 2
        IF(opts%lgradp) ircld_k%xpresave(lay, j) = ircld_k%xpresave(lay, j) + &
          &    100._jprb * ircld%wmixave(lay, ipf) * 1.e-6_jprb * 0.622_jprb * aux_k%relhum(lay, j) /      &
          &   (ircld%ppv(lay, ipf) * (0.622_jprb + ircld%wmixave(lay, ipf) * 1.e-6_jprb * 0.622_jprb))
!---------Compute vater vapour partial pressure-------------------------------------------
        ircld_k%PPV(lay, j)     = ircld_k%PPV(lay, j) / 100._jprb
        IF (ircld%TAVE(lay, ipf) > T00) THEN
          ircld_k%ESW(lay, j) = ircld_k%ESW(lay, j) + ircld_k%PPV(lay, j)
        ELSE IF (ircld%TAVE(lay, ipf) > TI .AND. ircld%TAVE(lay, ipf) <= T00) THEN
          ircld_k%ESI(lay, j)  = ircld_k%ESI(lay, j) + ircld_k%PPV(lay, j)
          ircld_k%ESW(lay, j)  =      &
            & ircld_k%ESW(lay, j) + ircld_k%PPV(lay, j) * ((ircld%TAVE(lay, ipf) - TI) / (T00 - TI)) ** 2
          ircld_k%ESI(lay, j)  =      &
            & ircld_k%ESI(lay, j) - ircld_k%PPV(lay, j) * ((ircld%TAVE(lay, ipf) - TI) / (T00 - TI)) ** 2
          ircld_k%TAVE(lay, j) = ircld_k%TAVE(lay, j) +                                &
            & ircld_k%PPV(lay, j) * (ircld%ESW(lay, ipf) - ircld%ESI(lay, ipf)) * 2 *  &
            & ((ircld%TAVE(lay, ipf) - TI) / (T00 - TI) ** 2)
        ELSE IF (ircld%TAVE(lay, ipf) <= TI) THEN
          ircld_k%ESI(lay, j) = ircld_k%ESI(lay, j) + ircld_k%PPV(lay, j)
        ENDIF
        ircld_k%TAVE(lay, j) = ircld_k%TAVE(lay, j) +                                       &
          & ircld_k%ESW(lay, j) * ircld%ESW(lay, ipf) * 17.502_jprb * (T00 - 32.19_jprb) /  &
          & (ircld%TAVE(lay, ipf) - 32.19_jprb) ** 2
        ircld_k%TAVE(lay, j) = ircld_k%TAVE(lay, j) +                                     &
          & ircld_k%ESI(lay, j) * ircld%ESI(lay, ipf) * 22.587_jprb * (T00 + 0.7_jprb) /  &
          & (ircld%TAVE(lay, ipf) + 0.7_jprb) ** 2
!-----------------------------------------------------------------------------------------
        profiles_k(j)%q(lev - 1) = profiles_k(j)%q(lev - 1) + ircld_k%wmixave(lay, j) / 2._jprb
        profiles_k(j)%q(lev)     = profiles_k(j)%q(lev) + ircld_k%wmixave(lay, j) / 2._jprb
        profiles_k(j)%t(lev - 1) = profiles_k(j)%t(lev - 1) + ircld_k%tave(lay, j) / 2._jprb
        profiles_k(j)%t(lev)     = profiles_k(j)%t(lev) + ircld_k%tave(lay, j) / 2._jprb
        IF (opts%lgradp) THEN
          profiles_k(j)%p(lev - 1) = profiles_k(j)%p(lev - 1) + ircld_k%xpresave(lay, j) / 2._jprb
          profiles_k(j)%p(lev)     = profiles_k(j)%p(lev) + ircld_k%xpresave(lay, j) / 2._jprb
        ENDIF
      ENDDO
    ENDDO
  ENDIF
!-----------------------------------------------------------------------------------------
  opdpcldl    = 0._jprb
  opdpaerl    = 0._jprb
  opdpcldlsun = 0._jprb
  opdpaerlsun = 0._jprb
  IF (LHOOK) CALL DR_HOOK('RTTOV_OPDPSCATTIR_K', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_opdpscattir_k
