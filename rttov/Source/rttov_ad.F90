
SUBROUTINE rttov_ad( &
            & errorstatus,       &
            & chanprof,          &
            & opts,              &
            & profiles,          &
            & profiles_ad,       &
            & coefs,             &
            & calcemis,          &
            & emissivity,        &
            & emissivity_ad,     &
            & emissivity_out,    &
            & emissivity_out_ad, &
            & transmission,      &
            & transmission_ad,   &
            & radiancedata,      &
            & radiancedata_ad,   &
            & traj,              &
            & traj_ad,           &
            & pccomp,            &
            & pccomp_ad,         &
            & channels_rec)
!
! Description:
! Adjoint of rttov_direct
! to compute multi-channel level to space transmittances,
! top of atmosphere and level to space radiances and brightness
! temperatures and optionally surface emissivities, for many
! profiles in a single call, for satellite
! infrared or microwave sensors. The code requires a coefficient file
! for each sensor for which simulated radiances are requested.
!
! Note that radiancedata_ad can be used for all its structure elements
! In normal case the element total or bt is the only one initialised but
! for some particular cases like for rttov_cld_ad some other elements
! have been already init.
! According to the argument opts%switchrad the main input total or bt is used
! opts%switchrad == true    bt is the input, brightness temperature
! opts%switchrad == false   total is the input, radiance
!
! The AD outputs profiles_ad and emissivity_ad should be allocated and
! initialised before calling the subroutine
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
! Method: The methodology is described in the following:
!
! Eyre J.R. and H.M. Woolf  1988 Transmittance of atmospheric gases
! in the microwave region: a fast model. Applied Optics 27  3244-3249
!
! Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
! systems.  ECMWF Research Dept. Tech. Memo. 176 (available from the
! librarian at ECMWF).
!
! Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
! Transfer Model for Assimilation of Satellite Radiance Observations.
! QJRMS, 125, 1407-1425.
!
! Matricardi, M., F. Chevallier and S. Tjemkes 2001 An improved general
! fast radiative transfer model for the assimilation of radiance
! observations. ECMWF Research Dept. Tech. Memo. 345
! (available from the librarian at ECMWF).
!
! Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
! transfer model for the infrared atmospheric sounding interferometer.
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.1   01/12/2002  New F90 code with structures (P Brunel A Smith)
!  1.2   02/01/2003  More comments added (R Saunders)
!  1.3   24/01/2003  Error return code by input profile (P Brunel)
!  1.4               Add WV Continuum and CO2 capability
!  1.5   04/12/2003  Optimisation (J Hague and D Salmond ECMWF)
!  1.6   02/06/2004  Change tests on id_comp_lvl == 7 by tests on fmv_model_ver (P. Brunel)
!  1.7   08/12/2005  Added surface humidity to lowest mean layer q (R Saunders)
!  1.8   01/06/2005  Marco Matricardi (ECMWF):
!              --    IASI capability added.
!              --    Variable trace gases added for IASI and AIRS.
!              --    Solar radiation added for IASI and AIRS
!              --    Altitude dependent local zenith angle added.
!              --    Linear in tau approximation for RT equation added.
!  1.9    06/12/2005 Minor mods to initialise of traj0_ad%sunglint (R Saunders)
!  1.10   08/02/2007 Removed polarisation indexing (R Saunders)
!         29/06/2007 Introduce interpolator (P.J.Rayer)
!              --    Profiles on USER levs -> COEF levs
!              --    Duplicate profaux and setgeometry (USER and COEF levs)
!              --    Predict atm/sol optical depths on COEF levs -> USER levs
!              --    Calculate transmittance and integrate RT on USER levs
!  1.11   11/10/2007 Move iaernum & iaertyp profile members to traj0%aux_prof P.Marguinaud
!  1.12   16/01/2008 Facility to apply regression limits (N. Bormann)
!  1.13   26/09/2008 Disable vectorization of the pointer association
!                    loops and nullify temp pointers; this is for the NEC
!                    under -Chopt. P. Marguinaud
!
!  1.13   03/10/2008 Fix bug in checking addsolar (A. Geer)
!  1.14   15/08/2009 User defined ToA. Layers distinct from levels (P.Rayer)
!  1.15   03/11/2009 Transmittances / optical depths on levels (A Geer)
!  1.16   02/12/2009 Introduced principal component capability. Marco Matricardi. ECMWF
!  1.17   17/06/2010 Introduced optional spacetop flag to zero opdeps at user's
!                    model-top in Zeeman channels (P Rayer)
!  1.28   05/07/2010 Remove addsolar flag from profiles structure (J Hocking)
!  1.29   14/10/2010 Remove rt8_mode (J Hocking)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
! A user guide and technical documentation is available at
! http://www.metoffice.com/research/interproj/nwpsaf/rtm/index.html
!
! Declarations:
! Modules used:
! Imported Parameters:
!INTF_OFF
#include "throw.h"
!INTF_ON
! Imported Type Definitions:
  USE rttov_types, ONLY :  &
       & rttov_coefs,       &
       & rttov_pccomp,      &
       & profile_Type,      &
       & transmission_Type, &
       & radiance_Type,     &
       & rttov_options,     &
       & rttov_chanprof,    &
       & rttov_traj
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY :  &
       & sensor_id_mw, &
       & sensor_id_ir, &
       & sensor_id_hi, &
       & sensor_id_po, &
       & ncldtyp
  USE yomhook, ONLY : LHOOK, DR_HOOK
  USE rttov_types, ONLY :  &
       & rttov_traj_sta,    &
       & rttov_traj_dyn      
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(profile_Type  )   , INTENT(IN)                           :: profiles(:)
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof(:)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(profile_Type  )   , INTENT(INOUT)                        :: profiles_ad(size(profiles))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis         (size(chanprof))
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus
  REAL   (KIND=jprb)     , INTENT(IN)                           :: emissivity       (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_ad    (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(OUT)                          :: emissivity_out   (size(chanprof))
  REAL   (KIND=jprb)     , INTENT(INOUT)                        :: emissivity_out_ad(size(chanprof))
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission                     ! in because of meme allocation
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission_ad                  ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                     ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata_ad                  ! in because of meme allocation
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL     , TARGET :: traj           , traj_ad         ! target is needed here, because
                                                                                                    ! of rttov_check_temp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp_ad
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
!INTF_END
#include "rttov_checkinput_ad.h"
#include "rttov_errorreport.h"
#include "rttov_setgeometry_ad.h"
#include "rttov_calcemis_ir_ad.h"
#include "rttov_profaux_ad.h"
#include "rttov_setpredictors_7_ad.h"
#include "rttov_setpredictors_8_ad.h"
#include "rttov_setpredictors_9_ad.h"
#include "rttov_opdep_ad.h"
#include "rttov_opdep_9_ad.h"
#include "rttov_opdep_9_solar_ad.h"
#include "rttov_transmit_ad.h"
#include "rttov_transmit_9_solar_ad.h"
#include "rttov_calcemis_mw_ad.h"
#include "rttov_integrate_ad.h"
#include "rttov_intavg_prof_ad.h"
#include "rttov_intavg_chan_ad.h"
#include "rttov_refsun_ad.h"
#include "rttov_init_prof.h"
#include "rttov_check_traj.h"
#include "rttov_init_predictor.h"
#include "rttov_init_raytracing.h"
#include "rttov_setpredictors_9_solar_ad.h"
#include "rttov_cldstr_ad.h"
#include "rttov_opdpscattir_ad.h"
#include "rttov_fresnel_ad.h"
#include "rttov_pcscores_ad.h"
#include "rttov_reconstruct_ad.h"
#include "rttov_calcbt_pc_ad.h"
#include "rttov_add_raytracing.h"
#include "rttov_add_aux_prof.h"
#include "rttov_add_prof.h"
#include "rttov_init_sunglint.h"
#include "rttov_init_ircld.h"
#include "rttov_init_opdp_path.h"
#include "rttov_init_transmission_aux.h"
#include "rttov_add_opdp_path.h"
#include "rttov_init_aux_prof.h"
#include "rttov_init_auxrad_stream.h"
#include "rttov_init_trans_scatt_ir.h"
#include "rttov_alloc_traj_dyn.h"
#include "rttov_alloc_traj_sta.h"
#include "rttov_direct.h"
!local variables:
  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_ad_dyn

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_ad
  TYPE(rttov_traj), POINTER :: traj0_ad

  TYPE(rttov_options) :: opts_COEF
  INTEGER(KIND=jpim)  :: nlevels
  INTEGER(KIND=jpim)  :: nprofiles                    ! Number of profiles
  INTEGER(KIND=jpim)  :: nchannels                    ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim)  :: npcscores
  INTEGER(KIND=jpim)  :: ERR
  REAL(KIND=JPRB)     :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  TRY
  IF (LHOOK) CALL DR_HOOK('RTTOV_AD', 0_jpim, ZHOOK_HANDLE)
  nprofiles           = size(profiles)
  nchannels           = size(chanprof)
  opts_COEF           = opts
  opts_COEF%addaerosl = .FALSE.
  opts_COEF%addclouds = .FALSE.
  nlevels             = profiles(1)%nlevels
!-------------
!0. initialize
!-------------
  errorstatus         = errorstatus_success

  IF (opts%addpc) THEN
    npcscores = size(pccomp%pcscores)
  ENDIF

  NULLIFY (traj0)
  NULLIFY (traj0_ad)
!-----------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!-----------------------------------------------------------------------
! done immediately in direct code - for quick escape if target values are out of bounds
  CALL rttov_check_traj( &
        & ERR,                 &
        & nprofiles,           &
        & nchannels,           &
        & opts,                &
        & nlevels,             &
        & coefs,               &
        & 1_jpim,              &
        & traj0 = traj0,       &
        & traj0_ad = traj0_ad, &
        & traj1 = traj1,       &
        & traj1_ad = traj1_ad, &
        & traj2 = traj,        &
        & traj2_ad = traj_ad)
  THROWM( ERR .NE. 0 , "rttov_check_traj fatal error")

  CALL rttov_direct( &
            & errorstatus,               &
            & chanprof,                  &
            & opts,                      &
            & profiles,                  &
            & coefs,                     &
            & calcemis,                  &
            & emissivity,                &
            & emissivity_out,            &
            & transmission,              &
            & radiancedata,              &
            & traj         = traj0,      &
            & traj_dyn     = traj0_dyn,  &
            & traj_sta     = traj0_sta,  &
            & pccomp       = pccomp,     &
            & channels_rec = channels_rec)
  
  IF (errorstatus == errorstatus_fatal) THEN
    err = errorstatus_fatal
    THROW(err.ne.0)
  ENDIF

  CALL rttov_alloc_traj_dyn (err, traj0_ad_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(err.ne.0)
  
!
! Adjoint
!----------------
!
!---------------------------------------------
!0. allocate and initialize local AD variables
!---------------------------------------------
!0.1 profile AD-variables on COEF levels (USER levels done in tstrad_ad)
!------------------------------------------------------------------------
  CALL rttov_init_opdp_path(traj0_ad%opdp_path)
  CALL rttov_init_opdp_path(traj0_ad%opdp_path_COEF)

  CALL rttov_init_ircld(traj0_ad%ircld)
  IF (opts%addaerosl .OR. opts%addclouds) THEN
    CALL rttov_init_trans_scatt_ir(traj0_ad%transmission_scatt_ir)
    CALL rttov_init_trans_scatt_ir(traj0_ad_dyn%transmission_scatt_ir_stream)
  ENDIF

  CALL rttov_init_auxrad_stream (traj0_ad_dyn%auxrad_stream)

  CALL rttov_init_transmission_aux (traj0_ad_dyn%transmission_aux)

  CALL rttov_init_prof(traj0_ad%profiles_COEF)
!0.2 auxillary profile AD-variables
!----------------------------------
! on USER levels
  CALL rttov_init_aux_prof(traj0_ad%aux_prof)
! on COEF levels
  CALL rttov_init_aux_prof(traj0_ad%aux_prof_COEF)
!0.3 predictor AD-variables on COEF levels
!------------------------------------------
  CALL rttov_init_predictor(traj0_ad%predictors)
!0.4 raytracing AD-variables
!---------------------------
  CALL rttov_init_raytracing(traj0_ad%raytracing)
  CALL rttov_init_raytracing(traj0_ad%raytracing_COEF)
!0.5 path opdep AD-variables on COEF levels
!--------------------------------------------------
! allocated with the direct arrays
! on USER levels
!0.6 transmission AD-variables on USER levels
!--------------------------------------------
! allocated with the direct arrays
!0.7 emissivity and traj0%reflectivity AD-variables - SURFACE
!-------------------------------------------------------
  traj0_ad%fresnrefl(:)    = 0._JPRB
! emissivity_ad is init before calling
  traj0_ad%reflectivity(:) = 0._JPRB
!--------------------------------------------------
!1. AD of radiative transfer integration - USER levs
!--------------------------------------------------

  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)

  IF (opts%addpc) THEN

    IF (opts%addradrec) THEN

      IF (opts%switchrad) THEN
        CALL rttov_calcbt_pc_ad(         &
              & traj0_sta%chanprof_in,   &
              & coefs%coef_pccomp,       &
              & pccomp,                  &
              & pccomp_ad)
      ENDIF

      CALL rttov_reconstruct_ad( &
            & traj0_sta%chanprof_in,   &
            & traj0_sta%chanprof_pc,   &
            & pccomp,                  &
            & pccomp_ad,               &
            & coefs%coef_pccomp)
    ENDIF

    CALL rttov_pcscores_ad( &
          & opts,                   &
          & chanprof,               &
          & traj0_sta%chanprof_pc,  &
          & pccomp,                 &
          & pccomp_ad,              &
          & coefs%coef_pccomp,      &
          & radiancedata_ad)
  ENDIF


  IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
    CALL rttov_init_sunglint(traj0_ad%sunglint)
  ENDIF

  CALL rttov_integrate_ad( &
        & addcosmic,                                 &
        & opts,                                      &
        & traj0_dyn%nstreams,                        &
        & chanprof,                                  &
        & emissivity_out,                            &
        & emissivity_out_ad,                         &
        & traj0%reflectivity,                        &
        & traj0_ad%reflectivity,                     &
        & traj0%fresnrefl,                           &
        & traj0_ad%fresnrefl,                        &
        & traj0%sunglint,                            &
        & traj0_ad%sunglint,                         &
        & traj0_sta%sun,                             &
        & traj0_dyn%transmission_aux,                &
        & traj0_ad_dyn%transmission_aux,             &
        & traj0_dyn%transmission_scatt_ir_stream,    &
        & traj0_ad_dyn%transmission_scatt_ir_stream, &
        & profiles,                                  &
        & profiles_ad,                               &
        & traj0%aux_prof,                            &
        & traj0_ad%aux_prof,                         &
        & coefs%coef,                                &
        & traj0%raytracing,                          &
        & traj0_ad%raytracing,                       &
        & traj0%ircld,                               &
        & traj0_ad%ircld,                            &
        & radiancedata,                              &
        & traj0_sta%auxrad,                          &
        & traj0_dyn%auxrad_stream,                   &
        & traj0_ad_dyn%auxrad_stream,                &
        & radiancedata_ad)! inout  (output if converstion Bt -> rad)
!---------------------------------------
!2. AD of channel emissivities - SURFACE
!---------------------------------------

  IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
    CALL rttov_fresnel_ad( &
          & chanprof,           &
          & profiles,           &
          & coefs%coef,         &
          & traj0%sunglint,     &
          & traj0_ad%sunglint,  &
          & traj0%fresnrefl,    &
          & traj0_ad%fresnrefl)
    CALL rttov_refsun_ad( &
          & profiles,            &
          & profiles_ad,         &
          & coefs%coef,          &
          & traj0%aux_prof,      &
          & traj0%sunglint,      &
          & traj0_ad%sunglint,   &
          & traj0%raytracing,    &
          & traj0_ad%raytracing)
  ENDIF


  IF (Any(calcemis)) THEN
! calculate surface emissivity for selected channels
! and traj0%reflectivity

    IF (coefs%coef%id_sensor == sensor_id_ir) THEN
!Infrared
      emissivity_out_ad(:) =  - traj0_ad%reflectivity(:) + emissivity_out_ad(:)
    ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
!Microwave
      CALL rttov_calcemis_mw_ad( &
            & opts,                               &
            & profiles,                           &
            & profiles_ad,                        &
            & traj0_sta%angles,                   &
            & coefs%coef,                         &
            & chanprof,                           &
            & traj0_dyn%transmission_aux,         &
            & traj0_ad_dyn%transmission_aux,      &
            & calcemis,                           &
            & emissivity_out_ad,                  &
            & traj0_ad%reflectivity)
    ELSE
! Hires
      emissivity_out_ad(:) =  - traj0_ad%reflectivity(:) + emissivity_out_ad(:)
      CALL rttov_calcemis_ir_ad( &
            & profiles,          &
            & profiles_ad,       &
            & coefs%coef,        &
            & opts%addpc,        &
            & coefs%coef_pccomp, &
            & chanprof,          &
            & calcemis,          &
            & emissivity_out_ad)
    ENDIF

  ENDIF


  WHERE (.NOT. calcemis)
    emissivity_out_ad = emissivity_out_ad - traj0_ad%reflectivity
  ENDWHERE

! Where( .not. calcemis )
!   emissivity_ad = emissivity_ad + emissivity_out_ad
! End Where
! Where( .not. calcemis )
!   emissivity_out_ad = 0._jprb
! End Where
  emissivity_ad = emissivity_ad + emissivity_out_ad
!------------------------------------
!3. AD of transmittances - USER levs
!------------------------------------

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_transmit_9_solar_ad( &
            & opts%addaerosl,                            &
            & opts%addclouds,                            &
            & profiles(1)%nlayers,                       &
            & chanprof,                                  &
            & profiles,                                  &
            & traj0_sta%sun,                             &
            & traj0%aux_prof,                            &
            & traj0_ad%aux_prof,                         &
            & coefs%coef,                                &
            & traj0%raytracing,                          &
            & traj0_ad%raytracing,                       &
            & traj0%ircld,                               &
            & traj0%opdp_path,                           &
            & traj0_ad%opdp_path,                        &
            & traj0_sta%odsun_level,                     &
            & traj0_sta%odsun_singlelayer,               &
            & traj0_sta%od_frac,                         &
            & traj0_dyn%transmission_aux,                &
            & traj0_ad_dyn%transmission_aux,             &
            & traj0_dyn%transmission_scatt_ir_stream,    &
            & traj0_ad_dyn%transmission_scatt_ir_stream, &
            & traj0_sta%tausun_ref,                      &
            & traj0_sta%tausun_ref_surf,                 &
            & traj0_sta%tausun_level,                    &
            & traj0_sta%tausun_surf)
    ENDIF
  ENDIF

  CALL rttov_transmit_ad( &
        & opts%addaerosl,                            &
        & opts%addclouds,                            &
        & profiles(1)%nlayers,                       &
        & chanprof,                                  &
        & traj0%aux_prof,                            &
        & traj0_ad%aux_prof,                         &
        & coefs%coef,                                &
        & traj0%ircld,                               &
        & traj0%opdp_path,                           &
        & traj0_ad%opdp_path,                        &
        & traj0_sta%od_level,                        &
        & transmission,                              &
        & transmission_ad,                           &
        & traj0_dyn%transmission_aux,                &
        & traj0_ad_dyn%transmission_aux,             &
        & traj0_dyn%transmission_scatt_ir_stream,    &
        & traj0_ad_dyn%transmission_scatt_ir_stream, &
        & traj0_sta%tau_ref,                         &
        & traj0_sta%tau_ref_surf,                    &
        & traj0_sta%tau_surf,                        &
        & traj0_sta%tau_level)
!-------------------------------------------------------------
!4. AD of optical depths of aerosols and/or clouds - USER levs
!-------------------------------------------------------------

  IF (opts%addaerosl .OR. opts%addclouds) THEN
    CALL rttov_opdpscattir_ad( &
          & profiles(1)%nlayers,                       &
          & chanprof,                                  &
          & opts,                                      &
          & traj0%aux_prof,                            &
          & traj0_ad%aux_prof,                         &
          & profiles,                                  &
          & profiles_ad,                               &
          & traj0_sta%sun,                             &
          & coefs%coef,                                &
          & coefs%coef_scatt_ir,                       &
          & traj0%raytracing,                          &
          & traj0_ad%raytracing,                       &
          & traj0%transmission_scatt_ir,               &
          & traj0_ad%transmission_scatt_ir,            &
          & traj0_dyn%transmission_scatt_ir_stream,    &
          & traj0_ad_dyn%transmission_scatt_ir_stream, &
          & coefs%optp,                                &
          & traj0%ircld,                               &
          & traj0_ad%ircld)
  ENDIF

!----------------------------------------------------
!5. AD of cloud streams and distributions - USER levs
!----------------------------------------------------

  IF (opts%addclouds) THEN
    CALL rttov_cldstr_ad( &
          & profiles,       &
          & profiles_ad,    &
          & traj0%ircld,    &
          & traj0_ad%ircld)
  ENDIF

!------------------------------------
!6. AD of optical depth interpolation
!------------------------------------
 
  if(opts%SpaceTop) traj0_ad%opdp_path%atm_level(1,:) = 0._jprb
  if(opts%SpaceTop) traj0_ad%opdp_path%sun_level(1,:) = 0._jprb

  IF (opts%addinterp) THEN
    CALL rttov_intavg_chan_ad( &
          & opts%lgradp,                    &
          & traj0_sta%sun,                  &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_COEF,            &
          & profiles,                       &
          & profiles_ad,                    &
          & traj0%opdp_path_COEF,           &
          & traj0_ad%opdp_path_COEF,        &
          & traj0%opdp_path,                &
          & traj0_ad%opdp_path)
  ELSE
    CALL rttov_add_opdp_path(traj0_ad%opdp_path_COEF, traj0_ad%opdp_path_COEF, traj0_ad%opdp_path)
  ENDIF

!! testop omitted
!---------------------------------------------------------
!7. AD of atmospheric and solar optical depths - COEF levs
!---------------------------------------------------------
! optical depth arrays allocated at start of AD-code

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_opdep_9_solar_ad( &
            & coefs%coef%nlayers,      &
            & chanprof,                &
            & profiles,                &
            & traj0_sta%sun,           &
            & traj0%predictors,        &
            & traj0_ad%predictors,     &
            & coefs%coef,              &
            & traj0%opdp_path_COEF,    &
            & traj0_ad%opdp_path_COEF, &
            & traj0_sta%opdpsun_ref_COEF)
    ENDIF
    CALL rttov_opdep_9_ad( &
          & coefs%coef%nlayers,      &
          & chanprof,                &
          & traj0%predictors,        &
          & traj0_ad%predictors,     &
          & traj0%aux_prof_COEF,     &
          & traj0_ad%aux_prof_COEF,  &
          & coefs%coef,              &
          & traj0%opdp_path_COEF,    &
          & traj0_ad%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
  ELSE
    CALL rttov_opdep_ad( &
          & coefs%coef%nlayers,      &
          & chanprof,                &
          & traj0%predictors,        &
          & traj0_ad%predictors,     &
          & traj0%aux_prof_COEF,     &
          & traj0_ad%aux_prof_COEF,  &
          & coefs%coef,              &
          & traj0%opdp_path_COEF,    &
          & traj0_ad%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
  ENDIF

!-------------------------------------------------------
!8. AD of RTTOV-7 RTTOV-8 RTTOV-9 predictors - COEF levs
!-------------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7_ad( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_ad%profiles_COEF,   &
          & traj0_sta%angles_COEF,    &
          & coefs%coef,               &
          & traj0%aux_prof_COEF,      &
          & traj0%predictors,         &
          & traj0_ad%predictors,      &
          & traj0%raytracing_COEF,    &
          & traj0_ad%raytracing_COEF)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8_ad( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_ad%profiles_COEF,   &
          & traj0_sta%angles_COEF,    &
          & coefs%coef,               &
          & traj0%aux_prof_COEF,      &
          & traj0%predictors,         &
          & traj0_ad%predictors,      &
          & traj0%raytracing_COEF,    &
          & traj0_ad%raytracing_COEF)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_setpredictors_9_solar_ad( &
            & opts,                     &
            & traj0%profiles_COEF,      &
            & traj0_ad%profiles_COEF,   &
            & coefs%coef,               &
            & traj0%predictors,         &
            & traj0_ad%predictors,      &
            & traj0%raytracing_COEF,    &
            & traj0_ad%raytracing_COEF)
    ENDIF
    CALL rttov_setpredictors_9_ad( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_ad%profiles_COEF,   &
          & traj0_sta%angles_COEF,    &
          & coefs%coef_pccomp,        &
          & coefs%coef,               &
          & traj0%predictors,         &
          & traj0_ad%predictors,      &
          & traj0%raytracing_COEF,    &
          & traj0_ad%raytracing_COEF)
  ENDIF

!-------------------------------------------------------
!9. AD of common geometric set-up for the RT integration
!-------------------------------------------------------
!9.1 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_setgeometry_ad( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_ad%profiles_COEF,   &
          & traj0%aux_prof_COEF,      &
          & coefs%coef,               &
          & traj0_sta%angles_COEF,    &
          & traj0%raytracing_COEF,    &
          & traj0_ad%raytracing_COEF)
  ELSE
    traj0_sta%angles_COEF = traj0_sta%angles
    CALL rttov_add_raytracing(traj0_ad%raytracing, traj0_ad%raytracing, traj0_ad%raytracing_COEF)
  ENDIF

!9.1 USER levs
!-------------
  CALL rttov_setgeometry_ad( &
        & opts,                &
        & profiles,            &
        & profiles_ad,         &
        & traj0%aux_prof,      &
        & coefs%coef,          &
        & traj0_sta%angles,    &
        & traj0%raytracing,    &
        & traj0_ad%raytracing)
!---------------------------------------------------------
!10. AD of cloud top, surface levels, ice cloud parameters
!---------------------------------------------------------
!10.1 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_profaux_ad( &
          & opts_COEF,              &
          & traj0%profiles_COEF,    &
          & traj0_ad%profiles_COEF, &
          & coefs%coef,             &
          & traj0%aux_prof_COEF,    &
          & traj0_ad%aux_prof_COEF)
  ELSE
    CALL rttov_add_aux_prof(traj0_ad%aux_prof, traj0_ad%aux_prof, traj0_ad%aux_prof_COEF)
  ENDIF

!10.2 USER levs
!-------------
  CALL rttov_profaux_ad( &
        & opts,              &
        & profiles,          &
        & profiles_ad,       &
        & coefs%coef,        &
        & traj0%aux_prof,    &
        & traj0_ad%aux_prof)
!------------------------------------------------------------------
!11 AD of check input data is within suitable physical limits - COEF levs
!------------------------------------------------------------------

  IF (opts%apply_reg_limits) THEN

    CALL rttov_checkinput_ad( &
          & opts,                        &
          & traj0_sta%profiles_COEF_ref, &
          & traj0_ad%profiles_COEF,      &
          & coefs%coef,                  &
          & coefs%coef_pccomp)

  ENDIF

!-----------------------------------------------------------------
!12. AD of profile variable interpolation - COEF levs to USER levs
!-----------------------------------------------------------------
! intvar logical array set by direct code

  IF (opts%addinterp) THEN
    CALL rttov_intavg_prof_ad( &
          & opts,                           &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & profiles_ad,                    &
          & traj0%profiles_COEF,            &
          & traj0_ad%profiles_COEF)!inout  target variables
  ELSE
! COEF levs same as USER levs
    CALL rttov_add_prof( &
          & profiles_ad,            &
          & profiles_ad,            &
          & traj0_ad%profiles_COEF, &
          & lair = .TRUE._jplm,     &
          & lground = .FALSE._jplm)
  ENDIF

  CALL rttov_add_prof( &
        & profiles_ad,            &
        & profiles_ad,            &
        & traj0_ad%profiles_COEF, &
        & lair = .FALSE._jplm,    &
        & lground = .TRUE._jplm)
!-----------------------------------------
!13. deallocate memory for local variables
!-----------------------------------------

  CALL rttov_alloc_traj_dyn (err, traj0_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
  THROW(err.ne.0)
!
  CALL rttov_alloc_traj_dyn (err, traj0_ad_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
  THROW(err.ne.0)
!
  CALL rttov_alloc_traj_sta (err, traj0_sta, opts, coefs%coef, nlevels, nchannels, nprofiles, 0_jpim, npcscores, channels_rec)
  THROW(err.ne.0)
!

  CALL rttov_check_traj( &
        & ERR,                 &
        & nprofiles,           &
        & nchannels,           &
        & opts,                &
        & nlevels,             &
        & coefs,               &
        & 0_jpim,              &
        & traj0 = traj0,       &
        & traj0_ad = traj0_ad, &
        & traj1 = traj1,       &
        & traj1_ad = traj1_ad, &
        & traj2 = traj,        &
        & traj2_ad = traj_ad)

  THROWM( ERR .NE. 0 , "deallocation check_traj")

  IF (LHOOK) CALL DR_HOOK('RTTOV_AD', 1_jpim, ZHOOK_HANDLE)
  CATCH
  errorstatus = ERR
  IF (LHOOK) CALL DR_HOOK('RTTOV_AD', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_ad
