!
SUBROUTINE rttov_tl( &
            & errorstatus,       &
            & chanprof,          &
            & opts,              &
            & profiles,          &
            & profiles_tl,       &
            & coefs,             &
            & calcemis,          &
            & emissivity,        &
            & emissivity_tl,     &
            & emissivity_out,    &
            & emissivity_out_tl, &
            & transmission,      &
            & transmission_tl,   &
            & radiancedata,      &
            & radiancedata_tl,   &
            & traj,              &
            & traj_tl,           &
            & pccomp,            &
            & pccomp_tl,         &
            & channels_rec)
!
! Description:
! Tangent Linear of rttov_direct
! to compute multi-channel level to space transmittances,
! top of atmosphere and level to space radiances and brightness
! temperatures and optionally surface emissivities, for many
! profiles in a single call, for satellite
! infrared or microwave sensors. The code requires a coefficient file
! for each sensor for which simulated radiances are requested.
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
! Rochon Y.J., L. Garand, D.S. Turner and S. Polavarapu Jacobian mapping
! between  vertical coordinate systems in data assimilation (submitted QJRMS
! June 2006)
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
!  1.7   08/12/2005  Added surface humidity to lowest mean layer q  (R Saunders)
!  1.8   23/01/2006  Marco Matricardi (ECMWF):
!              --    IASI capability added.
!              --    Variable trace gases added for IASI and AIRS.
!              --    Solar radiation added for IASI and AIRS
!              --    Altitude dependent local zenith angle added.
!              --    Linear in tau approximation for RT equation added.
!  1.9   23/01/2006 Changed intent of profiles_tl to inout (R Saunders)
!  1.10  31/03/2006  Initialised zlev (R Saunders)
!        06/02/2007  Removed polarisation indexing R. Saunders.
!        29/06/2007  Introduce interpolator (P.J.Rayer)
!              --    Profiles on USER levs -> COEF levs
!              --    Duplicate profaux and setgeometry (USER and COEF levs)
!              --    Predict atm/sol optical depths on COEF levs -> USER levs
!              --    Calculate transmittance and integrate RT on USER levs
!  1.11  11/10/2007  Move iaernum & iaertyp to profile_aux P.Marguinaud
!                    Initialize weights_calc array
!  1.12   16/01/2008 Facility to apply regression limits (N. Bormann)
!  1.13   26/09/2008 Disable vectorization of the pointer association
!                    loop and nullify temp pointers; this is for the NEC
!                    under -Chopt. P. Marguinaud
!  1.12   03/10/2008 Fix bug in checking addsolar
!  1.13   03/07/2009 Profile levels to include ToA. Distinguish between layer
!                    arrays  and level arrays - size, index labels, looping
!                    (P. Rayer)
!  1.14   03/11/2009 Transmittances / optical depths on levels (A Geer)
!  1.15   02/12/2009 Introduced principal component capability. Marco Matricardi. ECMWF
!  1.16   17/06/2010 Introduced optional spacetop flag to zero opdeps at user's
!                    model-top in Zeeman channels (P Rayer)
!  1.17   05/07/2010 Remove addsolar flag from profiles structure (J Hocking)
!  1.18    14/10/2010 Remove rt8_mode (J Hocking)
!
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
       & rttov_options,     &
       & radiance_Type,     &
       & transmission_Type, &
       & profile_Type,      &
       & rttov_coefs,       &
       & rttov_pccomp,      &
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
  USE rttov_types, ONLY : &
       & rttov_traj_dyn,    &
       & rttov_traj_sta
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(profile_Type)     , INTENT(IN)                           :: profiles   (:)
  INTEGER(KIND=jpim)     , INTENT(OUT)                          :: errorstatus
  TYPE(rttov_chanprof)   , INTENT(IN)                           :: chanprof   (:)
  TYPE(rttov_options )   , INTENT(IN)                           :: opts
  TYPE(profile_Type  )   , INTENT(INOUT)                        :: profiles_tl(size(profiles))
  TYPE(rttov_coefs   )   , INTENT(IN)   , TARGET                :: coefs
  LOGICAL(KIND=jplm)     , INTENT(IN)                           :: calcemis         (size(chanprof))
  REAL(KIND=jprb)        , INTENT(IN)                           :: emissivity       (size(chanprof))
  REAL(KIND=jprb)        , INTENT(IN)                           :: emissivity_tl    (size(chanprof))
  REAL(KIND=jprb)        , INTENT(OUT)                          :: emissivity_out   (size(chanprof))
  REAL(KIND=jprb)        , INTENT(OUT)                          :: emissivity_out_tl(size(chanprof))
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission                     ! in because of meme allocation
  TYPE(transmission_Type), INTENT(INOUT)                        :: transmission_tl                  ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata                     ! in because of meme allocation
  TYPE(radiance_Type    ), INTENT(INOUT)                        :: radiancedata_tl                  ! in because of meme allocation
  TYPE(rttov_traj       ), INTENT(INOUT), OPTIONAL     , TARGET :: traj           , traj_tl         ! target is needed here
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp
  TYPE(rttov_pccomp     ), OPTIONAL     , INTENT(INOUT)         :: pccomp_tl
  INTEGER(KIND=jpim)     , OPTIONAL     , INTENT(IN)            :: channels_rec(:)
!INTF_END
#include "rttov_errorreport.h"
#include "rttov_checkinput_tl.h"
#include "rttov_calcemis_ir_tl.h"
#include "rttov_profaux_tl.h"
#include "rttov_setgeometry_tl.h"
#include "rttov_setpredictors_7_tl.h"
#include "rttov_setpredictors_8_tl.h"
#include "rttov_setpredictors_9_tl.h"
#include "rttov_opdep_tl.h"
#include "rttov_opdep_9_tl.h"
#include "rttov_opdep_9_solar_tl.h"
#include "rttov_transmit_tl.h"
#include "rttov_transmit_9_solar_tl.h"
#include "rttov_calcemis_mw_tl.h"
#include "rttov_integrate_tl.h"
#include "rttov_intavg_prof_tl.h"
#include "rttov_intavg_chan_tl.h"
#include "rttov_refsun_tl.h"
#include "rttov_check_traj.h"
#include "rttov_init_prof.h"
#include "rttov_setpredictors_9_solar_tl.h"
#include "rttov_cldstr_tl.h"
#include "rttov_opdpscattir_tl.h"
#include "rttov_fresnel_tl.h"
#include "rttov_init_raytracing.h"
#include "rttov_pcscores_tl.h"
#include "rttov_reconstruct_tl.h"
#include "rttov_calcbt_pc_tl.h"
#include "rttov_copy_raytracing.h"
#include "rttov_copy_prof.h"
#include "rttov_copy_aux_prof.h"
#include "rttov_copy_opdp_path.h"
#include "rttov_init_opdp_path.h"
#include "rttov_init_aux_prof.h"
#include "rttov_alloc_traj_dyn.h"
#include "rttov_alloc_traj_sta.h"
#include "rttov_direct.h"
!local variables:
  INTEGER(KIND=jpim)               :: nlevels
  LOGICAL(KIND=jplm) :: addcosmic! switch for adding temp of cosmic background

  TYPE(rttov_traj_dyn) :: traj0_dyn
  TYPE(rttov_traj_sta) :: traj0_sta
  TYPE(rttov_traj_dyn) :: traj0_tl_dyn

  TYPE(rttov_traj), TARGET  :: traj1
  TYPE(rttov_traj), POINTER :: traj0
  TYPE(rttov_traj), TARGET  :: traj1_tl
  TYPE(rttov_traj), POINTER :: traj0_tl
  INTEGER(KIND=jpim) :: nprofiles                    ! Number of profiles
  INTEGER(KIND=jpim) :: nchannels                    ! Number of radiances computed (channels used * profiles)
  INTEGER(KIND=jpim)  :: npcscores
  TYPE(rttov_options) :: opts_COEF
  INTEGER(KIND=jpim)  :: ERR
  REAL   (KIND=JPRB)  :: ZHOOK_HANDLE

!- End of header --------------------------------------------------------
  TRY
!-------------
!0. initialize
!-------------
  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 0_jpim, ZHOOK_HANDLE)
  nprofiles           = size(profiles)
  nchannels           = size(chanprof)
  nlevels             = profiles(1)%nlevels
  opts_COEF           = opts
  opts_COEF%addclouds = .FALSE.
  opts_COEF%addaerosl = .FALSE.
  errorstatus         = errorstatus_success

  IF (opts%addpc) THEN
    npcscores = size(pccomp%pcscores)
  ENDIF

  NULLIFY (traj0)
  NULLIFY (traj0_tl)
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
        & traj0_tl = traj0_tl, &
        & traj1 = traj1,       &
        & traj1_tl = traj1_tl, &
        & traj2 = traj,        &
        & traj2_tl = traj_tl)
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

  CALL rttov_alloc_traj_dyn (err, traj0_tl_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 1_jpim)
  THROW(err.ne.0)
!
! Tangent Linear
!----------------
!
!---------------------------------------------------------------------------------------
!1. interpolator first call - input profiles from USER levs to COEF levs
!---------------------------------------------------------------------------------------
  CALL rttov_init_prof(traj0_tl%profiles_COEF)
! intvar logical array set by direct code

  IF (opts%addinterp) THEN
    CALL rttov_intavg_prof_tl( &
          & opts,                           &
          & nlevels,                        &
          & coefs%coef%nlevels,             &
          & profiles,                       &
          & profiles_tl,                    &
          & traj0%profiles_COEF,            &
          & traj0_tl%profiles_COEF)!inout  target variables
  ELSE
    CALL rttov_copy_prof( &
          & traj0_tl%profiles_COEF, &
          & profiles_tl,            &
          & larray = .TRUE._jplm,   &
          & lscalar = .FALSE._jplm)
  ENDIF

! complete profiles on COEF levs for surf, skin, cloud, aerosol, solar, angle
  CALL rttov_copy_prof( &
        & traj0_tl%profiles_COEF, &
        & profiles_tl,            &
        & larray = .FALSE._jplm,  &
        & lscalar = .TRUE._jplm)
  CALL rttov_init_raytracing(traj0_tl%raytracing)
!------------------------------------------------------------------
!2. check input data is within suitable physical limits - COEF levs
!------------------------------------------------------------------

  IF (opts%apply_reg_limits) THEN

    CALL rttov_checkinput_tl( &
          & opts,                       &
          & traj0_sta%profiles_COEF_ref,&
          & traj0_tl%profiles_COEF,     &
          & coefs%coef,                 &
          & coefs%coef_pccomp)

  ENDIF

!------------------------------------------------------------------------
!3. determine cloud top, surface levels, ice cloud parameters
!------------------------------------------------------------------------
!3.1 USER levs
!-------------
  CALL rttov_init_aux_prof(traj0_tl%aux_prof)
  CALL rttov_profaux_tl( &
        & opts,              &
        & profiles,          &
        & profiles_tl,       &
        & coefs%coef,        &
        & traj0%aux_prof,    &
        & traj0_tl%aux_prof)
  CALL rttov_init_aux_prof(traj0_tl%aux_prof_COEF)

  THROW(err.ne.0)

!3.2 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_profaux_tl( &
          & opts_COEF,              &
          & traj0%profiles_COEF,    &
          & traj0_tl%profiles_COEF, &
          & coefs%coef,             &
          & traj0%aux_prof_COEF,    &
          & traj0_tl%aux_prof_COEF)
  ELSE
    CALL rttov_copy_aux_prof(traj0_tl%aux_prof_COEF, traj0_tl%aux_prof)
  ENDIF

! TL on geometry
!------------------------------------------------------------------------
!4. set up common geometric variables for the RT integration
!------------------------------------------------------------------------
!4.1 USER levs
!-------------
  CALL rttov_setgeometry_tl( &
        & opts,                &
        & profiles,            &
        & profiles_tl,         &
        & traj0%aux_prof,      &
        & coefs%coef,          &
        & traj0_sta%angles,    &
        & traj0%raytracing,    &
        & traj0_tl%raytracing)
!4.2 COEF levs
!-------------

  IF (opts%addinterp) THEN
    CALL rttov_setgeometry_tl( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_tl%profiles_COEF,   &
          & traj0%aux_prof_COEF,      &
          & coefs%coef,               &
          & traj0_sta%angles_COEF,    &
          & traj0%raytracing_COEF,    &
          & traj0_tl%raytracing_COEF)
  ELSE
    traj0_sta%angles_COEF = traj0_sta%angles
    CALL rttov_copy_raytracing(traj0_tl%raytracing_COEF, traj0_tl%raytracing)
  ENDIF

! TL of predictors
!---------------------------------------------------
!5. calculate atm/sol predictors - atm/sol COEF levs
!---------------------------------------------------

  IF (coefs%coef%fmv_model_ver == 7) THEN
    CALL rttov_setpredictors_7_tl( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_tl%profiles_COEF,   &
          & traj0_sta%angles_COEF,    &
          & coefs%coef,               &
          & traj0%aux_prof_COEF,      &
          & traj0%predictors,         &
          & traj0_tl%predictors,      &
          & traj0%raytracing_COEF,    &
          & traj0_tl%raytracing_COEF)
  ELSE IF (coefs%coef%fmv_model_ver == 8) THEN
    CALL rttov_setpredictors_8_tl( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_tl%profiles_COEF,   &
          & traj0_sta%angles_COEF,    &
          & coefs%coef,               &
          & traj0%aux_prof_COEF,      &
          & traj0%predictors,         &
          & traj0_tl%predictors,      &
          & traj0%raytracing_COEF,    &
          & traj0_tl%raytracing_COEF)
  ELSE IF (coefs%coef%fmv_model_ver == 9) THEN
    CALL rttov_setpredictors_9_tl( &
          & opts,                     &
          & traj0%profiles_COEF,      &
          & traj0_tl%profiles_COEF,   &
          & traj0_sta%angles_COEF,    &
          & coefs%coef_pccomp,        &
          & coefs%coef,               &
          & traj0%predictors,         &
          & traj0_tl%predictors,      &
          & traj0%raytracing_COEF,    &
          & traj0_tl%raytracing_COEF)
    IF (opts%addsolar) THEN
      CALL rttov_setpredictors_9_solar_tl( &
            & opts,                     &
            & traj0%profiles_COEF,      &
            & traj0_tl%profiles_COEF,   &
            & coefs%coef,               &
            & traj0%predictors,         &
            & traj0_tl%predictors,      &
            & traj0%raytracing_COEF,    &
            & traj0_tl%raytracing_COEF)
    ENDIF
  ENDIF

!TL of optical depths
!--------------------------------------------------------------------
!6. predict atmospheric and solar optical depths - atm/sol COEF levs
!--------------------------------------------------------------------
  CALL rttov_init_opdp_path(traj0_tl%opdp_path)
  CALL rttov_init_opdp_path(traj0_tl%opdp_path_COEF)

  IF (coefs%coef%fmv_model_ver == 9) THEN
    CALL rttov_opdep_9_tl( &
          & coefs%coef%nlayers,      &
          & chanprof,                &
          & traj0%predictors,        &
          & traj0_tl%predictors,     &
          & traj0%aux_prof_COEF,     &
          & traj0_tl%aux_prof_COEF,  &
          & coefs%coef,              &
          & traj0%opdp_path_COEF,    &
          & traj0_tl%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
    IF (opts%addsolar) THEN
      CALL rttov_opdep_9_solar_tl( &
            & coefs%coef%nlayers,      &
            & chanprof,                &
            & profiles,                &
            & traj0_sta%sun,           &
            & traj0%predictors,        &
            & traj0_tl%predictors,     &
            & coefs%coef,              &
            & traj0%opdp_path_COEF,    &
            & traj0_tl%opdp_path_COEF, &
            & traj0_sta%opdpsun_ref_COEF)
    ENDIF
  ELSE
    CALL rttov_opdep_tl( &
          & coefs%coef%nlayers,      &
          & chanprof,                &
          & traj0%predictors,        &
          & traj0_tl%predictors,     &
          & traj0%aux_prof_COEF,     &
          & traj0_tl%aux_prof_COEF,  &
          & coefs%coef,              &
          & traj0%opdp_path_COEF,    &
          & traj0_tl%opdp_path_COEF, &
          & traj0_sta%opdp_ref_COEF)
  ENDIF

!--------------------------------------------------------------------------
!7. interpolator second  call - optical depths from COEF levs to USER levs
!--------------------------------------------------------------------------

  IF (opts%addinterp) THEN
    CALL rttov_intavg_chan_tl( &
          & opts%lgradp,                    &
          & traj0_sta%sun,                  &
          & coefs%coef%nlevels,             &
          & nlevels,                        &
          & chanprof,                       &
          & traj0%profiles_COEF,            &
          & profiles,                       &
          & profiles_tl,                    &
          & traj0%opdp_path_COEF,           &
          & traj0_tl%opdp_path_COEF,        &
          & traj0%opdp_path,                &
          & traj0_tl%opdp_path)
  ELSE
    CALL rttov_copy_opdp_path(traj0_tl%opdp_path, traj0_tl%opdp_path_COEF)
  ENDIF

  if(opts%SpaceTop) traj0_tl%opdp_path%atm_level(1,:) = 0._jprb
  if(opts%SpaceTop) traj0_tl%opdp_path%sun_level(1,:) = 0._jprb

!! testop omitted
! TL of ircld
!--------------------------------------------------------------------
!8. If clouds are present,calculate the number of streams and
!   the cloud distribution  in each stream - remains on USER levs
!--------------------------------------------------------------------

  IF (opts%addclouds) THEN
    CALL rttov_cldstr_tl( &
          & profiles,       &
          & profiles_tl,    &
          & traj0%ircld,    &
          & traj0_tl%ircld)
  ELSE
    traj0_tl%ircld%XSTRCLR = 0._jprb
    traj0_tl%ircld%NSTREAM = 0_jpim
  ENDIF

! TL of opdpscattir
!----------------------------------------------------------------------------
!9. Calculate optical depths of aerosols and/or clouds - USER levs
!----------------------------------------------------------------------------

  IF (opts%addaerosl .OR. opts%addclouds) THEN

    CALL rttov_opdpscattir_tl( &
          & profiles(1)%nlayers,                           &
          & chanprof,                                      &
          & opts,                                          &
          & traj0%aux_prof,                                &
          & traj0_tl%aux_prof,                             &
          & profiles,                                      &
          & profiles_tl,                                   &
          & traj0_sta%sun,                                 &
          & coefs%coef,                                    &
          & coefs%coef_scatt_ir,                           &
          & traj0%raytracing,                              &
          & traj0_tl%raytracing,                           &
          & traj0%transmission_scatt_ir,                   &
          & traj0_tl%transmission_scatt_ir,                &
          & traj0_dyn%transmission_scatt_ir_stream,        &
          & traj0_tl_dyn%transmission_scatt_ir_stream,     &
          & coefs%optp,                                    &
          & traj0%ircld,                                   &
          & traj0_tl%ircld)
  ENDIF

! TL of transmit
!----------------------------------------
!10. calculate transmittances - USER levs
!----------------------------------------

!TL of optical depths and transmittances
  CALL rttov_transmit_tl( &
        & opts%addaerosl,                                 &
        & opts%addclouds,                                 &
        & profiles(1)%nlayers,                            &
        & chanprof,                                       &
        & traj0%aux_prof,                                 &
        & traj0_tl%aux_prof,                              &
        & coefs%coef,                                     &
        & traj0%ircld,                                    &
        & traj0%opdp_path,                                &
        & traj0_tl%opdp_path,                             &
        & traj0_sta%od_level,                             &
        & transmission,                                   &
        & transmission_tl,                                &
        & traj0_dyn%transmission_aux,                     &
        & traj0_tl_dyn%transmission_aux,                  &
        & traj0_dyn%transmission_scatt_ir_stream,         &
        & traj0_tl_dyn%transmission_scatt_ir_stream,      &
        & traj0_sta%tau_ref,                              &
        & traj0_sta%tau_ref_surf,                         &
        & traj0_sta%tau_surf,                             &
        & traj0_sta%tau_level)

  IF (coefs%coef%fmv_model_ver == 9) THEN
    IF (opts%addsolar) THEN
      CALL rttov_transmit_9_solar_tl( &
            & opts%addaerosl,                                &
            & opts%addclouds,                                &
            & profiles(1)%nlayers,                           &
            & chanprof,                                      &
            & profiles,                                      &
            & traj0_sta%sun,                                 &
            & traj0%aux_prof,                                &
            & traj0_tl%aux_prof,                             &
            & coefs%coef,                                    &
            & traj0%raytracing,                              &
            & traj0_tl%raytracing,                           &
            & traj0%ircld,                                   &
            & traj0%opdp_path,                               &
            & traj0_tl%opdp_path,                            &
            & traj0_sta%odsun_level,                         &
            & traj0_sta%odsun_singlelayer,                   &
            & traj0_sta%od_frac,                             &
            & traj0_dyn%transmission_aux,                    &
            & traj0_tl_dyn%transmission_aux,                 &
            & traj0_dyn%transmission_scatt_ir_stream,        &
            & traj0_tl_dyn%transmission_scatt_ir_stream,     &
            & traj0_sta%tausun_ref,                          &
            & traj0_sta%tausun_ref_surf,                     &
            & traj0_sta%tausun_level,                        &
            & traj0_sta%tausun_surf)
    ENDIF
  ENDIF

! TL of emissivity
!-------------------------------------------------
!11. calculate channel emissivity values - SURFACE
!-------------------------------------------------
! Where( .not. calcemis )
!   emissivity_out_tl = emissivity_tl
! ElseWhere
!   emissivity_out_tl = 0._jprb
! End Where
  emissivity_out_tl = emissivity_tl

  WHERE (.NOT. calcemis)
    traj0_tl%reflectivity =  - emissivity_out_tl
  ENDWHERE


  IF (Any(calcemis)) THEN
! calculate surface emissivity for selected channels
! and traj0%reflectivity

    IF (coefs%coef%id_sensor == sensor_id_ir) THEN
!Infrared
! nothing to do
      traj0_tl%reflectivity(:) =  - emissivity_out_tl(:)
    ELSE IF (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po) THEN
!Microwave
      CALL rttov_calcemis_mw_tl( &
            & opts,                             &
            & profiles,                         &
            & profiles_tl,                      &
            & traj0_sta%angles,                 &
            & coefs%coef,                       &
            & chanprof,                         &
            & traj0_dyn%transmission_aux,       &
            & traj0_tl_dyn%transmission_aux,    &
            & calcemis,                         &
            & emissivity_out_tl,                &
            & traj0_tl%reflectivity)
    ELSE
! Hires
      CALL rttov_calcemis_ir_tl( &
            & profiles,          &
            & profiles_tl,       &
            & coefs%coef,        &
            & opts%addpc,        &
            & coefs%coef_pccomp, &
            & chanprof,          &
            & calcemis,          &
            & emissivity_out_tl)
      traj0_tl%reflectivity(:) =  - emissivity_out_tl(:)
    ENDIF

  ENDIF

! TL of reflectance and sun glint
!-------------------------------------------------------
!12. Compute Fresnel reflectance and sun glint - SURFACE
!-------------------------------------------------------

  IF (coefs%coef%id_sensor == sensor_id_hi .AND. coefs%coef%fmv_model_ver == 9 .AND. opts%addsolar) THEN
    CALL rttov_refsun_tl( &
          & profiles,            &
          & profiles_tl,         &
          & coefs%coef,          &
          & traj0%aux_prof,      &
          & traj0%sunglint,      &
          & traj0_tl%sunglint,   &
          & traj0%raytracing,    &
          & traj0_tl%raytracing)
    CALL rttov_fresnel_tl( &
          & chanprof,           &
          & profiles,           &
          & coefs%coef,         &
          & traj0%sunglint,     &
          & traj0_tl%sunglint,  &
          & traj0%fresnrefl,    &
          & traj0_tl%fresnrefl)
  ENDIF

! TL of RTE
!---------------------------------------------------------
!13. integrate the radiative transfer equation - USER levs
!---------------------------------------------------------

  addcosmic = (coefs%coef%id_sensor == sensor_id_mw .OR. coefs%coef%id_sensor == sensor_id_po)
  CALL rttov_integrate_tl( &
        & addcosmic,                                    &
        & opts,                                         &
        & traj0_dyn%nstreams,                           &
        & chanprof,                                     &
        & emissivity_out,                               &
        & emissivity_out_tl,                            &
        & traj0%reflectivity,                           &
        & traj0_tl%reflectivity,                        &
        & traj0%fresnrefl,                              &
        & traj0_tl%fresnrefl,                           &
        & traj0%sunglint,                               &
        & traj0_tl%sunglint,                            &
        & traj0_sta%sun,                                &
        & traj0_dyn%transmission_aux,                   &
        & traj0_tl_dyn%transmission_aux,                &
        & traj0_dyn%transmission_scatt_ir_stream,       &
        & traj0_tl_dyn%transmission_scatt_ir_stream,    &
        & profiles,                                     &
        & profiles_tl,                                  &
        & traj0%aux_prof,                               &
        & traj0_tl%aux_prof,                            &
        & coefs%coef,                                   &
        & traj0%raytracing,                             &
        & traj0_tl%raytracing,                          &
        & traj0%ircld,                                  &
        & traj0_tl%ircld,                               &
        & radiancedata,                                 &
        & traj0_sta%auxrad,                             &
        & traj0_dyn%auxrad_stream,                      &
        & traj0_tl_dyn%auxrad_stream,                   &
        & radiancedata_tl)

  IF (opts%addpc) THEN
    CALL rttov_pcscores_tl(          &
          & opts,                    &
          & chanprof,                &
          & traj0_sta%chanprof_pc,   &
          & pccomp,                  &
          & pccomp_tl,               &
          & coefs%coef_pccomp,       &
          & radiancedata_tl)

    IF (opts%addradrec) THEN
      CALL rttov_reconstruct_tl(     &
            & traj0_sta%chanprof_in, &
            & traj0_sta%chanprof_pc, &
            & pccomp,                &
            & pccomp_tl,             &
            & coefs%coef_pccomp)
      CALL rttov_calcbt_pc_tl(       &
            & traj0_sta%chanprof_in, &
            & coefs%coef_pccomp,     &
            & pccomp,                &
            & pccomp_tl)
    ENDIF

  ENDIF

!--------------------
!14. deallocate memory
!--------------------

  CALL rttov_alloc_traj_dyn (err, traj0_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
  THROW(err.ne.0)
!
  CALL rttov_alloc_traj_dyn (err, traj0_tl_dyn, opts, nchannels, profiles(1)%nlayers, traj0_dyn%nstreams, ncldtyp, 0_jpim)
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
        & traj0_tl = traj0_tl, &
        & traj1 = traj1,       &
        & traj1_tl = traj1_tl, &
        & traj2 = traj,        &
        & traj2_tl = traj_tl)

  THROWM( ERR .NE. 0 , "deallocation check_traj")

  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 1_jpim, ZHOOK_HANDLE)
  CATCH
  errorstatus = ERR
  IF (LHOOK) CALL DR_HOOK('RTTOV_TL', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_tl
