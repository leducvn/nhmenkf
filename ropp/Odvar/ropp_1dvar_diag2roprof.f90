! $Id: ropp_1dvar_diag2roprof.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* 1DVar/ropp_1dvar_diag2roprof *
!
! NAME
!    ropp_1dvar_diag2roprof - Add diagnostic information gathered during a
!                             1DVar retrieval to an ROProf data structure.
!
! SYNOPSIS
!    call ropp_1dvar_diag2roprof(obs, diag, ro_data, config)
!
! INPUTS
!    diag        - Diagnostic information structure
!    ro_data     - RO data file structure
!    config      - configuration options
!
! OUTPUT
!    ro_data     - Updated RO data file structure containing diagnostics
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

!-------------------------------------------------------------------------------
! 1. Bending angles
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_diag2roprof_bangle(obs, diag, ro_data, config)

! 1.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof, PCD_met
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_1dvar_copy, not_this => ropp_1dvar_diag2roprof_bangle

  IMPLICIT NONE

  TYPE(Obs1DBangle), INTENT(in)    :: obs
  TYPE(VarDiag),     INTENT(in)    :: diag
  TYPE(ROprof),      INTENT(inout) :: ro_data
  TYPE(VarConfig),   INTENT(in)    :: config

  LOGICAL :: dummy
  dummy = obs%obs_ok !! fix nag 'unused dummy variable'

! 1.2 1DVar diagnostics
! ---------------------

! 1.2.1 Overall quality of retrieval

  IF (diag % ok) THEN
     ro_data % overall_qual = 100.0_wp
  ELSE
     ro_data % overall_qual = 0.0_wp
     ro_data % PCD = IBSET(ro_data%PCD, PCD_met)
  ENDIF

! 1.2.2 Cost function

  CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J",                                         &
                      long_name = "Cost function value at convergence",        &
                      units     = "",                                          &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J)

! 1.2.3 Scaled cost function

  CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J_scaled",                                  &
                      long_name = "Scaled cost function value at convergence ('2J/m')", &
                      units     = "",                                          &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J_scaled)

! 1.2.4 Original cost function

  CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J_init",                                    &
                      long_name = "Initial cost function value",        &
                      units     = "",                                          &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J_init)

! 1.3 Extended 1DVar diagnostics
! ------------------------------

  IF (config % extended_1dvar_diag) THEN


    ! 1.3.1 Original cost function

    IF (ASSOCIATED(diag%J_bgr)) &
      CALL ropp_io_addvar(ro_data,                                           &
                      name      = "J_bgr",                                   &
                      long_name = "Background contribution to cost function profile",        &
                      units     = "",                                        &
                      range     = (/ 0.0_wp, 9999.0_wp/),                    &
                      DATA      = diag % J_bgr)

    IF (ASSOCIATED(diag%J_obs)) &
       CALL ropp_io_addvar(ro_data,                                          &
                      name      = "J_obs",                                   &
                      long_name = "Observation contribution to cost function profile",        &
                      units     = "",                                        &
                      range     = (/ 0.0_wp, 9999.0_wp/),                    &
                      DATA      = diag % J_obs)

     ! 1.3.1 Number of iterations

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "n_iter",                                 &
                         long_name = "Number of iterations",                   &
                         units     = "",                                       &
                         range     = (/ 0.0_wp, 9999.0_wp/),                   &
                         DATA      = REAL(diag % n_iter, wp))

     ! 1.3.2 Exit mode of minimiser

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "minropp_exit_mode",                      &
                         long_name = "Exit mode of minimiser (minROPP)",       &
                         units     = "",                                       &
                         range     = (/ 0.0_wp, 10.0_wp/),                     &
                         DATA      = REAL(diag % min_mode, wp))

     ! 1.3.3 O - B

     IF (ASSOCIATED(diag%OmB)) &
        CALL ropp_io_addvar(ro_data,                                           &
                         name      = "OmB",                                    &
                         long_name = "Observation - Background bending angle", &
                         units     = "radians",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmB)

     IF (ASSOCIATED(diag%OmB_sigma)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "OmB_sigma",                                    &
                         long_name = "Expected (Observation - Background bending angle) standard deviation", &
                         units     = "radians",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmB_sigma)

     ! 1.3.4 O - A

     IF (ASSOCIATED(diag%OmA)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "OmA",                                    &
                         long_name = "Observation - Analysis bending angle",   &
                         units     = "radians",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmA)

     IF (ASSOCIATED(diag%OmA_sigma)) &
        CALL ropp_io_addvar(ro_data,                                          &
                         name      = "OmA_sigma",                             &
                         long_name = "Expected (Observation - Analysis bending angle) standard deviation",   &
                         units     = "radians",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmA_sigma)

     IF (ASSOCIATED(diag%B_sigma)) &
        CALL ropp_io_addvar(ro_data,                                          &
                         name      = "B_sigma",                               &
                         long_name = "Background bending angle standard deviation", &
                         units     = "radians",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % B_sigma)

     ! 1.3.5 Probability of Gross Error

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "pge",                                    &
                         long_name = "Probability of Gross Error",             &
                         units     = "",                                       &
                         range     = (/ 0.0_wp, 1.0_wp/),                      &
                         DATA      = diag % pge)

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "pge_gamma",                                    &
                         long_name = "Probability of Gross Error 'gamma' value",             &
                         units     = "",                                       &
                         range     = (/ 0.0_wp, 1.0_wp/),                      &
                         DATA      = diag % pge_gamma)

  ENDIF

END SUBROUTINE ropp_1dvar_diag2roprof_bangle


!-------------------------------------------------------------------------------
! 2. Refractivity
!-------------------------------------------------------------------------------

SUBROUTINE ropp_1dvar_diag2roprof_refrac(obs, diag, ro_data, config)

! 2.1 Declarations
! ----------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof, PCD_met
  USE ropp_fm
  USE ropp_1dvar
! USE ropp_1dvar_copy, not_this => ropp_1dvar_diag2roprof_refrac

  IMPLICIT NONE

  TYPE(Obs1DRefrac), INTENT(in)    :: obs
  TYPE(VarDiag),     INTENT(in)    :: diag
  TYPE(ROprof),      INTENT(inout) :: ro_data
  TYPE(VarConfig),   INTENT(in)    :: config

  LOGICAL :: dummy
  dummy = obs%obs_ok !! fix nag 'unused dummy variable'

! 2.2 1DVar diagnostics
! ---------------------

! 2.2.1 Overall quality of retrieval

  IF (diag % ok) THEN
     ro_data % overall_qual = 100.0_wp
  ELSE
     ro_data % overall_qual = 0.0_wp
     ro_data % PCD = IBSET(ro_data%PCD, PCD_met)
  ENDIF

! 2.2.2 Cost function

  CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J",                                         &
                      long_name = "Cost function value at convergence",        &
                      units     = "1",                                         &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J)

! 2.2.3 Scaled cost function

  CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J_scaled",                                  &
                      long_name = "Scaled cost function value at convergence ('2J/m')", &
                      units     = "1",                                         &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J_scaled)

! 1.2.4 Original cost function

  CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J_init",                                    &
                      long_name = "Initial Cost function value",        &
                      units     = "",                                          &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J_init)

! 2.3 Extended 1DVar diagnostics
! ------------------------------

  IF (config % extended_1dvar_diag) THEN

    ! 2.3.1 Original cost function

     IF (ASSOCIATED(diag%J_bgr)) &
        CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J_bgr",                                    &
                      long_name = "Background contribution to cost function profile",        &
                      units     = "",                                          &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J_bgr)

    IF (ASSOCIATED(diag%J_obs)) &
       CALL ropp_io_addvar(ro_data,                                                 &
                      name      = "J_obs",                                    &
                      long_name = "Observation contribution to cost function profile",        &
                      units     = "",                                          &
                      range     = (/ 0.0_wp, 9999.0_wp/),                      &
                      DATA      = diag % J_obs)

     ! 2.3.1 Number of iterations

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "n_iter",                                 &
                         long_name = "Number of iterations",                   &
                         units     = "1",                                      &
                         range     = (/ 0.0_wp, 9999.0_wp/),                   &
                         DATA      = REAL(diag % n_iter, wp))

     ! 2.3.2 Exit mode of minimiser

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "minropp_exit_mode",                      &
                         long_name = "Exit mode of minimiser (minROPP)",       &
                         units     = "1",                                      &
                         range     = (/ 0.0_wp, 10.0_wp/),                     &
                         DATA      = REAL(diag % min_mode, wp))

     ! 2.3.3 O - B

     IF (ASSOCIATED(diag%OmB)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "OmB",                                    &
                         long_name = "Observation - Background refractivity",  &
                         units     = "1",                                      &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmB)

     IF (ASSOCIATED(diag%OmB_sigma)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "OmB_sigma",                                    &
                         long_name = "Expected (Observation - Background refractivity) standard deviation", &
                         units     = "1",                                      &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmB_sigma)

     ! 2.3.4 O - A

     IF (ASSOCIATED(diag%OmA)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "OmA",                                    &
                         long_name = "Observation - Analysis refractivity",    &
                         units     = "1",                                      &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmA)

     IF (ASSOCIATED(diag%OmA_sigma)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "OmA_sigma",                                    &
                         long_name = "Expected (Observation - Analysis refractivity) standard deviation", &
                         units     = "1",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % OmA_sigma)

     IF (ASSOCIATED(diag%B_sigma)) &
        CALL ropp_io_addvar(ro_data,                                              &
                         name      = "B_sigma",                                    &
                         long_name = "Background refractivity standard deviation", &
                         units     = "1",                                &
                         range     = (/ -999.9_wp, 999.9_wp /),                &
                         DATA      = diag % B_sigma)

     ! 2.3.5 Probability of Gross Error

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "pge",                                    &
                         long_name = "Probability of Gross Error",             &
                         units     = "1",                                      &
                         range     = (/ 0.0_wp, 1.0_wp/),                      &
                         DATA      = diag % pge)

     CALL ropp_io_addvar(ro_data,                                              &
                         name      = "pge_gamma",                                    &
                         long_name = "Probability of Gross Error 'gamma' value",             &
                         units     = "",                                       &
                         range     = (/ 0.0_wp, 1.0_wp/),                      &
                         DATA      = diag % pge_gamma)

  ENDIF

END SUBROUTINE ropp_1dvar_diag2roprof_refrac


