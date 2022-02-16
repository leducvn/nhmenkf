! $Id: ropp_pp_diag2roprof.f90 2369 2009-11-23 11:12:05Z frhl $

!****s* PP/ropp_pp_diag2roprof *
!
! NAME
!    ropp_pp_diag2roprof - Add diagnostic information gathered during 
!                          processing an occultation to a ROprof data structure
!
! SYNOPSIS
!    call ropp_pp_diag2roprof(diag, ro_data)
!
! INPUTS
!    type(PPdiag) :: diag       Diagnostic information structure
!    type(ROprof) :: ro_data    RO data file structure
!
! OUTPUT
!    type(ROprof) :: ro_data    Updated RO data structure with diagnostics
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_pp_diag2roprof(diag, ro_data)


  ! 1. Declarations
  ! ---------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
  USE ropp_pp_types, ONLY: ppDiag

  IMPLICIT NONE

  TYPE(ppDiag), INTENT(in)     :: diag
  TYPE(ROprof), INTENT(inout)  :: ro_data


  ! 2. Badness scores
  ! -----------------

  CALL ropp_io_addvar(ro_data,                                   &
                      name      = "L2_badness_score",            &
                      long_name = "L2 correction badness score", &
                      units     = " ",                           &
                      range     = (/0.0_wp,1000.0_wp/),          &
                      data      = diag%L2_badness)

  CALL ropp_io_addvar(ro_data,                                              &
                      name      = "SO_badness_score",                       &
                      long_name = "Statistical optimisation badness score", &
                      units     = " ",                                      &
                      range     = (/0.0_wp,1000.0_wp/),                     &
                      data      = diag%sq)

  CALL ropp_io_addvar(ro_data,                             &
                      name      = "Badness_score",         &
                      long_name = "L2 + SO badness score", &
                      units     = " ",                     &
                      range     = (/0.0_wp,1000.0_wp/),    &
                      data      = diag%L2_badness+diag%sq)


  ! 3. CT amplitude processing
  ! --------------------------
  
  IF (ASSOCIATED(diag%CTimpact))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "CT_impact",          &
                      long_name = "CT impact parameter grid", &
                      units     = "m ",                   &
                      range     = (/0.0_wp,50000.0_wp/),                &
                      data      = diag%CTimpact)

  IF (ASSOCIATED(diag%CTamplitude))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "CT_amplitude",          &
                      long_name = "CT amplitude as function of CT_impact", &
                      units     = " ",                   &
                      range     = (/0.0_wp,1000.0_wp/),                &
                      data      = diag%CTamplitude)

  IF (ASSOCIATED(diag%CTamplitude_smt))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "CT_amplitude_smt",          &
                      long_name = "Smoothed CT amplitude as function of CT_impact", &
                      units     = " ",                   &
                      range     = (/0.0_wp,1000.0_wp/),                &
                      data      = diag%CTamplitude_smt)
  
  IF (ASSOCIATED(diag%CTimpactL2))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "CT_impactL2",          &
                      long_name = "CT impact parameter grid L2", &
                      units     = "m ",                   &
                      range     = (/0.0_wp,50000.0_wp/),                &
                      data      = diag%CTimpactL2)

  IF (ASSOCIATED(diag%CTamplitudeL2))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "CT_amplitudeL2",          &
                      long_name = "CT amplitude as function of CT_impactL2", &
                      units     = " ",                   &
                      range     = (/0.0_wp,1000.0_wp/),                &
                      data      = diag%CTamplitudeL2)

  IF (ASSOCIATED(diag%CTamplitudeL2_smt))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "CT_amplitudeL2_smt",          &
                      long_name = "Smoothed CT amplitude as function of CT_impactL2", &
                      units     = " ",                   &
                      range     = (/0.0_wp,1000.0_wp/),                &
                      data      = diag%CTamplitudeL2_smt)


  ! 4. Ionospheric correction bending angle
  ! ---------------------------------------

  IF (ASSOCIATED(diag%ba_ion))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "IC_ion_bangle",                &
                      long_name = "Ionospheric bending angle", &
                      units     = "radians",                   &
                      range     =     &
          (/ro_data%lev1b%range%bangle(1), ro_data%lev1b%range%bangle(2)/), &
                      data      = diag%ba_ion)

  IF (ASSOCIATED(diag%err_ion))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "IC_cov_bangle_ion",                &
                      long_name = "Ionospheric bending angle error covariance", &
                      units     = "radians",                   &
                      range     =     &
          (/ro_data%lev1b%range%bangle(1), ro_data%lev1b%range%bangle(2)/), &
                      data      = diag%err_ion)

  IF (ASSOCIATED(diag%err_neut))  &
     CALL ropp_io_addvar(ro_data,                                 &
                      name      = "IC_cov_bangle_neut",                &
                      long_name = "LC bending angle error covariance", &
                      units     = "radians",                   &
                      range     =     &
          (/ro_data%lev1b%range%bangle(1), ro_data%lev1b%range%bangle(2)/), &
                      data      = diag%err_neut)

  IF (ASSOCIATED(diag%wt_data))  &
     CALL ropp_io_addvar(ro_data,                                  &
                      name      = "LC_weight",                              &
                      long_name = "SO-weight of observed LC bending angle", &
                      units     = " ",                                      &
                      range     = (/0.0_wp,1.0_wp/),                        &
                      data      = diag%wt_data)

END SUBROUTINE ropp_pp_diag2roprof
