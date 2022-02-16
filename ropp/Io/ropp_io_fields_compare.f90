! $Id: ropp_io_fields_compare.f90 2197 2009-06-23 09:11:17Z idculv $

!****si* Test/ropp_io_fields_compare *
!
! NAME
!    ropp_io_fields_compare - Compares two ROPP profiles for meaningful differences.
!
! SYNOPSIS
!    CALL ropp_io_fields_compare (prof1, prof2, ndiff, &
!         bg, GEOref, Lev1a, Lev1b, Lev2a, Lev2b, Lev2c, Lev2d, &
!         L1L2=L1L2, onedvar=onedvar, spectra=spectra)
!
! DESCRIPTION
!    Compares two ROPP profiles.
!    Used by ropp_<module>/tests/ropp_<module>_compare.f90 to compare output.
!
! INPUTS
!    ROprofs prof1 and prof2, and a set of logicals which define the 
!    'blocks' of data to be compared.
!
! OUTPUT
!    The (updated) number of 'significantly' different elements of the two profiles. 
!
! NOTES
!    The various thresholds for defining differences as significant
!    are subject to change, depending on experience with this routine.
!
! SEE ALSO
!    ropp_<module>/tests/ropp_<module>_compare.f90
!
! AUTHOR
!    Met Office, Exeter, UK.
!    Any comments on this software should be given via the ROM SAF
!    Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   (c) EUMETSAT. All rights reserved.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_io_fields_compare(prof1, prof2, ndiff, &
           bg, GEOref, Lev1a, Lev1b, Lev2a, Lev2b, Lev2c, Lev2d, &
           L1L2, onedvar, spectra, tdry)

  USE typesizes,     wp => EightByteReal
  USE ropp_io_types, ONLY: ROprof
  USE messages

! Input/output variables
  TYPE(ROprof), INTENT(in)        :: prof1            ! RO profile from file1.nc
  TYPE(ROprof), INTENT(in)        :: prof2            ! RO profile from file2.nc
  INTEGER, INTENT(inout)          :: ndiff            ! Number of differences
  LOGICAL, INTENT(in)             :: bg, GEOref, &    ! The substructures of
                                     Lev1a, Lev1b, &  ! ROprof to be compared
                                     Lev2a, Lev2b, &
                                     Lev2c, Lev2d
  LOGICAL, OPTIONAL, INTENT(in)   :: L1L2             ! False when comparing -direct_ion and
                                                      ! neutral L1 and L2 fields
  LOGICAL, OPTIONAL, INTENT(in)   :: onedvar          ! True to relax difference
                                                      ! tolerance for retrievals
  LOGICAL, OPTIONAL, INTENT(in)   :: spectra          ! True to define variables for spectra tests
  LOGICAL, OPTIONAL, INTENT(in)   :: tdry             ! True to relax tdry tolerances for occ and rs tests

! Local variables
! REAL(wp)                        :: tol              ! Difference tolerance ! Commented at 20 July, 2016
  LOGICAL                         :: L1L2_local       ! Local version of L1L2
  LOGICAL                         :: onedvar_local    ! Local version of onedvar
  LOGICAL                         :: spectra_local    ! Local version of spectra
  LOGICAL                         :: tdry_local       ! Local version of tdry

  CHARACTER(LEN=256)              :: routine          ! Routine name for diagnostics


! ------------------------------------------------------------------------------
! 0. Initialise
! ------------------------------------------------------------------------------

! 0.1 Set routine for error messages
! ----------------------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_fields_compare')

! 0.2 Define local logical variables
! ----------------------------------

  L1L2_local = .TRUE.
  IF ( PRESENT(L1L2) ) L1L2_local = L1L2

  onedvar_local = .FALSE.
  IF ( PRESENT(onedvar) ) onedvar_local = onedvar

  spectra_local = .FALSE.
  IF ( PRESENT(spectra) ) spectra_local = spectra

  tdry_local = .FALSE.
  IF ( PRESENT(tdry) ) tdry_local = tdry


! ------------------------------------------------------------------------------
! 1. Compare the elements of the profiles
! ------------------------------------------------------------------------------

! 1.1 Background
! --------------

  IF (bg) THEN

    CALL difference_char_0d(prof1%bg%source, prof2%bg%source, ndiff, 'bg%source')

    CALL difference_int_0d(prof1%bg%year, prof2%bg%year, 0, ndiff, 'bg%year')
    CALL difference_int_0d(prof1%bg%month, prof2%bg%month, 0, ndiff, 'bg%month')
    CALL difference_int_0d(prof1%bg%day, prof2%bg%day, 0, ndiff, 'bg%day')
    CALL difference_int_0d(prof1%bg%hour, prof2%bg%hour, 0, ndiff, 'bg%hour')
    CALL difference_int_0d(prof1%bg%minute, prof2%bg%minute, 0, ndiff, 'bg%minute')

    CALL difference_real_0d(prof1%bg%Fcperiod, prof2%bg%Fcperiod, 1.0e-6_wp, ndiff, 'bg%Fcperiod')

  ENDIF ! bg

! 1.2 Georeferencing
! ------------------

  IF (GEOref) THEN

    CALL difference_real_0d(prof1%GEOref%lat, prof2%GEOref%lat, 1.0e-6_wp, ndiff, 'GEOref%lat')
    CALL difference_real_0d(prof1%GEOref%lon, prof2%GEOref%lon, 1.0e-6_wp, ndiff, 'GEOref%lon')
    CALL difference_real_0d(prof1%GEOref%roc, prof2%GEOref%roc, 1.0_wp, ndiff, 'GEOref%roc')
    CALL difference_real_0d(prof1%GEOref%azimuth, prof2%GEOref%azimuth, 1.0e-6_wp, ndiff, 'GEOref%azimuth')
    CALL difference_real_0d(prof1%GEOref%undulation, prof2%GEOref%undulation, 1.0e-6_wp, ndiff, 'GEOref%undulation')
    IF (ABS(prof1%GEOref%undulation - prof2%GEOref%undulation) > 1.0e-6_wp) &
      CALL message ( msg_info, 'Suggest checking that GEOPOT_COEF and GEOPOT_CORR ' // &
                               'have been set correctly. (See grib2bgrasc man page.) \n')

    CALL difference_real_1d(prof1%GEOref%r_coc, prof2%GEOref%r_coc, 1.0_wp, ndiff, 'GEOref%r_coc')

  ENDIF ! GEOref

! 1.3 Level 1a
! ------------

  IF (Lev1a) THEN

    CALL difference_int_0d(prof1%Lev1a%Npoints, prof2%Lev1a%Npoints, 0, ndiff, 'Lev1a%Npoints')

    CALL difference_real_1d(prof1%Lev1a%dtime, prof2%Lev1a%dtime, 1.0e-6_wp, ndiff, 'Lev1a%dtime')
    CALL difference_real_1d(prof1%Lev1a%snr_L1ca, prof2%Lev1a%snr_L1ca, 1.0e-6_wp, ndiff, 'Lev1a%snr_L1ca')
    CALL difference_real_1d(prof1%Lev1a%snr_L1p, prof2%Lev1a%snr_L1p, 1.0e-6_wp, ndiff, 'Lev1a%snr_L1p')
    CALL difference_real_1d(prof1%Lev1a%snr_L2p, prof2%Lev1a%snr_L2p, 1.0e-6_wp, ndiff, 'Lev1a%snr_L2p')
! Relax tolerance here to +/- 5*2pi, because different compilers can introduce
! such phase jumps when calculating (eg) MODULO functions in accumulate_phase.
!    CALL difference_real_1d(prof1%Lev1a%phase_L1, prof2%Lev1a%phase_L1, 1.0e-6_wp, ndiff, 'Lev1a%phase_L1')
    CALL difference_real_1d(prof1%Lev1a%phase_L1, prof2%Lev1a%phase_L1, 1.0_wp, ndiff, 'Lev1a%phase_L1')
    CALL difference_real_1d(prof1%Lev1a%phase_L2, prof2%Lev1a%phase_L2, 1.0e-6_wp, ndiff, 'Lev1a%phase_L2')
    CALL difference_real_1d(prof1%Lev1a%phase_qual, prof2%Lev1a%phase_qual, 1.0e-6_wp, ndiff, 'Lev1a%phase_qual')

!  Open_loop_lcf is not (yet) part of standard Lev1a substructure.
!  Instead, if present, it is (probably) held in ROprof%vlist%VlistD1d%DATA.
    IF (ASSOCIATED(prof1%vlist%VlistD1d) .AND. ASSOCIATED(prof1%vlist%VlistD1d)) THEN
      IF (TRIM(ADJUSTL(prof1%vlist%VlistD1d%name)) == 'open_loop_lcf' .AND. &
          TRIM(ADJUSTL(prof2%vlist%VlistD1d%name)) == 'open_loop_lcf') THEN
        CALL difference_real_1d(prof1%vlist%VlistD1d%DATA, prof2%vlist%VlistD1d%DATA, 1.0e-6_wp, ndiff, 'Lev1a%open_loop_lcf')
      ENDIF
    ENDIF

    CALL difference_real_2d(prof1%Lev1a%r_gns, prof2%Lev1a%r_gns, 1.0_wp, ndiff, 'Lev1a%r_gns')
    CALL difference_real_2d(prof1%Lev1a%v_gns, prof2%Lev1a%v_gns, 1.0e-3_wp, ndiff, 'Lev1a%v_gns')
    CALL difference_real_2d(prof1%Lev1a%r_leo, prof2%Lev1a%r_leo, 1.0_wp, ndiff, 'Lev1a%r_leo')
    CALL difference_real_2d(prof1%Lev1a%v_leo, prof2%Lev1a%v_leo, 1.0e-3_wp, ndiff, 'Lev1a%v_leo')

  ENDIF ! Lev1a

! 1.4 Level 1b
! ------------

  IF (Lev1b) THEN

    CALL difference_int_0d(prof1%Lev1b%Npoints, prof2%Lev1b%Npoints, 0, ndiff, 'Lev1b%Npoints')

    IF ( onedvar_local ) THEN  ! Relax difference thresholds

      CALL difference_real_1d(prof1%Lev1b%lat_tp, prof2%Lev1b%lat_tp, 1.0e-4_wp, ndiff, 'Lev1b%lat_tp')
      CALL difference_real_1d(prof1%Lev1b%lon_tp, prof2%Lev1b%lon_tp, 1.0e-4_wp, ndiff, 'Lev1b%lon_tp')
      CALL difference_real_1d(prof1%Lev1b%azimuth_tp, prof2%Lev1b%azimuth_tp, 1.0e-2_wp, ndiff, 'Lev1b%azimuth_tp')
      CALL difference_real_1d(prof1%Lev1b%impact, prof2%Lev1b%impact, 1.0_wp, ndiff, 'Lev1b%impact')
      CALL difference_real_1d(prof1%Lev1b%impact_opt, prof2%Lev1b%impact_opt, 1.0_wp, ndiff, 'Lev1b%impact_opt')
      CALL difference_real_1d(prof1%Lev1b%bangle, prof2%Lev1b%bangle, 1.0e-6_wp, ndiff, 'Lev1b%bangle')
      CALL difference_real_1d(prof1%Lev1b%bangle_sigma, prof2%Lev1b%bangle_sigma, 1.0e-6_wp, ndiff, 'Lev1b%bangle_sigma')
      CALL difference_real_1d(prof1%Lev1b%bangle_opt, prof2%Lev1b%bangle_opt, 1.0e-6_wp, ndiff, 'Lev1b%bangle_opt')
      CALL difference_real_1d(prof1%Lev1b%bangle_opt_sigma, prof2%Lev1b%bangle_opt_sigma,1.0e-6_wp,ndiff,'Lev1b%bangle_opt_sigma')
      CALL difference_real_1d(prof1%Lev1b%bangle_qual, prof2%Lev1b%bangle_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_qual')
      CALL difference_real_1d(prof1%Lev1b%bangle_opt_qual, prof2%Lev1b%bangle_opt_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_opt_qual')

      IF (L1L2_local) THEN  ! These tests need to be switched off when comparing -direct_ion results

        CALL difference_real_1d(prof1%Lev1b%impact_L1, prof2%Lev1b%impact_L1, 1.0_wp, ndiff, 'Lev1b%impact_L1')
        CALL difference_real_1d(prof1%Lev1b%bangle_L1, prof2%Lev1b%bangle_L1, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L1')
        CALL difference_real_1d(prof1%Lev1b%bangle_L1_sigma, prof2%Lev1b%bangle_L1_sigma, 1.0e-6_wp, ndiff,'Lev1b%bangle_L1_sigma')
        CALL difference_real_1d(prof1%Lev1b%bangle_L1_qual, prof2%Lev1b%bangle_L1_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L1_qual')
        CALL difference_real_1d(prof1%Lev1b%impact_L2, prof2%Lev1b%impact_L2, 1.0_wp, ndiff, 'Lev1b%impact_L2')
        CALL difference_real_1d(prof1%Lev1b%bangle_L2, prof2%Lev1b%bangle_L2, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L2')
        CALL difference_real_1d(prof1%Lev1b%bangle_L2_sigma, prof2%Lev1b%bangle_L2_sigma, 1.0e-6_wp, ndiff,'Lev1b%bangle_L2_sigma')
        CALL difference_real_1d(prof1%Lev1b%bangle_L2_qual, prof2%Lev1b%bangle_L2_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L2_qual')

      ENDIF ! L1 and L2 testing

    ELSE

      CALL difference_real_1d(prof1%Lev1b%lat_tp, prof2%Lev1b%lat_tp, 1.0e-4_wp, ndiff, 'Lev1b%lat_tp')
      CALL difference_real_1d(prof1%Lev1b%lon_tp, prof2%Lev1b%lon_tp, 1.0e-4_wp, ndiff, 'Lev1b%lon_tp')
      CALL difference_real_1d(prof1%Lev1b%azimuth_tp, prof2%Lev1b%azimuth_tp, 1.0e-2_wp, ndiff, 'Lev1b%azimuth_tp')
      CALL difference_real_1d(prof1%Lev1b%impact, prof2%Lev1b%impact, 1.0_wp, ndiff, 'Lev1b%impact')
      CALL difference_real_1d(prof1%Lev1b%impact_opt, prof2%Lev1b%impact_opt, 1.0_wp, ndiff, 'Lev1b%impact_opt')
      CALL difference_real_1d(prof1%Lev1b%bangle, prof2%Lev1b%bangle, 1.0e-6_wp, ndiff, 'Lev1b%bangle')
      CALL difference_real_1d(prof1%Lev1b%bangle_sigma, prof2%Lev1b%bangle_sigma, 1.0e-6_wp, ndiff, 'Lev1b%bangle_sigma')
      CALL difference_real_1d(prof1%Lev1b%bangle_opt, prof2%Lev1b%bangle_opt, 1.0e-6_wp, ndiff, 'Lev1b%bangle_opt')
      CALL difference_real_1d(prof1%Lev1b%bangle_opt_sigma, prof2%Lev1b%bangle_opt_sigma,1.0e-6_wp,ndiff,'Lev1b%bangle_opt_sigma')
      CALL difference_real_1d(prof1%Lev1b%bangle_qual, prof2%Lev1b%bangle_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_qual')
      CALL difference_real_1d(prof1%Lev1b%bangle_opt_qual, prof2%Lev1b%bangle_opt_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_opt_qual')

      IF (L1L2_local) THEN  ! These tests need to be switched off when comparing -direct_ion results

        CALL difference_real_1d(prof1%Lev1b%impact_L1, prof2%Lev1b%impact_L1, 1.0_wp, ndiff, 'Lev1b%impact_L1')
        CALL difference_real_1d(prof1%Lev1b%bangle_L1, prof2%Lev1b%bangle_L1, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L1')
        CALL difference_real_1d(prof1%Lev1b%bangle_L1_sigma, prof2%Lev1b%bangle_L1_sigma, 1.0e-6_wp, ndiff,'Lev1b%bangle_L1_sigma')
        CALL difference_real_1d(prof1%Lev1b%bangle_L1_qual, prof2%Lev1b%bangle_L1_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L1_qual')
        CALL difference_real_1d(prof1%Lev1b%impact_L2, prof2%Lev1b%impact_L2, 1.0_wp, ndiff, 'Lev1b%impact_L2')
        CALL difference_real_1d(prof1%Lev1b%bangle_L2, prof2%Lev1b%bangle_L2, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L2')
        CALL difference_real_1d(prof1%Lev1b%bangle_L2_sigma, prof2%Lev1b%bangle_L2_sigma, 1.0e-6_wp, ndiff,'Lev1b%bangle_L2_sigma')
        CALL difference_real_1d(prof1%Lev1b%bangle_L2_qual, prof2%Lev1b%bangle_L2_qual, 1.0e-6_wp, ndiff, 'Lev1b%bangle_L2_qual')

      ENDIF ! L1 and L2 testing

    ENDIF

  ENDIF ! Lev1b

! 1.5 Level 2a
! ------------

  IF (Lev2a) THEN

    CALL difference_int_0d(prof1%Lev2a%Npoints, prof2%Lev2a%Npoints, 0, ndiff, 'Lev2a%Npoints')

    IF ( onedvar_local ) THEN  ! Relax difference thresholds
      CALL difference_real_1d(prof1%Lev2a%refrac, prof2%Lev2a%refrac, 0.5_wp, ndiff, 'Lev2a%refrac') ! Very high, but necessary
      CALL difference_real_1d(prof1%Lev2a%refrac_sigma, prof2%Lev2a%refrac_sigma, 1.0e-3_wp, ndiff, 'Lev2a%refrac_sigma')
    ELSE
      CALL difference_real_1d(prof1%Lev2a%refrac, prof2%Lev2a%refrac, 1.0e-6_wp, ndiff, 'Lev2a%refrac')
      CALL difference_real_1d(prof1%Lev2a%refrac_sigma, prof2%Lev2a%refrac_sigma, 1.0e-6_wp, ndiff, 'Lev2a%refrac_sigma')
    ENDIF

    IF ( tdry_local ) THEN ! Relax tdry tols for ropp_pp occ and rs tests
      CALL difference_real_1d(prof1%Lev2a%dry_temp, prof2%Lev2a%dry_temp, 1.0e-4_wp, ndiff, 'Lev2a%dry_temp')
    ELSE
      CALL difference_real_1d(prof1%Lev2a%dry_temp, prof2%Lev2a%dry_temp, 1.0e-6_wp, ndiff, 'Lev2a%dry_temp')
    ENDIF

    CALL difference_real_1d(prof1%Lev2a%alt_refrac, prof2%Lev2a%alt_refrac, 1.0_wp, ndiff, 'Lev2a%alt_refrac')
    CALL difference_real_1d(prof1%Lev2a%geop_refrac, prof2%Lev2a%geop_refrac, 1.0_wp, ndiff, 'Lev2a%geop_refrac')
    CALL difference_real_1d(prof1%Lev2a%refrac_qual, prof2%Lev2a%refrac_qual, 1.0e-6_wp, ndiff, 'Lev2a%refrac_qual')
    CALL difference_real_1d(prof1%Lev2a%dry_temp_sigma, prof2%Lev2a%dry_temp_sigma, 1.0e-6_wp, ndiff, 'Lev2a%dry_temp_sigma')
    CALL difference_real_1d(prof1%Lev2a%dry_temp_qual, prof2%Lev2a%dry_temp_qual, 1.0e-6_wp, ndiff, 'Lev2a%dry_temp_qual')

  ENDIF ! Lev2a

! 1.6 Level 2b
! ------------

  IF (Lev2b) THEN

    CALL difference_int_0d(prof1%Lev2b%Npoints, prof2%Lev2b%Npoints, 0, ndiff, 'Lev2b%Npoints')

    IF ( onedvar_local ) THEN  ! Relax difference thresholds

      CALL difference_real_1d(prof1%Lev2b%press, prof2%Lev2b%press, 1.0e-1_wp, ndiff, 'Lev2b%press')
      CALL difference_real_1d(prof1%Lev2b%temp, prof2%Lev2b%temp, 1.0e-1_wp, ndiff, 'Lev2b%temp')
      CALL difference_real_1d(prof1%Lev2b%shum, prof2%Lev2b%shum, 1.0e-1_wp, ndiff, 'Lev2b%shum')
      CALL difference_real_1d(prof1%Lev2b%geop, prof2%Lev2b%geop, 1.0e0_wp, ndiff, 'Lev2b%geop')
      CALL difference_real_1d(prof1%Lev2b%meteo_qual, prof2%Lev2b%meteo_qual, 1.0e-3_wp, ndiff, 'Lev2b%meteo_qual')

    ELSE

      CALL difference_real_1d(prof1%Lev2b%press, prof2%Lev2b%press, 1.0e-6_wp, ndiff, 'Lev2b%press')
      CALL difference_real_1d(prof1%Lev2b%temp, prof2%Lev2b%temp, 1.0e-6_wp, ndiff, 'Lev2b%temp')
      CALL difference_real_1d(prof1%Lev2b%shum, prof2%Lev2b%shum, 1.0e-6_wp, ndiff, 'Lev2b%shum')
      CALL difference_real_1d(prof1%Lev2b%geop, prof2%Lev2b%geop, 1.0e-3_wp, ndiff, 'Lev2b%geop')
      CALL difference_real_1d(prof1%Lev2b%meteo_qual, prof2%Lev2b%meteo_qual, 1.0e-6_wp, ndiff, 'Lev2b%meteo_qual')

    ENDIF

  ENDIF ! Lev2b

! 1.7 Level 2c
! ------------

  IF (Lev2c) THEN

    CALL  difference_int_0d(prof1%Lev2c%Npoints, prof2%Lev2c%Npoints, 0, ndiff, 'Lev2c%Npoints')

    CALL difference_real_0d(prof1%Lev2c%geop_sfc, prof2%Lev2c%geop_sfc, 1.0e-6_wp, ndiff, 'Lev2c%geop_sfc')
    CALL difference_real_0d(prof1%Lev2c%press_sfc, prof2%Lev2c%press_sfc, 1.0e-6_wp, ndiff, 'Lev2c%press_sfc')
    CALL difference_real_0d(prof1%Lev2c%press_sfc_qual, prof2%Lev2c%press_sfc_qual, 1.0e-6_wp, ndiff, 'Lev2c%press_sfc_qual')
    CALL difference_real_0d(prof1%Lev2c%ne_max, prof2%Lev2c%ne_max, 1.0e-6_wp, ndiff, 'Lev2c%ne_max')
    CALL difference_real_0d(prof1%Lev2c%h_peak, prof2%Lev2c%h_peak, 1.0e-6_wp, ndiff, 'Lev2c%h_peak')
    CALL difference_real_0d(prof1%Lev2c%h_width, prof2%Lev2c%h_width, 1.0e-6_wp, ndiff, 'Lev2c%h_width')

    CALL difference_real_0d(prof1%Lev2c%tph_bangle, prof2%Lev2c%tph_bangle, 1.0e0_wp, ndiff, 'Lev2c%tph_bangle')
    CALL difference_real_0d(prof1%Lev2c%tpa_bangle, prof2%Lev2c%tpa_bangle, 1.0e-6_wp, ndiff, 'Lev2c%tpa_bangle')
    CALL  difference_int_0d(prof1%Lev2c%tph_bangle_flag, prof2%Lev2c%tph_bangle_flag, 1, ndiff, 'Lev2c%tph_bangle_flag')

    CALL difference_real_0d(prof1%Lev2c%tph_refrac, prof2%Lev2c%tph_refrac, 1.0e0_wp, ndiff, 'Lev2c%tph_refrac')
    CALL difference_real_0d(prof1%Lev2c%tpn_refrac, prof2%Lev2c%tpn_refrac, 1.0e-6_wp, ndiff, 'Lev2c%tpn_refrac')
    CALL  difference_int_0d(prof1%Lev2c%tph_refrac_flag, prof2%Lev2c%tph_refrac_flag, 1, ndiff, 'Lev2c%tph_refrac_flag')

    CALL difference_real_0d(prof1%Lev2c%tph_tdry_lrt, prof2%Lev2c%tph_tdry_lrt, 1.0e0_wp, ndiff, 'Lev2c%tph_tdry_lrt')
    CALL difference_real_0d(prof1%Lev2c%tpt_tdry_lrt, prof2%Lev2c%tpt_tdry_lrt, 1.0e-6_wp, ndiff, 'Lev2c%tpt_tdry_lrt')
    CALL  difference_int_0d(prof1%Lev2c%tph_tdry_lrt_flag, prof2%Lev2c%tph_tdry_lrt_flag, 1, ndiff, 'Lev2c%tph_tdry_lrt_flag')
    CALL difference_real_0d(prof1%Lev2c%tph_tdry_cpt, prof2%Lev2c%tph_tdry_cpt, 1.0e0_wp, ndiff, 'Lev2c%tph_tdry_cpt')
    CALL difference_real_0d(prof1%Lev2c%tpt_tdry_cpt, prof2%Lev2c%tpt_tdry_cpt, 1.0e-6_wp, ndiff, 'Lev2c%tpt_tdry_cpt')
    CALL  difference_int_0d(prof1%Lev2c%tph_tdry_cpt_flag, prof2%Lev2c%tph_tdry_cpt_flag, 1, ndiff, 'Lev2c%tph_tdry_cpt_flag')
    CALL difference_real_0d(prof1%Lev2c%prh_tdry_cpt, prof2%Lev2c%prh_tdry_cpt, 1.0e0_wp, ndiff, 'Lev2c%prh_tdry_cpt')
    CALL difference_real_0d(prof1%Lev2c%prt_tdry_cpt, prof2%Lev2c%prt_tdry_cpt, 1.0e-6_wp, ndiff, 'Lev2c%prt_tdry_cpt')
    CALL  difference_int_0d(prof1%Lev2c%prh_tdry_cpt_flag, prof2%Lev2c%prh_tdry_cpt_flag, 1, ndiff, 'Lev2c%prh_tdry_cpt_flag')

    CALL difference_real_0d(prof1%Lev2c%tph_temp_lrt, prof2%Lev2c%tph_temp_lrt, 1.0e0_wp, ndiff, 'Lev2c%tph_temp_lrt')
    CALL difference_real_0d(prof1%Lev2c%tpt_temp_lrt, prof2%Lev2c%tpt_temp_lrt, 1.0e-6_wp, ndiff, 'Lev2c%tpt_temp_lrt')
    CALL  difference_int_0d(prof1%Lev2c%tph_temp_lrt_flag, prof2%Lev2c%tph_temp_lrt_flag, 1, ndiff, 'Lev2c%tph_temp_lrt_flag')
    CALL difference_real_0d(prof1%Lev2c%tph_temp_cpt, prof2%Lev2c%tph_temp_cpt, 1.0e0_wp, ndiff, 'Lev2c%tph_temp_cpt')
    CALL difference_real_0d(prof1%Lev2c%tpt_temp_cpt, prof2%Lev2c%tpt_temp_cpt, 1.0e-6_wp, ndiff, 'Lev2c%tpt_temp_cpt')
    CALL  difference_int_0d(prof1%Lev2c%tph_temp_cpt_flag, prof2%Lev2c%tph_temp_cpt_flag, 1, ndiff, 'Lev2c%tph_temp_cpt_flag')
    CALL difference_real_0d(prof1%Lev2c%prh_temp_cpt, prof2%Lev2c%prh_temp_cpt, 1.0e0_wp, ndiff, 'Lev2c%prh_temp_cpt')
    CALL difference_real_0d(prof1%Lev2c%prt_temp_cpt, prof2%Lev2c%prt_temp_cpt, 1.0e-6_wp, ndiff, 'Lev2c%prt_temp_cpt')
    CALL  difference_int_0d(prof1%Lev2c%prh_temp_cpt_flag, prof2%Lev2c%prh_temp_cpt_flag, 1, ndiff, 'Lev2c%prh_temp_cpt_flag')

  ENDIF ! Lev2c

! 1.8 Level 2d
! ------------

  IF (Lev2d) THEN

    CALL difference_int_0d(prof1%Lev2d%Npoints, prof2%Lev2d%Npoints, 0, ndiff, 'Lev2d%Npoints')

    CALL difference_char_0d(prof1%Lev2d%level_type, prof2%Lev2d%level_type, ndiff, 'Lev2d%level_type')

    CALL difference_real_1d(prof1%Lev2d%level_coeff_a, prof2%Lev2d%level_coeff_a, 1.0e-6_wp, ndiff, 'Lev2d%level_coeff_a')
    CALL difference_real_1d(prof1%Lev2d%level_coeff_b, prof2%Lev2d%level_coeff_b, 1.0e-6_wp, ndiff, 'Lev2d%level_coeff_b')

  ENDIF ! Lev2d

! 1.9 Spectra
! -----------

  IF (spectra_local) THEN

! NB: need to nest these extra data searches

! 1.9.1 1D fields

    IF (ASSOCIATED(prof1%vlist%VlistD1d) .AND. ASSOCIATED(prof2%vlist%VlistD1d)) THEN

      IF (TRIM(ADJUSTL(prof1%vlist%VlistD1d%name)) == 'freq' .AND. &
          TRIM(ADJUSTL(prof2%vlist%VlistD1d%name)) == 'freq') THEN
        CALL difference_real_1d(prof1%vlist%VlistD1d%DATA, prof2%vlist%VlistD1d%DATA, &
                                1.0e-6_wp, ndiff, 'spectra_freq')
      ENDIF
      IF (TRIM(ADJUSTL(prof1%vlist%VlistD1d%name)) == 'stime' .AND. &
          TRIM(ADJUSTL(prof2%vlist%VlistD1d%name)) == 'stime') THEN
        CALL difference_real_1d(prof1%vlist%VlistD1d%DATA, prof2%vlist%VlistD1d%DATA, &
                                1.0e-6_wp, ndiff, 'spectra_time')
      ENDIF

      IF (ASSOCIATED(prof1%vlist%VlistD1d%next) .AND. ASSOCIATED(prof2%vlist%VlistD1d%next)) THEN

        IF (TRIM(ADJUSTL(prof1%vlist%VlistD1d%next%name)) == 'freq' .AND. &
            TRIM(ADJUSTL(prof2%vlist%VlistD1d%next%name)) == 'freq') THEN
          CALL difference_real_1d(prof1%vlist%VlistD1d%next%DATA, prof2%vlist%VlistD1d%next%DATA, &
                                  1.0e-6_wp, ndiff, 'spectra_freq')
        ENDIF
        IF (TRIM(ADJUSTL(prof1%vlist%VlistD1d%next%name)) == 'stime' .AND. &
            TRIM(ADJUSTL(prof2%vlist%VlistD1d%next%name)) == 'stime') THEN
          CALL difference_real_1d(prof1%vlist%VlistD1d%next%DATA, prof2%vlist%VlistD1d%next%DATA, &
                                  1.0e-6_wp, ndiff, 'spectra_time')
        ENDIF

      ENDIF ! Second 1d extra variable

    ENDIF ! First 1d extra variable

! 1.9.2 2D fields

    IF (ASSOCIATED(prof1%vlist%VlistD2d) .AND. ASSOCIATED(prof2%vlist%VlistD2d)) THEN

      IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%name)) == 'amp' .AND. &
          TRIM(ADJUSTL(prof2%vlist%VlistD2d%name)) == 'amp') THEN
        CALL difference_real_2d(prof1%vlist%VlistD2d%DATA, prof2%vlist%VlistD2d%DATA, &
                                1.0e-6_wp, ndiff, 'spectra_amp')
      ENDIF
      IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%name)) == 'bangle' .AND. &
          TRIM(ADJUSTL(prof2%vlist%VlistD2d%name)) == 'bangle') THEN
        CALL difference_real_2d(prof1%vlist%VlistD2d%DATA, prof2%vlist%VlistD2d%DATA, &
                                1.0e-6_wp, ndiff, 'spectra_bangle')
      ENDIF
      IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%name)) == 'impact' .AND. &
          TRIM(ADJUSTL(prof2%vlist%VlistD2d%name)) == 'impact') THEN
        CALL difference_real_2d(prof1%vlist%VlistD2d%DATA, prof2%vlist%VlistD2d%DATA, &
                                1.0e-6_wp, ndiff, 'spectra_impact')
      ENDIF

      IF (ASSOCIATED(prof1%vlist%VlistD2d%next) .AND. ASSOCIATED(prof2%vlist%VlistD2d%next)) THEN

        IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%next%name)) == 'amp' .AND. &
            TRIM(ADJUSTL(prof2%vlist%VlistD2d%next%name)) == 'amp') THEN
          CALL difference_real_2d(prof1%vlist%VlistD2d%next%DATA, prof2%vlist%VlistD2d%next%DATA, &
                                  1.0e-6_wp, ndiff, 'spectra_amp')
        ENDIF
        IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%next%name)) == 'bangle' .AND. &
            TRIM(ADJUSTL(prof2%vlist%VlistD2d%next%name)) == 'bangle') THEN
          CALL difference_real_2d(prof1%vlist%VlistD2d%next%DATA, prof2%vlist%VlistD2d%next%DATA, &
                                  1.0e-6_wp, ndiff, 'spectra_bangle')
        ENDIF
        IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%next%name)) == 'impact' .AND. &
            TRIM(ADJUSTL(prof2%vlist%VlistD2d%next%name)) == 'impact') THEN
          CALL difference_real_2d(prof1%vlist%VlistD2d%next%DATA, prof2%vlist%VlistD2d%next%DATA, &
                                  1.0e-6_wp, ndiff, 'spectra_impact')
        ENDIF

        IF (ASSOCIATED(prof1%vlist%VlistD2d%next%next) .AND. ASSOCIATED(prof2%vlist%VlistD2d%next%next)) THEN

          IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%next%next%name)) == 'amp' .AND. &
              TRIM(ADJUSTL(prof2%vlist%VlistD2d%next%next%name)) == 'amp') THEN
            CALL difference_real_2d(prof1%vlist%VlistD2d%next%next%DATA, prof2%vlist%VlistD2d%next%next%DATA, &
                                    1.0e-6_wp, ndiff, 'spectra_amp')
          ENDIF
          IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%next%next%name)) == 'bangle' .AND. &
              TRIM(ADJUSTL(prof2%vlist%VlistD2d%next%next%name)) == 'bangle') THEN
            CALL difference_real_2d(prof1%vlist%VlistD2d%next%next%DATA, prof2%vlist%VlistD2d%next%next%DATA, &
                                    1.0e-6_wp, ndiff, 'spectra_bangle')
          ENDIF
          IF (TRIM(ADJUSTL(prof1%vlist%VlistD2d%next%next%name)) == 'impact' .AND. &
              TRIM(ADJUSTL(prof2%vlist%VlistD2d%next%next%name)) == 'impact') THEN
            CALL difference_real_2d(prof1%vlist%VlistD2d%next%next%DATA, prof2%vlist%VlistD2d%next%next%DATA, &
                                    1.0e-6_wp, ndiff, 'spectra_impact')
          ENDIF

        ENDIF ! Third 2d extra variable

      ENDIF ! Second 2d extra variable

    ENDIF ! First 2d extra variable


  ENDIF ! spectra


! ------------------------------------------------------------------------------
! 2. Clean up
! ------------------------------------------------------------------------------

    CALL message_set_routine(routine)


CONTAINS


! ------------------------------------------------------------------------------
  SUBROUTINE difference_char_0d(var1, var2, ndiff, varnam)
!
! Compare character scalars
! -------------------------

! Input/output variables
    CHARACTER(LEN=*),       INTENT(in)     :: var1   ! 1st variable
    CHARACTER(LEN=*),       INTENT(in)     :: var2   ! 2nd variable
    INTEGER,                INTENT(inout)  :: ndiff  ! number of differences
    CHARACTER(LEN=*),       INTENT(in)     :: varnam ! variable name

! Local variables
    LOGICAL                                :: diff   ! maximum difference


    diff = TRIM(ADJUSTL(var1)) /= TRIM(ADJUSTL(var2))

    IF ( diff ) THEN
      CALL message ( msg_error, 'prof1%' // TRIM(ADJUSTL(varnam)) // ' differs from ' // &
                                'prof2%' // TRIM(ADJUSTL(varnam)) // '\n' // &
                                '(prof1%' // TRIM(ADJUSTL(varnam)) // ' = ' // TRIM(ADJUSTL(var1)) // ';\n' // &
                                ' prof2%' // TRIM(ADJUSTL(varnam)) // ' = ' // TRIM(ADJUSTL(var2)) // ')' )
      ndiff = ndiff + 1
    END IF

  END SUBROUTINE difference_char_0d
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
  SUBROUTINE difference_int_0d(var1, var2, tol, ndiff, varnam)
!
! Compare integer scalars
! -----------------------

! Input/output variables
    INTEGER,                INTENT(in)     :: var1   ! 1st variable
    INTEGER,                INTENT(in)     :: var2   ! 2nd variable
    INTEGER,                INTENT(in)     :: tol    ! tolerance
    INTEGER,                INTENT(inout)  :: ndiff  ! number of differences
    CHARACTER(LEN=*),       INTENT(in)     :: varnam ! variable name

! Local variables
    INTEGER                                :: diff   ! maximum difference
    CHARACTER(LEN=5)                       :: sdiff  ! maximum difference
    CHARACTER(LEN=5)                       :: stol   ! maximum difference


    diff = ABS(var1 - var2)

    IF ( diff > tol ) THEN
      WRITE (sdiff, '(i5.5)') diff
      WRITE (stol,  '(i5.5)') tol
      CALL message ( msg_error, 'prof1%' // TRIM(ADJUSTL(varnam)) // ' differs from ' // &
                                'prof2%' // TRIM(ADJUSTL(varnam)) // '\n' // &
                                '(|diff| = ' // sdiff // ' > ' // stol // ')' )
      ndiff = ndiff + 1
    END IF

  END SUBROUTINE difference_int_0d
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
  SUBROUTINE difference_real_0d(var1, var2, tol, ndiff, varnam)
!
! Compare real scalars
! --------------------

! Input/output variables
    REAL(wp),               INTENT(in)     :: var1   ! 1st variable
    REAL(wp),               INTENT(in)     :: var2   ! 2nd variable
    REAL(wp),               INTENT(in)     :: tol    ! tolerance
    INTEGER,                INTENT(inout)  :: ndiff  ! number of differences
    CHARACTER(LEN=*),       INTENT(in)     :: varnam ! variable name

! Local variables
    REAL(wp)                               :: diff   ! maximum difference
    CHARACTER(LEN=12)                      :: sdiff  ! maximum difference
    CHARACTER(LEN=12)                      :: stol   ! maximum difference

    diff = ABS(var1 - var2)

    IF ( diff > tol ) THEN
      WRITE (sdiff, '(e12.5)') diff
      WRITE (stol,  '(e12.5)') tol
      CALL message ( msg_error, 'prof1%' // TRIM(ADJUSTL(varnam)) // ' differs from ' // &
                                'prof2%' // TRIM(ADJUSTL(varnam)) // '\n' // &
                                '(|diff| = ' // sdiff // ' > ' // stol // ')' )
      ndiff = ndiff + 1
    END IF

  END SUBROUTINE difference_real_0d
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
  SUBROUTINE difference_real_1d(var1, var2, tol, ndiff, varnam)
!
! Compare real vectors
! --------------------

! Input/output variables
    REAL(wp), DIMENSION(:), INTENT(in)     :: var1   ! 1st variable
    REAL(wp), DIMENSION(:), INTENT(in)     :: var2   ! 2nd variable
    REAL(wp),               INTENT(in)     :: tol    ! tolerance
    INTEGER,                INTENT(inout)  :: ndiff  ! number of differences
    CHARACTER(LEN=*),       INTENT(in)     :: varnam ! variable name

! Local variables
    REAL(wp)                               :: diff   ! maximum difference
    CHARACTER(LEN=12)                      :: sdiff  ! maximum difference
    CHARACTER(LEN=12)                      :: stol   ! maximum difference

    diff = MAXVAL(ABS(var1 - var2))

    IF ( diff > tol ) THEN
      WRITE (sdiff, '(e12.5)') diff
      WRITE (stol,  '(e12.5)') tol
      CALL message ( msg_error, 'prof1%' // TRIM(ADJUSTL(varnam)) // ' differs from ' // &
                                'prof2%' // TRIM(ADJUSTL(varnam)) // '\n' // &
                                '(max|diff| = ' // sdiff // ' > ' // stol // ')' )
      ndiff = ndiff + 1
    END IF

  END SUBROUTINE difference_real_1d
! ------------------------------------------------------------------------------


! ------------------------------------------------------------------------------
  SUBROUTINE difference_real_2d(var1, var2, tol, ndiff, varnam)
!
! Compare real arrays
! -------------------

! Input/output variables
    REAL(wp), DIMENSION(:, :), INTENT(in)  :: var1   ! 1st variable
    REAL(wp), DIMENSION(:, :), INTENT(in)  :: var2   ! 2nd variable
    REAL(wp),               INTENT(in)     :: tol    ! tolerance
    INTEGER,                INTENT(inout)  :: ndiff  ! number of differences
    CHARACTER(LEN=*),       INTENT(in)     :: varnam ! variable name

! Local variables
    REAL(wp)                               :: diff   ! maximum difference
    CHARACTER(LEN=12)                      :: sdiff  ! maximum difference
    CHARACTER(LEN=12)                      :: stol   ! maximum difference

    diff = MAXVAL(ABS(var1 - var2))

    IF ( diff > tol ) THEN
      WRITE (sdiff, '(e12.5)') diff
      WRITE (stol,  '(e12.5)') tol
      CALL message ( msg_error, 'prof1%' // TRIM(ADJUSTL(varnam)) // ' differs from ' // &
                                'prof2%' // TRIM(ADJUSTL(varnam)) // '\n' // &
                                '(max|diff| = ' // sdiff // ' > ' // stol // ')' )
      ndiff = ndiff + 1
    END IF

  END SUBROUTINE difference_real_2d
! ------------------------------------------------------------------------------


END SUBROUTINE ropp_io_fields_compare
