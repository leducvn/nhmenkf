! $Id: ropp_pp_tph_refrac.f90 3491 2013-02-06 12:43:43Z idculv $

SUBROUTINE ropp_pp_tph_refrac(ro_data, diag)

!****s* TropoPauseHeight/ropp_pp_tph_refrac *
!
! NAME
!   ropp_pp_tph_refrac
!
! SYNOPSIS
!   Tropopause height diagnostic based on refractivity
!
!   CALL ropp_pp_tph_refrac(ro_data, diag)
!
! DESCRIPTION
!   Diagnose tropopause height from the kinks in the log(refractivity) profile,
!   using the covariance transform method of Lewis (GRL 2009).
!
! INPUTS
!   TYPE(ROprof), INTENT(INOUT)    :: ro_data   ! input RO profile containing lev2a data
!   LOGICAL, OPTIONAL, INTENT(IN)  :: diag       ! extra diagnostics required
!
! OUTPUTS
!   TYPE(ROprof), INTENT(INOUT)    :: ro_data   ! output RO profile containing lev2a data
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
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
! USE ropp_pp, not_this => ropp_pp_tph_refrac
  USE ropp_pp

  IMPLICIT NONE

  TYPE(ROprof), INTENT(INOUT)          :: ro_data
  LOGICAL, OPTIONAL, INTENT(IN)        :: diag       ! extra diagnostics required

! Parameters in the refractivity TPH diagnosis algorithm

  REAL(wp), PARAMETER                  :: data_min_height=15.0e3_wp
  REAL(wp), PARAMETER                  :: data_max_height=30.0e3_wp

  REAL(wp), PARAMETER                  :: tph_min_const1=7.5e3_wp
  REAL(wp), PARAMETER                  :: tph_min_const2=2.5e3_wp
  REAL(wp)                             :: tph_min

  REAL(wp), PARAMETER                  :: tph_max_const1=17.5e3_wp
  REAL(wp), PARAMETER                  :: tph_max_const2=2.5e3_wp
  REAL(wp)                             :: tph_max

  REAL(wp), PARAMETER                  :: tph_dbl_range=2.0e3_wp
  REAL(wp), PARAMETER                  :: tph_cov_width=25.0e3_wp
  REAL(wp), PARAMETER                  :: tph_sharpness=1.05_wp
  REAL(wp), PARAMETER                  :: tph_sharpness1=1.10_wp
  REAL(wp), PARAMETER                  :: tph_sharp_width=5.0e3_wp
  REAL(wp), PARAMETER                  :: tph_low_check=10.0e3_wp
  REAL(wp), PARAMETER                  :: refrac_norm=1.0e3_wp
  REAL(wp), PARAMETER                  :: pi1=3.141592653589793238_wp
  REAL(wp), PARAMETER                  :: dtor=pi1/180.0_wp
  INTEGER,  PARAMETER                  :: min_no_points=2

! Local variables

  REAL(wp)                             :: max_cov,max_cov1
  REAL(wp)                             :: avg_cov,avg_cov1
  REAL(wp)                             :: tph,tpn
  REAL(wp)                             :: lat,cos2lat

  REAL(wp), DIMENSION(:), ALLOCATABLE  :: var_x,var_y,var_cov

  INTEGER                              :: tph_index,tph_index1
  INTEGER                              :: tph_qc_flag
  INTEGER                              :: n_points,n_points1
  INTEGER                              :: i,i1

  LOGICAL                              :: tph_possible
  LOGICAL                              :: l_diag

  CHARACTER(len=7)                     :: str_tph_min,str_tph_max
  CHARACTER(len=256)                   :: routine

  LOGICAL, DIMENSION(:), ALLOCATABLE   :: lregion,lregion1


!-------------------------------------------------------------------------------
! 2. Calculate TPH based on covariance transform of log(refractivity)
!-------------------------------------------------------------------------------

! 2.0 Initial messages
! --------------------

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_tph_refrac')

  CALL message(msg_info, "Calculating refractivity-based tropopause height \n")

! 2.1 Initialise variables
! ------------------------

  n_points = ro_data%lev2a%npoints

  IF ((n_points == 0) .OR. ro_data%lev2a%missing) THEN
    CALL message(msg_info, "No (valid) Level 2a data in profile ... " //    &
                           "will not calculate TPH based on refractivity \n")
    CALL message_set_routine(routine)
    RETURN
  END IF

  tph = ropp_MDFV  ;  tpn = ropp_MDFV  ;  tph_qc_flag = 0

  tph_possible = .TRUE.

  IF (PRESENT(diag)) THEN
    l_diag = diag
  ELSE
    l_diag = .FALSE.
  ENDIF

! 2.2 Check for numerically invalid data
! --------------------------------------

  IF (ANY((ro_data%lev2a%refrac <= ropp_ZERO) .AND. &
          (ro_data%lev2a%refrac >= ropp_MDTV))) THEN
    CALL message(msg_diag, "Negative refractivities encountered ... " // &
                           "using log of absolute value \n")
  END IF

  IF (ALL((ro_data%lev2a%refrac     < ropp_MDTV) .OR. &
          (ro_data%lev2a%alt_refrac < ropp_MDTV))) THEN
    CALL message(msg_error, "No common non-missing refractivities and " // &
                            "refractivity altitudes ... cannot calculate TPH\n")
    tph_possible = .FALSE.
    tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_data_invalid)
  END IF

  IF (ro_data%GEOref%lat < ropp_MDTV) THEN
    CALL message(msg_error, "Missing tangent point latitude ... " // &
                            "cannot calculate latitude-dependent parameters \n")
    tph_possible = .FALSE.
    tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_data_invalid)
  ELSE
    lat = ro_data%GEOref%lat
  END IF

  ALLOCATE(lregion(n_points))  ;  lregion = .FALSE.

  WHERE ((ro_data%lev2a%refrac     > ropp_MDTV) .AND. &
         (ro_data%lev2a%alt_refrac > ropp_MDTV)) lregion = .TRUE.
  n_points1 = COUNT(lregion)

  IF (n_points1 < min_no_points) THEN
    CALL message(msg_error, "Not enough non-missing points for " // &
                            "algorithm to proceed ... cannot calculate TPH\n")
    tph_possible = .FALSE.
    tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_data_invalid)
  END IF

  IF (BTEST(tph_qc_flag, TPH_QC_data_invalid)) THEN
    CALL message(msg_info, "Numerically invalid data")
  END IF

! 2.3 Check for scientifically invalid data
! -----------------------------------------

  IF (ro_data%GEOref%lat < ropp_MDTV) THEN
    tph_min = tph_min_const1 - tph_min_const2  ! Value at pole
    tph_max = tph_max_const1 + tph_max_const2  ! Value at equator
  ELSE
    cos2lat = COS(2.0_wp*dtor*lat)
    tph_min = tph_min_const1 + tph_min_const2*cos2lat
    tph_max = tph_max_const1 + tph_max_const2*cos2lat
  END IF

  IF (MINVAL(ro_data%lev2a%alt_refrac, mask=lregion) > data_min_height) THEN
    CALL message(msg_info, "Refractivity altitudes do not start deep enough")
    tph_possible = .FALSE.
    tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_prof_depth)
  END IF

  IF (MAXVAL(ro_data%lev2a%alt_refrac, mask=lregion) < data_max_height) THEN
    CALL message(msg_info, "Impact altitudes do not reach high enough")
    tph_possible = .FALSE.
    tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_prof_height)
  END IF

  IF (BTEST(tph_qc_flag,  TPH_QC_prof_depth) .OR. &
      BTEST(tph_qc_flag, TPH_QC_prof_height)) THEN
    CALL message(msg_info, "Scientifically invalid data")
  END IF

! 2.4 Calculate {alt refrac, log(refrac)} for valid data points
! -------------------------------------------------------------

  IF (tph_possible) THEN

    ALLOCATE(var_y(n_points1), var_x(n_points1))  ! Fill these with the non-missing data

    i1 = 1
    DO i=1,n_points
      IF (lregion(i)) THEN
        var_y(i1) = ro_data%lev2a%alt_refrac(i)
        var_x(i1) = LOG(ABS(ro_data%lev2a%refrac(i)/refrac_norm)+TINY(1.0_wp)) ! Assume refrac in N-units
        i1        = i1 + 1
      END IF
    END DO

    DEALLOCATE(lregion)

! Below this point, the length of all arrays should be n_points1

    ALLOCATE(lregion(n_points1), lregion1(n_points1))  ! Masks

! 2.5 Compute covariance transform of log(refrac) and find its maximum
! --------------------------------------------------------------------

    ALLOCATE(var_cov(n_points1)) ! Will hold covariance transform of var_x

    lregion = .FALSE.
    WHERE ((var_y < tph_max) .AND. &
           (var_y > tph_min)) lregion = .TRUE.

    IF (ANY(lregion)) THEN

      CALL message(msg_diag, "Calculating covariance transform")

      CALL ropp_pp_cov_transform(var_y, var_x, tph_cov_width, var_cov)

      IF (l_diag) THEN

! Append covariance transform of refractivity and its altitude to the 
! ROprof structure and thence the output file.

        CALL ropp_io_addvar_rodataD1d( &
             ro_data, &
             name      = "CT_refrac", &
             long_name = "Covariance transform of refractivity", &
             units     = "", &
             range     = (/-1000000.0_wp, 1000000.0_wp/), &
             DATA      = var_cov)

        CALL ropp_io_addvar_rodataD1d( &
             ro_data,   &
             name      = "CT_refrac_alt", &
             long_name = "Covariance transform of refractivity altitude", &
             units     = "m", &
             range     = (/-1000000.0_wp, 1000000.0_wp/), &
             DATA      = var_y)

      END IF

      max_cov = MAXVAL(var_cov, mask=lregion)

      tph_index = SUM(MINLOC(ABS(var_cov-max_cov), mask=lregion))

      tph = var_y(tph_index)

      tpn = refrac_norm*EXP(var_x(tph_index))  ! Refractivity at tropopause

! 2.6 Check for over-smoothness above TPH
! ---------------------------------------

      CALL message(msg_diag, "Checking over-smoothness above TPH")

      lregion1 = .FALSE.
      WHERE (((var_y-tph) >        ropp_ZERO) .AND. &
             ((var_y-tph) <= tph_sharp_width)) lregion1 = .TRUE.

      IF (ANY(lregion1)) THEN

        avg_cov = SUM(var_cov, mask=lregion1) / COUNT(lregion1)
        IF (max_cov <= (tph_sharpness*avg_cov)) THEN
          CALL message(msg_info, "Derived TPH not sharp enough above")
          tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_CT_smooth_above)
        END IF

      END IF

! 2.7 Check for over-smoothness below TPH
! ---------------------------------------

      CALL message(msg_diag, "Checking over-smoothness below TPH")

      lregion1 = .FALSE.
      WHERE (((var_y-tph) >  (-1.0_wp*tph_sharp_width)) .AND. &
             ((var_y-tph) <=                 ropp_ZERO)) lregion1 = .TRUE.

      IF (ANY(lregion1)) THEN

        avg_cov = SUM(var_cov, mask=lregion1) / COUNT(lregion1)
        IF (max_cov <= (tph_sharpness*avg_cov)) THEN
          CALL message(msg_info, "Derived TPH not sharp enough below")
          tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_CT_smooth_below)
        END IF

      END IF

    ELSE  ! Not possible to calculate CT

      tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_data_invalid)

    END IF

! 2.8 Check low TPH for possible double tropopause
! ------------------------------------------------

    IF ((.NOT. BTEST(tph_qc_flag, TPH_QC_data_invalid)) .AND. &
        (tph < tph_low_check)) THEN

      CALL message(msg_diag,"Checking low TPH for possible double tropopause")

      lregion1 = .FALSE.
      WHERE ((var_y < tph_max) .AND. &
             (var_y > tph_min) .AND. &
             (var_y > (tph+tph_dbl_range))) lregion1 = .TRUE.

      IF (ANY(lregion1)) THEN

        max_cov1 = MAXVAL(var_cov, mask=lregion1)

        tph_index1 = SUM(MINLOC(ABS(var_cov-max_cov1), mask=lregion1))

        IF (COUNT(ABS(var_y-var_y(tph_index1))<tph_dbl_range) > 0) THEN

          avg_cov1 = SUM(var_cov, &
                     mask=ABS(var_y-var_y(tph_index1))<tph_dbl_range) / &
                     COUNT(ABS(var_y-var_y(tph_index1))<tph_dbl_range)

          IF (var_cov(tph_index1) > (tph_sharpness*avg_cov1)) THEN  ! It's a reasonably well defined max

            CALL message(msg_info, "Possible double tropopause detected " // &
                                   "... comparing to lower one")
            IF (var_cov(tph_index) < (tph_sharpness1*var_cov(tph_index1))) THEN

              CALL message(msg_info, "Double tropopause detected")
              tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_double_trop)

            END IF

          END IF

        END IF

      END IF

    END IF

! 2.9 Check that TPH is not too low
! ---------------------------------

    IF (.NOT. BTEST(tph_qc_flag, TPH_QC_data_invalid)) THEN

      CALL message(msg_diag, "Checking that the TPH is not too low")

      IF ((tph < tph_min) .AND. (tph > ropp_MDTV)) THEN
        WRITE(str_tph_min, '(F7.3)') 1.0e-3_wp*tph_min
        CALL message(msg_info, "Derived TPH below min of " // str_tph_min // " km")
        tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_too_low)
      END IF

    END IF

! 2.10 Check that TPH is not too high
! -----------------------------------

    IF (.NOT. BTEST(tph_qc_flag, TPH_QC_data_invalid)) THEN

      CALL message(msg_diag, "Checking that the TPH is not too high")

      IF ((tph > tph_max) .AND. (tph > ropp_MDTV)) THEN
        WRITE(str_tph_max, '(F7.3)') 1.0e-3_wp*tph_max
        CALL message(msg_info, "Derived TPH above max of " // str_tph_max // " km")
        tph_qc_flag  = IBSET(tph_qc_flag, TPH_QC_too_high)
      END IF

    END IF

! 2.11 Clean up
! -------------

    DEALLOCATE(lregion, lregion1, var_y, var_x, var_cov)

  ELSE

    DEALLOCATE(lregion)

  END IF ! tph_possible = .TRUE.


!-------------------------------------------------------------------------------
! 3. Copy TPH variables to ROprof structure
!-------------------------------------------------------------------------------

  ro_data%lev2c%tph_refrac      = tph  ! Return a refractivity altitude
  ro_data%lev2c%tpn_refrac      = tpn  ! Return a refractivity
  ro_data%lev2c%tph_refrac_flag = tph_qc_flag


!-------------------------------------------------------------------------------
! 4. Reset routine name
!-------------------------------------------------------------------------------

  CALL message_set_routine(routine)


END SUBROUTINE ropp_pp_tph_refrac
