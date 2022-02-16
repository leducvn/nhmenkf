! $Id: ropp_pp_tph_tdry.f90 3491 2013-02-06 12:43:43Z idculv $

SUBROUTINE ropp_pp_tph_tdry(ro_data, diag)

!****s* TropoPauseHeight/ropp_pp_tph_tdry *
!
! NAME
!   ropp_pp_tph_tdry
!
! SYNOPSIS
!   Tropopause height diagnostic based on dry temperature
!
!   CALL ropp_pp_tph_tdry(ro_data, diag)
!
! DESCRIPTION
!   Diagnose tropopause height from the kinks in the dry temperature profile,
!   using the lapse rate and cold point criteria (Reichler et al, GRL, 2003).
!   This algorithm is written in terms of pressure.  We use the dry pressure, 
!   which is calculable from the dry temperature and refractivity.  The 
!   latter should therefore also be present in the file, but if it isn't,
!   a crude estimate is made (and a warning is broadcast).
!
! INPUTS
!   TYPE(ROprof), INTENT(INOUT)    :: ro_data    ! input RO profile containing lev2a data
!   LOGICAL, OPTIONAL, INTENT(IN)  :: diag       ! extra diagnostics required
!
!
! OUTPUTS
!   TYPE(ROprof), INTENT(INOUT)    :: ro_data    ! output RO profile containing lev2a data
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
! USE ropp_pp, not_this => ropp_pp_tph_tdry
  USE ropp_pp
  USE ropp_pp_constants, ONLY: R_dry, C_p, g_wmo

  IMPLICIT NONE

  TYPE(ROprof), INTENT(INOUT)          :: ro_data    ! input RO profile containing lev2a data
  LOGICAL, OPTIONAL, INTENT(IN)        :: diag       ! extra diagnostics required

! Parameters in the dry temperature-based TPH diagnosis algorithms

  REAL(wp), PARAMETER                  :: kappa=R_dry/C_p
  REAL(wp), PARAMETER                  :: p_ref=1000.0_wp ! hPa

  REAL(wp), PARAMETER                  :: tph_min_const1=7.5e3_wp
  REAL(wp), PARAMETER                  :: tph_min_const2=2.5e3_wp
  REAL(wp)                             :: tph_min

  REAL(wp), PARAMETER                  :: tph_max_const1=17.5e3_wp
  REAL(wp), PARAMETER                  :: tph_max_const2=2.5e3_wp
  REAL(wp)                             :: tph_max

  REAL(wp), PARAMETER                  :: wmo_lapse=-2.0e-3_wp
  REAL(wp), PARAMETER                  :: tph_2km=2.0e3_wp

  REAL(wp), PARAMETER                  :: tph_diff_max=2.0e3_wp
  REAL(wp), PARAMETER                  :: tropic_lat_bdy=30.0_wp

  REAL(wp), PARAMETER                  :: pi1=3.141592653589793238_wp
  REAL(wp), PARAMETER                  :: dtor=pi1/180.0_wp
  INTEGER,  PARAMETER                  :: min_no_points=3

! Local variables

  REAL(wp)                             :: tpt_lrt,tpt_cpt,prt_cpt
  REAL(wp)                             :: tph_lrt,tph_cpt,prh_cpt
  REAL(wp)                             :: p_tph,pi_tph
  REAL(wp)                             :: var_yi
  REAL(wp)                             :: avg_lapse
  REAL(wp)                             :: lat,cos2lat,interp_frac
  REAL(wp)                             :: var_1,var_n

  REAL(wp), DIMENSION(:), ALLOCATABLE  :: var_x,var_y,var_p,&
                                          var_lapse,var_pi,&
                                          var_y1,var_p1,&
                                          var_x_raw

  LOGICAL                              :: tph_possible
  LOGICAL                              :: l_diag

  INTEGER                              :: lrt_index,tph_lrt_qc_flag
  INTEGER                              :: cpt_index,tph_cpt_qc_flag
  INTEGER                              :: prf_index,prh_cpt_qc_flag
  INTEGER                              :: n_points,n_points1
  INTEGER                              :: i,i1

  CHARACTER(len=7)                     :: str_tph_min,str_tph_max
  CHARACTER(len=256)                   :: routine

  LOGICAL, DIMENSION(:), ALLOCATABLE   :: lregion,lregion1

  LOGICAL                              :: l_lapse


!-------------------------------------------------------------------------------
! 2. Calculate TPH based on lapse rate (wrt dry pressure) and cold point of dry temperature
!-------------------------------------------------------------------------------

! 2.0 Initial messages
! --------------------

  CALL message_get_routine(routine)

  CALL message_set_routine('ropp_pp_tph_tdry')

  CALL message(msg_info, "Calculating dry temperature-based tropopause height \n")

! 2.1 Initialise variables
! ------------------------

  n_points = ro_data%lev2a%npoints

  IF ((n_points == 0) .OR. ro_data%lev2a%missing) THEN
    CALL message(msg_info, "No (valid) Level 2a data in profile ... " //    &
                           "will not calculate TPH based on dry temperature \n")
    CALL message_set_routine(routine)
    RETURN
  END IF

  tph_lrt = ropp_MDFV  ;  tpt_lrt = ropp_MDFV  ;  tph_lrt_qc_flag = 0
  tph_cpt = ropp_MDFV  ;  tpt_cpt = ropp_MDFV  ;  tph_cpt_qc_flag = 0
  prh_cpt = ropp_MDFV  ;  prt_cpt = ropp_MDFV  ;  prh_cpt_qc_flag = 0

  tph_possible = .TRUE.

  IF (PRESENT(diag)) THEN
    l_diag = diag
  ELSE
    l_diag = .FALSE.
  ENDIF

! 2.2 Check for numerically invalid data
! --------------------------------------

  IF (ALL((ro_data%lev2a%dry_temp   < ropp_MDTV) .OR. &
          (ro_data%lev2a%alt_refrac < ropp_MDTV))) THEN
    CALL message(msg_error, "No common non-missing dry temperatures and " // &
                            "refractivity altitudes ... cannot calculate TPH \n")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_data_invalid)
  END IF

  IF (ro_data%GEOref%lat < ropp_MDTV) THEN
    CALL message(msg_error, "Missing tangent point latitude ... " // &
                            "cannot calculate latitude-dependent parameters \n")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_data_invalid)
  ELSE
    lat = ro_data%GEOref%lat
  END IF

  ALLOCATE(lregion(n_points))  ;  lregion = .FALSE.

  WHERE ((ro_data%lev2a%dry_temp   > ropp_MDTV) .AND. &
         (ro_data%lev2a%alt_refrac > ropp_MDTV)) lregion = .TRUE.
  n_points1 = COUNT(lregion)

  IF (n_points1 < min_no_points) THEN
    CALL message(msg_error, "Not enough non-missing points for algorithm " // &
                            "to proceed ... cannot calculate lapse rate TPH\n")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_data_invalid)
  END IF

  IF (BTEST(tph_lrt_qc_flag, TPH_QC_data_invalid)) THEN
    CALL message(msg_info, "Numerically invalid data")
  END IF

! 2.3 Check for scientifically invalid data
! -----------------------------------------

  IF (ro_data%GEOref%lat < ropp_MDTV) THEN
    tph_min = tph_min_const1 - tph_min_const2  ! Value at pole
    tph_max = tph_max_const1 + tph_max_const2  ! Value on equator
  ELSE
    cos2lat = COS(2.0_wp*dtor*lat)
    tph_min = tph_min_const1 + tph_min_const2*cos2lat
    tph_max = tph_max_const1 + tph_max_const2*cos2lat
  END IF

  IF (MINVAL(ro_data%lev2a%alt_refrac, mask=lregion) > tph_min) THEN
    CALL message(msg_info, "Refractivity altitudes do not start deep enough")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_prof_depth)
  END IF

  IF (MAXVAL(ro_data%lev2a%alt_refrac, mask=lregion) < tph_max) THEN
    CALL message(msg_info, "Refractivity altitudes do not reach high enough")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_prof_height)
  END IF

  IF (BTEST(tph_lrt_qc_flag,  TPH_QC_prof_depth) .OR. &
      BTEST(tph_lrt_qc_flag, TPH_QC_prof_height)) THEN
    CALL message(msg_info, "Scientifically invalid data")
  END IF

! 2.4 Prune out missing data
! --------------------------

  IF (tph_possible) THEN

    ALLOCATE(var_y(n_points1), var_p(n_points1), var_x(n_points1), var_x_raw(n_points1))  ! Fill these with the non-missing data

    i1 = 1
    DO i=1,n_points
      IF (lregion(i)) THEN
        var_y(i1) = ro_data%lev2a%alt_refrac(i)
        var_x(i1) = ro_data%lev2a%dry_temp(i)
        i1        = i1 + 1
      END IF
    END DO

! 2.5 Calculate dry pressure pdry = N Tdry / kappa1, if possible; 
!     otherwise estimated from Tdry profile
! -----------------------------------------

    CALL calc_pdry(ro_data%lev2a%refrac, lregion, var_x, var_y, var_p)

    DEALLOCATE(lregion)

! 2.6 Smooth Tdry and pdry
! ------------------------

    var_x_raw = var_x  ! Needed for CPT and PRH calculations later

    var_1 = var_x(1)   ;  var_n = var_x(n_points1)
    var_x = (CSHIFT(var_x, -1) + var_x + CSHIFT(var_x, 1)) / 3.0_wp
    var_x(1) = var_1   ;   var_x(n_points1) = var_n

    var_1 = var_p(1)   ;  var_n = var_p(n_points1)
    var_p = (CSHIFT(var_p, -1) + var_p + CSHIFT(var_p, 1)) / 3.0_wp
    var_p(1) = var_1   ;   var_p(n_points1) = var_n

! Below this point, the length of all arrays should be n_points1

    ALLOCATE(lregion(n_points1), lregion1(n_points1))  ! Masks

! 2.7 Calculate Exner pressure pi = (p/p_ref)**kappa
! --------------------------------------------------

    ALLOCATE(var_pi(n_points1))  !  Will hold Exner pressure Pi

    var_pi = (var_p / p_ref) ** kappa

! 2.8 Compute lapse rate at half points (i+1/2), assuming T varies linearly with pi
! ---------------------------------------------------------------------------------

    ALLOCATE(var_lapse(n_points1)) !  Will hold dT/dz (ie, -lapse rate)

    var_lapse = (CSHIFT(var_x,  1) - var_x) * &
                (CSHIFT(var_pi, 1) + var_pi) / &
                (CSHIFT(var_pi, 1) - var_pi + TINY(1.0_wp)) / &
                (CSHIFT(var_x,  1) + var_x  + TINY(1.0_wp))

    var_lapse(n_points1) = var_lapse(n_points1-1)

    var_lapse = -1.0_wp * var_lapse * (g_wmo/C_p)   ! var_lapse(i) = +dT/dz(i+1/2)

    IF (l_diag) THEN

! Append lapse rate, pressure and altitude to the ROprof structure and thence the output file.

      CALL ropp_io_addvar_rodataD1d( &
           ro_data, &
           name      = "LR_tdry", &
           long_name = "Dry temperature lapse rate", &
           units     = "K/km", &
           range     = (/-1000000.0_wp, 1000000.0_wp/), &
           DATA      = var_lapse)

      ALLOCATE(var_p1(n_points1))  ! This is where lapse rate is defined (approximately)
      var_p1 = (var_p + CSHIFT(var_p, 1)) / 2.0_wp
      var_p1(n_points1) = var_p(n_points1)
      CALL ropp_io_addvar_rodataD1d( &
           ro_data, &
           name      = "LR_tdry_press", &
           long_name = "Dry temperature lapse rate pressure", &
           units     = "hPa", &
           range     = (/-1000000.0_wp, 1000000.0_wp/), &
           DATA      = var_p1)
      DEALLOCATE(var_p1)

      ALLOCATE(var_y1(n_points1))  ! This is where lapse rate is defined (approximately)
      var_y1 = (var_y + CSHIFT(var_y, 1)) / 2.0_wp
      var_y1(n_points1) = var_y(n_points1)
      CALL ropp_io_addvar_rodataD1d( &
           ro_data, &
           name      = "LR_tdry_alt", &
           long_name = "Dry temperature lapse rate altitude", &
           units     = "m", &
           range     = (/-1000000.0_wp, 1000000.0_wp/), &
           DATA      = var_y1)
      DEALLOCATE(var_y1)

    END IF

! 2.9 Calculate TPH according to Reichler's method
! ------------------------------------------------

    l_lapse = .FALSE.

    lrt_index = ropp_MIFV

    DO i=2,n_points1-1 ! Need lrt_index > 1 and < npoints1 in next section

      IF ((var_lapse(i) < wmo_lapse  ) .OR. &
          (var_y(i)     < tph_min) .OR. &
          l_lapse) CYCLE

        var_yi = var_y(i)

        lregion = .FALSE.
        WHERE ((var_y >            var_yi)  .AND. &
               (var_y <= (var_yi+tph_2km))) lregion = .TRUE.

        IF (ANY(lregion)) THEN ! Calculate avg lapse rate in 2km above point

          avg_lapse = SUM(var_lapse, mask=lregion) / COUNT(lregion)

          IF ((avg_lapse > wmo_lapse) .AND. (var_lapse(i-1) < wmo_lapse)) THEN

            lrt_index = i
            l_lapse   = .TRUE.

          END IF

        END IF

    END DO

    IF (.NOT. l_lapse) THEN

      CALL message(msg_info, "Could not calculate lapse rate TPH ... " // &
                             "leaving as missing data \n")
      tph_lrt_qc_flag = IBSET(tph_lrt_qc_flag, TPH_QC_data_invalid)

    END IF

! 2.10 Calculate precise location of TPH by interpolation
! -------------------------------------------------------

    IF (.NOT. BTEST(tph_lrt_qc_flag, TPH_QC_data_invalid)) THEN

      IF ((var_lapse(lrt_index-1)-wmo_lapse) * &               ! We know lrt_index >= 2
          (var_lapse(lrt_index  )-wmo_lapse) < ropp_ZERO) THEN ! Interpolation should work

        CALL message(msg_diag, "Calculating TPH by interpolation")

        pi_tph =  0.5_wp*(var_pi(lrt_index-1) + var_pi(lrt_index  )) + &
                 (0.5_wp*(var_pi(lrt_index+1) - var_pi(lrt_index-1)) * &
                 (wmo_lapse            - var_lapse(lrt_index-1))     / &
                 (var_lapse(lrt_index) - var_lapse(lrt_index-1)))

        p_tph = p_ref * (pi_tph)**(1.0_wp/kappa)

! Around the tropopause, the temperature is by definition slowly varying,
! so we interpolate z linearly with ln(p), as for an isothermal atmosphere

        interp_frac = LOG(p_tph           /var_p(lrt_index-1)) / &
                      LOG(var_p(lrt_index)/var_p(lrt_index-1))

        tph_lrt = var_y(lrt_index-1) + &
                  (var_y(lrt_index) - var_y(lrt_index-1))*interp_frac

! Interpolate T linearly with z.  Since we're assuming z varies linearly with 
! log p, this implies T also varies linearly with log p, so we can use the same 
! interpolation coefficient.

        tpt_lrt = var_x(lrt_index-1) + &
                  (var_x(lrt_index) - var_x(lrt_index-1))*interp_frac

      ELSE

        CALL message(msg_info, "Could not interpolate to TPH ... " // &
                               "using nearest point \n")
        tph_lrt = var_y(lrt_index)
        tpt_lrt = var_x(lrt_index)

      END IF

    END IF

! 2.11 Calculate cold point tropopause if in tropics
! --------------------------------------------------

! Note that we need tph_lrt in this algorithm, so we cannot bypass the earlier calculation

! Cold point TPH should be defined in terms of the unsmoothed temperature var_x_raw

    IF ((.NOT. BTEST(tph_lrt_qc_flag, TPH_QC_data_invalid)) .AND. &
        (.NOT. BTEST(tph_cpt_qc_flag, TPH_QC_data_invalid))) THEN

      IF (ABS(lat) <= tropic_lat_bdy) THEN  ! In tropics

        lregion = .FALSE.
        WHERE ((var_y > tph_min) .AND. (var_y < tph_max)) lregion = .TRUE.

        IF (ANY(lregion)) THEN

          tpt_cpt   = MINVAL(var_x_raw, mask=lregion)
          cpt_index = SUM(MINLOC(ABS(var_x_raw-tpt_cpt), mask=lregion))
          tph_cpt   = var_y(cpt_index)

          IF (ABS(tph_cpt - tph_lrt) > tph_diff_max) THEN

            lregion1 = .FALSE.
            WHERE (ABS(var_y - tph_lrt) < tph_diff_max) lregion1 = .TRUE.

            IF (ANY(lregion1)) THEN

              tpt_cpt   = MINVAL(var_x_raw, mask=lregion1)
              cpt_index = SUM(MINLOC(ABS(var_x_raw-tpt_cpt), mask=lregion1))
              tph_cpt   = var_y(cpt_index)

            ELSE

              CALL message(msg_info, "Could not calculate cold point TPH ... " // &
                                     "leaving as missing data \n")
              tph_cpt_qc_flag = IBSET(tph_cpt_qc_flag, TPH_QC_data_invalid)

            END IF

          END IF

        ELSE

          CALL message(msg_info, "Could not calculate cold point TPH ... " // &
                                 "leaving as missing data \n")
          tph_cpt_qc_flag = IBSET(tph_cpt_qc_flag, TPH_QC_data_invalid)

        END IF

      ELSE

        CALL message(msg_info, "Outside tropics ... " // &
                               "not calculating cold point TPH \n")
        tph_cpt_qc_flag = IBSET(tph_cpt_qc_flag, TPH_QC_data_invalid)

      END IF ! In tropics

    ELSE

      CALL message(msg_info, "Could not calculate cold point TPH ... " // &
                                 "leaving as missing data \n")
      tph_cpt_qc_flag = IBSET(tph_cpt_qc_flag, TPH_QC_data_invalid)

    END IF

! 2.12 Check that TPH is not too low
! ----------------------------------

    IF (.NOT. BTEST(tph_lrt_qc_flag, TPH_QC_data_invalid)) THEN

      CALL message(msg_diag, "Checking that the lapse rate TPH is not too low")

      IF ((tph_lrt < tph_min) .AND. (tph_lrt > ropp_MDTV)) THEN
        WRITE(str_tph_min, '(F7.3)') 1.0e-3_wp*tph_min
        CALL message(msg_info, "Derived lapse rate TPH below min of " // str_tph_min // " km")
        tph_lrt_qc_flag = IBSET(tph_lrt_qc_flag, TPH_QC_too_low)
      END IF

    END IF

    IF (.NOT. BTEST(tph_cpt_qc_flag, TPH_QC_data_invalid)) THEN

      CALL message(msg_diag, "Checking that the cold point TPH is not too low")

      IF ((tph_cpt < tph_min) .AND. (tph_cpt > ropp_MDTV)) THEN
        WRITE(str_tph_min, '(F7.3)') 1.0e-3_wp*tph_min
        CALL message(msg_info, "Derived cold point TPH below min of " // str_tph_min // " km")
        tph_cpt_qc_flag = IBSET(tph_cpt_qc_flag, TPH_QC_too_low)
      END IF

    END IF

! 2.13 Check that TPH is not too high
! -----------------------------------

    IF (.NOT. BTEST(tph_lrt_qc_flag, TPH_QC_data_invalid)) THEN

      CALL message(msg_diag, "Checking that the lapse rate TPH is not too high")

      IF ((tph_lrt > tph_max) .AND. (tph_lrt > ropp_MDTV)) THEN
        WRITE(str_tph_max, '(F7.3)') 1.0e-3_wp*tph_max
        CALL message(msg_info, "Derived lapse rate TPH above max of " // str_tph_max // " km")
        tph_lrt_qc_flag = IBSET(tph_lrt_qc_flag, TPH_QC_too_high)
      END IF

    END IF

    IF (.NOT. BTEST(tph_cpt_qc_flag, TPH_QC_data_invalid)) THEN

      CALL message(msg_diag, "Checking that the cold point TPH is not too high")

      IF ((tph_cpt > tph_max) .AND. (tph_cpt > ropp_MDTV)) THEN
        WRITE(str_tph_max, '(F7.3)') 1.0e-3_wp*tph_max
        CALL message(msg_info, "Derived cold point TPH above max of " // str_tph_max // " km")
        tph_cpt_qc_flag = IBSET(tph_cpt_qc_flag, TPH_QC_too_high)
      END IF

    END IF

! 2.14 Calculate and store profile minimum temp and its height
! ------------------------------------------------------------

! Profile cold point should be defined in terms of the unsmoothed temperature var_x_raw

    IF (.NOT. BTEST(prh_cpt_qc_flag, TPH_QC_data_invalid)) THEN

      lregion = .FALSE.
      WHERE (var_x_raw > ropp_MDTV) lregion = .TRUE.

      IF (ANY(lregion)) THEN

        prt_cpt   = MINVAL(var_x_raw, mask=lregion)
        prf_index = SUM(MINLOC(ABS(var_x_raw-prt_cpt), mask=lregion))
        prh_cpt   = var_y(prf_index)

      ELSE

        CALL message(msg_info, "Could not calculate entire profile cold point ... " // &
                               "leaving as missing data \n")
        prh_cpt_qc_flag = IBSET(prh_cpt_qc_flag, TPH_QC_data_invalid)

      END IF

    END IF

! 2.15 Clean up
! -------------

    DEALLOCATE(lregion, lregion1, var_y, var_p, var_x, var_lapse, var_pi, var_x_raw)

  ELSE  ! Ensure we don't return valid cpt and prh flags if they cannot be calculated

    tph_cpt_qc_flag = ropp_MIFV
    prh_cpt_qc_flag = ropp_MIFV

    DEALLOCATE(lregion)

  END IF ! tph_possible = .TRUE.


!-------------------------------------------------------------------------------
! 3. Copy TPH variables to ROprof structure
!-------------------------------------------------------------------------------

  ro_data%lev2c%tph_tdry_lrt      = tph_lrt     ! Return a refractivity altitude
  ro_data%lev2c%tpt_tdry_lrt      = tpt_lrt     ! Return a dry temperature
  ro_data%lev2c%tph_tdry_lrt_flag = tph_lrt_qc_flag

  ro_data%lev2c%tph_tdry_cpt      = tph_cpt     ! Return a refractivity altitude
  ro_data%lev2c%tpt_tdry_cpt      = tpt_cpt     ! Return a dry temperature
  ro_data%lev2c%tph_tdry_cpt_flag = tph_cpt_qc_flag

  ro_data%lev2c%prh_tdry_cpt      = prh_cpt     ! Return a refractivity altitude
  ro_data%lev2c%prt_tdry_cpt      = prt_cpt     ! Return a dry temperature
  ro_data%lev2c%prh_tdry_cpt_flag = prh_cpt_qc_flag


!-------------------------------------------------------------------------------
! 4. Reset routine name
!-------------------------------------------------------------------------------

  CALL message_set_routine(routine)


CONTAINS


!-------------------------------------------------------------------------------
! 5. Dry pressure estimation
!-------------------------------------------------------------------------------

  SUBROUTINE calc_pdry(refrac, lregion, tdry, zdry, pdry)
!
! Estimate dry pressure from refractivity and dry temperature via
!
! pdry = N * Tdry / kappa1
!
! if possible. If not, estimate pdry from the dry temperature profile, 
! by assuming Tdry varies linearly between levels.  The hydrostatic equation 
! then implies
!
! (pdry(i+1)/pdry(i)) = (Tdry(i+1)/Tdry(i))**(-g/R beta(i))
!
! where
!
! beta(i) = (Tdry(i+1) - Tdry(i)) / (z(i+1) - z(i)).
!
! This method requires an estimate of p(1), which we crudely take to be
!
! pdry(1) = p_ref * EXP(-g*z(1)/R Tdry(1)).
!

! Definitions

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_utils, ONLY: ropp_MDTV, ropp_MDFV

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(IN)     :: refrac
    LOGICAL,  DIMENSION(:), INTENT(IN)     :: lregion
    REAL(wp), DIMENSION(:), INTENT(IN)     :: zdry, tdry
    REAL(wp), DIMENSION(:), INTENT(INOUT)  :: pdry

    REAL(wp)                               :: p1, beta

    INTEGER                                :: n_points, n_points1
    INTEGER                                :: i, i1

! Calculation

    n_points  = SIZE(lregion)

    n_points1 = COUNT(lregion)  ! We know this is >/= min_no_points from earlier check

    i1 = 1
    DO i=1,n_points
      IF (lregion(i)) THEN
        IF (refrac(i) > ropp_MDTV) THEN
          pdry(i1) = refrac(i) * tdry(i1) * 1.0e-2_wp / kappa1  !NB: pdry in hPa
        ELSE
          pdry(i1) = ropp_MDFV
        END IF
        i1 = i1 + 1
      END IF
    END DO

    IF (pdry(1) < ropp_MDTV) THEN
      CALL message(msg_warn, "Estimating level 1 pdry roughly")
      p1      = p_ref * EXP(-g_wmo*zdry(1)/R_dry/tdry(1))  ! Rough estimate, if necessary
      pdry(1) = p1
    END IF

    DO i1=2,n_points1
      IF (pdry(i1) < ropp_MDTV) THEN
        beta = (tdry(i1) - tdry(i1-1)) / (zdry(i1) - zdry(i1-1) + TINY(1.0_wp))
        pdry(i1) = pdry(i1-1) * (tdry(i1)/tdry(i1-1))**(-g_wmo/R_dry/(beta+TINY(1.0_wp)))
      END IF
    END DO

  END SUBROUTINE calc_pdry


END SUBROUTINE ropp_pp_tph_tdry
