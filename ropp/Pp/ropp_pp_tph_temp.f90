! $Id: ropp_pp_tph_temp.f90 3491 2013-02-06 12:43:43Z idculv $

SUBROUTINE ropp_pp_tph_temp(ro_data, diag)

!****s* TropoPauseHeight/ropp_pp_tph_temp *
!
! NAME
!   ropp_pp_tph_temp
!
! SYNOPSIS
!   Tropopause height diagnostic based on temperature
!
!   CALL ropp_pp_tph_temp(ro_data, diag)
!
! DESCRIPTION
!   Diagnose tropopause height from the kinks in the temperature profile,
!   using the lapse rate and cold point criteria (Reichler et al, GRL, 2003).
!
! INPUTS
!   TYPE(ROprof), INTENT(INOUT)    :: ro_data    ! input RO profile containing lev2b data
!   LOGICAL, OPTIONAL, INTENT(IN)  :: diag       ! extra diagnostics required
!
!
! OUTPUTS
!   TYPE(ROprof), INTENT(INOUT)    :: ro_data    ! output RO profile containing lev2b data
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
! USE ropp_pp, not_this => ropp_pp_tph_temp
  USE ropp_pp
  USE ropp_pp_constants, ONLY: R_dry, C_p, g_wmo

  IMPLICIT NONE

  TYPE(ROprof), INTENT(INOUT)          :: ro_data    ! input RO profile containing lev2b data
  LOGICAL, OPTIONAL, INTENT(IN)        :: diag       ! extra diagnostics required

! Parameters in the temperature-based TPH diagnosis algorithms

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

  REAL(wp), DIMENSION(:), ALLOCATABLE  :: var_x,var_y,var_p,&
                                          var_lapse,var_pi,&
                                          var_y1,var_p1

  REAL(wp), DIMENSION(:), ALLOCATABLE  :: geop

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
! 2. Calculate TPH based on lapse rate (wrt pressure) and cold point of temperature
!-------------------------------------------------------------------------------

! 2.0 Initial messages
! --------------------

  CALL message_get_routine(routine)

  CALL message_set_routine('ropp_pp_tph_temp')

  CALL message(msg_info, "Calculating temperature-based tropopause height \n")

! 2.1 Initialise variables
! ------------------------

  n_points = ro_data%lev2b%npoints

  IF ((n_points == 0) .OR. ro_data%lev2b%missing) THEN
    CALL message(msg_info, "No (valid) Level 2b data in profile ... " //    &
                           "will not calculate TPH based on temperature \n")
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

  ALLOCATE(geop(n_points))  ;  geop = ro_data%lev2b%geop

! 2.2 Check for numerically invalid data
! --------------------------------------

  IF (ALL((ro_data%lev2b%temp  < ropp_MDTV) .OR. &
          (ro_data%lev2b%press < ropp_MDTV) .OR. &
          (geop                < ropp_ZERO))) THEN

    CALL message(msg_diag, "No common non-missing temperatures, " // &
                           "pressures and positive geopotentials " // &
                           "... trying to generate geopotentials \n")

    CALL calc_geop(ro_data, geop)

  END IF

  IF (ALL((ro_data%lev2b%temp  < ropp_MDTV) .OR. &
          (ro_data%lev2b%press < ropp_MDTV) .OR. &
          (geop                < ropp_ZERO))) THEN
    CALL message(msg_error, "Still no common non-missing temperatures, " // &
                            "pressures and positive geopotentials " // &
                            "... cannot calculate TPH \n")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_data_invalid)
  END IF

  IF (ANY((ro_data%lev2b%press <= ropp_ZERO) .AND. &
          (ro_data%lev2b%press >= ropp_MDTV))) THEN
    CALL message(msg_error, "Non-positive pressures encountered ... " // &
                            "... cannot calculate TPH \n")
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

  WHERE ((ro_data%lev2b%temp  > ropp_MDTV) .AND. &
         (ro_data%lev2b%press > ropp_MDTV) .AND. &
         (              geop  > ropp_ZERO)) lregion = .TRUE.
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

  IF (MINVAL(geop, mask=lregion) > tph_min) THEN
    CALL message(msg_info, "Geopotential heights do not start deep enough")
    tph_possible = .FALSE.
    tph_lrt_qc_flag  = IBSET(tph_lrt_qc_flag, TPH_QC_prof_depth)
  END IF

  IF (MAXVAL(geop, mask=lregion) < tph_max) THEN
    CALL message(msg_info, "Geopotential heights do not reach high enough")
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

    ALLOCATE(var_y(n_points1), var_p(n_points1), var_x(n_points1))  ! Fill these with the non-missing data

    i1 = 1
    DO i=1,n_points
      IF (lregion(i)) THEN
        var_y(i1) = geop(i)
        var_p(i1) = ro_data%lev2b%press(i)
        var_x(i1) = ro_data%lev2b%temp(i)
        i1        = i1 + 1
      END IF
    END DO

    DEALLOCATE(lregion)

! Below this point, the length of all arrays should be n_points1

    ALLOCATE(lregion(n_points1), lregion1(n_points1))  ! Masks

! 2.5 Calculate Exner pressure pi = (p/p_ref)**kappa
! --------------------------------------------------

    ALLOCATE(var_pi(n_points1))  !  Will hold Exner pressure Pi

    var_pi = (var_p / p_ref) ** kappa

! 2.6 Compute lapse rate at half points (i+1/2), assuming T varies linearly with pi
! ---------------------------------------------------------------------------------

    ALLOCATE(var_lapse(n_points1)) !  Will hold dT/dz (ie, -lapse rate)

    var_lapse = (CSHIFT(var_x,  1) - var_x) * &
                (CSHIFT(var_pi, 1) + var_pi) / &
                (CSHIFT(var_pi, 1) - var_pi + TINY(1.0_wp)) / &
                (CSHIFT(var_x,  1) + var_x  + TINY(1.0_wp))

    var_lapse(n_points1) = var_lapse(n_points1-1)

    var_lapse = -1.0_wp * var_lapse * (g_wmo/C_p)   ! var_lapse(i) = +dT/dz(i+1/2)

    IF (l_diag) THEN

! Append lapse rate, pressure and geopotential height to the ROprof structure and thence the output file.

      CALL ropp_io_addvar_rodataD1d( &
           ro_data, &
           name      = "LR_temp", &
           long_name = "Temperature lapse rate", &
           units     = "K/km", &
           range     = (/-1000000.0_wp, 1000000.0_wp/), &
           DATA      = var_lapse)

      ALLOCATE(var_p1(n_points1))  ! This is where lapse rate is defined (approximately)
      var_p1 = (var_p + CSHIFT(var_p, 1)) / 2.0_wp
      var_p1(n_points1) = var_p(n_points1)
      CALL ropp_io_addvar_rodataD1d( &
           ro_data, &
           name      = "LR_temp_press", &
           long_name = "Temperature lapse rate pressure", &
           units     = "hPa", &
           range     = (/-1000000.0_wp, 1000000.0_wp/), &
           DATA      = var_p1)
      DEALLOCATE(var_p1)

      ALLOCATE(var_y1(n_points1))  ! This is where lapse rate is defined (approximately)
      var_y1 = (var_y + CSHIFT(var_y, 1)) / 2.0_wp
      var_y1(n_points1) = var_y(n_points1)
      CALL ropp_io_addvar_rodataD1d( &
           ro_data, &
           name      = "LR_temp_geop", &
           long_name = "Temperature lapse rate geopotential", &
           units     = "m", &
           range     = (/-1000000.0_wp, 1000000.0_wp/), &
           DATA      = var_y1)
      DEALLOCATE(var_y1)

    END IF

! 2.7 Calculate TPH according to Reichler's method
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

! 2.8 Calculate precise location of TPH by interpolation
! ------------------------------------------------------
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
! so we interpolate z linearly with ln(p), as for an isothermal atmosphere.

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

! 2.9 Calculate cold point tropopause if in tropics
! -------------------------------------------------

! Note that we need tph_lrt in this algorithm, so we cannot bypass the earlier calculation

    IF ((.NOT. BTEST(tph_lrt_qc_flag, TPH_QC_data_invalid)) .AND. &
        (.NOT. BTEST(tph_cpt_qc_flag, TPH_QC_data_invalid))) THEN

      IF (ABS(lat) <= tropic_lat_bdy) THEN ! In tropics

        lregion = .FALSE.
        WHERE ((var_y > tph_min) .AND. (var_y < tph_max)) lregion = .TRUE.

        IF (ANY(lregion)) THEN

          tpt_cpt   = MINVAL(var_x, mask=lregion)
          cpt_index = SUM(MINLOC(ABS(var_x-tpt_cpt), mask=lregion))
          tph_cpt   = var_y(cpt_index)

          IF (ABS(tph_cpt - tph_lrt) > tph_diff_max) THEN

            lregion1 = .FALSE.
            WHERE (ABS(var_y - tph_lrt) < tph_diff_max) lregion1 = .TRUE.

            IF (ANY(lregion1)) THEN

              tpt_cpt   = MINVAL(var_x, mask=lregion1)
              cpt_index = SUM(MINLOC(ABS(var_x-tpt_cpt), mask=lregion1))
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

! 2.10 Check that TPH is not too low
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

! 2.11 Check that TPH is not too high
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

! 2.12 Calculate and store profile minimum temp and its height
! ------------------------------------------------------------

    IF (.NOT. BTEST(prh_cpt_qc_flag, TPH_QC_data_invalid)) THEN

      lregion = .FALSE.
      WHERE (var_x > ropp_MDTV) lregion = .TRUE.

      IF (ANY(lregion)) THEN

        prt_cpt   = MINVAL(var_x, mask=lregion)
        prf_index = SUM(MINLOC(ABS(var_x-prt_cpt), mask=lregion))
        prh_cpt   = var_y(prf_index)

      ELSE

        CALL message(msg_info, "Could not calculate entire profile cold point ... " // &
                               "leaving as missing data \n")
        prh_cpt_qc_flag = IBSET(prh_cpt_qc_flag, TPH_QC_data_invalid)

      END IF

    END IF

! 2.13 Clean up
! -------------

    DEALLOCATE(lregion, lregion1, var_y, var_p, var_x, var_lapse, var_pi)

  ELSE  ! Ensure we don't return valid cpt and prh flags if they cannot be calculated

    tph_cpt_qc_flag = ropp_MIFV
    prh_cpt_qc_flag = ropp_MIFV

    DEALLOCATE(lregion)

  END IF ! tph_possible = .TRUE.

  DEALLOCATE(geop)


!-------------------------------------------------------------------------------
! 3. Copy TPH variables to ROprof structure
!-------------------------------------------------------------------------------

  ro_data%lev2c%tph_temp_lrt      = tph_lrt     ! Return a geopotential height
  ro_data%lev2c%tpt_temp_lrt      = tpt_lrt     ! Return a temperature
  ro_data%lev2c%tph_temp_lrt_flag = tph_lrt_qc_flag

  ro_data%lev2c%tph_temp_cpt      = tph_cpt     ! Return a geopotential height
  ro_data%lev2c%tpt_temp_cpt      = tpt_cpt     ! Return a temperature
  ro_data%lev2c%tph_temp_cpt_flag = tph_cpt_qc_flag

  ro_data%lev2c%prh_temp_cpt      = prh_cpt     ! Return a geopotential height
  ro_data%lev2c%prt_temp_cpt      = prt_cpt     ! Return a temperature
  ro_data%lev2c%prh_temp_cpt_flag = prh_cpt_qc_flag


!-------------------------------------------------------------------------------
! 4. Reset routine name
!-------------------------------------------------------------------------------

  CALL message_set_routine(routine)


CONTAINS


!-------------------------------------------------------------------------------
! 5. Geopotential calculation
!-------------------------------------------------------------------------------

  SUBROUTINE calc_geop(ro_data, geop)
!
! Calculate geopotential for ECMWF-like bgr files from {T, q, Ak, Bk}.
! Code lifted from relevant section of ropp_fm_roprof2state.f90,
! to which the reader is referred for explanatory comments.

! Definitions

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_utils
    USE ropp_io_types, ONLY: ROprof
    USE ropp_pp_constants, ONLY: R_dry, g_wmo

    IMPLICIT NONE

    TYPE(ROprof), INTENT(IN)              :: ro_data ! RO profile containing lev2b data
    REAL(wp), DIMENSION(:), INTENT(INOUT) :: geop

    REAL(wp), DIMENSION(:), ALLOCATABLE  :: p_hlv, geop_hlv
    REAL(wp), DIMENSION(:), ALLOCATABLE  :: Tvflv, akk, bkk, ln_prflv
    REAL(wp), DIMENSION(:), ALLOCATABLE  :: pres, del_p, del_geop, alpha
    REAL(wp)                             :: psfc

    INTEGER                              :: n_hlv, n_flv, lvl

    CHARACTER(len=256)                   :: level_type

! Calculation

    level_type = ro_data%Lev2d%level_type
    IF ( INDEX(level_type,'UNKNOWN') > 0) level_type = ro_data%bg%source

    CALL To_Upper(level_type)

    IF ((INDEX(level_type,'HYBRID') > 0) .OR. &
        (INDEX(level_type,'ECMWF')  > 0)) THEN  ! ECMWF background

      n_hlv = SIZE(ro_data%lev2d%level_coeff_a)
      n_flv = SIZE(ro_data%lev2b%press)

      ALLOCATE(p_hlv(n_hlv), geop_hlv(n_hlv))
      ALLOCATE(Tvflv(n_flv), akk(n_flv), bkk(n_flv), ln_prflv(n_flv))
      ALLOCATE(pres(n_flv), del_p(n_flv), del_geop(n_flv), alpha(n_flv))

      psfc = ro_data%lev2c%press_sfc
      p_hlv = ro_data%lev2d%level_coeff_a + ro_data%lev2d%level_coeff_b * psfc

      akk = 0.5_wp * (ro_data%lev2d%level_coeff_a(1:n_hlv-1) + &
                      ro_data%lev2d%level_coeff_a(2:n_hlv))

      bkk = 0.5_wp * (ro_data%lev2d%level_coeff_b(1:n_hlv-1) + &
                      ro_data%lev2d%level_coeff_b(2:n_hlv))

      pres = akk + bkk * psfc

      WHERE (ro_data%lev2b%shum > ropp_MDTV) ! very occasionally they're not
        Tvflv  = (1.0_wp + 0.61e-3_wp * ro_data%lev2b%shum) * ro_data%lev2b%temp !NB: shum in g/kg
      ELSEWHERE
        Tvflv  = (1.0_wp + 0.61e-3_wp *             0.0_wp) * ro_data%lev2b%temp !NB: shum in g/kg
      ENDWHERE

      del_p = p_hlv(1:n_hlv-1) - p_hlv(2:n_hlv)

      ln_prflv = LOG(p_hlv(1:n_hlv-1)/p_hlv(2:n_hlv))

      alpha    = 1.0_wp - p_hlv(2:n_hlv)/del_p * ln_prflv
      alpha(n_hlv-1) = LOG(2.0_wp)

      del_geop =   R_dry * Tvflv * ln_prflv / g_wmo

      DO lvl = 1, n_hlv
        geop_hlv(lvl) = ro_data%lev2c%geop_sfc + SUM(del_geop(1:lvl-1)) !! N.B:sum(X(1:0))=0
      END DO

      geop = geop_hlv(1:n_hlv-1) + alpha * R_dry * Tvflv / g_wmo

      DEALLOCATE(p_hlv, geop_hlv)
      DEALLOCATE(Tvflv, akk, bkk, ln_prflv)
      DEALLOCATE(pres, del_p, del_geop, alpha)

    ELSE

      CALL message(msg_warn, "Background profile not in ECMWF format ... " // &
                             "cannot not calculate geopotential \n")
      geop = ropp_MDFV

    END IF

  END SUBROUTINE calc_geop
 
END SUBROUTINE ropp_pp_tph_temp
