! $Id: ropp_fm_set_units.f90 4452 2015-01-29 14:42:02Z idculv $

!****s* Units/ropp_fm_set_units *
!
! NAME
!    ropp_fm_set_units - Set internally used units and valid ranges for ROPP
!                        forward models.
!
! SYNOPSIS
!    use ropp_io
!    use ropp_fm
!      ...
!    type(ROprof) :: rodata
!      ...
!    call ropp_fm_set_units(rodata)
!
! DESCRIPTION
!    This subroutine sets the units and valid range attributes within an ROprof
!    data structure to the units used internally in the forward models of the
!    ropp_fm package. For each variable to be defined a call to
!      call ropp_fm_set_units_range(old_units, new_units, range)
!    is made where old_units is the current unit, new_units is the unit assumed
!    in the forward model routines and range is the valid range array for that
!    variable, converted within the routine.
!
! INPUTS
!    rodata  Radio occultation profile data structure
!
! OUTPUT
!    rodata  As above, but with units and valid range attributes modified to
!             reflect standard physical units as used within ropp_fm and
!             ropp_1dvar.
!
! NOTES
!    Only the units and valid ranges for data levels which are used within the
!    forward models provided by ropp_fm are taken care of; thus, units for
!    level 1a (raw amplitude and phase) are left untouched.
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
! 1. ROprof, scalar arguments
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_set_units_roprof_sca(rodata)

! 1.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
! USE ropp_fm_copy,  not_this => ropp_fm_set_units_roprof_sca

  IMPLICIT NONE

  TYPE(ROprof), INTENT(inout) :: rodata

! 1.2 Header
! ----------

  CALL ropp_fm_set_units_range(rodata%units%overall_qual, "percent",     &
                               rodata%range%overall_qual)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%lat, "degrees_north", &
                               rodata%GEOref%range%lat)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%lon, "degrees_east",  &
                               rodata%GEOref%range%lon)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%roc, "metres",        &
                               rodata%GEOref%range%roc)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%r_coc, "kilometres",  &
                               rodata%GEOref%range%r_coc)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%azimuth, "radians",   &
                               rodata%GEOref%range%azimuth)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%undulation, "metres", &
                               rodata%GEOref%range%undulation)

! 1.3 Level 1b
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev1b%units%lat_tp, "degrees_north", &
                               rodata%Lev1b%range%lat_tp)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%lon_tp, "degrees_east",  &
                               rodata%Lev1b%range%lon_tp)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%azimuth_tp, "radians",   &
                               rodata%Lev1b%range%azimuth_tp)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%impact, "metres",        &
                               rodata%Lev1b%range%impact)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%bangle, "radians",       &
                               rodata%Lev1b%range%bangle)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%bangle_sigma, "radians", &
                               rodata%Lev1b%range%bangle_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%bangle_qual, "percent",  &
                               rodata%Lev1b%range%bangle_qual)

! 1.4 Level 2a
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2a%units%alt_refrac, "metres",    &
                               rodata%Lev2a%range%alt_refrac)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%geop_refrac,             &
                               "geopotential metres",                      &
                               rodata%Lev2a%range%geop_refrac)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%refrac, "N-units",       &
                               rodata%Lev2a%range%refrac)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%refrac_sigma, "N-units", &
                               rodata%Lev2a%range%refrac_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%refrac_qual, "percent",  &
                               rodata%Lev2a%range%refrac_qual)

! 1.5 Level 2b
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2b%units%geop, "geopotential metres", &
                               rodata%Lev2b%range%geop)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%geop_sigma,                  &
                               "geopotential metres",                          &
                               rodata%Lev2b%range%geop_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%press, "Pa",                 &
                               rodata%Lev2b%range%press)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%press_sigma, "Pa",           &
                               rodata%Lev2b%range%press_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%temp, "kelvin",              &
                               rodata%Lev2b%range%temp)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%temp_sigma, "kelvin",        &
                               rodata%Lev2b%range%temp_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%shum, "kg/kg",               &
                               rodata%Lev2b%range%shum)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%shum_sigma, "kg/kg",         &
                               rodata%Lev2b%range%shum_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%meteo_qual, "percent",       &
                               rodata%Lev2b%range%meteo_qual)

! 1.6 Level 2c
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%geop_sfc,                  &
                               "geopotential metres",                        &
                               rodata%Lev2c%range%geop_sfc)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%press_sfc, "Pa",           & 
                               rodata%Lev2c%range%press_sfc)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%press_sfc_sigma, "Pa",     &
                               rodata%Lev2c%range%press_sfc_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%press_sfc_qual, "percent", &
                               rodata%Lev2c%range%press_sfc_qual)

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%Ne_max, "m-3",             &
                               rodata%Lev2c%range%Ne_max)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%Ne_max_sigma, "m-3",       &
                               rodata%Lev2c%range%Ne_max_sigma)

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_peak, "m",               &
                               rodata%Lev2c%range%H_peak)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_peak_sigma, "m",         &
                               rodata%Lev2c%range%H_peak_sigma)

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_width, "m",              &
                               rodata%Lev2c%range%H_width)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_width_sigma, "m",        &
                               rodata%Lev2c%range%H_width_sigma)

! 1.7 Level 2d
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2d%units%level_coeff_a, "Pa", &
                               rodata%Lev2d%range%level_coeff_a)

END SUBROUTINE ropp_fm_set_units_roprof_sca


!-------------------------------------------------------------------------------
! 1b. ROprof 2D, scalar arguments
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_set_units_roprof_sca2d(rodata)

! 1.1 Declarations
! ----------------

  USE ropp_utils
  USE ropp_io
  USE ropp_io_types, ONLY: ROprof2d
! USE ropp_fm_copy,  not_this => ropp_fm_set_units_roprof_sca2d

  IMPLICIT NONE

  TYPE(ROprof2d), INTENT(inout) :: rodata

! 1.2 Header
! ----------

  CALL ropp_fm_set_units_range(rodata%units%overall_qual, "percent",     &
                               rodata%range%overall_qual)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%lat, "degrees_north", &
                               rodata%GEOref%range%lat)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%lon, "degrees_east",  &
                               rodata%GEOref%range%lon)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%roc, "metres",        &
                               rodata%GEOref%range%roc)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%r_coc, "kilometres",  &
                               rodata%GEOref%range%r_coc)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%azimuth, "radians",   &
                               rodata%GEOref%range%azimuth)
  CALL ropp_fm_set_units_range(rodata%GEOref%units%undulation, "metres", &
                               rodata%GEOref%range%undulation)

! 1.3 Level 1b
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev1b%units%lat_tp, "degrees_north", &
                               rodata%Lev1b%range%lat_tp)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%lon_tp, "degrees_east",  &
                               rodata%Lev1b%range%lon_tp)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%azimuth_tp, "radians",   &
                               rodata%Lev1b%range%azimuth_tp)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%impact, "metres",        &
                               rodata%Lev1b%range%impact)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%bangle, "radians",       &
                               rodata%Lev1b%range%bangle)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%bangle_sigma, "radians", &
                               rodata%Lev1b%range%bangle_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev1b%units%bangle_qual, "percent",  &
                               rodata%Lev1b%range%bangle_qual)

! 1.4 Level 2a
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2a%units%alt_refrac, "metres",    &
                               rodata%Lev2a%range%alt_refrac)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%geop_refrac,             &
                               "geopotential metres",                      &
                               rodata%Lev2a%range%geop_refrac)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%refrac, "N-units",       &
                               rodata%Lev2a%range%refrac)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%refrac_sigma, "N-units", &
                               rodata%Lev2a%range%refrac_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2a%units%refrac_qual, "percent",  &
                               rodata%Lev2a%range%refrac_qual)

! 1.5 Level 2b
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2b%units%geop, "geopotential metres", &
                               rodata%Lev2b%range%geop)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%geop_sigma,                  &
                               "geopotential metres",                          &
                               rodata%Lev2b%range%geop_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%press, "Pa",                 &
                               rodata%Lev2b%range%press)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%press_sigma, "Pa",           &
                               rodata%Lev2b%range%press_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%temp, "kelvin",              &
                               rodata%Lev2b%range%temp)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%temp_sigma, "kelvin",        &
                               rodata%Lev2b%range%temp_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%shum, "kg/kg",               &
                               rodata%Lev2b%range%shum)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%shum_sigma, "kg/kg",         &
                               rodata%Lev2b%range%shum_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2b%units%meteo_qual, "percent",       &
                               rodata%Lev2b%range%meteo_qual)

! 1.6 Level 2c
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%geop_sfc,                  &
                               "geopotential metres",                        &
                               rodata%Lev2c%range%geop_sfc)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%press_sfc, "Pa",           & 
                               rodata%Lev2c%range%press_sfc)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%press_sfc_sigma, "Pa",     &
                               rodata%Lev2c%range%press_sfc_sigma)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%press_sfc_qual, "percent", &
                               rodata%Lev2c%range%press_sfc_qual)

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%Ne_max, "m-3",             &
                               rodata%Lev2c%range%Ne_max)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%Ne_max_sigma, "m-3",       &
                               rodata%Lev2c%range%Ne_max_sigma)

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_peak, "m",               &
                               rodata%Lev2c%range%H_peak)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_peak_sigma, "m",         &
                               rodata%Lev2c%range%H_peak_sigma)

  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_width, "m",              &
                               rodata%Lev2c%range%H_width)
  CALL ropp_fm_set_units_range(rodata%Lev2c%units%H_width_sigma, "m",        &
                               rodata%Lev2c%range%H_width_sigma)

! 1.7 Level 2d
! ------------

  CALL ropp_fm_set_units_range(rodata%Lev2d%units%level_coeff_a, "Pa", &
                               rodata%Lev2d%range%level_coeff_a)

END SUBROUTINE ropp_fm_set_units_roprof_sca2d

!-------------------------------------------------------------------------------
! 2. ROprof, array arguments
!-------------------------------------------------------------------------------

SUBROUTINE ropp_fm_set_units_roprof_arr(rodata_set)

! 2.1 Declarations
! ----------------

  USE ropp_io
  USE ropp_io_types, ONLY: ROprof
! USE ropp_fm_copy,  not_this => ropp_fm_set_units_roprof_arr

  IMPLICIT NONE

  TYPE(ROprof), DIMENSION(:), INTENT(inout) :: rodata_set

  INTEGER                                   :: i

! 2.2 Loop over all elements
! --------------------------

  DO i = 1, SIZE(rodata_set)
     CALL ropp_fm_set_units(rodata_set(i))
  ENDDO

END SUBROUTINE ropp_fm_set_units_roprof_arr


!------------------------------------------------------------------------------
! 3. Convert range values to correct units if required
!------------------------------------------------------------------------------

SUBROUTINE ropp_fm_set_units_range(units_old, units_new, range)

  USE typesizes, ONLY: wp => EightByteReal
  USE unitconvert

  IMPLICIT NONE

  CHARACTER(len = *), INTENT(inout)     :: units_old
  CHARACTER(len = *), INTENT(in)        :: units_new
  REAL(wp), DIMENSION(2), INTENT(inout) :: range

  CHARACTER(len=20) :: units

  units = TRIM(units_new)

! 3.1 Check if new data units differ from current definition

  IF ( units_new /= "1" .AND. &
       units_old /= "1" .AND. &
       units_new /= units_old ) THEN

! 3.2 Convert valid range values into new data units
    
     CALL ut_convert(range, units_old, range, units)

! 3.3 Update current units with new unit definition

     units_old = units_new
  ENDIF

END SUBROUTINE
