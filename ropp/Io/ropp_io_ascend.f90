! $Id: ropp_io_ascend.f90 4060 2014-02-03 11:13:00Z idculv $
!
!****f* Initialisation/ropp_io_ascend
!
! NAME
!   ropp_io_ascend - ensure all profiles are in ascending altitude order
!
! SYNOPSIS
!   CALL ropp_io_ascend ( ROdata )
!
! INPUTS
!   ROdata  struc  RO data structure (original)
!
! OUTPUTS
!   ROdata  struc  RO data structure (inverted as necessary)
!
! CALLS
!    reverse
!
! USES
!   typsesizes
!   ropp_io
!
! DESCRIPTION
!   This subroutine tests the height-based parameter arrays for the
!   numerical difference between the first and last valid altitude values.
!   If the last is greater than the first, then the profile is assumed
!   to be in ascending order and we just return to the caller. Otherwise,
!   all parameters in the profile are inverted. The procedure is applied
!   to Level 1b, Level 2a, and Level 2b parts of the profiles 
!   separately. As Level 1a data is time-based and Level 2c is surface data
!   these parts are not processed.
!   Note: There is no checking that altitudes are monotonically increasing.
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

SUBROUTINE ropp_io_ascend ( ROdata ) ! (inout)

! 1. Declarations
!================

  USE messages
  USE arrays,        ONLY: reverse
  USE ropp_utils,    ONLY: ropp_MDTV
  USE ropp_io_types, ONLY: ROprof

  IMPLICIT NONE

! Argument list parameters

  TYPE (ROprof)        :: ROdata     ! RO data

! Local variables

  INTEGER              :: nLev     ! No. of levels
  INTEGER              :: i1, I2   ! Indices to first/last valid heights
  CHARACTER(len = 256) :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_ascend')

! 2. Level 1b profiles
!=====================

  nLev = ROdata%Lev1b%Npoints

  i1 = 1
  DO WHILE ( i1 < nLev )
    IF ( ROdata%Lev1b%Impact(i1) > ropp_MDTV ) EXIT
    i1 = i1 + 1
  END DO
  i2 = nLev
  DO WHILE ( i2 > i1 )
    IF ( ROdata%Lev1b%Impact(i2) > ropp_MDTV ) EXIT
    i2 = i2 - 1
  END DO

  IF ( nLev > 1) THEN
    IF (ROdata%Lev1b%Impact(i2) <    &
        ROdata%Lev1b%Impact(i1) ) THEN

      CALL message(msg_diag, &
         "Level1b (bending angle) data found to be in reverse order "      // &
         "(descending height order) \n Reversing data direction. "         // &
         "Reversing data direction." )
      CALL message(msg_diag, &
         "Note 1dVar users must ensure that the error description is " // &
         "appropriate for profiles in increasing height order. \n" // &
         "See User Guide for more details \n")

! 2.1 Latitude, Longitude & Azimuth
!----------------------------------

      ROdata%Lev1b%Lat_tp(1:nLev)     = reverse (ROdata%Lev1b%Lat_tp(1:nLev))

      ROdata%Lev1b%Lon_tp(1:nLev)     = reverse (ROdata%Lev1b%Lon_tp(1:nLev))

      ROdata%Lev1b%Azimuth_tp(1:nLev) = reverse (ROdata%Lev1b%Azimuth_tp(1:nLev))

! 2.2 L1 BA
!----------

      ROdata%Lev1b%Impact_L1(1:nLev)  = reverse (ROdata%Lev1b%Impact_L1(1:nLev))

      ROdata%Lev1b%Bangle_L1(1:nLev)  = reverse (ROdata%Lev1b%Bangle_L1(1:nLev))

      ROdata%Lev1b%Bangle_L1_sigma(1:nLev) = reverse (ROdata%Lev1b%Bangle_L1_sigma(1:nLev))

      ROdata%Lev1b%Bangle_L1_qual(1:nLev)  = reverse (ROdata%Lev1b%Bangle_L1_qual(1:nLev))

! 2.3 L2 BA
!----------

      ROdata%Lev1b%Impact_L2(1:nLev)  = reverse (ROdata%Lev1b%Impact_L2(1:nLev))

      ROdata%Lev1b%Bangle_L2(1:nLev)  = reverse (ROdata%Lev1b%Bangle_L2(1:nLev))

      ROdata%Lev1b%Bangle_L2_sigma(1:nLev) = reverse (ROdata%Lev1b%Bangle_L2_sigma(1:nLev))

      ROdata%Lev1b%Bangle_L2_qual(1:nLev)  = reverse (ROdata%Lev1b%Bangle_L2_qual(1:nLev))

! 2.4 Corrected BA
!-----------------

      ROdata%Lev1b%Impact(1:nLev)       = reverse (ROdata%Lev1b%Impact(1:nLev))

      ROdata%Lev1b%Bangle(1:nLev)       = reverse (ROdata%Lev1b%Bangle(1:nLev))

      ROdata%Lev1b%Bangle_sigma(1:nLev) = reverse (ROdata%Lev1b%Bangle_sigma(1:nLev))

      ROdata%Lev1b%Bangle_qual(1:nLev)  = reverse (ROdata%Lev1b%Bangle_qual(1:nLev))

! 2.5 Optimised BA
!-----------------

      ROdata%Lev1b%Impact_Opt(1:nLev)       = reverse ( &
                                                ROdata%Lev1b%Impact_Opt(1:nLev) )

      ROdata%Lev1b%Bangle_Opt(1:nLev)       = reverse ( &
                                                ROdata%Lev1b%Bangle_Opt(1:nLev) )

      ROdata%Lev1b%Bangle_Opt_sigma(1:nLev) = reverse ( &
                                          ROdata%Lev1b%Bangle_Opt_sigma(1:nLev) )


      ROdata%Lev1b%Bangle_Opt_qual(1:nLev)  = reverse ( &
                                           ROdata%Lev1b%Bangle_Opt_qual(1:nLev) )
    END IF
  END IF

! 3. Level 2a profiles
!=====================

  nLev = ROdata%Lev2a%Npoints

  i1 = 1
  DO WHILE ( i1 < nLev )
    IF ( ROdata%Lev2a%alt_refrac(i1) > ropp_MDTV ) EXIT
    i1 = i1 + 1
  END DO
  i2 = nLev
  DO WHILE ( i2 > i1 )
    IF ( ROdata%Lev2a%alt_refrac(i2) > ropp_MDTV ) EXIT
    i2 = i2 - 1
  END DO

  IF ( nLev > 1) THEN
    IF (ROdata%Lev2a%alt_refrac(i2) < &
        ROdata%Lev2a%alt_refrac(i1) ) THEN

      CALL message(msg_diag, &
         "Level2a (refractivity) data found to be in reverse order "        // &
         "(descending height order) \n Reversing data direction. ")
      CALL message(msg_diag, &
         "Note 1dVar users must ensure that the error description is " // &
         "appropriate for profiles in increasing height order. \n" // &
         "See User Guide for more details \n")

! 3.1 Altitude
!-------------

      ROdata%Lev2a%alt_refrac(1:nLev)  = reverse (ROdata%Lev2a%alt_refrac(1:nLev))

! 3.2 Geopotential
!-----------------

      ROdata%Lev2a%geop_refrac(1:nLev) = reverse (ROdata%Lev2a%geop_refrac(1:nLev))

! 3.3 Refractivity
!-----------------

      ROdata%Lev2a%Refrac(1:nLev)       = reverse (ROdata%Lev2a%Refrac(1:nLev))

      ROdata%Lev2a%Refrac_sigma(1:nLev) = reverse (ROdata%Lev2a%Refrac_sigma(1:nLev))

      ROdata%Lev2a%Refrac_qual(1:nLev) = reverse (ROdata%Lev2a%Refrac_qual(1:nLev))

! 3.4 Dry temperature
!--------------------

      ROdata%Lev2a%Dry_Temp(1:nLev)      = reverse (ROdata%Lev2a%Dry_Temp(1:nLev))

      ROdata%Lev2a%Dry_Temp_sigma(1:nLev) = reverse (ROdata%Lev2a%Dry_Temp_sigma(1:nLev))

      ROdata%Lev2a%Dry_Temp_qual(1:nLev) = reverse (ROdata%Lev2a%Dry_Temp_qual(1:nLev))

    END IF
  END IF

! 4. Level 2b profiles
!=====================

  nLev = ROdata%Lev2b%Npoints

  i1 = 1
  DO WHILE ( i1 < nLev )
    IF ( ROdata%Lev2b%geop(i1) > ropp_MDTV ) EXIT
    i1 = i1 + 1
  END DO
  i2 = nLev
  DO WHILE ( i2 > i1 )
    IF ( ROdata%Lev2b%geop(i2) > ropp_MDTV ) EXIT
    i2 = i2 - 1
  END DO

  IF ( nLev > 1) THEN
    IF (ROdata%Lev2b%geop(i2) < &
         ROdata%Lev2b%geop(i1) ) THEN

      CALL message(msg_diag, &
         "Level2b (geop,temp,shum,pres) data found to be in reverse order " // &
         "(descending height order) \n Reversing data direction.")      
      CALL message(msg_diag, &
         "Note 1dVar users must ensure that the error description is " // &
         "appropriate for profiles in increasing height order. \n" // &
         "See User Guide for more details \n")

      ROdata%Lev2b%geop(1:nLev)       = reverse (ROdata%Lev2b%geop(1:nLev))

! 4.1 Temperature
!-----------------

      ROdata%Lev2b%Temp(1:nLev)       = reverse (ROdata%Lev2b%Temp(1:nLev))

      ROdata%Lev2b%Temp_sigma(1:nLev) = reverse (ROdata%Lev2b%Temp_sigma(1:nLev))

! 4.2 Specific humidity
!-----------------------

      ROdata%Lev2b%SHum(1:nLev)       = reverse (ROdata%Lev2b%SHum(1:nLev))

      ROdata%Lev2b%SHum_sigma(1:nLev) = reverse (ROdata%Lev2b%SHum_sigma(1:nLev))

! 4.3 Pressure
!--------------

      ROdata%Lev2b%Press(1:nLev)       = reverse (ROdata%Lev2b%Press(1:nLev))

      ROdata%Lev2b%Press_sigma(1:nLev) = reverse (ROdata%Lev2b%Press_sigma(1:nLev))

! 4.4 Meteo quality
!------------------

      ROdata%Lev2b%Meteo_qual(1:nLev)  = reverse (ROdata%Lev2b%Meteo_qual(1:nLev))


! 5. Level 2d profiles - invert if Level 2b data is inverted
!=====================

      nLev = ROdata%Lev2d%Npoints
    
      IF (nLev > 1) THEN

        ROdata%Lev2d%Level_coeff_a(1:nLev) = reverse (ROdata%Lev2d%Level_coeff_a(1:nLev))
      
        ROdata%Lev2d%Level_coeff_b(1:nLev) = reverse (ROdata%Lev2d%Level_coeff_b(1:nLev))
      
      ENDIF


    END IF
  
  END IF

  CALL message_set_routine(routine)

END SUBROUTINE ropp_io_ascend
