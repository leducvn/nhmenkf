! $Id: ropp_io_test.f90 2197 2009-06-23 09:11:17Z idculv $

!****mi* Test/ropp_io_test *
!
! NAME
!    ropp_io_test - Interface module for the ropp_<module> test routines.
!
! SYNOPSIS
!    USE ropp_io_test
!
! DESCRIPTION
!    This module provides interfaces for some of the routines used to 
!    examine the results of the the 'make test' tests. 
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

MODULE ropp_io_test


  INTERFACE

    SUBROUTINE ropp_io_success(pass, test_name, comment)

      LOGICAL,          INTENT(IN)  :: pass
      CHARACTER(LEN=*), INTENT(IN)  :: test_name
      CHARACTER(LEN=*), INTENT(IN)  :: comment

    END SUBROUTINE ropp_io_success

  END INTERFACE


  INTERFACE

    SUBROUTINE ropp_io_fields_compare(prof1, prof2, ndiff, &
               bg, GEOref, Lev1a, Lev1b, Lev2a, Lev2b, Lev2c, Lev2d, &
               L1L2, onedvar, spectra, tdry)

      USE ropp_io_types, ONLY: ROprof

      TYPE(ROprof), INTENT(in)        :: prof1               ! RO profile from file1.nc
      TYPE(ROprof), INTENT(in)        :: prof2               ! RO profile from file2.nc
      INTEGER, INTENT(inout)          :: ndiff               ! Number of differences
      LOGICAL, INTENT(in)             :: bg, GEOref, &       ! The substructures of
                                         Lev1a, Lev1b, &     ! ROprof to be compared
                                         Lev2a, Lev2b, &
                                         Lev2c, Lev2d
      LOGICAL, OPTIONAL, INTENT(in)   :: L1L2
      LOGICAL, OPTIONAL, INTENT(in)   :: onedvar
      LOGICAL, OPTIONAL, INTENT(in)   :: spectra
      LOGICAL, OPTIONAL, INTENT(in)   :: tdry

    END SUBROUTINE ropp_io_fields_compare

  END INTERFACE


END MODULE ropp_io_test
