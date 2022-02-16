! $Id: ropp_io_isinrange.f90 3551 2013-02-25 09:51:28Z idculv $

!****f* Datatypes/isroppinrange *
!
! NAME
!    isroppinrange - Check elements of a data structure to be in a valid
!                     range.
!
! SYNOPSIS
!    true_or_not = isroppinrange(<data>)
!
! DESCRIPTION
!    This function checks if a given data item is within a valid range; it
!    compliments the isinrange() function from the ropp_utils package to
!    add data types used within ropp_io. As the latter carry their own
!    valid_range information as part of their structure, no valid range needs
!    to be specified.
!
! INPUTS
!    <data>       Data to be checked; currently limited to type(DT7type).
!
! OUTPUT
!    true_or_not  .true. if <data> is within the valid_range;
!                 .false. otherwise.
!
! CALLS
!    isinrange   from ropp_utils library
!
! NOTES
!    In the current version of ropp_io, isroppinrange only provides
!    range checking for the DT7 data type, i.e. for dates and times.
!
! SEE ALSO
!    DT7type
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

FUNCTION ropp_io_isinrangeDT7(dt) RESULT(inrange)

  USE ropp_utils,    ONLY: isinrange
! USE ropp_io,       not_this => ropp_io_isinrangeDT7
  USE ropp_io_types, ONLY: DT7type

  IMPLICIT NONE

  TYPE(DT7type), INTENT(in) :: dt
  LOGICAL                   :: inrange

!-------------------------------------------------------------------------------
! 3. Range checks
!-------------------------------------------------------------------------------

  inrange = isinrange(dt%year,   dt%range%year)   .AND. &
            isinrange(dt%month,  dt%range%month)  .AND. &
            isinrange(dt%day,    dt%range%day)    .AND. &
            isinrange(dt%hour,   dt%range%hour)   .AND. &
            isinrange(dt%minute, dt%range%minute) .AND. &
            isinrange(dt%second, dt%range%second) .AND. &
            isinrange(dt%msec,   dt%range%msec)

END FUNCTION ropp_io_isinrangeDT7
