! $Id: ropp_io_thin_skip.f90 3551 2013-02-25 09:51:28Z idculv $
!
!****s* Thin/ropp_io_thin_skip *
!
! NAME
!   ropp_io_thin_skip - skip profile levels to simply sub-sample data
!
! ARGUMENTS
!   nLev     (in)   int  No. of levels in full profile
!   nThinLev (in)   int  Max. no. of thinned levels required
!   skip1    (out)  int  Skip factor to select levels
!   skip2    (out)  int  Skip factor to reject levels
!   nSamp    (out)  int  Actual no. of thinned samples
!
! SYNPOSIS
!   call ropp_io_thin_skip ( ROdata%Lev1b%Npoints, nThinLev, &
!                            skip1, skip2, nSamp, DEBUG )
!   do i = 1, ROdata%Lev1b%Npoints, skip1 ! select every skip1-th sample
!     if ( mod(in,skip2) == 0 ) cycle     ! reject every skip2-th sample
!     ....
!   enddo
!
! CALLED BY
!   ropp_io_thin_select
!
! DESCRIPTION
!  Given the actual number of levels in a profile, and the
!  maximum desired number, this subroutine calculates skip
!  factors allowing appropriate thinning by
!    - selecting 1-in-N levels when levels > maxlevels * 2
!    - rejecting 1-in-N levels when max levels < levels < maxlevels * 2
!    - no thinning when levels < maxlevels
!  The routine can cater for thinning in cases where the input is
!  high resolution (raw sampling) or low resolution (pre-thinned) but
!  still too many levels for BUFR (either in terms of message length
!  or user-application array size limits - e.g. MetDB is limited to 375
!  levels).
!  Note that this approach effectively sub-samples the original profile
!  and there is no fixed no. of thinned levels (just a maximum).
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

SUBROUTINE ropp_io_thin_skip ( nLev,     & ! (in)
                               nThinLev, & ! (in)
                               skip1,    & ! (out)
                               skip2,    & ! (out)
                               nSamp )     ! (out)

  USE messages

! Argument list paremeters

  INTEGER, INTENT(IN)           :: nLev          ! no. of full samples
  INTEGER, INTENT(IN)           :: nThinLev      ! max. no. thinned samples
  INTEGER, INTENT(OUT)          :: skip1, skip2  ! skip factors
  INTEGER, INTENT(OUT)          :: nSamp         ! no. thinned samples
  CHARACTER(len = 256)          :: routine

! Local variables

  CHARACTER (LEN=10) :: strnum1, strnum2, strnum3        ! string numeric values

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_thin_skip')

! 1. Thinning by a factor of 2 or more (select 1-in-N)
!=====================================================

  IF ( nLev > nThinLev*2 ) THEN
    skip1 = 1 + ( nLev - 1 ) / nThinLev
    skip2 = 1 +   nLev
    nSamp = 1 + ( nLev - 1 ) / skip1
    IF ( skip1 > 1) THEN
      WRITE ( strnum1, "(I10)" ) skip1
      WRITE ( strnum2, "(I10)" ) nSamp
      WRITE ( strnum3, "(I10)" ) nLev
      CALL message(msg_diag, " Selecting by a factor of : 1-in-" // &
                   TRIM(ADJUSTL(strnum1))//" to "//           &
                   TRIM(ADJUSTL(strnum2))//" samples from "// &
                   TRIM(ADJUSTL(strnum3)))
    END IF

! 2. Thinning by a factor between 1 & 2 (reject 1-in-N)
!======================================================

  ELSE IF ( nLev > nThinLev ) THEN
    skip1 = 1
    skip2 = MAX ( nLev / ( nLev - nThinLev ), 2 )
    nSamp = nLev - nLev / skip2
    IF ( skip2 > 1 ) THEN
      WRITE ( strnum1, "(I10)" ) skip2
      WRITE ( strnum2, "(I10)" ) nSamp
      WRITE ( strnum3, "(I10)" ) nLev
      CALL message(msg_diag, " Rejecting by a factor of : 1-in-" // &
                   TRIM(ADJUSTL(strnum1))//" to "//           &
                   TRIM(ADJUSTL(strnum2))//" samples from "// &
                   TRIM(ADJUSTL(strnum3)))
    END IF

! 3. No thinning
!===============

  ELSE
    skip1 = 1
    skip2 = nLev + 1
    nSamp = nLev

  ENDIF

  CALL message_set_routine(routine)

END SUBROUTINE ropp_io_thin_skip
