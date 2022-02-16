! $Id: ropp_pp_version.f90 4452 2015-01-29 14:42:02Z idculv $

FUNCTION ropp_pp_version() RESULT (version)

!****f* common/ropp_pp_version *
!
! NAME
!   ropp_pp_version
!
! SYNOPSIS
!   Return ROPP_PP version ID string
!
!   USE ropp_pp
!   version = ropp_pp_version()
!
! DESCRIPTION
!   This function returns the (common) version string for the ROPP_PP
!   module. By default, this function should be called by all ROPP_PP
!   tools to display a version ID when the '-v' command-line switch is
!   used.
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

  CHARACTER (LEN=40) :: version

! Edit this string when ROPP_PP module version is updated

  version = "v8.0 31-Dec-2014"

END FUNCTION ropp_pp_version


