! $Id: ropp_pp_copy.f90 2107 2009-05-21 16:30:01Z frhl $

!****m* Modules/ropp_pp_copy *
!
! NAME
!    ropp_pp_copy - Interface module for the ROPP PP copying function to ROprof
!
! SYNOPSIS
!    use ropp_pp_copy
!
! DESCRIPTION
!    Data type/structure copying functions using ROprof structures used by the
!    ROPP preprocessing module.
!
! SEE ALSO
!    ropp_pp_diag2roprof
!
! AUTHOR
!   M Gorbunov, Russian Academy of Sciences, Russia.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998-2010 Michael Gorbunov <michael.gorbunov@zmaw.de>
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

MODULE ropp_pp_copy

! 1. Add diagnostic information to ROprof structure
! --------------------------------------------------

  INTERFACE ropp_pp_diag2roprof
     SUBROUTINE ropp_pp_diag2roprof(diag, ro_data)
       USE ropp_io_types
       USE ropp_pp_types
       TYPE(ppDiag),     INTENT(in)    :: diag
       TYPE(ROprof),     INTENT(inout) :: ro_data
     END SUBROUTINE ropp_pp_diag2roprof
  END INTERFACE

END MODULE ropp_pp_copy
