! $Id: ropp_1dvar_copy.f90 3551 2013-02-25 09:51:28Z idculv $

!****m* Modules/ropp_1dvar_copy *
!
! NAME
!    ropp_1dvar_copy - Interface module for the ROPP 1dVar
!
! SYNOPSIS
!    use ropp_1dVar_copy
!
! DESCRIPTION
!    Data type/structure copying functions using ROprof structures used by the
!    forward models of ROPP and converting units within the ROprof structure.
!
! SEE ALSO
!    ropp_1dvar_diag2roprof
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

MODULE ropp_1dvar_copy

! 1. Add diagnostic information to ROprof structure
! --------------------------------------------------

  INTERFACE ropp_1dvar_diag2roprof
     SUBROUTINE ropp_1dvar_diag2roprof_bangle(obs, diag, ro_data, config)
       USE ropp_io_types
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle), INTENT(in)    :: obs
       TYPE(VarDiag),     INTENT(in)    :: diag
       TYPE(ROprof),      INTENT(inout) :: ro_data
       TYPE(VarConfig),   INTENT(in)    :: config
     END SUBROUTINE ropp_1dvar_diag2roprof_bangle
     SUBROUTINE ropp_1dvar_diag2roprof_refrac(obs, diag, ro_data, config)
       USE ropp_io_types
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac), INTENT(in)    :: obs
       TYPE(VarDiag),     INTENT(in)    :: diag
       TYPE(ROprof),      INTENT(inout) :: ro_data
       TYPE(VarConfig),   INTENT(in)    :: config
     END SUBROUTINE ropp_1dvar_diag2roprof_refrac
  END INTERFACE

END MODULE ropp_1dvar_copy
