! $Id: ropp_1dvar.f90 3551 2013-02-25 09:51:28Z idculv $

!****t* Interface/Modules
!
! SYNOPSIS
!    use ropp_1dvar
!
!    use ropp_fm_types
!    use ropp_fm_constants
!
! DESCRIPTION
!    Access to the routines in the ROPP 1dVar library ropp_1dvar is
!    provided by the single Fortran 90 module ropp_1dvar. The ropp_1dvar module
!    also includes the module ropp_1dvar_types, which contains derived type
!    (structure) declarations used by the 1dVar routines.
!
!    Meteorological and physical constants are provided in the module
!    ropp_1dvar_constants; this module is loaded with the ropp_1dvar module.
!
! SEE ALSO
!    ropp_1dvar
!    ropp_1dvar_types
!    ropp_1dvar_constants
!    ropp_1dvar_copy
!
!****

!****m* Modules/ropp_1dvar *
!
! NAME
!    ropp_1dvar - Interface module for the ROPP 1DVar implementations.
!
! SYNOPSIS
!    use ropp_1dvar
!
! DESCRIPTION
!    This module provides interfaces for all 1DVar functions and routines
!    in the ROPP 1DVar library.
!
! NOTES
!
! SEE ALSO
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

MODULE ropp_1dvar

!-------------------------------------------------------------------------------
! 1. Other modules
!-------------------------------------------------------------------------------

  USE ropp_1dvar_types

!-------------------------------------------------------------------------------
! 2. 1DVar solvers and cost function calculators
!-------------------------------------------------------------------------------

! 2.1 1DVar solver
! ----------------

  INTERFACE ropp_1dvar_solve
     SUBROUTINE ropp_1dvar_solve_bangle(obs, bg, state, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1dBangle), INTENT(inout) :: obs
       TYPE(State1dFM),   INTENT(inout) :: bg
       TYPE(State1dFM),   INTENT(inout) :: state
       TYPE(VarConfig),   INTENT(in)    :: config
       TYPE(VarDiag),     INTENT(inout) :: diag
     END SUBROUTINE ropp_1dvar_solve_bangle
     SUBROUTINE ropp_1dvar_solve_refrac(obs, bg, state, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1dRefrac), INTENT(inout) :: obs
       TYPE(State1dFM),   INTENT(inout) :: bg
       TYPE(State1dFM),   INTENT(inout) :: state
       TYPE(VarConfig),   INTENT(in)    :: config
       TYPE(VarDiag),     INTENT(inout) :: diag
     END SUBROUTINE ropp_1dvar_solve_refrac
  END INTERFACE

! 2.2 1DVar cost functions
! ------------------------

  INTERFACE ropp_1dvar_cost
     SUBROUTINE ropp_1dvar_cost_bangle(yo, bg, control, precon, J, control_ad, &
                                       config, indic)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_1dvar_types
       USE matrix_types
       TYPE(Obs1dBangle),      INTENT(inout) :: yo
       TYPE(State1dFM),        INTENT(inout) :: bg
       TYPE(State1dFM),        INTENT(in)    :: control
       TYPE(matrix_sq),        INTENT(in)    :: precon
       REAL(wp),               INTENT(out)   :: J
       REAL(wp), DIMENSION(:), INTENT(out)   :: control_ad
       TYPE(VarConfig),        INTENT(in)    :: config
       INTEGER,                INTENT(inout) :: indic
     END SUBROUTINE ropp_1dvar_cost_bangle
     SUBROUTINE ropp_1dvar_cost_refrac(yo, bg, control, precon, J, control_ad, &
                                       config, indic)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_1dvar_types
       USE matrix_types
       TYPE(Obs1dRefrac),      INTENT(inout) :: yo
       TYPE(State1dFM),        INTENT(inout) :: bg
       TYPE(State1dFM),        INTENT(in)    :: control
       TYPE(matrix_sq),        INTENT(in)    :: precon
       REAL(wp),               INTENT(out)   :: J
       REAL(wp), DIMENSION(:), INTENT(out)   :: control_ad
       TYPE(VarConfig),        INTENT(in)    :: config
       INTEGER,                INTENT(inout) :: indic
     END SUBROUTINE ropp_1dvar_cost_refrac
  END INTERFACE

! 2.3 Minimiser
! -------------

  INTERFACE ropp_minropp
     SUBROUTINE ropp_1dvar_minropp(x, g, p, dJ, gconv, niter, indic, miter,   &
                                   maxstore)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(inout) :: x
       REAL(wp), DIMENSION(:), INTENT(inout) :: g
       REAL(wp), DIMENSION(:), INTENT(inout) :: p
       REAL(wp),               INTENT(in)    :: dJ
       REAL(wp),               INTENT(in)    :: gconv
       INTEGER,                INTENT(inout) :: niter
       INTEGER,                INTENT(inout) :: indic
       INTEGER,                INTENT(in)    :: miter
       INTEGER,                INTENT(in)    :: maxstore
     END SUBROUTINE ropp_1dvar_minropp
  END INTERFACE

  INTERFACE ropp_1dvar_levmarq
     SUBROUTINE ropp_1dvar_levmarq_bangle(obs, bg, state, config, diag)
        USE ropp_fm_types
        USE ropp_1dvar_types
        TYPE(Obs1dBangle), INTENT(inout) :: obs
        TYPE(State1dFM),   INTENT(inout) :: bg
        TYPE(State1dFM),   INTENT(inout) :: state
        TYPE(VarConfig),   INTENT(in)    :: config
        TYPE(VarDiag),     INTENT(inout) :: diag
      END SUBROUTINE ropp_1dvar_levmarq_bangle
      SUBROUTINE ropp_1dvar_levmarq_refrac(obs, bg, state, config, diag)
        USE ropp_fm_types
        USE ropp_1dvar_types
        TYPE(Obs1dRefrac), INTENT(inout) :: obs
        TYPE(State1dFM),   INTENT(inout) :: bg
        TYPE(State1dFM),   INTENT(inout) :: state
        TYPE(VarConfig),   INTENT(in)    :: config
        TYPE(VarDiag),     INTENT(inout) :: diag
      END SUBROUTINE ropp_1dvar_levmarq_refrac

   END INTERFACE

!-------------------------------------------------------------------------------
! 3. Preconditioning
!-------------------------------------------------------------------------------

  INTERFACE control2state
     SUBROUTINE ropp_control2state(precon, control, state)
       USE typesizes, ONLY: wp => EightByteReal
       USE matrix_types
       TYPE(matrix_sq),        INTENT(in)  :: precon
       REAL(wp), DIMENSION(:), INTENT(in)  :: control
       REAL(wp), DIMENSION(:), INTENT(out) :: state
     END SUBROUTINE ropp_control2state
  END INTERFACE

  INTERFACE control2state_ad
     SUBROUTINE ropp_control2state_ad(precon, control_ad, state_ad)
       USE typesizes, ONLY: wp => EightByteReal
       USE matrix_types
       TYPE(matrix_sq),        INTENT(in)    :: precon
       REAL(wp), DIMENSION(:), INTENT(inout) :: control_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: state_ad
     END SUBROUTINE ropp_control2state_ad
  END INTERFACE

  INTERFACE state2control
     SUBROUTINE ropp_state2control(precon, state, control)
       USE typesizes, ONLY: wp => EightByteReal
       USE matrix_types
       TYPE(matrix_sq),        INTENT(in)  :: precon
       REAL(wp), DIMENSION(:), INTENT(in)  :: state
       REAL(wp), DIMENSION(:), INTENT(out) :: control
     END SUBROUTINE ropp_state2control
  END INTERFACE

!-------------------------------------------------------------------------------
! 4. Error covariances
!-------------------------------------------------------------------------------

  INTERFACE ropp_1dvar_covar
     SUBROUTINE ropp_1dvar_covar_bg(bg, config)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(State1dFM), INTENT(inout) :: bg
       TYPE(VarConfig), INTENT(in)    :: config
     END SUBROUTINE ropp_1dvar_covar_bg
     SUBROUTINE ropp_1dvar_covar_bangle(obs, config)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle), INTENT(inout) :: obs
       TYPE(VarConfig),   INTENT(in)    :: config
     END SUBROUTINE ropp_1dvar_covar_bangle
     SUBROUTINE ropp_1dvar_covar_refrac(obs, config)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac), INTENT(inout) :: obs
       TYPE(VarConfig),   INTENT(in)    :: config
     END SUBROUTINE ropp_1dvar_covar_refrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 5. Diagnostics
!-------------------------------------------------------------------------------

! 5.1 1DVar diagnostics
! ---------------------

  INTERFACE ropp_1dvar_diagnostics
     SUBROUTINE ropp_1dvar_diag_1dbangle(obs, state, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle), INTENT(inout) :: obs
       TYPE(State1DFM),   INTENT(inout) :: state
       TYPE(VarConfig),   INTENT(in)    :: config
       TYPE(VarDiag),     INTENT(inout) :: diag
     END SUBROUTINE ropp_1dvar_diag_1dbangle
     SUBROUTINE ropp_1dvar_diag_1drefrac(obs, state, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac), INTENT(inout) :: obs
       TYPE(State1DFM),   INTENT(inout) :: state
       TYPE(VarConfig),   INTENT(in)    :: config
       TYPE(VarDiag),     INTENT(inout) :: diag
     END SUBROUTINE ropp_1dvar_diag_1drefrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 6. Common utilities
!-------------------------------------------------------------------------------

  INTERFACE ropp_1dvar_version
    FUNCTION ropp_1dvar_version() RESULT (version)
      CHARACTER (LEN=40) :: version
    END FUNCTION ropp_1dvar_version
  END INTERFACE ropp_1dvar_version

END MODULE ropp_1dvar
