! $Id: ropp_qc.f90 4010 2014-01-10 11:07:40Z idculv $

!****m* Modules/ropp_qc *
!
! NAME
!    ropp_qc - Interface module for the ROPP Quality Control implementations.
!
! SYNOPSIS
!    use ropp_qc
! 
! DESCRIPTION
!    This module provides interfaces for all quality control routines contained
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

MODULE ropp_qc

!-------------------------------------------------------------------------------
! 1. Generic quality control
!-------------------------------------------------------------------------------

  INTERFACE ropp_qc_genqc
     SUBROUTINE ropp_qc_genqc_bangle(obs, bg, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle)                     :: obs
       TYPE(State1DFM)                       :: bg
       TYPE(VarConfig)                       :: config
       TYPE(VarDiag)                         :: diag
     END SUBROUTINE ropp_qc_genqc_bangle
     SUBROUTINE ropp_qc_genqc_refrac(obs, bg, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac)                     :: obs
       TYPE(State1DFM)                       :: bg
       TYPE(VarConfig)                       :: config
       TYPE(VarDiag)                         :: diag
     END SUBROUTINE ropp_qc_genqc_refrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 2. O - B differences
!-------------------------------------------------------------------------------

  INTERFACE ropp_qc_OmB
     SUBROUTINE ropp_qc_OmB_bangle(obs, bg, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle)                     :: obs
       TYPE(State1DFM)                       :: bg
       TYPE(VarConfig)                       :: config
       TYPE(VarDiag)                         :: diag
     END SUBROUTINE ropp_qc_OmB_bangle
     SUBROUTINE ropp_qc_OmB_refrac(obs, bg, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac)                     :: obs
       TYPE(State1DFM)                       :: bg
       TYPE(VarConfig)                       :: config
       TYPE(VarDiag)                         :: diag
     END SUBROUTINE ropp_qc_OmB_refrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 3. Background quality control
!-------------------------------------------------------------------------------

  INTERFACE ropp_qc_bgqc
     SUBROUTINE ropp_qc_bgqc_1dbangle(obs, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle)              :: obs
       TYPE(VarConfig)                :: config
       TYPE(VarDiag)                  :: diag
     END SUBROUTINE ropp_qc_bgqc_1dbangle
     SUBROUTINE ropp_qc_bgqc_1drefrac(obs, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac)              :: obs
       TYPE(VarConfig)                :: config
       TYPE(VarDiag)                  :: diag
     END SUBROUTINE ropp_qc_bgqc_1drefrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 4. Probability of gross error
!-------------------------------------------------------------------------------

  INTERFACE ropp_qc_pge
     SUBROUTINE ropp_qc_pge_1dbangle(obs, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle)              :: obs
       TYPE(VarConfig)                :: config
       TYPE(VarDiag)                  :: diag
     END SUBROUTINE ropp_qc_pge_1dbangle
     SUBROUTINE ropp_qc_pge_1drefrac(obs, config, diag)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac)              :: obs
       TYPE(VarConfig)                :: config
       TYPE(VarDiag)                  :: diag
     END SUBROUTINE ropp_qc_pge_1drefrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 5. Valid data height cutoff
!-------------------------------------------------------------------------------

  INTERFACE ropp_qc_cutoff
     SUBROUTINE ropp_qc_cutoff_bangle(obs, config)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DBangle)              :: obs
       TYPE(VarConfig)                :: config
     END SUBROUTINE ropp_qc_cutoff_bangle
     SUBROUTINE ropp_qc_cutoff_refrac(obs, config)
       USE ropp_fm_types
       USE ropp_1dvar_types
       TYPE(Obs1DRefrac)              :: obs
       TYPE(VarConfig)                :: config
     END SUBROUTINE ropp_qc_cutoff_refrac
  END INTERFACE


!-------------------------------------------------------------------------------
! 6. Check ionospheric parameters
!-------------------------------------------------------------------------------

  INTERFACE ropp_qc_ion
    SUBROUTINE ropp_qc_ion_bangle(obs, config, diag)
      USE ropp_1dvar_types
      TYPE(Obs1dBangle), INTENT(inout)  :: obs
      TYPE(VarConfig), INTENT(in)       :: config
      TYPE(VarDiag), INTENT(inout)      :: diag
    END SUBROUTINE ropp_qc_ion_bangle
    SUBROUTINE ropp_qc_ion_bg(bg, obs, config, diag)
      USE ropp_1dvar_types
      TYPE(State1dFM), INTENT(inout)    :: bg
      TYPE(Obs1dBangle), INTENT(in)     :: obs
      TYPE(VarConfig), INTENT(in)       :: config
      TYPE(VarDiag), INTENT(inout)      :: diag
    END SUBROUTINE ropp_qc_ion_bg
  END INTERFACE


END MODULE ropp_qc
