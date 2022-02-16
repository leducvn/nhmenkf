! $Id: ropp_pp_preproc.f90 1961 2008-11-13 16:41:31Z frhl $

!****m* Modules/ropp_pp_preproc *
!
! NAME
!    ropp_pp_preproc - Interface module for the ROPP pre-processing module.
!
! SYNOPSIS
!    use ropp_pp_preproc
! 
! DESCRIPTION
!    Data type/structure copying functions using ROprof structures used by the
!    ropp_pp module.
!
! SEE ALSO
!    ropp_pp_set_coordinates
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

MODULE ropp_pp_preproc

!-------------------------------------------------------------------------------
! 1. Standard coordinate frames for radio occultation processing
!-------------------------------------------------------------------------------

!****t* Tools/Coordinates
!
! DESCRIPTION
!    Routines for dealing with coordinate frames
!
! SEE ALSO
!    ropp_pp_set_coordinates
!
!****

  INTERFACE ropp_pp_set_coordinates
     SUBROUTINE ropp_pp_set_coordinates_sca(ro_data)
       USE ropp_io_types
       TYPE(ROprof), INTENT(inout) :: ro_data
     END SUBROUTINE ropp_pp_set_coordinates_sca
     SUBROUTINE ropp_pp_set_coordinates_arr(ro_data_set)
       USE ropp_io_types
       TYPE(ROprof), DIMENSION(:), INTENT(inout) :: ro_data_set
     END SUBROUTINE ropp_pp_set_coordinates_arr
  END INTERFACE

!-------------------------------------------------------------------------------
! 2. Mission-specific pre-processing stages to update ROprof strucutre
!-------------------------------------------------------------------------------

!****t* Preprocessing/Preprocess
!
! DESCRIPTION
!    Routines for mission-specific pre-processing steps
!
! SEE ALSO
!    ropp_pp_preprocess
!
!****

  INTERFACE
     SUBROUTINE ropp_pp_preprocess(ro_data, config, diag)
       USE ropp_io_types
       USE ropp_pp_types
       TYPE(ROprof),   INTENT(inout) :: ro_data   ! Radio occultation data
       TYPE(PPConfig), INTENT(inout) :: config    ! Configuration options
       TYPE(PPDiag),   INTENT(inout) :: diag      ! Diagnostic output
     END SUBROUTINE ropp_pp_preprocess
  END INTERFACE
  
  INTERFACE
     SUBROUTINE ropp_pp_preprocess_COSMIC(ro_data, config, LCF)
       USE ropp_io_types
       USE ropp_pp_types
       TYPE(ROprof),   INTENT(inout) :: ro_data    ! Radio occultation data
       TYPE(PPConfig), INTENT(inout) :: config     ! Configuration options
       INTEGER, DIMENSION(:), INTENT(inout) :: LCF ! Lost carrier flag 
     END SUBROUTINE ropp_pp_preprocess_COSMIC
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_pp_preprocess_GRASRS(ro_data, config, LCF)
       USE ropp_io_types
       USE ropp_pp_types
       TYPE(ROprof),   INTENT(inout) :: ro_data    ! Radio occultation data
       TYPE(PPConfig), INTENT(inout) :: config     ! Configuration options
       INTEGER, DIMENSION(:), POINTER :: LCF ! Lost carrier flag
     END SUBROUTINE ropp_pp_preprocess_GRASRS
  END INTERFACE


!-------------------------------------------------------------------------------
! 3. Data cut-off
!-------------------------------------------------------------------------------
!****t* Preprocessing/Cutoff
!
! DESCRIPTION
!    Data cutoff based on amplitude, smoothed bending angle and impact parameter
!
! SEE ALSO
!
!****

   INTERFACE
     SUBROUTINE ropp_pp_cutoff(ro_data, config, var1, var2, var3)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_io_types
       USE ropp_pp_types
       TYPE(ROprof),       INTENT(inout) :: ro_data ! Radio occultation data
       TYPE(PPConfig),     INTENT(inout) :: config  ! Configuration options
       REAL(wp), DIMENSION(:), POINTER, OPTIONAL :: var1 ! Extra variable
       REAL(wp), DIMENSION(:), POINTER, OPTIONAL :: var2 ! Extra variable
       INTEGER,  DIMENSION(:), POINTER, OPTIONAL :: var3 ! Extra variable
     END SUBROUTINE ropp_pp_cutoff
  END INTERFACE 


   INTERFACE
     SUBROUTINE ropp_pp_cutoff_amplitude(ro_data, lcf, config)
       USE ropp_io_types
       USE ropp_pp_types
       TYPE(ROprof),       INTENT(inout) :: ro_data ! Radio occultation data
       INTEGER,  DIMENSION(:), POINTER   :: lcf     ! Lost carrier flag
       TYPE(PPConfig),     INTENT(inout) :: config  ! Configuration options
     END SUBROUTINE ropp_pp_cutoff_amplitude
  END INTERFACE
  

END MODULE ropp_pp_preproc
