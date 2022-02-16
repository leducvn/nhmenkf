! $Id: ropp_fm_copy.f90 4452 2015-01-29 14:42:02Z idculv $

!****m* Modules/ropp_fm_copy *
!
! NAME
!    ropp_fm_copy - Interface module for the ROPP forward models.
!
! SYNOPSIS
!    use ropp_fm_copy
!
! DESCRIPTION
!    Data type/structure copying functions using ROprof structures used by the
!    forward models of ROPP and converting units within the ROprof structure.
!
! SEE ALSO
!    ropp_fm_set_units
!    ropp_fm_roprof2state
!    ropp_fm_state2roprof
!    ropp_fm_roprof2obs
!    ropp_fm_obs2roprof
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

MODULE ropp_fm_copy

!-------------------------------------------------------------------------------
! 1. Standard internal units for foward models
!-------------------------------------------------------------------------------

!****t* Tools/Units
!
! DESCRIPTION
!    Routines for dealing with units.
!
! SEE ALSO
!    ropp_fm_set_units
!
!****

  INTERFACE ropp_fm_set_units
     SUBROUTINE ropp_fm_set_units_roprof_sca(ro_data)
       USE ropp_io_types
       TYPE(ROprof), INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_set_units_roprof_sca
     SUBROUTINE ropp_fm_set_units_roprof_sca2d(ro_data)
       USE ropp_io_types
       TYPE(ROprof2d), INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_set_units_roprof_sca2d
     SUBROUTINE ropp_fm_set_units_roprof_arr(ro_data_set)
       USE ropp_io_types
       TYPE(ROprof), DIMENSION(:), INTENT(inout) :: ro_data_set
     END SUBROUTINE ropp_fm_set_units_roprof_arr
     SUBROUTINE ropp_fm_set_units_roprof_arr2d(ro_data_set)
       USE ropp_io_types
       TYPE(ROprof2d), DIMENSION(:), INTENT(inout) :: ro_data_set
     END SUBROUTINE ropp_fm_set_units_roprof_arr2d
  END INTERFACE

!-------------------------------------------------------------------------------
! 2. Copy functions
!-------------------------------------------------------------------------------

!****t* Tools/Copying2
!
! DESCRIPTION
!    Routines for copying ROprof structures to and from state vectors and
!    observation vectors.
!
! SEE ALSO
!   ropp_fm_roprof2state
!   ropp_fm_state2roprof
!   ropp_fm_roprof2obs
!   ropp_fm_obs2roprof
!
!****

! 2.1 RO data (level 2b,c,d) -> state vector
! ------------------------------------------

  INTERFACE ropp_fm_roprof2state
     SUBROUTINE ropp_fm_roprof2state1d(ro_data, x)
       USE ropp_io_types, ONLY: ROprof
       USE ropp_fm_types, ONLY: State1dFM
       TYPE(ROprof),    INTENT(in)    :: ro_data
       TYPE(State1dFM), INTENT(inout) :: x
     END SUBROUTINE ropp_fm_roprof2state1d
     SUBROUTINE ropp_fm_roprof2state2d(ro_data, x)
       USE ropp_io_types, ONLY: ROprof2d
       USE ropp_fm_types, ONLY: State2dFM
       TYPE(ROprof2d),    INTENT(in)    :: ro_data
       TYPE(State2dFM), INTENT(inout) :: x
     END SUBROUTINE ropp_fm_roprof2state2d
     SUBROUTINE ropp_fm_roprof2d2state1d(ro_data, x)
       USE ropp_io_types, ONLY: ROprof2d
       USE ropp_fm_types, ONLY: State1dFM
       TYPE(ROprof2d),    INTENT(in)    :: ro_data
       TYPE(State1dFM), INTENT(inout) :: x
     END SUBROUTINE ropp_fm_roprof2d2state1d
  END INTERFACE

! 2.2 State vector -> RO data (level 2b,c,d)
! ------------------------------------------

  INTERFACE ropp_fm_state2roprof
     SUBROUTINE ropp_fm_state1d2roprof(x, ro_data)
       USE ropp_io_types, ONLY: ROprof
       USE ropp_fm_types, ONLY: State1dFM
       TYPE(State1dFM), INTENT(in)    :: x
       TYPE(ROprof),    INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_state1d2roprof
     SUBROUTINE ropp_fm_state2d2roprof(x, ro_data)
       USE ropp_io_types, ONLY: ROprof2d
       USE ropp_fm_types, ONLY: State2dFM
       TYPE(State2dFM), INTENT(in)    :: x
       TYPE(ROprof2d),  INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_state2d2roprof
  END INTERFACE

! 2.3 RO data (level 1b, 2a) -> observation vector
! ------------------------------------------------

  INTERFACE ropp_fm_roprof2obs
     SUBROUTINE ropp_fm_roprof2obs1drefrac(ro_data, obs)
       USE ropp_io_types, ONLY: ROprof
       USE ropp_fm_types, ONLY: Obs1dRefrac
       TYPE(ROprof),      INTENT(in)    :: ro_data
       TYPE(Obs1dRefrac), INTENT(inout) :: obs
     END SUBROUTINE ropp_fm_roprof2obs1drefrac
     SUBROUTINE ropp_fm_roprof2obs1dbangle(ro_data, obs)
       USE ropp_io_types, ONLY: ROprof
       USE ropp_fm_types, ONLY: Obs1dBangle
       TYPE(ROprof),      INTENT(in)    :: ro_data
       TYPE(Obs1dBangle), INTENT(inout) :: obs
     END SUBROUTINE ropp_fm_roprof2obs1dbangle
     SUBROUTINE ropp_fm_roprof2obs2dbangle(ro_data, obs)
       USE ropp_io_types
       USE ropp_fm_types
       TYPE(ROprof2d),      INTENT(in)    :: ro_data
       TYPE(Obs1dBangle), INTENT(inout) :: obs
     END SUBROUTINE ropp_fm_roprof2obs2dbangle
  END INTERFACE

! 2.4 Observation vector -> RO data (level 1b, 2a)
! ------------------------------------------------

  INTERFACE ropp_fm_obs2roprof
     SUBROUTINE ropp_fm_obs1drefrac2roprof(obs, ro_data)
       USE ropp_io_types, ONLY: ROprof
       USE ropp_fm_types, ONLY: Obs1dRefrac
       TYPE(Obs1dRefrac), INTENT(in)    :: obs
       TYPE(ROprof),      INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_obs1drefrac2roprof
     SUBROUTINE ropp_fm_obs1dbangle2roprof(obs, ro_data)
       USE ropp_io_types, ONLY: ROprof
       USE ropp_fm_types, ONLY: Obs1dBangle
       TYPE(Obs1dBangle), INTENT(in)    :: obs
       TYPE(ROprof),      INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_obs1dbangle2roprof
     SUBROUTINE ropp_fm_obs2dbangle2roprof(obs, ro_data)
       USE ropp_io_types
       USE ropp_fm_types
       TYPE(Obs1dBangle), INTENT(in)    :: obs
       TYPE(ROprof2d),   INTENT(inout) :: ro_data
     END SUBROUTINE ropp_fm_obs2dbangle2roprof
  END INTERFACE

END MODULE ropp_fm_copy
