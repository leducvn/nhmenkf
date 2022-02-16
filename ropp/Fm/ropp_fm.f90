! $Id: ropp_fm.f90 4452 2015-01-29 14:42:02Z idculv $

!****t* Interface/Modules
!
! SYNOPSIS
!    use ropp_fm
!
!    use ropp_fm_types
!    use ropp_fm_constants
!
! DESCRIPTION
!    Access to the routines in the ROPP Forward Model library ropp_fm is
!    provided by the single Fortran 90 module ropp_fm. The ropp_fm module
!    also includes the module ropp_fm_types, which contains derived type
!    (structure) declarations used by the forward models.
!
!    Meteorological and physical constants are provided in the module
!    ropp_fm_constants; this module is also loaded with the ropp_fm module.
!
! SEE ALSO
!    ropp_fm
!    ropp_fm_types
!    ropp_fm_constants
!
!****

!****m* Modules/ropp_fm *
!
! NAME
!    ropp_fm - Interface module for the ROPP forward models.
!
! SYNOPSIS
!    use ropp_fm
!
! DESCRIPTION
!    This module provides interfaces for all forward model routines in the
!    ROPP Forward Model library.
!
! NOTES
!    At present, the ROPP Forward Model library consists of two one-dimensional
!    forward models calculating and refractivity (as function of geopotential
!    height) and bending angle (as function of impact parameter), respectively.
!
! SEE ALSO
!    ropp_fm_types
!    ropp_fm_refrac_1d
!    ropp_fm_bangle_1d
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

MODULE ropp_fm

!-------------------------------------------------------------------------------
! 1. Other modules
!-------------------------------------------------------------------------------

  USE ropp_fm_types
  USE ropp_fm_constants

!-------------------------------------------------------------------------------
! 1. 1d bending angle forward model
!-------------------------------------------------------------------------------

!****t* ForwardModels/BendingAngle
!
! DESCRIPTION
!    Forward model calculating bending angles
!    (as function of impact parameters).
!
! SEE ALSO
!    ropp_fm_bangle_1d
!    ropp_fm_bangle_1d_ad
!    ropp_fm_bangle_1d_tl
!
!    ropp_fm_abel
!    ropp_fm_abel_ad
!    ropp_fm_abel_tl
!
!****

! 1.1 The forward model
! ---------------------

  INTERFACE ropp_fm_bangle_1d
     SUBROUTINE ropp_fm_bangle_1d(x, y)
       USE ropp_fm_types
       TYPE(State1dFM),   INTENT(in)    :: x
       TYPE(Obs1dBangle), INTENT(inout) :: y
     END SUBROUTINE ropp_fm_bangle_1d
  END INTERFACE

  INTERFACE ropp_fm_bangle_1d_tl
     SUBROUTINE ropp_fm_bangle_1d_tl(x, x_tl, y, y_tl)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),        INTENT(in)  :: x
       TYPE(State1dFM),        INTENT(in)  :: x_tl
       TYPE(Obs1dBangle),      INTENT(in)  :: y
       REAL(wp), DIMENSION(SIZE(y%bangle)), INTENT(out) :: y_tl
     END SUBROUTINE ropp_fm_bangle_1d_tl
  END INTERFACE

  INTERFACE ropp_fm_bangle_1d_ad
     SUBROUTINE ropp_fm_bangle_1d_ad(x, x_ad, y, y_ad)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),        INTENT(in)    :: x
       TYPE(State1dFM),        INTENT(inout) :: x_ad
       TYPE(Obs1dBangle),      INTENT(in)    :: y
       REAL(wp), DIMENSION(:), INTENT(inout) :: y_ad
     END SUBROUTINE ropp_fm_bangle_1d_ad
  END INTERFACE

  INTERFACE ropp_fm_bangle_1d_grad
     SUBROUTINE ropp_fm_bangle_1d_grad(x, y, K)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),          INTENT(in)  :: x
       TYPE(Obs1dBangle),        INTENT(in)  :: y
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: K
     END SUBROUTINE ropp_fm_bangle_1d_grad
  END INTERFACE


! 1.2 Abel transform
! ------------------

  INTERFACE
     SUBROUTINE ropp_fm_abel(nr, refrac, temp, roc, Tgrad_oper, impact, bangle)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: nr
       REAL(wp), DIMENSION(:), INTENT(in)  :: refrac
       REAL(wp), DIMENSION(:), INTENT(in)  :: temp
       REAL(wp),               INTENT(in)  :: roc
       LOGICAL ,               INTENT(in)  :: Tgrad_oper
       REAL(wp), DIMENSION(:), INTENT(in)  :: impact
       REAL(wp), DIMENSION(:), INTENT(out) :: bangle
     END SUBROUTINE ropp_fm_abel
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_abel_tl(nr, refrac, temp, temp_tl, roc, Tgrad_oper, impact, nr_tl, refrac_tl, bangle_tl)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)    :: nr
       REAL(wp), DIMENSION(:), INTENT(in)    :: refrac
       REAL(wp), DIMENSION(:), INTENT(in)    :: temp
       REAL(wp), DIMENSION(:), INTENT(in)    :: temp_tl
       REAL(wp),               INTENT(in)    :: roc 
       LOGICAL ,               INTENT(in)    :: Tgrad_oper
       REAL(wp), DIMENSION(:), INTENT(in)    :: impact
       REAL(wp), DIMENSION(:), INTENT(inout) :: nr_tl
       REAL(wp), DIMENSION(:), INTENT(inout) :: refrac_tl
       REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_tl
     END SUBROUTINE ropp_fm_abel_tl
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_abel_ad(nr, refrac, temp, temp_ad, roc, Tgrad_oper, impact, nr_ad, refrac_ad, bangle_ad)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)    :: nr
       REAL(wp), DIMENSION(:), INTENT(in)    :: refrac
       REAL(wp), DIMENSION(:), INTENT(in)    :: temp
       REAL(wp), DIMENSION(:), INTENT(inout) :: temp_ad
       REAL(wp),               INTENT(in)    :: roc
       LOGICAL ,               INTENT(in)    :: Tgrad_oper
       REAL(wp), DIMENSION(:), INTENT(in)    :: impact
       REAL(wp), DIMENSION(:), INTENT(inout) :: nr_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: refrac_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_ad
     END SUBROUTINE ropp_fm_abel_ad
  END INTERFACE

!-------------------------------------------------------------------------------
! 2. 1d refractivity forward model
!-------------------------------------------------------------------------------

!****t* ForwardModels/Refractivity
!
! DESCRIPTION
!    Forward model calculating refractivity (as function of geopotential).
!
! SEE ALSO
!    ropp_fm_refrac_1d
!    ropp_fm_refrac_1d_new
!    ropp_fm_refrac_1d_ad
!    ropp_fm_refrac_1d_new_ad
!    ropp_fm_refrac_1d_tl
!    ropp_fm_refrac_1d_new_tl
!    ropp_fm_refrac_1d_grad
!
!****

  INTERFACE ropp_fm_refrac_1d
     SUBROUTINE ropp_fm_refrac_1d(x, y, nonull)
       USE ropp_fm_types
       TYPE(State1dFM),   INTENT(in)    :: x
       TYPE(Obs1dRefrac), INTENT(inout) :: y
       LOGICAL, OPTIONAL, INTENT(in)    :: nonull
     END SUBROUTINE ropp_fm_refrac_1d
  END INTERFACE

  INTERFACE ropp_fm_refrac_1d_new
     SUBROUTINE ropp_fm_refrac_1d_new(x, y)
       USE ropp_fm_types
       TYPE(State1dFM),   INTENT(in)    :: x
       TYPE(Obs1dRefrac), INTENT(inout) :: y
     END SUBROUTINE ropp_fm_refrac_1d_new
  END INTERFACE

  INTERFACE ropp_fm_refrac_1d_tl
     SUBROUTINE ropp_fm_refrac_1d_tl(x, x_tl, y, y_tl)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),        INTENT(in)    :: x
       TYPE(State1dFM),        INTENT(in)    :: x_tl
       TYPE(Obs1dRefrac),      INTENT(inout) :: y
       REAL(wp), DIMENSION(:), INTENT(out)   :: y_tl
     END SUBROUTINE ropp_fm_refrac_1d_tl
  END INTERFACE

  INTERFACE ropp_fm_refrac_1d_new_tl
     SUBROUTINE ropp_fm_refrac_1d_new_tl(x, x_tl, y, y_tl)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),        INTENT(in)  :: x
       TYPE(State1dFM),        INTENT(in)  :: x_tl
       TYPE(Obs1dRefrac),      INTENT(in)  :: y
       REAL(wp), DIMENSION(:), INTENT(out) :: y_tl
     END SUBROUTINE ropp_fm_refrac_1d_new_tl
  END INTERFACE

  INTERFACE ropp_fm_refrac_1d_ad
     SUBROUTINE ropp_fm_refrac_1d_ad(x, x_ad, y, y_ad)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),        INTENT(in)    :: x
       TYPE(State1dFM),        INTENT(inout) :: x_ad
       TYPE(Obs1dRefrac),      INTENT(in)    :: y
       REAL(wp), DIMENSION(:), INTENT(inout) :: y_ad
     END SUBROUTINE ropp_fm_refrac_1d_ad
  END INTERFACE

  INTERFACE ropp_fm_refrac_1d_new_ad
     SUBROUTINE ropp_fm_refrac_1d_new_ad(x, x_ad, y, y_ad)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),        INTENT(in)    :: x
       TYPE(State1dFM),        INTENT(inout) :: x_ad
       TYPE(Obs1dRefrac),      INTENT(in)    :: y
       REAL(wp), DIMENSION(:), INTENT(inout) :: y_ad
     END SUBROUTINE ropp_fm_refrac_1d_new_ad
  END INTERFACE

  INTERFACE ropp_fm_refrac_1d_grad
     SUBROUTINE ropp_fm_refrac_1d_grad(x, y, K)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       TYPE(State1dFM),          INTENT(in)    :: x
       TYPE(Obs1dRefrac),        INTENT(inout) :: y
       REAL(wp), DIMENSION(:,:), INTENT(inout) :: K
     END SUBROUTINE ropp_fm_refrac_1d_grad
  END INTERFACE

!------------------------------------------------------------------------------
! 3. Model-specific state conversion routines - ECMWF hybrid-sigma coordinate
!------------------------------------------------------------------------------

!****t* ModelConversions/Model_ecmwf
!
! DESCRIPTION
!    Conversion routines for background data on ECWMF hybrid-sigma vertical
!    levels to p,q,T as function of geopotential height and state vector.
!
! SEE ALSO
!    ropp_fm_state2state_ecmwf
!    ropp_fm_state2state_ecmwf_ad
!    ropp_fm_state2state_ecmwf_tl
!
!****

  INTERFACE ropp_fm_state2state_ecmwf
     SUBROUTINE ropp_fm_state2state_ecmwf_1d(x)
       USE ropp_fm_types, ONLY: State1dFM
       TYPE(State1dFM), INTENT(inout) :: x
     END SUBROUTINE ropp_fm_state2state_ecmwf_1d
     SUBROUTINE ropp_fm_state2state_ecmwf_2d(x)
       USE ropp_fm_types, ONLY: State2dFM
       TYPE(State2dFM), INTENT(inout) :: x
     END SUBROUTINE ropp_fm_state2state_ecmwf_2d
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_state2state_ecmwf_tl(x, x_tl)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(in)    :: x
       TYPE(State1dFM), INTENT(inout) :: x_tl
     END SUBROUTINE ropp_fm_state2state_ecmwf_tl
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_state2state_ecmwf_ad(x, x_ad)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(in)    :: x
       TYPE(State1dFM), INTENT(inout) :: x_ad
     END SUBROUTINE ropp_fm_state2state_ecmwf_ad
  END INTERFACE

!------------------------------------------------------------------------------
! 4. Model-specific state conversion routines - MetOffice geopotential levels
!------------------------------------------------------------------------------

!****t* ModelConversions/Model_meto
!
! DESCRIPTION
!    Conversion routines for background data on Met Office geopotential height
!    based vertical levels to p,q,T as function of geopotential height and
!    state vector.
!
! SEE ALSO
!    ropp_fm_state2state_meto
!    ropp_fm_state2state_meto_ad
!    ropp_fm_state2state_meto_tl
!
!****

  INTERFACE
     SUBROUTINE ropp_fm_state2state_meto(x)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(inout) :: x
     END SUBROUTINE ropp_fm_state2state_meto
  END INTERFACE

 INTERFACE
     SUBROUTINE ropp_fm_state2state_meto_tl(x, x_tl)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(in)    :: x
       TYPE(State1dFM), INTENT(inout) :: x_tl
     END SUBROUTINE ropp_fm_state2state_meto_tl
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_state2state_meto_ad(x, x_ad)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(in)    :: x
       TYPE(State1dFM), INTENT(inout) :: x_ad
     END SUBROUTINE ropp_fm_state2state_meto_ad
  END INTERFACE

!-------------------------------------------------------------------------------
! 5. Vertical interpolation
!-------------------------------------------------------------------------------

!****t* VerticalLevels/Interpolation
!
! DESCRIPTION
!    Forward model routines dealing with vertical interpolation.
!
! SEE ALSO
!    ropp_fm_interpol          - Linear interpolation
!    ropp_fm_interpol_ad
!    ropp_fm_interpol_tl
!
!    ropp_fm_interpol_log      - Logarithmic interpolation
!    ropp_fm_interpol_log_ad
!    ropp_fm_interpol_log_tl
!
!****

! 5.1 Vertical interpolation
! --------------------------

  INTERFACE
     SUBROUTINE ropp_fm_interpol(x, newx, array, interp)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp), DIMENSION(:), INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp), DIMENSION(:), INTENT(out) :: interp
     END SUBROUTINE ropp_fm_interpol
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_interpol_tl(x, newx, array, x_tl, array_tl, interp_tl)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp), DIMENSION(:), INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp), DIMENSION(:), INTENT(in)  :: x_tl
       REAL(wp), DIMENSION(:), INTENT(in)  :: array_tl
       REAL(wp), DIMENSION(:), INTENT(out) :: interp_tl
     END SUBROUTINE ropp_fm_interpol_tl
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_interpol_ad(x, newx, array, x_ad, array_ad, interp_ad)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)    :: x
       REAL(wp), DIMENSION(:), INTENT(in)    :: newx
       REAL(wp), DIMENSION(:), INTENT(in)    :: array
       REAL(wp), DIMENSION(:), INTENT(inout) :: x_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: array_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: interp_ad
     END SUBROUTINE ropp_fm_interpol_ad
  END INTERFACE

! 5.2 Vertical interpolation of log()
! -----------------------------------

  INTERFACE
     SUBROUTINE ropp_fm_interpol_log(x, newx, array, interp)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp), DIMENSION(:), INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp), DIMENSION(:), INTENT(out) :: interp
     END SUBROUTINE ropp_fm_interpol_log
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_interpol_log_tl(x,newx,array, x_tl,array_tl,interp_tl)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)  :: x
       REAL(wp), DIMENSION(:), INTENT(in)  :: newx
       REAL(wp), DIMENSION(:), INTENT(in)  :: array
       REAL(wp), DIMENSION(:), INTENT(in)  :: x_tl
       REAL(wp), DIMENSION(:), INTENT(in)  :: array_tl
       REAL(wp), DIMENSION(:), INTENT(out) :: interp_tl
     END SUBROUTINE ropp_fm_interpol_log_tl
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_interpol_log_ad(x,newx,array, x_ad,array_ad,interp_ad)
       USE typesizes, ONLY: wp => EightByteReal
       REAL(wp), DIMENSION(:), INTENT(in)    :: x
       REAL(wp), DIMENSION(:), INTENT(in)    :: newx
       REAL(wp), DIMENSION(:), INTENT(in)    :: array
       REAL(wp), DIMENSION(:), INTENT(inout) :: x_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: array_ad
       REAL(wp), DIMENSION(:), INTENT(inout) :: interp_ad
     END SUBROUTINE ropp_fm_interpol_log_ad
  END INTERFACE

!-------------------------------------------------------------------------------
! 6. Copy functions
!-------------------------------------------------------------------------------

!****t* Tools/Copying
!
! DESCRIPTION
!    Routines for copying (and checking) data in and out of state and
!    observation vectors.
!
! SEE ALSO
!   ropp_fm_state2state
!   ropp_fm_obs2obs
!
!****

! 6.1 State vector -> state vector
! --------------------------------

  INTERFACE ropp_fm_state2state
     SUBROUTINE ropp_fm_state2state_1d(from_state, to_state)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(in)    :: from_state
       TYPE(State1dFM), INTENT(inout) :: to_state
     END SUBROUTINE ropp_fm_state2state_1d
  END INTERFACE

! 6.2 State vector -> state vector (assignment)
! ---------------------------------------------

  INTERFACE ASSIGNMENT(=)
     SUBROUTINE ropp_fm_assign_state_1d(to_state, from_state)
       USE ropp_fm_types
       TYPE(State1dFM), INTENT(inout) :: to_state
       TYPE(State1dFM), INTENT(in)    :: from_state
     END SUBROUTINE ropp_fm_assign_state_1d
  END INTERFACE

! 6.3 Observation vector -> observation vector
! --------------------------------------------

  INTERFACE ropp_fm_obs2obs
     SUBROUTINE ropp_fm_obs2obs_1dbangle(from_obs, to_obs)
       USE ropp_fm_types
       TYPE(Obs1DBangle), INTENT(in)    :: from_obs
       TYPE(Obs1DBangle), INTENT(inout) :: to_obs
     END SUBROUTINE ropp_fm_obs2obs_1dbangle
     SUBROUTINE ropp_fm_obs2obs_1drefrac(from_obs, to_obs)
       USE ropp_fm_types
       TYPE(Obs1DRefrac), INTENT(in)    :: from_obs
       TYPE(Obs1DRefrac), INTENT(inout) :: to_obs
     END SUBROUTINE ropp_fm_obs2obs_1drefrac
  END INTERFACE

! 6.4 Observation vector -> observation vector (assignment)
! ---------------------------------------------------------

  INTERFACE ASSIGNMENT(=)
     SUBROUTINE ropp_fm_assign_obs_1dbangle(to_obs, from_obs)
       USE ropp_fm_types
       TYPE(Obs1DBangle), INTENT(inout) :: to_obs
       TYPE(Obs1DBangle), INTENT(in)    :: from_obs
     END SUBROUTINE ropp_fm_assign_obs_1dbangle
     SUBROUTINE ropp_fm_assign_obs_1drefrac(to_obs, from_obs)
       USE ropp_fm_types
       TYPE(Obs1DRefrac), INTENT(inout) :: to_obs
       TYPE(Obs1DRefrac), INTENT(in)    :: from_obs
     END SUBROUTINE ropp_fm_assign_obs_1drefrac
  END INTERFACE

!-------------------------------------------------------------------------------
! 7. Memory freeing functions
!-------------------------------------------------------------------------------

!****t* Tools/Memory
!
! DESCRIPTION
!    Routines for freeing memory
!
! SEE ALSO
!   ropp_fm_free
!
!****

  INTERFACE ropp_fm_free
     SUBROUTINE ropp_fm_free_obs1drefrac(var)
       USE ropp_fm_types, ONLY: Obs1dRefrac
       TYPE(Obs1dRefrac), INTENT(inout) :: var
     END SUBROUTINE ropp_fm_free_obs1drefrac
     SUBROUTINE ropp_fm_free_obs1dbangle(var)
       USE ropp_fm_types, ONLY: Obs1dBangle
       TYPE(Obs1dBangle), INTENT(inout) :: var
     END SUBROUTINE ropp_fm_free_obs1dbangle
     SUBROUTINE ropp_fm_free_state1dfm(var)
       USE ropp_fm_types, ONLY: State1dFM
       TYPE(State1dFM), INTENT(inout) :: var
     END SUBROUTINE ropp_fm_free_state1dfm
     SUBROUTINE ropp_fm_free_state2dfm(var)
       USE ropp_fm_types, ONLY: State2dFM
       TYPE(State2dFM), INTENT(inout) :: var
     END SUBROUTINE ropp_fm_free_state2dfm
  END INTERFACE

!-------------------------------------------------------------------------------
! 8. 2d Operator Routines
!-------------------------------------------------------------------------------

!****t* ForwardModels/BendingAngle2d
!
! DESCRIPTION
!    Routines for 2d Bending Angle Operator
!
! SEE ALSO
!
!****

  INTERFACE ropp_fm_bangle_2d
     SUBROUTINE ropp_fm_bangle_2d(x, y)
       USE ropp_fm_types
       TYPE(State2dFM),   INTENT(inout) :: x
       TYPE(Obs1dBangle), INTENT(inout) :: y
     END SUBROUTINE ropp_fm_bangle_2d
  END INTERFACE


  INTERFACE ropp_fm_bangle_2d_tl
     SUBROUTINE ropp_fm_bangle_2d_tl(x, x_tl, y, y_tl)
       USE ropp_fm_types
       TYPE(State2dFM),   INTENT(in)    :: x, x_tl
       TYPE(Obs1dBangle), INTENT(inout) :: y
       TYPE(Obs1dBangle), INTENT(inout) :: y_tl
     END SUBROUTINE ropp_fm_bangle_2d_tl
  END INTERFACE


  INTERFACE ropp_fm_bangle_2d_ad
     SUBROUTINE ropp_fm_bangle_2d_ad(x, x_ad, y, y_ad)
       USE ropp_fm_types
       TYPE(State2dFM),   INTENT(in)    :: x
       TYPE(State2dFM),   INTENT(inout) :: x_ad
       TYPE(Obs1dBangle), INTENT(in)    :: y
       TYPE(Obs1dBangle), INTENT(inout) :: y_ad
     END SUBROUTINE ropp_fm_bangle_2d_ad
  END INTERFACE


  INTERFACE ropp_fm_alpha2drk
     SUBROUTINE ropp_fm_alpha2drk(kobs,klev,khoriz,ksplit,pdsep,pa,prefrac,  &
                                  pradius,pnr,proc,pz_2d,pa_path,palpha,prtan)
       USE typesizes, ONLY: wp => EightByteReal
       INTEGER,  INTENT(IN)  :: kobs
       INTEGER,  INTENT(IN)  :: klev
       INTEGER,  INTENT(IN)  :: khoriz
       INTEGER,  INTENT(IN)  :: ksplit
       REAL(wp), INTENT(IN)  :: pdsep
       REAL(wp), INTENT(IN)  :: pa(kobs)
       REAL(wp), INTENT(IN)  :: prefrac(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pradius(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pnr(klev,khoriz)
       REAL(wp), INTENT(IN)  :: proc
       REAL(wp), INTENT(IN)  :: pz_2d
       REAL(wp), INTENT(OUT) :: pa_path(kobs,2)
       REAL(wp), INTENT(OUT) :: palpha(kobs)
       REAL(wp), INTENT(OUT) :: prtan(kobs)
     END SUBROUTINE ropp_fm_alpha2drk
  END INTERFACE


  INTERFACE ropp_fm_alpha2drk_tl
     SUBROUTINE ropp_fm_alpha2drk_tl(kobs,klev,khoriz,ksplit,pdsep,pa,    &
                          prefrac,prefrac_prime,pradius,pradius_prime,    &
                          pnr,pnr_prime,proc,pz_2d,pa_path,pa_path_prime, &
                          palpha,palpha_prime)
       USE typesizes, ONLY: wp => EightByteReal
       INTEGER,  INTENT(IN)  :: kobs
       INTEGER,  INTENT(IN)  :: klev
       INTEGER,  INTENT(IN)  :: khoriz
       INTEGER,  INTENT(IN)  :: ksplit
       REAL(wp), INTENT(IN)  :: pdsep
       REAL(wp), INTENT(IN)  :: pa(kobs)
       REAL(wp), INTENT(IN)  :: prefrac(klev,khoriz)
       REAL(wp), INTENT(IN)  :: prefrac_prime(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pradius(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pradius_prime(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pnr(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pnr_prime(klev,khoriz)
       REAL(wp), INTENT(IN)  :: proc
       REAL(wp), INTENT(IN)  :: pz_2d
       REAL(wp), INTENT(OUT) :: pa_path(kobs,2)
       REAL(wp), INTENT(OUT) :: pa_path_prime(kobs,2)
       REAL(wp), INTENT(OUT) :: palpha(kobs)
       REAL(wp), INTENT(OUT) :: palpha_prime(kobs)
     END SUBROUTINE ropp_fm_alpha2drk_tl
  END INTERFACE


  INTERFACE ropp_fm_alpha2drk_ad
     SUBROUTINE ropp_fm_alpha2drk_ad(kobs,klev,khoriz,ksplit,pdsep,pa,   &
                              prefrac,prefrac_hat,pradius,pradius_hat,   &
                              pnr,pnr_hat,proc,pz_2d,pa_path_hat,palpha_hat)
       USE typesizes, ONLY: wp => EightByteReal
       INTEGER,  INTENT(IN)    :: kobs
       INTEGER,  INTENT(IN)    :: klev
       INTEGER,  INTENT(IN)    :: khoriz
       INTEGER,  INTENT(IN)    :: ksplit
       REAL(wp), INTENT(IN)    :: pdsep
       REAL(wp), INTENT(IN)    :: pa(kobs)
       REAL(wp), INTENT(IN)    :: prefrac(klev,khoriz)
       REAL(wp), INTENT(INOUT) :: prefrac_hat(klev,khoriz)
       REAL(wp), INTENT(IN)    :: pradius(klev,khoriz)
       REAL(wp), INTENT(INOUT) :: pradius_hat(klev,khoriz)
       REAL(wp), INTENT(IN)    :: pnr(klev,khoriz)
       REAL(wp), INTENT(INOUT) :: pnr_hat(klev,khoriz)
       REAL(wp), INTENT(IN)    :: proc
       REAL(wp), INTENT(IN)    :: pz_2d
       REAL(wp), INTENT(INOUT) :: pa_path_hat(kobs,2)
       REAL(wp), INTENT(INOUT) :: palpha_hat(kobs)
     END SUBROUTINE ropp_fm_alpha2drk_ad
  END INTERFACE


  INTERFACE ropp_fm_gpspderivs
     SUBROUTINE ropp_fm_gpspderivs(klev,khoriz,kk,pdsep,ptheta_min,ptheta_max,&
                           ptheta_tan,prtan,pamult,prefrac,pradius,py,pdydh)
       USE typesizes, ONLY: wp => EightByteReal
       INTEGER,  INTENT(IN)  :: klev           ! no. of refractivity levels
       INTEGER,  INTENT(IN)  :: khoriz         ! no. of horizontal locations
       INTEGER,  INTENT(IN)  :: kk
       REAL(wp), INTENT(IN)  :: pdsep
       REAL(wp), INTENT(IN)  :: ptheta_min
       REAL(wp), INTENT(IN)  :: ptheta_max
       REAL(wp), INTENT(IN)  :: ptheta_tan
       REAL(wp), INTENT(IN)  :: prtan
       REAL(wp), INTENT(IN)  :: pamult
       REAL(wp), INTENT(IN)  :: prefrac(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pradius(klev,khoriz)
       REAL(wp), INTENT(IN)  :: py(4)
       REAL(wp), INTENT(OUT) :: pdydh(4)
     END SUBROUTINE ropp_fm_gpspderivs
  END INTERFACE


  INTERFACE ropp_fm_gpspderivs_tl
     SUBROUTINE ropp_fm_gpspderivs_tl(klev,khoriz,kk,pdsep,ptheta_min,ptheta_max,   &
                           ptheta_tan,prtan,prtan_prime,pamult,prefrac,     &
                           prefrac_prime,pradius,pradius_prime,             &
                           py,py_prime,pdydh,pdydh_prime)
       USE typesizes, ONLY: wp => EightByteReal
       INTEGER,  INTENT(IN)  :: klev           ! no. of refractivity levels
       INTEGER,  INTENT(IN)  :: khoriz         ! no. of horizontal locations
       INTEGER,  INTENT(IN)  :: kk
       REAL(wp), INTENT(IN)  :: pdsep
       REAL(wp), INTENT(IN)  :: ptheta_min
       REAL(wp), INTENT(IN)  :: ptheta_max
       REAL(wp), INTENT(IN)  :: ptheta_tan
       REAL(wp), INTENT(IN)  :: prtan
       REAL(wp), INTENT(IN)  :: prtan_prime
       REAL(wp), INTENT(IN)  :: pamult
       REAL(wp), INTENT(IN)  :: prefrac(klev,khoriz)
       REAL(wp), INTENT(IN)  :: prefrac_prime(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pradius(klev,khoriz)
       REAL(wp), INTENT(IN)  :: pradius_prime(klev,khoriz)
       REAL(wp), INTENT(IN)  :: py(4)
       REAL(wp), INTENT(IN)  :: py_prime(4)
       REAL(wp), INTENT(OUT) :: pdydh(4)
       REAL(wp), INTENT(OUT) :: pdydh_prime(4)
     END SUBROUTINE ropp_fm_gpspderivs_tl
  END INTERFACE


  INTERFACE ropp_fm_gpspderivs_ad
     SUBROUTINE ropp_fm_gpspderivs_ad(klev,khoriz,kk,pdsep,ptheta_min,ptheta_max,   &
                           ptheta_tan,prtan,prtan_hat,pamult,prefrac,       &
                           prefrac_hat,pradius,pradius_hat,py,py_hat,pdydh_hat)
       USE typesizes, ONLY: wp => EightByteReal
       INTEGER,  INTENT(IN)    :: klev           ! no. of refractivity levels
       INTEGER,  INTENT(IN)    :: khoriz         ! no. of horizontal locations
       INTEGER,  INTENT(IN)    :: kk
       REAL(wp), INTENT(IN)    :: pdsep
       REAL(wp), INTENT(IN)    :: ptheta_min
       REAL(wp), INTENT(IN)    :: ptheta_max
       REAL(wp), INTENT(IN)    :: ptheta_tan
       REAL(wp), INTENT(IN)    :: prtan
       REAL(wp), INTENT(INOUT) :: prtan_hat
       REAL(wp), INTENT(IN)    :: pamult
       REAL(wp), INTENT(IN)    :: prefrac(klev,khoriz)
       REAL(wp), INTENT(INOUT) :: prefrac_hat(klev,khoriz)
       REAL(wp), INTENT(IN)    :: pradius(klev,khoriz)
       REAL(wp), INTENT(INOUT) :: pradius_hat(klev,khoriz)
       REAL(wp), INTENT(IN)    :: py(4)
       REAL(wp), INTENT(INOUT) :: py_hat(4)
       REAL(wp), INTENT(INOUT) :: pdydh_hat(4)
     END SUBROUTINE ropp_fm_gpspderivs_ad
  END INTERFACE

!-------------------------------------------------------------------------------
! 9. Compressibility Routines
!-------------------------------------------------------------------------------

!****t* ForwardModels/Compressibility
!
! DESCRIPTION
!    Routines for calculating geopotential height and wet and compressibilities
!    in 1D and 2D Forward Models, accounting for non-ideal gas effects.
!
! SEE ALSO
!   ropp_fm_refrac_1d
!   ropp_fm_refrac_new_1d
!   ropp_fm_refrac_1d_ad
!   ropp_fm_refrac_1d_new_ad
!   ropp_fm_refrac_1d_tl
!   ropp_fm_refrac_1d_new_tl
!   ropp_fm_bangle_1d
!   ropp_fm_bangle_1d_ad
!   ropp_fm_bangle_1d_tl
!   ropp_fm_bangle_2d
!   ropp_fm_bangle_2d_ad
!   ropp_fm_bangle_2d_tl
!
!****

  INTERFACE
     SUBROUTINE ropp_fm_compress(x, z_geop, zcomp_dry_inv, zcomp_wet_inv)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_fm_constants
       IMPLICIT NONE
       TYPE(State1dFM),              INTENT(in)  :: x
       REAL(wp), DIMENSION(x%n_lev), INTENT(out) :: z_geop
       REAL(wp), DIMENSION(x%n_lev), INTENT(out) :: zcomp_dry_inv
       REAL(wp), DIMENSION(x%n_lev), INTENT(out) :: zcomp_wet_inv
     END SUBROUTINE ropp_fm_compress
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_compress_tl&
     (x, x_tl, z_geop, z_geop_tl,zcomp_dry_inv, zcomp_dry_inv_tl,zcomp_wet_inv, zcomp_wet_inv_tl)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_fm_constants
       IMPLICIT NONE
       TYPE(State1dFM),              INTENT(in)    :: x
       TYPE(State1dFM),              INTENT(in)    :: x_tl
       REAL(wp), DIMENSION(x%n_lev), INTENT(out)   :: z_geop
       REAL(wp), DIMENSION(x%n_lev), INTENT(out)   :: z_geop_tl
       REAL(wp), DIMENSION(x%n_lev), INTENT(out)   :: zcomp_dry_inv
       REAL(wp), DIMENSION(x%n_lev), INTENT(out)   :: zcomp_dry_inv_tl
       REAL(wp), DIMENSION(x%n_lev), INTENT(out)   :: zcomp_wet_inv
       REAL(wp), DIMENSION(x%n_lev), INTENT(out)   :: zcomp_wet_inv_tl
     END SUBROUTINE ropp_fm_compress_tl
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_compress_ad&
     (x, x_ad, z_geop_ad,zcomp_dry_inv_ad,zcomp_wet_inv_ad)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_fm_constants
       IMPLICIT NONE
       TYPE(State1dFM),              INTENT(in)      :: x
       TYPE(State1dFM),              INTENT(inout)   :: x_ad
       REAL(wp), DIMENSION(x%n_lev), INTENT(inout)   :: z_geop_ad
       REAL(wp), DIMENSION(x%n_lev), INTENT(inout)   :: zcomp_dry_inv_ad
       REAL(wp), DIMENSION(x%n_lev), INTENT(inout)   :: zcomp_wet_inv_ad
     END SUBROUTINE ropp_fm_compress_ad
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_compress_single&
     (temp,pres,shum, zcomp_dry_inv, zcomp_wet_inv, zcomp1, zcomp2, zcomp3)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_fm_constants
       IMPLICIT NONE
       REAL(wp),           INTENT(in)  :: temp           ! Temperature
       REAL(wp),           INTENT(in)  :: pres           ! Pressure
       REAL(wp),           INTENT(in)  :: shum           ! Specific humidity
       REAL(wp),           INTENT(out) :: zcomp_dry_inv  ! Inverse of dry comp
       REAL(wp),           INTENT(out) :: zcomp_wet_inv  ! Inverse of wet comp
       REAL(wp), OPTIONAL, INTENT(out) :: zcomp1         ! Intermediate comp factor
       REAL(wp), OPTIONAL, INTENT(out) :: zcomp2         ! Intermediate comp factor
       REAL(wp), OPTIONAL, INTENT(out) :: zcomp3         ! Intermediate comp factor
     END SUBROUTINE ropp_fm_compress_single
  END INTERFACE

  INTERFACE
     SUBROUTINE ropp_fm_compress_single_tl&
     (temp, pres, shum, temp_tl, pres_tl, shum_tl, zcomp_dry_inv, zcomp_dry_inv_tl, &
     zcomp_wet_inv, zcomp_wet_inv_tl, zcomp1, zcomp1_tl, zcomp2, zcomp2_tl, zcomp3, zcomp3_tl)
       USE typesizes, ONLY: wp => EightByteReal
       USE ropp_fm_types
       USE ropp_fm_constants
       IMPLICIT NONE
       REAL(wp), INTENT(in)    :: temp              ! Temperature
       REAL(wp), INTENT(in)    :: pres              ! Pressure
       REAL(wp), INTENT(in)    :: shum              ! Specific humidity
       REAL(wp), INTENT(in)    :: temp_tl           ! Temperature TL
       REAL(wp), INTENT(in)    :: pres_tl           ! Pressure TL
       REAL(wp), INTENT(in)    :: shum_tl           ! Specific humidity TL
       REAL(wp), INTENT(out)   :: zcomp_dry_inv     ! inverse of dry comp
       REAL(wp), INTENT(out)   :: zcomp_dry_inv_tl  ! inverse of dry comp
       REAL(wp), INTENT(out)   :: zcomp_wet_inv     ! inverse of wet comp
       REAL(wp), INTENT(out)   :: zcomp_wet_inv_tl  ! inverse of wet comp
       REAL(wp),OPTIONAL, INTENT(out)  :: zcomp1    ! comp. factor 1
       REAL(wp),OPTIONAL, INTENT(out)  :: zcomp1_tl ! comp. factor 1 TL
       REAL(wp),OPTIONAL, INTENT(out)  :: zcomp2    ! comp. factor 2
       REAL(wp),OPTIONAL, INTENT(out)  :: zcomp2_tl ! comp. factor 2 TL
       REAL(wp),OPTIONAL, INTENT(out)  :: zcomp3    ! comp. factor 3
       REAL(wp),OPTIONAL, INTENT(out)  :: zcomp3_tl ! comp. factor 3 TL
     END SUBROUTINE ropp_fm_compress_single_tl
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_fm_compress_single_ad&
    (temp,temp_ad,pres,pres_ad,shum,shum_ad,zcomp_dry_inv_ad,zcomp_wet_inv_ad,&
    zcomp1_opt_ad,zcomp2_opt_ad,zcomp3_opt_ad)
      USE typesizes, ONLY: wp => EightByteReal
      USE ropp_fm_types
      USE ropp_fm_constants
      IMPLICIT NONE
      REAL(wp), INTENT(in)               :: temp              ! Temperature
      REAL(wp), INTENT(in)               :: pres              ! Pressure
      REAL(wp), INTENT(in)               :: shum              ! Specific humidity
      REAL(wp), INTENT(inout)            :: temp_ad           ! Temperature AD
      REAL(wp), INTENT(inout)            :: pres_ad           ! Pressure AD
      REAL(wp), INTENT(inout)            :: shum_ad           ! Specific humidity AD
      REAL(wp), INTENT(inout)            :: zcomp_dry_inv_ad  ! inverse of dry comp
      REAL(wp), INTENT(inout)            :: zcomp_wet_inv_ad  ! inverse of wet comp
      REAL(wp), OPTIONAL, INTENT(inout)  :: zcomp1_opt_ad     ! Intermediate comp. factor AD
      REAL(wp), OPTIONAL, INTENT(inout)  :: zcomp2_opt_ad     ! Intermediate comp. factor AD
      REAL(wp), OPTIONAL, INTENT(inout)  :: zcomp3_opt_ad     ! Intermediate comp. factor AD
    END SUBROUTINE ropp_fm_compress_single_ad
  END INTERFACE

  INTERFACE ropp_fm_compress_2d
    SUBROUTINE ropp_fm_compress_2d(x, z_geop, zcomp_dry_inv, zcomp_wet_inv)
      USE typesizes, ONLY: wp => EightByteReal
      USE ropp_fm_types
      USE ropp_fm_constants
      IMPLICIT NONE
      TYPE(State2dFM), INTENT(in)                         :: x              ! State vector
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: z_geop         ! adjusted geop height
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: zcomp_dry_inv  ! inverse of dry comp
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: zcomp_wet_inv  ! inverse of wet comp
    END SUBROUTINE ropp_fm_compress_2d
  END INTERFACE

  INTERFACE ropp_fm_compress_2d_tl
    SUBROUTINE ropp_fm_compress_2d_tl&
    (x, x_tl, z_geop, z_geop_tl,zcomp_dry_inv, zcomp_dry_inv_tl,zcomp_wet_inv, zcomp_wet_inv_tl)
      USE typesizes, ONLY: wp => EightByteReal
      USE ropp_fm_types
      USE ropp_fm_constants
      IMPLICIT NONE
      TYPE(State2dFM), INTENT(in)                         :: x                 ! State vector
      TYPE(State2dFM), INTENT(in)                         :: x_tl              ! State vector
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: z_geop            ! adjusted geop height
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: z_geop_tl         ! adjusted geop height
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: zcomp_dry_inv     ! inverse of dry comp
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: zcomp_dry_inv_tl  ! inverse of dry comp
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: zcomp_wet_inv     ! inverse of wet comp
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(out) :: zcomp_wet_inv_tl  ! inverse of wet comp
    END SUBROUTINE ropp_fm_compress_2d_tl
  END INTERFACE

  INTERFACE ropp_fm_compress_2d_ad
    SUBROUTINE ropp_fm_compress_2d_ad&
    (x, x_ad, z_geop_ad,zcomp_dry_inv_ad,zcomp_wet_inv_ad)
      USE typesizes, ONLY: wp => EightByteReal
      USE ropp_fm_types
      USE ropp_fm_constants
      IMPLICIT NONE
      TYPE(State2dFM), INTENT(in)                            :: x                 ! State vector
      TYPE(State2dFM), INTENT(inout)                         :: x_ad              ! State vector
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(inout)  :: z_geop_ad         ! adjusted geop height
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(inout)  :: zcomp_dry_inv_ad  ! inverse of dry comp
      REAL(wp), DIMENSION(x%n_lev,x%n_horiz), INTENT(inout)  :: zcomp_wet_inv_ad  ! inverse of wet comp
    END SUBROUTINE ropp_fm_compress_2d_ad
  END INTERFACE

!-------------------------------------------------------------------------------
! 10. Ionospheric Routines
!-------------------------------------------------------------------------------

!****t* ForwardModels/Ionosphere
!
! DESCRIPTION
!   Routines for calculating bending angle produced by model Chapman layer 
!   ionosphere in 1D Forward Model.
!
! SEE ALSO
!   ropp_fm_bangle_1d
!   ropp_fm_bangle_1d_ad
!   ropp_fm_bangle_1d_tl
!
!****

  INTERFACE
    SUBROUTINE ropp_fm_iono_bangle( &
      Ne_max, R_peak, H_width, n_L1, impact, bangle)
      USE typesizes, ONLY: wp => EightByteReal
      IMPLICIT NONE
      REAL(wp), INTENT(in)                  :: Ne_max
      REAL(wp), INTENT(in)                  :: R_peak
      REAL(wp), INTENT(in)                  :: H_width
      INTEGER                               :: n_L1
      REAL(wp), DIMENSION(:), INTENT(in)    :: impact
      REAL(wp), DIMENSION(:), INTENT(inout) :: bangle
    END SUBROUTINE ropp_fm_iono_bangle
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_fm_iono_bangle_tl( &
               Ne_max,    R_peak,    H_width, &
               Ne_max_tl, R_peak_tl, H_width_tl, &
               n_L1, impact, bangle, bangle_tl)
      USE typesizes, ONLY: wp => EightByteReal
      IMPLICIT NONE
      REAL(wp), INTENT(in)                  :: Ne_max
      REAL(wp), INTENT(in)                  :: R_peak
      REAL(wp), INTENT(in)                  :: H_width
      REAL(wp), INTENT(in)                  :: Ne_max_tl
      REAL(wp), INTENT(in)                  :: R_peak_tl
      REAL(wp), INTENT(in)                  :: H_width_tl
      INTEGER                               :: n_L1
      REAL(wp), DIMENSION(:), INTENT(in)    :: impact
      REAL(wp), DIMENSION(:), INTENT(in)    :: bangle
      REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_tl
    END SUBROUTINE ropp_fm_iono_bangle_tl
  END INTERFACE

  INTERFACE
    SUBROUTINE ropp_fm_iono_bangle_ad( &
               Ne_max,    R_peak,    H_width, &
               Ne_max_ad, R_peak_ad, H_width_ad, &
               n_L1, impact, bangle, bangle_ad)
      USE typesizes, ONLY: wp => EightByteReal
      IMPLICIT NONE
      REAL(wp), INTENT(in)                  :: Ne_max
      REAL(wp), INTENT(in)                  :: R_peak
      REAL(wp), INTENT(in)                  :: H_width
      REAL(wp), INTENT(inout)               :: Ne_max_ad
      REAL(wp), INTENT(inout)               :: R_peak_ad
      REAL(wp), INTENT(inout)               :: H_width_ad
      INTEGER                               :: n_L1
      REAL(wp), DIMENSION(:), INTENT(in)    :: impact
      REAL(wp), DIMENSION(:), INTENT(in)    :: bangle
      REAL(wp), DIMENSION(:), INTENT(inout) :: bangle_ad
    END SUBROUTINE ropp_fm_iono_bangle_ad
  END INTERFACE

  INTERFACE
    FUNCTION ropp_fm_zorro(lg) RESULT(z)
      USE typesizes, ONLY: wp => EightByteReal
      IMPLICIT NONE
      REAL(wp), DIMENSION(:)                :: lg
      REAL(wp), DIMENSION(SIZE(lg))         :: z
    END FUNCTION ropp_fm_zorro
  END INTERFACE

  INTERFACE
    FUNCTION ropp_fm_dzorro_dlg(lg) RESULT(dz_by_dl)
      USE typesizes, ONLY: wp => EightByteReal
      IMPLICIT NONE
      REAL(wp), DIMENSION(:)                :: lg
      REAL(wp), DIMENSION(SIZE(lg))         :: dz_by_dl
    END FUNCTION ropp_fm_dzorro_dlg
  END INTERFACE

  INTERFACE
    FUNCTION ropp_fm_asinh(x) RESULT(asinh)
      USE typesizes, ONLY: wp => EightByteReal
      IMPLICIT NONE
      REAL(wp), DIMENSION(:)                :: x
      REAL(wp), DIMENSION(SIZE(x))          :: asinh
    END FUNCTION ropp_fm_asinh
  END INTERFACE

!-------------------------------------------------------------------------------
! 11. Common utilities
!-------------------------------------------------------------------------------

  INTERFACE ropp_fm_version
    FUNCTION ropp_fm_version() RESULT (version)
      CHARACTER (LEN=40) :: version
    END FUNCTION ropp_fm_version
  END INTERFACE ropp_fm_version

END MODULE ropp_fm


