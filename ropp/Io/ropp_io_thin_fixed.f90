! $Id: ropp_io_thin_fixed.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* Thin/ropp_io_thin_fixed *
!
! NAME
!   ropp_io_thin_fixed - Thin data by interpolating onto fixed altitude
!                        levels.
!
! SYNOPSIS
!   call ropp_io_thin_fixed ( nLev,     FullLev, FullVal, &
!                             nThinLev, ThinLev, ThinVal, &
!                             Method,   sigma )
!
!
! INPUTS
!   nLev      int   Number of original levels
!   FullLev   dflt  Original level values of the data
!   FullVal   dflt  Original data  values
!   nThinLev  int   Number of thinned levels
!   ThinLev   dflt  Thinned level values
!   Method    chr   Required interpolation method. One of:
!                     SGLOG  : Savitzky-Golay smoothing + log interp.
!                     SGLIN  : Savitzky-Golay smoothing + lin interp.
!                     ASGLOG : Adaptive S-G   smoothing + log interp.
!                     ASGLIN : Adaptive S-G   smoothing + lin interp.
!                     LOG    : logarithmic interpolation
!                     LIN    : linear interpolation
!                    Default: LIN
!   sigma     log   .T. indicates error scaling. Applicable to ASGLIN
!                    and ASGLOG only (optional, default: .F.)
!
! OUTPUTS
!   ThinVal   dlft  Thinned data values on ThinLev levels
!
! CALLS
!   ropp_io_thin_sg
!
! CALLED BY
!   ropp_io_thin_select
!
! USES
!   typesizes
!   ropp_io
!   ropp_utils
!
! DESCRIPTION
!   This subroutine thins input data arrays to fixed number of levels.
!   The levels are assumed to be altitudes with data values on those
!   levels, but the algorithm is generic.
!   Several interpolation methods are supported, including:
!    a) Savitzky-Golay pre-smoothing & logarithmic interplation
!       with optional adaptive filter parameters
!    b) Logarithmic interpolation (no pre-smoothing)
!    c) Linear      interpolation (no pre-smoothing)
!   (For LOG interpolations, if either of the pair of values being interpolated
!   are negative, linear interpolation is performed for that fixed height 
!   point.)
!   The input data values are interpolated to the given fixed (thinned)
!   levels in the returned thinned array.
!   Notes:
!   1) The pairs of input (full) levels/values arrays and output (thinned)
!      levels/values arrays must be of the same order, rank and size.
!   2) The Adaptive SG takes account of the relative height 'resolutions'
!      (maximum level separations) of the original and thinned levels, and
!      modifies the width of the filter accordingly and the fitting polynomial
!      order is set to 2 (1 for standard SG). Thus, thining of hi-res data
!      applies more smoothing than for low-res data while still minimising
!      vertical correlations or loss of information.
!   3) If the full (input) height range does not fill the range of thinned
!      (fixed) heights, the out-of-range data values will be set missing.
!   4) If any input data values are flagged with 'missing' values
!      (e.g. -99999) these will be included in the smoothing and/or
!      interpolation which may not preserve any special missing data value for
!      the thinned data values (but ought to remain invalid for that
!      parameter in practice).
!   5) Because of the offsetting for LOG interpolation, data with
!      large negative values (e.g. missing data) may reduce the accuracy
!      of the interpolation at all levels.
!
! TODO
!   Possible future enhancements:
!   1) Replace (or add as option) the LOG interpolation post-SG
!      with a spline interpolation (recommendation in Ref.)
!   2) Add (unsmoothed) spline interpolation as a new option
!   3) Other methods we haven't yet thought of.
!
! REFERENCE
!   Mono-dimensional data thinning for GPS radio occulation data
!   SAF/GRAS/METO/ALG/ROPP/001
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

SUBROUTINE ropp_io_thin_fixed ( nLev,     FullLev, FullVal, &
                                nThinLev, ThinLev, ThinVal, &
                                Method,   sigma )

! Modules

  USE typeSizes,     ONLY: wp => EightByteReal
  USE messages
  USE ropp_io,       ONLY: ropp_io_thin_sg
  USE ropp_utils,    ONLY: sort, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Fixed values

  INTEGER,  PARAMETER  :: nl       = 1         ! SG: No. of points left  of center
  INTEGER,  PARAMETER  :: nr       = 1         ! SG: No. of points right of center
  INTEGER,  PARAMETER  :: m        = 1         ! SG: Order of polynomial to be fitted
  INTEGER,  PARAMETER  :: lder     = 0         ! SG: Order of derivative

! Argument list parameters

  CHARACTER (LEN=*), INTENT(IN)    :: Method
  INTEGER,           INTENT(IN)    :: nLev
  INTEGER,           INTENT(IN)    :: nThinLev
  REAL(wp),          INTENT(IN)    :: FullLev(1:nLev)
  REAL(wp),          INTENT(IN)    :: FullVal(1:nLev)
  REAL(wp),          INTENT(IN)    :: ThinLev(1:nThinLev)
  REAL(wp),          INTENT(OUT)   :: ThinVal(1:nThinLev)
  LOGICAL, OPTIONAL, INTENT(IN)    :: sigma

! Local variables

  CHARACTER (LEN=80) :: number
  REAL(wp) :: Val(1:nLev)
  REAL(wp) :: Lev(1:nLev)
  REAL(wp) :: SortVal(1:nLev)
  REAL(wp) :: SortLev(1:nLev)
  REAL(wp) :: ThinInt
  REAL(wp) :: minLev, maxLev, fullres, thinres
  REAL(wp) :: sigmafac, Val1, Val2
  INTEGER  :: i, j, k,   nsmooth
  LOGICAL  :: ascending, lsigma, logint
  INTEGER  :: idx(1:nLev)
  CHARACTER(len = 256) :: routine


  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_fixed')

!----------------------------------------------------------------------------
! 0. Initialise 
!----------------------------------------------------------------------------

  ThinVal(:)   = ropp_MDFV
  Val(:) = FullVal(:)
  Lev(:) = FullLev(:)

  IF ( .NOT. PRESENT(sigma) ) THEN
    lsigma = .FALSE.
  ELSE
    lsigma = sigma
  END IF
  sigmafac = 1.0_wp

!----------------------------------------------------------------------------
! 1. Sovitzky-Golay pre-smoothing (optional)
!    - with or without adaptive smoothing control
!----------------------------------------------------------------------------

  IF ( INDEX ( Method, "SGLOG" ) == 1 .OR. &
       INDEX ( Method, "SGLIN" ) == 1 ) THEN

    CALL ropp_io_thin_sg ( nLev, Val )

  ELSE IF ( INDEX ( Method, "ASGLOG" ) == 1 .OR. &
            INDEX ( Method, "ASGLIN" ) == 1 ) THEN

! Get resolutions & smoothing distance (half-width)

    fullres = ABS ( MAXVAL ( Lev(2:nLev)         - Lev(1:nLev-1)         ) )
    thinres = ABS ( MINVAL ( ThinLev(2:nThinlev) - ThinLev(1:nThinlev-1) ) )
    nsmooth = MAX ( 2, FLOOR ( thinres / MAX(fullres,0.1_wp) ) )
    nsmooth = 2 * ( nsmooth / 2 )

    IF ( lsigma ) sigmafac = SQRT(nsmooth+1.0_wp)

    WRITE(number, FMT="(I3)") nsmooth+1
    CALL message(msg_diag, "Adaptive SG smoothing over " // TRIM(number) // &
                           " points.")
    IF ( lsigma ) THEN
      WRITE(number, FMT="(F5.2)") sigmafac
      CALL message(msg_diag, "Uncorrelated error reduction by factor " // &
                             TRIM(number))
    END IF

    CALL ropp_io_thin_sg ( nLev, Val, npoints=INT(nsmooth/2), order=2 )
  END IF


!----------------------------------------------------------------------------
! 2. Sort impact parameters. The input impact parameter data will be a
!    monotonic function of height 
!----------------------------------------------------------------------------

  idx = sort(Lev)
  SortVal = Val(idx)
  SortLev = Lev(idx)

  Val = SortVal
  Lev = SortLev

!----------------------------------------------------------------------------
! 3. Interpolation (Log or Lin). Allow for an ascending or descending
!    input profile. The output (thinned) profile will be in the order
!    of the given thinned levels. Do not extrapolate!
!----------------------------------------------------------------------------

  IF ( Lev(1) < Lev(nLev) ) THEN
    ascending = .TRUE.
  ELSE
    ascending = .FALSE.
  END IF

  minLev = MINVAL(Lev)
  maxLev = MAXVAL(Lev)

!----------------------------------------------------------------------------
! 3.1 Linear interpolation
!----------------------------------------------------------------------------

  IF ( INDEX ( Method, "ASGLIN" ) == 1 .OR. &
       INDEX ( Method, "SGLIN"  ) == 1 .OR. &
       INDEX ( Method, "LIN"    ) == 1 ) THEN

      DO k = 1, nThinLev

! Skip upper/lower parts of target profile not covered by source profile

        IF ( ThinLev(k) < minLev .OR. &
             ThinLev(k) > maxLev) CYCLE

! Search for the pair of source levels just above & below the
! target level, allowing for ascending or descending profile.

        j = 2
        IF ( ascending ) THEN
          DO WHILE ( j < nLev .AND. Lev(j) < ThinLev(k) )
            j = j + 1
          END DO
          i = j - 1
        ELSE
          DO WHILE ( j < nLev-1 .AND. Lev(j) > ThinLev(k) )
            j = j + 1
          END DO
          i = j
          j = j - 1
        END IF

! If either source value is missing, either level value is missing, both source levels 
! are the same, or the gap between levels is larger than interval between thinned 
! levels, skip this target level. 

        IF ( k == nThinLev ) THEN
           ThinInt = ABS ( ThinLev(k) - ThinLev(k-1) )
        ELSE
           ThinInt = ABS ( ThinLev(k+1) - ThinLev(k) )
        ENDIF 

        IF ( Val(i) < ropp_MDTV .OR. Val(j) < ropp_MDTV .OR.       &
             Lev(i) < ropp_MDTV .OR. Lev(j) < ropp_MDTV .OR.       &
             Lev(i) == Lev(j) .OR.                                       &
             ABS(Lev(i)-Lev(j)) < 0.001_wp .OR.                          &
             ABS(Lev(i)-Lev(j)) > ThinInt ) CYCLE

! Simple linear interpolation

        ThinVal(k) = ( Val(j) + ( (ThinLev(k) - Lev(j)) &
                   / (Lev(i) - Lev(j)) * (Val(i) - Val(j)) ) ) &
                   * sigmafac
        
      END DO

!----------------------------------------------------------------------------
! 3.2 Logarithic interpolation.
!----------------------------------------------------------------------------

   ELSE IF ( INDEX ( Method, "ASGLOG" ) == 1 .OR. &
             INDEX ( Method, "SGLOG"  ) == 1 .OR. &
             INDEX ( Method, "LOG"    ) == 1 ) THEN
    
      DO k = 1, nThinLev

! Skip upper/lower parts of target profile not covered by source profile

        IF ( ThinLev(k) < minLev .OR. &
             ThinLev(k) > maxLev) CYCLE

! Search for the pair of source levels just above & below the
! target level, allowing for ascending or descending profile.

        j = 2
        IF ( ascending ) THEN
          DO WHILE ( j < nLev .AND. Lev(j) < ThinLev(k) )
            j = j + 1
          END DO
          i = j - 1
        ELSE
          DO WHILE ( j < nLev-1 .AND. Lev(j) > ThinLev(k) )
            j = j + 1
          END DO
          i = j
          j = j - 1
        END IF

! If either source value is missing, either level value is missing, both source levels 
! are the same, or the gap between levels is larger than interval between thinned 
! levels, skip this target level. 

        IF ( k == nThinLev ) THEN
           ThinInt = ABS ( ThinLev(k) - ThinLev(k-1) )
        ELSE
           ThinInt = ABS ( ThinLev(k+1) - ThinLev(k) )
        ENDIF
        
        IF ( Val(i) < ropp_MDTV .OR. Val(j) < ropp_MDTV .OR.       &
             Lev(i) < ropp_MDTV .OR. Lev(j) < ropp_MDTV .OR.       &
             Lev(i) == Lev(j) .OR.                                       &
             ABS(Lev(i)-Lev(j)) < 0.001_wp .OR.                          &
             ABS(Lev(i)-Lev(j)) > ThinInt ) CYCLE

! If both source values are larger than zero do interpolation in log
! space (and then invert target value), else fall back on linear
! interpolation

        IF ( Val(i) > 0.0_wp .AND. &
             Val(j) > 0.0_wp ) THEN
          Val1   = LOG(Val(i))
          Val2   = LOG(Val(j))
          logint = .TRUE.
        ELSE
          Val1   = Val(i)
          Val2   = Val(j)
          logint = .FALSE.
        END IF

        ThinVal(k) = ( Val2 + ( (ThinLev(k) - Lev(j)) &
                   / (Lev(i) - Lev(j)) * (Val1 - Val2) ) )

        IF ( logint ) ThinVal(k) = EXP(ThinVal(k))

        ThinVal(k)  = ThinVal(k) * sigmafac
      END DO

!----------------------------------------------------------------------------
! 3.3 Unknown option
!----------------------------------------------------------------------------

    ELSE
       
      CALL message(msg_diag, "*** Unknown thinning option" // &
                             "The available options for the fixed level " // &
                             "thinning are:")
      CALL message(msg_diag, "     ASGLOG, ASGLIN, SGLOG, SGLIN, LIN and LOG")

    END IF

    CALL message_set_routine(routine)

END SUBROUTINE ropp_io_thin_fixed
