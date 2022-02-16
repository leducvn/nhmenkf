! $Id: ropp_io_thin_select.f90 3551 2013-02-25 09:51:28Z idculv $
!
!****s* Thin/ropp_io_thin_select *
!
! NAME
!   ropp_io_thin_select - select & apply a thinning method
!
! SYNOPSIS
!
!   CALL ropp_io_thin_select ( nLev,     Lev,      Val,     &
!                              nThinLev, ThinLev,  ThinVal, &
!                              Method,   nSampLev, sigma )
!
! INPUTS
!   nLev      int   no. of full levels
!   nThinLev  int   (max) no. of thinned levels
!   Lev       dflt  array of full levels
!   Val       dflt  array of full values
!   ThinLev   dflt  array of thinned levels
!                   (with values set before input for SGLOG,
!                    SGLIN, ASGLOG, ASGLIN, LOG or LIN only)
!   ThinVal   dflt  array of thinned values
!   Method    chr   Method string. One of: SGLOG, SGLIN,
!                    ASGLOG, ASGLIN, LOG, LIN, SAMPLE, NONE
!   sigma     log   .T. indicates error smoothing
!                   (used for ASG* & SG* methods only)
!
! OUTPUTS
!   nSampLev  int   actual no. of thinned levels
!   ThinLev   dflt  array of thinned levels
!                   (with values set on output for SAMPLE only)
!   ThinVal   dflt  array of thinned values
!   Lev       dflt  array of full levels; on exit first nSamplev elements
!                   is a copy of ThinLev, remainder are set 'missing'
!   Val       dflt  array of full values; on exit first nSamplev elements
!                   is a copy of ThinLev, remainder are set 'missing'
!
! CALLS
!   ropp_io_thin_fixed
!   ropp_io_thin_skip
!   where
!
! CALLED BY
!   ropp_io_thin
!
! USES
!   typesizes
!   ropp_io
!   ropp_utils
!
! DESCRIPTION
!   This subroutine selects a thinning method from the input Method
!   string and calls the appropriate routine(s) to apply that
!   method to the input full levels and required output (thinned)
!   levels values to return the thinned values. The supported methods
!   include:
!     ASGLOG : Adaptive S-G smoothing filter with log interpolation
!              to fixed number of thinned levels (see Reference)
!     ASGLIN : Adaptive S-G smoothing filter with linear interpolation
!              to fixed number of thinned levels (see Reference)
!     SGLOG  : Savitzky-Golay smoothing filter with log interpolation
!              to fixed number of thinned levels (see Reference)
!     SGLIN  : Savitzky-Golay smoothing filter with linear interpolation
!              to fixed number of thinned levels (see Reference)
!     LOG    : Logarithmic interpolation to fixed levels (no smoothing)
!     LIN    : Linear      interpolation to fixed levels (no smoothing)
!     SAMPLE : Simple sub-sampling (select or reject 1-in-N) to no more
!              than the array size of the thinned levels.
!     NONE   : explicitly do nothing.
!
! REFERENCE
!   Mono-dimensional data thinning for GPS radio ccultations.
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

SUBROUTINE ropp_io_thin_select ( nLev,     Lev,      Val,     &
                                 nThinLev, ThinLev,  ThinVal, &
                                 Method,   nSampLev, sigma )

  USE typesizes,     ONLY: wp => EightByteReal
  USE messages
  USE ropp_io,       ONLY: ropp_io_thin_fixed, &
                           ropp_io_thin_skip
  USE ropp_utils,    ONLY: WHERE, ropp_MDFV, ropp_MDTV


  IMPLICIT NONE

! Argument list parameters

  CHARACTER (LEN=*), INTENT(IN)    :: Method            ! thinning method
  INTEGER,           INTENT(IN)    :: nLev              ! no. full levels
  INTEGER,           INTENT(IN)    :: nThinLev          ! (max) no thinned levels
  REAL(wp),         INTENT(INOUT) :: Lev(nLev)         ! input (full) levels
  REAL(wp),         INTENT(INOUT) :: Val(nLev)         ! input (full) values
  REAL(wp),          INTENT(INOUT) :: ThinLev(nThinLev) ! output (thinned) levels
  REAL(wp),          INTENT(OUT)   :: ThinVal(nThinLev) ! output (thinned) values
  INTEGER,           INTENT(OUT)   :: nSampLev          ! actual no. of thinned levels
  LOGICAL, OPTIONAL, INTENT(IN)    :: sigma             ! .T. smoothing errrors

! Local variables

  INTEGER                         :: skip1, skip2    ! Skip factors
  INTEGER                         :: i, j            ! Loop counters
  INTEGER                         :: nGood           ! No. of valid values
  INTEGER,  DIMENSION(:), POINTER :: idx  => NULL()  ! Hold WHERE output
  REAL(wp)                        :: LevGood(nLev)   ! Valid input levels
  REAL(wp)                        :: ValGood(nLev)   ! Valid input values
  CHARACTER(len = 256)            :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_thin_select')

! 1. Check for valid levels and values data
! -----------------------------------------

  idx => WHERE ( Lev > ropp_MDTV .AND. Val > ropp_MDTV, nGood )
  
  IF ( nGood > 0 ) THEN
     LevGood(1:nGood) = Lev(idx)
     ValGood(1:nGood) = Val(idx)
  ENDIF
  
! Thinning by one of several methods

! 2. Interpolation
!------------------

    IF ( INDEX ( Method, "SGLOG"  ) == 1 .OR. &
         INDEX ( Method, "SGLIN"  ) == 1 .OR. &
         INDEX ( Method, "ASGLOG" ) == 1 .OR. &
         INDEX ( Method, "ASGLIN" ) == 1 .OR. &
         INDEX ( Method, "LOG"    ) == 1 .OR. &
         INDEX ( Method, "LIN"    ) == 1 ) THEN

! Count valid values; only interpolate if at least two good values on
! valid levels, else just return all thinned values as missing on required
! thinned levels. Interpolation routine will handle profiles containing
! some missing data.

      IF ( nGood >= 2 ) THEN
        CALL ropp_io_thin_fixed ( nGood, LevGood(1:nGood), ValGood(1:nGood), &
                                  nThinLev, ThinLev,  ThinVal, &
                                  Method,   sigma )
      ELSE
        CALL message(msg_diag, "Warning: No valid data, only reducing array")
         ThinVal(1:nThinLev) = ropp_MDFV
      END IF
      nSampLev = nThinLev

      Val(1:nSampLev)      = ThinVal(1:nSampLev)
      Lev(1:nSampLev)      = ThinLev(1:nSampLev)
      Val(nSampLev+1:nlev) = ropp_MDFV
      Lev(nSampLev+1:nlev) = ropp_MDFV

! 3. Simple sub-sampling
!-----------------------

    ELSE IF ( INDEX ( Method, "SAMPLE" ) == 1 ) THEN

      IF ( nGood > nThinLev ) THEN
        CALL ropp_io_thin_skip ( nGood,  nThinLev,        &
                                 skip1, skip2, nSampLev )
        ThinVal(:) = ropp_MDFV
        ThinLev(:) = ropp_MDFV

        j = 0
        DO i = 1, nGood, skip1               ! select every skip1-th sample
          IF ( MOD(i,skip2) == 0 ) CYCLE     ! reject every skip2-th sample
          j = j + 1
          ThinVal(j) = ValGood(i)
          ThinLev(j) = LevGood(i)
        END DO
        nSampLev = j

        Val(1:nSampLev)      = ThinVal(1:nSampLev)
        Lev(1:nSampLev)      = ThinLev(1:nSampLev)
        Val(nSampLev+1:nLev) = ropp_MDFV
        Lev(nSampLev+1:nLev) = ropp_MDFV
      ELSE
        CALL message(msg_diag, "*** Thinning not required")
        nSamplev = nLev
      END IF

! 4. No sampling at all - do nothing
!------------------------------------

    ELSE IF ( INDEX ( Method, "NONE") == 1 ) THEN
      CALL message(msg_diag, "*** Thinning not required")
      nSampLev = nLev

! 5. Unknown method - ignore
!----------------------------

    ELSE 
      CALL message(msg_diag, "*** Thinning method " // TRIM(Method) // &
                                      " not supported")
      nSampLev = nGood

   END IF

   CALL message_set_routine(routine)

END SUBROUTINE ropp_io_thin_select
