! $Id: ropp_io_thin.f90 4452 2015-01-29 14:42:02Z idculv $
!
!****s* Thin/ropp_io_thin *
!
! NAME
!   ropp_io_thin - Thin RO profiles
!
! SYNOPSIS
!   CALL ropp_io_thin ( ROdata, ThinFile[, impactalt=..][, ranchk=..] )
!
! INPUT
!   ROdata     struc  RO profile structure (full profiles)
!   ThinFile   chr    path/name of a file containing thinning
!                     control information or a non-negative integer
!                     representing the max. no. of points to thin
!                     by sub-sampling. If input name is blank
!                     or "0", no thinning is attempted
!   impactalt   log   flag for thinning on impact altitudes
!                     .T. to thin on impact altitudes
!                     (= impact parameter - radius of curvature - undulation)
!                     .F. to thin on impact heights (default)
!                     (= impact parameter - radius of curvature)
!   ranchk     log    range check flag (optional) Valid values are:
!                     .T. to perform a range check before & after thinning
!                         (default) or
!                     .F. to skip range checking
!                     This flag should not normally be used or if present, be
!                     set to default .T.. It is intended for use only when
!                     deliberately invalid data is explicitly required to
!                     survive thinning, such as when testing the thinner
!                     itself or downstream Q/C.
!                     WARNING: the correct operation of these thinning
!                     procedures is not guaranteed if range checking is
!                     disabled!!
!
! OUTPUT
!   ROdata     struc  RO profile structure (thinned profiles)
!
! FILES
!   'ThinFile' is a plain-text file which should be formatted as follows:
!   Title=<title>
!   Method=<method>
!   Nlevels=<nlev>
!   Hlevels=
!   <level_1>
!   <level_2>
!   ...
!   <level_nlev>
!
!   where:
!    <title>  is a free-text description up to 70 characters
!    <method> specified the thinning method, and must be one of:
!      SGLOG  : Savitzky-Golay smoothing filter with log interpolation
!               to fixed number of thinned levels (see Reference)
!      SGLIN  : Savitzky-Golay smoothing filter with linear interpolation
!               to fixed number of thinned levels (see Reference)
!      ASGLOG : Adaptive S-G smoothing filter with log interpolation
!               to fixed number of thinned levels (see Reference)
!      ASGLIN : Adaptive S-G smoothing filter with linear interpolation
!               to fixed number of thinned levels (see Reference)
!      LOG    : Logarithmic interpolation to fixed levels (no smoothing)
!      LIN    : Linear      interpolation to fixed levels (no smoothing)
!      SAMPLE : Simple sub-sampling (select or reject 1-in-N) to no more
!               than the array size of the thinned levels.
!      NONE   : no thinning (default if file cannot be read or its name
!               blank or "0")
!    <nlev>     is the (maximum or fixed) number of levels for output.
!               If zero, disables thinning.
!    <levels_n> are the required set of fixed levels (for SGLOG, SGLIN,
!               LOG & LIN methods) as impact heights (impact parameter
!               minus radius of curvature). Values should be monotonically
!               increasing. May be omitted for SAMPLE and NONE but the key
!               'Hlevels=' must always be present.
!
! CALLS
!   ropp_io_ascend
!   ropp_io_rangecheck
!   ropp_io_thin_select
!   geometric2geopotential
!   message
!   message_get_routine
!   message_set_routine
!
! USES
!   typesizes
!   ropp_io
!   ropp_io_types
!   geodesy
!
! DESCRIPTION
!   This subroutine thins a complete RO structure profiles - ie Level 1b
!   L1+L2+C+O BA vs IP, Level 2a N vs ht, Level 2b T,q,H vs geopot ht
!   (NB it does not thin Level 1a (SNR/Phase vs time) profiles.)
!   The routine first reads & parses a thinning control file to determine
!   the required thinning method (several are supported)
!   and the number (maximum or fixed, depending on the method) and
!   if appropriate the values for the fixed levels (as impact heights).
!   The same thinning details are applied to all types of profiles.
!   If the file cannot be read sucessfully, no thinning will be performed.
!   If a logaritmic interpolation method is selected (SGLOG, ASGLOG, LOG),
!   note that this is only applied to profile parameters that behave
!   logarithmically with height (viz. bending angle, refractivity, pressure and
!   humidity); all other parameters are interploated in linear mode
!   (as if SGLIN, ASGLIN or LIN had been specified).
!
! REFERENCE
!   Monodimensional data thinning for GPS radio occultations.
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

SUBROUTINE ropp_io_thin ( ROdata,    & ! (inout)
                          ThinFile,  & ! (inout)
                          impactalt, & ! (in/optional)
                          ranchk )     ! (in/optional)

! Modules

  USE messages
  USE typesizes,     ONLY: wp => EightByteReal
  USE ropp_io_types, ONLY: ROprof
! USE ropp_io,       not_this => ropp_io_thin
  USE geodesy,       ONLY: geometric2geopotential
  USE ropp_utils,    ONLY: WHERE, ropp_MDFV, ropp_MDTV

  IMPLICIT NONE

! Argument list parameters

  TYPE(ROprof)                  :: ROdata
  CHARACTER (LEN=*), INTENT(IN) :: ThinFile
  LOGICAL, OPTIONAL, INTENT(IN) :: impactalt
  LOGICAL, OPTIONAL, INTENT(IN) :: ranchk

! Fixed parameters

  CHARACTER (LEN=*), PARAMETER :: ThinVer = "v3.1"    ! Thinner software version

  INTEGER,  PARAMETER :: ThinUnit = 22                ! Thinning control file stream
  LOGICAL,  PARAMETER :: sigma    = .TRUE.            ! Error parameter flag
  REAL(wp), PARAMETER :: dtor = 0.017453292519943295_wp ! Pi/180

! Local variables

  CHARACTER (LEN=80)     :: line       ! Line input from thin control file
  CHARACTER (LEN=65)     :: Title      ! Thinning title
  CHARACTER (LEN=10)     :: Method     ! Thinning method
  CHARACTER (LEN=10)     :: number     ! Number as chr string
  CHARACTER (LEN=10)     :: LocMethod  ! Local thinning method for parameters not linear in log space
  INTEGER                :: ierr       ! I/O error code
  INTEGER                :: i          ! Loop counter
  INTEGER                :: len_line   ! Length of NLEVELS= line in thinning file
  INTEGER                :: i1, i2, i3 ! Loop control
  INTEGER                :: nLev       ! No. of full input levels
  INTEGER                :: nThinLev   ! (Max) no. thinned output levels
  INTEGER                :: nMSampLev  ! no. thinned output levels of actually thinned profile, disregard missing
  INTEGER                :: nSampLev   ! Actual no. thinned output levels (for missing data just maximum)
  INTEGER                :: iPole      ! Which pole we are near, if either
  REAL(wp), ALLOCATABLE  :: Lev(:)     ! Full levels
  REAL(wp), ALLOCATABLE  :: Levtmp(:)  ! Temporary copy of impact parameters
  REAL(wp), ALLOCATABLE  :: Val(:)     ! Full values
  REAL(wp), ALLOCATABLE  :: ThinLev(:) ! Thinned levels
  REAL(wp), ALLOCATABLE  :: ThinVal(:) ! Thinned values
  REAL(wp), ALLOCATABLE  :: x_coord(:) ! Polar x_coord for occultations near pole
  REAL(wp), ALLOCATABLE  :: y_coord(:) ! Polar y_coord for occultations near pole
  REAL(wp), ALLOCATABLE  :: lon_temp(:)! Work space for lon
  REAL(wp), ALLOCATABLE  :: lat_temp(:)! Work space for lat
  REAL(wp), ALLOCATABLE  :: lon_save(:)! Store space for lon
  REAL(wp), ALLOCATABLE  :: lat_save(:)! Store space for lat
  REAL(wp)               :: N          ! Refractivity N=(n-1).10^6
  REAL(wp)               :: H          ! Height
  REAL(wp)               :: Ngood      ! Last 'good; refractivity value
  REAL(wp)               :: Hgood      ! Height for Ngood
  REAL(wp)               :: N300       ! (log of) Refractivity for 300 N-units
  REAL(wp)               :: roc        ! Radius of curvature
  REAL(wp)               :: undulation ! Undulation
  LOGICAL                :: fixlev     ! .T. if thinning on fixed levels (interpolation)
  LOGICAL                :: refrac_thinned   ! .T. when refractivity thinned
  LOGICAL                :: lontrans   ! Longitude transformation done yes/no
  LOGICAL                :: iranchk    ! range check internally
  LOGICAL                :: iimpactalt ! impact alt check internally

! holds WHERE output

  INTEGER, DIMENSION(:), POINTER :: idx  => NULL()
  INTEGER                        :: nidx
  CHARACTER(len = 256)           :: routine

! 0. Initialise/defaults
!=======================

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_io_thin')

  nThinLev       = 0           ! default: no thinning
  Method         = "SG"        ! default: SG algorithm
  Title          = "no title"
  fixlev         = .FALSE.
  refrac_thinned = .FALSE.
  ierr           = 0
  nidx           = 0

  IF (PRESENT(impactalt)) THEN
    iimpactalt = impactalt
  ELSE
      iimpactalt = .FALSE.
  ENDIF

  IF (PRESENT(ranchk)) THEN
    iranchk = ranchk
  ELSE
    iranchk = .TRUE.
  ENDIF

! (option) ensure standard missing data values are used; filter out
! any levels with missing vertical coordinates

  IF (iranchk) CALL ropp_io_rangecheck ( ROdata )

! 1. Load required thinning options
!==================================

! 1.1 If 'file' is a valid integer, assume SAMPLE method to that max
!     (if 1 or greater) or NONE (if 0 or negative). A blank name is
!     also treated as NONE.
! -------------------------------------------------------------------

  IF ( ThinFile /= " " ) READ ( ThinFile, "(BN,I10)", IOSTAT=ierr) nThinLev
  IF ( ierr == 0 ) THEN
    IF ( nThinLev > 0 ) THEN
      WRITE ( number, "(I10)" ) nThinLev
      Title  = "Sampling to a maximum of "//TRIM(ADJUSTL(number))//" levels"
      Method = "SAMPLE"
    ELSE
      Title  = "Thinning disabled"
      Method = "NONE"
    ENDIF

! 1.2 From named file
! -------------------

  ELSE

! Open levels file

    OPEN ( UNIT=ThinUnit, &
           FILE=ThinFile, &
           STATUS="OLD",  &
           ACTION="READ", &
           IOSTAT=ierr )

! Loop over header lines looking for keywords for
! No. levels, thinning method & levels

    IF ( ierr == 0 ) THEN
      DO
        READ ( UNIT=ThinUnit, &
                FMT="(A)",    &
             IOSTAT=ierr ) line

        IF ( ierr < 0 ) THEN
          CALL message(msg_error, "*** Unexpected EOF reading " // &
                                  TRIM(ThinFile))
          EXIT
        ELSE IF ( ierr /= 0 ) THEN
          CALL message(msg_error, "*** Error reading thinning keys from " // &
                                  TRIM(ThinFile))
          EXIT
        END IF
        line = ADJUSTL(line)
        CALL To_Upper ( line(1:7) )
        i = INDEX ( line, "TITLE=" )
        IF ( i > 0 ) Title = line(7:66)
        i = INDEX ( line, "METHOD=" )
        IF ( i > 0 ) THEN
          Method = TRIM(ADJUSTL(line(8:17)))
          CALL To_Upper ( Method )

          IF ( INDEX ( Method, "EUM"    ) == 1 ) Method = "LOG"

          IF ( INDEX ( Method, "SGLOG"  ) == 0 .AND. &
               INDEX ( Method, "SGLIN"  ) == 0 .AND. &
               INDEX ( Method, "ASGLIN" ) == 0 .AND. &
               INDEX ( Method, "ASGLOG" ) == 0 .AND. &
               INDEX ( Method, "LOG"    ) == 0 .AND. &
               INDEX ( Method, "LIN"    ) == 0 .AND. &
               INDEX ( Method, "SAMPLE" ) == 0 ) Method = "NONE"
        END IF
        i = INDEX ( line, "NLEVELS=" )
        IF ( i > 0 ) THEN      
           len_line = LEN(TRIM(line))
           IF (IACHAR(line(len_line:len_line)) == 13) len_line = len_line - 1 ! Strip off carriage return
           SELECT CASE (len_line)
             CASE (9)
             READ ( line(9:len_line), FMT='(i1)', IOSTAT=ierr ) i
             CASE (10)
             READ ( line(9:len_line), FMT='(i2)', IOSTAT=ierr ) i
             CASE (11)
             READ ( line(9:len_line), FMT='(i3)', IOSTAT=ierr ) i
             CASE DEFAULT
             READ ( line(9:len_line), FMT='(i3)', IOSTAT=ierr ) i
           END SELECT
           IF ( ierr == 0 .AND. i >= 0 .AND. i <= 5000) THEN
             nThinLev = i
           ELSE
             CALL message(msg_error,  &
              "** Error reading number of levels from " // TRIM(ThinFile))
             Title = "No thinning"
             Method = "NONE"
           ENDIF
        END IF
        i = INDEX ( line, "HLEVELS=" )
        IF ( i > 0 ) EXIT
     END DO
     IF ( INDEX ( Method, "NONE") == 1) THEN
       Title    = "No thinning"
       nThinLev = 0
     END IF
  ELSE
    CALL message(msg_error, "*** Failed to open thinning control file " // TRIM(ThinFile))
  END IF

END IF

! Save current thinning method if previous method was 'UNKNOWN'
! but preserve it if not, and current is 'NONE'; add thinner s/w version,
! but only if the thinning has been done by ROPP.

  IF ( (TRIM(ROdata%Thin_Method) /= "EUMETSAT thinned") .AND. &
       (TRIM(ROdata%Thin_Method) /= "Unthinned data") ) THEN  ! It hasn't come from EUMETSAT

    IF ( ROdata%Thin_Method(1:8) == "UNKNOWN" ) THEN
      ROdata%Thin_Method = TRIM(Method) // " (" // TRIM(Title) // ")"
    ELSE IF ( Method /= "NONE" ) THEN
      ROdata%Thin_Method = TRIM(Method) // " (" // TRIM(Title) // ")"
    END IF
    i = INDEX ( ROdata%Thin_Method, "[" )
    IF ( i == 0 ) i = LEN_TRIM(ROdata%Thin_Method) + 2
    ROdata%Thin_Method(i:) = "[" // TRIM(ThinVer) // "]"

  END IF

  CALL message(msg_diag, "Thinning method: " //  TRIM(ROdata%Thin_Method) )

! Allocate working memory for thinned arrays &
! read thinned impact heights for appropriate methods

  IF ( nThinLev > 0 ) THEN
    ALLOCATE ( ThinVal(nThinLev), &
               ThinLev(nThinLev), &
               STAT=ierr )
    IF ( ierr == 0 ) THEN
      IF ( INDEX ( Method(1:2), "SG"  ) == 1  .OR. &
           INDEX ( Method(1:3), "ASG" ) == 1 .OR. &
           INDEX ( Method(1:3), "LOG" ) == 1 .OR. &
           INDEX ( Method(1:3), "LIN" ) == 1 ) THEN

        READ ( ThinUnit, FMT='(F8.3)', IOSTAT=ierr ) ThinLev

        IF ( ierr == 0 ) THEN
          fixlev = .TRUE.
        ELSE
          WRITE(number, "(I4)") nThinLev
          CALL message(msg_error, "*** Error reading " // TRIM(number) // "thinning levels")
          nThinLev = 0
        END IF
      ELSE
        ThinLev = ropp_MDFV
      END IF

    ELSE
      CALL message(msg_fatal, "*** Failed to allocate memory for thinned arrays")
    END IF

  END IF
  CLOSE ( UNIT=ThinUnit, IOSTAT=ierr )

  IF ( nThinLev < 1 ) THEN
    CALL message_set_routine(routine)
    RETURN
  ENDIF

! Local thinning method for parameters not linear in log space

  LocMethod = Method
  IF ( INDEX ( Method(1:2), "SG" ) == 1 ) THEN
    LocMethod = "SGLIN"
  ELSE IF ( INDEX ( Method(1:3), "ASG" ) == 1 ) THEN
    LocMethod = "ASGLIN"
  ELSE IF ( INDEX ( Method(1:3), "LOG" ) == 1 ) THEN
    LocMethod = "LIN"
  END IF
  IF ( LocMethod /= Method ) CALL message(msg_diag, "Local thinning method: "//TRIM(LocMethod))

! 2. Thin bending angles as fn of impact parameter
!=================================================

! Convert input thinned levels from impact height to impact parameter

  IF ( ROdata%GEOref%RoC > ropp_MDTV ) THEN

    roc = ROdata%GEOref%RoC

  ELSE

    roc = 6378137.0_wp
    CALL message(msg_warn, "Error reading radius of curvature, using " // &
                           "default value 6378.137 km for computing output thinned levels")

  ENDIF

  IF (iimpactalt) THEN
     IF (.NOT. ( ROdata%GEOref%undulation > ropp_MDTV ) ) THEN
        CALL message(msg_fatal, "Undulation not defined for thinning!")
     ENDIF
     idx => WHERE (ThinLev > ropp_MDTV, nidx) 
     IF (nidx > 0) ThinLev(idx) = ThinLev(idx) + roc + ROdata%GEOref%undulation
  ELSE
     idx => WHERE (ThinLev > ropp_MDTV, nidx) 
     IF (nidx > 0) ThinLev(idx) = ThinLev(idx) + roc
  ENDIF

  nLev = ROdata%Lev1b%Npoints

  nMSampLev = nLev
  IF ( nLev >= nThinLev ) THEN
    ALLOCATE ( Val(nLev),    &
               Lev(nLev),    &
               Levtmp(nLev), &
               STAT=ierr )
    IF ( ierr /= 0 ) THEN
      CALL message(msg_fatal, "Failed to allocate memory for Level 1b arrays")
    ENDIF

    nSampLev = nLev

! 2.1 L1 BA, error and quality
!-----------------------------

    CALL message(msg_diag, "Thinning L1 Bending Angles...")
    Levtmp(:) = ROdata%Lev1b%Impact_L1(1:nlev)

    Lev(:) = ROdata%Lev1b%Impact_L1(1:nlev)
    Val(:) = ROdata%Lev1b%Bangle_L1(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Lev,      Val,      &
                               nThinLev, ThinLev,  ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev1b%Bangle_L1(1:nLev) = Val(:)
    ROdata%Lev1b%Impact_L1(1:nLev) = Lev(:)

    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning L1 Bending Angle errors...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_L1_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, sigma )

    ROdata%Lev1b%Bangle_L1_sigma(1:nLev) = Val(:)
!--
    CALL message(msg_diag, "Thinning L1 Bending Angle quality...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_L1_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev1b%Bangle_L1_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 2.2 L2 BA, error and quality
!-----------------------------

    CALL message(msg_diag, "Thinning L2 Bending Angles...")
    Levtmp(:) = ROdata%Lev1b%Impact_L2(1:nlev)

    Lev(:) = ROdata%Lev1b%Impact_L2(1:nlev)
    Val(:) = ROdata%Lev1b%Bangle_L2(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Lev,      Val,      &
                               nThinLev, ThinLev,  ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev1b%Bangle_L2(1:nLev) = Val(:)
    ROdata%Lev1b%Impact_L2(1:nLev) = Lev(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning L2 Bending Angle errors...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_L2_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, sigma )

    ROdata%Lev1b%Bangle_L2_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning L2 Bending Angle quality...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_L2_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev1b%Bangle_L2_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 2.3 Optimised BA, error and quality
!------------------------------------

    CALL message(msg_diag, "Thinning Optimised Bending Angles...")

    Levtmp(:) = ROdata%Lev1b%Impact_Opt(1:nlev)

    Lev(:) = ROdata%Lev1b%Impact_Opt(1:nlev)
    Val(:) = ROdata%Lev1b%Bangle_Opt(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Lev,      Val,      &
                               nThinLev, ThinLev,  ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev1b%Bangle_Opt(1:nLev) = Val(:)
    ROdata%Lev1b%Impact_Opt(1:nLev) = Lev(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Optimised Bending Angle errors...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_Opt_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, sigma )

    ROdata%Lev1b%Bangle_Opt_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Optimised Bending Angle quality...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_Opt_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev1b%Bangle_Opt_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 2.4 Corrected BA, error and quality
!------------------------------------

    CALL message(msg_diag, "Thinning Corrected Bending Angles...")

    Levtmp(:) = ROdata%Lev1b%Impact(1:nlev) ! used for all remaining L1b profile parameters

    Lev(:) = ROdata%Lev1b%Impact(1:nlev)
    Val(:) = ROdata%Lev1b%Bangle(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Lev,       Val,      &
                               nThinLev, ThinLev,   ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev1b%Bangle(1:nLev) = Val(:)
    ROdata%Lev1b%Impact(1:nLev) = Lev(:)
    nMSampLev = MIN ( nMSampLev, nSampLev)
!--
    CALL message(msg_diag, "Thinning Corrected Bending Angle errors...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, sigma )

    ROdata%Lev1b%Bangle_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Corrected Bending Angle quality...")

    Lev(:) = Levtmp(:)
    Val(:) = ROdata%Lev1b%Bangle_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev1b%Bangle_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 2.5 Latitude, Longitude & Azimuths
!-----------------------------------

! If we're within 5deg of a pole, transform lat and lon to local polar
! stereographic coordinates before interpolating, then transform back

    IF (MAXVAL(ROdata%Lev1b%Lat_tp, mask=(ROdata%Lev1b%Lat_tp > ropp_MDTV)) .GT. 85.0) THEN

      iPole =  1  ! Near north pole

    ELSE IF (MINVAL(ROdata%Lev1b%Lat_tp, mask=(ROdata%Lev1b%Lat_tp > ropp_MDTV)) .LT. -85.0) THEN

      iPole = -1  ! Near south pole

    ELSE

      iPole= 0    ! Not near either pole

    ENDIF

    IF (iPole .NE. 0) THEN  ! Make a polar stereographic projection

      ALLOCATE(x_coord(nLev), y_coord(nLev), lon_temp(nLev), lat_temp(nLev), lon_save(nLev), lat_save(nLev), STAT=ierr)
      IF ( ierr /= 0 ) THEN
        CALL message(msg_fatal, "Failed to allocate memory for x_coord, y_coord, lon_temp, lat_temp, lon_save and lat_save.")
      ENDIF

! Local x- and y-coordinates of tangent point

      x_coord(:) = ropp_MDFV
      y_coord(:) = ropp_MDFV

      idx => WHERE ((ROdata%Lev1b%Lat_tp > ropp_MDTV) .AND. (ROdata%Lev1b%Lon_tp > ropp_MDTV), nidx)

      IF (nidx > 0) THEN
        x_coord(idx) = &
                       tan(dtor*(45.0_wp-iPole*0.5_wp*ROdata%Lev1b%Lat_tp(idx)))* &
                       cos(dtor*ROdata%Lev1b%Lon_tp(idx))
        y_coord(idx) = &
                       tan(dtor*(45.0_wp-iPole*0.5_wp*ROdata%Lev1b%Lat_tp(idx)))* &
                       sin(dtor*ROdata%Lev1b%Lon_tp(idx))
      ENDIF

      Lev(:) = Levtmp(:)
      Val(:) = x_coord

      CALL message(msg_diag, "Thinning Level 1b x-cmpt of polar projection...")

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

      x_coord(1:nLev) = Val(:) ! Store

      Lev(:) = Levtmp(:)
      Val(:) = y_coord

      CALL message(msg_diag, "Thinning Level 1b y-cmpt of polar projection...")

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

      y_coord(1:nLev) = Val(:) ! Store

! Invert to generate thinned lats and lons

      lat_save(:) = ropp_MDFV
      lon_save(:) = ropp_MDFV

      idx => WHERE ((x_coord > ropp_MDTV) .AND. (y_coord > ropp_MDTV), nidx)

      IF (nidx > 0) THEN

        lat_save(idx) = &
                       (2.0_wp*iPole/dtor)*atan( &
                       (1.0_wp-sqrt(x_coord(idx)**2+y_coord(idx)**2))/ &
                       (1.0_wp+sqrt(x_coord(idx)**2+y_coord(idx)**2)))

        lon_save(idx) = atan2(y_coord(idx), x_coord(idx)+tiny(1.0_wp))/dtor

      ENDIF

! Wait until azimuths processed before setting 
! ROdata%Lev1b%Lat_tp(1:nLev) and ROdata%Lev1b%Lon_tp(1:nLev)
! to these saved values lat_save and lon_save.

! Azimuth also needs special treatment, as it changes by 180deg if the
! occultation passes right over the pole.
!
! Handle this by working with the sin and cos of the azimuth,
! which change sign as you pass over the pole, and inverting.
! (Because they both change sign, the tangent stays the same.)
!
! Note that ROPP azimuths are between 0 and 360 deg, and we need psi to
! be between -180 and 180 if atan2(sin(psi), cos(psi)) is to equal psi.

      idx => WHERE ((ROdata%Lev1b%Azimuth_tp > ropp_MDTV) .AND. &
                    (ROdata%Lev1b%Azimuth_tp > 180), nidx)

      IF (nidx > 0) ROdata%Lev1b%Azimuth_tp(idx) = &
                    ROdata%Lev1b%Azimuth_tp(idx) - 360.0_wp

! Calculate components in polar projection space of unit vectors in the 
! direction of the local azimuth, by rotating the local northward vector by psi.
! Strictly, these equations should involve iPole, but in the final expression 
! for azimuth these factors cancel out, so we omit them from the start.

      x_coord(:) = ropp_MDFV
      y_coord(:) = ropp_MDFV

      idx => WHERE ((ROdata%Lev1b%Azimuth_tp > ropp_MDTV) .AND. &
                    (ROdata%Lev1b%Lon_tp     > ropp_MDTV), nidx)

      IF (nidx > 0) THEN

        x_coord(idx) = cos(dtor*(ROdata%Lev1b%Lon_tp(idx) - &
                                 ROdata%Lev1b%Azimuth_tp(idx)))

        y_coord(idx) = sin(dtor*(ROdata%Lev1b%Lon_tp(idx) - &
                                 ROdata%Lev1b%Azimuth_tp(idx)))

      ENDIF

      Lev(:) = Levtmp(:)
      Val(:) = x_coord

      CALL message(msg_diag, "Thinning Level 1b x-cmpt of projected azimuths ...")

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

      x_coord(1:nLev) = Val(:) ! Store

      Lev(:) = Levtmp(:)
      Val(:) = y_coord

      CALL message(msg_diag, "Thinning Level 1b y-cmpt of projected azimuths ...")

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

      y_coord(1:nLev) = Val(:) ! Store

! Invert to generate thinned azimuths

      lon_temp(:) = ropp_MDFV
      lat_temp(:) = ropp_MDFV

      idx => WHERE ((x_coord > ropp_MDTV) .AND. &
                    (y_coord > ropp_MDTV) .AND. &
                    (lon_save > ropp_MDTV), nidx)

      IF (nidx > 0) THEN

        lon_temp(idx) = (cos(dtor*lon_save(idx))*x_coord(idx)  + &
                         sin(dtor*lon_save(idx))*y_coord(idx)) ! cos(azimuth) (no need to normalise)

        lat_temp(idx) = (sin(dtor*lon_save(idx))*x_coord(idx)  - &
                         cos(dtor*lon_save(idx))*y_coord(idx)) ! sin(azimuth) (no need to normalise)

      ENDIF
      
      IF (nidx > 0) THEN
        x_coord(idx) = atan2(lat_temp(idx), lon_temp(idx)+tiny(1.0_wp))/dtor
      ENDIF

      ROdata%Lev1b%Azimuth_tp(1:nLev) = x_coord(1:nLev)

      idx => WHERE ((ROdata%Lev1b%Azimuth_tp > ropp_MDTV) .AND. &
                    (ROdata%Lev1b%Azimuth_tp < 0), nidx)

      IF (nidx > 0) ROdata%Lev1b%Azimuth_tp(idx) = &
                    ROdata%Lev1b%Azimuth_tp(idx) + 360.0_wp

      ROdata%Lev1b%Lat_tp(1:nLev) = lat_save(1:nLev)
      ROdata%Lev1b%Lon_tp(1:nLev) = lon_save(1:nLev)

      idx => WHERE ((ROdata%Lev1b%Lat_tp < ropp_MDTV) .OR. &
                    (ROdata%Lev1b%Lon_tp < ropp_MDTV), nidx)

      IF (nidx > 0) ROdata%Lev1b%Azimuth_tp(idx) = ropp_MDFV

      nMSampLev = MIN ( nMSampLev, nSampLev )

      ROdata%Lev1b%Npoints = nMSampLev

      DEALLOCATE(x_coord, y_coord, lon_temp, lat_temp, lon_save, lat_save)

      DEALLOCATE ( Val, Lev )


    ELSE  ! No need to make a polar stereographic projection


      CALL message(msg_diag, "Thinning Level 1b Latitudes...")

      Lev(:) = Levtmp(:)
      Val(:) = ROdata%Lev1b%Lat_tp(1:nlev)

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                                 nThinLev,  ThinLev,  ThinVal,  &
                                 LocMethod, nSampLev )

      ROdata%Lev1b%Lat_tp(1:nLev) = Val(:)
      nMSampLev = MIN ( nMSampLev, nSampLev )

!--
      CALL message(msg_diag, "Thinning Level 1b Longitudes...")

      Lev(:) = Levtmp(:)
      Val(:) = ROdata%Lev1b%Lon_tp(1:nlev)

! transform Val if required from -180 to 180 to 0 to 360

      lontrans = .FALSE.
      DO i=2, nLev
         IF ((Val(i) .GT. ropp_MDTV) .AND. (Val(i-1) .GT. ropp_MDTV)) THEN
            IF (ABS(Val(i) - Val(i-1)) .GT. 300) lontrans = .TRUE.
         END IF
      END DO
      IF (lontrans) THEN
        idx => WHERE ( (Val .LT. 0) .AND. (Val .GT. ropp_MDTV), nidx)
        IF ( nidx > 0 ) Val(idx) = Val(idx) + 360.0_wp
      ENDIF

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                                 nThinLev,  ThinLev,  ThinVal,  &
                                 LocMethod, nSampLev )

! transformed from 0-360 to Val within -180 to +180
      idx => WHERE ( Val .GT. 180.0, nidx)
      IF ( nidx > 0 ) Val(idx) = Val(idx) - 360.0_wp

      ROdata%Lev1b%Lon_tp(1:nLev) = Val(:)

!--
      CALL message(msg_diag, "Thinning Level 1b Azimuths...")

      Lev(:) = Levtmp(:)
      Val(:) = ROdata%Lev1b%Azimuth_tp(1:nlev)

      CALL ropp_io_thin_select ( nlev,      Lev,      Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

      ROdata%Lev1b%Azimuth_tp(1:nLev) = Val(:)
      nMSampLev = MIN ( nMSampLev, nSampLev )

      ROdata%Lev1b%Npoints = nMSampLev

      DEALLOCATE ( Val, Lev )

      ! correct the impact altitudes back
      IF (iimpactalt) THEN
        idx => WHERE (ThinLev > ropp_MDTV, nidx) 
        IF (nidx > 0) ThinLev(idx) = ThinLev(idx) - ROdata%GEOref%undulation
      ENDIF

    ENDIF

  ELSE

    CALL message(msg_diag, "Level 1b data (BA) not thinned")

!-- output required number of thinned levels even if insufficient data to thin
    IF ( method /= "SAMPLE" ) THEN
      CALL ropp_io_init(ROdata%Lev1b, nThinLev)
      ROdata%Lev1b%Impact(1:nThinLev) = ThinLev
    ENDIF


  END IF

! 3. Thin refractivity and dry temperature as fn of geometric height (wrt geoid)
!===============================================================================

  nLev = ROdata%Lev2a%Npoints
  IF ( nLev > nThinLev ) THEN
    IF (ALLOCATED(Levtmp)) DEALLOCATE(Levtmp)
    ALLOCATE(Levtmp(nLev))
      Levtmp = ROdata%Lev2a%geop_refrac(1:nlev) + roc
    ALLOCATE ( Val(nLev),    &
               Lev(nLev),    &
               STAT=ierr )
    IF ( ierr /= 0 ) THEN
      CALL message(msg_fatal, "Failed to allocate memory for Level 2a arrays")
    ENDIF

    nSampLev  = nLev
    nMSampLev = nLev
!--
    CALL message(msg_diag, "Thinning Refractivity and Dry Temperature...")

! Actually we thin on Impact Parameter levels (ThinLev left from BA thinning)
! as this is independent of Level 1b profiles; converting ThinLev to
! geometric height here would need Impact Parameters which may not exist in
! the ROdata structure, or could be on different levels to Refractivity levels.
! After thinning, Impact Parameter levels (still in ThinLev) can then be
! converted to geometric height using the thinned refractivites.
!
! Convert full refractivity levels in geometric ht (wrt geoid) to IP using:
!  1) h = r - Rc - U
!  2) a = n.r
!  3) n = 1+N.10^-6
!  => a = (1+N.10^-6).(h+Rc+U)
!
! If N is invalid (which may be the case at high and low fixed levels),
! use a log interpolation between height of last good value (initialised
! with N=10^-6 at highest level) and a reference value (N=300) at the surface
! (h=0). This is just to give sensible impact parameter values for all of
! the fixed levels which can be inverted back to geometric heights later;
! they do not need to be accurate, since by definition, the refractvity is
! missing so these levels won't be useful anyway.

    IF ( fixlev ) THEN
      IF ( Levtmp(nLev) > Levtmp(1) ) THEN    ! ascending fixed levels
        i1 = nLev
        i2 =  1
        i3 = -1
      ELSE                                    ! descending fixed levels
        i1 = 1
        i2 = nLev
        i3 = 1
      ENDIF

      ! Log of reference refractivity at surface
      N300  = LOG(300.0_wp)
      ! Initialise with a small N at highest level
      Ngood = LOG(1e-6_wp)
      Hgood = Levtmp(i1) - roc

      ! Undulation value (if valid)
      IF ( ROdata%GEOref%Undulation > ropp_MDTV ) THEN

        undulation = ROdata%GEOref%Undulation

      ELSE

        undulation = 0.0_wp
        CALL message(msg_warn, "Error reading undulation, using " // &
           "default value 0 km for computing output thinned levels")

      ENDIF

      DO i = i1, i2, i3                       ! always scan downwards
        IF ( ROdata%Lev2a%Refrac(i) > 0.0_wp ) THEN ! use good N
          N     = ROdata%Lev2a%Refrac(i)
          Ngood = LOG(N)
          Hgood = Levtmp(i) - roc
        ELSE                                   ! interpolate N at H
          IF ( Levtmp(i) > ropp_MDTV ) THEN
            H = Levtmp(i) - roc
            N = EXP ( Ngood + ( Hgood - H ) * (N300 - Ngood) / Hgood )
            N = max(N, 1e-60_wp)
          ELSE
            N = 0.0_wp
          ENDIF
        ENDIF

        Lev(i) = ( 1.0_wp + 1e-6_wp * N ) &
               * ( ROdata%Lev2a%Alt_refrac(i) + roc + undulation )
      END DO
    ELSE
      Lev = ROdata%Lev2a%Alt_refrac
    END IF

!--
    CALL message(msg_diag, "Thinning Refractivity Geopotential Heights...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%geop_refrac(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev2a%geop_refrac(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Refractivity values...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%Refrac(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Levtmp,    Val,      &
                               nThinLev, ThinLev,   ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev2a%Refrac(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Refractivity errors...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%Refrac_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, .TRUE. )

    ROdata%Lev2a%Refrac_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Refractivity quality...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%Refrac_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev2a%Refrac_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 3.1 Dry temperature
!--------------------

    CALL message(msg_diag, "Thinning Dry Temperature values...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%Dry_Temp(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Levtmp,    Val,      &
                               nThinLev, ThinLev,   ThinVal,  &
                               LocMethod,   nSampLev )

    ROdata%Lev2a%Dry_Temp(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Dry Temperature errors...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%Dry_Temp_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, .TRUE. )

    ROdata%Lev2a%Dry_Temp_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )
!--
    CALL message(msg_diag, "Thinning Dry Temperature quality...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2a%Dry_Temp_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev2a%Dry_Temp_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! Convert thinned Impact Parameter levels (still in ThinLev for fixed
! level interpolation) back to geometric height (wrt geoid)
! using the just thinned refractivity values:
!  1) h = r - Rc - u
!  2) a = r . n
!  3) n = ( 1 + N.10^6 )
!  => h = a/(1+N.10^6) - Rc - u
! If N not valid (e.g. occulation, does not span range of fixed levels)
! use a log interpolation assuming N=300 at h=0 and the last good value.
! NB: this is only to provide realistic height values when thinned
! refractivity is missing; the interpolated N is used for this purpose ONLY
! - the returned N values will still be 'missing'.

    IF ( fixlev ) THEN
      N300  = LOG(300.0_wp) ! Log of reference refractivity
      Ngood = LOG(1e-6_wp)
      Hgood = ThinLev(nMSampLev) - roc
      DO i = nMSampLev, 1, -1
        IF ( ROdata%Lev2a%refrac(i) > 0.0_wp ) THEN
          N = ROdata%Lev2a%refrac(i)
          Ngood = LOG(N)
          Hgood = ThinLev(i) - roc
        ELSE
          H = ThinLev(i) - roc
          N = EXP ( Ngood + ( Hgood - H ) * ( N300 - Ngood ) / Hgood )
        ENDIF
        ThinLev(i) = ( ThinLev(i) / ( 1.0_wp + 1e-6_wp * N ) ) &
                   - roc - undulation
      END DO
    END IF
    ROdata%Lev2a%Alt_refrac(1:nMSampLev)       = ThinLev(1:nMSampLev)
    ROdata%Lev2a%Alt_refrac(nMSampLev+1:nLev)  = ropp_MDFV

    ROdata%Lev2a%Npoints = nMSampLev

    DEALLOCATE ( Val, Lev, Levtmp )

    refrac_thinned = .TRUE.   ! OK to thin T,q,p if needed

  ELSE
    CALL message(msg_diag, "Level 2a data not thinned")

!-- output required number of thinned levels even if insufficient data to thin
    IF ( method /= "SAMPLE" ) THEN
      CALL ropp_io_init(ROdata%Lev2a, nThinLev)
      ROdata%Lev2a%Alt_refrac(1:nThinLev) = ThinLev(:) - roc
    ENDIF

  END IF

! 4. Thin T,q,P as fn geopotential height (wrt geoid)
!====================================================

! WARNING: If refractivities were not thinned (because there were no
! Level 2a data present, or there were not enough full points to be
! thinned) we have a potential problem as we can't convert ThinLev
! Impact Parameters to geometric heights and hence to geopotential heights.
! In practice, if Level 2a profile has too few levels to thin, so
! will the Level 2b profile. In the case of ROdata containing *only*
! Level 2b data and enough points to require thinning, there will
! be an irreconcilable inconsistency between full (geopotential ht)
! and thinned levels (still in IP) - so best not to thin at all in
! this case.

  nLev = ROdata%Lev2b%Npoints
  IF ( nLev > nThinLev .AND. refrac_thinned ) THEN
    ALLOCATE ( Val(nLev),       &
               Lev(nLev),       &
               Levtmp(nLev),    &
               STAT=ierr )
    IF ( ierr /= 0 ) THEN
      CALL message(msg_fatal, "Failed to allocate memory for Level 2b arrays")
    END IF

    nSampLev  = nLev
    nMSampLev = nLev

    Lev(:) = ROdata%Lev2b%geop(1:nLev)
    Levtmp(:) = Lev(:)

! Starting from ThinLev in geometric height (wrt geoid)
! calculated during refractivity thinning in (3) above,
! convert to geopotential height (wrt geoid).

    ThinLev = geometric2geopotential ( ROdata%GEOref%Lat, ThinLev )

! 4.0.1 Geopotentials (always set)
!--------------------

    ROdata%Lev2b%geop(1:nThinLev)      = ThinLev(1:nThinLev)
    ROdata%Lev2b%geop(nThinLev+1:nLev) = ropp_MDFV

! 4.1 Temperature
!-----------------

    CALL message(msg_diag, "Thinning Temperature...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2b%Temp(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,    Val,      &
                               nThinLev,  ThinLev,   ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev2b%Temp(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 4.1.1 Temperature Errors (not always available)
!------------------------------------------------

    CALL message(msg_diag, "Thinning Temperature errors...")

    Levtmp(:) = Lev(:)
    Val(:)    = ROdata%Lev2b%Temp_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, sigma )

    ROdata%Lev2b%Temp_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )


! 4.2 Specific humidity
!-----------------------

    CALL message(msg_diag, "Thinning Humidity...")

    Levtmp(:) = Lev(:)
    Val(:)    = ROdata%Lev2b%SHum(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Levtmp,   Val,      &
                               nThinLev, ThinLev,  ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev2b%SHum(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev)

! 4.2.1 Specific humidity error (not always available)
!-----------------------------------------------------

    CALL message(msg_diag, "Thinning Humidity errors...")

    Levtmp(:) = Lev(:)
    Val(:)    = ROdata%Lev2b%SHum_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Levtmp,   Val,      &
                               nThinLev, ThinLev,  ThinVal,  &
                               Method,   nSampLev, sigma )

    ROdata%Lev2b%SHum_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 4.3 Pressure
!--------------

    CALL message(msg_diag, "Thinning Pressure...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2b%Press(1:nlev)

    CALL ropp_io_thin_select ( nlev,     Levtmp,    Val,      &
                               nThinLev, ThinLev,   ThinVal,  &
                               Method,   nSampLev )

    ROdata%Lev2b%Press(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev)

! 4.3.1 Pressure errors (not always avaiable)
!--------------------------------------------

    CALL message(msg_diag, "Thinning Pressure errors...")

    Levtmp(:) = Lev(:)
    Val(:) = ROdata%Lev2b%Press_sigma(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,   Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev, sigma )

    ROdata%Lev2b%Press_sigma(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

! 4.4 Meteo quality
!------------------

    CALL message(msg_diag, "Thinning Meteo quality...")

    Levtmp(:) = Lev(:)
    Val(:)    = ROdata%Lev2b%Meteo_qual(1:nlev)

    CALL ropp_io_thin_select ( nlev,      Levtmp,    Val,      &
                               nThinLev,  ThinLev,  ThinVal,  &
                               LocMethod, nSampLev )

    ROdata%Lev2b%Meteo_qual(1:nLev) = Val(:)
    nMSampLev = MIN ( nMSampLev, nSampLev )

    ROdata%Lev2b%Npoints = nMSampLev

    DEALLOCATE ( Val, Lev, Levtmp )

  ELSE
    IF ( refrac_thinned .AND. nLev > nThinLev) THEN
      CALL message(msg_diag, "Could not thin Level 2b as Level 2a data" // &
                             " not present or not thinned" )
    ELSE
      CALL message(msg_diag, "Level 2b data (T,q,P) not thinned")

 !-- output required number of thinned levels even if insufficient data to thin
      IF ( method /= "SAMPLE" ) THEN
        CALL ropp_io_init(ROdata%Lev2b, nThinLev)
        ROdata%Lev2b%geop(1:nThinLev) = ThinLev(:)
      ENDIF

    END IF

  END IF

! 5. Tidy up
!============

  IF ( ASSOCIATED ( idx ) ) DEALLOCATE ( idx )
  DEALLOCATE ( ThinVal, ThinLev )

  IF (iranchk) CALL ropp_io_rangecheck ( ROdata )

  CALL message_set_routine(routine)

END SUBROUTINE ropp_io_thin
