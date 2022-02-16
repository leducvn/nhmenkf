! $Id: ropp_pp_cutoff_amplitude.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_cutoff_amplitude(ro_data, LCF, config)

!****s* Preprocessing/ropp_pp_cutoff_amplitude *
!
! NAME
!    ropp_pp_cutoff_amplitude - Cut off occultation data based on amplitude
!                               and missing data flag
!
! SYNOPSIS
!    call ropp_pp_cutoff_amplitude(ro_data, LCF, config)
!
! DESCRIPTION
!    Cut off from amplitude (based on config%Acut parameter)
!
! INPUTS
!    type(ROprof)                      :: ro_data ! RO data strucuture
!    integer,  dimension(:)            :: LCF     ! Lost carrier flag
!    type(PPConfig)                    :: config  ! Configuration options

!
! OUTPUT
!    type(ROprof)                      :: ro_data ! Shrunk RO data strucuture
!    integer,  dimension(:)            :: LCF     ! Shrunk LCF
!    type(PPConfig)                    :: config  ! Configuration options
!
! NOTES
!   Requires ROprof data structure type, defined in ropp_io module. This
!   routine therefore requires that the ropp_io module is pre-installed before
!   compilation.
!
! REFERENCES
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

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
  USE ropp_utils, ONLY: impact_parameter
  USE messages
  USE ropp_io_types, ONLY: ROprof
  USE ropp_io, ONLY: ropp_io_shrink
! USE ropp_pp_preproc, not_this => ropp_pp_cutoff_amplitude
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPconfig

  IMPLICIT NONE

  TYPE(ROprof),       INTENT(inout) :: ro_data  ! RO data strucutre
  INTEGER,  DIMENSION(:), POINTER   :: LCF      ! Lost carrier flag
  TYPE(PPconfig),     INTENT(inout) :: config   ! Configuration options

  REAL(wp)                            :: p1, pN
  INTEGER                             :: i, n
  INTEGER                             :: imax, imin
  INTEGER                             :: w_smooth
  INTEGER                             :: ocd ! Occultation direction (1=rising)
  CHARACTER(len=8)                    :: imin_str, imax_str
  CHARACTER(len = 256)                :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_cutoff_amplitude')

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  n = ro_data%Lev1a%Npoints

  p1 = impact_parameter( ro_data%Lev1a%r_leo(1,:)-ro_data%georef%r_coc, &
                         ro_data%Lev1a%r_gns(1,:)-ro_data%georef%r_coc )
  pN = impact_parameter( ro_data%Lev1a%r_leo(n,:)-ro_data%georef%r_coc, &
                         ro_data%Lev1a%r_gns(n,:)-ro_data%georef%r_coc )

  ocd = NINT(SIGN(1.0_wp, pN-p1))

  config%Pmax = MAX(p1, pN)
  config%Pmin = MIN(p1, pN)

  w_smooth = CEILING( config%fw_go_smooth*(n-1)/ABS(config%Pmax-config%Pmin) )

!-------------------------------------------------------------------------------
! 3. Cut-off data from amplitude (snr_L1p) and missing data flag (LCF)
!-------------------------------------------------------------------------------

  ! 3.1 Determine cut-off limits

  imin = 1
  imax = n

  IF ( ocd < 0 ) THEN     ! setting occultation

     imin = 1
     DO i=n,1,-1
        IF (ro_data%Lev1a%snr_L1ca(i) >                               &
                 MAXVAL(ro_data%Lev1a%snr_L1ca(:))*config%Acut) THEN
           imax = MIN(i + 2*w_smooth, n)
           EXIT
        ENDIF
     ENDDO
     DO i=1,n
       IF (BTEST(LCF(i),3)) THEN
         imax = MAX(1, MIN(i-1, imax))
         EXIT
       ENDIF
     ENDDO

  ELSE                    ! rising occultation

     imax = n
     DO i=1,n
        IF (ro_data%Lev1a%snr_L1ca(i) >                               &
                 MAXVAL(ro_data%Lev1a%snr_L1ca(:))*config%Acut) THEN
           imin = MAX(i - 2*w_smooth, 1)
           EXIT
        ENDIF
     ENDDO
     DO i=n,1,-1
       IF (BTEST(LCF(i),3)) THEN
         imin = MIN(n, MAX(i+1, imin))
         EXIT
       ENDIF
     ENDDO

  ENDIF

  ! 3.2 Select data subset

  IF ((imin > 1 .OR. imax < n) .AND. (imax >= imin)) THEN
     WRITE(imin_str,  '(i8)') imin
     WRITE(imax_str,  '(i8)') imax
     CALL message(msg_info,"Cut-off (amplitude/LCF criterion). Keep data " // &
                           imin_str  // " to " // imax_str)

     CALL ropp_io_shrink(ro_data%Lev1a, imin, imax, 1)
     CALL shrink_varint(lcf, imin, imax, 1)
     n = ro_data%Lev1a%Npoints

  ENDIF

  CALL message_set_routine(routine)

CONTAINS

!-------------------------------------------------------------------------------
! 6. Select data subset for additional variables
!-------------------------------------------------------------------------------

  SUBROUTINE shrink_var(var, imin, imax, stride)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_utils, ONLY: copy_alloc

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), POINTER       :: var
    INTEGER,                INTENT(in)    :: imin
    INTEGER,                INTENT(in)    :: imax
    INTEGER,                INTENT(in)    :: stride
    REAL(wp), DIMENSION(:), POINTER       :: tmp => null()

    CALL copy_alloc(var(imin:imax:stride), tmp)
    DEALLOCATE(var)

    CALL copy_alloc(tmp, var)
    DEALLOCATE(tmp)

  END SUBROUTINE shrink_var

  SUBROUTINE shrink_varint(var, imin, imax, stride)

    USE ropp_utils, ONLY: copy_alloc

    IMPLICIT NONE

    INTEGER, DIMENSION(:), POINTER      :: var
    INTEGER,               INTENT(in)   :: imin
    INTEGER,               INTENT(in)   :: imax
    INTEGER,               INTENT(in)   :: stride
    INTEGER, DIMENSION(:), POINTER      :: tmp => null()

    CALL copy_alloc(var(imin:imax:stride), tmp)
    DEALLOCATE(var)

    CALL copy_alloc(tmp, var)
    DEALLOCATE(tmp)

  END SUBROUTINE shrink_varint


END SUBROUTINE ropp_pp_cutoff_amplitude
