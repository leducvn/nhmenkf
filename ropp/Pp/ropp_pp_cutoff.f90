! $Id: ropp_pp_cutoff.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_cutoff(ro_data, config, var1, var2, var3)

!****s* Preprocessing/ropp_pp_cutoff *
!
! NAME
!    ropp_pp_cutoff - Cut off occultation data based on amplitude, smoothed
!                     bending angle and impact parameter values
!
! SYNOPSIS
!    call ropp_pp_cutoff(ro_data, config, var1, var2, var3)
!
! DESCRIPTION
!    1) Compute smoothed bending angle profile
!    2) Cut off from amplitude (based on config%Acut parameter)
!    3) Cut off from bending angle profile (based on config%Bcut, config%Pcut)
!
! INPUTS
!    type(ROprof)                      :: ro_data ! RO data strucuture
!    type(PPConfig)                    :: config  ! Configuration options
!    real(wp), dimension(:), optional  :: var1    ! Extra variables to shrink
!    real(wp), dimension(:), optional  :: var2    ! Extra variables to shrink
!    integer,  dimension(:), optional  :: var3    ! Extra variables to shrink
!
! OUTPUT
!    type(ROprof)                      :: ro_data ! Shrunk RO data strucuture
!    type(PPConfig)                    :: config  ! Configuration options
!    real(wp), dimension(:), optional  :: var1    ! Shrunk extra variables
!    real(wp), dimension(:), optional  :: var2    ! Shrunk extra variables
!    integer,  dimension(:), optional  :: var3    ! Shrunk extra variables
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
  USE ropp_io, ONLY: ropp_io_init, &
                     ropp_io_shrink, ropp_io_free
! USE ropp_pp_preproc, not_this => ropp_pp_cutoff
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPconfig

  IMPLICIT NONE

  TYPE(ROprof),       INTENT(inout) :: ro_data  ! RO data strucutre
  TYPE(PPconfig),     INTENT(inout) :: config   ! Configuration options
  REAL(wp), DIMENSION(:), POINTER, OPTIONAL :: var1
  REAL(wp), DIMENSION(:), POINTER, OPTIONAL :: var2
  INTEGER,  DIMENSION(:), POINTER, OPTIONAL :: var3

  REAL(wp), DIMENSION(:), ALLOCATABLE :: impact_ht_L1
  REAL(wp), DIMENSION(:), ALLOCATABLE :: slta
  REAL(wp)                            :: sr0
  INTEGER                             :: stride
  INTEGER                             :: w_smooth
  INTEGER                             :: i, n
  INTEGER                             :: jmin, jmax
  INTEGER                             :: ocd ! Occultation direction (1=rising)
  INTEGER                             :: ocd1
  CHARACTER(len=8)                    :: imin_str, jmax_str, jmin_str
  CHARACTER(len = 256)                :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_cutoff')

!-------------------------------------------------------------------------------
! 2. Downsampling
!-------------------------------------------------------------------------------

  n = ro_data%Lev1a%Npoints
  sr0 = 1.0_wp/ABS(MINVAL(ro_data%Lev1a%dtime(2:n)-ro_data%Lev1a%dtime(1:n-1)))
!  stride = MAX(1, NINT(sr0/50.0_wp))
  stride = 1

  IF (stride > 1) THEN
    WRITE(imin_str,  '(f8.0)') sr0
    CALL message(msg_diag, "Downsampling from " // imin_str // " to 50.0 Hz")

    CALL ropp_io_shrink(ro_data%Lev1a, 1, n, stride)
    IF (PRESENT(var1)) CALL shrink_var(var1, 1, n, stride)
    IF (PRESENT(var2)) CALL shrink_var(var2, 1, n, stride)
    IF (PRESENT(var3)) CALL shrink_varint(var3, 1, n, stride)

    WRITE(imin_str,  '(i8)') ro_data%Lev1a%Npoints
    CALL message(msg_diag, imin_str // " downsampled data points")

  ENDIF

  CALL ropp_io_init(ro_data%Lev1b, n)

  ! 2.1 Smoothing window determination

  DO i=1,n
     ro_data%Lev1b%impact(i) =                                           &
        impact_parameter( ro_data%Lev1a%r_leo(i,:)-ro_data%georef%r_coc, &
                          ro_data%Lev1a%r_gns(i,:)-ro_data%georef%r_coc )
  ENDDO

  ocd = NINT(SIGN(1.0_wp, ro_data%Lev1b%impact(n)-ro_data%Lev1b%impact(1)))

  config%Pmax = MAXVAL(ro_data%Lev1b%impact(:))
  config%Pmin = MAX(MINVAL(ro_data%Lev1b%impact(:)), &
                    ro_data%georef%roc+2000.0_wp)

  w_smooth = CEILING( config%fw_go_smooth*(n-1) / ABS(config%Pmax-config%Pmin) )

!-------------------------------------------------------------------------------
! 2. Compute smoothed bending angle from phase - fill ro_data%Lev1b structure
!-------------------------------------------------------------------------------

  n = ro_data%Lev1a%Npoints


  ! 2.2 Calculate smoothed bending angle profile

  CALL ropp_pp_bending_angle_go(ro_data%Lev1a%dtime, ro_data%Lev1a%r_leo,  &
                                ro_data%Lev1a%r_gns, ro_data%georef%r_coc, &
                                ro_data%Lev1a%phase_L1,                    &
                                ro_data%Lev1a%phase_L2,                    &
                                w_smooth,                                  &
                                config%filter_method,                      &
                                ro_data%Lev1b%impact_L1,                   &
                                ro_data%Lev1b%bangle_L1,                   &
                                ro_data%Lev1b%impact_L2,                   &
                                ro_data%Lev1b%bangle_L2  )

!-------------------------------------------------------------------------------
! 4. Cut-off data from smoothed bending angle profile
!-------------------------------------------------------------------------------

  ! 4.1 Calculate impact height

  ALLOCATE(impact_ht_L1(n))
  ALLOCATE(slta(n))

  ocd1 = NINT(SIGN(1.0_wp, ro_data%Lev1b%impact_L1(n)-ro_data%Lev1b%impact_L1(1)))
  IF ( ocd < ocd1 ) ro_data%Lev1b%impact_L1(n) = ro_data%Lev1b%impact(n)
  IF ( ocd > ocd1 ) ro_data%Lev1b%impact_L1(1) = ro_data%Lev1b%impact(1)

  CALL ropp_pp_monotonous(ro_data%Lev1b%impact_L1, -1)

  impact_ht_L1(:) = ro_data%Lev1b%impact_L1 - ro_data%georef%roc
  slta(:) =  ro_data%Lev1b%impact - ro_data%georef%roc

  ! 4.2 Index limit determination

  Jmin = SUM(MAXLOC(ro_data%Lev1b%impact_L1, &
                    Mask = ((ro_data%Lev1b%bangle_L1 < config%Bcut) .AND.  &
                    (impact_ht_L1 > config%Pcut) .AND.                     &
                    (slta > config%Hcut))))

  Jmax = SUM(MINLOC(ro_data%Lev1b%impact_L1, &
                    Mask = ((ro_data%Lev1b%bangle_L1 < config%Bcut) .AND.  &
                    (impact_ht_L1 > config%Pcut) .AND.                     &
                    (slta > config%Hcut))))

  IF ((Jmin < 1) .OR. (Jmin > N) .OR. &
       (Jmax < 1) .OR. (Jmax > N)) THEN
     Jmin = 1
     Jmax = N
     RETURN
  END IF

  IF (Jmax < Jmin) THEN
     i    = Jmax
     Jmax = Jmin
     Jmin = i
  END IF

  ! 4.3 Select data subset

  IF (jmin > 1 .OR. jmax < n) THEN
     WRITE(jmin_str,  '(i8)') jmin
     WRITE(jmax_str,  '(i8)') jmax
     CALL message(msg_info,"Cut-off (bangle/impact criterion). Keep data " //  &
                            jmin_str  // " to " // jmax_str)

     CALL ropp_io_shrink(ro_data%Lev1a, jmin, jmax, 1)
     n = ro_data%Lev1a%Npoints

     IF (PRESENT(var1)) CALL shrink_var(var1, jmin, jmax, 1)
     IF (PRESENT(var2)) CALL shrink_var(var2, jmin, jmax, 1)
     IF (PRESENT(var3)) CALL shrink_varint(var3, jmin, jmax, 1)

   ENDIF

!-------------------------------------------------------------------------------
! 5. Clean up
!-------------------------------------------------------------------------------

  CALL ropp_io_free(ro_data%Lev1b)
  DEALLOCATE(impact_ht_L1)
  DEALLOCATE(slta)

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

    INTEGER, DIMENSION(:), POINTER       :: var
    INTEGER,               INTENT(in)    :: imin
    INTEGER,               INTENT(in)    :: imax
    INTEGER,               INTENT(in)    :: stride
    INTEGER, DIMENSION(:), POINTER       :: tmp => null()

    CALL copy_alloc(var(imin:imax:stride), tmp)
    DEALLOCATE(var)

    CALL copy_alloc(tmp, var)
    DEALLOCATE(tmp)

  END SUBROUTINE shrink_varint


END SUBROUTINE ropp_pp_cutoff
