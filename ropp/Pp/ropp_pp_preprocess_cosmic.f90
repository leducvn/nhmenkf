! $Id: ropp_pp_preprocess_COSMIC.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_preprocess_COSMIC(ro_data, config, LCF)

!****s* Preprocessing/ropp_pp_preprocess_COSMIC *
!
! NAME
!    ropp_pp_preprocess_COSMIC - Mission-specific Level1a data preprocessing
!                                for COSMIC 
!                   
! SYNOPSIS
!    call ropp_pp_preprocess_COSMIC(ro_data, config, LCF)
! 
! DESCRIPTION
!    Removal of phase discontinuity - updates phase_L1 and phase_L2
!
! INPUTS
!    type(ROprof)   :: ro_data      ! Radio occultation data strucuture
!    type(PPConfig) :: config       ! Configuration options
!    integer        :: LCF          ! Lost carrier flag
!
! OUTPUT
!    type(ROprof)   :: ro_data      ! Corrected radio occultation data
!    type(PPConfig) :: config       ! Configuration options
!    integer        :: LCF          ! Lost carrier flag
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
  USE ropp_utils
  USE ropp_io_types, ONLY: ROprof
! USE ropp_pp_preproc, not_this => ropp_pp_preprocess_COSMIC
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPConfig
  
  IMPLICIT NONE

  TYPE(ROprof),   INTENT(inout) :: ro_data      ! Radio occultation data struct
  TYPE(PPconfig), INTENT(inout) :: config       ! Configuration options
  INTEGER, DIMENSION(:), INTENT(inout) :: LCF   ! lost carrier flag

  INTEGER                       :: i, ic, igb   ! Index counters   
  INTEGER                       :: n            ! Number of data points
  INTEGER                       :: iol          ! Open loop discontinuity index
  INTEGER,  DIMENSION(:), ALLOCATABLE :: gb     ! Navigation bits
  REAL(wp), DIMENSION(:), ALLOCATABLE :: tg     ! Navigation bit GPS times
  REAL(wp), DIMENSION(:), ALLOCATABLE :: phase  ! phase data array
  REAL(wp), DIMENSION(:), ALLOCATABLE :: tmp   
  REAL(wp)                            :: sb     ! Phase base
  REAL(wp)                            :: start_time ! Occultation start time
  REAL(wp)                            :: s0     ! Straight-line path
  REAL(wp)                            :: tgl    ! Observation GPS times
  REAL(wp)                            :: tg0, tgN   ! Time limits
  INTEGER                             :: ngb        ! No. nav bits
  INTEGER                             :: unit
  INTEGER, DIMENSION(8)               :: DT8     ! date/time array
  LOGICAL                             :: exist

!-------------------------------------------------------------------------------
! 2. Remove phase discontinuity
!-------------------------------------------------------------------------------

  n = ro_data%Lev1a%Npoints
  ALLOCATE(phase(n))

  DO ic = 1, 2
     
    IF (ic == 1) phase(:) = ro_data%Lev1a%phase_L1(:)
    IF (ic == 2) phase(:) = ro_data%Lev1a%phase_L2(:)
      
    iol = SUM(MAXLOC(ABS(phase(2:n) - phase(1:n-1)))) 

    IF (iol > 1) THEN
     
       sb = phase(iol+1) - (phase(iol) + (phase(iol) - phase(iol-1)))
       phase(iol+1:n) = phase(iol+1:n) - sb
        
    ENDIF
     
    sb       = MINVAL(phase(:))
    phase(:) = phase(:) - sb
     
    IF (ic == 1) ro_data%Lev1a%phase_L1(:) = phase(:)
    IF (ic == 2) ro_data%Lev1a%phase_L2(:) = phase(:)

  ENDDO
  DEALLOCATE(phase)
 
!-------------------------------------------------------------------------------
! 3. Retrieve lost carrier flag information from input data file
!-------------------------------------------------------------------------------

    IF (ASSOCIATED(ro_data%vlist%VlistD1d)) THEN
      CALL message(msg_info, 'Reading lost carrier flag data from input file')
      ALLOCATE(tmp(n))
      tmp = ABS(ro_data%vlist%VlistD1d%data(1:n) + 999.0_wp)
      WHERE(tmp(:) > 1e-2_wp)
        LCF(:) = IBSET(LCF(:),0)
      ENDWHERE
      DEALLOCATE(tmp)
      DEALLOCATE(ro_data%vlist%VlistD1d)
    ELSE
      CALL message(msg_info, 'No lost carrier flag data in input file')
    ENDIF

!-------------------------------------------------------------------------------
! 4. Read external navigation bits (if available)
!-------------------------------------------------------------------------------
     
    ! 4.1 Read navigation bit file

    ngb = 0
    IF (config%navbit_file /= " ") THEN
      INQUIRE(file = config%navbit_file, exist=exist)
       
      IF (exist) THEN 
        unit = get_io_unit()
        OPEN(unit, file = config%navbit_file)
        CALL message(msg_info, 'Reading external NDM data from file ' //  &
                                   config%navbit_file)
         
        READ(unit, *) tg0, tgN, ngb
        IF (tg0 == -999.0 .or. tgN == -999.0 .or. ngb == 0) THEN 
          ngb = 0
          CLOSE(unit)
        ELSE
          ALLOCATE(gb(ngb))
          ALLOCATE(tg(ngb))
          DO i=1, ngb
            READ(unit, *) gb(i)
            tg(i) = tg0 + REAL((i-1),wp)*(tgN - tg0)/REAL(ngb,wp)
          ENDDO
          CLOSE(unit)
        ENDIF
         
      ELSE
        CALL message(msg_info, 'No external NDM data read')
      ENDIF
       
    ELSE
      CALL message(msg_info, 'No external NDM data read')
    ENDIF
     
    IF (ngb > 0) THEN
       
      ! 4.2 GPS transmitting time
      DT8 = (/ro_data%dtocc%year,   ro_data%dtocc%month,  &
              ro_data%dtocc%day,    0,                    &
              ro_data%dtocc%hour,   ro_data%dtocc%minute, &
              ro_data%dtocc%second, ro_data%dtocc%msec/)
      CALL TimeSince ( DT8, start_time, 1, Base="GPSSEC" )
     
      ! 4.3 Interpolate navigation bits to observation GPS times
      DO i=1,n
        s0 = SQRT(SUM((ro_data%Lev1a%r_leo(i,:) - ro_data%Lev1a%r_gns(i,:))**2))
        tgl = start_time + ro_data%Lev1a%dtime(i) -    &
                (ro_data%Lev1a%phase_L1(i) + s0)/c_light
        
        igb = FLOOR(50*(tgl - tg(1))) + 1
        IF ((igb >= 1) .and. (igb <= size(tg))) THEN
          IF (gb(igb)/10 == 1) THEN
            LCF(i)  = IBSET(LCF(i),1)
          ENDIF
          IF (MODULO(gb(igb),10) == 1) THEN
            LCF(i) = IBSET(LCF(i),2)
          ENDIF
        ENDIF
      ENDDO
       
    ENDIF
     
    IF (ALLOCATED(tg)) DEALLOCATE(tg)
    IF (ALLOCATED(gb)) DEALLOCATE(gb)

!-------------------------------------------------------------------------------
! 5. Scaling amplitude
!-------------------------------------------------------------------------------
     
    ro_data%Lev1a%snr_L1ca = ro_data%Lev1a%snr_L1ca * 0.1_wp
    ro_data%Lev1a%snr_L1p  = ro_data%Lev1a%snr_L1p  * 0.1_wp
    ro_data%Lev1a%snr_L2p  = ro_data%Lev1a%snr_L2p  * 0.1_wp

END SUBROUTINE ropp_pp_preprocess_COSMIC
