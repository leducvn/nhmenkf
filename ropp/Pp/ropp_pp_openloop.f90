! $Id: ropp_pp_openloop.f90 2021 2009-01-16 10:49:04Z frhl $

SUBROUTINE ropp_pp_openloop(time, phase_L1, phase_L2, phase_LM,    &
                            r_leo, r_gns, r_coc, LCF)

!****s* Preprocessing/ropp_pp_openloop *
!
! NAME
!    ropp_pp_openloop - Transform open-loop data into excess phase
!                   
! SYNOPSIS
!    call ropp_pp_openloop(time, phase_L1, phase_L2, phase_LM, 
!                             r_leo, r_gns, r_coc, LCF)
! 
! DESCRIPTION
!    This routine transforms open-loop data into excess phase.
!      1. Frequency down-conversion by subtracting phase mode.
!      2. Navigation bits removal.
!      3. Phase accumulation.
!      4. Restoring phase variation.
!
! INPUTS
!    real(wp), dimension(:)   :: time      ! time of samples (s)
!    real(wp), dimension(:)   :: phase_L1  ! excess phase L1 (m)
!    real(wp), dimension(:)   :: phase_L2  ! excess phase L2 (m)
!    real(wp), dimension(:)   :: phase_LM  ! model excess phase (m)
!    real(wp), dimension(:,:) :: r_leo     ! cartesian LEO coordinates (ECF)
!    real(wp), dimension(:,:) :: r_gns     ! cartesian GPS coordinates (ECF)
!    real(wp), dimension(:)   :: r_coc     ! cartesian centre curvature (ECF)
!    integer,  dimension(:)   :: lcf       ! lost carrier flag
!
! OUTPUT
!    real(wp), dimension(:)   :: phase_L1  ! excess phase L1 (m)
!    real(wp), dimension(:)   :: phase_L2  ! excess phase L2 (m)
!    integer,  dimension(:)   :: lcf       ! lost carrier flag
!
! NOTES
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
! USE ropp_pp, not_this => ropp_pp_openloop
  USE ropp_pp_constants, ONLY: pi, f_L1, f_L2, c_light
  USE ropp_utils, ONLY: impact_parameter
  USE messages

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(in)    :: time      ! time of samples (s)
  REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L1  ! excess phase L1 (m)
  REAL(wp), DIMENSION(:), INTENT(inout) :: phase_L2  ! excess phase L2 (m)
  REAL(wp), DIMENSION(:), INTENT(in)    :: phase_LM  ! model excess phase (m)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_leo     ! LEO coordinates (ECF)
  REAL(wp), DIMENSION(:,:), INTENT(in)  :: r_gns     ! GPS coordinates (ECF)
  REAL(wp), DIMENSION(:),   INTENT(in)  :: r_coc     ! centre curvature (ECF)
  INTEGER,  DIMENSION(:), INTENT(inout) :: LCF       ! Lost carrier flag

  INTEGER, PARAMETER  :: lb = 10  ! Maximum lag for navbits correlation
  INTEGER             :: n        ! Number of data
  INTEGER             :: i        ! Work index
  INTEGER             :: ic       ! Channel index
  REAL(wp)            :: SB       ! Phase base
  REAL(wp)            :: PS1      ! Start straight-line impact (m)
  REAL(wp)            :: PSn      ! End straight-line impact (m)
  INTEGER             :: imin     ! Start of open loop
  INTEGER             :: imax     ! End of open loop
! INTEGER             :: iol      ! Open loop discontinuity index ! Commented at 21 July, 2016
  INTEGER             :: kmin     ! Start of correlation interval
  INTEGER             :: kmax     ! End of correlation interval
  INTEGER             :: dk       ! Correlation interval width
  INTEGER             :: ni       ! Sum of squared vavbits variations
  INTEGER             :: inb      ! Correlation maximum position
  REAL(wp)            :: CNBmax   ! Correlation maximum
  INTEGER             :: bit_NB   ! Navbit position in LCF
  INTEGER             :: bit_NBQ  ! Navbit quality position in LCF
  INTEGER             :: ocd      ! Occultation direction:
                                  !   -1 - setting
                                  !    1 - rising
  
  REAL(wp), DIMENSION(2)                :: k   ! wave vectors
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: TG  ! GPS time (s)
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: DS  ! k(S - S_mod)
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: DDS ! OL delta phase difference
  INTEGER,  DIMENSION(:),   ALLOCATABLE :: DNB ! Navbits phase difference
  REAL(wp), DIMENSION(:),   ALLOCATABLE :: CNB ! Phase-navbits correlation

  CHARACTER(len=34) :: outstr
  LOGICAL           :: restore_pll = .TRUE.

  CHARACTER(len = 256)                :: routine

  CALL message_get_routine(routine)
  CALL message_set_routine('ropp_pp_openloop')

!-------------------------------------------------------------------------------
! 2. Initialisation
!-------------------------------------------------------------------------------

  bit_NB = 1
  bit_NBQ = 2

  n = SIZE(time)

  ALLOCATE(TG(n))
  ALLOCATE(DS(2,n))
  ALLOCATE(DDS(n))
  ALLOCATE(DnB(n))
  ALLOCATE(CNB(-lb:lb))

  ! 2.1 Subtract propagation time from observation time

  DO i=1,n
    TG(i) = time(i) - SQRT(SUM((r_gns(i,:) - r_leo(i,:))**2))/ C_light
  ENDDO

  ! 2.2 Determine occultation direction
  PS1 = Impact_Parameter(r_leo(1,:)-r_coc(:), r_gns(1,:)-r_coc(:))
  PSn = Impact_Parameter(r_leo(n,:)-r_coc(:), r_gns(n,:)-r_coc(:))
  
  ocd = NINT(SIGN(1.0_wp, PSn-PS1))

!-------------------------------------------------------------------------------
! 3. Frequency down-conversion and navigation bits removal
!-------------------------------------------------------------------------------

  ! 3.1 Determine open-loop record limits
  
  imin = n
  imax = 1
  
  DO i=1,n
     IF (BTest(LCF(i),0)) THEN
        imin = MIN(imin,i)
        imax = MAX(imax,i)
     END IF
  END DO
  
  ! 3.2 Determine open-loop record fragment for navbit correlation
  
  dk = NINT(4.0/ABS(TG(2)-TG(1)))
 
  SELECT CASE (ocd)
  CASE (-1)
!    iol  = imin ! Commented at 21 July, 2016
     kmin = imin
     kmax = MIN(imax, kmin + (dk-1))
  CASE (+1)
!    iol  = imax ! Commented at 21 July, 2016
     kmax = imax
     kmin = MAX(imin, kmax - (dk-1))
  END SELECT
  
  ! 3.3 Phase model subtraction

  k(1) = 2.0_wp*pi*f_L1/c_light
  k(2) = 2.0_wp*pi*f_L2/c_light

  DS(1,:) = k(1)*(phase_L1(:) - phase_LM(:))
  DS(2,:) = k(2)*(phase_L2(:) - phase_LM(:))
  
  ! 3.4 Navbits-phase correlation

  DO i=2,n
     DDS(i) = ABS(MODULO(pi + DS(1,i) - DS(1,i-1), 2.0_wp*pi) - pi) / pi
     DNB(i) = ABS(IBITS(LCF(i),bit_NB,1)-IBITS(LCF(i-1),bit_NB,1))
  END DO

  DDS(1) = 0.0_wp
  DNB(1) = 0.0_wp

  DO i=-lb,-1
    ni = SUM(DNB(kmin:kmax+i)**2)
    IF (ni > 0) THEN
      CNB(i) = SUM(DDS(kmin-i:kmax)*DNB(kmin:kmax+i)) /    &
                   SQRT(SUM(DDS(kmin-i:kmax)**2)*ni)
    ELSE
      CNB(i) = 0.0_wp
    ENDIF
  END DO
  
  DO i=0,lb
    ni = SUM(DNB(kmin+i:kmax)**2)
    IF (ni > 0) THEN
      CNB(i) = SUM(DDS(kmin:kmax-i)*DNB(kmin+i:kmax)) /  &
                   SQRT(SUM(DDS(kmin:kmax-i)**2)*ni)
    ELSE
      CNB(i) = 0.0_wp
    ENDIF
  END DO
  
  inb    = SUM(MAXLOC(CNB(:))) - 1 - lb
  CNBmax = MAXVAL(CNB(:))
  
  WRITE( outstr, '(F14.1,1X,F14.1)') TG(Kmin), TG(Kmax)
  CALL message(msg_diag,'Navbit correllation area: ' // TRIM(outstr))
  WRITE( outstr, '(F8.3,A,I4)') CNBmax, ' at ', INB
  CALL message(msg_diag,'Max. navbit correlation: ' // TRIM(outstr))

  ! 3.5 Shifting navigation bits to correlation maximum

  IF (inb > 0) THEN
     LCF(1:n-inb) = LCF(1+inb:n)
  ELSE IF (inb < 0) THEN
     LCF(1-inb:n) = LCF(1:n+inb)
  END IF

  ! 3.6 Discard navigation bits for low correlation
  
  IF (CNBmax < 0.6) THEN 
    bit_NB = 7
    bit_NBQ = 8
    CALL ropp_pp_internal_navbit(TG(imin:imax),ds(1,imin:imax),lcf(imin:imax))
  ENDIF
  
! 3.7 Navigation data bits removal or 2-quadrant correction
!     of missing bits and phase accumulation

  DO ic=1,1

   WHERE (BTEST(LCF(imin:imax),bit_NBQ))
      DS(ic,imin:imax) = DS(ic,imin:imax) + pi*IBITS(LCF(imin:imax),bit_NB,1)
   END WHERE
   
   DO i=2,n
      IF (BTEST(LCF(i),0)) THEN  ! In open loop mode
         IF (BTEST(LCF(i-1),bit_NBQ) .AND. BTEST(LCF(i),bit_NBQ)) THEN
            DS(ic,i) = DS(ic,i-1) +      &
               MODULO(DS(ic,i)-DS(ic,i-1)+pi, 2.0_wp*pi)-pi
         ELSE
            DS(ic,i) = DS(ic,i-1) +      &
               MODULO(DS(ic,i)-DS(ic,i-1)+pi/2.0_wp, pi)-pi/2.0_wp
         END IF
       ELSE                               ! In phase locked loop mode
         IF (Restore_PLL) THEN
            DS(ic,i) = DS(ic,i-1) +      &
               MODULO(DS(ic,i)-DS(ic,i-1)+pi/2.0_wp, pi) - pi/2.0_wp
         ELSE
            DS(ic,i) = DS(ic,i-1) +      &
               MODULO(DS(ic,i)-DS(ic,i-1)+pi, 2.0_wp*pi) - pi
         END IF
      END IF
   END DO

 END DO

!-------------------------------------------------------------------------------
! 4. Computation of excess phase
!-------------------------------------------------------------------------------

! 4.1 Phase accumulation

DO ic=1,2
   DO i=2,n
      DS(ic,i) = DS(ic,i-1) + MODULO(DS(ic,i)-DS(ic,i-1)+pi, 2.0_wp*pi) - pi
   END DO
END DO

! 4.2 Discontinuity removal

DO ic=1,2

   SELECT CASE (ocd)
      CASE (-1)
         IF (imin > 1) THEN
            SB = DS(ic,imin-1) - DS(ic,imin)
         ELSE
            SB = 0.0_wp
         END IF
      CASE (+1)
         IF (imax < n) THEN
            SB = DS(ic,imax+1) - DS(ic,imax)
         ELSE
            SB = 0.0_wp
         END IF
   END SELECT
   DS(ic,imin:imax) = DS(ic,imin:imax) + SB
   
END DO

! 4.3 Excess phase computation

phase_L1(:) = phase_LM(:) + DS(1,:)/k(1)
phase_L2(:) = phase_LM(:) + DS(2,:)/k(2)

!-------------------------------------------------------------------------------
! 5. Clean up
!-------------------------------------------------------------------------------

DEALLOCATE(TG)
DEALLOCATE(DS)
DEALLOCATE(DDS)
DEALLOCATE(DNB)
DEALLOCATE(CNB)

  CALL message_set_routine(routine)

CONTAINS

!****s* Preprocessing/ropp_pp_internal_navbit *
!
! NAME
!    ropp_pp_internal_navbit - get internal navigation bits
!
! SYNOPSIS
!    call ropp_pp_internal_navbit(time, ds, lcf)
!
! DESCRIPTION
!    Determine navigation bits frame position, polynomial regression of
!    phase in frames and detection of phase jumps by Pi
!
!****

  SUBROUTINE ropp_pp_internal_navbit(time, ds, lcf)

    IMPLICIT NONE
    
    REAL(wp), DIMENSION(:), INTENT(in)    :: time  ! OL time grid
    REAL(wp), DIMENSION(:), INTENT(in)    :: ds    ! phase deviation from model
    INTEGER,  DIMENSION(:), INTENT(inout) :: lcf   ! lost carrier flag

    INTEGER          :: n        ! Number of data
    INTEGER          :: i        ! Array index
    REAL(wp)         :: DSF      ! Phase change between frames

    REAL(wp), DIMENSION(:), ALLOCATABLE :: TGM  ! Time from navbit start
    INTEGER,  DIMENSION(:), ALLOCATABLE :: TGI  ! Navbit number
    INTEGER,  DIMENSION(:), ALLOCATABLE :: NB   ! Navigation bits

    REAL(wp), PARAMETER  :: gpsfl = 0.02   ! GPS navbit length
    
    ! 6.1 Memory allocation 
    
    n = SIZE(time)
    ALLOCATE(TGM(n))
    ALLOCATE(TGI(n))
    ALLOCATE(NB(n))
 
    ! 6.2 Positioning navbit frames
    
    TGM(:) = MODULO(time(:), GPSFL)
    TGI(:) = FLOOR(time(:)/GPSFL)
        
    ! 6.3 Comparison with external navbits positioning
        
    NB(:) = 0
    
    DO i=2,n
      If ((TGM(i-1) < TGM(i)) .and. (IBits(LCF(i-1),1,1) /= IBits(LCF(i),1,1))) then
        NB(i) = 1
      Else
        NB(i) = 0
      End If
    End Do

    
    Do i=2,N
      If ((TGI(i-1) == TGI(i)) .and. (IBits(LCF(i-1),1,1) /= IBits(LCF(i),1,1))) then
        NB(i) = 1
      Else
        NB(i) = 0
      End If
    End Do
    
    ! 6.4 Retrieval of navbits

    DO i=2,N
      IF (TGI(i-1) /= TGI(i)) THEN
        DSF = ABS(MODULO(DS(i) - DS(i-1) + Pi, 2*Pi) - Pi)
        IF (DSF > Pi/2) THEN
          NB(i) = MODULO(NB(i-1) + 1,2)
       ELSE
         NB(i) = NB(i-1)
       ENDIF
     ELSE
       NB(i) = NB(i-1)
     ENDIF
   ENDDO
   
   ! 6.5 Setting LCF
   
   WHERE (NB(:) == 1)
     LCF(:) = IBSET(LCF(:),7)
   ENDWHERE
   
   LCF(:) = IBSET(LCF(:),8)
   
   ! 6.6 Clean Up

   DEALLOCATE(TGM)
   DEALLOCATE(TGI)
   DEALLOCATE(NB)
   
 END SUBROUTINE ropp_pp_internal_navbit

END SUBROUTINE ropp_pp_openloop
