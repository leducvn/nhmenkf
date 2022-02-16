! $Id: ropp_pp_refrac_msis.f90 3551 2013-02-25 09:51:28Z idculv $

!****s* ModelRefraction/ropp_pp_refrac_MSIS *
!
! NAME
!    ropp_pp_refrac_MSIS - Compute refractivity from MSIS spherical harmonics
!
! SYNOPSIS
!    call ropp_pp_refrac_MSIS(file, month, lat, lon, alt, refrac, grad)
!
! DESCRIPTION
!    This subroutine calculates a climatological refractivity profile for 
!    a given month, latitude and longitude from MSIS data. Profiles are 
!    computed using Chebyshev polynomials and spherical harmonics with 
!    coefficients read from file.
!
! INPUTS
!    character(len=*) :: file     Model coefficients filename
!    integer,         :: month    Month of year
!    real(wp),        :: lat      Latitude
!    real(wp),        :: lon      Longitude
!    real(wp), dim(:) :: alt      Altitude levels on which to find N
!
! OUTPUT
!    character(len=*) :: file     Model coefficients filename
!    real(wp), dim(:) :: refrac   Refractivity field
!    real(wp), dim(:) :: grad     Refractivity field gradient (optional)
!
! AUTHOR
!   Met Office, Exeter, UK.
!   Any comments on this software should be given via the ROM SAF
!   Helpdesk at http://www.romsaf.org
!
! COPYRIGHT
!   Copyright (c) 1998 Stig Syndergaard <ssy@dmi.dk>.
!   For further details please refer to the file COPYRIGHT
!   which you should have received as part of this distribution.
!
!****

SUBROUTINE ropp_pp_refrac_MSIS(mfile, month, lat, lon, alt, refrac, grad_refrac)

!-------------------------------------------------------------------------------
! 1. Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal
! USE ropp_pp_MSIS, not_this => ropp_pp_refrac_MSIS
  USE ropp_pp_MSIS

  IMPLICIT NONE

  CHARACTER(len=*),    INTENT(inout)  :: mfile  ! Coefficients file
  INTEGER,                INTENT(in)  :: month    ! Month of year
  REAL(wp),               INTENT(in)  :: lat      ! Latitude
  REAL(wp),               INTENT(in)  :: lon      ! Longitude
  REAL(wp), DIMENSION(:), INTENT(in)  :: alt      ! Altitude (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: refrac   ! Refractivity (N-units)
  REAL(wp), DIMENSION(:), OPTIONAL, INTENT(out) :: grad_refrac  ! dN/dZ
  
  TYPE(MSIScoeff_ref), SAVE           :: coeff    ! Spherical harmonic coeffs
  
  REAL(wp), DIMENSION(:), ALLOCATABLE :: n10      ! Constants
  REAL(wp), DIMENSION(:), ALLOCATABLE :: n11
  REAL(wp), DIMENSION(:), ALLOCATABLE :: n20
  REAL(wp), DIMENSION(:), ALLOCATABLE :: n21
  REAL(wp), DIMENSION(:), ALLOCATABLE :: n10x
  REAL(wp), DIMENSION(:), ALLOCATABLE :: n11x
  REAL(wp), DIMENSION(:), ALLOCATABLE :: c_tmp

  REAL(wp), DIMENSION(:), ALLOCATABLE :: refrac_high
  REAL(wp), DIMENSION(:), ALLOCATABLE :: refrac_low

  REAL(wp) :: phi, theta, cosphi, sinphi, costheta, sintheta     
  REAL(wp) :: alpha, beta, refrac0
  
  REAL(wp), PARAMETER  :: dh = 10.0_wp           ! finite difference step (m)
  REAL(wp), PARAMETER  :: deg2rad = 0.0174532925_wp   ! degrees to radians
  INTEGER              :: i, k, ncoeff, ntmp
  LOGICAL              :: found
  CHARACTER(len=*), PARAMETER ::  &
       filepath(9) =  (/ 'data/MSIS_coeff.nc             ', &
                         '../data/MSIS_coeff.nc          ', &
                         '../../data/MSIS_coeff.nc       ', &
                         '../../../data/MSIS_coeff.nc    ', &
                         '../../../../data/MSIS_coeff.nc ', &
                         '*/data/MSIS_coeff.nc           ', &
                         '*/*/data/MSIS_coeff.nc         ', &
                         '*/*/*/data/MSIS_coeff.nc       ', &
                         '*/*/*/*/data/MSIS_coeff.nc     ' /)

!-------------------------------------------------------------------------------
! 2. Read coefficients file
!-------------------------------------------------------------------------------
  IF (.not. MSIS_read) THEN

    INQUIRE(File=mfile, Exist=Found)
    IF (.NOT. Found) THEN
      DO i=1,9
        INQUIRE(File=filepath(i), Exist=Found)
        IF (Found) THEN
          mfile = filepath(i)
          EXIT
        ENDIF
      ENDDO
    ENDIF
        
    CALL ropp_pp_read_MSIS(mfile, month, coeff)

    MSIS_read = .true.

  ENDIF
  
!-------------------------------------------------------------------------------
! 3. Constants for Clenshaw's recurrence formula
!-------------------------------------------------------------------------------
    
  ncoeff = SIZE(coeff%Aa1)
  ntmp   = SIZE(coeff%Ac0,2)
  
  ALLOCATE(n10(ncoeff))
  ALLOCATE(n11(ncoeff))
  ALLOCATE(n20(ncoeff))
  ALLOCATE(n21(ncoeff))
  ALLOCATE(n10x(ncoeff))
  ALLOCATE(n11x(ncoeff))
  ALLOCATE(c_tmp(ntmp))

  DO i=1,ncoeff
     n10(i)=(2.0_wp*i + 1.0_wp)/(i + 1.0_wp)
     n11(i)=(2.0_wp*i + 1.0_wp)/(i*1.0_wp)
     n20(i)=(i + 1.0_wp)/(i + 2.0_wp)
     n21(i)=(i + 2.0_wp)/(i + 1.0_wp)
  ENDDO
  
!-------------------------------------------------------------------------------
! 4. Calculate lat/lon dependent parameters
!-------------------------------------------------------------------------------

  phi = lon * deg2rad
  theta = (90.0_wp - lat) * deg2rad
  cosphi = COS(phi)
  sinphi = SIN(phi)
  costheta = COS(theta)
  sintheta = SIN(theta)

  n10x(:) = costheta * n10(:)
  n11x(:) = costheta * n11(:)
  
  DO i=1,ntmp
     c_tmp(i) = Spheric(coeff%Ac0(:,i), coeff%Ac1(:,i), coeff%Bc1(:,i))
  ENDDO

  alpha   = Spheric(coeff%Aa0, coeff%Aa1, coeff%Ba1)
  beta    = Spheric(coeff%Ab0, coeff%Ab1, coeff%Bb1)
  refrac0 = Spheric(coeff%Ad0, coeff%Ad1, coeff%Bd1)

!-------------------------------------------------------------------------------
! 5. Compute refractivity (and gradient)
!-------------------------------------------------------------------------------

  CALL calc_refrac(alt, refrac, alpha, beta, refrac0)

  IF(PRESENT(grad_refrac))THEN
     ALLOCATE(refrac_high(SIZE(refrac)))
     ALLOCATE(refrac_low(SIZE(refrac)))
     CALL calc_refrac(alt(:)+dh, refrac_high, alpha, beta, refrac0)
     CALL calc_refrac(alt(:)-dh, refrac_low, alpha, beta, refrac0)
     grad_refrac = (refrac_high - refrac_low)/(2.0_wp*dh)            ! N/m
     DEALLOCATE(refrac_high)
     DEALLOCATE(refrac_low)
  ENDIF

!-------------------------------------------------------------------------------
! 6. Clean up
!-------------------------------------------------------------------------------
  
  DEALLOCATE(c_tmp)
  DEALLOCATE(n10)
  DEALLOCATE(n11)
  DEALLOCATE(n20)
  DEALLOCATE(n21)
  DEALLOCATE(n10x)
  DEALLOCATE(n11x)

!  DEALLOCATE(coeff%Ac0)
!  DEALLOCATE(coeff%Ac1)
!  DEALLOCATE(coeff%Bc1)
!  DEALLOCATE(coeff%Aa0)
!  DEALLOCATE(coeff%Aa1)
!  DEALLOCATE(coeff%Ba1)
!  DEALLOCATE(coeff%Ab0)
!  DEALLOCATE(coeff%Ab1)
!  DEALLOCATE(coeff%Bb1)  
!  DEALLOCATE(coeff%Ad0)
!  DEALLOCATE(coeff%Ad1)
!  DEALLOCATE(coeff%Bd1)  

CONTAINS

!-------------------------------------------------------------------------------
! 7. Calculation of refractivity - Chebyshev polynomial expansion
!-------------------------------------------------------------------------------
  SUBROUTINE calc_refrac(alt, refrac, alpha, beta, refrac0)

    USE typesizes, ONLY: wp => EightByteReal

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(in)  :: alt
    REAL(wp), DIMENSION(:), INTENT(out) :: refrac
    REAL(wp),               INTENT(in)  :: alpha
    REAL(wp),               INTENT(in)  :: beta
    REAL(wp),               INTENT(in)  :: refrac0

    REAL(wp) :: x, z
    REAL(wp) :: dk, dk1, dk2
    REAL(wp) :: fx, hs
    REAL(wp) :: a0
    REAL(wp) :: alt_km
    REAL(wp), PARAMETER  :: zref = 100.0_wp   
    REAL(wp), PARAMETER  :: h0 = 100.0_wp   
    INTEGER  :: i

    a0 = 0.0_wp
    
    DO i=1,SIZE(alt)
    
       ! 6.0 Change units (m to km)
       
       alt_km = alt(i)/1.e3_wp
   
       ! 6.1 Change of variables
       
       IF(alt_km < a0) THEN 
          z = TANH(alt_km - a0)/ h0
       ELSE
          z = (alt_km - a0) / h0
     ENDIF
     
     x = 1.0_wp - 2.0_wp*EXP(-z)

     ! 5.2 Clenshaw's recurrence formula (Chebychev polynomials)

     dk1 = 0.0_wp
     dk2 = 0.0_wp

     DO k=ntmp,2,-1
        dk = 2.0_wp*x*dk1 - dk2 + c_tmp(k)
        dk2 = dk1
        dk1 = dk
     ENDDO
     
     fx = x*dk1 - dk2 + 0.5_wp*c_tmp(1)

     ! 5.3 Compute refractivity scale height

     hs = fx * EXP(-(z/zref)**2) + (alpha*z + beta)

     ! 5.4 Compute refractivity

     refrac(i) = 1.0e6_wp*refrac0*EXP(-(alt_km - a0)/hs)

  ENDDO

END SUBROUTINE calc_refrac


!-------------------------------------------------------------------------------
! 8. Calculation of spherical harmonics - Clenshaw's recurrence formula
!-------------------------------------------------------------------------------
  FUNCTION Spheric(A0, A1, B1) RESULT(sharm)
    
    USE typesizes, ONLY: wp => EightByteReal

    IMPLICIT NONE

    REAL(wp), DIMENSION(:), INTENT(in) :: A0
    REAL(wp), DIMENSION(:), INTENT(in) :: A1
    REAL(wp), DIMENSION(:), INTENT(in) :: B1
    REAL(wp)                           :: sharm

    REAL(wp) :: dn10, dn20, dn11, dn21, dn0, dn1
    INTEGER  :: n
    
    dn10 = 0.0_wp
    dn20 = 0.0_wp
    dn11 = 0.0_wp
    dn21 = 0.0_wp

    DO n=SIZE(A1),1,-1
       dn0 = n10x(n)*dn10 - n20(n)*dn20 + A0(n+1)
       dn1 = n11x(n)*dn11 - n21(n)*dn21 + A1(n)*cosphi + B1(n)*sinphi
       dn20 = dn10
       dn21 = dn11
       dn10 = dn0
       dn11 = dn1
    ENDDO
    
    sharm = costheta*dn10 + sintheta*dn11 - 0.5_wp*dn20 + A0(1)

  END FUNCTION Spheric
    

END SUBROUTINE ropp_pp_refrac_MSIS

