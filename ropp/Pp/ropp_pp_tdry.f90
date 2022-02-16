! $Id: ropp_pp_tdry.f90 1960 2008-11-13 16:28:24Z frhl $

SUBROUTINE ropp_pp_tdry(lat, alt, refrac, shum, t_dry, p_dry, Zmax)

!****s* Meteo/ropp_pp_tdry *
!
! NAME
!    ropp_pp_tdry - Compute temperature (and pressure profile) from 
!                   refractivity assuming zero humidity.
!
! SYNOPSIS
!    call ropp_pp_tdry(lat, alt, refrac, shum, t_dry, p_dry, Zmax)
!
! DESCRIPTION
!    Calculate P(z) from numerical integration of barometric formula
!
!            d ln(P(z))                    g(z)
!            ---------- = - -------------------------------------
!               dz           Rd T(N(z),P(z),Q(z)) (1 + eps(Q(z))
!
!    given N(z) and Q(z) and known dependence T(N, P, Q).
!    Calculate T(z) = T(N(z), P(z), Q(z)). 
!
! INPUTS
!    real(wp)                 :: lat       ! latitude
!    real(wp), dimension(:)   :: alt       ! altitude (z)
!    real(wp), dimension(:)   :: refrac    ! refraction (N(z))
!    real(wp), dimension(:)   :: shum      ! specific humidity (Q(z))
!    real(wp), optional       :: zmax      ! upper integration height 
!
! OUTPUT
!    real(wp), dimension(:)   :: t_dry     ! dry temperature (T(z))
!    real(wp), dimension(:)   :: p_dry     ! dry pressure (P(z))
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
! USE ropp_pp, not_this => ropp_pp_tdry
  USE ropp_pp
  
  IMPLICIT NONE
  
  REAL(wp),                       INTENT(in)  :: lat    ! Latitude
  REAL(wp), DIMENSION(:), TARGET, INTENT(in)  :: alt    ! Altitude (m)
  REAL(wp), DIMENSION(:),         INTENT(in)  :: refrac ! Refractivity 
  REAL(wp), DIMENSION(:), TARGET, INTENT(in)  :: shum   ! Spec humidity (kg/kg)
  REAL(wp), DIMENSION(:),         INTENT(out) :: t_dry  ! Dry temperature (K)
  REAL(wp), DIMENSION(:),         INTENT(out) :: p_dry  ! Dry pressure (hPa)
  REAL(wp), OPTIONAL,             INTENT(in)  :: Zmax   ! Maximum altitude (m)

  REAL(wp), TARGET, DIMENSION(:), ALLOCATABLE :: d2N    ! 2nd derivative ln(N)
  REAL(wp), TARGET, DIMENSION(:), ALLOCATABLE :: d2Q    ! 2nd derivative shum
  REAL(wp), TARGET, DIMENSION(:), ALLOCATABLE :: lnN    ! ln(N)

  REAL(wp), DIMENSION(:), ALLOCATABLE         :: ZI     ! Altitude grid
  REAL(wp), DIMENSION(:), ALLOCATABLE         :: TZ     ! Temperature on ZI
  REAL(wp), DIMENSION(:), ALLOCATABLE         :: D2T    ! 2nd derivative TZ
  REAL(wp), DIMENSION(:), ALLOCATABLE         :: LnPZ   ! ln(P) on ZI
  REAL(wp), DIMENSION(:), ALLOCATABLE         :: D2P    ! 2nd derivative LnPZ

  REAL(wp), DIMENSION(:), POINTER             :: PZ     ! Pointer to alt
  REAL(wp), DIMENSION(:), POINTER             :: PN     ! Pointer to refrac 
  REAL(wp), DIMENSION(:), POINTER             :: P2N    ! Pointer to d2N
  REAL(wp), DIMENSION(:), POINTER             :: PQ     ! Pointer to shum
  REAL(wp), DIMENSION(:), POINTER             :: P2Q    ! Pointer to d2Q
      
  REAL(wp)                                    :: GCLat  ! Geocentric latitude
  REAL(wp)                                    :: LnP    ! ln(p)
  REAL(wp)                                    :: LnNZ   ! ln(NZ) on ZI
  REAL(wp)                                    :: DLnNZ  ! 1st derivative LnNZ
  REAL(wp)                                    :: QZ     ! shum on ZI
  REAL(wp)                                    :: ZP     ! Altitude value
  REAL(wp)                                    :: Z1     ! Maximum altitude
  REAL(wp)                                    :: Zmin   ! Minimum altitude
  REAL(wp)                                    :: dz     ! Integration step
  INTEGER                                     :: KZ     ! No. integration steps
  INTEGER                                     :: i      ! Index

!  REAL(wp), PARAMETER :: dzi = 20.0_wp   ! Integration step size (m)
  REAL(wp), PARAMETER :: dzi = 15.0_wp   ! Integration step size (m)
  REAL(wp), PARAMETER :: Qmin = 1e-7_wp  ! Minimum value for QZ


  IF (ANY( refrac <= 0.0_wp )) THEN

    CALL message(msg_warn,   &
         "Cannot generate Tdry from non-positive refractivities \n")
    t_dry = refrac*0.0_wp + ropp_MDFV
    p_dry = refrac*0.0_wp + ropp_MDFV

  ELSE

!-------------------------------------------------------------------------------
! 2. Initialization
!-------------------------------------------------------------------------------

    ALLOCATE(d2N(SIZE(alt)))
    ALLOCATE(d2Q(SIZE(alt)))
    ALLOCATE(lnN(SIZE(alt)))

!-------------------------------------------------------------------------------
! 3. Calculate spline coefficients
!-------------------------------------------------------------------------------

    lnN = LOG(refrac)
    CALL ropp_pp_init_spline(alt(:), lnN(:), d2N)
    CALL ropp_pp_init_spline(alt(:), shum(:), d2Q)

!-------------------------------------------------------------------------------
! 4. Set up global variables for FTZP
!-------------------------------------------------------------------------------

    PZ   => alt
    PN   => lnN
    P2N  => d2N
    PQ   => shum
    P2Q  => d2Q

!-------------------------------------------------------------------------------
! 5. Compute number of steps and step size
!-------------------------------------------------------------------------------

    IF (PRESENT(Zmax)) THEN 
      Z1 = Zmax
    ELSE
      Z1 = MAXVAL(alt(:))
    ENDIF
    Zmin = MINVAL(alt(:))
    KZ = CEILING((Z1 - Zmin)/ABS(DZI))
    dz = -(Z1 - Zmin)/KZ

!-------------------------------------------------------------------------------
! 6. Array allocation
!-------------------------------------------------------------------------------

    ALLOCATE(ZI(0:KZ))
    ALLOCATE(TZ(0:KZ))
    ALLOCATE(D2T(0:KZ))
    ALLOCATE(LnPZ(0:KZ))
    ALLOCATE(D2P(0:KZ))

!-------------------------------------------------------------------------------
! 7. Numerical integration
!-------------------------------------------------------------------------------

    ! 7.1 Initial conditions

    ZP = Z1
    ZI(KZ) = ZP

    CALL ropp_pp_interpol_spline(alt, lnN, d2N, Z1, LnNZ, DLnNZ)
  
    GCLat = GCLat_from_GDLat(lat)    ! Geocentric latitude
    TZ(KZ) = - GravityGC(GCLat, ZP) / (R_dry*DLnNZ)
  
    LnP = LnNZ + LOG(TZ(KZ)/kappa1)
    LnPZ(KZ) = LnP

    ! 7.2 Integration

    DO i= kz-1, 0, -1
      CALL ropp_pp_runge_kutta(dz, zp, lnP)
      CALL ropp_pp_interpol_spline(pz, pn, p2n, zp, lnNz)
      CALL ropp_pp_interpol_spline(pz, pq, p2q, zp, QZ)

      QZ = MAX(Qmin, QZ)
      ZI(i) = ZP
      TZ(i) = T_from_NPQ(EXP(LnNZ), EXP(LnP), QZ)
      LNPZ(i) = LnP
   ENDDO

!-------------------------------------------------------------------------------
! 8. Interpolation of calculated T(Z) and P(Z)
!-------------------------------------------------------------------------------

    ! 8.1 Temperature
    CALL ropp_pp_init_spline(ZI, TZ, D2T)
    DO i=1,SIZE(alt)
      CALL ropp_pp_interpol_spline(zi, tz, d2t, alt(i), t_dry(i))
    ENDDO

    ! 8.2 Pressure
    CALL ropp_pp_init_spline(ZI, LnPZ, D2P)
    DO i=1,SIZE(alt)
      CALL ropp_pp_interpol_spline(zi, lnPZ, d2p, alt(i), p_dry(i))
    ENDDO

    p_dry(:) = EXP(p_dry(:))

!-------------------------------------------------------------------------------
! 9. Clean up
!-------------------------------------------------------------------------------

    DEALLOCATE(ZI)
    DEALLOCATE(TZ) 
    DEALLOCATE(D2T)
    DEALLOCATE(LnPZ)
    DEALLOCATE(D2P)

    DEALLOCATE(d2N)
    DEALLOCATE(d2Q)
    DEALLOCATE(lnN)


  ENDIF  ! Refrac > 0

CONTAINS

!-------------------------------------------------------------------------------
! 10. Conversion Geodetic to Geocentric latitude
!-------------------------------------------------------------------------------
  FUNCTION GCLat_from_GDLat(GDLat) RESULT(GCLat)
    
    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_constants, ONLY: Pi
    IMPLICIT NONE

    REAL(wp), INTENT(in) :: GDLat
    REAL(wp) :: GCLat

    REAL(wp) :: S, X(3), gd, flatfn, funsq
    REAL(wp), PARAMETER :: Re  = 6378137.0_wp
    REAL(wp), PARAMETER :: dtr = Pi/180.0_wp
    REAL(wp), PARAMETER :: f_ell = 1.0_wp / 298.257223563_wp

    ! 10.1 Convert geodetic to cartesian coordinates

    flatfn = (2.0_wp - f_ell)*f_ell
    funsq = (1.0_wp - f_ell)**2
    gd = Re / Sqrt(1.0_wp - flatfn*Sin(GDLat*dtr)**2)

    X(2) = COS(GDLat*dtr)*SIN(0.0_wp*dtr)*gd
    X(1) = COS(GDLat*dtr)*COS(0.0_wp*dtr)*gd
    X(3) = SIN(GDLat*dtr)*(gd*funsq)

    ! 10.2 Convert geodetic to spherical coordinates
    
    S = ACos(X(3)/(Sqrt(Sum(X(:)**2))))

    ! 10.3 Compute geodetic latitude
    
    GCLat = Pi/2.0_wp - S

  END FUNCTION GCLat_from_GDLat

!-------------------------------------------------------------------------------
! 11. Gravity calculation
!-------------------------------------------------------------------------------
  FUNCTION GravityGC(Lat, Alt)

    USE typesizes, ONLY: wp => EightByteReal
    IMPLICIT NONE

    REAL(wp), INTENT(in) :: Lat
    REAL(wp), INTENT(in) :: Alt
    REAL(wp) :: GravityGC

    REAL(wp) :: g_sur, R0, ge, f2, f4
    REAL(wp), PARAMETER :: Re  = 6378137.0_wp
    REAL(wp), PARAMETER :: m = 0.00345_wp
    REAL(wp), PARAMETER :: f_ell = 1.0_wp / 298.257223563_wp
    
    f2 = -f_ell + 5.0_wp*m/2.0_wp - 17.0_wp*f_ell*m/14.0_wp +   &
            15.0_wp*m**2/4.0_wp
    f4 = -f_ell**2/2.0_wp + 5.0_wp*f_ell*m/2.0_wp

    ge = 3.986004415e14_wp / (Re**2 * &
       (1.0_wp - f_ell + 3.0_wp*m/2.0_wp - 15.0_wp*m*f_ell/14.0_wp))

    g_sur = ge*(1.0_wp + f2*(Sin(Lat))**2 - (f4*(Sin(2.0_wp*Lat))**2)/4.0_wp)
    
    R0 = (g_sur/ge)*Re / (1.0_wp + f_ell + m +      &
            (-3.0_wp*f_ell + 5.0_wp*m/2.0_wp)*(Sin(Lat))**2)

    GravityGC = g_sur*(R0/(R0 + Alt))**2

  END FUNCTION GravityGC

!-------------------------------------------------------------------------------
! 12. Right-hand size of barometric formula
!-------------------------------------------------------------------------------
  
  FUNCTION FTZP(ZP, LnP)   
    
    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_constants, ONLY: epsilon_water, R_dry
    IMPLICIT NONE
    
    REAL(wp), INTENT(in) :: zp   ! altitude above reference ellipsoid
    REAL(wp), INTENT(in) :: lnP  ! logarithm of pressure
    REAL(wp)             :: ftzp ! right part of barometric formula
    
    REAL(wp)             :: lnNZ ! interpolated ln(N(z))
    REAL(wp)             :: QZ   ! interpolated Q(Z)
    REAL(wp)             :: TZP  ! T(z, P(z))
    REAL(wp)             :: GZ   ! gravity acceleration g(lat, z)

    CALL ropp_pp_interpol_spline(pz, pn, p2n, zp, lnNz)
    CALL ropp_pp_interpol_spline(pz, pq, p2q, zp, qz)

    qz = MAX(Qmin, Qz)
    TzP = T_from_NPQ(EXP(LnNZ), EXP(LnP), QZ)
    GZ = GravityGC(GClat, zp)

    ftzp = -GZ  / (R_dry * Tzp * (1.0_wp + (1.0_wp/epsilon_water-1.0_wp)*QZ))

  END FUNCTION FTZP
  
!-------------------------------------------------------------------------------
! 13. Calculate temperature from refractivity, pressure and humidity
!-------------------------------------------------------------------------------
 
  FUNCTION T_from_NPQ(N, P, Q) RESULT(T)

    USE typesizes, ONLY: wp => EightByteReal
    USE ropp_pp_constants, ONLY: R_dry, R_vap, kappa1, kappa2
    IMPLICIT NONE

    REAL(wp), INTENT(in) :: N
    REAL(wp), INTENT(in) :: P
    REAL(wp), INTENT(in) :: Q
    REAL(wp)             :: T

    REAL(wp), PARAMETER :: aq =  R_dry/R_vap
    REAL(wp), PARAMETER :: bq =  1.0_wp - aq

    ! 11.1 Analytical solution

    T = (kappa1*P + SQRT((kappa1*P)**2 + 4.0_wp*kappa2*N*P*Q/(aq + bq*Q))) /   &
         (2.0_wp*N)

  END FUNCTION T_from_NPQ
  

!-------------------------------------------------------------------------------
! 14. Runge_kutta numerical integration of X' = F(X,t)
!-------------------------------------------------------------------------------

  SUBROUTINE ropp_pp_runge_kutta(Dt, t, X)


    USE typesizes, ONLY: wp => EightByteReal
    
    IMPLICIT NONE

    REAL(wp),  INTENT(in)    :: dt   ! time integration step
    REAL(wp),  INTENT(inout) :: t    ! time variable
    REAL(wp),  INTENT(inout) :: X    ! dynamic variable vector
    
    REAL(wp) :: K1, K2, K3, K4
    
    K1 = FTZP(t,        X)
    K2 = FTZP(t + Dt/2, X + K1*Dt/2)
    K3 = FTZP(t + Dt/2, X + K2*Dt/2)
    K4 = FTZP(t + Dt,   X + K3*Dt)
    
    X  = X + (K1 + 2*K2 + 2*K3 + K4)*Dt/6
    t     = t + Dt
    
  END SUBROUTINE ropp_pp_runge_kutta
  
END SUBROUTINE ropp_pp_tdry
