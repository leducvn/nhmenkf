! $Id: ropp_pp_search_model_refraction.f90 2048 2009-04-07 15:45:10Z frhl $

SUBROUTINE ropp_pp_search_model_refraction(mfile, in_month, in_lat, in_lon, &
                                           in_impact, in_bangle,            & 
                                           impact, bangle_MSIS, config)

!****s* ModelRefraction/ropp_pp_search_model_refraction *
!
! NAME
!    ropp_pp_search_model_refraction - Calculate best-fit bending angle profile
!                                      for MSIS refractivity
!
! SYNOPSIS
!    call ropp_pp_model_refraction(mfile, in_month, in_lat, in_lon, 
!                                  in_impact, in_bangle, 
!                                  impact, bangle_MSIS, config)
! 
! DESCRIPTION
!    This subroutine calculates the best-fit bending angle profile from the
!    MSIS refractivity field by regression to an input bending angle profile 
!
! INPUT
!    character(len=*) :: mfile          Model coefficients file
!    integer,         :: in_month       Month of year
!    real(wp)         :: in_lat         Latitude  (deg)
!    real(wp)         :: in_lon         Longitude (deg)
!    real(wp), dim(:) :: in_impact      Input impact parameter (m)
!    real(wp), dim(:) :: in_bangle      Input bending angle (rad)
!    real(wp), dim(:) :: impact         Impact parameter (m)
!    type(PPConfig)   :: config         Configuration parameters
!
! OUTPUT
!    real(wp), dim(:) :: bangle_MSIS    Best-fit MSIS bending angles 
!                                       (on input impact parameter levels).
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
  USE messages
  USE ropp_utils, ONLY: WHERE
  USE ropp_pp_utils, ONLY: ropp_pp_regression
! USE ropp_pp, not_this => ropp_pp_search_model_refraction
  USE ropp_pp
  USE ropp_pp_types, ONLY: PPConfig

  IMPLICIT NONE

  CHARACTER(len=*),    INTENT(inout)  :: mfile       ! MSIS file path
  INTEGER,                INTENT(in)  :: in_month    ! Month of year
  REAL(wp),               INTENT(in)  :: in_lat      ! Latitude
  REAL(wp),               INTENT(in)  :: in_lon      ! Longitude
  REAL(wp), DIMENSION(:), INTENT(in)  :: in_impact   ! Input impact parameter (m)
  REAL(wp), DIMENSION(:), INTENT(in)  :: in_bangle   ! Input bending angle (rad)
  REAL(wp), DIMENSION(:), INTENT(in)  :: impact      ! MSIS impact parameter (m)
  REAL(wp), DIMENSION(:), INTENT(out) :: bangle_MSIS ! MSIS bending angle (rad)
  TYPE(PPConfig),         INTENT(in)  :: config      ! Configuration options

  INTEGER                             :: i, imonth, ilat, ilon   ! Index
  INTEGER                             :: n           ! No. of hi-res grid points

  REAL(wp), DIMENSION(:),ALLOCATABLE  :: ba          ! MSIS Bending angle 
  REAL(wp), DIMENSION(:),ALLOCATABLE  :: bangle_nh   ! Interpolated input bangle
  INTEGER,  DIMENSION(:), POINTER :: idx => NULL()   ! Array indices
  
  INTEGER                             :: month      ! Best-fit month
  REAL(wp)                            :: msislat    ! Best-fit latitude
  REAL(wp)                            :: msislon    ! Best-fit longitude
  REAL(wp)                            :: scanlat    ! Scanned latitude
  REAL(wp)                            :: scanlon    ! Scanned longitude
  INTEGER,  PARAMETER :: Nlat = 19        ! Number of latitudes to scan
  INTEGER,  PARAMETER :: Nlon = 18        ! Number of longitudes to scan
  INTEGER             :: n_reg, nrmax     ! Number of regression points
  REAL(wp)            :: mrnr, drnr
  REAL(wp)            :: d, dmin
  REAL(wp), DIMENSION(2)                :: rf           ! Regression factor
  REAL(wp), DIMENSION(:,:), ALLOCATABLE :: K            ! matrix for regression
  REAL(wp), DIMENSION(:), ALLOCATABLE   :: Y            ! regression data

  INTEGER            :: di
  INTEGER            :: IFmin    ! Start index of fitting area
  INTEGER            :: IFmax    ! End index of fitting area
  INTEGER            :: Imin     ! Lower index of fitting area
  INTEGER            :: Imax     ! Upper index of fitting area
  INTEGER            :: Jmax     ! Upper index of fitting area
  CHARACTER(len=247) :: outstr

!-------------------------------------------------------------------------------
! 2. Initialise 
!-------------------------------------------------------------------------------

  n = SIZE(impact)
  ALLOCATE(bangle_nh(n))
  ALLOCATE(ba(n))

  month = in_month
  msislat = in_lat
  msislon = in_lon

!-------------------------------------------------------------------------------
! 3. Interpolate
!------------------------------------------------------------------------------- 

  CALL ropp_pp_interpol(in_impact, impact, in_bangle, bangle_nh)

!-------------------------------------------------------------------------------
! 4. Global MSIS Search
!-------------------------------------------------------------------------------
  
  ! 4.1 Set regression interval mask

  idx => WHERE ( impact-config%r_curve >=  config%hmin_fit .AND.   &
                     impact-config%r_curve <= config%hmax_fit, n_reg)

  IF (n_reg < 100) CALL message(msg_warn, 'Too few data for global search')

  ! 4.2 Loop over all month/lat/lon combinations

  dmin = Huge(dmin)
  nrmax = 0

  DO imonth = 1,12
    
    MSIS_read = .false.     ! Reset MSIS read flag

    DO ilat = 1, nlat

      scanlat = -90.0_wp + REAL(ilat-1)*180.0_wp/REAL(nlat-1.0_wp)
      
      DO ilon = 1, nlon

        scanlon = REAL(ilon - 1.0_wp)*360.0_wp/REAL(nlon,wp)
        
        CALL ropp_pp_bangle_MSIS(mfile,imonth,scanlat,scanlon,   &
                                 impact-config%r_curve,ba)

        ! 4.3 1-parameter fit

        IF (config%nparm_fit == 1) THEN
          
          nrmax = MAX(n_reg, nrmax)
          mrnr = SUM(bangle_nh(idx)*ba(idx))/n_reg
          drnr = SUM(ba(idx)**2)/n_reg
          
          d = ABS(mrnr/drnr - 1.0_wp)

          IF (d < dmin) THEN
            dmin = d
            bangle_MSIS(:) = ba(:)
            month = imonth
            msislat = scanlat
            msislon = scanlon
          ENDIF
          
        ! 4.4 2-parameter fit
          
        ELSE IF (config%nparm_fit == 2) THEN
          
            imin = MINVAL(idx) 
            imax = MAXVAL(idx)
            di = SIGN(1, imax-imin)

            jmax = imin

            DO i=imin,imax,di
              IF (ABS(bangle_nh(i) - ba(i)) < 0.3*ba(i)) THEN
                jmax = i
              ELSE
                EXIT
              END IF
            END DO
            
            imax = jmax
            
            ifmin = MIN(imin, imax)
            ifmax = MAX(imin, imax)
            n_reg = SIZE(impact(ifmin:ifmax))
            
            IF (n_reg > 50) THEN

              ALLOCATE(K(n_reg, 2))
              ALLOCATE(Y(n_reg))

              K(:,1) = 1.0_wp
              K(:,2) = LOG(ba(ifmin:ifmax))
              Y(:)   = LOG(bangle_nh(ifmin:ifmax))

              CALL ropp_pp_regression(K, Y, rf)

              d = rf(1)**2 + (rf(2)-1.0_wp)**2

              IF (d < dmin) THEN
                dmin = d
                bangle_MSIS(:) = ba(:)
                month = imonth
                msislat = scanlat
                msislon = scanlon
              ENDIF
              
              DEALLOCATE(K)
              DEALLOCATE(Y)

            ENDIF

          ENDIF

        ENDDO
      ENDDO

      MSIS_read = .false.     ! Reset MSIS read flag

    ENDDO
    
    CALL ropp_pp_bangle_MSIS(mfile,month,msislat,msislon,impact-config%r_curve,bangle_MSIS)

    WRITE(outstr,'(2X,A,I2,A,F5.0,A,F5.0)') 'Month = ', month, ' Lat = ',   &
                                            msislat, ' Lon = ', msislon
        
    IF (dmin == Huge(dmin)) THEN
      CALL message(msg_warn, "Global MSIS search failed to find best fit profile.")
      CALL message(msg_warn, 'Local MSIS search: ' // outstr) 
    ELSE
      CALL message(msg_diag, 'Global MSIS search: ' // outstr)   
    ENDIF

    DEALLOCATE(ba)
    DEALLOCATE(bangle_nh)
    DEALLOCATE(idx)
    MSIS_read = .false.     ! Reset MSIS read flag

  END SUBROUTINE ropp_pp_search_model_refraction

  
