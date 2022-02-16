! $Id: ropp_pp_fourier_filter.f90 2048 2009-04-07 15:45:10Z frhl $

SUBROUTINE ropp_pp_fourier_filter(data, window)

!****s* FFT/ropp_pp_fourier_filter *
!
! NAME
!    ropp_pp_fourier_filter - Gaussian filter in spectral space
!
! SYNOPSIS
!    call ropp_pp_fourier_filter(data, window)
!
! DESCRIPTION
!    This routine filters an input complex signal by computing its spectrum,
!    filtering the signal in spectral space, and inverse transforming back
!    to the data space
!
! INPUTS
!    complex(wp), dim(:)  :: data      Input complex data signal
!    integer              :: window    Filtering window width (npoints)
!
! OUTPUT
!    complex(wp), dim(:)  :: data      Filtered sequence
!
! REFERENCES
!   W.H. Press, S.A. Teukolsjy, W.T. Vetterling and B.P. Flannery,
!   Numerical Recipes in C - The Art of Scientific Computing.
!   2nd Ed., Cambridge University Press, 1992.
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
  USE ropp_pp, ONLY: ropp_pp_FFT
  USE ropp_pp_constants, ONLY: Pi

  IMPLICIT NONE

  COMPLEX(wp), DIMENSION(:), INTENT(INOUT) :: data     ! Complex signal
  REAL(wp),                  INTENT(IN)    :: window   ! Window width (points)
  
  INTEGER  :: n     ! Number of data
  INTEGER  :: k     ! Array index
  REAL(wp) :: FFk   ! Filter 

!-------------------------------------------------------------------------------
! 2. Calculate spectrum
!-------------------------------------------------------------------------------

  n = SIZE(data)
  CALL ropp_pp_FFT(data, -1)

!-------------------------------------------------------------------------------
! 3. Filter in spectral space
!-------------------------------------------------------------------------------

  DO k=2,n/2
    
   FFk = 0.0_wp
   
   if ( (Pi*(k-1)*window/N)**2/4 < 500.0_wp) then
     FFk    = Exp(-(Pi*(k-1)*window/N)**2/4)
   endif
   data(k)   = data(k)*FFk
   data(N-k+2) = data(N-k+2)*FFk
 End Do
 
 FFk = 0.0_wp
 if ( (Pi*window/2)**2/4 < 500.0_wp ) then
   FFk    = Exp(-(Pi*window/2)**2/4)
 endif
 data(N/2+1) = data(N/2+1)*FFk
 
!-------------------------------------------------------------------------------
! 4. Inverse Fourier transform
!-------------------------------------------------------------------------------

 do k=1,N
   data(k) = data(k)/N
 enddo
 CALL ropp_pp_FFT(data, 1)

END SUBROUTINE ropp_pp_fourier_filter
