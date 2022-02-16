! $Id: ropp_pp_FFT.f90 2048 2009-04-07 15:45:10Z frhl $

!****s* FFT/ropp_pp_FFT_complex *
!
! NAME
!    ropp_pp_FFT - Compute Fast Fourier Transform of complex data
!                  using the Danielson-Lanczos Lemma
!
! SYNOPSIS
!    call ropp_pp_FFT(data, isign)
!
! DESCRIPTION
!    This subroutine computes the Fast Fourier transform (or its inverse)
!    of complex data using the Danielson-Lanczos Lemma. Based on dfour.f
!    routine provided in Numerical Recipes.
!
! INPUTS
!    complex(wp), dim(:)  :: data      Input complex data signal
!    integer           :: isign     FFT direction
!                                             > 0 - forward discrete FT
!                                             < 0 - inverse discrete FT
!
! OUTPUT
!    complex(wp), dim(:)  :: data      Transformed sequence
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

SUBROUTINE ropp_pp_FFT_complex(data, isign)

!-------------------------------------------------------------------------------
! 1.1 Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal

  IMPLICIT NONE

  COMPLEX(wp), DIMENSION(:), INTENT(INOUT) :: data     ! Complex signal
  INTEGER,                   INTENT(IN)    :: isign    ! FFT direction

  COMPLEX(wp)                              :: temp     ! Temporary storage
  REAL(wp)                                 :: ar, ai   ! Temporary storage
  INTEGER                                  :: nn       ! Array size
  INTEGER                                  :: i,istep,j,m,mmax,n  ! Indices
  REAL(wp)                                 :: theta,wi,wpi,wpr,wr,wtemp
  REAL(wp), PARAMETER                      :: pi = 3.141592653589793238_wp

!-------------------------------------------------------------------------------
! 1.2 Bit-reversal
!-------------------------------------------------------------------------------

  nn = SIZE(data)

  j=1
  DO i=1,nn

    IF (j > i) THEN
      temp = data(j)
      data(j) = data(i)
      data(i) = temp
    ENDIF
    m=nn/2

    DO WHILE((m >= 2) .AND. (j > m))
      j=j-m
      m=m/2
    ENDDO
    j=j+m

  ENDDO

!-------------------------------------------------------------------------------
! 1.3 Danielson-Lanczos algorithm
!-------------------------------------------------------------------------------

  n = 2*nn
  mmax=2

  DO WHILE(n > mmax)

    istep=2*mmax
    theta=6.28318530717959_wp/(isign*mmax)
    wpr=-2.0_wp*sin(0.5_wp*theta)**2
    wpi=sin(theta)
    wr=1.0_wp
    wi=0.0_wp

    DO m=1,mmax,2
      DO i=m,n,istep
        j=i+mmax

        ar = REAL(data((j-1)/2+1))
        ai = AIMAG(data((j-1)/2+1))
        temp = CMPLX(wr*ar-wi*ai, wr*ai+wi*ar, KIND=wp)

        data((j-1)/2+1) = data((i-1)/2+1) - temp
        data((i-1)/2+1) = data((i-1)/2+1) + temp

      ENDDO

      wtemp=wr
      wr=wr*wpr-wi*wpi+wr
      wi=wi*wpr+wtemp*wpi+wi
    ENDDO
    mmax=istep

  ENDDO

END SUBROUTINE ropp_pp_FFT_complex

!****s* FFT/ropp_pp_FFT_real *
!
! NAME
!    ropp_pp_FFT - Compute Fast Fourier Transform of complex data
!                  (real representation) using the Danielson-Lanczos Lemma
!
! SYNOPSIS
!    call ropp_pp_FFT(data, isign)
!
! DESCRIPTION
!    This subroutine computes the Fast Fourier transform (or its inverse)
!    of complex data using the Danielson-Lanczos Lemma. Based on dfour.f
!    routine provided in Numerical Recipes.
!
! INPUTS
!    real(wp), dim(:)  :: data      Input complex data signal
!    integer           :: isign     FFT direction
!                                             > 0 - forward discrete FT
!                                             < 0 - inverse discrete FT
!
! OUTPUT
!    real(wp), dim(:)  :: data      Transformed sequence
!
! NOTES
!   The signal array contains complex data stored as
!     data(1) = real part f(1)
!     data(2) = imag part f(1)
!     data(3) = real part f(2)
!     data(4) = imag part f(2)
!      ...
!     data(2*nn-1) = real part f(nn)
!     data(2*nn)   = imag part f(nn)
!   etc....
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
! 2. Real signal
!-------------------------------------------------------------------------------

SUBROUTINE ropp_pp_FFT_real(data, isign)

!-------------------------------------------------------------------------------
! 2.1 Declarations
!-------------------------------------------------------------------------------

  USE typesizes, ONLY: wp => EightByteReal

  IMPLICIT NONE

  REAL(wp), DIMENSION(:), INTENT(INOUT) :: data        ! Complex signal
  INTEGER,                INTENT(IN)    :: isign       ! FFT direction

  INTEGER                             :: i,istep,j,m,mmax,n   ! Indices
  REAL(wp)                            :: tempi,tempr          ! Temp storage
  REAL(wp)                            :: theta,wi,wpi,wpr,wr,wtemp

!-------------------------------------------------------------------------------
! 2.2 Bit-reversal
!-------------------------------------------------------------------------------

  n=SIZE(data)

  j=1
  DO i=1,n,2

    IF (j > i) THEN
      tempr=data(j)
      tempi=data(j+1)
      data(j)=data(i)
      data(j+1)=data(i+1)
      data(i)=tempr
      data(i+1)=tempi
    ENDIF
    m=n/2

    DO WHILE((m >= 2) .AND. (j > m))
      j=j-m
      m=m/2
    ENDDO
    j=j+m

  ENDDO

!-------------------------------------------------------------------------------
! 2.3 Danielson-Lanczos algorithm
!-------------------------------------------------------------------------------

    mmax=2

    DO WHILE(n > mmax)

      istep=2*mmax
      theta=6.28318530717959_wp/(isign*mmax)
      wpr=-2.0_wp*sin(0.5_wp*theta)**2
      wpi=sin(theta)
      wr=1.0_wp
      wi=0.0_wp

      DO m=1,mmax,2
        DO i=m,n,istep
          j=i+mmax
          tempr=sngl(wr)*data(j)-sngl(wi)*data(j+1)
          tempi=sngl(wr)*data(j+1)+sngl(wi)*data(j)
          data(j)=data(i)-tempr
          data(j+1)=data(i+1)-tempi
          data(i)=data(i)+tempr
          data(i+1)=data(i+1)+tempi
        ENDDO

        wtemp=wr
        wr=wr*wpr-wi*wpi+wr
        wi=wi*wpr+wtemp*wpi+wi
      ENDDO
      mmax=istep

    ENDDO

  END SUBROUTINE ropp_pp_FFT_real




