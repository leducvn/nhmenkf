!
SUBROUTINE rttov_calcrad( &
            & addcosmic,   &
            & chanprof,    &
            & profiles,    &
            & coeffs,      &
            & rad_cosmic,  &
            & rad_skin,    &
            & rad_surfair, &
            & rad_air)
! Description:
! To convert an array of atmospheric temperatures
!   to planck radiances in many channels
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2002, EUMETSAT, All Rights Reserved.
!
! Method: Uses band correction factors to convert T to radiance
!         which have been precomputed for each channel and are read from
!         the RT coefficient file.
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
! 1.0    01/12/2002  New F90 code with structures (P Brunel A Smith)
!                      based on PLNCX of previous RTTOV versions
! 1.1    02/01/2003  Comments added (R Saunders)
! 1.2    26/09/2003  Multiple polarisations (S English)
! 1.3    11/02/2005  Code vectorisation improved (D Dent)
! 1.4    26/01/2007  Removed polarisation (R Saunders)
! 1.5    15/04/2009  User defined ToA. Layers distinct from levels (P.Rayer)
! 1.6    03/11/2009  Transmittances on levels (A Geer)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!
! Declarations:
! Modules used:
! Imported Parameters:
! Imported Type Definitions:
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_Type
  USE parkind1, ONLY : jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : tcosmic
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(profile_Type  ), INTENT(IN)  :: profiles(:)                                     ! profiles
  TYPE(rttov_coef    ), INTENT(IN)  :: coeffs                                          ! coefficients (Planck)
  TYPE(rttov_chanprof), INTENT(IN)  :: chanprof(:)                                     ! Array of channel indices.
  LOGICAL(KIND=jplm)  , INTENT(IN)  :: addcosmic                                       ! switch for adding cosmic background
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_cosmic (size(chanprof)                     )! cosmic background radiance
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_skin   (size(chanprof)                     )! surface skin radiance
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_surfair(size(chanprof)                     )! 2m air radiance
  REAL(KIND=jprb)     , INTENT(OUT) :: rad_air    (profiles(1)%nlevels, size(chanprof))! level radiances
!INTF_END
! radiances are expressed in mw/cm-1/ster/sq.m
! and temperatures in Kelvin
!local variables:
  REAL   (KIND=jprb) :: t_effective                                       ! effective temperature
  REAL   (KIND=jprb) :: t_effective_1(size(chanprof)                     )! effective temperature
  REAL   (KIND=jprb) :: t_effective_2(size(chanprof)                     )! effective temperature
  REAL   (KIND=jprb) :: t_effective_3(profiles(1)%nlevels, size(chanprof))! effective temperature
  INTEGER(KIND=jpim) :: chan     , i, lev                                 ! loop indices
  INTEGER(KIND=jpim) :: nchannels                                         ! Number of radiances computed (channels used * profiles)
!- End of header ------------------------------------------------------
  nchannels = size(chanprof)
  IF (addcosmic) THEN
    DO i = 1, nchannels
      chan          = chanprof(i)%chan
      t_effective   = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * tcosmic
      rad_cosmic(i) = coeffs%planck1(chan) / (Exp(coeffs%planck2(chan) / t_effective) - 1.0_JPRB)
    ENDDO
  ELSE
    rad_cosmic(:) = 0.0_JPRB
  ENDIF
  DO i = 1, nchannels
    chan = chanprof(i)%chan
! point to corresponding profile (temp. for pressure levels, 2m, skin)
! NB on levels
    DO lev = 1, profiles(1)%nlevels
      t_effective_3(lev, i) = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * profiles(chanprof(i)%prof)%t(lev)
    ENDDO
  ENDDO
  DO i = 1, nchannels
    chan             = chanprof(i)%chan
    t_effective_1(i) = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * profiles(chanprof(i)%prof)%skin%t
    t_effective_2(i) = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * profiles(chanprof(i)%prof)%s2m%t
    rad_skin(i)      = coeffs%planck1(chan) / (Exp(coeffs%planck2(chan) / t_effective_1(i)) - 1.0_JPRB)
    rad_surfair(i)   = coeffs%planck1(chan) / (Exp(coeffs%planck2(chan) / t_effective_2(i)) - 1.0_JPRB)
  ENDDO
  DO i = 1, nchannels
    chan = chanprof(i)%chan
    DO lev = 1, profiles(1)%nlevels
      rad_air(lev, i) = coeffs%planck1(chan) / (Exp(coeffs%planck2(chan) / t_effective_3(lev, i)) - 1.0_JPRB)
    ENDDO
  ENDDO
END SUBROUTINE rttov_calcrad
