SUBROUTINE rttov_calcrad_k( &
            & chanprof,      &
            & profiles,      &
            & profiles_k,    &
            & coeffs,        &
            & rad_skin,      &
            & rad_surfair,   &
            & rad_air,       &
            & rad_skin_k,    &
            & rad_surfair_k, &
            & rad_air_k)
! Description:
! K code to convert an array of atmospheric temperatures
!   to planck radiances in many channels
! No K on Rad Cosmic, k = 0.
!
! derivative of Planck function with respect to temperature is
!
!                                     C2 * Nu
!              C1 * C2 * Nu**4 * Exp( ------- )
!                                        T
! B'(T,Nu) = ------------------------------------- dT
!                     (      C2 * Nu       )**2
!               T**2 *( Exp( ------- ) - 1 )
!                     (         T          )
!
!
! which can be reduced to the following, with
!  C1 = C1 * Nu**3
!  C2 = C2 * Nu
!
!              C2 * B(T,Nu) * (C1 + B(T,Nu))
!  B'(T,Nu) =  ----------------------------- dT
!                        C1 * T**2
!
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
! 1.3    29/03/2005  Add end of header comment (J. Cameron)
! 1.4    05/02/2007  Remove polarisation (R Saunders)
! 1.5    15/09/2009  User defined ToA. Layers distinct from levels (P.Rayer)
! 1.6    03/11/2009  Transmittances / optical depths on levels (A Geer)
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
!
  USE rttov_types, ONLY : rttov_chanprof, rttov_coef, profile_type
  USE parkind1, ONLY : jprb
!INTF_OFF
  USE parkind1, ONLY : jpim
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  TYPE(rttov_chanprof), INTENT(IN)    :: chanprof  (:)
  TYPE(profile_Type  ), INTENT(IN)    :: profiles  (:)
  TYPE(profile_Type  ), INTENT(INOUT) :: profiles_k(size(chanprof))
  TYPE(rttov_coef    ), INTENT(IN)    :: coeffs
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin     (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair  (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_air      (profiles(1)%nlevels, size(chanprof))
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_skin_k   (size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_surfair_k(size(chanprof)                     )
  REAL(KIND=jprb)     , INTENT(IN)    :: rad_air_k    (profiles(1)%nlevels, size(chanprof))
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: t_effective
  REAL   (KIND=jprb) :: t_effective_k
  INTEGER(KIND=jpim) :: chan         , i, lev! loop indices
  INTEGER(KIND=jpim) :: nchannels            ! Number of radiances computed (channels used * profiles)
!- End of header --------------------------------------------------------
  nchannels = size(chanprof)
!cdir nodep
  DO i = 1, nchannels
    chan = chanprof(i)%chan
! point to corresponding profile and geometry structures
    t_effective          = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * profiles(chanprof(i)%prof)%skin%t
    t_effective_k        = (coeffs%planck2(chan) * rad_skin(i) * (coeffs%planck1(chan) + rad_skin(i)) /      &
      & (coeffs%planck1(chan) * t_effective ** 2)) * rad_skin_k(i)
    profiles_k(i)%skin%t = profiles_k(i)%skin%t + coeffs%ff_bcs(chan) * t_effective_k
    t_effective          = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * profiles(chanprof(i)%prof)%s2m%t
    t_effective_k        = (coeffs%planck2(chan) * rad_surfair(i) * (coeffs%planck1(chan) + rad_surfair(i)) /      &
      & (coeffs%planck1(chan) * t_effective ** 2)) * rad_surfair_k(i)
    profiles_k(i)%s2m%t  = profiles_k(i)%s2m%t + coeffs%ff_bcs(chan) * t_effective_k
  ENDDO
  DO i = 1, nchannels
    chan = chanprof(i)%chan
    DO lev = 1, profiles(1)%nlevels
      t_effective          = coeffs%ff_bco(chan) + coeffs%ff_bcs(chan) * profiles(chanprof(i)%prof)%t(lev)
      t_effective_k        = (coeffs%planck2(chan) * rad_air(lev, i) * (coeffs%planck1(chan) + rad_air(lev, i)) /      &
        & (coeffs%planck1(chan) * t_effective ** 2)) * rad_air_k(lev, i)
      profiles_k(i)%t(lev) = profiles_k(i)%t(lev) + coeffs%ff_bcs(chan) * t_effective_k
    ENDDO
  ENDDO
END SUBROUTINE rttov_calcrad_k
