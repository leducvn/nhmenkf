!
SUBROUTINE rttov_init_transmission_aux(transmission_aux)
! Description:
! allocation/deallocation of a profile structure
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
!    Copyright 2007, EUMETSAT, All Rights Reserved.
!
! Method:
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!  1.0       01/12/2002  New F90 code
!  1.1       11/10/2007  Add addclouds, addaerosl, init logicals
!                        nullify unused pointers P.Marguinaud
!  1.2       03/11/2009  Transmittances / optical depths on levels (A Geer)
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
  USE rttov_types, ONLY : transmission_Type_aux
  USE parkind1, ONLY : jpim
!INTF_OFF
  USE parkind1, ONLY : jprb
!INTF_ON
  IMPLICIT NONE
! subroutine arguments
! scalar arguments with intent(out):
  TYPE(transmission_type_aux), INTENT(INOUT) :: transmission_aux
!INTF_END
! Local Arrays and Scalars:
!- End of header --------------------------------------------------------

  IF (Associated (transmission_aux%tau_surf))           transmission_aux%tau_surf           = 0._jprb
! IF (Associated (transmission_aux%tau_surf_r))         transmission_aux%tau_surf_r         = 0._jprb ! don't need to initialise!
! IF (Associated (transmission_aux%fac))                transmission_aux%fac                = 0._jprb ! don't need to initialise!
  IF (Associated (transmission_aux%tau_level))          transmission_aux%tau_level          = 0._jprb
  IF (Associated (transmission_aux%od_level))           transmission_aux%od_level           = 0._jprb
  IF (Associated (transmission_aux%od_singlelayer))     transmission_aux%od_singlelayer     = 0._jprb
  IF (Associated (transmission_aux%odsun_singlelayer))  transmission_aux%odsun_singlelayer  = 0._jprb
  IF (Associated (transmission_aux%tausun_surf))        transmission_aux%tausun_surf        = 0._jprb
  IF (Associated (transmission_aux%tausun_level))       transmission_aux%tausun_level       = 0._jprb
  IF (Associated (transmission_aux%od_sfrac))           transmission_aux%od_sfrac           = 0._jprb
! IF (Associated (transmission_aux%od_sfrac_r))         transmission_aux%od_sfrac_r         = 0._jprb
  IF (Associated (transmission_aux%odsun_sfrac))        transmission_aux%odsun_sfrac        = 0._jprb
  IF (Associated (transmission_aux%od_frac_ac))         transmission_aux%od_frac_ac         = 0._jprb
  IF (Associated (transmission_aux%odsun_frac_ac))      transmission_aux%odsun_frac_ac      = 0._jprb
  IF (Associated (transmission_aux%tau_surf_ac))        transmission_aux%tau_surf_ac        = 0._jprb
  IF (Associated (transmission_aux%tau_surf_acsun))     transmission_aux%tau_surf_acsun     = 0._jprb
  IF (Associated (transmission_aux%tau_ref_surf_ac))    transmission_aux%tau_ref_surf_ac    = 0._jprb
  IF (Associated (transmission_aux%tau_ref_surf_acsun)) transmission_aux%tau_ref_surf_acsun = 0._jprb
  IF (Associated (transmission_aux%od_frac_t))          transmission_aux%od_frac_t          = 0._jprb
  IF (Associated (transmission_aux%odsun_frac_t))       transmission_aux%odsun_frac_t       = 0._jprb
  IF (Associated (transmission_aux%tau_surf_t))         transmission_aux%tau_surf_t         = 0._jprb
  IF (Associated (transmission_aux%tausun_surf_t))      transmission_aux%tausun_surf_t      = 0._jprb
  IF (Associated (transmission_aux%tau_ref_surf_t))     transmission_aux%tau_ref_surf_t     = 0._jprb
  IF (Associated (transmission_aux%tausun_ref_surf_t))  transmission_aux%tausun_ref_surf_t  = 0._jprb
! IF (Associated (transmission_aux%refl_norm))          transmission_aux%refl_norm          = 0._jprb


END SUBROUTINE 
