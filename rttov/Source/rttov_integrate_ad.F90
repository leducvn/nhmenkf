Subroutine rttov_integrate_ad( &
   addcosmic, opts, maxnstreams, chanprof,                       &
   emissivity,                   emissivity_ad,                  &
   reflectivity,                 reflectivity_ad,                &
   fresnrefl,                    fresnrefl_ad,                   &
   sunglint,                     sunglint_ad,                    &
   sun,                                                          &
   transmission_aux,             transmission_aux_ad,            &
   transmission_scatt_ir_stream, transmission_scatt_ir_stream_ad,&
   profiles,                     profiles_ad,                    &
   aux_prof,                     aux_prof_ad,                    &
   coef,                                                         &
   raytracing,                   raytracing_ad,                  &
   ircld,                        ircld_ad,                       &
   rad,                                                          &
   auxrad,                                                       &
                                 auxrad_stream,                  &
                                 auxrad_stream_ad,               &
                                 rad_ad)

  ! in,                           !inout
!
! Description:
! To perform AD of integration of radiative transfer equation
! in rttov suite, calculating radiances and brightness temperature.
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
! Method:
! Eyre J.R. 1991 A fast radiative transfer model for satellite sounding
! systems.  ECMWF Research Dept. Tech. Memo. 176 (available from the
! librarian at ECMWF).
!
! Saunders R.W., M. Matricardi and P. Brunel 1999 An Improved Fast Radiative
! Transfer Model for Assimilation of Satellite Radiance Observations.
! QJRMS, 125, 1407-1425.
!
! Matricardi, M. 2003 RTIASI-4, a new version of the ECMWF fast radiative
! transfer model for the infrared atmospheric sounding interferometer.
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)  !
!
! Current Code Owner: SAF NWP
!
! History:
! Version   Date     Comment
! -------   ----     -------
!          25/06/91.    Original code.  J.R.EYRE    *ECMWF*
!          21/08/00.    Emissivity and reflectivity handled separately. Steve English
!          31/01/01.    More cloud computations. F. Chevallier
!          23/03/01     New coef. format, new channel numbers (P. Brunel)
!          31/01/01.    More cloud computations. F. Chevallier
!          28/09/01     Cosmic background temp added G.Deblonde
!          18/01/2002   Thread safe (D.Salmond)
!  1.0     01/12/2002   New F90 code with structures (P Brunel A Smith)
!  1.1     02/01/2003   Added comments (R Saunders)
!  1.2     06/05/2003   Init rad%downcld to 0 in section 1 (P  Brunel)
!  1.3     26/09/2003   Modified to allow for multiple polarisations (S English)
!  1.4     10/05/2004   Fixed bug in call to rttov_calcpolarisation (D.Salmond)
!  1.5     06/09/2004   Mods. for Vectorisation (D Salmond ECMWF & B  Carruthers, Cray)
!  1.6     28/02/2005   More improvements to vectorisation (D Dent)
!  1.7     29/03/2005   Add end of header comment (J. Cameron)
!  1.8     03/03/2006   Marco Matricardi (ECMWF):
!             --        IASI capability added.
!             --        Linear in tau approximation for RT equation introduced.
!             --        Solar radiation introduced for IASI and AIRS.
!  1.9     02/01/2007   Corrected some bugs in the old code (R Saunders)
!  1.10    08/02/2007   Removed polarisation indeces (R Saunders)
!  1.11    11/03/2007   Reintroduced overcast radiance (R Saunders)
!  1.12    27/11/2007   Optimised for IBM/NEC (D Salmond)
!  1.13    17/11/2007   Allowed adjoint/K sensitivity to clear BT/radiance (A Geer)
!  1.14    21/12/2007   Added polarimetric option (R. Saunders)
!  1.15    15/08/2009   User defined ToA. Layers distinct from levels. Top-layer
!                       brought into layer looping to shorten code (P.Rayer)
!  1.16    03/11/2009   Transmittances on levels (A Geer)
!  1.17    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
!                       related quantities are now layer arrays (Marco Matricardi).
!  1.18    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  1.19    14/10/2010   Remove rt8_mode (J Hocking)
!  1.20    14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                       should be performed (J Hocking)
!  2.0     14/12/2011   Re-written (D Rundle)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

Use rttov_types, Only : &
   rttov_chanprof, rttov_coef, profile_Type, profile_aux, transmission_type_aux, transmission_scatt_ir_type, sunglint_type, &
   radiance_Type, rttov_options, ircld_type, raytracing_type, radiance_aux

Use parkind1, Only : jpim, jprb, jplm

!INTF_OFF
Use rttov_const, Only : sensor_id_po

Use yomhook, Only : LHOOK, DR_HOOK
!INTF_ON

Implicit None

!subroutine arguments:
Logical(jplm), Intent(in)                       :: addcosmic
Type(rttov_options), Intent(in)                 :: opts
Integer(jpim), Intent(in)                       :: maxnstreams
Type(rttov_chanprof), Intent(in)                :: chanprof(:)
Type(profile_Type), Intent(in)                  :: profiles(:)
Real(jprb), Intent(in)                          :: emissivity(size(chanprof))
Real(jprb), Intent(in)                          :: reflectivity(size(chanprof))
Real(jprb), Intent(in)                          :: fresnrefl(size(chanprof))
Logical(jplm), Intent(in)                       :: sun(size(chanprof))
Type(rttov_coef), Intent(in)                    :: coef
Type(profile_aux) , Intent(in)                  :: aux_prof
Type(transmission_Type_aux), Intent(in)         :: transmission_aux
type(transmission_scatt_ir_type), intent(in)    :: transmission_scatt_ir_stream
Type(ircld_type), intent(in)                    :: ircld
Type(raytracing_type), intent(in)               :: raytracing
Type(sunglint_type), Intent(in)                 :: sunglint
Type(radiance_Type), Intent(in)                 :: rad
Type(radiance_aux), Intent(in)                  :: auxrad
Type(radiance_aux), Intent(in)                  :: auxrad_stream

Type(sunglint_type), Intent(inout)              :: sunglint_ad
Type(radiance_aux), Intent(inout)               :: auxrad_stream_ad
Real(jprb), Intent(inout)                       :: emissivity_ad(size(chanprof))
Real(jprb), Intent(inout)                       :: reflectivity_ad(size(chanprof))
Real(jprb), Intent(inout)                       :: fresnrefl_ad(size(chanprof))
Type(profile_Type), Intent(inout)               :: profiles_ad(size(profiles))
Type(profile_aux), Intent(inout)                :: aux_prof_ad
Type(transmission_Type_aux), Intent(inout)      :: transmission_aux_ad
type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_stream_ad
Type(ircld_type), intent(inout)                 :: ircld_ad
Type(raytracing_type), intent(inout)            :: raytracing_ad
Type(radiance_Type), Intent(inout)              :: rad_ad

!INTF_END

#include "rttov_calcbt_ad.h"
#include "rttov_calcrad_ad.h"

!local variables:
Real(jprb)    :: cfraction(size(chanprof))
Real(jprb)    :: cfraction_ad(size(chanprof))

Real(jprb)    :: rad_air_ad(profiles(1)%nlevels,size(chanprof))
Real(jprb)    :: rad_surfair_ad(size(chanprof))
Real(jprb)    :: rad_skin_ad(size(chanprof))

Real(jprb)    :: meanrad_up_ad(0:maxnstreams,size(chanprof))
Real(jprb)    :: zup_ad(profiles(1)%nlayers,0:maxnstreams,size(chanprof))
Real(jprb)    :: zdown_ad(profiles(1)%nlayers,0:maxnstreams,size(chanprof))
Real(jprb)    :: meanrad_down_ad(0:maxnstreams,size(chanprof))

real(jprb)    :: pfraction(size(chanprof))        ! cloud fraction
logical(jplm) :: sateqsun(profiles(1)%nlayers,size(profiles(:))) ! True where the solar zenith angle equal to observation angle
integer(jpim) :: pol_id(size(chanprof))       ! polarisation index  

integer(jpim) :: iv2lay(size(chanprof)), iv2lev(size(chanprof)), iv3lay(size(chanprof)), iv3lev(size(chanprof))

Integer(jpim) :: i, lay, iprof, ist, nlayers, nlevels, nchannels

Real(jprb)    :: p

Integer(jpim) :: narray(3)
logical(jplm) :: anysun ! Are *any* channels affected by sun (only need to make check once if negative)

REAL(JPRB) :: ZHOOK_HANDLE

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_AD',0_jpim,ZHOOK_HANDLE)
!- End of header --------------------------------------------------------

#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan

!X.  initialisation of local variables

cfraction_ad(:)   = 0._JPRB
rad_surfair_ad(:) = 0._JPRB
rad_skin_ad(:)    = 0._JPRB
rad_air_ad(:,:)   = 0._jprb

zdown_ad(:,:,:)   = 0._jprb
If(opts%addaerosl .OR. opts%addclouds) zup_ad(:,:,:) = 0._jprb

!---------------------------
!0. Initialise useful variables
!---------------------------

nchannels = size(chanprof)
nlayers = profiles(1)%nlayers
nlevels = nlayers + 1

! DAR: narray is a more compact way of passing nchannels, nlevels and maxnstreams to the internal subroutines
! Again module variables would be a more elegant way of solving this problem and would results in less duplicated code
narray(1) = nlevels; narray(2) = maxnstreams; narray(3) = nchannels

do i = 1, nchannels
   cfraction(i) = aux_prof%s(prof)%cfraction
   pfraction(i) = aux_prof%s(prof)%pfraction_surf
enddo

do iprof = 1, size(profiles)
   do lay = 1, nlayers
      if(raytracing%pathsat(lay,iprof) == raytracing%pathsun(lay,iprof)) then
         sateqsun(lay,iprof) = .true.
      else
         sateqsun(lay,iprof) = .false.
      endif
   enddo
enddo
Do i = 1, nchannels 
! case-1: surf lies above lev=nlevels
   iv3lev(i) = aux_prof%s(prof)%nearestlev_surf - 1   ! lowest lev above surf
! case-2: surf lies below lev=nlevels
   if (pfraction(i) < 0.0_JPRB) iv3lev(i) = iv3lev(i) + 1  ! iv3lev=iv2lev=lowest lev above surf

   iv2lev(i) = aux_prof%s(prof)%nearestlev_surf       ! highest lev below surf
   iv2lay(i) = iv2lev(i) - 1                          ! same layer as that numbered by  iv2 in RTTOV-9
   iv3lay(i) = iv3lev(i) - 1                          ! same layer as that numbered by  iv3 in RTTOV-9
End Do

if(coef%id_sensor == sensor_id_po) then
   do i = 1, nchannels
      pol_id(i) = coef%fastem_polar(chan) + 1_jpim
   enddo
Else
   pol_id(:) = 0_jpim
Endif 

anysun = any(sun)

!-----------------------------------------------
!9. convert radiances to brightness temperatures
!-----------------------------------------------
if (opts%switchrad)  then
   Call rttov_calcbt_ad(chanprof, coef, rad, &! in
                        rad_ad)    ! inout
endif

!---------------------------
!8. calculate total radiance (cloudy case)
!---------------------------
If(opts%addaerosl .OR. opts%addclouds) then
   rad_ad%cloudy(1:nchannels) = rad_ad%total(1:nchannels)
!----------------------------------------
!7. calculate complex cloudy radiances
!----------------------------------------
   Do i = 1, nchannels
      rad_ad%clear(i)        = rad_ad%clear(i) + rad_ad%cloudy(i) * ircld%XSTRCLR(prof)
      ircld_ad%XSTRCLR(prof) = ircld_ad%XSTRCLR(prof) + rad_ad%cloudy(i) * rad%clear(i)

      do ist = ircld%nstream(prof), 1, -1 !reverse loop nec. due to data dependence
         auxrad_stream_ad%cloudy(ist,i) = auxrad_stream_ad%cloudy(ist,i) + &
                                          rad_ad%cloudy(i) * (ircld%xstr(ist+1,prof) - ircld%xstr(ist,prof))

         ircld_ad%xstr(ist+1,prof)      = ircld_ad%xstr(ist+1,prof) + &
                                          rad_ad%cloudy(i) * auxrad_stream%cloudy(ist,i)

         ircld_ad%xstr(ist,prof)        = ircld_ad%xstr(ist,prof) - &
                                          rad_ad%cloudy(i) * auxrad_stream%cloudy(ist,i)
      Enddo
   End do  
Else
   if(.not. opts%addpc) then
      rad_ad%clear(:)  = rad_ad%clear(:) + (1.0_jprb - cfraction(:)) * rad_ad%total(:)
      rad_ad%cloudy(:) = cfraction(:) * rad_ad%total(:)
      cfraction_ad(:)  = (rad%cloudy(:) - rad%clear(:)) * rad_ad%total(:)
   else
      rad_ad%clear(:)  = rad_ad%clear(:) + rad_ad%total(:)
   endif

!6. Calculate "simple" cloudy radiances
!--------------------------------------------
!6.2 Interpolate to given cloud-top pressures
!--------------------------------------------
   do i = 1, nchannels
      lay = aux_prof%s(prof)%nearestlev_ctp - 1
      p   = aux_prof%s(prof)%pfraction_ctp

      rad_ad%overcast(lay,i)            = rad_ad%overcast(lay,i) + & !first use of rad_ad%overcast
                                         (1_jprb - p) * rad_ad%cloudy(i)
     
      rad_ad%overcast(lay-1,i)          = rad_ad%overcast(lay-1,i) + & !first use of rad_ad%overcast
                                          p * rad_ad%cloudy(i)
     
      aux_prof_ad%s(prof)%pfraction_ctp = aux_prof_ad%s(prof)%pfraction_ctp + &
                                         (rad%overcast(lay-1,i) - rad%overcast(lay,i)) * rad_ad%cloudy(i)
   End Do
! rad_ad%cloudy done

!---------------------------------------------------
!6. calculate overcast radiances
!---------------------------------------------------
!---------------
!6.1 Upward part
!---------------

ist = 0_jpim
!cdir nodep
Do i = 1, nchannels
   do lay = 1, nlayers
! overcast radiances at given cloud top
      zup_ad(lay,ist,i) = &!zup_ad(lay,ist,i) + &
                          rad_ad%overcast(lay,i) !first use of zup_ad for ist = 0
      rad_air_ad(lay+1,i) = &!rad_air_ad(lay+1,i) + &
                            rad_ad%overcast(lay,i) * transmission_aux%tau_level(lay+1,ist,i) ! first use of rad_air_ad
      transmission_aux_ad%tau_level(lay+1,ist,i) = transmission_aux_ad%tau_level(lay+1,ist,i) + &
                                                   auxrad%air(lay+1,i) * rad_ad%overcast(lay,i) 
   end do
enddo
! rad_ad%overcast done      
! rad_ad%clear done

 ! aux_prof_ad%cfraction write out now
 Do i = 1, nchannels
    aux_prof_ad%s(prof)%cfraction = aux_prof_ad%s(prof)%cfraction + cfraction_ad(i)
 End Do
Endif

!--------------------------------
!5. cosmic temperature correction
!--------------------------------

!calculate planck function corresponding to tcosmic=2.7k
!deblonde tcosmic for microwave sensors only

If (addcosmic) Then
   Do i = 1, nchannels
      ist = 0_jpim
      reflectivity_ad(i) = reflectivity_ad(i) + &
                           rad_ad%clear(i) * auxrad%cosmic(i) * &
                           transmission_aux%tau_surf(ist,i)**2_jpim
            
      transmission_aux_ad%tau_surf(ist,i) = & !transmission_aux_ad%tau_surf(ist,i) + & ! first use of transmission_aux_ad%tau_surf
                           rad_ad%clear(i) * 2.0_jprb * reflectivity(i) * auxrad%cosmic(i) * &
                           transmission_aux%tau_surf(ist,i)
            
      do ist = 1, ircld%nstream(prof)!rev loop not nec.
         reflectivity_ad(i) = reflectivity_ad(i) + &
                              auxrad_stream_ad%cloudy(ist,i) * auxrad%cosmic(i) * transmission_aux%tau_surf(ist,i)**2
            
         transmission_aux_ad%tau_surf(ist,i) = &! transmission_aux_ad%tau_surf(ist,i) + & !first use of transmission_aux_ad%tau_surf
                                               auxrad_stream_ad%cloudy(ist,i) * 2.0_jprb * &
                                               reflectivity(i) * auxrad%cosmic(i) * transmission_aux%tau_surf(ist,i)
      enddo
   enddo
Endif

!----------------------------------------------------------------------------------------
! 4. solar contribution
!----------------------------------------------------------------------------------------
if(coef%fmv_model_ver == 9) then
   if(anysun) call solar_contribution_ad(transmission_aux, transmission_aux_ad, &
                                         ircld, sunglint, sunglint_ad, fresnrefl, fresnrefl_ad, &
                                         reflectivity, reflectivity_ad, &
                                         coef, profiles, chanprof, sun, narray, auxrad_stream_ad, rad_ad)
endif

!-----------------------
!3. surface contribution
!-----------------------

Do i = 1, nchannels
  ist = 0_jpim
  emissivity_ad(i) = emissivity_ad(i) + &
                     rad_ad%clear(i) * auxrad%skin(i) * transmission_aux%tau_surf(ist,i)

  transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) + &
                                        rad_ad%clear(i) * auxrad%skin(i) * emissivity(i)

  rad_skin_ad(i) = rad_skin_ad(i) + &
                   rad_ad%clear(i) * emissivity(i) * transmission_aux%tau_surf(ist,i)

  do ist = 1, ircld%nstream(prof) !rev loop not nec.
     emissivity_ad(i) = emissivity_ad(i) + &
                        auxrad_stream_ad%cloudy(ist,i) * auxrad%skin(i) * transmission_aux%tau_surf(ist,i)

     transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) + &
                                           auxrad_stream_ad%cloudy(ist,i) * auxrad%skin(i) * emissivity(i)

     rad_skin_ad(i) = rad_skin_ad(i) + &
                      auxrad_stream_ad%cloudy(ist,i) * emissivity(i) * transmission_aux%tau_surf(ist,i)
  enddo


Enddo

!-------------------------------------------------------------------------------
!2a calculate near-surface layer contribution
!-------------------------------------------------------------------------------
!cdir nodep
Do i = 1, nchannels
   ist = 0_jpim
! DAR - upclear not added here because it is never set.
   meanrad_up_ad(ist,i) = &! meanrad_up_ad(ist,i) + &
                          rad_ad%clear(i)
           
   meanrad_down_ad(ist,i) = &! meanrad_down_ad(ist,i) + &
                            rad_ad%clear(i) * reflectivity(i) * transmission_aux%tau_surf(ist,i)**2

   transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) + &
                                         rad_ad%clear(i) * 2._JPRB * reflectivity(i) * &
                                         transmission_aux%tau_surf(ist,i) * auxrad_stream%meanrad_down(ist,i)

   reflectivity_ad(i) = reflectivity_ad(i) + &
                        rad_ad%clear(i) * transmission_aux%tau_surf(ist,i)**2 * auxrad_stream%meanrad_down(ist,i)

   do ist = 1, ircld%nstream(prof)
      ! add upward and downward parts

      meanrad_up_ad(ist,i) = & ! first use of meanrad_up_ad
                             zup_ad(iv2lay(i),ist,i) + &
                             auxrad_stream_ad%cloudy(ist,i)

      meanrad_down_ad(ist,i) = & ! first use of meanrad_down_ad
                               auxrad_stream_ad%cloudy(ist,i) * transmission_aux%tau_surf(ist,i)**2 * reflectivity(i)

      transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) + &
                                            auxrad_stream_ad%cloudy(ist,i) * &
                                            2._JPRB * transmission_aux%tau_surf(ist,i) * reflectivity(i) * &
                                            auxrad_stream%meanrad_down(ist,i)

      reflectivity_ad(i) = reflectivity_ad(i) + &
                           auxrad_stream_ad%cloudy(ist,i) * transmission_aux%tau_surf(ist,i)**2 * &
                           auxrad_stream%meanrad_down(ist,i)
   enddo
enddo

if((opts%addaerosl .OR. opts%addclouds) .and. anysun) &
     call solar_scattering_near_surf_ad(transmission_aux, transmission_aux_ad, &
                                        ircld, raytracing, raytracing_ad, coef, transmission_scatt_ir_stream_ad, &
                                        sun, sateqsun, iv3lay, iv3lev, pol_id, chanprof, narray, pfraction, &
                                        auxrad_stream, meanrad_up_ad, meanrad_down_ad)
 
call calc_near_surf_contribution_ad(transmission_aux, transmission_aux_ad, auxrad, ircld, & 
     chanprof, narray, iv3lay, iv3lev, pol_id, & 
     rad_surfair_ad, rad_air_ad, zup_ad, zdown_ad, &
     meanrad_up_ad, meanrad_down_ad)

 if(opts%addaerosl .OR. opts%addclouds) then
    Do i = 1, nchannels
       if(sun(i))then
          do ist = 0, ircld%nstream(prof)
             Do lay = 2, nlayers !rev loop not nec
                if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) then
                   zdown_ad(lay,ist,i) = 0._jprb
                endif
             enddo
          enddo
       endif
    enddo
 endif

 Do i = 1, nchannels
    do ist = ircld%nstream(prof),0,-1!rev loop nec due to data dependence
       Do lay = profiles(1) % nlayers, 2, -1
          zup_ad(lay-1,ist,i) = zup_ad(lay-1,ist,i) + zup_ad(lay,ist,i)
          zdown_ad(lay-1,ist,i) = zdown_ad(lay-1,ist,i) + zdown_ad(lay,ist,i)
       enddo
    enddo
 enddo

!-------------------------------------
!2. calculate atmospheric contribution
!-------------------------------------
if((opts%addaerosl .OR. opts%addclouds) .and. anysun) &
     call solar_scattering_air_ad(transmission_aux, transmission_aux_ad, &
                                      ircld, raytracing, raytracing_ad, coef, transmission_scatt_ir_stream_ad, &
                                      sun, sateqsun, chanprof, narray, auxrad_stream, & 
                                      zup_ad, zdown_ad)

call calc_atmospheric_radiance_ad(transmission_aux, transmission_aux_ad, auxrad, ircld, chanprof, narray, &
                                  rad_air_ad, zup_ad, zdown_ad)

Call rttov_calcrad_ad(chanprof, profiles, &! in
     profiles_ad, &! inout
     coef, auxrad%skin, auxrad%surfair, auxrad%air, rad_skin_ad, rad_surfair_ad, rad_air_ad)  ! in

IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_AD',1_jpim,ZHOOK_HANDLE)

contains

subroutine calc_atmospheric_radiance_ad(transmission_aux, transmission_aux_ad, auxrad, ircld, chanprof, narray, &
     rad_air_ad, zup_ad, zdown_ad)
  
  use rttov_types, Only : transmission_type_aux, radiance_aux, ircld_type, rttov_chanprof
  use rttov_const, only : min_tau
  use parkind1, only : jpim, jprb

  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(radiance_aux), intent(in)                :: auxrad
  type(ircld_type),     intent(in)              :: ircld
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  integer(jpim), intent(in)                     :: narray(:)

  Real(jprb), intent(inout)                     :: rad_air_ad(:,:)
  Real(jprb), intent(inout)                     :: zup_ad(:,0:,:), zdown_ad(:,0:,:)
  type(transmission_type_aux), intent(inout)    :: transmission_aux_ad
  
  integer(jpim) :: nlevels, nlayers, nchannels
  integer(jpim) :: i, ist, lay, lev, levm1

  Real(jprb) :: tau_prod
  Real(jprb) :: b1_ad, b2_ad, b3_ad
  Real(jprb) :: temp

  !unpack narray
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)

#define tau_lev transmission_aux%Tau_level(lay+1,ist,i)
#define tau_levm1 transmission_aux%Tau_level(lay,ist,i)
#define tau_lev_r transmission_aux%Tau_level_r(lay+1,ist,i)
#define tau_levm1_r transmission_aux%Tau_level_r(lay,ist,i)
#define tausun_lev transmission_aux%Tausun_level(lay+1,ist,i)
#define tausun_levm1 transmission_aux%Tausun_level(lay,ist,i)

#define dtau_lay (tau_levm1 - tau_lev)
#define daux_lay (auxrad%air(lev,i) - auxrad%air(levm1,i))
#define B1_2 auxrad%air(levm1,i) * dtau_lay
#define B2_2 daux_lay * tau_lev
#define B3_2 daux_lay * dtau_lay * transmission_aux%od_singlelayer_r(lay,ist,i)

  Do i = 1, nchannels
     do ist = 0, ircld%nstream(prof)
        Do lay = 1, nlayers
           levm1 = lay
           lev = lay + 1
           if(transmission_aux%anynegtau > 0_jprb) then

              if(lay > 1 .and. tau_levm1 >= 0.and. tau_lev < 0) then

                 If (tau_lev > min_tau) Then
                    tau_prod = tau_lev_r * tau_levm1_r
                    zup_ad(lay,ist,i) = zup_ad(lay,ist,i) + &
                         zdown_ad(lay,ist,i) * tau_prod

                    transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) - &
                         zdown_ad(lay,ist,i) * &
                         tau_levm1 * (tau_prod*tau_prod) * &
                         0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i)) * (tau_levm1 - tau_lev )

                    transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) - &
                         zdown_ad(lay,ist,i) * &
                         tau_lev * (tau_prod*tau_prod) * &
                         0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i)) * (tau_levm1 - tau_lev)
                 End If

                 rad_air_ad(lev,i) = rad_air_ad(lev,i) + &
                      zup_ad(lay,ist,i) * 0.5_JPRB * (tau_levm1 - tau_lev)

                 rad_air_ad(levm1,i) = rad_air_ad(levm1,i) + &
                      zup_ad(lay,ist,i) * 0.5_JPRB * (tau_levm1 - tau_lev)

                 transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) + &
                      zup_ad(lay,ist,i) * 0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i))

                 transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) - &
                      zup_ad(lay,ist,i) * 0.5_JPRB * (auxrad%air(lev,i) + auxrad%air(levm1,i))
              else
                 temp = transmission_aux%fac(1,lay,ist,i) * transmission_aux%fac(2,lay+1,ist,i) * zdown_ad(lay,ist,i)

                 B1_AD = temp * (tau_levm1_r * tau_lev_r)
                 B2_AD = temp * (tau_lev_r)**2
                 B3_AD = -temp * (tau_levm1_r * tau_lev_r)

                 B1_AD = B1_AD + zup_ad(lay,ist,i)
                 B2_AD = B2_AD - zup_ad(lay,ist,i)
                 B3_AD = B3_AD + zup_ad(lay,ist,i)

                 transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) + &
                      transmission_aux%fac(1,lay,ist,i) * (&
                      B1_AD * auxrad%air(levm1,i) - &
                      temp * &
                     (B1_2 - &
                      B3_2) * &
                     (tau_lev_r * tau_levm1_r**2) + &
                      B3_AD * daux_lay * transmission_aux%od_singlelayer_r(lay,ist,i))

                 transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) + &
                      transmission_aux%fac(1,lay,ist,i) * (&
                      -temp * &
                     (B1_2 - &
                      B3_2) * &
                     (tau_lev_r**2 * tau_levm1_r) - &
                      B3_AD * daux_lay * transmission_aux%od_singlelayer_r(lay,ist,i) - &
                      temp * 2.0_jprb * B2_2 * (tau_lev_r**3) - &
                      B1_AD * auxrad%air(levm1,i) + &
                      B2_AD * daux_lay)

                 transmission_aux_ad%od_singlelayer(lay,ist,i) = transmission_aux_ad%od_singlelayer(lay,ist,i) - &
                      transmission_aux%fac(1,lay,ist,i) * &
                      B3_AD * B3_2 * &
                      transmission_aux%od_singlelayer_r(lay,ist,i) 

                 rad_air_ad(levm1,i) = rad_air_ad(levm1,i) + transmission_aux%fac(1,lay,ist,i) * (&
                     -B3_AD * dtau_lay * transmission_aux%od_singlelayer_r(lay,ist,i) + &
                      B1_AD * dtau_lay - &
                      B2_AD * tau_lev)

                 rad_air_ad(lev,i) = rad_air_ad(lev,i) + transmission_aux%fac(1,lay,ist,i) * (&
                      B3_AD * dtau_lay * transmission_aux%od_singlelayer_r(lay,ist,i) + &
                      B2_AD * tau_lev)
              endif
           else
              temp = transmission_aux%fac(1,lay,ist,i) * transmission_aux%fac(2,lay+1,ist,i) * zdown_ad(lay,ist,i)

              B1_AD = temp * (tau_levm1_r * tau_lev_r)
              B2_AD = temp * (tau_lev_r)**2
              B3_AD = -temp * (tau_levm1_r * tau_lev_r)

              B1_AD = B1_AD + zup_ad(lay,ist,i)
              B2_AD = B2_AD - zup_ad(lay,ist,i)
              B3_AD = B3_AD + zup_ad(lay,ist,i)

              transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) + &
                   transmission_aux%fac(1,lay,ist,i) * (&
                   B1_AD * auxrad%air(levm1,i) - &
                   temp * &
                  (B1_2 - &
                   B3_2) * &
                  (tau_lev_r * tau_levm1_r**2) + &
                   B3_AD * daux_lay * transmission_aux%od_singlelayer_r(lay,ist,i))

              transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) + &
                   transmission_aux%fac(1,lay,ist,i) * (&
                   -temp * &
                  (B1_2 - &
                   B3_2) * &
                  (tau_lev_r**2 * tau_levm1_r) - &
                   B3_AD * daux_lay * transmission_aux%od_singlelayer_r(lay,ist,i) - &
                   temp * 2.0_jprb * B2_2 * (tau_lev_r**3) - &
                   B1_AD * auxrad%air(levm1,i) + &
                   B2_AD * daux_lay)
              
              transmission_aux_ad%od_singlelayer(lay,ist,i) = transmission_aux_ad%od_singlelayer(lay,ist,i) - &
                   transmission_aux%fac(1,lay,ist,i) * B3_AD * &
                   B3_2 * &
                   transmission_aux%od_singlelayer_r(lay,ist,i) 
              
              rad_air_ad(levm1,i) = rad_air_ad(levm1,i) + transmission_aux%fac(1,lay,ist,i) * (&
                   -B3_AD * dtau_lay * transmission_aux%od_singlelayer_r(lay,ist,i) + &
                   B1_AD * dtau_lay - &
                   B2_AD * tau_lev)
              
              rad_air_ad(lev,i) = rad_air_ad(lev,i) + transmission_aux%fac(1,lay,ist,i) * (&
                   B3_AD * dtau_lay * transmission_aux%od_singlelayer_r(lay,ist,i) + &
                   B2_AD * tau_lev)
           endif
        Enddo
     End Do
  End Do
end subroutine calc_atmospheric_radiance_ad

subroutine calc_near_surf_contribution_ad(transmission_aux, transmission_aux_ad, auxrad, ircld, & 
                                          chanprof, narray, iv3lay, iv3lev, pol_id, & 
                                          rad_surfair_ad, rad_air_ad, zup_ad, zdown_ad, &
                                          meanrad_up_ad, meanrad_down_ad)

  use rttov_const, only : min_tau, min_od
  use rttov_types, Only : rttov_chanprof, ircld_type, transmission_type_aux, radiance_aux
  use parkind1, only : jpim, jprb

  Implicit None

  type(transmission_type_aux), intent(in)    :: transmission_aux
  type(radiance_aux), intent(in)             :: auxrad
  type(ircld_type), intent(in)               :: ircld
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  integer(jpim), intent(in)                  :: narray(:)
  integer(jpim), intent(in)                  :: iv3lay(:), iv3lev(:)
  integer(jpim), intent(in)                  :: pol_id(:)

  real(jprb), intent(inout)                  :: meanrad_up_ad(0:,:), meanrad_down_ad(0:,:)

  type(transmission_type_aux), intent(inout) :: transmission_aux_ad
  real(jprb), intent(inout)                  :: rad_surfair_ad(:), rad_air_ad(:,:)
  real(jprb), intent(inout)                  :: zup_ad(1:,0:,:), zdown_ad(1:,0:,:)

!local
  integer(jpim) :: nlayers, nlevels, nstreams, nchannels

  real(jprb)    :: b1_3, b2_3, b3_3
  real(jprb)    :: b1_ad, b2_ad, b3_ad

  integer(jpim) :: i, lev, lay, ist

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)  
 
  Do i = 1, nchannels
     lay = iv3lay(i)
     lev = iv3lev(i)
   
     do ist = 0, ircld%nstream(prof) !rev loop not nec
        ! assume there is no atmospheric source term for 3rd/4th stokes vector elements
        if (pol_id(i) >= 6_jpim) meanrad_up_ad(ist,i) = 0.0_jprb

        zup_ad(lay,ist,i) = zup_ad(lay,ist,i) + &
                            meanrad_up_ad(ist,i)

        zdown_ad(lay,ist,i) = zdown_ad(lay,ist,i) + &
                              meanrad_down_ad(ist,i)

        if(transmission_aux%tau_surf(ist,i) < 0._JPRB) then
           if(transmission_aux%tau_level(lev,ist,i) >= 0._JPRB) then ! nearly always true?
              
              rad_surfair_ad(i) = rad_surfair_ad(i) + &
                                  0.5_JPRB * meanrad_up_ad(ist,i) * &
                                 (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i))

              rad_air_ad(lev,i) = rad_air_ad(lev,i) + &
                                  0.5_JPRB * meanrad_up_ad(ist,i)* &
                                 (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i))

              transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) + &
                                                         0.5_JPRB * meanrad_up_ad(ist,i) * &
                                                        (auxrad%surfair(i) + auxrad%air(lev,i))

              transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) - &
                                                    0.5_JPRB * meanrad_up_ad(ist,i) * &
                                                   (auxrad%surfair(i) + auxrad%air(lev,i))
           endif
        else
           if(transmission_aux%od_sfrac(ist,i) < min_od .OR. &
              ((transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i)) < min_od)) THEN
           else
              IF (transmission_aux%tau_surf(ist,i) > min_tau) THEN
!direct calc
                 B1_3 = auxrad%air(lev,i) * & 
                       (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i))

                 B2_3 = (auxrad%surfair(i) - auxrad%air(lev,i)) * transmission_aux%tau_surf(ist,i)

                 B3_3 = (auxrad%surfair(i) - auxrad%air(lev,i)) * &
                        (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i)) * &
                        (transmission_aux%od_sfrac_r(ist,i))
                 !first use B1,b2,b3c
                 B1_AD = meanrad_down_ad(ist,i) * (transmission_aux%tau_level_r(lev,ist,i) * &
                   transmission_aux%tau_surf_r(ist,i)) 
                 B3_AD = -meanrad_down_ad(ist,i) * (transmission_aux%tau_level_r(lev,ist,i) * &
                   transmission_aux%tau_surf_r(ist,i))
                 B2_AD = meanrad_down_ad(ist,i) * (transmission_aux%tau_surf_r(ist,i))**2

                 transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) - &
                       meanrad_down_ad(ist,i) * 2_jprb * B2_3 * transmission_aux%tau_surf_r(ist,i)**3_jpim

                 transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) - &
                       meanrad_down_ad(ist,i) * (B1_3 - B3_3) * &
                       transmission_aux%tau_level_r(lev,ist,i) * transmission_aux%tau_surf_r(ist,i)**2_jpim

                 transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) - &
                       meanrad_down_ad(ist,i) * (B1_3 - B3_3) * &
                       transmission_aux%tau_level_r(lev,ist,i)**2_jpim * transmission_aux%tau_surf_r(ist,i)
              ELSE
                 B1_AD = 0.0_JPRB
                 B2_AD = 0.0_JPRB
                 B3_AD = 0.0_JPRB
              ENDIF

              B1_AD = B1_AD + meanrad_up_ad(ist,i)
              B2_AD = B2_AD - meanrad_up_ad(ist,i)
              B3_AD = B3_AD + meanrad_up_ad(ist,i)

              transmission_aux_ad%od_sfrac(ist,i) = transmission_aux_ad%od_sfrac(ist,i) - &
                                                    B3_AD * (auxrad%surfair(i)-auxrad%air(lev,i)) * &
                                                   (transmission_aux%tau_level(lev,ist,i) - &
                                                   transmission_aux%tau_surf(ist,i)) * &
                                                   (transmission_aux%od_sfrac_r(ist,i)**2)

              transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) - &
                                                    B3_AD * (auxrad%surfair(i)-auxrad%air(lev,i)) * &
                                                    transmission_aux%od_sfrac_r(ist,i) + &
                                                    B2_AD * (auxrad%surfair(i)-auxrad%air(lev,i)) - &
                                                    B1_AD * auxrad%air(lev,i)
              
              transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) + &
                                                         B3_AD * (auxrad%surfair(i) - auxrad%air(lev,i)) * &
                                                         transmission_aux%od_sfrac_r(ist,i) + &
                                                         B1_AD * auxrad%air(lev,i)

              rad_air_ad(lev,i) = rad_air_ad(lev,i) - &
                   B3_AD * (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i)) * &
                   transmission_aux%od_sfrac_r(ist,i) - &
                   B2_AD * transmission_aux%tau_surf(ist,i) + &
                   B1_AD * (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i))

              rad_surfair_ad(i) = rad_surfair_ad(i) + &
                   B3_AD * (transmission_aux%tau_level(lev,ist,i) - transmission_aux%tau_surf(ist,i)) * &
                   transmission_aux%od_sfrac_r(ist,i) + &
                   B2_AD * transmission_aux%tau_surf(ist,i)
           endif
        endif
     Enddo
  End do
end subroutine calc_near_surf_contribution_ad

subroutine solar_scattering_air_ad(transmission_aux, transmission_aux_ad, &
                                       ircld, raytracing, raytracing_ad, coef, transmission_scatt_ir_stream_ad, &
                                       sun, sateqsun, chanprof, narray, auxrad_stream, & 
                                       zup_ad, zdown_ad)

  use rttov_const, only : min_tau, z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, &
                          rttov_chanprof, rttov_coef, transmission_scatt_ir_type, raytracing_type
  use parkind1, only : jpim, jprb, jplm 

  Implicit None
  
  type(transmission_type_aux), intent(in)         :: transmission_aux
  type(ircld_type),     intent(in)                :: ircld
  type(raytracing_type), intent(in)               :: raytracing
  type(rttov_coef),     intent(in)                :: coef
  logical(jplm), intent(in)                       :: sun(:)
  logical(jplm), intent(in)                       :: sateqsun(:,:)
  type(rttov_chanprof), intent(in)                :: chanprof(:)
  integer(jpim), intent(in)                       :: narray(:)
  type(radiance_aux), intent(in)                  :: auxrad_stream

  type(transmission_type_aux), intent(inout)      :: transmission_aux_ad
  type(raytracing_type), intent(inout)            :: raytracing_ad
  type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_stream_ad
  real(jprb),intent(inout)                        :: zup_ad(:,0:,:), zdown_ad(:,0:,:)

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist, lay, lev, levm1
  
  Real(jprb) :: fac_2_ad(7)
  Real(jprb) :: temp(3,narray(1)-1)

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)  

  Do i = 1, nchannels
     if(sun(i)) then

        do ist = 0, ircld%nstream(prof)!rev loop not nec

           lay = 1
           if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) then
              zdown_ad(lay,ist,i) = 0._jprb
           endif

           temp(1,1:nlayers) = auxrad_stream%fac6_2(1:nlayers,ist,i) * auxrad_stream%fac2_2(1:nlayers,ist,i) * &
                               zdown_ad(1:nlayers,ist,i)        
           temp(2,1:nlayers) = transmission_aux%tau_level_r(1:nlayers,ist,i) * &
             transmission_aux%tau_level_r(2:nlayers+1,ist,i)
           
           Do lay = 1, nlayers !rev loop not nec

              levm1 = lay
              lev   = lay + 1
             
              fac_2_ad(:) = 0.0_jprb ! DAR: if errors then may need to zero everything. Not just 2 and 6

              if(lay == 1 .or. tau_lev > min_tau) then
                 !-------------------Downward single scattering of the solar beam--------------------------
                 if(.not. sateqsun(lay,prof)) then
                    temp(3,lay) = (auxrad_stream%fac5_2(lay,ist,i) - auxrad_stream%fac4_2(lay,ist,i)) * &
                                   auxrad_stream%fac7_2(lay,ist,i)

                    if(lay > 1_jpim) then
                       transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) - &
                            temp(1,lay) * temp(2,lay)**2 * temp(3,lay) * tausun_levm1 * tau_lev

                       transmission_aux_ad%tausun_level(levm1,ist,i) = transmission_aux_ad%tausun_level(levm1,ist,i) + &
                            temp(1,lay) * temp(2,lay) * temp(3,lay)
                    endif
                   
                    fac_2_ad(4) = &!fac4_2_ad 
                                - temp(1,lay) * temp(2,lay) * &
                                auxrad_stream%fac7_2(lay,ist,i) * transmission_aux%tausun_level(lay,ist,i)
                    
                    fac_2_ad(5) = &!fac5_2_ad + &
                                temp(1,lay) * temp(2,lay) * &
                                auxrad_stream%fac7_2(lay,ist,i) * transmission_aux%tausun_level(lay,ist,i)
                    
                    fac_2_ad(7) = &!fac7_2_ad + &
                                temp(1,lay) * temp(2,lay) * &
                               (auxrad_stream%fac5_2(lay,ist,i) - auxrad_stream%fac4_2(lay,ist,i)) * &
                                transmission_aux%tausun_level(lay,ist,i)
                    
                    raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
                         fac_2_ad(7) * raytracing%pathsun(lay,prof) / &
                        (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2

                    raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) - &
                         fac_2_ad(7) * raytracing%pathsat(lay,prof) / &
                        (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2
                 else
                    temp(3,lay) = auxrad_stream%fac4_2(lay,ist,i) * transmission_aux%odsun_singlelayer(lay,ist,i)

                    if(lay > 1) then
!                       transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) - &
!                            temp(1,lay) * temp(2,lay) * temp(3,lay) !DAR 3108 think this should be temp22 * temp3 * tau_lev as above from integrate_tl
                       transmission_aux_ad%tau_level(levm1,ist,i) = transmission_aux_ad%tau_level(levm1,ist,i) - &
                            temp(1,lay) * temp(2,lay)**2_jpim * temp(3,lay) * tau_lev 
                       
                       transmission_aux_ad%tausun_level(levm1,ist,i) = transmission_aux_ad%tausun_level(levm1,ist,i) + &
                            temp(1,lay) * temp(2,lay) * temp(3,lay)
                    endif
                    
                    transmission_aux_ad%odsun_singlelayer(lay,ist,i) = transmission_aux_ad%odsun_singlelayer(lay,ist,i) + &
                         temp(1,lay) * temp(2,lay) * auxrad_stream%fac4_2(lay,ist,i) * tausun_levm1
                    
                    fac_2_ad(4) = &!fac4_2_ad + &
                                temp(1,lay) * temp(2,lay) * &
                                transmission_aux%tausun_level(lay,ist,i) * transmission_aux%odsun_singlelayer(lay,ist,i)
                 endif
                 
                 transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) - &
                      temp(1,lay) * temp(3,lay) * temp(2,lay)**2_jpim * tausun_levm1 * tau_levm1

                 fac_2_ad(2) = &!fac2_2_ad + &
                             zdown_ad(lay,ist,i) * auxrad_stream%fac6_2(lay,ist,i) * &
                             temp(2,lay) * temp(3,lay) * tausun_levm1
                 
                 fac_2_ad(6) = &!fac6_2_ad + &
                             zdown_ad(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                             temp(2,lay) * temp(3,lay) * tausun_levm1
              else
                 zdown_ad(lay,ist,i) = 0._JPRB
              endif  ! min_tau

       !----------------Upward single scattering of the solar beam-------------------------------
              fac_2_ad(1) = &!fac1_2_ad + &
                          zup_ad(lay,ist,i) * &
                          auxrad_stream%fac2_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i) * &
                         (tausun_levm1 - tausun_lev)

              fac_2_ad(2) = fac_2_ad(2) + &
                          zup_ad(lay,ist,i) * &
                          auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i) * &
                         (tausun_levm1 - tausun_lev)

              fac_2_ad(3) = &!fac3_2_ad + &
                          zup_ad(lay,ist,i) * &
                          auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * &
                         (tausun_levm1 - tausun_lev)

              if(lay > 1) then
                 transmission_aux_ad%tausun_level(levm1,ist,i) = transmission_aux_ad%tausun_level(levm1,ist,i) + &
                      zup_ad(lay,ist,i) * &
                      auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i)
              endif

              transmission_aux_ad%tausun_level(lev,ist,i) = transmission_aux_ad%tausun_level(lev,ist,i) - &
                  zup_ad(lay,ist,i) * &
                  auxrad_stream%fac1_2(lay,ist,i) * auxrad_stream%fac2_2(lay,ist,i) * auxrad_stream%fac3_2(lay,ist,i)

              transmission_aux_ad%od_singlelayer(lay,ist,i) = transmission_aux_ad%od_singlelayer(lay,ist,i) - &
                  fac_2_ad(5) * auxrad_stream%fac5_2(lay,ist,i) !exp(od_singlelayer)

              transmission_aux_ad%odsun_singlelayer(lay,ist,i) = transmission_aux_ad%odsun_singlelayer(lay,ist,i) - &
                   fac_2_ad(4) * auxrad_stream%fac4_2(lay,ist,i) !exp(od_sunsinglelayer)

             raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
                  fac_2_ad(3) * raytracing%pathsun(lay,prof) / &
                  (raytracing%pathsat(lay,prof)+raytracing%pathsun(lay,prof))**2

             raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) - &
                  fac_2_ad(3) * raytracing%pathsat(lay,prof) / &
                  (raytracing%pathsat(lay,prof)+raytracing%pathsun(lay,prof))**2

             transmission_scatt_ir_stream_ad%ssa(ist,i,lay) = transmission_scatt_ir_stream_ad%ssa(ist,i,lay) + &
                  fac_2_ad(2)

             transmission_scatt_ir_stream_ad%azphacup(ist,i,lay) = transmission_scatt_ir_stream_ad%azphacup(ist,i,lay) + &
                  fac_2_ad(1) * coef%ss_solar_spectrum(chan) * z4pi_r

              transmission_scatt_ir_stream_ad%azphacdo(ist,i,lay) = transmission_scatt_ir_stream_ad%azphacdo(ist,i,lay) + &
                  fac_2_ad(6) * coef%ss_solar_spectrum(chan) * z4pi_r
          enddo
       enddo
    endif
 enddo
end subroutine solar_scattering_air_ad

subroutine solar_scattering_near_surf_ad(transmission_aux, transmission_aux_ad, &
                                         ircld, raytracing, raytracing_ad, coef, transmission_scatt_ir_stream_ad, &
                                         sun, sateqsun, iv3lay, iv3lev, pol_id, chanprof, narray, pfraction, &
                                         auxrad_stream, meanrad_up_ad, meanrad_down_ad)

  use rttov_const, only : min_tau, z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, &
                          rttov_chanprof, rttov_coef, transmission_scatt_ir_type, raytracing_type
  use parkind1, only : jpim, jprb, jplm 

  Implicit None
  
  type(transmission_type_aux), intent(in)         :: transmission_aux
  type(ircld_type),     intent(in)                :: ircld
  type(raytracing_type), intent(in)               :: raytracing
  type(rttov_coef),     intent(in)                :: coef
  logical(jplm), intent(in)                       :: sun(:)
  logical(jplm), intent(in)                       :: sateqsun(:,:)
  integer(jpim), intent(in)                       :: iv3lay(:), iv3lev(:)
  integer(jpim), intent(in)                       :: pol_id(:)       ! polarisation index  
  type(rttov_chanprof), intent(in)                :: chanprof(:)
  integer(jpim), intent(in)                       :: narray(:)
  real(jprb), intent(in)                          :: pfraction(:)
  type(radiance_aux), intent(in)                  :: auxrad_stream

  type(transmission_type_aux), intent(inout)      :: transmission_aux_ad
  type(raytracing_type), intent(inout)            :: raytracing_ad
  type(transmission_scatt_ir_type), intent(inout) :: transmission_scatt_ir_stream_ad
  real(jprb),intent(inout)                        :: meanrad_up_ad(0:,:), meanrad_down_ad(0:,:)

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist, lay, lev, lay1
  
  Real(jprb) :: fac_3_ad(7)
  Real(jprb) :: temp(4,0:narray(2))

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)  
  
  Do i = 1, nchannels
     if(sun(i))then
        lay = iv3lay(i)
        lev = iv3lev(i)

        if(pfraction(i) < 0.0_JPRB )then
           lay1 = iv3lay(i)
        else
           lay1 = iv3lay(i) + 1
        endif

       do ist = 0, ircld%nstream(prof) !rev loop not nec.
          if(auxrad_stream%meanrad_down(ist,i) < 0._jprb) then
             meanrad_down_ad(ist,i) = 0._jprb
          endif
          
          fac_3_ad(:) = 0.0_jprb
          temp(1,ist) = auxrad_stream%fac6_3(ist,i) * auxrad_stream%fac2_3(ist,i) * meanrad_down_ad(ist,i)
          temp(2,ist) = transmission_aux%tau_level_r(lev,ist,i) * transmission_aux%tau_surf_r(ist,i)
          temp(3,ist) = temp(2,ist)**2_jpim

          if(.not. sateqsun(lay,prof)) then                    
             temp(4,ist) = (auxrad_stream%fac5_3(ist,i) - auxrad_stream%fac4_3(ist,i)) * auxrad_stream%fac7_3(ist,i)

             fac_3_ad(4) = &!fac4_3_ad &
                         - temp(1,ist) * temp(2,ist) * auxrad_stream%fac7_3(ist,i) * &
                         transmission_aux%tausun_level(lev,ist,i)
          
             fac_3_ad(5) = &!fac5_3_ad + &
                         + temp(1,ist) * temp(2,ist) * auxrad_stream%fac7_3(ist,i) * &
                         transmission_aux%tausun_level(lev,ist,i)
          
             fac_3_ad(7) = &!fac7_3_ad + &
                         temp(1,ist) * temp(2,ist) * &
                        (auxrad_stream%fac5_3(ist,i) - auxrad_stream%fac4_3(ist,i)) * &
                        transmission_aux%tausun_level(lev,ist,i)

             raytracing_ad%pathsat(lay1,prof) = raytracing_ad%pathsat(lay1,prof) + &
                                                fac_3_ad(7) * raytracing%pathsun(lay1,prof) / &
                                               (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2

             raytracing_ad%pathsun(lay1,prof) = raytracing_ad%pathsun(lay1,prof) - &
                                                fac_3_ad(7) * raytracing%pathsat(lay1,prof) / &
                                               (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2
          else
             temp(4,ist) = auxrad_stream%fac4_3(ist,i) * transmission_aux%odsun_sfrac(ist,i)
             
             transmission_aux_ad%odsun_sfrac(ist,i) = transmission_aux_ad%odsun_sfrac(ist,i) + &
                                                      temp(1,ist) * temp(2,ist) * &
                                                      auxrad_stream%fac4_3(ist,i) * &
                                                      transmission_aux%tausun_level(lev,ist,i)

             fac_3_ad(4) = &!fac4_3_ad + &
                         temp(1,ist) * temp(2,ist) * &
                         transmission_aux%tausun_level(lev,ist,i) * transmission_aux%odsun_sfrac(ist,i)
          endif

          if(transmission_aux%tau_level(lev,ist,i) > min_tau) then

             transmission_aux_ad%tau_surf(ist,i) = transmission_aux_ad%tau_surf(ist,i) - &
                                                   temp(1,ist) * temp(3,ist) * temp(4,ist) * &
                                                   transmission_aux%tausun_level(lev,ist,i) * &
                                                   transmission_aux%tau_level(lev,ist,i)

             transmission_aux_ad%tau_level(lev,ist,i) = transmission_aux_ad%tau_level(lev,ist,i) - &
                                                        temp(1,ist) * temp(3,ist) * temp(4,ist) * &
                                                        transmission_aux%tausun_level(lev,ist,i) * &
                                                        transmission_aux%tau_surf(ist,i)
 
             transmission_aux_ad%tausun_level(lev,ist,i) = transmission_aux_ad%tausun_level(lev,ist,i) + &
                                                           temp(1,ist) * temp(2,ist) * temp(4,ist)

             fac_3_ad(2) = &!fac2_3_ad + &
                         meanrad_down_ad(ist,i) * auxrad_stream%fac6_3(ist,i) * &
                         temp(2,ist) * temp(4,ist) * transmission_aux%tausun_level(lev,ist,i)

             fac_3_ad(6) = &!fac6_3_ad + &
                         meanrad_down_ad(ist,i) * auxrad_stream%fac2_3(ist,i) * &
                         temp(2,ist) * temp(4,ist) * transmission_aux%tausun_level(lev,ist,i)
          endif

       !--------------Upward single scattering of the solar beam---------------------------------

          if (pol_id(i) >= 6_jpim) meanrad_up_ad(ist,i) = 0.0_jprb

          transmission_aux_ad%tausun_level(lev,ist,i) = transmission_aux_ad%tausun_level(lev,ist,i) + &
                                                        meanrad_up_ad(ist,i) * &
                          auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i)

          transmission_aux_ad%tausun_surf(ist,i) = transmission_aux_ad%tausun_surf(ist,i) - &
                                                   meanrad_up_ad(ist,i) * &
                     auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i)

          fac_3_ad(3) = &!fac3_3_ad + &
                      meanrad_up_ad(ist,i) * &
                      auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac2_3(ist,i) * &
                     (transmission_aux%tausun_level(lev,ist,i) - transmission_aux%tausun_surf(ist,i))

          fac_3_ad(2) = fac_3_ad(2) + &
                      meanrad_up_ad(ist,i) * &
                      auxrad_stream%fac1_3(ist,i) * auxrad_stream%fac3_3(ist,i) * &
                     (transmission_aux%tausun_level(lev,ist,i) - transmission_aux%tausun_surf(ist,i))

          fac_3_ad(1) = &!fac1_3_ad + &
                      meanrad_up_ad(ist,i) * &
                      auxrad_stream%fac2_3(ist,i) * auxrad_stream%fac3_3(ist,i) * &
                     (transmission_aux%tausun_level(lev,ist,i) - transmission_aux%tausun_surf(ist,i))

          transmission_scatt_ir_stream_ad%azphacdo(ist,i,lay) = transmission_scatt_ir_stream_ad%azphacdo(ist,i,lay) + &
               fac_3_ad(6) * coef%ss_solar_spectrum(chan) * z4pi_r

          transmission_aux_ad%od_sfrac(ist,i) = transmission_aux_ad%od_sfrac(ist,i) - &
               fac_3_ad(5) * auxrad_stream%fac5_3(ist,i)

          transmission_aux_ad%odsun_sfrac(ist,i) = transmission_aux_ad%odsun_sfrac(ist,i) - &
               fac_3_ad(4) * auxrad_stream%fac4_3(ist,i)

          raytracing_ad%pathsat(lay,prof) = raytracing_ad%pathsat(lay,prof) + &
               fac_3_ad(3) * &
               raytracing%pathsun(lay,prof) / (raytracing%pathsat(lay,prof) + raytracing%pathsun(lay,prof))**2

          raytracing_ad%pathsun(lay,prof) = raytracing_ad%pathsun(lay,prof) - &
               fac_3_ad(3) * &
               raytracing%pathsat(lay,prof) / (raytracing%pathsat(lay,prof) + raytracing%pathsun(lay,prof))**2

          transmission_scatt_ir_stream_ad%ssa(ist,i,lay) = transmission_scatt_ir_stream_ad%ssa(ist,i,lay) + &
               fac_3_ad(2)

          transmission_scatt_ir_stream_ad%azphacup(ist,i,lay) = transmission_scatt_ir_stream_ad%azphacup(ist,i,lay) + &
               fac_3_ad(1) * coef%ss_solar_spectrum(chan) * z4pi_r

       enddo
    endif
 enddo
end subroutine solar_scattering_near_surf_ad

! Currently nothing is being done for sea ice.       
subroutine solar_contribution_ad(transmission_aux, transmission_aux_ad, &
                                 ircld, sunglint, sunglint_ad, fresnrefl, fresnrefl_ad, reflectivity, reflectivity_ad, &
                                 coef, profiles, chanprof, sun, narray, auxrad_stream_ad, rad_ad)

  use rttov_const, only : surftype_land, surftype_sea, surftype_seaice
  use rttov_types, Only : rttov_chanprof, ircld_type, rttov_coef, profile_type, &
                          transmission_type_aux, sunglint_type, radiance_type, radiance_aux
  use parkind1, only : jpim, jprb, jplm 

  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(ircld_type),     intent(in)              :: ircld
  type(sunglint_type), intent(in)               :: sunglint
  real(jprb), intent(in)                        :: fresnrefl(:)
  real(jprb), intent(in)                        :: reflectivity(:)
  type(rttov_coef), intent(in)                  :: coef
  type(profile_type), intent(in)                :: profiles(:)
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  logical(jplm), intent(in)                     :: sun(:)
  integer(jpim), intent(in)                     :: narray(:)

  type(radiance_type), intent(in)               :: rad_ad
  type(radiance_aux), intent(in)                :: auxrad_stream_ad

  type(transmission_type_aux), intent(inout)    :: transmission_aux_ad
  type(sunglint_type), intent(inout)            :: sunglint_ad
  real(jprb), intent(inout)                     :: fresnrefl_ad(:)
  real(jprb), intent(inout)                     :: reflectivity_ad(:)
  
  integer(jpim) :: nlevels, nlayers, nstreams, nchannels, i, ist
  real(jprb)    :: temp

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)  

  Do i = 1, nchannels
     if(sun(i)) then
        nstreams = ircld%nstream(prof)
        if(profiles(prof)%skin%surftype == surftype_sea) then
           temp = fresnrefl(i) * coef%ss_solar_spectrum(chan) * &
                 (transmission_aux%tausun_surf(0,i) * rad_ad%clear(i) + &
                  sum(transmission_aux%tausun_surf(1:nstreams,i) * auxrad_stream_ad%cloudy(1:nstreams,i)))

           if(sunglint%s(prof)%windsp .eq. 0._jprb) then
              ! note: must increment because of degeneracy of channels onto same prof
              sunglint_ad%s(prof)%glint = sunglint_ad%s(prof)%glint + temp * transmission_aux%refl_norm(i) 
           else
              sunglint_ad%s(prof)%glint = sunglint_ad%s(prof)%glint + temp
           endif

           ist = 0_jpim

           transmission_aux_ad%tausun_surf(ist,i) = transmission_aux_ad%tausun_surf(ist,i) + &
                                                    sunglint%s(prof)%glint * fresnrefl(i) * &
                                                    rad_ad%clear(i) * coef%ss_solar_spectrum(chan)
           fresnrefl_ad(i) = fresnrefl_ad(i) + & !first use but zeroed outside rttov_integrate so left alone.
                             sunglint%s(prof)%glint * transmission_aux%tausun_surf(ist,i) * &
                             rad_ad%clear(i) * coef%ss_solar_spectrum(chan)
         
           do ist=1, ircld%nstream(prof)
              transmission_aux_ad%tausun_surf(ist,i) = transmission_aux_ad%tausun_surf(ist,i) + &
                                                      sunglint%s(prof)%glint * fresnrefl(i) * &
                                                      auxrad_stream_ad%cloudy(ist,i) * coef%ss_solar_spectrum(chan)
              fresnrefl_ad(i) = fresnrefl_ad(i) + &
                                sunglint%s(prof)%glint * transmission_aux%tausun_surf(ist,i) * &
                                auxrad_stream_ad%cloudy(ist,i) * coef%ss_solar_spectrum(chan)
           enddo

           if(sunglint%s(prof)%windsp .eq. 0._jprb) then
              transmission_aux_ad%tausun_surf(:,i) = transmission_aux_ad%tausun_surf(:,i) * transmission_aux%refl_norm(i)
              fresnrefl_ad(i) = fresnrefl_ad(i) * transmission_aux%refl_norm(i)
           endif

        else ! sea-ice = land
           ist = 0_jpim
           reflectivity_ad(i) = reflectivity_ad(i) + &
                                rad_ad%clear(i) * transmission_aux%tausun_surf(ist,i) * &
                                coef%ss_solar_spectrum(chan) * transmission_aux%refl_norm(i)
           transmission_aux_ad%tausun_surf(ist,i) = transmission_aux_ad%tausun_surf(ist,i) + &
                                                    rad_ad%clear(i) * reflectivity(i) * &
                                                    coef%ss_solar_spectrum(chan) * transmission_aux%refl_norm(i)
           do ist=1, ircld%nstream(prof)
              reflectivity_ad(i) = reflectivity_ad(i) + &
                                   auxrad_stream_ad%cloudy(ist,i) * transmission_aux%tausun_surf(ist,i) * &
                                   coef%ss_solar_spectrum(chan) * transmission_aux%refl_norm(i)

              transmission_aux_ad%tausun_surf(ist,i) = transmission_aux_ad%tausun_surf(ist,i) + &
                                                       auxrad_stream_ad%cloudy(ist,i) * reflectivity(i) * &
                                                       coef%ss_solar_spectrum(chan) * transmission_aux%refl_norm(i)
           enddo
        endif
     endif
  enddo
end subroutine solar_contribution_ad

End Subroutine rttov_integrate_ad
