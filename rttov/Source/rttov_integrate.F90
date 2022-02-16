subroutine rttov_integrate(addcosmic, opts, maxnstreams , chanprof, &! in
                           emissivity, reflectivity, fresnrefl, sunglint, &! in
                           sun, transmission_aux, transmission_scatt_ir_stream, &!in
                           profiles, aux_prof, coef, raytracing, ircld, &! in
                           rad, auxrad, auxrad_stream)      ! inout

! Description:
! To perform integration of radiative transfer equation
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
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF).
!
! Matricardi, M. 2005 The inclusion of aerosols and clouds in RTIASI,the
! ECMWF radiative transfer model for the infrared atmospheric sounding
! interferometer.
! ECMWF Research Dept. Tech. Memo. 474 (available from the librarian at ECMWF).
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
!  1.4     03/09/2004   Mods. for Vectorisation (D Salmond ECMWF & BCarruthers, Cray)
!  1.5     28/02/2005   Further mods to vectorisation (D Dent)
!  1.6     01/06/2005   Marco Matricardi (ECMWF):
!             --        IASI capability added.
!             --        Linear in tau approximation for RT equation introduced.
!             --        Solar radiation introduced for IASI and AIRS.
!  1.7     01/06/2006   Marco Matricardi (ECMWF):
!                       Multiple scattering in the infrared introduced for
!                       water clouds,cirrus clouds and aerosols.
!  1.8     26/01/2007   Removed polarisation (R Saunders)
!  1.9     08/03/2007   Reintroduced overcast cloud (R Saunders)
!  1.10    13/07/2007   Added extra variables requested by P. Watts(RSaunders)
!  1.11    15/11/2007   Changed skin to surfair for overcast radiance calc (R Saunders)
!  1.12    27/11/2007   Optimised for NEC/IBM (D. Salmond)
!  1.13    21/12/2007   Added polarimetric option (R. Saunders)
!  1.14    27/02/2009   Profile levels to include ToA. Distinguish between
!                       layer arrays and level arrays - size, index labels,
!                       looping (P. Rayer)
!  1.15    17/04/2009   Model-top brought into calculation of source function (P.Rayer)
!  1.16    29/06/2009   Top-layer brought into layer looping to shorten code (P.Rayer)
!  1.17    03/11/2009   Transmittances on levels (A Geer)
!  1.18    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
!                       related quantities are now layer arrays (Marco Matricardi).
!  1.19    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  1.20    14/10/2010   Remove rt8_mode (J Hocking)
!  1.21    14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                       should be performed (J Hocking)
!  2.0     14/12/2011   Re-written (D Rundle)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".

  use parkind1, Only : jpim, jprb, jplm

  use rttov_types, Only : rttov_chanprof, rttov_coef, rttov_options, profile_type, profile_aux, transmission_type_aux, &
                          transmission_scatt_ir_type, sunglint_type, radiance_type, ircld_type, raytracing_type, radiance_aux

!INTF_OFF
  use rttov_const, Only : sensor_id_po, min_od, min_tau, z4pi_r, surftype_land, surftype_sea, surftype_seaice, deg2rad, pi_r

  use yomhook, Only : LHOOK, DR_HOOK
!INTF_ON

  Implicit None

!subroutine arguments:

  logical(jplm),   intent(in)                   :: addcosmic   ! switch for adding cosmic background
  type(rttov_options),  intent(in)              :: opts        ! options structure
  integer(jpim),   intent(in)                   :: maxnstreams                !
  type(rttov_chanprof), intent(in)              :: chanprof(:)     ! Channel indices
  type(profile_type),   intent(in)              :: profiles(:) ! Profiles
  real(jprb),      intent(in)                   :: emissivity(size(chanprof))   ! surface emissivity
  real(jprb),      intent(in)                   :: reflectivity(size(chanprof)) ! surface reflectivity
  real(jprb),      intent(in)                   :: fresnrefl(size(chanprof))
  logical(jplm),   intent(in)                   :: sun(size(chanprof))
  type(ircld_type),     intent(in)              :: ircld
  type(raytracing_type),intent(in)              :: raytracing
  type(sunglint_type),  intent(in)              :: sunglint
  type(transmission_type_aux), intent(inout)    :: transmission_aux            ! transmittances and single-layer od
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir_stream
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef
  type(radiance_aux),   intent(inout)           :: auxrad_stream
  type(radiance_type),  intent(inout)           :: rad    ! radiances (mw/cm-1/ster/sq.m) and BTs
  type(radiance_aux),   intent(inout)           :: auxrad ! auxillary radiances
!INTF_END

!module variables: 

  integer(jpim)              :: i, lev, ist, lay, nchannels, nlayers, nlevels, prof, iprof, chan, totnstreams ! counter variables

  integer(jpim), allocatable :: nstreams(:)
  integer(jpim), allocatable :: iv2lev(:), iv2lay(:), iv3lev(:), iv3lay(:)
  integer(jpim), allocatable :: pol_id(:)       ! polarisation index  

  real(jprb), allocatable    :: cfraction(:), pfraction(:)        ! cloud fraction
  logical(jplm), allocatable :: sateqsun(:,:) ! True where the solar zenith angle equal to observation angle

  logical(jplm)              :: anysun ! Are *any* channels affected by sun (only need to make check once if negative

! Define macros for commonly used variables
#define prof chanprof(i)%Prof
#define chan chanprof(i)%Chan

#define tau_surf_r transmission_aux%Tau_surf_r(ist,i)
#define tau_surf transmission_aux%Tau_surf(ist,i)
#define tau_layer_r transmission_aux%Tau_level_r(lay,ist,i)
#define tau_layer transmission_aux%Tau_level(lay,ist,i)
#define tau_level_r transmission_aux%Tau_level_r(lay+1,ist,i)
#define tau_level transmission_aux%Tau_level(lay+1,ist,i)
#define daux_lay (auxrad%air(lay+1,i) - auxrad%air(lay,i))
#define dtau_lay (tau_layer - tau_level)
#define od_singlelayer_r transmission_aux%Od_singlelayer_r(lay,ist,i)
#define od_singlelayer transmission_aux%Od_singlelayer(lay,ist,i)
#define fac1 transmission_aux%fac(1,lay,ist,i)
#define fac2 transmission_aux%fac(2,lay+1,ist,i)
#define surf_fac transmission_aux%Surf_fac(ist,i)

#include "rttov_calcbt.h"
#include "rttov_calcrad.h"

  REAL(JPRB) :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE',0_jpim,ZHOOK_HANDLE)

!---------------------------
!0. Initialise useful variables
!---------------------------

call init_rttov_integrate_mod_vars(chanprof, profiles, ircld, aux_prof, coef, raytracing)

totnstreams = maxnstreams
anysun = any(sun)
rad%down = 0.0_jprb

!----------------------------
!1. calculate layer radiances
!----------------------------

call rttov_calcrad(addcosmic, chanprof, profiles,  coef, &!in
                   auxrad%cosmic, auxrad%skin, auxrad%surfair, auxrad%air) !out

!-------------------------------------
!2. calculate atmospheric contribution
!-------------------------------------
 
call calc_atmospheric_radiance(transmission_aux, auxrad, &! in
                               auxrad_stream) ! out

!---Scattering of the solar beam--------------------------------------------------
if((opts%addaerosl .OR. opts%addclouds) .and. anysun) &
     call  solar_scattering_air(transmission_aux, sun, raytracing, coef, chanprof, transmission_scatt_ir_stream, & !in
                                    auxrad_stream) !out

#ifndef RTTOV_KEYRADONLY
ist = 0_jpim
! Phil Watts outputs
Do i = 1, nchannels
  Do lay = 2, nlayers
     If (tau_level > min_tau ) Then
        rad%down(lay,i) = auxrad_stream%down(lay,0,i) * tau_level
     Else
        rad%down(lay,i) = rad%down(lay-1,i)
     Endif
     rad%surf(lay,i) = auxrad%air(lay+1,i)
  End Do
enddo
#endif

!-------------------------------------------------------------------------------
!2a calculate near-surface layer contribution
!-------------------------------------------------------------------------------

call calc_near_surf_contribution(transmission_aux, auxrad, &!in
                                 auxrad_stream) !inout

!---Scattering of the solar beam---------------------------
if((opts%addaerosl .OR. opts%addclouds) .and. anysun) &
     call solar_scattering_near_surf(transmission_aux, sun, chanprof, & !in
                                         auxrad_stream) ! out

do i = 1, nchannels
   ist = 0_jpim  

#ifndef RTTOV_KEYRADONLY
! clear sky radiance without reflection term
   rad%upclear(i) = auxrad_stream%meanrad_up(ist,i)

!  Phil Watts outputs
  If (tau_surf > min_tau ) Then
    rad%down(iv2lay(i),i) = auxrad_stream%meanrad_down(ist,i) * tau_surf
  Else
    rad%down(iv2lay(i),i) = rad%down(iv3lay(i),i)
  Endif
  rad%surf(iv2lay(i),i) = auxrad%skin(i)

  ! clear sky downwelling radiance
  rad%dnclear(i) = auxrad_stream%meanrad_down(ist,i) * tau_surf
  ! reflected clear sky downwelling radiance
  rad%reflclear(i) = rad%dnclear(i) * reflectivity(i) * tau_surf
#endif

  rad%clear(i) = auxrad_stream%meanrad_up(ist,i) + auxrad_stream%meanrad_down(ist,i) * reflectivity(i) * tau_surf**2_jpim

  auxrad_stream%up(iv2lay(i),ist,i) = auxrad_stream%meanrad_up(ist,i)

   do ist = 1, nstreams(i)
      auxrad_stream%up(iv2lay(i),ist,i) = auxrad_stream%meanrad_up(ist,i)
      auxrad_stream%cloudy(ist,i) = auxrad_stream%meanrad_up(ist,i) + &
                                    auxrad_stream%meanrad_down(ist,i) * reflectivity(i) * tau_surf**2_jpim
   enddo
enddo

!-----------------------
!3. calculate surface contribution
!-----------------------

!cdir nodep
Do i = 1, nchannels
   ist = 0_jpim
! clear sky radiance without reflection term

rad%clear(i) = rad%clear(i) + auxrad%skin(i) * emissivity(i) * tau_surf

#ifndef RTTOV_KEYRADONLY
rad%upclear(i) = rad%upclear(i) + auxrad%skin(i) * emissivity(i) * tau_surf
#endif

!cdir nodep
   do ist = 1, nstreams(i)
      auxrad_stream%cloudy(ist,i) = auxrad_stream%cloudy(ist,i) + auxrad%skin(i) * emissivity(i) * tau_surf
   End Do
End Do

!--------------------------------
!4. Add solar contribution
!--------------------------------

if(coef%fmv_model_ver == 9) then
   if(anysun) call solar_contribution(transmission_aux, sun, sunglint, fresnrefl, reflectivity, & !in
                                          coef, profiles, chanprof, & !in
                                          auxrad_stream, rad) ! out
endif

!--------------------------------
!5. cosmic temperature correction
!--------------------------------
if (addcosmic) then
   ist = 0
   rad%clear(1:nchannels) = rad%clear(1:nchannels) + reflectivity(1:nchannels) * &
                            transmission_aux%Tau_surf(ist,1:nchannels)**2_jpim * auxrad%cosmic(1:nchannels)
   do i = 1,nchannels 
      do ist = 1, nstreams(i)
         auxrad_stream%cloudy(ist,i) = auxrad_stream%cloudy(ist,i) + &
                                       auxrad%cosmic(i) * reflectivity(i) * tau_surf**2_jpim
      enddo
   enddo
Endif

!---------------------------------------------------
!6. calculate overcast radiances
!---------------------------------------------------
!---------------
!6.1 Upward part
!---------------

rad%up(:,1:nchannels)= auxrad_stream%up(:,0,1:nchannels)  ! Phil Watts output

ist = 0
Do i = 1, nchannels
   do lay = 1,nlayers
      lev = lay + 1
! overcast radiances at given cloud top
      rad%overcast(lay,i) = auxrad_stream%up(lay,ist,i) + auxrad%air(lev,i) * tau_level
   end do
enddo

! Add surface component to overcast radiances
do i = 1, nchannels
   lay = iv2lay(i)
   rad%overcast(lay,i) = auxrad_stream%up(lay,ist,i) + tau_surf * auxrad%surfair(i)
enddo

! DAR: can this be just addclouds? Also, is it necessary to calculate overcast radiances when using RTTOV-CLOUDY? 
if(opts%addaerosl .OR. opts%addclouds) then 
!----------------------------------------
!7. calculate complex cloudy radiances
!----------------------------------------
   rad%cloudy = 0._jprb ! nchannels
   Do i = 1, nchannels
      do ist = 1, nstreams(i)
         rad%cloudy(i) = rad%cloudy(i) + &
                         auxrad_stream%cloudy(ist,i) * &
                         (ircld%xstr(ist+1,prof) - ircld%xstr(ist,prof))
      enddo
      rad%cloudy(i) = rad%cloudy(i) + rad%clear(i) * ircld%XSTRCLR(prof)
   End do
   
!---------------------------
!8. calculate total radiance (cloudy case)
!---------------------------
   rad%total(1:nchannels) = rad%cloudy(1:nchannels)
else
!6. Calculate "simple" cloudy radiances
!--------------------------------------------
!6.2 Interpolate to given cloud-top pressures
!--------------------------------------------
do i = 1, nchannels
   
   lay = aux_prof%s(prof)%nearestlev_ctp - 1
   rad%cloudy(i) = rad%overcast(lay,i) * &
                  (1.0_JPRB - aux_prof%s(prof)%pfraction_ctp) + &
                   rad%overcast(lay-1,i) * (aux_prof%s(prof)%pfraction_ctp)
End Do
!---------------------------
!8. calculate total radiance (clear case)
!---------------------------
   if(opts%addpc) then
      rad%total(1:nchannels) = rad%clear(1:nchannels)
   else     
      rad%total(1:nchannels) = rad%clear(1:nchannels) + &
                               cfraction(1:nchannels) * &
                              (rad%cloudy(1:nchannels) - rad%clear(1:nchannels) )
   endif
Endif

!-----------------------------------------------
!9. convert radiances to brightness temperatures
!-----------------------------------------------
Call rttov_calcbt(chanprof, coef, rad)

! deallocate module arrays - DAR: can move these to a more useful array if they are commonly used?
deallocate(nstreams, iv2lev, iv2lay, iv3lev, iv3lay, pol_id, cfraction, pfraction, sateqsun) 

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE',1_jpim,ZHOOK_HANDLE)

contains

! DAR: This subroutine calculates the individual layer clear-sky radiances and then does a cumulative sum to determine the 
!      radiance observed from the top of atmosphere to a particular layer.
! DAR: As expected, this routine consumes most of the time (loops over channels, streams and levels)
!      and has been most heavily optimised (for IBM only so far, Intel shows neutral impact - will look into this)
!      I have taken out the code that switches the order of the loops on the NEC and will let MF test this impact       
!      See subroutine comments for more details of individual changes.
subroutine calc_atmospheric_radiance(transmission_aux, auxrad, auxrad_stream)

  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(radiance_aux), intent(in)                :: auxrad
  type(radiance_aux), intent(inout)             :: auxrad_stream


! DAR: fac and fac2 contain real 1 or 0 depending on whether a calculation should be performed or not.
!      IBM performance was suffering as a result of doing lots of branching (mispredicts?) so these arrays are populated in advance
!      and the calculation is performed regardless.

! DAR: I've removed a lot of the 'temporary' variables that were eventually summed together because they were taking up a signficant
!      amount of time on the IBM. Now the big sum is done with all the variables in full because the 16 prefetch streams can handle 
!      this (POWER6) - maybe will have to split this up for POWER7 (only 12). Will do more testing on Intel with vtune to check 
!      performance impact.

! DAR: I've also removed od_singlelayer_r and added it to transmission_aux so it can be reused the TL/AD/K code when ready 
!      and doesn't have to be recalculated - This should be faster...
  
! DAR: Slow code path - may be faster on intel?
  if(transmission_aux%anynegtau .gt. 0.0_jprb) then
     do i = 1, nchannels
      do ist = 0, nstreams(i)
         do lay = 1, nlayers
! DAR: As explained above, every difference is calculated
!      explicitly as this is apparently faster on the IBM. Also, changed indices on auxrad_stream%up/down to be consistent with
!      every other variable.
         auxrad_stream%up(lay,ist,i) = fac1 * &
              (dtau_lay * &
              (auxrad%air(lay,i) + daux_lay * od_singlelayer_r) &
              - daux_lay * tau_level)         

         auxrad_stream%down(lay,ist,i) = fac1 * fac2 * &
                                      ((dtau_lay * (auxrad%air(1,i) - daux_lay * od_singlelayer_r) * &
                                        (tau_level_r * tau_layer_r)) + &
                                        daux_lay * tau_level_r)
         enddo
      enddo
   enddo

! DAR: This is slow. Hopefully it'll never run and can be deleted
   do  i = 1, nchannels
      do ist = 0, nstreams(i)
         do lay = 1,nlayers
            if(tau_level < 0.0_jprb .and. tau_layer >= 0._jprb) then
               auxrad_stream%up(lay,ist,i) = (0.5_JPRB * (auxrad%air(lay+1,i) + auxrad%air(lay,i))) * dtau_lay
               auxrad_stream%down(lay,ist, i) = 0.0_jprb
            endif
         enddo
      enddo
   enddo

   do i = 1,nchannels
      do ist = 0, nstreams(i)
         do lay = 2, nlayers
            auxrad_stream%up(lay,ist,i) = auxrad_stream%up(lay, ist,i) + auxrad_stream%up(lay-1,ist,i)
            auxrad_stream%down(lay,ist,i) = auxrad_stream%down(lay,ist,i) + auxrad_stream%down(lay-1,ist,i)
         enddo
      enddo
   enddo
   
else ! DAR: fast code (on IBM) when no special cases. Note cumulative sum is done at same time as it's significantly quicker than 
   !      doing it later. This code is still very hard to read though.
   do i = 1, nchannels
      do ist = 0, nstreams(i)
         lay = 1
         auxrad_stream%up(lay,ist,i) = fac1 * &
                                      (dtau_lay * &
                                      (auxrad%air(lay,i) + daux_lay * od_singlelayer_r) - &
                                       daux_lay * tau_level)

         auxrad_stream%down(lay,ist,i) = fac1 * fac2 * &
                                       ((dtau_lay * (auxrad%air(lay,i) - daux_lay * od_singlelayer_r) * &
                                        (tau_level_r * tau_layer_r)) + &
                                         daux_lay * tau_level_r)

         do lay = 2, nlayers
            auxrad_stream%up(lay,ist,i) = auxrad_stream%up(lay-1, ist, i) + fac1 * &
                 (dtau_lay * &
                 (auxrad%air(lay,i) + daux_lay * od_singlelayer_r) - &
                 daux_lay * tau_level)         

            auxrad_stream%down(lay,ist,i) = auxrad_stream%down(lay-1, ist, i) + &
                 fac1 * fac2 * &
                 ((dtau_lay * (auxrad%air(lay,i) - daux_lay * od_singlelayer_r) * &
                 (tau_level_r * tau_layer_r)) + &
                 daux_lay * tau_level_r)
         enddo
      enddo
   enddo
endif

end subroutine calc_atmospheric_radiance

! DAR: This is the next biggest consumer of CPU time and should be looked at next. The trouble seems to be that you have to change a
!      lot of non-consecutive data. Maybe the RTTOV 9 comments can be removed
subroutine calc_near_surf_contribution(transmission_aux, auxrad, & !in
                                       auxrad_stream) !inout
  Implicit None

  type(transmission_type_aux), intent(in) :: transmission_aux
  type(radiance_aux), intent(in)          :: auxrad
  type(radiance_aux),   intent(inout)     :: auxrad_stream

  real(jprb) :: rad_tmp(0:totnstreams,nchannels) 

#define B1_3  auxrad%air(lev,i) * (tau_level  - tau_surf)
#define B2_3 (auxrad%surfair(i) - auxrad%air(lev,i)) * tau_surf
#define B3_3 (auxrad%surfair(i) - auxrad%air(lev,i)) * (tau_level - tau_surf) * (transmission_aux%od_sfrac_r(ist,i))             
 
  Do i = 1, nchannels
     lay = iv3lay(i)
     lev = iv3lev(i)
     do ist = 0, nstreams(i)

        if(tau_surf < 0._JPRB) then 
           rad_tmp(ist,i) = 0.0_JPRB
           
           if(tau_level >= 0._JPRB) then
              auxrad_stream%meanrad_up(ist,i) = (0.5_JPRB * (auxrad%surfair(i) + auxrad%air(lev,i))) * (tau_level - tau_surf)
           endif

        else !tau_surf >= 0

           if(transmission_aux%od_sfrac(ist,i) < min_od .OR. ((tau_level - tau_surf) < min_od)) THEN
           rad_tmp(ist,i) = 0.0_JPRB

! small optical depth or optical depth change set radiance to zero
              auxrad_stream%meanrad_up(ist,i) = 0.0_JPRB
           else
              ! compute linear in tau radiances
              auxrad_stream%meanrad_up(ist,i) = B1_3 - &
                                                B2_3 + &
                                                B3_3
              rad_tmp(ist,i) = surf_fac * & 
                              (B1_3 - &
                               B3_3) * &
                               tau_level_r * tau_surf_r + &
                               B2_3 * tau_surf_r**2                           
           endif
        endif
        
        if (pol_id(i) >= 6_jpim ) then
           auxrad_stream%meanrad_up(ist,i) = 0.0_jprb ! auxrad_stream%meanrad_down not zeroed DAR
        else
           auxrad_stream%meanrad_up(ist,i) = auxrad_stream%meanrad_up(ist,i) + auxrad_stream%up(lay,ist,i)
        endif

        auxrad_stream%meanrad_down(ist,i) = auxrad_stream%down(lay,ist,i) + rad_tmp(ist,i)
     enddo
  enddo

end subroutine calc_near_surf_contribution

! DAR: I can't find any explanation for this code and it only affects a few channels so I've mainly left it alone for the time being
subroutine solar_scattering_air(transmission_aux, sun, raytracing, coef, chanprof, &
                         transmission_scatt_ir_stream, auxrad_stream)
  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux
  type(raytracing_type), intent(in)             :: raytracing
  type(rttov_coef),     intent(in)              :: coef
  type(radiance_aux), intent(inout)             :: auxrad_stream
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir_stream

  logical(jplm), intent(in)                     :: sun(:)
  real(jprb)                                    :: temp(nlayers,0:totnstreams)
    
  do i=1,nchannels
     if(sun(i)) then

        do ist = 0, nstreams(i)
           auxrad_stream%fac6_2(:,ist,i) = coef%ss_solar_spectrum(chan) * z4pi_r * &
                                           transmission_scatt_ir_stream%azphacdo(ist,i,:)
           
           auxrad_stream%fac1_2(:,ist,i) = coef%ss_solar_spectrum(chan) * z4pi_r * &
                                           transmission_scatt_ir_stream%azphacup(ist,i,:)
 
           auxrad_stream%fac2_2(:,ist,i) = transmission_scatt_ir_stream%ssa(ist,i,:)
        enddo

        do lay = 1, nlayers
           auxrad_stream%fac3_2(lay,:,i) = raytracing%pathsat(lay,prof) / &
                                          (raytracing%pathsun(lay,prof) + raytracing%pathsat(lay,prof))
        enddo

        do ist = 0, nstreams(i)
          auxrad_stream%fac4_2(:,ist,i) = exp(-transmission_aux%odsun_singlelayer(:,ist,i)) 
          auxrad_stream%fac5_2(:,ist,i) = exp(-transmission_aux%Od_singlelayer(:,ist,i)) ! capital O to avoid fpp
        enddo
        
        do lay = 1, nlayers
           if(.not. sateqsun(lay,prof)) then
              auxrad_stream%fac7_2(lay,:,i) = raytracing%pathsat(lay,prof) / &
                                             (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))
           endif
        enddo

#define fac1_2 auxrad_stream%Fac1_2(lay,ist,i)
#define fac2_2 auxrad_stream%Fac2_2(lay,ist,i)
#define fac3_2 auxrad_stream%Fac3_2(lay,ist,i)
#define fac4_2 auxrad_stream%Fac4_2(lay,ist,i)
#define fac5_2 auxrad_stream%Fac5_2(lay,ist,i)
#define fac6_2 auxrad_stream%Fac6_2(lay,ist,i)
#define fac7_2 auxrad_stream%Fac7_2(lay,ist,i)
#define tausun_layer transmission_aux%tausun_level(lay,ist,i)

#define dfac54_2 (fac5_2 - fac4_2)

  !----------------Upward single scattering of the solar beam-----------------------
  !
  ! auxrad_stream has already been layer integrated so must add integrated (not single layer) flux to each layer
        do ist=0, nstreams(i)
          temp(:,ist) = (auxrad_stream%Fac1_2(:,ist,i) * &
                        auxrad_stream%Fac2_2(:,ist,i) * &
                        auxrad_stream%Fac3_2(:,ist,i) * &
                        (transmission_aux%tausun_level(1:nlevels-1,ist,i) - &
                        transmission_aux%tausun_level(2:nlevels,ist,i)))
        enddo

        do ist=0, nstreams(i)
           do lay=2, nlayers
              temp(lay,ist) = temp(lay,ist) + temp(lay-1,ist) 
           enddo
        enddo

        do ist=0, nstreams(i)
          auxrad_stream%up(:,ist,i) = auxrad_stream%up(:,ist,i) + temp(:,ist)
        enddo

  !-------------------Downward single scattering of the solar beam------------------
        do ist = 0, nstreams(i)
           do lay = 1, nlayers

                 temp(lay,ist) = fac2 * &
                                 fac6_2 * fac2_2 * &
                                 tausun_layer * tau_layer_r * tau_level_r

                 if(.not. sateqsun(lay,prof)) then ! prof is fn of channel
                    temp(lay,ist) = temp(lay,ist) * fac7_2 * dfac54_2
                 else
                    temp(lay,ist) = temp(lay,ist) * fac4_2 * transmission_aux%odsun_singlelayer(lay,ist,i)
                 endif
           enddo
        enddo

      do ist = 0, nstreams(i)     
         lay = 1_jpim
         auxrad_stream%down_ref(lay,ist,i) = auxrad_stream%down(lay,ist,i) + temp(lay,ist)
         do lay = 2, nlayers
               temp(lay,ist) = temp(lay-1,ist) + temp(lay,ist) 
               auxrad_stream%down_ref(lay,ist,i) = auxrad_stream%down(lay,ist,i) + temp(lay,ist)
               auxrad_stream%down(lay,ist,i) = max(auxrad_stream%down_ref(lay,ist,i), 0.0_jprb)
         enddo
      enddo

   end if
enddo

end subroutine solar_scattering_air

!DAR: As above
subroutine solar_scattering_near_surf(transmission_aux, sun, chanprof, &
                                          auxrad_stream)
  Implicit None
  
  type(transmission_type_aux), intent(in)    :: transmission_aux
  logical(jplm), intent(in)                  :: sun(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)

  type(radiance_aux), intent(inout)          :: auxrad_stream

  integer(jpim) :: lay1 
  real(jprb) :: temp(0:totnstreams)

#define fac1_3 auxrad_stream%Fac1_3(ist,i)
#define fac2_3 auxrad_stream%Fac2_3(ist,i)
#define fac3_3 auxrad_stream%Fac3_3(ist,i)
#define fac4_3 auxrad_stream%Fac4_3(ist,i)
#define fac5_3 auxrad_stream%Fac5_3(ist,i)
#define fac6_3 auxrad_stream%Fac6_3(ist,i)
#define fac7_3 auxrad_stream%Fac7_3(ist,i)
#define dfac54_3 (fac5_3 - fac4_3)
#define dtausun_surf (transmission_aux%Tausun_level(lev,ist,i) - transmission_aux%Tausun_surf(ist,i))
#define tausun_level transmission_aux%Tausun_level(lev,ist,i) 

  do i=1, nchannels
     if(sun(i)) then     

        lay = iv3lay(i)
        lev = iv3lay(i) + 1

        if(pfraction(i) < 0.0_JPRB )then
           lay1 = iv3lay(i)
        else
           lay1 = iv3lay(i) + 1
        endif

        auxrad_stream%Fac4_3(0:nstreams(i),i) = exp(-transmission_aux%odsun_sfrac(0:nstreams(i),i))
        auxrad_stream%Fac5_3(0:nstreams(i),i) = exp(-transmission_aux%od_sfrac(0:nstreams(i),i))   
        
        do ist = 0, nstreams(i)

           fac1_3 = fac1_2
           fac2_3 = fac2_2
           fac3_3 = fac3_2
           fac6_3 = fac6_2
        
          !--------------Upward single scattering of the solar beam-------------------------
           auxrad_stream%meanrad_up(ist,i) = auxrad_stream%meanrad_up(ist,i) + &
                               fac1_3 * fac2_3 * fac3_3 * dtausun_surf                              


! assume there is no atmospheric source term for 3rd/4th stokes vector elements
           if (pol_id(i) >= 6_jpim ) auxrad_stream%meanrad_up(ist,i) = 0.0_jprb
        enddo

          !--------------Downward single scattering of the solar beam-----------------------
        do ist = 0, nstreams(i)
           temp(ist) = fac6_3 * fac2_3 * tausun_level * (tau_level_r * tau_surf_r)
        enddo

        if(.not. sateqsun(lay1,prof))then
           do ist = 0, nstreams(i)
              auxrad_stream%Fac7_3(ist,i) = auxrad_stream%Fac7_2(lay1,ist,i)
                
              auxrad_stream%meanrad_down(ist,i) = auxrad_stream%meanrad_down(ist,i) + &
                                                  transmission_aux%fac(2,lev,ist,i) * temp(ist) * fac7_3 * dfac54_3
           enddo
        else
           do ist = 0, nstreams(i)
              auxrad_stream%meanrad_down(ist,i) = auxrad_stream%meanrad_down(ist,i) + &
                                  transmission_aux%fac(2,lev,ist,i) * temp(ist) * fac4_3 * &
                                  transmission_aux%odsun_sfrac(ist,i)
           enddo
        endif
        do ist = 0, nstreams(i)
          auxrad_stream%meanrad_down(ist,i) = max(auxrad_stream%meanrad_down(ist,i), 0._jprb)
        enddo
     endif
  end do

end subroutine solar_scattering_near_surf

subroutine solar_contribution(transmission_aux, sun, sunglint, fresnrefl, reflectivity, &
                                  coef, profiles, chanprof, auxrad_stream, rad)
  Implicit None

  type(transmission_type_aux), intent(inout) :: transmission_aux
  type(sunglint_type), intent(in)            :: sunglint
  real(jprb), intent(in)                     :: fresnrefl(:)
  real(jprb), intent(in)                     :: reflectivity(:) ! surface reflectivity
  type(rttov_coef), intent(in)               :: coef
  type(profile_type), intent(in)             :: profiles(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  logical(jplm), intent(in)                  :: sun(:)
  type(radiance_aux), intent(inout)          :: auxrad_stream
  type(radiance_type), intent(inout)         :: rad
  
  real(jprb)    :: temp(0:totnstreams) ! ist

  Do i = 1, nchannels
     if(sun(i)) then
                
        temp(0:nstreams(i)) = coef%ss_solar_spectrum(chan) * transmission_aux%tausun_surf(0:nstreams(i),i)
        transmission_aux%refl_norm(i) = cos(profiles(prof)%sunzenangle * deg2rad) * pi_r

        if(profiles(prof)%skin%surftype == surftype_sea) then
           if(sunglint%s(prof)%windsp == 0._jprb) then
              do ist = 0, nstreams(i)
                temp(ist) = temp(ist) * fresnrefl(i) * sunglint%s(prof)%glint * transmission_aux%refl_norm(i)
              enddo
           else
              do ist = 0, nstreams(i)
                temp(ist) = temp(ist) * fresnrefl(i) * sunglint%s(prof)%glint
              enddo
           endif
        else ! land or sea ice
          do ist = 0, nstreams(i)
            temp(ist) = temp(ist) * reflectivity(i) * transmission_aux%refl_norm(i)
          enddo
        endif
        rad%clear(i) = rad%clear(i) + temp(0)
        auxrad_stream%cloudy(1:nstreams(i),i)= auxrad_stream%cloudy(1:nstreams(i),i) + temp(1:nstreams(i))
     endif
  enddo
end subroutine solar_contribution

subroutine init_rttov_integrate_mod_vars(chanprof, profiles, ircld, aux_prof, coef, raytracing)

  type(rttov_chanprof), intent(in)              :: chanprof(:)     ! Channel indices
  type(profile_type),   intent(in)              :: profiles(:) ! Profiles
  type(ircld_type),     intent(in)              :: ircld
  type(raytracing_type),intent(in)              :: raytracing
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef

!store calculation size
nchannels = size(chanprof)
nlayers = profiles(1)%nlayers
nlevels = nlayers + 1

! allocate arrays
! DAR: Will this be a problem for malloc locking for parallel users?
allocate(nstreams(nchannels), iv2lev(nchannels), iv2lay(nchannels), iv3lev(nchannels), iv3lay(nchannels), pol_id(nchannels), &
         cfraction(nchannels), pfraction(nchannels), sateqsun(nlayers,size(profiles)))

! populate module arrays
Do i = 1, nchannels
   nstreams(i) = ircld%nstream(prof)
enddo

do i = 1, nchannels
   cfraction(i) = aux_prof%s(prof)%cfraction
   pfraction(i) = aux_prof%s(prof)%pfraction_surf
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

do iprof = 1, size(profiles)
   do lay = 1, nlayers
      sateqsun(lay,iprof) = .false.
      if(raytracing%pathsat(lay,iprof) == raytracing%pathsun(lay,iprof)) sateqsun(lay,iprof) = .true.
   enddo
enddo

end subroutine init_rttov_integrate_mod_vars

end subroutine rttov_integrate
