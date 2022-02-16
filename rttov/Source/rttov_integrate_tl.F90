Subroutine rttov_integrate_tl( &
   addcosmic, opts, maxnstreams, chanprof,      &! in
   emissivity, emissivity_tl,    &! in
   reflectivity, reflectivity_tl,  &! in
   fresnrefl, fresnrefl_tl,     &! in
   sunglint, sunglint_tl,      &! in
   sun,              &! in
   transmission_aux, transmission_aux_tl,  &! in
   transmission_scatt_ir_stream, transmission_scatt_ir_stream_tl,&
   profiles, profiles_tl,      &! in
   aux_prof, aux_prof_tl,      &! in
   coef,     &! in
   raytracing, raytracing_tl,    &! in
   ircld, ircld_tl,         &! in
   rad , &! in
   auxrad ,          &! in
   auxrad_stream, auxrad_stream_tl, &! in
   rad_tl           ) ! inout
!
! Description:
! To perform TL of integration of radiative transfer equation
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
! ECMWF Research Dept. Tech. Memo. 425 (available from the librarian at ECMWF)
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
!  1.4     28/02/2005   Improved vectorisation (D Dent)
!  1.5     29/03/2005   Add end of header comment (J. Cameron)
!  1.6     01/06/2005   Marco Matricardi (ECMWF):
!             --        IASI capability added.
!             --        Linear in tau approximation for RT equation introduced.
!             --        Solar radiation introduced for IASI and AIRS.
!  1.7     03/01/2007   Corrected bugs in'old code' R Saunders
!  1.8     08/02/2007   Removed polarisation index (R Saunders)
!  1.9     11/03/2007   Reintroduced overcast radiance (R Saunders)
!  1.10    15/11/2007   Changed skin to surfair for overcast radiances(RSaunders)
!  1.11    27/11/2007   Optimised for NEC/IBM (D Salmond)
!  1.12    21/12/2007   Added polarimetric option (R. Saunders)
!  1.13    15/07/2009   User defined ToA. Layers distinct from levels. Top-layer
!                       brought into layer looping to shorten code (P.Rayer)
!  1.14    03/11/2009   Transmittances on levels (A Geer)
!  1.15    02/12/2009   Introduced principal component capability. Pathsat, Pathsun and
!                       related quantities are now layer arrays (Marco Matricardi).
!  1.16    05/07/2010   Remove addsolar flag from profiles structure (J Hocking)
!  1.17    14/10/2010   Remove rt8_mode (J Hocking)
!  1.18    14/12/2010   Use traj0_sta%sun array to flag channels for which solar calculations
!                       should be performed (J Hocking)
!  2.0     14/12/2011   Re-written (D Rundle)
!
! 2010/03 Code cleaning Pascal Brunel, Philippe Marguinaud
!
!
!
! Code Description:
!   Language:           Fortran 90.
!   Software Standards: "European Standards for Writing and
!     Documenting Exchangeable Fortran 90 Code".
!

  use rttov_types, Only : rttov_chanprof, rttov_coef, rttov_options, profile_type, profile_aux, transmission_type_aux, & 
       transmission_scatt_ir_type, sunglint_type, radiance_type, ircld_type, raytracing_type, radiance_aux

  use parkind1, Only : jpim, jprb, jplm
!INTF_OFF  

  use rttov_const, Only : sensor_id_po

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
  type(transmission_type_aux), intent(in)       :: transmission_aux            ! transmittances and single-layer od
  type(transmission_scatt_ir_type), intent(in)  :: transmission_scatt_ir_stream
  type(profile_aux) ,   intent(in)              :: aux_prof ! auxillary profiles info.
  type(rttov_coef),     intent(in)              :: coef
  type(radiance_aux),   intent(in)              :: auxrad_stream
  type(radiance_type),  intent(in)              :: rad    ! radiances (mw/cm-1/ster/sq.m) and BTs
  type(radiance_aux),   intent(in)              :: auxrad ! auxillary radiances

  Real(jprb),                Intent(in)         :: emissivity_tl(size(chanprof))
  Real(jprb),                Intent(in)         :: reflectivity_tl(size(chanprof))
  Real(jprb),                Intent(in)         :: fresnrefl_tl(size(chanprof))
  Type(profile_Type),  Intent(in)               :: profiles_tl(size(profiles))
  Type(ircld_type)                 ,intent(in)  :: ircld_tl
  Type(raytracing_type),intent(in)              :: raytracing_tl
  Type(sunglint_Type), Intent(in)               :: sunglint_tl
  Type(transmission_Type_aux), Intent(in)       :: transmission_aux_tl
  type(transmission_scatt_ir_type) ,intent(in)  :: transmission_scatt_ir_stream_tl
  Type(profile_aux) ,  Intent(in)               :: aux_prof_tl
  Type(radiance_Type), Intent(inout)            :: rad_tl ! in because of mem allocation
  Type(radiance_aux),  Intent(inout)            :: auxrad_stream_tl

!INTF_END

#include "rttov_calcbt_tl.h"
#include "rttov_calcrad_tl.h"

!local variables: 
  integer(jpim) :: i, lev, ist, lay, nchannels, nlayers, nlevels, narray(3), iprof, prof, chan, nstreams ! counter variables
  real(jprb)    :: cfraction(size(chanprof)), cfraction_tl(size(chanprof)), pfraction(size(chanprof))        ! cloud fraction
  logical(jplm) :: sateqsun(profiles(1)%nlayers,size(profiles(:))) ! True where the solar zenith angle equal to observation angle
  logical(jplm) :: anysun ! Are *any* channels affected by sun (only need to make check once if negative)
  real(jprb)    :: meanrad_up_tl(0:maxnstreams,size(chanprof)), meanrad_down_tl(0:maxnstreams,size(chanprof))
  integer(jpim) :: pol_id(size(chanprof))       ! polarisation index  

  Real(jprb)    :: rad_air_tl(profiles(1) % nlevels, size(chanprof)), rad_surfair_tl(size(chanprof)), rad_skin_tl(size(chanprof))
  Real(jprb)    :: zup_tl(0:profiles(1)%nlayers,0:maxnstreams,size(chanprof)), &
                   zdown_tl(0:profiles(1)%nlayers,0:maxnstreams,size(chanprof))

  integer(jpim) :: iv2lay(size(chanprof)), iv2lev(size(chanprof)), iv3lay(size(chanprof)), iv3lev(size(chanprof))
  real(jprb)    :: p1, p2

  REAL(JPRB)    :: ZHOOK_HANDLE

!- End of header ------------------------------------------------------

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_TL',0_jpim,ZHOOK_HANDLE)

#define tau_surf_tl transmission_aux_tl%Tau_surf(ist,i)
#define tau_layer_tl_r transmission_aux_tl%Tau_level_r(lay,ist,i)
#define tau_layer_tl transmission_aux_tl%Tau_level(lay,ist,i)
#define tau_level_tl_r transmission_aux_tl%Tau_level_r(lay+1,ist,i)
#define tau_level_tl transmission_aux_tl%Tau_level(lay+1,ist,i)
#define daux_lay_tl (rad_air_tl(lay+1,i) - rad_air_tl(lay,i))
#define dtau_lay_tl (tau_layer_tl - tau_level_tl)
#define od_singlelayer_tl transmission_aux_tl%Od_singlelayer(lay,ist,i)

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
   prof = chanprof(i)%prof
   cfraction(i) = aux_prof%s(prof)%cfraction
   cfraction_tl(i) = aux_prof_tl%s(prof)%cfraction
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
   prof = chanprof(i)%prof
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
      chan = chanprof(i)%chan        
      pol_id(i) = coef%fastem_polar(chan) + 1_jpim
   enddo
Else
   pol_id(:) = 0_jpim
Endif 

anysun = any(sun)

!----------------------------
!1. calculate layer radiances
!----------------------------
Call rttov_calcrad_tl( chanprof, &! in
     profiles, profiles_tl,  &! in
     coef,         &! in
     auxrad%skin, auxrad%surfair, auxrad%air,    &! in
     rad_skin_tl, rad_surfair_tl, rad_air_tl) ! out

!-------------------------------------
!2. calculate atmospheric contribution
!-------------------------------------

call calc_atmospheric_radiance_tl(transmission_aux, transmission_aux_tl, auxrad, ircld, chanprof, narray, &
                                        rad_air_tl, &
                                        zup_tl, zdown_tl)

!---Scattering of the solar beam--------------------------------------------------
   if((opts%addaerosl .OR. opts%addclouds) .and. anysun) &
        call  solar_scattering_air_tl(transmission_aux, transmission_aux_tl, ircld, raytracing, raytracing_tl, &
                                       coef, chanprof, narray, transmission_scatt_ir_stream_tl, &
                                       sun, sateqsun, auxrad_stream, zup_tl, zdown_tl)

!-------------------------------------------------------------------------------
!2a calculate near-surface layer contribution
!-------------------------------------------------------------------------------
call calc_near_surf_contribution_tl(transmission_aux, transmission_aux_tl, auxrad, ircld, & !in
                                             chanprof, narray, iv3lay, iv3lev, pol_id,& !in
                                             rad_surfair_tl, rad_air_tl, zup_tl, zdown_tl, & !inout
                                             meanrad_up_tl, meanrad_down_tl) ! inout

!---Scattering of the solar beam---------------------------
  if((opts%addaerosl .OR. opts%addclouds) .and. anysun) &
     call solar_scattering_near_surf_tl(transmission_aux, transmission_aux_tl, &
                                        ircld, raytracing, raytracing_tl, coef, transmission_scatt_ir_stream_tl, &
                                        sun, sateqsun, iv3lay, iv3lev, pol_id, chanprof, narray, pfraction, &
                                        auxrad_stream, meanrad_up_tl, meanrad_down_tl)
  do i = 1, nchannels
     prof = chanprof(i)%prof
     nstreams = ircld%nstream(prof)
   ! clear sky radiance without reflection term
   ! without surface contribution at this line  
     ist = 0_jpim
     rad_tl%reflclear(i) = tau_surf * &
                          (reflectivity(i) * &
                          (meanrad_down_tl(ist,i) * tau_surf + &
                           2._JPRB * tau_surf_tl * auxrad_stream%meanrad_down(ist,i)) + &
                           reflectivity_tl(i) * auxrad_stream%meanrad_down(ist,i) * tau_surf)

     rad_tl%clear(i) = meanrad_up_tl(ist,i) + rad_tl%reflclear(i)
     rad_tl%upclear(i) = meanrad_up_tl(ist,i)
   ! clear sky downwelling radiance

   ! reflected clear sky downwelling radiance
     zup_tl(iv2lay(i),ist,i) = meanrad_up_tl(ist,i)

     do ist = 1, nstreams
        auxrad_stream_tl%cloudy(ist,i) = meanrad_up_tl(ist,i) + &
                                         tau_surf * (reflectivity(i) * &
                                        (meanrad_down_tl(ist,i) * tau_surf + &
                                         auxrad_stream%meanrad_down(ist,i) * tau_surf_tl * 2._JPRB) + &
                                         reflectivity_tl(i) * auxrad_stream%meanrad_down(ist,i) * tau_surf)
        zup_tl(iv2lay(i),ist,i) = meanrad_up_tl(ist,i)
     enddo
  enddo

!-----------------------
!3. calculate surface contribution
!-----------------------
Do i = 1, nchannels
! clear sky radiance without reflection term
   ist = 0_jpim

  rad_tl%upclear(i) = rad_tl%upclear(i) + &
                      auxrad%skin(i) * (emissivity_tl(i) * tau_surf + emissivity(i) * tau_surf_tl) + &
                      rad_skin_tl(i) *  emissivity(i) * tau_surf

  rad_tl%clear(i) = rad_tl%clear(i) + &
                    auxrad%skin(i) * (emissivity_tl(i) * tau_surf + emissivity(i) * tau_surf_tl) + &
                    rad_skin_tl(i) *  emissivity(i) * tau_surf

  prof = chanprof(i)%prof
  do ist=1, ircld%nstream(prof)
     auxrad_stream_tl%cloudy(ist,i) = auxrad_stream_tl%cloudy(ist,i) + auxrad%skin(i) * &
                                     (emissivity_tl(i) * tau_surf     + &
                                      emissivity(i)    * tau_surf_tl) + &
                                      rad_skin_tl(i) * emissivity(i) * tau_surf ! DAR replace with rskinetau
  End Do
End Do

!--------------------------------
!4. Add solar contribution
!--------------------------------

if(coef%fmv_model_ver == 9) then
   if(anysun) call solar_contribution_tl(transmission_aux, transmission_aux_tl, &
                                         ircld, sunglint, sunglint_tl, fresnrefl, fresnrefl_tl, &
                                         reflectivity, reflectivity_tl, &
                                         coef, profiles, chanprof, sun, narray, auxrad_stream_tl, rad_tl)
endif

!--------------------------------
!5. cosmic temperature correction
!--------------------------------

!calculate planck function corresponding to tcosmic = 2.7k
!deblonde tcosmic for microwave sensors only
!DAR: addcosmic should probably go into opts structure or should be determined from sensor type/frequency?
if (addcosmic) then                          
   do i = 1,nchannels 
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)

      ist = 0_jpim

      rad_tl%clear(i) = rad_tl%clear(i) + &
                        auxrad%cosmic(i) * tau_surf * &
                       (reflectivity_tl(i) * tau_surf + &
                        reflectivity(i) * 2.0_jprb * tau_surf_tl)

      do ist = 1, nstreams
         auxrad_stream_tl%cloudy(ist,i) =  auxrad_stream_tl%cloudy(ist,i) + &
                                           auxrad%cosmic(i) * tau_surf * &
                                          (reflectivity_tl(i) * tau_surf + &
                                           reflectivity(i) * 2.0_jprb * tau_surf_tl)
      enddo
   enddo
Endif

!---------------------------------------------------
!6. calculate overcast radiances
!---------------------------------------------------
!---------------
!6.1 Upward part
!---------------


ist = 0_jpim
!cdir nodep
Do i = 1, nchannels
   do lay = 1,nlayers
      lev = lay + 1
! overcast radiances at given cloud top
      rad_tl%overcast(lay,i) = zup_tl(lay,ist,i) + &
                               rad_air_tl(lev,i) * tau_level + &
                               auxrad%air(lev,i) * tau_level_tl
   end do
enddo

! Add surface component to overcast radiances
do i = 1, nchannels
   lay = iv2lay(i)

   rad_tl%overcast(lay,i) = zup_tl(lay,ist,i) + &
                            tau_surf_tl * auxrad%surfair(i) + &
                            tau_surf    * rad_surfair_tl(i)
enddo

! DAR: can this be just addclouds? Also, is it necessary to calculate overcast radiances when using RTTOV-CLOUDY? 
if(opts%addaerosl .OR. opts%addclouds) then 
!----------------------------------------
!7. calculate complex cloudy radiances
!----------------------------------------
   rad_tl%cloudy = 0._jprb ! nchannels
   Do i = 1, nchannels
      prof = chanprof(i)%prof

      do ist = 1, ircld%nstream(prof)
         rad_tl%cloudy(i) = rad_tl%cloudy(i) + & 
                            auxrad_stream_tl%cloudy(ist,i) * (ircld%xstr(ist+1,prof)    - ircld%xstr(ist,prof)   ) + &
                            auxrad_stream%cloudy(ist,i)    * (ircld_tl%xstr(ist+1,prof) - ircld_tl%xstr(ist,prof))
      enddo

      rad_tl%cloudy(i)= rad_tl%cloudy(i) + &
                        rad_tl%clear(i) * ircld%XSTRCLR(prof) + &
                        rad%clear(i)    * ircld_tl%XSTRCLR(prof)
   End do
   
!---------------------------
!8. calculate total radiance (cloudy case)
!---------------------------
   rad_tl%total(1:nchannels) = rad_tl%cloudy(1:nchannels)
else
!6. Calculate "simple" cloudy radiances
!--------------------------------------------
!6.2 Interpolate to given cloud-top pressures
!--------------------------------------------
do i = 1, nchannels
   prof = chanprof(i)%prof
   lay = aux_prof%s(prof)%nearestlev_ctp - 1

   p1   = aux_prof%s(chanprof(i)%prof)%pfraction_ctp
   p2   = aux_prof_tl%s(chanprof(i)%prof)%pfraction_ctp

   rad_tl%cloudy(i) = rad_tl%overcast(lay,i) + &
                     (rad_tl%overcast(lay-1,i) - rad_tl%overcast(lay,i)) * p1 + &
                     (rad%overcast(lay-1,i)    - rad%overcast(lay,i)   ) * p2

End Do

!---------------------------
!8. calculate total radiance (clear case)
!---------------------------
   if(opts%addpc) then
      rad_tl%total(1:nchannels) = rad_tl%clear(1:nchannels)
   else     
      rad_tl%total(1:nchannels) = rad_tl%clear(1:nchannels) + &
                                  cfraction(1:nchannels)    * (rad_tl%cloudy(1:nchannels) - rad_tl%clear(1:nchannels) ) + &
                                  cfraction_tl(1:nchannels) * (rad%cloudy(1:nchannels)    - rad%clear(1:nchannels)    )
   endif
Endif

!-----------------------------------------------
!9. convert radiances to brightness temperatures
!-----------------------------------------------
Call rttov_calcbt_tl(chanprof, coef, rad, &! in
                     rad_tl)    ! inout

if (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATE_TL',1_jpim,ZHOOK_HANDLE)

contains

subroutine calc_atmospheric_radiance_tl(transmission_aux, transmission_aux_tl, auxrad, ircld, chanprof, narray, &
                                        rad_air_tl, zup_tl, zdown_tl)

  use rttov_types, Only : transmission_type_aux, radiance_aux, ircld_type, rttov_chanprof
  use parkind1, only : jpim, jprb

  Implicit None

  type(transmission_type_aux), intent(in)       :: transmission_aux, transmission_aux_tl
  type(radiance_aux), intent(in)                :: auxrad
  type(ircld_type),     intent(in)              :: ircld
  type(rttov_chanprof), intent(in)              :: chanprof(:)
  integer(jpim), intent(in)                     :: narray(:)
  Real(jprb), intent(in) :: rad_air_tl(:,:)
  Real(jprb), intent(out) :: zup_tl(0:,0:,:), zdown_tl(0:,0:,:)

  integer(jpim) :: nlevels, nlayers, nstreams, nchannels, prof
  integer(jpim) :: i, ist, lay

  Real(jprb) :: tau_lev_r_tl,tau_lev_rm1_tl,od_singlelayer_r_tl

  Real(jprb) :: b1_2,b2_2,b3_2
  Real(jprb) :: b1_tl, b2_tl, b3_tl

!unpack narray
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)

  if(transmission_aux%anynegtau .gt. 0.0_jprb) then
     
     do i = 1, nchannels
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)
      do ist = 0, nstreams
         do lay = 1, nlayers

            B1_2 = auxrad%air(lay,i) * dtau_lay
            B2_2 = daux_lay * tau_level
            B3_2 = daux_lay * dtau_lay * od_singlelayer_r

            od_singlelayer_r_tl = -od_singlelayer_tl * od_singlelayer_r**2
            tau_lev_r_tl   = -tau_level_tl * tau_level_r**2
            tau_lev_rm1_tl = -tau_layer_tl * tau_layer_r**2

            B1_TL = rad_air_tl(lay,i) * dtau_lay + auxrad%air(lay,i) * dtau_lay_tl
            B2_TL = daux_lay_tl       * tau_level  + &
                    daux_lay          * tau_level_tl
            B3_TL = od_singlelayer_r  * &
                   (daux_lay_tl * dtau_lay + &
                    daux_lay * dtau_lay_tl) + &
                    od_singlelayer_r_tl  * (daux_lay * dtau_lay)

            zup_tl(lay,ist,i)   = transmission_aux%fac(1,lay,ist,i) * (B1_TL - B2_TL + B3_TL)
            zdown_tl(lay,ist,i) = transmission_aux%fac(1,lay,ist,i) * transmission_aux%fac(2,lay+1,ist,i) * &
                                 (tau_level_r * ( &
                                 (B1_TL - B3_TL) * (tau_layer_r) + &
                                 (B2_TL * tau_level_r          + &
                                  B2_2 * (2._jprb * tau_lev_r_tl)))                          + &
                                 (B1_2 - B3_2)  * (tau_lev_rm1_tl * tau_level_r + &
                                  tau_lev_r_tl * tau_layer_r))
         enddo
      enddo
   enddo

! DAR: This is slow. Hopefully it'll never run.
   if(transmission_aux%anynegtau .gt. 0_jprb) then
      do  i = 1,nchannels
         prof = chanprof(i)%prof
         nstreams = ircld%nstream(prof)
         do ist = 0, nstreams
            do lay = 1,nlayers
               if(tau_level < 0.0_jprb .and. tau_layer >= 0._jprb) then 
                  zup_tl(lay, ist, i) = 0.5_JPRB  * &
                                      ((auxrad%air(lay+1,i) + auxrad%air(lay,i)) * &
                                       (tau_layer_tl - tau_level_tl) + &
                                       (rad_air_tl(lay+1,i) + rad_air_tl(lay,i)) * &
                                       (tau_layer - tau_level))
                  zdown_tl(lay, ist, i) = 0.0_jprb
               endif
            enddo
         enddo
      enddo
   endif

!Cumulative sum
   do i = 1,nchannels
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)
      do ist = 0, nstreams
         do lay = 2, nlayers
            zup_tl(lay,ist,i) = zup_tl(lay,ist,i) + zup_tl(lay-1,ist,i) 
            zdown_tl(lay,ist,i) = zdown_tl(lay,ist,i) + zdown_tl(lay-1,ist,i) 
         enddo
      enddo
   enddo
else ! DAR: fast code (on IBM) when no special cases
   do i = 1, nchannels
      prof = chanprof(i)%prof
      do ist = 0, ircld%nstream(prof)
         zup_tl(0,ist,i) = 0._jprb
         zdown_tl(0,ist,i) = 0._jprb
         do lay = 1, nlayers

            B1_2 = auxrad%air(lay,i) * dtau_lay
            B2_2 = daux_lay * tau_level
            B3_2 = daux_lay * dtau_lay * od_singlelayer_r
            
            od_singlelayer_r_tl = -od_singlelayer_tl * od_singlelayer_r**2
            tau_lev_r_tl   = -tau_level_tl * tau_level_r**2
            tau_lev_rm1_tl = -tau_layer_tl   * tau_layer_r**2

            B1_TL = rad_air_tl(lay,i) * dtau_lay + auxrad%air(lay,i) * dtau_lay_tl
            B2_TL = daux_lay_tl       * tau_level  + &
                    daux_lay          * tau_level_tl
            B3_TL = od_singlelayer_r  * &
                   (daux_lay_tl * dtau_lay + &
                    daux_lay * dtau_lay_tl) + &
                    od_singlelayer_r_tl  * (daux_lay * dtau_lay)
            ! zup/down_tl(0...) = 0
            zup_tl(lay,ist,i)   = zup_tl(lay-1,ist,i) + transmission_aux%fac(1,lay,ist,i) * (B1_TL - B2_TL + B3_TL) 
            zdown_tl(lay,ist,i) = zdown_tl(lay-1,ist,i)+ transmission_aux%fac(1,lay,ist,i) * &
                 transmission_aux%fac(2,lay+1,ist,i) * &
                 (tau_level_r * ( &
                 (B1_TL - B3_TL) * (tau_layer_r) + &
                 (B2_TL * tau_level_r          + &
                 B2_2 * (2._jprb * tau_lev_r_tl)))                          + &
                 (B1_2 - B3_2)  * (tau_lev_rm1_tl * tau_level_r + &
                 tau_lev_r_tl * tau_layer_r))
         enddo
      enddo
   enddo
endif

end subroutine calc_atmospheric_radiance_tl

! DAR: This is the next biggest consumer of CPU time and should be looked at next. The trouble seems to be that you have to change a
!      lot of non-consecutive data. Maybe the RTTOV 9 comments can be removed
subroutine calc_near_surf_contribution_tl(transmission_aux, transmission_aux_tl, auxrad, ircld, & !in
                                       chanprof, narray, iv3lay, iv3lev, pol_id, & !in
                                       rad_surfair_tl, rad_air_tl, zup_tl, zdown_tl, &
                                       meanrad_up_tl, meanrad_down_tl) ! out

  use rttov_const, only : min_od
  use rttov_types, Only : rttov_chanprof, ircld_type, transmission_type_aux, radiance_aux
  use parkind1, only : jpim, jprb

  Implicit None

  type(transmission_type_aux), intent(in) :: transmission_aux, transmission_aux_tl
  type(radiance_aux), intent(in)          :: auxrad
  type(ircld_type), intent(in)            :: ircld
  type(rttov_chanprof), intent(in)        :: chanprof(:)
  integer(jpim), intent(in)               :: narray(:)
  integer(jpim), intent(in)               :: iv3lay(:), iv3lev(:)
  integer(jpim), intent(in)               :: pol_id(:)
  real(jprb), intent(in)                  :: rad_surfair_tl(:), rad_air_tl(:,:)
  real(jprb), intent(in)                  :: zup_tl(0:,0:,:), zdown_tl(0:,0:,:)

  real(jprb), intent(inout)               :: meanrad_up_tl(0:,:), meanrad_down_tl(0:,:)

  integer(jpim) :: nlayers, nlevels, nstreams, nchannels
  integer(jpim) :: prof

  integer(jpim) :: i, lev, lay, ist
 
  real(jprb) ::  rad_tmp, rad_tmp_tl
  
#define B1_3 (auxrad%air(lev,i) * (tau_level  - tau_surf))
#define B2_3 ((auxrad%surfair(i) - auxrad%air(lev,i)) * tau_surf)
#define B3_3 ((auxrad%surfair(i) - auxrad%air(lev,i)) * (tau_level - tau_surf) * (transmission_aux%od_sfrac_r(ist,i)))             

#define B1_TL1 (rad_air_tl(lev,i) * (tau_level  - tau_surf))
#define B1_TL2 (auxrad%air(lev,i) * (tau_level_tl - tau_surf_tl))
#define B2_TL1 (tau_surf * (rad_surfair_tl(i) - rad_air_tl(lev,i)))
#define B2_TL2 (tau_surf_tl * (auxrad%surfair(i) - auxrad%air(lev,i)))
#define B3_TL1 ((auxrad%surfair(i) - auxrad%air(lev,i)) * (tau_level_tl - tau_surf_tl))
#define B3_TL2 ((rad_surfair_tl(i) - rad_air_tl(lev,i)) * (tau_level - tau_surf))
#define B3_TL3 (-B3_3 * transmission_aux_tl%od_sfrac(ist,i))
#define B3_TL4 (transmission_aux%od_sfrac_r(ist,i))

! z=x/y -=> z' = x'-zy' 

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)  
  
  Do i = 1, nchannels
     prof = chanprof(i)%prof
     nstreams = ircld%nstream(prof)

     lay = iv3lay(i)
     lev = iv3lev(i)
  
     do ist = 0, nstreams !ircld%nstreams(prof)

        if(tau_surf < 0._JPRB) then ! DAR 0 < min_tau
           if(tau_level >= 0._JPRB) then
   
              meanrad_up_tl(ist,i) = 0.5_JPRB * (((rad_surfair_tl(i) + rad_air_tl(lev,i)) * &
                                    (tau_level - tau_surf)) + &
                                   ((auxrad%surfair(i)     + auxrad%air(lev,i)) * &
                                    (tau_level_tl - tau_surf_tl)))

              if (pol_id(i) >= 6_jpim ) meanrad_up_tl(ist,i) = 0.0_jprb
           endif
           
           meanrad_down_tl(ist,i) = zdown_tl(lay,ist,i)
           meanrad_up_tl(ist,i)   = zup_tl(lay,ist,i) + meanrad_up_tl(ist,i)
           
        else !tau_level<0 .or. transmission_aux%tau_surf(ist,i) >= 0

           if(transmission_aux%od_sfrac(ist,i) < min_od .OR. &
                ((tau_level - tau_surf) < min_od)) THEN

              ! small optical depth or optical depth change set radiance to zero

              meanrad_down_tl(ist,i) = zdown_tl(lay,ist,i)
              meanrad_up_tl(ist,i)   = zup_tl(lay,ist,i)

           else
              ! compute linear in tau radiances
              ! PGF90 compatibility for b3_tl macros
              meanrad_up_tl(ist,i) = &
                   B1_TL1 + &
                   B1_TL2 - &
                  (B2_TL1 + &
                   B2_TL2) + &
                   (B3_TL1 + & 
                   B3_TL2 + &
                   B3_TL3) * &
                   B3_TL4
                 
              rad_tmp = transmission_aux%surf_fac(ist,i) * &
                       (B1_3 - &
                        B3_3) * &
                        (tau_level_r * tau_surf_r) + &
                            B2_3 * (tau_surf_r)**2_jpim
              
              rad_tmp_tl = tau_surf_r * ( &                                
                  (tau_level * ((B1_TL1 + &
                   B1_TL2) - &
                   (B3_TL1 + &
                   B3_TL2 + &
                   B3_TL3) * &
                   B3_TL4) - &
                   (B1_3 - &
                   B3_3) * tau_level_tl) * &
                   tau_level_r**2_jpim + &
                   ((B2_TL1 + &
                     B2_TL2) * &
                     tau_surf - B2_3 * tau_surf_tl) * &
                   tau_surf_r**2_jpim - &
                   rad_tmp * tau_surf_tl) 
                 ! f = 1/z * (x/y + w/z) => f' = 1/z * ((yx'-xy')/y^2 + (wz'-zw')/z^2 - f)
                    
              if (pol_id(i) >= 6_jpim ) meanrad_up_tl(ist,i) = 0.0_jprb

              meanrad_down_tl(ist,i) = zdown_tl(lay,ist,i) + rad_tmp_tl
              meanrad_up_tl(ist,i)   = zup_tl(lay,ist,i)   + meanrad_up_tl(ist,i)

           endif
        endif
        if (pol_id(i) >= 6_jpim) meanrad_up_tl(ist,i) = 0.0_jprb ! meanrad_down not zeroed DAR
     enddo
  enddo
  
end subroutine calc_near_surf_contribution_tl

subroutine solar_scattering_air_tl(transmission_aux, transmission_aux_tl, ircld, raytracing, raytracing_tl, &
                                       coef, chanprof, narray, transmission_scatt_ir_stream_tl, &
                                       sun, sateqsun, auxrad_stream, zup_tl, zdown_tl)

  use rttov_const, only : z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, raytracing_type, &
                          rttov_coef, rttov_chanprof, transmission_scatt_ir_type
  use parkind1, only : jpim, jplm, jprb

  Implicit None

  type(transmission_type_aux), intent(in)      :: transmission_aux, transmission_aux_tl
  type(ircld_type),     intent(in)             :: ircld
  type(raytracing_type), intent(in)            :: raytracing, raytracing_tl
  type(rttov_coef),     intent(in)             :: coef
  type(radiance_aux), intent(in)               :: auxrad_stream
  type(rttov_chanprof), intent(in)             :: chanprof(:)
  type(transmission_scatt_ir_type), intent(in) :: transmission_scatt_ir_stream_tl
  logical(jplm), intent(in)                    :: sun(:)
  logical(jplm), intent(in)                    :: sateqsun(:,:)    
  integer(jpim), intent(in)                    :: narray(:)
  real(jprb), intent(inout)                    :: zup_tl(0:,0:,:), zdown_tl(0:,0:,:)

!local variables
  integer(jpim) :: chan, prof
  integer(jpim) :: nlevels, nlayers, nstreams, maxnstreams, nchannels
  integer(jpim) :: i, ist, lay

  real(jprb)    :: temp_tl(1:narray(1)-1_jpim,0:narray(2))

  Real(jprb) :: fac1_2_tl(1:narray(1)-1_jpim,0:narray(2)), fac2_2_tl(1:narray(1)-1_jpim,0:narray(2)), &
                fac3_2_tl(1:narray(1)-1_jpim), fac4_2_tl(1:narray(1)-1_jpim,0:narray(2)), &
                fac5_2_tl(1:narray(1)-1_jpim,0:narray(2)), fac6_2_tl(1:narray(1)-1_jpim,0:narray(2)), &
                fac7_2_tl(1:narray(1)-1_jpim)

! macro definitions for different factors

#define fac1_2 auxrad_stream%Fac1_2(lay,ist,i)
#define fac2_2 auxrad_stream%Fac2_2(lay,ist,i)
#define fac3_2 auxrad_stream%Fac3_2(lay,ist,i)
#define fac4_2 auxrad_stream%Fac4_2(lay,ist,i)
#define fac5_2 auxrad_stream%Fac5_2(lay,ist,i)
#define fac6_2 auxrad_stream%Fac6_2(lay,ist,i)
#define fac7_2 auxrad_stream%Fac7_2(lay,ist,i)

#define fac1_2_tl Fac1_2_tl(lay,ist)
#define fac2_2_tl Fac2_2_tl(lay,ist)
#define fac3_2_tl Fac3_2_tl(lay)
#define fac4_2_tl Fac4_2_tl(lay,ist)
#define fac5_2_tl Fac5_2_tl(lay,ist)
#define fac6_2_tl Fac6_2_tl(lay,ist)
#define fac7_2_tl Fac7_2_tl(lay)

#define tausun_layer_tl transmission_aux_tl%Tausun_level(lay,ist,i)
#define tausun_level_tl transmission_aux_tl%Tausun_level(lay+1,ist,i)
#define tausun_layer transmission_aux%Tausun_level(lay,ist,i)
#define tausun_level transmission_aux%Tausun_level(lay+1,ist,i)

#define dfac54_2_tl (fac5_2_tl - fac4_2_tl)
#define dfac54_2 (fac5_2 - fac4_2)
    
! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  maxnstreams = narray(2)
  nchannels = narray(3)
  
  do i = 1, nchannels
    if(sun(i)) then

      chan = chanprof(i)%chan
      prof = chanprof(i)%prof
      nstreams = ircld%nstream(prof)
      
      do ist = 0, nstreams
        do lay = 1, nlayers
          fac6_2_tl = coef%ss_solar_spectrum(chan) * z4pi_r * &
                      transmission_scatt_ir_stream_tl%azphacdo(ist,i,lay)
          
          fac1_2_tl = coef%ss_solar_spectrum(chan) * z4pi_r * &
                      transmission_scatt_ir_stream_tl%azphacup(ist,i,lay)

          fac2_2_tl = transmission_scatt_ir_stream_tl%ssa(ist,i,lay)

          fac4_2_tl = -transmission_aux_tl%odsun_singlelayer(lay,ist,i) * fac4_2
          fac5_2_tl = -transmission_aux_tl%Od_singlelayer(lay,ist,i) * fac5_2
        enddo
      enddo

      do lay = 1,nlayers
        fac3_2_tl = (raytracing_tl%pathsat(lay,prof) * raytracing%pathsun(lay,prof) - &
                     raytracing%pathsat(lay,prof) * raytracing_tl%pathsun(lay,prof)) / &
                    (raytracing%pathsun(lay,prof) + raytracing%pathsat(lay,prof))**2_jpim
      enddo

  !----------------Upward single scattering of the solar beam-----------------------
      do ist = 0, nstreams
        do lay = 1, nlayers
          temp_tl(lay,ist) = (tausun_layer - tausun_level) * &
            (fac1_2_tl * fac2_2 * fac3_2 + & 
              fac2_2_tl * fac3_2 * fac1_2 + & 
              fac3_2_tl * fac1_2 * fac2_2) + &
            (fac1_2 * fac2_2 * fac3_2 * &
            (tausun_layer_tl - tausun_level_tl))
        enddo
      enddo
        
      do ist = 0, nstreams
        do lay = 2, nlayers
          temp_tl(lay,ist) = temp_tl(lay,ist) + temp_tl(lay-1,ist) 
        enddo
      enddo

      do ist = 0, nstreams
        zup_tl(1:nlayers,ist,i) = zup_tl(1:nlayers,ist,i) + temp_tl(1:nlayers,ist) ! because zup_tl goes from 0
      enddo
        
  !-------------------Downward single scattering of the solar beam------------------

      do lay = 1, nlayers
        if(.not. sateqsun(lay,prof)) then
          fac7_2_tl = (raytracing%pathsun(lay,prof) * raytracing_tl%pathsat(lay,prof) - &
                       raytracing_tl%pathsun(lay,prof) * raytracing%pathsat(lay,prof)) / &
                      (raytracing%pathsun(lay,prof) - raytracing%pathsat(lay,prof))**2_jpim
        endif
      enddo

      do ist = 0, nstreams
        do lay = 1, nlayers
          if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb)then
            zdown_tl(lay,ist,i) = 0._jprb
            cycle
          endif
! nested multiplication of derivative
          if(.not. sateqsun(lay,prof)) then
            temp_tl(lay,ist) = &
              transmission_aux%fac(2,lay+1,ist,i) * &
              (1.0_jprb / (tau_layer * tau_level)) * (&
                fac6_2_tl * &
                fac2_2 * &
                fac7_2 * &
                dfac54_2 * &
                tausun_layer + &
                fac6_2 * (&
                  fac2_2_tl * &
                  fac7_2 * &
                  dfac54_2 * &
                  tausun_layer + &
                  fac2_2 * (&
                    fac7_2_tl * &
                    dfac54_2 * &
                    tausun_layer + &
                    fac7_2 * (&
                      dfac54_2_tl * &
                      tausun_layer + &
                      dfac54_2 * (&
                        (tausun_layer_tl - &
                          (tausun_layer/&
                          (tau_layer * tau_level))*&
                            (tau_layer_tl * tau_level + &
                            tau_layer * tau_level_tl)))))))
          else
            temp_tl(lay,ist) = &
              transmission_aux%fac(2,lay+1,ist,i) * &
              (1.0_jprb / (tau_layer * tau_level)) * (&
               fac6_2_tl * &
               fac2_2 * &
               fac4_2 * &
               transmission_aux%odsun_singlelayer(lay,ist,i) * &
               tausun_layer + &
               fac6_2 * (&
                fac2_2_tl * &
                fac4_2 * &
                transmission_aux%odsun_singlelayer(lay,ist,i) * &
                tausun_layer + &
                fac2_2 * (&
                 fac4_2_tl * &
                 transmission_aux%odsun_singlelayer(lay,ist,i) * &
                 tausun_layer + &
                 fac4_2 * (&
                  transmission_aux_tl%odsun_singlelayer(lay,ist,i) * &
                  tausun_layer + &
                  dfac54_2 * (&
                   (tausun_layer_tl - &
                   (tausun_layer/&
                    (tau_layer * tau_level))*&
                    (tau_layer_tl * tau_level + &
                     tau_layer * tau_level_tl)))))))
          endif
        enddo
      enddo

      do ist = 0, nstreams     
        lay = 1_jpim
        if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) then
          zdown_tl(lay,ist,i) = 0._jprb
        else
          zdown_tl(lay,ist,i) = zdown_tl(lay,ist,i) + temp_tl(lay,ist)              
        endif
            
        do lay = 2, nlayers
          temp_tl(lay,ist) = temp_tl(lay-1,ist) + temp_tl(lay,ist) 
          if(auxrad_stream%down_ref(lay,ist,i) < 0._jprb) then
            zdown_tl(lay,ist,i) = max(zdown_tl(lay,ist,i), 0.0_jprb)
          else
            zdown_tl(lay,ist,i) = zdown_tl(lay,ist,i) + temp_tl(lay,ist)
          endif
        enddo
      enddo
        
    endif
  enddo
   
end subroutine solar_scattering_air_tl

subroutine solar_scattering_near_surf_tl(transmission_aux, transmission_aux_tl, &
                                         ircld, raytracing, raytracing_tl, coef, transmission_scatt_ir_stream_tl, &
                                         sun, sateqsun, iv3lay, iv3lev, pol_id, chanprof, narray, pfraction, &
                                         auxrad_stream, meanrad_up_tl, meanrad_down_tl)

  use rttov_const, only : z4pi_r
  use rttov_types, Only : transmission_type_aux, ircld_type, radiance_aux, &
                          rttov_chanprof, rttov_coef, transmission_scatt_ir_type, raytracing_type
  use parkind1, only : jpim, jprb, jplm 

  Implicit None
  
  type(transmission_type_aux), intent(in)      :: transmission_aux, transmission_aux_tl
  type(ircld_type),     intent(in)             :: ircld
  type(raytracing_type), intent(in)            :: raytracing, raytracing_tl
  type(rttov_coef),     intent(in)             :: coef
  type(transmission_scatt_ir_type), intent(in) :: transmission_scatt_ir_stream_tl
  logical(jplm), intent(in)                    :: sun(:)
  logical(jplm), intent(in)                    :: sateqsun(:,:)
  integer(jpim), intent(in)                    :: iv3lay(:), iv3lev(:)
  integer(jpim), intent(in)                    :: pol_id(:)       ! polarisation index  
  type(rttov_chanprof), intent(in)             :: chanprof(:)
  integer(jpim), intent(in)                    :: narray(:)
  real(jprb), intent(in)                       :: pfraction(:)

  type(radiance_aux), intent(in)               :: auxrad_stream
  real(jprb),intent(inout)                     :: meanrad_up_tl(0:,:), meanrad_down_tl(0:,:)

  integer(jpim) :: prof, chan
  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist, lay, lev, lay1
  
  Real(jprb) :: fac1_3_tl(0:narray(2)), fac2_3_tl(0:narray(2)), fac3_3_tl, fac4_3_tl(0:narray(2)), &
                fac5_3_tl(0:narray(2)), fac6_3_tl(0:narray(2)), fac7_3_tl(0:narray(2))

#define fac1_3 auxrad_stream%Fac1_3(ist,i)
#define fac2_3 auxrad_stream%Fac2_3(ist,i)
#define fac3_3 auxrad_stream%Fac3_3(ist,i)
#define fac4_3 auxrad_stream%Fac4_3(ist,i)
#define fac5_3 auxrad_stream%Fac5_3(ist,i)
#define fac6_3 auxrad_stream%Fac6_3(ist,i)
#define fac7_3 auxrad_stream%Fac7_3(ist,i)

#define dfac54_3_tl (fac5_3_tl(ist) - fac4_3_tl(ist))
#define dfac54_3 (fac5_3 - fac4_3)

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nstreams = narray(2)
  nchannels = narray(3)  

  do i=1, nchannels
     if(sun(i)) then
        prof = chanprof(i)%prof
        chan = chanprof(i)%chan
        nstreams = ircld%nstream(prof)
        lay = iv3lay(i)
        lev = iv3lev(i)

        if(pfraction(i) < 0.0_JPRB) then
           lay1 = iv3lay(i)
        else
           lay1 = iv3lay(i) + 1_jpim
        endif
        
        do ist = 0, nstreams
          fac4_3_tl(ist) = -transmission_aux_tl%odsun_sfrac(ist,i) * auxrad_stream%Fac4_3(ist,i) !(exp(-x))'=-exp(-x)
          fac5_3_tl(ist) = -transmission_aux_tl%od_sfrac(ist,i) * auxrad_stream%Fac5_3(ist,i)

          fac6_3_tl(ist) = coef%ss_solar_spectrum(chan) * z4pi_r * transmission_scatt_ir_stream_tl%azphacdo(ist,i,lay)
            
          fac1_3_tl(ist) = coef%ss_solar_spectrum(chan) * z4pi_r * transmission_scatt_ir_stream_tl%azphacup(ist,i,lay)
  
          fac2_3_tl(ist) = transmission_scatt_ir_stream_tl%ssa(ist,i,lay)
        enddo
        
        fac3_3_tl = (raytracing_tl%pathsat(lay,prof) * raytracing%pathsun(lay,prof) - &
                     raytracing%pathsat(lay,prof)    * raytracing_tl%pathsun(lay,prof)) / &
                    (raytracing%pathsun(lay,prof) + raytracing%pathsat(lay,prof))**2_jpim
 
!--------------Upward single scattering of the solar beam-------------------------      
        do ist = 0, nstreams       
           meanrad_up_tl(ist,i) = meanrad_up_tl(ist,i) + &
                                  fac1_3_tl(ist) * fac2_3 * fac3_3 * &
                                 (tausun_level - transmission_aux%tausun_surf(ist,i)) + &
                                  fac1_3 * (fac2_3_tl(ist) * fac3_3 * &
                                 (tausun_level - transmission_aux%tausun_surf(ist,i)) + &
                                  fac2_3 * (fac3_3_tl  * &
                                 (tausun_level - transmission_aux%tausun_surf(ist,i)) + &
                                  fac3_3 * (&
                                 (tausun_level_tl - transmission_aux_tl%tausun_surf(ist,i)))))
        enddo

        if (pol_id(i) >= 6_jpim ) meanrad_up_tl(ist,i) = 0.0_jprb

!--------------Downward single scattering of the solar beam-----------------------
        if(.not. sateqsun(lay1,prof))then
           do ist = 0, nstreams
              fac7_3_tl(ist) = (raytracing%pathsun(lay1,prof)    * raytracing_tl%pathsat(lay1,prof)  - &
                                raytracing_tl%pathsun(lay1,prof) * raytracing%pathsat(lay1,prof)   ) / &
                               (raytracing%pathsun(lay1,prof) - raytracing%pathsat(lay1,prof))**2

              meanrad_down_tl(ist,i) = meanrad_down_tl(ist,i) + transmission_aux%fac(2,lev,ist,i) * &
                                      (tau_level_r * tau_surf_r) * (&
                                       fac6_3_tl(ist) * fac2_3 * fac7_3 * dfac54_3 * tausun_level + &
                                       fac6_3 * (&
                                                fac2_3_tl(ist) * fac7_3 * dfac54_3 * tausun_level + &
                                                        fac2_3 * (&
                                                         fac7_3_tl(ist) * dfac54_3 * tausun_level + &
                                                                 fac7_3 * (&
                                                                       dfac54_3_tl * tausun_level + &
                                                                          dfac54_3 * (&
                                      (tausun_level_tl - &
                                      (tausun_level * (tau_level_r * tau_surf_r))*&
                                      (tau_level_tl * tau_surf + &
                                       tau_level * tau_surf_tl)))))))
           enddo
        else
           do ist = 0, nstreams
              meanrad_down_tl(ist,i) = meanrad_down_tl(ist,i) + transmission_aux%fac(2,lev,ist,i) * &
                                      (tau_level_r * tau_surf_r) * (&
                                       fac6_3_tl(ist) * fac2_3 * fac4_3 * transmission_aux%odsun_sfrac(ist,i) * &
                                       tausun_level + &
                                               fac6_3 * (&
                                                fac2_3_tl(ist) * fac4_3 * transmission_aux%odsun_sfrac(ist,i) * &
                                                tausun_level + &
                                                        fac2_3 * (&
                                                         fac4_3_tl(ist) * transmission_aux%odsun_sfrac(ist,i) * &
                                                         tausun_level + &
                                                                 fac4_3 * (&
                                                                       transmission_aux_tl%odsun_sfrac(ist,i) * &
                                                                       tausun_level + &
                                                                       transmission_aux%odsun_sfrac(ist,i) * (&
                                      (tausun_level_tl - &
                                      (tausun_level/(tau_level * tau_surf)) * &
                                      (tau_level_tl * tau_surf + &
                                      tau_level * tau_surf_tl)))))))

           enddo
        endif
        do ist = 0, nstreams
            if(auxrad_stream%meanrad_down(ist,i) < 0._jprb)then
              meanrad_down_tl(ist,i) = 0._jprb
            endif
        enddo

     endif

  end do
  
end subroutine solar_scattering_near_surf_tl

! Currently nothing is being done for sea ice.       
subroutine solar_contribution_tl(transmission_aux, transmission_aux_tl, &
                                     ircld, sunglint, sunglint_tl, fresnrefl, fresnrefl_tl, reflectivity, reflectivity_tl, &
                                     coef, profiles, chanprof, sun, narray, auxrad_stream_tl, rad_tl)

  use rttov_const, only : surftype_land, surftype_sea, surftype_seaice, deg2rad, pi_r
  use rttov_types, Only : rttov_chanprof, ircld_type, rttov_coef, profile_type, &
                          transmission_type_aux, sunglint_type, radiance_type, radiance_aux
  use parkind1, only : jpim, jprb, jplm 

  Implicit None

  type(transmission_type_aux), intent(in)    :: transmission_aux, transmission_aux_tl
  type(ircld_type),     intent(in)           :: ircld
  type(sunglint_type), intent(in)            :: sunglint, sunglint_tl
  real(jprb), intent(in)                     :: fresnrefl(:), fresnrefl_tl(:)
  real(jprb), intent(in)                     :: reflectivity(:), reflectivity_tl(:) ! surface reflectivity
  type(rttov_coef), intent(in)               :: coef
  type(profile_type), intent(in)             :: profiles(:)
  type(rttov_chanprof), intent(in)           :: chanprof(:)
  logical(jplm), intent(in)                  :: sun(:)
  integer(jpim), intent(in)                  :: narray(:)

  type(radiance_type), intent(inout)         :: rad_tl
  type(radiance_aux), intent(inout)          :: auxrad_stream_tl

  
  integer(jpim) :: chan, prof
  integer(jpim) :: nlevels, nlayers, nstreams, nchannels
  integer(jpim) :: i, ist
  real(jprb)    :: temp_tl(0:narray(2)) ! ist

! unpack
  nlevels = narray(1)
  nlayers = nlevels - 1_jpim
  nchannels = narray(3)  

  Do i = 1, nchannels
     if(sun(i)) then
        prof = chanprof(i)%prof
        chan = chanprof(i)%chan
        nstreams = ircld%nstream(prof)

        if(profiles(prof)%skin%surftype == surftype_sea) then
           do ist = 0, nstreams
              temp_tl(ist) = coef%ss_solar_spectrum(chan) * &
                            (fresnrefl_tl(i) * sunglint%s(prof)%glint * transmission_aux%tausun_surf(ist,i)  + &
                            sunglint_tl%s(prof)%glint * fresnrefl(i) * transmission_aux%tausun_surf(ist,i)  + &
                            transmission_aux_tl%tausun_surf(ist,i) * fresnrefl(i) * sunglint%s(prof)%glint)
           enddo
           
           if(sunglint%s(prof)%windsp .eq. 0._jprb) then
              do ist = 0, nstreams
                temp_tl(ist) = temp_tl(ist) * transmission_aux%refl_norm(i)
              enddo
           endif
        else !sea-ice = land
           do ist = 0, nstreams
              temp_tl(ist) = coef%ss_solar_spectrum(chan) * cos(profiles(prof)%sunzenangle * deg2rad) * pi_r * &
                           (transmission_aux_tl%tausun_surf(ist,i) * reflectivity(i) + &
                           transmission_aux%tausun_surf(ist,i) * reflectivity_tl(i))
           enddo
        endif
        rad_tl%clear(i) = rad_tl%clear(i) + temp_tl(0)
        auxrad_stream_tl%cloudy(1:nstreams,i) = auxrad_stream_tl%cloudy(1:nstreams,i) + temp_tl(1:nstreams)
     endif
  enddo
end subroutine solar_contribution_tl

end subroutine rttov_integrate_tl
