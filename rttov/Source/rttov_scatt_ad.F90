!
Subroutine rttov_scatt_ad( & 
     & errorstatus,        &! out
     & nlevels,            &! in
     & chanprof,           &! in
     & frequencies,        &! in
     & profiles,           &! in
     & cld_profiles,       &! in
     & coef_rttov,         &! in
     & coef_scatt,         &! in
     & calcemiss,          &! in
     & emissivity_in,      &! in
     & profiles_ad,        &! inout
     & cld_profiles_ad,    &! inout
     & emissivity_in_ad,   &! inout
     & radiance,           &! inout
     & radiance_ad,        &! inout 
     & lnewcld,            &! in, optional
     & lusercfrac)          ! in, optional
     
  ! Description:
  ! AD of subroutine
  ! to compute microwave multi-channel radiances and brightness
  ! temperatures for many profiles per call in a cloudy and/or rainy sky.
  !
  ! NOTE ON ADJOINT / K: in RTTOV_SCATT the adjoint and "K" routines are 
  ! combined. Normal adjoint behaviour is that the number of profiles in 
  ! the adjoint variables ("nprofilesad") is the same as in the direct 
  ! variables ("nprofiles") but to get "K" behaviour, the adjoint increments 
  ! for each channel are kept separate and "nprofilesad" should be the same 
  ! as the number of channels, "nchannels"       
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
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and 
  !     precipitation. Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 21 pp.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative 
  !     transfer modeling in clouds and precipitation. Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans:comparison with SSM/I observations. 
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Smith, E. A., P. Bauer, F. S. Marzano, C. D. Kummerow, D. McKague, 
  !     A. Mugnai, G. Panegrossi, 2002:
  !     Intercomparison of microwave radiative transfer models for precipitating clouds.
  !     IEEE Trans. Geosci. Remote Sens. 40, 541-549.
  ! - Bauer, P., Moreau, E., Chevallier, F., O'Keeffe, U., 2006:
  !     Multiple-scattering microwave radiative transfer for data assimilation applications
  !     Quart. J. R. Meteorol. Soc. 132, 1259-1281
  ! - Geer, A.J., Bauer, P. and O'Dell, C.W., 2009: 
  !     A Revised Cloud Overlap Scheme for Fast Microwave Radiative Transfer in Rain and Cloud.
  !     Journal of Applied Meteorology and Climatology. 48, 2257-2270
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    07/2005      Polarimetry fixes   (U. O'Keeffe)
  !   1.6    11/2005      Add errorstatus to iniscatt arguments and use a temporary
  !                       radiance type for the calcpolarisation call (J Cameron)
  !   1.7    11/2007      RTTOV9 version      (A. Geer)
  !   1.8    07/2008      Clear sky speed-ups (A. Geer)
  !   1.9    10/2008      Revised cloud partitioning (A. Geer)
  !   1.10   11/2009      User may supply the average cloud fraction (A. Geer)
  !   1.11   11/2009      Adapted for RTTOV10 (A. Geer)
  !   1.12   04/2010      Tidied up after code cleaning (A. Geer)
  !   1.13   02/2010      Revised cloud partitioning is now the default (A. Geer)
  !   1.14   08/2010      Fix for polarimetric channels until they can be 
  !                       hanled properly by RTTOV_SCATT (W. Bell)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Parameters:
  
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  
  Use rttov_const, Only :   &
       & sensor_id_mw,      &
       & errorstatus_success, &
       & errorstatus_fatal,   &  
       & adk_adjoint,         &
       & adk_k,               &
       & ccthres,             &  
       & sensor_id_po

  Use rttov_types, Only :    &
       & rttov_coefs          ,&
       & rttov_scatt_coef     ,&
       & geometry_Type        ,&
       & profile_Type         ,&
       & profile_cloud_Type   ,&
       & profile_scatt_aux    ,&
       & transmission_Type    ,&
       & radiance_Type        ,&
       & rttov_chanprof       ,&
       & rttov_options   


  Use parkind1, Only : jpim, jprb, jplm
  
  Implicit None


!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in)  :: nlevels ! Number of levels
  Type(rttov_chanprof),Intent (in)  :: chanprof(:)                  ! Indices
  Type (profile_Type), Intent (in)  :: profiles(:)                  ! Atmospheric profiles
  Type (profile_Type), Intent (inout) :: profiles_ad(:)    
  Integer (Kind=jpim), Intent (in)  :: frequencies (size(chanprof)) ! Frequency indices
  Integer (Kind=jpim), Intent (out) :: errorstatus                  ! Error return flag

  Logical (Kind=jplm), Intent (in)    :: calcemiss        (size(chanprof))          ! Switch for emmissivity calculation
  Real    (Kind=jprb), Intent (in)    :: emissivity_in    (size(chanprof))          ! Surface emissivity 
  Real    (Kind=jprb), Intent (inout) :: emissivity_in_ad (size(chanprof))          ! Surface emissivity 
  

  Type (rttov_coefs),         Intent (in)    :: coef_rttov                     ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                     ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (size(profiles))    ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (inout) :: cld_profiles_ad (size(profiles_ad))   
  Type (radiance_Type),       Intent (inout) :: radiance         ! Radiances
  Type (radiance_Type),       Intent (inout) :: radiance_ad          

  Logical (Kind=jplm), optional, Intent(in) :: lnewcld ! T = revised cloud/rain treatment; F = old treatment; Default T
  Logical (Kind=jplm), optional, Intent(in) :: lusercfrac ! T = take av. cloud fraction from that supplied; Default F

!INTF_END

#include "rttov_direct.h"
#include "rttov_iniscatt.h"
#include "rttov_eddington.h"
#include "rttov_ad.h"
#include "rttov_k.h"
#include "rttov_iniscatt_ad.h"
#include "rttov_eddington_ad.h"
#include "rttov_errorreport.h"
   
  Integer (Kind=jpim), target :: sa__mclayer    (size(chanprof))
  Integer (Kind=jpim), target :: sa_ad__mclayer (size(chanprof))
      
  Real (Kind=jprb), target :: t__tau_total     (size(chanprof))
  Real (Kind=jprb), target :: t__tau_levels    (nlevels,size(chanprof))
  Real (Kind=jprb), target :: t_ad__tau_total  (size(chanprof))
  Real (Kind=jprb), target :: t_ad__tau_levels (nlevels,size(chanprof))
  
  Real (Kind=jprb), target :: sa__cfrac   (size(profiles))  
  Real (Kind=jprb), target :: sa__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa__ref_cld (size(chanprof))
  
  Real (Kind=jprb), target :: sa__tbd (size(profiles),nlevels+1)
  
  Real (Kind=jprb), target :: sa__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa__h      (size(chanprof),nlevels)
  
  Real (Kind=jprb), target :: sa__b0     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__b1     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__bn     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__dz     (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__clw    (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__ciw    (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__totalice (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__rain   (size(profiles),nlevels)
  Real (Kind=jprb), target :: sa__sp     (size(profiles),nlevels)

  Real (Kind=jprb), target :: sa_ad__cfrac   (size(profiles_ad))
  Real (Kind=jprb), target :: sa_ad__ems_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_ad__ref_bnd (size(chanprof))
  Real (Kind=jprb), target :: sa_ad__ems_cld (size(chanprof))
  Real (Kind=jprb), target :: sa_ad__ref_cld (size(chanprof))
  
  Real (Kind=jprb), target :: sa_ad__tbd (size(profiles_ad),nlevels+1)
  
  Real (Kind=jprb), target :: sa_ad__delta  (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__tau    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__ext    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__ssa    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__asm    (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__lambda (size(chanprof),nlevels)
  Real (Kind=jprb), target :: sa_ad__h      (size(chanprof),nlevels)
  
  Real (Kind=jprb), target :: sa_ad__b0     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__b1     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__bn     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__dz     (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__clw    (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__ciw    (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__totalice (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__rain   (size(profiles_ad),nlevels)
  Real (Kind=jprb), target :: sa_ad__sp     (size(profiles_ad),nlevels)

!* Local variables:
  Integer (Kind=jpim) :: nprofiles, nprofilesad, nchannels
  Logical (Kind=jplm) :: usenewcld, usercfrac
  Integer (Kind=jpim) :: iprof, ichan
  Integer (Kind=jpim) :: iprofad, adk  
  Integer (Kind=jpim) :: pol_id, chan

  Real (Kind=jprb), Dimension (size(chanprof)) :: cld_bt, cld_bt_ad
       
  Type (transmission_Type)   :: transmission, transmission_ad
  Type (geometry_Type)       :: angles (size(profiles))
  Type (profile_scatt_aux)   :: scatt_aux, scatt_aux_ad   

  Type(rttov_options)  :: opts
  Real    (Kind=jprb)  :: emissivity_out (size(chanprof))
  Real    (Kind=jprb)  :: emissivity_out_ad (size(chanprof))

  Character (len=80) :: errMessage
  Character (len=15) :: NameOfRoutine = 'rttov_scatt_ad '

  REAL(KIND=JPRB) :: ZHOOK_HANDLE
  
  !- End of header --------------------------------------------------------
  
  if (lhook) call dr_hook('RTTOV_SCATT_AD',0_jpim,zhook_handle)        

  nprofiles   = size(profiles)
  nprofilesad = size(profiles_ad)
  nchannels   = size(chanprof)
  
  usenewcld = .true.
  if(present(lnewcld)) then
    usenewcld = lnewcld
  endif

  usercfrac      = .false.
  if(present(lusercfrac)) then
    if(lusercfrac) usercfrac = .true.
  endif
  
  errorstatus = errorstatus_success

  if (nprofilesad == nprofiles) then 
     adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
     adk = adk_k         ! K mode
  else
     Write( errMessage, '( "incorrect number of profiles in adjoint/K")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  endif 

  transmission % tau_total  => t__tau_total
  transmission % tau_levels => t__tau_levels

  transmission_ad % tau_total  => t_ad__tau_total
  transmission_ad % tau_levels => t_ad__tau_levels

  scatt_aux % cfrac    => sa__cfrac
  scatt_aux % ems_bnd  => sa__ems_bnd
  scatt_aux % ref_bnd  => sa__ref_bnd
  scatt_aux % ems_cld  => sa__ems_cld
  scatt_aux % ref_cld  => sa__ref_cld
  scatt_aux % tbd      => sa__tbd
  scatt_aux % mclayer  => sa__mclayer
  scatt_aux % delta    => sa__delta
  scatt_aux % tau      => sa__tau
  scatt_aux % ext      => sa__ext
  scatt_aux % ssa      => sa__ssa
  scatt_aux % asm      => sa__asm
  scatt_aux % lambda   => sa__lambda
  scatt_aux % h        => sa__h
  scatt_aux % b0       => sa__b0
  scatt_aux % b1       => sa__b1
  scatt_aux % bn       => sa__bn
  scatt_aux % dz       => sa__dz
  scatt_aux % clw      => sa__clw
  scatt_aux % ciw      => sa__ciw
  scatt_aux % totalice => sa__totalice
  scatt_aux % rain     => sa__rain
  scatt_aux % sp       => sa__sp

  scatt_aux_ad % cfrac    => sa_ad__cfrac
  scatt_aux_ad % ems_bnd  => sa_ad__ems_bnd
  scatt_aux_ad % ref_bnd  => sa_ad__ref_bnd
  scatt_aux_ad % ems_cld  => sa_ad__ems_cld
  scatt_aux_ad % ref_cld  => sa_ad__ref_cld
  scatt_aux_ad % tbd      => sa_ad__tbd
  scatt_aux_ad % mclayer  => sa_ad__mclayer
  scatt_aux_ad % delta    => sa_ad__delta
  scatt_aux_ad % tau      => sa_ad__tau
  scatt_aux_ad % ext      => sa_ad__ext
  scatt_aux_ad % ssa      => sa_ad__ssa
  scatt_aux_ad % asm      => sa_ad__asm
  scatt_aux_ad % lambda   => sa_ad__lambda
  scatt_aux_ad % h        => sa_ad__h
  scatt_aux_ad % b0       => sa_ad__b0
  scatt_aux_ad % b1       => sa_ad__b1
  scatt_aux_ad % bn       => sa_ad__bn
  scatt_aux_ad % dz       => sa_ad__dz
  scatt_aux_ad % clw      => sa_ad__clw
  scatt_aux_ad % ciw      => sa_ad__ciw
  scatt_aux_ad % totalice => sa_ad__totalice
  scatt_aux_ad % rain     => sa_ad__rain
  scatt_aux_ad % sp       => sa_ad__sp

  ! Check inputs
  ! ------------
  Do iprof = 1, nprofiles
    If (  profiles(iprof) % s2m % p /= cld_profiles(iprof) % ph(nlevels+1)  ) Then
      errorstatus = errorstatus_fatal
      Write( errMessage, '( "Surface pressure and lowest half level should be identical")' )
      Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
      IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT',1_jpim,ZHOOK_HANDLE)
      Return
    End If
  End Do

!*         1.   Gas absorption

  ! Profiles will be interpolated from model/RTTOV-SCATT levels to 
  ! RTTOV coefficient levels within RTTOV itself.
  opts%addinterp  = .true.  
  opts%addclouds  = .false.
  opts%addsolar   = .false.
  opts%addaerosl  = .false.
  opts%addpc      = .false.
  opts%addradrec  = .false.
  opts%switchrad  = .true.   ! Required, in AD only, and for some historical reason
 
  Call rttov_direct(          &   
        & errorstatus,        &! out
        & chanprof,           &! in
        & opts,               &! in
        & profiles,           &! in
        & coef_rttov,         &! in
        & calcemiss,          &! in
        & emissivity_in,      &! inout
        & emissivity_out,     &
        & transmission,       &! inout
        & radiance     )       ! inout 

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_direct")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If

  scatt_aux % ems_cld (:) = emissivity_in (:)      
  scatt_aux % ref_cld (:) = 1.0_JPRB - emissivity_in (:)      
  
!*  2.   Initialisations for Eddington
  Call rttov_iniscatt(       &
        & errorstatus,       &! out
        & nlevels,           &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & chanprof,          &! in
        & frequencies,       &! in
        & profiles,          &! in  
        & cld_profiles,      &! in
        & coef_rttov%coef,   &! in  
        & coef_scatt,        &! in  
        & transmission,      &! in
        & calcemiss,         &! in
        & usenewcld,         &! in
        & usercfrac,         &! in
        & angles,            &! out
        & scatt_aux   )       ! inout   

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If

!* 3.   Eddington (in temperature space)
   Call rttov_eddington(     &
        & nlevels,           &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & chanprof%prof,     &! in
        & angles,            &! in
        & profiles,          &! in  
        & scatt_aux,         &! in
        & cld_bt)             ! out   

!*  4.   Combine clear and cloudy parts
  Do ichan = 1, nchannels

    iprof = chanprof(ichan) % prof
    if (scatt_aux % cfrac (iprof) > ccthres ) then 

      radiance % bt (ichan) = cld_bt (ichan) * scatt_aux % cfrac (iprof) & 
        & + radiance % bt_clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))

    else

      radiance % bt (ichan) = radiance % bt_clear (ichan)

    endif
  End Do
    
!* ADJOINT PART
  
  scatt_aux_ad % cfrac   (:)   = 0.0_JPRB
  scatt_aux_ad % ems_bnd (:)   = 0.0_JPRB
  scatt_aux_ad % ref_bnd (:)   = 0.0_JPRB
  scatt_aux_ad % ems_cld (:)   = 0.0_JPRB
  scatt_aux_ad % ref_cld (:)   = 0.0_JPRB
  scatt_aux_ad % tbd     (:,:) = 0.0_JPRB
  scatt_aux_ad % delta   (:,:) = 0.0_JPRB
  scatt_aux_ad % tau     (:,:) = 0.0_JPRB
  scatt_aux_ad % ext     (:,:) = 0.0_JPRB
  scatt_aux_ad % ssa     (:,:) = 0.0_JPRB
  scatt_aux_ad % asm     (:,:) = 0.0_JPRB
  scatt_aux_ad % lambda  (:,:) = 0.0_JPRB
  scatt_aux_ad % h       (:,:) = 0.0_JPRB
  scatt_aux_ad % b0      (:,:) = 0.0_JPRB
  scatt_aux_ad % b1      (:,:) = 0.0_JPRB
  scatt_aux_ad % bn      (:,:) = 0.0_JPRB
  scatt_aux_ad % dz      (:,:) = 0.0_JPRB
  scatt_aux_ad % clw     (:,:) = 0.0_JPRB
  scatt_aux_ad % ciw     (:,:) = 0.0_JPRB
  scatt_aux_ad % totalice (:,:) = 0.0_JPRB
  scatt_aux_ad % rain    (:,:) = 0.0_JPRB
  scatt_aux_ad % sp      (:,:) = 0.0_JPRB

  transmission_ad % tau_total       (:)   = 0.0_JPRB
  transmission_ad % tau_levels      (:,:) = 0.0_JPRB 

  emissivity_out_ad(:)=0.0_JPRB 

  cld_bt_ad (:) = 0.0_JPRB

!* 4.   Combine clear and cloudy parts  
  Do ichan = 1, nchannels
    iprof = chanprof (ichan) % prof

    if (adk == adk_adjoint) then
      iprofad = iprof
    else if (adk == adk_k) then
      iprofad = ichan
    endif

! Use only clear sky part of RT calculation for polarimetric channels

    chan=chanprof (ichan) % chan

    If(coef_rttov % coef % id_sensor == sensor_id_po) Then

      pol_id = coef_rttov % coef % fastem_polar(chan) + 1_jpim

    Else

      pol_id = 0_jpim

    Endif

    If (pol_id >= 6_jpim ) Then
      radiance_ad % bt_clear (ichan) = radiance_ad % bt_clear (ichan) + &
       & radiance_ad % bt (ichan)

      radiance_ad % bt (ichan) = 0.0_JPRB
    End If

! Usual case

    if (scatt_aux % cfrac (iprof) > ccthres ) then 

      cld_bt_ad (ichan) = cld_bt_ad (ichan) + &
       & scatt_aux % cfrac (iprof) * radiance_ad % bt (ichan) 
      
      scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad) + &
       & (cld_bt (ichan) - radiance % bt_clear (ichan)) * radiance_ad % bt (ichan) 

      radiance_ad  % bt_clear (ichan) = radiance_ad % bt_clear (ichan) + &
       & (1.0_JPRB - scatt_aux % cfrac (iprof)) * radiance_ad % bt (ichan)

      radiance_ad % bt (ichan) = 0.0_JPRB

    else

      radiance_ad % bt_clear (ichan) = radiance_ad % bt_clear (ichan) + &
       & radiance_ad % bt (ichan)

      radiance_ad % bt (ichan) = 0.0_JPRB

    endif

  End Do

!* 3.   Eddington (in temperature space)
  Call rttov_eddington_ad(   &
        & nlevels,           &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & nprofilesad,       &! in
        & chanprof % prof,   &! in
        & angles,            &! in
        & profiles,          &! in  
        & profiles_ad,       &! inout  
        & scatt_aux,         &! in
        & scatt_aux_ad,      &! inout
        & cld_bt,            &! out  
        & cld_bt_ad)          ! inout   

!*  2.   Initialisations for Eddington  
  Call rttov_iniscatt_ad(    &
        & errorstatus,       &! out
        & nlevels,           &! in
        & nchannels,         &! in
        & nprofiles,         &! in
        & nprofilesad,       &! in
        & chanprof,          &! in
        & frequencies,       &! in
        & profiles,          &! in 
        & profiles_ad,       &! inout
        & cld_profiles,      &! in
        & cld_profiles_ad,   &! inout
        & coef_rttov%coef,   &! in 
        & coef_scatt,        &! in  
        & transmission,      &! in
        & transmission_ad,   &! inout
        & calcemiss,         &! in
        & usenewcld,         &! in
        & usercfrac,         &! in
        & angles,            &! out
        & scatt_aux,         &! inout
        & scatt_aux_ad)       ! inout   

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_iniscatt_ad")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If
  
  emissivity_in_ad (:) = emissivity_in_ad (:) - scatt_aux_ad % ref_cld (:)
  scatt_aux_ad % ref_cld (:) = 0.0_JPRB      
  emissivity_in_ad (:) = emissivity_in_ad (:) + scatt_aux_ad % ems_cld (:)  
  scatt_aux_ad % ems_cld (:) = 0.0_JPRB   
  
!*         1.   Gas absorption

  if (adk == adk_adjoint) then  
    Call rttov_ad(        &
       & errorstatus,     &! out
       & chanprof,        &! in
       & opts,            &
       & profiles,        &! in
       & profiles_ad,     &! inout
       & coef_rttov,      &! in
       & calcemiss,       &! in
       & emissivity_in,   &! inout
       & emissivity_in_ad,&! inout
       & emissivity_out,  &
       & emissivity_out_ad,&
       & transmission,    &! inout
       & transmission_ad, &! inout
       & radiance,        &! inout
       & radiance_ad     ) ! inout 
  else if (adk == adk_k) then 
    Call rttov_k(        &
       & errorstatus,     &! out
       & chanprof,        &! in
       & opts,            &! in
       & profiles,        &! in
       & profiles_ad,     &! inout
       & coef_rttov,      &! in
       & calcemiss,       &! in
       & emissivity_in,   &! inout
       & emissivity_in_ad,&! inout
       & emissivity_out,  &
       & emissivity_out_ad,&
       & transmission,    &! inout
       & transmission_ad, &! inout
       & radiance,        &! inout
       & radiance_ad     ) ! inout 
  endif

  If ( errorstatus == errorstatus_fatal ) Then
     Write( errMessage, '( "error in rttov_ad")' )
     Call rttov_errorreport (errorstatus_fatal, errMessage, NameOfRoutine)
     IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_AD',1_jpim,ZHOOK_HANDLE)
     Return
  End If

!*  4b.   Combine clear and cloudy parts - again - because rttov_ad/k has
!*        overwritten the bt array with clear-sky values.
  Do ichan = 1, nchannels
    iprof = chanprof (ichan) % prof
    if (scatt_aux % cfrac (iprof) > ccthres ) then 

      radiance % bt (ichan) = cld_bt (ichan) * scatt_aux % cfrac (iprof) &
       & + radiance % bt_clear (ichan) * (1.0_JPRB - scatt_aux % cfrac (iprof))

    else

      radiance % bt (ichan) = radiance % bt_clear (ichan)

    endif

! Use only clear sky part of RT calculation for polarimetric channels

    chan=chanprof(ichan) % chan

    If(coef_rttov % coef % id_sensor == sensor_id_po) Then

      pol_id = coef_rttov % coef % fastem_polar(chan) + 1_jpim

    Else

      pol_id = 0_jpim

    Endif

    If (pol_id >= 6_jpim ) radiance % bt (ichan) = radiance % bt_clear (ichan)

  End Do

  if (lhook) call dr_hook('RTTOV_SCATT_AD',1_jpim,zhook_handle)

End Subroutine rttov_scatt_ad
