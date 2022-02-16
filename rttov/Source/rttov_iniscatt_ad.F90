!
Subroutine rttov_iniscatt_ad (&
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
      & coef_rttov,        &! in
      & coef_scatt,        &! in
      & transmission,      &! in
      & transmission_ad,   &! inout
      & calcemiss,         &! in
      & usenewcld,         &! in
      & usercfrac,         &! in
      & angles,            &! out
      & scatt_aux,         &! inout
      & scatt_aux_ad)       ! inout 

  !
  ! Description:
  ! AD of routine to
  ! Calculate some variables related to the input precipitation profile
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
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !   1.0    09/2002      Initial version     (F. Chevallier)
  !   1.1    05/2003      RTTOV7.3 compatible (F. Chevallier)
  !   1.2    03/2004      Added polarimetry   (R. Saunders)
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keefe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    10/2005      Fixes for rttov8 indexing   (U. O'Keeffe)
  !   1.6    11/2005      Limit lines to 132 characters
  !                       add errorstatus to arguments
  !                       change stop to return (J Cameron)
  !   1.7    09/2006      Add if loop to stop use of iccmax index
  !                       if = 0 (A. Doherty)
  !   1.8    11/2007      RTTOV9 / cleanup (A. Geer)   
  !   1.9    03/2008      Revised cloud partitioning (A. Geer)  
  !   1.10   03/2009      Safety check on cloud fraction (A. Geer) 
  !   1.11   11/2009      RTTOV transmittances / optical depths come on levels (A Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:
  
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
    
  Use rttov_types, Only :    &
       & rttov_coef           ,&
       & rttov_scatt_coef     ,&
       & transmission_type    ,&
       & transmission_type_aux,&
       & geometry_Type        ,&
       & profile_scatt_aux    ,&
       & profile_Type         ,&
       & profile_cloud_Type   ,&
       & rttov_chanprof      ,&
       & rttov_options  

  Use rttov_const, Only:    &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & gravity,             &
       & pressure_top,        &
       & rgp,                 &
       & rm,                  &
       & rho_rain,            &
       & rho_snow,            &
       & ccthres,             & 
       & adk_adjoint,         &
       & adk_k

  Use parkind1, Only : jpim, jprb, jplm
  
  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles  ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nprofilesad  ! Number of profiles in adjoint
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent(out) :: errorstatus             ! Error return code
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Logical (Kind=jplm), Intent (in) :: calcemiss   (nchannels) ! Emissivity flags
  Logical (Kind=jplm), Intent (in) :: usenewcld               ! New or old cloud partition
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction

  Type (profile_Type),        Intent (in)    :: profiles        (nprofiles)   ! Atmospheric profiles
  Type (profile_Type),        Intent (inout) :: profiles_ad     (nprofilesad) ! Atmospheric profiles
  Type (rttov_coef),          Intent (in)    :: coef_rttov                    ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                    ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (nprofiles)   ! Cloud profiles
  Type (profile_cloud_Type),  Intent (inout) :: cld_profiles_ad (nprofilesad) ! Cloud profiles
  Type (transmission_Type),   Intent (in)    :: transmission                  ! Transmittances and optical depths
  Type (transmission_Type),   Intent (inout) :: transmission_ad               ! Transmittances and optical depths
  Type (geometry_Type),       Intent (out)   :: angles          (nprofiles)   ! Zenith angles
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux                     ! Auxiliary profile variables
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux_ad                  ! Auxiliary profile variables

!INTF_END
 
!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, ichan, iccmax(nprofiles)
  Integer (Kind=jpim) :: iprofad, adk
  Real    (Kind=jprb) :: p1, p2, pm, p1_ad, p2_ad, pm_ad, dp2dz, de2mr, zccmax

  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf   ! Pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh  ! Half-level pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: od
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop  ! Thickness of top half of level [km]
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot  ! Thickness of bottom half of level [km]
  Real    (Kind=jprb), Dimension (2:nlevels)           :: dzr    ! Thickness of RTTOV level [km]
  Real    (Kind=jprb), Dimension (nchannels,nlevels+1) :: od_rttov ! Single layer optical depths on RTTOV levels

  Real    (Kind=jprb), Dimension (nchannels)           :: zod_up_cld   ! Optical depth from top of the atmosphere 
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels)   :: presf_ad     ! Pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels+1) :: presfh_ad    ! Half-level pressure levels [hPa]
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: od_ad
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels)   :: dztop_ad 
  Real    (Kind=jprb), Dimension (nprofilesad,nlevels)   :: dzbot_ad  
  Real    (Kind=jprb), Dimension (2:nlevels)           :: dzr_ad    
  Real    (Kind=jprb), Dimension (nchannels,nlevels+1) :: od_rttov_ad 
  Real    (Kind=jprb), Dimension (nchannels)           :: zod_up_cld_ad    ! Optical depth from top of the atmosphere 
   
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_0
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_1, ssa_1, asm_1
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_2, ssa_2, asm_2
  Real    (Kind=jprb), Dimension (nchannels,nlevels)   :: ext_3, ssa_3, asm_3
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: clw_scale, ciw_scale, &
   & totalice_scale, rain_scale, sp_scale
  Real    (Kind=jprb), Dimension (nprofiles,nlevels)   :: clw_precf, ciw_precf, &
   & rain_precf, sp_precf, totalice_precf

  Real    (Kind=jprb) :: hydro_weights(nprofiles,nlevels), hydro_column(nprofiles)   ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: cfrac(nprofiles)
  Real    (Kind=jprb) :: hydro_weights_ad(nlevels), hydro_column_ad        

  Type (transmission_Type_aux) :: transmissioncld     ! Clear+cloud transmittances with cloud
  Type (transmission_Type_aux) :: transmissioncld_ad  ! Clear+cloud transmittances with cloud
  Type(rttov_options)  :: opts

  Character (len=80) :: errMessage
  Character (len=18) :: NameOfRoutine = 'rttov_iniscatt_ad '

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "rttov_mieproc.h"
#include "rttov_iniedd.h"
#include "rttov_calcemis_mw.h"
#include "rttov_mieproc_ad.h"
#include "rttov_iniedd_ad.h"
#include "rttov_calcemis_mw_ad.h"
#include "rttov_setgeometry.h"
#include "rttov_errorreport.h"
#include "rttov_calcemis_mw_k.h"

  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INISCATT_AD',0_jpim,zhook_handle)

  if (nprofilesad == nprofiles) then 
    adk = adk_adjoint   ! Adjoint mode
  else if (nprofilesad == nchannels) then
    adk = adk_k         ! K mode
  endif 

  ! NB this only truly needed if using FASTEM 3 
  allocate (transmissioncld % tau_surf (0:0,nchannels))
  allocate (transmissioncld_ad % tau_surf (0:0,nchannels))

  errorstatus = errorstatus_success


  de2mr =  1.0E+05_JPRB * rm / rgp
  dp2dz = -1.0E-03_JPRB * rgp / gravity / rm 

  scatt_aux % ext (:,:) = 0.0_JPRB
  scatt_aux % ssa (:,:) = 0.0_JPRB
  scatt_aux % asm (:,:) = 0.0_JPRB

!* Security on user-defined pressures
  Do iprof = 1, nprofiles
     Do ilayer = 1, nlevels
        If (profiles (iprof) % p (ilayer) >= pressure_top) Then
            presf (iprof,ilayer) = profiles (iprof) % p (ilayer)
       else
            presf (iprof,ilayer) = pressure_top
       Endif
     Enddo
     Do ilayer = 1, nlevels + 1
        If (cld_profiles (iprof) % ph (ilayer) >= pressure_top ) Then
            presfh (iprof,ilayer) = cld_profiles (iprof) % ph (ilayer)
        else
            presfh (iprof,ilayer) = pressure_top
        Endif
     Enddo
  Enddo

!* Set up geometric variables
  Call rttov_setgeometry ( &
    & opts,       & ! in
    & profiles,        & ! in
    & coef=coef_rttov, & ! in
    & angles=angles)     ! out

!* Compute temperature at layer boundaries (K)
  Do iprof = 1, nprofiles
     scatt_aux % tbd (iprof,nlevels+1) = profiles (iprof) % s2m % t
     scatt_aux % tbd (iprof,1)         = profiles (iprof) % t(1)
  Enddo

  Do ilayer = 1, nlevels-1
     Do iprof = 1, nprofiles     
        p1 = presf  (iprof,ilayer+1)
        p2 = presf  (iprof,ilayer  )
        pm = presfh (iprof,ilayer+1)

        scatt_aux % tbd (iprof,ilayer+1) =  profiles (iprof) % t (ilayer+1) &
                                       & + (profiles (iprof) % t (ilayer)   & 
                                       & -  profiles (iprof) % t (ilayer+1)) &
                                       & / log(p2/p1) * log(pm/p1)        
     Enddo
  Enddo

!* Initialise cloud and rain properties of the cloudy/rainy column
  clw_scale   (:,:) = 0.0_JPRB
  ciw_scale   (:,:) = 0.0_JPRB
  totalice_scale   (:,:) = 0.0_JPRB
  rain_scale  (:,:) = 0.0_JPRB
  sp_scale    (:,:) = 0.0_JPRB

  if(usenewcld) then

    ! Weighted average partitioning between "cloudy" and "clear" columns
    Do ilayer=1,nlevels
      Do iprof = 1, nprofiles  
        clw_scale  (iprof,ilayer) = cld_profiles (iprof) % clw  (ilayer) 
        rain_scale (iprof,ilayer) = cld_profiles (iprof) % rain (ilayer) 
        if ( cld_profiles (iprof) % use_totalice ) then
          totalice_scale (iprof,ilayer) = cld_profiles (iprof) % totalice (ilayer)
        else
          ciw_scale  (iprof,ilayer) = cld_profiles (iprof) % ciw  (ilayer)
          sp_scale   (iprof,ilayer) = cld_profiles (iprof) % sp   (ilayer)
        endif 
      Enddo
    Enddo

  else

    ! Maximum cloud fraction partitioning between "cloudy" and "clear" columns
    iccmax (:) = 0

    Do iprof = 1, nprofiles  
      zccmax = 0.0_JPRB

      Do ilayer = 1, nlevels
        if (cld_profiles (iprof) % cc (ilayer) > zccmax) then
          zccmax = cld_profiles (iprof) % cc (ilayer)
          iccmax(iprof) = ilayer
        end if 
      end do 
      scatt_aux % cfrac (iprof) = zccmax
    Enddo

    do ilayer=1,nlevels
      do iprof = 1, nprofiles  
        If (scatt_aux % cfrac (iprof) > ccthres) Then
          clw_scale (iprof,ilayer) = cld_profiles (iprof) % clw  (ilayer) / scatt_aux % cfrac (iprof)
          rain_scale (iprof,ilayer) = cld_profiles (iprof) % rain (ilayer) / scatt_aux % cfrac (iprof)
          if ( cld_profiles (iprof) % use_totalice ) then
            totalice_scale (iprof,ilayer) = cld_profiles (iprof) % totalice (ilayer) / scatt_aux % cfrac (iprof)
          else
            ciw_scale (iprof,ilayer) = cld_profiles (iprof) % ciw  (ilayer) / scatt_aux % cfrac (iprof)
            sp_scale (iprof,ilayer) = cld_profiles (iprof) % sp   (ilayer) / scatt_aux % cfrac (iprof)
          endif
        Endif
      enddo
    enddo

  endif
 
!* Nadir heights (km)
  Do ilayer = nlevels, 1, -1
     Do iprof = 1, nprofiles
        p1 = presfh (iprof,ilayer+1)
        p2 = presfh (iprof,ilayer  )
        pm = presf  (iprof,ilayer  )

        If (p1 <= p2) then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "iniscatt : problem with user-defined pressure layering")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
           if (lhook) call dr_hook('RTTOV_INISCATT',1_jpim,zhook_handle)
           Return
        End If

        scatt_aux % dz (iprof,ilayer) = dp2dz * Log(p2/p1) * profiles (iprof) % t (ilayer)

        dzbot (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles (iprof) % t (ilayer)
        dztop (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles (iprof) % t (ilayer)

     Enddo
  Enddo

!* Get single-layer optical depths (at nadir and in hPa-1) and put onto model half levels
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof

! Top RTTOV level to space    
     od_rttov (ichan,1)         = -1.0_jprb * log( transmission % tau_levels (1,ichan) )     
     Do ilayer = 2, nlevels 
        od_rttov (ichan,ilayer) = log( transmission % tau_levels (ilayer-1,ichan) ) &
                              & - log( transmission % tau_levels (ilayer,ichan) )
     Enddo
! Surface to bottom RTTOV (full pressure) level
     od_rttov (ichan,nlevels+1) = log( transmission % tau_levels (nlevels,ichan) ) &
                              & - log( transmission % tau_total (ichan) )

! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzr(ilayer) = dzbot(iprof,ilayer-1) + dztop(iprof,ilayer)
     Enddo

! Re-allocate optical depths between half pressure levels        
     od (ichan,1)         = od_rttov(ichan,1) &
                        & + od_rttov(ichan,2) * dzbot(iprof,1) / dzr(2)     
     Do ilayer = 2, nlevels - 1  
        od (ichan,ilayer) = od_rttov(ichan,ilayer)   * dztop(iprof,ilayer) / dzr(ilayer) &
                        & + od_rttov(ichan,ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) 
     Enddo
     od (ichan,nlevels)   = od_rttov(ichan,nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + od_rttov(ichan,nlevels+1) 
      
  Enddo

!* Change units
  Do ilayer = 1, nlevels

!* Optical depths in km-1 and at nadir
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
     
        scatt_aux % ext (ichan,ilayer) = od (ichan,ilayer) &
          & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen 
     
        ext_0 (ichan,ilayer) = scatt_aux % ext (ichan,ilayer)     

        if (scatt_aux % ext (ichan,ilayer) < 1.0E-10_JPRB) scatt_aux % ext (ichan,ilayer) = 1.0E-10_JPRB
     Enddo

!* Condensate from g/g to g/m^3
     Do iprof = 1, nprofiles

        scatt_aux % clw (iprof,ilayer) = clw_scale (iprof,ilayer) &
          & * presf (iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer)
        scatt_aux % ciw (iprof,ilayer) = ciw_scale (iprof,ilayer) &
          & * presf (iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer)
        scatt_aux % totalice (iprof,ilayer) = totalice_scale (iprof,ilayer) &
          & * presf (iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer)
    
!* Rates from kg/m^2/s to g/m^3
        rain_scale (iprof,ilayer) =  rain_scale (iprof,ilayer) / rho_rain
        sp_scale   (iprof,ilayer) =  sp_scale   (iprof,ilayer) / rho_snow

        rain_scale (iprof,ilayer) =  rain_scale (iprof,ilayer) * 3600.0_JPRB 
        sp_scale   (iprof,ilayer) =  sp_scale   (iprof,ilayer) * 3600.0_JPRB

        if (rain_scale (iprof,ilayer) > 0.0_JPRB) scatt_aux % rain (iprof,ilayer) = &
          & (rain_scale (iprof,ilayer) * coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
        if (sp_scale   (iprof,ilayer) > 0.0_JPRB) scatt_aux % sp   (iprof,ilayer) = &
          & (sp_scale   (iprof,ilayer) * coef_scatt % conv_sp   (1))**(coef_scatt % conv_sp   (2))
      Enddo
 Enddo

!* Store clear-sky absorption/scattering parameters
  ext_1 (:,:) = scatt_aux % ext (:,:)
  ssa_1 (:,:) = scatt_aux % ssa (:,:)
  asm_1 (:,:) = scatt_aux % asm (:,:)

!* Store hydrometeor amounts in g/m^3
  clw_precf (:,:) = scatt_aux % clw (:,:)
  ciw_precf (:,:) = scatt_aux % ciw (:,:)
  totalice_precf (:,:) = scatt_aux % totalice (:,:)
  rain_precf (:,:) = scatt_aux % rain (:,:)
  sp_precf (:,:) = scatt_aux % sp (:,:)

  if (usenewcld) then

    ! Calculate a hydrometeor-weighted average cloudy sky fraction 
    scatt_aux % cfrac (:)   = 0.0_JPRB

    Do iprof = 1, nprofiles

      if( usercfrac ) then

        !* User-supplied cloud fraction
        scatt_aux % cfrac (iprof) = cld_profiles (iprof) % cfrac

      else

        !* Partial column of hydrometeors in g/m^2
        hydro_weights(iprof,:) = &
          & ( scatt_aux % rain (iprof,:) + scatt_aux % sp (iprof,:) &
          &   + scatt_aux % ciw (iprof,:)  + scatt_aux % clw (iprof,:)  &
          &   + scatt_aux % totalice (iprof,:) ) &
          & * scatt_aux % dz (iprof,:)
        hydro_column(iprof) = sum( hydro_weights(iprof,:) )

        !* Weighted mean cloud fraction
        if (hydro_column(iprof) > 1e-10_JPRB) then
          scatt_aux % cfrac (iprof) = &
            & sum(hydro_weights(iprof,:) * cld_profiles (iprof) % cc (:)) / hydro_column(iprof)
          cfrac(iprof) = scatt_aux % cfrac (iprof) ! Store for adjoint use
          if ( cfrac(iprof) < 0.0_JPRB ) scatt_aux % cfrac (iprof) = 0.0_JPRB
          if ( cfrac(iprof) > 1.0_JPRB ) scatt_aux % cfrac (iprof) = 1.0_JPRB
        else
          scatt_aux % cfrac (iprof) = 0.0_JPRB
        endif

      endif

      !* Partition all cloud and rain into the cloudy column
      If (scatt_aux % cfrac (iprof) > ccthres) Then
        scatt_aux % clw  (iprof,:) = scatt_aux % clw  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % ciw  (iprof,:) = scatt_aux % ciw  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % totalice (iprof,:) = scatt_aux % totalice (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % rain (iprof,:) = scatt_aux % rain (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % sp   (iprof,:) = scatt_aux % sp   (iprof,:) / scatt_aux % cfrac (iprof)
      else
        scatt_aux % clw  (iprof,:) = 0.0_JPRB
        scatt_aux % ciw  (iprof,:) = 0.0_JPRB
        scatt_aux % totalice (iprof,:) = 0.0_JPRB
        scatt_aux % rain (iprof,:) = 0.0_JPRB
        scatt_aux % sp   (iprof,:) = 0.0_JPRB
      Endif
    Enddo
  endif

!* Cloud/rain absorption/scattering parameters
  Call rttov_mieproc (      &
       & nlevels,           &! in
       & nchannels,         &! in
       & nprofiles,         &! in
       & frequencies,       &! in
       & chanprof%prof,     &! in
       & profiles,          &! in
       & coef_scatt,        &! in
       & scatt_aux)          ! inout 
       
!* Store clear+cloud+rain absorption/scattering parameters
  ext_2 (:,:) = scatt_aux % ext (:,:)
  ssa_2 (:,:) = scatt_aux % ssa (:,:)
  asm_2 (:,:) = scatt_aux % asm (:,:)

!* Scattering parameters for Eddington RT
  Call rttov_iniedd(        &
       & nlevels,           &! in
       & nchannels ,        &! in
       & nprofiles ,        &! in
       & chanprof%prof,     &! in
       & angles    ,        &! in
       & scatt_aux)          ! inout 

!* Store delta-scaled clear+cloud+rain absorption/scattering parameters
  ext_3 (:,:) = scatt_aux % ext (:,:)
  ssa_3 (:,:) = scatt_aux % ssa (:,:)
  asm_3 (:,:) = scatt_aux % asm (:,:)

!* Surface emissivities
  zod_up_cld (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     Do ilayer = 1, nlevels     
        zod_up_cld (ichan) = zod_up_cld (ichan) + scatt_aux % ext (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) 
     Enddo
     if (zod_up_cld (ichan) >= 30.0_JPRB) zod_up_cld (ichan) = 30.0_JPRB
     transmissioncld % tau_surf (0,ichan) = Exp(-1.0_JPRB * zod_up_cld (ichan) / angles (iprof) % coszen)
  Enddo
  
  Call rttov_calcemis_mw(      &
       & opts,                 &! in
       & profiles,             &! in
       & angles,               &! in
       & coef_rttov,           &! in
       & chanprof,            &! in
       & transmissioncld,      &! in
       & calcemiss,            &! in
       & scatt_aux % ems_cld,  &! inout
       & scatt_aux % ref_cld,  &! out
       & errorstatus          ) ! inout 

!* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux % ems_bnd (:) = scatt_aux % ems_cld (:)
  scatt_aux % ref_bnd (:) = scatt_aux % ref_cld (:)

!* ADJOINT PART
!* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux_ad % ems_cld (:) = scatt_aux_ad % ems_cld (:) + scatt_aux_ad % ems_bnd (:)
  scatt_aux_ad % ems_bnd (:) = 0.0_JPRB
  
  scatt_aux_ad % ref_cld (:) = scatt_aux_ad % ref_cld (:) + scatt_aux_ad % ref_bnd (:)  
  scatt_aux_ad % ref_bnd (:) = 0.0_JPRB
  
  transmissioncld_ad % tau_surf (0,:) = 0.0_JPRB

  if (adk == adk_adjoint) then   
    Call rttov_calcemis_mw_ad(           &
            & opts,                      &! in
            & profiles,                  &! in
            & profiles_ad,               &! inout
            & angles,                    &! in
            & coef_rttov,                &! in
            & chanprof,                 &! in
            & transmissioncld    ,       &! in
            & transmissioncld_ad,        &! in
            & calcemiss,                 &! in
            & scatt_aux_ad % ems_cld,    &! inout
            & scatt_aux_ad % ref_cld)     ! inout 
  else if (adk == adk_k) then 
    Call rttov_calcemis_mw_k(            &
            & opts,                      &! in
            & profiles,                  &! in
            & profiles_ad,               &! inout
            & angles,                    &! in
            & coef_rttov,                &! in
            & chanprof,                 &! in
            & transmissioncld    ,       &! in
            & transmissioncld_ad,        &! in
            & calcemiss,                 &! in
            & scatt_aux_ad % ems_cld,    &! inout
            & scatt_aux_ad % ref_cld)     ! inout 
  endif
 
  zod_up_cld_ad (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     zod_up_cld_ad (ichan) = zod_up_cld_ad (ichan) - transmissioncld_ad % tau_surf (0,ichan) &
                         & * transmissioncld % tau_surf (0,ichan) / angles (iprof) % coszen
     transmissioncld_ad % tau_surf (0,ichan) = 0.0_JPRB

     if (zod_up_cld (ichan) == 30.0_JPRB) zod_up_cld_ad (ichan) = 0.0_JPRB

     Do ilayer = 1, nlevels
        iprof = chanprof(ichan) % prof
        if (adk == adk_adjoint) then
          iprofad = iprof  
        else if (adk == adk_k) then
          iprofad = ichan  
        endif
   
        scatt_aux_ad % ext (ichan,ilayer) = scatt_aux_ad % ext (ichan,ilayer) & 
                                        & + scatt_aux % dz  (iprof,ilayer) * zod_up_cld_ad (ichan)
        scatt_aux_ad % dz  (iprofad,ilayer) = scatt_aux_ad % dz  (iprofad,ilayer) & 
                                        & + scatt_aux % ext (ichan,ilayer) * zod_up_cld_ad (ichan)
     Enddo
  Enddo
  zod_up_cld_ad (:) = 0.0_JPRB

  scatt_aux % ext (:,:) = ext_2 (:,:) 
  scatt_aux % ssa (:,:) = ssa_2 (:,:) 
  scatt_aux % asm (:,:) = asm_2 (:,:) 

!* Scattering parameters for Eddington RT
  Call rttov_iniedd_ad(     &
       & nlevels,           &! in
       & nchannels ,        &! in
       & nprofiles ,        &! in
       & nprofilesad,       &! in
       & chanprof%prof,     &! in
       & angles    ,        &! in
       & scatt_aux ,        &! inout
       & scatt_aux_ad)       ! inout 
       
!* Cloud/rain absorption/scattering parameters
  scatt_aux % ext (:,:) = ext_1 (:,:) 
  scatt_aux % ssa (:,:) = ssa_1 (:,:) 
  scatt_aux % asm (:,:) = asm_1 (:,:) 

  Call rttov_mieproc_ad (&
       & nlevels,           &! in
       & nchannels,         &! in
       & nprofiles,         &! in
       & nprofilesad,       &! in
       & frequencies,       &! in
       & chanprof%prof,     &! in
       & profiles,          &! in
       & coef_scatt,        &! in
       & scatt_aux,         &! inout
       & scatt_aux_ad)       ! inout 
       
  if(usenewcld) then

    ! Calculate a hydrometeor-weighted average cloudy sky fraction 
    Do iprofad = 1, nprofilesad

      if (adk == adk_adjoint) then
        iprof = iprofad
      else if (adk == adk_k) then
        iprof = chanprof(iprofad) % prof
      endif

      !* Partition all cloud and rain into the cloudy column
      If (scatt_aux % cfrac (iprof) > ccthres) Then

        scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad)        &
          & - sum ( scatt_aux_ad % clw (iprofad,:)  * clw_precf  (iprof,:)   &
          & +       scatt_aux_ad % ciw (iprofad,:)  * ciw_precf  (iprof,:)   &
          & +       scatt_aux_ad % totalice (iprofad,:)  * totalice_precf  (iprof,:)   &
          & +       scatt_aux_ad % rain (iprofad,:) * rain_precf (iprof,:)   &
          & +       scatt_aux_ad % sp (iprofad,:)   * sp_precf   (iprof,:) ) &
          & / ( scatt_aux % cfrac (iprof) ** 2 )

        scatt_aux_ad % clw  (iprofad,:) = scatt_aux_ad % clw  (iprofad,:) / scatt_aux % cfrac (iprof)
        scatt_aux_ad % ciw  (iprofad,:) = scatt_aux_ad % ciw  (iprofad,:) / scatt_aux % cfrac (iprof)
        scatt_aux_ad % totalice (iprofad,:) = scatt_aux_ad %totalice (iprofad,:) / scatt_aux % cfrac (iprof)
        scatt_aux_ad % rain (iprofad,:) = scatt_aux_ad % rain (iprofad,:) / scatt_aux % cfrac (iprof) 
        scatt_aux_ad % sp   (iprofad,:) = scatt_aux_ad % sp   (iprofad,:) / scatt_aux % cfrac (iprof) 

      Else
        scatt_aux_ad % clw  (iprofad,:) = 0.0_JPRB
        scatt_aux_ad % ciw  (iprofad,:) = 0.0_JPRB
        scatt_aux_ad % totalice  (iprofad,:) = 0.0_JPRB
        scatt_aux_ad % rain (iprofad,:) = 0.0_JPRB
        scatt_aux_ad % sp   (iprofad,:) = 0.0_JPRB
      Endif

      if( usercfrac ) then

        !* User-supplied cloud fraction
        cld_profiles_ad (iprofad) % cfrac = scatt_aux_ad % cfrac (iprofad)  
        scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

      else

        !* Weighted mean cloud fraction
        if (hydro_column(iprof) > 1e-10_JPRB) then

          if ( cfrac(iprof) < 0.0_JPRB ) scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB
          if ( cfrac(iprof) > 1.0_JPRB ) scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

          cld_profiles_ad (iprofad) % cc (:) = cld_profiles_ad (iprofad) % cc (:) &
            & + scatt_aux_ad % cfrac (iprofad) * hydro_weights(iprof,:) / hydro_column(iprof)

          hydro_weights_ad(:) = scatt_aux_ad % cfrac (iprofad) * cld_profiles (iprof) % cc (:) &
            & / hydro_column(iprof)

          hydro_column_ad = -1.0_JPRB * scatt_aux_ad % cfrac (iprofad) &
            & * sum(hydro_weights(iprof,:) * cld_profiles (iprof) % cc (:)) &
            & / ( hydro_column(iprof) ** 2 )

        else
          hydro_weights_ad(:) = 0.0_JPRB
          hydro_column_ad     = 0.0_JPRB
        endif

        scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB

        !* Partial column of hydrometeors in g/m^2

        hydro_weights_ad(:) = hydro_weights_ad(:) + hydro_column_ad
        hydro_column_ad = 0.0_JPRB

        scatt_aux_ad % rain (iprofad,:) = scatt_aux_ad % rain (iprofad,:) &
                 & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
        scatt_aux_ad % sp (iprofad,:)   = scatt_aux_ad % sp (iprofad,:)   &
                 & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
        scatt_aux_ad % ciw (iprofad,:)  = scatt_aux_ad % ciw (iprofad,:)  &
                 & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
        scatt_aux_ad % totalice (iprofad,:)  = scatt_aux_ad % totalice (iprofad,:)  &
                 & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)
        scatt_aux_ad % clw (iprofad,:)  = scatt_aux_ad % clw (iprofad,:)  &
                 & + hydro_weights_ad(:) * scatt_aux % dz (iprof,:)

        scatt_aux_ad % dz (iprofad,:) = scatt_aux_ad % dz (iprofad,:) &
          & + hydro_weights_ad(:) * ( clw_precf (iprof,:) + ciw_precf (iprof,:) &
          & + rain_precf (iprof,:) + sp_precf (iprof,:) + totalice_precf (iprof,:))

        hydro_weights_ad(:) = 0.0_JPRB

      endif

    Enddo

  Endif

!* Change units
  presfh_ad (:,:) = 0.0_JPRB
  presf_ad  (:,:) = 0.0_JPRB
  
  Do ilayer = 1,nlevels
     Do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
          iprof = iprofad  
        else if (adk == adk_k) then
          iprof = chanprof(iprofad) % prof  
        endif

!* Rates from kg/m^2/s to g/m^3
        if (sp_scale   (iprof,ilayer) > 0.0_JPRB) then 
           scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) &
             & * (coef_scatt % conv_sp (2)) * (sp_scale (iprof,ilayer)**(coef_scatt % conv_sp (2) - 1.0_JPRB)) &
             & * (coef_scatt % conv_sp (1))**(coef_scatt % conv_sp (2)) 
        else
           scatt_aux_ad % sp   (iprofad,ilayer) = 0.0_JPRB
        endif

        if (rain_scale (iprof,ilayer) > 0.0_JPRB) then 
           scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) &
             & * (coef_scatt % conv_rain (2)) * (rain_scale (iprof,ilayer)**(coef_scatt % conv_rain (2) - 1.0_JPRB)) &
             & * (coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
        else
            scatt_aux_ad % rain (iprofad,ilayer) = 0.0_JPRB
        endif

        scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) * 3600.0_JPRB
        scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) * 3600.0_JPRB 

        scatt_aux_ad % sp   (iprofad,ilayer) = scatt_aux_ad % sp   (iprofad,ilayer) / rho_snow
        scatt_aux_ad % rain (iprofad,ilayer) = scatt_aux_ad % rain (iprofad,ilayer) / rho_rain

!* Condensate from g/g to g/m^3

        presf_ad (iprofad,ilayer) = presf_ad (iprofad,ilayer) + ciw_scale (iprof,ilayer) & 
                                           & * de2mr / profiles (iprof) % t (ilayer) * scatt_aux_ad % ciw (iprofad,ilayer) 
        profiles_ad (iprofad) % t (ilayer) =  profiles_ad (iprofad) % t (ilayer) & 
                                           & - ciw_scale (iprof,ilayer) * presf (iprof,ilayer) * de2mr &
                                           & / (profiles    (iprof) % t (ilayer) & 
                                           & * profiles (iprof) % t (ilayer)) * scatt_aux_ad % ciw (iprofad,ilayer) 
        scatt_aux_ad % ciw (iprofad,ilayer) = presf(iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer) & 
                                           & * scatt_aux_ad % ciw (iprofad,ilayer) 

        presf_ad (iprofad,ilayer) = presf_ad (iprofad,ilayer) + totalice_scale (iprof,ilayer) & 
                                           & * de2mr / profiles (iprof) % t (ilayer) * scatt_aux_ad % totalice (iprofad,ilayer) 
        profiles_ad (iprofad) % t (ilayer) =  profiles_ad (iprofad) % t (ilayer) & 
                                           & - totalice_scale (iprof,ilayer) * presf (iprof,ilayer) * de2mr &
                                           & / (profiles    (iprof) % t (ilayer) & 
                                           & * profiles (iprof) % t (ilayer)) * scatt_aux_ad % totalice (iprofad,ilayer) 
        scatt_aux_ad % totalice (iprofad,ilayer) = presf(iprof,ilayer) * de2mr / profiles (iprof) % t (ilayer) & 
                                           & * scatt_aux_ad % totalice (iprofad,ilayer) 

        presf_ad (iprofad,ilayer) = presf_ad (iprofad,ilayer) + clw_scale (iprof,ilayer) & 
                                           & * de2mr / profiles (iprof) % t (ilayer) * scatt_aux_ad % clw (iprofad,ilayer)   
        profiles_ad (iprofad) % t (ilayer) =  profiles_ad (iprofad) % t (ilayer) & 
                                           & - clw_scale (iprof,ilayer) * presf (iprof,ilayer) * de2mr &
                                           & / (profiles    (iprof) % t (ilayer) & 
                                           & * profiles (iprof) % t (ilayer)) * scatt_aux_ad % clw (iprofad,ilayer)    
        scatt_aux_ad % clw (iprofad,ilayer) = presf (iprof,ilayer) * de2mr / profiles(iprof) % t (ilayer) & 
                                           & * scatt_aux_ad % clw (iprofad,ilayer)
  
     Enddo

!* Optical depths in km-1 and at nadir
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
        if (adk == adk_adjoint) then
          iprofad = iprof  
        else if (adk == adk_k) then
          iprofad = ichan  
        endif
     
        If (ext_0 (ichan,ilayer) < 1.0E-10_JPRB) scatt_aux_ad % ext (ichan,ilayer) = 0.0_JPRB

        od_ad (ichan,ilayer) = scatt_aux_ad % ext (ichan,ilayer) &
          & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen   
        scatt_aux_ad % dz (iprofad,ilayer) = scatt_aux_ad % dz (iprofad,ilayer) &
          & - od (ichan,ilayer) * angles (iprof) % coszen &
          & / (scatt_aux % dz (iprof,ilayer) * scatt_aux % dz (iprof,ilayer)) &
          & * scatt_aux_ad % ext (ichan,ilayer) 
        scatt_aux_ad % ext (ichan,ilayer) = 0.0_JPRB
     Enddo
  Enddo

  dzbot_ad(:,:)    = 0.0_JPRB
  dztop_ad(:,:)    = 0.0_JPRB
  od_rttov_ad(:,:) = 0.0_JPRB

  Do ichan = 1, nchannels

     iprof = chanprof(ichan) % prof
     if (adk == adk_adjoint) then
        iprofad = iprof  
     else if (adk == adk_k) then
        iprofad = ichan  
     endif

     dzr_ad(:)         = 0.0_JPRB

     Do ilayer = 2, nlevels  
       dzr(ilayer) = dzbot(iprof,ilayer-1) + dztop(iprof,ilayer)
     Enddo

! Re-allocate optical depths between half pressure levels 
     od_rttov_ad(ichan,1) = od_rttov_ad(ichan,1) + od_ad (ichan,1)  
     od_rttov_ad(ichan,2) = od_rttov_ad(ichan,2) + od_ad (ichan,1) &
                        & * dzbot(iprof,1) / dzr(2) 
     dzbot_ad(iprofad,1)  = dzbot_ad(iprofad,1)  + od_ad (ichan,1) &
                        & * od_rttov(ichan,2) / dzr(2) 
     dzr_ad(2)            = dzr_ad(2) - od_ad (ichan,1) * od_rttov(ichan,2) &
                        & * dzbot(iprof,1) / dzr(2)**2   
   
     Do ilayer = 2, nlevels - 1  
        od_rttov_ad(ichan,ilayer)   = od_rttov_ad(ichan,ilayer)                     & 
          & + od_ad (ichan,ilayer) * dztop(iprof,ilayer) / dzr(ilayer)
        dztop_ad(iprofad,ilayer)    = dztop_ad(iprofad,ilayer)                      &  
          & + od_ad (ichan,ilayer) * od_rttov(ichan,ilayer) / dzr(ilayer)
        dzr_ad(ilayer)              = dzr_ad(ilayer)                                &  
          & - od_ad (ichan,ilayer) * od_rttov(ichan,ilayer) * dztop(iprof,ilayer)   &
          & / dzr(ilayer)**2
        od_rttov_ad(ichan,ilayer+1) = od_rttov_ad(ichan,ilayer+1)                   &
          & + od_ad (ichan,ilayer) * dzbot(iprof,ilayer) / dzr(ilayer+1)
        dzbot_ad(iprofad,ilayer)    = dzbot_ad(iprofad,ilayer)                      &
          & + od_ad (ichan,ilayer) * od_rttov(ichan,ilayer+1) / dzr(ilayer+1)     
        dzr_ad(ilayer+1)            = dzr_ad(ilayer+1)                              &         
          & - od_ad (ichan,ilayer) * od_rttov(ichan,ilayer+1) * dzbot(iprof,ilayer) &
          & / dzr(ilayer+1)**2
     Enddo

     od_rttov_ad(ichan,nlevels)   = od_rttov_ad(ichan,nlevels)                    &
       & + od_ad (ichan,nlevels) * dztop(iprof,nlevels) / dzr(nlevels)
     dztop_ad(iprofad,nlevels)    = dztop_ad(iprofad,nlevels)                     &
       & + od_ad (ichan,nlevels) * od_rttov(ichan,nlevels) / dzr(nlevels)
     dzr_ad(nlevels)              = dzr_ad(nlevels)                               &             
       & - od_ad (ichan,nlevels) * od_rttov(ichan,nlevels) * dztop(iprof,nlevels) &
       & / dzr(nlevels)**2
     od_rttov_ad(ichan,nlevels+1) = od_rttov_ad(ichan,nlevels+1) + od_ad (ichan,nlevels)

     od_ad (ichan,:) = 0.0_JPRB

! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzbot_ad(iprofad,ilayer-1) = dzbot_ad(iprofad,ilayer-1) + dzr_ad(ilayer)
       dztop_ad(iprofad,ilayer)   = dztop_ad(iprofad,ilayer)   + dzr_ad(ilayer)
       dzr_ad(ilayer)           = 0.0_JPRB
     Enddo

     transmission_ad % tau_levels (nlevels,ichan) =     &
          & + transmission_ad % tau_levels (nlevels,ichan) &
          & + od_rttov_ad(ichan,nlevels+1) / transmission % tau_levels (nlevels,ichan) 

     transmission_ad % tau_total (ichan) =     &
          & + transmission_ad % tau_total (ichan) &
          & - od_rttov_ad(ichan,nlevels+1) / transmission % tau_total (ichan) 
     
     Do ilayer = nlevels, 2, -1

        transmission_ad % tau_levels (ilayer-1,ichan) =     &
          & + transmission_ad % tau_levels (ilayer-1,ichan) &
          & + od_rttov_ad(ichan,ilayer) / transmission % tau_levels (ilayer-1,ichan) 

        transmission_ad % tau_levels (ilayer,ichan) =     &
          & + transmission_ad % tau_levels (ilayer,ichan) &
          & - od_rttov_ad(ichan,ilayer) / transmission % tau_levels (ilayer,ichan) 
       
     Enddo

     transmission_ad % tau_levels (1,ichan) = transmission_ad % tau_levels (1,ichan) &
                   & - od_rttov_ad (ichan,1) / transmission % tau_levels (1,ichan)

  Enddo

  od_rttov_ad (:,:) = 0.0_JPRB  
 
!* Nadir heights (km)
  Do ilayer = 1, nlevels
     Do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
          iprof = iprofad  
        else if (adk == adk_k) then
          iprof = chanprof(iprofad) % prof  
        endif

        p1 = presfh (iprof,ilayer+1)
        p2 = presfh (iprof,ilayer  )
        pm = presf  (iprof,ilayer  )

        p1_ad = 0.0_JPRB
        p2_ad = 0.0_JPRB
        pm_ad = 0.0_JPRB

        p2_ad = p2_ad + dp2dz / p2 * profiles (iprof) % t (ilayer) * dztop_ad (iprofad,ilayer)
        pm_ad = pm_ad - dp2dz / pm * profiles (iprof) % t (ilayer) * dztop_ad (iprofad,ilayer)
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                       & + dp2dz * Log(p2/pm) * dztop_ad (iprofad,ilayer) 
        dztop_ad (iprofad,ilayer) = 0.0_JPRB

        pm_ad = pm_ad + dp2dz / pm * profiles (iprof) % t (ilayer) * dzbot_ad (iprofad,ilayer)
        p1_ad = p1_ad - dp2dz / p1 * profiles (iprof) % t (ilayer) * dzbot_ad (iprofad,ilayer)
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                       & + dp2dz * Log(pm/p1) * dzbot_ad (iprofad,ilayer) 
        dzbot_ad (iprofad,ilayer) = 0.0_JPRB

        p2_ad = p2_ad + dp2dz / p2 * profiles (iprof) % t (ilayer) * scatt_aux_ad % dz (iprofad,ilayer)
        p1_ad = p1_ad - dp2dz / p1 * profiles (iprof) % t (ilayer) * scatt_aux_ad % dz (iprofad,ilayer)
        profiles_ad (iprofad) % t (ilayer) = profiles_ad (iprofad) % t (ilayer) & 
                                       & + dp2dz * Log(p2/p1) * scatt_aux_ad % dz (iprofad,ilayer) 
        scatt_aux_ad % dz (iprofad,ilayer) = 0.0_JPRB

        presf_ad  (iprofad,ilayer)   = presf_ad  (iprofad,ilayer)   + pm_ad
        presfh_ad (iprofad,ilayer)   = presfh_ad (iprofad,ilayer)   + p2_ad
        presfh_ad (iprofad,ilayer+1) = presfh_ad (iprofad,ilayer+1) + p1_ad
     Enddo
  Enddo

  !* Initialise cloud and rain properties of the cloudy/rainy column
  if(usenewcld) then

    ! Weighted average partitioning between "cloudy" and "clear" columns
    Do ilayer=1,nlevels
      Do iprofad = 1, nprofilesad
        if (adk == adk_adjoint) then
          iprof = iprofad  
        else if (adk == adk_k) then
          iprof = chanprof(iprofad) % prof  
        endif

        cld_profiles_ad (iprofad) % clw  (ilayer) = cld_profiles_ad (iprofad) % clw  (ilayer) & 
                                         & + scatt_aux_ad % clw  (iprofad,ilayer) 
        cld_profiles_ad (iprofad) % rain (ilayer) = cld_profiles_ad (iprofad) % rain (ilayer) & 
                                         & + scatt_aux_ad % rain (iprofad,ilayer) 
        if ( cld_profiles (iprof) % use_totalice ) then
          cld_profiles_ad (iprofad) % totalice  (ilayer) = cld_profiles_ad (iprofad) % totalice  (ilayer) & 
                                         & + scatt_aux_ad % totalice  (iprofad,ilayer) 
        else
          cld_profiles_ad (iprofad) % ciw  (ilayer) = cld_profiles_ad (iprofad) % ciw  (ilayer) & 
                                         & + scatt_aux_ad % ciw  (iprofad,ilayer) 
          cld_profiles_ad (iprofad) % sp   (ilayer) = cld_profiles_ad (iprofad) % sp   (ilayer) & 
                                         & + scatt_aux_ad % sp   (iprofad,ilayer) 
        endif
      Enddo
    Enddo

  else

    ! Maximum cloud fraction partitioning between "cloudy" and "clear" columns
    Do ilayer=1,nlevels
      Do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
          iprof = iprofad  
        else if (adk == adk_k) then
          iprof = chanprof(iprofad) % prof  
        endif

        If (scatt_aux % cfrac (iprof) > ccthres) Then

          cld_profiles_ad (iprofad) % clw (ilayer) = cld_profiles_ad (iprofad) % clw (ilayer) & 
            & + scatt_aux_ad % clw (iprofad,ilayer) / scatt_aux % cfrac (iprof)
          cld_profiles_ad (iprofad) % rain (ilayer) = cld_profiles_ad (iprofad) % rain (ilayer) & 
            & + scatt_aux_ad % rain (iprofad,ilayer) / scatt_aux % cfrac (iprof)

          if ( cld_profiles (iprof) % use_totalice ) then

            cld_profiles_ad (iprofad) % totalice (ilayer) = cld_profiles_ad (iprofad) % totalice (ilayer) &
              & + scatt_aux_ad % totalice  (iprofad,ilayer) / scatt_aux % cfrac (iprof)

          else

            cld_profiles_ad (iprofad) % ciw  (ilayer) = cld_profiles_ad (iprofad) % ciw (ilayer) & 
              & + scatt_aux_ad % ciw  (iprofad,ilayer) / scatt_aux % cfrac (iprof)
            cld_profiles_ad (iprofad) % sp   (ilayer) = cld_profiles_ad (iprofad) % sp  (ilayer) & 
              & + scatt_aux_ad % sp   (iprofad,ilayer) / scatt_aux % cfrac (iprof)

          endif

          scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad) &
            & - (cld_profiles (iprof) % clw  (ilayer) * scatt_aux_ad % clw  (iprofad,ilayer)  &
            & +  cld_profiles (iprof) % rain (ilayer) * scatt_aux_ad % rain (iprofad,ilayer)) &
            & / (scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof))

          if ( cld_profiles (iprof) % use_totalice ) then

            scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad) &
              & - (cld_profiles (iprof) % totalice (ilayer) * scatt_aux_ad % totalice (iprofad,ilayer)) &
              & / (scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof))

          else

            scatt_aux_ad % cfrac (iprofad) = scatt_aux_ad % cfrac (iprofad) &
              & - (cld_profiles (iprof) % ciw  (ilayer) * scatt_aux_ad % ciw  (iprofad,ilayer) &
              & +  cld_profiles (iprof) % sp   (ilayer) * scatt_aux_ad % sp   (iprofad,ilayer)) &
              & / (scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof))

          endif
       
        Endif
     
        if (iccmax(iprof) >0) then
          cld_profiles_ad (iprofad) % cc (iccmax (iprof)) =   &
            & cld_profiles_ad (iprofad) % cc (iccmax (iprof)) &
            & + scatt_aux_ad % cfrac (iprofad)
        endif 
        scatt_aux_ad % cfrac (iprofad) = 0.0_JPRB    
      Enddo
    Enddo

  endif
  
  scatt_aux_ad % clw   (:,:) = 0.0_JPRB
  scatt_aux_ad % ciw   (:,:) = 0.0_JPRB
  scatt_aux_ad % totalice  (:,:) = 0.0_JPRB
  scatt_aux_ad % rain  (:,:) = 0.0_JPRB
  scatt_aux_ad % sp    (:,:) = 0.0_JPRB

!* Temperature at layer boundaries (K)
  Do ilayer = nlevels - 1, 1, -1
     Do iprofad = 1, nprofilesad

        if (adk == adk_adjoint) then
           iprof = iprofad  
        else if (adk == adk_k) then
           iprof = chanprof(iprofad) % prof  
        endif

        p1 = presf  (iprof,ilayer+1)
        p2 = presf  (iprof,ilayer  )
        pm = presfh (iprof,ilayer+1)

        profiles_ad (iprofad) % t (ilayer+1) = profiles_ad (iprofad) % t (ilayer+1) &
         & + scatt_aux_ad % tbd (iprofad,ilayer+1)                                  &
         & - 1.0_JPRB / Log(p2/p1) * Log(pm/p1) * scatt_aux_ad % tbd (iprofad,ilayer+1) 
        profiles_ad (iprofad) % t (ilayer)   = profiles_ad (iprofad) % t (ilayer)   &
         & + 1.0_JPRB / Log(p2/p1) * Log(pm/p1) * scatt_aux_ad % tbd (iprofad,ilayer+1) 
     
        p1_ad = (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & / (Log(p2/p1) * Log(p2/p1)) / p1 * Log(pm/p1) * scatt_aux_ad % tbd (iprofad,ilayer+1) &
            & - (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & /  Log(p2/p1) / p1 * scatt_aux_ad % tbd (iprofad,ilayer+1)     
        p2_ad = -1.0_JPRB &
            & * (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & *  Log(pm/p1) / (Log(p2/p1) * Log(p2/p1)) / p2 * scatt_aux_ad % tbd (iprofad,ilayer+1) 
        pm_ad = (profiles (iprof) % t (ilayer) - profiles (iprof) % t (ilayer+1)) &
            & /  Log(p2/p1) / pm * scatt_aux_ad % tbd (iprofad,ilayer+1) 
        scatt_aux_ad % tbd (iprofad,ilayer+1) = 0.0_JPRB

        presf_ad  (iprofad,ilayer+1) = presf_ad  (iprofad,ilayer+1) + p1_ad
        presf_ad  (iprofad,ilayer  ) = presf_ad  (iprofad,ilayer  ) + p2_ad
        presfh_ad (iprofad,ilayer+1) = presfh_ad (iprofad,ilayer+1) + pm_ad
     Enddo
  Enddo

  Do iprofad = 1, nprofilesad
     profiles_ad (iprofad) % s2m % t = profiles_ad (iprofad) % s2m % t + scatt_aux_ad % tbd (iprofad,nlevels+1) 
     profiles_ad (iprofad) % t (1)   = profiles_ad (iprofad) % t (1)   + scatt_aux_ad % tbd (iprofad,1) 
  Enddo
  scatt_aux_ad % tbd (:,:) = 0.0_JPRB

!* Security on user-defined pressures
  Do iprofad = 1, nprofilesad

     if (adk == adk_adjoint) then
        iprof = iprofad  
     else if (adk == adk_k) then
        iprof = chanprof(iprofad) % prof 
     endif

     Do ilayer = 1, nlevels
        If (profiles    (iprof) % p (ilayer) >= pressure_top) &
         &  profiles_ad (iprofad) % p (ilayer) = profiles_ad (iprofad) % p (ilayer) + presf_ad (iprofad,ilayer) 
        presf_ad  (iprofad,ilayer) = 0.0_JPRB
     Enddo
     Do ilayer = 1, nlevels + 1
        If (cld_profiles    (iprof) % ph (ilayer) >= pressure_top) &
         &  cld_profiles_ad (iprofad) % ph (ilayer) = cld_profiles_ad (iprofad) % ph (ilayer) + presfh_ad (iprofad,ilayer) 
        presfh_ad (iprofad,ilayer) = 0.0_JPRB
     Enddo
  Enddo

  scatt_aux % ext (:,:) = ext_3 (:,:) 
  scatt_aux % ssa (:,:) = ssa_3 (:,:) 
  scatt_aux % asm (:,:) = asm_3 (:,:) 

!* Deallocate
  deallocate (transmissioncld    % tau_surf)
  deallocate (transmissioncld_ad % tau_surf)

  if (lhook) call dr_hook('RTTOV_INISCATT_AD',1_jpim,zhook_handle)
 
End Subroutine rttov_iniscatt_ad
