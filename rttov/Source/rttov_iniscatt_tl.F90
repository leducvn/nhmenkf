!
Subroutine rttov_iniscatt_tl (&
      & errorstatus,       &! out
      & nlevels,           &! in
      & nchannels,         &! in
      & nprofiles,         &! in
      & chanprof,          &! in
      & frequencies,       &! in
      & profiles,          &! in  
      & profiles_tl,       &! in  
      & cld_profiles,      &! in 
      & cld_profiles_tl,   &! in 
      & coef_rttov,        &! in
      & coef_scatt,        &! in
      & transmission,      &! in
      & transmission_tl,   &! in
      & calcemiss,         &! in
      & usenewcld,         &! in
      & usercfrac,         &! in
      & angles,            &! out
      & scatt_aux,         &! inout
      & scatt_aux_tl)       ! inout 

  !
  ! Description:
  ! Calculates some variables related to the input precipitation profile
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
  !   1.3    08/2004      Polarimetry fixes   (U. O'Keeffe)
  !   1.4    11/2004      Clean-up            (P. Bauer)
  !   1.5    10/2005      Fixes for rttov8 indexing   (U. O'Keeffe)
  !   1.6    11/2005      Add errorstatus to arguments (J. Cameron)   
  !   1.7    09/2006      Use zccmax_tl instead of iccmax index (A. Doherty)
  !   1.8    11/2007      RTTOV9 / cleanup    (A. Geer) 
  !   1.9    03/2008      Revised cloud partitioning (A. Geer)  
  !   1.10   03/2009      Safety check on cloud fraction (A. Geer)
  !   1.11   11/2009      User may supply average cloud fraction (A. Geer)
  !   1.12   11/2009      RTTOV transmittances / optical depths come on levels (A Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declaratiochannelsns:
  ! Modules used:
  ! Imported Type Definitions:
  
  USE YOMHOOK, ONLY: LHOOK , DR_HOOK
  
  Use rttov_types, Only :     &
       & rttov_coef          ,&
       & rttov_scatt_coef    ,&
       & transmission_type   ,&
       & transmission_type_aux  ,&
       & geometry_Type       ,&
       & profile_scatt_aux   ,&
       & profile_Type        ,&
       & profile_cloud_Type  ,&
       & rttov_chanprof      ,&
       & rttov_options 

  Use rttov_const, Only:      &
       & errorstatus_success, &
       & errorstatus_fatal,   &
       & gravity,             &
       & pressure_top,        &
       & rgp,                 &
       & rm,                  &
       & rho_rain,            &
       & rho_snow,            &
       & ccthres 

  Use parkind1, Only : jpim, jprb, jplm

  Implicit None

!* Subroutine arguments:

  Integer (Kind=jpim), Intent (in) :: nlevels    ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles  ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels  ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent(out) :: errorstatus             ! Error return code
  Integer (Kind=jpim), Intent (in) :: frequencies (nchannels) ! Frequency indices
  Type(rttov_chanprof), Intent(in) :: chanprof    (nchannels) ! Channel and profile indices
  Logical (Kind=jplm), Intent (in) :: calcemiss   (nchannels) ! Emissivity flags
  Logical (Kind=jplm), Intent (in) :: usenewcld               ! New or old cloud partition
  Logical (Kind=jplm), Intent (in) :: usercfrac               ! User has supplied cloud fraction

  Type (profile_Type),        Intent (in)    :: profiles        (nprofiles)   ! Atmospheric profiles
  Type (profile_Type),        Intent (in)    :: profiles_tl     (nprofiles)   ! Atmospheric profiles
  Type (rttov_coef),          Intent (in)    :: coef_rttov                    ! RTTOV Coefficients
  Type (rttov_scatt_coef),    Intent (in)    :: coef_scatt                    ! RTTOV_SCATT Coefficients
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles    (nprofiles)   ! Cloud profiles 
  Type (profile_cloud_Type),  Intent (in)    :: cld_profiles_tl (nprofiles)   ! Cloud profiles 
  Type (transmission_Type),   Intent (in)    :: transmission                  ! Transmittances and optical depths
  Type (transmission_Type),   Intent (in)    :: transmission_tl               ! Transmittances and optical depths
  Type (geometry_Type),       Intent (out)   :: angles          (nprofiles)   ! Zenith angles
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux                     ! Auxiliary profile variables
  Type (profile_scatt_aux),   Intent (inout) :: scatt_aux_tl                  ! Auxiliary profile variables

!INTF_END

!* Local variables
  Integer (Kind=jpim) :: ilayer, iprof, ichan
  Real    (Kind=jprb) :: p1, p2, pm, p1_tl, p2_tl, pm_tl, dp2dz, de2mr, zccmax, zccmax_tl

  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf  ! Pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh ! Half-level pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop  ! Thickness of top half of level [m]
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot  ! Thickness of bottom half of level [m]  
  Real (Kind=jprb), Dimension (2:nlevels)           :: dzr    ! Thickness of RTTOV level [m]
  Real (Kind=jprb), Dimension (nchannels,nlevels)   :: od     ! Single layer optical depths on our levels
  Real (Kind=jprb), Dimension (nlevels+1)           :: od_rttov ! Single layer optical depths on RTTOV levels
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dztop_tl    
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: dzbot_tl  
  Real (Kind=jprb), Dimension (2:nlevels)           :: dzr_tl  
  Real (Kind=jprb), Dimension (nchannels,nlevels)   :: od_tl       
  Real (Kind=jprb), Dimension (nlevels+1)           :: od_rttov_tl 
  Real (Kind=jprb), Dimension (nchannels)           :: zod_up_cld   ! Optical depth from top of the atmosphere 
  Real (Kind=jprb), Dimension (nprofiles,nlevels)   :: presf_tl     ! Pressure levels [hPa]
  Real (Kind=jprb), Dimension (nprofiles,nlevels+1) :: presfh_tl    ! Half-level pressure levels [hPa]
  Real (Kind=jprb), Dimension (nchannels)           :: zod_up_cld_tl    ! Optical depth from top of the atmosphere 

  Real    (Kind=jprb) :: hydro_weights(nlevels), hydro_column              ! Total hydrometeor amounts [g m^-2]
  Real    (Kind=jprb) :: hydro_weights_tl(nlevels), hydro_column_tl        

  Type (transmission_Type_aux) :: transmissioncld, transmissioncld_tl             ! Clear+cloud transmittances with cloud
  Type(rttov_options)  :: opts

  Character (len=80) :: errMessage
  Character (len=18) :: NameOfRoutine = 'rttov_iniscatt_tl '
  
  REAL(KIND=JPRB) :: ZHOOK_HANDLE

#include "rttov_mieproc_tl.h"
#include "rttov_iniedd_tl.h"
#include "rttov_calcemis_mw.h"
#include "rttov_calcemis_mw_tl.h"
#include "rttov_setgeometry.h"
#include "rttov_errorreport.h"
  !- End of header --------------------------------------------------------

  if (lhook) call dr_hook('RTTOV_INISCATT_TL',0_jpim,zhook_handle)

  errorstatus = errorstatus_success

  allocate (transmissioncld    % tau_surf (0:0,nchannels))
  allocate (transmissioncld_tl % tau_surf (0:0,nchannels))

  de2mr =  1.0E+05_JPRB * rm / rgp
  dp2dz = -1.0E-03_JPRB * rgp / gravity / rm 

  scatt_aux    % ext (:,:) = 0.0_JPRB
  scatt_aux    % ssa (:,:) = 0.0_JPRB
  scatt_aux    % asm (:,:) = 0.0_JPRB
  scatt_aux_tl % ext (:,:) = 0.0_JPRB
  scatt_aux_tl % ssa (:,:) = 0.0_JPRB
  scatt_aux_tl % asm (:,:) = 0.0_JPRB
  
!* Security on user-defined pressures
  Do iprof = 1, nprofiles
     Do ilayer = 1, nlevels
        If (profiles (iprof) % p (ilayer) >= pressure_top) Then
            presf_tl (iprof,ilayer) = profiles_tl (iprof) % p (ilayer)
            presf    (iprof,ilayer) = profiles    (iprof) % p (ilayer)
        else
            presf_tl (iprof,ilayer) = 0.0_JPRB
            presf    (iprof,ilayer) = pressure_top
        Endif
     Enddo
     Do ilayer = 1, nlevels + 1
        If (cld_profiles(iprof) % ph (ilayer) >= pressure_top) Then
            presfh_tl (iprof,ilayer) = cld_profiles_tl (iprof) % ph (ilayer)
            presfh    (iprof,ilayer) = cld_profiles    (iprof) % ph (ilayer)
        else
            presfh_tl (iprof,ilayer) = 0.0_JPRB
            presfh    (iprof,ilayer) = pressure_top
        Endif
     Enddo
  Enddo

!* Geometric variables
  Call rttov_setgeometry ( &
    & opts,       & ! in
    & profiles,        & ! in
    & coef=coef_rttov, & ! in
    & angles=angles)     ! out

!* Temperature at layer boundaries (K)
  Do iprof = 1, nprofiles
     scatt_aux_tl % tbd (iprof,nlevels+1) = profiles_tl (iprof) % s2m % t
     scatt_aux    % tbd (iprof,nlevels+1) = profiles    (iprof) % s2m % t
     scatt_aux_tl % tbd (iprof,1)         = profiles_tl (iprof) % t (1)
     scatt_aux    % tbd (iprof,1)         = profiles    (iprof) % t (1)
  Enddo

  Do ilayer = 1, nlevels-1
     Do iprof = 1, nprofiles
        p1_tl = presf_tl  (iprof,ilayer+1)
        p1    = presf     (iprof,ilayer+1)
        p2_tl = presf_tl  (iprof,ilayer  )
        p2    = presf     (iprof,ilayer  )
        pm_tl = presfh_tl (iprof,ilayer+1)
        pm    = presfh    (iprof,ilayer+1)

        scatt_aux_tl % tbd (iprof,ilayer+1) =  profiles_tl (iprof) % t (ilayer+1)   &
                                          & + (profiles_tl (iprof) % t (ilayer)     &
                                          & -  profiles_tl (iprof) % t (ilayer+1))  &
                                          & / log(p2/p1) * log(pm/p1)                   &
                                          & + (profiles    (iprof) % t (ilayer)     &
                                          & -  profiles    (iprof) % t (ilayer+1))  &
                                          & / (-1.0_JPRB *  log(p2/p1) * log(p2/p1) )   &
                                          & * (p2_tl / p2 - p1_tl / p1) * log(pm/p1)    & 
                                          & + (profiles    (iprof) % t  (ilayer)    &
                                          & -  profiles    (iprof) % t  (ilayer+1)) &
                                          & / log(p2/p1) * (pm_tl / pm - p1_tl / p1) 
        scatt_aux    % tbd (iprof,ilayer+1) =  profiles    (iprof) % t  (ilayer+1)  &
                                          & + (profiles    (iprof) % t  (ilayer)    &
                                          & -  profiles    (iprof) % t  (ilayer+1)) &
                                          & / log(p2/p1) * log(pm/p1) 
     Enddo
  Enddo

!* Initialise cloud and rain properties of the cloudy/rainy column
  scatt_aux % clw   (:,:) = 0.0_JPRB
  scatt_aux % ciw   (:,:) = 0.0_JPRB
  scatt_aux % totalice   (:,:) = 0.0_JPRB
  scatt_aux % rain  (:,:) = 0.0_JPRB
  scatt_aux % sp    (:,:) = 0.0_JPRB
  scatt_aux_tl % clw   (:,:) = 0.0_JPRB
  scatt_aux_tl % ciw   (:,:) = 0.0_JPRB
  scatt_aux_tl % totalice   (:,:) = 0.0_JPRB
  scatt_aux_tl % rain  (:,:) = 0.0_JPRB
  scatt_aux_tl % sp    (:,:) = 0.0_JPRB

  if(usenewcld) then

    ! Weighted average partitioning between "cloudy" and "clear" columns
    Do ilayer=1,nlevels
      Do iprof = 1, nprofiles  
        scatt_aux % clw  (iprof,ilayer) = cld_profiles (iprof) % clw  (ilayer) 
        scatt_aux_tl % clw  (iprof,ilayer) = cld_profiles_tl (iprof) % clw  (ilayer) 
        if ( cld_profiles (iprof) % use_totalice ) then
          scatt_aux % totalice (iprof,ilayer) = cld_profiles (iprof) % totalice (ilayer) 
          scatt_aux_tl % totalice (iprof,ilayer) = cld_profiles_tl (iprof) % totalice (ilayer) 
        else
          scatt_aux % ciw  (iprof,ilayer) = cld_profiles (iprof) % ciw  (ilayer) 
          scatt_aux % sp   (iprof,ilayer) = cld_profiles (iprof) % sp   (ilayer) 
          scatt_aux_tl % ciw  (iprof,ilayer) = cld_profiles_tl (iprof) % ciw  (ilayer) 
          scatt_aux_tl % sp   (iprof,ilayer) = cld_profiles_tl (iprof) % sp   (ilayer) 
        endif 
        scatt_aux % rain (iprof,ilayer) = cld_profiles (iprof) % rain (ilayer) 
        scatt_aux_tl % rain (iprof,ilayer) = cld_profiles_tl (iprof) % rain (ilayer) 
      Enddo
    Enddo

  else

    ! Maximum cloud fraction partitioning between "cloudy" and "clear" columns
    Do iprof = 1, nprofiles  
      zccmax    = 0.0_JPRB
      zccmax_tl = 0.0_JPRB

      Do ilayer = 1, nlevels
        if (cld_profiles (iprof) % cc (ilayer) > zccmax) then
          zccmax    = cld_profiles    (iprof) % cc (ilayer)
          zccmax_tl = cld_profiles_tl (iprof) % cc (ilayer) 
        end if 
      end do 
      scatt_aux % cfrac (iprof) = zccmax
      scatt_aux_tl % cfrac (iprof) = zccmax_tl
    Enddo

    do ilayer=1,nlevels
      do iprof = 1, nprofiles  
        If (scatt_aux % cfrac (iprof) > ccthres) Then
          scatt_aux_tl % clw  (iprof,ilayer) = &
            &  (cld_profiles_tl (iprof) % clw (ilayer) * scatt_aux % cfrac (iprof)  &
            & - scatt_aux_tl % cfrac (iprof) * cld_profiles (iprof) % clw (ilayer)) &
            & /(scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof)) 
          scatt_aux % clw  (iprof,ilayer) = cld_profiles (iprof) % clw  (ilayer) / scatt_aux % cfrac (iprof)
          scatt_aux_tl % rain (iprof,ilayer) = &
            &  (cld_profiles_tl (iprof) % rain (ilayer) * scatt_aux % cfrac (iprof)  &
            & - scatt_aux_tl % cfrac (iprof) * cld_profiles (iprof) % rain (ilayer)) &
            & /(scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof)) 
          scatt_aux % rain (iprof,ilayer) = cld_profiles (iprof) % rain (ilayer) / scatt_aux % cfrac (iprof)
          if ( cld_profiles (iprof) % use_totalice ) then
            scatt_aux_tl % totalice (iprof,ilayer) = &
              &  (cld_profiles_tl (iprof) % totalice (ilayer) * scatt_aux % cfrac (iprof) & 
              & - scatt_aux_tl % cfrac (iprof) * cld_profiles (iprof) %totalice (ilayer)) &
              & /(scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof)) 
            scatt_aux % totalice (iprof,ilayer) = cld_profiles (iprof) % totalice (ilayer) / scatt_aux % cfrac (iprof)
          else
            scatt_aux_tl % ciw  (iprof,ilayer) = &
              &  (cld_profiles_tl (iprof) % ciw (ilayer) * scatt_aux % cfrac (iprof)  &
              & - scatt_aux_tl % cfrac (iprof) * cld_profiles (iprof) % ciw (ilayer)) &
              & / (scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof)) 
            scatt_aux % ciw (iprof,ilayer) = cld_profiles (iprof) % ciw  (ilayer) / scatt_aux % cfrac (iprof)
            scatt_aux_tl % sp (iprof,ilayer) = &
              &  (cld_profiles_tl (iprof) % sp (ilayer) * scatt_aux % cfrac (iprof)  &
              & - scatt_aux_tl % cfrac (iprof) * cld_profiles (iprof) % sp (ilayer)) &
              & / (scatt_aux % cfrac (iprof) * scatt_aux % cfrac (iprof)) 
            scatt_aux % sp   (iprof,ilayer) = cld_profiles (iprof) % sp   (ilayer) / scatt_aux % cfrac (iprof)
          endif
        Endif
      enddo
    enddo

  endif

!* Nadir heights (km)
  Do ilayer = nlevels, 1, -1
     Do iprof = 1, nprofiles
        p1_tl = presfh_tl (iprof,ilayer+1)
        p1    = presfh    (iprof,ilayer+1)
        p2_tl = presfh_tl (iprof,ilayer  )
        p2    = presfh    (iprof,ilayer  )
        pm_tl = presf_tl  (iprof,ilayer  )
        pm    = presf     (iprof,ilayer  )

        If (p1 <= p2) then
           errorstatus = errorstatus_fatal
           Write( errMessage, '( "iniscatt : problem with user-defined pressure layering")' )
           Call rttov_errorreport (errorstatus, errMessage, NameOfRoutine)
           if (lhook) call dr_hook('RTTOV_INISCATT',1_jpim,zhook_handle)
           Return
        End If

        scatt_aux_tl % dz (iprof,ilayer) = dp2dz * (Log(p2/p1) * profiles_tl (iprof) % t (ilayer) &
                         & + (p2_tl / p2 - p1_tl / p1) * profiles (iprof) % t (ilayer))        
        scatt_aux    % dz (iprof,ilayer) = dp2dz * Log(p2/p1) * profiles (iprof) % t (ilayer) 

        dzbot_tl (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles_tl (iprof) % t (ilayer) &
                              & + dp2dz * (pm_tl / pm - p1_tl / p1) * profiles (iprof) % t (ilayer)
        dztop_tl (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles_tl (iprof) % t (ilayer) &
                              & + dp2dz * (p2_tl / p2 - pm_tl / pm) * profiles (iprof) % t (ilayer)

        dzbot (iprof,ilayer) = dp2dz * Log(pm/p1) * profiles (iprof) % t (ilayer)
        dztop (iprof,ilayer) = dp2dz * Log(p2/pm) * profiles (iprof) % t (ilayer)

     Enddo
  Enddo
 
!* Get single layer optical depths (at nadir and in hPa-1) and put onto model half levels
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
        
! Top RTTOV level to space   
     od_rttov_tl (1)      = -1.0_jprb / transmission % tau_levels (1,ichan)  &
                        & * transmission_tl % tau_levels (1,ichan) 
     od_rttov (1)         = -1.0_jprb * log( transmission % tau_levels (1,ichan) )     
     Do ilayer = 2, nlevels 
        od_rttov_tl (ilayer) = transmission_tl % tau_levels (ilayer-1,ichan) &
                        & / transmission % tau_levels (ilayer-1,ichan)    &
                        & - transmission_tl % tau_levels (ilayer,ichan)   &
                        & / transmission % tau_levels (ilayer,ichan) 
        od_rttov (ilayer) = log( transmission % tau_levels (ilayer-1,ichan) ) &
                        & - log( transmission % tau_levels (ilayer,ichan) )
     Enddo
! Surface to bottom RTTOV (full pressure) level
     od_rttov_tl (nlevels+1) = transmission_tl % tau_levels (nlevels,ichan) &
                        & / transmission % tau_levels (nlevels,ichan)    &
                        & - transmission_tl % tau_total (ichan)   &
                        & / transmission % tau_total (ichan) 
     od_rttov (nlevels+1) = log( transmission % tau_levels (nlevels,ichan) ) &
                        & - log( transmission % tau_total (ichan) )

! RTTOV layer thickness
     Do ilayer = 2, nlevels  
       dzr(ilayer)    = dzbot(iprof,ilayer-1)    + dztop(iprof,ilayer)
       dzr_tl(ilayer) = dzbot_tl(iprof,ilayer-1) + dztop_tl(iprof,ilayer)
     Enddo
 
! Re-allocate optical depths between half pressure levels       
     od_tl (ichan,1)      = od_rttov_tl(1) & 
                        & + od_rttov_tl(2) * dzbot(iprof,1) / dzr(2) &
                        & + dzbot_tl(iprof,1) * od_rttov(2) / dzr(2) &
                        & - dzr_tl(2) * od_rttov(2) * dzbot(iprof,1) / dzr(2)**2
     od (ichan,1)         = od_rttov(1) &
                        & + od_rttov(2) * dzbot(iprof,1) / dzr(2)     
     Do ilayer = 2, nlevels - 1  
        od_tl (ichan,ilayer) = od_rttov_tl(ilayer) * dztop(iprof,ilayer) / dzr(ilayer)  &
                        & + dztop_tl(iprof,ilayer) * od_rttov(ilayer)    / dzr(ilayer)  &
                        & - dzr_tl(ilayer) * od_rttov(ilayer) * dztop(iprof,ilayer)     &
                        & / dzr(ilayer)**2                                              & 
                        & + od_rttov_tl(ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) &
                        & + dzbot_tl(iprof,ilayer) * od_rttov(ilayer+1) / dzr(ilayer+1) &
                        & - dzr_tl(ilayer+1) * od_rttov(ilayer+1) * dzbot(iprof,ilayer) &
                        & / dzr(ilayer+1)**2
        od (ichan,ilayer) = od_rttov(ilayer)   * dztop(iprof,ilayer) / dzr(ilayer)   &
                        & + od_rttov(ilayer+1) * dzbot(iprof,ilayer) / dzr(ilayer+1) 
     Enddo
     od_tl (ichan,nlevels) = od_rttov_tl(nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + dztop_tl(iprof,nlevels) * od_rttov(nlevels)   / dzr(nlevels) &
                        & - dzr_tl(nlevels) * od_rttov(nlevels) * dztop(iprof,nlevels)   &
                        & / dzr(nlevels)**2                                              &
                        & + od_rttov_tl(nlevels+1)
     od (ichan,nlevels)   = od_rttov(nlevels)  * dztop(iprof,nlevels) / dzr(nlevels) &
                        & + od_rttov(nlevels+1) 


  Enddo

!* Change units
  Do ilayer = 1,nlevels

!* Optical depths in km-1 and at nadir
     Do ichan = 1, nchannels
        iprof = chanprof(ichan) % prof
     
        scatt_aux_tl % ext (ichan,ilayer) = od_tl (ichan,ilayer) &
                                        & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen &
                                        & - od    (ichan,ilayer) &
                                        & * scatt_aux_tl % dz(iprof,ilayer) / (scatt_aux % dz (iprof,ilayer) &
                                        & * scatt_aux % dz (iprof,ilayer)) &
                                        & * angles (iprof) % coszen 
        scatt_aux    % ext (ichan,ilayer) = od (ichan,ilayer)  &
                                       & / scatt_aux % dz (iprof,ilayer) * angles (iprof) % coszen 

        If (scatt_aux    % ext (ichan,ilayer) < 1.0e-10_JPRB) Then
            scatt_aux_tl % ext (ichan,ilayer) = 0.0_JPRB
            scatt_aux    % ext (ichan,ilayer) = 1.0e-10_JPRB
        Endif
     Enddo

!* Condensate from g/g to g/m^3
      Do iprof = 1, nprofiles
         scatt_aux_tl % clw (iprof,ilayer) = (scatt_aux_tl % clw (iprof,ilayer) * presf    (iprof,ilayer) &
                                         & / profiles    (iprof) % t (ilayer) &
                                         & +  scatt_aux    % clw (iprof,ilayer) * presf_tl (iprof,ilayer) &
                                         & / profiles    (iprof) % t (ilayer) &
                                         & -  scatt_aux    % clw (iprof,ilayer) * presf    (iprof,ilayer) &
                                         & * profiles_tl (iprof) % t (ilayer) &
                                         & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer))) * de2mr

         scatt_aux    % clw (iprof,ilayer) = scatt_aux % clw (iprof,ilayer) * presf (iprof,ilayer) * de2mr & 
                                         & / profiles (iprof) % t (ilayer)
     
         scatt_aux_tl % ciw (iprof,ilayer) = (scatt_aux_tl % ciw (iprof,ilayer) * presf    (iprof,ilayer) & 
                                         & / profiles    (iprof) % t (ilayer) &
                                         & +  scatt_aux    % ciw (iprof,ilayer) * presf_tl (iprof,ilayer) &
                                         & / profiles    (iprof) % t (ilayer) &
                                         & -  scatt_aux    % ciw (iprof,ilayer) * presf    (iprof,ilayer) &
                                         & * profiles_tl (iprof) % t (ilayer) &
                                         & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer))) * de2mr

         scatt_aux    % ciw (iprof,ilayer) = scatt_aux % ciw (iprof,ilayer) * presf (iprof,ilayer) * de2mr & 
                                         & / profiles (iprof) % t (ilayer) 

         scatt_aux_tl % totalice (iprof,ilayer) = (scatt_aux_tl % totalice (iprof,ilayer) * presf (iprof,ilayer) & 
                                         & / profiles    (iprof) % t (ilayer) &
                                         & +  scatt_aux    % totalice (iprof,ilayer) * presf_tl (iprof,ilayer) &
                                         & / profiles    (iprof) % t (ilayer) &
                                         & -  scatt_aux    % totalice (iprof,ilayer) * presf    (iprof,ilayer) &
                                         & * profiles_tl (iprof) % t (ilayer) &
                                         & / (profiles (iprof) % t (ilayer) * profiles (iprof) % t (ilayer))) * de2mr

         scatt_aux    % totalice (iprof,ilayer) = scatt_aux % totalice (iprof,ilayer) * presf (iprof,ilayer) * de2mr & 
                                         & / profiles (iprof) % t (ilayer) 
     
!* Rates from kg/m^2/s to g/m^3
         scatt_aux_tl % rain (iprof,ilayer) = scatt_aux_tl % rain (iprof,ilayer) / rho_rain
         scatt_aux    % rain (iprof,ilayer) = scatt_aux    % rain (iprof,ilayer) / rho_rain
         scatt_aux_tl % sp   (iprof,ilayer) = scatt_aux_tl % sp   (iprof,ilayer) / rho_snow
         scatt_aux    % sp   (iprof,ilayer) = scatt_aux    % sp   (iprof,ilayer) / rho_snow

         scatt_aux_tl % rain (iprof,ilayer) = scatt_aux_tl % rain (iprof,ilayer) * 3600.0_JPRB 
         scatt_aux    % rain (iprof,ilayer) = scatt_aux    % rain (iprof,ilayer) * 3600.0_JPRB 
         scatt_aux_tl % sp   (iprof,ilayer) = scatt_aux_tl % sp   (iprof,ilayer) * 3600.0_JPRB
         scatt_aux    % sp   (iprof,ilayer) = scatt_aux    % sp   (iprof,ilayer) * 3600.0_JPRB

         if (scatt_aux    % rain (iprof,ilayer) > 0.0_JPRB) then
             scatt_aux_tl % rain (iprof,ilayer) =  scatt_aux_tl % rain (iprof,ilayer) &
               & * (coef_scatt % conv_rain (2)) * (scatt_aux % rain (iprof,ilayer)**(coef_scatt % conv_rain (2) - 1.0_JPRB)) &
               & * (coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
             scatt_aux    % rain (iprof,ilayer) = (scatt_aux    % rain (iprof,ilayer) &
                                          & *  coef_scatt % conv_rain (1))**(coef_scatt % conv_rain (2)) 
         else
             scatt_aux_tl % rain (iprof,ilayer) = 0.0_JPRB
         endif     
         if (scatt_aux    % sp   (iprof,ilayer) > 0.0_JPRB) then
             scatt_aux_tl % sp   (iprof,ilayer) = scatt_aux_tl % sp    (iprof,ilayer) &
              & * (coef_scatt % conv_sp   (2)) * (scatt_aux % sp   (iprof,ilayer)**(coef_scatt % conv_sp   (2) - 1.0_JPRB)) &
              & * (coef_scatt % conv_sp   (1))**(coef_scatt % conv_sp   (2)) 
             scatt_aux    % sp   (iprof,ilayer) = (scatt_aux   % sp    (iprof,ilayer) &
              & *  coef_scatt % conv_sp   (1))**(coef_scatt%conv_sp     (2))
         else
             scatt_aux_tl % sp   (iprof,ilayer) = 0.0_JPRB
         endif     
     end do   
  Enddo

  if(usenewcld) then 

    !* Calculate a hydrometeor-weighted average cloudy sky fraction 
    scatt_aux % cfrac (:)   = 0.0_JPRB
    scatt_aux_tl % cfrac (:)   = 0.0_JPRB

    Do iprof = 1, nprofiles  

      if( usercfrac ) then

        !* User-supplied cloud fraction
        scatt_aux    % cfrac (iprof) = cld_profiles    (iprof) % cfrac
        scatt_aux_tl % cfrac (iprof) = cld_profiles_tl (iprof) % cfrac

      else

        !* Partial column of hydrometeors in g/m^2
        hydro_weights(:) = &
          & ( scatt_aux % rain (iprof,:) + scatt_aux % sp (iprof,:) &
          & + scatt_aux % ciw (iprof,:)  + scatt_aux % clw (iprof,:) &
          & + scatt_aux % totalice (iprof,:) ) &
          & * scatt_aux % dz (iprof,:)
        hydro_column = sum(hydro_weights)

        hydro_weights_tl(:) = &
          & ( scatt_aux_tl % rain (iprof,:) + scatt_aux_tl % sp (iprof,:) &
          & + scatt_aux_tl % ciw (iprof,:)  + scatt_aux_tl % clw (iprof,:) &
          & + scatt_aux_tl % totalice (iprof,:) ) &
          & * scatt_aux % dz (iprof,:) &
          & + ( scatt_aux % rain (iprof,:) + scatt_aux % sp (iprof,:) &
          & + scatt_aux % ciw (iprof,:)  + scatt_aux % clw (iprof,:) &
          & + scatt_aux % totalice (iprof,:) ) &
          & * scatt_aux_tl % dz (iprof,:) 
        hydro_column_tl = sum(hydro_weights_tl)

        !* Weighted mean cloud fraction    
        if (hydro_column > 1e-10_JPRB) then

          scatt_aux % cfrac (iprof) = &
            & sum(hydro_weights(:) * cld_profiles (iprof) % cc (:)) / hydro_column 

          scatt_aux_tl % cfrac (iprof) = &
            & + sum(hydro_weights_tl(:) * cld_profiles (iprof) % cc (:)) / hydro_column &
            & + sum(hydro_weights(:) * cld_profiles_tl (iprof) % cc (:)) / hydro_column &
            & - hydro_column_tl * sum(hydro_weights(:) * cld_profiles (iprof) % cc (:)) &
            & / ( hydro_column**2 )

          if ( scatt_aux % cfrac (iprof) < 0.0_JPRB ) then
            scatt_aux    % cfrac (iprof) = 0.0_JPRB
            scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
          endif

          if ( scatt_aux % cfrac (iprof) > 1.0_JPRB ) then
            scatt_aux    % cfrac (iprof) = 1.0_JPRB
            scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
          endif

        else
          scatt_aux % cfrac (iprof)    = 0.0_JPRB
          scatt_aux_tl % cfrac (iprof) = 0.0_JPRB
        endif

      endif
 
      !* Partition all cloud and rain into the cloudy column
      If (scatt_aux % cfrac (iprof) > ccthres) Then

        scatt_aux_tl % clw  (iprof,:) = &
          &   scatt_aux_tl % clw  (iprof,:) / scatt_aux % cfrac (iprof) &
          & - scatt_aux_tl % cfrac (iprof) * scatt_aux % clw  (iprof,:) &
          & / (scatt_aux % cfrac (iprof)**2)
        scatt_aux_tl % ciw  (iprof,:) = &
          &   scatt_aux_tl % ciw  (iprof,:) / scatt_aux % cfrac (iprof) &
          & - scatt_aux_tl % cfrac (iprof) * scatt_aux % ciw  (iprof,:) &
          & / (scatt_aux % cfrac (iprof)**2)
        scatt_aux_tl % totalice  (iprof,:) = &
          &   scatt_aux_tl % totalice  (iprof,:) / scatt_aux % cfrac (iprof) &
          & - scatt_aux_tl % cfrac (iprof) * scatt_aux % totalice  (iprof,:) &
          & / (scatt_aux % cfrac (iprof)**2)
        scatt_aux_tl % rain  (iprof,:) = &
          &   scatt_aux_tl % rain  (iprof,:) / scatt_aux % cfrac (iprof) &
          & - scatt_aux_tl % cfrac (iprof) * scatt_aux % rain  (iprof,:) &
          & / (scatt_aux % cfrac (iprof)**2)
        scatt_aux_tl % sp  (iprof,:) = &
          &   scatt_aux_tl % sp  (iprof,:) / scatt_aux % cfrac (iprof) &
          & - scatt_aux_tl % cfrac (iprof) * scatt_aux % sp  (iprof,:) &
          & / (scatt_aux % cfrac (iprof)**2)

        scatt_aux % clw  (iprof,:) = scatt_aux % clw  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % ciw  (iprof,:) = scatt_aux % ciw  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % totalice  (iprof,:) = scatt_aux % totalice  (iprof,:) / scatt_aux % cfrac (iprof)
        scatt_aux % rain (iprof,:) = scatt_aux % rain (iprof,:) / scatt_aux % cfrac (iprof) 
        scatt_aux % sp   (iprof,:) = scatt_aux % sp   (iprof,:) / scatt_aux % cfrac (iprof) 

      else 
        scatt_aux % clw  (iprof,:) = 0.0_JPRB
        scatt_aux % ciw  (iprof,:) = 0.0_JPRB
        scatt_aux % totalice  (iprof,:) = 0.0_JPRB
        scatt_aux % rain (iprof,:) = 0.0_JPRB
        scatt_aux % sp   (iprof,:) = 0.0_JPRB
        scatt_aux_tl % clw  (iprof,:) = 0.0_JPRB
        scatt_aux_tl % ciw  (iprof,:) = 0.0_JPRB
        scatt_aux_tl % totalice  (iprof,:) = 0.0_JPRB
        scatt_aux_tl % rain (iprof,:) = 0.0_JPRB
        scatt_aux_tl % sp   (iprof,:) = 0.0_JPRB
      Endif
    Enddo
  endif 
  
!* Cloud/rain absorption/scattering parameters
  Call rttov_mieproc_tl (   &
       & nlevels,           &! in
       & nchannels,         &! in
       & nprofiles,         &! in
       & frequencies,       &! in
       & chanprof%prof,     &! in
       & profiles,          &! in
       & coef_scatt,        &! in
       & scatt_aux,         &! inout
       & scatt_aux_tl)       ! inout 
       
!* Scattering parameters for Eddington RT
  Call rttov_iniedd_tl(     &
       & nlevels,           &! in
       & nchannels ,        &! in
       & nprofiles ,        &! in
       & chanprof%prof,     &! in
       & angles    ,        &! in
       & scatt_aux ,        &! inout
       & scatt_aux_tl)       ! inout 
       
!* Surface emissivities
  zod_up_cld_tl (:) = 0.0_JPRB
  zod_up_cld    (:) = 0.0_JPRB
  
  Do ichan = 1, nchannels
     iprof = chanprof(ichan) % prof
     
     Do ilayer = 1, nlevels     
        zod_up_cld_tl (ichan) = zod_up_cld_tl (ichan) &
                            & + scatt_aux_tl % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer)  &
                            & + scatt_aux    % ext (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer)
        zod_up_cld    (ichan) = zod_up_cld    (ichan) &
                            & + scatt_aux    % ext (ichan,ilayer) * scatt_aux    % dz (iprof,ilayer)  
     Enddo
     if (zod_up_cld (ichan) >= 30.0_JPRB) then
         zod_up_cld    (ichan) = 30.0_JPRB
         zod_up_cld_tl (ichan) =  0.0_JPRB              
     endif
     
     transmissioncld    % tau_surf (0,ichan) = Exp(-1.0_JPRB * zod_up_cld (ichan) / angles (iprof) % coszen)
     transmissioncld_tl % tau_surf (0,ichan) = -1.0_JPRB * zod_up_cld_tl (ichan)  / angles (iprof) % coszen &
                                         & * transmissioncld % tau_surf (0,ichan)
  Enddo

  Call rttov_calcemis_mw(         &
       & opts,                    &! in
       & profiles,                &! in
       & angles,                  &! in
       & coef_rttov,              &! in
       & chanprof,               &! in
       & transmissioncld,         &! in
       & calcemiss,               &! in
       & scatt_aux % ems_cld,     &! inout
       & scatt_aux % ref_cld,     &! out
       & errorstatus          )    ! inout 
       
  Call rttov_calcemis_mw_tl(      &
       & opts,                    &! in
       & profiles,                &! in
       & profiles_tl,             &! in
       & angles,                  &! in
       & coef_rttov,              &! in
       & chanprof,               &! in
       & transmissioncld,         &! in
       & transmissioncld_tl,      &! in
       & calcemiss,               &! in
       & scatt_aux_tl % ems_cld,  &! inout
       & scatt_aux_tl % ref_cld)   ! out 

!* Hemispheric emissivity (= Fastem's effective emissivity)
  scatt_aux_tl % ems_bnd (:) = scatt_aux_tl % ems_cld (:)
  scatt_aux    % ems_bnd (:) = scatt_aux    % ems_cld (:)
  scatt_aux_tl % ref_bnd (:) = scatt_aux_tl % ref_cld (:)
  scatt_aux    % ref_bnd (:) = scatt_aux    % ref_cld (:)

!* Deallocate
  Deallocate (transmissioncld    % tau_surf)
  Deallocate (transmissioncld_tl % tau_surf)

  if (lhook) call dr_hook('RTTOV_INISCATT_TL',1_jpim,zhook_handle)

End Subroutine rttov_iniscatt_tl
