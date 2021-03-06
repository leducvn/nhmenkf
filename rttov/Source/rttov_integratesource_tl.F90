!      
Subroutine rttov_integratesource_tl (&        
     & nlevels,       &! in
     & nchannels,     &! in
     & nprofiles,     &! in
     & lprofiles,     &! in
     & angles,        &! in
     & scatt_aux,     &! in
     & scatt_aux_tl,  &! in
     & dp,            &! in
     & dp_tl,         &! in
     & dm,            &! in
     & dm_tl,         &! in
     & j_do,          &! inout
     & j_do_tl,       &! inout
     & j_up,          &! inout
     & j_up_tl)        ! inout 

  ! Description:
  ! integrate source in Eddington
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
  ! - Bauer, P., 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part I: Model description.
  !     NWP SAF Report No. NWPSAF-EC-TR-005, 27 pp.
  ! - Chevallier, F., and P. Bauer, 2003:
  !     Model rain and clouds over oceans: comparison with SSM/I observations.
  !     Mon. Wea. Rev., 131, 1240-1255.
  ! - Moreau, E., P. Bauer and F. Chevallier, 2002: Microwave radiative transfer modeling in clouds and precipitation.
  !     Part II: Model evaluation.
  !     NWP SAF Report No. NWPSAF-EC-TR-006, 27 pp.
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0       09/2002   Initial version     (E. Moreau)
  !  1.1       05/2003   RTTOV7.3 compatible (F. Chevallier)
  !  1.2       03/2004   Added polarimetry   (R. Saunders)
  !  1.3       08/2004   Polarimetry fixes   (U. O'Keefe)
  !  1.4       11/2004   Clean-up            (P. Bauer)
  !  1.5       11/2007   RTTOV9 version      (A. Geer)
  !  1.6       07/2008   Speed-ups / tidied  (A. Geer)
  !  1.7       03/2010   Use mclayer rather than min_ssa (A. Geer)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !   Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  ! Imported Type Definitions:
  
  Use rttov_types, Only :    &
       & profile_scatt_aux    ,&
       & geometry_Type 

  Use rttov_const, Only: min_ssa, ccthres

  USE YOMHOOK, ONLY: LHOOK , DR_HOOK

  Use parkind1, Only : jpim     ,jprb

  Implicit None

!* Subroutine arguments:
  Integer (Kind=jpim), Intent (in) :: nlevels               ! Number of levels
  Integer (Kind=jpim), Intent (in) :: nprofiles             ! Number of profiles
  Integer (Kind=jpim), Intent (in) :: nchannels             ! Number of channels*profiles=radiances
  Integer (Kind=jpim), Intent (in) :: lprofiles (nchannels) ! Profile indices

  Type (profile_scatt_aux), Intent(in) :: scatt_aux         ! Auxiliary profile variables for RTTOV_SCATT
  Type (profile_scatt_aux), Intent(in) :: scatt_aux_tl      ! Auxiliary profile variables for RTTOV_SCATT
  Type (geometry_Type),     Intent(in) :: angles (nprofiles)! Zenith angles  

  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dp       ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dm       ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do     ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up     ! Upward source terms
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dp_tl    ! D+ for boundary conditions
  Real (Kind=jprb), Intent (in)   , dimension (nchannels,nlevels) :: dm_tl    ! D- for boundary conditions
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_do_tl  ! Downward source terms 
  Real (Kind=jprb), Intent (inout), dimension (nchannels,nlevels) :: j_up_tl  ! Upward source terms

!INTF_END

!* Local variables
  Real    (Kind=jprb) :: ja   , jb   , jc   , jd   , aa   , bb   , cp   , cm   , ztmp
  Real    (Kind=jprb) :: ja_tl, jb_tl, jc_tl, jd_tl, aa_tl, bb_tl, cp_tl, cm_tl, ztmp_tl
  Integer (Kind=jpim) :: iprof, ichan, ilayer

  REAL(KIND=JPRB) :: ZHOOK_HANDLE

  !- End of header --------------------------------------------------------

  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATESOURCE_TL',0_jpim,ZHOOK_HANDLE)

!* Channels * Profiles      
 do ilayer=1,nlevels
  do ichan = 1, nchannels
     iprof = lprofiles (ichan)
    
     if (ilayer >= scatt_aux % mclayer(ichan) .and. & 
       & scatt_aux % cfrac (iprof) > ccthres ) then 

!* Coefficients
     aa_tl  =  scatt_aux_tl % b0  (iprof,ilayer) - 1.5_JPRB * angles( iprof) % coszen &
             & * (scatt_aux_tl % asm (ichan,ilayer) * scatt_aux    % ssa (ichan,ilayer) * scatt_aux    % b1 (iprof,ilayer)  &
             & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer) * scatt_aux    % b1 (iprof,ilayer)  &
             & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux    % ssa (ichan,ilayer) * scatt_aux_tl % b1 (iprof,ilayer)) &
             & /  scatt_aux % h (ichan,ilayer) &
             & + 1.5_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) &
             & *  angles (iprof) % coszen * scatt_aux % b1 (iprof,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
             & / (scatt_aux % h (ichan,ilayer) * scatt_aux % h (ichan,ilayer))             
     aa    = scatt_aux % b0 (iprof,ilayer) - 1.5_JPRB * scatt_aux % asm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) &
            & * angles (iprof) % coszen * scatt_aux % b1 (iprof,ilayer) / scatt_aux % h (ichan,ilayer) 
            
     bb_tl  = scatt_aux_tl % b1 (iprof,ilayer)
     bb     = scatt_aux    % b1 (iprof,ilayer)

     cp_tl  = (dp_tl (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) + dp (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer)) &
             & * (1.0_JPRB - 1.5_JPRB * angles (iprof) % coszen * scatt_aux % asm (ichan,ilayer) * &
             &scatt_aux % lambda (ichan,ilayer) &
             & / scatt_aux % h (ichan,ilayer) ) &
             & - 1.5_JPRB * dp (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * angles (iprof) % coszen  &
             & * (scatt_aux_tl % asm (ichan,ilayer) * scatt_aux    % lambda (ichan,ilayer) / scatt_aux    % h( ichan,ilayer) &
             & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux_tl % lambda (ichan,ilayer) / scatt_aux    % h (ichan,ilayer) &
             & -  scatt_aux    % asm (ichan,ilayer) * scatt_aux    % lambda (ichan,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
             & / (scatt_aux    % h   (ichan,ilayer) * scatt_aux    % h      (ichan,ilayer))) 
     cp     = dp (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * (1.0_JPRB - 1.5_JPRB * scatt_aux % asm (ichan,ilayer) &
             & * angles (iprof) % coszen * scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer)) 
          
     cm_tl  = (dm_tl (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) + dm (ichan,ilayer) * scatt_aux_tl % ssa (ichan,ilayer)) &
             & * (1.0_JPRB + 1.5_JPRB * angles (iprof) % coszen * scatt_aux % asm (ichan,ilayer) * &
             &scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer))&
             & + 1.5_JPRB * dm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * angles (iprof) % coszen &
             & * (scatt_aux_tl % asm (ichan,ilayer) * scatt_aux    % lambda (ichan,ilayer) / scatt_aux    % h (ichan,ilayer) &
             & +  scatt_aux    % asm (ichan,ilayer) * scatt_aux_tl % lambda (ichan,ilayer) / scatt_aux    % h (ichan,ilayer) &
             & -  scatt_aux    % asm (ichan,ilayer) * scatt_aux    % lambda (ichan,ilayer) * scatt_aux_tl % h (ichan,ilayer) &
             & / (scatt_aux    % h   (ichan,ilayer) * scatt_aux    % h      (ichan,ilayer))) 
     cm     = dm (ichan,ilayer) * scatt_aux % ssa (ichan,ilayer) * (1.0_JPRB + 1.5_JPRB * scatt_aux % asm (ichan,ilayer) &
             & * angles (iprof) % coszen * scatt_aux % lambda (ichan,ilayer) / scatt_aux % h (ichan,ilayer)) 

!* Downward radiance source terms    
        ja_tl  = -1.0_JPRB * scatt_aux_tl % tau (ichan,ilayer)
        ja     =  1.0_JPRB - scatt_aux    % tau (ichan,ilayer)
 
        jb_tl  = -1.0_JPRB * angles(iprof) % coszen &
                & * (scatt_aux_tl % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) * scatt_aux % ext (ichan,ilayer))*&
                &(1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
                & +  scatt_aux_tl % tau (ichan,ilayer) /  scatt_aux % ext (ichan,ilayer)) &
                & -  scatt_aux_tl % tau (ichan,ilayer) *  scatt_aux % dz  (iprof,ilayer) - &
                & scatt_aux % tau (ichan,ilayer) * scatt_aux_tl % dz (iprof,ilayer) 
        jb     = angles (iprof) % coszen / scatt_aux % ext (ichan,ilayer) * (1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
                & - scatt_aux % tau (ichan,ilayer) * scatt_aux % dz (iprof,ilayer) 

        ztmp     = exp (scatt_aux % dz (iprof,ilayer) * (scatt_aux % lambda (ichan,ilayer) - scatt_aux % ext (ichan,ilayer) / &
         &angles (iprof) % coszen))
        ztmp_tl  = ztmp * (scatt_aux_tl % dz (iprof,ilayer) * (scatt_aux    % lambda (ichan,ilayer) - &
         &scatt_aux    % ext (ichan,ilayer) / angles (iprof) % coszen) &
                &              + scatt_aux    % dz (iprof,ilayer) * (scatt_aux_tl % lambda (ichan,ilayer) - &
                &scatt_aux_tl % ext (ichan,ilayer) / angles (iprof) % coszen)) 

        jc_tl = (scatt_aux_tl % ext (ichan,ilayer)    * scatt_aux % lambda (ichan,ilayer) &
            & -  scatt_aux_tl % lambda (ichan,ilayer) * scatt_aux % ext (ichan,ilayer)) &
            & * ( ztmp - 1.0_JPRB ) * angles (iprof) % coszen &
            & / ((scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen - scatt_aux % ext (ichan,ilayer)) **2 )&
            & + ztmp_tl * scatt_aux % ext (ichan,ilayer) &
            & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen - scatt_aux % ext (ichan,ilayer))     
        jc    =  scatt_aux % ext (ichan,ilayer) * (ztmp - 1.0_JPRB) &
            & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen - scatt_aux % ext (ichan,ilayer)) 

        ztmp    = exp (scatt_aux % dz (iprof,ilayer) * (scatt_aux % lambda (ichan,ilayer) + scatt_aux % ext (ichan,ilayer) / &
                & angles (iprof) % coszen))
        ztmp_tl = ztmp * ( scatt_aux_tl % dz (iprof,ilayer) & 
              & * (scatt_aux % lambda (ichan,ilayer) + scatt_aux % ext (ichan,ilayer) / angles (iprof) % coszen) &
              & + scatt_aux % dz (iprof,ilayer) * (scatt_aux_tl % lambda (ichan,ilayer) &
              & + scatt_aux_tl % ext (ichan,ilayer) / angles (iprof) % coszen)) 

        jd_tl = (scatt_aux_tl % ext (ichan,ilayer)    * scatt_aux % lambda (ichan,ilayer) &
            & -  scatt_aux_tl % lambda (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) )  &
            & * (1.0_JPRB - 1.0_JPRB / ztmp) * angles (iprof) % coszen &
            & / ((scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen + scatt_aux % ext (ichan,ilayer)) ** 2) &
            & +  scatt_aux % ext    (ichan,ilayer) * ztmp_tl  / ztmp  / ztmp &
            & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen + scatt_aux % ext (ichan,ilayer))    
        jd    =    scatt_aux % ext  (ichan,ilayer) * (1.0_JPRB - 1.0_JPRB / ztmp ) &
            & / (scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen + scatt_aux % ext (ichan,ilayer)) 

        j_do_tl (ichan,ilayer) = ja_tl  * aa  + ja  * aa_tl  &
                        & + jb_tl  * bb  + jb  * bb_tl  &
                        & + jc_tl  * cp  + jc  * cp_tl  &
                        & + jd_tl  * cm  + jd  * cm_tl  
        j_do    (ichan,ilayer) = ja  * aa  + jb  * bb  + jc  * cp  + jd  * cm 

!* Upward radiance source terms    

        ja_tl  = -1.0_JPRB * scatt_aux_tl % tau (ichan,ilayer)
        ja     =  1.0_JPRB - scatt_aux    % tau (ichan,ilayer)
       
        jb_tl  = angles (iprof) % coszen  &
                & * (scatt_aux_tl % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) * scatt_aux    % ext (ichan,ilayer)) * &
                &(1.0_JPRB - scatt_aux % tau (ichan,ilayer)) &
                & +  scatt_aux_tl % tau (ichan,ilayer) /  scatt_aux % ext (ichan,ilayer)) + scatt_aux_tl % dz  (iprof,ilayer) 
        jb     =  scatt_aux    % dz  (iprof,ilayer) - angles (iprof) % coszen / scatt_aux % ext (ichan,ilayer) * &
         &(1.0_JPRB - scatt_aux % tau (ichan,ilayer)) 

        ztmp     = exp (scatt_aux % dz (iprof,ilayer) * scatt_aux % lambda (ichan,ilayer))
        ztmp_tl  = (scatt_aux_tl % dz (iprof,ilayer) * scatt_aux    % lambda (ichan,ilayer) &
             &   +  scatt_aux    % dz (iprof,ilayer) * scatt_aux_tl % lambda (ichan,ilayer)) * ztmp 

        jc_tl  =  (  scatt_aux_tl % ext   (ichan,ilayer) * scatt_aux % lambda (ichan,ilayer)   &
                & -  scatt_aux_tl % lambda(ichan,ilayer) * scatt_aux % ext    (ichan,ilayer) ) &
                & * ( ztmp - scatt_aux % tau (ichan,ilayer)) * angles(iprof) % coszen &
                & / ( (scatt_aux % ext (ichan,ilayer) + scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) ** 2) &
                & + (ztmp_tl - scatt_aux_tl % tau (ichan,ilayer)) * scatt_aux % ext (ichan,ilayer) &
                & / (scatt_aux % ext (ichan,ilayer) + scatt_aux % lambda (ichan,ilayer) * angles(iprof) % coszen ) 
        jc     = scatt_aux % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) + scatt_aux % lambda (ichan,ilayer) &
                & * angles (iprof) % coszen) * (ztmp - scatt_aux % tau (ichan,ilayer)) 

        jd_tl  =  ( (scatt_aux_tl % ext (ichan,ilayer) * (-1.0_JPRB) * scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen &
                & + scatt_aux_tl % lambda (ichan,ilayer) * scatt_aux % ext (ichan,ilayer) * angles (iprof) % coszen ) &
                & * (1.0_JPRB/ztmp - scatt_aux % tau (ichan,ilayer) ) &
                & / (scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) &
                & +  scatt_aux % ext (ichan,ilayer) * (-1.0_JPRB * ztmp_tl / ztmp / ztmp  - scatt_aux_tl % tau (ichan,ilayer))) &
                & / (scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) * angles (iprof) % coszen) 
        jd     = scatt_aux    % ext (ichan,ilayer) / (scatt_aux % ext (ichan,ilayer) - scatt_aux % lambda (ichan,ilayer) &
                & * angles (iprof) % coszen) * (1.0_JPRB / ztmp  - scatt_aux % tau (ichan,ilayer)) 

        j_up_tl (ichan,ilayer) = ja_tl  * aa  + ja  * aa_tl  &
                        & + jb_tl  * bb  + jb  * bb_tl &
                        & + jc_tl  * cp  + jc  * cp_tl &
                        & + jd_tl  * cm  + jd  * cm_tl  
        j_up    (ichan,ilayer) = ja  * aa  + jb  * bb  + jc  * cp  + jd  * cm     
     end if 
  end do
 end do
  
  IF (LHOOK) CALL DR_HOOK('RTTOV_INTEGRATESOURCE_TL',1_jpim,ZHOOK_HANDLE)

End subroutine rttov_integratesource_tl
