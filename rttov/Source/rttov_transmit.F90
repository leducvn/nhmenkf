!
SUBROUTINE rttov_transmit( &
            & addaerosl,                    &
            & addclouds,                    &
            & nlayers,                      &
            & chanprof,                     &
            & aux,                          &
            & coef,                         &
            & ircld,                        &
            & opdp_path,                    &
            & od_level,                     &
            & transmission,                 &
            & transmission_aux,             &
            & transmission_scatt_ir_stream, &
            & tau_ref,                      &
            & tau_ref_surf,                 &
            & tau_surf,                     &
            & tau_level)
!
! Description:
!  To calculate optical depths for a number of channels
!  and profiles from every pressure level to space.
! To interpolate optical depths on to levels of radiative transfer model
! (which, at present, entails only surface transmittance, as
! other optical depths are on *rt* levels) and to convert
! optical depths to transmittances.
!
! Code based on OPDEP and RTTAU from previous versions of RTTOV
! Only one profile per call
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
!  1.0    01/06/2005  Marco Matricardi (ECMWF):
!            --       New routine based on rttov_transmit.F90.
!                     Variable trace gases, CO2, N2O,CO and CH4
!                     have been introduced for AIRS and IASI
!  1.1    29/01/2007  Removed polarisation R Saunders
!  1.2    22/08/2007  Optimised (D Salmond)
!  1.3    04/06/2008  Fix od_frac and od_frac_ac calculation near surface level (PB PM)
!  1.4    27/02/2009  Profile levels to include ToA. Distinguish between
!                     layer arrays and level arrays - size, index
!                     labels, looping (P. Rayer)
!  1.5    03/11/2009  Transmittances on levels (A Geer)
!  1.6    08/03/2010  Missing levsurf assignment (N Bormann)
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
  USE rttov_types, ONLY :  &
       & rttov_chanprof,             &
       & rttov_coef,                 &
       & opdp_path_Type,             &
       & transmission_Type,          &
       & transmission_Type_aux,      &
       & transmission_scatt_ir_type, &
       & profile_aux,                &
       & ircld_type
  USE parkind1, ONLY : jpim, jprb, jplm
!INTF_OFF
  USE rttov_const, ONLY : sensor_id_hi, max_optical_depth, min_od, min_tau
  USE yomhook, ONLY : LHOOK, DR_HOOK
!INTF_ON
  IMPLICIT NONE
!subroutine arguments:
  LOGICAL(KIND=jplm), INTENT(IN)    :: addaerosl
  LOGICAL(KIND=jplm), INTENT(IN)    :: addclouds
  INTEGER(KIND=jpim)              , INTENT(IN)    :: nlayers                                     ! Number of pressure levels
  TYPE(rttov_chanprof            ), INTENT(IN)    :: chanprof(:)                                 ! Channel indices
  TYPE(rttov_coef                ), INTENT(IN)    :: coef                                        ! Coefficients
  TYPE(opdp_path_Type            ), INTENT(INOUT) :: opdp_path
  TYPE(transmission_Type         ), INTENT(INOUT) :: transmission                                ! Transmittances and single-layer od
  TYPE(transmission_Type_aux     ), INTENT(INOUT) :: transmission_aux                            ! Transmittances and single-layer od
  TYPE(transmission_scatt_ir_type), INTENT(INOUT) :: transmission_scatt_ir_stream
  TYPE(profile_aux               ), INTENT(IN)    :: aux                                         ! auxillary profiles informations
  TYPE(ircld_type                ), INTENT(IN)    :: ircld
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_ref     (nlayers + 1   , size(chanprof))
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_ref_surf(size(chanprof)                )
  REAL(KIND=jprb)                 , INTENT(OUT)   :: od_level    (nlayers + 1   , size(chanprof))! sat to level optical depth
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_level   (nlayers + 1   , size(chanprof))! sat to level transmittance
  REAL(KIND=jprb)                 , INTENT(OUT)   :: tau_surf    (size(chanprof)                )! sat to surfacetransmittance
!INTF_END
!local variables:
  REAL   (KIND=jprb) :: od_surf(size(chanprof))                                       ! sat to surface optical depth
  REAL   (KIND=jprb) :: od_surf_ac
  REAL   (KIND=jprb) :: od_frac       (size(chanprof)                )
  REAL   (KIND=jprb) :: od_singlelayer(nlayers       , size(chanprof))                ! sat to layer optical depth
  REAL   (KIND=jprb) :: small_val 
  
  INTEGER(KIND=jpim) :: lev         , lay    , chan   , i, j, prof, ist, klevels, klayers! loop variables
  INTEGER(KIND=jpim) :: nlevels
! cloud liquid water local variables
  INTEGER(KIND=jpim) :: levsurf
  INTEGER(KIND=jpim) :: nchannels                                                     ! Number of radiances computed (channels used * profiles)
  REAL   (KIND=JPRB) :: ZHOOK_HANDLE

#ifdef RTTOV_INTEL_MKL

  integer(kind=8) :: oldmode, mode

include '/home/h03/sa_app/intel/mkl/include/mkl_vml.f90'
!    oldmode = vmlsetmode( VML_LA )
    oldmode = vmlsetmode(IOR(VML_EP, VML_ERRMODE_ERRNO))
!    call vmlsetmode(mode)
#endif

!- End of header --------------------------------------------------------


!--------------------------------------------------------------
!1. Assemble layer optical depths and convert to transmittances
!--------------------------------------------------------------
! - optical depths here are negative except for od_singlelayer
! - in rttov_opdep, already checked that value of opticaldepth is sensible
  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT', 0_jpim, ZHOOK_HANDLE)
  small_val = (tiny(1._jprb)) ** (0.333333_jprb) ! XLF doesn't like 1/3
  nchannels = size(chanprof)
  nlevels   = nlayers + 1
! single layer optical depths local variable
  DO j = 1, nchannels
    DO lay = 1, nlayers
      od_singlelayer(lay, j) =  - (opdp_path%atm_level(lay + 1, j) - opdp_path%atm_level(lay, j))
    ENDDO
  ENDDO
! level to space optical depths - local variable
! Introduce klevels/layers to stop NEC compiler merging next 2 loops
  klevels = nlevels
  klayers = nlayers
! gamma correction of local variables
  if(any(coef%ff_gam(:) .ne. 1.0_jprb)) then
!cdir NODEP
!cdir COLLAPSE
     DO j = 1, nchannels
        chan = chanprof(j)%chan
        DO lev = 1, klevels
           od_level(lev, j) = max(coef%ff_gam(chan) * opdp_path%atm_level(lev, j), -max_optical_depth)
        ENDDO
        od_singlelayer(:, j) = coef%ff_gam(chan) * od_singlelayer(:, j)
     ENDDO
  else
!cdir NODEP
!cdir COLLAPSE
     do j=1,nchannels
        DO lev = 1, nlevels
           od_level(lev, j) = max(opdp_path%atm_level(lev, j), -max_optical_depth)
        enddo
     enddo
  endif
! On some computers when optical depth is too thick
! there is an underflow during the conversion in
! transmittances. In that case uncomment following line
! and the declaration statement of max_optical_depth

#ifdef RTTOV_INTEL_MKL
  call vdexp(nlevels*nchannels, od_level, tau_ref)
#else
  tau_ref(:,:)               = Exp(od_level(:,:))
#endif
  DO j = 1, nchannels
    DO lev = 1, nlevels
       tau_level(lev, j) = tau_ref(lev,j)
       transmission%tau_levels(lev, j) = tau_ref(lev, j)
    ENDDO
  ENDDO

  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (coef%tt_val_chn(chan) == 1) THEN
        DO lev = 1, nlevels
          IF (tau_level(lev, j) < coef%tt_a0(chan)) THEN
            tau_level(lev, j)               = coef%tt_a1(chan)
            transmission%tau_levels(lev, j) = tau_level(lev, j)
          ENDIF
        ENDDO
      ENDIF
    ENDDO
  ENDIF
!---------------------------------------------------------------------------------------
!2. Compute optical depth and transmittance at surface
!---------------------------------------------------------------------------------------
  DO j = 1, nchannels
    prof       = chanprof(j)%prof
! as defined in rttov_profaux
    levsurf    = aux%s(prof)%nearestlev_surf
! layer above this
! NB all od_level -ve
! if surface below nlevels, pfraction -ve
    od_surf(j) = od_level(levsurf, j) + aux%s(prof)%pfraction_surf * (od_level(levsurf - 1, j) - od_level(levsurf, j))
    IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
      od_frac(j) = od_surf(j) - od_level(levsurf - 1, j)
    ELSE
      od_frac(j) = od_surf(j) - od_level(levsurf, j)
    ENDIF
    tau_surf(j)     = Exp(od_surf(j))
    tau_ref_surf(j) = tau_surf(j)
  ENDDO
  IF (coef%id_sensor == sensor_id_hi .AND. coef%fmv_model_ver >= 9) THEN
    DO j = 1, nchannels
      chan = chanprof(j)%chan
      IF (coef%tt_val_chn(chan) == 1) THEN
        IF (tau_surf(j) < coef%tt_a0(chan)) THEN
          tau_surf(j) = coef%tt_a1(chan)
        ENDIF
      ENDIF
    ENDDO
  ENDIF
!---Loop over the streams-----------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      prof = chanprof(j)%prof
      levsurf    = aux%s(prof)%nearestlev_surf
! recall previous local definition of laysurf
! all arrays here based on layers, not levels
      DO ist = 0, ircld%nstream(prof)
        od_surf_ac = transmission_scatt_ir_stream%opdpac(ist, j, levsurf) + aux%s(prof)%pfraction_surf * (     &
          & transmission_scatt_ir_stream%opdpac(ist, j, levsurf - 1) - transmission_scatt_ir_stream%opdpac(ist, j, levsurf))
        IF (aux%s(prof)%pfraction_surf >= 0._jprb) THEN
          transmission_aux%od_frac_ac(ist, j) = od_surf_ac - transmission_scatt_ir_stream%opdpac(ist, j, levsurf - 1)
        ELSE
          transmission_aux%od_frac_ac(ist, j) = od_surf_ac - transmission_scatt_ir_stream%opdpac(ist, j, levsurf)
        ENDIF
        transmission_aux%tau_surf_ac(ist, j)     = Exp(-od_surf_ac)
        transmission_aux%tau_ref_surf_ac(ist, j) = transmission_aux%tau_surf_ac(ist, j)
        transmission_aux%od_frac_t(ist, j)       =  - od_frac(j) + transmission_aux%od_frac_ac(ist, j)
        IF (tau_surf(j) >= 0) THEN
          transmission_aux%tau_surf_t(ist, j) = tau_surf(j) * transmission_aux%tau_surf_ac(ist, j)
        ELSE
          transmission_aux%tau_surf_t(ist, j) = tau_surf(j)
        ENDIF
        transmission_aux%tau_ref_surf_t(ist, j) = transmission_aux%tau_surf_t(ist, j)
      ENDDO
    ENDDO
  ENDIF
!-----------------------------------------------------------------------------------------
!---------------------------------------------------------------------------------------
!3. Store transmittances for other streams and single stream for o/p
!---------------------------------------------------------------------------------------
  IF (addaerosl .OR. addclouds) THEN
    DO j = 1, nchannels
      transmission%tau_total(j) = transmission_aux%tau_surf_t(0, j)
    ENDDO
    DO j = 1, nchannels
      prof = chanprof(j)%prof! Profile index
      DO ist = 0, ircld%nstream(prof)
        DO lev = 1, nlevels
           IF (tau_level(lev, j) >= 0) THEN
            transmission_aux%tau_level(lev, ist, j) = tau_level(lev, j) * &
                 exp(-min(transmission_scatt_ir_stream%opdpac(ist, j, lev), max_optical_depth))
            transmission_aux%tau_level_r(lev, ist, j) = 1.0_jprb / transmission_aux%tau_level(lev, ist, j)
          ELSE
            transmission_aux%tau_level(lev, ist, j) = tau_level(lev, j)
            transmission_aux%tau_level_r(lev, ist, j) = 1.0_jprb / transmission_aux%tau_level(lev, ist, j)
            transmission_aux%od_level(lev, ist, j)  = od_level(lev, j)
          ENDIF

        ENDDO
        DO lay = 1, nlayers
          transmission_aux%od_singlelayer(lay, ist, j) = max(small_val,   &
            & od_singlelayer(lay, j) + transmission_scatt_ir_stream%opdpacl(ist, j, lay))
          transmission_aux%od_singlelayer_r(lay, ist, j) = 1.0_jprb / transmission_aux%od_singlelayer(lay, ist, j)
        ENDDO
      ENDDO
      DO ist = 0, ircld%nstream(prof)
        transmission_aux%od_sfrac(ist, j) = max(small_val, transmission_aux%od_frac_t(ist, j))
        transmission_aux%tau_surf(ist, j) = max(small_val, transmission_aux%tau_surf_t(ist, j))
      ENDDO
    ENDDO
    DO j = 1, nchannels
      DO lev = 1, nlevels
        IF (tau_level(lev, j) >= 0) THEN
          transmission%tau_levels(lev, j) = transmission_aux%tau_level(lev, 0, j)
        ENDIF
      ENDDO
    ENDDO
  ELSE
     DO j = 1, nchannels
        transmission%tau_total(j) = tau_surf(j)
     ENDDO
     
#ifdef RTTOV_INTEL_MKL
     call vdinv(nlevels*ircld%nstream(1)*nchannels,transmission_aux%tau_level,transmission_aux%tau_level_r)
     call vdinv(nlayers*ircld%nstream(1)*nchannels,transmission_aux%od_singlelayer,transmission_aux%od_singlelayer_r)
#endif
     
     DO j = 1, nchannels
        prof = chanprof(j)%prof! Profile index
        DO ist = 0, ircld%nstream(prof)
           transmission_aux%od_sfrac(ist, j) = max(small_val, -od_frac(j))
!           transmission_aux%od_sfrac(ist, j) = -od_frac(j)
           transmission_aux%tau_surf(ist, j) = max(small_val, tau_surf(j)) ! DAR stop changes in emissivity_ad/k for v. small tau 
                                                                           ! and stop nag from complaining for under/overflows
!           transmission_aux%tau_surf(ist, j) = tau_surf(j)                                                        

           DO lay = 1, nlayers
            transmission_aux%tau_level(lay, ist, j) = tau_level(lay, j)
#ifndef RTTOV_INTEL_MKL
            transmission_aux%tau_level_r(lay, ist, j) = 1.0_jprb / transmission_aux%tau_level(lay, ist, j)
#endif
            transmission_aux%od_level(lay, ist, j)  = od_level(lay, j)
            transmission_aux%od_singlelayer(lay, ist, j) = max(small_val, od_singlelayer(lay, j))
!            transmission_aux%od_singlelayer(lay, ist, j) = od_singlelayer(lay, j)
#ifndef RTTOV_INTEL_MKL
            transmission_aux%od_singlelayer_r(lay, ist, j) = 1.0_jprb / transmission_aux%od_singlelayer(lay, ist, j)
#endif
         ENDDO
         transmission_aux%tau_level(nlevels, ist, j) = tau_level(nlevels, j)
#ifndef RTTOV_INTEL_MKL
         transmission_aux%tau_level_r(nlevels, ist, j) = 1.0_jprb / tau_level(nlevels, j)
#endif
         transmission_aux%od_level(nlevels, ist, j)  = od_level(nlevels, j)
      ENDDO
   ENDDO
ENDIF

!transmission_aux%tau_surf_r = min(1e100_jprb, 1.0_jprb / transmission_aux%tau_surf) ! initi to zero means this is not allowed.
do i=1,nchannels
   prof = chanprof(i)%prof! Profile index
   do ist=0,ircld%nstream(prof)
      transmission_aux%tau_surf_r(ist,i) = 1.0_jprb / transmission_aux%tau_surf(ist,i)
      transmission_aux%od_sfrac_r(ist,i) = 1.0_jprb / transmission_aux%od_sfrac(ist,i)
   enddo
enddo
! Moved all of this code from rttov_integrate because it isn't to do with integrating! Because it's common in TL/AD/K code it's saved for later use. 

transmission_aux%anynegtau = -0.5_jprb

do i =1,nchannels
   prof = chanprof(i)%prof
   do ist = 0, ircld%nstream(prof)
      do lay = 1, nlayers
#ifdef RTTOV_XLF 
! IBM codepath uses fsel intrinsic that doesn't generated branched code
! DAR: if anynegtau > 0 then select a slower code path that does more tests (because something has gone wrong earlier?)
           transmission_aux%anynegtau = transmission_aux%anynegtau + fsel(transmission_aux%tau_level(lay,ist,i), 0.0_jprb, 1.0_jprb)
! DAR: tau_level differences smaller than min optical depth. Is this condition necessary?
           transmission_aux%fac(1,lay,ist,i) = fsel(transmission_aux%tau_level(lay, ist, i) - transmission_aux%tau_level(lay+1, ist, i) - min_od, &
                1.0_jprb, 0.0_jprb)
! DAR: layer ods smaller than min optical depth. This should be done (is done?) in rttov_transmit as well.
           transmission_aux%fac(1,lay,ist,i) = fsel(transmission_aux%od_singlelayer(lay,ist,i) - min_od, transmission_aux%fac(1,lay,ist,i), 0.0_jprb)
! DAR: is tau_level less than min_tau? THIS should probably be done elsewhere (rttov_transmit?)
!      is it a problem if tau_level is smaller than min_tau? Can we let this underflow gracefully? It's not being used as a divsor
!      consider removing this condition.
           transmission_aux%fac(2,lay,ist,i) = fsel(transmission_aux%tau_level(lay,ist,i) - min_tau, 1.0_jprb, 0.0_jprb)
#else
! Non-IBM codepath         
           lev = lay + 1

           if(transmission_aux%tau_level(lay,ist,i) < 0.0_jprb) then 
              transmission_aux%anynegtau = 1.0_jprb
           endif

! DAR: separated fac1 and fac2 becuase fac2 can be used separately from fac1 in Phil Watts calcs

           if(transmission_aux%od_singlelayer(lay,ist,i) < min_od) then
              transmission_aux%fac(1,lay,ist,i) = 0.0_jprb
           else
              if(transmission_aux%tau_level(lay, ist, i) - transmission_aux%tau_level(lev, ist, i) < min_od) then
                 transmission_aux%fac(1,lay,ist,i) = 0.0_jprb
              else
                 transmission_aux%fac(1,lay,ist,i) = 1.0_jprb
              endif
           endif
          
           if(transmission_aux%tau_level(lay,ist,i) < min_tau) then
              transmission_aux%fac(2,lay,ist,i) = 0.0_jprb
           else
              transmission_aux%fac(2,lay,ist,i) = 1.0_jprb
           endif
#endif  
        enddo
        if(transmission_aux%tau_level(nlevels,ist,i) < min_tau) then
           transmission_aux%fac(2,nlevels,ist,i) = 0.0_jprb
        else
           transmission_aux%fac(2,nlevels,ist,i) = 1.0_jprb
        endif
     enddo
  enddo

do i =1,nchannels
   prof = chanprof(i)%prof
   do ist = 0, ircld%nstream(prof)
#ifdef RTTOV_XLF 
      transmission_aux%surf_fac(ist,i) = fsel(transmission_aux%tau_surf(ist, i) - min_tau, 1.0_jprb, 0.0_jprb)
#else
! Non-IBM codepath         
      if(transmission_aux%tau_surf(ist, i) > min_tau) then
         transmission_aux%surf_fac(ist,i) = 1.0_jprb
      else
         transmission_aux%surf_fac(ist,i) = 0.0_jprb
      endif
#endif
   enddo
enddo

  IF (LHOOK) CALL DR_HOOK('RTTOV_TRANSMIT', 1_jpim, ZHOOK_HANDLE)
END SUBROUTINE rttov_transmit
