subroutine rttovscatt_test_one ( nchannels, coef_rttov, coef_scatt, &
                        & chanprof, &
                        & frequencies,  &
                        & emissivity, &
                        & use_totalice, &
                        & use_newcld, &
                        & use_cfrac)        
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2010, EUMETSAT, All Rights Reserved.

!
! Main part of the RTTOV-SCATT test routines
!
!       2005   Peter Bauer        First version
! 13/12/2007   Alan Geer          RTTOV9 version; streamlined
!
!-------------------------------------------------------

Use mod_rttovscatt_test, only: kflevg, kproma, fastem_land_coeff, ioin, zenangle, ioout, ioout2

USE YOMHOOK   ,ONLY : LHOOK,   DR_HOOK

Use rttov_const, only :   &
  & errorstatus_fatal,    &
  & errorstatus_success,  &
  & default_err_unit,     &
  & fastem_sp            ,&    
  & q_mixratio_to_ppmv 

Use rttov_types, only :   &
  & geometry_type        ,&
  & rttov_coefs          ,&
  & rttov_scatt_coef     ,&
  & rttov_chanprof       ,&
  & rttov_options        ,&
  & profile_type         ,&
  & profile_cloud_type   ,&
  & transmission_type    ,&
  & profile_scatt_aux    ,&
  & radiance_type   


Use parkind1, only: jpim, jprb, jplm

    IMPLICIT NONE
  

integer (kind=jpim),  intent (in) :: nchannels
real    (kind=jprb),  intent (in) , dimension (nchannels) :: emissivity    
integer (kind=jpim),  intent (in) , dimension (nchannels) :: frequencies
Type(rttov_chanprof), Intent (in) , dimension (nchannels) :: chanprof 
logical (kind=jplm), intent (in) :: use_totalice 
logical (kind=jplm), intent (in) :: use_newcld
logical (kind=jplm), intent (in) :: use_cfrac

type (rttov_coefs     ), intent (inout) :: coef_rttov        
type (rttov_scatt_coef), intent (inout) :: coef_scatt  

!INTF_END

#include "rttov_setpressure.h"
#include "rttov_scatt.h"
#include "rttov_scatt_tl.h"
#include "rttov_scatt_ad.h"
#include "rttov_alloc_rad.h"
#include "rttov_alloc_prof.h"
#include "rttov_alloc_scatt_prof.h"   

!* FORWARD   
type (profile_type)         :: profiles_d1     (kproma)
type (profile_cloud_type)   :: cld_profiles_d1 (kproma)
type (radiance_type)        :: radiance_d1  

real (kind=jprb), dimension (nchannels) :: emissivity_d1 

!* TL   
type (profile_type)        :: profiles_d2     (kproma)
type (profile_type)        :: profiles_tl     (kproma)
type (profile_type)        :: profiles_tl2     (kproma)
type (profile_cloud_type)  :: cld_profiles_d2 (kproma)
type (profile_cloud_type)  :: cld_profiles_tl (kproma)
type (profile_cloud_type)  :: cld_profiles_tl2 (kproma)
type (radiance_type)       :: radiance_d2             
type (radiance_type)       :: radiance_d3             
type (radiance_type)       :: radiance_tl  
type (radiance_type)       :: radiance_tl2  

real (kind=jprb), dimension (nchannels) :: emissivity_d2 
real (kind=jprb), dimension (nchannels) :: emissivity_tl 
real (kind=jprb), dimension (nchannels) :: emissivity_tl2 
           
!* AD   
type (profile_type)        :: profiles_ad     (kproma)
type (profile_cloud_type)  :: cld_profiles_ad (kproma)
type (radiance_type)       :: radiance_ad  

type (profile_type)        :: profiles_ad2     (kproma)
type (profile_cloud_type)  :: cld_profiles_ad2 (kproma)
type (radiance_type)       :: radiance_ad2  

real (kind=jprb), dimension (nchannels) :: emissivity_ad 
    
real (kind=jprb), dimension (nchannels) :: emissivity_ad2 
       
!* K   
type (profile_type)        :: profiles_k     (nchannels)
type (profile_cloud_type)  :: cld_profiles_k (nchannels)
type (radiance_type)       :: radiance_k  

real (kind=jprb), dimension (nchannels) :: emissivity_k

!* OTHER                      
logical (kind=jplm) :: calcemiss (nchannels)

integer (kind=jpim) :: errorstatus, erroralloc
integer (kind=jpim) :: i_lev, i_proma, i_chan, i_btout, i_lambda, i_fast 

real    (kind=jprb) :: t_2m, q_2m, u_10m, v_10m, ls, p_sfc, t_sfc, rlat, rlon
real    (kind=jprb) :: lambda, zeps, zdelta1, zdelta2, threshold, z
Real    (Kind=jprb)       :: ratio(2)
Real(Kind=jprb),    Allocatable :: radiance_total_ref (:)
REAL(KIND=JPRB) :: ZHOOK_HANDLE
real    (kind=jprb), parameter :: imposed_cfrac    = 0.74_JPRB  
real    (kind=jprb), dimension (kproma) :: cfrac  
real    (kind=jprb), dimension (kflevg) :: snow  

integer (kind=jpim), parameter :: iallocate = 1, ideallocate = 0              

!- End of header ------------------------------------------------------

IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TEST',0_jpim,ZHOOK_HANDLE)
threshold = 1.0E-8_JPRB 
  
!* FORWARD-MODEL TEST ***********************************************************************************
!* Set-up
errorstatus = errorstatus_success
emissivity_d1 (1:nchannels) = emissivity    (1:nchannels)
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) < 0.01_JPRB

call allocate_profs( kflevg, kproma, profiles_d1, cld_profiles_d1, use_totalice, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d1, kflevg-1_jpim, iallocate)   
Allocate ( radiance_total_ref  ( nchannels ) )

!* Read profiles
open (ioin,  file = 'profiles2_fmt', status = 'old')
  
do i_proma = 1, kproma

  !* Surface
  read (ioin,'(10e16.6)') &
    & rlon,            &! longitude (deg)
    & rlat,            &! latitude (deg)
    & ls,              &! land-sea mask (1=land)
    & t_sfc,           &! surface temperature (K)
    & p_sfc,           &! surface pressure (Pa)
    & t_2m,            &! 2-meter temperature (K)
    & q_2m,            &! 2-meter specific humidity (kg/kg)
    & u_10m,           &! 10-meter wind u (m/s)
    & v_10m             ! 10-meter wind u (m/s) 

  if (kflevg > 60) then  

    ! Use a 61 level profile for comparing RTTOV10 results to RTTOV9. This
    ! takes the 60 level profile and extrapolates at constant value to TOA 
  
    !* Profile    
    read (ioin,'(10e16.6)') profiles_d1 (i_proma) % t    (2:kflevg)     ! temperature (K)
    read (ioin,'(10e16.6)') profiles_d1 (i_proma) % q    (2:kflevg)     ! specific humidity (kg/kg)
    read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % cc   (2:kflevg) ! cloud cover
    read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % clw  (2:kflevg) ! liquid water (kg/kg)

    ! Cloud ice water (kg/kg)
    if (use_totalice) then 
      ! Simple test of the total ice functionality - assume cloud ice water only and ignore snow.
      read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % totalice  (2:kflevg) 
    else
      read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % ciw  (2:kflevg) 
    endif

    read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % rain (2:kflevg) ! rain (kg/m2/s)

    ! solid precipitation (kg/m2/s) 
    if (use_totalice) then 
      ! Ignore for total ice - it's in the wrong units!
      read (ioin,'(10e16.6)') snow (2:kflevg)
    else
      read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % sp   (2:kflevg) 
    endif  

    !* Convert kg/kg to ppmv
    profiles_d1 (i_proma) % q (2:kflevg) = profiles_d1 (i_proma) % q (2:kflevg) * q_mixratio_to_ppmv

    !* Populate full and half pressure arrays
    call rttov_setpressure (p_sfc, profiles_d1 (i_proma) % p(2:kflevg), & 
      & cld_profiles_d1 (i_proma) % ph(2:kflevg+1))                       
    profiles_d1     (i_proma) % p  (2:kflevg)   = profiles_d1     (i_proma) % p  (2:kflevg) * 0.01_JPRB
    cld_profiles_d1 (i_proma) % ph (2:kflevg+1) = cld_profiles_d1 (i_proma) % ph (2:kflevg+1) * 0.01_JPRB

    !* Add an explicit top level to replicate RTTOV-9 behaviour for our tests. 
    profiles_d1 (i_proma) % p (1) = 0.00498_jprb ! hPa

    !* Replicate RTTOV-9: T and q are extrapolated at constant value.
    profiles_d1 (i_proma) % q (1) = profiles_d1 (i_proma) % q (2) * 0.5 ! Avoid coef limit exceeded warnings JAH
    profiles_d1 (i_proma) % t (1) = profiles_d1 (i_proma) % t (2) 

    !* This also requires a new top level for SCATT use which didn't exist before
    !* and changes to top 2 half level pressures but should not have any major effects.
    cld_profiles_d1 (i_proma) % ph  (1) = 0.0_jprb 
    cld_profiles_d1 (i_proma) % ph  (2) = 0.01_jprb 
    cld_profiles_d1 (i_proma) % cc  (1) = 0.0_jprb
    cld_profiles_d1 (i_proma) % clw (1) = 0.0_jprb
    cld_profiles_d1 (i_proma) % rain(1) = 0.0_jprb 
    if (use_totalice) then 
      cld_profiles_d1 (i_proma) % totalice(1) = 0.0_jprb  
    else
      cld_profiles_d1 (i_proma) % ciw(1) = 0.0_jprb  
      cld_profiles_d1 (i_proma) % sp (1) = 0.0_jprb
    endif
 
  else

    read (ioin,'(10e16.6)') profiles_d1 (i_proma) % t    (1:kflevg)     ! temperature (K)
    read (ioin,'(10e16.6)') profiles_d1 (i_proma) % q    (1:kflevg)     ! specific humidity (kg/kg)
    read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % cc   (1:kflevg) ! cloud cover
    read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % clw  (1:kflevg) ! liquid water (kg/kg)

    ! Cloud ice water (kg/kg)
    if (use_totalice) then 
      ! Simple test of the total ice functionality - assume cloud ice water only and ignore snow.
      read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % totalice  (1:kflevg) 
    else
      read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % ciw  (1:kflevg) 
    endif

    read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % rain (1:kflevg) ! rain (kg/m2/s)

    ! solid precipitation (kg/m2/s) 
    if (use_totalice) then 
      ! Ignore for total ice - it's in the wrong units!
      read (ioin,'(10e16.6)') snow 
    else
      read (ioin,'(10e16.6)') cld_profiles_d1 (i_proma) % sp   (1:kflevg) 
    endif

    !* Convert kg/kg to ppmv
    profiles_d1 (i_proma) % q (1:kflevg) = profiles_d1 (i_proma) % q (1:kflevg) * q_mixratio_to_ppmv

    !* Populate full and half pressure arrays
    call rttov_setpressure (p_sfc, profiles_d1 (i_proma) % p, cld_profiles_d1 (i_proma) % ph)
    profiles_d1     (i_proma) % p  (:) = profiles_d1     (i_proma) % p  (:) * 0.01_JPRB
    cld_profiles_d1 (i_proma) % ph (:) = cld_profiles_d1 (i_proma) % ph (:) * 0.01_JPRB

  endif

  q_2m = q_2m / (1.0_JPRB - q_2m)       
  profiles_d1 (i_proma) % s2m % p = cld_profiles_d1 (i_proma) % ph(kflevg+1)
  profiles_d1 (i_proma) % s2m % q = q_2m * q_mixratio_to_ppmv
  profiles_d1 (i_proma) % s2m % o = 0.0_JPRB
  profiles_d1 (i_proma) % s2m % t = t_2m  
  profiles_d1 (i_proma) % s2m % u = u_10m 
  profiles_d1 (i_proma) % s2m % v = v_10m 
  profiles_d1 (i_proma) % s2m % wfetc = 0.0_JPRB ! Unused (solar computations only)
      
  profiles_d1 (i_proma) % skin % surftype  = int (1.0_JPRB - ls)   
  profiles_d1 (i_proma) % skin % watertype = 0  ! Ocean water   
  profiles_d1 (i_proma) % skin % t         = t_sfc 
  profiles_d1 (i_proma) % skin % fastem(:) = fastem_land_coeff (:)

  profiles_d1 (i_proma) % zenangle   = zenangle
  profiles_d1 (i_proma) % azangle    = 0.0_JPRB   ! default value ! FASTEM-3 will require an actual value here
  profiles_d1 (i_proma) % ctp        = 500.0_JPRB ! default value
  profiles_d1 (i_proma) % cfraction  = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % elevation  = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % sunzenangle  = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % sunazangle   = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % latitude     = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % longitude    = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % Be           = 0.0_JPRB   ! default value
  profiles_d1 (i_proma) % cosbk        = 0.0_JPRB   ! default value
 
  if(use_cfrac) cld_profiles_d1 (i_proma) % cfrac = imposed_cfrac
   
end do  
close (ioin)  

!* Reference forward model run
call rttov_scatt ( &
  & errorstatus,         &! out 
  & kflevg,              &! in
  & chanprof,            &! in
  & frequencies,         &! in
  & profiles_d1,         &! inout  
  & cld_profiles_d1,     &! in
  & coef_rttov,          &! in
  & coef_scatt,          &! in
  & calcemiss,           &! in
  & emissivity_d1,       &! inout
  & radiance_d1,         &! inout
  & cfrac,               &! out, diagnostic only  
  & lnewcld = use_newcld, lusercfrac = use_cfrac)

! main output:
! radiance_d1%tb        = cloud-affected Tbs
! radiance_d1%tb_clear  = clear-sky Tbs
  
! Full output
write(ioout,'(A10,i4)') 'nchan ', nchannels

write(ioout,*) '--------------------------------------------------------------------------'
write(ioout,*) '--------------------------------------------------------------------------'
write(ioout,*) 'This dataset is made of ',kproma,' ECMWF model profiles'
write(ioout,*)
if (use_totalice) then
  write(ioout,*) 'Single "total ice" content hydrometeor type used (Met Office style).'
else
  write(ioout,*) 'Separate snow flux and ice content hydrometeor types used (ECMWF style).'
endif
if (use_newcld) then
  write(ioout,*) 'New average cloud fraction approach being used.'
  if (use_cfrac) write(ioout,*) 'User imposing average cloud fraction'
else
  write(ioout,*) 'Original maximum cloud fraction approach being used.'
endif
write(ioout,*)
write(ioout,*) 'Call to RTTOV_SCATT'
write(ioout,*) '-------------------'
write(ioout,*)
write(ioout,*) 'Channel  cloudy Tb '

do i_chan = 1, nchannels
  write (ioout,'(i4,3x,30f24.16)') i_chan, radiance_d1 % bt (i_chan)
enddo

write(ioout,*)
write(ioout,*) 'Channel  clear Tb '

do i_chan = 1, nchannels
  write (ioout,'(i4,3x,30f24.16)') i_chan, radiance_d1 % bt_clear (i_chan)
enddo
  
! Partial output for comparison to reference output
write(ioout2,'(A10,i4)') 'nchan ', nchannels

write(ioout2,*) '--------------------------------------------------------------------------'
write(ioout2,*) '--------------------------------------------------------------------------'
write(ioout2,*) 'This dataset is made of ',kproma,' ECMWF model profiles'
write(ioout2,*)
if (use_totalice) then
  write(ioout2,*) 'Single "total ice" content hydrometeor type used (Met Office style).'
else
  write(ioout2,*) 'Separate snow flux and ice content hydrometeor types used (ECMWF style).'
endif
if (use_newcld) then
  write(ioout2,*) 'New average cloud fraction approach being used.'
  if (use_cfrac) write(ioout2,*) 'User imposing average cloud fraction'
else
  write(ioout2,*) 'Original maximum cloud fraction approach being used.'
endif
write(ioout2,*)
write(ioout2,*) 'Call to RTTOV_SCATT'
write(ioout2,*) '-------------------'
write(ioout2,*)
write(ioout2,*) 'Channel  cloudy Tb '

do i_chan = 1, nchannels
  write (ioout2,'(i4,3x,30f14.6)') i_chan, radiance_d1 % bt (i_chan)
enddo

write(ioout2,*)
write(ioout2,*) 'Channel  clear Tb '

do i_chan = 1, nchannels
  write (ioout2,'(i4,3x,30f14.6)') i_chan, radiance_d1 % bt_clear (i_chan)
enddo


!* TANGENT-LINEAR TEST ***********************************************************************************

zeps = 0.01_JPRB

call allocate_profs( kflevg, kproma, profiles_tl, cld_profiles_tl, use_totalice, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d3, kflevg-1_jpim, iallocate)   
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl, kflevg-1_jpim, iallocate)  

call set_perturbation ( kproma, profiles_tl, cld_profiles_tl, profiles_d1, cld_profiles_d1, zeps, use_totalice)

emissivity_tl (1:nchannels) = emissivity_d1 (1:nchannels) * zeps 
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) < 0.01_JPRB

call rttov_scatt_tl ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_tl,        &! in
  & cld_profiles_tl,    &! in
  & emissivity_tl,      &! inout
  & radiance_d3,        &! inout
  & radiance_tl,        &! inout  
  & lnewcld = use_newcld, lusercfrac = use_cfrac)
      
! Save radiance as a reference for the trajectory
! TL is used instead of rttov_scatt because
! calcemis = F and reflectivities have not been saved
radiance_total_ref(:) = radiance_d1%bt(:)

! Full output
write(ioout,*)
write(ioout,*) 'Call to RTTOV_SCATT_TL'
write(ioout,*) '----------------------'
write(ioout,*)
write(ioout,*) 'Channel  cloudy Tb TL '

do i_chan = 1, nchannels
  write (ioout,'(i4,3x,30f24.16)') i_chan, radiance_tl % bt (i_chan)
enddo

write(ioout,*)
write(ioout,*) 'Channel  clear Tb TL '

do i_chan = 1, nchannels
  write (ioout,'(i4,3x,30f24.16)') i_chan, radiance_tl % bt_clear (i_chan)
enddo

! Partial output for comparison to reference output!
write(ioout2,*)
write(ioout2,*) 'Call to RTTOV_SCATT_TL'
write(ioout2,*) '----------------------'
write(ioout2,*)
write(ioout2,*) 'Channel  cloudy Tb TL '

do i_chan = 1, nchannels
  write (ioout2,'(i4,3x,30f14.6)') i_chan, radiance_tl % bt (i_chan)
enddo

write(ioout2,*)
write(ioout2,*) 'Channel  clear Tb TL '

do i_chan = 1, nchannels
  write (ioout2,'(i4,3x,30f14.6)') i_chan, radiance_tl % bt_clear (i_chan)
enddo


write(ioout,*)
write(ioout,*) 'Test TL'
write(ioout,*) '-------'
write(ioout,*)

!---------------------------
! second run of TL
!---------------------------
lambda = 0.5_JPRB
  
call allocate_profs( kflevg, kproma, profiles_tl2, cld_profiles_tl2, use_totalice, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl2, kflevg-1_jpim, iallocate)  

call set_perturbation ( kproma, profiles_tl2, cld_profiles_tl2, profiles_tl, cld_profiles_tl, lambda, use_totalice)
  
emissivity_tl2 (1:nchannels) = emissivity_tl (1:nchannels) *  lambda
calcemiss     (1:nchannels) = emissivity_tl (1:nchannels) < 0.01_JPRB

call rttov_scatt_tl ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_tl2,       &! in
  & cld_profiles_tl2,   &! in
  & emissivity_tl2,     &! inout
  & radiance_d1,        &! inout
  & radiance_tl2,       &! inout 
  & lnewcld = use_newcld, lusercfrac = use_cfrac)
  
!---------------------------

do i_chan = 1, nchannels
  if( abs(lambda * radiance_tl%bt_clear(i_chan) - radiance_tl2%bt_clear(i_chan)) > threshold ) &
    call error_stop( 'TL test fails for radiance_tl%bt_clear for channel ', i_chan )
  if( abs(lambda * radiance_tl%bt(i_chan) - radiance_tl2%bt(i_chan)) > threshold ) &
    call error_stop( 'TL test fails for radiance_tl%bt for channel ', i_chan )
end do
 

! Now run the Taylor test
!-------------------------

call allocate_profs( kflevg, kproma, profiles_d2, cld_profiles_d2, use_totalice, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d2, kflevg-1_jpim, iallocate)  

write(ioout,*) '  Until the RTTOV internal interpolation has TL/adjoint sensitivity '
write(ioout,*) ' to pressure changes, the following test will show very non-linear  '
write(ioout,*) ' behaviour from the operator, i.e (h*-h)/H very different from 1.'

do i_lambda = 10, 1, -1
  lambda = 10.0_JPRB ** (-1.0_JPRB * i_lambda) 

  errorstatus = errorstatus_success

  emissivity_d2 (1:nchannels) = emissivity_d1 (1:nchannels) + emissivity_tl (1:nchannels) * lambda
  calcemiss     (1:nchannels) = emissivity_d2 (1:nchannels) < 0.01_JPRB

  do i_proma = 1, kproma    
  
    !* Add perturbations
    profiles_d2 (i_proma) % p (:) = profiles_d1 (i_proma) % p (:) + profiles_tl (i_proma) % p (:) * lambda
    profiles_d2 (i_proma) % t (:) = profiles_d1 (i_proma) % t (:) + profiles_tl (i_proma) % t (:) * lambda 
    profiles_d2 (i_proma) % q (:) = profiles_d1 (i_proma) % q (:) + profiles_tl (i_proma) % q (:) * lambda 

    cld_profiles_d2 (i_proma) % cfrac   = cld_profiles_d1 (i_proma) % cfrac   + cld_profiles_tl (i_proma) % cfrac   * lambda
    cld_profiles_d2 (i_proma) % ph  (:) = cld_profiles_d1 (i_proma) % ph  (:) + cld_profiles_tl (i_proma) % ph  (:) * lambda
    cld_profiles_d2 (i_proma) % cc  (:) = cld_profiles_d1 (i_proma) % cc  (:) + cld_profiles_tl (i_proma) % cc  (:) * lambda 
    cld_profiles_d2 (i_proma) % clw (:) = cld_profiles_d1 (i_proma) % clw (:) + cld_profiles_tl (i_proma) % clw (:) * lambda 
    cld_profiles_d2 (i_proma) % rain (:) = cld_profiles_d1 (i_proma) % rain (:) + cld_profiles_tl (i_proma) % rain (:)*lambda 
    if( use_totalice) then 
      cld_profiles_d2 (i_proma) % totalice (:) = cld_profiles_d1 (i_proma) % totalice (:)  + &
                                               & cld_profiles_tl (i_proma) % totalice (:) * lambda 
    else
      cld_profiles_d2 (i_proma) % ciw (:) = cld_profiles_d1 (i_proma) % ciw (:)  + cld_profiles_tl (i_proma) % ciw(:)*lambda 
      cld_profiles_d2 (i_proma) % sp   (:) = cld_profiles_d1 (i_proma) % sp (:)  + cld_profiles_tl (i_proma) % sp (:)*lambda  
    endif

!* Fill in RTTOV/RTTOVSCATT arrays once per profile
    profiles_d2 (i_proma) % s2m % p = profiles_d1 (i_proma) % s2m % p + profiles_tl (i_proma) % s2m % p * lambda
    profiles_d2 (i_proma) % s2m % q = profiles_d1 (i_proma) % s2m % q + profiles_tl (i_proma) % s2m % q * lambda
    profiles_d2 (i_proma) % s2m % o = profiles_d1 (i_proma) % s2m % o + profiles_tl (i_proma) % s2m % o * lambda
    profiles_d2 (i_proma) % s2m % t = profiles_d1 (i_proma) % s2m % t + profiles_tl (i_proma) % s2m % t * lambda
    profiles_d2 (i_proma) % s2m % u = profiles_d1 (i_proma) % s2m % u + profiles_tl (i_proma) % s2m % u * lambda
    profiles_d2 (i_proma) % s2m % v = profiles_d1 (i_proma) % s2m % v + profiles_tl (i_proma) % s2m % v * lambda
    profiles_d2 (i_proma) % s2m % wfetc = 0.0_JPRB
    profiles_d2 (i_proma) % skin % surftype   = profiles_d1 (i_proma) % skin % surftype
    profiles_d2 (i_proma) % skin % t = profiles_d1 (i_proma) % skin % t + profiles_tl (i_proma) % skin % t * lambda 
    profiles_d2 (i_proma) % skin % fastem (:) = profiles_d1 (i_proma) % skin % fastem (:) + &
                                              & profiles_tl (i_proma) % skin % fastem (:) * lambda 
    profiles_d2 (i_proma) % skin % watertype = 0.0_JPRB  ! Ocean water   
    
    profiles_d2 (i_proma) % zenangle   = zenangle
    profiles_d2 (i_proma) % azangle    = 0.0_JPRB    ! default value
    profiles_d2 (i_proma) % ctp        = 500.0_JPRB  ! default value
    profiles_d2 (i_proma) % cfraction  =   0.0_JPRB  ! default value
    profiles_d2 (i_proma) % elevation  = 0.0_JPRB   ! default value
    profiles_d2 (i_proma) % sunzenangle  = 0.0_JPRB   ! default value
    profiles_d2 (i_proma) % sunazangle   = 0.0_JPRB   ! default value
    profiles_d2 (i_proma) % latitude     = 0.0_JPRB   ! default value
    profiles_d2 (i_proma) % longitude    = 0.0_JPRB   ! default value
    profiles_d2 (i_proma) % Be           = 0.0_JPRB   ! default value
    profiles_d2 (i_proma) % cosbk        = 0.0_JPRB   ! default value

  end do    

  !* Reference forward model run
  call rttov_scatt ( &
    & errorstatus,         &! out 
    & kflevg,              &! in
    & chanprof,            &! in
    & frequencies,         &! in
    & profiles_d2,         &! inout  
    & cld_profiles_d2,     &! in
    & coef_rttov,          &! in
    & coef_scatt,          &! in
    & calcemiss,           &! in
    & emissivity_d2,       &! inout
    & radiance_d2,         &! inout  
    & cfrac,               &! out, diagnostic only           
    & lnewcld = use_newcld, lusercfrac = use_cfrac)
     
  write(ioout,*)
  write(ioout,*) 'Chan      Lambda         Cloudy  [h(x*)-h(x)]/H(x*-x)  Clear '

  do i_chan = 1, nchannels
    ratio(1) = (radiance_d2 % bt(i_chan) - radiance_d1 % bt(i_chan)) / (lambda * radiance_tl % bt(i_chan))
    ratio(2) = (radiance_d2 % bt_clear(i_chan) - radiance_d1 % bt_clear(i_chan)) / (lambda * radiance_tl % bt_clear(i_chan))
    write (ioout,'(i4,3x,1e8.2,2f25.16)') i_chan,  lambda, ratio(1), ratio(2)    
  enddo
enddo
  
!* ADJOINT TEST ***********************************************************************************

write(ioout,*)
write(ioout,*) 'Test AD'
write(ioout,*) '-------'
write(ioout,*)

write(ioout,*) '1 - Test Linearity'
write(ioout,*)

write(ioout2,*)
write(ioout2,*) 'Test AD'
write(ioout2,*) '-------'
write(ioout2,*)

write(ioout2,*) '1 - Test Linearity'
write(ioout2,*)

call allocate_profs( kflevg, kproma, profiles_ad, cld_profiles_ad, use_totalice, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad, kflevg-1_jpim, iallocate)  

call zero_profs( kproma, profiles_ad, cld_profiles_ad, use_totalice) 

emissivity_ad (1:nchannels) = 0.0_JPRB 
   
! Set perturbations
!
radiance_ad % bt_clear(:) = 0.05_JPRB * radiance_d1 % bt_clear(:)
radiance_ad % bt(:)       = 0.05_JPRB * radiance_d1 % bt(:)
radiance_ad % clear     (:)   = 0.0_JPRB
radiance_ad % cloudy    (:)   = 0.0_JPRB
radiance_ad % total     (:)   = 0.0_JPRB
radiance_ad % upclear   (:)   = 0.0_JPRB
radiance_ad % dnclear   (:)   = 0.0_JPRB
radiance_ad % reflclear (:)   = 0.0_JPRB
radiance_ad % overcast  (:,:) = 0.0_JPRB
radiance_ad % up        (:,:) = 0.0_JPRB
radiance_ad % down      (:,:) = 0.0_JPRB
radiance_ad % surf      (:,:) = 0.0_JPRB
    
call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_ad,        &! in
  & cld_profiles_ad,    &! in
  & emissivity_ad,      &! inout
  & radiance_d2,        &! inout
  & radiance_ad,        &! inout  
  & lnewcld = use_newcld, lusercfrac = use_cfrac)
   
If ( errorstatus == errorstatus_fatal ) Then
  write ( ioout, '(A30,i4)' ) 'rttov_scatt_ad error'
  write ( ioout2, '(A30,i4)' ) 'rttov_scatt_ad error'
  Stop
End If

If ( Any( abs(radiance_total_ref(:) - radiance_d2%bt(:)) > threshold * radiance_total_ref(:)  ))  Then 
  write(default_err_unit,*) 'wrong forward model in AD'
  write(default_err_unit,*) radiance_total_ref(:)
  write(default_err_unit,*) abs(radiance_total_ref(:)-radiance_d2%bt(:)) / (threshold * radiance_total_ref(:))
  Stop
Endif

!---------------------------
! Second run of AD

call allocate_profs( kflevg, kproma, profiles_ad2, cld_profiles_ad2, use_totalice, iallocate)
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad2, kflevg-1_jpim, iallocate)  

call zero_profs( kproma, profiles_ad2, cld_profiles_ad2, use_totalice) 
   
emissivity_ad2 (1:nchannels) = 0.0_JPRB 
     
! Set perturbations
!
radiance_ad2 % bt_clear(:) = 0.05_JPRB * radiance_d1 % bt_clear(:) * lambda
radiance_ad2 % bt(:)       = 0.05_JPRB * radiance_d1 % bt(:) * lambda
radiance_ad2 % clear     (:)   = 0.0_JPRB
radiance_ad2 % cloudy    (:)   = 0.0_JPRB
radiance_ad2 % total     (:)   = 0.0_JPRB
radiance_ad2 % upclear   (:)   = 0.0_JPRB
radiance_ad2 % dnclear   (:)   = 0.0_JPRB
radiance_ad2 % reflclear (:)   = 0.0_JPRB
radiance_ad2 % overcast  (:,:) = 0.0_JPRB
radiance_ad2 % up        (:,:) = 0.0_JPRB
radiance_ad2 % down      (:,:) = 0.0_JPRB
radiance_ad2 % surf      (:,:) = 0.0_JPRB

call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_ad2,       &! in
  & cld_profiles_ad2,   &! in
  & emissivity_ad2,     &! inout
  & radiance_d2,        &! inout
  & radiance_ad2,       &! inout  
  & lnewcld = use_newcld, lusercfrac = use_cfrac)

If ( errorstatus == errorstatus_fatal ) Then
  write ( ioout, * ) 'rttov_scatt_ad error'
  write ( ioout2, * ) 'rttov_scatt_ad error'
  Stop
End If

do i_proma = 1, kproma
  do i_lev = 1, profiles_ad (i_proma) % nlevels
    if ( abs(lambda * profiles_ad (i_proma) % t (i_lev) - profiles_ad2 (i_proma) % t (i_lev)) > threshold ) &
      call error_stop( 'test AD 1 fails', i_lev ) 
    if ( abs(lambda * profiles_ad (i_proma) % q (i_lev) - profiles_ad2 (i_proma) % q (i_lev)) > threshold ) &
      call error_stop( 'test AD 2 fails', i_lev )
    if ( abs(lambda * profiles_ad (i_proma) % p (i_lev) - profiles_ad2 (i_proma) % p (i_lev)) > threshold ) &
       call error_stop( 'test AD 3 fails', i_lev )
  enddo
enddo


do i_proma = 1, kproma
  if ( use_cfrac .and. ((abs(lambda * cld_profiles_ad (i_proma) % cfrac - cld_profiles_ad2 (i_proma) % cfrac)) > threshold)) &
    call error_stop( 'test AD 3b fails', i_lev )
  do i_lev = 1, cld_profiles_ad (i_proma) % nlevels
    if ( abs(lambda * cld_profiles_ad (i_proma) % ph (i_lev) - cld_profiles_ad2 (i_proma) % ph (i_lev)) > threshold ) &
      call error_stop( 'test AD 4 fails', i_lev )
    if ( abs(lambda * cld_profiles_ad (i_proma) % cc (i_lev) - cld_profiles_ad2 (i_proma) % cc (i_lev)) > threshold ) &
      call error_stop( 'test AD 5 fails', i_lev )
    if ( abs(lambda * cld_profiles_ad (i_proma) % clw (i_lev) - cld_profiles_ad2 (i_proma) % clw (i_lev)) > threshold ) &
      call error_stop( 'test AD 6 fails', i_lev )
    if ( abs(lambda * cld_profiles_ad (i_proma) % rain (i_lev) - cld_profiles_ad2 (i_proma) % rain (i_lev))>threshold ) &
      call error_stop( 'test AD 8 fails', i_lev )
    if (use_totalice) then 
      if(abs(lambda * cld_profiles_ad (i_proma) % totalice(i_lev)-cld_profiles_ad2(i_proma)%totalice(i_lev))>threshold ) &
        call error_stop( 'test AD 7 fails', i_lev )
    else
      if ( abs(lambda * cld_profiles_ad (i_proma) % ciw (i_lev) - cld_profiles_ad2 (i_proma) % ciw (i_lev)) > threshold ) &
        call error_stop( 'test AD 7 fails', i_lev )
      if ( abs(lambda * cld_profiles_ad (i_proma) % sp (i_lev) - cld_profiles_ad2 (i_proma) % sp (i_lev)) > threshold ) &
        call error_stop( 'test AD 9 fails', i_lev )
    endif
  enddo
enddo

do i_proma = 1, kproma
  if ( abs(lambda * profiles_ad (i_proma) % s2m % t - profiles_ad2 (i_proma) % s2m % t) > threshold ) &
    call error_stop( 'test AD 12 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % q - profiles_ad2 (i_proma) % s2m % q) > threshold ) &
    call error_stop( 'test AD 13 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % p - profiles_ad2 (i_proma) % s2m % p) > threshold ) &
    call error_stop( 'test AD 14 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % u - profiles_ad2 (i_proma) % s2m % u) > threshold ) &
    call error_stop( 'test AD 15 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % s2m % v - profiles_ad2 (i_proma) % s2m % v) > threshold ) &
    call error_stop( 'test AD 16 fails', i_proma )
  if ( abs(lambda * profiles_ad (i_proma) % skin % t - profiles_ad2 (i_proma) % skin % t) > threshold ) &
    call error_stop( 'test AD 17 fails', i_proma ) 
enddo

do i_chan = 1, nchannels
  if ( abs(lambda * emissivity_ad (i_chan) - emissivity_ad2 (i_chan) ) > threshold ) &
    call error_stop( 'test AD 18 fails', i_chan )
enddo

write(ioout,*) '2 - Test Equality of Norms'
write(ioout,*)

write(ioout2,*) '2 - Test Equality of Norms'
write(ioout2,*)

call set_perturbation ( kproma, profiles_tl, cld_profiles_tl, profiles_d1, cld_profiles_d1, zeps, use_totalice) 
call zero_profs( kproma, profiles_ad, cld_profiles_ad, use_totalice) 
     
emissivity_d1 (1:nchannels) = 0.0_JPRB 
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) < 0.01_JPRB

emissivity_tl (1:nchannels) = emissivity_d1 (1:nchannels) * zeps 
emissivity_ad (1:nchannels) = 0.0_JPRB 

radiance_tl % bt_clear(:) = 0._JPRB
radiance_tl % bt(:)       = 0._JPRB

call rttov_scatt_tl ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_tl,        &! in
  & cld_profiles_tl,    &! in
  & emissivity_tl,      &! inout
  & radiance_d3,        &! inout
  & radiance_tl,        &! inout  
  & lnewcld = use_newcld, lusercfrac = use_cfrac)

If ( errorstatus == errorstatus_fatal ) Then
  write (ioout, * ) 'rttov_scatt_tl error'
  write (ioout2, * ) 'rttov_scatt_tl error'
  Stop
End If

!* compute <subtl(delta_x),delta_z>

zdelta1 = 0.0_JPRB

do i_chan = 1, nchannels
  zdelta1 = zdelta1 + (radiance_tl % bt(i_chan)) ** 2.0_JPRB 
  zdelta1 = zdelta1 + (radiance_tl % bt_clear(i_chan)) ** 2.0_JPRB  
enddo
   
!* Initialize     
radiance_ad % bt_clear(:) = radiance_tl % bt_clear (:) 
radiance_ad % bt(:)       = radiance_tl % bt (:) 
radiance_ad % clear     (:)   = 0.0_JPRB
radiance_ad % cloudy    (:)   = 0.0_JPRB
radiance_ad % total     (:)   = 0.0_JPRB
radiance_ad % upclear   (:)   = 0.0_JPRB
radiance_ad % dnclear   (:)   = 0.0_JPRB
radiance_ad % reflclear (:)   = 0.0_JPRB
radiance_ad % overcast  (:,:) = 0.0_JPRB
radiance_ad % up        (:,:) = 0.0_JPRB
radiance_ad % down      (:,:) = 0.0_JPRB
radiance_ad % surf      (:,:) = 0.0_JPRB

radiance_d1 % bt_clear(:) = 0.0_JPRB
radiance_d1 % bt(:)       = 0.0_JPRB
radiance_d1 % clear     (:)   = 0.0_JPRB
radiance_d1 % cloudy    (:)   = 0.0_JPRB
radiance_d1 % total     (:)   = 0.0_JPRB
radiance_d1 % upclear   (:)   = 0.0_JPRB
radiance_d1 % dnclear   (:)   = 0.0_JPRB
radiance_d1 % reflclear (:)   = 0.0_JPRB
radiance_d1 % overcast  (:,:) = 0.0_JPRB
radiance_d1 % up        (:,:) = 0.0_JPRB
radiance_d1 % down      (:,:) = 0.0_JPRB
radiance_d1 % surf      (:,:) = 0.0_JPRB   

!---------------------------
! Now run AD code with TL radiances in input
! move TL results to AD radiance increments
  
call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_ad,        &! in
  & cld_profiles_ad,    &! in
  & emissivity_ad,      &! inout
  & radiance_d2,        &! inout
  & radiance_ad,        &! inout  
  & lnewcld = use_newcld, lusercfrac = use_cfrac)

If ( errorstatus == errorstatus_fatal ) Then
  write ( ioout, * ) 'rttov_scatt_ad error'
  write ( ioout2, * ) 'rttov_scatt_ad error'
  Stop
End If

!* compute <delta_x,subad(delta_z)>

zdelta2 = 0.0_JPRB

do i_proma = 1, kproma

  if(use_cfrac) zdelta2 = zdelta2 &
    & + cld_profiles_tl (i_proma) % cfrac * cld_profiles_ad (i_proma) % cfrac 

  do i_lev = 1, kflevg
    zdelta2 = zdelta2 &
      & + cld_profiles_tl (i_proma) % cc   (i_lev) * cld_profiles_ad (i_proma) % cc   (i_lev) &
      & + cld_profiles_tl (i_proma) % clw  (i_lev) * cld_profiles_ad (i_proma) % clw  (i_lev) &
      & + cld_profiles_tl (i_proma) % rain (i_lev) * cld_profiles_ad (i_proma) % rain (i_lev) 
    if (use_totalice) then 
      zdelta2 = zdelta2 &
        & + cld_profiles_tl (i_proma) % totalice  (i_lev) * cld_profiles_ad (i_proma) % totalice  (i_lev) 
    else  
      zdelta2 = zdelta2 &
        & + cld_profiles_tl (i_proma) % ciw  (i_lev) * cld_profiles_ad (i_proma) % ciw  (i_lev) &
        & + cld_profiles_tl (i_proma) % sp   (i_lev) * cld_profiles_ad (i_proma) % sp   (i_lev) 
    endif
  enddo

  do i_lev = 1, kflevg + 1
    zdelta2 = zdelta2 &
      & + cld_profiles_tl (i_proma) % ph   (i_lev) * cld_profiles_ad (i_proma) % ph   (i_lev) 
  enddo

  do i_lev = 1, kflevg  
    zdelta2 = zdelta2  &
      & + profiles_tl (i_proma) % p   (i_lev) * profiles_ad (i_proma) % p   (i_lev) &
      & + profiles_tl (i_proma) % t   (i_lev) * profiles_ad (i_proma) % t   (i_lev) &
      & + profiles_tl (i_proma) % q   (i_lev) * profiles_ad (i_proma) % q   (i_lev)  
  enddo

  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % p * profiles_ad (i_proma) % s2m % p
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % q * profiles_ad (i_proma) % s2m % q
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % o * profiles_ad (i_proma) % s2m % o
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % t * profiles_ad (i_proma) % s2m % t
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % u * profiles_ad (i_proma) % s2m % u
  zdelta2 = zdelta2 + profiles_tl (i_proma) % s2m % v * profiles_ad (i_proma) % s2m % v

  zdelta2 = zdelta2 + profiles_tl (i_proma) % skin % t * profiles_ad (i_proma) % skin % t  

  do i_fast = 1, fastem_sp     
    zdelta2 = zdelta2 + profiles_tl (i_proma) % skin % fastem(i_fast)*profiles_ad(i_proma)%skin%fastem(i_fast)
  enddo         

  do i_chan = 1, nchannels 
    zdelta2 = zdelta2 + emissivity_tl (i_chan) * emissivity_ad (i_chan)
  enddo
enddo 

if (zdelta2 == 0._JPRB) then
  z = 1._JPRB
else
  z = zdelta2
endif
  
write (ioout,'(A10,f23.16)') 'delta1 = ', zdelta1
write (ioout,'(A10,f23.16)') 'delta2 = ', zdelta2

write (ioout,fmt= &
  & '('' The difference is '',f22.1, '' times the zero of the machine '')') &
  & abs((zdelta2-zdelta1)/epsilon(z)/z) 

if( abs((zdelta2-zdelta1)/epsilon(z)/z) < 50) then
  write (ioout,*) 'AD is OK'
  write (ioout2,*) 'AD is OK'
else
  call error_stop( 'Adjoint test fails', 0_jpim)
endif

!* K-TEST ***********************************************************************************

write(ioout,*)
write(ioout,*) 'Test K'
write(ioout,*) '------'
write(ioout,*)
write(ioout2,*)
write(ioout2,*) 'Test K'
write(ioout2,*) '------'
write(ioout2,*)

call allocate_profs( kflevg, nchannels, profiles_k, cld_profiles_k, use_totalice, iallocate)
call zero_profs( nchannels, profiles_k, cld_profiles_k, use_totalice) 
call rttov_alloc_rad ( erroralloc, nchannels, radiance_k, kflevg-1_jpim, iallocate)  

emissivity_d1 (1:nchannels) = 0.0_JPRB 
calcemiss     (1:nchannels) = emissivity_d1 (1:nchannels) < 0.01_JPRB      
emissivity_k  (1:nchannels) = 0.0_JPRB 

radiance_k % bt_clear(:) = 0.0_JPRB
radiance_k % bt(:)       = 1.0_JPRB ! (ACTIVATES K behaviour in adjoint code)
radiance_k % clear     (:)   = 0.0_JPRB
radiance_k % cloudy    (:)   = 0.0_JPRB
radiance_k % total     (:)   = 0.0_JPRB
radiance_k % upclear   (:)   = 0.0_JPRB
radiance_k % dnclear   (:)   = 0.0_JPRB
radiance_k % reflclear (:)   = 0.0_JPRB
radiance_k % overcast  (:,:) = 0.0_JPRB
radiance_k % up        (:,:) = 0.0_JPRB
radiance_k % down      (:,:) = 0.0_JPRB
radiance_k % surf      (:,:) = 0.0_JPRB   

call rttov_scatt_ad ( &
  & errorstatus,        &! out
  & kflevg,             &! in
  & chanprof,           &! in
  & frequencies,        &! in
  & profiles_d1,        &! inout  
  & cld_profiles_d1,    &! in
  & coef_rttov,         &! in
  & coef_scatt,         &! in
  & calcemiss,          &! in
  & emissivity_d1,      &! inout
  & profiles_k,         &! in (ACTIVATES K behaviour in adjoint code: note dimension)
  & cld_profiles_k,     &! in
  & emissivity_k,       &! inout
  & radiance_d1,        &! inout
  & radiance_k,         &! inout 
  & lnewcld = use_newcld, lusercfrac = use_cfrac)

If ( Any( abs(radiance_total_ref(:) - radiance_d1 % bt(:)) > threshold * radiance_total_ref(:) ))  Then
  write(ioout,*) 'wrong forward model in K'
  write(ioout,*) radiance_total_ref(:)
  write(ioout,*) radiance_d1 % bt(:)
  write(ioout,*) abs(radiance_total_ref(:)-radiance_d1 %bt(:)) / ( threshold * radiance_total_ref(:))
  write(ioout2,*) 'wrong forward model in K'
  write(ioout2,*) radiance_total_ref(:)
  write(ioout2,*) radiance_d1 % bt(:)
  write(ioout2,*) abs(radiance_total_ref(:)-radiance_d1 %bt(:)) / ( threshold * radiance_total_ref(:))
  Stop
Endif

!---------------------------
! Compares K to AD

do i_btout = 1, nchannels   

  radiance_ad % bt_clear(:) = 0.0_JPRB
  radiance_ad % bt(:)       = 0.0_JPRB
  radiance_ad % clear     (:)   = 0.0_JPRB
  radiance_ad % cloudy    (:)   = 0.0_JPRB
  radiance_ad % total     (:)   = 0.0_JPRB
  radiance_ad % upclear   (:)   = 0.0_JPRB
  radiance_ad % dnclear   (:)   = 0.0_JPRB
  radiance_ad % reflclear (:)   = 0.0_JPRB
  radiance_ad % overcast  (:,:) = 0.0_JPRB
  radiance_ad % up        (:,:) = 0.0_JPRB
  radiance_ad % down      (:,:) = 0.0_JPRB
  radiance_ad % surf      (:,:) = 0.0_JPRB   
  
  radiance_ad % bt (i_btout) = 1.0_JPRB
  
  call zero_profs( kproma, profiles_ad, cld_profiles_ad, use_totalice) 
              
  emissivity_ad (1:nchannels) = 0.0_JPRB    
      
  call rttov_scatt_ad ( &
    & errorstatus,    &! out
    & kflevg,             &! in
    & chanprof,           &! in
    & frequencies,        &! in
    & profiles_d1,        &! inout  
    & cld_profiles_d1,    &! in
    & coef_rttov,         &! in
    & coef_scatt,         &! in
    & calcemiss,          &! in
    & emissivity_d1,      &! inout
    & profiles_ad,        &! in
    & cld_profiles_ad,    &! in
    & emissivity_ad,      &! inout
    & radiance_d2,        &! inout
    & radiance_ad,        &! inout  
    & lnewcld = use_newcld, lusercfrac = use_cfrac)

  i_proma = chanprof(i_btout)%prof
                           
  do i_lev = 1, profiles_ad(i_proma) % nlevels
    if ( abs (profiles_ad (i_proma) % p   (i_lev) - profiles_k (i_btout) % p   (i_lev)) > threshold ) &
      call error_stop( 'test K 1 fails',i_lev)
    if ( abs (profiles_ad (i_proma) % t   (i_lev) - profiles_k (i_btout) % t   (i_lev)) > threshold ) &
      call error_stop( 'test K 2 fails',i_lev)
    if ( abs (profiles_ad (i_proma) % q   (i_lev) - profiles_k (i_btout) % q  (i_lev)) > threshold ) &
      call error_stop( 'test K 3 fails',i_lev)    
  End Do

  if( use_cfrac .and. (abs (cld_profiles_ad (i_proma) % cfrac - cld_profiles_k (i_btout) % cfrac) > threshold))&
     call error_stop( 'test K 4 fails',i_lev)  

  do i_lev = 1, cld_profiles_ad(i_proma) % nlevels
    if ( abs (cld_profiles_ad (i_proma) % ph   (i_lev) - cld_profiles_k (i_btout) % ph   (i_lev)) > threshold)&
      call error_stop( 'test K 6 fails',i_lev)
    if ( abs (cld_profiles_ad (i_proma) % cc   (i_lev) - cld_profiles_k (i_btout) % cc   (i_lev)) > threshold)&
      call error_stop( 'test K 8 fails',i_lev)
    if ( abs (cld_profiles_ad (i_proma) % clw   (i_lev) - cld_profiles_k (i_btout) % clw   (i_lev))>threshold)&
      call error_stop( 'test K 9 fails', i_lev)
    if ( abs (cld_profiles_ad (i_proma) % rain   (i_lev) - cld_profiles_k (i_btout) % rain   (i_lev))>threshold)&
      call error_stop( 'test K 11 fails',i_lev)
    if (use_totalice) then
      if (abs(cld_profiles_ad(i_proma)%totalice(i_lev)-cld_profiles_k(i_btout)%totalice(i_lev))>threshold)&
        call error_stop( 'test K 10 fails',i_lev)
    else
      if ( abs (cld_profiles_ad (i_proma) % ciw   (i_lev) - cld_profiles_k (i_btout) % ciw(i_lev))>threshold)&
        call error_stop( 'test K 10 fails',i_lev)
      if ( abs (cld_profiles_ad (i_proma) % sp   (i_lev) - cld_profiles_k (i_btout) % sp   (i_lev))>threshold)&
        call error_stop( 'test K 12 fails',i_lev)
    endif
  End Do
     
  if ( abs (profiles_ad (i_proma)  % s2m % p - profiles_k (i_btout) %  s2m % p) > threshold ) &
    call error_stop( 'test K 13 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % q - profiles_k (i_btout) %  s2m % q) > threshold ) &
    call error_stop( 'test K 14 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % o - profiles_k (i_btout) %  s2m % o) > threshold ) &
    call error_stop( 'test K 15 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % t - profiles_k (i_btout) %  s2m % t) > threshold ) &
    call error_stop( 'test K 16 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % u - profiles_k (i_btout) %  s2m % u) > threshold ) &
    call error_stop( 'test K 17 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % s2m % v - profiles_k (i_btout) %  s2m % v) > threshold ) &
    call error_stop( 'test K 18 fails',i_lev )
  if ( abs (profiles_ad (i_proma)  % skin % t - profiles_k (i_btout) %  skin % t) > threshold ) &
    call error_stop( 'test K 19 fails',i_lev )  
  if ( abs (emissivity_ad (i_btout) - emissivity_k (i_btout) ) > threshold ) &
    call error_stop('test K 23 fails',i_lev )
          
enddo

write(ioout,*) 'K is ok'
write(ioout,*)
write(ioout,*) 'End of RTTOVSCATT tests'
write(ioout,*)

write(ioout2,*) 'K is ok'
write(ioout2,*)
write(ioout2,*) 'End of RTTOVSCATT tests'
write(ioout2,*)

call allocate_profs( kflevg, kproma, profiles_d1, cld_profiles_d1, use_totalice, ideallocate)
call allocate_profs( kflevg, kproma, profiles_tl, cld_profiles_tl, use_totalice, ideallocate)
call allocate_profs( kflevg, kproma, profiles_tl2, cld_profiles_tl2, use_totalice, ideallocate)
call allocate_profs( kflevg, kproma, profiles_ad, cld_profiles_ad, use_totalice, ideallocate)
call allocate_profs( kflevg, kproma, profiles_d2, cld_profiles_d2, use_totalice, ideallocate)
call allocate_profs( kflevg, kproma, profiles_ad2, cld_profiles_ad2, use_totalice, ideallocate)
call allocate_profs( kflevg, nchannels, profiles_k, cld_profiles_k, use_totalice, ideallocate)

call rttov_alloc_rad ( erroralloc, nchannels, radiance_d1, kflevg-1_jpim, ideallocate)   
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d3, kflevg-1_jpim, ideallocate)   
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl, kflevg-1_jpim, ideallocate)  
call rttov_alloc_rad ( erroralloc, nchannels, radiance_tl2, kflevg-1_jpim, ideallocate)  
call rttov_alloc_rad ( erroralloc, nchannels, radiance_d2, kflevg-1_jpim, ideallocate)  
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad, kflevg-1_jpim, ideallocate)  
call rttov_alloc_rad ( erroralloc, nchannels, radiance_ad2, kflevg-1_jpim, ideallocate)  
call rttov_alloc_rad ( erroralloc, nchannels, radiance_k, kflevg-1_jpim, ideallocate)  

IF (LHOOK) CALL DR_HOOK('RTTOV_SCATT_TEST',1_jpim,ZHOOK_HANDLE)

!*******

contains

!*******

!* Allocate and prepare / deallocate profile structures
subroutine allocate_profs( kflevg, kproma, profiles, cld_profiles, use_totalice, asw)

Integer(Kind=jpim), Intent(in) :: kflevg       ! number of levels
Integer(Kind=jpim), Intent(in) :: kproma       ! number of profiles
Integer(Kind=jpim), Intent(in) :: asw          ! 1=allocate, 0=deallocate
Logical(Kind=jplm), Intent(in) :: use_totalice ! Choose separate ciw and snow, or totalice
Type(profile_type), Intent (inout)       :: profiles (kproma) 
Type(profile_cloud_type), Intent (inout) :: cld_profiles (kproma)

Type(rttov_options) :: opts ! Defaults to everything switched off
integer(kind=jpim)  :: err
   
call rttov_alloc_prof( err, kproma, profiles, kflevg, opts, asw, init = .true._jplm)
call rttov_alloc_scatt_prof ( kproma, cld_profiles, kflevg, use_totalice, asw, init = .true._jplm)

end subroutine

!*******

subroutine set_perturbation( kproma, profiles_tl, cld_profiles_tl, &
  & profiles_d, cld_profiles_d, zeps, use_totalice)

Integer(Kind=jpim), Intent(in) :: kproma       ! number of profiles (or channels in case of K/AD)
Type(profile_type),       Intent (inout) :: profiles_tl  (kproma)
Type(profile_cloud_type), Intent (inout) :: cld_profiles_tl (kproma)
Type(profile_type),       Intent (in)    :: profiles_d   (kproma)
Type(profile_cloud_type), Intent (in)    :: cld_profiles_d (kproma)
Real(Kind=jprb),          Intent (in)    :: zeps
Logical(Kind=jplm), Intent(in) :: use_totalice ! Choose separate ciw and snow, or totalice 

do i_proma = 1, kproma

  !* Set perturbation
  
  profiles_tl (i_proma) % p (:)        = profiles_d (i_proma) % p (:)        * zeps
  profiles_tl (i_proma) % t (1:kflevg) = profiles_d (i_proma) % t (1:kflevg) * zeps
  profiles_tl (i_proma) % q (1:kflevg) = profiles_d (i_proma) % q (1:kflevg) * zeps

  cld_profiles_tl (i_proma) % cfrac           = cld_profiles_d (i_proma) % cfrac  * zeps

  cld_profiles_tl (i_proma) % ph (:)          = cld_profiles_d (i_proma) % ph (:) * zeps
  cld_profiles_tl (i_proma) % cc   (1:kflevg) = cld_profiles_d (i_proma) % cc   (1:kflevg) * zeps
  cld_profiles_tl (i_proma) % clw  (1:kflevg) = cld_profiles_d (i_proma) % clw  (1:kflevg) * zeps
  cld_profiles_tl (i_proma) % rain (1:kflevg) = cld_profiles_d (i_proma) % rain (1:kflevg) * zeps
  if (use_totalice) then 
    cld_profiles_tl (i_proma) % totalice  (1:kflevg) = cld_profiles_d (i_proma) % totalice(1:kflevg)*zeps
  else
    cld_profiles_tl (i_proma) % ciw  (1:kflevg) = cld_profiles_d (i_proma) % ciw  (1:kflevg) * zeps
    cld_profiles_tl (i_proma) % sp   (1:kflevg) = cld_profiles_d (i_proma) % sp   (1:kflevg) * zeps
  endif 

    
  profiles_tl (i_proma) % s2m % p = profiles_d (i_proma) % s2m % p * zeps
  profiles_tl (i_proma) % s2m % q = profiles_d (i_proma) % s2m % q * zeps
  profiles_tl (i_proma) % s2m % o = profiles_d (i_proma) % s2m % o * zeps
  profiles_tl (i_proma) % s2m % t = profiles_d (i_proma) % s2m % t * zeps
  profiles_tl (i_proma) % s2m % u = profiles_d (i_proma) % s2m % u * zeps
  profiles_tl (i_proma) % s2m % v = profiles_d (i_proma) % s2m % v * zeps


  profiles_tl (i_proma) % s2m % wfetc = 0.0_JPRB ! Unused (solar computations only)

  profiles_tl (i_proma) % skin % surftype = -1 
  profiles_tl (i_proma) % skin % watertype = -1   

  profiles_tl (i_proma) % skin % t          = profiles_d (i_proma) % skin % t          * zeps
  profiles_tl (i_proma) % skin % fastem (:) = profiles_d (i_proma) % skin % fastem (:) * zeps

  profiles_tl (i_proma) % zenangle   = -1._jprb
  profiles_tl (i_proma) % azangle    = -1._jprb
  profiles_tl (i_proma) % ctp    = 0._jprb
  profiles_tl (i_proma) % cfraction    = 0._jprb
  profiles_tl (i_proma) % sunzenangle   = 0._jprb
  profiles_tl (i_proma) % sunazangle    = 0._jprb
  profiles_tl (i_proma) % latitude      = 0._jprb
  profiles_tl (i_proma) % longitude     = 0._jprb
  profiles_tl (i_proma) % elevation     = 0._jprb
  profiles_tl (i_proma) % Be            = 0._jprb
  profiles_tl (i_proma) % cosbk         = 0._jprb
  profiles_tl (i_proma) % idg    = 0._jprb
  profiles_tl (i_proma) % ish    = 0._jprb

enddo 

end subroutine

!*******

subroutine zero_profs( kproma, profiles_ad, cld_profiles_ad, use_totalice)

Integer(Kind=jpim), Intent(in) :: kproma       ! number of profiles (or channels in case of K/AD)
Type(profile_type),       Intent (inout) :: profiles_ad (kproma) 
Type(profile_cloud_type), Intent (inout) :: cld_profiles_ad (kproma)
Logical(Kind=jplm), Intent(in) :: use_totalice ! Choose separate ciw and snow, or totalice 

do i_proma = 1, kproma

  cld_profiles_ad (i_proma) % cfrac  = 0.0_JPRB  

  cld_profiles_ad (i_proma) % ph (:) = 0.0_JPRB  
  cld_profiles_ad (i_proma) % cc   (:) = 0.0_JPRB  
  cld_profiles_ad (i_proma) % clw  (:) = 0.0_JPRB  

  cld_profiles_ad (i_proma) % rain (:) = 0.0_JPRB  
  if (use_totalice) then 
    cld_profiles_ad (i_proma) % totalice (:) = 0.0_JPRB  
  else
    cld_profiles_ad (i_proma) % sp   (:) = 0.0_JPRB  
    cld_profiles_ad (i_proma) % ciw  (:) = 0.0_JPRB  
  endif
    
  profiles_ad (i_proma) % p   (:) = 0.0_JPRB
  profiles_ad (i_proma) % t   (:) = 0.0_JPRB 
  profiles_ad (i_proma) % q   (:) = 0.0_JPRB 

  profiles_ad (i_proma) % s2m % p = 0.0_JPRB 
  profiles_ad (i_proma) % s2m % q = 0.0_JPRB 
  profiles_ad (i_proma) % s2m % o = 0.0_JPRB 
  profiles_ad (i_proma) % s2m % t = 0.0_JPRB 
  profiles_ad (i_proma) % s2m % u = 0.0_JPRB 
  profiles_ad (i_proma) % s2m % v = 0.0_JPRB 
  profiles_ad (i_proma) % s2m % wfetc = 0.0_JPRB
  
  profiles_ad (i_proma) % skin % surftype  = -1 
  profiles_ad (i_proma) % skin % watertype = -1
  profiles_ad (i_proma) % skin % t          = 0.0_JPRB 
  profiles_ad (i_proma) % skin % fastem (:) = 0.0_JPRB 

  profiles_ad (i_proma) % zenangle   = -1
  profiles_ad (i_proma) % azangle    = -1
  profiles_ad (i_proma) % ctp        = 0.0_JPRB 
  profiles_ad (i_proma) % cfraction  = 0.0_JPRB 
  profiles_ad (i_proma) % sunzenangle   = 0._jprb
  profiles_ad (i_proma) % sunazangle    = 0._jprb
  profiles_ad (i_proma) % latitude      = 0._jprb
  profiles_ad (i_proma) % longitude     = 0._jprb
  profiles_ad (i_proma) % elevation     = 0._jprb
  profiles_ad (i_proma) % Be            = 0._jprb
  profiles_ad (i_proma) % cosbk         = 0._jprb
  profiles_ad (i_proma) % idg    = 0._jprb
  profiles_ad (i_proma) % ish    = 0._jprb

enddo 

end subroutine

!*******
subroutine error_stop( text, inum)

character (len=*), intent(in)   :: text
integer( kind=jpim), intent(in) :: inum

write(default_err_unit,*) ' -------------------------------------------------- '
write(default_err_unit,*) ' ------------------- TEST HAS FAILED -------------- '
write(default_err_unit,*)
write(default_err_unit,*) text, inum
write(default_err_unit,*)
write(default_err_unit,*) ' -------------------------------------------------- '
write(default_err_unit,*) ' -------------------------------------------------- '

stop

end subroutine
       
End subroutine rttovscatt_test_one
