subroutine rttov_scatt_setupindex (&
  & nprofiles,  &  ! number of profiles
  & n_chan,     &  ! number of channels 
  & coef,       &  ! coef structure read in from rttov coef file
  & nchannels,  &  ! number of calculated channels
  & chanprof,   &  ! channels and profile numbers
  & frequencies, & ! array, frequency number for each "channel"
  & lchannel_subset) ! OPTIONAL array of logical flags to indicate a subset of channels
                           
! See RTTOV User guide for definitions of all these variables.

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

! P. Bauer, P. Lopez, E. Moreau, D. Salmond   ECMWF    May 2004

! 1. Set up indices for profiles/channels

! RTTOV_SCATT_SETUPINDEX is called from ONEDVAR_OBSOP_RTTOV, ONEDVAR_OBSOP_RTTOV_GRAD

! Modifications:

! 10/10/2006 Chris O'Dell : Made to work for any sensor (removed SSM/I hardcoding).
!                         : Changed the way lsprofiles2 is assigned to be slightly 
!                         : more general. Changed order of passed variables for 
!                         : style consistency.
! 11/11/2007 Alan Geer    : RTTOV9 version, cleaned.
! 06/06/2008 Alan Geer    : Added optional channel subsetting 
! 26/02/2009 W. Bell      : Changes for Windsat.

!* KIND     
use parkind1,    only : jpim, jprb, jplm
use rttov_const, only : npolar_compute
use rttov_types, only : rttov_coef, rttov_chanprof   

implicit none

integer (kind=jpim), intent ( in) :: nprofiles
integer (kind=jpim), intent ( in) :: n_chan 
type   (rttov_coef), intent ( in) :: coef 
integer (kind=jpim), intent ( in) :: nchannels
logical (kind=jplm), optional, intent ( in) :: lchannel_subset(nprofiles, n_chan)
 
integer  (kind=jpim), intent (out), dimension (nchannels) :: frequencies
type(rttov_chanprof), Intent (out), dimension (nchannels) :: chanprof ! Channel and profile indices

!INTF_END

integer (kind=jpim) :: i_prof, i_chan, j_chan, i_freq, polid
integer (kind=jpim) :: old_polid
real (kind=jprb)    :: cwn, old_cwn
logical (kind=jplm) :: luse(nprofiles, n_chan)
!- End of header ------------------------------------------------------

luse(:,:) = .true.
if (present(lchannel_subset)) luse = lchannel_subset

!* Set index arrays
j_chan = 0  ! counter to store calculated channels
old_cwn = 0.0_jprb
old_polid = -1

do i_prof = 1, nprofiles
  i_freq = 0
  do i_chan = 1, n_chan
   
    polid = coef % fastem_polar (i_chan) + 1 ! polarisation ID of this channel
    cwn = coef % ff_cwn(i_chan)

    ! Below is a test to see if this channel represents a new frequency.
    ! The following 3 conditions must all hold in order to be same frequency as 
    ! previous channel:
    !   1) Same central wave number as previous channel
    !      (NB check the absolute diff because some SSMI/S channel pairs have 
    !          slightly different central wavenumbers: 1E-3 cm-1 == 0.03 GHz)
    !   2) Polarization ID eq 4,5,6 or 7 (single V or H polarisation only)
    !   3) Polarization ID different from last channel
    if ( (abs(cwn - old_cwn) > 1.E-3_jprb) .or. &
      & ((polid /= 4) .and. (polid /= 5) .and. (polid /= 6) .and. (polid /= 7)) &
      & .or. (polid == old_polid) ) then  

      i_freq = i_freq + 1

    endif

    if( luse(i_prof,i_chan) ) then

      ! Profile and frequency number for each calculated channel
      j_chan = j_chan + 1
      chanprof   (j_chan)%chan  = i_chan
      frequencies(j_chan)       = i_freq
      chanprof   (j_chan)%prof  = i_prof

    endif

    old_polid = polid
    old_cwn = cwn

  end do
end do 

end subroutine rttov_scatt_setupindex
