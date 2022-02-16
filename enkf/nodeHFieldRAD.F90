module NodeHFieldRAD_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl
   use interpolate, only : PosGrid, convert_LatLon_to_GridPos, convert_GridPos_offset, &
                         & FlagOverTop, interpolation2D
   use interpolate_TLAD, only : TL_interpolation2D, AD_interpolation2D
   use parkind1, only : jpim, jprb
   use rttovlib
   use NodeInfo_class
   use NodeObsSpaceFieldRAD_class
   use NodeObsValidField_class
   use NodeObsField_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeProfileControl_class
   use NodeMPI
   implicit none
   !
   integer(4), parameter :: uslevels = 47
   real(kind=r_size), parameter :: KGKG2PPMV = 1.60771704d6
   !real(kind=r_size), dimension(uslevels) :: us_p, us_t, us_q, us_o3
   real(kind=r_size), dimension(uslevels) :: us_p, us_t, us_o3
   !real(kind=r_size), dimension(uslevels) :: ref_p, ref_t, ref_q
   real(kind=r_size), dimension(uslevels) :: ref_q
   integer(4) :: ilev
   !
   DATA (us_p(ilev),ilev=1,uslevels)/ &
   0.005000E+00,0.013100E+00,0.030400E+00,0.064400E+00, &
   0.100000E+00,0.290000E+00,0.690000E+00,0.142000E+01,0.261100E+01, &
   0.440700E+01,0.695000E+01,0.103700E+02,0.148100E+02,0.204000E+02, &
   0.272600E+02,0.355100E+02,0.452900E+02,0.567300E+02,0.699700E+02, &
   0.851800E+02,0.102050E+03,0.122040E+03,0.143840E+03,0.167950E+03, &
   0.194360E+03,0.222940E+03,0.253710E+03,0.286600E+03,0.321500E+03, &
   0.358280E+03,0.396810E+03,0.436950E+03,0.478540E+03,0.521460E+03, &
   0.565540E+03,0.610600E+03,0.656430E+03,0.702730E+03,0.749120E+03, &
   0.795090E+03,0.839950E+03,0.882800E+03,0.922460E+03,0.957440E+03, &
   0.985880E+03,0.100543E+04,0.101325E+04 /
   !
   DATA (us_t(ilev),ilev=1,uslevels)/ &
   0.191313E+03,0.197579E+03,0.208936E+03,0.224429E+03, &
   0.231635E+03,0.252872E+03,0.269142E+03,0.265346E+03,0.252376E+03, &
   0.241730E+03,0.232520E+03,0.227413E+03,0.225093E+03,0.222998E+03, &
   0.221165E+03,0.219444E+03,0.217888E+03,0.216647E+03,0.216718E+03, &
   0.216697E+03,0.216701E+03,0.216701E+03,0.216694E+03,0.216720E+03, &
   0.216690E+03,0.216482E+03,0.221171E+03,0.226731E+03,0.231547E+03, &
   0.236427E+03,0.241075E+03,0.245536E+03,0.249835E+03,0.253959E+03, &
   0.257916E+03,0.261712E+03,0.265345E+03,0.268812E+03,0.272104E+03, &
   0.275206E+03,0.278094E+03,0.280738E+03,0.283098E+03,0.285114E+03, &
   0.286711E+03,0.287788E+03,0.288214E+03 /
   !
!  DATA (us_q(ilev),ilev=1,uslevels)/ &
!  0.256269E-02,0.306150E-02,0.323828E-02,0.325424E-02,0.314645E-02, &
!  0.307367E-02,0.302166E-02,0.296545E-02,0.288478E-02,0.281744E-02, &
!  0.271834E-02,0.260025E-02,0.248594E-02,0.241979E-02,0.238108E-02, &
!  0.239867E-02,0.242665E-02,0.313469E-02,0.383665E-02,0.710243E-02, &
!  0.119370E-01,0.208427E-01,0.357383E-01,0.648740E-01,0.128563E+00, &
!  0.232985E+00,0.321762E+00,0.436822E+00,0.601353E+00,0.779304E+00, &
!  0.100684E+01,0.130068E+01,0.162530E+01,0.199245E+01,0.242870E+01, &
!  0.288133E+01,0.328199E+01,0.364070E+01,0.398304E+01,0.430064E+01, &
!  0.456475E+01,0.474702E+01,0.481971E+01 /
   !
   DATA (us_o3(ilev),ilev=1,uslevels)/ &
   0.308392E+00,0.308392E+00,0.308392E+00,0.308392E+00, &
   0.639447E+00,0.133264E+01,0.272510E+01,0.508174E+01,0.701138E+01, &
   0.783251E+01,0.761789E+01,0.685889E+01,0.609873E+01,0.557377E+01, &
   0.490937E+01,0.409623E+01,0.318304E+01,0.249272E+01,0.179687E+01, &
   0.127884E+01,0.894568E+00,0.642060E+00,0.490098E+00,0.377517E+00, &
   0.308549E+00,0.226650E+00,0.149429E+00,0.109116E+00,0.804815E-01, &
   0.590528E-01,0.519773E-01,0.458820E-01,0.405801E-01,0.384802E-01, &
   0.363139E-01,0.340682E-01,0.332325E-01,0.331802E-01,0.330332E-01, &
   0.323680E-01,0.311295E-01,0.297832E-01,0.286691E-01,0.278203E-01, &
   0.271840E-01,0.267617E-01,0.265947E-01 /
   !
   DATA (ref_q(ilev),ilev=1,uslevels)/ &
   0.616179E-02,0.576365E-02,0.541568E-02,0.510539E-02, &
   0.487837E-02,0.449743E-02,0.417236E-02,0.383499E-02,0.364302E-02, &
   0.326617E-02,0.291474E-02,0.276716E-02,0.264687E-02,0.255156E-02, &
   0.248366E-02,0.244131E-02,0.240914E-02,0.238102E-02,0.235978E-02, &
   0.237089E-02,0.272652E-02,0.380070E-02,0.611720E-02,0.106726E-01, &
   0.194080E-01,0.386951E-01,0.775005E-01,0.142208E+00,0.251103E+00, &
   0.403920E+00,0.615826E+00,0.895425E+00,0.125042E+01,0.167965E+01, &
   0.217088E+01,0.271999E+01,0.333378E+01,0.399480E+01,0.468579E+01, &
   0.541720E+01,0.618229E+01,0.696163E+01,0.774010E+01,0.851887E+01, &
   0.958083E+01,0.940562E+01,0.933204E+01 /
   !
   type NodeHFieldRAD
      character(len=10) :: name
      integer :: nobs
      type(PosGrid), dimension(:), allocatable :: pg
      ! For ensemble linear tangent H
      integer :: ne
      integer, dimension(:), allocatable :: nrank
      real(r_size), dimension(:,:), allocatable :: eigval
      real(r_size), dimension(:,:,:), allocatable :: Y
      real(r_size), dimension(:,:,:,:,:), allocatable :: X, eigvec
   end type NodeHFieldRAD
   !
   interface new
      module procedure new_NodeHFieldRAD
   end interface
   interface destroy
      module procedure destroy_NodeHFieldRAD
   end interface
   interface display
      module procedure display_NodeHFieldRAD
   end interface
   interface get_name
      module procedure get_name_NodeHFieldRAD
   end interface
   interface get_nobs
      module procedure get_nobs_NodeHFieldRAD
   end interface
   interface get_xyloc
      module procedure get_xyloc_NodeHFieldRAD1
      module procedure get_xyloc_NodeHFieldRAD2
   end interface
!   interface apply_Hlogp
!      module procedure apply_Hlogp_NodeHFieldRAD
!   end interface
   interface apply_H
      module procedure apply_H_NodeHFieldRAD1
      module procedure apply_H_NodeHFieldRAD2
      module procedure apply_H_NodeHFieldRAD3
   end interface
   interface initialize_DH
      module procedure initialize_DH_NodeHFieldRAD
   end interface
   interface apply_DH
      module procedure apply_DH_NodeHFieldRAD
   end interface
   interface apply_DHT
      module procedure apply_DHT_NodeHFieldRAD
   end interface
   !
contains
   !
   subroutine new_NodeHFieldRAD(self, info, obsspace)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      if (self%nobs > 0) then
         allocate(self%pg(self%nobs))
         do iobs = 1, self%nobs
            call convert_LatLon_to_GridPos(self%pg(iobs), obsspace%lat(iobs), obsspace%lon(iobs), 1)
	    call convert_GridPos_offset(self%pg(iobs), is-1, js-1)
         end do
      end if
      !
      return
   end subroutine new_NodeHFieldRAD
   !
   !
   !
   subroutine destroy_NodeHFieldRAD(self)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      !
      if (allocated(self%pg)) deallocate(self%pg)
      if (allocated(self%nrank)) then
	 deallocate(self%Y, self%X, self%nrank, self%eigval, self%eigvec)
      end if
      !
      return
   end subroutine destroy_NodeHFieldRAD
   !
   !
   !
   subroutine display_NodeHFieldRAD(self)
      implicit none
      type(NodeHFieldRAD), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs
      !
      return
   end subroutine display_NodeHFieldRAD
   !
   !
   !
   subroutine get_name_NodeHFieldRAD(self, name)
      implicit none
      type(NodeHFieldRAD), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeHFieldRAD
   !
   !
   !
   subroutine get_nobs_NodeHFieldRAD(self, nobs)
      implicit none
      type(NodeHFieldRAD), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeHFieldRAD
   !
   !
   !
   subroutine get_xyloc_NodeHFieldRAD1(self, xyloc)
      implicit none
      type(NodeHFieldRAD), intent(in) :: self
      type(NodeObsField), intent(inout) :: xyloc
      integer :: iobs
      !
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px
	 xyloc%field(iobs,2) = self%pg(iobs)%py
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldRAD1
   !
   !
   !
   subroutine get_xyloc_NodeHFieldRAD2(self, info, xyloc)
      implicit none
      type(NodeHFieldRAD), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsField), intent(inout) :: xyloc
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px - is + 1
	 xyloc%field(iobs,2) = self%pg(iobs)%py - js + 1
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldRAD2
   !
   !
   !
   subroutine extend_profile(glevels, nlevels, p, t, q, o3, p_n, t_n, q_n, o3_n)
      implicit none
      !
      integer(4) ,intent(in) :: glevels, nlevels
      real(kind=r_size), dimension(glevels), intent(in) :: p, t, q, o3
      real(kind=r_size), dimension(nlevels), intent(out) :: p_n, t_n, q_n, o3_n
      integer(4) :: jlev, klev, ntop
      real(kind=r_size) :: lrate(uslevels)
      !
      do jlev = 1, uslevels
         if (p(2) < us_p(jlev)) then
            ntop = jlev
            exit
         end if
      end do
      do jlev = 2, glevels, 1
         p_n(jlev+nlevels-glevels) = p(jlev)
         t_n(jlev+nlevels-glevels) = t(jlev)
         q_n(jlev+nlevels-glevels) = q(jlev)
         o3_n(jlev+nlevels-glevels) = o3(jlev)
      end do
      !
      lrate(ntop) = (us_t(ntop-1)-us_t(ntop))/(us_p(ntop)-us_p(ntop-1))
      t_n( nlevels-glevels+1) = t(2) + lrate(ntop)*(p(2)-us_p(ntop-1)) ! T is extrapolated with us standard T laps rate.
      p_n( nlevels-glevels+1) = us_p( ntop-1) ! p is defined in rt levels
      q_n( nlevels-glevels+1) = ref_q(ntop-1)
      o3_n(nlevels-glevels+1) = us_o3(ntop-1)
      do klev = 2, ntop-1
         if (nlevels-glevels-klev+2<1) exit
         t_n( nlevels-glevels-klev+2) = t_n(nlevels-glevels+1-klev+2) + us_t(ntop-klev)-us_t(ntop-klev+1)
         p_n( nlevels-glevels-klev+2) = us_p(ntop-klev)
         q_n( nlevels-glevels-klev+2) = ref_q(ntop-klev)
         o3_n(nlevels-glevels-klev+2) = us_o3(ntop-klev)
      enddo
      !
      return
   end subroutine extend_profile
   !
   !
   !
   subroutine extend_Dprofile(glevels, nlevels, p, p_tl, t_tl, q_tl, o3_tl, p_n_tl, t_n_tl, q_n_tl, o3_n_tl)
      implicit none
      !
      integer(4) ,intent(in) :: glevels, nlevels
      real(kind=r_size), dimension(glevels), intent(in) :: p, p_tl, t_tl, q_tl, o3_tl
      real(kind=r_size), dimension(nlevels), intent(out) :: p_n_tl, t_n_tl, q_n_tl, o3_n_tl
      integer(4) :: jlev, klev, ntop
      real(kind=r_size) :: lrate(uslevels)
      !
      do jlev = 1, uslevels
         if (p(2) < us_p(jlev)) then
            ntop = jlev
            exit
         end if
      end do
      do jlev = 2, glevels,1
         p_n_tl(jlev+nlevels-glevels) = p_tl(jlev)
         t_n_tl(jlev+nlevels-glevels) = t_tl(jlev)
         q_n_tl(jlev+nlevels-glevels) = q_tl(jlev)
         o3_n_tl(jlev+nlevels-glevels) = o3_tl(jlev)
      end do
      !
      lrate(ntop) = (us_t(ntop-1)-us_t(ntop))/(us_p(ntop)-us_p(ntop-1))
      t_n_tl( nlevels-glevels+1) = t_tl(2) + lrate(ntop)*(p_tl(2)) ! T is extrapolated with us standard T laps rate.
      p_n_tl( nlevels-glevels+1) = 0.0d0
      q_n_tl( nlevels-glevels+1) = 0.0d0
      o3_n_tl(nlevels-glevels+1) = 0.0d0 
      do klev = 2, ntop-1
         if (nlevels-glevels-klev+2 < 1) exit
         t_n_tl( nlevels-glevels-klev+2) = t_n_tl(nlevels-glevels+1-klev+2)
         p_n_tl( nlevels-glevels-klev+2) = 0.0d0
         q_n_tl( nlevels-glevels-klev+2) = 0.0d0
         o3_n_tl(nlevels-glevels-klev+2) = 0.0d0
      enddo
      !
      return
   end subroutine extend_Dprofile
   !
   !
   !
   subroutine extend_DprofileT(glevels, nlevels, p, p_n_ad, t_n_ad, q_n_ad, o3_n_ad, p_ad, t_ad, q_ad, o3_ad)
      implicit none
      !
      integer(4) ,intent(in) :: glevels, nlevels
      real(kind=r_size), dimension(glevels), intent(in) :: p
      real(kind=r_size), dimension(nlevels), intent(inout) :: p_n_ad, t_n_ad, q_n_ad, o3_n_ad
      real(kind=r_size),dimension(glevels),intent(out) :: p_ad, t_ad, q_ad, o3_ad
      integer(4) :: jlev, klev, ntop
      real(kind=r_size) :: lrate(uslevels)
      !
      do jlev = 1, uslevels
         if (p(2) < us_p(jlev)) then
            ntop = jlev
            exit
         end if
      end do
      lrate(ntop) = (us_t(ntop-1)-us_t(ntop))/(us_p(ntop)-us_p(ntop-1))
      !
      p_ad(:)=0.0d0
      t_ad(:)=0.0d0
      q_ad(:)=0.0d0
      o3_ad(:)=0.0d0
      do klev = 2, ntop-1
         if (nlevels-glevels-klev+2<1) exit
         t_n_ad(nlevels-glevels+1-klev+2) = t_n_ad(nlevels-glevels-klev+2)
      enddo
      p_ad(2) = p_ad(2) + lrate(ntop)*t_n_ad(nlevels-glevels+1)
      t_ad(2) = t_ad(2) +             t_n_ad(nlevels-glevels+1)
      !
      do jlev = 2, glevels,1
         p_ad(jlev) = p_ad(jlev) + p_n_ad(jlev+nlevels-glevels)
         t_ad(jlev) = t_ad(jlev) + t_n_ad(jlev+nlevels-glevels)
         q_ad(jlev) = q_ad(jlev) + q_n_ad(jlev+nlevels-glevels)
         o3_ad(jlev) = o3_ad(jlev) + o3_n_ad(jlev+nlevels-glevels)
      end do
      !
      return
   end subroutine extend_DprofileT
   !
   !
   !
   subroutine apply_rttov(iobs, nlevels, nchannels, obsspace, usrf, vsrf, pprf_rt, tprf_rt, qprf_rt, ierrorstatus, y)
      use rttov_const, only :         &
     & errorstatus_success,           &
     & errorstatus_fatal,             &
     & nplatforms,ninst,sensor_id_mw, &
     & inst_id_goesim,inst_id_gmsim,  &
     & platform_name,inst_name,       &
     & mh2o, mair
      use rttov_types, only : &
     & rttov_coefs,           &
     & rttov_chanprof,        &
     & profile_Type,          &
     & transmission_Type,     &
     & radiance_Type
      implicit none
      !
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      integer, intent(in) :: iobs, nlevels
      integer(kind=jpim), intent(in) :: nchannels               ! Number of radiances computed (channels used * pro
      real(kind=r_size), intent(in) :: usrf, vsrf
      real(kind=r_size), intent(in), dimension(nlevels) :: pprf_rt, tprf_rt, qprf_rt
      integer, intent(out) :: ierrorstatus
      type(NodeObsField), intent(inout) :: y
      !
      integer(kind=jpim), parameter :: nprofiles = 1                ! Number of profiles
      type(rttov_chanprof), allocatable :: chanprof(:)
      type(profile_Type) :: profiles(nprofiles)
      type(rttov_coefs), pointer :: coefs
      !type(rttov_coef_scatt), pointer :: coef_scatt
      logical :: calcemis(nchannels)
      !
      ! Forward model outputs
      type(transmission_Type) :: transmission
      type(radiance_Type) :: radiancedata
      real(kind=jprb) :: emissivity(nchannels)
      real(kind=jprb) :: emissivity_out(nchannels)
      integer(kind=jpim) :: errorstatus  ! rttov error return code
      integer(kind=jpim), allocatable :: setup_errorstatus(:) ! setup return code
      !
      ! set-up
      integer, parameter :: ERR_UNIT=6
      integer, parameter :: VERBOSITY_LEVEL=1 ! 0: no error messages output
                                              ! 1: fatal errors only printed
                                              ! 2: warning erros only printed
                                              ! 3: information message
      integer:: num_sensor      ! number of the sensors used
      integer:: r_satn          ! satn : satellite number
      integer:: r_nidx          ! nidx : number of channel index
      character(len=8) :: r_plat ! plat : platform
      character(len=8) :: r_inst ! inst : instrument
      integer :: tvsinst(3), r_cidx(10)
      integer, allocatable :: tvschan(:)
      integer :: ich,jch
      integer :: i_coef, c_plat, c_inst, c_chan, ctchan, ierr
      integer(kind=jpim) :: alloc_status(60)
      integer :: i, iz
      logical :: hyperch
      integer :: kchanidx
      !
#include "rttov_alloc_prof.h"
#include "rttov_errorhandling.h"
#include "rttov_direct.h"
#include "rttov_alloc_rad.h"
      !
      if (rttv_chidx(jpchus+1,obsspace%satid(iobs)) == 0) then   !non.hyperCH
         hyperch = .false.
      else if(rttv_chidx(jpchus+1,obsspace%satid(iobs)) /= 0) then  !hyperCH
         hyperch = .true.
      end if
      errorstatus = 0
      coefs => rt%coefs(obsspace%satid(iobs))
      !
      call  rttov_alloc_prof(errorstatus, nprofiles, profiles, nlevels, rt%opts, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Profiles Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata, nlevels-1_jpim, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata Memory allocation error ", errorstatus
         stop 99
      endif
      !
      profiles(1)%nlevels = nlevels
      profiles(1)%p(1:profiles(1)%nlevels) = pprf_rt(1:nlevels)
      profiles(1)%t(1:profiles(1)%nlevels) = tprf_rt(1:nlevels)
      profiles(1)%q(1:profiles(1)%nlevels) = qprf_rt(1:nlevels)
      profiles(1)%s2m%u = usrf !obsspace%us(iobs)
      profiles(1)%s2m%v = vsrf !obsspace%vs(iobs)
      profiles(1)%s2m%p = profiles(1)%p(nlevels) !obsspace%ps(iobs)
      profiles(1)%s2m%t = profiles(1)%t(nlevels) !obsspace%ts(iobs)
      profiles(1)%s2m%q = profiles(1)%q(nlevels)
      profiles(1)%s2m%wfetc = 100000
      !
      profiles(1)%zenangle = obsspace%sat_zenith(iobs)
      profiles(1)%azangle  = obsspace%sat_azimuth(iobs)
      !profiles(1)%azangle = 0.0d0
      profiles(1)%sunzenangle = 0.0d0
      profiles(1)%sunazangle  = 0.0d0
      profiles(1)%elevation   = 0.0d0
      profiles(1)%latitude    = 0.0d0
      profiles(1)%ctp         = 500.0d0
      profiles(1)%cfraction   =   0.0d0
      profiles(1)%skin%t      = obsspace%sst(iobs)
      profiles(1)%skin%surftype = obsspace%surftype(iobs)
      profiles(1)%skin%watertype = 1   ! 0=freshwaterm, 1=ocean water
      profiles(1)%skin%fastem(:) = 0.d0  ! If these are used, RTTOV will stop!
      if(.not. associated(coefs)) write(6,*) ".not. associated( coefs )"
      !
      calcemis(1:nchannels) = .true.
      emissivity(:) = 0.d0
      emissivity_out(:) = 0.d0
      if (profiles(1)%skin%surftype == 0 .or. profiles(1)%skin%surftype ==2 ) then
         emissivity(:) = 0.9d0
         calcemis(:)= .false.
      end if
      allocate(transmission%tau_levels(profiles(1)%nlevels,nchannels))
      allocate(transmission%tau_total(nchannels))
      transmission%tau_levels = 0._jprb
      transmission%tau_total = 0._jprb
      allocate(chanprof(nchannels))
      do i = 1, nchannels
         chanprof(i)%chan = i
         chanprof(i)%prof = 1
      enddo
      profiles(1)%idg = 0._jprb
      profiles(1)%ish = 0._jprb
      !
      errorstatus=0
      call rttov_direct(errorstatus, chanprof, rt%opts, profiles, coefs, &
                      & calcemis, emissivity, emissivity_out, transmission, radiancedata)    
      ierrorstatus = 0
      if (errorstatus > 1) then
         write(6,*) "@@@@@@ This is errorstatus from rttov @@@@@@", errorstatus 
         ierrorstatus = errorstatus
      endif 
      if (any(emissivity_out >= 1.0d0))then
         write(6,*) "@@@@@@ This is error rttov emissivity estimation @@@@@@"
         ierrorstatus = errorstatus_fatal
      endif 
      !
      do ich = 1, obsspace%nchannel(iobs)
         if (hyperch) then  !hyperCH
            kchanidx = rttv_ch2378tidx(obsspace%channel(iobs,ich),obsspace%satid(iobs))
         else
            kchanidx = obsspace%channel(iobs,ich)
         end if
         if (kchanidx <= 0) cycle ! process only ch defined in parm_rttv
         jch = kchanidx
         y%field(iobs,ich) = radiancedata%bt(jch)
      enddo
      !
      deallocate(chanprof)
      deallocate(transmission%tau_levels)
      deallocate(transmission%tau_total)
      call rttov_alloc_prof(errorstatus, nprofiles, profiles, nlevels, rt%opts, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Profiles Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata, nlevels-1_jpim, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata Memory deallocation error ", errorstatus
         stop 99
      endif
      !
      return
   end subroutine apply_rttov
   !
   !
   !
   subroutine apply_Drttov(iobs, nlevels, nchannels, obsspace, usrf9, vsrf9, usrf, vsrf, &
                         & pprf_rt9, tprf_rt9, qprf_rt9, pprf_rt, tprf_rt, qprf_rt, ierrorstatus, y)
      use rttov_const, only :         &
     & errorstatus_success,           &
     & errorstatus_fatal,             &
     & nplatforms,ninst,sensor_id_mw, &
     & inst_id_goesim,inst_id_gmsim,  &
     & platform_name,inst_name,       &
     & mh2o, mair
      use rttov_types, only : &
     & rttov_coefs,           &
     & rttov_chanprof,        &
     & profile_Type,          &
     & transmission_Type,     &
     & radiance_Type
      implicit none
      !
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      integer, intent(in) :: iobs, nlevels
      integer(kind=jpim), intent(in) :: nchannels               ! Number of radiances computed (channels used * pro
      real(kind=r_size), intent(in) :: usrf9, vsrf9, usrf, vsrf
      real(kind=r_size), intent(in), dimension(nlevels) :: pprf_rt9, tprf_rt9, qprf_rt9
      real(kind=r_size), intent(in), dimension(nlevels) :: pprf_rt, tprf_rt, qprf_rt
      integer, intent(out) :: ierrorstatus
      type(NodeObsField), intent(inout) :: y
      !
      integer(kind=jpim), parameter :: nprofiles = 1                ! Number of profiles
      type(rttov_chanprof), allocatable :: chanprof(:)
      type(profile_Type) :: profiles(nprofiles)
      type(rttov_coefs), pointer :: coefs
      !type(rttov_coef_scatt), pointer :: coef_scatt
      logical :: calcemis(nchannels)
      !
      ! Forward model outputs
      type(transmission_Type) :: transmission
      type(radiance_Type) :: radiancedata
      type(transmission_Type) :: transmission_tl
      type(radiance_Type) :: radiancedata_tl
      real(kind=jprb) :: emissivity(nchannels)
      real(kind=jprb) :: emissivity_out(nchannels)
      ! TL variables for rttov_tl calls
      type(profile_Type) :: profiles_tl(nprofiles)
      real(Kind=jprb) :: emissivity_tl(nchannels)
      real(Kind=jprb) :: emissivity_out_tl(nchannels)
      !
      integer(Kind=jpim) :: errorstatus  ! rttov error return code
      integer(Kind=jpim), allocatable :: setup_errorstatus(:) ! setup return code
      integer, parameter :: ERR_UNIT=6
      integer, parameter :: VERBOSITY_LEVEL=1 ! 0: no error messages output
                                              ! 1: fatal errors only printed
                                              ! 2: warning erros only printed
                                              ! 3: information message
      integer:: num_sensor      ! number of the sensors used
      integer:: r_satn          ! satn : satellite number
      integer:: r_nidx          ! nidx : number of channel index
      character(len=8):: r_plat ! plat : platform
      character(len=8):: r_inst ! inst : instrument
      integer:: tvsinst(3),r_cidx(10)
      integer,allocatable:: tvschan(:)
      integer :: ich,jch
      integer :: i_coef,c_plat,c_inst,c_chan,ctchan,ierr
      integer(Kind=jpim) :: alloc_status(60)
      integer :: i,iz
      logical :: hyperch
      integer :: kchanidx
      !
#include "rttov_alloc_prof.h"
#include "rttov_errorhandling.h"
#include "rttov_tl.h"
#include "rttov_alloc_rad.h"
      !
      if (rttv_chidx(jpchus+1,obsspace%satid(iobs)) == 0) then !non.hyperCH
         hyperch = .false.
      else if (rttv_chidx(jpchus+1,obsspace%satid(iobs)) /= 0) then  !hyperCH
         hyperch = .true.
      end if
      errorstatus=0
      coefs => rt%coefs(obsspace%satid(iobs))
      !
      rt%opts%addinterp = .true.
      rt%opts%addsolar   = .false.
      rt%opts%addaerosl  = .false.
      rt%opts%addclouds  = .false.
      rt%opts%switchrad=.true.
      rt%opts%spacetop=.true.
      rt%opts%lgradp=.false.
      rt%opts%use_q2m   = .false.
      rt%opts%apply_reg_limits    = .true.
      rt%opts%verbose_checkinput_warnings    = .true.
      rt%opts%ozone_data    = .false.
      rt%opts%co2_data      = .false.
      rt%opts%n2o_data      = .false.
      rt%opts%co_data       = .false.
      rt%opts%ch4_data      = .false.
      rt%opts%clw_data      = .true.
      rt%opts%addrefrac     = .false.
      rt%opts%do_checkinput = .true.
      call rttov_alloc_prof(errorstatus, nprofiles, profiles, nlevels, rt%opts, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Profiles Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_prof(errorstatus, nprofiles, profiles_tl, nlevels, rt%opts, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Profiles_tl Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata, nlevels-1_jpim, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata_tl, nlevels-1_jpim, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata_tl Memory allocation error ", errorstatus
         stop 99
      endif
      !
      profiles(1)%nlevels=nlevels
      profiles(1)%p(1:profiles(1)%nlevels) = pprf_rt9(1:nlevels)
      profiles(1)%t(1:profiles(1)%nlevels) = tprf_rt9(1:nlevels)
      profiles(1)%q(1:profiles(1)%nlevels) = qprf_rt9(1:nlevels)
      profiles(1)%s2m%u = usrf9 !obsspace%us(iobs)
      profiles(1)%s2m%v = vsrf9 !obsspace%vs(iobs) 
      profiles(1)%s2m%p = profiles(1)%p(nlevels) !obsspace%ps(iobs)
      profiles(1)%s2m%t = profiles(1)%t(nlevels) !obsspace%ts(iobs)
      profiles(1)%s2m%q = profiles(1)%q(nlevels)
      profiles(1)%s2m%wfetc = 100000
      profiles(1)%zenangle = obsspace%sat_zenith(iobs)
      profiles(1)%azangle  = obsspace%sat_azimuth(iobs)
      !profiles(1)%azangle = 0.0d0
      profiles(1)%sunzenangle = 0.0d0
      profiles(1)%sunazangle  = 0.0d0
      profiles(1)%elevation   = 0.0d0
      profiles(1)%latitude    = 0.0d0
      profiles(1)%ctp         = 500.0d0
      profiles(1)%cfraction   =   0.0d0
      profiles(1)%skin%t = obsspace%sst(iobs)
      profiles(1)%skin%surftype = obsspace%surftype(iobs)
      profiles(1)%skin%watertype = 1   ! 0=freshwaterm, 1=ocean water
      profiles(1)%skin%fastem(:) = 0.d0  ! If these are used, RTTOV will stop!
      if(.not. associated( coefs ) ) write(6,*) ".not. associated( coefs )"
      !
      calcemis(1:nchannels) = .true.
      emissivity(:) = 0.d0
      emissivity_out(:) = 0.d0
      if (profiles(1)%skin%surftype == 0 .or. profiles(1)%skin%surftype == 2) then
         emissivity(:) = 0.9d0
         calcemis(:)= .false.
      end if
      emissivity_tl(:) = 0.d0
      emissivity_out_tl(:) = 0.d0
      allocate(transmission%tau_levels(profiles(1)%nlevels,nchannels))
      allocate(transmission%tau_total(nchannels))
      allocate(transmission_tl%tau_levels(profiles(1)%nlevels,nchannels))
      allocate(transmission_tl%tau_total(nchannels))
      transmission%tau_levels = 0._jprb
      transmission%tau_total = 0._jprb
      transmission_tl%tau_levels = 0._jprb
      transmission_tl%tau_total = 0._jprb
      !
      allocate(chanprof(nchannels))
      do i = 1,nchannels
         chanprof(i)%chan = i
         chanprof(i)%prof = 1
      enddo
      profiles(1)%idg = 0._jprb
      profiles(1)%ish = 0._jprb
      ! initialise TL output profile variables
      do i = 1, nprofiles
         profiles_tl(i)%zenangle   = -1       ! no meaning
         profiles_tl(i)%azangle    = -1       ! no meaning
         profiles_tl(i)%sunzenangle= -1       ! no meaning
         profiles_tl(i)%sunazangle = -1       ! no meaning
         profiles_tl(i)%latitude   = -1       ! no meaning
         profiles_tl(i)%elevation  = -1       ! no meaning
         profiles_tl(i)%p(:)   = 0._JPRB ! pressure levels
         profiles_tl(i)%t(:)   = 0._JPRB ! temperarure
         if (rt%opts%ozone_data) profiles_tl(i) % o3(:)  = 0._JPRB ! O3 ppmv
         if (rt%opts%co2_data) profiles_tl(i) % co2(:) = 0._JPRB ! CO2 ppmv
         if (rt%opts%n2o_data) profiles_tl(i) % n2o(:) = 0._JPRB ! N2O ppmv
         if (rt%opts%co_data)  profiles_tl(i) % co(:)  = 0._JPRB ! CO ppmv
         if (rt%opts%ch4_data) profiles_tl(i) % ch4(:) = 0._JPRB ! CH4 ppmv
         if (rt%opts%addaerosl) profiles_tl(i) % aerosols(:,:) = 0._JPRB 
         if (rt%opts%addclouds)then
            profiles_tl(i) % cloud(:,:)=0._JPRB
            profiles_tl(i) % cfrac(:,:)=0._JPRB
         endif
         profiles_tl(i) % clw(:) = 0._JPRB ! clw
	 profiles_tl(i)%nlevels = nlevels
	 profiles_tl(i)%p(1:profiles_tl(i)%nlevels) = pprf_rt(1:nlevels)
         profiles_tl(i)%t(1:profiles_tl(i)%nlevels) = tprf_rt(1:nlevels)
         profiles_tl(i)%q(1:profiles_tl(i)%nlevels) = qprf_rt(1:nlevels)
         profiles_tl(i) % s2m % u = usrf! wind components
         profiles_tl(i) % s2m % v = vsrf! wind components
	 profiles_tl(i) % s2m % p = profiles_tl(i)%p(nlevels) ! pressure
         profiles_tl(i) % s2m % t = profiles_tl(i)%t(nlevels) ! temperarure
         profiles_tl(i) % s2m % q = profiles_tl(i)%q(nlevels) ! WV
	 profiles_tl(i) % s2m % o = 0._JPRB! O3
	 profiles_tl(i) % s2m % wfetc = 0._JPRB! wind fetc 
         profiles_tl(i) % skin % surftype = -1  ! no meaning
         profiles_tl(i) % skin % t        = 0._JPRB  ! on temperarure
         profiles_tl(i) % skin % fastem   = 0._JPRB  ! Fastem
         profiles_tl(i) % ctp       = 0._JPRB  ! cloud top pressure
         profiles_tl(i) % cfraction = 0._JPRB  ! cloud fraction
      end do
      !
      radiancedata_tl % clear(:)              = 0._JPRB
      radiancedata_tl % cloudy(:)             = 0._JPRB
      radiancedata_tl % total(:)              = 0._JPRB
      radiancedata_tl % bt(:)                 = 0._JPRB
      radiancedata_tl % bt_clear(:)           = 0._JPRB
      radiancedata_tl % upclear(:)            = 0._JPRB
      radiancedata_tl % reflclear(:)          = 0._JPRB
      radiancedata_tl % overcast(:,:)         = 0._JPRB
      radiancedata_tl%bt(:)=0.0d0
      do ich = 1, obsspace%nchannel(iobs)
         if (hyperch) then  !hyperCH
            kchanidx = rttv_ch2378tidx(obsspace%channel(iobs,ich),obsspace%satid(iobs))
         else
            kchanidx = obsspace%channel(iobs,ich)
         end if
         if (kchanidx <= 0) cycle ! process only ch defined in parm_rttv
         jch = kchanidx
      end do
      !
      errorstatus = 0
      call rttov_tl(errorstatus, chanprof, rt%opts, profiles, profiles_tl, coefs, calcemis, &
                  & emissivity, emissivity_tl, emissivity_out, emissivity_out_tl, &
                  & transmission, transmission_tl, radiancedata, radiancedata_tl)    
      ierrorstatus = 0
      if (errorstatus > 1) then
         write(6,*) "@@@@@@ This is errorstatus from rttov_tl @@@@@@", errorstatus 
         ierrorstatus = errorstatus
      endif 
      if (any(emissivity_out >= 1.0d0)) then
         write(6,*) "@@@@@@ This is error rttov emissivity estimation @@@@@@"
         ierrorstatus = errorstatus_fatal
      endif 
      do ich = 1, obsspace%nchannel(iobs)
         if (hyperch) then  !hyperCH
            kchanidx = rttv_ch2378tidx(obsspace%channel(iobs,ich),obsspace%satid(iobs))
         else
            kchanidx = obsspace%channel(iobs,ich)
         end if
         if (kchanidx <= 0) cycle !process only ch defined in parm_rttv
         jch = kchanidx
	 y%field(iobs,ich) = radiancedata_tl%bt(jch) 
      enddo
      !
      deallocate(chanprof)
      deallocate( transmission % tau_levels )
      deallocate( transmission % tau_total )
      deallocate( transmission_tl % tau_levels )
      deallocate( transmission_tl % tau_total )
      call rttov_alloc_prof(errorstatus, nprofiles, profiles, nlevels, rt%opts, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Profiles Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_prof(errorstatus, nprofiles, profiles_tl, nlevels, rt%opts, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Profiles_tl Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata, nlevels-1_jpim, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata_tl, nlevels-1_jpim, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata_tl Memory deallocation error ", errorstatus
         stop 99
      endif
      !
      return
   end subroutine apply_Drttov
   !
   !
   !
   subroutine apply_DrttovT(iobs, nlevels, nchannels, obsspace, usrf9, vsrf9, &
                          & pprf_rt9, tprf_rt9, qprf_rt9, y, ierrorstatus, usrf, vsrf, pprf_rt, tprf_rt, qprf_rt)
      use rttov_const, only :         &
     & errorstatus_success,           &
     & errorstatus_fatal,             &
     & nplatforms,ninst,sensor_id_mw, &
     & inst_id_goesim,inst_id_gmsim,  &
     & platform_name,inst_name,       &
     & mh2o, mair
      use rttov_types, only : &
     & rttov_coefs,           &
     & rttov_chanprof,        &
     & profile_Type,          &
     & transmission_Type,     &
     & radiance_Type      
      implicit none
      !
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      integer, intent(in) :: iobs, nlevels
      integer(kind=jpim), intent(in) :: nchannels               ! Number of radiances computed (channels used * pro
      type(NodeObsField), intent(in) :: y
      real(kind=r_size), intent(in) :: usrf9, vsrf9
      real(kind=r_size), intent(in), dimension(nlevels) :: pprf_rt9, tprf_rt9, qprf_rt9
      integer, intent(out) :: ierrorstatus
      real(kind=r_size), intent(out) :: usrf, vsrf
      real(kind=r_size), intent(out), dimension(nlevels) :: pprf_rt, tprf_rt, qprf_rt
      !
      integer(kind=jpim), parameter :: nprofiles = 1                ! Number of profiles
      type(rttov_chanprof), allocatable :: chanprof(:)
      type(profile_Type) :: profiles(nprofiles)
      type(rttov_coefs), pointer :: coefs
      !type(rttov_coef_scatt), pointer :: coef_scatt
      logical :: calcemis(nchannels)
      !
      ! Forward model outputs
      type(transmission_Type) :: transmission
      type(radiance_Type) :: radiancedata
      type(transmission_Type) :: transmission_ad
      type(radiance_Type) :: radiancedata_ad
      real(kind=jprb) :: emissivity(nchannels)
      real(kind=jprb) :: emissivity_out(nchannels)
      ! AD variables for rttov_ad calls
      type(profile_Type) :: profiles_ad(nprofiles)
      real(Kind=jprb) :: emissivity_ad(nchannels)
      real(Kind=jprb) :: emissivity_out_ad(nchannels)
      !
      integer(Kind=jpim) :: errorstatus  ! rttov error return code
      integer(Kind=jpim), allocatable :: setup_errorstatus(:) ! setup return code
      integer, parameter :: ERR_UNIT=6
      integer, parameter :: VERBOSITY_LEVEL=1 ! 0: no error messages output
                                              ! 1: fatal errors only printed
                                              ! 2: warning erros only printed
                                              ! 3: information message
      integer:: num_sensor      ! number of the sensors used
      integer:: r_satn          ! satn : satellite number
      integer:: r_nidx          ! nidx : number of channel index
      character(len=8):: r_plat ! plat : platform
      character(len=8):: r_inst ! inst : instrument
      integer:: tvsinst(3),r_cidx(10)
      integer,allocatable:: tvschan(:)
      integer :: ich,jch
      integer :: i_coef,c_plat,c_inst,c_chan,ctchan,ierr
      integer(Kind=jpim) :: alloc_status(60)
      integer :: i,iz
      logical :: hyperch
      integer :: kchanidx
      !
#include "rttov_alloc_prof.h"
#include "rttov_errorhandling.h"
#include "rttov_ad.h"
#include "rttov_alloc_rad.h"
      !
      if (rttv_chidx(jpchus+1,obsspace%satid(iobs)) == 0) then !non.hyperCH
         hyperch = .false.
      else if (rttv_chidx(jpchus+1,obsspace%satid(iobs)) /= 0) then  !hyperCH
         hyperch = .true.
      end if
      errorstatus=0
      coefs => rt%coefs(obsspace%satid(iobs))
      !
      rt%opts%addinterp = .true.
      rt%opts%addsolar   = .false.
      rt%opts%addaerosl  = .false.
      rt%opts%addclouds  = .false.
      rt%opts%switchrad=.true.
      rt%opts%spacetop=.true.
      rt%opts%lgradp=.false.
      rt%opts%use_q2m   = .false.
      rt%opts%apply_reg_limits    = .true.
      rt%opts%verbose_checkinput_warnings    = .true.
      rt%opts%ozone_data    = .false.
      rt%opts%co2_data      = .false.
      rt%opts%n2o_data      = .false.
      rt%opts%co_data       = .false.
      rt%opts%ch4_data      = .false.
      rt%opts%clw_data      = .true.
      rt%opts%addrefrac     = .false.
      rt%opts%do_checkinput = .true.
      call rttov_alloc_prof(errorstatus, nprofiles, profiles, nlevels, rt%opts, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Profiles Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_prof(errorstatus, nprofiles, profiles_ad, nlevels, rt%opts, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Profiles_ad Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata, nlevels-1_jpim, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata Memory allocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata_ad, nlevels-1_jpim, 1_jpim, init=.true.)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata_ad Memory allocation error ", errorstatus
         stop 99
      endif
      !
      profiles(1)%nlevels=nlevels
      profiles(1)%p(1:profiles(1)%nlevels) = pprf_rt9(1:nlevels)
      profiles(1)%t(1:profiles(1)%nlevels) = tprf_rt9(1:nlevels)
      profiles(1)%q(1:profiles(1)%nlevels) = qprf_rt9(1:nlevels)
      profiles(1)%s2m%u = usrf9 !obsspace%us(iobs)
      profiles(1)%s2m%v = vsrf9 !obsspace%vs(iobs) 
      profiles(1)%s2m%p = profiles(1)%p(nlevels) !obsspace%ps(iobs)
      profiles(1)%s2m%t = profiles(1)%t(nlevels) !obsspace%ts(iobs)
      profiles(1)%s2m%q = profiles(1)%q(nlevels) 
      profiles(1)%s2m%wfetc = 100000
      profiles(1)%zenangle = obsspace%sat_zenith(iobs)
      profiles(1)%azangle  = obsspace%sat_azimuth(iobs)
      !profiles(1)%azangle = 0.0d0
      profiles(1)%sunzenangle = 0.0d0
      profiles(1)%sunazangle  = 0.0d0
      profiles(1)%elevation   = 0.0d0
      profiles(1)%latitude    = 0.0d0
      profiles(1)%ctp         = 500.0d0
      profiles(1)%cfraction   =   0.0d0
      profiles(1)%skin%t = obsspace%sst(iobs)
      profiles(1)%skin%surftype = obsspace%surftype(iobs)
      profiles(1)%skin%watertype = 1   ! 0=freshwaterm, 1=ocean water
      profiles(1)%skin%fastem(:) = 0.d0  ! If these are used, RTTOV will stop!
      if(.not. associated( coefs ) ) write(6,*) ".not. associated( coefs )"
      !
      calcemis(1:nchannels) = .true.
      emissivity(:) = 0.d0
      emissivity_out(:) = 0.d0
      if (profiles(1)%skin%surftype == 0 .or. profiles(1)%skin%surftype == 2) then
         emissivity(:) = 0.9d0
         calcemis(:)= .false.
      end if
      emissivity_ad(:) = 0.d0
      emissivity_out_ad(:) = 0.d0
      allocate(transmission%tau_levels(profiles(1)%nlevels,nchannels))
      allocate(transmission%tau_total(nchannels))
      allocate(transmission_ad%tau_levels(profiles(1)%nlevels,nchannels))
      allocate(transmission_ad%tau_total(nchannels))
      transmission%tau_levels = 0._jprb
      transmission%tau_total = 0._jprb
      transmission_ad%tau_levels = 0._jprb
      transmission_ad%tau_total = 0._jprb
      !
      allocate(chanprof(nchannels))
      do i = 1,nchannels
         chanprof(i)%chan = i
         chanprof(i)%prof = 1
      enddo
      profiles(1)%idg = 0._jprb
      profiles(1)%ish = 0._jprb
      ! initialise AD output profile variables
      do i = 1, nprofiles
         profiles_ad(i)%zenangle   = -1       ! no meaning
         profiles_ad(i)%azangle    = -1       ! no meaning
         profiles_ad(i)%sunzenangle= -1       ! no meaning
         profiles_ad(i)%sunazangle = -1       ! no meaning
         profiles_ad(i)%latitude   = -1       ! no meaning
         profiles_ad(i)%elevation  = -1       ! no meaning
         if (rt%opts%ozone_data) profiles_ad(i) % o3(:)  = 0._JPRB ! O3 ppmv
         if (rt%opts%co2_data) profiles_ad(i) % co2(:) = 0._JPRB ! CO2 ppmv
         if (rt%opts%n2o_data) profiles_ad(i) % n2o(:) = 0._JPRB ! N2O ppmv
         if (rt%opts%co_data)  profiles_ad(i) % co(:)  = 0._JPRB ! CO ppmv
         if (rt%opts%ch4_data) profiles_ad(i) % ch4(:) = 0._JPRB ! CH4 ppmv
         if (rt%opts%addaerosl) profiles_ad(i) % aerosols(:,:) = 0._JPRB 
         if (rt%opts%addclouds)then
            profiles_ad(i) % cloud(:,:)=0._JPRB
            profiles_ad(i) % cfrac(:,:)=0._JPRB
         endif
         profiles_ad(i) % clw(:) = 0._JPRB ! clw
	 profiles_ad(i)%nlevels = nlevels
	 profiles_ad(i) % p(:)   = 0._JPRB ! pressure levels
         profiles_ad(i) % t(:)   = 0._JPRB ! temperarure
	 profiles_ad(i) % q(:)   = 0._JPRB ! WV 
         profiles_ad(i) % s2m % t = 0._JPRB! temperarure
         profiles_ad(i) % s2m % q = 0._JPRB! WV
         profiles_ad(i) % s2m % o = 0._JPRB! O3
         profiles_ad(i) % s2m % p = 0._JPRB! pressure
         profiles_ad(i) % s2m % u = 0._JPRB! wind components
         profiles_ad(i) % s2m % v = 0._JPRB! wind components
         profiles_ad(i) % s2m % wfetc = 0._JPRB! wind fetc 
         profiles_ad(i) % skin % surftype = -1  ! no meaning
         profiles_ad(i) % skin % t        = 0._JPRB  ! on temperarure
         profiles_ad(i) % skin % fastem   = 0._JPRB  ! Fastem
         profiles_ad(i) % ctp       = 0._JPRB  ! cloud top pressure
         profiles_ad(i) % cfraction = 0._JPRB  ! cloud fraction
      end do
      !
      radiancedata_ad % clear(:)              = 0._JPRB
      radiancedata_ad % cloudy(:)             = 0._JPRB
      radiancedata_ad % total(:)              = 0._JPRB
      radiancedata_ad % bt(:)                 = 0._JPRB
      radiancedata_ad % bt_clear(:)           = 0._JPRB
      radiancedata_ad % upclear(:)            = 0._JPRB
      radiancedata_ad % reflclear(:)          = 0._JPRB
      radiancedata_ad % overcast(:,:)         = 0._JPRB
      radiancedata_ad%bt(:)=0.0d0
      do ich = 1, obsspace%nchannel(iobs)
         if (hyperch) then  !hyperCH
            kchanidx = rttv_ch2378tidx(obsspace%channel(iobs,ich),obsspace%satid(iobs))
         else
            kchanidx = obsspace%channel(iobs,ich)
         end if
         if (kchanidx <= 0) cycle ! process only ch defined in parm_rttv
         jch = kchanidx
         radiancedata_ad%bt(jch) = radiancedata_ad%bt(jch) + y%field(iobs,ich)
      end do
      !
      errorstatus = 0
      call rttov_ad(errorstatus, chanprof, rt%opts, profiles, profiles_ad, coefs, calcemis, &
                  & emissivity, emissivity_ad, emissivity_out, emissivity_out_ad, &
                  & transmission, transmission_ad, radiancedata, radiancedata_ad)    
      ierrorstatus = 0
      if (errorstatus > 1) then
         write(6,*) "@@@@@@ This is errorstatus from rttov_ad @@@@@@", errorstatus 
         ierrorstatus = errorstatus
      endif 
      if (any(emissivity_out >= 1.0d0)) then
         write(6,*) "@@@@@@ This is error rttov emissivity estimation @@@@@@"
         ierrorstatus = errorstatus_fatal
      endif 
      do ich = 1, obsspace%nchannel(iobs)
         if (hyperch) then  !hyperCH
            kchanidx = rttv_ch2378tidx(obsspace%channel(iobs,ich),obsspace%satid(iobs))
         else
            kchanidx = obsspace%channel(iobs,ich)
         end if
         if (kchanidx <= 0) cycle !process only ch defined in parm_rttv
         jch = kchanidx
      enddo
      !
      pprf_rt(1:nlevels) = 0.0d0
      tprf_rt(1:nlevels) = 0.0d0
      qprf_rt(1:nlevels) = 0.0d0
      pprf_rt(1:nlevels) = pprf_rt(1:nlevels) + profiles_ad(1)%p(1:nlevels)
      tprf_rt(1:nlevels) = tprf_rt(1:nlevels) + profiles_ad(1)%t(1:nlevels)
      qprf_rt(1:nlevels) = qprf_rt(1:nlevels) + profiles_ad(1)%q(1:nlevels)
      pprf_rt(nlevels) = pprf_rt(nlevels) + profiles_ad(1)%s2m%p
      tprf_rt(nlevels) = tprf_rt(nlevels) + profiles_ad(1)%s2m%t
      qprf_rt(nlevels) = qprf_rt(nlevels) + profiles_ad(1)%s2m%q
      usrf = 0.0d0
      usrf = usrf + profiles_ad(1)%s2m%u
      vsrf = 0.0d0
      vsrf = vsrf + profiles_ad(1)%s2m%v
      !
      deallocate(chanprof)
      deallocate( transmission % tau_levels )
      deallocate( transmission % tau_total )
      deallocate( transmission_ad % tau_levels )
      deallocate( transmission_ad % tau_total )
      call rttov_alloc_prof(errorstatus, nprofiles, profiles, nlevels, rt%opts, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Profiles Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_prof(errorstatus, nprofiles, profiles_ad, nlevels, rt%opts, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Profiles_ad Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata, nlevels-1_jpim, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata Memory deallocation error ", errorstatus
         stop 99
      endif
      call rttov_alloc_rad(errorstatus, nchannels, radiancedata_ad, nlevels-1_jpim, 0_jpim)
      if (errorstatus /= 0) then
         write(6,*) " Radiancedata_ad Memory deallocation error ", errorstatus
         stop 99
      endif
      !
      return
   end subroutine apply_DrttovT
   !
   !
   !
   subroutine interpolate_rad1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(x%nx,x%ny,x%nlev) :: u, v, p, t, qv
      !
      integer :: nze, ierror
      integer(kind=jpim) :: nchannel
      real(r_size) :: us, vs
      real(r_size), dimension(x%nlev) :: p1D, t1D, qv1D, o31D
      real(r_size), dimension(x%nlev+10) :: p1De, t1De, qv1De, o31De
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      nze = nz + 10
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(x, 'u', it, u, nz)
	 call get_field(x, 'v', it, v, nz)
	 call get_field(x, 'p', it, p, nz)
	 call get_field(x, 't', it, t, nz)
         call get_field(x, 'qv', it, qv, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          call interpolation2D(us, self%pg(iobs), u(:,:,1))
		  call interpolation2D(vs, self%pg(iobs), v(:,:,1))
                  do k = 1, nz 
                     call interpolation2D(p1D(nz+1-k), self%pg(iobs), p(:,:,k))
		     call interpolation2D(t1D(nz+1-k), self%pg(iobs), t(:,:,k))
		     call interpolation2D(qv1D(nz+1-k), self%pg(iobs), qv(:,:,k))
		     p1D(nz+1-k) = p1D(nz+1-k)*0.01d0
		     qv1D(nz+1-k) = qv1D(nz+1-k)*1000.0d0
		     qv1D(nz+1-k) = max(qv1D(nz+1-k),0.001d0)
                  end do
                  o31D(:) = 0.0d0              
		  call extend_profile(nz, nze, p1D, t1D, qv1D, o31D, p1De, t1De, qv1De, o31De)
		  qv1De(:) = qv1De(:)*0.001d0*KGKG2PPMV ! g/Kg --> Kg/Kg --> PPMV
                  ierror = 0
		  nchannel = rt%coefs(obsspace%satid(iobs))%coef%fmv_chn
		  call apply_rttov(iobs, nze, nchannel, obsspace, us, vs, p1De, t1De, qv1De, ierror, y)
                  if (ierror >= 2) then
                     valid%field(iobs,1:obsspace%nchannel(iobs)) = 0
		     y%field(iobs,1:obsspace%nchannel(iobs)) = 0.d0
                  else
		     valid%field(iobs,1:obsspace%nchannel(iobs)) = 1
		  end if          
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_rad1
   !
   !
   !
   subroutine interpolate_rad2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      logical :: update
      integer :: nz, nt, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it, imin, jmin, ii, jj
      real(r_size), dimension(:,:,:), allocatable, save :: u, v, p, t, qv
      !
      integer :: nze, ierror
      integer(kind=jpim) :: nchannel
      real(r_size) :: us, vs
      real(r_size), dimension(x%nlev) :: p1D, t1D, qv1D, o31D
      real(r_size), dimension(x%nlev+10) :: p1De, t1De, qv1De, o31De
      !
      nz = x%nlev; nt = obsspace%nt
      nze = nz + 10
      if (.not. allocated(u)) allocate(u(nx,ny,nz))
      if (.not. allocated(v)) allocate(v(nx,ny,nz))
      if (.not. allocated(p)) allocate(p(nx,ny,nz))
      if (.not. allocated(t)) allocate(t(nx,ny,nz))
      if (.not. allocated(qv)) allocate(qv(nx,ny,nz))
      !
      do it = 1, nt
	 do ixyt = 1, nxyt
	    if (it /= k2ijt(ixyt,3)) cycle
	    i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	    update = .False.
	    imin = max(1,i-1); jmin = max(1,j-1)
	    do ii = imin, i
	    do jj = jmin, j
	       if (.not. processed(ii,jj)) update = .True.
	    end do
	    end do
	    if (.not. update) cycle
	    !
	    call get_field(x, 'u', ixyt, u(i,j,1:nz), nz)
	    call get_field(x, 'v', ixyt, v(i,j,1:nz), nz)
	    call get_field(x, 'p', ixyt, p(i,j,1:nz), nz)
	    call get_field(x, 't', ixyt, t(i,j,1:nz), nz)
	    call get_field(x, 'qv', ixyt, qv(i,j,1:nz), nz)
	 end do
	 !
	 do ixyt = 1, nxyt
	    if (it /= k2ijt(ixyt,3)) cycle
	    i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	    if (processed(i,j)) cycle
	    if (obsspace%mobs(i,j,it) == 0) cycle
	    !
	    iobs1 = obsspace%iobs(i,j,it)
	    iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	    do iobs = iobs1, iobs2
	       call interpolation2D(us, self%pg(iobs), u(:,:,1))
	       call interpolation2D(vs, self%pg(iobs), v(:,:,1))
               do k = 1, nz 
                  call interpolation2D(p1D(nz+1-k), self%pg(iobs), p(:,:,k))
		  call interpolation2D(t1D(nz+1-k), self%pg(iobs), t(:,:,k))
		  call interpolation2D(qv1D(nz+1-k), self%pg(iobs), qv(:,:,k))
		  p1D(nz+1-k) = p1D(nz+1-k)*0.01d0
		  qv1D(nz+1-k) = qv1D(nz+1-k)*1000.0d0
		  qv1D(nz+1-k) = max(qv1D(nz+1-k),0.001d0)
               end do
               o31D(:) = 0.0d0           
	       call extend_profile(nz, nze, p1D, t1D, qv1D, o31D, p1De, t1De, qv1De, o31De)
	       qv1De(:) = qv1De(:)*0.001d0*KGKG2PPMV ! g/Kg --> Kg/Kg --> PPMV
               ierror = 0
	       nchannel = rt%coefs(obsspace%satid(iobs))%coef%fmv_chn
	       call apply_rttov(iobs, nze, nchannel, obsspace, us, vs, p1De, t1De, qv1De, ierror, y)
               if (ierror >= 2) then
                  valid%field(iobs,1:obsspace%nchannel(iobs)) = 0
		  y%field(iobs,1:obsspace%nchannel(iobs)) = 0.d0
               else
		  valid%field(iobs,1:obsspace%nchannel(iobs)) = 1
	       end if
	    end do
         end do
      end do
      !
      return
   end subroutine interpolate_rad2
   !
   !
   !
   subroutine interpolate_rad3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nlev
      integer :: ixyt, iobs, iobs1, iobs2, i, j, k, it
      real(r_size), dimension(:,:,:), allocatable, save :: u, v, p, t, qv
      !
      integer :: nze, ierror
      integer(kind=jpim) :: nchannel
      real(r_size) :: us, vs
      real(r_size), dimension(x%nlev) :: p1D, t1D, qv1D, o31D
      real(r_size), dimension(x%nlev+10) :: p1De, t1De, qv1De, o31De
      !
      nx = obsspace%nx; ny = obsspace%ny; nz = x%nlev
      nze = nz + 10
      if (.not. allocated(u)) allocate(u(nx,ny,nz))
      if (.not. allocated(v)) allocate(v(nx,ny,nz))
      if (.not. allocated(p)) allocate(p(nx,ny,nz))
      if (.not. allocated(t)) allocate(t(nx,ny,nz))
      if (.not. allocated(qv)) allocate(qv(nx,ny,nz))
      !
      do it = 1, nt
         do i = ip, ip+1
	 do j = jp, jp+1
	    ixyt = ijt2k(i-ip+1,j-jp+1,it)
	    if (ixyt == 0) cycle
	    !
	    call get_field(x, 'u', ixyt, u(i,j,1:nz), nz)
	    call get_field(x, 'v', ixyt, v(i,j,1:nz), nz)
	    call get_field(x, 'p', ixyt, p(i,j,1:nz), nz)
	    call get_field(x, 't', ixyt, t(i,j,1:nz), nz)
	    call get_field(x, 'qv', ixyt, qv(i,j,1:nz), nz)
	 end do
	 end do
	 !
	 i = ip; j = jp
	 ixyt = ijt2k(1,1,it)
	 if (ixyt == 0) cycle
	 if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    call interpolation2D(us, self%pg(iobs), u(:,:,1))
	    call interpolation2D(vs, self%pg(iobs), v(:,:,1))
            do k = 1, nz 
               call interpolation2D(p1D(nz+1-k), self%pg(iobs), p(:,:,k))
	       call interpolation2D(t1D(nz+1-k), self%pg(iobs), t(:,:,k))
	       call interpolation2D(qv1D(nz+1-k), self%pg(iobs), qv(:,:,k))
	       p1D(nz+1-k) = p1D(nz+1-k)*0.01d0
	       qv1D(nz+1-k) = qv1D(nz+1-k)*1000.0d0
	       qv1D(nz+1-k) = max(qv1D(nz+1-k),0.001d0)
            end do
            o31D(:) = 0.0d0        
	    call extend_profile(nz, nze, p1D, t1D, qv1D, o31D, p1De, t1De, qv1De, o31De)
	    qv1De(:) = qv1De(:)*0.001d0*KGKG2PPMV ! g/Kg --> Kg/Kg --> PPMV
            ierror = 0
	    nchannel = rt%coefs(obsspace%satid(iobs))%coef%fmv_chn
	    call apply_rttov(iobs, nze, nchannel, obsspace, us, vs, p1De, t1De, qv1De, ierror, y)
            if (ierror >= 2) then
               valid%field(iobs,1:obsspace%nchannel(iobs)) = 0
	       y%field(iobs,1:obsspace%nchannel(iobs)) = 0.d0
            else
	       valid%field(iobs,1:obsspace%nchannel(iobs)) = 1
	    end if
	 end do
      end do
      !
      return
   end subroutine interpolate_rad3
   !
   !
   !
   subroutine interpolate_Drad(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: u, v, p, t, qv, du, dv, dp, dt, dqv
      !
      integer :: nze, ierror
      integer(kind=jpim) :: nchannel
      real(r_size) :: us, vs, dus, dvs
      real(r_size), dimension(x%nlev) :: p1D, t1D, qv1D, o31D, dp1D, dt1D, dqv1D, do31D
      real(r_size), dimension(x%nlev+10) :: p1De, t1De, qv1De, o31De, dp1De, dt1De, dqv1De, do31De
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      nze = nz + 10
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      if (y%nobs > 0) y%field(:,:) = 0.d0
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) == 0) cycle
	 call get_field(xbck, 'u', 1, u, nz)
	 call get_field(xbck, 'v', 1, v, nz)
	 call get_field(xbck, 'p', it, p, nz)
	 call get_field(xbck, 't', it, t, nz)
         call get_field(xbck, 'qv', it, qv, nz)
	 call get_field(x, 'u', 1, du, nz)
	 call get_field(x, 'v', 1, dv, nz)
	 call get_field(x, 'p', it, dp, nz)
	 call get_field(x, 't', it, dt, nz)
         call get_field(x, 'qv', it, dqv, nz)
	 !
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          call interpolation2D(us, self%pg(iobs), u(:,:,1))
		  call interpolation2D(vs, self%pg(iobs), v(:,:,1))
		  call TL_interpolation2D(dus, du(:,:,1), self%pg(iobs))
		  call TL_interpolation2D(dvs, dv(:,:,1), self%pg(iobs))
                  do k = 1, nz 
                     call interpolation2D(p1D(nz+1-k), self%pg(iobs), p(:,:,k))
		     call interpolation2D(t1D(nz+1-k), self%pg(iobs), t(:,:,k))
		     call interpolation2D(qv1D(nz+1-k), self%pg(iobs), qv(:,:,k))
		     p1D(nz+1-k) = p1D(nz+1-k)*0.01d0
		     qv1D(nz+1-k) = qv1D(nz+1-k)*1000.0d0
		     qv1D(nz+1-k) = max(qv1D(nz+1-k),0.001d0)
		     call TL_interpolation2D(dp1D(nz+1-k), dp(:,:,k), self%pg(iobs))
		     call TL_interpolation2D(dt1D(nz+1-k), dt(:,:,k), self%pg(iobs))
		     call TL_interpolation2D(dqv1D(nz+1-k), dqv(:,:,k), self%pg(iobs))
		     dp1D(nz+1-k) = dp1D(nz+1-k)*0.01d0
		     dqv1D(nz+1-k) = dqv1D(nz+1-k)*1000.0d0
                  end do
		  o31D(:) = 0.0d0
		  do31D(:) = 0.0d0              
		  call extend_profile(nz, nze, p1D, t1D, qv1D, o31D, p1De, t1De, qv1De, o31De)
		  call extend_Dprofile(nz, nze, p1D, dp1D, dt1D, dqv1D, do31D, dp1De, dt1De, dqv1De, do31De)
		  qv1De(:) = qv1De(:)*0.001d0*KGKG2PPMV ! g/Kg --> Kg/Kg --> PPMV
		  dqv1De(:) = dqv1De(:)*0.001d0*KGKG2PPMV ! g/Kg --> Kg/Kg --> PPMV
		  !
		  ierror = 0
		  nchannel = rt%coefs(obsspace%satid(iobs))%coef%fmv_chn
		  call apply_Drttov(iobs, nze, nchannel, obsspace, us, vs, dus, dvs, &
                                   & p1De, t1De, qv1De, dp1De, dt1De, dqv1De, ierror, y)
	          if (ierror >= 2) then
                     valid%field(iobs,1:obsspace%nchannel(iobs)) = 0
		     y%field(iobs,1:obsspace%nchannel(iobs)) = 0.d0
		  else
	             valid%field(iobs,1:obsspace%nchannel(iobs)) = 1
		  end if
	       end do
	    end do
	 end do
      end do
      !
      return
   end subroutine interpolate_Drad
   !
   !
   !
   subroutine interpolate_DradT(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      integer :: nx, ny, nz, nt, nlev
      integer :: i, j, k, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      real(r_size), dimension(xbck%nx,xbck%ny,xbck%nlev) :: u, v, p, t, qv, du, dv, dp, dt, dqv
      !
      integer :: nze, ierror
      integer(kind=jpim) :: nchannel
      real(r_size) :: us, vs, dus, dvs
      real(r_size), dimension(x%nlev) :: p1D, t1D, qv1D, o31D, dp1D, dt1D, dqv1D, do31D
      real(r_size), dimension(x%nlev+10) :: p1De, t1De, qv1De, o31De, dp1De, dt1De, dqv1De, do31De
      !
      nx = x%nx; ny = x%ny; nz = x%nlev; nt = x%nt
      nze = nz + 10
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      do it = 1, nt
         du = 0.0d0
	 dv = 0.0d0
	 dp = 0.0d0
         dt = 0.0d0
	 dqv = 0.d0
	 !
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
	    call get_field(xbck, 'u', it, u, nz)
	    call get_field(xbck, 'v', it, v, nz)
	    call get_field(xbck, 'p', it, p, nz)
	    call get_field(xbck, 't', it, t, nz)
            call get_field(xbck, 'qv', it, qv, nz)
	    !
	    do j = js, je
	       do i = is, ie
		  iobs1 = obsspace%iobs(i,j,it)
		  iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
		  do iobs = iobs1, iobs2
		     if (valid%field(iobs,1) == 0) cycle
		     call interpolation2D(us, self%pg(iobs), u(:,:,1))
		     call interpolation2D(vs, self%pg(iobs), v(:,:,1))
                     do k = 1, nz 
                        call interpolation2D(p1D(nz+1-k), self%pg(iobs), p(:,:,k))
		        call interpolation2D(t1D(nz+1-k), self%pg(iobs), t(:,:,k))
		        call interpolation2D(qv1D(nz+1-k), self%pg(iobs), qv(:,:,k))
		        p1D(nz+1-k) = p1D(nz+1-k)*0.01d0
			qv1D(nz+1-k) = qv1D(nz+1-k)*1000.0d0
		        qv1D(nz+1-k) = max(qv1D(nz+1-k),0.001d0)
                     end do
                     o31D(:) = 0.0d0                 
		     call extend_profile(nz, nze, p1D, t1D, qv1D, o31D, p1De, t1De, qv1De, o31De)
		     qv1De(:) = qv1De(:)*0.001d0*KGKG2PPMV
		     !
		     ierror = 0
		     nchannel = rt%coefs(obsspace%satid(iobs))%coef%fmv_chn
		     call apply_DrttovT(iobs, nze, nchannel, obsspace, us, vs, &
                                      & p1De, t1De, qv1De, y, ierror, dus, dvs, dp1De, dt1De, dqv1De)
		     if (ierror >= 2) cycle
		     dqv1De(:) = dqv1De(:)*0.001d0*KGKG2PPMV
		     do31De(:) = 0.0d0
		     call extend_DprofileT(nz, nze, p1D, dp1De, dt1De, dqv1De, do31De, dp1D, dt1D, dqv1D, do31D)
		     dp1D(:) = dp1D(:)*0.01d0
		     dqv1D(:) = dqv1D(:)*1000.0d0
		     do k = 1, nz
                        call AD_interpolation2D(dp1D(nz+1-k), dp(:,:,k), self%pg(iobs))
                        call AD_interpolation2D(dt1D(nz+1-k), dt(:,:,k), self%pg(iobs))
			call AD_interpolation2D(dqv1D(nz+1-k), dqv(:,:,k), self%pg(iobs))
                     end do
		     call AD_interpolation2D(dus, du(:,:,1), self%pg(iobs))
		     call AD_interpolation2D(dvs, dv(:,:,1), self%pg(iobs))
		  end do
	       end do
	    end do
	 end if
	 !
	 call add_halo(info, du(:,:,1:1), nx, ny, 1)
	 call add_halo(info, dv(:,:,1:1), nx, ny, 1)
	 call add_halo(info, dp, nx, ny, nz)
	 call add_halo(info, dt, nx, ny, nz)
	 call add_halo(info, dqv, nx, ny, nz)
	 call add_field(x, 'u', it, du, nz)
	 call add_field(x, 'v', it, dv, nz)
	 call add_field(x, 'p', it, dp, nz)
	 call add_field(x, 't', it, dt, nz)
	 call add_field(x, 'qv', it, dqv, nz)
      end do
      !
      return
   end subroutine interpolate_DradT
   !
   !
   !
   subroutine apply_H_NodeHFieldRAD1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      call interpolate_rad1(self, info, x, obsspace, valid, y)
      !
      return
   end subroutine apply_H_NodeHFieldRAD1
   !
   !
   !
   subroutine apply_H_NodeHFieldRAD2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      call interpolate_rad2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      !
      return
   end subroutine apply_H_NodeHFieldRAD2
   !
   !
   !
   subroutine apply_H_NodeHFieldRAD3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      call interpolate_rad3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      !
      return
   end subroutine apply_H_NodeHFieldRAD3
   !
   !
   !
   subroutine initialize_DH_NodeHFieldRAD(self, info, x, obsspace, valid, y, ne)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: y
      !
      !call initialize_Drad(self, info, x, obsspace, valid, y)
      !
      return
   end subroutine initialize_DH_NodeHFieldRAD
   !
   !
   !
   subroutine apply_DH_NodeHFieldRAD(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      call interpolate_Drad(self, info, xbck, x, obsspace, valid, y)
      !
      return
   end subroutine apply_DH_NodeHFieldRAD
   !
   !
   !
   subroutine apply_DHT_NodeHFieldRAD(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldRAD), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldRAD), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      !
      call interpolate_DradT(self, info, xbck, obsspace, valid, y, x)
      !
      return
   end subroutine apply_DHT_NodeHFieldRAD
   !
   !
   !
end module NodeHFieldRAD_class
