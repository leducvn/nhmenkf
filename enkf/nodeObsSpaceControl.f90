module NodeObsSpaceControl_class
! Author: Le Duc
! Created date: 15 Mar 2016
   use NodeInfo_class
   use NodeObsSpaceField_class
   implicit none
   !
   type NodeObsSpaceFieldPointer
      type(NodeObsSpaceField), pointer :: p
   end type NodeObsSpaceFieldPointer
   !
   type NodeObsSpaceControl
      integer :: ncontrol, nx, ny, nt
      type(NodeObsSpaceField), pointer :: u, v, t, p, rh, pwv, rvl, tc, re, ba, rad
      type(NodeObsSpaceFieldPointer), dimension(:), allocatable :: control
   end type NodeObsSpaceControl
   !
   interface new_conv
      module procedure new_conv_NodeObsSpaceControl
   end interface
   interface new_gnss
      module procedure new_gnss_NodeObsSpaceControl
   end interface
   interface new_rad
      module procedure new_rad_NodeObsSpaceControl
   end interface
   interface new_convgnss
      module procedure new_convgnss_NodeObsSpaceControl
   end interface
   interface new_convrad
      module procedure new_convrad_NodeObsSpaceControl
   end interface
   interface destroy
      module procedure destroy_NodeObsSpaceControl
   end interface
   interface display
      module procedure display_NodeObsSpaceControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeObsSpaceControl
   end interface
   interface get_name
      module procedure get_name_NodeObsSpaceControl
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsSpaceControl
   end interface
   interface get_nsubobs
      module procedure get_nsubobs_NodeObsSpaceControl
   end interface
   interface get_mobs
      module procedure get_mobs_NodeObsSpaceControl1
      module procedure get_mobs_NodeObsSpaceControl2
   end interface
   interface read_obs
      module procedure read_obs_NodeObsSpaceControl
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsSpaceControl
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsSpaceControl
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_conv_NodeObsSpaceControl(self, nx, ny, nt)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer, intent(in) :: nx, ny, nt
      integer :: ivar
      !
      self%ncontrol = 7
      self%nx = nx
      self%ny = ny
      self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      ! pwv
      !allocate(self%pwv)
      !call new(self%pwv, 'pwv', 'cnv', nx, ny, nt)
      !ivar = ivar + 1
      !self%control(ivar)%p => self%pwv
      !return
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! p
      allocate(self%p)
      call new(self%p, 'p', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ! rh
      allocate(self%rh)
      call new(self%rh, 'rh', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rh
      ! pwv
      allocate(self%pwv)
      call new(self%pwv, 'pwv', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%pwv
      ! rvl
      allocate(self%rvl)
      call new(self%rvl, 'rvl', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rvl
      !
      return
   end subroutine new_conv_NodeObsSpaceControl
   !
   !
   !
   subroutine new_gnss_NodeObsSpaceControl(self, nx, ny, nt)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer, intent(in) :: nx, ny, nt
      integer :: ivar
      !
      self%ncontrol = 1
      self%nx = nx
      self%ny = ny
      self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! re
      allocate(self%re)
      call new(self%re, 're', 'gnss', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%re
      ! ba
      !allocate(self%ba)
      !call new(self%ba, 'ba', 'gnss', nx, ny, nt)
      !ivar = ivar + 1
      !self%control(ivar)%p => self%ba
      !
      return
   end subroutine new_gnss_NodeObsSpaceControl
   !
   !
   !
   subroutine new_rad_NodeObsSpaceControl(self, nx, ny, nt)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer, intent(in) :: nx, ny, nt
      integer :: ivar
      !
      self%ncontrol = 1
      self%nx = nx
      self%ny = ny
      self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! rad
      allocate(self%rad)
      call new(self%rad, 'rad', 'rad', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rad
      !
      return
   end subroutine new_rad_NodeObsSpaceControl
   !
   !
   !
   subroutine new_convgnss_NodeObsSpaceControl(self, nx, ny, nt)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer, intent(in) :: nx, ny, nt
      integer :: ivar
      !
      self%ncontrol = 8
      self%nx = nx
      self%ny = ny
      self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      ! tc
      !allocate(self%tc)
      !call new(self%tc, 'tc', 'tc', nx, ny, nt)
      !ivar = ivar + 1
      !self%control(ivar)%p => self%tc
      !return
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! p
      allocate(self%p)
      call new(self%p, 'p', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ! rh
      allocate(self%rh)
      call new(self%rh, 'rh', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rh
      ! pwv
      allocate(self%pwv)
      call new(self%pwv, 'pwv', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%pwv
      ! rvl
      allocate(self%rvl)
      call new(self%rvl, 'rvl', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rvl
      ! re
      allocate(self%re)
      call new(self%re, 're', 'gnss', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%re
      ! ba
      !allocate(self%ba)
      !call new(self%ba, 'ba', 'gnss', nx, ny, nt)
      !ivar = ivar + 1
      !self%control(ivar)%p => self%ba
      !
      return
   end subroutine new_convgnss_NodeObsSpaceControl
   !
   !
   !
   subroutine new_convrad_NodeObsSpaceControl(self, nx, ny, nt)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer, intent(in) :: nx, ny, nt
      integer :: ivar
      !
      self%ncontrol = 9
      self%nx = nx
      self%ny = ny
      self%nt = nt
      allocate(self%control(self%ncontrol))
      ivar = 0
      !
      ! u
      allocate(self%u)
      call new(self%u, 'u', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%u
      ! v
      allocate(self%v)
      call new(self%v, 'v', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%v
      ! t
      allocate(self%t)
      call new(self%t, 't', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%t
      ! p
      allocate(self%p)
      call new(self%p, 'p', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%p
      ! rh
      allocate(self%rh)
      call new(self%rh, 'rh', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rh
      ! pwv
      allocate(self%pwv)
      call new(self%pwv, 'pwv', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%pwv
      ! rvl
      allocate(self%rvl)
      call new(self%rvl, 'rvl', 'cnv', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rvl
      ! re
      allocate(self%re)
      call new(self%re, 're', 'gnss', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%re
      ! ba
      !allocate(self%ba)
      !call new(self%ba, 'ba', 'gnss', nx, ny, nt)
      !ivar = ivar + 1
      !self%control(ivar)%p => self%ba
      ! rad
      allocate(self%rad)
      call new(self%rad, 'rad', 'rad', nx, ny, nt)
      ivar = ivar + 1
      self%control(ivar)%p => self%rad
      !
      return
   end subroutine new_convrad_NodeObsSpaceControl
   !
   !
   !
   subroutine destroy_NodeObsSpaceControl(self)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeObsSpaceControl
   !
   !
   !
   subroutine display_NodeObsSpaceControl(self)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer :: ivar
      !
      print*, 'Control: ', self%ncontrol
      do ivar = 1, self%ncontrol
         call display(self%control(ivar)%p)
      end do
      !
      return
   end subroutine display_NodeObsSpaceControl
   !
   !
   !
   subroutine get_ncontrol_NodeObsSpaceControl(self, ncontrol)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeObsSpaceControl
   !
   !
   !
   subroutine get_name_NodeObsSpaceControl(self, ivar, name)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeObsSpaceControl
   !
   !
   !
   subroutine get_nobs_NodeObsSpaceControl(self, ivar, nobs)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nobs
      !
      call get_nobs(self%control(ivar)%p, nobs)
      !
      return
   end subroutine get_nobs_NodeObsSpaceControl
   !
   !
   !
   subroutine get_nsubobs_NodeObsSpaceControl(self, ivar, nsubobs)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nsubobs
      !
      call get_nsubobs(self%control(ivar)%p, nsubobs)
      !
      return
   end subroutine get_nsubobs_NodeObsSpaceControl
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceControl1(self, i, j, mobs)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer, intent(in) :: i, j
      integer, intent(out) :: mobs
      integer :: ivar, m
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         call get_mobs(self%control(ivar)%p, i, j, m)
	 mobs = mobs + m
      end do
      !
      return
   end subroutine get_mobs_NodeObsSpaceControl1
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceControl2(self, i, j, it, mobs)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      integer, intent(in) :: i, j, it
      integer, intent(out) :: mobs
      integer :: ivar, m
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         call get_mobs(self%control(ivar)%p, i, j, it, m)
	 mobs = mobs + m
      end do
      !
      return
   end subroutine get_mobs_NodeObsSpaceControl2
   !
   !
   !
   subroutine read_obs_NodeObsSpaceControl(self, myid)
      implicit none
      type(NodeObsSpaceControl), intent(inout) :: self
      integer, intent(in) :: myid
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call read_obs(self%control(ivar)%p, myid)
      end do
      !
      return
   end subroutine read_obs_NodeObsSpaceControl
   !
   !
   !
   subroutine scatter_obs_NodeObsSpaceControl(self, info, local_object)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(inout) :: local_object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call scatter_obs(self%control(ivar)%p, info, local_object%control(ivar)%p)
      end do
      !
      return
   end subroutine scatter_obs_NodeObsSpaceControl
   !
   !
   !
   subroutine gather_obs_NodeObsSpaceControl(self, info, global_object)
      implicit none
      type(NodeObsSpaceControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(inout) :: global_object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call gather_obs(self%control(ivar)%p, info, global_object%control(ivar)%p)
      end do
      !
      return
   end subroutine gather_obs_NodeObsSpaceControl
   !
   !
   !
end module NodeObsSpaceControl_class
