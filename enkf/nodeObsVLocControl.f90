module NodeObsVLocControl_class
! Author: Le Duc
! Created date: 15 Mar 2016
   use variable, only : r_size
   use NodeInfo_class
   use NodeObsVLocField_class
   use NodeObsSpaceControl_class
   use NodeObsControl_class
   implicit none
   !
   type NodeObsVLocFieldPointer
      type(NodeObsVLocField), pointer :: p
   end type NodeObsVLocFieldPointer
   !
   type NodeObsVLocControl
      integer :: ncontrol
      type(NodeObsVLocField), pointer :: u, v, t, p, rh, pwv, rvl, tc, re, ba, rad
      type(NodeObsVLocFieldPointer), dimension(:), allocatable :: control
   end type NodeObsVLocControl
   !
   interface new
      module procedure new_NodeObsVLocControl1
      module procedure new_NodeObsVLocControl2
   end interface
   interface destroy
      module procedure destroy_NodeObsVLocControl
   end interface
   interface display
      module procedure display_NodeObsVLocControl
   end interface
   interface assignment(=)
      module procedure copy_NodeObsVLocControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeObsVLocControl
   end interface
   interface get_name
      module procedure get_name_NodeObsVLocControl
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsVLocControl
   end interface
   interface get_field
      module procedure get_field_NodeObsVLocControl
   end interface
   interface set_field
      module procedure set_field_NodeObsVLocControl1
      module procedure set_field_NodeObsVLocControl2
   end interface
   interface add_vloc
      module procedure add_vloc_NodeObsVLocControl
   end interface
   interface subtract_vloc
      module procedure subtract_vloc_NodeObsVLocControl
   end interface
   interface power_vloc
      module procedure power_vloc_NodeObsVLocControl
   end interface
   interface multiply_vloc
      module procedure multiply_vloc_NodeObsVLocControl1
      module procedure multiply_vloc_NodeObsVLocControl2
      module procedure multiply_vloc_NodeObsVLocControl3
   end interface
   interface divide_vloc
      module procedure divide_vloc_NodeObsVLocControl1
      module procedure divide_vloc_NodeObsVLocControl2
   end interface
   interface set_logp
      module procedure set_logp_NodeObsVLocControl1
      module procedure set_logp_NodeObsVLocControl2
   end interface
   interface set_profile
      module procedure set_profile_NodeObsVLocControl1
      module procedure set_profile_NodeObsVLocControl2
   end interface
   interface ecorap
      module procedure ecorap_NodeObsVLocControl
   end interface
   interface allreduce_ensvloc
      module procedure allreduce_ensvloc_NodeObsVLocControl
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodeObsVLocControl1(self, obsspace, nz)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer, intent(in) :: nz
      character(len=10) :: name
      integer :: nvar, ivar
      !
      call get_ncontrol(obsspace, nvar)
      self%ncontrol = nvar
      allocate(self%control(self%ncontrol))
      do ivar = 1, nvar
         allocate(self%control(ivar)%p)
         call new(self%control(ivar)%p, obsspace%control(ivar)%p, nz)
	 !
	 call get_name(obsspace, ivar, name)
	 if (trim(name) == 'u') then
	    self%u => self%control(ivar)%p
	 else if (trim(name) == 'v') then
	    self%v => self%control(ivar)%p
	 else if (trim(name) == 't') then
	    self%t => self%control(ivar)%p
	 else if (trim(name) == 'p') then
	    self%p => self%control(ivar)%p
	 else if (trim(name) == 'rh') then
	    self%rh => self%control(ivar)%p
	 else if (trim(name) == 'pwv') then
	    self%pwv => self%control(ivar)%p
	 else if (trim(name) == 'rvl') then
	    self%rvl => self%control(ivar)%p
	 else if (trim(name) == 'tc') then
	    self%tc => self%control(ivar)%p
	 else if (trim(name) == 're') then
	    self%re => self%control(ivar)%p
	 else if (trim(name) == 'ba') then
	    self%ba => self%control(ivar)%p
	 else if (trim(name) == 'rad') then
	    self%rad => self%control(ivar)%p
	 end if
      end do
      !
      return
   end subroutine new_NodeObsVLocControl1
   !
   !
   !
   subroutine new_NodeObsVLocControl2(self, obsspace, nz, nvar)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer, intent(in) :: nz, nvar
      character(len=10) :: name
      integer :: ncontrol, ivar
      !
      call get_ncontrol(obsspace, ncontrol)
      self%ncontrol = ncontrol
      allocate(self%control(self%ncontrol))
      do ivar = 1, ncontrol
         allocate(self%control(ivar)%p)
         call new(self%control(ivar)%p, obsspace%control(ivar)%p, nz, nvar)
	 !
	 call get_name(obsspace, ivar, name)
	 if (trim(name) == 'u') then
	    self%u => self%control(ivar)%p
	 else if (trim(name) == 'v') then
	    self%v => self%control(ivar)%p
	 else if (trim(name) == 't') then
	    self%t => self%control(ivar)%p
	 else if (trim(name) == 'p') then
	    self%p => self%control(ivar)%p
	 else if (trim(name) == 'rh') then
	    self%rh => self%control(ivar)%p
	 else if (trim(name) == 'pwv') then
	    self%pwv => self%control(ivar)%p
	 else if (trim(name) == 'rvl') then
	    self%rvl => self%control(ivar)%p
	 else if (trim(name) == 'tc') then
	    self%tc => self%control(ivar)%p
	 else if (trim(name) == 're') then
	    self%re => self%control(ivar)%p
	 else if (trim(name) == 'ba') then
	    self%ba => self%control(ivar)%p
	 else if (trim(name) == 'rad') then
	    self%rad => self%control(ivar)%p
	 end if
      end do
      !
      return
   end subroutine new_NodeObsVLocControl2
   !
   !
   !
   subroutine destroy_NodeObsVLocControl(self)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeObsVLocControl
   !
   !
   !
   subroutine display_NodeObsVLocControl(self)
      implicit none
      type(NodeObsVLocControl), intent(in) :: self
      integer :: ivar
      !
      print*, 'Control: ', self%ncontrol
      do ivar = 1, self%ncontrol
         call display(self%control(ivar)%p)
      end do
      !
      return
   end subroutine display_NodeObsVLocControl
   !
   !
   !
   subroutine copy_NodeObsVLocControl(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsVLocControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         self%control(ivar)%p = object%control(ivar)%p
      end do
      !
      return
   end subroutine copy_NodeObsVLocControl
   !
   !
   !
   subroutine get_ncontrol_NodeObsVLocControl(self, ncontrol)
      implicit none
      type(NodeObsVLocControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeObsVLocControl
   !
   !
   !
   subroutine get_name_NodeObsVLocControl(self, ivar, name)
      implicit none
      type(NodeObsVLocControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeObsVLocControl
   !
   !
   !
   subroutine get_nobs_NodeObsVLocControl(self, ivar, nobs)
      implicit none
      type(NodeObsVLocControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nobs
      !
      call get_nobs(self%control(ivar)%p, nobs)
      !
      return
   end subroutine get_nobs_NodeObsVLocControl
   !
   !
   !
   subroutine get_field_NodeObsVLocControl(self, k, object)
      implicit none
      type(NodeObsVLocControl), intent(in) :: self
      integer, intent(in) :: k
      type(NodeObsControl), intent(inout) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call get_field(self%control(ivar)%p, k, object%control(ivar)%p)
      end do
      !
      return
   end subroutine get_field_NodeObsVLocControl
   !
   !
   !
   subroutine set_field_NodeObsVLocControl1(self, const)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine set_field_NodeObsVLocControl1
   !
   !
   !
   subroutine set_field_NodeObsVLocControl2(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsVLocControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_field_NodeObsVLocControl2
   !
   !
   !
   subroutine add_vloc_NodeObsVLocControl(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsVLocControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add_vloc(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine add_vloc_NodeObsVLocControl
   !
   !
   !
   subroutine subtract_vloc_NodeObsVLocControl(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsVLocControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call subtract_vloc(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine subtract_vloc_NodeObsVLocControl
   !
   !
   !
   subroutine power_vloc_NodeObsVLocControl(self, const)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call power_vloc(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine power_vloc_NodeObsVLocControl
   !
   !
   !
   subroutine multiply_vloc_NodeObsVLocControl1(self, const)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply_vloc(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine multiply_vloc_NodeObsVLocControl1
   !
   !
   !
   subroutine multiply_vloc_NodeObsVLocControl2(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsVLocControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply_vloc(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine multiply_vloc_NodeObsVLocControl2
   !
   !
   !
   subroutine multiply_vloc_NodeObsVLocControl3(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply_vloc(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine multiply_vloc_NodeObsVLocControl3
   !
   !
   !
   subroutine divide_vloc_NodeObsVLocControl1(self, const)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_vloc(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine divide_vloc_NodeObsVLocControl1
   !
   !
   !
   subroutine divide_vloc_NodeObsVLocControl2(self, object)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeObsVLocControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_vloc(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine divide_vloc_NodeObsVLocControl2
   !
   !
   !
   subroutine set_logp_NodeObsVLocControl1(self, info, obsspace, x)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeControl), intent(in) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_logp(self%control(ivar)%p, info, obsspace%control(ivar)%p, x)
      end do
      !
      return
   end subroutine set_logp_NodeObsVLocControl1
   !
   !
   !
   subroutine set_logp_NodeObsVLocControl2(self, info, processed, k2ijt, obsspace, x, nx, ny, nxyt)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeProfileControl), intent(in) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_logp(self%control(ivar)%p, info, processed, k2ijt, obsspace%control(ivar)%p, x, nx, ny, nxyt)
      end do
      !
      return
   end subroutine set_logp_NodeObsVLocControl2
   !
   !
   !
   subroutine set_profile_NodeObsVLocControl1(self, info, obsspace, x)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeControl), intent(in) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_profile(self%control(ivar)%p, info, obsspace%control(ivar)%p, x)
      end do
      !
      return
   end subroutine set_profile_NodeObsVLocControl1
   !
   !
   !
   subroutine set_profile_NodeObsVLocControl2(self, info, processed, k2ijt, obsspace, x, nx, ny, nxyt)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeProfileControl), intent(in) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_profile(self%control(ivar)%p, info, processed, k2ijt, obsspace%control(ivar)%p, x, nx, ny, nxyt)
      end do
      !
      return
   end subroutine set_profile_NodeObsVLocControl2
   !
   !
   !
   subroutine ecorap_NodeObsVLocControl(self, scale, logp, vloc)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      real(r_size), intent(in) :: scale
      type(NodeObsVLocControl), intent(in) :: logp
      type(NodeObsVLocControl), intent(inout) :: vloc
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call ecorap(self%control(ivar)%p, scale, logp%control(ivar)%p, vloc%control(ivar)%p)
      end do
      !
      return
   end subroutine ecorap_NodeObsVLocControl
   !
   !
   !
   subroutine allreduce_ensvloc_NodeObsVLocControl(self)
      implicit none
      type(NodeObsVLocControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensvloc(self%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ensvloc_NodeObsVLocControl
   !
   !
   !
end module NodeObsVLocControl_class
