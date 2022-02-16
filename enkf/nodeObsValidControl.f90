module NodeObsValidControl_class
! Author: Le Duc
! Created date: 15 Mar 2016
   use NodeInfo_class
   use NodeObsValidField_class
   use NodeObsSpaceControl_class
   implicit none
   !
   type NodeObsValidFieldPointer
      type(NodeObsValidField), pointer :: p
   end type NodeObsValidFieldPointer
   !
   type NodeObsValidControl
      integer :: ncontrol, nx, ny
      type(NodeObsValidField), pointer :: u, v, t, p, rh, pwv, rvl, tc, re, ba, rad
      type(NodeObsValidFieldPointer), dimension(:), allocatable :: control
   end type NodeObsValidControl
   !
   interface new
      module procedure new_NodeObsValidControl1
      module procedure new_NodeObsValidControl2
   end interface
   interface destroy
      module procedure destroy_NodeObsValidControl
   end interface
   interface display
      module procedure display_NodeObsValidControl
   end interface
   interface assignment(=)
      module procedure copy_NodeObsValidControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeObsValidControl
   end interface
   interface get_name
      module procedure get_name_NodeObsValidControl
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsValidControl
   end interface
   interface set_field
      module procedure set_field_NodeObsValidControl1
      module procedure set_field_NodeObsValidControl2
      module procedure set_field_NodeObsValidControl3
      module procedure set_field_NodeObsValidControl4
   end interface
   interface and
      module procedure and_NodeObsValidControl1
      module procedure and_NodeObsValidControl2
      module procedure and_NodeObsValidControl3
      module procedure and_NodeObsValidControl4
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsValidControl
   end interface
   interface broadcast_ensvalid
      module procedure broadcast_ensvalid_NodeObsValidControl
   end interface
   interface allreduce_ensvalid
      module procedure allreduce_ensvalid_NodeObsValidControl1
      module procedure allreduce_ensvalid_NodeObsValidControl2
      module procedure allreduce_ensvalid_NodeObsValidControl3
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodeObsValidControl1(self, obsspace)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=10) :: name
      integer :: nvar, ivar
      !
      call get_ncontrol(obsspace, nvar)
      self%ncontrol = nvar
      self%nx = obsspace%nx
      self%ny = obsspace%ny
      allocate(self%control(self%ncontrol))
      do ivar = 1, nvar
         allocate(self%control(ivar)%p)
         call new(self%control(ivar)%p, obsspace%control(ivar)%p)
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
   end subroutine new_NodeObsValidControl1
   !
   !
   !
   subroutine new_NodeObsValidControl2(self, obsspace, obsname)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=*), intent(in) :: obsname
      character(len=10) :: name
      integer :: nvar, ivar
      !
      self%ncontrol = 1
      self%nx = obsspace%nx
      self%ny = obsspace%ny
      allocate(self%control(self%ncontrol))
      allocate(self%control(1)%p)
      call get_ncontrol(obsspace, nvar)
      do ivar = 1, nvar
	 call get_name(obsspace, ivar, name)
	 if (trim(name) /= trim(obsname)) cycle
	 call new(self%control(1)%p, obsspace%control(ivar)%p)
	 if (trim(name) == 'u') then
	    self%u => self%control(1)%p
	 else if (trim(name) == 'v') then
	    self%v => self%control(1)%p
	 else if (trim(name) == 't') then
	    self%t => self%control(1)%p
	 else if (trim(name) == 'p') then
	    self%p => self%control(1)%p
	 else if (trim(name) == 'rh') then
	    self%rh => self%control(1)%p
	 else if (trim(name) == 'pwv') then
	    self%pwv => self%control(1)%p
	 else if (trim(name) == 'rvl') then
	    self%rvl => self%control(1)%p
	 else if (trim(name) == 'tc') then
	    self%tc => self%control(1)%p
	 else if (trim(name) == 're') then
	    self%re => self%control(1)%p
	 else if (trim(name) == 'ba') then
	    self%ba => self%control(1)%p
	 else if (trim(name) == 'rad') then
	    self%rad => self%control(1)%p
	 end if
      end do
      !
      return
   end subroutine new_NodeObsValidControl2
   !
   !
   !
   subroutine destroy_NodeObsValidControl(self)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeObsValidControl
   !
   !
   !
   subroutine display_NodeObsValidControl(self)
      implicit none
      type(NodeObsValidControl), intent(in) :: self
      integer :: ivar
      !
      print*, 'Control: ', self%ncontrol
      do ivar = 1, self%ncontrol
         call display(self%control(ivar)%p)
      end do
      !
      return
   end subroutine display_NodeObsValidControl
   !
   !
   !
   subroutine copy_NodeObsValidControl(self, object)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      type(NodeObsValidControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         self%control(ivar)%p = object%control(ivar)%p
      end do
      !
      return
   end subroutine copy_NodeObsValidControl
   !
   !
   !
   subroutine get_ncontrol_NodeObsValidControl(self, ncontrol)
      implicit none
      type(NodeObsValidControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeObsValidControl
   !
   !
   !
   subroutine get_name_NodeObsValidControl(self, ivar, name)
      implicit none
      type(NodeObsValidControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeObsValidControl
   !
   !
   !
   subroutine get_nobs_NodeObsValidControl(self, ivar, nobs)
      implicit none
      type(NodeObsValidControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nobs
      !
      call get_nobs(self%control(ivar)%p, nobs)
      !
      return
   end subroutine get_nobs_NodeObsValidControl
   !
   !
   !
   subroutine set_field_NodeObsValidControl1(self, const)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine set_field_NodeObsValidControl1
   !
   !
   !
   subroutine set_field_NodeObsValidControl2(self, object)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      type(NodeObsValidControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_field_NodeObsValidControl2
   !
   !
   !
   subroutine set_field_NodeObsValidControl3(self, processed, k2ijt, obsspace, const, nxyt)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer, intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, const, nxyt)
      end do
      !
      return
   end subroutine set_field_NodeObsValidControl3
   !
   !
   !
   subroutine set_field_NodeObsValidControl4(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine set_field_NodeObsValidControl4
   !
   !
   !
   subroutine and_NodeObsValidControl1(self, object)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      type(NodeObsValidControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call and(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine and_NodeObsValidControl1
   !
   !
   !
   subroutine and_NodeObsValidControl2(self, object, obsname)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      type(NodeObsValidControl), intent(in) :: object
      character(len=*), intent(in) :: obsname
      integer :: ivar
      character(len=10) :: name
      !
      do ivar = 1, self%ncontrol
         call get_name(self, ivar, name)
	 if (trim(name) /= trim(obsname)) cycle
	 call and(self%control(ivar)%p, object%control(1)%p)
      end do
      !
      return
   end subroutine and_NodeObsValidControl2
   !
   !
   !
   subroutine and_NodeObsValidControl3(self, ip, jp, obsspace, object)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call and(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine and_NodeObsValidControl3
   !
   !
   !
   subroutine and_NodeObsValidControl4(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call and(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine and_NodeObsValidControl4
   !
   !
   !
   subroutine gather_obs_NodeObsValidControl(self, info, obsspace, global_obsspace, global_object)
      implicit none
      type(NodeObsValidControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(in) :: obsspace, global_obsspace
      type(NodeObsValidControl), intent(inout) :: global_object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call gather_obs(self%control(ivar)%p, info, obsspace%control(ivar)%p, &
	               & global_obsspace%control(ivar)%p, global_object%control(ivar)%p)
      end do
      !
      return
   end subroutine gather_obs_NodeObsValidControl
   !
   !
   !
   subroutine broadcast_ensvalid_NodeObsValidControl(self, source)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: source
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call broadcast_ensvalid(self%control(ivar)%p, source)
      end do
      !
      return
   end subroutine broadcast_ensvalid_NodeObsValidControl
   !
   !
   !
   subroutine allreduce_ensvalid_NodeObsValidControl1(self)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensvalid(self%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ensvalid_NodeObsValidControl1
   !
   !
   !
   subroutine allreduce_ensvalid_NodeObsValidControl2(self, ip, jp, obsspace)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensvalid(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ensvalid_NodeObsValidControl2
   !
   !
   !
   subroutine allreduce_ensvalid_NodeObsValidControl3(self, processed, k2ijt, obsspace, nxyt)
      implicit none
      type(NodeObsValidControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensvalid(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine allreduce_ensvalid_NodeObsValidControl3
   !
   !
   !
end module NodeObsValidControl_class
