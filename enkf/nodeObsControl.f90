module NodeObsControl_class
! Author: Le Duc
! Created date: 15 Mar 2016
   use variable, only : r_size
   use NodeInfo_class
   use NodeObsField_class
   use NodeObsSpaceControl_class
   use NodeObsValidControl_class
   implicit none
   !
   type NodeObsFieldPointer
      type(NodeObsField), pointer :: p
   end type NodeObsFieldPointer
   !
   type NodeObsControl
      integer :: ncontrol, nx, ny
      type(NodeObsField), pointer :: u, v, t, p, rh, pwv, rvl, tc, re, ba, rad
      type(NodeObsFieldPointer), dimension(:), allocatable :: control
   end type NodeObsControl
   !
   interface new
      module procedure new_NodeObsControl1
      module procedure new_NodeObsControl2
      module procedure new_NodeObsControl3
   end interface
   interface destroy
      module procedure destroy_NodeObsControl
   end interface
   interface display
      module procedure display_NodeObsControl
   end interface
   interface assignment(=)
      module procedure copy_NodeObsControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeObsControl
   end interface
   interface get_name
      module procedure get_name_NodeObsControl
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsControl
   end interface
   interface extract_y
      module procedure extract_y_NodeObsControl0
      module procedure extract_y_NodeObsControl1
      module procedure extract_y_NodeObsControl2
      module procedure extract_y_NodeObsControl3
      module procedure extract_y_NodeObsControl4
   end interface
   interface set_obs
      module procedure set_obs_NodeObsControl
   end interface
   interface set_error
      module procedure set_error_NodeObsControl
   end interface
   interface set_field
      module procedure set_field_NodeObsControl1
      module procedure set_field_NodeObsControl2
      module procedure set_field_NodeObsControl3
      module procedure set_field_NodeObsControl4
      module procedure set_field_NodeObsControl5
      module procedure set_field_NodeObsControl6
      module procedure set_field_NodeObsControl7
   end interface
   interface add_obs
      module procedure add_obs_NodeObsControl1
      module procedure add_obs_NodeObsControl2
      module procedure add_obs_NodeObsControl3
   end interface
   interface subtract_obs
      module procedure subtract_obs_NodeObsControl1
      module procedure subtract_obs_NodeObsControl2
      module procedure subtract_obs_NodeObsControl3
   end interface
   interface power_obs
      module procedure power_obs_NodeObsControl
   end interface
   interface loga10_obs
      module procedure loga10_obs_NodeObsControl
   end interface
   interface multiply_obs
      module procedure multiply_obs_NodeObsControl1
      module procedure multiply_obs_NodeObsControl2
   end interface
   interface divide_obs
      module procedure divide_obs_NodeObsControl1
      module procedure divide_obs_NodeObsControl2
      module procedure divide_obs_NodeObsControl3
      module procedure divide_obs_NodeObsControl4
      module procedure divide_obs_NodeObsControl5
      module procedure divide_obs_NodeObsControl6
   end interface
   interface compute_normsquare
      module procedure compute_normsquare_NodeObsControl
   end interface
   interface innerproduct
      module procedure innerproduct_NodeObsControl1
      module procedure innerproduct_NodeObsControl2
      module procedure innerproduct_NodeObsControl3
   end interface
   interface qccheck
      module procedure qccheck_NodeObsControl1
      module procedure qccheck_NodeObsControl2
   end interface
   interface stdcheck
      module procedure stdcheck_NodeObsControl
   end interface
   interface random_obs
      module procedure random_obs_NodeObsControl1
      module procedure random_obs_NodeObsControl2
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsControl
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsControl
   end interface
   interface broadcast_ensobs
      module procedure broadcast_ensobs_NodeObsControl
   end interface
   interface allreduce_ensobs
      module procedure allreduce_ensobs_NodeObsControl1
      module procedure allreduce_ensobs_NodeObsControl2
      module procedure allreduce_ensobs_NodeObsControl3
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodeObsControl1(self, obsspace)
      implicit none
      type(NodeObsControl), intent(inout) :: self
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
   end subroutine new_NodeObsControl1
   !
   !
   !
   subroutine new_NodeObsControl2(self, obsspace, obsname)
      implicit none
      type(NodeObsControl), intent(inout) :: self
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
   end subroutine new_NodeObsControl2
   !
   !
   !
   subroutine new_NodeObsControl3(self, obsspace, nsubobs)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer, intent(in) :: nsubobs
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
         call new(self%control(ivar)%p, obsspace%control(ivar)%p, nsubobs)
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
   end subroutine new_NodeObsControl3
   !
   !
   !
   subroutine destroy_NodeObsControl(self)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeObsControl
   !
   !
   !
   subroutine display_NodeObsControl(self)
      implicit none
      type(NodeObsControl), intent(in) :: self
      integer :: ivar
      !
      print*, 'Control: ', self%ncontrol
      do ivar = 1, self%ncontrol
         call display(self%control(ivar)%p)
      end do
      !
      return
   end subroutine display_NodeObsControl
   !
   !
   !
   subroutine copy_NodeObsControl(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         self%control(ivar)%p = object%control(ivar)%p
      end do
      !
      return
   end subroutine copy_NodeObsControl
   !
   !
   !
   subroutine get_ncontrol_NodeObsControl(self, ncontrol)
      implicit none
      type(NodeObsControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeObsControl
   !
   !
   !
   subroutine get_name_NodeObsControl(self, ivar, name)
      implicit none
      type(NodeObsControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeObsControl
   !
   !
   !
   subroutine get_nobs_NodeObsControl(self, ivar, nobs)
      implicit none
      type(NodeObsControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nobs
      !
      call get_nobs(self%control(ivar)%p, nobs)
      !
      return
   end subroutine get_nobs_NodeObsControl
   !
   !
   !
   subroutine extract_y_NodeObsControl0(self, obsspace, valid, mobs, y, mmax)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: mmax
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      integer, intent(out) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: ivar
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         call extract_y(self%control(ivar)%p,  obsspace%control(ivar)%p, valid%control(ivar)%p, mobs, y, mmax)
      end do
      !
      return
   end subroutine extract_y_NodeObsControl0
   !
   !
   !
   subroutine extract_y_NodeObsControl1(self, ip, jp, correlated, obsspace, valid, mobs, y, mmax)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp, mmax
      logical, dimension(self%ncontrol), intent(in) :: correlated
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      integer, intent(out) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: ivar
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         if (.not. correlated(ivar)) cycle
	 call extract_y(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, valid%control(ivar)%p, mobs, y, mmax)
      end do
      !
      return
   end subroutine extract_y_NodeObsControl1
   !
   !
   !
   subroutine extract_y_NodeObsControl2(self, i, j, ip, jp, correlated, obsspace, valid, mobs, y, np, mmax)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: i, j, np, mmax
      integer, dimension(np), intent(in) :: ip, jp
      logical, dimension(self%ncontrol), intent(in) :: correlated
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      integer, intent(out) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: ivar
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         if (.not. correlated(ivar)) cycle
	 call extract_y(self%control(ivar)%p, i, j, ip, jp, obsspace%control(ivar)%p, valid%control(ivar)%p, mobs, y, np, mmax)
      end do
      !
      return
   end subroutine extract_y_NodeObsControl2
   !
   !
   !
   subroutine extract_y_NodeObsControl3(self, i, j, ip, jp, correlated, obsspace, valid, xyloc, mobs, y, np, mmax)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: i, j, np, mmax
      integer, dimension(np), intent(in) :: ip, jp
      logical, dimension(self%ncontrol), intent(in) :: correlated
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      type(NodeObsControl), intent(in) :: xyloc
      integer, intent(out) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: ivar
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         if (.not. correlated(ivar)) cycle
	 call extract_y(self%control(ivar)%p, i, j, ip, jp, obsspace%control(ivar)%p, valid%control(ivar)%p, xyloc%control(ivar)%p, mobs, y, np, mmax)
      end do
      !
      return
   end subroutine extract_y_NodeObsControl3
   !
   !
   !
   subroutine extract_y_NodeObsControl4(self, i, j, ip, jp, correlated, obsspace, valid, xyloc, vloc, mobs, y, np, mmax)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: i, j, np, mmax
      integer, dimension(np), intent(in) :: ip, jp
      logical, dimension(self%ncontrol), intent(in) :: correlated
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      type(NodeObsControl), intent(in) :: xyloc
      type(NodeObsControl), intent(in) :: vloc
      integer, intent(out) :: mobs
      real(r_size), dimension(mmax), intent(out) :: y
      integer :: ivar
      !
      mobs = 0
      do ivar = 1, self%ncontrol
         if (.not. correlated(ivar)) cycle
	 call extract_y(self%control(ivar)%p, i, j, ip, jp, obsspace%control(ivar)%p, valid%control(ivar)%p, xyloc%control(ivar)%p, vloc%control(ivar)%p, mobs, y, np, mmax)
      end do
      !
      return
   end subroutine extract_y_NodeObsControl4
   !
   !
   !
   subroutine set_obs_NodeObsControl(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call set_obs(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_obs_NodeObsControl
   !
   !
   !
   subroutine set_error_NodeObsControl(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsSpaceControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call set_error(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_error_NodeObsControl
   !
   !
   !
   subroutine set_field_NodeObsControl1(self, const)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine set_field_NodeObsControl1
   !
   !
   !
   subroutine set_field_NodeObsControl2(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_field_NodeObsControl2
   !
   !
   !
   subroutine set_field_NodeObsControl3(self, object, obsname)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      character(len=*), intent(in) :: obsname
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(self, ivar, name)
         if (trim(name) /= trim(obsname)) cycle
	 call set_field(self%control(ivar)%p, object%control(1)%p)
      end do
      !
      return
   end subroutine set_field_NodeObsControl3
   !
   !
   !
   subroutine set_field_NodeObsControl4(self, ip, jp, obsspace, const)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, const)
      end do
      !
      return
   end subroutine set_field_NodeObsControl4
   !
   !
   !
   subroutine set_field_NodeObsControl5(self, ip, jp, obsspace, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine set_field_NodeObsControl5
   !
   !
   !
   subroutine set_field_NodeObsControl6(self, processed, k2ijt, obsspace, const, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, const, nxyt)
      end do
      !
      return
   end subroutine set_field_NodeObsControl6
   !
   !
   !
   subroutine set_field_NodeObsControl7(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call set_field(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine set_field_NodeObsControl7
   !
   !
   !
   subroutine add_obs_NodeObsControl1(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add_obs(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine add_obs_NodeObsControl1
   !
   !
   !
   subroutine add_obs_NodeObsControl2(self, ip, jp, obsspace, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add_obs(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine add_obs_NodeObsControl2
   !
   !
   !
   subroutine add_obs_NodeObsControl3(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call add_obs(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine add_obs_NodeObsControl3
   !
   !
   !
   subroutine subtract_obs_NodeObsControl1(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call subtract_obs(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine subtract_obs_NodeObsControl1
   !
   !
   !
   subroutine subtract_obs_NodeObsControl2(self, ip, jp, obsspace, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call subtract_obs(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine subtract_obs_NodeObsControl2
   !
   !
   !
   subroutine subtract_obs_NodeObsControl3(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call subtract_obs(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine subtract_obs_NodeObsControl3
   !
   !
   !
   subroutine power_obs_NodeObsControl(self, const)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call power_obs(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine power_obs_NodeObsControl
   !
   !
   !
   subroutine loga10_obs_NodeObsControl(self)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call loga10_obs(self%control(ivar)%p)
      end do
      !
      return
   end subroutine loga10_obs_NodeObsControl
   !
   !
   !
   subroutine multiply_obs_NodeObsControl1(self, const)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply_obs(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine multiply_obs_NodeObsControl1
   !
   !
   !
   subroutine multiply_obs_NodeObsControl2(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call multiply_obs(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine multiply_obs_NodeObsControl2
   !
   !
   !
   subroutine divide_obs_NodeObsControl1(self, const)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_obs(self%control(ivar)%p, const)
      end do
      !
      return
   end subroutine divide_obs_NodeObsControl1
   !
   !
   !
   subroutine divide_obs_NodeObsControl2(self, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_obs(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine divide_obs_NodeObsControl2
   !
   !
   !
   subroutine divide_obs_NodeObsControl3(self, ip, jp, obsspace, const)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_obs(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, const)
      end do
      !
      return
   end subroutine divide_obs_NodeObsControl3
   !
   !
   !
   subroutine divide_obs_NodeObsControl4(self, ip, jp, obsspace, object)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_obs(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine divide_obs_NodeObsControl4
   !
   !
   !
   subroutine divide_obs_NodeObsControl5(self, processed, k2ijt, obsspace, const, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      real(r_size), intent(in) :: const
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_obs(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, const, nxyt)
      end do
      !
      return
   end subroutine divide_obs_NodeObsControl5
   !
   !
   !
   subroutine divide_obs_NodeObsControl6(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsControl), intent(in) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call divide_obs(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, object%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine divide_obs_NodeObsControl6
   !
   !
   !
   subroutine compute_normsquare_NodeObsControl(self, normsquare)
      implicit none
      type(NodeObsControl), intent(in) :: self
      !type(NodeObsValidControl), intent(in) :: valid
      real(r_size), intent(out) :: normsquare
      integer :: ivar
      real(r_size) :: tmp
      !
      normsquare = 0.d0
      do ivar = 1, self%ncontrol
	 !call compute_normsquare(self%control(ivar)%p, valid%control(ivar)%p, tmp)
	 call compute_normsquare(self%control(ivar)%p, tmp)
	 normsquare = normsquare + tmp
      end do
      !
      return
   end subroutine compute_normsquare_NodeObsControl
   !
   !
   !
   subroutine innerproduct_NodeObsControl1(self, object, product)
      implicit none
      type(NodeObsControl), intent(in) :: self
      type(NodeObsControl), intent(in) :: object
      real(r_size), intent(out) :: product
      integer :: ivar
      real(r_size) :: tmp
      !
      product = 0.d0
      do ivar = 1, self%ncontrol
	 call innerproduct(self%control(ivar)%p, object%control(ivar)%p, tmp)
	 product = product + tmp
      end do
      !
      return
   end subroutine innerproduct_NodeObsControl1
   !
   !
   !
   subroutine innerproduct_NodeObsControl2(self, object, valid, product)
      implicit none
      type(NodeObsControl), intent(in) :: self
      type(NodeObsControl), intent(in) :: object
      type(NodeObsValidControl), intent(in) :: valid
      real(r_size), intent(out) :: product
      integer :: ivar
      real(r_size) :: tmp
      !
      product = 0.d0
      do ivar = 1, self%ncontrol
	 call innerproduct(self%control(ivar)%p, object%control(ivar)%p, valid%control(ivar)%p, tmp)
	 product = product + tmp
      end do
      !
      return
   end subroutine innerproduct_NodeObsControl2
   !
   !
   !
   subroutine innerproduct_NodeObsControl3(self, object, obsspace, valid, product)
      implicit none
      type(NodeObsControl), intent(in) :: self
      type(NodeObsControl), intent(in) :: object
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      real(r_size), intent(out) :: product
      integer :: ivar
      real(r_size) :: tmp
      !
      product = 0.d0
      do ivar = 1, self%ncontrol
	 call innerproduct(self%control(ivar)%p, object%control(ivar)%p, obsspace%control(ivar)%p, valid%control(ivar)%p, tmp)
	 product = product + tmp
      end do
      !
      return
   end subroutine innerproduct_NodeObsControl3
   !
   !
   !
   subroutine qccheck_NodeObsControl1(self, threshold, valid)
      implicit none
      type(NodeObsControl), intent(in) :: self
      real(r_size), intent(in) :: threshold
      type(NodeObsValidControl), intent(inout) :: valid
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         !call get_name(self, ivar, name)
         !if (trim(name) == 'p') then
            !call qccheck(self%control(ivar)%p, 8.d0, valid%control(ivar)%p)
         !else
	    call qccheck(self%control(ivar)%p, threshold, valid%control(ivar)%p)
	 !end if
      end do
      !
      return
   end subroutine qccheck_NodeObsControl1
   !
   !
   !
   subroutine qccheck_NodeObsControl2(self, threshold, ip, jp, obsspace, valid)
      implicit none
      type(NodeObsControl), intent(in) :: self
      real(r_size), intent(in) :: threshold
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call qccheck(self%control(ivar)%p, threshold, ip, jp, obsspace%control(ivar)%p, valid%control(ivar)%p)
      end do
      !
      return
   end subroutine qccheck_NodeObsControl2
   !
   !
   !
   subroutine stdcheck_NodeObsControl(self, valid)
      implicit none
      type(NodeObsControl), intent(in) :: self
      type(NodeObsValidControl), intent(inout) :: valid
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call stdcheck(self%control(ivar)%p, valid%control(ivar)%p)
      end do
      !
      return
   end subroutine stdcheck_NodeObsControl
   !
   !
   !
   subroutine random_obs_NodeObsControl1(self)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call random_obs(self%control(ivar)%p)
      end do
      !
      return
   end subroutine random_obs_NodeObsControl1
   !
   !
   !
   subroutine random_obs_NodeObsControl2(self, valid)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      type(NodeObsValidControl), intent(in) :: valid
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
	 call random_obs(self%control(ivar)%p, valid%control(ivar)%p)
      end do
      !
      return
   end subroutine random_obs_NodeObsControl2
   !
   !
   !
   subroutine scatter_obs_NodeObsControl(self, info, obsspace, local_obsspace, local_object)
      implicit none
      type(NodeObsControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(in) :: obsspace, local_obsspace
      type(NodeObsControl), intent(inout) :: local_object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call scatter_obs(self%control(ivar)%p, info, obsspace%control(ivar)%p, &
	                & local_obsspace%control(ivar)%p, local_object%control(ivar)%p)
      end do
      !
      return
   end subroutine scatter_obs_NodeObsControl
   !
   !
   !
   subroutine gather_obs_NodeObsControl(self, info, obsspace, global_obsspace, global_object)
      implicit none
      type(NodeObsControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(in) :: obsspace, global_obsspace
      type(NodeObsControl), intent(inout) :: global_object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call gather_obs(self%control(ivar)%p, info, obsspace%control(ivar)%p, &
	               & global_obsspace%control(ivar)%p, global_object%control(ivar)%p)
      end do
      !
      return
   end subroutine gather_obs_NodeObsControl
   !
   !
   !
   subroutine broadcast_ensobs_NodeObsControl(self, source)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: source
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call broadcast_ensobs(self%control(ivar)%p, source)
      end do
      !
      return
   end subroutine broadcast_ensobs_NodeObsControl
   !
   !
   !
   subroutine allreduce_ensobs_NodeObsControl1(self)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensobs(self%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ensobs_NodeObsControl1
   !
   !
   !
   subroutine allreduce_ensobs_NodeObsControl2(self, ip, jp, obsspace)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: ip, jp
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensobs(self%control(ivar)%p, ip, jp, obsspace%control(ivar)%p)
      end do
      !
      return
   end subroutine allreduce_ensobs_NodeObsControl2
   !
   !
   !
   subroutine allreduce_ensobs_NodeObsControl3(self, processed, k2ijt, obsspace, nxyt)
      implicit none
      type(NodeObsControl), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceControl), intent(in) :: obsspace
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call allreduce_ensobs(self%control(ivar)%p, processed, k2ijt, obsspace%control(ivar)%p, nxyt)
      end do
      !
      return
   end subroutine allreduce_ensobs_NodeObsControl3
   !
   !
   !
end module NodeObsControl_class
