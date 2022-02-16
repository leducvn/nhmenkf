module NodeHControl_class
! Author: Le Duc
! Created date: 15 Mar 2016
   use NodeInfo_class
   use NodeHField_class
   use NodeObsSpaceControl_class
   use NodeObsValidControl_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeProfileControl_class
   implicit none
   !
   type NodeHFieldPointer
      type(NodeHField), pointer :: p
   end type NodeHFieldPointer
   !
   type NodeHControl
      integer :: ncontrol
      type(NodeHField), pointer :: u, v, t, p, rh, pwv, rvl, tc, re, ba, rad
      type(NodeHFieldPointer), dimension(:), allocatable :: control
   end type NodeHControl
   !
   interface new
      module procedure new_NodeHControl
   end interface
   interface destroy
      module procedure destroy_NodeHControl
   end interface
   interface display
      module procedure display_NodeHControl
   end interface
   interface get_ncontrol
      module procedure get_ncontrol_NodeHControl
   end interface
   interface get_name
      module procedure get_name_NodeHControl
   end interface
   interface get_nobs
      module procedure get_nobs_NodeHControl
   end interface
   interface get_xyloc
      module procedure get_xyloc_NodeHControl1
      module procedure get_xyloc_NodeHControl2
   end interface
   interface apply_Hlogp
      module procedure apply_Hlogp_NodeHControl
   end interface
   interface apply_H
      module procedure apply_H_NodeHControl1
      module procedure apply_H_NodeHControl2
      module procedure apply_H_NodeHControl3
      module procedure apply_H_NodeHControl4
   end interface
   interface initialize_DH
      module procedure initialize_DH_NodeHControl
   end interface
   interface apply_DH
      module procedure apply_DH_NodeHControl1
      module procedure apply_DH_NodeHControl2
   end interface
   interface apply_DHlocal
      module procedure apply_DHlocal_NodeHControl
   end interface
   interface apply_DHT
      module procedure apply_DHT_NodeHControl1
      module procedure apply_DHT_NodeHControl2
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodeHControl(self, info, obsspace)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=10) :: name
      integer :: nvar, ivar
      !
      call get_ncontrol(obsspace, nvar)
      self%ncontrol = nvar
      allocate(self%control(self%ncontrol))
      do ivar = 1, nvar
         allocate(self%control(ivar)%p)
         call new(self%control(ivar)%p, info, obsspace%control(ivar)%p)
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
   end subroutine new_NodeHControl
   !
   !
   !
   subroutine destroy_NodeHControl(self)
      implicit none
      type(NodeHControl), intent(inout) :: self
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call destroy(self%control(ivar)%p)
      end do
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeHControl
   !
   !
   !
   subroutine display_NodeHControl(self)
      implicit none
      type(NodeHControl), intent(in) :: self
      integer :: ivar
      !
      print*, 'Control: ', self%ncontrol
      do ivar = 1, self%ncontrol
         call display(self%control(ivar)%p)
      end do
      !
      return
   end subroutine display_NodeHControl
   !
   !
   !
   subroutine get_ncontrol_NodeHControl(self, ncontrol)
      implicit none
      type(NodeHControl), intent(in) :: self
      integer, intent(out) :: ncontrol
      !
      ncontrol = self%ncontrol
      !
      return
   end subroutine get_ncontrol_NodeHControl
   !
   !
   !
   subroutine get_name_NodeHControl(self, ivar, name)
      implicit none
      type(NodeHControl), intent(in) :: self
      integer, intent(in) :: ivar
      character(len=10), intent(out) :: name
      !
      call get_name(self%control(ivar)%p, name)
      !
      return
   end subroutine get_name_NodeHControl
   !
   !
   !
   subroutine get_nobs_NodeHControl(self, ivar, nobs)
      implicit none
      type(NodeHControl), intent(in) :: self
      integer, intent(in) :: ivar
      integer, intent(out) :: nobs
      !
      call get_nobs(self%control(ivar)%p, nobs)
      !
      return
   end subroutine get_nobs_NodeHControl
   !
   !
   !
   subroutine get_xyloc_NodeHControl1(self, object)
      implicit none
      type(NodeHControl), intent(in) :: self
      type(NodeObsControl), intent(inout) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_xyloc(self%control(ivar)%p, object%control(ivar)%p)
      end do
      !
      return
   end subroutine get_xyloc_NodeHControl1
   !
   !
   !
   subroutine get_xyloc_NodeHControl2(self, info, object)
      implicit none
      type(NodeHControl), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsControl), intent(inout) :: object
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_xyloc(self%control(ivar)%p, info, object%control(ivar)%p)
      end do
      !
      return
   end subroutine get_xyloc_NodeHControl2
   !
   !
   !
   subroutine apply_Hlogp_NodeHControl(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_Hlogp(self%control(ivar)%p, info, x, obsspace%control(ivar)%p, &
	                & valid%control(ivar)%p, y%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_Hlogp_NodeHControl
   !
   !
   !
   subroutine apply_H_NodeHControl1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_H(self%control(ivar)%p, info, x, obsspace%control(ivar)%p, &
	            & valid%control(ivar)%p, y%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_H_NodeHControl1
   !
   !
   !
   subroutine apply_H_NodeHControl2(self, info, x, obsspace, obsname, valid, y)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=*), intent(in) :: obsname
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(obsspace, ivar, name)
	 if (trim(name) /= trim(obsname)) cycle
         call apply_H(self%control(ivar)%p, info, x, obsspace%control(ivar)%p, &
	            & valid%control(1)%p, y%control(1)%p)
      end do
      !
      return
   end subroutine apply_H_NodeHControl2
   !
   !
   !
   subroutine apply_H_NodeHControl3(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_H(self%control(ivar)%p, info, processed, k2ijt, x, obsspace%control(ivar)%p, &
	            & valid%control(ivar)%p, y%control(ivar)%p, nx, ny, nxyt)
      end do
      !
      return
   end subroutine apply_H_NodeHControl3
   !
   !
   !
   subroutine apply_H_NodeHControl4(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_H(self%control(ivar)%p, info, ip, jp, ijt2k, x, obsspace%control(ivar)%p, &
	            & valid%control(ivar)%p, y%control(ivar)%p, nt)
      end do
      !
      return
   end subroutine apply_H_NodeHControl4
   !
   !
   !
   subroutine initialize_DH_NodeHControl(self, info, x, obsspace, obsname, valid, y, ne)
      implicit none
      type(NodeHControl), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=*), intent(in) :: obsname
      type(NodeObsValidControl), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: y
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(obsspace, ivar, name)
	 if (trim(name) /= trim(obsname)) cycle
         call initialize_DH(self%control(ivar)%p, info, x, obsspace%control(ivar)%p, &
                          & valid%control(1)%p, y, ne)
      end do
      !
      return
   end subroutine initialize_DH_NodeHControl
   !
   !
   !
   subroutine apply_DH_NodeHControl1(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_DH(self%control(ivar)%p, info, xbck, x, obsspace%control(ivar)%p, &
	             & valid%control(ivar)%p, y%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_DH_NodeHControl1
   !
   !
   !
   subroutine apply_DH_NodeHControl2(self, info, xbck, x, obsspace, obsname, valid, y)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=*), intent(in) :: obsname
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(obsspace, ivar, name)
	 if (trim(name) /= trim(obsname)) cycle
         call apply_DH(self%control(ivar)%p, info, xbck, x, obsspace%control(ivar)%p, &
	             & valid%control(1)%p, y%control(1)%p)
      end do
      !
      return
   end subroutine apply_DH_NodeHControl2
   !
   !
   !
   subroutine apply_DHlocal_NodeHControl(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(inout) :: valid
      type(NodeObsControl), intent(inout) :: y
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(obsspace, ivar, name)
	 if (trim(name) == 'tc') cycle
         call apply_DH(self%control(ivar)%p, info, xbck, x, obsspace%control(ivar)%p, &
	             & valid%control(ivar)%p, y%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_DHlocal_NodeHControl
   !
   !
   !
   subroutine apply_DHT_NodeHControl1(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceControl), intent(in) :: obsspace
      type(NodeObsValidControl), intent(in) :: valid
      type(NodeObsControl), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_DHT(self%control(ivar)%p, info, xbck, obsspace%control(ivar)%p, &
	              & valid%control(ivar)%p, y%control(ivar)%p, x)
      end do
      !
      return
   end subroutine apply_DHT_NodeHControl1
   !
   !
   !
   subroutine apply_DHT_NodeHControl2(self, info, xbck, obsspace, obsname, valid, y, x)
      implicit none
      type(NodeHControl), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceControl), intent(in) :: obsspace
      character(len=*), intent(in) :: obsname
      type(NodeObsValidControl), intent(in) :: valid
      type(NodeObsControl), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      character(len=10) :: name
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call get_name(obsspace, ivar, name)
	 if (trim(name) /= trim(obsname)) cycle
         call apply_DHT(self%control(ivar)%p, info, xbck, obsspace%control(ivar)%p, &
	              & valid%control(1)%p, y%control(1)%p, x)
      end do
      !
      return
   end subroutine apply_DHT_NodeHControl2
   !
   !
   !
end module NodeHControl_class
