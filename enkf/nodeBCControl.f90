module NodeBCControl_class
! Author: Le Duc
! Created date: 21 Jul 2014
   use variable, only : r_size
   use NodeBCField_class
   use NodeControl_class
   use NodeProfileControl_class
   implicit none
   !
   type NodeBCFieldPointer
      type(NodeBCField), pointer :: p
   end type NodeBCFieldPointer
   !
   type NodeBCControl
      integer :: ncontrol
      type(NodeBCFieldPointer), dimension(:), allocatable :: control
   end type NodeBCControl
   !
   interface new
      module procedure new_NodeBCControl
   end interface
   interface destroy
      module procedure destroy_NodeBCControl
   end interface
   interface display
      module procedure display_NodeBCControl
   end interface
   interface apply_bc
      module procedure apply_bc_NodeBCControl1
      module procedure apply_bc_NodeBCControl2
   end interface
   interface apply_bcT
      module procedure apply_bcT_NodeBCControl
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodeBCControl(self, control)
      implicit none
      type(NodeBCControl), intent(inout) :: self
      type(NodeControl), intent(in) :: control
      character(len=10) :: name
      integer :: ivar, nvar, nlev
      !
      call get_ncontrol(control, nvar)
      self%ncontrol = nvar
      allocate(self%control(self%ncontrol))
      do ivar = 1, nvar
         call get_name(control, ivar, name)
	 call get_nlev(control, ivar, nlev)
	 allocate(self%control(ivar)%p)
         call new(self%control(ivar)%p, name, nlev)
      end do
      !
      return
   end subroutine new_NodeBCControl
   !
   !
   !
   subroutine destroy_NodeBCControl(self)
      implicit none
      type(NodeBCControl), intent(inout) :: self
      !
      if (allocated(self%control)) deallocate(self%control)
      !
      return
   end subroutine destroy_NodeBCControl
   !
   !
   !
   subroutine display_NodeBCControl(self)
      implicit none
      type(NodeBCControl), intent(in) :: self
      integer :: ivar
      !
      print*, 'Control: ', self%ncontrol
      do ivar = 1, self%ncontrol
         call display(self%control(ivar)%p)
      end do
      !
      return
   end subroutine display_NodeBCControl
   !
   !
   !
   subroutine apply_bc_NodeBCControl1(self, x)
      implicit none
      type(NodeBCControl), intent(in) :: self
      type(NodeControl), intent(inout) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_bc(self%control(ivar)%p, x%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_bc_NodeBCControl1
   !
   !
   !
   subroutine apply_bc_NodeBCControl2(self, x)
      implicit none
      type(NodeBCControl), intent(in) :: self
      type(NodeProfileControl), intent(inout) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_bc(self%control(ivar)%p, x%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_bc_NodeBCControl2
   !
   !
   !
   subroutine apply_bcT_NodeBCControl(self, x)
      implicit none
      type(NodeBCControl), intent(in) :: self
      type(NodeControl), intent(inout) :: x
      integer :: ivar
      !
      do ivar = 1, self%ncontrol
         call apply_bcT(self%control(ivar)%p, x%control(ivar)%p)
      end do
      !
      return
   end subroutine apply_bcT_NodeBCControl
   !
   !
   !
end module NodeBCControl_class
