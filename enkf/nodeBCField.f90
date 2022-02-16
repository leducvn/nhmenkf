module NodeBCField_class
! Author: Le Duc
! Created date: 05 Apr 2016
   use variable, only : r_size
   use NodeField_class
   use NodeProfileField_class
   implicit none
   !
   type NodeBCField
      character(len=10) :: name
      integer :: nlev, top
   end type NodeBCField
   integer, parameter :: dirichlet = 1, neumann = 2
   !
   interface new
      module procedure new_NodeBCField
   end interface
   interface display
      module procedure display_NodeBCField
   end interface
   interface apply_bc
      module procedure apply_bc_NodeBCField1
      module procedure apply_bc_NodeBCField2
   end interface
   interface apply_bcT
      module procedure apply_bcT_NodeBCField
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodeBCField(self, name, nlev)
      implicit none
      type(NodeBCField), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: nlev
      !
      self%name = name
      self%nlev = nlev
      self%top = neumann
      !
      return
   end subroutine new_NodeBCField
   !
   !
   !
   subroutine display_NodeBCField(self)
      implicit none
      type(NodeBCField), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Boundary condition: ', self%nlev, self%top
      !
      return
   end subroutine display_NodeBCField
   !
   !
   !
   subroutine apply_bc_NodeBCField1(self, x)
      implicit none
      type(NodeBCField), intent(in) :: self
      type(NodeField), intent(inout) :: x
      integer :: nz
      !
      nz = self%nlev
      if (nz == 1) return
      if (self%top == neumann) then
         x%field(:,:,nz,:) = x%field(:,:,nz-1,:)
      end if
      !
      return
   end subroutine apply_bc_NodeBCField1
   !
   !
   !
   subroutine apply_bc_NodeBCField2(self, x)
      implicit none
      type(NodeBCField), intent(in) :: self
      type(NodeProfileField), intent(inout) :: x
      integer :: nz
      !
      if (x%nxyt == 0) return
      nz = self%nlev
      if (nz == 1) return
      if (self%top == neumann) then
         x%field(:,nz) = x%field(:,nz-1)
      end if
      !
      return
   end subroutine apply_bc_NodeBCField2
   !
   !
   !
   subroutine apply_bcT_NodeBCField(self, x)
      implicit none
      type(NodeBCField), intent(in) :: self
      type(NodeField), intent(inout) :: x
      integer :: nz
      !
      nz = self%nlev
      if (nz == 1) return
      if (self%top == neumann) then
         x%field(:,:,nz-1,:) = x%field(:,:,nz-1,:) + x%field(:,:,nz,:)
         x%field(:,:,nz,:) = 0.d0
      end if
      !
      return
   end subroutine apply_bcT_NodeBCField
   !
   !
   !
end module NodeBCField_class
