module NodeProfileField_class
! Author: Le Duc
! Created date: 13 Jul 2014
   use variable, only : r_size
   use NodeInfo_class
   use NodeObsSpaceField_class
   use NodeMPI
   implicit none
   !
   type NodeProfileField
      character(len=10) :: name
      integer :: nxyt, nlev, nx, ny
      real(r_size), dimension(:,:), allocatable :: field
   end type NodeProfileField
   !
   interface new
      module procedure new_NodeProfileField
   end interface
   interface destroy
      module procedure destroy_NodeProfileField
   end interface
   interface display
      module procedure display_NodeProfileField
   end interface
   interface assignment(=)
      module procedure copy_NodeProfileField
   end interface
   interface get_name
      module procedure get_name_NodeProfileField
   end interface
   interface get_nxyt
      module procedure get_nxyt_NodeProfileField
   end interface
   interface get_nlev
      module procedure get_nlev_NodeProfileField
   end interface
   interface get_field
      module procedure get_field_NodeProfileField1
      module procedure get_field_NodeProfileField2
      module procedure get_field_NodeProfileField3
      module procedure get_field_NodeProfileField4
   end interface
   interface set_field
      module procedure set_field_NodeProfileField1
      module procedure set_field_NodeProfileField2
      module procedure set_field_NodeProfileField3
      module procedure set_field_NodeProfileField4
      module procedure set_field_NodeProfileField5
   end interface
   interface add
      module procedure add_NodeProfileField
   end interface
   interface subtract
      module procedure subtract_NodeProfileField
   end interface
   interface multiply
      module procedure multiply_NodeProfileField1
      module procedure multiply_NodeProfileField2
   end interface
   interface divide
      module procedure divide_NodeProfileField1
      module procedure divide_NodeProfileField2
   end interface
   interface allreduce_ens
      module procedure allreduce_ens_NodeProfileField
   end interface
   !
contains
   !
   subroutine new_NodeProfileField(self, name, nxyt, nlev, nx, ny)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      character(len=*), intent(in) :: name
      integer, intent(in) :: nxyt, nlev, nx, ny
      !
      self%name = name
      self%nxyt = nxyt
      self%nlev = nlev
      self%nx = nx
      self%ny = ny
      if (nxyt > 0) then
         allocate(self%field(nxyt,nlev))
         self%field(:,:) = 0.d0
      end if
      !
      return
   end subroutine new_NodeProfileField
   !
   !
   !
   subroutine destroy_NodeProfileField(self)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      !
      if (allocated(self%field)) deallocate(self%field)
      !
      return
   end subroutine destroy_NodeProfileField
   !
   !
   !
   subroutine display_NodeProfileField(self)
      implicit none
      type(NodeProfileField), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nxyt, self%nlev
      !
      return
   end subroutine display_NodeProfileField
   !
   !
   !
   subroutine copy_NodeProfileField(self, f)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      type(NodeProfileField), intent(in) :: f
      !
      if (self%nxyt == 0) return
      self%field(:,:) = f%field(:,:)
      !
      return
   end subroutine copy_NodeProfileField
   !
   !
   !
   subroutine get_name_NodeProfileField(self, name)
      implicit none
      type(NodeProfileField), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeProfileField
   !
   !
   !
   subroutine get_nxyt_NodeProfileField(self, nxyt)
      implicit none
      type(NodeProfileField), intent(in) :: self
      integer, intent(out) :: nxyt
      !
      nxyt = self%nxyt
      !
      return
   end subroutine get_nxyt_NodeProfileField
   !
   !
   !
   subroutine get_nlev_NodeProfileField(self, nlev)
      implicit none
      type(NodeProfileField), intent(in) :: self
      integer, intent(out) :: nlev
      !
      nlev = self%nlev
      !
      return
   end subroutine get_nlev_NodeProfileField
   !
   !
   !
   subroutine get_field_NodeProfileField1(self, it, k2ijt, field, nx, ny)
      implicit none
      type(NodeProfileField), intent(in) :: self
      integer, intent(in) :: it, nx, ny
      integer, dimension(self%nxyt,3), intent(in) :: k2ijt
      real(r_size), dimension(nx,ny,self%nlev), intent(inout) :: field
      integer :: i, j, k
      !
      do k = 1, self%nxyt
         if (it /= k2ijt(k,3)) cycle
	 i = k2ijt(k,1); j = k2ijt(k,2)
         field(i,j,:) = self%field(k,:)
      end do
      !
      return
   end subroutine get_field_NodeProfileField1
   !
   !
   !
   subroutine get_field_NodeProfileField2(self, ixyt, field)
      implicit none
      type(NodeProfileField), intent(in) :: self
      integer, intent(in) :: ixyt
      real(r_size), dimension(self%nlev), intent(inout) :: field
      !
      field(:) = self%field(ixyt,:)
      !
      return
   end subroutine get_field_NodeProfileField2
   !
   !
   !
   subroutine get_field_NodeProfileField3(self, field)
      implicit none
      type(NodeProfileField), intent(in) :: self
      real(r_size), dimension(self%nxyt,self%nlev), intent(inout) :: field
      !
      field(:,:) = self%field(:,:)
      !
      return
   end subroutine get_field_NodeProfileField3
   !
   !
   !
   subroutine get_field_NodeProfileField4(self, ixyt, i, j, field, nx, ny)
      implicit none
      type(NodeProfileField), intent(in) :: self
      integer, intent(in) :: ixyt, i, j, nx, ny
      real(r_size), dimension(nx,ny,self%nlev), intent(inout) :: field
      !
      field(i,j,:) = self%field(ixyt,:)
      !
      return
   end subroutine get_field_NodeProfileField4
   !
   !
   !
   subroutine set_field_NodeProfileField1(self, field)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      real(r_size), intent(in) :: field
      !
      if (self%nxyt == 0) return
      self%field(:,:) = field
      !
      return
   end subroutine set_field_NodeProfileField1
   !
   !
   !
   subroutine set_field_NodeProfileField2(self, ixyt, field)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      integer, intent(in) :: ixyt
      real(r_size), dimension(self%nlev), intent(in) :: field
      !
      self%field(ixyt,:) = field(:)
      !
      return
   end subroutine set_field_NodeProfileField2
   !
   !
   !
   subroutine set_field_NodeProfileField3(self, field)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      real(r_size), dimension(self%nxyt,self%nlev), intent(in) :: field
      !
      self%field(:,:) = field(:,:)
      !
      return
   end subroutine set_field_NodeProfileField3
   !
   !
   !
   subroutine set_field_NodeProfileField4(self, it, k2ijt, field, nx, ny)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      integer, intent(in) :: it, nx, ny
      integer, dimension(self%nxyt,3), intent(in) :: k2ijt
      real(r_size), dimension(nx,ny,self%nlev), intent(in) :: field
      integer :: i, j, ixyt
      !
      do ixyt = 1, self%nxyt
         if (it /= k2ijt(ixyt,3)) cycle
	 i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
         self%field(ixyt,:) = field(i,j,:)
      end do
      !
      return
   end subroutine set_field_NodeProfileField4
   !
   !
   !
   subroutine set_field_NodeProfileField5(self, processed, k2ijt, obsspace, object)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(self%nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeProfileField), intent(in) :: object
      logical :: update
      integer :: ixyt, i, j, imin, jmin, ii, jj
      !
      do ixyt = 1, self%nxyt
	 i = k2ijt(ixyt,1); j = k2ijt(ixyt,2)
	 update = .False.
	 imin = max(1,i-1); jmin = max(1,j-1)
	 do ii = imin, i
	 do jj = jmin, j
	    if (.not. processed(ii,jj)) update = .True.
	 end do
	 end do
	 if (.not. update) cycle
         self%field(ixyt,:) = object%field(ixyt,:)
      end do
      !
      return
   end subroutine set_field_NodeProfileField5
   !
   !
   !
   subroutine add_NodeProfileField(self, object)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      type(NodeProfileField), intent(in) :: object
      !
      if (self%nxyt == 0) return
      self%field(:,:) = self%field(:,:) + object%field(:,:)
      !
      return
   end subroutine add_NodeProfileField
   !
   !
   !
   subroutine subtract_NodeProfileField(self, object)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      type(NodeProfileField), intent(in) :: object
      !
      if (self%nxyt == 0) return
      self%field(:,:) = self%field(:,:) - object%field(:,:)
      !
      return
   end subroutine subtract_NodeProfileField
   !
   !
   !
   subroutine multiply_NodeProfileField1(self, const)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nxyt == 0) return
      self%field(:,:) = const*self%field(:,:)
      !
      return
   end subroutine multiply_NodeProfileField1
   !
   !
   !
   subroutine multiply_NodeProfileField2(self, object)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      type(NodeProfileField), intent(in) :: object
      !
      if (self%nxyt == 0) return
      self%field(:,:) = object%field(:,:)*self%field(:,:)
      !
      return
   end subroutine multiply_NodeProfileField2
   !
   !
   !
   subroutine divide_NodeProfileField1(self, const)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      real(r_size), intent(in) :: const
      !
      if (self%nxyt == 0) return
      self%field(:,:) = self%field(:,:)/const
      !
      return
   end subroutine divide_NodeProfileField1
   !
   !
   !
   subroutine divide_NodeProfileField2(self, object)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      type(NodeProfileField), intent(in) :: object
      !
      if (self%nxyt == 0) return
      where (object%field(:,:) < 1.e-12)
	 self%field(:,:) = 0.d0
      else where
	 self%field(:,:) = self%field(:,:)/object%field(:,:)
      end where
      !
      return
   end subroutine divide_NodeProfileField2
   !
   !
   !
   subroutine allreduce_ens_NodeProfileField(self)
      implicit none
      type(NodeProfileField), intent(inout) :: self
      !
      if (self%nxyt == 0) return
      call allreduce2D('e', self%field(:,:), self%nxyt, self%nlev)
      !
      return
   end subroutine allreduce_ens_NodeProfileField
   !
   !
   !
end module NodeProfileField_class
   
