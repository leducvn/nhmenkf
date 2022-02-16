module NodeSigma_class
! Author: Le Duc
! Created date: 22 Dec 2014
   use variable, only : r_size, r_sngl
   use NodeInfo_class
   use NodeMPI
   implicit none
   !
   type NodeSigma
      integer :: nx, ny, nlev
      real(r_size), dimension(:,:,:), allocatable :: sigma
   end type NodeSigma
   !
   interface new
      module procedure new_NodeSigma
   end interface
   interface destroy
      module procedure destroy_NodeSigma
   end interface
   interface display
      module procedure display_NodeSigma
   end interface
   interface get_nlev
      module procedure get_nlev_NodeSigma
   end interface
   interface get_sigma
      module procedure get_sigma_NodeSigma
   end interface
   interface set_sigma1D
      module procedure set_sigma1D_NodeSigma
   end interface
   interface set_sigma3D
      module procedure set_sigma3D_NodeSigma
   end interface
   interface apply_sigma
      module procedure apply_sigma_NodeSigma
   end interface
   !
contains
   !
   subroutine new_NodeSigma(self, nx, ny, nlev)
      implicit none
      type(NodeSigma), intent(inout) :: self
      integer, intent(in) :: nx, ny, nlev
      !
      self%nx = nx
      self%ny = ny
      self%nlev = nlev
      allocate(self%sigma(nx,ny,nlev))
      !
      return
   end subroutine new_NodeSigma
   !
   !
   !
   subroutine destroy_NodeSigma(self)
      implicit none
      type(NodeSigma), intent(inout) :: self
      !
      if (allocated(self%sigma)) deallocate(self%sigma)
      !
      return
   end subroutine destroy_NodeSigma
   !
   !
   !
   subroutine display_NodeSigma(self)
      implicit none
      type(NodeSigma), intent(in) :: self
      !
      print*, 'Dimension: ', self%nx, self%ny, self%nlev
      !
      return
   end subroutine display_NodeSigma
   !
   !
   !
   subroutine get_nlev_NodeSigma(self, nlev)
      implicit none
      type(NodeSigma), intent(in) :: self
      integer, intent(out) :: nlev
      !
      nlev = self%nlev
      !
      return
   end subroutine get_nlev_NodeSigma
   !
   !
   !
   subroutine get_sigma_NodeSigma(self, field)
      implicit none
      type(NodeSigma), intent(in) :: self
      real(r_size), dimension(self%nx,self%ny,self%nlev), intent(out) :: field
      !
      field(:,:,:) = self%sigma
      !
      return
   end subroutine get_sigma_NodeSigma
   !
   !
   !
   subroutine set_sigma1D_NodeSigma(self, info)
      implicit none
      type(NodeSigma), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      character(5) :: infile='sigma'
      integer :: nproc, myid, nlev, k
      real(r_sngl), dimension(:), allocatable :: sigma
      !
      nlev = self%nlev
      call get_myid(info, myid, nproc)
      if (myid == 0) then
         allocate(sigma(nlev))
	 open(90, file=infile, form='unformatted')
	 read(90) sigma(:)
	 close(90)
	 !print*, 'Sigma: ', nlev, minval(sigma), maxval(sigma)
	 self%sigma(1,1,:) = sigma(:)
	 deallocate(sigma)
      end if
      call broadcast1D('all', self%sigma(1,1,:), 0, nlev)
      do k = 1, nlev
	 self%sigma(:,:,k) = self%sigma(1,1,k)
      end do
      !
      return
   end subroutine set_sigma1D_NodeSigma
   !
   !
   !
   subroutine set_sigma3D_NodeSigma(self, info)
      implicit none
      type(NodeSigma), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      character(5) :: infile='sigma'
      integer :: nproc, myid, myide, nx0, ny0, nlev, nx, ny
      real(r_sngl), dimension(:,:,:), allocatable :: sigma
      real(r_size), dimension(:,:,:), allocatable :: datum
      !
      nlev = self%nlev
      call get_myid(info, myid, nproc)
      call get_myide(info, myide)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      call get_nx(info, nx)
      call get_ny(info, ny)
      allocate(datum(nx0,ny0,nlev))
      !
      if (myid == 0) then
	 allocate(sigma(nx0,ny0,nlev))
	 print*, 'Read ', infile
	 open(90, file=infile, form='unformatted')
	 read(90) sigma
	 close(90)
	 datum(:,:,:) = sigma
	 deallocate(sigma)
      end if
      if (myide == 0) call scatter(info, 'xy', datum, 0, nx0, ny0, nlev)
      call broadcast3D('e', datum(1:nx,1:ny,:), 0, nx, ny, nlev)
      self%sigma(1:nx,1:ny,:) = datum(1:nx,1:ny,:)
      deallocate(datum)
      !
      return
   end subroutine set_sigma3D_NodeSigma
   !
   !
   !
   subroutine apply_sigma_NodeSigma(self, info, field, nx, ny, nz0)
      implicit none
      type(NodeSigma), intent(in) :: self
      integer, intent(in) :: nx, ny, nz0
      type(NodeInfo), intent(in) :: info
      real(r_size), dimension(nx,ny,nz0), intent(inout) :: field
      !
      field = self%sigma*field
      !
      return
   end subroutine apply_sigma_NodeSigma
   !
   !
   !
end module NodeSigma_class
