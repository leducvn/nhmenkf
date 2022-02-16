module NodeVModulation_class
! Author: Le Duc
! Created date: 17 Aug 2018
   use variable, only : r_size, r_sngl
   use NodeInfo_class
   use NodeMPI
   implicit none
   !
   type NodeVModulation
      integer :: ncorr, nrank
      real(r_size), dimension(:), allocatable :: dcoef, lcoef
      real(r_size), dimension(:,:), allocatable :: lvector
   end type NodeVModulation
   !
   interface new
      module procedure new_NodeVModulation
   end interface
   interface destroy
      module procedure destroy_NodeVModulation
   end interface
   interface display
      module procedure display_NodeVModulation
   end interface
   interface set_VModulation
      module procedure set_VModulation_NodeVModulation
   end interface
   interface modulate
      module procedure modulate_NodeVModulation1
      module procedure modulate_NodeVModulation2
      module procedure modulate_NodeVModulation3
   end interface
   interface multiply
      module procedure multiply_NodeVModulation
   end interface
   !
contains
   !
   subroutine new_NodeVModulation(self, ncorr)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      integer, intent(in) :: ncorr
      !
      self%ncorr = ncorr
      !
      return
   end subroutine new_NodeVModulation
   !
   !
   !
   subroutine destroy_NodeVModulation(self)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      !
      if (allocated(self%dcoef)) deallocate(self%dcoef)
      if (allocated(self%lcoef)) deallocate(self%lcoef)
      if (allocated(self%lvector)) deallocate(self%lvector)
      !
      return
   end subroutine destroy_NodeVModulation
   !
   !
   !
   subroutine display_NodeVModulation(self)
      implicit none
      type(NodeVModulation), intent(in) :: self
      !
      print*, 'Dimension: ', self%ncorr, self%nrank
      !
      return
   end subroutine display_NodeVModulation
   !
   !
   !
   subroutine set_VModulation_NodeVModulation(self, info, filename)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: filename
      integer :: nproc, myid, ncorr, nrank
      real(r_sngl), dimension(:), allocatable :: dcoef, lcoef
      real(r_sngl), dimension(:,:), allocatable :: lvector
      !
      ncorr = self%ncorr
      call get_myid(info, myid, nproc)
      if (myid == 0) then
	 allocate(dcoef(ncorr), lcoef(ncorr))
         allocate(lvector(ncorr,ncorr))
	 open(90, file=trim(filename), form='unformatted')
	 read(90) nrank
	 read(90) dcoef(:)
	 read(90) lcoef(1:nrank), lvector(:,1:nrank)
	 close(90)
      end if
      !
      call int_broadcast0D('all', nrank, 0)
      self%nrank = nrank
      if (allocated(self%dcoef)) deallocate(self%dcoef)
      if (allocated(self%lcoef)) deallocate(self%lcoef)
      if (allocated(self%lvector)) deallocate(self%lvector)
      allocate(self%dcoef(ncorr), self%lcoef(nrank), self%lvector(ncorr,nrank))
      if (myid == 0) then
         self%dcoef(:) = dcoef(:)
	 self%lcoef(:) = lcoef(1:nrank)
	 self%lvector(:,:) = lvector(:,1:nrank)
	 deallocate(dcoef, lcoef, lvector)
      end if
      call broadcast1D('all', self%dcoef, 0, ncorr)
      call broadcast1D('all', self%lcoef, 0, nrank)
      call broadcast2D('all', self%lvector, 0, ncorr, nrank)
      !
      return
   end subroutine set_VModulation_NodeVModulation
   !
   !
   !
   subroutine modulate_NodeVModulation1(self, imode, field, nx, ny, nlev)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      integer, intent(in) :: imode, nx, ny, nlev
      real(r_size), dimension(nx,ny,nlev), intent(inout) :: field
      integer :: i, j
      !
      do i = 1, nx
	 do j = 1, ny
	    if (imode == 0) then
	       field(i,j,:) = self%dcoef(:)*field(i,j,:)
	    else
	       field(i,j,:) = self%lvector(:,imode)*field(i,j,:)
	    end if
	 end do
      end do
      !
      return
   end subroutine modulate_NodeVModulation1
   !
   !
   !
   subroutine modulate_NodeVModulation2(self, info, imode, field, nx, ny, nlev)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      integer, intent(in) :: imode, nx, ny, nlev
      type(NodeInfo), intent(in) :: info
      real(r_size), dimension(nx,ny,nlev), intent(inout) :: field
      integer :: dis, die, djs, dje
      integer :: is, ie, js, je, i, j
      !
      call get_di(info, dis, die)
      is = 1 + dis; ie = nx - die
      call get_dj(info, djs, dje)
      js = 1 + djs; je = ny - dje
      do i = is, ie
	 do j = js, je
	    if (imode == 0) then
	       field(i,j,:) = self%dcoef(:)*field(i,j,:)
	    else
	       field(i,j,:) = self%lvector(:,imode)*field(i,j,:)
	    end if
	 end do
      end do
      !
      return
   end subroutine modulate_NodeVModulation2
   !
   !
   !
   subroutine modulate_NodeVModulation3(self, imode, field, nxyt, nlev)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      integer, intent(in) :: imode, nxyt, nlev
      real(r_size), dimension(nxyt,nlev), intent(inout) :: field
      integer :: ixyt
      !
      do ixyt = 1, nxyt
	 if (imode == 0) then
	    field(ixyt,:) = self%dcoef(:)*field(ixyt,:)
	 else
	    field(ixyt,:) = self%lvector(:,imode)*field(ixyt,:)
	 end if
      end do
      !
      return
   end subroutine modulate_NodeVModulation3
   !
   !
   !
   subroutine multiply_NodeVModulation(self, info, imode, field, nx, ny, nlev)
      implicit none
      type(NodeVModulation), intent(inout) :: self
      integer, intent(in) :: imode, nx, ny, nlev
      type(NodeInfo), intent(in) :: info
      real(r_size), dimension(nx,ny,nlev), intent(inout) :: field
      integer :: dis, die, djs, dje
      integer :: is, ie, js, je
      !
      call get_di(info, dis, die)
      is = 1 + dis; ie = nx - die
      call get_dj(info, djs, dje)
      js = 1 + djs; je = ny - dje
      field(is:ie,js:je,:) = self%lcoef(imode)*field(is:ie,js:je,:)
      !
      return
   end subroutine multiply_NodeVModulation
   !
   !
   !
end module NodeVModulation_class
