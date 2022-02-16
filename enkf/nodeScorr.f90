module NodeScorr_class
! Author: Le Duc
! Created date: 22 Dec 2014
   use variable, only : r_size, r_sngl
   use NodeInfo_class
   use NodeMPI
   implicit none
   !
   type NodeScorr
      character(1) :: header
      integer :: ncorr, nrank, nlev
      real(r_size), dimension(:,:), allocatable :: scale, eigval
      real(r_size), dimension(:,:,:), allocatable :: eigvec
      real(r_size), dimension(:,:,:), allocatable :: xtmp, ytmp, ztmp
   end type NodeScorr
   !
   interface new
      module procedure new_NodeScorr
   end interface
   interface destroy
      module procedure destroy_NodeScorr
   end interface
   interface display
      module procedure display_NodeScorr
   end interface
   interface get_nlev
      module procedure get_nlev_NodeScorr
   end interface
   interface set_Scorr
      module procedure set_Scorr_NodeScorr
   end interface
   interface apply_Scorr
      module procedure apply_Scorr_NodeScorr
   end interface
   !
contains
   !
   subroutine new_NodeScorr(self, ncorr, nlev, header)
      implicit none
      type(NodeScorr), intent(inout) :: self
      integer, intent(in) :: ncorr, nlev
      character(1), intent(in) :: header
      !
      self%ncorr = ncorr
      self%nlev = nlev
      self%header = header
      !
      return
   end subroutine new_NodeScorr
   !
   !
   !
   subroutine destroy_NodeScorr(self)
      implicit none
      type(NodeScorr), intent(inout) :: self
      !
      if (allocated(self%eigval)) deallocate(self%eigval)
      if (allocated(self%eigvec)) deallocate(self%eigvec)
      if (allocated(self%xtmp)) deallocate(self%xtmp)
      if (allocated(self%ytmp)) deallocate(self%ytmp)
      if (allocated(self%ztmp)) deallocate(self%ztmp)
      !
      return
   end subroutine destroy_NodeScorr
   !
   !
   !
   subroutine display_NodeScorr(self)
      implicit none
      type(NodeScorr), intent(in) :: self
      !
      print*, 'Dimension: ', self%ncorr, self%nrank, self%nlev
      !
      return
   end subroutine display_NodeScorr
   !
   !
   !
   subroutine get_nlev_NodeScorr(self, nlev)
      implicit none
      type(NodeScorr), intent(in) :: self
      integer, intent(out) :: nlev
      !
      nlev = self%nlev
      !
      return
   end subroutine get_nlev_NodeScorr
   !
   !
   !
   subroutine set_Scorr_NodeScorr(self, info, mode)
      implicit none
      type(NodeScorr), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      character(*), intent(in) :: mode
      integer :: nproc, myid, ncorr, nrank, nlev
      real(r_sngl), dimension(:,:), allocatable :: scale, eigval
      real(r_sngl), dimension(:,:,:), allocatable :: eigvec
      !
      ncorr = self%ncorr
      nlev = self%nlev
      call get_myid(info, myid, nproc)
      if (myid == 0) then
	 allocate(scale(ncorr,nlev),eigval(ncorr,nlev))
         allocate(eigvec(ncorr,ncorr,nlev))
	 open(90, file=self%header//trim(mode), form='unformatted')
	 read(90) nrank
	 read(90) scale(:,:), eigval(1:nrank,:), eigvec(:,1:nrank,:)
	 !print*, 'Scale: ', ncorr, nlev, minval(scale), maxval(scale)
	 close(90)
      end if
      !
      call int_broadcast0D('all', nrank, 0)
      self%nrank = nrank
      if (allocated(self%scale)) deallocate(self%scale)
      if (allocated(self%eigval)) deallocate(self%eigval)
      if (allocated(self%eigvec)) deallocate(self%eigvec)
      allocate(self%scale(ncorr,nlev), self%eigval(nrank,nlev), self%eigvec(ncorr,nrank,nlev))
      if (myid == 0) then
         self%scale(:,:) = scale(:,:)
	 self%eigval(:,:) = eigval(1:nrank,:)
	 self%eigvec(:,:,:) = eigvec(:,1:nrank,:)
	 deallocate(scale, eigval, eigvec)
      end if
      call broadcast2D('all', self%scale, 0, ncorr, nlev)
      call broadcast2D('all', self%eigval, 0, nrank, nlev)
      call broadcast3D('all', self%eigvec, 0, ncorr, nrank, nlev)
      !
      return
   end subroutine set_Scorr_NodeScorr
   !
   !
   !
   subroutine apply_Scorr_NodeScorr(self, adjoint, info, field, nx, ny, nz0)
      ! no temporal localization
      implicit none
      type(NodeScorr), intent(inout) :: self
      logical, intent(in) :: adjoint
      integer, intent(in) :: nx, ny, nz0
      type(NodeInfo), intent(in) :: info
      real(r_size), dimension(nx,ny,nz0), intent(inout) :: field
      integer :: nrank, nx0, ny0, dis, die, djs, dje
      integer :: is0, ie0, js0, je0, is, ie, js, je, i, j, k
      real(r_size), dimension(self%nrank) :: newu
      !
      nrank = self%nrank
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      call get_xindex(info, 1, is0)
      call get_xindex(info, nx, ie0)
      call get_di(info, dis, die)
      is0 = is0 + dis; ie0 = ie0 - die
      is = 1 + dis; ie = nx - die
      call get_yindex(info, 1, js0)
      call get_yindex(info, ny, je0)
      call get_dj(info, djs, dje)
      js0 = js0 + djs; je0 = je0 - dje
      js = 1 + djs; je = ny - dje
      !
      ! x-dim Sb
      if (self%header == 'x') then
	 if (.not. allocated(self%xtmp)) allocate(self%xtmp(nx0,ny,nz0))
	 self%xtmp = 0.d0
	 do k = 1, nz0
	    do j = js, je
	       if (adjoint) field(is:ie,j,k) = field(is:ie,j,k)/self%scale(is0:ie0,k)
	       do i = 1, nrank
		  newu(i) = sum(self%eigvec(is0:ie0,i,k)*field(is:ie,j,k))
		  newu(i) = self%eigval(i,k)*newu(i)
	       end do
	       do i = 1, nrank
		  self%xtmp(:,j,k) = self%xtmp(:,j,k) + newu(i)*self%eigvec(:,i,k)
	       end do
	       if (.not. adjoint) self%xtmp(:,j,k) = self%xtmp(:,j,k)/self%scale(:,k)
	    end do
	 end do
	 call allreduce3D(self%header, self%xtmp, nx0, ny, nz0)
	 field(is:ie,js:je,:) = self%xtmp(is0:ie0,js:je,:)
      !
      ! y-dim Sb
      else if (self%header == 'y') then
	 if (.not. allocated(self%ytmp)) allocate(self%ytmp(nx,ny0,nz0))
         self%ytmp = 0.d0
	 do k = 1, nz0
	    do i = is, ie
	       if (adjoint) field(i,js:je,k) = field(i,js:je,k)/self%scale(js0:je0,k)
	       do j = 1, nrank
		  newu(j) = sum(self%eigvec(js0:je0,j,k)*field(i,js:je,k))
		  newu(j) = self%eigval(j,k)*newu(j)
	       end do
	       do j = 1, nrank
		  self%ytmp(i,:,k) = self%ytmp(i,:,k) + newu(j)*self%eigvec(:,j,k)
	       end do
	       if (.not. adjoint) self%ytmp(i,:,k) = self%ytmp(i,:,k)/self%scale(:,k)
	    end do
	 end do
	 call allreduce3D(self%header, self%ytmp, nx, ny0, nz0)
	 field(is:ie,js:je,:) = self%ytmp(is:ie,js0:je0,:)
      !
      ! z-dim Sb
      else if (self%header == 'z') then
	 if (.not. allocated(self%ztmp)) allocate(self%ztmp(nx,ny,nz0))
         self%ztmp = 0.d0
	 do i = is, ie
	    do j = js, je
	       if (adjoint) field(i,j,:) = field(i,j,:)/self%scale(:,1)
	       do k = 1, nrank
		  newu(k) = sum(self%eigvec(:,k,1)*field(i,j,:))
		  newu(k) = self%eigval(k,1)*newu(k)
	       end do
	       do k = 1, nrank
		  self%ztmp(i,j,:) = self%ztmp(i,j,:) + newu(k)*self%eigvec(:,k,1)
	       end do
	       if (.not. adjoint) self%ztmp(i,j,:) = self%ztmp(i,j,:)/self%scale(:,1)
	    end do
	 end do
	 field(is:ie,js:je,:) = self%ztmp(is:ie,js:je,:)
      end if
      !
      !do k = 1, nz0
	 !do j = js, je
	    !do i = is, ie
	       !if (field(i,j,k) /= field(i,j,k)) x%field(i,j,k) = 0.d0
	    !end do
	 !end do
      !end do
      !
      return
   end subroutine apply_Scorr_NodeScorr
   !
   !
   !
end module NodeScorr_class
