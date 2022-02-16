module NodeObsValidField_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, nepe
   use NodeInfo_class
   use NodeObsSpaceField_class
   use NodeMPI
   implicit none
   !
   type NodeObsValidField
      character(len=10) :: name
      integer :: nobs, nsubobs, nx, ny
      integer, dimension(:,:), allocatable :: field
   end type NodeObsValidField
   !
   interface new
      module procedure new_NodeObsValidField
   end interface
   interface destroy
      module procedure destroy_NodeObsValidField
   end interface
   interface display
      module procedure display_NodeObsValidField
   end interface
   interface assignment(=)
      module procedure copy_NodeObsValidField
   end interface
   interface get_name
      module procedure get_name_NodeObsValidField
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsValidField
   end interface
   interface set_field
      module procedure set_field_NodeObsValidField1
      module procedure set_field_NodeObsValidField2
      module procedure set_field_NodeObsValidField3
      module procedure set_field_NodeObsValidField4
   end interface
   interface and
      module procedure and_NodeObsValidField1
      module procedure and_NodeObsValidField2
      module procedure and_NodeObsValidField3
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsValidField
   end interface
   interface broadcast_ensvalid
      module procedure broadcast_ensvalid_NodeObsValidField1
      module procedure broadcast_ensvalid_NodeObsValidField2
   end interface
   interface allreduce_ensvalid
      module procedure allreduce_ensvalid_NodeObsValidField1
      module procedure allreduce_ensvalid_NodeObsValidField2
      module procedure allreduce_ensvalid_NodeObsValidField3
   end interface
   !
contains
   !
   !
   !
   subroutine new_NodeObsValidField(self, obsspace)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer :: iobs
      !
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      self%nsubobs = obsspace%nsubobs
      self%nx = obsspace%nx
      self%ny = obsspace%ny
      if (self%nobs > 0) then
	 allocate(self%field(self%nobs,self%nsubobs))
	 if (trim(obsspace%obstype) == 'RAD') then
	    self%field(:,:) = 0
	    do iobs = 1, self%nobs
	       self%field(iobs,1:obsspace%rad%nchannel(iobs)) = 1
	    end do
	 else
	    self%field(:,:) = 1
	 end if
      end if
      !
      return
   end subroutine new_NodeObsValidField
   !
   !
   !
   subroutine destroy_NodeObsValidField(self)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      !
      if (allocated(self%field)) deallocate(self%field)
      !
      return
   end subroutine destroy_NodeObsValidField
   !
   !
   !
   subroutine display_NodeObsValidField(self)
      implicit none
      type(NodeObsValidField), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs, self%nsubobs
      !
      return
   end subroutine display_NodeObsValidField
   !
   !
   !
   subroutine copy_NodeObsValidField(self, object)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      type(NodeObsValidField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = object%field(:,:)
      !
      return
   end subroutine copy_NodeObsValidField
   !
   !
   !
   subroutine get_name_NodeObsValidField(self, name)
      implicit none
      type(NodeObsValidField), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsValidField
   !
   !
   !
   subroutine get_nobs_NodeObsValidField(self, nobs)
      implicit none
      type(NodeObsValidField), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsValidField
   !
   !
   !
   subroutine set_field_NodeObsValidField1(self, const)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: const
      !
      if (self%nobs == 0) return
      self%field(:,:) = const
      !
      return
   end subroutine set_field_NodeObsValidField1
   !
   !
   !
   subroutine set_field_NodeObsValidField2(self, object)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      type(NodeObsValidField), intent(in) :: object
      !
      if (self%nobs == 0) return
      self%field(:,:) = object%field(:,:)
      !
      return
   end subroutine set_field_NodeObsValidField2
   !
   !
   !
   subroutine set_field_NodeObsValidField3(self, processed, k2ijt, obsspace, const, nxyt)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer, intent(in) :: const
      integer :: ixyt, i, j, it, iobs, iobs1, iobs2, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = const
	 end do
      end do
      !
      return
   end subroutine set_field_NodeObsValidField3
   !
   !
   !
   subroutine set_field_NodeObsValidField4(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: object
      integer :: ixyt, i, j, it, iobs, iobs1, iobs2, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    self%field(iobs,1:nsubobs) = object%field(iobs,1:nsubobs)
	 end do
      end do
      !
      return
   end subroutine set_field_NodeObsValidField4
   !
   !
   !
   subroutine and_NodeObsValidField1(self, object)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      type(NodeObsValidField), intent(in) :: object
      integer :: iobs, j
      !
      do iobs = 1, self%nobs
         do j = 1, self%nsubobs
            if (self%field(iobs,j) == 1 .and. object%field(iobs,j) == 1) then
	       self%field(iobs,j) = 1
	    else
	       self%field(iobs,j) = 0
	    end if
	 end do
      end do
      !
      return
   end subroutine and_NodeObsValidField1
   !
   !
   !
   subroutine and_NodeObsValidField2(self, i, j, obsspace, object)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: object
      integer :: it, iobs1, iobs2, iobs, jobs, nsubobs
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = self%nsubobs
	    end if
	    do jobs = 1, nsubobs
	       if (self%field(iobs,jobs) == 1 .and. object%field(iobs,jobs) == 1) then
		  self%field(iobs,jobs) = 1
	       else
		  self%field(iobs,jobs) = 0
	       end if
	    end do
	 end do
      end do
	 
      !
      return
   end subroutine and_NodeObsValidField2
   !
   !
   !
   subroutine and_NodeObsValidField3(self, processed, k2ijt, obsspace, object, nxyt)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: object
      integer :: ixyt, i, j, it, iobs1, iobs2, iobs, jobs, nsubobs
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 do iobs = iobs1, iobs2
	    if (trim(obsspace%obstype) == 'RAD') then
	       nsubobs = obsspace%rad%nchannel(iobs)
	    else
	       nsubobs = obsspace%nsubobs
	    end if
	    do jobs = 1, nsubobs
	       if (self%field(iobs,jobs) == 1 .and. object%field(iobs,jobs) == 1) then
		  self%field(iobs,jobs) = 1
	       else
		  self%field(iobs,jobs) = 0
	       end if
	    end do
	 end do
      end do
      !
      return
   end subroutine and_NodeObsValidField3
   !
   !
   !
   subroutine gather_obs_NodeObsValidField(self, info, obsspace, global_obsspace, global_object)
      ! gather all obs from all nodes into one big obs. This big obs exists in all nodes.
      implicit none
      type(NodeObsValidField), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceField), intent(in) :: obsspace, global_obsspace
      type(NodeObsValidField), intent(inout) :: global_object
      !
      integer :: istat(MPI_STATUS_SIZE)
      integer :: nproc, myid, nx0, ny0, nx, ny, nt, nobs, nmax
      integer :: id, ierror, iobs1, iobs2, iobs3, iobs4, i, j, it
      integer :: idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: myidx, myidy, myide, nxpe, nype, nepe, ide, source, destination
      integer, dimension(:,:,:), allocatable :: iobs
      integer, dimension(:,:), allocatable :: field
      !
      call get_myid(info, myid, nproc)
      call get_myidx(info, myidx); call get_nxpe(info, nxpe)
      call get_myidy(info, myidy); call get_nype(info, nype)
      call get_myide(info, myide); call get_nepe(info, nepe)
      nx0 = global_obsspace%nx
      ny0 = global_obsspace%ny
      nx = obsspace%nx
      ny = obsspace%ny
      nt = obsspace%nt
      !
      allocate(iobs(nx0,ny0,nt))
      iobs(1:nx,1:ny,:) = obsspace%iobs
      call int_gather(info, 'xy', iobs, 0, nx0, ny0, nt)
      if (myidx == 0 .and. myidy == 0) then
         if (myide == 0) then
            do ide = 1, nepe-1
               destination = ide*nxpe*nype
               call MPI_SEND(iobs, nx0*ny0*nt, MPI_INTEGER, destination, 0, MPI_COMM_WORLD, ierror)
            end do
         else
            source = 0
            call MPI_RECV(iobs, nx0*ny0*nt, MPI_INTEGER, source, 0, MPI_COMM_WORLD, istat, ierror)
         end if
      end if
      !
      if (myidx == 0 .and. myidy == 0) then
	 nmax = 0
	 do id = 0, nxpe*nype-1
	    call get_info(info, id, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    nmax = max(nmax, sum(global_obsspace%mobs(is+dis:ie-die,js+djs:je-dje,:)))
	 end do
	 allocate(field(nmax,self%nsubobs))
	 !
	 do id = 0, nxpe*nype-1
	    source = myide*nxpe*nype + id
	    if (id == 0) then
	       nobs = self%nobs
	    else
	       call MPI_RECV(nobs, 1, MPI_INTEGER, source, 0, MPI_COMM_WORLD, istat, ierror)
	    end if
	    if (nobs == 0) cycle
	    if (id == 0) then
	       field(1:nobs,:) = self%field
	    else
	       call MPI_RECV(field(1:nobs,:), nobs*self%nsubobs, MPI_INTEGER, source, 0, MPI_COMM_WORLD, istat, ierror)
	    end if
	    !
	    call get_info(info, id, idx, idy, is, ie, js, je, dis, die, djs, dje)
	    do it = 1, nt
	       do j = js+djs, je-dje
		  do i = is+dis, ie-die
		     if (global_obsspace%mobs(i,j,it) > 0) then
			iobs1 = global_obsspace%iobs(i,j,it)
		        iobs2 = iobs1 + global_obsspace%mobs(i,j,it) - 1
			iobs3 = iobs(i,j,it)
			iobs4 = iobs3 + global_obsspace%mobs(i,j,it) - 1
			global_object%field(iobs1:iobs2,:) = field(iobs3:iobs4,:)
		     end if
		  end do
	       end do
	    end do
         end do
	 deallocate(field)
      else
	 destination = myide*nxpe*nype
	 call MPI_SEND(self%nobs, 1, MPI_INTEGER, destination, 0, MPI_COMM_WORLD, ierror)
	 if (self%nobs > 0) call MPI_SEND(self%field, self%nobs*self%nsubobs, MPI_INTEGER, destination, 0, MPI_COMM_WORLD, ierror)
      end if
      !
      if (myidx == 0 .and. myidy == 0) then
         do id = 1, nxpe*nype-1
            destination = myide*nxpe*nype + id
            call MPI_SEND(global_object%field, global_object%nobs*global_object%nsubobs, MPI_INTEGER, destination, 0, MPI_COMM_WORLD, ierror)
         end do
      else
         source = myide*nxpe*nype
         call MPI_RECV(global_object%field, global_object%nobs*global_object%nsubobs, MPI_INTEGER, source, 0, MPI_COMM_WORLD, istat, ierror)
      end if
      deallocate(iobs)
      !
      return
   end subroutine gather_obs_NodeObsValidField
   !
   !
   !
   subroutine broadcast_ensvalid_NodeObsValidField1(self, source)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: source
      !
      if (self%nobs == 0) return
      call int_broadcast2D('e', self%field, source, self%nobs, self%nsubobs)
      !
      return
   end subroutine broadcast_ensvalid_NodeObsValidField1
   !
   !
   !
   subroutine broadcast_ensvalid_NodeObsValidField2(self, i, j, obsspace, source)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer, intent(in) :: source
      integer :: it, iobs1, iobs2
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 call int_broadcast2D('e', self%field(iobs1:iobs2,:), source, obsspace%mobs(i,j,it), self%nsubobs)
      end do
      !
      return
   end subroutine broadcast_ensvalid_NodeObsValidField2
   !
   !
   !
   subroutine allreduce_ensvalid_NodeObsValidField1(self)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      !
      if (self%nobs == 0) return
      call int_allreduce2D('e', self%field, self%nobs, self%nsubobs)
      where (self%field == nepe)
         self%field = 1
      else where
         self%field = 0
      end where
      !
      return
   end subroutine allreduce_ensvalid_NodeObsValidField1
   !
   !
   !
   subroutine allreduce_ensvalid_NodeObsValidField2(self, i, j, obsspace)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: i, j
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer :: it, iobs1, iobs2
      !
      do it = 1, obsspace%nt
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 call int_allreduce2D('e', self%field(iobs1:iobs2,:), obsspace%mobs(i,j,it), self%nsubobs)
	 where (self%field(iobs1:iobs2,:) == nepe)
	    self%field(iobs1:iobs2,:) = 1
	 else where
	    self%field(iobs1:iobs2,:) = 0
	 end where
      end do
      !
      return
   end subroutine allreduce_ensvalid_NodeObsValidField2
   !
   !
   !
   subroutine allreduce_ensvalid_NodeObsValidField3(self, processed, k2ijt, obsspace, nxyt)
      implicit none
      type(NodeObsValidField), intent(inout) :: self
      integer, intent(in) :: nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeObsSpaceField), intent(in) :: obsspace
      integer :: ixyt, i, j, it, iobs1, iobs2
      !
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 if (processed(i,j)) cycle
         if (obsspace%mobs(i,j,it) == 0) cycle
	 iobs1 = obsspace%iobs(i,j,it)
	 iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	 call int_allreduce2D('e', self%field(iobs1:iobs2,:), obsspace%mobs(i,j,it), self%nsubobs)
	 where (self%field(iobs1:iobs2,:) == nepe)
	    self%field(iobs1:iobs2,:) = 1
	 else where
	    self%field(iobs1:iobs2,:) = 0
	 end where
      end do
      !
      return
   end subroutine allreduce_ensvalid_NodeObsValidField3
   !
   !
   !
end module NodeObsValidField_class
