module NodeObsSpaceFieldTC_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl, iflg_incremental, tcobs_format
   use NodeInfo_class
   use NodeMPI
   implicit none
   !
   ! Obs are stored as a linked list at each horizontal grid point. To access
   ! obs fast, all link lists are concatenated into a one-dimensional array. The
   ! iobs and mobs are the indices fo accessing: iobs is the index of each link
   ! list in this array (=0 if no link list) and mobs is the number of obs in each link list.
   ! We do not put all information like lon, lat, ... into a structure because these fields
   ! depend on observation types.
   type NodeObsSpaceFieldTC
      character(len=10) :: name, filename
      ! For TC Vital, nsubobs = 3 (tclon, tclat, tcpmsl)
      integer :: nobs, nsubobs, nx, ny, nt
      real(r_sngl), dimension(:), allocatable :: lon, lat
      real(r_size), dimension(:,:), allocatable :: obs, error
      integer, dimension(:,:,:), allocatable :: iobs, mobs
   end type NodeObsSpaceFieldTC
   !
   interface new
      module procedure new_NodeObsSpaceFieldTC
   end interface
   interface destroy
      module procedure destroy_NodeObsSpaceFieldTC
   end interface
   interface display
      module procedure display_NodeObsSpaceFieldTC
   end interface
   interface get_name
      module procedure get_name_NodeObsSpaceFieldTC
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsSpaceFieldTC
   end interface
   interface get_nsubobs
      module procedure get_nsubobs_NodeObsSpaceFieldTC
   end interface
   interface get_mobs
      module procedure get_mobs_NodeObsSpaceFieldTC1
      module procedure get_mobs_NodeObsSpaceFieldTC2
   end interface
   interface read_obs
      module procedure read_obs_NodeObsSpaceFieldTC
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsSpaceFieldTC
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsSpaceFieldTC
   end interface
   !
contains
   !
   subroutine new_NodeObsSpaceFieldTC(self, name, filename, nx, ny, nt)
      implicit none
      type(NodeObsSpaceFieldTC), intent(inout) :: self
      character(len=*), intent(in) :: name, filename
      integer, intent(in) :: nx, ny, nt
      !
      self%name = name
      self%filename = filename
      self%nx = nx
      self%ny = ny
      self%nt = nt
      allocate(self%iobs(nx,ny,nt), self%mobs(nx,ny,nt))
      self%iobs = 0
      self%mobs = 0
      self%nobs = 0
      self%nsubobs = 0
      !
      return
   end subroutine new_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine destroy_NodeObsSpaceFieldTC(self)
      implicit none
      type(NodeObsSpaceFieldTC), intent(inout) :: self
      integer :: iobs
      !
      if (allocated(self%iobs)) deallocate(self%iobs)
      if (allocated(self%mobs)) deallocate(self%mobs)
      if (allocated(self%lon)) deallocate(self%lon)
      if (allocated(self%lat)) deallocate(self%lat)
      if (allocated(self%obs)) deallocate(self%obs)
      if (allocated(self%error)) deallocate(self%error)
      !
      return
   end subroutine destroy_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine display_NodeObsSpaceFieldTC(self)
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      !
      print*, 'Observation name: ', self%name
      print*, 'File name: ', self%filename
      print*, 'Number: ', self%nobs, self%nsubobs
      !
      return
   end subroutine display_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine get_name_NodeObsSpaceFieldTC(self, name)
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine get_nobs_NodeObsSpaceFieldTC(self, nobs)
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine get_nsubobs_NodeObsSpaceFieldTC(self, nsubobs)
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      integer, intent(out) :: nsubobs
      !
      nsubobs = self%nsubobs
      !
      return
   end subroutine get_nsubobs_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldTC1(self, i, j, mobs)
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      integer, intent(in) :: i, j
      integer, intent(out) :: mobs
      !
      mobs = sum(self%mobs(i,j,:))*self%nsubobs
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldTC1
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldTC2(self, i, j, it, mobs)
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      integer, intent(in) :: i, j, it
      integer, intent(out) :: mobs
      !
      mobs = self%mobs(i,j,it)*self%nsubobs
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldTC2
   !
   !
   !
   subroutine txtread_nobs(self, nobs)
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      integer, intent(out) :: nobs
      !
      character(2) :: infile = '00'
      integer :: nx, ny, nt
      integer :: it, i, j, ierror
      real(kind=r_sngl) :: lat, lon  ! lat, lon must be single
      real(kind=r_size) :: obslon, obslat, tclon, tclat, tcpmsl, locerror, pmslerror
      type(PosGrid) :: pg
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      !
      nobs = 0
      LOOP_SLOT: do it = 1, nt
         write(infile(1:2),'(I2.2)') it
	 open(90, file=trim(self%filename)//infile, iostat=ierror)
         if (ierror > 0) cycle
	 !
	 LOOP_RECORD: do
	    read(90,'(7f8.3)',end=10) obslon, obslat, tclon, tclat, tcpmsl, locerror, pmslerror
	    goto 11
 10	    exit
 11	    continue
	    !
	    ! convert lat,lon to grid coordinates
	    lon = real(obslon, r_sngl)
	    lat = real(obslat, r_sngl)
	    call convert_LatLon_to_GridPos(pg, lat, lon, 1)
	    i = pg%i1; j = pg%j1
	    if (i < 2 .or. i > nx-1 .or. j < 2 .or. j > ny-1) cycle
	    nobs = nobs + 1
         end do LOOP_RECORD
         close(90)
      end do LOOP_SLOT
      !
      return
   end subroutine txtread_nobs
   !
   !
   !
   subroutine txtread_obs(self, myid)
      ! self is a global ObsSpace
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos
      implicit none
      type(NodeObsSpaceFieldTC), intent(inout) :: self
      integer, intent(in) :: myid
      !
      character(2) :: infile = '00'
      integer :: nx, ny, nt, nobs, nsubobs, ierror
      integer :: it, i, j, iobs, iobs1, iobs2, m, itmp, jtmp
      real(kind=r_sngl) :: lat, lon
      real(r_size) :: obslon, obslat, tclon, tclat, tcpmsl, locerror, pmslerror
      real(r_size), dimension(3) :: obs, error
      type(PosGrid) :: pg
      integer, dimension(:), allocatable :: i0, j0
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      !
      if (myid == 0) call txtread_nobs(self, nobs)
      call int_broadcast0D('all', nobs, 0)
      self%nobs = nobs
      nsubobs = 3
      self%nsubobs = nsubobs
      if (nobs == 0) return
      allocate(self%lon(nobs), self%lat(nobs))
      allocate(self%obs(nobs,nsubobs), self%error(nobs,nsubobs))
      !
      if (myid == 0) then
	 allocate(i0(nobs), j0(nobs))
	 !
	 iobs = 0
	 LOOP_SLOT: do it = 1, nt
	    write(infile(1:2),'(I2.2)') it
	    open(90, file=trim(self%filename)//infile, iostat=ierror)
	    if (ierror > 0) cycle
	    !
	    LOOP_RECORD: do
	       read(90,'(7f8.3)',end=10) obslon, obslat, tclon, tclat, tcpmsl, locerror, pmslerror
	       goto 11
 10	       exit
 11	       continue
	       !
	       ! convert lat,lon to grid coordinates
	       lon = real(obslon, r_sngl)
	       lat = real(obslat, r_sngl)
	       call convert_LatLon_to_GridPos(pg, lat, lon, 1)
	       i = pg%i1; j = pg%j1
	       if (i < 2 .or. i > nx-1 .or. j < 2 .or. j > ny-1) cycle
	       ! pressure: [hPa] -> [Pa]
	       tcpmsl = tcpmsl*100.0d0
	       pmslerror = pmslerror*100.0d0
	       !
	       iobs = iobs + 1
	       self%mobs(i,j,it) = self%mobs(i,j,it) + 1
	       i0(iobs) = i; j0(iobs) = j
	       self%lon(iobs) = lon; self%lat(iobs) = lat
	       self%obs(iobs,1) = tclon; self%error(iobs,1) = locerror
	       self%obs(iobs,2) = tclat; self%error(iobs,2) = locerror
	       self%obs(iobs,3) = tcpmsl; self%error(iobs,3) = pmslerror
	    end do LOOP_RECORD
	    close(90)
	 end do LOOP_SLOT
	 !
	 ! shell sorting
	 do it = 1, nt
	    iobs1 = sum(self%mobs(:,:,1:it-1))
	    if (it == 1) iobs1 = 0
	    iobs2 = sum(self%mobs(:,:,1:it))
	    if (iobs2 > iobs1) then
	       m = iobs2-iobs1
	       do while (m > 1)
		  m = (m+2)/3
		  do i = m+1, iobs2-iobs1
		     do j = i, m+1, -m
			if (j0(iobs1+j-m) > j0(iobs1+j) .or. &
			 & (j0(iobs1+j-m) == j0(iobs1+j) .and. i0(iobs1+j-m) > i0(iobs1+j))) then
			   itmp = i0(iobs1+j); i0(iobs1+j) = i0(iobs1+j-m); i0(iobs1+j-m) = itmp
			   jtmp = j0(iobs1+j); j0(iobs1+j) = j0(iobs1+j-m); j0(iobs1+j-m) = jtmp
			   lon = self%lon(iobs1+j); self%lon(iobs1+j) = self%lon(iobs1+j-m); self%lon(iobs1+j-m) = lon
			   lat = self%lat(iobs1+j); self%lat(iobs1+j) = self%lat(iobs1+j-m); self%lat(iobs1+j-m) = lat
			   obs(:) = self%obs(iobs1+j,:); self%obs(iobs1+j,:) = self%obs(iobs1+j-m,:); self%obs(iobs1+j-m,:) = obs(:)
			   error(:) = self%error(iobs1+j,:); self%error(iobs1+j,:) = self%error(iobs1+j-m,:); self%error(iobs1+j-m,:) = error(:)
			end if
		     end do
		  end do
	       end do
	    end if
	 end do
	 !
	 iobs = 1
	 do it = 1, nt
	    do j = 1, ny
	       do i = 1, nx
	          if (self%mobs(i,j,it) == 0) then
		     self%iobs(i,j,it) = 0
		  else
		     self%iobs(i,j,it) = iobs
		  end if
		  iobs = iobs + self%mobs(i,j,it)
	       end do
	    end do
	 end do
	 deallocate(i0, j0)
      end if
      !
      call int_broadcast3D('all', self%iobs, 0, nx, ny, nt)
      call int_broadcast3D('all', self%mobs, 0, nx, ny, nt)
      call real_broadcast1D('all', self%lon, 0, nobs)
      call real_broadcast1D('all', self%lat, 0, nobs)
      call broadcast2D('all', self%obs, 0, nobs, nsubobs)
      call broadcast2D('all', self%error, 0, nobs, nsubobs)
      !
      return
   end subroutine txtread_obs
   !
   !
   !
   subroutine read_obs_NodeObsSpaceFieldTC(self, myid)
      implicit none
      type(NodeObsSpaceFieldTC), intent(inout) :: self
      integer, intent(in) :: myid
      !
      if (trim(tcobs_format) == 'txt') then
	 call txtread_obs(self, myid)
      end if
      !
      return
   end subroutine read_obs_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine scatter_obs_NodeObsSpaceFieldTC(self, info, local_object)
      ! scatter a big obsspace from Node 0 to all nodes.
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldTC), intent(inout) :: local_object
      !
      integer :: myid, myidx, myidy, nproc, nxpe, nype, nx0, ny0, nx, ny, nt, nobs, nsubobs
      integer :: kobs, iobs1, iobs2, iobs3, iobs4, i, j, it
      integer :: is, ie, js, je, dis, die, djs, dje
      integer, dimension(:,:,:), allocatable :: iobs, mobs
      !
      call get_nxpe(info, nxpe); call get_nype(info, nype)
      call get_myidx(info, myidx); call get_myidy(info, myidy)
      myid = myidy*nxpe+myidx
      nx0 = self%nx
      ny0 = self%ny
      nx = local_object%nx
      ny = local_object%ny
      nt = local_object%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      if (info%indexmode == 0 .or. info%indexmode == 1) then ! EnSRF
         is = 1; ie = nx - 1
         js = 1; je = ny - 1
      else
         is = 1 + dis; ie = nx - die
         js = 1 + djs; je = ny - dje
      end if
      allocate(iobs(nx0,ny0,nt), mobs(nx0,ny0,nt))
      !
      if (myid == 0) then
	 iobs(:,:,:) = self%iobs
         mobs(:,:,:) = self%mobs
      end if
      call int_scatter(info, 'xy', iobs, 0, nx0, ny0, nt)
      call int_scatter(info, 'xy', mobs, 0, nx0, ny0, nt)
      !
      nobs = sum(mobs(is:ie,js:je,:))
      local_object%nobs = nobs
      nsubobs = self%nsubobs
      local_object%nsubobs = nsubobs
      if (nobs > 0) then
	 local_object%mobs(is:ie,js:je,:) = mobs(is:ie,js:je,:)
	 allocate(local_object%lon(nobs), local_object%lat(nobs))
	 allocate(local_object%obs(nobs,nsubobs), local_object%error(nobs,nsubobs))
	 !
	 kobs = 1
	 do it = 1, nt
	    do j = js, je
	       do i = is, ie
	          if (mobs(i,j,it) > 0) then
		     local_object%iobs(i,j,it) = kobs
		     iobs1 = kobs
		     iobs2 = iobs1 + mobs(i,j,it) - 1
		     iobs3 = iobs(i,j,it)
		     iobs4 = iobs3 + mobs(i,j,it) - 1
		     local_object%lon(iobs1:iobs2) = self%lon(iobs3:iobs4)
		     local_object%lat(iobs1:iobs2) = self%lat(iobs3:iobs4)
		     local_object%obs(iobs1:iobs2,:) = self%obs(iobs3:iobs4,:)
		     local_object%error(iobs1:iobs2,:) = self%error(iobs3:iobs4,:)
		  else
		     local_object%iobs(i,j,it) = 0
		  end if
		  kobs = kobs + mobs(i,j,it)
	       end do
	    end do
	 end do
      end if
      deallocate(iobs, mobs)
      !
      return
   end subroutine scatter_obs_NodeObsSpaceFieldTC
   !
   !
   !
   subroutine gather_obs_NodeObsSpaceFieldTC(self, info, global_object)
      ! this subroutine is only for completion and we realy do not use this
      implicit none
      type(NodeObsSpaceFieldTC), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldTC), intent(inout) :: global_object
      !
      integer :: nproc, myid, nx0, ny0, nx, ny, nt, nobs, i, j, it, iobs
      !
      call get_myid(info, myid, nproc)
      nx0 = global_object%nx
      ny0 = global_object%ny
      nx = self%nx
      ny = self%ny
      nt = self%nt
      !
      global_object%mobs(1:nx,1:ny,:) = self%mobs
      call int_gather(info, 'xy', global_object%mobs, 0, nx0, ny0, nt)
      call int_broadcast3D('all', global_object%mobs, 0, nx0, ny0, nt)
      if (myid == 0) nobs = sum(global_object%mobs)
      call int_broadcast0D('all', nobs, 0)
      global_object%nobs = nobs
      !
      if (myid == 0) then
	 iobs = 1
	 do it = 1, nt
	    do j = 1, ny0
	       do i = 1, nx0
		  if (global_object%mobs(i,j,it) == 0) then
		     global_object%iobs(i,j,it) = 0
		  else
		     global_object%iobs(i,j,it) = iobs
		  end if
		  iobs = iobs + global_object%mobs(i,j,it)
	       end do
	    end do
	 end do
      end if
      call int_broadcast3D('all', global_object%iobs, 0, nx0, ny0, nt)
      !
      return
   end subroutine gather_obs_NodeObsSpaceFieldTC
   !
   !
   !
end module NodeObsSpaceFieldTC_class
