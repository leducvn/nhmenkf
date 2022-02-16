module NodeHFieldTC_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl
   use interpolate, only : PosGrid, convert_LatLon_to_GridPos, convert_GridPos_offset
   use NodeInfo_class
   use NodeObsSpaceFieldTC_class
   use NodeObsValidField_class
   use NodeObsField_class
   use NodeObsControl_class
   use NodeControl_class
   use NodeProfileControl_class
   use NodeMPI
   implicit none
   !
   type NodeHFieldTC
      character(len=10) :: name
      integer :: nobs
      type(PosGrid), dimension(:), allocatable :: pg
      ! For ensemble linear tangent H
      integer :: ne
      integer, dimension(:), allocatable :: nrank
      real(r_size), dimension(:,:), allocatable :: eigval
      real(r_size), dimension(:,:,:), allocatable :: Y
      real(r_size), dimension(:,:,:,:,:), allocatable :: X, eigvec
   end type NodeHFieldTC
   !
   interface new
      module procedure new_NodeHFieldTC
   end interface
   interface destroy
      module procedure destroy_NodeHFieldTC
   end interface
   interface display
      module procedure display_NodeHFieldTC
   end interface
   interface get_name
      module procedure get_name_NodeHFieldTC
   end interface
   interface get_nobs
      module procedure get_nobs_NodeHFieldTC
   end interface
   interface get_xyloc
      module procedure get_xyloc_NodeHFieldTC1
      module procedure get_xyloc_NodeHFieldTC2
   end interface
   interface apply_H
      module procedure apply_H_NodeHFieldTC1
      module procedure apply_H_NodeHFieldTC2
      module procedure apply_H_NodeHFieldTC3
   end interface
   interface initialize_DH
      module procedure initialize_DH_NodeHFieldTC
   end interface
   interface apply_DH
      module procedure apply_DH_NodeHFieldTC
   end interface
   interface apply_DHT
      module procedure apply_DHT_NodeHFieldTC
   end interface
   !
contains
   !
   subroutine new_NodeHFieldTC(self, info, obsspace)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      self%name = obsspace%name
      self%nobs = obsspace%nobs
      if (self%nobs > 0) then
         allocate(self%pg(self%nobs))
         do iobs = 1, self%nobs
            call convert_LatLon_to_GridPos(self%pg(iobs), obsspace%lat(iobs), obsspace%lon(iobs), 1)
         end do
      end if
      !
      return
   end subroutine new_NodeHFieldTC
   !
   !
   !
   subroutine destroy_NodeHFieldTC(self)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      !
      if (allocated(self%pg)) deallocate(self%pg)
      if (allocated(self%nrank)) then
	 deallocate(self%Y, self%X, self%nrank, self%eigval, self%eigvec)
      end if
      !
      return
   end subroutine destroy_NodeHFieldTC
   !
   !
   !
   subroutine display_NodeHFieldTC(self)
      implicit none
      type(NodeHFieldTC), intent(in) :: self
      !
      print*, 'Field name: ', self%name
      print*, 'Dimension: ', self%nobs
      !
      return
   end subroutine display_NodeHFieldTC
   !
   !
   !
   subroutine get_name_NodeHFieldTC(self, name)
      implicit none
      type(NodeHFieldTC), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeHFieldTC
   !
   !
   !
   subroutine get_nobs_NodeHFieldTC(self, nobs)
      implicit none
      type(NodeHFieldTC), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeHFieldTC
   !
   !
   !
   subroutine get_xyloc_NodeHFieldTC1(self, xyloc)
      implicit none
      type(NodeHFieldTC), intent(in) :: self
      type(NodeObsField), intent(inout) :: xyloc
      integer :: iobs
      !
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px
	 xyloc%field(iobs,2) = self%pg(iobs)%py
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldTC1
   !
   !
   !
   subroutine get_xyloc_NodeHFieldTC2(self, info, xyloc)
      implicit none
      type(NodeHFieldTC), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsField), intent(inout) :: xyloc
      integer :: is, js, iobs
      !
      call get_xindex(info, 1, is)
      call get_yindex(info, 1, js)
      do iobs = 1, self%nobs
         xyloc%field(iobs,1) = self%pg(iobs)%px - is + 1
	 xyloc%field(iobs,2) = self%pg(iobs)%py - js + 1
      end do
      !
      return
   end subroutine get_xyloc_NodeHFieldTC2
   !
   !
   !
   subroutine detect_tcproc(info, tcid, ntcproc, tcproc)
      implicit none
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: tcid, ntcproc
      integer, dimension(ntcproc), intent(inout) :: tcproc
      integer :: myid, nproc, myide, nxpe, nype, ierror, iproc, itcproc
      integer, dimension(:), allocatable :: tcidr
      !
      call get_myid(info, myid, nproc); call get_myide(info, myide)
      call get_nxpe(info, nxpe); call get_nype(info, nype)
      nproc = nxpe*nype
      if (myid == myide*nproc) allocate(tcidr(nproc))
      !
      call MPI_GATHER(tcid, 1, MPI_INTEGER, tcidr, 1, MPI_INTEGER, 0, MPI_COMM_XYWORLD, ierror)
      if (myid == myide*nproc) then
         itcproc = 0
         do iproc = 1, nproc
            if (tcidr(iproc) < 0) cycle
            itcproc = itcproc + 1
            tcproc(itcproc) = tcidr(iproc)
         end do
         deallocate(tcidr)
      end if
      call MPI_BCAST(tcproc, ntcproc, MPI_INTEGER, 0, MPI_COMM_XYWORLD, ierror)
      return
      !
   end subroutine detect_tcproc
   !
   !
   !
   subroutine interpolate_tc1(self, info, x, obsspace, valid, y)
      use variable, only : proj, resolution, xi, xj, xlon, xlat, slon, slat
      use nhmlib, only : ij2lonlat
      use tclib, only : dbuffer, search_max
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      integer :: myid, nproc, nx0, ny0, nx, ny, nt, nobs, ntcproc
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      integer :: tcid, itcproc, ibuffer, fi, fj, imin, imax, jmin, jmax, mobs, jobs
      integer, dimension(100) :: tcproc
      real(r_size) :: xmax, ymax, fmax, lonmax, latmax
      real(r_size), dimension(:), allocatable :: xx, yy
      real(r_size), dimension(:,:,:), allocatable :: data, pmsl
      !
      call get_myid(info, myid, nproc)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nx = x%nx
      ny = x%ny
      nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      ! allocate
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      call int_allreduce0D('xy', nobs)
      if (nobs == 0) return
      allocate(data(nx0,ny0,1))
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      if (nobs > 0) then
         ibuffer = nint(1000*dbuffer/resolution)
         allocate(xx(nx0), yy(ny0))
         allocate(pmsl(nx0,ny0,1))
         do i = 1, nx0
            xx(i) = i*1.d0
         end do
         do j = 1, ny0
            yy(j) = j*1.d0
         end do
      end if
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            ntcproc = 1
         else
            ntcproc = 0
         end if
         call int_allreduce0D('xy', ntcproc)
	 if (ntcproc == 0) cycle
	 if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            tcid = myid
         else
            tcid = -1
         end if
	 call detect_tcproc(info, tcid, ntcproc, tcproc(1:ntcproc))
	 !
	 ! get pmsl
         do itcproc = 1, ntcproc
            call get_field(x, 'pmsl', it, data(1:nx,1:ny,1:1), 1)
	    call gather(info, 'xy', data, tcproc(itcproc), nx0, ny0, 1)
            if (myid == tcproc(itcproc)) pmsl(:,:,1) = data(:,:,1)
         end do
         tcid = 0
         do itcproc = 1, ntcproc
            if (myid == tcproc(itcproc)) tcid = 1
         end do
	 if (tcid == 0) cycle
	 !
	 mobs = sum(obsspace%mobs(is:ie,js:je,it))
	 jobs = 0
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          fi = int(self%pg(iobs)%px); fj = int(self%pg(iobs)%py)
	          imin = max(1,fi-ibuffer); imax = min(nx0, fi+ibuffer)
	          jmin = max(1,fj-ibuffer); jmax = min(ny0, fj+ibuffer)
		  call search_max(xx(imin:imax), yy(jmin:jmax), -pmsl(imin:imax,jmin:jmax,1), imax-imin+1, jmax-jmin+1, xmax, ymax, fmax)
		  if (xmax < -1000.) then
		     valid%field(iobs,:) = 0
		     y%field(iobs,:) = 0.d0
		  else
		     valid%field(iobs,:) = 1
		     xmax = xmax
		     ymax = ny0 + 1 - ymax
		     call ij2lonlat(lonmax, latmax, xmax, ymax, proj, 1.d0*resolution, 1.d0*xi, 1.d0*xj, 1.d0*xlon, 1.d0*xlat, 1.d0*slon, 1.d0*slat)
		     y%field(iobs,1) = lonmax
		     y%field(iobs,2) = latmax
		     y%field(iobs,3) = -fmax
		  end if
		  !if (iobs == nobs) print*, 'TC:', myid, y%field(iobs,:)
		  jobs = jobs + 1
		  if (jobs == mobs) exit
	       end do
	       if (jobs == mobs) exit
	    end do
	    if (jobs == mobs) exit
	 end do
      end do
      if (nobs > 0) deallocate(xx, yy, pmsl)
      deallocate(data)
      !
      return
   end subroutine interpolate_tc1
   !
   !
   !
   subroutine initialize_Dtc(self, info, xpert, obsspace, valid, ypert, ne)
      use variable, only : resolution
      use tclib, only : tcloc, dbuffer, dloc
      use enkflib, only : mtx_eigen
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: xpert
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: ypert
      !
      integer :: myid, nproc, nx0, ny0, nx, ny, nt, nobs, nsubobs, ntcproc
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      integer :: nxtc0, nytc0, nrank0, nxtc, nytc, nrank, mobs, jobs
      integer :: tcid, itcproc, ibuffer, imember, fi, fj, gi, gj, imin, imax, jmin, jmax, i1, j1, k1, i2, j2, k2
      integer, dimension(100) :: tcproc
      real(r_size) :: Lloc, pi, pj
      real(r_size), dimension(:), allocatable :: eigval
      real(r_size), dimension(:,:), allocatable :: eigvec
      real(r_size), dimension(:,:,:), allocatable :: data, pmsl, cov, loc, obsloc
      !
      call get_myid(info, myid, nproc)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nx = xpert(1)%nx
      ny = xpert(1)%ny
      nt = xpert(1)%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      ! allocate
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      call int_allreduce0D('xy', nobs)
      if (nobs == 0) return
      allocate(data(nx0,ny0,1))
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      if (nobs > 0) then
         nsubobs = obsspace%nsubobs
	 self%ne = ne
         ibuffer = nint(1000.d0*dbuffer/resolution)
	 Lloc = 1000.d0*dloc/resolution
	 nxtc0 = 2*ibuffer + 1; nytc0 = 2*ibuffer + 1
	 nrank0 = nxtc0*nytc0
	 allocate(pmsl(nx0,ny0,1))
	 allocate(cov(nrank0,nrank0,nobs))
	 if (tcloc == 1) allocate(loc(nrank0,nrank0,nobs), obsloc(nxtc0,nytc0,nobs))
	 allocate(eigval(nrank0), eigvec(nrank0,nrank0))
	 allocate(self%Y(nsubobs,ne,nobs), self%X(nxtc0,nytc0,1,ne,nobs))
	 allocate(self%nrank(nobs), self%eigval(nrank0,nobs), self%eigvec(nxtc0,nytc0,1,nrank0,nobs))
	 cov(:,:,:) = 0.d0
	 self%Y(:,:,:) = 0.d0; self%X(:,:,:,:,:) = 0.d0
	 self%nrank(:) = 0; self%eigval(:,:) = 0.d0; self%eigvec(:,:,:,:,:) = 0.d0
      end if
      !
      ! calculate covariance
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            ntcproc = 1
         else
            ntcproc = 0
         end if
         call int_allreduce0D('xy', ntcproc)
	 if (ntcproc == 0) cycle
	 if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            tcid = myid
         else
            tcid = -1
         end if
	 call detect_tcproc(info, tcid, ntcproc, tcproc(1:ntcproc))
	 !
	 do imember = 1, ne
	    ! get pmsl
            do itcproc = 1, ntcproc
               call get_field(xpert(imember), 'pmsl', it, data(1:nx,1:ny,1:1), 1)
               call gather(info, 'xy', data, tcproc(itcproc), nx0, ny0, 1)
               if (myid == tcproc(itcproc)) pmsl(:,:,1) = data(:,:,1)
            end do
            tcid = 0
            do itcproc = 1, ntcproc
               if (myid == tcproc(itcproc)) tcid = 1
            end do
	    if (tcid == 0) cycle
	    !
	    mobs = sum(obsspace%mobs(is:ie,js:je,it))
	    jobs = 0
	    do j = js, je
	       do i = is, ie
	          iobs1 = obsspace%iobs(i,j,it)
	          iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	          do iobs = iobs1, iobs2
	             if (valid%field(iobs,1) == 0) then
	                jobs = jobs + 1
	                cycle
	             end if
	             pi = self%pg(iobs)%px; pj = self%pg(iobs)%py
	             fi = int(pi); fj = int(pj)
	             imin = max(1,fi-ibuffer); imax = min(nx0, fi+ibuffer)
	             jmin = max(1,fj-ibuffer); jmax = min(ny0, fj+ibuffer)
	             nxtc = imax - imin + 1; nytc = jmax - jmin + 1
	             do j1 = jmin, jmax
	                gj = j1 - jmin + 1
	                do i1 = imin, imax
	                   gi = i1 - imin + 1
		           k1 = (gj-1)*nytc + gi-1 + 1
		           if (tcloc == 1 .and. imember == 1) obsloc(gi,gj,iobs) = exp(-0.5d0*((pi-i1)**2+(pj-j1)**2)/Lloc**2)
		           self%X(gi,gj,1,imember,iobs) = pmsl(i1,j1,1)
	                   do j2 = jmin, jmax
		              fj = j2 - jmin + 1
		              do i2 = imin, imax
		                 fi = i2 - imin + 1
			         k2 = (fj-1)*nytc + fi-1 + 1
			         if (tcloc == 1 .and. imember == 1) loc(k1,k2,iobs) = exp(-0.5d0*((i2-i1)**2+(j2-j1)**2)/Lloc**2)
			         cov(k1,k2,iobs) = cov(k1,k2,iobs) + pmsl(i1,j1,1)*pmsl(i2,j2,1)
		              end do
		           end do
	                end do
	             end do
	             self%Y(:,imember,iobs) = ypert(imember)%tc%field(iobs,:)
	             if (tcloc == 1) self%X(1:nxtc,1:nytc,1,imember,iobs) = self%X(1:nxtc,1:nytc,1,imember,iobs)*obsloc(1:nxtc,1:nytc,iobs)
	             jobs = jobs + 1
	             if (jobs == mobs) exit
	          end do
	          if (jobs == mobs) exit
	       end do
	       if (jobs == mobs) exit
	    end do
	 end do
      end do
      if (nobs == 0) then
         deallocate(data)
         return
      end if
      call allreduce3D('e', cov, nrank0, nrank0, nobs)
      if (tcloc == 1) cov = cov*loc
      !
      ! eigen-decomposition
      do it = 1, nt
         mobs = sum(obsspace%mobs(is:ie,js:je,it))
	 jobs = 0
         do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          if (valid%field(iobs,1) == 0) then
	             jobs = jobs + 1
	             cycle
	          end if
	          pi = self%pg(iobs)%px; pj = self%pg(iobs)%py
	          fi = int(pi); fj = int(pj)
	          imin = max(1,fi-ibuffer); imax = min(nx0, fi+ibuffer)
	          jmin = max(1,fj-ibuffer); jmax = min(ny0, fj+ibuffer)
	          nxtc = imax - imin + 1; nytc = jmax - jmin + 1
	          nrank = nxtc*nytc
	          call mtx_eigen(1, nrank, cov(1:nrank,1:nrank,iobs), eigval(1:nrank), eigvec(1:nrank,1:nrank), self%nrank(iobs))
	          do k2 = 1, self%nrank(iobs)
	             if (eigval(k2) > 0.) then
	                self%eigval(k2,iobs) = 1.d0/eigval(k2)
	             else
	                self%eigval(k2,iobs) = 0.d0
	             end if
	             do gj = 1, nytc
	                do gi = 1, nxtc
		           k1 = (gj-1)*nytc + gi-1 + 1
		           self%eigvec(gi,gj,1,k2,iobs) = eigvec(k1,k2)
	                end do
	             end do
	          end do
	          !if (iobs == nobs) print*, 'TCINI:', myid, self%nrank(iobs), self%eigval(1,iobs), self%eigval(self%nrank(iobs),iobs)
	          jobs = jobs + 1
	          if (jobs == mobs) exit
	       end do
	       if (jobs == mobs) exit
	    end do
	    if (jobs == mobs) exit
         end do
      end do
      if (tcloc == 1) deallocate(obsloc, loc)
      deallocate(pmsl, cov, eigval, eigvec)
      !
      return
   end subroutine initialize_Dtc
   !
   !
   !
   subroutine interpolate_Dtc(self, info, x, obsspace, valid, y)
      use variable, only : resolution
      use tclib, only : dbuffer
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      integer :: myid, nproc, nx0, ny0, nx, ny, nt, ne, nobs, nsubobs, ntcproc
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      integer :: nxtc0, nytc0, nrank0, nxtc, nytc, mobs, jobs
      integer :: tcid, itcproc, ibuffer, imember, fi, fj, imin, imax, jmin, jmax, ivec
      integer, dimension(100) :: tcproc
      real(r_size), dimension(:), allocatable :: coordinate, weight
      real(r_size), dimension(:,:,:), allocatable :: data, dpmsl
      !
      call get_myid(info, myid, nproc)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nx = x%nx
      ny = x%ny
      nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      ! allocate
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      call int_allreduce0D('xy', nobs)
      if (nobs == 0) return
      allocate(data(nx0,ny0,1))
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      if (nobs > 0) then
         nsubobs = obsspace%nsubobs
         ne = self%ne
         ibuffer = nint(1000*dbuffer/resolution)
         nxtc0 = 2*ibuffer + 1; nytc0 = 2*ibuffer + 1
	 nrank0 = nxtc0*nytc0
         allocate(dpmsl(nx0,ny0,1))
         allocate(coordinate(nrank0), weight(ne))
      end if
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            ntcproc = 1
         else
            ntcproc = 0
         end if
         call int_allreduce0D('xy', ntcproc)
	 if (ntcproc == 0) cycle
	 if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            tcid = myid
         else
            tcid = -1
         end if
	 call detect_tcproc(info, tcid, ntcproc, tcproc(1:ntcproc))
	 !
	 ! get pmsl
         do itcproc = 1, ntcproc
            call get_field(x, 'pmsl', it, data(1:nx,1:ny,1:1), 1)
            call gather(info, 'xy', data, tcproc(itcproc), nx0, ny0, 1)
            if (myid == tcproc(itcproc)) dpmsl(:,:,1) = data(:,:,1)
         end do
         tcid = 0
         do itcproc = 1, ntcproc
            if (myid == tcproc(itcproc)) tcid = 1
         end do
	 if (tcid == 0) cycle
	 !
	 mobs = sum(obsspace%mobs(is:ie,js:je,it))
	 jobs = 0
	 do j = js, je
	    do i = is, ie
	       iobs1 = obsspace%iobs(i,j,it)
	       iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	       do iobs = iobs1, iobs2
	          fi = int(self%pg(iobs)%px); fj = int(self%pg(iobs)%py)
	          imin = max(1,fi-ibuffer); imax = min(nx0, fi+ibuffer)
	          jmin = max(1,fj-ibuffer); jmax = min(ny0, fj+ibuffer)
	          nxtc = imax - imin + 1; nytc = jmax - jmin + 1
	          do ivec = 1, self%nrank(iobs)
	             coordinate(ivec) = sum(self%eigvec(1:nxtc,1:nytc,1,ivec,iobs)*dpmsl(imin:imax,jmin:jmax,1))
	             coordinate(ivec) = coordinate(ivec)*self%eigval(ivec,iobs)
	          end do
	          dpmsl(1:nxtc,1:nytc,1) = 0.d0
	          do ivec = 1, self%nrank(iobs)
	             dpmsl(1:nxtc,1:nytc,1) = dpmsl(1:nxtc,1:nytc,1) + coordinate(ivec)*self%eigvec(1:nxtc,1:nytc,1,ivec,iobs)
	          end do
	          valid%field(iobs,:) = 1
	          y%field(iobs,:) = 0.d0
	          do imember = 1, ne
	             weight(imember) = sum(self%X(1:nxtc,1:nytc,1,imember,iobs)*dpmsl(1:nxtc,1:nytc,1))
	             y%field(iobs,:) = y%field(iobs,:) + weight(imember)*self%Y(:,imember,iobs)
	          end do
	          jobs = jobs + 1
	          if (jobs == mobs) exit
	       end do
	       if (jobs == mobs) exit
	    end do
	    if (jobs == mobs) exit
	 end do
      end do
      if (nobs == 0) then
         deallocate(data)
         return
      end if
      call allreduce2D('e', y%field, nobs, nsubobs)
      !print*, 'DTC:', myid, y%field(nobs,:)
      deallocate(coordinate, weight, dpmsl)
      !
      return
   end subroutine interpolate_Dtc
   !
   !
   !
   subroutine interpolate_DtcT(self, info, obsspace, valid, y, x)
      use variable, only : resolution
      use tclib, only : dbuffer
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      !
      integer :: myid, nproc, nx0, ny0, nx, ny, nt, ne, nobs, ntcproc
      integer :: i, j, it, iobs, iobs1, iobs2, dis, die, djs, dje, is, ie, js, je
      integer :: nxtc0, nytc0, nrank0, nxtc, nytc, mobs, jobs
      integer :: tcid, itcproc, ibuffer, imember, fi, fj, imin, imax, jmin, jmax, ivec
      integer, dimension(100) :: tcproc
      real(r_size), dimension(:), allocatable :: coordinate, weight
      real(r_size), dimension(:,:,:), allocatable :: data, dpmsl, dpmslall
      !
      call get_myid(info, myid, nproc)
      call get_nx0(info, nx0)
      call get_ny0(info, ny0)
      nx = x%nx
      ny = x%ny
      nt = x%nt
      call get_di(info, dis, die)
      call get_dj(info, djs, dje)
      is = 1 + dis; ie = nx - die
      js = 1 + djs; je = ny - dje
      !
      ! allocate
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      call int_allreduce0D('xy', nobs)
      if (nobs == 0) return
      allocate(data(nx0,ny0,1))
      nobs = sum(obsspace%mobs(is:ie,js:je,1:nt))
      if (nobs > 0) then
         ne = self%ne
         ibuffer = nint(1000*dbuffer/resolution)
         nxtc0 = 2*ibuffer + 1; nytc0 = 2*ibuffer + 1
	 nrank0 = nxtc0*nytc0
         allocate(dpmsl(nx0,ny0,1), dpmslall(nx0,ny0,1))
         allocate(coordinate(nrank0), weight(ne))
      end if
      !
      do it = 1, nt
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            ntcproc = 1
         else
            ntcproc = 0
         end if
         call int_allreduce0D('xy', ntcproc)
	 if (ntcproc == 0) cycle
	 if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            tcid = myid
         else
            tcid = -1
         end if
	 call detect_tcproc(info, tcid, ntcproc, tcproc(1:ntcproc))
	 !
         if (sum(obsspace%mobs(is:ie,js:je,it)) > 0) then
            dpmslall(:,:,1) = 0.d0
            mobs = sum(obsspace%mobs(is:ie,js:je,it))
            jobs = 0
            do j = js, je
	       do i = is, ie
	          iobs1 = obsspace%iobs(i,j,it)
	          iobs2 = obsspace%iobs(i,j,it) + obsspace%mobs(i,j,it) - 1
	          do iobs = iobs1, iobs2
	             if (valid%field(iobs,1) == 0) then
	                jobs = jobs + 1
	                cycle
	             end if
	             fi = int(self%pg(iobs)%px); fj = int(self%pg(iobs)%py)
	             imin = max(1,fi-ibuffer); imax = min(nx0, fi+ibuffer)
	             jmin = max(1,fj-ibuffer); jmax = min(ny0, fj+ibuffer)
	             nxtc = imax - imin + 1; nytc = jmax - jmin + 1
	             dpmsl(1:nxtc,1:nytc,1) = 0.d0
	             do imember = 1, ne
	                weight(imember) = sum(self%Y(1:2,imember,iobs)*y%field(iobs,1:2))
	                if (valid%field(iobs,3) == 1) weight(imember) = weight(imember) + self%Y(3,imember,iobs)*y%field(iobs,3)
	                dpmsl(1:nxtc,1:nytc,1) = dpmsl(1:nxtc,1:nytc,1) + weight(imember)*self%X(1:nxtc,1:nytc,1,imember,iobs)
	             end do
	             call allreduce3D('e', dpmsl(1:nxtc,1:nytc,:), nxtc, nytc, 1)
	             !
	             do ivec = 1, self%nrank(iobs)
	                coordinate(ivec) = sum(self%eigvec(1:nxtc,1:nytc,1,ivec,iobs)*dpmsl(1:nxtc,1:nytc,1))
	                coordinate(ivec) = coordinate(ivec)*self%eigval(ivec,iobs)
	             end do
	             dpmsl(:,:,1) = 0.d0
	             do ivec = 1, self%nrank(iobs)
	                dpmsl(imin:imax,jmin:jmax,1) = dpmsl(imin:imax,jmin:jmax,1) + coordinate(ivec)*self%eigvec(1:nxtc,1:nytc,1,ivec,iobs)
	             end do
	             dpmslall(imin:imax,jmin:jmax,1) = dpmslall(imin:imax,jmin:jmax,1) + dpmsl(imin:imax,jmin:jmax,1)
	             jobs = jobs + 1
	             if (jobs == mobs) exit
	          end do
	          if (jobs == mobs) exit
	       end do
	       if (jobs == mobs) exit
	    end do
	 end if
	 !
	 ! add pmsl
         do itcproc = 1, ntcproc
            if (myid == tcproc(itcproc)) data(:,:,1) = dpmslall(:,:,1)
            call scatter(info, 'xy', data, tcproc(itcproc), nx0, ny0, 1)
            call add_field(x, 'pmsl', it, data(1:nx,1:ny,1:1), 1)
         end do
      end do
      if (nobs > 0) deallocate(coordinate, weight, dpmsl, dpmslall)
      deallocate(data)
      !
      return
   end subroutine interpolate_DtcT
   !
   !
   !
   subroutine apply_H_NodeHFieldTC1(self, info, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      call interpolate_tc1(self, info, x, obsspace, valid, y)
      !
      return
   end subroutine apply_H_NodeHFieldTC1
   !
   !
   !
   subroutine apply_H_NodeHFieldTC2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: nx, ny, nxyt
      logical, dimension(nx,ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      !call interpolate_tc2(self, info, processed, k2ijt, x, obsspace, valid, y, nx, ny, nxyt)
      !
      return
   end subroutine apply_H_NodeHFieldTC2
   !
   !
   !
   subroutine apply_H_NodeHFieldTC3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      integer, intent(in) :: ip, jp, nt
      integer, dimension(2,2,nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      !call interpolate_tc3(self, info, ip, jp, ijt2k, x, obsspace, valid, y, nt)
      !
      return
   end subroutine apply_H_NodeHFieldTC3
   !
   !
   !
   subroutine initialize_DH_NodeHFieldTC(self, info, x, obsspace, valid, y, ne)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      integer, intent(in) :: ne
      type(NodeInfo), intent(in) :: info
      type(NodeControl), dimension(ne), intent(in) :: x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsControl), dimension(ne), intent(in) :: y
      !
      call initialize_Dtc(self, info, x, obsspace, valid, y, ne)
      !
      return
   end subroutine initialize_DH_NodeHFieldTC
   !
   !
   !
   subroutine apply_DH_NodeHFieldTC(self, info, xbck, x, obsspace, valid, y)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck, x
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(inout) :: valid
      type(NodeObsField), intent(inout) :: y
      !
      call interpolate_Dtc(self, info, x, obsspace, valid, y)
      !
      return
   end subroutine apply_DH_NodeHFieldTC
   !
   !
   !
   subroutine apply_DHT_NodeHFieldTC(self, info, xbck, obsspace, valid, y, x)
      implicit none
      type(NodeHFieldTC), intent(inout) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeControl), intent(in) :: xbck
      type(NodeObsSpaceFieldTC), intent(in) :: obsspace
      type(NodeObsValidField), intent(in) :: valid
      type(NodeObsField), intent(in) :: y
      type(NodeControl), intent(inout) :: x
      !
      call interpolate_DtcT(self, info, obsspace, valid, y, x)
      !
      return
   end subroutine apply_DHT_NodeHFieldTC
   !
   !
   !
end module NodeHFieldTC_class
