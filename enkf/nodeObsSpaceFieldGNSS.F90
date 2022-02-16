module NodeObsSpaceFieldGNSS_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl, iflg_incremental, gnssobs_format
   use NodeInfo_class
   use NodeMPI
   implicit none
   !
   ! --- CDA parameters
   integer, parameter :: LX1=14, LX2=25, LX3=50, LX4=65200 ! array size of each part of CDA
   integer, parameter :: NXHT =    250 ! array size of level block in one record
   integer, parameter :: NXIH =     20 ! 
   integer, parameter :: NXEL =     60 ! array size of elem block in one level
   integer, parameter :: NXIE =     40 ! 
   integer, parameter :: NXQC =     50
   integer, parameter :: MISS = -2**15 ! missing value
   !
   integer, parameter :: LBSRF = 6        ! surface
   integer, parameter :: KHCD_IVNM = 100  ! level name
   integer, parameter :: KHCD_IVUT = 200  ! unit of height
   integer, parameter :: KHCD_RPHT = 500  ! reported height
   !
   integer, parameter :: KECD_OBQC = 200  ! QCed observation value
   integer, parameter :: KECD_OBER = 400  ! observation error
   integer, parameter :: KECD_DVAL = 500  ! D-value
   !
   integer, parameter :: MREFR = 61*256+ 3 ! refractivity (gnss ro)
   integer, parameter :: MBANG = 61*256+ 4 ! bending angle (gnss ro)
   integer, parameter :: MIMPA =  7*256+40 ! impact parameter (need for gnss ro)
   !
   !integer, parameter :: IVITP_HEIGHT = 0  ! with height(m)
   !
   ! Obs are stored as a linked list at each horizontal grid point. To access
   ! obs fast, all link lists are concatenated into a one-dimensional array. The
   ! iobs and mobs are the indices fo accessing: iobs is the index of each link
   ! list in this array (=0 if no link list) and mobs is the number of obs in each link list.
   ! We do not put all information like lon, lat, ... into a structure because these fields
   ! depend on observation types.
   type NodeObsSpaceFieldGNSS
      character(len=10) :: name, filename
      ! For GNSS obs, nsubobs = 1
      integer :: nobs, nsubobs, nx, ny, nt
      integer, dimension(:), allocatable :: cdakind
      real(r_size), dimension(:), allocatable :: lon, lat, height
      real(r_size), dimension(:), allocatable :: rcurve, azimuth, undulation, impact
      real(r_size), dimension(:,:), allocatable :: obs, error
      integer, dimension(:,:,:), allocatable :: iobs, mobs
   end type NodeObsSpaceFieldGNSS
   !
   interface new
      module procedure new_NodeObsSpaceFieldGNSS
   end interface
   interface destroy
      module procedure destroy_NodeObsSpaceFieldGNSS
   end interface
   interface display
      module procedure display_NodeObsSpaceFieldGNSS
   end interface
   interface get_name
      module procedure get_name_NodeObsSpaceFieldGNSS
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsSpaceFieldGNSS
   end interface
   interface get_nsubobs
      module procedure get_nsubobs_NodeObsSpaceFieldGNSS
   end interface
   interface get_mobs
      module procedure get_mobs_NodeObsSpaceFieldGNSS1
      module procedure get_mobs_NodeObsSpaceFieldGNSS2
   end interface
   interface read_obs
      module procedure read_obs_NodeObsSpaceFieldGNSS
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsSpaceFieldGNSS
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsSpaceFieldGNSS
   end interface
   !
contains
   !
   subroutine new_NodeObsSpaceFieldGNSS(self, name, filename, nx, ny, nt)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(inout) :: self
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
   end subroutine new_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine destroy_NodeObsSpaceFieldGNSS(self)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(inout) :: self
      integer :: iobs
      !
      if (allocated(self%iobs)) deallocate(self%iobs)
      if (allocated(self%mobs)) deallocate(self%mobs)
      if (allocated(self%cdakind)) deallocate(self%cdakind)
      if (allocated(self%lon)) deallocate(self%lon)
      if (allocated(self%lat)) deallocate(self%lat)
      if (allocated(self%height)) deallocate(self%height)
      if (allocated(self%rcurve)) deallocate(self%rcurve)
      if (allocated(self%azimuth)) deallocate(self%azimuth)
      if (allocated(self%undulation)) deallocate(self%undulation)
      if (allocated(self%impact)) deallocate(self%impact)
      if (allocated(self%obs)) deallocate(self%obs)
      if (allocated(self%error)) deallocate(self%error)
      !
      return
   end subroutine destroy_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine display_NodeObsSpaceFieldGNSS(self)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      !
      print*, 'Observation name: ', self%name
      print*, 'File name: ', self%filename
      print*, 'Number: ', self%nobs, self%nsubobs
      !
      return
   end subroutine display_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine get_name_NodeObsSpaceFieldGNSS(self, name)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine get_nobs_NodeObsSpaceFieldGNSS(self, nobs)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine get_nsubobs_NodeObsSpaceFieldGNSS(self, nsubobs)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      integer, intent(out) :: nsubobs
      !
      nsubobs = self%nsubobs
      !
      return
   end subroutine get_nsubobs_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldGNSS1(self, i, j, mobs)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      integer, intent(in) :: i, j
      integer, intent(out) :: mobs
      !
      mobs = sum(self%mobs(i,j,:))*self%nsubobs
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldGNSS1
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldGNSS2(self, i, j, it, mobs)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      integer, intent(in) :: i, j, it
      integer, intent(out) :: mobs
      !
      mobs = self%mobs(i,j,it)*self%nsubobs
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldGNSS2
   !
   !
   !
   subroutine cdaread_nobs(self, nobs)
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos, proj
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      integer, intent(out) :: nobs
      !
      integer(4) :: NHT, IVNM, IEND
      integer(2) :: IC1(1:LX1)
      integer(2) :: IC2(1:LX2)
      integer(2) :: IC3(1:LX3)
      integer(2) :: IC4(1:LX4)
      !
      real   (4) :: HGT (1:NXHT)
      integer(4) :: NEL (1:NXHT)
      integer(4) :: NHCD(1:NXHT)
      integer(4) :: KHCD(1:NXIH,1:NXHT)
      real   (4) :: AHCD(1:NXIH,1:NXHT)
      integer(4) :: MEL (1:NXEL,1:NXHT)
      integer(4) :: NECD(1:NXEL,1:NXHT)
      integer(4) :: JGQC(1:NXEL,1:NXHT)
      integer(4) :: IUSE(1:NXEL,1:NXHT)
      integer(4) :: NOQC(1:NXEL,1:NXHT)
      integer(4) :: KECD(1:NXIE,1:NXEL,1:NXHT)
      real   (4) :: AECD(1:NXIE,1:NXEL,1:NXHT)
      integer(4) :: JEQC(1:NXQC,1:NXEL,1:NXHT)
      integer(4) :: KEQC(1:NXQC,1:NXEL,1:NXHT)
      !
      logical :: limpact
      character(2) :: infile = '00'
      integer :: nx, ny, nt, cdacode, cdakind, ztype
      integer :: it, nh, kh, ne, ke, i, j, ierror
      real(kind=r_sngl) :: lat, lon, height  ! lat, lon, rhgt must be single
      type(PosGrid) :: pg
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      if (trim(self%name) == 're') then
	 cdacode = MREFR
      else if (trim(self%name) == 'ba') then
	 cdacode = MBANG
      end if
      !
      nobs = 0
      LOOP_SLOT: do it = 1, nt
         write(infile(1:2),'(I2.2)') it
	 open(90, file=trim(self%filename)//infile, form='unformatted', access='sequential', action='read', iostat = ierror)
         if (ierror > 0) cycle
	 !
	 LOOP_RECORD: do
	    call CDARD(IC1, IC2, IC3, IC4, IEND, 90)
	    if (IEND == 1) exit
	    !--- in case of DATE record
	    if (IC1(7) == 0) cycle
	    !
	    ! convert lat,lon to grid coordinates
	    cdakind = IC1(8)
	    ztype = -1
	    lon = real(IC1(10),r_sngl)/100.e0
	    lat = real(IC1(9),r_sngl)/100.e0
	    call convert_LatLon_to_GridPos(pg, lat, lon, 1)
	    i = pg%i1; j = pg%j1
	    if (i < 2 .or. i > nx-1 .or. j < 2 .or. j > ny-1) cycle
	    !
	    ! unpack CDA part 4
	    call CDAUNPK(IC1, IC2, IC3, IC4, NHT, HGT, NEL, MEL, &
		       & NHCD, KHCD, AHCD, NECD, KECD, AECD, &
		       & JGQC, IUSE, NOQC, JEQC, KEQC, &
		       & MISS, NXHT, NXEL, NXIH, NXIE, NXQC)
            LOOP_LEVEL: do nh = 1, NHT
	       do kh = 1, NHCD(nh)
		  if (KHCD(kh,nh) == KHCD_RPHT) then           
		     height = AHCD(kh,nh)       ! reported height
		     ztype = 0 !IVITP_HEIGHT            ! interpolate with m
		     exit
		  end if
	       end do
	       ! height data is necessary for upper data
	       if (ztype < 0) cycle
	       if (int(height) == MISS) cycle
               !
	       ! find impact parameter
	       if (trim(self%name) == 'ba') then
		     LOOP_ELEM_PRE: do ne=1,NEL(nh)
			! variable selection
			limpact = .false.
			if( MEL(ne,nh) == MIMPA ) then  ! impact parameter
			   limpact = .true.
			   exit
			end if
		     end do LOOP_ELEM_PRE
		     if(.not. limpact) cycle
	       end if
	       !
	       ! set obsereved value, observation error
               do ne = 1, NEL(nh)
                  if (IUSE(ne,nh) /= 1) cycle             ! use only IUSE flag == 1
                  if (mel(ne,nh) /= cdacode) cycle
		  nobs = nobs + 1
	       end do
            end do LOOP_LEVEL
         end do LOOP_RECORD
         close(90)
      end do LOOP_SLOT
      print*, self%name, 'nobs = ', nobs
      !
      return
   end subroutine cdaread_nobs
   !
   !
   !
   subroutine cdaread_obs(self, myid)
      ! self is a global ObsSpace
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos, proj
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(inout) :: self
      integer, intent(in) :: myid
      !
      integer(2) :: IC1(1:LX1)
      integer(2) :: IC2(1:LX2)
      integer(2) :: IC3(1:LX3)
      integer(2) :: IC4(1:LX4)
      !
      integer(4) :: NHT, IVNM, IEND
      integer(4) :: ifound_obqc, ifound_ober, ifound_dval
      real(4)    :: OBQC, OBER, DVAL
      real   (4) :: HGT (1:NXHT)
      integer(4) :: NEL (1:NXHT)
      integer(4) :: NHCD(1:NXHT)
      integer(4) :: KHCD(1:NXIH,1:NXHT)
      real   (4) :: AHCD(1:NXIH,1:NXHT)
      integer(4) :: MEL (1:NXEL,1:NXHT)
      integer(4) :: NECD(1:NXEL,1:NXHT)
      integer(4) :: JGQC(1:NXEL,1:NXHT)
      integer(4) :: IUSE(1:NXEL,1:NXHT)
      integer(4) :: NOQC(1:NXEL,1:NXHT)
      integer(4) :: KECD(1:NXIE,1:NXEL,1:NXHT)
      real   (4) :: AECD(1:NXIE,1:NXEL,1:NXHT)
      integer(4) :: JEQC(1:NXQC,1:NXEL,1:NXHT)
      integer(4) :: KEQC(1:NXQC,1:NXEL,1:NXHT)
      !
      character(2) :: infile = '00'
      logical :: limpact
      integer :: cdacode, cdakind, ztype
      integer :: nx, ny, nt, nobs, nsubobs, ierror
      integer :: it, nh, kh, ne, ke, i, j, iobs, iobs1, iobs2, m, itmp, jtmp
      real(kind=r_sngl) :: lat, lon
      real(r_size) :: obs, error, lonobs, latobs, height, rcurve, azimuth, undulation, impact
      type(PosGrid) :: pg
      integer, dimension(:), allocatable :: i0, j0
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      !
      if (myid == 0) call cdaread_nobs(self, nobs)
      call int_broadcast0D('all', nobs, 0)
      self%nobs = nobs
      nsubobs = 1
      self%nsubobs = nsubobs
      if (nobs == 0) return
      allocate(self%cdakind(nobs))
      allocate(self%lon(nobs), self%lat(nobs), self%height(nobs))
      allocate(self%obs(nobs,nsubobs), self%error(nobs,nsubobs))
      if (trim(self%name) == 'ba') allocate(self%rcurve(nobs), self%azimuth(nobs), self%undulation(nobs), self%impact(nobs))
      !
      if (myid == 0) then
	 allocate(i0(nobs), j0(nobs))
	 if (trim(self%name) == 're') then
	    cdacode = MREFR
         else if (trim(self%name) == 'ba') then
	    cdacode = MBANG
         end if
	 !
	 iobs = 0
	 LOOP_SLOT: do it = 1, nt
	    write(infile(1:2),'(I2.2)') it
	    open(90, file=trim(self%filename)//infile, form='unformatted', access='sequential', action='read', iostat = ierror)
	    if (ierror > 0) cycle
	    !
	    LOOP_RECORD: do
	       call CDARD(IC1, IC2, IC3, IC4, IEND, 90)
	       if (IEND == 1) exit
	       !--- in case of DATE record
	       if (IC1(7) == 0) cycle
	       !
	       ! convert lat,lon to grid coordinates
	       cdakind = IC1(8)
	       ztype = -1
	       lon = real(IC1(10),r_sngl)/100.e0
	       lat = real(IC1(9),r_sngl)/100.e0
	       call convert_LatLon_to_GridPos(pg, lat, lon, 1)
	       i = pg%i1; j = pg%j1
	       if (i < 2 .or. i > nx-1 .or. j < 2 .or. j > ny-1) cycle
	       !
	       call movec(IEND, 1, ic2(13), 1, 4)
#ifdef SWP_ON
               call swp2swp4_l(IEND)
#endif
               latobs = real(IEND,r_size)*0.000001d0
               call movec(IEND, 1, ic2(15), 1, 4)
#ifdef SWP_ON
               call swp2swp4_l(IEND)
#endif
               lonobs = real(IEND,r_size)*0.000001d0
               if (trim(self%name) == 'ba') then
                  call movec(IEND, 1, ic3(9), 1, 4)
#ifdef SWP_ON
                  call swp2swp4_l(IEND)
#endif
                  rcurve = real(IEND,r_size)*0.1d0
		  azimuth = real(ic3(11),r_size)*0.01d0
                  undulation = real(ic3(12),r_size)*0.01d0
               end if
	       !
	       ! unpack CDA part 4
	       call CDAUNPK(IC1, IC2, IC3, IC4, NHT, HGT, NEL, MEL, &
			  & NHCD, KHCD, AHCD, NECD, KECD, AECD, &
			  & JGQC, IUSE, NOQC, JEQC, KEQC, &
			  & MISS, NXHT, NXEL, NXIH, NXIE, NXQC)
	       LOOP_LEVEL: do nh = 1, NHT
		  do kh = 1, NHCD(nh)
		     if (KHCD(kh,nh) == KHCD_RPHT) then           
			height = AHCD(kh,nh)       ! reported height
			ztype = 0 ! IVITP_HEIGHT            ! interpolate with m
		        exit
		     end if
		  end do
		  ! height data is necessary for upper data
		  if (ztype < 0) cycle
		  if (int(height) == MISS) cycle
		  !
		  ! find impact parameter
		  if (trim(self%name) == 'ba') then
		     LOOP_ELEM_PRE: do ne=1,NEL(nh)
			! variable selection
			limpact = .false.
			if( MEL(ne,nh) == MIMPA ) then  ! impact parameter
			   LOOP_SELECT_PRE: do ke=1,NECD(ne,nh)
			      ! QCed observataion value
			      if( KECD(ke,ne,nh) == KECD_OBQC ) then
				 impact = AECD(ke,ne,nh)
				 exit
			      end if
			   end do LOOP_SELECT_PRE
			   limpact = .true.
			   exit
			end if
		     end do LOOP_ELEM_PRE
		     if(.not. limpact) cycle
		  end if
				  !
				  ! set obsereved value, observation error
				  do ne = 1, NEL(nh)
					 if (IUSE(ne,nh) /= 1) cycle             ! use only IUSE flag == 1
					 if (mel(ne,nh) /= cdacode) cycle
					 !
					 ifound_obqc = 0; ifound_ober = 0; ifound_dval = 0
					 do ke = 1, NECD(ne,nh)
						! QCed observataion value
						if (KECD(ke,ne,nh) == KECD_OBQC .and. ifound_obqc == 0) then
						   obqc = AECD(ke,ne,nh)
						   ifound_obqc = 1
						end if
						! observation error
						if (KECD(ke,ne,nh) == KECD_OBER .and. ifound_ober == 0) then
						   ober = AECD(ke,ne,nh)
						   ifound_ober = 1
						end if
						! D-value
						if (KECD(ke,ne,nh) == KECD_DVAL .and. ifound_dval == 0) then
						   dval = AECD(ke,ne,nh)
						   ifound_dval = 1
						end if
						if (ifound_obqc*ifound_ober*ifound_dval == 1) exit
					 end do
					 !
					 if (iflg_incremental == 1) then ! INCREMENTAL APPROACH (2006/8/17)
						obs = dval
					 else
						obs = obqc
					 end if
					 error = ober
					 ! convert unit
					 ! microrad to rad
					 if (mel(ne,nh) == MBANG) then  
						obs = obs*0.000001d0
						error = error*0.000001d0
					 end if
					 !
					 iobs = iobs + 1
					 self%mobs(i,j,it) = self%mobs(i,j,it) + 1
					 self%cdakind(iobs) = cdakind
					 i0(iobs) = i; j0(iobs) = j
					 self%lon(iobs) = lonobs; self%lat(iobs) = latobs
					 self%height(iobs) = height
					 self%obs(iobs,1) = obs; self%error(iobs,1) = error
					 if (trim(self%name) == 'ba') then
					 	self%rcurve(iobs) = rcurve; self%azimuth(iobs) = azimuth
					    self%undulation(iobs) = undulation; self%impact(iobs) = impact
					 end if
					 !print*, cdakind, mel(ne,nh), lon, lat, height, obs, error
                  end do
               end do LOOP_LEVEL
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
						   cdakind = self%cdakind(iobs1+j); self%cdakind(iobs1+j) = self%cdakind(iobs1+j-m); self%cdakind(iobs1+j-m) = cdakind
						   lonobs = self%lon(iobs1+j); self%lon(iobs1+j) = self%lon(iobs1+j-m); self%lon(iobs1+j-m) = lonobs
						   latobs = self%lat(iobs1+j); self%lat(iobs1+j) = self%lat(iobs1+j-m); self%lat(iobs1+j-m) = latobs
						   height = self%height(iobs1+j); self%height(iobs1+j) = self%height(iobs1+j-m); self%height(iobs1+j-m) = height
						   obs = self%obs(iobs1+j,1); self%obs(iobs1+j,1) = self%obs(iobs1+j-m,1); self%obs(iobs1+j-m,1) = obs
						   error = self%error(iobs1+j,1); self%error(iobs1+j,1) = self%error(iobs1+j-m,1); self%error(iobs1+j-m,1) = error
						   if (trim(self%name) == 'ba') then
							  rcurve = self%rcurve(iobs1+j); self%rcurve(iobs1+j) = self%rcurve(iobs1+j-m); self%rcurve(iobs1+j-m) = rcurve
							  azimuth = self%azimuth(iobs1+j); self%azimuth(iobs1+j) = self%azimuth(iobs1+j-m); self%azimuth(iobs1+j-m) = azimuth
							  undulation = self%undulation(iobs1+j); self%undulation(iobs1+j) = self%undulation(iobs1+j-m); self%undulation(iobs1+j-m) = undulation
							  impact = self%impact(iobs1+j); self%impact(iobs1+j) = self%impact(iobs1+j-m); self%impact(iobs1+j-m) = impact
						   end if
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
      call int_broadcast1D('all', self%cdakind, 0, nobs)
      call broadcast1D('all', self%lon, 0, nobs)
      call broadcast1D('all', self%lat, 0, nobs)
      call broadcast1D('all', self%height, 0, nobs)
      call broadcast2D('all', self%obs, 0, nobs, nsubobs)
      call broadcast2D('all', self%error, 0, nobs, nsubobs)
      if (trim(self%name) == 'ba') then
      	 call broadcast1D('all', self%rcurve, 0, nobs)
	 call broadcast1D('all', self%azimuth, 0, nobs)
	 call broadcast1D('all', self%undulation, 0, nobs)
	 call broadcast1D('all', self%impact, 0, nobs)
      end if
      !
      return
   end subroutine cdaread_obs
   !
   !
   !
   subroutine read_obs_NodeObsSpaceFieldGNSS(self, myid)
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(inout) :: self
      integer, intent(in) :: myid
      !
      if (trim(gnssobs_format) == 'cda') then
	     call cdaread_obs(self, myid)
      else if (trim(gnssobs_format) == 'txt') then
	     !call txtread_obs(self, myid)
      end if
      !
      return
   end subroutine read_obs_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine scatter_obs_NodeObsSpaceFieldGNSS(self, info, local_object)
      ! scatter a big obsspace from Node 0 to all nodes.
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldGNSS), intent(inout) :: local_object
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
		 allocate(local_object%cdakind(nobs))
		 allocate(local_object%lon(nobs), local_object%lat(nobs), local_object%height(nobs))
		 allocate(local_object%obs(nobs,nsubobs), local_object%error(nobs,nsubobs))
		 if (trim(self%name) == 'ba') then
		 	allocate(local_object%rcurve(nobs), local_object%azimuth(nobs))
		 	allocate(local_object%undulation(nobs), local_object%impact(nobs))
		 end if
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
					 local_object%cdakind(iobs1:iobs2) = self%cdakind(iobs3:iobs4)
					 local_object%lon(iobs1:iobs2) = self%lon(iobs3:iobs4)
					 local_object%lat(iobs1:iobs2) = self%lat(iobs3:iobs4)
					 local_object%height(iobs1:iobs2) = self%height(iobs3:iobs4)
					 local_object%obs(iobs1:iobs2,:) = self%obs(iobs3:iobs4,:)
					 local_object%error(iobs1:iobs2,:) = self%error(iobs3:iobs4,:)
					 if (trim(self%name) == 'ba') then
					    local_object%rcurve(iobs1:iobs2) = self%rcurve(iobs3:iobs4)
					    local_object%azimuth(iobs1:iobs2) = self%azimuth(iobs3:iobs4)
					    local_object%undulation(iobs1:iobs2) = self%undulation(iobs3:iobs4)
					    local_object%impact(iobs1:iobs2) = self%impact(iobs3:iobs4)
				     end if
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
   end subroutine scatter_obs_NodeObsSpaceFieldGNSS
   !
   !
   !
   subroutine gather_obs_NodeObsSpaceFieldGNSS(self, info, global_object)
      ! this subroutine is only for completion and we realy do not use this
      implicit none
      type(NodeObsSpaceFieldGNSS), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldGNSS), intent(inout) :: global_object
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
   end subroutine gather_obs_NodeObsSpaceFieldGNSS
   !
   !
   !
end module NodeObsSpaceFieldGNSS_class
