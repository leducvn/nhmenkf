module NodeObsSpaceFieldCNV_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl, pi, iflg_incremental, iflg_surface_obs, iflg_skip_psobs, cnvobs_format
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
   integer, parameter :: MZHGT = 10*256+ 7 ! geopotential height
   integer, parameter :: MWUEW = 11*256+ 3 ! wind U (east-west)
   integer, parameter :: MWVSN = 11*256+ 4 ! wind V (south-north)
   integer, parameter :: MZSMD = 57*256+ 2 ! height over model surface
   integer, parameter :: MWULM = 57*256+ 3 ! wind U (on Lambert)
   integer, parameter :: MWVLM = 57*256+ 4 ! wind V (on Lambert)
   integer, parameter :: MPSMD = 57*256+ 1 ! model surface pressure
   integer, parameter :: MTEMP = 12*256+ 1 ! temperature
   integer, parameter :: MRHMD = 13*256+ 3 ! rerative humidity
   integer, parameter :: MPCWT = 13*256+16 ! precipitable water
   integer, parameter :: MWRAD = 21*256+14 ! dopper radar: radial velocity
   !
   integer, parameter :: IVITP_HEIGHT = 0  ! with height(m)
   integer, parameter :: IVITP_P      = 1  ! with pressure (not implemented)
   integer, parameter :: IVITP_LOGP   = 2  ! with log of pressure
   integer, parameter :: IVITP_SURF   = 3  ! surface 2D interpolation
   integer, parameter :: IVITP_GRID   = 4  ! with grid coordinate
   !
   ! Obs are stored as a linked list at each horizontal grid point. To access
   ! obs fast, all link lists are concatenated into a one-dimensional array. The
   ! iobs and mobs are the indices fo accessing: iobs is the index of each link
   ! list in this array (=0 if no link list) and mobs is the number of obs in each link list.
   ! We do not put all information like lon, lat, ... into a structure because these fields
   ! depend on observation types.
   type NodeObsSpaceFieldCNV
      character(len=10) :: name, filename
      ! For conventional obs, nsubobs = 1
      integer :: nobs, nsubobs, nx, ny, nt
      integer, dimension(:), allocatable :: cdakind, ztype
      real(r_sngl), dimension(:), allocatable :: lon, lat, height
      real(r_size), dimension(:), allocatable :: range, azimuth
      real(r_size), dimension(:,:), allocatable :: obs, error
      integer, dimension(:,:,:), allocatable :: iobs, mobs
   end type NodeObsSpaceFieldCNV
   !
   interface new
      module procedure new_NodeObsSpaceFieldCNV
   end interface
   interface destroy
      module procedure destroy_NodeObsSpaceFieldCNV
   end interface
   interface display
      module procedure display_NodeObsSpaceFieldCNV
   end interface
   interface get_name
      module procedure get_name_NodeObsSpaceFieldCNV
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsSpaceFieldCNV
   end interface
   interface get_nsubobs
      module procedure get_nsubobs_NodeObsSpaceFieldCNV
   end interface
   interface get_mobs
      module procedure get_mobs_NodeObsSpaceFieldCNV1
      module procedure get_mobs_NodeObsSpaceFieldCNV2
   end interface
   interface read_obs
      module procedure read_obs_NodeObsSpaceFieldCNV
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsSpaceFieldCNV
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsSpaceFieldCNV
   end interface
   !
contains
   !
   subroutine new_NodeObsSpaceFieldCNV(self, name, filename, nx, ny, nt)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(inout) :: self
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
   end subroutine new_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine destroy_NodeObsSpaceFieldCNV(self)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(inout) :: self
      integer :: iobs
      !
      if (allocated(self%iobs)) deallocate(self%iobs)
      if (allocated(self%mobs)) deallocate(self%mobs)
      if (allocated(self%cdakind)) deallocate(self%cdakind)
      if (allocated(self%ztype)) deallocate(self%ztype)
      if (allocated(self%lon)) deallocate(self%lon)
      if (allocated(self%lat)) deallocate(self%lat)
      if (allocated(self%height)) deallocate(self%height)
      if (allocated(self%range)) deallocate(self%range)
      if (allocated(self%azimuth)) deallocate(self%azimuth)
      if (allocated(self%obs)) deallocate(self%obs)
      if (allocated(self%error)) deallocate(self%error)
      !
      return
   end subroutine destroy_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine display_NodeObsSpaceFieldCNV(self)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      !
      print*, 'Observation name: ', self%name
      print*, 'File name: ', self%filename
      print*, 'Number: ', self%nobs, self%nsubobs
      !
      return
   end subroutine display_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine get_name_NodeObsSpaceFieldCNV(self, name)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine get_nobs_NodeObsSpaceFieldCNV(self, nobs)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine get_nsubobs_NodeObsSpaceFieldCNV(self, nsubobs)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      integer, intent(out) :: nsubobs
      !
      nsubobs = self%nsubobs
      !
      return
   end subroutine get_nsubobs_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldCNV1(self, i, j, mobs)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      integer, intent(in) :: i, j
      integer, intent(out) :: mobs
      !
      mobs = sum(self%mobs(i,j,:))*self%nsubobs
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldCNV1
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldCNV2(self, i, j, it, mobs)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      integer, intent(in) :: i, j, it
      integer, intent(out) :: mobs
      !
      mobs = self%mobs(i,j,it)*self%nsubobs
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldCNV2
   !
   !
   !
   subroutine cdaread_nobs(self, nobs)
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos, proj
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      integer, intent(out) :: nobs
      !
      integer :: MWUCP, MWVCP
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
      character(2) :: infile = '00'
      integer :: nx, ny, nt, cdacode, cdakind, ztype
      integer :: it, nh, kh, ne, ke, i, j, ierror
      real(kind=r_sngl) :: lat, lon, height  ! lat, lon, rhgt must be single
      type(PosGrid) :: pg
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      if (proj == 'LMN ' .or. proj == 'LMS ') then
         MWUCP = MWULM
         MWVCP = MWVLM
      else if (proj == 'MER ') then
         MWUCP = MWUEW
         MWVCP = MWVSN
      else
	 MWUCP = MWUEW
	 MWVCP = MWVSN
      end if
      if (trim(self%name) == 'u') then
	     cdacode = MWUCP
      else if (trim(self%name) == 'v') then
	     cdacode = MWVCP
      else if (trim(self%name) == 't') then
	     cdacode = MTEMP
      else if (trim(self%name) == 'p') then
	     cdacode = MPSMD
      else if (trim(self%name) == 'rh') then
	     cdacode = MRHMD
      else if (trim(self%name) == 'pwv') then
	     cdacode = MPCWT
      else if (trim(self%name) == 'rvl') then
	     cdacode = MWRAD
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
			if (cdakind == 15100 .or. cdakind == 15205) cycle  ! skip ATOVS and SATEM report
			if (cdakind == 17100 .or. cdakind == 17200 .or. &
		   &    cdakind == 17350 .or. cdakind == 17400 ) cycle ! skip Satellite TPW retrievals
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
			   IVNM = 0
			   do kh = 1, NHCD(nh)
				  if (KHCD(kh,nh) == KHCD_IVNM) then   ! vertical indicator
					 IVNM = AHCD(kh,nh)
					 exit
				  end if
			   end do
			   !
			   ! surface
			   if (btest(IVNM,LBSRF)) then
				  ztype = IVITP_SURF                ! interpolate with SURF-2D
				  if (cdakind/1000 == 2) then
				     height = 0.0e0                   ! height = 0m
				  elseif (cdakind/1000 == 1 .or. cdakind/1000 == 3 .or.  &
					    & cdakind/1000 == 14 .or. cdakind/1000 == 17.or. &
					    & cdakind == 5900) then
				     height = real(IC3(9),r_sngl)    ! model surface altitue
				  elseif (cdakind/1000 == 16 ) then
					 height = 10.0
					 ztype = IVITP_HEIGHT
				  else
				     cycle
				  end if
			   ! upper
			   else
				  if ((cdakind/1000) == 6 .or. (cdakind/1000) == 7 .or. (cdakind/1000) == 19) then
				     do kh = 1, NHCD(nh)
					    if (KHCD(kh,nh) == KHCD_RPHT) then           
					       height = AHCD(kh,nh)       ! reported height
					       exit
					    end if
				     end do
				     ztype = IVITP_HEIGHT            ! interpolate with m
				  elseif (cdakind/1000 == 16) then
					 height = 10.0
					 ztype = IVITP_HEIGHT
				  else
				     do kh = 1,NHCD(nh)
					    if (KHCD(kh,nh) == KHCD_IVUT) then   ! unit of height
						   if (AHCD(kh,nh) == 1) then    ! unit = hPa
							  height = HGT(nh)
							  ztype = IVITP_LOGP        ! interpolate with logP
							  exit
						   end if
					    end if
				     end do
				  end if
				  ! height data is necessary for upper data
				  if (ztype < 0) cycle
				  if (int(height) == MISS) cycle
               end if
               !
	           ! set obsereved value, observation error
               do ne = 1, NEL(nh)
                  if (IUSE(ne,nh) /= 1) cycle             ! use only IUSE flag == 1
                  if (mel(ne,nh) /= cdacode) cycle
		          if (mel(ne,nh) == MPSMD .and. iflg_skip_psobs == 1) cycle
		          ! select surface observation
	              select case (iflg_surface_obs)
                  case (1)
		             if (ztype == IVITP_SURF) then
						if (cdakind/100 == 164 .and. (mel(ne,nh) == MWUCP .or. mel(ne,nh) == MWVCP)) then
						else if((cdakind/1000 == 17 .or. cdakind/1000 == 14) .and. mel(ne,nh) == MPCWT) then
						else if((cdakind/1000 == 1 .or. cdakind/1000 == 2) .and. mel(ne,nh) == MPSMD) then ! model surface pressure
						else if (cdakind == 1400) then
						else if (cdakind == 5900 ) then
						else
						   cycle
						end if
		             end if
	              case default
                  end select
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
      type(NodeObsSpaceFieldCNV), intent(inout) :: self
      integer, intent(in) :: myid
      !
      integer :: MWUCP, MWVCP
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
      integer :: cdacode, cdakind, ztype
      integer :: nx, ny, nt, nobs, nsubobs, ierror
      integer :: it, nh, kh, ne, ke, i, j, iobs, iobs1, iobs2, m, itmp, jtmp
      real(kind=r_sngl) :: lat, lon, height
      real(r_size) :: obs, error, range, azimuth
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
      allocate(self%cdakind(nobs), self%ztype(nobs))
      allocate(self%lon(nobs), self%lat(nobs), self%height(nobs))
      allocate(self%obs(nobs,nsubobs), self%error(nobs,nsubobs))
      if (trim(self%name) == 'rvl') allocate(self%range(nobs), self%azimuth(nobs))
      !
      if (myid == 0) then
		 allocate(i0(nobs), j0(nobs))
		 if (proj == 'LMN ' .or. proj == 'LMS ') then
			MWUCP = MWULM
			MWVCP = MWVLM
		 else if (proj == 'MER ') then
			MWUCP = MWUEW
			MWVCP = MWVSN
		 else
			MWUCP = MWUEW
			MWVCP = MWVSN
	     end if
		 if (trim(self%name) == 'u') then
		    cdacode = MWUCP
		 else if (trim(self%name) == 'v') then
		    cdacode = MWVCP
		 else if (trim(self%name) == 't') then
		    cdacode = MTEMP
		 else if (trim(self%name) == 'p') then
		    cdacode = MPSMD
		 else if (trim(self%name) == 'rh') then
		    cdacode = MRHMD
		 else if (trim(self%name) == 'pwv') then
		    cdacode = MPCWT
		 else if (trim(self%name) == 'rvl') then
	        cdacode = MWRAD
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
			   if (cdakind == 15100 .or. cdakind == 15205) cycle  ! skip ATOVS and SATEM report
               if (cdakind == 17100 .or. cdakind == 17200 .or. &
              &    cdakind == 17350 .or. cdakind == 17400 ) cycle ! skip Satellite TPW retrievals
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
				  IVNM = 0
				  do kh = 1, NHCD(nh)
					 if (KHCD(kh,nh) == KHCD_IVNM) then   ! vertical indicator
						IVNM = AHCD(kh,nh)
						exit
					 end if
				  end do
				  !
				  ! surface
				  if (btest(IVNM,LBSRF)) then
					 ztype = IVITP_SURF                ! interpolate with SURF-2D
					 if (cdakind/1000 == 2) then
						height = 0.0e0                   ! height = 0m
					 elseif (cdakind/1000 == 1 .or. cdakind/1000 == 3 .or.  &
						   & cdakind/1000 == 14 .or. cdakind/1000 == 17.or. &
						   & cdakind == 5900) then
						height = real(IC3(9),r_sngl)    ! model surface altitue
					 elseif (cdakind/1000 == 16 ) then
						height = 10.0
						ztype = IVITP_HEIGHT
					 else
						cycle
					 end if
				  ! upper
				  else
					 if ((cdakind/1000) == 6 .or. (cdakind/1000) == 7 .or. (cdakind/1000) == 19) then
						do kh = 1, NHCD(nh)
						   if (KHCD(kh,nh) == KHCD_RPHT) then           
							  height = AHCD(kh,nh)       ! reported height
							  exit
						   end if
						end do
						ztype = IVITP_HEIGHT            ! interpolate with m
					 elseif (cdakind/1000 == 16) then
						height = 10.0
						ztype = IVITP_HEIGHT
					 else
						do kh = 1,NHCD(nh)
						   if (KHCD(kh,nh) == KHCD_IVUT) then   ! unit of height
							  if (AHCD(kh,nh) == 1) then    ! unit = hPa
								 height = HGT(nh)
								 ztype = IVITP_LOGP        ! interpolate with logP
								 exit
							  end if
						   end if
						end do
					 end if
					 ! height data is necessary for upper data
					 if (ztype < 0) cycle
					 if (int(height) == MISS) cycle
				  end if
				  select case(ztype)
				  case(IVITP_HEIGHT) ! vertical interpolation with height
					 height = height
				  case(IVITP_LOGP)   ! vertical interpolation with log_p
					 height = log(height*100.e0)
				  case(IVITP_GRID)   ! vertical interpolation with grid coordinate
					 height = 0.0d0
				  case(IVITP_SURF)   ! surface data: 2D interpolation
					 height = height
				  case default
				  end select
				  !
				  ! set obsereved value, observation error
				  do ne = 1, NEL(nh)
					 if (IUSE(ne,nh) /= 1) cycle             ! use only IUSE flag == 1
					 if (mel(ne,nh) /= cdacode) cycle
					 if (mel(ne,nh) == MPSMD .and. iflg_skip_psobs == 1) cycle
					 ! select surface observation
					 select case (iflg_surface_obs)
					 case (1)
						if (ztype == IVITP_SURF) then
						   if (cdakind/100 == 164 .and. (mel(ne,nh) == MWUCP .or. mel(ne,nh) == MWVCP)) then
						   else if((cdakind/1000 == 17 .or. cdakind/1000 == 14) .and. mel(ne,nh) == MPCWT) then
						   else if((cdakind/1000 == 1 .or. cdakind/1000 == 2) .and. mel(ne,nh) == MPSMD) then ! model surface pressure
						   else if (cdakind == 1400) then
						   else if (cdakind == 5900 ) then
						   else
							  cycle
						   end if
						end if
					 case default
					 end select
					 !
					 if (trim(self%name) == 'rvl') then
					 	call movec(IEND, 1, ic3(15), 1, 4)
#ifdef SWP_ON
                        call swp2swp4_l(IEND)
#endif
					    range = real(IEND,r_size)*0.3*pi/180.0d0
					    azimuth = real(ic3(17),r_size)*0.01d0
					    if (azimuth < 0.) azimuth = azimuth + 360.
					    azimuth = azimuth*pi/180.0d0
					 end if
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
					 ! relative humidity: [%] -> [no unit]
					 if (mel(ne,nh) == MRHMD) then  
						obs = obs/100.0d0
						error = error/100.0d0
						! pressure: [hPa] -> [Pa]
					 else if (mel(ne,nh) == MPSMD) then  ! model surface pressure
						obs = obs*100.0d0
						error = error*100.0d0
					 end if
					 !
					 iobs = iobs + 1
					 self%mobs(i,j,it) = self%mobs(i,j,it) + 1
					 self%cdakind(iobs) = cdakind; self%ztype(iobs) = ztype
					 i0(iobs) = i; j0(iobs) = j
					 self%lon(iobs) = lon; self%lat(iobs) = lat
					 self%height(iobs) = height
					 self%obs(iobs,1) = obs; self%error(iobs,1) = error
					 if (trim(self%name) == 'rvl') then
					    self%range(iobs) = range; self%azimuth(iobs) = azimuth
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
						   ztype = self%ztype(iobs1+j); self%ztype(iobs1+j) = self%ztype(iobs1+j-m); self%ztype(iobs1+j-m) = ztype
						   lon = self%lon(iobs1+j); self%lon(iobs1+j) = self%lon(iobs1+j-m); self%lon(iobs1+j-m) = lon
						   lat = self%lat(iobs1+j); self%lat(iobs1+j) = self%lat(iobs1+j-m); self%lat(iobs1+j-m) = lat
						   height = self%height(iobs1+j); self%height(iobs1+j) = self%height(iobs1+j-m); self%height(iobs1+j-m) = height
						   obs = self%obs(iobs1+j,1); self%obs(iobs1+j,1) = self%obs(iobs1+j-m,1); self%obs(iobs1+j-m,1) = obs
						   error = self%error(iobs1+j,1); self%error(iobs1+j,1) = self%error(iobs1+j-m,1); self%error(iobs1+j-m,1) = error
						   if (trim(self%name) == 'rvl') then
							  range = self%range(iobs1+j); self%range(iobs1+j) = self%range(iobs1+j-m); self%range(iobs1+j-m) = range
							  azimuth = self%azimuth(iobs1+j); self%azimuth(iobs1+j) = self%azimuth(iobs1+j-m); self%azimuth(iobs1+j-m) = azimuth
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
      call int_broadcast1D('all', self%ztype, 0, nobs)
      call real_broadcast1D('all', self%lon, 0, nobs)
      call real_broadcast1D('all', self%lat, 0, nobs)
      call real_broadcast1D('all', self%height, 0, nobs)
      call broadcast2D('all', self%obs, 0, nobs, nsubobs)
      call broadcast2D('all', self%error, 0, nobs, nsubobs)
      if (trim(self%name) == 'rvl') then
	     call broadcast1D('all', self%range, 0, nobs)
	     call broadcast1D('all', self%azimuth, 0, nobs)
      end if
      !
      return
   end subroutine cdaread_obs
   !
   !
   !
   subroutine read_obs_NodeObsSpaceFieldCNV(self, myid)
      implicit none
      type(NodeObsSpaceFieldCNV), intent(inout) :: self
      integer, intent(in) :: myid
      !
      if (trim(cnvobs_format) == 'cda') then
	     call cdaread_obs(self, myid)
      else if (trim(cnvobs_format) == 'txt') then
	     !call txtread_obs(self, myid)
      end if
      !
      return
   end subroutine read_obs_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine scatter_obs_NodeObsSpaceFieldCNV(self, info, local_object)
      ! scatter a big obsspace from Node 0 to all nodes.
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldCNV), intent(inout) :: local_object
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
		 allocate(local_object%cdakind(nobs), local_object%ztype(nobs))
		 allocate(local_object%lon(nobs), local_object%lat(nobs), local_object%height(nobs))
		 allocate(local_object%obs(nobs,nsubobs), local_object%error(nobs,nsubobs))
		 if (trim(self%name) == 'rvl') allocate(local_object%range(nobs), local_object%azimuth(nobs))
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
					 local_object%ztype(iobs1:iobs2) = self%ztype(iobs3:iobs4)
					 local_object%lon(iobs1:iobs2) = self%lon(iobs3:iobs4)
					 local_object%lat(iobs1:iobs2) = self%lat(iobs3:iobs4)
					 local_object%height(iobs1:iobs2) = self%height(iobs3:iobs4)
					 local_object%obs(iobs1:iobs2,:) = self%obs(iobs3:iobs4,:)
					 local_object%error(iobs1:iobs2,:) = self%error(iobs3:iobs4,:)
					 if (trim(self%name) == 'rvl') then
					    local_object%range(iobs1:iobs2) = self%range(iobs3:iobs4)
					    local_object%azimuth(iobs1:iobs2) = self%azimuth(iobs3:iobs4)
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
   end subroutine scatter_obs_NodeObsSpaceFieldCNV
   !
   !
   !
   subroutine gather_obs_NodeObsSpaceFieldCNV(self, info, global_object)
      ! this subroutine is only for completion and we realy do not use this
      implicit none
      type(NodeObsSpaceFieldCNV), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldCNV), intent(inout) :: global_object
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
   end subroutine gather_obs_NodeObsSpaceFieldCNV
   !
   !
   !
end module NodeObsSpaceFieldCNV_class
