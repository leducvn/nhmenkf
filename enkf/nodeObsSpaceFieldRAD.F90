module NodeObsSpaceFieldRAD_class
! Author: Le Duc
! Created date: 13 Mar 2016
   use variable, only : r_size, r_sngl, pi, iflg_incremental, radobs_format
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
   integer, parameter :: KECD_GUES = 600  ! first guess
   !
   integer, parameter :: MZHGT = 10*256+ 7 ! geopotential height
   integer, parameter :: MWUEW = 11*256+ 3 ! wind U (east-west)
   integer, parameter :: MWVSN = 11*256+ 4 ! wind V (south-north)
   integer, parameter :: MZSMD = 57*256+ 2 ! height over model surface
   integer, parameter :: MPSMD = 57*256+ 1 ! model surface pressure
   integer, parameter :: MTEMP = 12*256+ 1 ! temperature
   integer, parameter :: MSRAD = 12*256+62 ! satellite radiance
   integer, parameter :: MPRES = 10*256+ 4 ! Pressure
   integer, parameter :: MEMIS = 3825      ! emissivity
   !
   integer, parameter :: IVITP_SRAD   = 5  ! satellite radiance
   !
   ! Obs are stored as a linked list at each horizontal grid point. To access
   ! obs fast, all link lists are concatenated into a one-dimensional array. The
   ! iobs and mobs are the indices fo accessing: iobs is the index of each link
   ! list in this array (=0 if no link list) and mobs is the number of obs in each link list.
   ! We do not put all information like lon, lat, ... into a structure because these fields
   ! depend on observation types.
   type NodeObsSpaceFieldRAD
      character(len=10) :: name, filename
      ! For radiance obs, nsubobs = 100
      integer :: nobs, nsubobs, nx, ny, nt
      integer, dimension(:), allocatable :: cdakind, surftype, satid, nchannel
      integer, dimension(:,:), allocatable :: channel
      real(r_sngl), dimension(:), allocatable :: lon, lat
      real(r_size), dimension(:), allocatable :: sat_zenith, sat_azimuth, sst
      !real(r_size), dimension(:), allocatable :: us, vs, ps, ts
      real(r_size), dimension(:,:), allocatable :: obs, error
      integer, dimension(:,:,:), allocatable :: iobs, mobs
   end type NodeObsSpaceFieldRAD
   !
   interface new
      module procedure new_NodeObsSpaceFieldRAD
   end interface
   interface destroy
      module procedure destroy_NodeObsSpaceFieldRAD
   end interface
   interface display
      module procedure display_NodeObsSpaceFieldRAD
   end interface
   interface get_name
      module procedure get_name_NodeObsSpaceFieldRAD
   end interface
   interface get_nobs
      module procedure get_nobs_NodeObsSpaceFieldRAD
   end interface
   interface get_nsubobs
      module procedure get_nsubobs_NodeObsSpaceFieldRAD
   end interface
   interface get_mobs
      module procedure get_mobs_NodeObsSpaceFieldRAD1
      module procedure get_mobs_NodeObsSpaceFieldRAD2
   end interface
   interface read_obs
      module procedure read_obs_NodeObsSpaceFieldRAD
   end interface
   interface scatter_obs
      module procedure scatter_obs_NodeObsSpaceFieldRAD
   end interface
   interface gather_obs
      module procedure gather_obs_NodeObsSpaceFieldRAD
   end interface
   !
contains
   !
   subroutine new_NodeObsSpaceFieldRAD(self, name, filename, nx, ny, nt)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(inout) :: self
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
   end subroutine new_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine destroy_NodeObsSpaceFieldRAD(self)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(inout) :: self
      integer :: iobs
      !
      if (allocated(self%iobs)) deallocate(self%iobs)
      if (allocated(self%mobs)) deallocate(self%mobs)
      if (allocated(self%cdakind)) deallocate(self%cdakind)
      if (allocated(self%surftype)) deallocate(self%surftype)
      if (allocated(self%satid)) deallocate(self%satid)
      if (allocated(self%lon)) deallocate(self%lon)
      if (allocated(self%lat)) deallocate(self%lat)
      if (allocated(self%nchannel)) deallocate(self%nchannel)
      if (allocated(self%channel)) deallocate(self%channel)
      if (allocated(self%sat_zenith)) deallocate(self%sat_zenith)
      if (allocated(self%sat_azimuth)) deallocate(self%sat_azimuth)
      if (allocated(self%sst)) deallocate(self%sst)
      !if (allocated(self%us)) deallocate(self%us)
      !if (allocated(self%vs)) deallocate(self%vs)
      !if (allocated(self%ps)) deallocate(self%ps)
      !if (allocated(self%ts)) deallocate(self%ts)
      if (allocated(self%obs)) deallocate(self%obs)
      if (allocated(self%error)) deallocate(self%error)
      !
      return
   end subroutine destroy_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine display_NodeObsSpaceFieldRAD(self)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      !
      print*, 'Observation name: ', self%name
      print*, 'File name: ', self%filename
      print*, 'Number: ', self%nobs, self%nsubobs
      !
      return
   end subroutine display_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine get_name_NodeObsSpaceFieldRAD(self, name)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      character(len=10), intent(out) :: name
      !
      name = self%name
      !
      return
   end subroutine get_name_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine get_nobs_NodeObsSpaceFieldRAD(self, nobs)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      integer, intent(out) :: nobs
      !
      nobs = self%nobs
      !
      return
   end subroutine get_nobs_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine get_nsubobs_NodeObsSpaceFieldRAD(self, nsubobs)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      integer, intent(out) :: nsubobs
      !
      nsubobs = self%nsubobs
      !
      return
   end subroutine get_nsubobs_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldRAD1(self, i, j, mobs)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      integer, intent(in) :: i, j
      integer, intent(out) :: mobs
      integer :: it, iobs, iobs1, iobs2
      !
      !mobs = sum(self%mobs(i,j,:))*self%nsubobs
      mobs = 0
      do it = 1, self%nt
         iobs1 = self%iobs(i,j,it)
         iobs2 = self%iobs(i,j,it) + self%mobs(i,j,it) - 1
         do iobs = iobs1, iobs2
            mobs = mobs + self%nchannel(iobs)
         end do
      end do
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldRAD1
   !
   !
   !
   subroutine get_mobs_NodeObsSpaceFieldRAD2(self, i, j, it, mobs)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      integer, intent(in) :: i, j, it
      integer, intent(out) :: mobs
      integer :: iobs, iobs1, iobs2
      !
      !mobs = self%mobs(i,j,it)*self%nsubobs
      mobs = 0
      iobs1 = self%iobs(i,j,it)
      iobs2 = self%iobs(i,j,it) + self%mobs(i,j,it) - 1
      do iobs = iobs1, iobs2
         mobs = mobs + self%nchannel(iobs)
      end do
      !
      return
   end subroutine get_mobs_NodeObsSpaceFieldRAD2
   !
   !
   !
   subroutine cdaread_nobs(self, myid, nobs)
      use rttovlib, only : get_satid
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos, proj
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      integer, intent(in) :: myid
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
      character(2) :: infile = '00'
      character(len=8) :: satname
      integer :: nx, ny, nt, cdacode, cdakind, ztype, nsubobs, satid
      integer :: it, nh, kh, ne, ke, i, j, k, ierror, ichannel
      real(kind=r_sngl) :: lat, lon, height  ! lat, lon, rhgt must be single
      type(PosGrid) :: pg
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      cdacode = MSRAD
      !
      nobs = 0
      nsubobs = 100
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
!#ifdef SWP_ON
                        do k = 3, 6
	                   call swap_c2_l(IC3(k))
	                end do
!#endif
			call movec(satname, 1, IC3(3), 1, 8)
                        call get_satid(satname, cdakind, satid)
                        !if (satid == 0) print*, satname, cdakind
                        if (satid == 0) cycle
			!
			! unpack CDA part 4
			call CDAUNPK(IC1, IC2, IC3, IC4, NHT, HGT, NEL, MEL, &
						  & NHCD, KHCD, AHCD, NECD, KECD, AECD, &
						  & JGQC, IUSE, NOQC, JEQC, KEQC, &
						  & MISS, NXHT, NXEL, NXIH, NXIE, NXQC)
			ichannel = 0
	        LOOP_LEVEL: do nh = 1, NHT-1
	           IVNM = 0
	           do kh = 1, NHCD(nh)
		          if (KHCD(kh,nh) == KHCD_IVNM) then   ! vertical indicator
		             IVNM = AHCD(kh,nh)
		             exit
		          end if
	           end do
               !
               do kh = 1,NHCD(nh)
                  if (AHCD(kh,nh) == 3) then    ! unit = channel number (satellite)
                     height = HGT(nh)
                     ztype = IVITP_SRAD        ! temporary only 2D interpolation? !!!!!
                  end if
               end do
		       ! height data is necessary for radiance data
		       if (int(height) == MISS) cycle
               !
	           ! set obsereved value, observation error
               do ne = 1, NEL(nh)
                  if (IUSE(ne,nh) /= 1) cycle             ! use only IUSE flag == 1
                  if (mel(ne,nh) /= cdacode) cycle
                  if (ichannel == nsubobs) exit
                  ichannel = ichannel + 1
	           end do
            end do LOOP_LEVEL
            if (ichannel > 0) nobs = nobs + 1
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
      use rttovlib, only : get_satid, idtype_ssmi, idtype_tmi, idtype_amse, idtype_gmi, idtype_airs, idtype_iasi, idtype_cris
      use interpolate, only : PosGrid, convert_LatLon_to_GridPos, proj
      implicit none
      type(NodeObsSpaceFieldRAD), intent(inout) :: self
      integer, intent(in) :: myid
      !
      integer(2) :: IC1(1:LX1)
      integer(2) :: IC2(1:LX2)
      integer(2) :: IC3(1:LX3)
      integer(2) :: IC4(1:LX4)
      !
      integer(4) :: NHT, IVNM, IEND
      integer(4) :: ifound_obqc, ifound_ober, ifound_dval
      real   (4) :: OBQC, OBER, DVAL
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
      character(len=8) :: satname
      integer :: cdacode, cdakind, ztype, surftype, satid
      integer :: nx, ny, nt, nobs, nsubobs, nchannel, ierror
      integer :: it, nh, kh, ne, ke, i, j, k, iobs, iobs1, iobs2, m, itmp, jtmp, ichannel
      real(kind=r_sngl) :: lat, lon, height, tskin
      real(r_size) :: sat_zenith, sat_azimuth, sst, us, vs, ts, ps
      type(PosGrid) :: pg
      integer, dimension(:), allocatable :: i0, j0, channel
      real(r_size), dimension(:), allocatable :: obs, error
      !
      nx = self%nx
      ny = self%ny
      nt = self%nt
      !
      if (myid == 0) call cdaread_nobs(self, myid, nobs)
      call int_broadcast0D('all', nobs, 0)
      self%nobs = nobs
      nsubobs = 100
      self%nsubobs = nsubobs
      if (nobs == 0) return
      allocate(self%cdakind(nobs), self%surftype(nobs), self%satid(nobs))
      allocate(self%nchannel(nobs), self%channel(nobs,nsubobs))
      allocate(self%lon(nobs), self%lat(nobs))
      allocate(self%sat_zenith(nobs), self%sat_azimuth(nobs))
      allocate(self%sst(nobs))
      !allocate(self%us(nobs), self%vs(nobs), self%ps(nobs), self%ts(nobs))
      allocate(self%obs(nobs,nsubobs), self%error(nobs,nsubobs))
      self%obs(:,:) = 0.d0
      self%error(:,:) = 0.d0
      !
      if (myid == 0) then
         allocate(i0(nobs), j0(nobs), channel(nsubobs))
         allocate(obs(nsubobs), error(nsubobs))
         cdacode = MSRAD
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
               ! satellite information
!#ifdef SWP_ON
               do k = 3, 6
	          call swap_c2_l(IC3(k))
	       end do
!#endif
               call movec(satname, 1, IC3(3), 1, 8)
               call get_satid(satname, cdakind, satid)
               if (satid == 0) cycle
               surftype = IC3(7)
	       if (cdakind == idtype_ssmi .or. cdakind == idtype_tmi .or. &
	         & cdakind == idtype_amse .or. cdakind == idtype_gmi) surftype = 1 ! always sea
               !nscan = IC3(10)
               sat_zenith  = real(IC3(11))/100.d0
               sat_azimuth = real(IC3(12))/100.d0
               !sun_zenith  = 0.0d0
               !sun_azimuth = 0.0d0
               !call movec(tskin, 1, IC3(17), 1, 4)
!#ifdef SWP_ON
               !call swp2swp4_l(tskin)
!#endif
               !clw = tskin ! clw [g/m2]
               !call movec(tskin, 1, IC3(31), 1, 4)
!#ifdef SWP_ON
               !call swp2swp4_l(tskin)
!#endif
               !tpw = tskin ! guess tpw [kg/m2]
               call movec(tskin, 1, IC3(23), 1, 4)
#ifdef SWP_ON
               call swp2swp4_l(tskin)
#endif
               sst = tskin ! sst [K]
               !
               ! unpack CDA part 4
               call CDAUNPK(IC1, IC2, IC3, IC4, NHT, HGT, NEL, MEL, &
              &             NHCD, KHCD, AHCD, NECD, KECD, AECD, &
              &             JGQC, IUSE, NOQC, JEQC, KEQC, &
              &             MISS, NXHT, NXEL, NXIH, NXIE, NXQC)
               ichannel = 0
               LOOP_LEVEL: do nh = 1, NHT-1
                  IVNM = 0
                  do kh = 1, NHCD(nh)
                     if (KHCD(kh,nh) == KHCD_IVNM) then ! vertical indicator
                        IVNM = AHCD(kh,nh)
                        exit
                     end if
                  end do
                  !
                  do kh = 1, NHCD(nh)
		     if (AHCD(kh,nh) == 3) then ! unit = channel number (satellite)           
			height = HGT(nh)
			ztype = IVITP_SRAD ! temporary only 2D interpolation? !!!!!
		     end if
		  end do
		  ! height data is necessary for sattelite data
		  if (int(height) == MISS) cycle
		  !
		  ! set obsereved value, observation error
		  do ne = 1, NEL(nh)
		     if (IUSE(ne,nh) /= 1) cycle             ! use only IUSE flag == 1
		     if (mel(ne,nh) /= cdacode) cycle
		     if (ichannel == nsubobs) exit
		     ichannel = ichannel + 1
		     if (cdakind == idtype_airs .or. cdakind == idtype_iasi .or. cdakind == idtype_cris) then
		        channel(ichannel) = int(height)
		     else
		        channel(ichannel) = nh
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
			obs(ichannel) = dval
		     else
		        obs(ichannel) = obqc
		     end if
		     error(ichannel) = ober
		  end do
	       end do LOOP_LEVEL
	       nchannel = ichannel
	       if (nchannel == 0) cycle LOOP_RECORD
			!
			nh = NHT
			   LOOP_ELEM_SURF: do ne = 1, NEL(nh)
				  ! SURF_U
				  if (mel(ne,nh) == MWUEW) then
					 do ke = 1, NECD(ne,nh)
						if (kecd(ke,ne,nh) == KECD_OBQC) then
						   us = aecd(ke,ne,nh)
						end if
					 end do
				  ! SURF_V
				  else if (mel(ne,nh) == MWVSN) then
					 do ke = 1, NECD(ne,nh)
						if (kecd(ke,ne,nh) == KECD_OBQC) then
						   vs = aecd(ke,ne,nh)
						end if
					 end do
				  ! SURF_P
				  else if (mel(ne,nh) == MPSMD) then
					 do ke = 1, NECD(ne,nh)
						if (kecd(ke,ne,nh) == KECD_OBQC) then
						   ps = aecd(ke,ne,nh)
						end if 
					 end do
				  ! SURF_T
				  else if (mel(ne,nh) == MTEMP) then
					 do ke = 1, NECD(ne,nh)
					   if (kecd(ke,ne,nh) == KECD_OBQC) then
						 ts = aecd(ke,ne,nh)
					   end if 
					 end do
				  end if 
			   end do LOOP_ELEM_SURF
			   !
			   iobs = iobs + 1
			   self%mobs(i,j,it) = self%mobs(i,j,it) + 1
			   self%cdakind(iobs) = cdakind; self%surftype(iobs) = surftype; self%satid(iobs) = satid
			   self%nchannel(iobs) = nchannel; self%channel(iobs,1:nchannel) = channel(1:nchannel)
			   i0(iobs) = i; j0(iobs) = j
			   self%lon(iobs) = lon; self%lat(iobs) = lat
			   self%sat_zenith(iobs) = sat_zenith; self%sat_azimuth(iobs) = sat_azimuth
			   self%sst(iobs) = sst
			   !self%us(iobs) = us; self%vs(iobs) = vs; self%ps(iobs) = ps; self%ts(iobs) = ts
			   self%obs(iobs,1:nchannel) = obs(1:nchannel); self%error(iobs,1:nchannel) = error(1:nchannel)
			   !print*, cdakind, satid, surftype, nchannel, sat_zenith, sat_azimuth, sst, lon, lat, obs(1), error(1)
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
						   surftype = self%surftype(iobs1+j); self%surftype(iobs1+j) = self%surftype(iobs1+j-m); self%surftype(iobs1+j-m) = surftype
						   satid = self%satid(iobs1+j); self%satid(iobs1+j) = self%satid(iobs1+j-m); self%satid(iobs1+j-m) = satid
						   nchannel = self%nchannel(iobs1+j); self%nchannel(iobs1+j) = self%nchannel(iobs1+j-m); self%nchannel(iobs1+j-m) = nchannel
						   channel(:) = self%channel(iobs1+j,:); self%channel(iobs1+j,:) = self%channel(iobs1+j-m,:); self%channel(iobs1+j-m,:) = channel(:)
						   lon = self%lon(iobs1+j); self%lon(iobs1+j) = self%lon(iobs1+j-m); self%lon(iobs1+j-m) = lon
						   lat = self%lat(iobs1+j); self%lat(iobs1+j) = self%lat(iobs1+j-m); self%lat(iobs1+j-m) = lat
						   sat_zenith = self%sat_zenith(iobs1+j); self%sat_zenith(iobs1+j) = self%sat_zenith(iobs1+j-m); self%sat_zenith(iobs1+j-m) = sat_zenith
						   sat_azimuth = self%sat_azimuth(iobs1+j); self%sat_azimuth(iobs1+j) = self%sat_azimuth(iobs1+j-m); self%sat_azimuth(iobs1+j-m) = sat_azimuth
						   sst = self%sst(iobs1+j); self%sst(iobs1+j) = self%sst(iobs1+j-m); self%sst(iobs1+j-m) = sst
						   !us = self%us(iobs1+j); self%us(iobs1+j) = self%us(iobs1+j-m); self%us(iobs1+j-m) = us
						   !vs = self%vs(iobs1+j); self%vs(iobs1+j) = self%vs(iobs1+j-m); self%vs(iobs1+j-m) = vs
						   !ps = self%ps(iobs1+j); self%ps(iobs1+j) = self%ps(iobs1+j-m); self%ps(iobs1+j-m) = ps
						   !ts = self%ts(iobs1+j); self%ts(iobs1+j) = self%ts(iobs1+j-m); self%ts(iobs1+j-m) = ts
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
      call int_broadcast1D('all', self%cdakind, 0, nobs)
      call int_broadcast1D('all', self%surftype, 0, nobs)
      call int_broadcast1D('all', self%satid, 0, nobs)
      call int_broadcast1D('all', self%nchannel, 0, nobs)
      call int_broadcast2D('all', self%channel, 0, nobs, nsubobs)
      call real_broadcast1D('all', self%lon, 0, nobs)
      call real_broadcast1D('all', self%lat, 0, nobs)
      call broadcast1D('all', self%sat_zenith, 0, nobs)
      call broadcast1D('all', self%sat_azimuth, 0, nobs)
      call broadcast1D('all', self%sst, 0, nobs)
      !call broadcast1D('all', self%us, 0, nobs)
      !call broadcast1D('all', self%vs, 0, nobs)
      !call broadcast1D('all', self%ps, 0, nobs)
      !call broadcast1D('all', self%ts, 0, nobs)
      call broadcast2D('all', self%obs, 0, nobs, nsubobs)
      call broadcast2D('all', self%error, 0, nobs, nsubobs)
      !
      return
   end subroutine cdaread_obs
   !
   !
   !
   subroutine read_obs_NodeObsSpaceFieldRAD(self, myid)
      implicit none
      type(NodeObsSpaceFieldRAD), intent(inout) :: self
      integer, intent(in) :: myid
      !
      if (trim(radobs_format) == 'cda') then
	     call cdaread_obs(self, myid)
      else if (trim(radobs_format) == 'txt') then
	     !call txtread_obs(self, myid)
      end if
      !
      return
   end subroutine read_obs_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine scatter_obs_NodeObsSpaceFieldRAD(self, info, local_object)
      ! scatter a big obsspace from Node 0 to all nodes.
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldRAD), intent(inout) :: local_object
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
		 allocate(local_object%cdakind(nobs), local_object%surftype(nobs), local_object%satid(nobs))
		 allocate(local_object%nchannel(nobs), local_object%channel(nobs,nsubobs))
		 allocate(local_object%lon(nobs), local_object%lat(nobs))
		 allocate(local_object%sat_zenith(nobs), local_object%sat_azimuth(nobs), local_object%sst(nobs))
		 !allocate(local_object%us(nobs), local_object%vs(nobs), local_object%ps(nobs), local_object%ts(nobs))
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
					 local_object%cdakind(iobs1:iobs2) = self%cdakind(iobs3:iobs4)
					 local_object%surftype(iobs1:iobs2) = self%surftype(iobs3:iobs4)
					 local_object%satid(iobs1:iobs2) = self%satid(iobs3:iobs4)
					 local_object%nchannel(iobs1:iobs2) = self%nchannel(iobs3:iobs4)
					 local_object%channel(iobs1:iobs2,:) = self%channel(iobs3:iobs4,:)
					 local_object%lon(iobs1:iobs2) = self%lon(iobs3:iobs4)
					 local_object%lat(iobs1:iobs2) = self%lat(iobs3:iobs4)
					 local_object%sat_zenith(iobs1:iobs2) = self%sat_zenith(iobs3:iobs4)
					 local_object%sat_azimuth(iobs1:iobs2) = self%sat_azimuth(iobs3:iobs4)
					 local_object%sst(iobs1:iobs2) = self%sst(iobs3:iobs4)
					 !local_object%us(iobs1:iobs2) = self%us(iobs3:iobs4)
					 !local_object%vs(iobs1:iobs2) = self%vs(iobs3:iobs4)
					 !local_object%ps(iobs1:iobs2) = self%ps(iobs3:iobs4)
					 !local_object%ts(iobs1:iobs2) = self%ts(iobs3:iobs4)
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
   end subroutine scatter_obs_NodeObsSpaceFieldRAD
   !
   !
   !
   subroutine gather_obs_NodeObsSpaceFieldRAD(self, info, global_object)
      ! this subroutine is only for completion and we realy do not use this
      implicit none
      type(NodeObsSpaceFieldRAD), intent(in) :: self
      type(NodeInfo), intent(in) :: info
      type(NodeObsSpaceFieldRAD), intent(inout) :: global_object
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
   end subroutine gather_obs_NodeObsSpaceFieldRAD
   !
   !
   !
end module NodeObsSpaceFieldRAD_class
