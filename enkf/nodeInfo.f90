module NodeInfo_class
! Author: Le Duc
! Created date: 27 Jul 2015
   implicit none
   !
   type NodeInfo
      integer :: indexmode, myid, myidx, myidy, myide, nproc, nxpe, nype, nepe
      integer :: nx0, ny0, ne0, nx, ny, ne, dis, die, djs, dje, halo
      integer, dimension(:), allocatable :: xindex, yindex, eindex
      integer, dimension(:), allocatable :: istart, iend, jstart, jend
   end type NodeInfo
   !
   interface new
      module procedure new_NodeInfo1
      module procedure new_NodeInfo2
   end interface
   interface destroy
      module procedure destroy_NodeInfo
   end interface
   interface display
      module procedure display_NodeInfo
   end interface
   interface get_halo
      module procedure get_halo_NodeInfo
   end interface
   interface get_nx
      module procedure get_nx_NodeInfo
   end interface
   interface get_ny
      module procedure get_ny_NodeInfo
   end interface
   interface get_ne
      module procedure get_ne_NodeInfo
   end interface
   interface get_nx0
      module procedure get_nx0_NodeInfo
   end interface
   interface get_ny0
      module procedure get_ny0_NodeInfo
   end interface
   interface get_ne0
      module procedure get_ne0_NodeInfo
   end interface
   interface get_myid
      module procedure get_myid_NodeInfo
   end interface
   interface get_myidx
      module procedure get_myidx_NodeInfo
   end interface
   interface get_myidy
      module procedure get_myidy_NodeInfo
   end interface
   interface get_myide
      module procedure get_myide_NodeInfo
   end interface
   interface get_di
      module procedure get_di_NodeInfo
   end interface
   interface get_dj
      module procedure get_dj_NodeInfo
   end interface
   interface get_nxpe
      module procedure get_nxpe_NodeInfo
   end interface
   interface get_nype
      module procedure get_nype_NodeInfo
   end interface
   interface get_nepe
      module procedure get_nepe_NodeInfo
   end interface
   !
contains
   !
   subroutine new_NodeInfo1(self, indexmode, nxpe, nype, nepe, myid, nx0, ny0, ne0, halo)
      implicit none
      type(NodeInfo), intent(inout) :: self
      integer, intent(in) :: indexmode, nxpe, nype, nepe, myid, nx0, ny0, ne0, halo
      integer :: myidx, myidy, myide, nx, ny, ne
      integer :: is, ie, js, je, dis, die, djs, dje, ks, ke, ip
      !
      self%indexmode = indexmode
      self%nxpe = nxpe
      self%nype = nype
      self%nepe = nepe
      self%nproc = nxpe*nype*nepe
      self%nx0 = nx0
      self%ny0 = ny0
      self%ne0 = ne0
      self%halo = halo
      call get_info(self, myid, myidx, myidy, is, ie, js, je, dis, die, djs, dje)
      call get_info_nohalo(self, myid, myide, ks, ke)
      self%myid = myid
      self%myidx = myidx
      self%myidy = myidy
      self%myide = myide
      self%dis = dis
      self%die = die
      self%djs = djs
      self%dje = dje
      !
      nx = ie - is + 1
      self%nx = nx
      allocate(self%xindex(nx))
      do ip = is, ie
	 self%xindex(ip-is+1) = ip
      end do
      !
      ny = je - js + 1
      self%ny = ny
      allocate(self%yindex(ny))
      do ip = js, je
	 self%yindex(ip-js+1) = ip
      end do
      !
      ne = ke - ks + 1
      self%ne = ne
      allocate(self%eindex(ne))
      do ip = ks, ke
	 self%eindex(ip-ks+1) = ip
      end do
      !
      return
   end subroutine new_NodeInfo1
   !
   !
   !
   subroutine new_NodeInfo2(self, indexmode, nxpe, nype, nepe, myid, nx0, ny0, ne0, halo, istart, iend, jstart, jend)
      implicit none
      type(NodeInfo), intent(inout) :: self
      integer, intent(in) :: indexmode, nxpe, nype, nepe, myid, nx0, ny0, ne0, halo
      integer, dimension(nxpe), intent(in) :: istart, iend
      integer, dimension(nype), intent(in) :: jstart, jend
      integer :: myidx, myidy, myide, nx, ny, ne
      integer :: is, ie, js, je, dis, die, djs, dje, ks, ke, ip
      !
      self%indexmode = indexmode
      self%nxpe = nxpe; self%nype = nype; self%nepe = nepe
      self%nproc = nxpe*nype*nepe
      self%nx0 = nx0; self%ny0 = ny0; self%ne0 = ne0
      self%halo = halo
      allocate(self%istart(nxpe), self%iend(nxpe), self%jstart(nype), self%jend(nype))
      self%istart(:) = istart(:); self%iend(:) = iend(:)
      self%jstart(:) = jstart(:); self%jend(:) = jend(:)
      call get_info(self, myid, myidx, myidy, is, ie, js, je, dis, die, djs, dje)
      call get_info_nohalo(self, myid, myide, ks, ke)
      self%myid = myid
      self%myidx = myidx; self%myidy = myidy; self%myide = myide
      self%dis = dis; self%die = die
      self%djs = djs; self%dje = dje
      !
      nx = ie - is + 1
      self%nx = nx
      allocate(self%xindex(nx))
      do ip = is, ie
	 self%xindex(ip-is+1) = ip
      end do
      !
      ny = je - js + 1
      self%ny = ny
      allocate(self%yindex(ny))
      do ip = js, je
	 self%yindex(ip-js+1) = ip
      end do
      !
      ne = ke - ks + 1
      self%ne = ne
      allocate(self%eindex(ne))
      do ip = ks, ke
	 self%eindex(ip-ks+1) = ip
      end do
      !
      return
   end subroutine new_NodeInfo2
   !
   !
   !
   subroutine destroy_NodeInfo(self)
      implicit none
      type(NodeInfo), intent(inout) :: self
      !
      if (allocated(self%xindex)) deallocate(self%xindex)
      if (allocated(self%yindex)) deallocate(self%yindex)
      if (allocated(self%eindex)) deallocate(self%eindex)
      !
      return
   end subroutine destroy_NodeInfo
   !
   !
   !
   subroutine display_NodeInfo(self)
      implicit none
      type(NodeInfo), intent(in) :: self
      !
      print*, 'Node: ', self%nproc, self%myid
      print*, 'Dimension0: ', self%nx0, self%ny0, self%ne0
      print*, 'Dimension: ', self%nx, self%ny, self%ne, self%halo
      !
      return
   end subroutine display_NodeInfo
   !
   !
   !
   subroutine get_halo_NodeInfo(self, halo)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: halo
      !
      halo = self%halo
      !
      return
   end subroutine get_halo_NodeInfo
   !
   !
   !
   subroutine get_nxpe_NodeInfo(self, nxpe)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: nxpe
      !
      nxpe = self%nxpe
      !
      return
   end subroutine get_nxpe_NodeInfo
   !
   !
   !
   subroutine get_nype_NodeInfo(self, nype)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: nype
      !
      nype = self%nype
      !
      return
   end subroutine get_nype_NodeInfo
   !
   !
   !
   subroutine get_nepe_NodeInfo(self, nepe)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: nepe
      !
      nepe = self%nepe
      !
      return
   end subroutine get_nepe_NodeInfo
   !
   !
   !
   subroutine get_nx_NodeInfo(self, nx)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: nx
      !
      nx = self%nx
      !
      return
   end subroutine get_nx_NodeInfo
   !
   !
   !
   subroutine get_ny_NodeInfo(self, ny)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: ny
      !
      ny = self%ny
      !
      return
   end subroutine get_ny_NodeInfo
   !
   !
   !
   subroutine get_ne_NodeInfo(self, ne)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: ne
      !
      ne = self%ne
      !
      return
   end subroutine get_ne_NodeInfo
   !
   !
   !
   subroutine get_nx0_NodeInfo(self, nx0)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: nx0
      !
      nx0 = self%nx0
      !
      return
   end subroutine get_nx0_NodeInfo
   !
   !
   !
   subroutine get_ny0_NodeInfo(self, ny0)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: ny0
      !
      ny0 = self%ny0
      !
      return
   end subroutine get_ny0_NodeInfo
   !
   !
   !
   subroutine get_ne0_NodeInfo(self, ne0)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: ne0
      !
      ne0 = self%ne0
      !
      return
   end subroutine get_ne0_NodeInfo
   !
   !
   !
   subroutine get_myid_NodeInfo(self, myid, nproc)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: myid, nproc
      !
      myid = self%myid
      nproc = self%nproc
      !
      return
   end subroutine get_myid_NodeInfo
   !
   !
   !
   subroutine get_myidx_NodeInfo(self, myidx)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: myidx
      !
      myidx = self%myidx
      !
      return
   end subroutine get_myidx_NodeInfo
   !
   !
   !
   subroutine get_myidy_NodeInfo(self, myidy)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: myidy
      !
      myidy = self%myidy
      !
      return
   end subroutine get_myidy_NodeInfo
   !
   !
   !
   subroutine get_myide_NodeInfo(self, myide)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: myide
      !
      myide = self%myide
      !
      return
   end subroutine get_myide_NodeInfo
   !
   !
   !
   subroutine get_di_NodeInfo(self, dis, die)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: dis, die
      !
      dis = self%dis
      die = self%die
      !
      return
   end subroutine get_di_NodeInfo
   !
   !
   !
   subroutine get_dj_NodeInfo(self, djs, dje)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(out) :: djs, dje
      !
      djs = self%djs
      dje = self%dje
      !
      return
   end subroutine get_dj_NodeInfo
   !
   !
   !
   subroutine get_xindex(self, ip, ix)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(in) :: ip
      integer, intent(out) :: ix
      !
      ix = self%xindex(ip)
      !
      return
   end subroutine get_xindex
   !
   !
   !
   subroutine get_yindex(self, ip, jy)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(in) :: ip
      integer, intent(out) :: jy
      !
      jy = self%yindex(ip)
      !
      return
   end subroutine get_yindex
   !
   !
   !
   subroutine get_eindex(self, ip, ke)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(in) :: ip
      integer, intent(out) :: ke
      !
      ke = self%eindex(ip)
      !
      return
   end subroutine get_eindex
   !
   !
   !
   subroutine get_info(self, id, idx, idy, is, ie, js, je, dis, die, djs, dje)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(in) :: id
      integer, intent(out) :: idx, idy, is, ie, js, je, dis, die, djs, dje
      integer :: nxpe, nype, myid
      !
      nxpe = self%nxpe
      nype = self%nype
      myid = id
      idx = mod(myid,nxpe)
      myid = myid/nxpe
      idy = mod(myid,nype)
      !
      if (self%indexmode == 0) then
         call set_index0(self%nx0, nxpe, idx, self%halo, self%istart, self%iend, is, ie, dis, die)
         call set_index0(self%ny0, nype, idy, self%halo, self%jstart, self%jend, js, je, djs, dje)
      else if (self%indexmode == 1) then
         call set_index1(self%nx0, nxpe, idx, self%halo, is, ie, dis, die)
         call set_index1(self%ny0, nype, idy, self%halo, js, je, djs, dje)
      else if (self%indexmode == 2) then
         call set_index2(self%nx0, nxpe, idx, self%halo, is, ie, dis, die)
         call set_index2(self%ny0, nype, idy, self%halo, js, je, djs, dje)
      else
         call set_index3(self%nx0, nxpe, idx, self%halo, is, ie, dis, die)
         call set_index3(self%ny0, nype, idy, self%halo, js, je, djs, dje)
      end if
      !
      return
   end subroutine get_info
   !
   !
   !
   subroutine get_info_nohalo(self, id, ide, ks, ke)
      implicit none
      type(NodeInfo), intent(in) :: self
      integer, intent(in) :: id
      integer, intent(out) :: ide, ks, ke
      integer :: nxpe, nype, nepe, myid
      !
      nxpe = self%nxpe
      nype = self%nype
      nepe = self%nepe
      myid = id
      !idx = mod(myid,nxpe)
      myid = myid/nxpe
      !idy = mod(myid,nype)
      ide = myid/nype
      !
      call set_index_nohalo(self%ne0, nepe, ide, ks, ke)
      !
      return
   end subroutine get_info_nohalo
   !
   !
   !
   subroutine load_balance(nxpe, nype, halo, nobs, is, ie, js, je, nx, ny) ! halo: radius of influence
      implicit none
      integer, intent(in) :: nxpe, nype, halo, nx, ny
      integer, dimension(nx,ny), intent(in) :: nobs
      integer, dimension(nxpe), intent(out) :: is, ie
      integer, dimension(nype), intent(out) :: js, je
      logical :: found
      integer :: nproc, n, m, id, idx, idy, ds, de, imin, imax, jmin, jmax, idmin, idxmin, idymin, k
      real :: score, score0
      real, dimension(nxpe*nype) :: nflop, nflop0
      !
      nproc = nxpe*nype
      do idx = 1, nxpe
         call set_index1(nx, nxpe, idx-1, 0, is(idx), ie(idx), ds, de)
      end do
      do idy = 1, nype
         call set_index1(ny, nype, idy-1, 0, js(idy), je(idy), ds, de)
      end do
      do id = 1, nproc
         idx = mod(id-1,nxpe) + 1; idy = (id-1)/nxpe + 1
	 imin = max(1,is(idx)-halo); imax = min(nx,ie(idx)+halo)
	 jmin = max(1,js(idy)-halo); jmax = min(ny,je(idy)+halo)
	 n = (ie(idx)-is(idx)+1)*(je(idy)-js(idy)+1)
	 m = sum(nobs(imin:imax,jmin:jmax))
	 if (m == 0) m = 1
	 nflop(id) = log10(1.d0*n) + 2.d0*log10(1.d0*m)
	 !print*, idx, idy, ie(idx)-is(idx)+1, je(idy)-js(idy)+1, nflop(id)
      end do
      !
      !score = maxval(nflop) - minval(nflop)
      score = sqrt(sum(nflop**2)/nproc-sum(nflop)**2/nproc**2)
      do while (.True.)
         idmin = minloc(nflop,dim=1)
	 idxmin = mod(idmin-1,nxpe) + 1; idymin = (idmin-1)/nxpe + 1
	 !print*, idxmin, idymin, score
	 !
	 found = .False.
	 do k = 1, 4
	    ! Left
	    if (k == 1) then
	       if (idxmin == 1) then
		  score0 = score
	       else if (ie(idxmin-1) == is(idxmin-1)) then
		  score0 = score
	       else
		  is(idxmin) = is(idxmin) - 1; ie(idxmin-1) = ie(idxmin-1) - 1
		  do idx = idxmin-1, idxmin
		  do idy = 1, nype
		     id = (idy-1)*nxpe + idx-1 + 1
		     imin = max(1,is(idx)-halo); imax = min(nx,ie(idx)+halo)
		     jmin = max(1,js(idy)-halo); jmax = min(ny,je(idy)+halo)
		     n = (ie(idx)-is(idx)+1)*(je(idy)-js(idy)+1)
		     m = sum(nobs(imin:imax,jmin:jmax))
		     if (m == 0) m = 1
		     nflop0(id) = nflop(id)
		     nflop(id) = log10(1.d0*n) + 2.d0*log10(1.d0*m)
		  end do
		  end do
		  !score0 = maxval(nflop) - minval(nflop)
		  score0 = sqrt(sum(nflop**2)/nproc-sum(nflop)**2/nproc**2)
		  if (score0 >= score) then
		     is(idxmin) = is(idxmin) + 1; ie(idxmin-1) = ie(idxmin-1) + 1
	             do idx = idxmin-1, idxmin
		     do idy = 1, nype
		        id = (idy-1)*nxpe + idx-1 + 1
		        nflop(id) = nflop0(id)
		     end do
		     end do
		  end if
	       end if
	       if (score0 < score) then
		  found = .True.
		  score = score0
		  exit
	       end if
	    ! Right
	    else if (k == 2) then
	       if (idxmin == nxpe) then
		  score0 = score
	       else if (ie(idxmin+1) == is(idxmin+1)) then
		  score0 = score
	       else
		  ie(idxmin) = ie(idxmin) + 1; is(idxmin+1) = is(idxmin+1) + 1
	          do idx = idxmin, idxmin+1
		  do idy = 1, nype
		     id = (idy-1)*nxpe + idx-1 + 1
		     imin = max(1,is(idx)-halo); imax = min(nx,ie(idx)+halo)
		     jmin = max(1,js(idy)-halo); jmax = min(ny,je(idy)+halo)
		     n = (ie(idx)-is(idx)+1)*(je(idy)-js(idy)+1)
		     m = sum(nobs(imin:imax,jmin:jmax))
		     if (m == 0) m = 1
		     nflop0(id) = nflop(id)
		     nflop(id) = log10(1.d0*n) + 2.d0*log10(1.d0*m)
		  end do
		  end do
		  !score0 = maxval(nflop) - minval(nflop)
		  score0 = sqrt(sum(nflop**2)/nproc-sum(nflop)**2/nproc**2)
		  if (score0 >= score) then
		     ie(idxmin) = ie(idxmin) - 1; is(idxmin+1) = is(idxmin+1) - 1
	             do idx = idxmin, idxmin+1
		     do idy = 1, nype
		        id = (idy-1)*nxpe + idx-1 + 1
		        nflop(id) = nflop0(id)
		     end do
		     end do
		  end if
	       end if
	       if (score0 < score) then
		  found = .True.
		  score = score0
		  exit
	       end if
	    ! Bottom
	    else if (k == 3) then
	       if (idymin == 1) then
		  score0 = score
	       else if (je(idymin-1) == js(idymin-1)) then
		  score0 = score
	       else
		  js(idymin) = js(idymin) - 1; je(idymin-1) = je(idymin-1) - 1
	          do idy = idymin-1, idymin
		  do idx = 1, nxpe
		     id = (idy-1)*nxpe + idx-1 + 1
		     imin = max(1,is(idx)-halo); imax = min(nx,ie(idx)+halo)
		     jmin = max(1,js(idy)-halo); jmax = min(ny,je(idy)+halo)
		     n = (ie(idx)-is(idx)+1)*(je(idy)-js(idy)+1)
		     m = sum(nobs(imin:imax,jmin:jmax))
		     if (m == 0) m = 1
		     nflop0(id) = nflop(id)
		     nflop(id) = log10(1.d0*n) + 2.d0*log10(1.d0*m)
		  end do
		  end do
		  !score0 = maxval(nflop) - minval(nflop)
		  score0 = sqrt(sum(nflop**2)/nproc-sum(nflop)**2/nproc**2)
		  if (score0 >= score) then
		     js(idymin) = js(idymin) + 1; je(idymin-1) = je(idymin-1) + 1
	             do idy = idymin-1, idymin
		     do idx = 1, nxpe
		        id = (idy-1)*nxpe + idx-1 + 1
		        nflop(id) = nflop0(id)
		     end do
		     end do
		  end if
	       end if
	       if (score0 < score) then
		  found = .True.
		  score = score0
		  exit
	       end if
	    ! Top
	    else if (k == 4) then
	       if (idymin == nype) then
		  score0 = score
	       else if (je(idymin+1) == js(idymin+1)) then
		  score0 = score
	       else
		  je(idymin) = je(idymin) + 1; js(idymin+1) = js(idymin+1) + 1
	          do idy = idymin, idymin+1
		  do idx = 1, nxpe
		     id = (idy-1)*nxpe + idx-1 + 1
		     imin = max(1,is(idx)-halo); imax = min(nx,ie(idx)+halo)
		     jmin = max(1,js(idy)-halo); jmax = min(ny,je(idy)+halo)
		     n = (ie(idx)-is(idx)+1)*(je(idy)-js(idy)+1)
		     m = sum(nobs(imin:imax,jmin:jmax))
		     if (m == 0) m = 1
		     nflop0(id) = nflop(id)
		     nflop(id) = log10(1.d0*n) + 2.d0*log10(1.d0*m)
		  end do
		  end do
		  !score0 = maxval(nflop) - minval(nflop)
		  score0 = sqrt(sum(nflop**2)/nproc-sum(nflop)**2/nproc**2)
		  if (score0 >= score) then
		     je(idymin) = je(idymin) - 1; js(idymin+1) = js(idymin+1) - 1
	             do idy = idymin, idymin+1
		     do idx = 1, nxpe
		        id = (idy-1)*nxpe + idx-1 + 1
		        nflop(id) = nflop0(id)
		     end do
		     end do
		  end if
	       end if
	       if (score0 < score) then
		  found = .True.
		  score = score0
		  exit
	       end if
	    end if
	 end do
	 if (.not. found) exit
      end do
      do id = 1, nproc
         idx = mod(id-1,nxpe) + 1; idy = (id-1)/nxpe + 1
	 !print*, idx, idy, ie(idx)-is(idx)+1, je(idy)-js(idy)+1, nflop(id)
      end do
      !
      return
   end subroutine load_balance
   !
   !
   !
   subroutine set_index0(nx0, nxpe, idx, halo, istart, iend, is, ie, ds, de) ! halo: radius of influence
      ! This subroutine extends the computational domain a number of points (halo) for all directions
      implicit none
      integer, intent(in) :: nx0, nxpe, idx, halo
      integer, dimension(nxpe), intent(in) :: istart, iend
      integer, intent(out) :: is, ie, ds, de
      integer :: nx, nxmod
      !
      is = istart(idx+1)
      ie = iend(idx+1)
      ! now add halo
      if (idx > 0) then
	 ds = halo
	 if ((is-ds) < 1) ds = is - 1
	 is = is - ds
      else
	 ds = 0
      end if
      if (idx < nxpe-1) then
	 de = halo
	 if ((ie+de) > nx0) de = nx0 - ie
	 ie = ie + de
      else
	 de = 0
      end if
      !
      return
   end subroutine set_index0
   !
   !
   !
   subroutine set_index1(nx0, nxpe, idx, halo, is, ie, ds, de) ! halo: radius of influence
      ! This subroutine extends the computational domain a number of points (halo) for all directions
      implicit none
      integer, intent(in) :: nx0, nxpe, idx, halo
      integer, intent(out) :: is, ie, ds, de
      integer :: nx, nxmod
      !
      nx = nx0/nxpe
      nxmod = mod(nx0,nxpe)
      if (idx < nxmod) then
	 nx = nx + 1
         is = idx*nx + 1
      else
         is = (idx-nxmod)*nx + nxmod*(nx+1) + 1
      end if
      ie = is + nx - 1
      ! now add halo
      if (idx > 0) then
	 ds = halo
	 if ((is-ds) < 1) ds = is - 1
	 is = is - ds
      else
	 ds = 0
      end if
      if (idx < nxpe-1) then
	 de = halo
	 if ((ie+de) > nx0) de = nx0 - ie
	 ie = ie + de
      else
	 de = 0
      end if
      !
      return
   end subroutine set_index1
   !
   !
   !
   subroutine set_index2(nx0, nxpe, idx, halo, is, ie, ds, de) ! halo: halo size
      ! This subroutine extends the computational domain a number of points (halo) for all directions
      implicit none
      integer, intent(in) :: nx0, nxpe, idx, halo
      integer, intent(out) :: is, ie, ds, de
      integer :: nx, nxmod
      !
      nx = (nx0-2*halo)/nxpe ! 2*halo: two boundary regions
      nxmod = mod(nx0-2*halo,nxpe)
      if (idx < nxmod) then
	 nx = nx + 1
         is = idx*nx + 1
      else
         is = (idx-nxmod)*nx + nxmod*(nx+1) + 1
      end if
      if (idx > 0) is = is + halo
      if (idx == nxpe-1) nx = nx + halo
      ie = is + nx - 1
      ! now add halo
      if (idx > 0) then
	 ds = halo
	 is = is - ds
      else
	 ds = 0
      end if
      if (idx < nxpe-1) then
	 de = halo
	 ie = ie + de
      else
	 de = 0
      end if
      !
      return
   end subroutine set_index2
   !
   !
   !
   subroutine set_index3(nx0, nxpe, idx, halo, is, ie, ds, de) ! halo: halo size
      ! The halo region is only needed for interpolation, so it's not necessary to extend the computational domain
      ! in all direction like in set_index1. We only need to extend it in the north and east boundaries.
      implicit none
      integer, intent(in) :: nx0, nxpe, idx, halo
      integer, intent(out) :: is, ie, ds, de
      integer :: nx, nxmod
      !
      nx = (nx0-halo)/nxpe ! halo: the right boundary region
      nxmod = mod(nx0-halo,nxpe)
      if (idx < nxmod) then
	 nx = nx + 1
         is = idx*nx + 1
      else
         is = (idx-nxmod)*nx + nxmod*(nx+1) + 1
      end if
      if (idx == nxpe-1) nx = nx + halo
      ie = is + nx - 1
      ! now add halo
      ds = 0
      if (idx < nxpe-1) then
	 de = halo
	 ie = ie + de
      else
	 de = 0
      end if
      !
      return
   end subroutine set_index3
   !
   !
   !
   subroutine set_index_nohalo(nx0, nxpe, idx, is, ie)
      implicit none
      integer, intent(in) :: nx0, nxpe, idx
      integer, intent(out) :: is, ie
      integer :: nx, nxmod
      !
      nx = nx0/nxpe
      nxmod = mod(nx0,nxpe)
      if (idx < nxmod) then
	 nx = nx + 1
         is = idx*nx + 1
      else
         is = (idx-nxmod)*nx + nxmod*(nx+1) + 1
      end if
      ie = is + nx - 1
      !
      return
   end subroutine set_index_nohalo
   !
   !
   !
end module NodeInfo_class
