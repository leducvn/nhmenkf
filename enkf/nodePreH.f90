module NodePreH_class
! Author: Le Duc
! Created date: 5 Dec 2015
   use variable, only : nsoil, r_size
   use NodeControl_class
   use NodeProfileControl_class
   implicit none
   !
   type NodePreH
      logical :: derivative
      integer :: nx, ny, nz, nt
      real(r_size), dimension(:), allocatable :: vdz, zrp, dvtransp, dvtransw
      real(r_size), dimension(:,:), allocatable :: zs
      real(r_size), dimension(:,:,:), allocatable :: g2, deltaz, z, g, ps
      real(r_size), dimension(:,:,:,:), allocatable :: t, p
      real(r_size), dimension(:,:,:), allocatable :: dpmsldps, dpmsldts, dtsdps, dtsdp, dtsdt
      real(r_size), dimension(:,:,:,:), allocatable :: dphdph, drhdp, drhdqv, dpwvdp, dpwvdt, dpwvdqv
      real(r_size), dimension(:,:,:,:,:), allocatable :: dphdt, dphdqv, drhdt
   end type NodePreH
   !
   interface new
      module procedure new_NodePreH
   end interface
   interface destroy
      module procedure destroy_NodePreH
   end interface
   interface set_logp
      module procedure set_logp_NodePreH
   end interface
   interface set_tqv
      module procedure set_tqv_NodePreH
   end interface
   interface set_rain
      module procedure set_rain_NodePreH
   end interface
   interface set_pmsl
      module procedure set_pmsl_NodePreH
   end interface
   interface apply_PreH
      module procedure apply_PreH_NodePreH1
      module procedure apply_PreH_NodePreH2
      module procedure apply_PreH_NodePreH3
   end interface
   interface apply_DPreH
      module procedure apply_DPreH_NodePreH
   end interface
   interface apply_DPreHT
      module procedure apply_DPreHT_NodePreH
   end interface
   !
   !
   !
contains
   !
   !
   !
   subroutine new_NodePreH(self, derivative, vdz, zrp, vctransp, dvtransp, dvtransw, zs, nx, ny, nz, nt)
      implicit none
      integer, intent(in) :: nx, ny, nz, nt
      type(NodePreH), intent(inout) :: self
      logical, intent(in) :: derivative
      real(r_size), dimension(nz), intent(in) :: vdz, zrp, vctransp, dvtransp, dvtransw
      real(r_size), dimension(nx,ny), intent(in) :: zs
      integer :: k
      !
      self%derivative = derivative
      self%nx = nx; self%ny = ny; self%nz = nz; self%nt = nt
      !
      allocate(self%vdz(nz), self%zrp(nz), self%dvtransp(nz), self%dvtransw(nz))
      allocate(self%zs(nx,ny), self%g2(nx,ny,nz), self%deltaz(nx,ny,nz), self%z(nx,ny,nz), self%g(nx,ny,nz))
      if (derivative) then
	 allocate(self%t(nx,ny,nz,nt), self%p(nx,ny,nz,nt), self%ps(nx,ny,nt))
	 allocate(self%dpmsldps(nx,ny,nt), self%dpmsldts(nx,ny,nt))
	 allocate(self%dtsdps(nx,ny,nt), self%dtsdp(nx,ny,nt), self%dtsdt(nx,ny,nt))
	 allocate(self%dphdph(nx,ny,nz,nt), self%dphdt(nx,ny,nz,2,nt), self%dphdqv(nx,ny,nz,2,nt))
	 allocate(self%drhdp(nx,ny,nz,nt), self%drhdqv(nx,ny,nz,nt), self%drhdt(nx,ny,nz,3,nt))
	 allocate(self%dpwvdp(nx,ny,nz,nt), self%dpwvdt(nx,ny,nz,nt), self%dpwvdqv(nx,ny,nz,nt))
      end if
      !
      self%vdz = vdz
      self%zrp = zrp
      self%dvtransp = dvtransp
      self%dvtransw = dvtransw
      self%zs = zs
      do k = 1, nz-1
         self%g2(:,:,k) = 1. + dvtransw(k)*zs(:,:)
      end do
      self%deltaz(:,:,1) = zrp(1)*self%g2(:,:,1)
      do k = 2, nz
	 self%deltaz(:,:,k) = (zrp(k)-zrp(k-1))*self%g2(:,:,k-1)
      end do
      do k = 1, nz
         self%g2(:,:,k) = 1. + zs*dvtransp(k)
      end do
      self%z(:,:,1) = zs
      self%g(:,:,1) = 1.0d0
      do k = 1, nz
	 self%z(:,:,k) = zrp(k) + zs*vctransp(k)
	 self%g(:,:,k) = real(k,r_size)
      end do
      !
   end subroutine new_NodePreH
   !
   !
   !
   subroutine destroy_NodePreH(self)
      implicit none
      type(NodePreH), intent(inout) :: self
      !
      if (allocated(self%vdz)) deallocate(self%vdz)
      if (allocated(self%zrp)) deallocate(self%zrp)
      if (allocated(self%dvtransp)) deallocate(self%dvtransp)
      if (allocated(self%dvtransw)) deallocate(self%dvtransw)
      if (allocated(self%zs)) deallocate(self%zs)
      if (allocated(self%g2)) deallocate(self%g2)
      if (allocated(self%deltaz)) deallocate(self%deltaz)
      if (allocated(self%z)) deallocate(self%z)
      if (allocated(self%g)) deallocate(self%g)
      !
      if (allocated(self%t)) deallocate(self%t)
      if (allocated(self%p)) deallocate(self%p)
      if (allocated(self%ps)) deallocate(self%ps)
      !
      if (allocated(self%dphdph)) deallocate(self%dphdph)
      if (allocated(self%dphdt)) deallocate(self%dphdt)
      if (allocated(self%dphdqv)) deallocate(self%dphdqv)
      !
      if (allocated(self%drhdp)) deallocate(self%drhdp)
      if (allocated(self%drhdqv)) deallocate(self%drhdqv)
      if (allocated(self%drhdt)) deallocate(self%drhdt)
      !
      if (allocated(self%dpwvdp)) deallocate(self%dpwvdp)
      if (allocated(self%dpwvdt)) deallocate(self%dpwvdt)
      if (allocated(self%dpwvdqv)) deallocate(self%dpwvdqv)
      !
      return
   end subroutine destroy_NodePreH
   !
   !
   !
   subroutine set_logp_NodePreH(self, xobs, x)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      type(NodeControl), intent(in) :: xobs
      type(NodeControl), intent(inout) :: x
      !
      integer :: nz, nt, it
      real(r_size), dimension(self%nx,self%ny,self%nz) :: logp
      !
      nz = self%nz
      nt = self%nt
      do it = 1, nt
         call get_field(xobs, 'logp',  it, logp, nz)
	 call get_field(xobs, 'logps', it, logp(:,:,1:1), 1)
	 call set_field(x, 'logp', it, logp, nz)
      end do
      !
      return
   end subroutine set_logp_NodePreH
   !
   !
   !
   subroutine set_tqv_NodePreH(self, xobs, x)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      type(NodeControl), intent(in) :: xobs
      type(NodeControl), intent(inout) :: x
      !
      integer :: nz, nt, it
      real(r_size), dimension(self%nx,self%ny,self%nz) :: tmp
      !
      nz = self%nz
      nt = self%nt
      do it = 1, nt
         call get_field(xobs, 't',  it, tmp, nz)
	 call set_field(x, 't', it, tmp, nz)
	 call get_field(xobs, 'qv',  it, tmp, nz)
	 call set_field(x, 'qv', it, tmp, nz)
      end do
      !
      return
   end subroutine set_tqv_NodePreH
   !
   !
   !
   subroutine set_rain_NodePreH(self, x, xobs)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      type(NodeControl), intent(in) :: x
      type(NodeControl), intent(inout) :: xobs
      !
      integer :: nt, it
      real(r_size), dimension(self%nx,self%ny,1) :: rain
      !
      nt = self%nt
      do it = 1, nt
         call get_field(x, 'rain', it, rain, 1)
	 call set_field(xobs, 'rain', it, rain, 1)
      end do
      !
      return
   end subroutine set_rain_NodePreH
   !
   !
   !
   subroutine set_pmsl_NodePreH(self, xobs, x)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      type(NodeControl), intent(in) :: xobs
      type(NodeControl), intent(inout) :: x
      !
      integer :: nt, it
      real(r_size), dimension(self%nx,self%ny,1) :: pmsl
      !
      nt = self%nt
      do it = 1, nt
         call get_field(xobs, 'pmsl', it, pmsl, 1)
	 call set_field(x, 'pmsl', it, pmsl, 1)
      end do
      !
      return
   end subroutine set_pmsl_NodePreH
   !
   !
   !
   subroutine apply_Hph1(self, it, ps, t, qv, ph)
   ! modified so that t, qv at 1st level are t, qv at 2nd level
      use variable, only : gvrd
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny), intent(in) :: ps
      real(r_size), dimension(self%nx,self%ny,self%nz), intent(in) :: t, qv
      real(r_size), dimension(self%nx,self%ny,self%nz), intent(out) :: ph
      integer :: k
      real(r_size), dimension(self%nx,self%ny) :: tmean, factor
      real(r_size), dimension(self%nx,self%ny,self%nz) :: tv, tratio
      !
      tratio = 1.d0+0.608d0*qv
      tv = t*tratio
      tmean = tv(:,:,2)
      factor = gvrd*self%deltaz(:,:,1)/tmean
      ph(:,:,1) = ps*exp(-factor)
      if (self%derivative) then
	 self%dphdph(:,:,1,it) = exp(-factor)
	 self%dphdt(:,:,1,2,it) = ph(:,:,1)*factor/tmean*tratio(:,:,2)
	 self%dphdqv(:,:,1,2,it) = 0.608d0*ph(:,:,1)*factor/tmean*t(:,:,2)
      end if
      do k = 2, self%nz
         if (k == 2) then
	    tmean = tv(:,:,k)
	 else
	    tmean = 0.5d0*(tv(:,:,k-1)+tv(:,:,k))
	 end if
	 factor = gvrd*self%deltaz(:,:,k)/tmean
	 ph(:,:,k) = ph(:,:,k-1)*exp(-factor)
	 if (self%derivative) then
	    self%dphdph(:,:,k,it) = exp(-factor)
	    if (k == 2) then
	       self%dphdt(:,:,k,1,it) = 0.d0
	       self%dphdqv(:,:,k,1,it) = 0.d0
	       self%dphdt(:,:,k,2,it) = ph(:,:,k)*factor/tmean*tratio(:,:,k)
	       self%dphdqv(:,:,k,2,it) = 0.608d0*ph(:,:,k)*factor/tmean*t(:,:,k)
	    else
	       self%dphdt(:,:,k,1,it) = 0.5d0*ph(:,:,k)*factor/tmean*tratio(:,:,k-1)
	       self%dphdqv(:,:,k,1,it) = 0.5d0*0.608d0*ph(:,:,k)*factor/tmean*t(:,:,k-1)
	       self%dphdt(:,:,k,2,it) = 0.5d0*ph(:,:,k)*factor/tmean*tratio(:,:,k)
	       self%dphdqv(:,:,k,2,it) = 0.5d0*0.608d0*ph(:,:,k)*factor/tmean*t(:,:,k)
	    end if
	 end if
      end do
      !
      return
   end subroutine apply_Hph1
   !
   !
   !
   subroutine apply_Hph2(self, i, j, ps, t, qv, ph)
   ! modified so that t, qv at 1st level are t, qv at 2nd level
      use variable, only : gvrd
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: i, j
      real(r_size), intent(in) :: ps
      real(r_size), dimension(self%nz), intent(in) :: t, qv
      real(r_size), dimension(self%nz), intent(out) :: ph
      integer :: k
      real(r_size) :: tmean, factor
      real(r_size), dimension(:), allocatable, save :: tv, tratio
      !
      if (.not. allocated(tv)) allocate(tv(self%nz))
      if (.not. allocated(tratio)) allocate(tratio(self%nz))
      tratio(:) = 1.d0+0.608d0*qv(:)
      tv(:) = t(:)*tratio(:)
      !
      tmean = tv(2)
      factor = gvrd*self%deltaz(i,j,1)/tmean
      ph(1) = ps*exp(-factor)
      do k = 2, self%nz
         if (k == 2) then
	    tmean = tv(k)
	 else
	    tmean = 0.5d0*(tv(k-1)+tv(k))
	 end if
	 factor = gvrd*self%deltaz(i,j,k)/tmean
	 ph(k) = ph(k-1)*exp(-factor)
      end do
      !
      return
   end subroutine apply_Hph2
   !
   !
   !
   subroutine apply_DHph(dps, dt, dqv, dphdph, dphdt, dphdqv, dph, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny), intent(in) :: dps
      real(r_size), dimension(nx,ny,nz), intent(in) :: dt, dqv
      real(r_size), dimension(nx,ny,nz), intent(in) :: dphdph
      real(r_size), dimension(nx,ny,nz,2), intent(in) :: dphdt, dphdqv
      real(r_size), dimension(nx,ny,nz), intent(out) :: dph
      integer :: k
      !
      dph(:,:,1) = dphdph(:,:,1)*dps + dphdt(:,:,1,2)*dt(:,:,2) + dphdqv(:,:,1,2)*dqv(:,:,2)
      do k = 2, nz
         dph(:,:,k) = dphdph(:,:,k)*dph(:,:,k-1) + &
	            & dphdt(:,:,k,1)*dt(:,:,k-1) + dphdt(:,:,k,2)*dt(:,:,k) + &
	            & dphdqv(:,:,k,1)*dqv(:,:,k-1) + dphdqv(:,:,k,2)*dqv(:,:,k)
      end do
      !
      return
   end subroutine apply_DHph
   !
   !
   !
   subroutine apply_DHphT(dph, dphdph, dphdt, dphdqv, dps, dt, dqv, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(inout) :: dph
      real(r_size), dimension(nx,ny,nz), intent(in) :: dphdph
      real(r_size), dimension(nx,ny,nz,2), intent(in) :: dphdt, dphdqv
      real(r_size), dimension(nx,ny), intent(inout) :: dps
      real(r_size), dimension(nx,ny,nz), intent(inout) :: dt, dqv
      integer :: k
      !
      do k = nz, 2, -1
	 dph(:,:,k-1) = dph(:,:,k-1) + dphdph(:,:,k)*dph(:,:,k)
	 dt(:,:,k-1)  = dt(:,:,k-1)  + dphdt(:,:,k,1)*dph(:,:,k)
	 dt(:,:,k)    = dt(:,:,k)    + dphdt(:,:,k,2)*dph(:,:,k)
	 dqv(:,:,k-1) = dqv(:,:,k-1) + dphdqv(:,:,k,1)*dph(:,:,k)
	 dqv(:,:,k)   = dqv(:,:,k)   + dphdqv(:,:,k,2)*dph(:,:,k)
      end do
      dps = dps + dphdph(:,:,1)*dph(:,:,1)
      dt(:,:,2) = dt(:,:,2) + dphdt(:,:,1,2)*dph(:,:,1)
      dqv(:,:,2) = dqv(:,:,2) + dphdqv(:,:,1,2)*dph(:,:,1)
      !
      return
   end subroutine apply_DHphT
   !
   !
   !
   subroutine apply_Hrh1(self, it, p, t, qv, rh)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny,self%nz), intent(in) :: p, t, qv
      real(r_size), dimension(self%nx,self%ny,self%nz), intent(out) :: rh
      !
      real(r_size), parameter :: alfa = 0.378d0, beta = 0.622d0
      integer :: i, j, k
      real(r_size) :: e0c, al, bl, e0i, ai, bi, tc, es0, es1
      real(r_size), dimension(self%nx,self%ny,self%nz) :: es, e
      real(r_size), dimension(self%nx,self%ny,self%nz,3) :: desdt
      !
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      if (self%derivative) desdt = 0.d0
      do k = 1, self%nz
         do j = 1, self%ny
            do i = 1, self%nx
               tc = t(i,j,k) - tkelvn
               if (tc >= tcw) then
                  es(i,j,k) = e0c*exp(al*tc/(bl+tc))
	          if (self%derivative) desdt(i,j,k,3) = es(i,j,k)*al*bl/(bl+tc)**2
               else if(tc <= tci) then
                  es(i,j,k) = e0i*exp(ai*tc/(bi+tc))
	          if (self%derivative) desdt(i,j,k,1) = es(i,j,k)*ai*bi/(bi+tc)**2
               else
		  es0 = e0i*exp(ai*tc/(bi+tc))
	          es1 = e0c*exp(al*tc/(bl+tc))
                  es(i,j,k) = es1*(tc-tci)/(tcw-tci) + es0*(tcw-tc)/(tcw-tci)
		  if (self%derivative) desdt(i,j,k,2) = es1*al*bl/(bl+tc)**2*(tc-tci)/(tcw-tci) + es1/(tcw-tci) + &
		                                      & es0*ai*bi/(bi+tc)**2*(tcw-tc)/(tcw-tci) - es0/(tcw-tci)
               end if
            end do
         end do
      end do
      !
      e = qv*p/(alfa*qv+beta)
      rh = e/es
      if (self%derivative) then
	 self%drhdp(:,:,:,it) = qv/(alfa*qv+beta)/es
	 self%drhdqv(:,:,:,it) = beta*p/(alfa*qv+beta)**2/es
	 do k = 1, 3
	    self%drhdt(:,:,:,k,it) = -e/es**2*desdt(:,:,:,k)
	 end do
      end if
      !
      return
   end subroutine apply_Hrh1
   !
   !
   !
   subroutine apply_Hrh2(self, i, j, p, t, qv, rh)
      use variable, only : tkelvn, tcw, e0cw, tetn1w, tetn2w, tetn3w, &
                         &         tci, e0ci, tetn1i, tetn2i, tetn3i
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: i, j
      real(r_size), dimension(self%nz), intent(in) :: p, t, qv
      real(r_size), dimension(self%nz), intent(out) :: rh
      !
      real(r_size), parameter :: alfa = 0.378d0, beta = 0.622d0
      integer :: is, ie, js, je, ii, jj, k
      real(r_size) :: e0c, al, bl, e0i, ai, bi, tc, es0, es1
      real(r_size), dimension(:), allocatable, save :: es, e
      !
      if (.not. allocated(es)) allocate(es(self%nz))
      if (.not. allocated(e)) allocate(e(self%nz))
      e0c = e0cw
      al = tetn1w
      bl = tetn2w - tetn3w
      e0i = e0ci
      ai = tetn1i
      bi = tetn2i - tetn3i
      do k = 1, self%nz
	 tc = t(k) - tkelvn
	 if (tc >= tcw) then
	    es(k) = e0c*exp(al*tc/(bl+tc))
	 else if(tc <= tci) then
	    es(k) = e0i*exp(ai*tc/(bi+tc))
	 else
	    es0 = e0i*exp(ai*tc/(bi+tc))
	    es1 = e0c*exp(al*tc/(bl+tc))
	    es(k) = es1*(tc-tci)/(tcw-tci) + es0*(tcw-tc)/(tcw-tci)
	 end if
      end do
      !
      e(:) = qv(:)*p(:)/(alfa*qv(:)+beta)
      rh(:) = e(:)/es(:)
      !
      return
   end subroutine apply_Hrh2
   !
   !
   !
   subroutine apply_DHrh(dp, dt, dqv, t, drhdp, drhdqv, drhdt, drh, nx, ny ,nz)
      use variable, only : tkelvn, tcw, tci
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: dp, dt, dqv
      real(r_size), dimension(nx,ny,nz), intent(in) :: t, drhdp, drhdqv
      real(r_size), dimension(nx,ny,nz,3), intent(in) :: drhdt
      real(r_size), dimension(nx,ny,nz), intent(out) :: drh
      !
      integer :: i, j, k
      real(r_size) :: tc
      !
      drh = drhdp*dp + drhdqv*dqv
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tc = t(i,j,k) - tkelvn
               if (tc >= tcw) then
                  drh(i,j,k) = drh(i,j,k) + drhdt(i,j,k,3)*dt(i,j,k)
               else if(tc <= tci) then
                  drh(i,j,k) = drh(i,j,k) + drhdt(i,j,k,1)*dt(i,j,k)
               else
		  drh(i,j,k) = drh(i,j,k) + drhdt(i,j,k,2)*dt(i,j,k)
               end if
            end do
         end do
      end do
      !
      return
   end subroutine apply_DHrh
   !
   !
   !
   subroutine apply_DHrhT(drh, t, drhdp, drhdqv, drhdt, dp, dt, dqv, nx, ny ,nz)
      use variable, only : tkelvn, tcw, tci
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: drh
      real(r_size), dimension(nx,ny,nz), intent(in) :: t, drhdp, drhdqv
      real(r_size), dimension(nx,ny,nz,3), intent(in) :: drhdt
      real(r_size), dimension(nx,ny,nz), intent(inout) :: dp, dt, dqv
      !
      integer :: i, j, k
      real(r_size) :: tc
      !
      dp = dp + drhdp*drh
      dqv = dqv + drhdqv*drh
      do k = 1, nz
         do j = 1, ny
            do i = 1, nx
               tc = t(i,j,k) - tkelvn
               if (tc >= tcw) then
                  dt(i,j,k) = dt(i,j,k) + drhdt(i,j,k,3)*drh(i,j,k)
               else if(tc <= tci) then
                  dt(i,j,k) = dt(i,j,k) + drhdt(i,j,k,1)*drh(i,j,k)
               else
		  dt(i,j,k) = dt(i,j,k) + drhdt(i,j,k,2)*drh(i,j,k)
               end if
            end do
         end do
      end do
      !
      return
   end subroutine apply_DHrhT
   !
   !
   !
   subroutine apply_Hpwv1(self, it, p, t, qv, pwv)
      use variable, only : rd
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny,self%nz), intent(in) :: p, t, qv
      real(r_size), dimension(self%nx,self%ny), intent(out) :: pwv
      !
      integer :: k
      real(r_size), dimension(self%nx,self%ny) :: tmp
      real(r_size), dimension(self%nx,self%ny,self%nz) :: rhog2
      !
      rhog2 = p/(rd*t)*self%g2
      pwv = 0.d0
      do k = 2, self%nz-1
	 pwv = pwv + qv(:,:,k)*rhog2(:,:,k)*self%vdz(k)
      end do
      !
      if (self%derivative) then
	 self%dpwvdp(:,:,:,it) = 0.d0
	 self%dpwvdt(:,:,:,it) = 0.d0
	 self%dpwvdqv(:,:,:,it) = 0.d0
	 do k = 2, self%nz-1
	    tmp = self%g2(:,:,k)*self%vdz(k)/rd
	    self%dpwvdp(:,:,k,it) = tmp*qv(:,:,k)/t(:,:,k)
	    self%dpwvdt(:,:,k,it) = -tmp*p(:,:,k)*qv(:,:,k)/t(:,:,k)**2
	    self%dpwvdqv(:,:,k,it) = tmp*p(:,:,k)/t(:,:,k)
	 end do
      end if
      !
      return
   end subroutine apply_Hpwv1
   !
   !
   !
   subroutine apply_Hpwv2(self, i, j, p, t, qv, pwv)
      use variable, only : rd
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: i, j
      real(r_size), dimension(self%nz), intent(in) :: p, t, qv
      real(r_size), intent(out) :: pwv
      !
      integer :: k
      real(r_size), dimension(:), allocatable, save :: rhog2
      !
      if (.not. allocated(rhog2)) allocate(rhog2(self%nz))
      rhog2(:) = p(:)/(rd*t(:))*self%g2(i,j,:)
      pwv = 0.d0
      do k = 2, self%nz-1
	 pwv = pwv + qv(k)*rhog2(k)*self%vdz(k)
      end do
      !
      return
   end subroutine apply_Hpwv2
   !
   !
   !
   subroutine apply_DHpwv(dp, dt, dqv, dpwvdp, dpwvdt, dpwvdqv, dpwv, nx, ny, nz)
      use variable, only : rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny,nz), intent(in) :: dp, dt, dqv
      real(r_size), dimension(nx,ny,nz), intent(in) :: dpwvdp, dpwvdt, dpwvdqv
      real(r_size), dimension(nx,ny), intent(out) :: dpwv
      !
      integer :: k
      !
      dpwv = 0.d0
      do k = 2, nz-1
	 dpwv = dpwv + dpwvdp(:,:,k)*dp(:,:,k) + dpwvdt(:,:,k)*dt(:,:,k) + dpwvdqv(:,:,k)*dqv(:,:,k)
      end do
      !
      return
   end subroutine apply_DHpwv
   !
   !
   !
   subroutine apply_DHpwvT(dpwv, dpwvdp, dpwvdt, dpwvdqv, dp, dt, dqv, nx, ny, nz)
      use variable, only : rd
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny), intent(in) :: dpwv
      real(r_size), dimension(nx,ny,nz), intent(in) :: dpwvdp, dpwvdt, dpwvdqv
      real(r_size), dimension(nx,ny,nz), intent(inout) :: dp, dt, dqv
      !
      integer :: k
      !
      do k = nz-1, 2, -1
	 dp(:,:,k) = dp(:,:,k) + dpwvdp(:,:,k)*dpwv
	 dt(:,:,k) = dt(:,:,k) + dpwvdt(:,:,k)*dpwv
	 dqv(:,:,k) = dqv(:,:,k) + dpwvdqv(:,:,k)*dpwv
      end do
      !
      return
   end subroutine apply_DHpwvT
   !
   !
   !
   subroutine apply_Hpmsl1(self, it, ps, p, t, pmsl)
      use variable, only : gvrd, gamma, gmrgh, gmrgi
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: it
      real(r_size), dimension(self%nx,self%ny), intent(in) :: ps
      real(r_size), dimension(self%nx,self%ny,self%nz), intent(in) :: p, t
      real(r_size), dimension(self%nx,self%ny), intent(out) :: pmsl
      real(r_size), dimension(self%nx,self%ny) :: ts
      !
      ts = t(:,:,2)*(ps/p(:,:,2))**gmrgh
      pmsl = ps*(1.+gamma*self%zs/ts)**gmrgi
      if (self%derivative) then
	 self%dtsdps(:,:,it) = gmrgh*ts/ps
	 self%dtsdp(:,:,it) = -gmrgh*ts/p(:,:,2)
	 self%dtsdt(:,:,it) = ts/t(:,:,2)
	 self%dpmsldps(:,:,it) = pmsl/ps
	 self%dpmsldts(:,:,it) = -gvrd*pmsl*self%zs/(ts**2+gamma*self%zs*ts)
      end if
      !
      return
   end subroutine apply_Hpmsl1
   !
   !
   !
   subroutine apply_Hpmsl2(self, i, j, ps, p, t, pmsl)
      use variable, only : gvrd, gamma, gmrgh, gmrgi
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: i, j
      real(r_size), intent(in) :: ps
      real(r_size), dimension(self%nz), intent(in) :: p, t
      real(r_size), intent(out) :: pmsl
      real(r_size) :: ts
      !
      ts = t(2)*(ps/p(2))**gmrgh
      pmsl = ps*(1.+gamma*self%zs(i,j)/ts)**gmrgi
      !
      return
   end subroutine apply_Hpmsl2
   !
   !
   !
   subroutine apply_DHpmsl(dps, dp, dt, dpmsldps, dpmsldts, dtsdps, dtsdp, dtsdt, dpmsl, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny), intent(in) :: dps
      real(r_size), dimension(nx,ny,nz), intent(in) :: dp, dt
      real(r_size), dimension(nx,ny), intent(in) :: dpmsldps, dpmsldts, dtsdps, dtsdp, dtsdt
      real(r_size), dimension(nx,ny), intent(out) :: dpmsl
      real(r_size), dimension(nx,ny) :: dts
      !
      dts = dtsdps*dps + dtsdp*dp(:,:,2) + dtsdt*dt(:,:,2)
      dpmsl = dpmsldps*dps + dpmsldts*dts
      !
      return
   end subroutine apply_DHpmsl
   !
   !
   !
   subroutine apply_DHpmslT(dpmsl, dpmsldps, dpmsldts, dtsdps, dtsdp, dtsdt, dps, dp, dt, nx, ny, nz)
      implicit none
      !
      integer, intent(in) :: nx, ny, nz
      real(r_size), dimension(nx,ny), intent(in) :: dpmsl
      real(r_size), dimension(nx,ny), intent(in) :: dpmsldps, dpmsldts, dtsdps, dtsdp, dtsdt
      real(r_size), dimension(nx,ny), intent(inout) :: dps
      real(r_size), dimension(nx,ny,nz), intent(inout) :: dp, dt
      real(r_size), dimension(nx,ny) :: dts
      !
      dps = dps + dpmsldps*dpmsl
      dts = 0.
      dts = dts + dpmsldts*dpmsl
      dps = dps + dtsdps*dts
      dp(:,:,2) = dp(:,:,2) + dtsdp*dts
      dt(:,:,2) = dt(:,:,2) + dtsdt*dts
      !
      return
   end subroutine apply_DHpmslT
   !
   !
   !
   subroutine apply_PreH_NodePreH1(self, control_mode, xbck, x, xobs)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: control_mode
      type(NodeControl), intent(in) :: xbck
      type(NodeControl), intent(in) :: x
      type(NodeControl), intent(inout) :: xobs
      !
      integer :: nx, ny, nz, nt, it
      real(r_size), dimension(self%nx,self%ny) :: ps, dps, pwv, pmsl
      real(r_size), dimension(self%nx,self%ny,1) :: tmp
      real(r_size), dimension(self%nx,self%ny,self%nz) :: u, v, t, p, qv, pnh, ph, rh
      real(r_size), dimension(self%nx,self%ny,self%nz) :: du, dv, dt, dp, dqv
      !
      nx = self%nx
      ny = self%ny
      nz = self%nz
      nt = self%nt
      pnh(:,:,1) = 0.d0
      do it = 1, nt
         call get_field(xbck, 'u', it, u, nz)
	 call get_field(xbck, 'v', it, v, nz)
	 call get_field(xbck, 't', it, t, nz)
	 call get_field(xbck, 'p', it, p, nz)
	 call get_field(xbck, 'qv', it, qv, nz)
	 ps = p(:,:,1)
	 pnh(:,:,2:nz) = p(:,:,2:nz)
	 !
	 call get_field(x, 'u', it, du, nz)
	 call get_field(x, 'v', it, dv, nz)
	 call get_field(x, 't', it, dt, nz)
	 call get_field(x, 'p', it, dp, nz)
	 call get_field(x, 'qv', it, dqv, nz)
	 dps = dp(:,:,1)
	 u = u + du
	 v = v + dv
	 t = t + dt
	 pnh(:,:,2:nz) = pnh(:,:,2:nz) + dp(:,:,2:nz)
	 qv = qv + dqv
	 ps = ps + dps
	 !
	 call apply_Hph1(self, it, ps, t, qv, ph)
	 if (control_mode < 10) then
	    p(:,:,:) = ph(:,:,:)
	 else
	    p = ph + pnh
	 end if
	 call apply_Hrh1(self, it, p, t, qv, rh)
	 call apply_Hpwv1(self, it, p, t, qv, pwv)
	 call apply_Hpmsl1(self, it, ps, p, t, pmsl)
	 !
	 call set_field(xobs, 'u',  it, u, nz)
	 call set_field(xobs, 'v',  it, v, nz)
	 call set_field(xobs, 't',  it, t, nz)
	 call set_field(xobs, 'p',  it, p, nz)
	 call set_field(xobs, 'logp',  it, log(p), nz)
	 call set_field(xobs, 'qv', it, qv, nz)
	 call set_field(xobs, 'rh', it, rh, nz)
	 tmp(:,:,1) = pwv
	 call set_field(xobs, 'pwv', it, tmp, 1)
	 tmp(:,:,1) = ps
	 call set_field(xobs, 'ps', it, tmp, 1)
	 call set_field(xobs, 'logps', it, log(tmp), 1)
	 tmp(:,:,1) = pmsl
	 call set_field(xobs, 'pmsl', it, tmp, 1)
	 if (self%derivative) then
	    self%t(:,:,:,it) = t
	    self%p(:,:,:,it) = p
	    self%ps(:,:,it)  = ps
	 end if
      end do
      call set_field(xobs, 'z', 1, self%z, nz)
      call set_field(xobs, 'g', 1, self%g, nz)
      !
      return
   end subroutine apply_PreH_NodePreH1
   !
   !
   !
   subroutine apply_PreH_NodePreH2(self, control_mode, processed, k2ijt, xbck, x, xobs, nxyt)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: control_mode, nxyt
      logical, dimension(self%nx,self%ny), intent(in) :: processed
      integer, dimension(nxyt,3), intent(in) :: k2ijt
      type(NodeProfileControl), intent(in) :: xbck, x
      type(NodeProfileControl), intent(inout) :: xobs
      !
      logical :: update
      integer :: nz, ixyt, i, j, it, imin, jmin, ii, jj
      real(r_size) :: ps, dps, pwv, pmsl
      real(r_size), dimension(:), allocatable, save :: tmp
      real(r_size), dimension(:), allocatable, save :: u, v, t, p, qv, pnh, ph, rh
      real(r_size), dimension(:), allocatable, save :: du, dv, dt, dp, dqv
      !
      if (.not. allocated(tmp)) allocate(tmp(1))
      if (.not. allocated(u)) allocate(u(self%nz))
      if (.not. allocated(v)) allocate(v(self%nz))
      if (.not. allocated(t)) allocate(t(self%nz))
      if (.not. allocated(p)) allocate(p(self%nz))
      if (.not. allocated(qv)) allocate(qv(self%nz))
      if (.not. allocated(pnh)) allocate(pnh(self%nz))
      if (.not. allocated(ph)) allocate(ph(self%nz))
      if (.not. allocated(rh)) allocate(rh(self%nz))
      if (.not. allocated(du)) allocate(du(self%nz))
      if (.not. allocated(dv)) allocate(dv(self%nz))
      if (.not. allocated(dt)) allocate(dt(self%nz))
      if (.not. allocated(dp)) allocate(dp(self%nz))
      if (.not. allocated(dqv)) allocate(dqv(self%nz))
      !
      nz = self%nz
      do ixyt = 1, nxyt
         i = k2ijt(ixyt,1); j = k2ijt(ixyt,2); it = k2ijt(ixyt,3)
	 update = .False.
	 imin = max(1,i-1); jmin = max(1,j-1)
	 do ii = imin, i
	 do jj = jmin, j
	    if (.not. processed(ii,jj)) update = .True.
	 end do
	 end do
	 if (.not. update) cycle
	 !
	 pnh(1) = 0.d0
	 call get_field(xbck, 'u', ixyt, u, nz)
	 call get_field(xbck, 'v', ixyt, v, nz)
	 call get_field(xbck, 't', ixyt, t, nz)
	 call get_field(xbck, 'p', ixyt, p, nz)
	 call get_field(xbck, 'qv', ixyt, qv, nz)
	 ps = p(1)
	 pnh(2:nz) = p(2:nz)
	 !
	 call get_field(x, 'u', ixyt, du, nz)
	 call get_field(x, 'v', ixyt, dv, nz)
	 call get_field(x, 't', ixyt, dt, nz)
	 call get_field(x, 'p', ixyt, dp, nz)
	 call get_field(x, 'qv', ixyt, dqv, nz)
	 dps = dp(1)
	 u(:) = u(:) + du(:)
	 v(:) = v(:) + dv(:)
	 t(:) = t(:) + dt(:)
	 pnh(2:nz) = pnh(2:nz) + dp(2:nz)
	 qv(:) = qv(:) + dqv(:)
	 ps = ps + dps
	 !
	 call apply_Hph2(self, i, j, ps, t, qv, ph)
	 if (control_mode < 10) then
	    p(:) = ph(:)
	 else
	    p(:) = ph(:) + pnh(:)
	 end if
	 call apply_Hrh2(self, i, j, p, t, qv, rh)
	 call apply_Hpwv2(self, i, j, p, t, qv, pwv)
	 call apply_Hpmsl2(self, i, j, ps, p, t, pmsl)
	 !
	 call set_field(xobs, 'u', ixyt, u, nz)
	 call set_field(xobs, 'v', ixyt, v, nz)
	 call set_field(xobs, 't', ixyt, t, nz)
	 call set_field(xobs, 'p', ixyt, p, nz)
	 call set_field(xobs, 'logp', ixyt, log(p(:)), nz)
	 call set_field(xobs, 'qv', ixyt, qv, nz)
	 call set_field(xobs, 'rh', ixyt, rh, nz)
	 tmp(1) = pwv
	 call set_field(xobs, 'pwv', ixyt, tmp, 1)
	 tmp(1) = ps
	 call set_field(xobs, 'ps', ixyt, tmp, 1)
	 tmp(1) = log(ps)
	 call set_field(xobs, 'logps', ixyt, tmp, 1)
	 tmp(1) = pmsl
	 call set_field(xobs, 'pmsl', ixyt, tmp, 1)
	 call set_field(xobs, 'z', ixyt, self%z(i,j,:), nz)
	 call set_field(xobs, 'g', ixyt, self%g(i,j,:), nz)
      end do
      !
      return
   end subroutine apply_PreH_NodePreH2
   !
   !
   !
   subroutine apply_PreH_NodePreH3(self, control_mode, ip, jp, ijt2k, xbck, x, xobs)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: control_mode, ip, jp
      integer, dimension(2,2,self%nt), intent(in) :: ijt2k
      type(NodeProfileControl), intent(in) :: xbck, x
      type(NodeProfileControl), intent(inout) :: xobs
      !
      integer :: nz, nt, ixyt, i, j, it
      real(r_size) :: ps, dps, pwv, pmsl
      real(r_size), dimension(:), allocatable, save :: tmp
      real(r_size), dimension(:), allocatable, save :: u, v, t, p, qv, pnh, ph, rh
      real(r_size), dimension(:), allocatable, save :: du, dv, dt, dp, dqv
      !
      if (.not. allocated(tmp)) allocate(tmp(1))
      if (.not. allocated(u)) allocate(u(self%nz))
      if (.not. allocated(v)) allocate(v(self%nz))
      if (.not. allocated(t)) allocate(t(self%nz))
      if (.not. allocated(p)) allocate(p(self%nz))
      if (.not. allocated(qv)) allocate(qv(self%nz))
      if (.not. allocated(pnh)) allocate(pnh(self%nz))
      if (.not. allocated(ph)) allocate(ph(self%nz))
      if (.not. allocated(rh)) allocate(rh(self%nz))
      if (.not. allocated(du)) allocate(du(self%nz))
      if (.not. allocated(dv)) allocate(dv(self%nz))
      if (.not. allocated(dt)) allocate(dt(self%nz))
      if (.not. allocated(dp)) allocate(dp(self%nz))
      if (.not. allocated(dqv)) allocate(dqv(self%nz))
      !
      nz = self%nz; nt = self%nt
      do i = ip, ip+1
      do j = jp, jp+1
         do it = 1, nt
	    ixyt = ijt2k(i-ip+1,j-jp+1,it)
	    if (ixyt == 0) cycle
	    !
	    pnh(1) = 0.d0
	    call get_field(xbck, 'u', ixyt, u, nz)
	    call get_field(xbck, 'v', ixyt, v, nz)
	    call get_field(xbck, 't', ixyt, t, nz)
	    call get_field(xbck, 'p', ixyt, p, nz)
	    call get_field(xbck, 'qv', ixyt, qv, nz)
	    ps = p(1)
	    pnh(2:nz) = p(2:nz)
	    !
	    call get_field(x, 'u', ixyt, du, nz)
	    call get_field(x, 'v', ixyt, dv, nz)
	    call get_field(x, 't', ixyt, dt, nz)
	    call get_field(x, 'p', ixyt, dp, nz)
	    call get_field(x, 'qv', ixyt, dqv, nz)
	    dps = dp(1)
	    u(:) = u(:) + du(:)
	    v(:) = v(:) + dv(:)
	    t(:) = t(:) + dt(:)
	    pnh(2:nz) = pnh(2:nz) + dp(2:nz)
	    qv(:) = qv(:) + dqv(:)
	    ps = ps + dps
	    !
	    call apply_Hph2(self, i, j, ps, t, qv, ph)
	    if (control_mode < 10) then
	       p(:) = ph(:)
	    else
	       p(:) = ph(:) + pnh(:)
	    end if
	    call apply_Hrh2(self, i, j, p, t, qv, rh)
	    call apply_Hpwv2(self, i, j, p, t, qv, pwv)
	    call apply_Hpmsl2(self, i, j, ps, p, t, pmsl)
	    !
	    call set_field(xobs, 'u', ixyt, u, nz)
	    call set_field(xobs, 'v', ixyt, v, nz)
	    call set_field(xobs, 't', ixyt, t, nz)
	    call set_field(xobs, 'p', ixyt, p, nz)
	    call set_field(xobs, 'logp', ixyt, log(p(:)), nz)
	    call set_field(xobs, 'qv', ixyt, qv, nz)
	    call set_field(xobs, 'rh', ixyt, rh, nz)
	    tmp(1) = pwv
	    call set_field(xobs, 'pwv', ixyt, tmp, 1)
	    tmp(1) = ps
	    call set_field(xobs, 'ps', ixyt, tmp, 1)
	    tmp(1) = log(ps)
	    call set_field(xobs, 'logps', ixyt, tmp, 1)
	    tmp(1) = pmsl
	    call set_field(xobs, 'pmsl', ixyt, tmp, 1)
	    call set_field(xobs, 'z', ixyt, self%z(i,j,:), nz)
	    call set_field(xobs, 'g', ixyt, self%g(i,j,:), nz)
         end do
      end do
      end do
      !
      return
   end subroutine apply_PreH_NodePreH3
   !
   !
   !
   subroutine apply_DPreH_NodePreH(self, control_mode, x, xobs)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: control_mode
      type(NodeControl), intent(in) :: x
      type(NodeControl), intent(inout) :: xobs
      !
      integer :: nx, ny, nz, nt, it
      real(r_size), dimension(self%nx,self%ny) :: dps, dlogps, dpwv, dpmsl
      real(r_size), dimension(self%nx,self%ny,1) :: tmp
      real(r_size), dimension(self%nx,self%ny,self%nz) :: du, dv, dt, dp, dqv, dpnh, dph, dlogp, drh
      !
      nx = self%nx
      ny = self%ny
      nz = self%nz
      nt = self%nt
      dpnh(:,:,1) = 0.d0
      do it = 1, nt
	 call get_field(x, 'u', it, du, nz)
	 call get_field(x, 'v', it, dv, nz)
	 call get_field(x, 't', it, dt, nz)
	 call get_field(x, 'p', it, dp, nz)
	 call get_field(x, 'qv', it, dqv, nz)
	 dps = dp(:,:,1)
	 dpnh(:,:,2:nz) = dp(:,:,2:nz)
	 !
	 call apply_DHph(dps, dt, dqv, self%dphdph(:,:,:,it), self%dphdt(:,:,:,:,it), self%dphdqv(:,:,:,:,it), dph, nx, ny, nz)
	 if (control_mode < 10) then
	    dp(:,:,:) = dph(:,:,:)
	 else
	    dp = dph + dpnh
	 end if
	 call apply_DHrh(dp, dt, dqv, self%t(:,:,:,it), &
	               & self%drhdp(:,:,:,it), self%drhdqv(:,:,:,it), self%drhdt(:,:,:,:,it), drh, nx, ny, nz)
	 call apply_DHpwv(dp, dt, dqv, self%dpwvdp(:,:,:,it), self%dpwvdt(:,:,:,it), self%dpwvdqv(:,:,:,it), dpwv, nx, ny, nz)
	 call apply_DHpmsl(dps, dp, dt, self%dpmsldps(:,:,it), self%dpmsldts(:,:,it), &
	                 & self%dtsdps(:,:,it), self%dtsdp(:,:,it), self%dtsdt(:,:,it), dpmsl, nx, ny, nz)
	 dlogp = dp/self%p(:,:,:,it)
	 dlogps = dps/self%ps(:,:,it)
	 !
	 call set_field(xobs, 'u',  it, du, nz)
	 call set_field(xobs, 'v',  it, dv, nz)
	 call set_field(xobs, 't',  it, dt, nz)
	 call set_field(xobs, 'p',  it, dp, nz)
	 call set_field(xobs, 'logp',  it, dlogp, nz)
	 call set_field(xobs, 'qv', it, dqv, nz)
	 call set_field(xobs, 'rh', it, drh, nz)
	 tmp(:,:,1) = dpwv
	 call set_field(xobs, 'pwv', it, tmp, 1)
	 tmp(:,:,1) = dps
	 call set_field(xobs, 'ps', it, tmp, 1)
	 tmp(:,:,1) = dlogps
	 call set_field(xobs, 'logps', it, tmp, 1)
	 tmp(:,:,1) = dpmsl
	 call set_field(xobs, 'pmsl', it, tmp, 1)
	 
      end do
      dph(:,:,:) = 0.d0
      call set_field(xobs, 'z', 1, dph, nz)
      call set_field(xobs, 'g', 1, dph, nz)
      !
      return
   end subroutine apply_DPreH_NodePreH
   !
   !
   !
   subroutine apply_DPreHT_NodePreH(self, control_mode, xobs, x)
      implicit none
      !
      type(NodePreH), intent(inout) :: self
      integer, intent(in) :: control_mode
      type(NodeControl), intent(in) :: xobs
      type(NodeControl), intent(inout) :: x
      !
      integer :: nx, ny, nz, nt, it
      real(r_size), dimension(self%nx,self%ny) :: dps, dlogps, dpwv, dpmsl
      real(r_size), dimension(self%nx,self%ny,1) :: tmp
      real(r_size), dimension(self%nx,self%ny,4) :: dtsoil
      real(r_size), dimension(self%nx,self%ny,self%nz) :: du, dv, dt, dp, dqv, dpnh, dph, dlogp, drh
      !
      nx = self%nx
      ny = self%ny
      nz = self%nz
      nt = self%nt
      dtsoil = 0.d0
      do it = 1, nt
         call get_field(xobs, 'u',  it, du, nz)
	 call get_field(xobs, 'v',  it, dv, nz)
	 call get_field(xobs, 't',  it, dt, nz)
	 call get_field(xobs, 'p',  it, dp, nz)
	 call get_field(xobs, 'logp',  it, dlogp, nz)
	 call get_field(xobs, 'qv', it, dqv, nz)
	 call get_field(xobs, 'rh', it, drh, nz)
	 call get_field(xobs, 'pwv', it, tmp, 1)
	 dpwv = tmp(:,:,1)
	 call get_field(xobs, 'ps', it, tmp, 1)
	 dps = tmp(:,:,1)
	 call get_field(xobs, 'logps', it, tmp, 1)
	 dlogps = tmp(:,:,1)
	 call get_field(xobs, 'pmsl', it, tmp, 1)
	 dpmsl = tmp(:,:,1)
	 !
	 dps = dps + dlogps/self%ps(:,:,it)
	 dp = dp + dlogp/self%p(:,:,:,it)
	 call apply_DHpmslT(dpmsl, self%dpmsldps(:,:,it), self%dpmsldts(:,:,it), &
	                  & self%dtsdps(:,:,it), self%dtsdp(:,:,it), self%dtsdt(:,:,it), dps, dp, dt, nx, ny, nz)
	 call apply_DHpwvT(dpwv, self%dpwvdp(:,:,:,it), self%dpwvdt(:,:,:,it), self%dpwvdqv(:,:,:,it), dp, dt, dqv, nx, ny, nz)
	 call apply_DHrhT(drh, self%t(:,:,:,it), self%drhdp(:,:,:,it), self%drhdqv(:,:,:,it), self%drhdt(:,:,:,:,it), &
	                & dp, dt, dqv, nx, ny, nz)
	 dph = 0.d0
	 dpnh = 0.d0
	 if (control_mode < 10) then
	    dph = dph + dp
	 else
	    dph = dph + dp
	    dpnh = dpnh + dp
	 end if
	 call apply_DHphT(dph, self%dphdph(:,:,:,it), self%dphdt(:,:,:,:,it), self%dphdqv(:,:,:,:,it), dps, dt, dqv, nx, ny, nz)
	 !
	 dp(:,:,1) = dps(:,:)
	 dp(:,:,2:nz) = dpnh(:,:,2:nz)
	 call set_field(x, 'u', it, du, nz)
	 call set_field(x, 'v', it, dv, nz)
	 call set_field(x, 't', it, dt, nz)
	 call set_field(x, 'p', it, dp, nz)
	 call set_field(x, 'qv', it, dqv, nz)
	 call set_field(x, 'tsoil', it, dtsoil, nsoil)
      end do
      !
      return
   end subroutine apply_DPreHT_NodePreH
   !
   !
   !
end module NodePreH_class
