program outer2inner
   use variable
   use outer2innerlib
   use mpi
   implicit none
   !
   character(8) :: infile='xdriving', outfile='xnesting'
   integer :: nproc, myid, ierror
   integer :: k, idate(5), hsdiff = 0
   real :: dx0, dy0, pmsl_mean, ptop_mean, ps_mean
   !
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)
   !
   ! read namelist
   open(10, file='namelist.txt')
   read(10, control)
   read(10, model)
   read(10, driving)
   close(10)
   resolution = 1000.*resolution
   resolution0 = 1000.*resolution0
   dx0 = resolution0/resolution
   dy0 = dx0
   idate(1) = iyear
   idate(2) = imonth
   idate(3) = iday
   idate(4) = ihour
   idate(5) = iminute
   !
   ! model
   allocate(vdz(nz), vrdz(nz), vrdz2(nz))
   allocate(zrp(nz), vctransp(nz), dvtransp(nz))
   allocate(zrw(nz), vctransw(nz), dvtransw(nz))
   allocate(zs(nx,ny), landsea(nx,ny), fmap(nx,ny,6), lon(nx,ny), lat(nx,ny))
   ! horizontal
   call read_nhmcst(30, nx, ny, lon, lat, zs, landsea)
   call cfmap(proj, slat, lat, nx, ny, fmap)
   ! vertical
   call vrgdis(nz, nz, 1, nz, nz, 1, dz2, dz1, dz2, vdz, vrdz, vrdz2)
   call setzrp(vrdz2, nz, zrp)
   call calc_zcoordinate(vctrans_type, n_vctrans, ztop, zl_vctrans, zh_vctrans, zrp, nz, vctransp, dvtransp)
   call setzrw(vrdz, nz, zrw)
   call calc_zcoordinate(vctrans_type, n_vctrans, zrw(nz-1), zl_vctrans, zh_vctrans, zrw, nz, vctransw, dvtransw)
   deallocate(vdz, vrdz, vrdz2)
   deallocate(lon, lat)
   !
   ! driving
   allocate(vdz0(nz0), vrdz0(nz0), vrdz20(nz0))
   allocate(zrp0(nz0), vctransp0(nz0), dvtransp0(nz0))
   allocate(zrw0(nz0), vctransw0(nz0), dvtransw0(nz0))
   allocate(zs0(nx0,ny0), landsea0(nx0,ny0))
   allocate(lon0(nx0,ny0), lat0(nx0,ny0), lonu0(nx0,ny0), latu0(nx0,ny0), lonv0(nx0,ny0), latv0(nx0,ny0))
   ! horizontal
   call read_nhmcst(40, nx0, ny0, lon0, lat0, zs0, landsea0)
   call setorg(proj0, resolution0, slat0, slon0, xi0, xj0, xlat0, xlon0, lat0, lon0, nx0, ny0, latu0, lonu0, latv0, lonv0)
   ! vertical
   call vrgdis(nz0, nz0, 1, nz0, nz0, 1, dz20, dz10, dz20, vdz0, vrdz0, vrdz20)
   call setzrp(vrdz20, nz0, zrp0)
   call calc_zcoordinate(vctrans_type0, n_vctrans0, ztop0, zl_vctrans0, zh_vctrans0, zrp0, nz0, vctransp0, dvtransp0)
   call setzrw(vrdz0, nz0, zrw0)
   call calc_zcoordinate(vctrans_type0, n_vctrans0, zrw0(nz0-1), zl_vctrans0, zh_vctrans0, zrw0, nz0, vctransw0, dvtransw0)
   deallocate(vdz0, vrdz0, vrdz20)
   deallocate(landsea0)
   !
   ! coordinate
   allocate(fi(nx0,ny0), fj(nx0,ny0), fiu(nx0,ny0), fju(nx0,ny0), fiv(nx0,ny0), fjv(nx0,ny0))
   call setfifj(proj, nx, ny, resolution, resolution, slat, slon, xi, xj, xlat, xlon, &
              & latu0, lonu0, latv0, lonv0, lat0, lon0, nx0, ny0, fi, fj, fiu, fju, fiv, fjv)
   deallocate(lon0, lat0, lonu0, latu0, lonv0, latv0)
   !
   ! allocate
   ! large domain
   allocate(ps0(nx0,ny0), pb0(nx0,ny0), pmsl0(nx0,ny0), ptop0(nx0,ny0), ts0(nx0,ny0))
   allocate(pts0(nx0,ny0), us0(nx0,ny0), vs0(nx0,ny0), qvs0(nx0,ny0))
   allocate(tsoil0(nx0,ny0,ngm))
   allocate(u0(nx0,ny0,nz0), v0(nx0,ny0,nz0), w0(nx0,ny0,nz0), pt0(nx0,ny0,nz0), p0(nx0,ny0,nz0))
   allocate(qv0(nx0,ny0,nz0), qc0(nx0,ny0,nz0), qi0(nx0,ny0,nz0), qr0(nx0,ny0,nz0), qs0(nx0,ny0,nz0), qg0(nx0,ny0,nz0))
   allocate(h0(nx0,ny0,nz0), t0(nx0,ny0,nz0), pnh0(nx0,ny0,nz0))
   ! small domain
   allocate(deltazs(nx,ny), sst(nx,ny), dptsdt(nx,ny))
   allocate(ps(nx,ny), ts(nx,ny), ptop(nx,ny), pmsl(nx,ny), pts(nx,ny))
   allocate(us(nx,ny), vs(nx,ny), qvs(nx,ny))
   allocate(tsoil(nx,ny,ngm), p(nx,ny,nz), t(nx,ny,nz), rh(nx,ny,nz))
   allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), pt(nx,ny,nz), pnh(nx,ny,nz))
   allocate(qv(nx,ny,nz), qc(nx,ny,nz), qi(nx,ny,nz), qr(nx,ny,nz), qs(nx,ny,nz), qg(nx,ny,nz))
   !
   ! const fields
   call average2D(dx0, dy0, fi, fj, zs0, nx0, ny0, nx, ny, deltazs)
   deltazs = deltazs - zs
   call read_sst(20, nx, ny, sst)
   dptsdt = 0.
   !
   ! Read input fields
   if (input_mode == nhm) then
      call read_nhm(infile, nx0, ny0, nz0, u0, v0, w0, t0, pnh0, qv0, qc0, qi0, qr0, qs0, qg0, tsoil0)
      ps0(:,:) = pnh0(:,:,1)
      pnh0(:,:,1) = 0.d0
      call average3Duv(hsdiff, dx0, dy0, 0, 0, zrp0, vctransp0, dvtransp0, zrp, vctransp, dvtransp, zs0, zs, fiu, fju, u0, nx0, ny0, nz0, nx, ny, nz, u)
      call average3Duv(hsdiff, dx0, dy0, 0, 0, zrp0, vctransp0, dvtransp0, zrp, vctransp, dvtransp, zs0, zs, fiv, fjv, v0, nx0, ny0, nz0, nx, ny, nz, v)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, w0, nx0, ny0, nz0, nx, ny, nz, w)
      call average3D(1, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, t0, nx0, ny0, nz0, nx, ny, nz, t)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, pnh0, nx0, ny0, nz0, nx, ny, nz, pnh)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qv0, nx0, ny0, nz0, nx, ny, nz, qv)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qc0, nx0, ny0, nz0, nx, ny, nz, qc)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qi0, nx0, ny0, nz0, nx, ny, nz, qi)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qr0, nx0, ny0, nz0, nx, ny, nz, qr)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qs0, nx0, ny0, nz0, nx, ny, nz, qs)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qg0, nx0, ny0, nz0, nx, ny, nz, qg)
      !call ntransg(proj0, proj, slon0, slon, latu, lonu, latv, lonv, du, dv, nx, ny, nz)
      !do k = 1, nz
         !print*,'U:', k, minval(u(:,:,k)), maxval(u(:,:,k))
      !end do
      !
      ! tsoil
      do k = 1, ngm
         call average2D(dx0, dy0, fi, fj, tsoil0(:,:,k), nx0, ny0, nx, ny, tsoil(:,:,k))
         where (landsea < 0.5)
	    tsoil(:,:,k) = sst
         else where
	    tsoil(:,:,k) = tsoil(:,:,k) + gamma*deltazs
         end where
      end do
      ! ts
      where (landsea >= 0.5)
         t(:,:,1) = t(:,:,1) + gamma*deltazs
      end where
      ! ps
      call average2D(dx0, dy0, fi, fj, ps0, nx0, ny0, nx, ny, ps)
      ps = ps*(1+gamma*deltazs/tsoil(:,:,1))**gmrgi
      pnh(:,:,1) = ps
      !
      call write_nhm(outfile, u, v, w, t, pnh, qv, qc, qi, qr, qs, qg, tsoil, nx, ny, nz)
   else
      call read_nonhydro_icbc(infile, 1, 1, 1, nx0, ny0, nz0, pmsl_mean, ptop_mean, ps_mean, ps0, ps0, ps0, ps0, tsoil0, u0, v0, w0, pt0, p0, qv0, qc0, qi0, qr0, qs0, qg0)
      !ps0 = 100.d0*ps0
      p0 = 100.d0*p0
      pt0 = pt0 + ptrf
      ! pnh
      pb0 = p0(:,:,1)
      call pbpt2p(zrp0, dvtransw0, zs0, pb0, pt0, qv0, nx0, ny0, nz0, pnh0)
      pnh0 = p0 - pnh0
      pnh0(:,:,1) = 0.d0
      t0 = pt0*(p0/presrf)**rdvcp
      ps0 = 0.5*(pt0(:,:,1)*(1.+0.608*qv0(:,:,1)) + pt0(:,:,2)*(1.+0.608*qv0(:,:,2))) ! virpt
      ps0 = (pb0/presrf)**rdvcp + gvcp/ps0*zrp0(1)*(1.+dvtransw0(1)*zs0) ! pais
      ps0 = presrf*ps0**cpdvrd
      !
      u0(:,:,1) = u0(:,:,2)
      v0(:,:,1) = v0(:,:,2)
      w0(:,:,1) = w0(:,:,2)
      call average3Duv(hsdiff, dx0, dy0, 1, 1, zrp0, vctransp0, dvtransp0, zrp, vctransp, dvtransp, zs0, zs, fiu, fju, u0, nx0, ny0, nz0, nx, ny, nz, u)
      call average3Duv(hsdiff, dx0, dy0, 2, 2, zrp0, vctransp0, dvtransp0, zrp, vctransp, dvtransp, zs0, zs, fiv, fjv, v0, nx0, ny0, nz0, nx, ny, nz, v)
      call average3D(hsdiff, dx0, dy0, zrw0, vctransw0, zrw, vctransw, zs0, zs, fi, fj, w0, nx0, ny0, nz0, nx, ny, nz, w)
      call average3D(1, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, pt0, nx0, ny0, nz0, nx, ny, nz, pt)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, pnh0, nx0, ny0, nz0, nx, ny, nz, pnh)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qv0, nx0, ny0, nz0, nx, ny, nz, qv)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qc0, nx0, ny0, nz0, nx, ny, nz, qc)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qi0, nx0, ny0, nz0, nx, ny, nz, qi)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qr0, nx0, ny0, nz0, nx, ny, nz, qr)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qs0, nx0, ny0, nz0, nx, ny, nz, qs)
      call average3D(hsdiff, dx0, dy0, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, qg0, nx0, ny0, nz0, nx, ny, nz, qg)
      !call ntransg(proj0, proj, slon0, slon, latu, lonu, latv, lonv, du, dv, nx, ny, nz)
      ! boundary condition
      u(:,:,1) = 0.
      v(:,:,1) = 0.
      w(:,:,1) = 0.
      !
      ! tsoil
      do k = 1, ngm
         call average2D(dx0, dy0, fi, fj, tsoil0(:,:,k), nx0, ny0, nx, ny, tsoil(:,:,k))
         where (landsea < 0.5)
	    tsoil(:,:,k) = sst
         else where
	    tsoil(:,:,k) = tsoil(:,:,k) + gamma*deltazs
         end where
      end do
      where (landsea >= 0.5)
         sst = tsoil(:,:,1) + gamma*zs
      end where
      !
      ! ps, p
      call average2D(dx0, dy0, fi, fj, ps0, nx0, ny0, nx, ny, ps)
      ps = ps*(1+gamma*deltazs/tsoil(:,:,1))**gmrgi
      call pspt2p(zrp, dvtransw, zs, ps, pt, qv, nx, ny, nz, p)
      pnh(:,:,1) = 0.
      p = p + pnh
      ! saturation check
      t = pt*(p/presrf)**rdvcp
      call qv2rh(p, t, qv, nx, ny, nz, rh)
      call rh2qv(p, t, rh, nx, ny, nz, qv)
      !
      ! pmsl, ptop
      !ts0 = t0(:,:,2)*(ps0/p0(:,:,2))**gmrgh
      !pmsl0 = ps0*(1.+gamma*zs0/ts0)**gmrgi
      pmsl0 = ps0*(1.+gamma*zs0/tsoil0(:,:,1))**gmrgi
      do k = 1, nz0
	 h0(:,:,k) = zrp0(k) + zs0*vctransp0(k)
      end do
      call calptop(zrw0(nz0-1), h0, t0, p0, nx0, ny0, nz0, ptop0)
      call average2D(dx0, dy0, fi, fj, pmsl0, nx0, ny0, nx, ny, pmsl)
      call average2D(dx0, dy0, fi, fj, ptop0, nx0, ny0, nx, ny, ptop)
      call calpmean(resolution, resolution, pmsl, ps, ptop, fmap, nx, ny, pmsl_mean, ps_mean, ptop_mean)
      ! pts
      pts = tsoil(:,:,1)*(presrf/ps)**rdvcp
      !
      ! output
      p = 0.01*p
      ps = 0.01*ps
      ps_mean = 0.01*ps_mean
      ptop_mean = 0.01*ptop_mean
      pmsl_mean = 0.01*pmsl_mean
      pt = pt - ptrf
      pts = pts - ptrf
      call write_ic(outfile, 0, idate, vctrans_type, n_vctrans, 10, 20, 10, 20, nz, nz, &
		  & pmsl_mean, ptop_mean, ps_mean, 60., zl_vctrans, zh_vctrans, resolution, resolution, &
		  & resolution, resolution, resolution, resolution, dz2, dz1, dz2, &
		  & pts, dptsdt, sst, ps, tsoil, u, v, w, pt, p, qv, qc, qi, qr, qs, qg, nx, ny, nz)
  
   end if
   !
   deallocate(u0, v0, w0, pt0, p0, h0, t0, pnh0)
   deallocate(qv0, qc0, qi0, qr0, qs0, qg0)
   deallocate(zs0, pb0, pmsl0, ptop0, ts0, ps0, tsoil0)
   deallocate(pts0, us0, vs0, qvs0)
   deallocate(zrp0, zrw0, vctransp0, vctransw0, dvtransp0, dvtransw0)
   deallocate(zrp, zrw, vctransp, vctransw, dvtransp, dvtransw)
   deallocate(zs, landsea, fmap)
   deallocate(fi, fj, fiu, fju, fiv, fjv)
   deallocate(deltazs, ptop, pmsl, t, rh, pnh)
   deallocate(ps, pts, us, vs, ts, qvs, dptsdt, sst, tsoil)
   deallocate(u, v, w, pt, p, qv, qc, qi, qr, qs, qg)
   call MPI_FINALIZE(ierror)
   !
end program outer2inner
