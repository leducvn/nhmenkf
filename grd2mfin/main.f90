program grd2mfin
   use variable
   use grd2mfinlib
   implicit none
   !
   character(4) :: meanfile='mean', pertfile='pert'
   integer :: k, idate(5)
   real :: pmsl_mean, ps_mean, ptop_mean
   !
   ! read namelist
   open(10, file='namelist.txt')
   read(10, control)
   read(10, model)
   close(10)
   resolution = 1000.*resolution
   idate(1) = iyear
   idate(2) = imonth
   idate(3) = iday
   idate(4) = ihour
   idate(5) = iminute
   !
   ! inner
   allocate(vdz(nz), vrdz(nz), vrdz2(nz))
   allocate(zrp(nz), vctransp(nz), dvtransp(nz))
   allocate(zrw(nz), vctransw(nz), dvtransw(nz))
   allocate(lon(nx,ny), lat(nx,ny), zs(nx,ny), landsea(nx,ny), fmap(nx,ny,6))
   ! horizontal
   call read_nhmcst(30, nx, ny, lon, lat, zs, landsea)
   call cfmap(proj, slat, lat, nx, ny, fmap)
   ! vertical
   call vrgdis(nz, nz, 1, nz, nz, 1, dz2, dz1, dz2, vdz, vrdz, vrdz2)
   call setzrp(vrdz2, nz, zrp)
   call calc_zcoordinate(vctrans_type, n_vctrans, ztop, zl_vctrans, zh_vctrans, zrp, nz, vctransp, dvtransp)
   call setzrw(vrdz, nz, zrw)
   call calc_zcoordinate(vctrans_type, n_vctrans, zrw(nz-1), zl_vctrans, zh_vctrans, zrw, nz, vctransw, dvtransw)
   deallocate(vdz, vrdz2, lon, lat)
   !
   ! allocate
   allocate(pb(nx,ny), pts(nx,ny))
   allocate(tsoil(nx,ny,ngm))
   allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), pt(nx,ny,nz), pnh(nx,ny,nz))
   allocate(qv(nx,ny,nz), qc(nx,ny,nz), qi(nx,ny,nz), qr(nx,ny,nz), qs(nx,ny,nz), qg(nx,ny,nz))
   allocate(dps(nx,ny), dtsoil(nx,ny,ngm))
   allocate(dua(nx,ny,nz), dva(nx,ny,nz), dwa(nx,ny,nz), du(nx,ny,nz), dv(nx,ny,nz), dw(nx,ny,nz), dt(nx,ny,nz), dpnh(nx,ny,nz))
   allocate(dqv(nx,ny,nz), dqc(nx,ny,nz), dqi(nx,ny,nz), dqr(nx,ny,nz), dqs(nx,ny,nz), dqg(nx,ny,nz))
   allocate(ps(nx,ny), dptsdt(nx,ny), sst(nx,ny), ptop(nx,ny), pmsl(nx,ny))
   allocate(p(nx,ny,nz), ph(nx,ny,nz), t(nx,ny,nz), rh(nx,ny,nz), h(nx,ny,nz))
   !
   ! Read mean and perturbation fields
   open(90, file=meanfile, form='unformatted')
   if (hydrostatic == 1) then
      call read_hydro_icbc(90, nx, ny, nz, pmsl_mean, ptop_mean, ps_mean, pts, dptsdt, sst, ps, tsoil, u, v, pt, qv)
      ps = 100.d0*ps
      pt = pt + ptrf
      pb = 0.5*(pt(:,:,1)*(1.+0.608*qv(:,:,1)) + pt(:,:,2)*(1.+0.608*qv(:,:,2))) ! virpt
      pb = (ps/presrf)**rdvcp - gvcp/pb*zrp(1)*(1.+dvtransw(1)*zs) ! paib
      pb = presrf*pb**cpdvrd
      call pbpt2p(zrp, dvtransw, zs, pb, pt, qv, nx, ny, nz, p)
      pnh(:,:,:) = 0.d0
      t = pt*(p/presrf)**rdvcp
   else
      call read_nonhydro_icbc(90, 1, 1, 1, nx, ny, nz, pmsl_mean, ptop_mean, ps_mean, pts, dptsdt, sst, ps, tsoil, u, v, w, pt, p, qv, qc, qi, qr, qs, qg)
      !ps = 100.d0*ps
      pt = pt + ptrf
      p = 100.d0*p
      t = pt*(p/presrf)**rdvcp
      pb = p(:,:,1)
      call pbpt2p(zrp, dvtransw, zs, pb, pt, qv, nx, ny, nz, ph)
      pnh = p - ph
      pnh(:,:,1) = 0.d0
      ps = 0.5*(pt(:,:,1)*(1.+0.608*qv(:,:,1)) + pt(:,:,2)*(1.+0.608*qv(:,:,2))) ! virpt
      ps = (pb/presrf)**rdvcp + gvcp/ps*zrp(1)*(1.+dvtransw(1)*zs) ! pais
      ps = presrf*ps**cpdvrd
   end if
   !
   call read_nhm(pertfile, nx, ny, nz, dua, dva, dwa, dt, dpnh, dqv, dqc, dqi, dqr, dqs, dqg, dtsoil)
   dps(:,:) = dpnh(:,:,1)
   dua(:,:,1) = 0.d0
   dva(:,:,1) = 0.d0
   dwa(:,:,1) = 0.d0
   dpnh(:,:,1) = 0.d0
   dt(:,:,1) = dt(:,:,2)
   dqv(:,:,1) = dqv(:,:,2)
   dqc(:,:,1) = dqc(:,:,2)
   dqi(:,:,1) = dqi(:,:,2)
   dqr(:,:,1) = dqr(:,:,2)
   dqs(:,:,1) = dqs(:,:,2)
   dqg(:,:,1) = dqg(:,:,2)
   call uvw_a2c(vrdz, dua, dva, dwa, nx, ny, nz, du, dv, dw)
   !
   ! add perturbations
   u = u + du
   v = v + dv
   w = w + dw
   t = t + dt
   pnh = pnh + dpnh
   qv = qv + dqv
   qc = qc + dqc
   qi = qi + dqi
   qr = qr + dqr
   qs = qs + dqs
   qg = qg + dqg
   tsoil = tsoil + dtsoil
   ps = ps + dps
   where (qv < 0.)
      qv = 0.
   end where
   where (qc < 0.)
      qc = 0.
   end where
   where (qi < 0.)
      qi = 0.
   end where
   where (qr < 0.)
      qr = 0.
   end where
   where (qs < 0.)
      qs = 0.
   end where
   where (qg < 0.)
      qg = 0.
   end where
   !
   ! adjust
   call pst2p(zrp, dvtransw, zs, ps, t, qv, nx, ny, nz, p)
   p = p + pnh
   pt = t*(presrf/p)**rdvcp
   ! saturation check
   call qv2rh(p, t, qv, nx, ny, nz, rh)
   call rh2qv(p, t, rh, nx, ny, nz, qv)
   !
   ! surface fields
   do k = 1, ngm
      where (landsea < 0.5)
         tsoil(:,:,k) = sst
      end where
   end do
   pts = tsoil(:,:,1)*(presrf/ps)**rdvcp
   where (landsea >= 0.5)
      sst = tsoil(:,:,1) + gamma*zs
   end where
   pmsl = ps*(1.+gamma*zs/tsoil(:,:,1))**gmrgi
   do k = 1, nz
      h(:,:,k) = zrp(k) + zs*vctransp(k)
   enddo
   call calptop(zrw(nz-1), h, t, p, nx, ny, nz, ptop)
   call calpmean(resolution, resolution, pmsl, ps, ptop, fmap, nx, ny, pmsl_mean, ps_mean, ptop_mean)
   dptsdt = 0.d0
   !
   p = 0.01*p
   ps = 0.01*ps
   ps_mean = 0.01*ps_mean
   ptop_mean = 0.01*ptop_mean
   pmsl_mean = 0.01*pmsl_mean
   pt = pt - ptrf
   pts = pts - ptrf
   if (hydrostatic == 1) then
      call write_hydro_ic(90, idate, vctrans_type, n_vctrans, 10, 20, 10, 20, nz, nz, &
		  & pmsl_mean, ptop_mean, ps_mean, 60., zl_vctrans, zh_vctrans, resolution, resolution, &
		  & resolution, resolution, resolution, resolution, dz2, dz1, dz2, &
		  & pts, dptsdt, sst, ps, tsoil, u, v, pt, qv, nx, ny, nz)
   else
      call write_nonhydro_ic(90, idate, vctrans_type, n_vctrans, 10, 20, 10, 20, nz, nz, &
		  & pmsl_mean, ptop_mean, ps_mean, 60., zl_vctrans, zh_vctrans, resolution, resolution, &
		  & resolution, resolution, resolution, resolution, dz2, dz1, dz2, &
		  & pts, dptsdt, sst, ps, tsoil, u, v, w, pt, p, qv, qc, qi, qr, qs, qg, nx, ny, nz)
   end if
   deallocate(zrp, zrw, vctransp, vctransw, dvtransp, dvtransw)
   deallocate(vrdz, zs, landsea, fmap)
   deallocate(dua, dva, dwa, dps, du, dv, dw, dt, dpnh, dqv, dqc, dqi, dqr, dqs, dqg, dtsoil)
   deallocate(pb, ptop, pmsl, h, t, rh, ph, pnh)
   deallocate(ps, pts, dptsdt, sst, tsoil)
   deallocate(u, v, w, pt, p, qv, qc, qi, qr, qs, qg)
   !
   !stop
end program grd2mfin

