program nus2mfin
   use variable
   use nus2mfinlib
   implicit none
   !
   character(6) :: level, dummy
   character(4) :: outfile = 'mfin'
   integer :: idate(5), ibase, ivalid, k, ierror
   real :: pmsl_mean, ps_mean, ptop_mean
   include 'nusdas_fort.h'
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
   call nwp_ymdhm2seq(iyear, imonth, iday, ihour, iminute, ivalid)
   ibase = ivalid - tend
   !
   ! inner
   allocate(vdz(nz), vrdz(nz), vrdz2(nz))
   allocate(zrp(nz), vctransp(nz), dvtransp(nz))
   allocate(zrw(nz), vctransw(nz), dvtransw(nz))
   allocate(lon(nx,ny), lat(nx,ny), zs(nx,ny), landsea(nx,ny), fmap(nx,ny,6), g2(nx,ny,nz))
   ! horizontal
   call read_nhmcst(30, nx, ny, lon, lat, zs, landsea)
   call cfmap(proj, slat, lat, nx, ny, fmap)
   ! vertical
   call vrgdis(nz, nz, 1, nz, nz, 1, dz2, dz1, dz2, vdz, vrdz, vrdz2)
   call setzrp(vrdz2, nz, zrp)
   call calc_zcoordinate(vctrans_type, n_vctrans, ztop, zl_vctrans, zh_vctrans, zrp, nz, vctransp, dvtransp)
   call setzrw(vrdz, nz, zrw)
   call calc_zcoordinate(vctrans_type, n_vctrans, zrw(nz-1), zl_vctrans, zh_vctrans, zrw, nz, vctransw, dvtransw)
   do k = 1, nz
      g2(:,:,k) = 1.0 + zs(:,:)*dvtransp(k)
   enddo
   deallocate(vdz, vrdz, vrdz2, lon, lat)
   !
   ! allocate
   allocate(tsoil(nx,ny,ngm))
   allocate(pairf(nx,ny,nz), dnsg2(nx,ny,nz), rhou(nx,ny,nz), rhov(nx,ny,nz), rhow(nx,ny,nz))
   allocate(pt(nx,ny,nz), prs(nx,ny,nz))
   allocate(qv(nx,ny,nz), qc(nx,ny,nz), qi(nx,ny,nz), qr(nx,ny,nz), qs(nx,ny,nz), qg(nx,ny,nz))
   allocate(ps(nx,ny), pb(nx,ny), pts(nx,ny), dptsdt(nx,ny), sst(nx,ny), ts(nx,ny), ptop(nx,ny), pmsl(nx,ny))
   allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz))
   allocate(t(nx,ny,nz), p(nx,ny,nz), h(nx,ny,nz))
   !
   ! read
   do k = 1, nz
      write(dummy, '(i6)') k
      level = adjustl(dummy)
      call read_field(type1, type2, type3, ibase, member, ibase,  level, 'PAIRF ', dtype, nx, ny, ierror, pairf(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'DNSG2 ', dtype, nx, ny, ierror, dnsg2(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'RU    ', dtype, nx, ny, ierror, rhou(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'RV    ', dtype, nx, ny, ierror, rhov(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'RW    ', dtype, nx, ny, ierror, rhow(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'PT    ', dtype, nx, ny, ierror, pt(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'PRS   ', dtype, nx, ny, ierror, prs(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'QV    ', dtype, nx, ny, ierror, qv(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'QC    ', dtype, nx, ny, ierror, qc(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'QCI   ', dtype, nx, ny, ierror, qi(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'QR    ', dtype, nx, ny, ierror, qr(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'QS    ', dtype, nx, ny, ierror, qs(:,:,k))
      call read_field(type1, type2, type3, ibase, member, ivalid, level, 'QG    ', dtype, nx, ny, ierror, qg(:,:,k))
      if (k <= ngm) call read_field(type1, type2, type3, ibase, member, ivalid, level, 'TIN   ', dtype, nx, ny, ierror, tsoil(:,:,k))
   enddo
   call uvw_ru2u(fmap, dnsg2, rhou, rhov, rhow, nx, ny, nz, u, v, w)
   p = prs/g2 + presrf*pairf**cpdvrd
   pb = p(:,:,1)
   pt = pt + ptrf
   t = pt*(p/presrf)**rdvcp
   ps = 0.5*(t(:,:,1)*(1.+0.608*qv(:,:,1))+t(:,:,2)*(1.+0.608*qv(:,:,2))) ! virt
   ps = g0/rd/ps*zrp(1)*(1.+dvtransw(1)*zs)
   ps = pb*exp(ps)
   !
   ! surface fields
   call read_sst(20, nx, ny, sst)
   do k = 1, ngm
      where (landsea < 0.5)
         tsoil(:,:,k) = sst
      end where
   end do
   pts = tsoil(:,:,1)*(presrf/ps)**rdvcp
   where (landsea >= 0.5)
      sst = tsoil(:,:,1) + gamma*zs
   end where
   !ts = t(:,:,2)*(ps/p(:,:,2))**gmrgh
   !pmsl = ps*(1.+gamma*zs/ts)**gmrgi
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
   call write_ic(90, idate, vctrans_type, n_vctrans, 10, 20, 10, 20, nz, nz, &
               & pmsl_mean, ptop_mean, ps_mean, 60., zl_vctrans, zh_vctrans, resolution, resolution, &
	       & resolution, resolution, resolution, resolution, dz2, dz1, dz2, &
	       & pts, dptsdt, sst, ps, tsoil, u, v, w, pt, p, qv, qc, qi, qr, qs, qg, nx, ny, nz)
   deallocate(pb, g2, pairf, dnsg2, rhou, rhov, rhow, prs)
   deallocate(zrp, zrw, vctransp, vctransw, dvtransp, dvtransw)
   deallocate(zs, landsea, fmap)
   deallocate(ts, ptop, pmsl, h, t)
   deallocate(ps, pts, dptsdt, sst, tsoil)
   deallocate(u, v, w, pt, p, qv, qc, qi, qr, qs, qg)
   !
   !stop
end program nus2mfin
