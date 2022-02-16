program nus2grd
   use variable
   use nus2grdlib
   use mpi
   implicit none
   !
   character(6) :: level, dummy
   character(8) :: fgfile='fg000000'
   integer :: nproc, myid
   integer :: islot, ivalid, k, j, i, ngm = 4
   integer :: dnum, ierror
   integer :: gsize(2)
   real :: latlon(4)
   include 'nusdas_fort.h'
   !
   call MPI_INIT(ierror)
   call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierror)
   call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierror)
   !
   ! read namelist
   open(10, file='namelist.txt')
   read(10, control)
   !write(6, control)
   close(10)
   call nwp_ymdhm2seq(iyear, imonth, iday, ihour, iminute, ibase)
   !
   ! projection
   dnum = 1
   call nusdas_inq_def(type1, type2, type3, N_PROJECTION, proj, dnum, ierror)
   if (ierror /= dnum) then
      print*, 'nusdas_inq_def error: N_PROJECTION, ierror=', ierror
      stop
   end if
   dnum = 4
   call nusdas_inq_def(type1, type2, type3, N_STAND_LATLON, latlon, dnum, ierror)
   if (ierror /= dnum) then
      print*, 'nusdas_inq_def error: N_STAND_LATLON, ierror=', ierror
      stop
   end if
   slat = latlon(1)
   !
   ! read nx, ny, nz
   dnum = 2
   call nusdas_inq_def(type1, type2, type3, N_GRID_SIZE, gsize, dnum, ierror)
   if (ierror /= dnum) then
      print*, 'nusdas_inq_def error: N_GRID_SIZE, ierror=', ierror
      stop
   end if
   nx = gsize(1)
   ny = gsize(2)
   call nusdas_subc_eta_inq_nz(type1, type2, type3, ibase, member, ibase, 'ZHYB', nz, ierror)
   if (ierror <= 0) then
      print*, 'nusdas_subc_eta_inq_nz error, ierror=', ierror
      stop
   endif
   !print*, proj, slat, nx, ny, nz
   !
   ! vertical coordinate
   allocate(zrp(nz))
   allocate(zrw(nz))
   allocate(vctrans_p(nz))
   allocate(vctrans_w(nz))
   allocate(dvtrans_p(nz))
   allocate(dvtrans_w(nz))
   call nusdas_subc_zhyb(type1, type2, type3, ibase, member, ibase, nz, ptrf, presrf, zrp, zrw, vctrans_p, vctrans_w, dvtrans_p, dvtrans_w, 'get', ierror)
   if (ierror /= 0) then
      print*, 'nusdas_subc_zhyb error, error=', ierror
      stop
   end if
   !
   ! const fields
   allocate(vrdz2(nz), fmap(nx,ny,6), g2(nx,ny,nz))
   allocate(zs(nx,ny), lat(nx,ny), pairf(nx,ny,nz))
   do k = 1, nz-1
      vrdz2(k) = 1./(zrp(k+1)-zrp(k))
   end do
   vrdz2(nz) = vrdz2(nz-1)
   level = 'SURF  '
   call read_field(type1, type2, type3, ibase, member, ibase, level, 'ZS    ', dtype, nx, ny, ierror, zs)
   call read_field(type1, type2, type3, ibase, member, ibase, level, 'FLAT  ', dtype, nx, ny, ierror, lat)
   do k = 1, nz
      g2(:,:,k) = 1.0 + zs(:,:)*dvtrans_p(k)
   enddo
   call cfmap(proj, slat, lat, nx, ny, fmap)
   do k = 1, nz
      write(dummy, '(i6)') k
      level = adjustl(dummy)
      call read_field(type1, type2, type3, ibase, member, ibase, level, 'PAIRF ', dtype, nx, ny, ierror, pairf(:,:,k))
   end do
   !
   ! allocate
   allocate(pb(nx,ny), us(nx,ny), vs(nx,ny), ts(nx,ny), ps(nx,ny), tds(nx,ny), rhs(nx,ny), qvs(nx,ny))
   allocate(tsoil(nx,ny,ngm))
   allocate(dnsg2(nx,ny,nz), rhou(nx,ny,nz), rhov(nx,ny,nz), rhow(nx,ny,nz), uc(nx,ny,nz), vc(nx,ny,nz), wc(nx,ny,nz), prs(nx,ny,nz))
   allocate(u(nx,ny,nz), v(nx,ny,nz), w(nx,ny,nz), pt(nx,ny,nz), t(nx,ny,nz), p(nx,ny,nz), pnh(nx,ny,nz))
   allocate(qv(nx,ny,nz), qc(nx,ny,nz), qi(nx,ny,nz), qr(nx,ny,nz), qs(nx,ny,nz), qg(nx,ny,nz))
   !
   do islot = 1, nslot
      if (mod(islot-1,nproc) /= myid) cycle
      ivalid = ibase + (islot-1)*interval
      do k = 1, nz
         write(dummy, '(i6)') k
         level = adjustl(dummy)
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
      !
      call uvw_ru2u(fmap, dnsg2, rhou, rhov, rhow, nx, ny, nz, uc, vc, wc)
      call uvw_c2a(vrdz2, uc, vc, wc, nx, ny, nz, u, v, w)
      ! p, pai
      p = prs/g2 + presrf*pairf**cpdvrd
      pb = p(:,:,1)
      pt = pt + ptrf
      call compute_phydro(zrp, dvtrans_w, zs, pb, pt, qv, nx, ny, nz, pnh)
      pnh = p - pnh
      t = pt*(p/presrf)**rdvcp
      ps = 0.5*(pt(:,:,1)*(1.+0.608*qv(:,:,1)) + pt(:,:,2)*(1.+0.608*qv(:,:,2))) ! virpt
      ps = (pb/presrf)**rdvcp + gvcp/ps*zrp(1)*(1.+dvtrans_w(1)*zs) ! pais
      ps = presrf*ps**cpdvrd
      !
      ! surface fields
      level = 'SURF  '
      call read_field(type1s, type2s, type3s, ibase, member, ivalid, level, 'U     ', dtype, nx, ny, ierror, us)
      call read_field(type1s, type2s, type3s, ibase, member, ivalid, level, 'V     ', dtype, nx, ny, ierror, vs)
      call read_field(type1s, type2s, type3s, ibase, member, ivalid, level, 'T     ', dtype, nx, ny, ierror, ts)
      !call read_field(type1s, type2s, type3s, ibase, member, ivalid, level, 'P     ', dtype, nx, ny, ierror, ps)
      call read_field(type1s, type2s, type3s, ibase, member, ivalid, level, 'TTD   ', dtype, nx, ny, ierror, tds)
      !ps = 100.*ps
      tds = ts - tds
      call compute_rhs(ps, ts, tds, nx, ny, rhs)
      call rh2qv(ps, ts, rhs, nx, ny, qvs)
      !
      u(:,:,1) = us
      v(:,:,1) = vs
      t(:,:,1) = ts
      pnh(:,:,1) = ps
      qv(:,:,1) = qvs
      write(fgfile(3:4),'(I2.2)') islot
      call write_nhm(fgfile, u, v, w, t, pnh, qv, qc, qi, qr, qs, qg, tsoil, nx, ny, nz, ngm)
   enddo
   !
   deallocate(zrp, zrw, vctrans_p, vctrans_w, dvtrans_p, dvtrans_w)
   deallocate(vrdz2, fmap, g2)
   deallocate(zs, lat, pairf)
   deallocate(tsoil)
   deallocate(dnsg2, rhou, rhov, rhow, uc, vc, wc, prs)
   deallocate(pb, us, vs, ts, ps, tds, rhs, qvs)
   deallocate(u, v, w, pt, t, p, pnh, qv, qc, qi, qr, qs, qg)
   call MPI_FINALIZE(ierror)
   !
end program nus2grd

