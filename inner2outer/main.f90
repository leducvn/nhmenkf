program inner2outer
   use variable
   use inner2outerlib
   implicit none
   !
   character(8) :: infile='xdriving'
   character(6) :: outfile='xmodel'
   integer :: hsdiff = 0, k
   !
   ! read namelist
   open(10, file='namelist.txt')
   read(10, model)
   read(10, driving)
   close(10)
   resolution = 1000.*resolution
   resolution0 = 1000.*resolution0
   !
   ! model
   allocate(vdz(nz), vrdz(nz), vrdz2(nz))
   allocate(zrp(nz), vctransp(nz), dvtransp(nz))
   allocate(zrw(nz), vctransw(nz), dvtransw(nz))
   allocate(zs(nx,ny), landsea(nx,ny))
   allocate(lon(nx,ny), lat(nx,ny), lonu(nx,ny), latu(nx,ny), lonv(nx,ny), latv(nx,ny))
   ! horizontal
   call read_nhmcst(30, nx, ny, lon, lat, zs, landsea)
   call setorg(proj, resolution, slat, slon, xi, xj, xlat, xlon, lat, lon, nx, ny, latu, lonu, latv, lonv)
   ! vertical
   call vrgdis(nz, nz, 1, nz, nz, 1, dz2, dz1, dz2, vdz, vrdz, vrdz2)
   call setzrp(vrdz2, nz, zrp)
   call calc_zcoordinate(vctrans_type, n_vctrans, ztop, zl_vctrans, zh_vctrans, zrp, nz, vctransp, dvtransp)
   call setzrw(vrdz, nz, zrw)
   call calc_zcoordinate(vctrans_type, n_vctrans, zrw(nz-1), zl_vctrans, zh_vctrans, zrw, nz, vctransw, dvtransw)
   deallocate(vdz, vrdz, vrdz2)
   deallocate(landsea)
   !
   ! driving
   allocate(vdz0(nz0), vrdz0(nz0), vrdz20(nz0))
   allocate(zrp0(nz0), vctransp0(nz0), dvtransp0(nz0))
   allocate(zrw0(nz0), vctransw0(nz0), dvtransw0(nz0))
   allocate(zs0(nx0,ny0), landsea0(nx0,ny0), lon0(nx0,ny0), lat0(nx0,ny0))
   ! horizontal
   call read_nhmcst(40, nx0, ny0, lon0, lat0, zs0, landsea0)
   call vrgdis(nz0, nz0, 1, nz0, nz0, 1, dz20, dz10, dz20, vdz0, vrdz0, vrdz20)
   call setzrp(vrdz20, nz0, zrp0)
   call calc_zcoordinate(vctrans_type0, n_vctrans0, ztop0, zl_vctrans0, zh_vctrans0, zrp0, nz0, vctransp0, dvtransp0)
   call setzrw(vrdz0, nz0, zrw0)
   call calc_zcoordinate(vctrans_type0, n_vctrans0, zrw0(nz0-1), zl_vctrans0, zh_vctrans0, zrw0, nz0, vctransw0, dvtransw0)
   deallocate(vdz0, vrdz0, vrdz20)
   deallocate(landsea0, lon0, lat0)
   !
   ! coordinate
   allocate(fi(nx,ny), fj(nx,ny), fiu(nx,ny), fju(nx,ny), fiv(nx,ny), fjv(nx,ny))
   call setfifj(proj0, nx0, ny0, resolution0, resolution0, slat0, slon0, xi0, xj0, xlat0, xlon0, &
              & latu, lonu, latv, lonv, lat, lon, nx, ny, fi, fj, fiu, fju, fiv, fjv)
   deallocate(lon, lat, lonu, latu, lonv, latv)
   !
   ! allocate
   ! large domain
   allocate(dps0(nx0,ny0), dtsoil0(nx0,ny0,ngm))
   allocate(du0(nx0,ny0,nz0), dv0(nx0,ny0,nz0), dw0(nx0,ny0,nz0), dt0(nx0,ny0,nz0), dpnh0(nx0,ny0,nz0))
   allocate(dqv0(nx0,ny0,nz0), dqc0(nx0,ny0,nz0), dqi0(nx0,ny0,nz0), dqr0(nx0,ny0,nz0), dqs0(nx0,ny0,nz0), dqg0(nx0,ny0,nz0))
   ! small domain
   allocate(dps(nx,ny), dtsoil(nx,ny,ngm))
   allocate(du(nx,ny,nz), dv(nx,ny,nz), dw(nx,ny,nz), dt(nx,ny,nz), dpnh(nx,ny,nz))
   allocate(dqv(nx,ny,nz), dqc(nx,ny,nz), dqi(nx,ny,nz), dqr(nx,ny,nz), dqs(nx,ny,nz), dqg(nx,ny,nz))
   !
   ! Read increment fields
   call read_nhm(infile, nx0, ny0, nz0, du0, dv0, dw0, dt0, dpnh0, dqv0, dqc0, dqi0, dqr0, dqs0, dqg0, dtsoil0)
   dps0(:,:) = dpnh0(:,:,1)
   dpnh0(:,:,1) = 0.d0
   !
   ! 2D interpolation
   do k = 1, ngm
      call interpolate2D(fi, fj, dtsoil0(:,:,k), nx0, ny0, nx, ny, dtsoil(:,:,k))
   end do
   call interpolate2D(fi, fj, dps0, nx0, ny0, nx, ny, dps)
   !
   ! 3D interpolation
   call interpolate3Duv(hsdiff, 0, 0, zrp0, vctransp0, dvtransp0, zrp, vctransp, dvtransp, zs0, zs, fiu, fju, du0, nx0, ny0, nz0, nx, ny, nz, du)
   call interpolate3Duv(hsdiff, 0, 0, zrp0, vctransp0, dvtransp0, zrp, vctransp, dvtransp, zs0, zs, fiv, fjv, dv0, nx0, ny0, nz0, nx, ny, nz, dv)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dw0, nx0, ny0, nz0, nx, ny, nz, dw)
   call interpolate3D(1, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dt0, nx0, ny0, nz0, nx, ny, nz, dt)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dpnh0, nx0, ny0, nz0, nx, ny, nz, dpnh)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dqv0, nx0, ny0, nz0, nx, ny, nz, dqv) 
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dqc0, nx0, ny0, nz0, nx, ny, nz, dqc)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dqi0, nx0, ny0, nz0, nx, ny, nz, dqi)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dqr0, nx0, ny0, nz0, nx, ny, nz, dqr)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dqs0, nx0, ny0, nz0, nx, ny, nz, dqs)
   call interpolate3D(hsdiff, zrp0, vctransp0, zrp, vctransp, zs0, zs, fi, fj, dqg0, nx0, ny0, nz0, nx, ny, nz, dqg)
   !call ntransg(proj0, proj, slon0, slon, latu, lonu, latv, lonv, du, dv, nx, ny, nz)
   !
   dpnh(:,:,1) = dps(:,:)
   call write_nhm(outfile, du, dv, dw, dt, dpnh, dqv, dqc, dqi, dqr, dqs, dqg, dtsoil, nx, ny, nz)
   !
   ! deallocate
   deallocate(zrp0, zrw0, vctransp0, vctransw0, dvtransp0, dvtransw0)
   deallocate(zrp, zrw, vctransp, vctransw, dvtransp, dvtransw)
   deallocate(fi, fj, fiu, fju, fiv, fjv)
   deallocate(zs0, zs)
   deallocate(du0, dv0, dw0, dt0, dpnh0, dqv0, dqc0, dqi0, dqr0, dqs0, dqg0, dtsoil0, dps0)
   deallocate(du, dv, dw, dt, dpnh, dqv, dqc, dqi, dqr, dqs, dqg, dtsoil, dps)
   !
end program inner2outer
