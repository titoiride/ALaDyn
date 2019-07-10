 !*****************************************************************************************************!
 !                            Copyright 2008-2018  The ALaDyn Collaboration                            !
 !*****************************************************************************************************!

 !*****************************************************************************************************!
 !  This file is part of ALaDyn.                                                                       !
 !                                                                                                     !
 !  ALaDyn is free software: you can redistribute it and/or modify                                     !
 !  it under the terms of the GNU General Public License as published by                               !
 !  the Free Software Foundation, either version 3 of the License, or                                  !
 !  (at your option) any later version.                                                                !
 !                                                                                                     !
 !  ALaDyn is distributed in the hope that it will be useful,                                          !
 !  but WITHOUT ANY WARRANTY; without even the implied warranty of                                     !
 !  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                                      !
 !  GNU General Public License for more details.                                                       !
 !                                                                                                     !
 !  You should have received a copy of the GNU General Public License                                  !
 !  along with ALaDyn.  If not, see <http://www.gnu.org/licenses/>.                                    !
 !*****************************************************************************************************!

 module grid_part_lib

 use common_param
 use grid_param

 implicit none

 real(dp),parameter :: half=0.5,thr=0.75,two_third=2./3.,one_sixth=1./6.
 real(sp),parameter :: shx=3.,shy=3.,shz=3.
 integer(kind=2) :: err_ind

 contains
 !    Templates of spl=1,2,3 order shape functions for grid-particle  connection
 !==================================================
 ! Computes 
 subroutine set_int_pshape(spl,xx,ax,ind)
 integer,intent(in) :: spl
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: ax(0:3)
 integer,intent(out) :: ind
 real(dp) :: sx,sx2,sx3
 !To integer grid points
 ax(0:3)=0.0
 select case(spl)
 case(1)
  ind=int(xx)
  ax(1)=xx-real(ind,dp)
  ax(0)=1.-ax(1)
 case(2)
  ind=int(xx+0.5)
  sx=xx-real(ind,dp)
  sx2=sx*sx
  ax(1)=0.75-sx2
  ax(2)=0.5*(0.25+sx2+sx)
  ax(0)=1.-ax(1)-ax(2)
 case(3)
  ind=int(xx)
  sx=xx-real(ind,dp)
  sx2=sx*sx
  sx3=sx2*sx
  ax(1)=two_third-sx2+0.5*sx3
  ax(2)=one_sixth+0.5*(sx+sx2-sx3)
  ax(3)=one_sixth*sx3
  ax(0)=1.-ax(1)-ax(2)-ax(3)
 end select
 end subroutine set_int_pshape
!========================
 subroutine set_hint_pshape(spl,xx,ax,ind)
 integer,intent(in) :: spl
 real(dp),intent(in) :: xx
 integer,intent(out) :: ind
 real(dp),intent(out) :: ax(0:3)
 real(dp) :: sx,sx2,sx3
 !To half-integer grid points
 ax(0:3)=0.0
 select case(spl)
 case(1)
  sx=xx+0.5
  ind=int(sx)
  ax(1)=sx-real(ind,dp)
  ax(0)=1.-ax(1)
 case(2)
  ind=int(xx)
  sx=xx-0.5-real(ind,dp)
  sx2=sx*sx
  ax(1)=0.75-sx2
  ax(2)=0.5*(0.25+sx2+sx)
  ax(0)=1.-ax(1)-ax(2)
 case(3)
  ind=int(xx+0.5)
  sx=xx-real(ind,dp)
  sx2=sx*sx
  sx3=sx2*sx
  ax(1)=two_third-sx2+0.5*sx3
  ax(2)=one_sixth+0.5*(sx+sx2-sx3)
  ax(3)=one_sixth*sx3
  ax(0)=1.-ax(1)-ax(2)-ax(3) ! to (i-1/2,i+1/2,i+3/2,i+5/2) half-int
 end select
 end subroutine set_hint_pshape
 !=======================
 !============== Mapping for stretched grids==========

 subroutine map2dy_part_sind(np,sind,ic1,ym,pt)
 integer,intent(in) :: np,sind,ic1
 real(dp),intent(in) :: ym
 real(dp),intent(inout) :: pt(:,:)
 real(dp) :: yp,ys,ximn, yp_loc
 integer :: n
 !========================
 !  enter the y=part(ic1,n) particle position in stretched grid
 !            y=y(xi)
 !  exit      xi=part(ic1,n) the  particle position in uniform grid
 !               normalized to the Dxi cell size
 !==========================================
 select case(sind)
 case(1)       !y<0
  ys=str_ygrid%smin
  ximn=-dyi_inv*atan(sy_rat*(ym-ys))
  do n=1,np
   yp=pt(n,ic1)
   yp_loc=ximn+dy_inv*(yp-ys)
   if(yp <= ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   pt(n,ic1)= yp_loc
  end do
 case(2)       !y>0
  ys=str_ygrid%smax
  if(ym>ys)then
   ximn=dyi_inv*atan(sy_rat*(ys-ym))
  else
   ximn=dy_inv*(ys-ym)
  endif
  do n=1,np
   yp=pt(n,ic1)
   yp_loc=dy_inv*(yp-ym)
   if(yp > ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   pt(n,ic1)= yp_loc
  end do
 end select
 end subroutine map2dy_part_sind

 subroutine map2dz_part_sind(np,sind,ic1,zm,pt)
 integer,intent(in) :: np,sind,ic1
 real(dp),intent(in) :: zm
 real(dp),intent(inout) :: pt(:,:)
 real(dp) :: zp,zs,zimn,zp_loc
 integer :: n
 !========================
 !  enter the y=part(ic1,n) particle position in stretched grid
 !            y=y(xi)
 !  exit      xi=part(ic1,n) the  particle position in uniform grid
 !               normalized to the Dxi cell size
 !==========================================
 select case(sind)
 case(1)       !z<0
  zs=str_zgrid%smin
  zimn=-dzi_inv*atan(sz_rat*(zm-zs))
  do n=1,np
   zp=pt(n,ic1)
   zp_loc=zimn+dz_inv*(zp-zs)
   if(zp <= zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic1)= zp_loc
  end do
 case(2)       !z>0
  zs=str_zgrid%smax
  if(zm>zs)then
   zimn=dzi_inv*atan(sz_rat*(zs-zm))
  else
   zimn=dz_inv*(zs-zm)
  endif
  do n=1,np
   zp=pt(n,ic1)
   zp_loc=dz_inv*(zp-zm)
   if(zp > zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic1)= zp_loc
  end do
 end select
 end subroutine map2dz_part_sind
 !--------------------------

 subroutine map3d_part_sind(pt,np,sind,ic1,ic2,ym,zm)
 integer,intent(in) :: np,sind,ic1,ic2
 real(dp),intent(in) :: ym,zm
 real(dp),intent(inout) :: pt(:,:)
 real(dp) :: yp,zp,yp_loc,zp_loc,ys,zs,ximn,zimn
 integer :: n
 !========================
 !  enter the y=part(n,ic1) z=part(n,ic2) particle positions
 !        in stretched grids    y=y(xi), z(zi)
 !  exit   xi=part(n,ic1) zi=part(n,ic2)
 !    particle positions in uniform grid
 !    normalized to the (Dxi Dzi) cell sizes
 !==========================================

 select case(sind)
 case(1)       !y<0 z<0 corner ys>ymn
  ys=str_ygrid%smin
  ximn=-dyi_inv*atan(sy_rat*(ym-ys))
  zs=str_zgrid%smin
  zimn=-dzi_inv*atan(sz_rat*(zm-zs))
  do n=1,np
   yp=pt(n,ic1)
   zp=pt(n,ic2)
   yp_loc=ximn+dy_inv*(yp-ys)
   zp_loc=zimn+dz_inv*(zp-zs)
   if(yp < ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp < zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic1)=yp_loc
   pt(n,ic2)=zp_loc
  end do
 case(2)       !z<0
  zs=str_zgrid%smin
  zimn=-dzi_inv*atan(sz_rat*(zm-zs))
  do n=1,np
   yp=pt(n,ic1)
   pt(n,ic1)=dy_inv*(yp-ym)
   zp=pt(n,ic2)
   zp_loc=zimn+dz_inv*(zp-zs)
   if(zp < zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic2)=zp_loc
  end do
 case(3)       !y>0 z<0 corner
  ys=str_ygrid%smax
  if(ym>ys)then
   ximn=dyi_inv*atan(sy_rat*(ys-ym))
  else
   ximn=dy_inv*(ys-ym)
  endif
  zs=str_zgrid%smin
  zimn=-dzi_inv*atan(sz_rat*(zm-zs))
  do n=1,np
   yp=pt(n,ic1)
   zp=pt(n,ic2)
   yp_loc=dy_inv*(yp-ym)
   zp_loc=zimn+dz_inv*(zp-zs)
   if(yp > ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp < zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic1)=yp_loc
   pt(n,ic2)=zp_loc
  end do
 case(4)       !y>0
  ys=str_ygrid%smax
  if(ym>ys)then
   ximn=dyi_inv*atan(sy_rat*(ys-ym))
  else
   ximn=dy_inv*(ys-ym)
  endif
  do n=1,np
   yp=pt(n,ic1)
   zp=pt(n,ic2)
   yp_loc=dy_inv*(yp-ym)
   if(yp > ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   pt(n,ic1)=yp_loc
   pt(n,ic2)=dz_inv*(zp-zm)
  end do
 case(5)       !y>0 z>0 corner
  ys=str_ygrid%smax
  if(ym>ys)then
   ximn=dyi_inv*atan(sy_rat*(ys-ym))
  else
   ximn=dy_inv*(ys-ym)
  endif
  zs=str_zgrid%smax
  if(zm>zs)then
   zimn=dzi_inv*atan(sz_rat*(zs-zm))
  else
   zimn=dz_inv*(zs-zm)
  endif
  do n=1,np
   yp=pt(n,ic1)
   zp=pt(n,ic2)
   yp_loc=dy_inv*(yp-ym)
   zp_loc=dz_inv*(zp-zm)
   if(yp > ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp > zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic1)=yp_loc
   pt(n,ic2)=zp_loc
  end do
 case(6)       !z>0
  zs=str_zgrid%smax
  if(zm>zs)then
   zimn=dzi_inv*atan(sz_rat*(zs-zm))
  else
   zimn=dz_inv*(zs-zm)
  endif
  do n=1,np
   yp=pt(n,ic1)
   pt(n,ic1)=dy_inv*(yp-ym)
   zp=pt(n,ic2)
   zp_loc=dz_inv*(zp-zm)
   if(zp > zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic2)=zp_loc
  end do
 case(7)       !y<0 z>0 corner
  ys=str_ygrid%smin
  ximn=-dyi_inv*atan(sy_rat*(ym-ys))
  zs=str_zgrid%smax
  if(zm>zs)then
   zimn=dzi_inv*atan(sz_rat*(zs-zm))
  else
   zimn=dz_inv*(zs-zm)
  endif
  do n=1,np
   yp=pt(n,ic1)
   zp=pt(n,ic2)
   yp_loc=ximn+dy_inv*(yp-ys)
   zp_loc=dz_inv*(zp-zm)
   if(yp < ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp > zs)zp_loc=zimn+dzi_inv*atan(sz_rat*(zp-zs))
   pt(n,ic1)=yp_loc
   pt(n,ic2)=zp_loc
  end do
 case(8)       !y<0
  ys=str_ygrid%smin
  ximn=-dyi_inv*atan(sy_rat*(ym-ys))
  do n=1,np
   yp=pt(n,ic1)
   zp=pt(n,ic2)
   yp_loc=ximn+dy_inv*(yp-ys)
   if(yp < ys)yp_loc=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   pt(n,ic1)=yp_loc
   pt(n,ic2)=dz_inv*(zp-zm)
  end do
 end select
 end subroutine map3d_part_sind
 !========================================
 !===========================================
 !     Computes particle density charge on grid integer points
 !==============================================
 !DIR$ ATTRIBUTES INLINE :: ql_interpolate,qq_interpolate
!====================
  subroutine qlh_2d_spline(xp,ax,axh,ay,ayh,ix,ihx,iy,ihy)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),axh(0:1),ay(0:2),ayh(0:1)
   integer,intent(inout) :: ix,ihx,iy,ihy
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   axh(1)=sx+0.5
   axh(0)=1.-axh(1)

   xx=shy+xp(2)
   iy=int(xx+0.5)
   sx=xx-real(iy,dp)
   sx2=sx*sx
   ay(1)=0.75-sx2
   ay(2)=0.5*(0.25+sx2+sx)
   ay(0)=1.-ay(1)-ay(2)

   ayh(1)=sx+0.5
   ayh(0)=1.-ayh(1)

   ix=ix-1
   ihx=ix
   iy=iy-1
   ihy=iy

  end subroutine qlh_2d_spline
!====================
  subroutine qqh_1d_spline(xp,ax,axh,ix,ihx)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),axh(0:2)
   integer,intent(inout) :: ix,ihx
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   ihx=int(xx)
   sx=xx-real(ihx,dp)
   sx2=sx*sx
   axh(1)=0.75-sx2
   axh(2)=0.5*(0.25+sx2+sx)
   axh(0)=1.-axh(1)-axh(2)
   ix=ix-1
   ihx=ihx-1

  end subroutine qqh_1d_spline
!=======================
  subroutine qqh_2d_spline(xp,ax,axh,ay,ayh,ix,ihx,iy,ihy)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),axh(0:2),ay(0:2),ayh(0:2)
   integer,intent(inout) :: ix,ihx,iy,ihy
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   ihx=int(xx)
   sx=xx-real(ihx,dp)-0.5
   sx2=sx*sx
   axh(1)=0.75-sx2
   axh(2)=0.5*(0.25+sx2+sx)
   axh(0)=1.-axh(1)-axh(2)

   xx=shy+xp(2)
   iy=int(xx+0.5)
   sx=xx-real(iy,dp)
   sx2=sx*sx
   ay(1)=0.75-sx2
   ay(2)=0.5*(0.25+sx2+sx)
   ay(0)=1.-ay(1)-ay(2)

   ihy=int(xx)
   sx=xx-real(ihy,dp)-0.5
   sx2=sx*sx
   ayh(1)=0.75-sx2
   ayh(2)=0.5*(0.25+sx2+sx)
   ayh(0)=1.-ayh(1)-ayh(2)

   ix=ix-1
   ihx=ihx-1
   iy=iy-1
   ihy=ihy-1
  end subroutine qqh_2d_spline
!=======================================
  subroutine qlh_3d_spline(xp,ax,axh,ay,ayh,az,azh,ix,ihx,iy,ihy,iz,ihz)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),axh(0:1),ay(0:2),ayh(0:1),az(0:2),azh(0:1)
   integer,intent(inout) :: ix,ihx,iy,ihy,iz,ihz
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   axh(1)=sx+0.5
   axh(0)=1.-axh(1)

   xx=shy+xp(2)
   iy=int(xx+0.5)
   sx=xx-real(iy,dp)
   sx2=sx*sx
   ay(1)=0.75-sx2
   ay(2)=0.5*(0.25+sx2+sx)
   ay(0)=1.-ay(1)-ay(2)

   ayh(1)=sx+0.5
   ayh(0)=1.-ayh(1)

   xx=shz+xp(3)
   iz=int(xx+0.5)
   sx=xx-real(iz,dp)
   sx2=sx*sx
   az(1)=0.75-sx2
   az(2)=0.5*(0.25+sx2+sx)
   az(0)=1.-az(1)-az(2)

   azh(1)=sx+0.5
   azh(0)=1.-azh(1)

   ix=ix-1
   ihx=ix
   iy=iy-1
   ihy=iy
   iz=iz-1
   ihz=iz
  end subroutine qlh_3d_spline
!=================================
  subroutine qqh_3d_spline(xp,ax,axh,ay,ayh,az,azh,ix,ihx,iy,ihy,iz,ihz)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),axh(0:2),ay(0:2),ayh(0:2),az(0:2),azh(0:2)
   integer,intent(inout) :: ix,ihx,iy,ihy,iz,ihz
   real(sp) :: xx,sx,sx2

   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   ihx=int(xx)
   sx=xx-real(ihx,dp)-0.5
   sx2=sx*sx
   axh(1)=0.75-sx2
   axh(2)=0.5*(0.25+sx2+sx)
   axh(0)=1.-axh(1)-axh(2)

   xx=shy+xp(2)
   iy=int(xx+0.5)
   sx=xx-real(iy,dp)
   sx2=sx*sx
   ay(1)=0.75-sx2
   ay(2)=0.5*(0.25+sx2+sx)
   ay(0)=1.-ay(1)-ay(2)

   ihy=int(xx)
   sx=xx-real(ihy,dp)-0.5
   sx2=sx*sx
   ayh(1)=0.75-sx2
   ayh(2)=0.5*(0.25+sx2+sx)
   ayh(0)=1.-ayh(1)-ayh(2)

   xx=shz+xp(3)
   iz=int(xx+0.5)
   sx=xx-real(iz,dp)
   sx2=sx*sx
   az(1)=0.75-sx2
   az(2)=0.5*(0.25+sx2+sx)
   az(0)=1.-az(1)-az(2)

   ihz=int(xx)
   sx=xx-real(ihz,dp)-0.5
   sx2=sx*sx
   azh(1)=0.75-sx2
   azh(2)=0.5*(0.25+sx2+sx)
   azh(0)=1.-azh(1)-azh(2)

   ix=ix-1
   ihx=ihx-1
   iy=iy-1
   ihy=ihy-1
   iz=iz-1
   ihz=ihz-1
  end subroutine qqh_3d_spline

  subroutine qden_1d_wgh(xp,ax,ix)
   real(dp),intent(in) :: xp(:)
   integer,intent(inout) :: ix
   real(dp),intent(inout) :: ax(0:2)
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)
  end subroutine qden_1d_wgh

  subroutine qden_2d_wgh(xp,ax,ay,ix,iy)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),ay(0:2)
   integer,intent(inout) :: ix,iy
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   xx=shy+xp(2)
   iy=int(xx+0.5)
   sx=xx-real(iy,dp)
   sx2=sx*sx
   ay(1)=0.75-sx2
   ay(2)=0.5*(0.25+sx2+sx)
   ay(0)=1.-ay(1)-ay(2)
  end subroutine qden_2d_wgh
!------------------------------
  subroutine qden_3d_wgh(xp,ax,ay,az,ix,iy,iz)
   real(dp),intent(in) :: xp(:)
   real(dp),intent(inout) :: ax(0:2),ay(0:2),az(0:2)
   integer,intent(inout) :: ix,iy,iz
   real(sp) :: xx,sx,sx2
!======================
   xx=shx+xp(1)
   ix=int(xx+0.5)
   sx=xx-real(ix,dp)
   sx2=sx*sx
   ax(1)=0.75-sx2
   ax(2)=0.5*(0.25+sx2+sx)
   ax(0)=1.-ax(1)-ax(2)

   xx=shy+xp(2)
   iy=int(xx+0.5)
   sx=xx-real(iy,dp)
   sx2=sx*sx
   ay(1)=0.75-sx2
   ay(2)=0.5*(0.25+sx2+sx)
   ay(0)=1.-ay(1)-ay(2)

   xx=shz+xp(3)
   iz=int(xx+0.5)
   sx=xx-real(iz,dp)
   sx2=sx*sx
   az(1)=0.75-sx2
   az(2)=0.5*(0.25+sx2+sx)
   az(0)=1.-az(1)-az(2)
  end subroutine qden_3d_wgh
!====================================

 !DIR$ ATTRIBUTES INLINE :: ql_interpolate
 subroutine set_local_2d_positions(pt_loc,n1,np)

 real(dp),intent(inout) :: pt_loc(:,:)
 integer,intent(in) :: n1,np
 integer :: n
 !=========================
 do n=1,np
  pt_loc(n,1)=dx_inv*(pt_loc(n,1)-xmn)
 end do
 if(n1==0)return
 if(n_str==0)then
  do n=1,np
    pt_loc(n,2)=dy_inv*(pt_loc(n,2)-ymn)              !
  end do
 else
  call map2dy_part_sind(np,n_str,2,ymn,pt_loc)
 endif
 if(n1==1)return

 do n=1,np
  pt_loc(n,3)=dx_inv*(pt_loc(n,3)-xmn)
 end do
 if(n_str==0)then
  do n=1,np
    pt_loc(n,4)=dy_inv*(pt_loc(n,4)-ymn)              !
  end do
 else
  call map2dy_part_sind(np,n_str,4,ymn,pt_loc)
 endif

 end subroutine set_local_2d_positions
!======================
 subroutine set_local_3d_positions(pt_loc,n1,np)
  real(dp),intent(inout) :: pt_loc(:,:)
  integer,intent(in) :: n1,np
  integer :: n
 !=========================
  do n=1,np
   pt_loc(n,1)=dx_inv*(pt_loc(n,1)-xmn)
  end do
  if(n_str==0)then
   do n=1,np
    pt_loc(n,2)=dy_inv*(pt_loc(n,2)-ymn)
    pt_loc(n,3)=dz_inv*(pt_loc(n,3)-zmn)
   end do
  else
   call map3d_part_sind(pt_loc,np,n_str,2,3,ymn,zmn)
  endif
  if(n1==1)return
  do n=1,np
   pt_loc(n,4)=dx_inv*(pt_loc(n,4)-xmn)
  end do
  if(n_str==0)then
   do n=1,np
    pt_loc(n,5)=dy_inv*(pt_loc(n,5)-ymn)
    pt_loc(n,6)=dz_inv*(pt_loc(n,6)-zmn)
   end do
  else
   call map3d_part_sind(pt_loc,np,n_str,5,6,ymn,zmn)
  endif

 end subroutine set_local_3d_positions
 
 end module grid_part_lib
 !==========================
