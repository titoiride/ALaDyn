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

 use grid_param

 implicit none

 real(dp),parameter :: half=0.5,thr=0.75,two_third=2./3.,one_sixth=1./6.
 real(sp),parameter :: shx=3.,shy=3.,shz=3.
 integer(kind=2) :: err_ind

 contains
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
 !    Linear and quadratic shape functions for grid-particle  connection
 !==================================================
 subroutine set_hint_pshape(spl,xx,ax,ind)
 integer,intent(in) :: spl
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: ax(0:3)
 integer,intent(out) :: ind
 real(dp) :: sx,sx2,sx3
 !To half-integer grid points
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
 subroutine set_int_pshape(spl,xx,ax,ind)
 integer,intent(in) :: spl
 real(dp),intent(in) :: xx
 real(dp),intent(out) :: ax(0:3)
 integer,intent(out) :: ind
 real(dp) :: sx,sx2,sx3
 !To integer grid points
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
 !===========================================
 !     Computes particle density charge on grid integer points
 !==============================================
 !DIR$ ATTRIBUTES INLINE :: qq_interpolate
!====================
 subroutine qq_density_spline(sh,x0,x1,a0,a1,ix0,ix1)
 real(sp),intent(in)  :: sh
 real(dp),intent(in)  :: x0,x1
 real(dp),intent(inout) :: a0(3),a1(3)
 integer, intent(out) :: ix0,ix1
 real(dp) :: xx,sx,sx2

 xx=sh+x0
 ix0=int(xx+0.5)
 sx=xx-real(ix0,dp)
 sx2=sx*sx
 a0(2)=0.75-sx2
 a0(3)=0.5*(0.25+sx2+sx)
 a0(1)=1.-a0(2)-a0(3)
 xx=sh+x1
 ix1=int(xx+0.5)
 sx=xx-real(ix1,dp)
 sx2=sx*sx
 a1(2)=0.75-sx2
 a1(3)=0.5*(0.25+sx2+sx)
 a1(1)=1.-a1(2)-a1(3)
 ix0=ix0-2
 ix1=ix1-2
 end subroutine qq_density_spline
!==========================================
 subroutine qlql_density_spline(sh,x0,x1,a0,a1,ah0,ah1,ix0,ix1,ih0,ih1)
 real(sp),intent(in)  :: sh
 real(dp),intent(in)  :: x0,x1
 real(dp),intent(inout) :: a0(3),a1(3),ah0(2),ah1(2)
 integer, intent(out) :: ix0,ix1,ih0,ih1
 real(dp) :: xx,sx,sx2

 xx=sh+x0
 ix0=int(xx+0.5)
 sx=xx-real(ix0,dp)
 sx2=sx*sx
 a0(2)=0.75-sx2
 a0(3)=0.5*(0.25+sx2+sx)
 a0(1)=1.-a0(2)-a0(3)

 ah0(2)=sx+0.5
 ah0(1)=1.-ah0(2)

 xx=sh+x1
 ix1=int(xx+0.5)
 sx=xx-real(ix1,dp)
 sx2=sx*sx
 a1(2)=0.75-sx2
 a1(3)=0.5*(0.25+sx2+sx)
 a1(1)=1.-a1(2)-a1(3)

 ah1(2)=sx+0.5
 ah1(1)=1.-ah1(2)

 ix0=ix0-2
 ih0=ix0
 ix1=ix1-2
 ih1=ix1
 end subroutine qlql_density_spline
!==========================================
!
 subroutine qq_interpolate(xp,ax1,axh,ay1,ayh,ix,iy,ixh,iyh)
 real(dp),intent(in)  :: xp(:)
 real(dp),intent(inout) :: ax1(:),axh(:),ay1(:),ayh(:)
 integer, intent(out) :: ix,iy,ixh,iyh
 real(dp) :: xx,sx,sx2

 xx=shx+xp(1)
 ix=int(xx+0.5)
 sx=xx-real(ix,dp)
 sx2=sx*sx
 ax1(2)=0.75-sx2
 ax1(3)=0.5*(0.25+sx2+sx)
 ax1(1)=1.-ax1(2)-ax1(3)

 ixh=int(xx)
 sx=xx-0.5-real(ixh,dp)
 sx2=sx*sx
 axh(2)=0.75-sx2
 axh(3)=0.5*(0.25+sx2+sx)
 axh(1)=1.-axh(3)-axh(2)

 xx=shy+xp(2)

 iy=int(xx+0.5)
 sx=xx-real(iy,dp)
 sx2=sx*sx
 ay1(2)=0.75-sx2
 ay1(3)=0.5*(0.25+sx2+sx)
 ay1(1)=1.-ay1(3)-ay1(2)

 iyh=int(xx)
 sx=xx-0.5-real(iyh,dp)
 sx2=sx*sx
 ayh(2)=0.75-sx2
 ayh(3)=0.5*(0.25+sx2+sx)
 ayh(1)=1.-ayh(3)-ayh(2)

 ix=ix-2
 iy=iy-2
 ixh=ixh-2
 iyh=iyh-2
 end subroutine qq_interpolate

 !DIR$ ATTRIBUTES INLINE :: ql_interpolate
 subroutine ql_interpolate(xp,ax1,axh,ay1,ayh,ix,iy,ixh,iyh)
 real(dp),intent(in)  :: xp(:)
 real(dp),intent(inout) :: ax1(:),axh(:),ay1(:),ayh(:)
 integer, intent(out) :: ix,iy,ixh,iyh
 real(dp) :: xx,sx,sx2

 xx=shx+xp(1)
 ix=int(xx+0.5)
 sx=xx-real(ix,dp)
 sx2=sx*sx
 ax1(2)=0.75-sx2
 ax1(3)=0.5*(0.25+sx2+sx)
 ax1(1)=1.-ax1(2)-ax1(3)

 axh(2)=sx+0.5
 axh(1)=1.-axh(2)

 xx=shy+xp(2)

 iy=int(xx+0.5)
 sx=xx-real(iy,dp)
 sx2=sx*sx
 ay1(2)=0.75-sx2
 ay1(3)=0.5*(0.25+sx2+sx)
 ay1(1)=1.-ay1(3)-ay1(2)

 ayh(2)=sx+0.5
 ayh(1)=1.-ayh(2)

 ix=ix-2
 iy=iy-2
 ixh=ix
 iyh=iy
 end subroutine ql_interpolate
 
 end module grid_part_lib
 !==========================
