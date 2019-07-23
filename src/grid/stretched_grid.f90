!*****************************************************************************************************!
!                            Copyright 2008-2019  The ALaDyn Collaboration                            !
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

 module stretched_grid

  use common_param
  use grid_param

  implicit none
  private
  real(dp) :: const, smin, smax, xs
  real(dp), parameter :: SYMM_CENTER = zero_dp
  
  public :: map2dy_part_sind, map3d_part_sind
  
 contains
  !============== Mapping for stretched grids==========

  pure function invert_stretched_grid(yp_in) result(stretched)
   real(dp), intent(in) :: yp_in
   real(dp) :: stretched
   stretched = dyi_inv*atan(sy_rat*(yp_in+const-xs))+ny_stretch
  end function

  pure function invert_uniform_grid(yp_in) result(uniform)
   real(dp), intent(in) :: yp_in
   real(dp) :: uniform
   uniform = (yp_in + const)*dy_inv
  end function

  subroutine map2dy_part_sind(np, sind, ic1, ym, pt)
   integer, intent (in) :: np, sind, ic1
   real (dp), intent (in) :: ym
   real (dp), intent (inout) :: pt(:, :)
   real (dp) :: yp, yp_loc
   integer :: n
   !========================
   !  enter the y=part(ic1,n) particle position in stretched grid
   !            y=y(xi)
   !  exit      xi=part(ic1,n) the  particle position in uniform grid
   !               normalized to the Dxi cell size
   !==========================================
   const = one_dp*ny*dy/2
   smin = str_ygrid%smin
   smax = str_ygrid%smax
   xs = ny_stretch*dy
   do n = 1, np
    yp = pt(n, ic1)
    if (yp <= smin) then
     yp_loc = invert_stretched_grid(yp)
    else if (yp >= smax) then
     yp = 2*SYMM_CENTER - yp
     yp_loc = invert_stretched_grid(yp)
     yp_loc = ny - yp_loc
    else
     yp_loc = invert_uniform_grid(yp)
    end if
    pt(n, ic1) = yp_loc
   end do
  end subroutine

  subroutine map2dy_part_sind_old(np, sind, ic1, ym, pt)
   integer, intent (in) :: np, sind, ic1
   real (dp), intent (in) :: ym
   real (dp), intent (inout) :: pt(:, :)
   real (dp) :: yp, ys, ximn, yp_loc
   integer :: n
   !========================
   !  enter the y=part(ic1,n) particle position in stretched grid
   !            y=y(xi)
   !  exit      xi=part(ic1,n) the  particle position in uniform grid
   !               normalized to the Dxi cell size
   !==========================================
   const = one_dp*ny*dy/2
   smin = str_ygrid%smin
   smax = str_ygrid%smax
   xs = ny_stretch*dy
   select case (sind)
   case (1) !y<0
    ys = str_ygrid%smin
    ximn = -dyi_inv*atan(sy_rat*(ym-ys))
    do n = 1, np
     yp = pt(n, ic1)
     yp_loc = ximn + dy_inv*(yp-ys)
     if (yp<=ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     pt(n, ic1) = yp_loc
    end do
   case (2) !y>0
    ys = str_ygrid%smax
    if (ym>ys) then
     ximn = dyi_inv*atan(sy_rat*(ys-ym))
    else
     ximn = dy_inv*(ys-ym)
    end if
    do n = 1, np
     yp = pt(n, ic1)
     yp_loc = dy_inv*(yp-ym)
     if (yp>ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     pt(n, ic1) = yp_loc
    end do
   end select
  end subroutine

  subroutine map2dz_part_sind(np, sind, ic1, zm, pt)
   integer, intent (in) :: np, sind, ic1
   real (dp), intent (in) :: zm
   real (dp), intent (inout) :: pt(:, :)
   real (dp) :: zp, zs, zimn, zp_loc
   integer :: n
   !========================
   !  enter the y=part(ic1,n) particle position in stretched grid
   !            y=y(xi)
   !  exit      xi=part(ic1,n) the  particle position in uniform grid
   !               normalized to the Dxi cell size
   !==========================================
   const = one_dp*ny*dy/2
   smin = str_ygrid%smin
   smax = str_ygrid%smax
   xs = ny_stretch*dy
   select case (sind)
   case (1) !z<0
    zs = str_zgrid%smin
    zimn = -dzi_inv*atan(sz_rat*(zm-zs))
    do n = 1, np
     zp = pt(n, ic1)
     zp_loc = zimn + dz_inv*(zp-zs)
     if (zp<=zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic1) = zp_loc
    end do
   case (2) !z>0
    zs = str_zgrid%smax
    if (zm>zs) then
     zimn = dzi_inv*atan(sz_rat*(zs-zm))
    else
     zimn = dz_inv*(zs-zm)
    end if
    do n = 1, np
     zp = pt(n, ic1)
     zp_loc = dz_inv*(zp-zm)
     if (zp>zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic1) = zp_loc
    end do
   end select
  end subroutine

  subroutine map3d_part_sind(pt, np, sind, ic1, ic2, ym, zm)
   integer, intent (in) :: np, sind, ic1, ic2
   real (dp), intent (in) :: ym, zm
   real (dp), intent (inout) :: pt(:, :)
   real (dp) :: yp, zp, yp_loc, zp_loc, ys, zs, ximn, zimn
   integer :: n
   !========================
   !  enter the y=part(n,ic1) z=part(n,ic2) particle positions
   !        in stretched grids    y=y(xi), z(zi)
   !  exit   xi=part(n,ic1) zi=part(n,ic2)
   !    particle positions in uniform grid
   !    normalized to the (Dxi Dzi) cell sizes
   !==========================================
   const = one_dp*ny*dy/2
   smin = str_ygrid%smin
   smax = str_ygrid%smax
   xs = ny_stretch*dy
   select case (sind)
   case (1) !y<0 z<0 corner ys>ymn
    ys = str_ygrid%smin
    ximn = -dyi_inv*atan(sy_rat*(ym-ys))
    zs = str_zgrid%smin
    zimn = -dzi_inv*atan(sz_rat*(zm-zs))
    do n = 1, np
     yp = pt(n, ic1)
     zp = pt(n, ic2)
     yp_loc = ximn + dy_inv*(yp-ys)
     zp_loc = zimn + dz_inv*(zp-zs)
     if (yp<ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     if (zp<zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic1) = yp_loc
     pt(n, ic2) = zp_loc
    end do
   case (2) !z<0
    zs = str_zgrid%smin
    zimn = -dzi_inv*atan(sz_rat*(zm-zs))
    do n = 1, np
     yp = pt(n, ic1)
     pt(n, ic1) = dy_inv*(yp-ym)
     zp = pt(n, ic2)
     zp_loc = zimn + dz_inv*(zp-zs)
     if (zp<zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic2) = zp_loc
    end do
   case (3) !y>0 z<0 corner
    ys = str_ygrid%smax
    if (ym>ys) then
     ximn = dyi_inv*atan(sy_rat*(ys-ym))
    else
     ximn = dy_inv*(ys-ym)
    end if
    zs = str_zgrid%smin
    zimn = -dzi_inv*atan(sz_rat*(zm-zs))
    do n = 1, np
     yp = pt(n, ic1)
     zp = pt(n, ic2)
     yp_loc = dy_inv*(yp-ym)
     zp_loc = zimn + dz_inv*(zp-zs)
     if (yp>ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     if (zp<zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic1) = yp_loc
     pt(n, ic2) = zp_loc
    end do
   case (4) !y>0
    ys = str_ygrid%smax
    if (ym>ys) then
     ximn = dyi_inv*atan(sy_rat*(ys-ym))
    else
     ximn = dy_inv*(ys-ym)
    end if
    do n = 1, np
     yp = pt(n, ic1)
     zp = pt(n, ic2)
     yp_loc = dy_inv*(yp-ym)
     if (yp>ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     pt(n, ic1) = yp_loc
     pt(n, ic2) = dz_inv*(zp-zm)
    end do
   case (5) !y>0 z>0 corner
    ys = str_ygrid%smax
    if (ym>ys) then
     ximn = dyi_inv*atan(sy_rat*(ys-ym))
    else
     ximn = dy_inv*(ys-ym)
    end if
    zs = str_zgrid%smax
    if (zm>zs) then
     zimn = dzi_inv*atan(sz_rat*(zs-zm))
    else
     zimn = dz_inv*(zs-zm)
    end if
    do n = 1, np
     yp = pt(n, ic1)
     zp = pt(n, ic2)
     yp_loc = dy_inv*(yp-ym)
     zp_loc = dz_inv*(zp-zm)
     if (yp>ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     if (zp>zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic1) = yp_loc
     pt(n, ic2) = zp_loc
    end do
   case (6) !z>0
    zs = str_zgrid%smax
    if (zm>zs) then
     zimn = dzi_inv*atan(sz_rat*(zs-zm))
    else
     zimn = dz_inv*(zs-zm)
    end if
    do n = 1, np
     yp = pt(n, ic1)
     pt(n, ic1) = dy_inv*(yp-ym)
     zp = pt(n, ic2)
     zp_loc = dz_inv*(zp-zm)
     if (zp>zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic2) = zp_loc
    end do
   case (7) !y<0 z>0 corner
    ys = str_ygrid%smin
    ximn = -dyi_inv*atan(sy_rat*(ym-ys))
    zs = str_zgrid%smax
    if (zm>zs) then
     zimn = dzi_inv*atan(sz_rat*(zs-zm))
    else
     zimn = dz_inv*(zs-zm)
    end if
    do n = 1, np
     yp = pt(n, ic1)
     zp = pt(n, ic2)
     yp_loc = ximn + dy_inv*(yp-ys)
     zp_loc = dz_inv*(zp-zm)
     if (yp<ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     if (zp>zs) zp_loc = zimn + dzi_inv*atan(sz_rat*(zp-zs))
     pt(n, ic1) = yp_loc
     pt(n, ic2) = zp_loc
    end do
   case (8) !y<0
    ys = str_ygrid%smin
    ximn = -dyi_inv*atan(sy_rat*(ym-ys))
    do n = 1, np
     yp = pt(n, ic1)
     zp = pt(n, ic2)
     yp_loc = ximn + dy_inv*(yp-ys)
     if (yp<ys) yp_loc = ximn + dyi_inv*atan(sy_rat*(yp-ys))
     pt(n, ic1) = yp_loc
     pt(n, ic2) = dz_inv*(zp-zm)
    end do
   end select
  end subroutine
  !========================================
 end module