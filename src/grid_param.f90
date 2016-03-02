!*****************************************************************************************************!
!             Copyright 2008-2016 Pasquale Londrillo, Stefano Sinigardi, Andrea Sgattoni              !
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

 module grid_param

 use precision_def
 use struct_def

 implicit none

 type(grid),allocatable :: loc_rgrid(:),loc_ygrid(:),loc_zgrid(:),loc_xgrid(:)
 real(dp),allocatable :: loc_yg(:,:,:),loc_zg(:,:,:),loc_xg(:,:,:)
 real(dp),allocatable :: x(:),xw(:),y(:),z(:),dx1(:),dy1(:),dz1(:)
 real(dp),allocatable :: xh(:),yh(:),zh(:),dx1h(:),dy1h(:),dz1h(:)
 real(dp),allocatable :: loc_rg(:,:,:),r(:),rh(:),dr1(:),dr1h(:),dvr(:)
 integer,allocatable :: str_indx(:,:)
 real(dp),allocatable :: rpt(:),wgp(:)
 real(dp) :: xtot,xmax,xmin,ymax,ymin,zmax,zmin
 real(dp) :: dx,dx_inv,dxi_inv,dy,dz,dy_inv,dyi_inv,dz_inv,dzi_inv
 real(dp) :: djc(3)
 real(dp) :: aph,L_s,Lx_s,dxi,dyi,dzi,sy_rat,sz_rat,sx_rat
 type(sgrid) :: str_xgrid,str_ygrid,str_zgrid
 integer :: loc_ygr_max,loc_zgr_max,loc_xgr_max
 integer,parameter :: sh_ix=3

 !--------------------------

 contains

 !--------------------------

 subroutine grid_alloc(n1,n1_loc,n2,n2_loc,n3,n3_loc,&
  n2_targ,pml_sz,npey,npez,npex)
 integer,intent(in) :: n1,n1_loc,n2,n2_loc,n3,n3_loc,&
  n2_targ,pml_sz,npey,npez,npex
 integer :: n_loc,nyt_loc,nyt_res,nzt_loc,nzt_res,nyt_hres,nzt_hres
 integer :: p

 allocate(loc_ygrid(0:npey-1),loc_zgrid(0:npez-1))
 allocate(loc_rgrid(0:npey-1),loc_xgrid(0:npex-1))
 allocate(x(n1+1),xw(n1+1),dx1(n1+1),y(n2+1),z(n3+1),dy1(n2+1),dz1(n3+1))
 allocate(dx1h(n1+1),dy1h(n2+1),dz1h(n3+1))
 allocate(xh(n1+1),yh(n2+1),zh(n3+1))
 allocate(r(n2+1),rh(n2+1),dr1(n2+1),dr1h(n2+1),dvr(0:n2+1))

 loc_xgr_max=n1_loc
 do p=0,npex-1
  loc_xgrid(p)%ng=n1_loc
 end do
 if(pml_sz==0)then
  loc_ygr_max=n2_loc
  loc_zgr_max=n3_loc
  do p=0,npey-1
   loc_ygrid(p)%ng=n2_loc
   loc_rgrid(p)%ng=n2_loc
  end do
  do p=0,npez-1
   loc_zgrid(p)%ng=n3_loc
  end do
 else
  n_loc=1
  if(npey>1)n_loc=pml_sz
  nyt_loc=0
  if(npey>2)then
   nyt_loc=n2_targ/(npey-2)
   loc_ygrid(1:npey-2)%ng=nyt_loc
   nyt_res=n2_targ-nyt_loc*(npey-2)
   if(nyt_res >0)then
    nyt_hres=nyt_res/2
    n_loc=n_loc+nyt_hres
   endif
  endif
  loc_ygrid(0)%ng=n_loc
  loc_ygrid(npey-1)%ng=n_loc
  loc_ygr_max=max(n_loc,nyt_loc)

  n_loc=1
  if(npez>1)n_loc=pml_sz
  nzt_loc=0
  if(npez>2)then
   nzt_loc=n2_targ/(npez-2)
   loc_zgrid(1:npez-2)%ng=nzt_loc
   nzt_res=n2_targ-nzt_loc*(npez-2)
   if(nzt_res >0)then
    nzt_hres=nzt_res/2
    n_loc=n_loc+nzt_hres
   endif
  endif
  loc_zgrid(0)%ng=n_loc
  loc_zgrid(npez-1)%ng=n_loc
  loc_zgr_max=max(n_loc,nzt_loc)
 endif

 allocate(loc_yg(0:loc_ygr_max+1,4,0:npey-1))
 allocate(loc_zg(0:loc_zgr_max+1,4,0:npez-1))
 allocate(loc_xg(0:loc_xgr_max+1,4,0:npex-1))
 allocate(loc_rg(0:loc_ygr_max+1,5,0:npey-1)) !FIX SERVE SE NON SIAMO IN CILINDRICO???

 allocate(str_indx(0:npey-1,0:npez-1))
 str_indx(0:npey-1,0:npez-1)=0

 end subroutine grid_alloc

 !--------------------------

 subroutine set_fyzxgrid(npey,npez,npex,sh)
 integer,intent(in) :: npey,npez,npex,sh
 integer :: i,ii,p,ip,n_loc

 ! Defines initial local p-grid coordinate and loc n_cell
 ! y-grid decomposed on n2_loc uniform grid size
 loc_ygrid(0)%gmin=y(1)
 loc_rgrid(0)%gmin=r(1)
 ip=loc_ygrid(0)%ng
 n_loc=ip
 loc_ygrid(0)%gmax=y(ip+1)
 loc_rgrid(0)%gmax=r(ip+1)
 loc_ygrid(0)%p_ind(1)=min(sh,n_loc)
 loc_rgrid(0)%p_ind(1)=min(sh,n_loc)
 loc_ygrid(0)%p_ind(2)=n_loc+loc_ygrid(0)%p_ind(1)-1
 loc_rgrid(0)%p_ind(2)=n_loc+loc_rgrid(0)%p_ind(1)-1

 p=0
 do i=1,n_loc+1
  loc_yg(i,1,p)=y(i)
  loc_yg(i,2,p)=yh(i)
  loc_yg(i,3,p)=dy1(i)
  loc_yg(i,4,p)=dy1h(i)

  loc_rg(i,1,p)=r(i) !FIX SERVE SE NON SIAMO IN CILINDRICO???
  loc_rg(i,2,p)=rh(i)
  loc_rg(i,3,p)=dr1(i)
  loc_rg(i,4,p)=dr1h(i)
  loc_rg(i,5,p)=dvr(i)
 end do

 if(npey >1)then
  ip=loc_ygrid(0)%ng
  if(npey >2) then
   do p=1,npey-2
    n_loc=loc_ygrid(p-1)%ng
    loc_yg(0,1:4,p)=loc_yg(n_loc,1:4,p-1)
    loc_rg(0,1:4,p)=loc_rg(n_loc,1:4,p-1) !FIX SERVE SE NON SIAMO IN CILINDRICO???
    n_loc=loc_ygrid(p)%ng
    do i=1,n_loc+1
     ii=i+ip
     loc_yg(i,1,p)=y(ii)
     loc_yg(i,2,p)=yh(ii)
     loc_yg(i,3,p)=dy1(ii)
     loc_yg(i,4,p)=dy1h(ii)
     loc_rg(i,1,p)=r(ii) !FIX SERVE SE NON SIAMO IN CILINDRICO???
     loc_rg(i,2,p)=rh(ii)
     loc_rg(i,3,p)=dr1(ii)
     loc_rg(i,4,p)=dr1h(ii)
     loc_rg(i,5,p)=dvr(ii)
    end do
    loc_ygrid(p)%gmin=loc_ygrid(p-1)%gmax
    loc_rgrid(p)%gmin=loc_rgrid(p-1)%gmax !FIX SERVE SE NON SIAMO IN CILINDRICO???

    ip=ip+n_loc
    loc_ygrid(p)%gmax=y(ip+1)
    loc_rgrid(p)%gmax=r(ip+1) !FIX SERVE SE NON SIAMO IN CILINDRICO???

    loc_ygrid(p)%p_ind(1)=sh
    loc_ygrid(p)%p_ind(2)=n_loc+loc_ygrid(p)%p_ind(1)-1

    loc_rgrid(p)%p_ind(1)=sh !FIX SERVE SE NON SIAMO IN CILINDRICO???
    loc_rgrid(p)%p_ind(2)=n_loc+loc_rgrid(p)%p_ind(1)-1
   end do
  endif
  p=npey-1
  n_loc=loc_ygrid(p-1)%ng
  loc_yg(0,1:4,p)=loc_yg(n_loc,1:4,p-1)
  loc_rg(0,1:4,p)=loc_rg(n_loc,1:4,p-1) !FIX SERVE SE NON SIAMO IN CILINDRICO???
  n_loc=loc_ygrid(p)%ng
  do i=1,n_loc+1
   ii=i+ip
   loc_yg(i,1,p)=y(ii)
   loc_yg(i,2,p)=yh(ii)
   loc_yg(i,3,p)=dy1(ii)
   loc_yg(i,4,p)=dy1h(ii)

   loc_rg(i,1,p)=r(ii) !FIX SERVE SE NON SIAMO IN CILINDRICO???
   loc_rg(i,2,p)=rh(ii)
   loc_rg(i,3,p)=dr1(ii)
   loc_rg(i,4,p)=dr1h(ii)
   loc_rg(i,5,p)=dvr(ii)
  end do
  loc_ygrid(p)%gmin=loc_ygrid(p-1)%gmax
  loc_rgrid(p)%gmin=loc_rgrid(p-1)%gmax
  ip=ip+n_loc
  loc_ygrid(p)%gmax=y(ip+1)
  loc_ygrid(p)%p_ind(1)=sh
  loc_ygrid(p)%p_ind(2)=n_loc+loc_ygrid(p)%p_ind(1)-1

  loc_rgrid(p)%gmax=r(ip+1)
  loc_rgrid(p)%p_ind(1)=sh
  loc_rgrid(p)%p_ind(2)=n_loc+loc_rgrid(p)%p_ind(1)-1
 endif
 !=========================
 loc_zgrid(0)%gmin=z(1)
 ip=loc_zgrid(0)%ng
 n_loc=ip
 loc_zgrid(0)%gmax=z(ip+1)
 loc_zgrid(0)%p_ind(1)=min(sh,n_loc)
 loc_zgrid(0)%p_ind(2)=n_loc+loc_zgrid(0)%p_ind(1)-1

 p=0
 do i=1,n_loc+1
  loc_zg(i,1,p)=z(i)
  loc_zg(i,2,p)=zh(i)
  loc_zg(i,3,p)=dz1(i)
  loc_zg(i,4,p)=dz1h(i)
 end do

 if(npez >1)then
  ip=loc_zgrid(0)%ng
  if(npez >2) then
   do p=1,npez-2
    n_loc=loc_zgrid(p-1)%ng
    loc_zg(0,1:4,p)=loc_zg(n_loc,1:4,p-1)
    n_loc=loc_zgrid(p)%ng
    do i=1,n_loc+1
     ii=i+ip
     loc_zg(i,1,p)=z(ii)
     loc_zg(i,2,p)=zh(ii)
     loc_zg(i,3,p)=dz1(ii)
     loc_zg(i,4,p)=dz1h(ii)
    end do
    loc_zgrid(p)%gmin=loc_zgrid(p-1)%gmax
    ip=ip+n_loc
    loc_zgrid(p)%gmax=z(ip+1)
    loc_zgrid(p)%p_ind(1)=sh
    loc_zgrid(p)%p_ind(2)=n_loc+loc_zgrid(p)%p_ind(1)-1
   end do
  endif
  p=npez-1
  n_loc=loc_zgrid(p-1)%ng
  loc_zg(0,1:4,p)=loc_zg(n_loc,1:4,p-1)
  n_loc=loc_zgrid(p)%ng
  do i=1,n_loc+1
   ii=i+ip
   loc_zg(i,1,p)=z(ii)
   loc_zg(i,2,p)=zh(ii)
   loc_zg(i,3,p)=dz1(ii)
   loc_zg(i,4,p)=dz1h(ii)
  end do
  loc_zgrid(p)%gmin=loc_zgrid(p-1)%gmax
  ip=ip+n_loc
  loc_zgrid(p)%gmax=z(ip+1)
  loc_zgrid(p)%p_ind(1)=sh
  loc_zgrid(p)%p_ind(2)=n_loc+loc_zgrid(p)%p_ind(1)-1
 endif
 !======================
 loc_xgrid(0)%gmin=x(1)
 ip=loc_xgrid(0)%ng
 n_loc=ip
 loc_xgrid(0)%gmax=x(ip+1)
 loc_xgrid(0)%p_ind(1)=min(sh,n_loc)
 loc_xgrid(0)%p_ind(2)=n_loc+loc_xgrid(0)%p_ind(1)-1

 p=0
 do i=1,n_loc+1
  loc_xg(i,1,p)=x(i)
  loc_xg(i,2,p)=xh(i)
  loc_xg(i,3,p)=dx1(i)
  loc_xg(i,4,p)=dx1h(i)
 end do

 if(npex >1)then
  ip=loc_xgrid(0)%ng
  if(npex >2) then
   do p=1,npex-2
    n_loc=loc_xgrid(p-1)%ng
    loc_xg(0,1:4,p)=loc_xg(n_loc,1:4,p-1)
    n_loc=loc_xgrid(p)%ng
    do i=1,n_loc+1
     ii=i+ip
     loc_xg(i,1,p)=x(ii)
     loc_xg(i,2,p)=xh(ii)
     loc_xg(i,3,p)=dx1(ii)
     loc_xg(i,4,p)=dx1h(ii)
    end do
    loc_xgrid(p)%gmin=loc_xgrid(p-1)%gmax
    ip=ip+n_loc
    loc_xgrid(p)%gmax=x(ip+1)
    loc_xgrid(p)%p_ind(1)=sh
    loc_xgrid(p)%p_ind(2)=n_loc+loc_xgrid(p)%p_ind(1)-1
   end do
  endif
  p=npex-1
  n_loc=loc_xgrid(p-1)%ng
  loc_xg(0,1:4,p)=loc_xg(n_loc,1:4,p-1)
  n_loc=loc_xgrid(p)%ng
  do i=1,n_loc+1
   ii=i+ip
   loc_xg(i,1,p)=x(ii)
   loc_xg(i,2,p)=xh(ii)
   loc_xg(i,3,p)=dx1(ii)
   loc_xg(i,4,p)=dx1h(ii)
  end do
  loc_xgrid(p)%gmin=loc_xgrid(p-1)%gmax
  ip=ip+n_loc
  loc_xgrid(p)%gmax=x(ip+1)
  loc_xgrid(p)%p_ind(1)=sh
  loc_xgrid(p)%p_ind(2)=n_loc+loc_xgrid(p)%p_ind(1)-1
 endif
 end subroutine set_fyzxgrid
 !======================
 subroutine set_fxgrid(npex,sh)
 integer,intent(in) :: npex,sh
 integer :: i,ii,p,ip,n_loc

 loc_xgrid(0)%gmin=x(1)
 ip=loc_xgrid(0)%ng
 n_loc=ip
 loc_xgrid(0)%gmax=x(ip+1)
 loc_xgrid(0)%p_ind(1)=min(sh,n_loc)
 loc_xgrid(0)%p_ind(2)=n_loc+loc_xgrid(0)%p_ind(1)-1

 p=0
 do i=1,n_loc+1
  loc_xg(i,1,p)=x(i)
  loc_xg(i,2,p)=xh(i)
  loc_xg(i,3,p)=dx1(i)
  loc_xg(i,4,p)=dx1h(i)
 end do

 if(npex >1)then
  ip=loc_xgrid(0)%ng
  if(npex >2) then
   do p=1,npex-2
    n_loc=loc_xgrid(p-1)%ng
    loc_xg(0,1:4,p)=loc_xg(n_loc,1:4,p-1)
    n_loc=loc_xgrid(p)%ng
    do i=1,n_loc+1
     ii=i+ip
     loc_xg(i,1,p)=x(ii)
     loc_xg(i,2,p)=xh(ii)
     loc_xg(i,3,p)=dx1(ii)
     loc_xg(i,4,p)=dx1h(ii)
    end do
    loc_xgrid(p)%gmin=loc_xgrid(p-1)%gmax
    ip=ip+n_loc
    loc_xgrid(p)%gmax=x(ip+1)
    loc_xgrid(p)%p_ind(1)=sh
    loc_xgrid(p)%p_ind(2)=n_loc+loc_xgrid(p)%p_ind(1)-1
   end do
  endif
  p=npex-1
  n_loc=loc_xgrid(p-1)%ng
  loc_xg(0,1:4,p)=loc_xg(n_loc,1:4,p-1)
  n_loc=loc_xgrid(p)%ng
  do i=1,n_loc+1
   ii=i+ip
   loc_xg(i,1,p)=x(ii)
   loc_xg(i,2,p)=xh(ii)
   loc_xg(i,3,p)=dx1(ii)
   loc_xg(i,4,p)=dx1h(ii)
  end do
  loc_xgrid(p)%gmin=loc_xgrid(p-1)%gmax
  ip=ip+n_loc
  loc_xgrid(p)%gmax=x(ip+1)
  loc_xgrid(p)%p_ind(1)=sh
  loc_xgrid(p)%p_ind(2)=n_loc+loc_xgrid(p)%p_ind(1)-1
 endif
 end subroutine set_fxgrid

 subroutine set_str_ind(npey,npez,ndm)
 integer,intent(in) :: npey,npez,ndm
 integer :: p,q,ip(4)

 str_indx(0:npey-1,0:npez-1)=0
 ip=0
 if(ndm <3)then
  do p=0,npey-1
   if(str_ygrid%smin >loc_ygrid(p)%gmin)ip(1)=p
   if(str_ygrid%smax >= loc_ygrid(p)%gmin)ip(2)=p
  end do

  p=0
  do q=0,ip(1)
   str_indx(q,p)=1
  end do
  do q=ip(2),npey-1
   str_indx(q,p)=2
  end do
  return
 endif
 do p=0,npey-1
  if(str_ygrid%smin > loc_ygrid(p)%gmin)ip(1)=p
  if(str_ygrid%smax >= loc_ygrid(p)%gmin)ip(2)=p
 end do
 do p=0,npez-1
  if(str_zgrid%smin >loc_zgrid(p)%gmin)ip(3)=p
  if(str_zgrid%smax >=loc_zgrid(p)%gmin)ip(4)=p
 end do

 do p=0,ip(3)
  str_indx(0:npey-1,p)=2
  do q=0,ip(1)
   str_indx(q,p)=1
  end do
  do q=ip(2),npey-1
   str_indx(q,p)=3
  end do
 end do
 do p=ip(3)+1,ip(4)-1
  do q=0,ip(1)
   str_indx(q,p)=8
  end do
  do q=ip(2),npey-1
   str_indx(q,p)=4
  end do
 end do
 do p=ip(4),npez-1
  str_indx(0:npey-1,p)=6
  do q=0,ip(1)
   str_indx(q,p)=7
  end do
  do q=ip(2),npey-1
   str_indx(q,p)=5
  end do
 end do
 end subroutine set_str_ind

 !--------------------------

 subroutine map2dy_part_sind(np,sind,ic1,ym,pt)
 integer,intent(in) :: np,sind,ic1
 real(dp),intent(in) :: ym
 real(dp),intent(inout) :: pt(:,:)
 real(dp) :: yp,ys,ximn
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
   yp=pt(ic1,n)
   pt(ic1,n)=ximn+dy_inv*(yp-ys)
   if(yp <= ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
  end do
 case(2)       !y>0
  ys=str_ygrid%smax
  if(ym>ys)then
   ximn=dyi_inv*atan(sy_rat*(ys-ym))
  else
   ximn=dy_inv*(ys-ym)
  endif
  do n=1,np
   yp=pt(ic1,n)
   pt(ic1,n)=dy_inv*(yp-ym)
   if(yp > ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
  end do
 end select
 end subroutine map2dy_part_sind

 subroutine map2dz_part_sind(np,sind,ic1,zm,pt)
 integer,intent(in) :: np,sind,ic1
 real(dp),intent(in) :: zm
 real(dp),intent(inout) :: pt(:,:)
 real(dp) :: zp,zs,zimn
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
   zp=pt(ic1,n)
   pt(ic1,n)=zimn+dz_inv*(zp-zs)
   if(zp <= zs)pt(ic1,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
  end do
 case(2)       !z>0
  zs=str_zgrid%smax
  if(zm>zs)then
   zimn=dzi_inv*atan(sz_rat*(zs-zm))
  else
   zimn=dz_inv*(zs-zm)
  endif
  do n=1,np
   zp=pt(ic1,n)
   pt(ic1,n)=dz_inv*(zp-zm)
   if(zp > zs)pt(ic1,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
  end do
 end select
 end subroutine map2dz_part_sind

 !--------------------------

 subroutine map3d_part_sind(pt,np,sind,ic1,ic2,ym,zm)
 integer,intent(in) :: np,sind,ic1,ic2
 real(dp),intent(in) :: ym,zm
 real(dp),intent(inout) :: pt(:,:)
 real(dp) :: yp,zp,ys,zs,ximn,zimn
 integer :: n
 !========================
 !  enter the y=part(ic1,n) z=part(ic2,n) particle positions
 !        in stretched grids    y=y(xi), z(zi)
 !  exit   xi=part(ic1,n) zi=part(ic2,n)
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
   yp=pt(ic1,n)
   zp=pt(ic2,n)
   pt(ic1,n)=ximn+dy_inv*(yp-ys)
   pt(ic2,n)=zimn+dz_inv*(zp-zs)
   if(yp < ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp < zs)pt(ic2,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
  end do
 case(2)       !z<0
  zs=str_zgrid%smin
  zimn=-dzi_inv*atan(sz_rat*(zm-zs))
  do n=1,np
   yp=pt(ic1,n)
   pt(ic1,n)=dy_inv*(yp-ym)
   zp=pt(ic2,n)
   pt(ic2,n)=zimn+dz_inv*(zp-zs)
   if(zp < zs)pt(ic2,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
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
   yp=pt(ic1,n)
   zp=pt(ic2,n)
   pt(ic1,n)=dy_inv*(yp-ym)
   pt(ic2,n)=zimn+dz_inv*(zp-zs)
   if(yp > ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp < zs)pt(ic2,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
  end do
 case(4)       !y>0
  ys=str_ygrid%smax
  if(ym>ys)then
   ximn=dyi_inv*atan(sy_rat*(ys-ym))
  else
   ximn=dy_inv*(ys-ym)
  endif
  do n=1,np
   yp=pt(ic1,n)
   zp=pt(ic2,n)
   pt(ic2,n)=dz_inv*(zp-zm)
   pt(ic1,n)=dy_inv*(yp-ym)
   if(yp > ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
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
   yp=pt(ic1,n)
   zp=pt(ic2,n)
   pt(ic1,n)=dy_inv*(yp-ym)
   pt(ic2,n)=dz_inv*(zp-zm)
   if(yp > ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp > zs)pt(ic2,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
  end do
 case(6)       !z>0
  zs=str_zgrid%smax
  if(zm>zs)then
   zimn=dzi_inv*atan(sz_rat*(zs-zm))
  else
   zimn=dz_inv*(zs-zm)
  endif
  do n=1,np
   yp=pt(ic1,n)
   pt(ic1,n)=dy_inv*(yp-ym)
   zp=pt(ic2,n)
   pt(ic2,n)=dz_inv*(zp-zm)
   if(zp > zs)pt(ic2,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
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
   yp=pt(ic1,n)
   zp=pt(ic2,n)
   pt(ic1,n)=ximn+dy_inv*(yp-ys)
   pt(ic2,n)=dz_inv*(zp-zm)
   if(yp < ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
   if(zp > zs)pt(ic2,n)=zimn+dzi_inv*atan(sz_rat*(zp-zs))
  end do
 case(8)       !y<0
  ys=str_ygrid%smin
  ximn=-dyi_inv*atan(sy_rat*(ym-ys))
  do n=1,np
   zp=pt(ic2,n)
   pt(ic2,n)=dz_inv*(zp-zm)
   yp=pt(ic1,n)
   pt(ic1,n)=ximn+dy_inv*(yp-ys)
   if(yp < ys)pt(ic1,n)=ximn+dyi_inv*atan(sy_rat*(yp-ys))
  end do
 end select
 end subroutine map3d_part_sind

 !--------------------------

 end module grid_param
