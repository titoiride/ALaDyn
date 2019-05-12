
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

 module curr_and_fields_util
 use array_wspace
 use grid_param
 use mpi_curr_interface
 use mpi_field_interface
 use grid_part_connect
 use grid_fields

 implicit none
 !===============================
 ! MOVING WINDOW SECTION
 !=============================
 contains
 !============================
 subroutine set_lpf_acc(ef,sp_loc,apt,np,ndm,nf,nst,xmn,ymn,zmn)
 
 real(dp),intent(in) :: ef(:,:,:,:)
 type(species),intent(in) :: sp_loc
 real(dp),intent(inout) :: apt(:,:)
 integer,intent(in) :: np,ndm,nf,nst
 real(dp),intent(in) :: xmn,ymn,zmn
 ! Uses alternating order quadratic or linear shapes

 select case(ndm)
 case(1)
  call set_part1d_acc(ef,sp_loc,apt,np,nf,xmn)
 case(2)
  call set_part2d_hcell_acc(ef,sp_loc,apt,np,nf,nst,xmn,ymn)
 case(3)
  call set_part3d_hcell_acc(ef,sp_loc,apt,np,nst,xmn,ymn,zmn)
 end select
 end subroutine set_lpf_acc

 subroutine field_charge_multiply(sp_loc,apt,n0,np,nf)

  type(species),intent(in) :: sp_loc
  real(dp),intent(inout) :: apt(:,:)
  integer,intent(in) :: n0,np,nf
  integer :: p,ch
  ch=size(apt,2)
 !==========================
   do p=n0,np
    wgh_cmp=sp_loc%part(p,ch)
    apt(p,1:nf)=charge*apt(p,1:nf)
   end do
 ! EXIT p-assigned (E,B) fields multiplied by charge
 end subroutine field_charge_multiply

 subroutine curr_accumulate(sp_loc,pdata,curr,npt0,npt,f_ch,n_st,xb,yb,zb)
 type(species),intent(in) :: sp_loc
 real(dp),intent(inout) :: pdata(:,:),curr(:,:,:,:)
 integer,intent(in) :: npt0,npt,f_ch,n_st
 ! real(dp),intent(in) :: dtloc
 real(dp),intent(in) :: xb,yb,zb
 !=========================
 ! charge preservinga iform=0, 1
 !iform=0 => (D_xJx,D_yJy,D_zJz) inverted on each particle=> (Jx,Jy,Jz)
 !iform=1    as iform=0
 !iform=2    no charge preserving
 !=========================
 if(npt==0)return
  if(ndim < 3)then
   if(f_ch <2)then
    call esirkepov_2d_curr(sp_loc,pdata,curr,npt0,npt,n_st,curr_ndim,ndim,xb,yb)
   else
    call ncdef_2d_curr(sp_loc,pdata,curr,npt0,npt,n_st,curr_ndim,ndim,xb,yb)
   endif
   return
  endif
  if(f_ch <2)then
   call esirkepov_3d_curr(sp_loc,pdata,curr,npt0,npt,n_st,xb,yb,zb)
  else
   call ncdef_3d_curr(sp_loc,pdata,curr,npt0,npt,n_st,xb,yb,zb)
  endif
 !========================
 ! accumulates for each species currents on curr(i1:n1p,j1:n2p,k1:n3p,1:compnent)
 !============================
 end subroutine curr_accumulate
 !==============================
 subroutine curr_mpi_collect(curr,i1,n1p,j1,n2p,k1,n3p,njc)

 real(dp),intent(inout) :: curr(:,:,:,:)
 integer,intent(in) :: i1,n1p,j1,n2p,k1,n3p,njc
 integer :: i,j,k,jj,kk
 real(dp) :: dery,derhy,derz,derhz
 !============sums data on ghost points
 if(prl)then
  call fill_curr_yzxbdsdata(curr,i1,n1p,j1,n2p,k1,n3p,njc)
 endif
 call jc_xyzbd(curr,i1,n1p,j1,n2p,k1,n3p,njc)
 !=================
 if(iform <2)then
  do i=1,ndim
   curr(i1:n1p,j1:n2p,k1:n3p,i)=djc(i)*curr(i1:n1p,j1:n2p,k1:n3p,i)
  end do
 endif
 if(Stretch)then
  select case(curr_ndim)
  case(2)
   do k=k1,n3p
    do j=j1,n2p
     jj=j-2
     dery=loc_yg(jj,3,imody)
     derhy=loc_yg(jj,4,imody)
     do i=i1,n1p
      curr(i,j,k,1)=dery*curr(i,j,k,1)
      curr(i,j,k,2)=derhy*curr(i,j,k,2)
     end do
    end do
   end do
  case(3)
   do k=k1,n3p
    kk=k-2
    derz=loc_zg(kk,3,imodz)
    derhz=loc_zg(kk,4,imodz)
    do j=j1,n2p
     jj=j-2
     dery=loc_yg(jj,3,imody)
     derhy=loc_yg(jj,4,imody)
     do i=i1,n1p
      curr(i,j,k,1)=dery*derz*curr(i,j,k,1)
      curr(i,j,k,2)=derhy*derz*curr(i,j,k,2)
      curr(i,j,k,3)=dery*derhz*curr(i,j,k,3)
     end do
    end do
   end do
  end select
 endif
 end subroutine curr_mpi_collect
 !=================================
 ! END SECTION ON GRID DEFINED PARTICLE VARIABLES
 !================================================
 !==========================
 subroutine pfields_prepare(ef,i1,i2,j1,j2,k1,k2,nc,spr,spl)
 real(dp),intent(inout) :: ef(:,:,:,:)
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,nc,spr,spl
 !===================
 ! Enter fields ef[i1:i2,j1,j2,k1,k2,1:nc)
 ! Exit fields ef[i1:i2,j1,j2,k1,k2,1:nc)
 !===========================
 if(prl)then
  call fill_ebfield_yzxbdsdata(ef,i1,i2,j1,j2,k1,k2,1,nc,spr,spl)
  !======================
  ! Adds point data => spl to the left spr to the right  
  ! Sets periodic BCs for
  ! iby,ibz,ibx =2
  !==================
 endif
 call field_xyzbd(ef,i1,i2,j1,j2,k1,k2,nc,spl,spr)
 ! extends for one point at the box boundaries
 !==================================================
 end subroutine pfields_prepare
!================================
 subroutine advance_lpf_fields(ef,curr,dt_lp,v_b,&
                               i1,i2,j1,j2,k1,k2,ibd)
 real(dp),intent(inout) :: ef(:,:,:,:)
 real(dp),intent(in) :: curr(:,:,:,:)
 real(dp),intent(in) :: dt_lp,v_b
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,ibd
 integer :: str,stl,ik
 real(dp) :: dtx,dty,dtz,dth_lp,dthx,dthy,dthz

 dth_lp=0.5*dt_lp

 dtx=dx_inv*dt_lp
 dty=dy_inv*dt_lp
 dtz=dz_inv*dt_lp

 dthx=0.5*dtx
 dthy=0.5*dty
 dthz=0.5*dtz

 !=======================
 ! A LPF order Lpf with time-centered source term
 !=============================
 if(prl)then
  str=0
  stl=1
  call fill_ebfield_yzxbdsdata(ef,i1,i2,j1,j2,k1,k2,1,curr_ndim,str,stl)
  ! To fill electric field data
  ! sends stl to the left,
  ! recvs stl points from right at (nyp+stl), (nzp+stl)
 endif
 !============== first substep dt/2 advance of B-field
 call ef_bds(ef,i1,i2,j1,j2,k1,k2,zero_dp,ibd)
 ! Uses upper BCs of E fields: ibd=0 for inflow-outflow
 !                             ibd=1 for symmetric
 if(Comoving)then
  call field_xcomov_advect(ef,i1,i2,j1,j2,k1,k2,curr_ndim+1,nfield,dthx,-v_b,0)
 ! Solves for one-half step backward advection B field explicit
 ! B^{n} => B^{n}+v_b*Dth[Dx]B^n  v_b >0
 endif
 !==================================
 ! solves B^{n+1/2}= B^n -Dth[rot E]^n 
 !============================
 call rotE(ef,i1,i2,j1,j2,k1,k2,dthx,dthy,dthz)
 !=============================
 !============== central tep dt advance of E-field
 if(prl)then
  str=1
  stl=0
  call fill_ebfield_yzxbdsdata(ef,i1,i2,j1,j2,k1,k2,curr_ndim+1,nfield,str,stl)
  ! sends nyp+1-str to the right
  ! recvs str points from left at (1-str)
 endif
 !======================
 call bf_bds(ef,i1,i2,j1,j2,k1,k2,dt_lp,ibd)
 ! Uses lower BCs for B fields: ibd=0 for inflow-outflow
 !                             ibd=1 for symmetric
 if(Comoving)then
  call field_xcomov_advect(ef,i1,i2,j1,j2,k1,k2,1,curr_ndim,dthx,-v_b,0)
 ! Solves for half-step backward advection E field expplicit
 ! E^{n} => E^{n}+v_b*Dth[Dx]E^n
 endif
 !=======================
 ! solves E^{n+1}= E^n +DT[rot B]^{n+1/2}- ompe*DT*J^{n+1/2} 
 !==================
 call rotB(ef,i1,i2,j1,j2,k1,k2,dtx,dty,dtz)
 !===================
 ! adds currents
 do ik=1,curr_ndim
  ef(i1:i2,j1:j2,k1:k2,ik)=ef(i1:i2,j1:j2,k1:k2,ik)-&
                           ompe*curr(i1:i2,j1:j2,k1:k2,ik)
 end do
 if(Comoving)then
  call field_xcomov_advect(ef,i1,i2,j1,j2,k1,k2,1,curr_ndim,dthx,-v_b,2)
 ! Solves for backward advection E field implicit
 ! E^{n+1} => E^{n}+v_b*Dth[Dx]E^{n+1}
 endif
 !============== second substep dt/2 advance of B-field
 if(prl)then
  str=0
  stl=1
  call fill_ebfield_yzxbdsdata(ef,i1,i2,j1,j2,k1,k2,1,curr_ndim,str,stl)
 endif
 ! E field gets stl points from right (nyp+stl), (nzp+stl)
 call ef_bds(ef,i1,i2,j1,j2,k1,k2,dt_lp,ibd)
 ! solves B^{n+1}= B^{n+1/2} -Dth[rot E]^{n+1} 
 !===================
 call rotE(ef,i1,i2,j1,j2,k1,k2,dthx,dthy,dthz)
 !==============
 if(Comoving)then
  call field_xcomov_advect(ef,i1,i2,j1,j2,k1,k2,curr_ndim+1,nfield,dthx,-v_b,2)
 ! Solves for one-half step backward advection B field implicit
 ! B^{n+1} => B^{n+1/2}+v_b*Dth[Dx]B^{n+1}
 endif
 !===============
 end subroutine advance_lpf_fields
!===================================================
 subroutine advect_bunch_fields(fb,curr,&
  v_b,dt_lp,i1,i2,j1,j2,k1,k2,init_time)

 real(dp),intent(inout) :: fb(:,:,:,:),curr(:,:,:,:)
 real(dp),intent(in) :: v_b,dt_lp
 integer,intent(in) :: i1,i2,j1,j2,k1,k2
 logical,intent(in) :: init_time
 real(dp) :: dtx,dty,dtz,dthx,dthy,dthz
 integer :: ix,iy,iz

 dtx=dx_inv*dt_lp
 dty=dy_inv*dt_lp
 dtz=dz_inv*dt_lp

 dthx=0.5*dtx
 dthy=0.5*dty
 dthz=0.5*dtz
 !In 2D nfield=3 nbfield=4=nfield+1   in fb(4)=Jx[i+1/2,j,k] at t=0
 !fb=[Ex,Ey,Ez,Jx,By,Bz]
 !================================================
 call fill_ebfield_xbdsdata(&
                     fb,i1,i2,j1,j2,k1,k2,1,nbfield,1,1)
 if(init_time)then
  call field_xadvect(fb,i1,i2,j1,j2,k1,k2,4,4,dthx,-v_b)
  fb(i1:i2,j1:j2,k1:k2,4)=dt_lp*fb(i1:i2,j1:j2,k1:k2,4)
 endif
 call field_xadvect(fb,i1,i2,j1,j2,k1,k2,1,nbfield,dtx,v_b)
 do iz=k1,k2
  do iy=j1,j2
   do ix=i1,i2
    curr(ix,iy,iz,1)=curr(ix,iy,iz,1)-fb(ix,iy,iz,4)
  end do
  end do
 end do
 ! Subtracts from the ongitudinal bunch current the advected initial bunch
 ! current  
 end subroutine advect_bunch_fields
 !============================
 end module curr_and_fields_util
!==============================
