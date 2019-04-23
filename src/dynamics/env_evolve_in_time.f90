
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

 module env_evolve
 use window
 use boris_push
 use curr_and_fields_util
 use mpi_part_interface
 use ionize
 use fluid_density_momenta

 implicit none
 !===============================
 ! MOVING WINDOW SECTION
 !=============================
 contains
 !============================
 subroutine advance_lpf_envelope(curr,evf,dt_loc,omg,i1,nxp,j1,nyp,k1,nzp)

 real(dp),intent(inout) :: curr(:,:,:,:),evf(:,:,:,:)
 real(dp),intent(in) :: dt_loc,omg
 integer,intent(in) :: i1,nxp,j1,nyp,k1,nzp
 integer :: i,j,k
 integer :: nst,str,stl,cind,ib

 !====== enter env(3:4)=A^{n-1} and env(1:2)= A^{n}
 ! enters jc(3)=<wgh*n/gamp> >0

 nst=0
 if(Stretch)nst=1
 !ord=2
 cind=1                !cind=0 FFT cind=1 grid deriv
 ib=2          !ib=1 implicit ib=2 optimazid explicit
 !optimized advection scheme
 if(Comoving)ib=0
 if(prl)then
  str=1
  stl=1
  call fill_ebfield_yzxbdsdata(&
                evf,i1,nxp,j1,nyp,k1,nzp,1,2,str,stl)
 endif
 do k=k1,nzp
  do j=j1,nyp
   do i=i1,nxp
    curr(i,j,k,1)=-ompe*curr(i,j,k,3)*evf(i,j,k,1)
    curr(i,j,k,2)=-ompe*curr(i,j,k,3)*evf(i,j,k,2)
   end do
  end do
 end do
 !  jc(1:2)=-ompe*chi*env(1:2)  the J_{env} source term
!==================================
 if(ib==0)then
  !call env_comov_maxw_solve(jc,evf,i1,nxp,j1,nyp,k1,nzp,omg,&
                    !dx_inv,dy_inv,dz_inv,dt_loc)
  call env_lpf_solve(jc,evf,ib,i1,nxp,j1,nyp,k1,nzp,omg,&
                    dx_inv,dy_inv,dz_inv,dt_loc)
 !========================================================
 else
 !=================== second order in time full wave equation
  call env_maxw_solve(jc,evf,i1,nxp,j1,nyp,k1,nzp,omg,&
                    dx_inv,dy_inv,dz_inv,dt_loc)
 endif
 ! =================================
 end subroutine advance_lpf_envelope
 ! =================================
 ! ENVELOPE model in LP regime
 !============================
 subroutine env_den_collect(eden,i1,n1p,j1,n2p,k1,n3p)

  real(dp),intent(inout) :: eden(:,:,:,:)
  integer,intent(in) :: i1,n1p,j1,n2p,k1,n3p
  integer :: i,j,k,kk,jj,ic
  real(dp) :: dery,derz

  ic=1
  if(prl)then
   call fill_curr_yzxbdsdata(eden,i1,n1p,j1,n2p,k1,n3p,ic)
  endif
  call den_zyxbd(eden,i1,n1p,j1,n2p,k1,n3p,ic)
 !=======Enters normalized <w*n/gam> > 0
 !==================================
 if(Stretch)then
  ic=1
  if(ndim==2)then
   k=1
   do j=j1,n2p
    jj=j-2
    dery=loc_yg(jj,3,imody)
    do i=i1,n1p
     eden(i,j,k,ic)=dery*eden(i,j,k,ic)
    end do
   end do
   return
  endif
  do k=k1,n3p
   kk=k-2
   derz=loc_zg(kk,3,imodz)
   do j=j1,n2p
    jj=j-2
    dery=derz*loc_yg(jj,3,imody)
    do i=i1,n1p
     eden(i,j,k,ic)=dery*eden(i,j,k,ic)
    end do
   end do
  end do
 endif
 !=========================
 !exit in eden(ic)  the source terms chi() >0 of the envelope equation
 !-----------------------------------------------
 end subroutine env_den_collect
 !=========================
 subroutine env_two_fields_average(evf,ev1f,av,i1,i2,j1,j2,k1,k2,spl_in,spr_in)
 real(dp),intent(in) :: evf(:,:,:,:),ev1f(:,:,:,:)
 real(dp),intent(out) :: av(:,:,:,:)
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,spl_in,spr_in
 integer :: ix,iy,iz
 real(dp) :: ar,ai
 !===================
 do iz=k1,k2
  do iy=j1,j2
   do ix=i1,i2
    ar=0.5*(evf(ix,iy,iz,1)+evf(ix,iy,iz,3))   !A^{n+1/2}=(A^n+1+A^n)/2
    ai=0.5*(evf(ix,iy,iz,2)+evf(ix,iy,iz,4))
    av(ix,iy,iz,1)=0.5*(ar*ar+ai*ai)
    ar=0.5*(ev1f(ix,iy,iz,1)+ev1f(ix,iy,iz,3))   !A^{n+1/2}=(A^n+1+A^n)/2
    ai=0.5*(ev1f(ix,iy,iz,2)+ev1f(ix,iy,iz,4))
    av(ix,iy,iz,1)=av(ix,iy,iz,1)+0.5*(ar*ar+ai*ai)
    ! |A|^2/2 at t^{n+1/2}=> gamp^{n+1/2}
    !  NO overlap assumed
   end do
  end do
 end do
 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,1,1,spr_in,spl_in)
 !call field_xyzbd(av,i1,i2,j1,j2,k1,k2,1,spr_in,spl_in)
 !=====================
 end subroutine env_two_fields_average
!===========================
 subroutine env_fields_average(evf,av,i1,i2,j1,j2,k1,k2,spl_in,spr_in)
 real(dp),intent(in) :: evf(:,:,:,:)
 real(dp),intent(out) :: av(:,:,:,:)
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,spl_in,spr_in
 integer :: ix,iy,iz,ord
 real(dp) :: ar,ai
 !===================
 ord=2
 do iz=k1,k2
  do iy=j1,j2
   do ix=i1,i2
    ar=0.5*(evf(ix,iy,iz,1)+evf(ix,iy,iz,3))   !A^{n+1/2}=(A^n+1+A^n)/2
    ai=0.5*(evf(ix,iy,iz,2)+evf(ix,iy,iz,4))
    av(ix,iy,iz,1)=0.5*(ar*ar+ai*ai)
    ! |A|^2/2 at t^{n+1/2}=> gamp^{n+1/2}
   end do
  end do
 end do
 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,1,1,spr_in,spl_in)
 call env_grad(av,i1,i2,j1,j2,k1,k2,ord,dx_inv,dy_inv,dz_inv)
 !Exit staggered grad|A|^2/2 in jc(2:4) or jc(2:3) 

 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,2,curr_ndim+1,spr_in,spl_in)
 !=====================
 end subroutine env_fields_average
 !===========================
 subroutine env_amp_prepare(envf,av,i1,i2,j1,j2,k1,k2,ord,spl_in,spr_in)
 real(dp),intent(in) :: envf(:,:,:,:)
 real(dp),intent(out) :: av(:,:,:,:)
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,ord,spl_in,spr_in
 integer :: ix,iy,iz,spl,spr
 !real(dp) :: ar,ai
 !===================
 do iz=k1,k2
  do iy=j1,j2
   do ix=i1,i2
    av(ix,iy,iz,1)=0.5*(envf(ix,iy,iz,1)*envf(ix,iy,iz,1)+&
     envf(ix,iy,iz,2)*envf(ix,iy,iz,2))
    !|A|^2/2 at current t^n time level
   end do
  end do
 end do
 spl=spl_in
 spr=spr_in
 if(spl >2)spl=2
 if(spr >2)spr=2

 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,1,1,spr,spl)

 call env_grad(av,i1,i2,j1,j2,k1,k2,ord,dx_inv,dy_inv,dz_inv)
 !Exit staggered grad|A|^2/2 in jc(2:4) or jc(2:3) 

 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,1,curr_ndim+1,spr,spl)

 !call field_xyzbd(av,i1,i2,j1,j2,k1,k2,nj_dim,spr,spl)
 !=====================
 end subroutine env_amp_prepare
!=============================
 subroutine env_amp_two_fields_prepare(envf,env1f,av,i1,i2,j1,j2,k1,k2,ord,spl_in,spr_in)
 real(dp),intent(in) :: envf(:,:,:,:),env1f(:,:,:,:)
 real(dp),intent(out) :: av(:,:,:,:)
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,ord,spl_in,spr_in
 integer :: ix,iy,iz,spl,spr
 !real(dp) :: ar,ai
 !===================
 do iz=k1,k2
  do iy=j1,j2
   do ix=i1,i2
    av(ix,iy,iz,1)=0.5*(envf(ix,iy,iz,1)*envf(ix,iy,iz,1)+&
     envf(ix,iy,iz,2)*envf(ix,iy,iz,2))
    av(ix,iy,iz,1)=av(ix,iy,iz,1)+0.5*(env1f(ix,iy,iz,1)*env1f(ix,iy,iz,1)+&
     env1f(ix,iy,iz,2)*env1f(ix,iy,iz,2))
    !|A|^2/2 at current t^n time level
   end do
  end do
 end do
 spl=spl_in
 spr=spr_in
 if(spl >2)spl=2
 if(spr >2)spr=2

 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,1,1,spr,spl)

 call env_grad(av,i1,i2,j1,j2,k1,k2,ord,dx_inv,dy_inv,dz_inv)
 !Exit staggered grad|A|^2/2 in jc(2:4) or jc(2:3)

 if(prl)call fill_ebfield_yzxbdsdata(av,i1,i2,j1,j2,k1,k2,2,curr_ndim+1,spr,spl)

 !call field_xyzbd(av,i1,i2,j1,j2,k1,k2,nj_dim,spr,spl)
 !=====================
 end subroutine env_amp_two_fields_prepare
 !=======================================
 subroutine env_lpf2_evolve(dt_loc,it_loc,initial_time)

 real(dp),intent(in) :: dt_loc
 integer,intent(in) :: it_loc
 logical,intent(in) :: initial_time
 integer :: np,ic,nyf,nzf,n_st
 integer :: i1,j1,k1,i2,id_ch
 real(dp) :: xm,ym,zm,Ltz,ef2_ion,loc_ef2_ion(2)
 !============================
 xm=loc_xgrid(imodx)%gmin
 zm=loc_zgrid(imodz)%gmin
 ym=loc_ygrid(imody)%gmin

 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 j1=loc_ygrid(imody)%p_ind(1)
 nyf=loc_ygrid(imody)%p_ind(2)
 k1=loc_zgrid(imodz)%p_ind(1)
 nzf=loc_zgrid(imodz)%p_ind(2)
 n_st=0
 if(Stretch)n_st=str_indx(imody,imodz)
 !====================
 if(prl)call fill_ebfield_yzxbdsdata(&
                ebf,i1,i2,j1,nyf,k1,nzf,1,nfield,2,2)
!======================================
 if(Ionization)then
  if(it_loc==0)then
   call init_random_seed(mype)
  endif
  id_ch=nd2+1
  if(Two_color)then
   if(prl)call fill_ebfield_yzxbdsdata(&
                env1,i1,i2,j1,nyf,k1,nzf,1,2,2,2)
   do ic=2,nsp_ionz
    np=loc_npart(imody,imodz,imodx,ic)
    if(np>0)then
     call set_ion_env_field(env1,spec(ic),ebfp,np,ndim,n_st,xm,ym,zm,om1)
     loc_ef2_ion(1)=maxval(ebfp(1:np,id_ch))
     loc_ef2_ion(1)=sqrt(loc_ef2_ion(1))
     if(ef2_ion > lp_max)then
      lp_max=1.1*ef2_ion
      call set_field_ioniz_wfunction(&
       ion_min(ic-1),atomic_number(ic-1),ic,ionz_lev,ionz_model,lp_max,dt_loc)
     endif
     call ionization_cycle(spec(ic),ebfp,np,ic,it_loc,1,de_inv)
    endif
   end do
!================
  else
   if(prl)call fill_ebfield_yzxbdsdata(&
                env,i1,i2,j1,nyf,k1,nzf,1,2,2,2)
   do ic=2,nsp_ionz
    np=loc_npart(imody,imodz,imodx,ic)
    if(np>0)then
     call set_ion_env_field(env,spec(ic),ebfp,np,ndim,n_st,xm,ym,zm,oml)
      loc_ef2_ion(1)=maxval(ebfp(1:np,id_ch))
      loc_ef2_ion(1)=sqrt(loc_ef2_ion(1))
      ef2_ion=loc_ef2_ion(1)
      if(ef2_ion > lp_max)then
       write(6,'(a22,i6,2E11.4)')'reset high ionz field ',mype,ef2_ion,lp_max
       lp_max=1.1*ef2_ion
       call set_field_ioniz_wfunction(&
       ion_min(ic-1),atomic_number(ic-1),ic,ionz_lev,ionz_model,lp_max,dt_loc)
      endif
     call ionization_cycle(spec(ic),ebfp,np,ic,it_loc,1,de_inv)
    endif
   end do
  endif
 endif
!=================================
 ic=1
 Ltz=Lorentz_fact(ic)
 !Lorentz force on a particle => ebfp0(1:3), velocity ebfp0(3:4) stored
 ! all the time for t^{n-1} levels ==========
 !===========================
 jc(:,:,:,:)=0.0
 np=loc_npart(imody,imodz,imodx,ic)
 if(Two_color)then
  call env_amp_two_fields_prepare(env,env1,jc,i1,i2,j1,nyf,k1,nzf,2,2,2)
 else
  call env_amp_prepare(env,jc,i1,i2,j1,nyf,k1,nzf,2,2,2)
 endif
  !======================================
  ! exit jc(1)=|a|^2/2 at t^n
  !      jc(2:4)=grad|a|^2/2 at t^n
  ! For two-color |A|= |A_0|+|A_1|
  !======================================
  ebfp(:,:)=0.0
   call set_env_acc(ebf,jc,spec(ic),ebfp,np,curr_ndim,n_st,dt_loc,xm,ym,zm)
                            !call field_charge_multiply(spec(ic),ebfp,1,np,nfield)
  !exit ebfp(1:3)=q*[E+F] ebfp(4:6)=q*B/gamp, ebfp(7)=wgh/gamp at t^n
  !Lorentz force already multiplied by particle charge
  !jc(1:4) not modified
  !====================
   call lpf_env_momenta(spec(ic),ebfp,np,dt_loc,Ltz)
  ! Updates particle momenta P^{n-1/2} => P^{n+1/2}
  ! stores in ebfp(1:3)=old (x,y,z)^n ebfp(7)=wgh/gamp >0
  !======================
  if(Hybrid)then
 !+++++++++++++++++++++++++++++++++++++++++++++++++++++++
   call set_env_momentum_density_flux(up,ebf,jc,ebf0,flux,i1,i2,j1,nyf,k1,nzf)
    !exit jc(1)=q^2*n/gam, jc(2:4) ponderomotive force on a grid
    !ebf0= total fields flux(1:4)=(P,den)^n 
!============================
   call update_adam_bash_fluid_variables(up,up0,flux,ebf0,dt_loc,&
                                        i1,i2,j1,nyf,k1,nzf,it_loc,Ltz,initial_time)
   ! In up exit updated momenta-density variables u^{n+1}
   ! in  u0^{n} stores Dt*F(u^n), in flux(1:fdim)=(P,den)^{n+1/2}

   flux(i1:i2,j1:nyf,k1:nzf,curr_ndim+2)=jc(i1:i2,j1:nyf,k1:nzf,1)

   ! in flux(fdim+1) exit the fluid contribution of the sorce term q^2*n/gam
   ! for the envelope field solver
  endif
  jc(:,:,:,1)=0.0
   call set_env_density(ebfp,jc,np,curr_ndim,1,n_st,xm,ym,zm)
   call env_den_collect(jc,i1,i2,j1,nyf,k1,nzf)
  ! in jc(1)the particle contribution of the source term <q^2*n/gamp>
  ! to be added to the fluid contribution if (Hybrid)
  jc(i1:i2,j1:nyf,k1:nzf,3)=jc(i1:i2,j1:nyf,k1:nzf,1)
  if(Hybrid)then
   jc(i1:i2,j1:nyf,k1:nzf,3)=jc(i1:i2,j1:nyf,k1:nzf,3)+&
                             flux(i1:i2,j1:nyf,k1:nzf,curr_ndim+2)
  endif
!===================
 ! in the envelope equation (A^{n-1},A^n)==> (A^n,A^{n+1})
 ! jc(3) = <q^2n/gam>
 ! Jc(1:2)=-ompe*jc(3)*A at level t^n
 !==================================================
 !==================================================
 call advance_lpf_envelope(jc,env,dt_loc,oml,i1,i2,j1,nyf,k1,nzf)
 !advance (A^n, J^n) => A^{n+1}, A^{n-1}=> A^n
 ! jc(3) not modified
 if(Two_color)call advance_lpf_envelope(jc,env1,dt_loc,om1,i1,i2,j1,nyf,k1,nzf)
 !advance (A_1^n, J^n) => A_1^{n+1}, A_1^{n-1}=> A_1^n
!=======================
 if(Two_color)then
  call env_two_fields_average(env,env1,jc,i1,i2,j1,nyf,k1,nzf,2,2)
 else
  call env_fields_average(env,jc,i1,i2,j1,nyf,k1,nzf,2,2)
 endif
 ! In jc(1)= Phi= |A|^2/2 +|A_1|/2 at t^{n+1/2} 
 if(Hybrid)then
  flux(i1:i2,j1:nyf,k1:nzf,curr_ndim+2)=jc(i1:i2,j1:nyf,k1:nzf,1)
  !stores in flux()
 endif
 call set_env_grad_interp(jc,spec(ic),ebfp,np,curr_ndim,n_st,xm,ym,zm)
  !=============================
  ! Exit p-interpolated field variables
  ! at time level t^{n+1/2} and positions at time t^n
  ! in ebfp(1:3)=grad|A|^2/2 ebfp(4)=|A|^2/2 in 3D
  ! in ebfp(1:2)=grad|A|^2/2 ebfp(3)=|A|^2/2 in 2D
  !=====================================
   call lpf_env_positions(spec(ic),ebfp,np,dt_loc,vbeam)
  !===========================
  ! ebfp(1:3) dt*V^{n+1/2}  ebfp(4:6) old positions for curr J^{n+1/2}
  ! ebfp(7)=dt*gam_inv
  !=======collects in jc(1:curr_ndim) currents due to electrons
  jc(:,:,:,:)=0.0
  call curr_accumulate(spec(ic),ebfp,jc,1,np,iform,n_st,xm,ym,zm)
 !===========================
  call curr_mpi_collect(jc,i1,i2,j1,nyf,k1,nzf,curr_ndim)
 if(Hybrid)then
  !In flux(1:fdim+1) are stored fluid (P,den) at t^{n+1/2}
  !In flux(fdim+1) is stored |A|^2/2 at t^{n+1/2}
  call fluid_curr_accumulate(flux,jc,dt_loc,i1,i2,j1,nyf,k1,nzf)
  !Computes fluid contribution => J^{n+1/2} and adds to particle contribution
 endif
 !====================
 ! Jc(1:3) for total curr Dt*J^{n+1/2}
 lp_in(1)=lp_in(1)+dt_loc
 call advance_lpf_fields(ebf,jc,dt_loc,vbeam,&
                          i1,i2,j1,nyf,k1,nzf,0)
 ! (E,B) fields at time t^{n+1}
 !-----------------------------
 end subroutine env_lpf2_evolve
!=============== END ENV PIC SECTION
 subroutine ENV_run(t_loc,dt_loc,iter_loc)

 real(dp),intent(in) :: t_loc,dt_loc
 integer,intent(in) :: iter_loc
 logical,parameter :: mw=.false.
 logical :: init_time
 !+++++++++++++++++++++++++++++++++
 !for vbeam >0 uses the xw=(x+vbeam*t)
 !x=xi=(xw-vbeam*t) fixed
 !+++++++++++++++++++++++++++++++++
 init_time=.false.
 if(t_loc<1.e-08)init_time=.true.
 if(w_speed>0.0)then ! moves the computational box with w_speed>0.
  if(iter_loc==0)call LP_window_xshift(dt_loc,w_sh,iter_loc)
  if(t_loc>=wi_time.and.t_loc < wf_time)then
   if(mod(iter_loc,w_sh)==0)then
    call LP_window_xshift(dt_loc,w_sh,iter_loc)
   endif
  endif
 endif
 if(Comoving)then
  if(t_loc>=wi_time.and.t_loc< wf_time)then
   if(mod(iter_loc,w_sh)==0)then
    call comoving_coordinate(vbeam,dt_loc,w_sh,iter_loc)
   endif
  endif
 endif
  !=========================
  call env_lpf2_evolve(dt_loc,iter_loc,init_time)
!================================
 if(Part)call cell_part_dist(mw)
 !
 end subroutine ENV_run
 !============================
 ! END ENVELOPE MODULE
 !==============================
 end module env_evolve
!==============================
