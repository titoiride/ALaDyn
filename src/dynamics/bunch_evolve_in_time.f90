
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

 module bunch_evolve
  use window
  use boris_push
  use curr_and_fields_util
  use mpi_part_interface
  use ionize
  use fluid_density_momenta

 implicit none
 !
 contains
 ! MODULE for beam-driven wakefields
 !==============================
 ! Leap-frog integrators
 !===========================
 subroutine eb_fields_collect(ef,efb1,efb2,ef_tot,nfl)
  real(dp),intent(in) :: ef(:,:,:,:),efb1(:,:,:,:),efb2(:,:,:,:)
  real(dp),intent(inout) :: ef_tot(:,:,:,:)
  integer,intent(in) :: nfl

  ef_tot(:,:,:,1:nfl)=ef(:,:,:,1:nfl)+efb1(:,:,:,1:nfl)
  ef_tot(:,:,:,1:3)=ef_tot(:,:,:,1:3)+efb2(:,:,:,1:3)
  if(nfl ==6)then
   ef_tot(:,:,:,5:6)=ef_tot(:,:,:,5:6)+efb2(:,:,:,5:6)
  endif
 end subroutine eb_fields_collect 
!=============================
 subroutine lpf2_eb_evolve(dt_loc,iter_loc,initial_time)

 real(dp),intent(in) :: dt_loc
 integer,intent(in) :: iter_loc
 logical,intent(in) :: initial_time
 integer :: np,ic
 integer :: i1,i2,j1,j2,k1,k2,n_st,id_ch
 real(dp) :: xm,ym,zm,vb,Ltz,ef2_ion(1),loc_ef2_ion(1)
 !============================
 ! Fields are in ebf() (wake) and ebf_bunch() bunches
 ! particles are in spec(1)+ebfp plasma bunch(1:2)+ebfb (drive and witness)
 !==============================
 xm=loc_xgrid(imodx)%gmin
 ym=loc_ygrid(imody)%gmin
 zm=loc_zgrid(imodz)%gmin
 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 j1=loc_ygrid(imody)%p_ind(1)
 j2=loc_ygrid(imody)%p_ind(2)
 k1=loc_zgrid(imodz)%p_ind(1)
 k2=loc_zgrid(imodz)%p_ind(2)
 n_st=0
 if(Stretch)n_st=str_indx(imody,imodz)
 id_ch=nd2+1
 !====================
 !Ghost cell values for field assignement on particles
 vb=vbeam
 call eb_fields_collect(ebf,ebf1_bunch,ebf_bunch,ebf0,nfield)
!   in ebf0() the total fields  bunch+wake
 call pfields_prepare(ebf0,i1,i2,j1,j2,k1,k2,nfield,2,2)
 !call pfields_prepare(ebf_bunch,i1,i2,j1,j2,k1,k2,nfield,1,1)
 !call pfields_prepare(ebf1_bunch,i1,i2,j1,j2,k1,k2,nfield,1,1)
!======== first new plasma electrons are injected by ionization
 if(Ionization)then
  if(iter_loc==0)then
   call init_random_seed(mype)
  endif
  do ic=2,nsp_ionz
   np=loc_npart(imody,imodz,imodx,ic)
   if(np>0)then
    call set_ion_Efield(ebf0,spec(ic),ebfp,np,n_st,ndim,nsp_run,dt_loc,xm,ym,zm)
                       !call set_ion_two_Ebfield(ebf,ebf_bunch,ebf1_bunch,spec(ic),&
                       !                ebfp,np,n_st,ndim,xm,ym,zm)
    if(mod(iter_loc,50)==0)then     !refresh ionization tables
     loc_ef2_ion(1)=maxval(ebfp(1:np,id_ch))
     loc_ef2_ion(1)=sqrt(loc_ef2_ion(1))/514.   !In atomic units
     if(ef2_ion(1) > eb_max)then
      eb_max=1.1*ef2_ion(1)
      call set_field_ioniz_wfunction(ion_min(ic-1),atomic_number(ic-1),ic,&
                                     ionz_lev,ionz_model,eb_max,dt_loc)
     endif
    endif
   call ionization_cycle(&
                        spec(ic),ebfp,np,ic,iter_loc,0,deb_inv)
   endif
  end do
  !======== injects new electrons, with weights equal to ion weights 
 endif
!=======================
! STEP 1
!Advances momenta and position of plasma particles using total field
!==========================
  jc(:,:,:,:)=0.0
  do ic=1,nsp_run
   np=loc_npart(imody,imodz,imodx,ic)
   Ltz=Lorentz_fact(ic)
   if(np >0)then
    call set_lpf_acc(ebf0,spec(ic),ebfp,np,ndim,nfield,n_st,xm,ym,zm)
    call field_charge_multiply(spec(ic),ebfp,1,np,nfield)
!==================================
   !EXIT ebfp(1:6)  total=wake + bunch fields assigned to plasma particle position
!==================================
    if(initial_time)call init_lpf_momenta(spec(ic),ebfp,1,np,dt_loc,Ltz)
    call lpf_momenta_and_positions(spec(ic),ebfp,1,np,dt_loc,vb,Ltz)
!  EXIT p^{n+1/2}, v^{n+1/2}, x^{n+1}
!  in x^{n+1} are stored in ebfp(1:3) old x^n are stored in ebfp((4:6)  
!  in ebfp(7) is stored dt_loc/gamma
   call curr_accumulate(spec(ic),ebfp,jc,1,np,iform,n_st,xm,ym,zm)
   endif
  end do
  call curr_mpi_collect(jc,i1,i2,j1,j2,k1,k2,curr_ndim)
! EXIT current density jc(1:3) due to plasma particles density and velocity at
! time t^{n+1/2}
!======================================================
 if(Hybrid)then
  call set_momentum_density_flux(up,flux,i1,i2,j1,j2,k1,k2)
!============================
  ! in the flux() array exit: (px,py,pz,den,vx,vy,vz) at t^n
!============================
  call update_adam_bash_fluid_variables(&
                        up,up0,flux,ebf0,dt_loc,i1,i2,j1,j2,k1,k2,iter_loc,Ltz,initial_time)
   ! In up exit updated momenta-density variables u^{n+1}
   ! in  u0^{n} stores Dt*F(u^n), in flux(1:fdim)=(P,den)^{n+1/2}
   ! In flux(1:curr_ndim+1) are stored fluid (P,den) at t^{n+1/2}
   flux(i1:i2,j1:j2,k1:k2,curr_ndim+2)=0.0
  call fluid_curr_accumulate(flux,jc,dt_loc,i1,i2,j1,j2,k1,k2)
  !Computes fluid contribution => J_f^{n+1/2} and adds to particle contribution
!  ===============  END plasma fluid section
 endif
 call advance_lpf_fields(ebf,jc,dt_loc,vbeam,i1,i2,j1,j2,k1,k2,1) 
!==============================
! STEP 2
!Advances momenta and position of bunch particles
!==========================
 jc(:,:,:,:)=0.0
 do ic=1,nsb
  np=loc_nbpart(imody,imodz,imodx,ic)
  Ltz=Lorentz_bfact(ic)
  if(np >0)then
   call set_lpf_acc(ebf0,bunch(ic),ebfb,np,ndim,nfield,n_st,xm,ym,zm)
  endif
  call field_charge_multiply(bunch(ic),ebfb,1,np,nfield)
!==================================
   !EXIT ebfb(1:6)  total=wake + bunch fields assigned to bunch particle position
!==================================
  if(initial_time)call init_lpf_momenta(bunch(ic),ebfb,1,np,dt_loc,Ltz)
  call lpf_momenta_and_positions(bunch(ic),ebfb,1,np,dt_loc,vb,Ltz)
  call curr_accumulate(bunch(ic),ebfb,jc,1,np,iform,n_st,xm,ym,zm)
 enddo
 call curr_mpi_collect(jc,i1,i2,j1,j2,k1,k2,curr_ndim)
! EXIT current density jb(1:3) due to bunch particles density and velocity at
! time t^{n+1/2}
! STEP3  advances fields
!=======================
 !======================= boundary ibx as for Maxwell equation
 if(ibeam > 0)then
  call advect_bunch_fields(ebf_bunch,jc,&
                          bet0,dt_loc,i1,i2,j1,j2,k1,k2,initial_time)
  !ebf_bunch(1:6)=[Ex,Ey,Ez,Jbx,By,Bz], Bx=0
 endif
 call advance_lpf_fields(ebf1_bunch,jc,dt_loc,vbeam,i1,i2,j1,j2,k1,k2,1)  !here reflecting bds
 !=========================
 end subroutine lpf2_eb_evolve
 !================================
 subroutine BUNCH_run(t_loc,dt_loc,iter_loc)

 real(dp),intent(in) :: t_loc,dt_loc
 integer,intent(in) :: iter_loc
 real(dp) :: ts
 real(dp),parameter :: eps=1.e-08
 logical :: init_time
 logical,parameter :: mw=.false.
 !+++++++++++++++++++++++++++++++++
 !for vbeam >0 uses the xw=(x+vbeam*t)
 !x=xi=(xw-vbeam*t) fixed
 ts=t_loc
 init_time=.false.
 if(t_loc<eps)init_time=.true.
 if(w_speed>0.0)then ! moves the computational box with w_speed>0.
  if(iter_loc==0)call BUNCH_window_xshift(dt_loc,w_sh,iter_loc)
  if(ts>=wi_time.and.ts< wf_time)then
   if(mod(iter_loc,w_sh)==0)then
    call BUNCH_window_xshift(dt_loc,w_sh,iter_loc)
   endif
  endif
 endif
 if(Comoving)then
  if(ts>=wi_time.and.ts< wf_time)then
   if(mod(iter_loc,w_sh)==0)then
    call comoving_coordinate(vbeam,dt_loc,w_sh,iter_loc)
   endif
  endif
 endif
 !=========================
  call lpf2_eb_evolve(dt_loc,iter_loc,init_time)
 !
 call cell_part_dist(mw)
 call cell_bpart_dist(mw)
 ts=ts+dt_loc
 end subroutine BUNCH_run
 !================================
 end module bunch_evolve
!==============================
