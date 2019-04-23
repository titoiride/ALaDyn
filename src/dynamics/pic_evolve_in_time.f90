
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

 module pic_evolve

 use window
 use boris_push
 use curr_and_fields_util
 use mpi_part_interface
 use init_grid_field
 use ionize
 use fluid_density_momenta

 implicit none
 !===============================
 contains
 !=======================
 subroutine lpf2_evolve(t_loc,dt_loc,iter_loc,initial_time)
 real(dp),intent(in) :: t_loc,dt_loc
 integer,intent(in) :: iter_loc
 logical,intent(in) :: initial_time
 integer :: lp,ic,np,i1,i2,j1,j2,k1,k2,n_st,id_ch
 real(dp) :: xm,ym,zm,Ltz,ef2_ion(1),loc_ef2_ion(1)
 !============================
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
 !====================
 call pfields_prepare(ebf,i1,i2,j1,j2,k1,k2,nfield,2,2)
 if(Ionization)then
  if(iter_loc==0)then
   call init_random_seed(mype)
  endif
  id_ch=nd2+1
  do ic=2,nsp_ionz
   np=loc_npart(imody,imodz,imodx,ic)
   if(np>0)then
    call set_ion_Efield(ebf,spec(ic),ebfp,np,n_st,ndim,nsp_run,dt_loc,xm,ym,zm)
    if(mod(iter_loc,100)==0)then     !refresh ionization tables, if needed
     loc_ef2_ion(1)=maxval(ebfp(1:np,id_ch))
     loc_ef2_ion(1)=sqrt(loc_ef2_ion(1))
     ef2_ion(1)=loc_ef2_ion(1)
     !if(prl)call allreduce_dpreal(MAXV,loc_ef2_ion,ef2_ion,1)
     if(ef2_ion(1) > lp_max)then
      lp_max=1.1*ef2_ion(1)
      call set_field_ioniz_wfunction(&
                     ion_min(ic-1),atomic_number(ic-1),ic,ionz_lev,ionz_model,lp_max,dt_loc)
     endif
    endif
    call ionization_cycle(spec(ic),ebfp,np,ic,iter_loc,0,de_inv)
   endif
    !======== injects new electrons. 
  end do
 endif
 !===================END IONIZATION MODULE============
  !    ions enter with new ionization levels and new electrons
  !                   are injected
 !=============================================
 jc(:,:,:,:)=0.0
 !curr_clean
  do ic=1,nsp_run
   np=loc_npart(imody,imodz,imodx,ic)
   Ltz=Lorentz_fact(ic)
   if(np >0)then
    !==============
    !============
    call set_lpf_acc(ebf,spec(ic),ebfp,np,ndim,nfield,n_st,xm,ym,zm)
    call field_charge_multiply(spec(ic),ebfp,1,np,nfield)
   
    if(initial_time)call init_lpf_momenta(spec(ic),ebfp,1,np,dt_loc,Ltz)
    call lpf_momenta_and_positions(spec(ic),ebfp,1,np,dt_loc,vbeam,Ltz)
    ! For each species :
    ! ebfp(1:3) store (X^{n+1}-X_n)=V^{n+1/2}*dt
    ! ebfp(4:7) store old x^n positions and dt/gam at t^{n+1/2}
    !
    call curr_accumulate(spec(ic),ebfp,jc,1,np,iform,n_st,xm,ym,zm)
    !================= only old ion charge saved
   endif
  enddo
  !==========================================
  if(Part)call curr_mpi_collect(jc,i1,i2,j1,j2,k1,k2,curr_ndim)
  !================ sums and normalize currents
 if(Hybrid)then
  call set_momentum_density_flux(up,flux,i1,i2,j1,j2,k1,k2)

  call update_adam_bash_fluid_variables(&
                        up,up0,flux,ebf,dt_loc,i1,i2,j1,j2,k1,k2,iter_loc,Ltz,initial_time)
  ! In flux(1:curr_ndim+1) are stored fluid (P,den) at t^{n+1/2}
  call fluid_curr_accumulate(flux,jc,dt_loc,i1,i2,j1,j2,k1,k2)
!=====================================
  ! In jc(1:3) exit total current density array Dt*(Jx,Jy,Jz)^{n+1/2}
 endif
  !=======================
  ! Inject fields at i=i1-1  for inflow Lp_inject=T
  call wave_field_left_inject(xm)  !(Bz=Ey By=Ez are injected at i1-1 point
  call advance_lpf_fields(ebf,jc,dt_loc,vbeam,i1,i2,j1,j2,k1,k2,0)
 !============================
 contains
 subroutine wave_field_left_inject(x_left)
  real(dp),intent(in) :: x_left
  real(dp) :: tnew
  integer :: wmodel_id,ic

  wmodel_id=model_id
  if(Plane_wave)wmodel_id=0
  tnew=t_loc    !Set inflow values [B_z{n}(i1-1/2) E_y{n}(i-1}
  Lp_inject=.false.
  do ic=1,nb_laser
   if(lp_in(ic) < x_left.and.lp_end(ic)>= xm)then
    Lp_inject=.true.
    if(model_id<3) call inflow_lp_fields(&
              ebf,lp_amp,tnew,t0_lp,w0_x,w0_y,xf_loc(ic),oml,wmodel_id,i1,j1,j2,k1,k2)
    if(model_id==3)call inflow_cp_fields(&
              ebf,lp_amp,tnew,t0_lp,w0_x,w0_y,xf_loc(ic),wmodel_id,i1,j1,j2,k1,k2)
   endif
   lp_in(ic)=lp_in(ic)+dt_loc
   lp_end(ic)=lp_end(ic)+dt_loc
  end do
  if(Two_color)then
   if(lp_ionz_in < x_left.and.lp_ionz_end >=xm)then
    Lp_inject=.true.
    call inflow_lp_fields(&
       ebf,lp1_amp,tnew,t1_lp,w1_x,w1_y,xf1,om1,model_id,i1,j1,j2,k1,k2)
   endif
   lp_ionz_in=lp_ionz_in+dt_loc
   lp_ionz_end=lp_ionz_end+dt_loc
  endif
 end subroutine wave_field_left_inject
 !-----------------------------
 end subroutine lpf2_evolve
 !===============
 ! END SECTION for Leap-frog one-cycle integrator in LP regime
 !===============
 subroutine LP_run(t_loc,dt_loc,iter_loc)

 real(dp),intent(in) :: t_loc,dt_loc
 integer,intent(in) :: iter_loc
 real(dp) :: ts
 logical :: init_time
 logical,parameter :: mw=.false.
 !+++++++++++++++++++++++++++++++++
 ts=t_loc
 init_time=.false.
 if(t_loc<1.e-08)init_time=.true.
 if(w_speed>0.0)then ! moves the computational box with w_speed>0.
  if(iter_loc==0)call LP_window_xshift(dt_loc,w_sh,iter_loc)
  if(ts>=wi_time.and.ts< wf_time)then
   if(mod(iter_loc,w_sh)==0)then
    call LP_window_xshift(dt_loc,w_sh,iter_loc)
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
 !==============================
 !vbeam=-w_speed
 !vbeam >0 uses the xw=(x+vbeam*t)
 !x=xi=(xw-vbeam*t) fixed
 !==============================
 !=========================
  call lpf2_evolve(t_loc,dt_loc,iter_loc,init_time)
!===================================
 if(Part)call cell_part_dist(mw)
 ts=ts+dt_loc
 end subroutine LP_run

 end module pic_evolve

