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
 !Common parameters defined in
 !   1) input.nml list (=> START/read_input.f90)
 !   2) initial conditions (=> START/set_init_param.f90, init)
 !   3) during time evolution
 !================================

 module common_param
  use phys_param
  use control_bunch_input

  implicit none

 integer,parameter :: Ref_nlayer=6,Ref_nlas=8
                                               !namelist input parameters
 integer :: nx,ny,nz,ny_targ
 real(dp) :: k0,yx_rat,zx_rat
 integer :: ibx,iby,ibz,ibeam
 integer :: LPf_ord,der_ord,str_flag,iform,model_id,dmodel_id
 integer :: nsp,nsb,ionz_lev,ionz_model,ion_min(3),ion_max(3),atomic_number(3)
 integer :: nb_laser,nb_1,np_per_xc(Ref_nlayer),np_per_yc(Ref_nlayer)
 real(dp) :: mass_number(3),t0_pl(4)
 real (dp) :: lpx(7),lpy(2),n_over_nc,np1,np2,r_c
 real(dp) :: t0_lp,xc_lp,tau_fwhm,w0_y,a0,lam0,lp_delay(Ref_nlas)
 real(dp) :: lp_offset,t1_lp,tau1_fwhm,w1_y,a1,lam1,a_symm_rat
 real(dp) :: xc_1,gam_1,sxb_1,syb_1,epsy_1,epsz_1,dg_1, charge_1,t_inject
 integer :: nouts,iene,nvout,nden,npout,nbout,jump,pjump
 integer :: new_sim,id_new,dump
 real(dp) :: gam_min,xp0_out,xp1_out,yp_out,tmax,cfl
 real(dp) :: w_speed,wi_time,wf_time
 integer :: tkjump,nkjump,track_tot_nstep
 real(dp) :: txmin,txmax,tymin, tymax, tzmin, tzmax,t_in, t_out
 integer :: nprocx,nprocy,nprocz 

                                          ! parameters set by initial
                                          ! conditiona
 logical :: G_prof,P_tracking,Comoving 
 logical :: Beam, Hybrid, Wake,Envelope,Solid_target
 logical :: Ionization,Ions
 logical :: Part,Stretch,Channel,Inject_beam
 logical :: Lp_active,Lp_inject,Plane_wave,Lin_lp,Circ_lp,Relativistic,Two_color
 logical :: Enable_ionization(2),Symmetrization_pulse
 logical :: Charge_cons,High_gamma,Test

 integer :: nx_loc,ny_loc,nz_loc,npty,nptz,nptx_max,ncmp_max,nx_alloc
 integer :: loc_npty(8),loc_nptz(8),nptx(8),loc_nptx(8),sptx_max(4),nxf,npt_buffer(4)
 integer :: sh_targ,nx_stretch,ny_stretch,nz_stretch
 integer :: mp_per_cell(6),nref,np_per_zc(Ref_nlayer),ppc(Ref_nlayer)
 integer :: loc_nyc_max,loc_nzc_max,loc_nxc_max,ndim_max
 real(dp) :: djc(3),ratio_mpc(6),pavg_npart(4),wgh_ion,dt
 real(dp) :: mass(4),mass_rat(4),charge_to_mass(4),unit_charge(4),Lorentz_fact(4)
 real(dp) :: n0_ref,pmass,ompe,vbeam,curr_max(3),j0_norm,ratio_mpfluid,chann_fact
 real(dp) :: gam0,bet0,u0_b, nb_over_np,b_charge

 real(dp) :: oml,E0,lp_pow,ZR,lp_intensity,lp_xsize,P_c
 real(dp) :: w0_x,lp_amp,xf,lp_max,eb_max,lp_energy,lp_rad
 real(dp) :: xc1_lp,xf1,ZR1,lp1_rad,lp1_amp,om1,w1_x

 real(dp) :: t0_b,el_lp,el_D,lambda_p,omega_p,lpvol
 real(dp) :: nc0,ncrit,n1_over_n,n2_over_n
 real(dp) :: np_per_cell,nb_per_cell,np_per_nmacro,nmacro
 real(dp) :: targ_in,targ_end,lx_fwhm
 real(dp) :: lp_in(Ref_nlas),lp_end(Ref_nlas),lp_ionz_in,lp_ionz_end,xf_loc(Ref_nlas),xc_loc(Ref_nlas)
 real(dp) :: y0_cent(Ref_nlas),z0_cent(Ref_nlas),y1_cent,z1_cent
 real(dp) :: ymin_t,ymax_t,zmin_t,zmax_t,rmin_t,rmax_t
 ! tracking param
 integer  :: track_tot_part

 integer :: pot_ndim,nb_max,pe_nbmax,nb_min,pe_nbmin
 integer :: Tsc_ord,t_ord,spl_ord
 integer :: nsp_run,nsp_ionz
 integer :: ndim,curr_ndim,nj_dim,nd2,nfield,nbfield,nfcomp,mod_ord,w_sh
 real(dp) :: macro_charge

 real(dp) :: tnow,tscale
 real(dp) :: energy_in_targ
 integer(kind=8) :: nptot_global


 end module common_param


