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

 module set_init_param

  use code_util
  use common_param
  use set_grid_param
  use ionz_data
  use control_bunch_input, only: lorentz_bfact, gam, jb_norm, nb_tot, &
   bunch_volume, syb, rhob, sxb, bunch_shape, charge_right, ppc_bunch, &
   charge_left, l_particles, bunch_charge, reduced_charge, bunch_type

  implicit none

 contains

!--------------------------

  subroutine set_initial_param


! sets general parameters and grid depending on initial conditions
   integer :: i
   real (dp) :: bunch_charge_density, gvol, gvol_inv, nm_fact, ncell
   real (dp) :: aph_fwhm, c1_fact, c2_fact
!================== grid dimension definition ==============
   ndim = 1
   if (ny>1) ndim = 2
   if (nz>1) ndim = 3
!===========time integrator schemes
   t_ord = lpf_ord
   if (t_ord==4) der_ord = t_ord
!===================================
   nx_loc = nx/nprocx
   ny_loc = ny/nprocy
   nz_loc = nz/nprocz
   sh_targ = ny/2 - ny_targ/2
   nx_stretch = 0
   stretch = .false.
   ny_stretch = 0
   if (ndim<2) str_flag = 0
   if (model_id>4) str_flag = 0
!===================
   if (str_flag==1) then !size_of_stretch defined in phys_param module => common_param
    stretch = .true.
    ny_stretch = nint(real(ny,dp)*size_of_stretch_along_y) !set to ny/6 
   end if
   if (str_flag==2) then
    stretch = .true.
    ny_stretch = nint(real(ny,dp)*1.5*size_of_stretch_along_y) !set to ny/4
   end if
   loc_nyc_max = loc_ygr_max
   loc_nzc_max = loc_zgr_max
   loc_nxc_max = loc_xgr_max
! i========Sets grid data=============
   gvol_inv = 1.
   gvol = 1.
   dx = 1./k0
   call set_grid(nx, ny, nz, ibx, nx_stretch, ny_stretch, k0, yx_rat, &
     zx_rat)
   dt = cfl*dx
   select case (ndim)
   case (1)
    dt = cfl/sqrt(dx_inv*dx_inv)
    gvol_inv = dx_inv*dx_inv*dx_inv
   case (2)
    dt = cfl/sqrt(dx_inv*dx_inv+dy_inv*dy_inv)
    gvol_inv = dx_inv*dy_inv*dy_inv
   case (3)
    dt = cfl/sqrt(dx_inv*dx_inv+dy_inv*dy_inv+dz_inv*dz_inv)
    gvol_inv = dx_inv*dy_inv*dz_inv
   end select
   gvol = 1./gvol_inv
   djc(1) = dx
   djc(2:3) = 0.0
   if (ndim==2) then
    djc(1) = 0.5*dx
    djc(2) = 0.5*dy
   end if
   if (ndim==3) then
    djc(1) = dx/3.
    djc(2) = dy/3.
    djc(3) = dz/3.
   end if
   ymin_t = y(1)
   ymax_t = y(ny+1)
   zmin_t = z(1)
   zmax_t = z(nz+1)
   if (ndim>1) then
    ymin_t = y(sh_targ+1)
    ymax_t = y(ny+1-sh_targ)
   end if
   if (ndim>2) then
    zmin_t = z(1+sh_targ)
    zmax_t = z(nz+1-sh_targ)
   end if
!======================================
   hybrid = .false.
   high_gamma = .false.
   if (gam_min>1.0) high_gamma = .true.
   test = .false.
   part = .false.
   lp_active = .false.
   lp_inject = .false.
   ionization = .false.
   charge_cons = .false.
   if (iform>1) charge_cons = .true.
   inject_beam = .false.
   beam = .false.
   relativistic = .false.
   ions = .false.
   envelope = .false.
   circ_lp = .false.
   plane_wave = .false.
   Two_color = .false.
   wake = .false.
   solid_target = .false.
   channel = .false.
   nm_fact = 1.
   if (np_per_xc(1)>0 .and. np_per_yc(1)>0) part = .true.
   if (nsp>1) part = .true.
   if (nsp>1) ions = .true.
   if (model_id<5) then
    lp_active = .true.
   else
    beam = .true.
   end if
   if (ibeam==2) hybrid = .true.
!===============================
! Multispecies target with max 3 ionic species (nsp=4)
!=================
!========================== particle on cells
   do i = 1, 6
    mp_per_cell(i) = np_per_xc(i)*np_per_yc(i)
    if (ndim==3) mp_per_cell(i) = np_per_xc(i)*np_per_yc(i)*np_per_zc(i)
   end do
   j0_norm = 1.
   nref = mp_per_cell(1)
   ratio_mpc = 1.
   if (mp_per_cell(1)>0) then
    j0_norm = 1./real(nref, dp)
    do i = 1, 6
     ratio_mpc(i) = 0.0
     if (mp_per_cell(i)>0) ratio_mpc(i) = real(mp_per_cell(1), dp)/ &
       real(mp_per_cell(i), dp)
    end do
   end if
   ratio_mpfluid = 1.0
   if (hybrid) then
    ratio_mpfluid = 0.1
    if (mp_per_cell(1)==0) ratio_mpfluid = 0.0
    if (ny_targ==0) ratio_mpfluid = 0.0
   end if
   j0_norm = j0_norm*ratio_mpfluid
!========================== multispecies
! mass-charge parameters four species: three ion species+ electrons
! Ions charges defined by initial conditions.
   unit_charge(1) = electron_charge_norm
   unit_charge(2) = ion_min(1)*proton_charge_norm
   unit_charge(3) = ion_min(2)*proton_charge_norm
   unit_charge(4) = ion_min(3)*proton_charge_norm
!=================================
   mass(1) = electron_mass_norm
   mass(2) = mass_number(1)*proton_mass_norm
   mass(3) = mass_number(2)*proton_mass_norm
   mass(4) = mass_number(3)*proton_mass_norm
   do i = 1, 4
    mass_rat(i) = 1./mass(i)
    charge_to_mass(i) = unit_charge(i)*mass_rat(i)
   end do
   lorentz_fact(1:4) = mass_rat(1:4) !for (E,B) fields in mc2/e/l0(mu) unit
   e0 = electron_mass
!=====================================
! Check species molecular number
!=====================================
   if (nsp>1) then
    do i = 1, nsp - 1
     call set_atoms_per_molecule(atomic_number(i), n_mol_atoms(i))
    end do
   end if
!======================================
! Set ionization param
! WARNING  ion ordering in input data must set ionizing species first
!==========================================================
!wgh_ion=j0_norm
   nsp_ionz = 1
   if (nsp>1) then
    do i = 1, nsp - 1 !index of ionizing ion species 2,.. nsp_ionz <= nsp
     if (ion_min(i)<ion_max(i)) nsp_ionz = i + 1
    end do
    nsp_ionz = min(nsp_ionz, nsp)
    do i = 1, 3
     if (mass_number(i)<1.) call set_atomic_weight(atomic_number(i), &
       mass_number(i))
    end do
    if (ionz_lev>0) ionization = .true.
    if (ionization) call set_ionization_coeff(atomic_number, nsp_ionz)
!uses ion index 1,2,,,nsp-1
    wgh_ion = concentration(1)/(real(mp_per_cell(2),dp))
!if(ion_min(1)>1)wgh_ion=1./(real(ion_min(1),dp)*real(mp_per_cell(2),dp))
   end if
!=======Moving window
   if (w_speed<0.0) then
    vbeam = -w_speed
    comoving = .true.
   else
    comoving = .false.
    vbeam = 0.0
   end if
!==========================
! Target parameters  enter n_over_nc
   if (lp_active) then
    if (iby==2) plane_wave = .true.
    mod_ord = 1
    if (model_id<3) lin_lp = .true.
    if (model_id==3) circ_lp = .true.
    if (model_id==4) then
     mod_ord = 2
     envelope = .true.
    end if
    if (n_over_nc>1.) then
     solid_target = .true.
    else
     wake = .true.
    end if
    relativistic = .true.
    nfield = 3
    curr_ndim = 2
    nbfield = 4
    if (ndim>2) then
     nfield = 6
     curr_ndim = ndim
     nbfield = 6
    end if
    if (circ_lp) then
     nfield = 6
     curr_ndim = 3
    end if
    if (lp_offset>0.0) Two_color = .true.
!====================
!to be multiplied by the particle charge in the equation of motion
!E=in unit mc^2/(e*l0)=[TV/m]/E0; E0*E in [TV/m]=[MV/mu] unit
!E0(A,phi) in MV unit
!==========================================
    np_per_cell = 1
!====================================
!=======================
! Code Units for laser fields
    if (wake) then
     n_plasma = 0
     if (nsp>1) then
      do i = 1, nsp - 1
       n_plasma = n_plasma + concentration(i)*n_mol_atoms(i)*ion_min(i)
      end do
      if (n_plasma<epsilon) then
       np_per_xc(1) = 0
       np_per_yc(1) = 0
       n_plasma = zero_dp
      end if
     else if (nsp==1) then
      n_plasma = one_dp
     end if
    end if
    ncrit = pi/(rc0*lam0*lam0) !critical density in units n0=10^21/cm^3=10^9/mu^3
    n0_ref = 1.e03*ncrit*n_over_nc ! reference density in 10^18/cc=10^6/mu^3 unit
    nm_fact = ncrit*(1.e+9) ! critical density (1/mu^3)
!n_over_nc*nm_fact = n0_ref*1.e+06=  background density in [1/mu^3]
!===================================     
    oml = pi2/lam0 !laser frequency in unit c/l0
    om1 = pi2/lam1 !laser frequency in unit c/l0
    lp_amp = a0*oml
    lp1_amp = a1*om1 !field in unit 0.51 MV/m
    lp_max = 1.5*lp_amp
    if (Two_color) lp_max = max(lp_max, 1.5*lp1_amp)
!=============================
    nc0 = oml*oml !nc0=(2*pi/lam0)** 2
    ompe = nc0*n_over_nc !squared adimensional plasma frequency :
!===============================
! Parabolic plasma channel profile const  r_c=w0_y matched condition
    chann_fact = 0.0
    if (r_c>0.0) then
     channel = .true.
     c1_fact = w0_y*w0_y/(r_c*r_c)
     c2_fact = lam0*lam0/(r_c*r_c)
     chann_fact = c1_fact*c2_fact/(pi*pi*n_over_nc)
    end if
!========== Laser parameters
    lp_intensity = 1.37*(a0/lam0)*(a0/lam0) !in units 10^18 W/cm^2
    lp_rad = w0_y*sqrt(2.*log(2.)) !FWHM focal spot
    zr = pi*w0_y*w0_y/lam0
    lp1_rad = w1_y*sqrt(2.*log(2.)) !FWHM focal spot
    zr1 = pi*w1_y*w1_y/lam1
    lp_pow = 0.5*pi*lp_intensity*w0_y*w0_y !in units 10^10 W
    if (ndim==2) lp_pow = 0.5*dy*lp_intensity*w0_y !in units 10^10 W
    lp_pow = 0.01*lp_pow !in TW= 1.e-03[J/fs]  units
    if (plane_wave) then !plane LP wave
     lp_pow = 0.01*lp_intensity*ly_box*lz_box !in TW = 10^{-3}J/fs
     zr = 0.0
    end if
!======= Defines scales (w0_x,w1_x) in longitudinal (t-x) pulse profile
    if (g_prof) then
     aph_fwhm = sqrt(2.*log(2.))
    else
     aph_fwhm = 2.*acos(sqrt(0.5*sqrt(2.)))/pi
!aph_fwhm=2.*acos(sqrt(sqrt(0.5*sqrt(2.))))/pi this is only valid for cos^4 pulse shape
    end if
    lx_fwhm = tau_fwhm*speed_of_light ! In micron unit
    w0_x = lx_fwhm/aph_fwhm
    w1_x = speed_of_light*tau1_fwhm/aph_fwhm
!=================
    lp_energy = 1.e-03*tau_fwhm*lp_pow
    energy_in_targ = 0.0
    el_lp = lam0
    if (n_over_nc>0.0) then
     p_c = 0.0174/n_over_nc !Critical power in TW
     lambda_p = lam0/sqrt(n_over_nc)
     el_lp = lambda_p/pi2
     el_d = t0_pl(1)*el_lp !Debye length t0_pl(1)= V_T/c at t=0
     omega_p = 1./el_lp
    end if
    bet0 = 0.0
    lpvol = el_lp*el_lp*el_lp
    if (nsb>0) inject_beam = .true.
!=====================
    if (inject_beam) then
!ON input phase space coordinates, beam size, 
!         total macro-particle number nb_tot(1), total charge (pC)
!====================================
     gam0 = gam(1) !the initial gamma factor
     u0_b = sqrt(gam0*gam0-1.) !the beam x-momentum
     bet0 = u0_b/gam0 !the beam velocity
!==================
     b_charge = nb_tot(1)*e_charge !the charge of bunch macro-particle
     np_per_nmacro = bunch_charge(1)/b_charge !real particles/macro particles=real charge/macro charge
!===============
     if (ndim<3) then
      bunch_volume(1) = pi2*sxb(1)*syb(1)*dy !the bunch volume (mu^3) in 2D Gaussian
     else
      bunch_volume(1) = pi2*sqrt(pi2)*sxb(1)*syb(1)*syb(1) !the bunch volume (mu^3) in 3D Gussian bunch
     end if
     rhob(1) = bunch_charge(1)/(e_charge*bunch_volume(1)) !bunch density (1/mu^3)

     n0_ref = nm_fact*n_over_nc ! background plasma density (1/mu^3)
     rhob(1) = rhob(1)/n0_ref !ratio beam density/background plasma density
!----------------------------
     ncell = bunch_volume(1)*gvol_inv
     nb_per_cell = nb_tot(1)/ncell
     jb_norm(1) = rhob(1)/nb_per_cell !
    end if
   end if !EMD Laser-driven section
!===================================
   if (beam) then !  e-Beams section
!========================================
! uses l0=1mu
! for fields E_u=GV/m = 1.e-03*mc^2/(l0*e*E0)  (A,phi) in kVolt unit
! n0= 10^18/cc =10^6/mu^3
!========================================
! the electron radius r_e=rc0*10^{-8}mu
!the squared adimensional plasma frequency on n0 density
!   r_e*l0*l0*n0= 10^{-3}*rc0
!omp^2=   4*pi*l0*l0*rc*n0=4*pi*rc0*1.e-03*E0/
!=============================
    nm_fact = 1.e+06 !electron density[mu^{-3}] in the n0=10^18/cm^3 plasma
    ncrit = 1.0
    wake = .true.
    solid_target = .false.
    n0_ref = n_over_nc
!======================
    eb_max = 0.0
    nc0 = 2.*e0*pi2*rc0
!======================
    ompe = nc0*n_over_nc
!======================
    jb_norm = 1.
    b_charge = 1.
    gam0 = 1.
    lorentz_fact(1:4) = 1.e-03*lorentz_fact(1:4)/e0
    omega_p = 0.02*sqrt(10.*pi*rc0*n0_ref)
    el_lp = 1./omega_p
    nb_over_np = rhob(1)
    lam0 = 1. !Is the unit of length
    gam0 = gam(1)
    u0_b = sqrt(gam0*gam0-1.) !the beam x-momentum
    bet0 = u0_b/gam0 !the beam velocity
    lorentz_bfact(1:5) = lorentz_fact(1)
    lambda_p = pi2*el_lp
    lpvol = el_lp*el_lp*el_lp
    nfield = 3
    nbfield = 4
    if (ndim==3) then
     nfield = 6
     nbfield = nfield
    end if
    curr_ndim = ndim
    do i = 1, 5
     if (bunch_type(i)==3) lorentz_bfact(i) = lorentz_fact(1)/ &
       proton_mass_norm
    end do
    mod_ord = 3
    bunch_volume = 1.
    do i = 1, nsb
!bunch_volume(i)=pi2*sqrt(pi2)*sxb(i)*syb(i)*syb(i)    !the bunch volume (mu^3) !ciao
!bunch_charge_density=rhob(i)*n_over_nc*nm_fact*e_charge  !pC/mu^3
!bunch_charge(i)=bunch_charge_density*bunch_volume(i)  !bunch charge in [pC]
!reduced_charge(i)=rhob(i)*bunch_volume(i)/lpvol

!---------------------------!
     if (bunch_shape(i)==1) then !--- nomi da verificare e aggiungere
      if (ndim==3) then
       bunch_volume(i) = pi2*sqrt(pi2)*sxb(i)*syb(i)*syb(i) !the bunch volume (mu^3) in 3D Gussian
      else
       bunch_volume(i) = pi2*sxb(i)*syb(i)*dy !the bunch volume (mu^3) in 2D Gaussian
      end if
      bunch_charge_density = rhob(i)*n_over_nc*nm_fact*e_charge !pC/mu^3
      bunch_charge(i) = bunch_charge_density*bunch_volume(i)
      reduced_charge(i) = rhob(i)*bunch_volume(i)/lpvol
     end if

     if (bunch_shape(i)==2) then
      bunch_volume(i) = pi*syb(i)*syb(i)*sxb(i)
      bunch_charge_density = (charge_left(i)+charge_right(i))/2.0* &
        n_over_nc*nm_fact*e_charge !pC/mu^3
      bunch_charge(i) = bunch_charge_density*bunch_volume(i)
      reduced_charge(i) = (charge_left(i)+charge_right(i))/2.0* &
        bunch_volume(i)/lpvol
     end if

     if (bunch_shape(i)==3) then
      bunch_volume(i) = pi2*syb(i)*syb(i)*sxb(i)
      bunch_charge_density = (charge_left(i)+charge_right(i))/2.0* &
        n_over_nc*nm_fact*e_charge !pC/mu^3
      bunch_charge(i) = bunch_charge_density*bunch_volume(i)
      reduced_charge(i) = (charge_left(i)+charge_right(i))/2.0* &
        bunch_volume(i)/lpvol
     end if

     if (bunch_shape(i)==4) then
      bunch_volume(i) = pi*syb(i)*syb(i)*sxb(i)
      bunch_charge_density = rhob(i)*n_over_nc*nm_fact*e_charge !pC/mu^3
      bunch_charge(i) = bunch_charge_density*bunch_volume(i)
      reduced_charge(i) = rhob(i)*bunch_volume(i)/lpvol
     end if

    end do
    if (l_particles) then
     do i = 1, nsb
      nb_tot(i) = nint(gvol_inv*bunch_volume(i)/j0_norm)
     end do
    end if
    do i = 1, nsb
     if (bunch_shape(i)==1 .and. nb_tot(i)>0) then
      jb_norm(i) = rhob(i)*gvol_inv*bunch_volume(i)/(nb_tot(i)*j0_norm)
     end if
     if (bunch_shape(i)==1 .and. nb_tot(i)==-1) then
      jb_norm(i) = 1.0_dp/(product(ppc_bunch(i,:))*j0_norm)
     end if
     if (bunch_shape(i)==2 .and. nb_tot(i)>0) then
      jb_norm(i) = (charge_left(i)+charge_right(i))/2.0*gvol_inv* &
        bunch_volume(i)/(nb_tot(i)*j0_norm)
     end if
     if (bunch_shape(i)==1 .and. nb_tot(i)==-1) then
      jb_norm(i) = 1.0_dp/(product(ppc_bunch(i,:))*j0_norm)
     end if
     if (bunch_shape(i)==3 .and. nb_tot(i)>0) then
      jb_norm(i) = (charge_left(i)+charge_right(i))/2.0*gvol_inv* &
        bunch_volume(i)/(nb_tot(i)*j0_norm)
     end if
     if (bunch_shape(i)==3 .and. nb_tot(i)==-1) then
      jb_norm(i) = 1.0_dp/(product(ppc_bunch(i,:))*j0_norm)
     end if
     if (bunch_shape(i)==4 .and. nb_tot(i)>0) then
      jb_norm(i) = rhob(i)*gvol_inv*bunch_volume(i)/(nb_tot(i)*j0_norm)
     end if
     if (bunch_shape(i)==1 .and. nb_tot(i)==-1) then
      jb_norm(i) = 1.0_dp/(product(ppc_bunch(i,:))*j0_norm)
     end if
    end do

!--- I am forcing this part to be again with the correct input values ---!
    do i = 1, nsb
     if (ppc_bunch(i,1)>0) nb_tot(i) = -1
    end do
!--- *** ---!

    b_charge = bunch_charge(1)
   end if
!============================
!  SET PARAM all cases
   if (hybrid) nfcomp = curr_ndim + 1
   nsp_run = nsp
   if (wake) then
    nsp_run = 1 !only electrons running
    if (dmodel_id==3) wgh_ion = np1*wgh_ion
   end if
!==============  in 2D dz=dy mp_per_cell =mp_x*mp_y,  mp_z=1
   np_per_cell = nm_fact*n_over_nc*gvol !here the real particle number per cell
!===========================
   if (mp_per_cell(1)>0) then
    nmacro = nref*gvol_inv
    macro_charge = np_per_cell*e_charge !real charge [in pC units]/cell
!===========================================
!all parameters to be multiplied by the macro weight= j0_norm=1/mp_per_cell
!==========================
! j0_norm*np_percell = number of real particles for each
! macroparticle
!========================
    nmacro = nref*real(nx*ny*nz, dp)
! here the total initial nmacro particles
   end if
!========================= Memory allocation
   nx_alloc = nint(dx_inv*sum(lpx(1:5)))
   nx_alloc = min(nx_loc, nx_alloc)
   npt_buffer = 0
   do i = 1, nsp
    npt_buffer(i) = nx_alloc*ny_loc*nz_loc*mp_per_cell(i)
   end do
!===============================
   call set_ftgrid(nx, ny, nz)
!===============================================================
!density per macroparticle: np_per_nmacro=nm_fact*n_over_nc/nmacro

!np_per_nmacro= electron density over el_macro density
!multiplied by nb_over_np gives the bunch elecron density over bunch macro
!under the condition :np_per_cell is the same for plasma and bunch macro
!-----------------------------
!========================= driving beams parameters
   pot_ndim = 0
   nd2 = 2*curr_ndim
   nj_dim = curr_ndim
   if (envelope) nj_dim = curr_ndim + 1
   if (beam) then
    if (ibeam==2) then
     nj_dim = curr_ndim + 1
     pot_ndim = nj_dim + 1
    end if
   end if


  end subroutine
!--------------------------
 end module

!--------------------------
