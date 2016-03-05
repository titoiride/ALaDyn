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

 module all_param

 use precision_def
 use mpi_var
 use phys_param
 use grid_and_particles
 use code_util
 use fftg_param
 use grid_param
 use control_bunch_input
 use ionize

 implicit none



 contains


 subroutine lpf_sch_coeff

 a_lpf(1:4)=1.
 b_lpf(1:4)=1.
 lpf_form(1:4)=iform
 if(LPf_ord==3)then
  lpf_form(1:3)=iform
  lpf_form(4)=0
  b_lpf(1)=1./(2.-2**(1./3.))  !alpha
  b_lpf(2)=1.-2.*b_lpf(1)      !beta=1.-2*alpha
  b_lpf(3)=b_lpf(1)            !alpha
  b_lpf(4)=0.0
  a_lpf(1)=0.5*b_lpf(1)        !alpha/2
  a_lpf(2)=0.5*(1.-b_lpf(1))   !1/2(1-alpha)
  a_lpf(3)=a_lpf(2)
  a_lpf(4)=a_lpf(1)
  c_lpf(1)=0.0
  c_lpf(2)=a_lpf(1)
  c_lpf(3)=0.5*b_lpf(2)
  c_lpf(4)=c_lpf(2)
 endif

 end subroutine lpf_sch_coeff

 !--------------------------

 subroutine rk_tsch_coeff(rk)
 integer,intent(in) :: rk
 integer :: ip

 select case(rk)
 case(3)
  b_rk(1)=1.
  b_rk(2)=0.25
  b_rk(3)=2./3.
  a_rk(1)=0.0
  a_rk(2)=1.-b_rk(2)
  a_rk(3)=1.-b_rk(3)

 case(4)
  rk_form(1:4)=iform
  if(iform==2)rk_form(4)=0
  a_rk(1)=1./6.
  a_rk(2)=1./3.
  a_rk(3)=1./3.
  a_rk(4)=a_rk(1)
  c_rk(0)=-4./3.
  c_rk(1)=1./3.
  c_rk(2)=2./3.
  c_rk(3)=c_rk(1)
  c_rk(4)=1.
  b_rk(1)=0.5
  b_rk(2)=0.5
  b_rk(3)=1.0
  b_rk(4)=1.0/6.0

 case(5)      !4th order five-stage optimized
  rk_form(1:5)=iform
  a_rk(1)=0.4
  a_rk(2)=0.0581597561
  a_rk(3)=1.0746983979
  a_rk(4)=-0.5338616876
  a_rk(5)=0.0

  b_rk(1)=-0.01
  b_rk(2)=0.4818402439
  b_rk(3)=-0.8615395999
  b_rk(4)=1.057293281
  b_rk(5)=0.3324060746

  c_rk(1)=0.0
  c_rk(2)=0.39
  c_rk(3)=0.53
  c_rk(4)=0.6349990419
  c_rk(5)=0.1337322374

  do ip=2,5
   c_rk(ip)=b_rk(ip-1)+c_rk(ip-1)
  end do

  do ip=2,5
   c_rk(ip)=a_rk(ip-1)+c_rk(ip)
  end do

 case(6)          !4th-order six-stage optimized
  b_rk(1)=0.10893125722541
  b_rk(2)=0.13201701492152
  b_rk(3)=0.38911623225517
  b_rk(4)=-0.59203884581148
  b_rk(5)=0.47385028714844
  b_rk(6)=0.48812405426094
  a_rk(1)=0.17985400977138
  a_rk(2)=0.14081893152111
  a_rk(3)=0.08255631629428
  a_rk(4)=0.65804425034331
  a_rk(5)=0.31862993413251
  a_rk(6)=0.
  c_rk(1)=0.0
  do ip=2,6
   c_rk(ip)=b_rk(ip-1)+c_rk(ip-1)
  end do
  do ip=2,6
   c_rk(ip)=a_rk(ip-1)+c_rk(ip)
  end do

 end select

 end subroutine rk_tsch_coeff

 !--------------------------

 subroutine param
 ! sets general parameters and grid depending on initial conditions
 integer :: i, sh_t
 real(dp) :: b_charge_den,gvol,gvol_inv,a_ch,k_ch,nm_fact

 a_lpf(1:4)=1.
 b_lpf(1:4)=1.
 if(LPf_ord >2)call rk_tsch_coeff(LPf_ord)

 nxcell=nx
 nycell=ny
 nzcell=nz

 ndim=1
 if(ny>1)ndim=2
 if(nz>1)ndim=3

 nx_loc=nx/nprocx
 nx1_loc=nx/nprocy
 ny_loc=ny/nprocy
 nz_loc=nz/nprocz
 sh_t=ny/2-ny_targ/2
 Lp_inject=.false.
 Cyl_coord=.false.
 Ionization=.false.
 Impact_ioniz=.false.
 !===============================
 ! Multispecies target with max 3 ionic species
 !=================
 Charge_cons=.true.
 if(iform>1)Charge_cons=.false.
 !==========================
 nsp_run=nsp
 !===============================
 do i=1, nsp-1                  !ion species
  if(ion_min(i)< ion_max(i))nsp_ionz=i+1
 enddo
 nsp_ionz=min(nsp_ionz,nsp)
 if(ionz_lev >0) Ionization=.true.
 if(Ionization) call set_ionization_coeff(atomic_number,nsp_ionz)
 !==========================
 Wake=.false.
 Solid_target=.false.
 if(n_over_nc < 0.5)Wake=.true.
 if(n_over_nc >1.)Solid_target=.true.
 !=================================
 PML=.false. ! FIX SERVE O LO RIMUOVIAMO?
 Stretch=.false.
 pml_size=0
 n_stretch=0
 nx_stretch=0
 if(str_flag==1)then
  Stretch=.true.
  n_stretch=nint(real(ny,dp)*size_of_stretch_along_y)
 endif
 if(str_flag==2)then
  Stretch=.true.
  n_stretch=nint(real(ny,dp)*size_of_stretch_along_y)
  nx_stretch=nint(real(nx,dp)*size_of_stretch_along_x)
 endif
 call grid_alloc(nx,nx_loc,ny,ny_loc,nz,nz_loc,ny_targ,&
  pml_size,nprocy,nprocz,nprocx)
 loc_nyc_max=loc_ygr_max
 loc_nzc_max=loc_zgr_max
 loc_nxc_max=loc_xgr_max

 pi2=2.0*pi
 test=.false.
 Part=.false.

 if(np_per_xc(1) > 0)Part=.true.

 np_per_cell=1
 cmp=.false.
 Lin_lp=.false.
 Relativistic=.false.
 Ions=.false.
 Beam=.false.
 Pbeam=.false.
 Envelope=.false.
 Circ_lp=.false.
 Plane_wave=.false.
 if(iby==2)Plane_wave=.true.

 ifilt=0
 if(der_ord >3)cmp=.true.
 t_ord=LPf_ord
 RK_ord=LPf_ord
 Lp_active=.false.

 ! mass-charge parameters four species: three ion species+ electrons
 ! Ions charges defined by initial conditions
 unit_charge(1)=electron_charge_norm
 unit_charge(2)=ion_min(1)*proton_charge_norm
 unit_charge(3)=ion_min(2)*proton_charge_norm
 unit_charge(4)=ion_min(3)*proton_charge_norm
 !=================================
 mass(1)=electron_mass_norm
 mass(2)=mass_number(1)*proton_mass_norm
 mass(3)=mass_number(2)*proton_mass_norm
 mass(4)=mass_number(3)*proton_mass_norm
 if(dmodel_id==0)then
  mass(2)=mass(1)
 endif
 do i=1,4
  mass_rat(i)=1./mass(i)
  charge_to_mass(i)=unit_charge(i)*mass_rat(i)
 end do
 !======================================
 Lorentz_fact(1:4)=mass_rat(1:4)  !to be multiplied by the particle charge
 !==========================================
 nm_fact=1.
 gvol_inv=1.
 gvol=1.
 j0_norm=1.
 jb_norm=1.
 E0=electron_mass !E=in unit mc^2/(e*l0)=[TV/m]/E0; E0*E in [TV/m]=[MV/mu] unit

 do i=1,6
  mp_per_cell(i)=np_per_xc(i)*np_per_yc(i)
  if(ndim==3) mp_per_cell(i)=np_per_xc(i)*np_per_yc(i)*np_per_zc(i)
 end do
 ratio_mpc=1.
 nref=mp_per_cell(1)
 if(mp_per_cell(1) >0)then
  j0_norm=1./real(nref,dp)
  do i=1,6
   ratio_mpc(i)=real(mp_per_cell(1),dp)/real(mp_per_cell(i),dp)
  enddo
 endif
 !================================
 if(w_speed < 0.0)then
  vbeam=-w_speed
  Comoving=.true.
 else
  Comoving=.false.
  vbeam=0.0
 endif

 if(model_id < 5)then
  mod_ord=1
  ifilt=0
  Relativistic=.true.
  part_dcmp=.false.
  nfield=3
  curr_ndim=2
  if(ndim > 2)then
   nfield=6
   curr_ndim=ndim
  endif
  dx=1./k0
  call set_grid(nx,ny,nz)
  dt=cfl*dx
  select case(ndim)
  case(1)
   dt=cfl/sqrt(dx_inv*dx_inv)
   gvol_inv=dx_inv*dx_inv*dx_inv
  case(2)
   dt=cfl/sqrt(dx_inv*dx_inv+dy_inv*dy_inv)
   gvol_inv=dx_inv*dy_inv*dy_inv
  case(3)
   dt=cfl/sqrt(dx_inv*dx_inv+dy_inv*dy_inv+dz_inv*dz_inv)
   gvol_inv=dx_inv*dy_inv*dz_inv
  end select
  gvol=1./gvol_inv
  Lp_active=.true.

  if(model_id < 3) Lin_lp = .true.
  if(model_id == 3) Circ_lp = .true.
  if(model_id ==4) then
   mod_ord=2
   Envelope=.true.
  endif
  !=======================
  ! Code Units for laser fields
  nc=pi/(rc0*lam0*lam0)    !critical density in units n0=10^21/cm^3=10^9/mu^3
  nm_fact=nc*(1.e+9)       ! critical density (1/mu^3)
  oml=pi2/lam0             !laser frequency in unit c/l0
  lp_amp=a0*oml
  lp_max=lp_amp
  nc0=oml*oml              !nc0=(2*pi/lam0)** 2
  ompe=nc0*n_over_nc       !squared adimensional plasma frequency :

  lp_intensity=1.37*(a0/lam0)*(a0/lam0)  !in units 10^18 W/cm^2
  if(Plane_wave)then        !plane LP wave
   lp_pow=0.01*lp_intensity*Ly_box*Lz_box  !in TW = 10^{-3}J/fs
   ZR=0.0
  else
   lp_rad=w0_y*sqrt(2.*log(2.))           !FWHM focal spot
   ZR=pi*w0_y*w0_y/lam0
   lp_pow=0.5*pi*lp_intensity*w0_y*w0_y   !in units 10^10 W
   if(ndim==2)lp_pow=0.5*dy*lp_intensity*w0_y   !in units 10^10 W
   lp_pow=0.01*lp_pow                     !in TW= 1.e-03[J/fs]  units
  endif
  lx_FWHM=2.*w0_x*acos(sqrt(0.5*sqrt(2.)))/pi !FWHM llength in microns
  tau_FWHM=lx_FWHM/speed_of_light
  lp_energy=1.e-03*tau_FWHM*lp_pow
  energy_in_targ=0.0
  el_lp=lam0
  chann_rad=0.0
  if(n_over_nc > 0.0)then
   el_lp=lam0/(pi2*sqrt(n_over_nc))
   el_D=t0_pl(1)*el_lp         !Debye length t0_pl(1)= V_T/c at t=0
   omega_p=1./el_lp
   k_ch=8.*13.5*1000.*n_over_nc/(17.*pi*pi)
   a_ch=(k_ch*lp_pow)**(1./3.)
   if(n_over_nc >1.)chann_rad=pi*lpx(3)*n_over_nc/(a0*lam0)
  endif
  bet0=0.0
  lambda_p=pi2*el_lp
  lpvol=el_lp*el_lp*el_lp

 endif


 if(model_id> 4)then !  e-Beams section
  select case(model_id)
   !========================================
   ! uses l0=1mu
   ! for fields E_u=GV/m = 1.e-03*mc^2/(l0*e*E0)
   ! n0= 10^18/cc =10^6/mu^3
   !========================================
   ! the electron radius r_c=rc0*10^{-8}mu
   !the squared adimensional plasma frequency on n0 density
   !   r_c*l0*l0*n0= 10^{-3}*rc0
   !omp^2=   4*pi*l0*l0*rc*n0=4*pi*rc0*1.e-03*E0/1.e-03
   !=============================
  case(5)
   nm_fact=1.e+06      !electron density[mu^{-3}] in the n0=10^18/cm^3 plasma
   nc=1.0
   Wake=.true.
   Solid_target=.false.
   !======================
   eb_max=0.0
   nc0=2.*E0*pi2*rc0
   !======================
   ompe=nc0*n_over_nc
   Lorentz_fact(1:4)=1.e-03*Lorentz_fact(1:4)/E0
   omega_p=0.02*sqrt(10.*pi*rc0*n_over_nc)
   el_lp=1./omega_p
   nb_over_np=rhob(1)
   lam0=1.             !Is the unit of length
   gam0=gam(1)
   u0_b=sqrt(gam0*gam0-1.)   !the beam x-momentum
   bet0=u0_b/gam0            !the beam velocity
   Lorentz_bfact(1:5)=Lorentz_fact(1)
   lambda_p=pi2*el_lp
   lpvol=el_lp*el_lp*el_lp
   Beam=.true.
   nfield=3
   nbfield=4
   if(ndim ==3)then
    nfield=6
    nbfield=nfield
   endif
   curr_ndim=ndim
   do i=1,5
    if(bunch_type(i)==3)Lorentz_bfact(i)=Lorentz_fact(1)/proton_mass_norm
   end do
   mod_ord=3
   dx=1./k0
   jb_norm=1.
   bvol=1.
   call set_grid(nx,ny,nz)
   dt=cfl/sqrt(dx_inv*dx_inv+dy_inv*dy_inv+dz_inv*dz_inv)
   do i=1,nsb
    bvol(i)=pi2*sqrt(pi2)*sxb(i)*syb(i)*syb(i)    !the bunch volume (mu^3)
    b_charge_den=rhob(i)*n_over_nc*nm_fact*e_charge  !pC/mu^3
    bcharge(i)=b_charge_den*bvol(i)  !bunch charge in [pC]
    Qbch(i)=rhob(i)*bvol(i)/lpvol
   end do
   gvol_inv=dx_inv*dy_inv*dz_inv
   gvol=1./gvol_inv
   if(L_particles)then
    do i=1,nsb
     nb_tot(i) = nint(gvol_inv* bvol(i)/j0_norm)
    end do
   endif
   do i=1,nsb
    jb_norm(i)=rhob(i)*gvol_inv*bvol(i)/(nb_tot(i)*j0_norm)
   end do
   b_charge=bcharge(1)

   !==================
  end select
 endif
 if(Wake)nsp_run=1  !only electrons running
 !============================
 !==============  in 2D dz=dy mp_per_cell =mp_x*mp_y,  mp_z=1
 if(mp_per_cell(1) >0)then
  nmacro=nref*gvol_inv
  np_per_cell=nm_fact*n_over_nc*gvol !here the real particle number per cell
  macro_charge=np_per_cell*e_charge  !real charge [in pC units]/cell
  !===========================================
  !all parameters to be multiplied by the macro weight= j0_norm=1/mp_per_cell
  !==========================
  ! j0_norm*np_percell = number of real particles for each
  ! macroparticle
  !========================
  nmacro=nref*real(nx*ny*nz,dp)
  ! here the total initial nmacro particles
 endif

 !===============================================================
 !charge per macroparticle: np_per_nmacro=nm_fact*n_over_nc/nmacro
 !Qbch=bvol/lpvol

 !np_per_nmacro= electron density over el_macro density
 !multiplied by nb_over_np gives the bunch elecron density over bunch macro
 !under the condition :np_per_cell is the same for plasma and bunch macro
 !-----------------------------
 !========================= driving beams parameters

 pot_ndim=0
 nd2=2*curr_ndim
 nj_dim=curr_ndim
 if(Envelope)nj_dim=nj_dim+1
 if(ibeam==2)then
  pot_ndim=4
  if(curr_ndim >2)pot_ndim=5
  nj_dim=nj_dim+1
 endif

 djc(1)=dx
 djc(2:3)=0.0
 if(ndim==2)then
  djc(1)=0.5*dx
  djc(2)=0.5*dy
 endif
 if(ndim==3)then
  djc(1)=dx/3.
  djc(2)=dy/3.
  djc(3)=dz/3.
 endif
 ymin_t=y(1)
 ymax_t=y(ny+1)
 zmin_t=z(1)
 zmax_t=z(nz+1)
 if(ndim >1)then
  ymin_t=y(sh_t+1)
  ymax_t=y(ny+1-sh_t)
 endif
 if(ndim >2 )then
  zmin_t=z(1+sh_t)
  zmax_t=z(nz+1-sh_t)
 endif
 end subroutine param
 !--------------------------
 subroutine set_grid(n1,n2,n3)
 integer,intent(in) :: n1,n2,n3
 integer :: i,ns1
 real(dp) :: yy,yyh,sm,spec

 aph=pi*0.4
 dxi=1.
 dyi=1.
 dzi=1.
 sx_rat=1.
 sy_rat=1.
 sz_rat=1.
 sm=0.0
 spec=0.0
 dx_inv=1.0/dx
 do i=1,n1+1
  x(i)=dx*real(i-1,dp)    !xminx(1)=0,.....,xmax-dx=x(nx)
  xh(i)=x(i)+0.5*dx
  dx1(i)=1.
  dx1h(i)=1.
 end do
 dxi=dx
 dxi_inv=dx_inv
 ns1=n1+1-nx_stretch
 if(nx_stretch >0)then
  dxi=aph/real(nx_stretch,dp)
  dxi_inv=1./dxi
  Lx_s=dx*dxi_inv
  sx_rat=dxi*dx_inv
  spec=x(ns1)
  do i=ns1,n1+1
   yy=dxi*real(i-ns1,dp)
   yyh=yy+dxi*0.5
   x(i)=spec+Lx_s*tan(yy)
   xh(i)=spec+Lx_s*tan(yyh)
   dx1h(i)=cos(yyh)*cos(yyh)
   dx1(i)=cos(yy)*cos(yy)
  end do
 endif
 str_xgrid%sind(1)=nx_stretch
 str_xgrid%sind(2)=ns1
 str_xgrid%smin=x(1)
 str_xgrid%smax=x(ns1)
 xw=x
 xmax=x(n1)
 if(ibx==2)xmax=x(n1+1)
 xmin=x(1)

 dy=1.
 dy_inv=1./dy
 dyi=dy
 dyi_inv=1./dy
 ymin=0.0
 ymax=0.0
 y=0.0
 yh=0.0
 dy1=1.
 dy1h=1.
 Ly_box=1.
 if(n2 > 1)then
  dy=yx_rat*dx
  dy_inv=1./dy
  dyi=dy
  dyi_inv=dy_inv
  do i=1,n2+1
   y(i)=dy*real(i-1-n2/2,dp)
   yh(i)=y(i)+0.5*dy
   dy1(i)=1.
   dy1h(i)=1.
  end do
  ns1=n2+1-n_stretch
  if(Stretch)then
   dyi=aph/real(n_stretch,dp)
   dyi_inv=1./dyi
   L_s=dy*dyi_inv
   sy_rat=dyi*dy_inv
   sm=y(n_stretch+1)
   spec=y(ns1)
   do i=1,n_stretch
    yy=dyi*real(i-1-n_stretch,dp)
    yyh=yy+0.5*dyi
    y(i)=sm+L_s*tan(yy)
    yh(i)=sm+L_s*tan(yyh)
    dy1h(i)=cos(yyh)*cos(yyh)
    dy1(i)=cos(yy)*cos(yy)
   end do
   do i=ns1,n2+1
    yy=dyi*real(i-ns1,dp)
    yyh=yy+dyi*0.5
    y(i)=spec+L_s*tan(yy)
    yh(i)=spec+L_s*tan(yyh)
    dy1h(i)=cos(yyh)*cos(yyh)
    dy1(i)=cos(yy)*cos(yy)
   end do
  endif
  str_ygrid%sind(1)=n_stretch
  str_ygrid%sind(2)=ns1
  str_ygrid%smin=y(n_stretch+1)
  str_ygrid%smax=y(ns1)
  ymin=y(1)
  ymax=y(n2+1)
  Ly_box=ymax-ymin

  r=y
  rh=yh
  dr1=dy1
  dr1h=dy1h
  if(Cyl_coord)then
   do i=1,n2+1
    r(i)=dy*real(i-1,dp)
    rh(i)=r(i)+0.5*dy
    dr1(i)=1.
    dr1h(i)=1.
   end do
  endif
 endif

 dz=1.
 dz_inv=1./dz
 dzi=dz
 dzi_inv=1./dz
 zmin=0.0
 zmax=0.0
 z=0.0
 zh=0.0
 dz1=1.
 dz1h=1.
 Lz_box=1.
 if(n3 > 1)then
  dz=yx_rat*dx
  dz_inv=1./dz
  do i=1,n3+1
   z(i)=dz*real(i-1-n3/2,dp)
   zh(i)=z(i)+0.5*dz
   dz1(i)=1.
   dz1h(i)=1.
  end do
  ns1=n3+1-n_stretch
  if(Stretch)then
   dzi=aph/real(n_stretch,dp)
   dzi_inv=1./dzi
   L_s=dz*dzi_inv
   sz_rat=dzi*dz_inv
   sm=z(n_stretch+1)
   spec=z(ns1)
   do i=1,n_stretch
    yy=dzi*real(i-1-n_stretch,dp)
    yyh=yy+0.5*dzi
    z(i)=sm+L_s*tan(yy)
    zh(i)=sm+L_s*tan(yyh)
    dz1h(i)=cos(yyh)*cos(yyh)
    dz1(i)=cos(yy)*cos(yy)
   end do
   do i=ns1,n3+1
    yy=dzi*real(i-ns1,dp)
    yyh=yy+dzi*0.5
    z(i)=spec+L_s*tan(yy)
    zh(i)=spec+L_s*tan(yyh)
    dz1h(i)=cos(yyh)*cos(yyh)
    dz1(i)=cos(yy)*cos(yy)
   end do
  endif
  str_zgrid%sind(1)=n_stretch
  str_zgrid%sind(2)=ns1
  str_zgrid%smin=sm
  str_zgrid%smax=spec
  zmin=z(1)
  zmax=z(n3+1)
  Lz_box=zmax-zmin
 endif
 !================
 end subroutine set_grid

 !--------------------------

 subroutine set_ftgrid(n1,n2,n3,ksh,lxbox,lybox)
 integer,intent(in) :: n1,n2,n3,ksh
 real(dp),intent(in) :: lxbox,lybox
 integer :: i
 real(dp) :: wkx,wky,wkz

 wkx=2.*acos(-1.)/lxbox !lxbox=x(n1+1)-x(1)
 wky=2.*acos(-1.)/lybox !lybox=y(n2+1)-y(1)

 wkz=wky

 allocate(aky(n2+2),akz(n3+2))
 allocate(sky(n2+2),skz(n3+2))
 allocate(ak2y(n2+2),ak2z(n3+2),ak2x(n1+1))
 allocate(akx(1:n1+1),skx(1:n1+1))

 ak2x=0.0
 select case(ksh)
 case(0)  ! staggered k-grid
  akx=0.0
  do i=1,n1/2
   akx(i)=wkx*(real(i,dp)-0.5)
   skx(i)=2.*sin(0.5*dx*akx(i))/dx
  end do
  aky=0.0
  ak2y=0.0
  if(n2>1)then
   do i=1,n2/2
    aky(i)=wky*(real(i,dp)-0.5)
    aky(n2+1-i)=-aky(i)
   end do
   do i=1,n2
    sky(i)=2.*sin(0.5*dy*aky(i))/dy
   end do
   ak2y=aky*aky
  endif
  akz=0.0
  ak2z=0.0
  if(n3 >1)then
   do i=1,n3/2
    akz(i)=wkz*(real(i,dp)-0.5)
    akz(n3+1-i)=-akz(i)
   end do
   do i=1,n3
    skz(i)=2.*sin(0.5*dz*akz(i))/dz
   end do
   ak2z=akz*akz
  endif

 case(1)    !standard FT k-grid
  do i=1,n1/2
   akx(i)=wkx*real(i-1,dp)
   akx(n1+2-i)=-akx(i)
  end do
  do i=1,n1+1
   skx(i)=2.*sin(0.5*dx*akx(i))/dx
  end do
  aky=0.0
  ak2y=0.0
  if(n2 > 1)then
   do i=1,n2/2
    aky(i)=wky*real(i-1,dp)
    aky(n2+2-i)=-aky(i)
    sky(i)=2.*sin(0.5*dy*aky(i))/dy
   end do
   ak2y=aky*aky
  endif
  akz=0.0
  ak2z=0.0
  if(n3 > 1)then
   do i=1,n3/2
    akz(i)=wkz*real(i-1,dp)
    akz(n3+2-i)=-akz(i)
   end do
   do i=1,n3
    skz(i)=2.*sin(0.5*dz*akz(i))/dz
   end do
   ak2z=akz*akz
  endif

 case(2)  ! for the sine/cosine transform
  wkx=acos(-1.0)/lxbox
  wky=acos(-1.0)/lybox
  wkz=wky
  do i=1,n1+1
   akx(i)=wkx*real(i-1,dp)
   skx(i)=2.*sin(0.5*dx*akx(i))/dx
  end do
  aky=0.0
  ak2y=0.0
  if(n2>1)then
   do i=1,n2+1
    aky(i)=wky*real(i-1,dp)
    sky(i)=2.*sin(0.5*dy*aky(i))/dy
   end do
   ak2y=aky*aky
  endif
  akz=0.0
  ak2z=0.0
  if(n3 >1)then
   do i=1,n3+1
    akz(i)=wkz*real(i-1,dp)
    skz(i)=2.*sin(0.5*dz*akz(i))/dz
   end do
   ak2z=akz*akz
  endif
 end select
 end subroutine set_ftgrid

 !--------------------------


 end module all_param


 !--------------------------
