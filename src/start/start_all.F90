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

 module start_all

 use array_alloc
 use set_init_param
 use read_input
 use set_grid_param
 use ionize
 use pic_in
 use pic_dump
 use run_data_info
 use system_utilities

 implicit none

 contains

 subroutine start

 integer :: ic,nxp,nyp,nzp,ncmp


 !enable loop to attach with gdb only if really needed
 !WARNING if enabled with no need, the program sleeps at start without doing anything!
 !To enable the flag, uncomment the corresponding line in CMakeLists.txt
#ifdef ENABLE_GDB_ATTACH
 call gdbattach
#endif

 !Read parameters from input.nml file
 call read_main_input

 call check_grid_size

 call set_initial_param      !sets parameters related to initial condition
                             !call set_grid() to define global grid and grid
                             !parameters
!=============================
 call start_parallel(nd2,nsp,nsb)
 !
 if(mpi_err >0)then
  if(pe0)write(6,*)' ERROR in mpi domain decomposition'
  call end_parallel
  stop
 endif
  !=== Ascii art generated on http://patorjk.com/software/taag using the Star Wars font ===
 if(pe0)then
  write(6,*)"     ___       __           ___       _______  ____    ____ .__   __.  ____         _____" 
  write(6,*)"    /   \     |  |         /   \     |       \ \   \  /   / |  \ |  | |__  \       /  _  \ "
  write(6,*)"   /  ^  \    |  |        /  ^  \    |  .--.  | \   \/   /  |   \|  |    )  |     |  | |  |" 
  write(6,*)"  /  /_\  \   |  |       /  /_\  \   |  |  |  |  \_    _/   |  . `  |   /  /      |  | |  |"
  write(6,*)" /  _____  \  |  `----. /  _____  \  |  '--'  |    |  |     |  |\   |  /  /_   __ |  |_|  |" 
  write(6,*)"/__/     \__\ |_______|/__/     \__\ |_______/     |__|     |__| \__| |_____| (__) \_____/ "
  write(6,*)"                                                                                        "
 endif                                                                             
 if(pe0)then
  write(6,*)'==========================================================================================='
  write(6,'(a52,i1,a1,i2,a36)') ' =                                  Code version    ',major_version,'.'&
  &,minor_version,'                                   ='
  write(6,*)'==========================================================================================='
  call create_initial_folders
  call write_read_nml
 endif
!======================================
!                            
 call mpi_loc_grid(nx_loc,ny_loc,nz_loc,&
                     nprocx,nprocy,nprocz)
                                          !Exit
                                          !loc_xgrid(nprocx),loc_ygrid(nprocy),loc_ygrid(nprocz) local grid data
 call set_fyzxgrid(npe_yloc,npe_zloc,npe_xloc,sh_ix)  !local grid parameters and
                                                      !coordinate struct  loc_xg,loc_yg,loc_zg
 if (Stretch) call set_str_ind(npe_yloc,npe_zloc,ndim)
 !---------------------------------
 nyp=loc_ygrid(imody)%p_ind(2)  !Ny_loc+2
 nzp=loc_zgrid(imodz)%p_ind(2)  !Nz_loc+2
 nxp=loc_xgrid(imodx)%p_ind(2)  !Nx_loc+2
!======================================
 call set_output_grid(jump,nprocx,nprocy,nprocz)   
                        !defines (nhx(nprocx), nhy(nprocy),nhz(nprocz) arrays of grid points
                        ! for output data
                        !allocates wdata() and gwdata()
 !=====================
 ! Allocates basic arrays, defines grid parameters, boundary index etc
 call set_field_param              !local arrays and coefficients for space derivatives
 mem_size=0
 mem_psize=0
 if (nvout>nfield) nvout=nfield
 !============
 !     Extended local grid
 !====== Fields and current arrays allocated on [1: N_loc+5]
 ncmp=nfield
 !==========================
 call v_alloc(nxp,nyp,nzp,nfield,nj_dim,&
             ndim,ibeam,LPf_ord,der_ord,Envelope,Two_color,Comoving,mem_size)
 if(Hybrid)then
  call fluid_alloc(nxp,nyp,nzp,nfcomp,ndim,LPF_ord,mem_size)
  ncmp=max(ncmp,nfcomp)
 endif
 if (Beam) then
  call bv_alloc(nxp,nyp,nzp,nbfield,ndim,ibeam,mem_size)
 endif
 call mpi_buffer_alloc(nx_loc,ny_loc,nz_loc,ncmp)
 !====To activate/disactivate diagnostics in envar() en_data routines=========
 Diag=.true.
 if(iene==0)then
  Diag = .false.
  iene=1
 endif
 Tpart=.false.
 inject_ind=-1
 if(Inject_beam)then
  inject_ind=nint(t_inject/dt)
 endif
!========================
 if(Ionization)then
  do ic=2,nsp_ionz
   call set_field_ioniz_wfunction(ion_min(ic-1),atomic_number(ic-1),ic,ionz_lev,ionz_model,lp_max,dt)
  end do
  if(Pe0) call Ioniz_data(lp_max,ion_min,atomic_number,ionz_lev,ionz_model)
 endif

 select case(new_sim)


 case (0) ! Set initial conditions, allocate fields and particles
  iout=id_new
  ienout=0
  call init
  tstart=0.0
  last_iter=0
  tdia=tstart
  tout=tstart
  ! to count outputs in energy-data (iene+1 times)
  ! in general data (nouts+1 times)
  dt_loc=dt
  iter_max=1
  dtout=(tmax-tstart)/nouts
  dtdia=(tmax-tstart)/iene
  if(tmax >0.0)then
   iter_max=int(tmax/dt)
   dt_loc=tmax/float(iter_max)
  endif

 case (1) ! reads from dump evolved data
  if (.not.L_first_output_on_restart) then
   iout=id_new
   ienout=0
  else
   iout=id_new+1
   ienout=0
  endif
  call restart(last_iter,tstart)
  call mpi_barrier(comm,error)
!=============================
  call set_fxgrid(npe_xloc,sh_ix)
  if(tmax >0.0)then
   iter_max=int(tmax/dt)
   dt_loc=tmax/float(iter_max)
  endif
  dtout=tmax/nouts
  dtdia=tmax/iene
  tmax=tmax+tstart
  if(.not.L_first_output_on_restart) then
   tdia=tstart+dtdia
   tout=tstart+dtout
  else
   tdia=tstart
   tout=tstart
  endif

 end select
 contains 

  subroutine check_grid_size
  if(mod(nx,2)/=0)then
   write(6,*)' Wrong x dimension'
   stop
  endif
  if(ny==0)then
   write(6,*)' Wrong y dimension'
   stop
  endif
  if(ny>1)then
   if(mod(ny,2)/=0)then
    write(6,*)' Wrong y dimension'
    stop
   endif
  endif
  end subroutine check_grid_size
!==================================
 end subroutine start

!===================
 end module start_all
 !---------------------------

