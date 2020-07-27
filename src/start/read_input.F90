!*****************************************************************************************************!
!                            Copyright 2008-2020  The ALaDyn Collaboration                            !
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

 module read_input

  use read_namelist
#if defined(JSON_FORTRAN)
  use read_json
#endif
  use sim_params_types

  implicit none

 contains


  subroutine read_main_input(parameters_out, exist_json_out)
  !! Reads the input file
   type(parameters_t), intent(inout) :: parameters_out
   logical, intent(inout) :: exist_json_out
   logical :: exist_nml = .false.

   input_json_filename = 'input.json'
   input_namelist_filename = 'input.nml'
   inquire (file=input_namelist_filename, exist=exist_nml)

#if defined(JSON_FORTRAN)
   inquire (file=input_json_filename, exist=exist_json_out)   
   if (exist_json_out) then
    call read_input_json( input_json_filename, namelist, parameters_out)
   end if
#endif

   if (.not. exist_json_out .and. exist_nml) then
    write (6, *) '==========================================='
    write (6, *) 'No JSON input file has been found.'
    write (6, *) 'Instead, the old nml format will be used.'
    write (6, *) '==========================================='
    call read_input_nml
   else if ((.not. exist_json_out) .and. (.not. exist_nml)) then
    write (6, *) 'No input file (.json or .nml) has been found'
    stop
   end if
  end subroutine

#if defined(JSON_FORTRAN)
  subroutine assign_parameters( parameters_in, exist_json_in )
  !! Assignes all the parameters contained in parameters_in and read from the json input file
  !! to the global parameters used in the simulation
    type(parameters_t), intent(in) :: parameters_in
    logical, intent(in) :: exist_json_in

    if (.not. exist_json_in) then
      return
    end if
    !========================
    ! Assign grid parameters
    !========================
    if( .not. allocated(parameters_in%grid_params) ) then
      write(6, *) 'ERROR: grid parameters not found'
      stop
    end if
    nx = parameters_in%grid_params%nx
    ny = parameters_in%grid_params%ny
    nz = parameters_in%grid_params%nz
    ny_targ = parameters_in%grid_params%ny_targ
    k0 = parameters_in%grid_params%k0
    yx_rat = parameters_in%grid_params%yx_rat
    zx_rat = parameters_in%grid_params%zx_rat

    !==============================
    ! Assign simulation parameters
    !==============================
    if( .not. allocated(parameters_in%sim_params) ) then
      write(6, *) 'ERROR: simulation parameters not found'
      stop
    end if
    lpf_ord = parameters_in%sim_params%lpf_ord
    der_ord = parameters_in%sim_params%der_ord
    str_flag = parameters_in%sim_params%str_flag
    iform = parameters_in%sim_params%iform
    model_id = parameters_in%sim_params%model_id
    dmodel_id = parameters_in%sim_params%dmodel_id
    ibx = parameters_in%sim_params%ibx
    iby = parameters_in%sim_params%iby
    ibz = parameters_in%sim_params%ibz
    ibeam = parameters_in%sim_params%ibeam
    density_limiter = parameters_in%sim_params%density_limiter

    !==============================
    ! Assign target parameters
    !==============================
    if( allocated(parameters_in%targ_params) ) then
      nsp = parameters_in%targ_params%nsp
      nsb = parameters_in%targ_params%nsb
      ionz_lev = parameters_in%targ_params%ionz_lev
      ionz_model = parameters_in%targ_params%ionz_model

      if (allocated(parameters_in%targ_params%ion_min)) then
        ion_min(1:nsp - 1) = parameters_in%targ_params%ion_min(1:nsp - 1)
      else
        if (nsp > 1 .and. ionz_lev > 0) then
          write(6, *) "Please specify the ionization level 'ion_min'"
          stop
        end if
        ion_min(:) = 1
      end if

      if (allocated(parameters_in%targ_params%ion_max)) then
        ion_max(1:nsp - 1) = parameters_in%targ_params%ion_max(1:nsp - 1)
      else
        if (nsp > 1 .and. ionz_lev > 0) then
          write(6, *) "Please specify the ionization level 'ion_max'"
          stop
        end if
        ion_max(:) = 1
      end if

      if (allocated(parameters_in%targ_params%atomic_number)) then
        atomic_number(1:nsp - 1) = parameters_in%targ_params%atomic_number(1:nsp - 1)
      else
        if (nsp > 1) then
          write(6, *) "Please specify the species 'atomic_number'"
          stop
        end if
        atomic_number(:) = 1
      end if

      if (allocated(parameters_in%targ_params%mass_number)) then
        mass_number(1:nsp - 1) = parameters_in%targ_params%mass_number(1:nsp - 1)
      else
        if (nsp > 1) then
          write(6, *) "Please specify the species 'mass_number'"
          stop
        end if
        mass_number(:) = 1
      end if

      if (allocated(parameters_in%targ_params%t0_pl)) then
        t0_pl(1:nsp) = parameters_in%targ_params%t0_pl(1:nsp)
      else
        t0_pl(1:nsp) = zero_dp
      end if

      if (allocated(parameters_in%targ_params%np_per_xc)) then
        np_per_xc(1:nsp) = parameters_in%targ_params%np_per_xc(1:nsp)
      else
        np_per_xc(1:nsp) = -1
      end if

      if (allocated(parameters_in%targ_params%np_per_yc)) then
        np_per_yc(1:nsp) = parameters_in%targ_params%np_per_yc(1:nsp)
      else
        np_per_yc(1:nsp) = -1
      end if

      if (allocated(parameters_in%targ_params%np_per_zc)) then
        np_per_zc(1:nsp) = parameters_in%targ_params%np_per_zc(1:nsp)
      else
        np_per_zc(1:nsp) = -1
      end if

      if (allocated(parameters_in%targ_params%concentration)) then
        concentration(1:nsp) = parameters_in%targ_params%concentration(1:nsp)
      else
        concentration(1:nsp) = zero_dp
        concentration(1) = one_dp
      end if

      if (allocated(parameters_in%targ_params%lpx)) then
        if (size(parameters_in%targ_params%lpx) < 7) then
          write(6, *) "'lpx' size must be of 7 elements"
          stop
        end if
        lpx(1:7) = parameters_in%targ_params%lpx(1:7)
      else
        lpx(1:7) = zero_dp
      end if

      if (allocated(parameters_in%targ_params%lpy)) then
        if (size(parameters_in%targ_params%lpy) < 2) then
          write(6, *) "'lpy' size must be of 2 elements"
          stop
        end if
        lpy(1:2) = parameters_in%targ_params%lpy(1:2)
      else
        lpy(1:2) = zero_dp
      end if

      n0_ref = parameters_in%targ_params%n0_ref
      np1 = parameters_in%targ_params%np1
      np2 = parameters_in%targ_params%np2
      r_c = parameters_in%targ_params%r_c
      l_disable_rng_seed = parameters_in%targ_params%l_disable_rng_seed
    else
      nsp = 0
      np_per_xc(:) = 0
      np_per_yc(:) = 0
      np_per_zc(:) = 0
      lpx(:) = zero_dp
      lpy(:) = zero_dp
    end if

    !==============================
    ! Assign beam parameters
    !==============================
    if( allocated(parameters_in%beam_params) ) then
      nb_1 = parameters_in%beam_params%nb_1
      xc_1 = parameters_in%beam_params%xc_1
      gam_1 = parameters_in%beam_params%gam_1
      sxb_1 = parameters_in%beam_params%sxb_1
      syb_1 = parameters_in%beam_params%syb_1
      epsy_1 = parameters_in%beam_params%epsy_1
      epsz_1 = parameters_in%beam_params%epsz_1
      dg_1 = parameters_in%beam_params%dg_1
      charge_1 = parameters_in%beam_params%charge_1
      ap1_twiss = parameters_in%beam_params%ap1_twiss
      bt1_twiss = parameters_in%beam_params%bt1_twiss
      t_inject = parameters_in%beam_params%t_inject
    else
      nsb = 0
    end if

    !==============================
    ! Assign window parameters
    !==============================
    if( allocated(parameters_in%window_params) ) then
      w_sh = parameters_in%window_params%w_sh
      wi_time = parameters_in%window_params%wi_time
      wf_time = parameters_in%window_params%wf_time
      w_speed = parameters_in%window_params%w_speed
    else
      w_sh = 0
      wi_time = zero_dp
      wf_time = zero_dp
      w_speed = zero_dp
    end if

    !==============================
    ! Assign output parameters
    !==============================
    if( .not. allocated(parameters_in%output_params) ) then
      write(6, *) 'ERROR: output parameters not found'
      stop
    end if
    nouts = parameters_in%output_params%nouts
    iene = parameters_in%output_params%iene
    nvout = parameters_in%output_params%nvout
    nden = parameters_in%output_params%nden
    ncurr = parameters_in%output_params%ncurr
    npout = parameters_in%output_params%npout
    jump = parameters_in%output_params%jump
    pjump = parameters_in%output_params%pjump
    gam_min = parameters_in%output_params%gamma_min
    xp0_out = parameters_in%output_params%xp0_out
    xp1_out = parameters_in%output_params%xp1_out
    yp_out = parameters_in%output_params%yp_out
    tmax = parameters_in%output_params%tmax
    cfl = parameters_in%output_params%cfl
    new_sim = parameters_in%output_params%new_sim
    id_new = parameters_in%output_params%id_new
    dump = parameters_in%output_params%dump
    l_force_singlefile_output = parameters_in%output_params%l_force_single_output
    l_print_j_on_grid = parameters_in%output_params%l_print_j_on_grid
    l_first_output_on_restart = parameters_in%output_params%l_first_output_on_restart
    l_env_modulus = parameters_in%output_params%l_env_modulus
    time_interval_dumps = parameters_in%output_params%time_interval_dumps

    !==============================
    ! Assign MPI parameters
    !==============================
    if( allocated(parameters_in%mpi_params) ) then
      nprocx = parameters_in%mpi_params%nprocx
      nprocy = parameters_in%mpi_params%nprocy
      nprocz = parameters_in%mpi_params%nprocz
    else
      nprocx = 1
      nprocy = 1
      nprocz = 1
    end if

    
    !==============================
    ! Assign tracking parameters
    !==============================
    if( allocated(parameters_in%track_params) ) then
      p_tracking = parameters_in%track_params%p_tracking
      tkjump = parameters_in%track_params%tkjump
      nkjump = parameters_in%track_params%nkjump
      txmin = parameters_in%track_params%txmin
      txmax = parameters_in%track_params%txmax
      tymin = parameters_in%track_params%tymin
      tymax = parameters_in%track_params%tymax
      tzmin = parameters_in%track_params%tzmin
      tzmax = parameters_in%track_params%tzmax
      t_in = parameters_in%track_params%t_in
      t_out = parameters_in%track_params%t_out
    else
      p_tracking = .false.
    end if

  end subroutine
#else
  subroutine assign_parameters( parameters_in, exist_json_in )
    type(parameters_t), intent(in) :: parameters_in
    logical, intent(in) :: exist_json_in
    type(parameters_t) :: parameters_mock

    if (.not. exist_json_in) then
      return
    end if

    parameters_mock = parameters_in
  end subroutine
#endif
 end module
