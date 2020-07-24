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

module read_json

 use json_module
 use precision_def
 use sim_params_types

 implicit none

 type(json_file)    :: namelist
 type(parameters_t) :: params

 contains

 subroutine read_input_json( input_name, namelist_in )
  character(LEN=:), allocatable, intent(in) :: input_name
  type(json_file), intent(inout) :: namelist_in
  type(json_core) :: json
  type(json_value), pointer :: inp => null()
  type(json_value), pointer :: grid => null()
  type(json_value), pointer :: sim => null()
  type(json_value), pointer :: targ_desc => null()
  logical :: status_ok, found
  character(:), allocatable :: error
  character(13) :: filename_json_out

  ! Initializing the json file
  call namelist_in%initialize()
  call json%initialize()

  ! Loading the file
  call namelist_in%load_file(trim(input_name))
  if ( namelist_in%failed() ) then
   write(6, *) 'Error while reading ', trim(input_name)
   stop
  end if
  ! Check the namelist for errors
  call namelist_in%check_for_errors(status_ok, error_msg = error)
  if ( .not. status_ok ) then
   write(6, *) error
   stop
  end if

  !========================================
  ! Read grid parameters
  !========================================
  
  call namelist_in%get('grid', grid, found)
  if ( .not. found ) then
   write(6, *) "'Grid' not found in json input"
   stop
  end if
  call read_json_grid( grid, params%grid_params )  

  !========================================
  ! Read simulation parameters
  !========================================
  
  call namelist_in%get('simulation', sim, found)
  if ( .not. found ) then
   write(6, *) "'Simulation' not found in json input"
   stop
  end if
  call read_json_sim( sim, params%sim_params )  

  !========================================
  ! Read target parameters
  !========================================
  
  call namelist_in%get('target_description', targ_desc, found)
  if ( .not. found ) then
   write(6, *) "'Target_description' not found in json input"
   stop
  end if

  write(filename_json_out, '(a6,i2.2,a5)') 'input_', 0, '.json'
  call json%create_object(inp, '')
  call json%add(inp, grid)
  call json%add(inp, sim)
  call json%add(inp, targ_desc)
  call json%print(inp, filename_json_out)

  
  call namelist_in%destroy()
  call json%destroy()
  stop

 end subroutine

 subroutine read_json_grid( grid_json, grid_params)
  type(json_value), pointer, intent(in) :: grid_json
  type(grid_parameters_t), intent(inout) :: grid_params
  type(json_core) :: json
  type(grid_parameters_t) :: grid_params_init
  logical :: found

  call json%initialize()

  call json%get(grid_json, 'nx', grid_params%nx, found, &
   default=grid_params_init%nx)
  call json%get(grid_json, 'ny', grid_params%ny, found, &
   default=grid_params_init%ny)
  call json%get(grid_json, 'nz', grid_params%nz, found, &
   default=grid_params_init%nz)
  call json%get(grid_json, 'ny_targ', grid_params%ny_targ, found, &
   default=grid_params_init%ny_targ)
  call json%get(grid_json, 'k0', grid_params%k0, found, &
   default=grid_params_init%k0)
  call json%get(grid_json, 'yx_rat', grid_params%yx_rat, found, &
   default=grid_params_init%yx_rat)
  call json%get(grid_json, 'zx_rat', grid_params%zx_rat, found, &
   default=grid_params_init%zx_rat)

  call json%destroy()
 end subroutine

 subroutine read_json_sim( sim_json, sim_params)
  type(json_value), pointer, intent(in) :: sim_json
  type(simulation_parameters_t), intent(inout) :: sim_params
  type(json_core) :: json
  type(simulation_parameters_t) :: sim_params_init
  logical :: found

  call json%initialize()

  call json%get(sim_json, 'lpf_ord', sim_params%lpf_ord, found, &
   default=sim_params_init%lpf_ord)
  call json%get(sim_json, 'der_ord', sim_params%der_ord, found, &
   default=sim_params_init%der_ord)
  call json%get(sim_json, 'str_flag', sim_params%str_flag, found, &
   default=sim_params_init%str_flag)
  call json%get(sim_json, 'iform', sim_params%iform, found, &
   default=sim_params_init%iform)
  call json%get(sim_json, 'model_id', sim_params%model_id, found, &
   default=sim_params_init%model_id)
  call json%get(sim_json, 'dmodel_id', sim_params%dmodel_id, found, &
   default=sim_params_init%dmodel_id)
  call json%get(sim_json, 'ibx', sim_params%ibx, found, &
   default=sim_params_init%ibx)
  call json%get(sim_json, 'iby', sim_params%iby, found, &
   default=sim_params_init%iby)
  call json%get(sim_json, 'ibz', sim_params%ibz, found, &
   default=sim_params_init%ibz)
  call json%get(sim_json, 'ibeam', sim_params%ibeam, found, &
   default=sim_params_init%ibeam)
  call json%get(sim_json, 'density_limiter', sim_params%density_limiter, found, &
   default=sim_params_init%density_limiter)

  call json%destroy()
 end subroutine

 subroutine read_json_targ( targ_json, targ_params)
  type(json_value), pointer, intent(in) :: targ_json
  type(targ_description_parameters_t), intent(inout) :: targ_params
  type(json_core) :: json
  type(targ_description_parameters_t):: targ_params_init
  logical :: found

  call json%initialize()

  call json%get(targ_json, 'nsp', targ_params%nsp, found, &
   default=targ_params_init%nsp)
  call json%get(targ_json, 'nsb', targ_params%nsb, found, &
   default=targ_params_init%nsb)
  call json%get(targ_json, 'ionz_lev', targ_params%ionz_lev, found, &
   default=targ_params_init%ionz_lev)
  call json%get(targ_json, 'ion_min', targ_params%ion_min, found, &
   default=targ_params_init%ion_min)
  call json%get(targ_json, 'ion_max', targ_params%ion_max, found, &
   default=targ_params_init%ion_max)
  call json%get(targ_json, 'atomic_number', targ_params%atomic_number, found, &
   default=targ_params_init%atomic_number)
  call json%get(targ_json, 'mass_number', targ_params%mass_number, found, &
   default=targ_params_init%mass_number)
  call json%get(targ_json, 't0_pl', targ_params%t0_pl, found, &
   default=targ_params_init%t0_pl)
  call json%get(targ_json, 'np_per_xc', targ_params%np_per_xc, found, &
   default=targ_params_init%np_per_xc)
  call json%get(targ_json, 'np_per_yc', targ_params%np_per_yc, found, &
   default=targ_params_init%np_per_yc)
  call json%get(targ_json, 'np_per_zc', targ_params%np_per_zc, found, &
   default=targ_params_init%np_per_zc)
  call json%get(targ_json, 'concetration', targ_params%concentration, found, &
   default=targ_params_init%concentration)
  call json%get(targ_json, 'lpx', targ_params%lpx, found, &
   default=targ_params_init%lpx)
  call json%get(targ_json, 'lpy', targ_params%lpy, found, &
   default=targ_params_init%lpy)
  call json%get(targ_json, 'n0_ref', targ_params%n0_ref, found, &
   default=targ_params_init%n0_ref)
  call json%get(targ_json, 'np1', targ_params%np1, found, &
   default=targ_params_init%np1)
  call json%get(targ_json, 'np2', targ_params%np2, found, &
   default=targ_params_init%np2)
  call json%get(targ_json, 'r_c', targ_params%r_c, found, &
   default=targ_params_init%r_c)
  call json%get(targ_json, 'l_disable_rng_seed', targ_params%l_disable_rng_seed, found, &
   default=targ_params_init%l_disable_rng_seed)

  call json%destroy()
 end subroutine

end module