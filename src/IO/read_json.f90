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
  character(:), allocatable :: filename_json_out

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

  call json%destroy()
 end subroutine

end module