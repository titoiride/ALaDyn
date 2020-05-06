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

 module grid_tracking_part_connect

  use pstruct_data
  use fstruct_data
  use grid_part_lib
  use particles_def
  use mpi_var

  type(interp_coeff), private, allocatable, save :: interp
  real(dp), dimension(:, :), allocatable, private :: gtpc_xx 

  contains

  subroutine a_on_tracking_particles( field, spec_in, np )
   real (dp), intent (in) :: field(:, :, :)
   type (species_new), intent (inout) :: spec_in
   integer, intent(in) :: np
   real(dp), allocatable, dimension(:) :: ap
   real(dp) :: dvol
   integer :: npt, i1, i2, j1, j2, k1, k2, n, cc
   logical, dimension(:), allocatable :: track_mask
   real(dp), dimension(:), allocatable :: interpolated_field

   npt = spec_in%pick_tot_tracked_parts()

   if (npt == 0) return
   allocate( ap(npt), source=zero_dp )
   allocate( interpolated_field(np), source=zero_dp )
   call interp_realloc(interp, np, spec_in%pick_dimensions())
   call array_realloc_1d( track_mask, np )
   associate( inds => spec_in%call_component(INDEX_COMP, lb=1, ub=np) )

    track_mask(1:np) = (int(inds) > 0)
 
   end associate
   
   cc = COUNT(track_mask(1:np) )
   if (npt /= cc ) then
    call write_warning(text='Error in counting tracked particles', task=mype)
   end if
 
   select case (ndim)
   case(2)
    k2 = 1
    call xx_realloc(gtpc_xx, np, 2)
    gtpc_xx(1:np, 1) = set_local_positions( spec_in, X_COMP, track_mask )
    gtpc_xx(1:np, 2) = set_local_positions( spec_in, Y_COMP, track_mask )

    call qden_2d_wgh( gtpc_xx(1:np, 1:2), interp )

    associate( ax1 => interp%coeff_x_rank2, &
               ay1 => interp%coeff_y_rank2, &
               i => interp%ix_rank2, &
               j => interp%iy_rank2 )

     do n = 1, npt
      do j1 = 0, 2
       j2 = j(n) + j1
       dvol = ay1( n, j2 )
       do i1 = 0, 2
        i2 = i(n) + i1
        ap(n) = ap(n) + ax1(n, i2)*dvol*field(i2, j2, k2)
       end do
      end do
     end do

    end associate

    interpolated_field = UNPACK( ap(1:npt), track_mask(1:np), interpolated_field(1:np) )

    call spec_in%set_component( interpolated_field, A_PARTICLE, lb=1, ub=np)

   case(3)
    call xx_realloc(gtpc_xx, np, 3)
    gtpc_xx(1:np, 1) = set_local_positions( spec_in, X_COMP, track_mask )
    gtpc_xx(1:np, 2) = set_local_positions( spec_in, Y_COMP, track_mask )
    gtpc_xx(1:np, 3) = set_local_positions( spec_in, Z_COMP, track_mask )

    call qden_3d_wgh( gtpc_xx(1:np, 1:3), interp )

    associate( ax1 => interp%coeff_x_rank2, &
               ay1 => interp%coeff_y_rank2, &
               az1 => interp%coeff_z_rank2, &
               i => interp%ix_rank2, &
               j => interp%iy_rank2, &
               k => interp%iz_rank2 )

               do n = 1, npt
                do k1 = 0, 2
                 k2 = k(n) + k1
                 dvol = az1( n, k2 )
                 do j1 = 0, 2
                  j2 = j(n) + j1
                  dvol = dvol*ay1( n, j2 )
                  do i1 = 0, 2
                   i2 = i(n) + i1
                   ap(n) = ap(n) + ax1(n, i2)*dvol*field(i2, j2, k2)
                  end do
                 end do
                end do
               end do

    end associate

    interpolated_field = UNPACK( ap(1:npt), track_mask(1:np), interpolated_field(1:np) )

    call spec_in%set_component( interpolated_field, A_PARTICLE, lb=1, ub=np)

   end select
  end subroutine

 end module