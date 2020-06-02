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

  implicit none
  type(interp_coeff), private, allocatable, save :: interp
  real(dp), dimension(:, :), allocatable, private :: gtpc_xx
  integer, parameter :: Y_POLARIZATION = 1
  integer, parameter :: Z_POLARIZATION = 2

  contains

  subroutine a_on_tracking_particles( field, spec_in, spec_aux_in, np, order, polarization, mask_in )
   !! Subroutine that interpolates vector potential A on the tracking particles.

   real (dp), intent (in) :: field(:, :, :)
   type (species_new), intent (inout) :: spec_in
   type (species_aux), intent (in) :: spec_aux_in
   integer, intent(in) :: np, order, polarization
   logical, dimension(:), intent(in), optional :: mask_in
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
   if ( PRESENT( mask_in ) ) then
    track_mask(1:np) = mask_in(1:np)
   else
    associate( inds => spec_in%call_component(INDEX_COMP, lb=1, ub=np) )

     track_mask(1:np) = (int(inds) > 0)

    end associate
   end if

   cc = COUNT(track_mask(1:np) )
   if (npt /= cc ) then
    call write_warning(text='Error in counting tracked particles', task=mype)
   end if
   
   select case (ndim)
   case(2)

    k2 = 1
    call xx_realloc(gtpc_xx, npt, 2)
    gtpc_xx(1:npt, 1) = set_local_positions( spec_aux_in, X_COMP, mask_in=track_mask(1:np) )
    gtpc_xx(1:npt, 2) = set_local_positions( spec_aux_in, Y_COMP, mask_in=track_mask(1:np) )

    select case (order)

    case(0)
     ! Envelope case
     call qden_2d_wgh( gtpc_xx(1:npt, 1:2), interp )

     associate( ax1 => interp%coeff_x_rank2, &
                ay1 => interp%coeff_y_rank2, &
                i => interp%ix_rank2, &
                j => interp%iy_rank2 )

                do n = 1, npt
                 do j1 = 0, 2
                  j2 = j(n) + j1
                  dvol = ay1( n, j1 )
                  do i1 = 0, 2
                   i2 = i(n) + i1
                   ap(n) = ap(n) + ax1(n, i1)*dvol*field(i2, j2, k2)
                  end do
                 end do
                end do

     end associate

    case(1:2)
     ! Full Pic case

     call qqh_2d_spline( gtpc_xx(1:npt, 1:2), interp )

     select case(polarization)
     !Integration order used for computing a from Ey
     case(Y_POLARIZATION)


      associate( ax1 => interp%coeff_x_rank2, &
                 ayh1 => interp%h_coeff_y_rank2, &
                 i => interp%ix_rank2, &
                 jh => interp%ihy_rank2 )

       do n = 1, npt
        do j1 = 0, 2
         j2 = jh(n) + j1
         dvol = ayh1( n, j1 )
         do i1 = 0, 2
          i2 = i(n) + i1
          ap(n) = ap(n) + ax1(n, i1)*dvol*field(i2, j2, k2)
         end do
        end do
       end do

      end associate

     case(Z_POLARIZATION)

      associate( ax1 => interp%coeff_x_rank2, &
                 ay1 => interp%coeff_y_rank2, &
                 i => interp%ix_rank2, &
                 j => interp%iy_rank2 )

       do n = 1, npt
        do j1 = 0, 2
         j2 = j(n) + j1
         dvol = ay1( n, j1 )
         do i1 = 0, 2
          i2 = i(n) + i1
          ap(n) = ap(n) + ax1(n, i1)*dvol*field(i2, j2, k2)
         end do
        end do
       end do

      end associate

     end select
    end select

   case(3)

    call xx_realloc(gtpc_xx, npt, 3)
    gtpc_xx(1:npt, 1) = set_local_positions( spec_aux_in, X_COMP, mask_in=track_mask(1:np) )
    gtpc_xx(1:npt, 2) = set_local_positions( spec_aux_in, Y_COMP, mask_in=track_mask(1:np) )
    gtpc_xx(1:npt, 3) = set_local_positions( spec_aux_in, Z_COMP, mask_in=track_mask(1:np) )

    select case(order)
    case(0)
    ! Envelope case
     call qden_3d_wgh( gtpc_xx(1:npt, 1:3), interp )

     associate( ax1 => interp%coeff_x_rank2, &
                ay1 => interp%coeff_y_rank2, &
                az1 => interp%coeff_z_rank2, &
                i => interp%ix_rank2, &
                j => interp%iy_rank2, &
                k => interp%iz_rank2 )

                do n = 1, npt
                 do k1 = 0, 2
                  k2 = k(n) + k1
                  dvol = az1( n, k1 )
                  do j1 = 0, 2
                   j2 = j(n) + j1
                   dvol = dvol*ay1( n, j1 )
                   do i1 = 0, 2
                    i2 = i(n) + i1
                    ap(n) = ap(n) + ax1(n, i1)*dvol*field(i2, j2, k2)
                   end do
                  end do
                 end do
                end do

     end associate

    case(1)
     !Full Pic
     call qqh_3d_spline( gtpc_xx(1:npt, 1:3), interp )
     select case(polarization)
     !Integration order used for computing a from Ey
     case(Y_POLARIZATION)


      associate( ax1 => interp%coeff_x_rank2, &
                 ayh1 => interp%h_coeff_y_rank2, &
                 az1 => interp%coeff_z_rank2, &
                 i => interp%ix_rank2, &
                 jh => interp%ihy_rank2, &
                 k => interp%iz_rank2 )

                 do n = 1, npt
                  do k1 = 0, 2
                   k2 = k(n) + k1
                   dvol = az1( n, k1 )
                   do j1 = 0, 2
                    j2 = jh(n) + j1
                    dvol = dvol*ayh1( n, j1 )
                    do i1 = 0, 2
                     i2 = i(n) + i1
                     ap(n) = ap(n) + ax1(n, i1)*dvol*field(i2, j2, k2)
                    end do
                   end do
                  end do
                 end do

      end associate

     case(Z_POLARIZATION)

      associate( ax1 => interp%h_coeff_x_rank2, &
                 ay1 => interp%coeff_y_rank2, &
                 azh1 => interp%h_coeff_z_rank2, &
                 i => interp%ihx_rank2, &
                 j => interp%iy_rank2, &
                 kh => interp%ihz_rank2 )

                 do n = 1, npt
                  do k1 = 0, 2
                   k2 = kh(n) + k1
                   dvol = azh1( n, k1 )
                   do j1 = 0, 2
                    j2 = j(n) + j1
                    dvol = dvol*ay1( n, j1 )
                    do i1 = 0, 2
                     i2 = i(n) + i1
                     ap(n) = ap(n) + ax1(n, i1)*dvol*field(i2, j2, k2)
                    end do
                   end do
                  end do
                 end do

      end associate
     end select
    end select
   end select
   
   ! Field assignment on tracked particles
   interpolated_field = UNPACK( ap(1:npt), track_mask(1:np), interpolated_field(1:np) )
   call spec_in%set_component( interpolated_field, A_PARTICLE, lb=1, ub=np)

  end subroutine

 end module