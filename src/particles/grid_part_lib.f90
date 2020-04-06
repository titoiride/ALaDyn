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

 module grid_part_lib

  use common_param
  use grid_param
  use stretched_grid
  use particles_def
  use interpolation_lib
  use array_alloc, only: array_realloc_1d

  implicit none

  real (dp), parameter :: half = 0.5
  real (sp), parameter :: shx = 3., shy = 3., shz = 3.
  integer (kind=2) :: err_ind
  real(dp), allocatable, dimension(:), private, save :: gpl_xx
  real(dp), allocatable, dimension(:), private, save :: gpl_sx
  integer, parameter, private :: max_order = 2
  integer, parameter, private :: max_h_order = 2

  interface qqh_1d_spline
   module procedure :: qqh_1d_spline_real
   module procedure :: qqh_1d_spline_vector
  end interface
  
  interface qden_1d_wgh
   module procedure :: qden_1d_wgh_real
   module procedure :: qden_1d_wgh_vector
  end interface
  
  interface qlh_2d_spline
   module procedure :: qlh_2d_spline_real
   module procedure :: qlh_2d_spline_vector
  end interface
    
  interface qqh_2d_spline
   module procedure :: qqh_2d_spline_real
   module procedure :: qqh_2d_spline_vector
  end interface
    
    
  interface qden_2d_wgh
   module procedure :: qden_2d_wgh_real
   module procedure :: qden_2d_wgh_vector
  end interface
  
  interface cden_2d_wgh
   module procedure :: cden_2d_wgh_real
   module procedure :: cden_2d_wgh_vector
  end interface

  interface qlh_3d_spline
   module procedure :: qlh_3d_spline_real
   module procedure :: qlh_3d_spline_vector
  end interface

  interface qqh_3d_spline
   module procedure :: qqh_3d_spline_real
   module procedure :: qqh_3d_spline_vector
  end interface

  interface qden_3d_wgh
   module procedure :: qden_3d_wgh_real
   module procedure :: qden_3d_wgh_vector
  end interface

  interface cden_3d_wgh
   module procedure :: cden_3d_wgh_real
   module procedure :: cden_3d_wgh_vector
  end interface

  interface set_local_positions
   module procedure set_local_positions_sa
   module procedure set_local_positions_sn
  end interface

  contains

  !===========================================
  !     Computes particle density charge on grid integer points
  !==============================================
  !DIR$ ATTRIBUTES INLINE :: ql_interpolate,qq_interpolate
  !====================

  
  subroutine interp_realloc(interp_in, npart_new, dimensions)
   type(interp_coeff), allocatable, intent(inout) :: interp_in
   integer, intent(in) :: npart_new, dimensions

   if ( allocated(interp_in) ) then   
    if (interp_in%n_parts < npart_new ) then
     call interp_in%sweep()
     call interp_in%new_interp(npart_new, max_order, &
      max_h_order, dimensions)
    end if
   else
    allocate( interp_in )
    call interp_in%new_interp(npart_new, max_order, &
      max_h_order, dimensions)
   end if

  end subroutine
  !==========================================
  subroutine xx_realloc( xx_in, npart_new, dimensions )
   real(dp), allocatable, dimension(:, :), intent(inout) :: xx_in
   integer, intent(in) :: npart_new, dimensions
   integer :: allocstatus, deallocstatus

   if ( allocated(xx_in) ) then   
    if (SIZE(xx_in, DIM=1) < npart_new .or. &
     SIZE(xx_in, DIM=2) < dimensions ) then
     deallocate(xx_in, stat=deallocstatus)
     allocate(xx_in(npart_new, dimensions), stat=allocstatus)
    end if
   else
    allocate(xx_in(npart_new, dimensions), stat=allocstatus)
   end if
  end subroutine
  !==========================================

   !=======================
   !  1D
   !=======================
  subroutine qqh_1d_spline_real( xp, interp_in )
    !! Quadratic interpolation at both integer and half-integer
    !! cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   real (dp) :: xx, sx
   !======================

   xx = shx + xp(1)
   interp_in%ix = int(xx+half)
   sx = xx - real(interp_in%ix, dp)
   call second_order( sx, interp_in%coeff_x )

   interp_in%ihx = int(xx)
   sx = xx - real(interp_in%ihx, dp)
   call second_order( sx, interp_in%h_coeff_x )

   interp_in%ix = interp_in%ix - 1
   interp_in%ihx = interp_in%ihx - 1

  end subroutine
  !=======================

  subroutine qqh_1d_spline_vector( xp, interp_in )
    !! Quadratic interpolation at both integer and half-integer
    !! cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   call array_realloc_1d(gpl_xx, length)
   call array_realloc_1d(gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int( gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   interp_in%ihx_rank2(1:length) = int(gpl_xx(1:length))
   gpl_sx(1:length) = gpl_xx(1:length) - &
   real(interp_in%ihx_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%h_coeff_x_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%ihx_rank2(1:length) = interp_in%ihx_rank2(1:length) - 1

  end subroutine
  !=======================

  subroutine qden_1d_wgh_real( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix
   real (dp) :: xx, sx
   !======================

   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x)

   interp_in%ix = ix - 1
  end subroutine

  subroutine qden_1d_wgh_vector( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )

   call array_realloc_1d( gpl_xx, length)
   call array_realloc_1d( gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int( gpl_xx(1:length) + half )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
  end subroutine

  !=======================
  !  2D
  !=======================

  subroutine qlh_2d_spline_real( xp, interp_in )
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, iy
   real (dp) :: xx, sx
   !======================

   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x )
   call first_order( sx + half, interp_in%h_coeff_x )

   xx = shy + xp(2)
   iy = int(xx+half)
   sx = xx - real(iy, dp)

   call second_order( sx, interp_in%coeff_y )
   call first_order( sx + half, interp_in%h_coeff_y )


   interp_in%ix = ix - 1
   interp_in%ihx = ix
   interp_in%iy = iy - 1
   interp_in%ihy = iy
   
  end subroutine

  subroutine qlh_2d_spline_vector( xp, interp_in )
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   call array_realloc_1d(gpl_xx, length)
   call array_realloc_1d(gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )
   call first_order( gpl_sx(1:length) + half, interp_in%h_coeff_x_rank2 )

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )
   call first_order( gpl_sx(1:length) + half, interp_in%h_coeff_y_rank2 )

   interp_in%ihx_rank2(1:length) = interp_in%ix_rank2(1:length)
   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%ihy_rank2(1:length) = interp_in%iy_rank2(1:length)
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1
   
  end subroutine
  !====================

  subroutine qqh_2d_spline_real( xp, interp_in )
   !! Quadratic interpolation at both integer and half-integer
   !! cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, ihx, iy, ihy
   real (dp) :: xx, sx
   !======================

   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x )

   ihx = int(xx)
   sx = xx - real(ihx, dp) - half

   call second_order( sx, interp_in%h_coeff_x )

   xx = shy + xp(2)
   iy = int(xx+half)
   sx = xx - real(iy, dp)

   call second_order( sx, interp_in%coeff_y )

   ihy = int(xx)
   sx = xx - real(ihy, dp) - half

   call second_order( sx, interp_in%h_coeff_y )

   interp_in%ix = ix - 1
   interp_in%ihx = ihx - 1
   interp_in%iy = iy - 1
   interp_in%ihy = ihy - 1
  end subroutine
  !=======================================

  subroutine qqh_2d_spline_vector( xp, interp_in )
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   call array_realloc_1d(gpl_xx, length)
   call array_realloc_1d(gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   interp_in%ihx_rank2(1:length) = int(gpl_xx(1:length))
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ihx_rank2(1:length), dp) - half

   call second_order( gpl_sx(1:length), interp_in%h_coeff_x_rank2 )

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )

   interp_in%ihy_rank2(1:length) = int(gpl_xx(1:length))
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ihy_rank2(1:length), dp) - half

   call second_order( gpl_sx(1:length), interp_in%h_coeff_y_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%ihx_rank2(1:length) = interp_in%ihx_rank2(1:length) - 1
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1
   interp_in%ihy_rank2(1:length) = interp_in%ihy_rank2(1:length) - 1
   
  end subroutine
  !====================

  subroutine qden_2d_wgh_real( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, iy
   real (dp) :: xx, sx
   !======================
   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x )

   xx = shy + xp(2)
   iy = int(xx+half)
   sx = xx - real(iy, dp)

   call second_order( sx, interp_in%coeff_y )

   interp_in%ix = ix - 1
   interp_in%iy = iy - 1
  end subroutine
  !======================

  subroutine qden_2d_wgh_vector( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )
   call array_realloc_1d( gpl_xx, length )
   call array_realloc_1d( gpl_sx, length )

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int( gpl_xx(1:length) + half )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int( gpl_xx(1:length) + half )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1

  end subroutine

  subroutine cden_2d_wgh_real( xp, interp_in )
   !! Cubic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, iy
   real (dp) :: xx, sx
   !======================
   xx = shy + xp(1)
   ix = int(xx)
   sx = xx - real(ix, dp)

   call third_order( sx, interp_in%coeff_x )
   
   xx = shy + xp(2)
   iy = int(xx)
   sx = xx - real(iy, dp)
   
   call third_order( sx, interp_in%coeff_y )

   interp_in%ix = ix - 1
   interp_in%iy = iy - 1
  end subroutine
  !==========================

  subroutine cden_2d_wgh_vector( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )
   call array_realloc_1d( gpl_xx, length )
   call array_realloc_1d( gpl_sx, length )

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int( gpl_xx(1:length) )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call third_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int( gpl_xx(1:length) )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call third_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1

  end subroutine
  !=======================
  !  3D
  !=======================
  subroutine qlh_3d_spline_real( xp, interp_in )
   !! Quadratic interpolation at integer cells, linear at
   !! half-integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, iy, iz
   real (dp) :: xx, sx
   !======================

   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x )
   call first_order( sx + half, interp_in%h_coeff_x )

   xx = shy + xp(2)
   iy = int(xx+half)
   sx = xx - real(iy, dp)

   call second_order( sx, interp_in%coeff_y )
   call first_order( sx + half, interp_in%h_coeff_y )

   xx = shz + xp(3)
   iz = int(xx+half)
   sx = xx - real(iz, dp)

   call second_order( sx, interp_in%coeff_z )
   call first_order( sx + half, interp_in%h_coeff_z )

   interp_in%ix = ix - 1
   interp_in%ihx = ix
   interp_in%iy = iy - 1
   interp_in%ihy = iy
   interp_in%iz = iz - 1
   interp_in%ihz = iz
  end subroutine
  !=================================

  subroutine qlh_3d_spline_vector( xp, interp_in )
   !! Quadratic interpolation at integer cells, linear at
   !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )
   call array_realloc_1d( gpl_xx, length)
   call array_realloc_1d( gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )
   call first_order( gpl_sx(1:length), interp_in%h_coeff_x_rank2 )

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )
   call first_order( gpl_sx(1:length), interp_in%h_coeff_y_rank2 )

   gpl_xx(1:length) = shz + xp(1:length, 3)
   interp_in%iz_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iz_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_z_rank2 )
   call first_order( gpl_sx(1:length), interp_in%h_coeff_z_rank2 )

   interp_in%ihx_rank2(1:length) = interp_in%ix_rank2(1:length)
   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%ihy_rank2(1:length) = interp_in%iy_rank2(1:length)
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1
   interp_in%ihz_rank2(1:length) = interp_in%iz_rank2(1:length)
   interp_in%iz_rank2(1:length) = interp_in%iz_rank2(1:length) - 1
  end subroutine
  !=================================
  subroutine qqh_3d_spline_real( xp, interp_in )
   !! Quadratic interpolation at both integer and half-integer
   !! cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, ihx, iy, ihy, iz, ihz
   real (dp) :: xx, sx
   !=====================

   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x )

   ihx = int(xx)
   sx = xx - real(ihx, dp) - half

   call second_order( sx, interp_in%h_coeff_x )

   xx = shy + xp(2)
   iy = int(xx+half)
   sx = xx - real(iy, dp)

   call second_order( sx, interp_in%coeff_y )

   ihy = int(xx)
   sx = xx - real(ihy, dp) - half

   call second_order( sx, interp_in%h_coeff_y )

   xx = shz + xp(3)
   iz = int(xx+half)
   sx = xx - real(iz, dp)

   call second_order( sx, interp_in%coeff_z )

   ihz = int(xx)
   sx = xx - real(ihz, dp) - half

   call second_order( sx, interp_in%h_coeff_z )

   interp_in%ix = ix - 1
   interp_in%ihx = ihx - 1
   interp_in%iy = iy - 1
   interp_in%ihy = ihy - 1
   interp_in%iz = iz - 1
   interp_in%ihz = ihz - 1
  end subroutine

  subroutine qqh_3d_spline_vector( xp, interp_in )
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   call array_realloc_1d( gpl_xx, length)
   call array_realloc_1d( gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   interp_in%ihx_rank2(1:length) = int(gpl_xx(1:length))
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ihx_rank2(1:length), dp) - half

   call second_order( gpl_sx(1:length), interp_in%h_coeff_x_rank2 ) 

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )

   interp_in%ihy_rank2(1:length) = int(gpl_xx(1:length))
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ihy_rank2(1:length), dp) - half

   call second_order( gpl_sx(1:length), interp_in%h_coeff_y_rank2 )

   gpl_xx(1:length) = shz + xp(1:length, 3)
   interp_in%iz_rank2(1:length) = int(gpl_xx(1:length) + half)
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iz_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_z_rank2 )

   interp_in%ihz_rank2(1:length) = int(gpl_xx(1:length))
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ihz_rank2(1:length), dp) - half

    call second_order( gpl_sx(1:length), interp_in%h_coeff_z_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%ihx_rank2(1:length) = interp_in%ihx_rank2(1:length) - 1
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1
   interp_in%ihy_rank2(1:length) = interp_in%ihy_rank2(1:length) - 1
   interp_in%iz_rank2(1:length) = interp_in%iz_rank2(1:length) - 1
   interp_in%ihz_rank2(1:length) = interp_in%ihz_rank2(1:length) - 1
   
  end subroutine
  !====================

  subroutine qden_3d_wgh_real( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, iy, iz
   real (dp) :: xx, sx
   !======================
   xx = shx + xp(1)
   ix = int(xx+half)
   sx = xx - real(ix, dp)

   call second_order( sx, interp_in%coeff_x )

   xx = shy + xp(2)
   iy = int(xx+half)
   sx = xx - real(iy, dp)

   call second_order( sx, interp_in%coeff_y )

   xx = shz + xp(3)
   iz = int(xx+half)
   sx = xx - real(iz, dp)

   call second_order( sx, interp_in%coeff_z )

   interp_in%ix = ix - 1
   interp_in%iy = iy - 1
   interp_in%iz = iz - 1
  end subroutine
  !================================

  subroutine qden_3d_wgh_vector( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )
   call array_realloc_1d(gpl_xx, length)
   call array_realloc_1d(gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int( gpl_xx(1:length) + half )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_x_rank2 )

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int( gpl_xx(1:length) + half )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call second_order( gpl_sx(1:length), interp_in%coeff_y_rank2 )

   gpl_xx(1:length) = shz + xp(1:length, 3)
   interp_in%iz_rank2(1:length) = int( gpl_xx(1:length) + half )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iz_rank2(1:length), dp)

    call second_order( gpl_sx(1:length), interp_in%coeff_z_rank2 )

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1
   interp_in%iz_rank2(1:length) = interp_in%iz_rank2(1:length) - 1

  end subroutine

  subroutine cden_3d_wgh_real( xp, interp_in )
   !! Cubic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff), intent(inout) :: interp_in
   integer :: ix, iy, iz
   real (dp) :: xx, sx
   !======================
   !cubic spline to integer index
   xx = shx + xp(1)
   ix = int(xx)
   sx = xx - real(ix, dp)

   call third_order( sx, interp_in%coeff_x )

   xx = shy + xp(2)
   iy = int(xx)
   sx = xx - real(iy, dp)

   call third_order( sx, interp_in%coeff_y )

   xx = shz + xp(3)
   iz = int(xx)
   sx = xx - real(iz, dp)

   call third_order( sx, interp_in%coeff_z )

   interp_in%ix = ix - 1
   interp_in%iy = iy - 1
   interp_in%iz = iz - 1

  end subroutine

  subroutine cden_3d_wgh_vector( xp, interp_in )
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff), intent(inout) :: interp_in
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )
   call array_realloc_1d(gpl_xx, length)
   call array_realloc_1d(gpl_sx, length)

   gpl_xx(1:length) = shx + xp(1:length, 1)
   interp_in%ix_rank2(1:length) = int( gpl_xx(1:length) )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%ix_rank2(1:length), dp)

   call third_order( gpl_sx(1:length), interp_in%coeff_x_rank2)

   gpl_xx(1:length) = shy + xp(1:length, 2)
   interp_in%iy_rank2(1:length) = int( gpl_xx(1:length) )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iy_rank2(1:length), dp)

   call third_order( gpl_sx(1:length), interp_in%coeff_y_rank2)

   gpl_xx(1:length) = shz + xp(1:length, 3)
   interp_in%iz_rank2(1:length) = int( gpl_xx(1:length) )
   gpl_sx(1:length) = gpl_xx(1:length) - &
    real(interp_in%iz_rank2(1:length), dp)

   call third_order( gpl_sx(1:length), interp_in%coeff_y_rank2)

   interp_in%ix_rank2(1:length) = interp_in%ix_rank2(1:length) - 1
   interp_in%iy_rank2(1:length) = interp_in%iy_rank2(1:length) - 1
   interp_in%iz_rank2(1:length) = interp_in%iz_rank2(1:length) - 1

  end subroutine
  !====================================

  !DIR$ ATTRIBUTES INLINE :: ql_interpolate
  function set_local_positions_sn( pt_array, component ) result(position)
   type (species_new), intent(in) :: pt_array
   integer, intent(in) :: component
   integer :: np
   real(dp), allocatable, dimension(:) :: position

   np = pt_array%how_many()
   allocate( position(np) )
   position(1:np) = pt_array%call_component( component, lb=1, ub=np )

   select case(component)
   case(X_COMP)
    call map2dx_part_sind( np, position )
   case(Y_COMP)
    call map2dy_part_sind( np, position )
   case(Z_COMP)
    call map2dz_part_sind( np, position )
   end select

  end function

  function set_local_positions_sa( pt_array, component ) result(position)
   type (species_aux), intent(in) :: pt_array
   integer, intent(in) :: component
   integer :: np
   real(dp), allocatable, dimension(:) :: position

   np = pt_array%how_many()
   allocate( position(np) )
   position(1:np) = pt_array%call_component( component, lb=1, ub=np )

   select case(component)
   case(X_COMP)
    call map2dx_part_sind( np, position )
   case(Y_COMP)
    call map2dy_part_sind( np, position )
   case(Z_COMP)
    call map2dz_part_sind( np, position )
   end select

  end function

  subroutine set_local_2d_positions(pt_loc, n1, np)

   real (dp), intent (inout) :: pt_loc(:, :)
   integer, intent (in) :: n1, np
   integer :: n
   !=========================
   do n = 1, np
    pt_loc(n, 1) = dx_inv*(pt_loc(n,1)-xmn)
   end do
   if (n1==0) return
   if (n_str==0) then
    do n = 1, np
     pt_loc(n, 2) = dy_inv*(pt_loc(n,2)-ymn)
    end do
   else
    call map2dy_part_sind(np, 2, pt_loc)
   end if
   if (n1==1) return

   do n = 1, np
    pt_loc(n, 3) = dx_inv*(pt_loc(n,3)-xmn)
   end do
   if (n_str==0) then
    do n = 1, np
     pt_loc(n, 4) = dy_inv*(pt_loc(n,4)-ymn)
    end do
   else
    call map2dy_part_sind(np, 4, pt_loc)
   end if

  end subroutine
  !======================
  subroutine set_local_3d_positions(pt_loc, n1, np)
   real (dp), intent (inout) :: pt_loc(:, :)
   integer, intent (in) :: n1, np
   integer :: n
   !=========================
   do n = 1, np
    pt_loc(n, 1) = dx_inv*(pt_loc(n,1)-xmn)
   end do
   if (n_str==0) then
    do n = 1, np
     pt_loc(n, 2) = dy_inv*(pt_loc(n,2)-ymn)
     pt_loc(n, 3) = dz_inv*(pt_loc(n,3)-zmn)
    end do
   else
    call map3d_part_sind( pt_loc, np, 2, 3 )
   end if
   if (n1==1) return
   do n = 1, np
    pt_loc(n, 4) = dx_inv*(pt_loc(n,4)-xmn)
   end do
   if (n_str==0) then
    do n = 1, np
     pt_loc(n, 5) = dy_inv*(pt_loc(n,5)-ymn)
     pt_loc(n, 6) = dz_inv*(pt_loc(n,6)-zmn)
    end do
   else
    call map3d_part_sind( pt_loc, np, 5, 6 )
   end if

  end subroutine

 end module
