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

  implicit none

  real (dp), parameter :: half = 0.5
  real (sp), parameter :: shx = 3., shy = 3., shz = 3.
  integer (kind=2) :: err_ind

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

  contains

  !===========================================
  !     Computes particle density charge on grid integer points
  !==============================================
  !DIR$ ATTRIBUTES INLINE :: ql_interpolate,qq_interpolate
  !====================

   !=======================
   !  1D
   !=======================
  pure function qqh_1d_spline_real( xp ) result(interp)
    !! Quadratic interpolation at both integer and half-integer
    !! cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, ihx
   real (dp) :: xx, sx
   !======================

   interp%order = 2
   interp%h_order = 2

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   ihx = int(xx)
   sx = xx - real(ihx, dp)

   interp%h_coeff_x = second_order( sx )

   interp%ix = ix - 1
   interp%ihx = ihx - 1

  end function
  !=======================

  pure function qqh_1d_spline_vector( xp ) result(interp)
    !! Quadratic interpolation at both integer and half-integer
    !! cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   integer :: length
   real (dp), allocatable, dimension(:) :: xx, sx, ix, ihx
   !======================

   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( ihx(length) )
   allocate( sx(length) )

   interp%order = 2
   interp%h_order = 2

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int( xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   ihx(1:length) = int(xx(1:length))
   sx(1:length) = xx(1:length) - real(ihx(1:length), dp)

   interp%h_coeff_x_rank2 = second_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%ihx_rank2 = ihx(1:length) - 1

  end function
  !=======================

  pure function qden_1d_wgh_real( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix
   real (dp) :: xx, sx
   !======================
   interp%order = 2

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   interp%ix = ix - 1
  end function

  pure function qden_1d_wgh_vector( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   integer :: length
   real (dp), allocatable, dimension(:) :: xx, sx, ix
   !======================
   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( sx(length) )

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int( xx(1:length) + 0.5 )
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
  end function

  !=======================
  !  2D
  !=======================

  pure function qlh_2d_spline_real( xp ) result(interp)
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, ihx, iy, ihy
   real (dp) :: xx, sx
   !======================

   interp%order = 2
   interp%h_order = 1

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   interp%h_coeff_x = first_order( sx + half )

   xx = shy + xp(2)
   iy = int(xx+0.5)
   sx = xx - real(iy, dp)

   interp%coeff_y = second_order( sx )

   interp%h_coeff_y = first_order( sx + half )

   interp%ix = ix - 1
   interp%ihx = ix
   interp%iy = iy - 1
   interp%ihy = iy
   
  end function

  pure function qlh_2d_spline_vector( xp ) result(interp)
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   real (dp), allocatable, dimension(:) :: xx, sx, ix, iy
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( iy(length) )
   allocate( sx(length) )

   interp%order = 2
   interp%h_order = 1

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   interp%h_coeff_x_rank2 = first_order( sx(1:length) + half )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = second_order( sx(1:length) )

   interp%h_coeff_y_rank2 = first_order( sx(1:length) + half )

   interp%ix_rank2 = ix(1:length) - 1
   interp%ihx_rank2 = ix(1:length)
   interp%iy_rank2 = iy(1:length) - 1
   interp%ihy_rank2 = iy(1:length)
   
  end function
  !====================

  pure function qqh_2d_spline_real( xp ) result(interp)
   !! Quadratic interpolation at both integer and half-integer
   !! cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, ihx, iy, ihy
   real (dp) :: xx, sx
   !======================

   interp%order = 2
   interp%h_order = 2

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   ihx = int(xx)
   sx = xx - real(ihx, dp) - 0.5

   interp%h_coeff_x = second_order( sx )

   xx = shy + xp(2)
   iy = int(xx+0.5)
   sx = xx - real(iy, dp)

   interp%coeff_y = second_order( sx )

   ihy = int(xx)
   sx = xx - real(ihy, dp) - 0.5

   interp%h_coeff_y = second_order( sx )

   interp%ix = ix - 1
   interp%ihx = ihx - 1
   interp%iy = iy - 1
   interp%ihy = ihy - 1
  end function
  !=======================================

  pure function qqh_2d_spline_vector( xp ) result(interp)
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   real (dp), allocatable, dimension(:) :: xx, sx, ix, ihx, iy ,ihy
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( ihx(length) )
   allocate( iy(length) )
   allocate( ihy(length) )
   allocate( sx(length) )

   interp%order = 2
   interp%h_order = 2

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   ihx(1:length) = int(xx(1:length))
   sx(1:length) = xx(1:length) - real(ihx(1:length), dp) - 0.5

   interp%h_coeff_x_rank2 = second_order( sx(1:length) )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = second_order( sx(1:length) )

   ihy(1:length) = int(xx(1:length))
   sx(1:length) = xx(1:length) - real(ihy(1:length), dp) - 0.5

   interp%h_coeff_y_rank2 = second_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%ihx_rank2 = ihx(1:length) - 1
   interp%iy_rank2 = iy(1:length) - 1
   interp%ihy_rank2 = ihy(1:length) - 1
   
  end function
  !====================

  pure function qden_2d_wgh_real( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, iy
   real (dp) :: xx, sx
   !======================
   interp%order = 2

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   xx = shy + xp(2)
   iy = int(xx+0.5)
   sx = xx - real(iy, dp)

   interp%coeff_y = second_order( sx )

   interp%ix = ix - 1
   interp%iy = iy - 1
  end function
  !======================

  pure function qden_2d_wgh_vector( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   integer :: length
   real (dp), allocatable, dimension(:) :: xx, sx, ix, iy
   !======================
   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( iy(length) )
   allocate( sx(length) )

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int( xx(1:length) + 0.5 )
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int( xx(1:length) + 0.5 )
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = second_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%iy_rank2 = iy(1:length) - 1

  end function

  pure function cden_2d_wgh_real( xp ) result(interp)
   !! Cubic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, iy
   real (dp) :: xx, sx
   !======================
   interp%order = 3

   xx = shy + xp(1)
   ix = int(xx)
   sx = xx - real(ix, dp)

   interp%coeff_x = third_order( sx )

   xx = shy + xp(2)
   iy = int(xx)
   sx = xx - real(iy, dp)

   interp%coeff_y = third_order( sx )

   interp%ix = ix - 1
   interp%iy = iy - 1
  end function
  !==========================

  pure function cden_2d_wgh_vector( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   integer :: length
   real (dp), allocatable, dimension(:) :: xx, sx, ix, iy
   !======================
   length = SIZE( xp, DIM=1 )
   interp%order = 3

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( iy(length) )
   allocate( sx(length) )

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int( xx(1:length) )
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = third_order( sx(1:length) )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int( xx(1:length) )
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = third_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%iy_rank2 = iy(1:length) - 1

  end function
  !=======================
  !  3D
  !=======================
  pure function qlh_3d_spline_real( xp ) result(interp)
   !! Quadratic interpolation at integer cells, linear at
   !! half-integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, ihx, iy, ihy, iz, ihz
   real (dp) :: xx, sx
   !======================

   interp%order = 2
   interp%h_order = 1

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   interp%h_coeff_x = first_order( sx + half )

   xx = shy + xp(2)
   iy = int(xx+0.5)
   sx = xx - real(iy, dp)

   interp%coeff_y = second_order( sx )

   interp%h_coeff_y = first_order( sx + half )

   xx = shz + xp(3)
   iz = int(xx+0.5)
   sx = xx - real(iz, dp)

   interp%coeff_z = second_order( sx )

   interp%h_coeff_z = first_order( sx + half )

   interp%ix = ix - 1
   interp%ihx = ix
   interp%iy = iy - 1
   interp%ihy = iy
   interp%iz = iz - 1
   interp%ihz = iz
  end function
  !=================================

  pure function qlh_3d_spline_vector( xp ) result(interp)
   !! Quadratic interpolation at integer cells, linear at
   !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   real (dp), allocatable, dimension(:) :: xx, sx, ix, iy, iz
   integer :: length
   !======================
   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( iy(length) )
   allocate( iz(length) )
   allocate( sx(length) )

   interp%order = 2
   interp%h_order = 1

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   interp%h_coeff_x_rank2 = first_order( sx(1:length) + half )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = second_order( sx(1:length) )

   interp%h_coeff_y_rank2 = first_order( sx(1:length) + half )

   xx(1:length) = shz + xp(1:length, 3)
   iz(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(iz(1:length), dp)

   interp%coeff_z_rank2 = second_order( sx(1:length) )

   interp%h_coeff_z_rank2 = first_order( sx(1:length) + half )

   interp%ix_rank2 = ix(1:length) - 1
   interp%ihx_rank2 = ix(1:length)
   interp%iy_rank2 = iy(1:length) - 1
   interp%ihy_rank2 = iy(1:length)
   interp%iz_rank2 = iz(1:length) - 1
   interp%ihz_rank2 = iz(1:length)
  end function
  !=================================
  pure function qqh_3d_spline_real( xp ) result(interp)
   !! Quadratic interpolation at both integer and half-integer
   !! cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, ihx, iy, ihy, iz, ihz
   real (dp) :: xx, sx
   !=====================

   interp%order = 2
   interp%h_order = 2

   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   ihx = int(xx)
   sx = xx - real(ihx, dp) - 0.5

   interp%h_coeff_x = second_order( sx )

   xx = shy + xp(2)
   iy = int(xx+0.5)
   sx = xx - real(iy, dp)

   interp%coeff_y = second_order( sx )

   ihy = int(xx)
   sx = xx - real(ihy, dp) - 0.5

   interp%h_coeff_y = second_order( sx )

   xx = shz + xp(3)
   iz = int(xx+0.5)
   sx = xx - real(iz, dp)

   interp%coeff_z = second_order( sx )

   ihz = int(xx)
   sx = xx - real(ihz, dp) - 0.5

   interp%h_coeff_z = second_order( sx )

   interp%ix = ix - 1
   interp%ihx = ihx - 1
   interp%iy = iy - 1
   interp%ihy = ihy - 1
   interp%iz = iz - 1
   interp%ihz = ihz - 1
  end function

  pure function qqh_3d_spline_vector( xp ) result(interp)
  !! Quadratic interpolation at integer cells, linear at
  !! half-integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   real (dp), allocatable, dimension(:) :: xx, sx, ix, ihx, iy ,ihy, &
    iz, ihz
   integer :: length
   !======================

   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( ihx(length) )
   allocate( iy(length) )
   allocate( ihy(length) )
   allocate( iz(length) )
   allocate( ihz(length) )
   allocate( sx(length) )

   interp%order = 2
   interp%h_order = 2

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   ihx(1:length) = int(xx(1:length))
   sx(1:length) = xx(1:length) - real(ihx(1:length), dp) - 0.5

   interp%h_coeff_x_rank2 = second_order( sx(1:length) )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = second_order( sx(1:length) )

   ihy(1:length) = int(xx(1:length))
   sx(1:length) = xx(1:length) - real(ihy(1:length), dp) - 0.5

   interp%h_coeff_y_rank2 = second_order( sx(1:length) )

   xx(1:length) = shz + xp(1:length, 3)
   iz(1:length) = int(xx(1:length) + 0.5)
   sx(1:length) = xx(1:length) - real(iz(1:length), dp)

   interp%coeff_z_rank2 = second_order( sx(1:length) )

   ihz(1:length) = int(xx(1:length))
   sx(1:length) = xx(1:length) - real(ihz(1:length), dp) - 0.5

   interp%h_coeff_z_rank2 = second_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%ihx_rank2 = ihx(1:length) - 1
   interp%iy_rank2 = iy(1:length) - 1
   interp%ihy_rank2 = ihy(1:length) - 1
   interp%iz_rank2 = iz(1:length) - 1
   interp%ihz_rank2 = ihz(1:length) - 1
   
  end function
  !====================

  pure function qden_3d_wgh_real( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, iy, iz
   real (dp) :: xx, sx
   !======================
   xx = shx + xp(1)
   ix = int(xx+0.5)
   sx = xx - real(ix, dp)

   interp%coeff_x = second_order( sx )

   xx = shy + xp(2)
   iy = int(xx+0.5)
   sx = xx - real(iy, dp)

   interp%coeff_y = second_order( sx )

   xx = shz + xp(3)
   iz = int(xx+0.5)
   sx = xx - real(iz, dp)

   interp%coeff_z = second_order( sx )

   interp%ix = ix - 1
   interp%iy = iy - 1
   interp%iz = iz - 1
  end function
  !================================

  pure function qden_3d_wgh_vector( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   integer :: length
   real (dp), allocatable, dimension(:) :: xx, sx, ix, iy, iz
   !======================
   length = SIZE( xp, DIM=1 )
   interp%order = 2

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( iy(length) )
   allocate( iz(length) )
   allocate( sx(length) )

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int( xx(1:length) + 0.5 )
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = second_order( sx(1:length) )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int( xx(1:length) + 0.5 )
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = second_order( sx(1:length) )

   xx(1:length) = shz + xp(1:length, 3)
   iz(1:length) = int( xx(1:length) + 0.5 )
   sx(1:length) = xx(1:length) - real(iz(1:length), dp)

   interp%coeff_z_rank2 = second_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%iy_rank2 = iy(1:length) - 1
   interp%iz_rank2 = iz(1:length) - 1

  end function

  pure function cden_3d_wgh_real( xp ) result(interp)
   !! Cubic interpolation at integer cells
   real (dp), intent (in) :: xp(:)
   type(interp_coeff) :: interp
   integer :: ix, iy, iz
   real (dp) :: xx, sx
   !======================
   !cubic spline to integer index
   xx = shx + xp(1)
   ix = int(xx)
   sx = xx - real(ix, dp)

   interp%coeff_x = third_order( sx )

   xx = shy + xp(2)
   iy = int(xx)
   sx = xx - real(iy, dp)

   interp%coeff_y = third_order( sx )

   xx = shz + xp(3)
   iz = int(xx)
   sx = xx - real(iz, dp)

   interp%coeff_z = third_order( sx )

   interp%ix = ix - 1
   interp%iy = iy - 1
   interp%iz = iz - 1

  end function

  pure function cden_3d_wgh_vector( xp ) result(interp)
   !! Quadratic interpolation at integer cells
   real (dp), intent (in), dimension(:, :) :: xp
   type(interp_coeff) :: interp
   integer :: length
   real (dp), allocatable, dimension(:) :: xx, sx, ix, iy, iz
   !======================
   length = SIZE( xp, DIM=1 )
   interp%order = 3

   allocate( xx(length) )
   allocate( ix(length) )
   allocate( iy(length) )
   allocate( iz(length) )
   allocate( sx(length) )

   xx(1:length) = shx + xp(1:length, 1)
   ix(1:length) = int( xx(1:length) )
   sx(1:length) = xx(1:length) - real(ix(1:length), dp)

   interp%coeff_x_rank2 = third_order( sx(1:length) )

   xx(1:length) = shy + xp(1:length, 2)
   iy(1:length) = int( xx(1:length) )
   sx(1:length) = xx(1:length) - real(iy(1:length), dp)

   interp%coeff_y_rank2 = third_order( sx(1:length) )

   xx(1:length) = shz + xp(1:length, 3)
   iz(1:length) = int( xx(1:length) )
   sx(1:length) = xx(1:length) - real(iz(1:length), dp)

   interp%coeff_z_rank2 = third_order( sx(1:length) )

   interp%ix_rank2 = ix(1:length) - 1
   interp%iy_rank2 = iy(1:length) - 1
   interp%iz_rank2 = iz(1:length) - 1

  end function
  !====================================

  !DIR$ ATTRIBUTES INLINE :: ql_interpolate
  function set_local_positions( pt_array, component ) result(position)
   type (species_new), intent(in) :: pt_array
   integer, intent(in) :: component
   integer :: np
   real(dp), allocatable, dimension(:) :: position

   np = pt_array%how_many()
   position = pt_array%call_component( component )

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
