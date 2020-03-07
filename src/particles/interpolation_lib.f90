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

 module interpolation_lib

  use precision_def

  implicit none

  real (dp), parameter, private :: PP1 = 0.25, PP2 = 0.5, &
   PP3 = 0.75, PP4 = 2./3., PP5 = 1./6.

   type interp_coeff
   !! Type containing particle-grid relations during
   !! interpolation on grid
    integer :: order
    !! Interpolation order on integer cells
    integer :: h_order
    !! Interpolation order on half-integer cells
    real(dp), allocatable, dimension(:) :: coeff_x
    !! Interpolation coefficients along x (integer cells)
    real(dp), allocatable, dimension(:) :: coeff_y
    !! Interpolation coefficients along y (integer cells)
    real(dp), allocatable, dimension(:) :: coeff_z
    !! Interpolation coefficients along z (integer cells)
    real(dp), allocatable, dimension(:) :: h_coeff_x
    !! Interpolation coefficients along x (half-integer cells)
    real(dp), allocatable, dimension(:) :: h_coeff_y
    !! Interpolation coefficients along y (half-integer cells)
    real(dp), allocatable, dimension(:) :: h_coeff_z
    !! Interpolation coefficients along z (half-integer cells)
    real(dp), allocatable, dimension(:, :) :: coeff_x_rank2
    !! Interpolation coefficients along x for vectors (integer cells)
    real(dp), allocatable, dimension(:, :) :: coeff_y_rank2
    !! Interpolation coefficients along y for vectors (integer cells)
    real(dp), allocatable, dimension(:, :) :: coeff_z_rank2
    !! Interpolation coefficients along z for vectors (integer cells)
    real(dp), allocatable, dimension(:, :) :: h_coeff_x_rank2
    !! Interpolation coefficients along x for vectors (half-integer cells)
    real(dp), allocatable, dimension(:, :) :: h_coeff_y_rank2
    !! Interpolation coefficients along y for vectors (half-integer cells)
    real(dp), allocatable, dimension(:, :) :: h_coeff_z_rank2
    !! Interpolation coefficients along z for vectors (half-integer cells)
    integer :: ix
    !! Cell index along x
    integer :: ihx
    !! Half-cell index along x
    integer :: iy
    !! Cell index along y
    integer :: ihy
    !! Half-cell index along y
    integer :: iz
    !! Cell index along z
    integer :: ihz
    !! Half-cell index along z
    integer, allocatable, dimension(:) :: ix_rank2
    !! Cell index along x for vectors
    integer, allocatable, dimension(:) :: ihx_rank2
    !! Half-cell index along x for vectors
    integer, allocatable, dimension(:) :: iy_rank2
    !! Cell index along y for vectors
    integer, allocatable, dimension(:) :: ihy_rank2
    !! Half-cell index along y for vectors
    integer, allocatable, dimension(:) :: iz_rank2
    !! Cell index along z for vectors
    integer, allocatable, dimension(:) :: ihz_rank2
    !! Half-cell index along z for vectors
   end type

  interface zeroth_order
   module procedure :: zeroth_order_real
   module procedure :: zeroth_order_vector
  end interface

  interface first_order
   module procedure :: first_order_real
   module procedure :: first_order_vector
  end interface

  interface second_order
   module procedure :: second_order_real
   module procedure :: second_order_vector
  end interface

  interface third_order
   module procedure :: third_order_real
   module procedure :: third_order_vector
  end interface

  contains
  ! =======================================================
  !    Templates of spl=1,2,3 order shape functions for grid-particle  connection
  !==================================================
  ! Computes 
  subroutine set_int_pshape(spl, xx, ax, ind)
   integer, intent (in) :: spl
   real (dp), intent (in) :: xx
   real (dp), intent (out) :: ax(0:3)
   integer, intent (out) :: ind
   real (dp) :: sx, sx2, sx3
   !To integer grid points
   ax(0:3) = 0.0
   select case (spl)
   case (1)
    ind = int(xx)
    ax(1) = xx - real(ind, dp)
    ax(0) = one_dp - ax(1)
   case (2)
    ind = int(xx+0.5)
    sx = xx - real(ind, dp)
    sx2 = sx*sx
    ax(1) = PP3 - sx2
    ax(2) = PP2*(PP1+sx2+sx)
    ax(0) = one_dp - ax(1) - ax(2)
   case (3)
    ind = int(xx)
    sx = xx - real(ind, dp)
    sx2 = sx*sx
    sx3 = sx2*sx
    ax(1) = PP4 - sx2 + PP2*sx3
    ax(2) = PP5 + PP2*(sx+sx2-sx3)
    ax(3) = PP5*sx3
    ax(0) = one_dp - ax(1) - ax(2) - ax(3)
   end select
  end subroutine
  !========================
  subroutine set_hint_pshape(spl, xx, ax, ind)
   integer, intent (in) :: spl
   real (dp), intent (in) :: xx
   integer, intent (out) :: ind
   real (dp), intent (out) :: ax(0:3)
   real (dp) :: sx, sx2, sx3
   !To half-integer grid points
   ax(0:3) = 0.0
   select case (spl)
   case (1)
    sx = xx + 0.5
    ind = int(sx)
    ax(1) = sx - real(ind, dp)
    ax(0) = one_dp - ax(1)
   case (2)
    ind = int(xx)
    sx = xx - PP2 - real(ind, dp)
    sx2 = sx*sx
    ax(1) = PP3 - sx2
    ax(2) = PP2*(PP1+sx2+sx)
    ax(0) = one_dp - ax(1) - ax(2)
   case (3)
    ind = int(xx+0.5)
    sx = xx - real(ind, dp)
    sx2 = sx*sx
    sx3 = sx2*sx
    ax(1) = PP4 - sx2 + PP2*sx3
    ax(2) = PP5 + PP2*(sx+sx2-sx3)
    ax(3) = PP5*sx3
    ax(0) = one_dp - ax(1) - ax(2) - ax(3) ! to (i-1/2,i+1/2,i+3/2,i+5/2) half-int
   end select
  end subroutine
  !=======================
  ! End templates
  ! =======================================================
  pure function zeroth_order_real(deltax) result(ax)
   real(dp), intent(in) :: deltax
   real(dp) :: ax(1)
   ! Here dx should be computed as
   ! dx = x - int(x + 0.5)
   ax(1) = one_dp
  end function

  pure function zeroth_order_vector(deltax) result(ax)
   real(dp), intent(in), dimension(:) :: deltax
   real(dp), allocatable :: ax(:, :)
   integer :: length
   ! Here dx should be computed as
   ! dx = x - int(x + 0.5)
   length = SIZE( deltax )
   allocate(ax(length, 1))
   ax(1:length, 1) = one_dp
  end function

  pure function first_order_real(deltax) result(ax)
   real(dp), intent(in) :: deltax
   real(dp) :: ax(2)
   ! Here dx should be computed as
   ! dx = x - int(x)
   ax(2) = deltax
   ax(1) = one_dp - ax(2)
  end function

  pure function first_order_vector(deltax) result(ax)
   real(dp), intent(in), dimension(:) :: deltax
   real(dp), allocatable :: ax(:, :)
   integer :: length
   ! Here dx should be computed as
   ! dx = x - int(x)
   length = SIZE( deltax )
   allocate(ax(length, 2))
   ax(1:length, 2) = deltax
   ax(1:length, 1) = one_dp - ax(1:length, 2)
  end function

  pure function second_order_real(deltax) result(ax)
   real(dp), intent(in) :: deltax
   real(dp) :: ax(3)
   real(dp) :: dx2
   ! Here dx should be computed as
   ! dx = x - int(x + 0.5)
   dx2 = deltax*deltax
   ax(2) = PP3 - dx2
   ax(3) = PP2*(PP1 + deltax + dx2)
   ax(1) = one_dp - ax(2) - ax(3)
  end function

  pure function second_order_vector(deltax) result(ax)
   real(dp), intent(in), dimension(:) :: deltax
   real(dp), allocatable :: ax(:, :), dx2(:)
   integer :: length
   ! Here dx should be computed as
   ! dx = x - int(x)
   length = SIZE( deltax )
   allocate(ax(length, 3))
   allocate(dx2(length))
   dx2(1:length) = deltax(1:length)*deltax(1:length)
   ax(1:length, 2) = PP3 - dx2(1:length)
   ax(1:length, 3) = PP2*(PP1 + deltax(1:length) + dx2(1:length))
   ax(1:length, 1) = one_dp - ax(1:length, 2) - ax(1:length, 3)
  end function

  pure function third_order_real(deltax) result(ax)
   real(dp), intent(in) :: deltax
   real(dp) :: ax(4)
   real(dp) :: dx2, dx3
   ! Here dx should be computed as
   ! dx = x - int(x)
   dx2 = deltax*deltax
   dx3 = dx2*deltax

   ax(2) = PP4 - dx2 + PP2*dx3
   ax(3) = PP5 + PP2*(deltax + dx2 - dx3)
   ax(4) = PP5*dx3
   ax(1) = one_dp - ax(2) - ax(3) - ax(4)
  end function

  pure function third_order_vector(deltax) result(ax)
   real(dp), intent(in), dimension(:) :: deltax
   real(dp), allocatable :: ax(:, :), dx2(:), dx3(:)
   integer :: length
   ! Here dx should be computed as
   ! dx = x - int(x)
   length = SIZE( deltax )
   allocate(ax(length, 3))
   allocate(dx2(length))
   dx2(1:length) = deltax(1:length)*deltax(1:length)
   dx3(1:length) = dx2(1:length)*deltax(1:length)

   ax(1:length, 2) = PP4 - dx2(1:length) + PP2*dx3(1:length)
   ax(1:length, 3) = PP5 + PP2*(deltax(1:length) + dx2(1:length) - dx3(1:length))
   ax(1:length, 4) = PP5*dx3(1:length)
   ax(1:length, 1) = one_dp - ax(1:length, 2) - ax(1:length, 3) - ax(1:length, 4)
  end function

 end module
