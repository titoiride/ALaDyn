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

module particles_aux_def

 use particles_def
 implicit none
 public

 integer, parameter :: EX_COMP = 1
 integer, parameter :: EY_COMP = 2
 integer, parameter :: EZ_COMP = 3
 integer, parameter :: BX_COMP = 4
 integer, parameter :: BY_COMP = 5
 integer, parameter :: BZ_COMP = 6

 integer, parameter :: AUX1_COMP = 7
 integer, parameter :: AUX2_COMP = 8
 integer, parameter :: AUX3_COMP = 9
 integer, parameter :: AUX4_COMP = 10
 integer, parameter :: AUX5_COMP = 11
 integer, parameter :: AUX6_COMP = 12
 integer, parameter :: AUX7_COMP = 13

 integer, parameter :: OLD_X_COMP = 14
 integer, parameter :: OLD_Y_COMP = 15
 integer, parameter :: OLD_Z_COMP = 16
 integer, parameter :: OLD_PX_COMP = 17
 integer, parameter :: OLD_PY_COMP = 18
 integer, parameter :: OLD_PZ_COMP = 19
 integer, parameter :: OLD_GAMMA_COMP = 20

 type, extends(species_new) :: species_aux
  !! Auxiliary species for operations on species type

  real(dp), allocatable :: aux1(:)
  logical :: allocated_aux1
  !! True if aux1 array is allocated

  real(dp), allocatable :: aux2(:)
  logical :: allocated_aux2
  !! True if aux2 array is allocated

  real(dp), allocatable :: aux3(:)
  logical :: allocated_aux3
  !! True if aux3 array is allocated

  real(dp), allocatable :: aux4(:)
  !! Bx field interpolated on particles position
  logical :: allocated_aux4
  !! True if aux4 array is allocated

  real(dp), allocatable :: aux5(:)
  !! By field interpolated on particles position
  logical :: allocated_aux5
  !! True if aux5 array is allocated

  real(dp), allocatable :: aux6(:)
  !! Bz field interpolated on particles position
  logical :: allocated_aux6
  !! True if aux6 array is allocated

  real(dp), allocatable :: aux7(:)
  !! Bz field interpolated on particles position
  logical :: allocated_aux7
  !! True if aux6 array is allocated

  contains
   procedure, private :: set_component_aux
   procedure, public :: call_component => call_component_aux
   procedure, public :: copy_scalars_from => copy_scalars_from_aux

 end type species_aux

 contains

 subroutine copy_scalars_from_aux( this, other )
  !! Copies all the non-array values from a `species_new` to another
  class(species_aux), intent(out) :: this
  type(species_new), intent(in) :: other

  this%charge = other%charge
  this%dimensions = other%dimensions
  this%initialized = other%initialized
  this%n_part = other%n_part
  this%allocated_x = other%allocated_x
  this%allocated_y = other%allocated_y
  this%allocated_z = other%allocated_z
  this%allocated_px = other%allocated_px
  this%allocated_py = other%allocated_py
  this%allocated_pz = other%allocated_pz
  this%allocated_gamma = other%allocated_gamma
  this%allocated_weight = other%allocated_weight
  this%allocated_index = other%allocated_index

 end subroutine

 pure function call_component_aux( this, component, lb, ub ) result(comp)
 !! Function that hides the underlying array and calls the
 !! corresponding component from the particle structure.

  class(species_aux), intent(in) :: this
  integer, intent(in) :: component
  integer, intent(in), optional :: lb, ub
  real(dp), allocatable, dimension(:) :: comp
  integer :: lowb, upb, n_parts

  n_parts = this%n_part

  if ( present(lb) ) then
   lowb = lb
  else
   lowb = 1
  end if

  if ( present(ub) ) then
   upb = ub
  else
   upb = n_parts
  end if

  ! WARNING: allocation status should be checked
  select case(component)
  case(EX_COMP)
   comp = this%aux1(lowb:upb)
  case(EY_COMP)
   comp = this%aux2(lowb:upb)
  case(EZ_COMP)
   comp = this%aux3(lowb:upb)
  case(BX_COMP)
   comp = this%aux4(lowb:upb)
  case(BY_COMP)
   comp = this%aux5(lowb:upb)
  case(BZ_COMP)
   comp = this%aux6(lowb:upb)

  case(AUX1_COMP)
   comp = this%aux1(lowb:upb)
  case(AUX2_COMP)
   comp = this%aux2(lowb:upb)
  case(AUX3_COMP)
   comp = this%aux3(lowb:upb)
  case(AUX4_COMP)
   comp = this%aux4(lowb:upb)
  case(AUX5_COMP)
   comp = this%aux5(lowb:upb)
  case(AUX6_COMP)
   comp = this%aux6(lowb:upb)
  case(AUX7_COMP)
   comp = this%aux7(lowb:upb)

  case(OLD_X_COMP)
   comp = this%aux1(lowb:upb)
  case(OLD_Y_COMP)
   comp = this%aux2(lowb:upb)
  case(OLD_Z_COMP)
   comp = this%aux3(lowb:upb)
  case(OLD_PX_COMP)
   comp = this%aux4(lowb:upb)
  case(OLD_PY_COMP)
   comp = this%aux5(lowb:upb)
  case(OLD_PZ_COMP)
   comp = this%aux6(lowb:upb)
  case(OLD_GAMMA_COMP)
   comp = this%aux7(lowb:upb)

  end select

 end function

 subroutine set_component_aux( this, values, component, lb, ub )
  !! Assigns an array of real values to a given `species_new` component
  class(species_aux), intent(out) :: this
  real (dp), intent(in) :: values(:)
  integer, intent(in) :: component
  integer, intent(in), optional :: lb, ub
  integer :: lowb, upb

  if ( present(lb) ) then
   lowb = lb
  else
   lowb = lbound(this%call_component( component ), 1)
  end if

  if ( present(ub) ) then
   upb = ub
  else
   upb = ubound(this%call_component( component ), 1)
  end if

  select case(component)
  case(EX_COMP)
   this%aux1(lowb:upb) = values(:)
  case(EY_COMP)
   this%aux2(lowb:upb) = values(:)
  case(EZ_COMP)
   this%aux3(lowb:upb) = values(:)
  case(BX_COMP)
   this%aux4(lowb:upb) = values(:)
  case(BY_COMP)
   this%aux5(lowb:upb) = values(:)
  case(BZ_COMP)
   this%aux6(lowb:upb) = values(:)

  case(AUX1_COMP)
   this%aux1(lowb:upb) = values(:)
  case(AUX2_COMP)
   this%aux2(lowb:upb) = values(:)
  case(AUX3_COMP)
   this%aux3(lowb:upb) = values(:)
  case(AUX4_COMP)
   this%aux4(lowb:upb) = values(:)
  case(AUX5_COMP)
   this%aux5(lowb:upb) = values(:)
  case(AUX6_COMP)
   this%aux6(lowb:upb) = values(:)
  case(AUX7_COMP)
   this%aux7(lowb:upb) = values(:)

  case(OLD_X_COMP)
   this%aux1(lowb:upb) = values(:)
  case(OLD_Y_COMP)
   this%aux2(lowb:upb) = values(:)
  case(OLD_Z_COMP)
   this%aux3(lowb:upb) = values(:)
  case(OLD_PX_COMP)
   this%aux4(lowb:upb) = values(:)
  case(OLD_PY_COMP)
   this%aux5(lowb:upb) = values(:)
  case(OLD_PZ_COMP)
   this%aux6(lowb:upb) = values(:)
  case(OLD_GAMMA_COMP)
   this%aux7(lowb:upb) = values(:)

  end select

 end subroutine
end module