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

 type, extends(species_new) :: species_aux
  !! Auxiliary species for operations on species type

  real(dp), allocatable :: ex(:)
  !! Ex field interpolated on particles position
  logical :: allocated_ex
  !! True if ex array is allocated

  real(dp), allocatable :: ey(:)
  !! Ey field interpolated on particles position
  logical :: allocated_ey
  !! True if ey array is allocated

  real(dp), allocatable :: ez(:)
  !! Ez field interpolated on particles position
  logical :: allocated_ez
  !! True if ez array is allocated

  real(dp), allocatable :: bx(:)
  !! Bx field interpolated on particles position
  logical :: allocated_bx
  !! True if bx array is allocated

  real(dp), allocatable :: by(:)
  !! By field interpolated on particles position
  logical :: allocated_by
  !! True if by array is allocated

  real(dp), allocatable :: bz(:)
  !! Bz field interpolated on particles position
  logical :: allocated_bz
  !! True if bz array is allocated

  contains

   procedure, public :: call_component => call_component_aux
   procedure, private :: copy_scalars_from => copy_scalars_from_aux

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
 !! @warning
 !! This function gives back always an array of reals!
 !! When using for weights and particle indexes remember to
 !! cast it again to the right type.
 !! @endwarning

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
   comp = this%ex(lowb:upb)
  case(EY_COMP)
   comp = this%ey(lowb:upb)
  case(EZ_COMP)
   comp = this%ez(lowb:upb)
  case(BX_COMP)
   comp = this%bx(lowb:upb)
  case(BY_COMP)
   comp = this%by(lowb:upb)
  case(BZ_COMP)
   comp = this%bz(lowb:upb)
  end select

 end function
end module