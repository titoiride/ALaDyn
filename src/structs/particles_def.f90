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

module particles_def

 use base_species
 implicit none

 type species
  real (dp), allocatable :: part(:, :)
 end type
 
 type, extends(base_species_T) :: species_new
  contains
   procedure, public :: append => append_spec
   procedure, public :: call_component => call_component_spec
   procedure, public :: extend => extend_spec
   procedure, public :: new_species => new_species_spec
   procedure, public :: set_component_real
   procedure, public :: set_component_integer
   procedure, pass :: sel_particles_bounds => sel_particles_bounds_spec
   procedure, pass :: sel_particles_index => sel_particles_index_spec
   procedure, pass :: sweep => sweep_spec
 end type

 interface operator(*)
  module procedure :: multiply_number_spec
 end interface

 contains

 !========================================
 ! CONSTRUCTOR
 !========================================
 
 subroutine new_species_spec( this, n_particles, curr_ndims, tracked )
  !! Constructor for the `species_new` type
  class(species_new), intent(inout) :: this
  integer, intent(in) :: n_particles, curr_ndims
  logical, intent(in), optional :: tracked
  integer :: allocstatus
 
  if (n_particles < 0) then
   return
  end if
  if ( .not. allocated(this%initialized)) then
   allocate(this%initialized)
  end if
  this%initialized = .true.
  call this%set_part_number(n_particles)
  call this%set_dimensions(curr_ndims)
  if ( present(tracked) ) then
   call this%track( tracked )
  else
   call this%track( .false. )
  end if
  this%allocated_x = .false.
  this%allocated_y = .false.
  this%allocated_z = .false.
  this%allocated_px = .false.
  this%allocated_py = .false.
  this%allocated_pz = .false.
  this%allocated_gamma = .false.
  this%allocated_weight = .false.
  this%allocated_index = .false.
  if (n_particles == 0) then
   this%empty = .true.
   return
  end if
  this%empty = .false.
 
  select case(curr_ndims)
  
  case(1)
  
   allocate( this%x(n_particles), stat=allocstatus)
   this%allocated_x = .true.
   allocate( this%px(n_particles), stat=allocstatus)
   this%allocated_px = .true.
   allocate( this%gamma_inv(n_particles), stat=allocstatus)
   this%allocated_gamma = .true.
   allocate( this%weight(n_particles), stat=allocstatus)
   this%allocated_weight = .true.
   if (this%istracked()) then
    allocate( this%part_index(n_particles), stat=allocstatus)
    this%allocated_index = .true.
   end if
  case(2)
  
   allocate( this%x(n_particles), stat=allocstatus)
   this%allocated_x = .true.
   allocate( this%px(n_particles), stat=allocstatus)
   this%allocated_px = .true.
   allocate( this%y(n_particles), stat=allocstatus)
   this%allocated_y = .true.
   allocate( this%py(n_particles), stat=allocstatus)
   this%allocated_py = .true.
   allocate( this%gamma_inv(n_particles), stat=allocstatus)
   this%allocated_gamma = .true.
   allocate( this%weight(n_particles), stat=allocstatus)
   this%allocated_weight = .true.
   if (this%istracked()) then
    allocate( this%part_index(n_particles), stat=allocstatus)
    this%allocated_index = .true.
   end if
  
  case(3)
  
   allocate( this%x(n_particles), stat=allocstatus)
   this%allocated_x = .true.
   allocate( this%px(n_particles), stat=allocstatus)
   this%allocated_px = .true.
   allocate( this%y(n_particles), stat=allocstatus)
   this%allocated_y = .true.
   allocate( this%py(n_particles), stat=allocstatus)
   this%allocated_py = .true.
   allocate( this%z(n_particles), stat=allocstatus)
   this%allocated_z = .true.
   allocate( this%pz(n_particles), stat=allocstatus)
   this%allocated_pz = .true.
   allocate( this%gamma_inv(n_particles), stat=allocstatus)
   this%allocated_gamma = .true.
   allocate( this%weight(n_particles), stat=allocstatus)
   this%allocated_weight = .true.
   if (this%istracked()) then
    allocate( this%part_index(n_particles), stat=allocstatus)
    this%allocated_index = .true.
   end if
  end select
 end subroutine

 !========================================
 ! CLEANING PROCEDURES
 !========================================
 
 subroutine sweep_spec( this )
  class(species_new), intent(inout) :: this

  if ( .not. allocated(this%initialized) ) then
   return
  else
   deallocate( this%initialized )
  end if

  if ( this%allocated_x ) then
   this%allocated_x = .false.
   deallocate( this%x )
  end if
  if ( this%allocated_y ) then
   this%allocated_y = .false.
   deallocate( this%y )
  end if
  if ( this%allocated_z ) then
   this%allocated_z = .false.
   deallocate( this%z )
  end if
  if ( this%allocated_px ) then
   this%allocated_px = .false.
   deallocate( this%px )
  end if
  if ( this%allocated_py ) then
   this%allocated_py = .false.
   deallocate( this%py )
  end if
  if ( this%allocated_pz ) then
   this%allocated_pz = .false.
   deallocate( this%pz )
  end if
  if ( this%allocated_weight ) then
   this%allocated_weight = .false.
   deallocate( this%weight )
  end if
  if ( this%allocated_gamma ) then
   this%allocated_gamma = .false.
   deallocate( this%gamma_inv )
  end if
  if ( this%allocated_index ) then
   this%allocated_index = .false.
   deallocate( this%part_index )
  end if

 end subroutine
 !========================================
 ! TYPE BOUND PROCEDURES
 !========================================

 function append_spec( this, other ) result(spec)
  class(species_new), intent(inout) :: this
  class(base_species_T), intent(inout) :: other
  type(species_new) :: spec
  integer :: tot_size

  if ( .not. allocated(this%initialized) .and. .not. allocated(other%initialized) ) then
   return
  else if ( .not. allocated(this%initialized) ) then
   call this%new_species(0, other%pick_dimensions(), tracked=other%istracked())
  else if ( .not. allocated(other%initialized) ) then
   call other%new_species(0, this%pick_dimensions(), tracked=this%istracked())
  end if

  if ( allocated(this%initialized) .and. allocated(other%initialized) ) then
   if ( this%empty .and. .not. other%empty ) then
    call spec%copy(other, 1, other%how_many())
    return
   else if (other%empty .and. .not. this%empty) then
    call spec%copy(this, 1, this%how_many())
    return
   else if (other%empty .and. this%empty) then
    return
   end if
  end if
  
  tot_size = this%how_many()+other%how_many()
  call spec%new_species(tot_size, this%pick_dimensions(), tracked=this%istracked())
  
  if(other%allocated_x) then
   call assign(spec%x, [this%call_component(X_COMP), other%call_component(X_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_y) then
   call assign(spec%y, [this%call_component(Y_COMP), other%call_component(Y_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_z) then
   call assign(spec%z, [this%call_component(Z_COMP), other%call_component(Z_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_px) then
   call assign(spec%px, [this%call_component(PX_COMP), other%call_component(PX_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_py) then
   call assign(spec%py, [this%call_component(PY_COMP), other%call_component(PY_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_pz) then
   call assign(spec%pz, [this%call_component(PZ_COMP), other%call_component(PZ_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_gamma) then
   call assign(spec%gamma_inv, [this%call_component(INV_GAMMA_COMP), other%call_component(INV_GAMMA_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_weight) then
   call assign(spec%weight, [this%call_component(W_COMP), other%call_component(W_COMP)], &
    1, tot_size, tot_size)
  end if
  if(other%allocated_index) then
   call assign(spec%part_index, [this%call_component(INDEX_COMP), other%call_component(INDEX_COMP)], &
    1, tot_size, tot_size)
  end if
 end function

 pure function call_component_spec( this, component, lb, ub ) result(comp)
 !! Function that hides the underlying array and calls the
 !! corresponding component from the particle structure.
 !! @warning
 !! This function gives back always an array of reals!
 !! When using for weights and particle indexes remember to
 !! cast it again to the right type.
 !! @endwarning

  class(species_new), intent(in) :: this
  integer, intent(in) :: component
  integer, intent(in), optional :: lb, ub
  real(dp), allocatable, dimension(:) :: comp
  integer :: lowb, upb, n_parts

  n_parts = this%how_many()

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
  case(X_COMP)
   comp = this%x(lowb:upb)
  case(Y_COMP)
   comp = this%y(lowb:upb)
  case(Z_COMP)
   comp = this%z(lowb:upb)
  case(PX_COMP)
   comp = this%px(lowb:upb)
  case(PY_COMP)
   comp = this%py(lowb:upb)
  case(PZ_COMP)
   comp = this%pz(lowb:upb)
  case(INV_GAMMA_COMP)
   comp = this%gamma_inv(lowb:upb)
  case(W_COMP)
   comp = real( this%weight(lowb:upb), dp )
  case(INDEX_COMP)
   comp = real( this%part_index(lowb:upb), dp )
  end select

 end function

 subroutine extend_spec( this, new_number )
  class(species_new), intent(inout) :: this
  integer, intent(in) :: new_number
  type(species_new) :: temp
  integer :: n_parts, dimensions


  if ( .not. allocated(this%initialized) ) then
   write (6, *) 'Warning, cannot extend uninitialized species'
   return
  end if

  n_parts = this%how_many()
  if ( n_parts >= new_number ) then
   return
  end if
  call temp%copy(this)
  call this%sweep()
  call this%new_species(new_number, this%pick_dimensions(), tracked=this%istracked())
  call this%copy(temp)

 end subroutine

 subroutine sel_particles_bounds_spec( this, out_sp, lower_bound, upper_bound )
 !! Function that selects particles with respect to the given array boundaries
 !! (Memory position, NOT a particle index)
  class(species_new), intent(in) :: this
  class(base_species_T), intent(inout) :: out_sp
  integer, intent(in) :: lower_bound, upper_bound
  integer :: tot_len

  tot_len = upper_bound - lower_bound

  if ( allocated(out_sp%initialized) ) then
   call out_sp%sweep()
  end if

  call out_sp%new_species( tot_len, this%pick_dimensions(), tracked=this%istracked() )
  call out_sp%set_charge(this%pick_charge())

  if( this%allocated_x ) then
   out_sp%x = this%x(lower_bound:upper_bound)
  end if
  if( this%allocated_y ) then
   out_sp%y = this%y(lower_bound:upper_bound)
  end if
  if( this%allocated_z ) then
   out_sp%z = this%z(lower_bound:upper_bound)
  end if
  if( this%allocated_px ) then
   out_sp%px = this%px(lower_bound:upper_bound)
  end if
  if( this%allocated_py ) then
   out_sp%py = this%py(lower_bound:upper_bound)
  end if
  if( this%allocated_pz ) then
   out_sp%pz = this%pz(lower_bound:upper_bound)
  end if
  if( this%allocated_gamma ) then
   out_sp%gamma_inv = this%gamma_inv(lower_bound:upper_bound)
  end if
  if( this%allocated_weight ) then
   out_sp%weight = this%weight(lower_bound:upper_bound)
  end if
  if( this%allocated_index ) then
   out_sp%part_index = this%part_index(lower_bound:upper_bound)
  end if

  call out_sp%set_part_number(out_sp%array_size())

 end subroutine

 subroutine sel_particles_index_spec( this, out_sp, index_array )
  class(species_new), intent(in) :: this
  class(base_species_T), intent(inout) :: out_sp
  integer, dimension(:), intent(in) :: index_array
  integer :: i, tot_len, n

  tot_len = SIZE(index_array, DIM=1)

  if ( allocated(out_sp%initialized) ) then
   call out_sp%sweep()
  end if

  call out_sp%new_species( tot_len, this%pick_dimensions(), tracked=this%istracked() )
  call out_sp%set_charge(this%pick_charge())

  if( this%allocated_x ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%x(i) = this%x(n)
   end do
  end if
  if( this%allocated_y ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%y(i) = this%y(n)
   end do
  end if
  if( this%allocated_z ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%z(i) = this%z(n)
   end do
  end if

  if( this%allocated_px ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%px(i) = this%px(n)
   end do
  end if
  if( this%allocated_py ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%py(i) = this%py(n)
   end do
  end if
  if( this%allocated_pz ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%pz(i) = this%pz(n)
   end do
  end if
  if( this%allocated_gamma ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%gamma_inv(i) = this%gamma_inv(n)
   end do
  end if

  if( this%allocated_weight ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%weight(i) = this%weight(n)
   end do
  end if
  if( this%allocated_index ) then
   do i = 1, tot_len
    n = index_array(i)
    out_sp%part_index(i) = this%part_index(n)
   end do
  end if

 end subroutine

 subroutine set_component_real( this, values, component, lb, ub )
  !! Assigns an array of real values to a given `species_new` component
  class(species_new), intent(inout) :: this
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
  case(X_COMP)
   call assign(this%x, values, lowb, upb, this%how_many())
   this%allocated_x = .true.
  case(Y_COMP)
   call assign(this%y, values, lowb, upb, this%how_many())
   this%allocated_y = .true.
  case(Z_COMP)
   call assign(this%z, values, lowb, upb, this%how_many())
   this%allocated_z = .true.
  case(PX_COMP)
   call assign(this%px, values, lowb, upb, this%how_many())
   this%allocated_px = .true.
  case(PY_COMP)
   call assign(this%py, values, lowb, upb, this%how_many())
   this%allocated_py = .true.
  case(PZ_COMP)
   call assign(this%pz, values, lowb, upb, this%how_many())
   this%allocated_pz = .true.
  case(INV_GAMMA_COMP)
   call assign(this%gamma_inv, values, lowb, upb, this%how_many())
   this%allocated_gamma = .true.
  case(W_COMP)
   call assign(this%weight, values, lowb, upb, this%how_many())
   this%allocated_weight = .true.
  end select

 end subroutine

 subroutine set_component_integer( this, values, component, lb, ub )
  !! Assigns an array of integer values to a given `species_new` component
  class(species_new), intent(inout) :: this
  integer, intent(in) :: values(:)
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
  case(X_COMP)
   call assign(this%x, values, lowb, upb, this%how_many())
   this%allocated_x = .true.
  case(Y_COMP)
   call assign(this%y, values, lowb, upb, this%how_many())
   this%allocated_y = .true.
  case(Z_COMP)
   call assign(this%z, values, lowb, upb, this%how_many())
   this%allocated_z = .true.
  case(PX_COMP)
   call assign(this%px, values, lowb, upb, this%how_many())
   this%allocated_px = .true.
  case(PY_COMP)
   call assign(this%py, values, lowb, upb, this%how_many())
   this%allocated_py = .true.
  case(PZ_COMP)
   call assign(this%pz, values, lowb, upb, this%how_many())
   this%allocated_pz = .true.
  case(INV_GAMMA_COMP)
   call assign(this%gamma_inv, values, lowb, upb, this%how_many())
   this%allocated_gamma = .true.
  case(W_COMP)
   call assign(this%weight, values, lowb, upb, this%how_many())
   this%allocated_weight = .true.
  case(INDEX_COMP)
   call assign(this%part_index, values, lowb, upb, this%how_many())
   this%allocated_index = .true.
   call this%track( .true. )
  end select

 end subroutine
 !========================================
 ! NOT TYPE BOUND PROCEDURES
 !========================================

 function multiply_number_spec( this, number ) result(dot)
  type(species_new), intent(in) :: this
  real(dp), intent(in) :: number
  type(species_new) :: dot
  integer :: np
  !===========================================
  ! WARNING: THIS FUNCTION IS CLEARLY WRONG
  !===========================================
  np = this%how_many()
  call dot%new_species(this%how_many(), this%pick_dimensions(), tracked=this%istracked())
  call dot%set_charge(this%pick_charge())

  if( this%allocated_x ) then
   call assign(dot%x, number*this%call_component(X_COMP), 1, np)
  end if
  if( this%allocated_y ) then
   call assign(dot%y, number*this%call_component(Y_COMP), 1, np)
  end if
  if( this%allocated_z ) then
   call assign(dot%z, number*this%call_component(Z_COMP), 1, np)
  end if
  if( this%allocated_px ) then
   call assign(dot%px, number*this%call_component(PX_COMP), 1, np)
  end if
  if( this%allocated_py ) then
   call assign(dot%py, number*this%call_component(PY_COMP), 1, np)
  end if
  if( this%allocated_pz ) then
   call assign(dot%pz, number*this%call_component(PZ_COMP), 1, np)
  end if
  if( this%allocated_gamma ) then
   call assign(dot%gamma_inv, number*this%call_component(INV_GAMMA_COMP), 1, np)
  end if
  if( this%allocated_weight ) then
   call assign(dot%weight, number*this%call_component(W_COMP), 1, np)
  end if
  if( this%allocated_index ) then
   call assign(dot%part_index, number*this%call_component(INDEX_COMP), 1, np)
  end if

 end function
end module