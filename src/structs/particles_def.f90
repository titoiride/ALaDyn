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

 use precision_def
 implicit none
 public

 integer, parameter :: X_COMP = 1
 integer, parameter :: Y_COMP = 2
 integer, parameter :: Z_COMP = 3
 integer, parameter :: PX_COMP = 4
 integer, parameter :: PY_COMP = 5
 integer, parameter :: PZ_COMP = 6
 integer, parameter :: GAMMA_COMP = 7
 integer, parameter :: W_COMP = 8
 integer, parameter :: INDEX_COMP = -1

 type species
 real (dp), allocatable :: part(:, :)
 end type
 
 type species_new

  logical :: initialized
  !! Flag that states if the species has been initialized
  real(dp) :: charge
  !! Particle charge
  integer, private :: n_part
  !! Number of particles
  integer :: dimensions
  !! Number of dimensions in which particles live

  real(dp), allocatable :: x(:)
  !! Array containig the x particle positions
  logical :: allocated_x
  !! True if x array is allocated

  real(dp), allocatable :: y(:)
  !! Array containig the y particle positions
  logical :: allocated_y
  !! True if y array is allocated

  real(dp), allocatable :: z(:)
  !! Array containig the z particle positions
  logical :: allocated_z
  !! True if z array is allocated

  real(dp), allocatable :: px(:)
  !! Array containig the x particle momenta
  logical :: allocated_px
  !! True if px array is allocated

  real(dp), allocatable :: py(:)
  !! Array containig the y particle momenta
  logical :: allocated_py
  !! True if py array is allocated

  real(dp), allocatable :: pz(:)
  !! Array containig the z particle momenta
  logical :: allocated_pz
  !! True if pz array is allocated
  
  real(dp), allocatable :: gamma(:)
  !! Array containig the Lorentz gamma factor
  logical :: allocated_gamma
  !! True if gamma array is allocated
  
  real (sp), allocatable :: weight(:)
  !! Array containig the particle weights
  logical :: allocated_weight
  !! True if weight array is allocated
  
  integer, allocatable :: part_index(:)
  !! Array containig the particle index
  logical :: allocated_index
  !! True if index array is allocated
  contains
   procedure, public :: call_component
   procedure, public :: count_particles
   procedure, private :: copy_scalars_from
   procedure, public :: flatten
   procedure, public :: how_many
   procedure, public :: sel_particles
   procedure, public :: total_size
   procedure, private :: pack_species_logical
   procedure, private :: pack_species_array
   procedure, private :: set_charge_int
   procedure, private :: set_charge_real
   procedure, public :: set_part_number
   procedure, private :: set_component_real
   procedure, private :: set_component_integer
   generic :: pack_species => pack_species_array, pack_species_logical
   generic :: set_component => set_component_real, set_component_integer
   generic :: set_charge => set_charge_real, set_charge_int
 end type

 interface species_new
  module procedure :: new_species_new
 end interface

 contains

 pure function count_particles( this ) result( number )
  class(species_new), intent(in) :: this
  integer :: number

  if ( this%allocated_x ) then
   number = SIZE( this%x(:) )
  else if ( this%allocated_y ) then
   number = SIZE( this%y(:) )
  else if ( this%allocated_z ) then
   number = SIZE( this%z(:) )
  else if ( this%allocated_px ) then
   number = SIZE( this%px(:) )
  else if ( this%allocated_py ) then
   number = SIZE( this%py(:) )
  else if ( this%allocated_pz ) then
   number = SIZE( this%pz(:) )
  else if ( this%allocated_gamma ) then
   number = SIZE( this%gamma(:) )
  else if ( this%allocated_index ) then
   number = SIZE( this%part_index(:) )
  else if ( this%allocated_weight ) then
   number = SIZE( this%weight(:) )
  end if

 end function

 subroutine copy_scalars_from( this, other )
  !! Copies all the non-array values from a `species_new` to another
  class(species_new), intent(out) :: this
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

 pure function call_component( this, component, lb, ub ) result(comp)
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
  case(GAMMA_COMP)
   comp = this%gamma(lowb:upb)
  case(W_COMP)
   comp = real( this%weight(lowb:upb), dp )
  case(INDEX_COMP)
   comp = real( this%part_index(lowb:upb), dp )
  end select

 end function

 pure function flatten( this ) result(flat_array)
  class(species_new), intent(in) :: this
  integer :: array_size, num_comps, i
  real(dp), allocatable :: temp(:, :), flat_array(:)

  array_size = this%count_particles()
  num_comps = this%total_size()
  allocate(temp( array_size, num_comps ))

  i = 1
  if( this%allocated_x ) then
   temp( :, i ) = this%x(:)
   i = i + 1
  end if
  if( this%allocated_y ) then
   temp( :, i ) = this%y(:)
   i = i + 1
  end if
  if( this%allocated_z ) then
   temp( :, i ) = this%z(:)
   i = i + 1
  end if
  if( this%allocated_px ) then
   temp( :, i ) = this%px(:)
   i = i + 1
  end if
  if( this%allocated_py ) then
   temp( :, i ) = this%py(:)
   i = i + 1
  end if
  if( this%allocated_pz ) then
   temp( :, i ) = this%pz(:)
   i = i + 1
  end if
  if( this%allocated_gamma ) then
   temp( :, i ) = this%gamma(:)
   i = i + 1
  end if

  flat_array = PACK( temp(:, :), .true. )

 end function

 pure function how_many( this ) result(n_parts)
  !! Number of particles in the species
  class(species_new), intent(in) :: this
  integer :: n_parts

  n_parts = this%n_part

 end function

 function sel_particles( this, lower_bound, upper_bound ) result(sel)
 !! Function that selects particles with respect to the given array boundaries
 !! (Memory position, NOT a particle index)
  class(species_new), intent(in) :: this
  integer, intent(in) :: lower_bound, upper_bound
  type(species_new) :: sel


  call sel%copy_scalars_from(this)

  if( this%allocated_x ) then
   sel%x = this%x(lower_bound:upper_bound)
  end if
  if( this%allocated_y ) then
   sel%y = this%y(lower_bound:upper_bound)
  end if
  if( this%allocated_z ) then
   sel%z = this%z(lower_bound:upper_bound)
  end if
  if( this%allocated_px ) then
   sel%px = this%px(lower_bound:upper_bound)
  end if
  if( this%allocated_py ) then
   sel%py = this%py(lower_bound:upper_bound)
  end if
  if( this%allocated_pz ) then
   sel%pz = this%pz(lower_bound:upper_bound)
  end if
  if( this%allocated_gamma ) then
   sel%gamma = this%gamma(lower_bound:upper_bound)
  end if
  if( this%allocated_weight ) then
   sel%weight = this%weight(lower_bound:upper_bound)
  end if
  if( this%allocated_index ) then
   sel%part_index = this%part_index(lower_bound:upper_bound)
  end if

  sel%n_part = sel%count_particles()

 end function

 pure function total_size( this ) result(size)
  class(species_new), intent(in) :: this
  integer :: size, i
  i = 0
  if( this%allocated_x ) then
   i = i + 1
  end if
  if( this%allocated_y ) then
   i = i + 1
  end if
  if( this%allocated_z ) then
   i = i + 1
  end if
  if( this%allocated_px ) then
   i = i + 1
  end if
  if( this%allocated_py ) then
   i = i + 1
  end if
  if( this%allocated_pz ) then
   i = i + 1
  end if
  if( this%allocated_gamma ) then
   i = i + 1
  end if
  size = i

 end function

 function new_species_new( n_particles, curr_ndims ) result(this)
  !! Constructor for the `species_new` type
  integer, intent(in) :: n_particles, curr_ndims
  type(species_new) :: this
  integer :: allocstatus
 
  if (n_particles <= 0) then
   this%initialized = .false.
   this%n_part = 0
   return
  end if
 
  this%initialized = .true.
  this%n_part = n_particles
  this%dimensions = curr_ndims
  this%allocated_x = .false.
  this%allocated_y = .false.
  this%allocated_z = .false.
  this%allocated_px = .false.
  this%allocated_py = .false.
  this%allocated_pz = .false.
  this%allocated_gamma = .false.
  this%allocated_weight = .false.
  this%allocated_index = .false.
 
  select case(curr_ndims)
  
  case(1)
  
   allocate( this%x(n_particles), stat=allocstatus)
   this%allocated_x = .true.
   allocate( this%px(n_particles), stat=allocstatus)
   this%allocated_px = .true.
   allocate( this%gamma(n_particles), stat=allocstatus)
   this%allocated_gamma = .true.
   allocate( this%weight(n_particles), stat=allocstatus)
   this%allocated_weight = .true.
   allocate( this%part_index(n_particles), stat=allocstatus)
   this%allocated_index = .true.
  case(2)
  
   allocate( this%x(n_particles), stat=allocstatus)
   this%allocated_x = .true.
   allocate( this%px(n_particles), stat=allocstatus)
   this%allocated_px = .true.
   allocate( this%y(n_particles), stat=allocstatus)
   this%allocated_y = .true.
   allocate( this%py(n_particles), stat=allocstatus)
   this%allocated_py = .true.
   allocate( this%gamma(n_particles), stat=allocstatus)
   this%allocated_gamma = .true.
   allocate( this%weight(n_particles), stat=allocstatus)
   this%allocated_weight = .true.
   allocate( this%part_index(n_particles), stat=allocstatus)
   this%allocated_index = .true.
  
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
   allocate( this%gamma(n_particles), stat=allocstatus)
   this%allocated_gamma = .true.
   allocate( this%weight(n_particles), stat=allocstatus)
   this%allocated_weight = .true.
   allocate( this%part_index(n_particles), stat=allocstatus)
   this%allocated_index = .true.
  end select
 end function

 function pack_species_logical( this, mask ) result(packed)
  class(species_new), intent(in) :: this
  logical, intent(in) :: mask
  type(species_new) :: packed

  call packed%copy_scalars_from(this)

  if( this%allocated_x ) then
   packed%x = PACK( this%x(:), mask)
  end if
  if( this%allocated_y ) then
   packed%y = PACK( this%y(:), mask)
  end if
  if( this%allocated_z ) then
   packed%z = PACK( this%z(:), mask)
  end if
  if( this%allocated_px ) then
   packed%px = PACK( this%px(:), mask)
  end if
  if( this%allocated_py ) then
   packed%py = PACK( this%py(:), mask)
  end if
  if( this%allocated_pz ) then
   packed%pz = PACK( this%pz(:), mask)
  end if
  if( this%allocated_gamma ) then
   packed%gamma = PACK( this%gamma(:), mask)
  end if
  if( this%allocated_weight ) then
   packed%weight = PACK( this%weight(:), mask)
  end if
  if( this%allocated_index ) then
   packed%part_index = PACK( this%part_index(:), mask)
  end if

  packed%n_part = packed%count_particles()
 end function

 function pack_species_array( this, mask ) result(packed)
  class(species_new), intent(in) :: this
  logical, intent(in) :: mask(:)
  type(species_new) :: packed

  call packed%copy_scalars_from(this)

  if( this%allocated_x ) then
   packed%x = PACK( this%x(:), mask(:) )
  end if
  if( this%allocated_y ) then
   packed%y = PACK( this%y(:), mask(:) )
  end if
  if( this%allocated_z ) then
   packed%z = PACK( this%z(:), mask(:) )
  end if
  if( this%allocated_px ) then
   packed%px = PACK( this%px(:), mask(:) )
  end if
  if( this%allocated_py ) then
   packed%py = PACK( this%py(:), mask(:) )
  end if
  if( this%allocated_pz ) then
   packed%pz = PACK( this%pz(:), mask(:) )
  end if
  if( this%allocated_gamma ) then
   packed%gamma = PACK( this%gamma(:), mask(:) )
  end if
  if( this%allocated_weight ) then
   packed%weight = PACK( this%weight(:), mask(:) )
  end if
  if( this%allocated_index ) then
   packed%part_index = PACK( this%part_index(:), mask(:) )
  end if

  packed%n_part = packed%count_particles()

 end function

 subroutine set_component_real( this, values, component, lb, ub )
  !! Assigns an array of real values to a given `species_new` component
  class(species_new), intent(out) :: this
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
   this%x(lowb:upb) = values(:)
  case(Y_COMP)
   this%y(lowb:upb) = values(:)
  case(Z_COMP)
   this%z(lowb:upb) = values(:)
  case(PX_COMP)
   this%px(lowb:upb) = values(:)
  case(PY_COMP)
   this%py(lowb:upb) = values(:)
  case(PZ_COMP)
   this%pz(lowb:upb) = values(:)
  case(GAMMA_COMP)
   this%gamma(lowb:upb) = values(:)
  case(W_COMP)
   this%weight(lowb:upb) = real(values(:), sp)
  end select

 end subroutine

 subroutine set_component_integer( this, values, component, lb, ub )
  !! Assigns an array of integer values to a given `species_new` component
  class(species_new), intent(out) :: this
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
   this%x(lowb:upb) = real(values(:), dp)
  case(Y_COMP)
   this%y(lowb:upb) = real(values(:), dp)
  case(Z_COMP)
   this%z(lowb:upb) = real(values(:), dp)
  case(PX_COMP)
   this%px(lowb:upb) = real(values(:), dp)
  case(PY_COMP)
   this%py(lowb:upb) = real(values(:), dp)
  case(PZ_COMP)
   this%pz(lowb:upb) = real(values(:), dp)
  case(GAMMA_COMP)
   this%gamma(lowb:upb) = real(values(:), dp)
  case(W_COMP)
   this%weight(lowb:upb) = real(values(:), sp)
  case(INDEX_COMP)
   this%part_index(lowb:upb) = values(:)
  end select

 end subroutine

 subroutine set_charge_int( this, ch)
  class(species_new), intent(out) :: this
  integer, intent(in) :: ch

  this%charge = real(ch, dp)
 end subroutine

 subroutine set_charge_real( this, ch)
  class(species_new), intent(out) :: this
  real(dp), intent(in) :: ch

  this%charge = ch
 end subroutine

 subroutine set_part_number( this, n_parts)
  class(species_new), intent(out) :: this
  integer, intent(in) :: n_parts

  this%n_part = n_parts
 end subroutine

 pure function component_dictionary( component ) result(cdir)
 !! Dictionary wrapper for send_recieve routines in parallel.F90
 !! that need the cdir parameter
  integer, intent(in) :: component
  integer :: cdir
  select case(component)
  case(X_COMP)
   cdir = 3
  case(Y_COMP)
   cdir = 1
  case(Z_COMP)
   cdir = 2
  end select
 end function

 pure function link_position_momentum( component ) result (pm_comp)
 !! Dictionary that gives back the corresponding position (momentum)
 !! when a momentum (position) is given, e.g.
 !!
 !! ```
 !!     comp = X_COMP
 !!
 !!     link_position_momentum( comp ) => gives PX_COMP
 !! ```
  integer, intent(in) :: component
  integer :: pm_comp
  select case(component)
  case(X_COMP)
   pm_comp = PX_COMP
  case(Y_COMP)
   pm_comp = PY_COMP
  case(Z_COMP)
   pm_comp = PZ_COMP
  case(PX_COMP)
   pm_comp = X_COMP
  case(PY_COMP)
   pm_comp = Y_COMP
  case(PZ_COMP)
   pm_comp = Z_COMP
  end select
 end function

end module