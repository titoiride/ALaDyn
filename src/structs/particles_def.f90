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
 integer, parameter :: W_COMP = 7
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
  real(dp), allocatable :: x(:)
  !! Array containig the x particle positions
  real(dp), allocatable :: y(:)
  !! Array containig the y particle positions
  real(dp), allocatable :: z(:)
  !! Array containig the z particle positions
  real(dp), allocatable :: px(:)
  !! Array containig the x particle momenta
  real(dp), allocatable :: py(:)
  !! Array containig the y particle momenta
  real(dp), allocatable :: pz(:)
  !! Array containig the z particle momenta
  real (sp), allocatable :: weight(:)
  !! Array containig the particle weights
  integer, allocatable :: part_index(:)
  !! Array containig the particle index
  contains
   procedure, public :: call_component
   procedure, public :: how_many
   procedure, public :: total_size
   procedure, public :: set_charge_int
   procedure, public :: set_charge_real
   procedure, public :: set_part_number
   procedure, public :: set_component_real
   procedure, public :: set_component_integer
   generic :: set_component => set_component_real, set_component_integer
   generic :: set_charge => set_charge_real, set_charge_int
 end type

 interface species_new
  module procedure :: new_species_new
 end interface

 contains

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
  case(W_COMP)
   comp = real( this%weight(lowb:upb), dp )
  case(INDEX_COMP)
   comp = real( this%part_index(lowb:upb), dp )
  end select

 end function

 pure function how_many( this ) result(n_parts)
  class(species_new), intent(in) :: this
  integer :: n_parts

  n_parts = this%n_part

 end function

 pure function total_size( this ) result(size)
  class(species_new), intent(in) :: this
  integer :: size, i
  i = 0
  if(allocated(this%x)) then
   i = i + 1
  end if
  if(allocated(this%y)) then
   i = i + 1
  end if
  if(allocated(this%z)) then
   i = i + 1
  end if
  if(allocated(this%px)) then
   i = i + 1
  end if
  if(allocated(this%py)) then
   i = i + 1
  end if
  if(allocated(this%pz)) then
   i = i + 1
  end if
  size = i*this%how_many()

 end function

 function new_species_new( n_particles, curr_ndims ) result(this)
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
 
  select case(curr_ndims)
  
  case(1)
  
   allocate( this%x(n_particles), stat=allocstatus)
   allocate( this%px(n_particles), stat=allocstatus)
   allocate( this%weight(n_particles), stat=allocstatus)
   allocate( this%part_index(n_particles), stat=allocstatus)
  
  case(2)
  
   allocate( this%x(n_particles), stat=allocstatus)
   allocate( this%px(n_particles), stat=allocstatus)
   allocate( this%y(n_particles), stat=allocstatus)
   allocate( this%py(n_particles), stat=allocstatus)
   allocate( this%weight(n_particles), stat=allocstatus)
   allocate( this%part_index(n_particles), stat=allocstatus)
  
  case(3)
  
   allocate( this%x(n_particles), stat=allocstatus)
   allocate( this%px(n_particles), stat=allocstatus)
   allocate( this%y(n_particles), stat=allocstatus)
   allocate( this%py(n_particles), stat=allocstatus)
   allocate( this%z(n_particles), stat=allocstatus)
   allocate( this%pz(n_particles), stat=allocstatus)
   allocate( this%weight(n_particles), stat=allocstatus)
   allocate( this%part_index(n_particles), stat=allocstatus)
  end select
 end function

 subroutine set_component_real( this, values, component, lb, ub )
  class(species_new), intent(out) :: this
  real (dp), intent(in) :: values(:)
  integer, intent(in) :: component
  integer, intent(in), optional :: lb, ub
  integer :: lowb, upb

  if ( present(lb) ) then
   lowb = lb
  else
   lowb = lbound(this%call_component( component ))
  end if

  if ( present(ub) ) then
   upb = ub
  else
   upb = ubound(this%x)
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
  case(W_COMP)
   this%weight(lowb:upb) = values(:)
  end select

 end subroutine

 subroutine set_component_integer( this, values, component, lb, ub )
  class(species_new), intent(out) :: this
  integer, intent(in) :: values(:)
  integer, intent(in) :: component
  integer, intent(in), optional :: lb, ub
  integer :: lowb, upb

  if ( present(lb) ) then
   lowb = lb
  else
   lowb = lbound(this%call_component( component ))
  end if

  if ( present(ub) ) then
   upb = ub
  else
   upb = ubound(this%x)
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
 !! when a momentum (position) is given, e.g.7
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