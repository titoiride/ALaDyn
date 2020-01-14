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
   procedure, public :: initialize_component_real
   procedure, public :: initialize_component_integer
   generic :: initialize_component => initialize_component_real, initialize_component_integer
   generic :: set_charge => set_charge_real, set_charge_int
 end type

 interface species_new
  module procedure :: new_species_new
 end interface

 contains

 pure function call_component( this, component ) result(comp)
  class(species_new), intent(in) :: this
  integer, intent(in) :: component
  real(dp), allocatable, dimension(:) :: comp
  integer :: n_parts

  n_parts = this%n_part

  select case(component)
  case(X_COMP)
   comp = this%x(1:n_parts)
  case(Y_COMP)
   comp = this%y(1:n_parts)
  case(Z_COMP)
   comp = this%z(1:n_parts)
  case(PX_COMP)
   comp = this%px(1:n_parts)
  case(PY_COMP)
   comp = this%py(1:n_parts)
  case(PZ_COMP)
   comp = this%pz(1:n_parts)
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

 subroutine initialize_component_real( this, values, component )
  class(species_new), intent(out) :: this
  real (dp), intent(in) :: values(:)
  integer, intent(in) :: component

  select case(component)
  case(X_COMP)
   this%x(:) = values(:)
  case(Y_COMP)
   this%y(:) = values(:)
  case(Z_COMP)
   this%z(:) = values(:)
  case(PX_COMP)
   this%px(:) = values(:)
  case(PY_COMP)
   this%py(:) = values(:)
  case(PZ_COMP)
   this%pz(:) = values(:)
  case(W_COMP)
   this%weight(:) = values(:)
  case default
   stop
  end select

 end subroutine

 subroutine initialize_component_integer( this, values, component )
  class(species_new), intent(out) :: this
  integer, intent(in) :: values(:)
  integer, intent(in) :: component

  select case(component)
  case(X_COMP)
   this%x(:) = real(values(:), dp)
  case(Y_COMP)
   this%y(:) = real(values(:), dp)
  case(Z_COMP)
   this%z(:) = real(values(:), dp)
  case(PX_COMP)
   this%px(:) = real(values(:), dp)
  case(PY_COMP)
   this%py(:) = real(values(:), dp)
  case(PZ_COMP)
   this%pz(:) = real(values(:), dp)
  case(W_COMP)
   this%weight(:) = real(values(:), sp)
  case(INDEX_COMP)
   this%part_index(:) = values(:)
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
 !! Dictionary for send_recieve routines in parallel.F90
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

end module