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

 character(1), parameter :: X_COMP = 'x'
 character(1), parameter :: Y_COMP = 'y'
 character(1), parameter :: Z_COMP = 'z'
 character(2), parameter :: PX_COMP = 'px'
 character(2), parameter :: PY_COMP = 'py'
 character(2), parameter :: PZ_COMP = 'pz'
 character(6), parameter :: W_COMP = 'weight'
 character(5), parameter :: INDEX_COMP = 'index'

 type species
 real (dp), allocatable :: part(:, :)
 end type
 
 type species_new
  logical :: initialized
  !! Flag that states if the species has been initialized
  integer :: charge
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
   procedure, public :: how_many
   procedure, public :: total_size
   procedure, public :: initialize_component_real
   procedure, public :: initialize_component_integer
   generic :: initialize_component => initialize_component_real, initialize_component_integer
 end type

 interface species_new
  module procedure :: new_species_new
 end interface

 contains

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
  class(species_new), intent(inout) :: this
  real (dp), intent(in) :: values(:)
  character(*), intent(in) :: component

  select case(component)
  case('x')
   this%x(:) = values(:)
  case('y')
   this%y(:) = values(:)
  case('z')
   this%z(:) = values(:)
  case('px')
   this%px(:) = values(:)
  case('py')
   this%py(:) = values(:)
  case('pz')
   this%pz(:) = values(:)
  case('weight')
   this%weight(:) = values(:)
  case default
   stop
  end select

 end subroutine

 subroutine initialize_component_integer( this, values, component )
  class(species_new), intent(inout) :: this
  integer, intent(in) :: values(:)
  character(*), intent(in) :: component

  select case(component)
  case('x')
   this%x(:) = real(values(:), dp)
  case('y')
   this%y(:) = real(values(:), dp)
  case('z')
   this%z(:) = real(values(:), dp)
  case('px')
   this%px(:) = real(values(:), dp)
  case('py')
   this%py(:) = real(values(:), dp)
  case('pz')
   this%pz(:) = real(values(:), dp)
  case('weight')
   this%weight(:) = real(values(:), sp)
  case('index')
   this%part_index(:) = values(:)
  end select

 end subroutine

end module