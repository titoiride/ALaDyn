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

module base_species

 use precision_def
 implicit none
 public

 integer, parameter :: X_COMP = 1
 integer, parameter :: Y_COMP = 2
 integer, parameter :: Z_COMP = 3
 integer, parameter :: PX_COMP = 4
 integer, parameter :: PY_COMP = 5
 integer, parameter :: PZ_COMP = 6
 integer, parameter :: INV_GAMMA_COMP = 7
 integer, parameter :: W_COMP = 8
 integer, parameter :: INDEX_COMP = -1
 
 type, abstract :: base_species_T

  logical :: initialized
  !! Flag that states if the species has been initialized
  real(dp) :: charge
  !! Particle charge
  integer, public :: n_part
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
  
  real(dp), allocatable :: gamma_inv(:)
  !! Array containig the inverse of Lorentz gamma factor
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
   procedure, pass :: compute_gamma
   procedure, pass :: count_particles
   procedure, pass :: how_many
   procedure, pass :: new_species => new_species_abstract
   procedure, pass :: set_charge_int
   procedure, pass :: set_charge_real
   procedure, pass :: set_part_number
   procedure(call_component_abstract), deferred, pass :: call_component
   procedure(copy_scalars_abstract), deferred, pass :: copy_scalars_from
   procedure(set_component_abstract_real), deferred, pass :: set_component_real
   procedure(set_component_abstract_integer), deferred, pass :: set_component_integer
   generic :: set_component => set_component_real, set_component_integer
   generic :: set_charge => set_charge_int, set_charge_real
  end type

  
  abstract interface
   pure function call_component_abstract( this, component, lb, ub ) result(comp)
    import :: base_species_T, dp
    implicit none
    class(base_species_T), intent(in)                                     :: this
    integer,       intent(in)                                     :: component
    integer,       intent(in), optional                           :: lb, ub
    real(dp),                           allocatable, dimension(:) :: comp
   end function
  end interface

  abstract interface
   subroutine set_component_abstract_real( this, values, component, lb, ub )
    import :: base_species_T, dp
    implicit none
    class(base_species_T), intent(inout)                     :: this
    real(dp),      intent(in), dimension(:)          :: values
    integer,       intent(in)                        :: component
    integer,       intent(in),              optional :: lb, ub
   end subroutine

   subroutine set_component_abstract_integer( this, values, component, lb, ub )
    import :: base_species_T, dp
    implicit none
    class(base_species_T), intent(inout)                     :: this
    integer,       intent(in), dimension(:)          :: values
    integer,       intent(in)                        :: component
    integer,       intent(in),              optional :: lb, ub
   end subroutine
  end interface

  abstract interface
   subroutine copy_scalars_abstract( this, other )
    import base_species_T, dp
    implicit none
    class(base_species_T), intent(inout) :: this
    class(base_species_T),  intent(in)  :: other
   end subroutine
  end interface

  contains
  
  !==== Constructor ===
  subroutine new_species_abstract( this, n_particles, curr_ndims )
   !! Constructor for the `species_new` type
   class(base_species_T), intent(inout) :: this
   integer, intent(in) :: n_particles, curr_ndims
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
    allocate( this%gamma_inv(n_particles), stat=allocstatus)
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
    allocate( this%gamma_inv(n_particles), stat=allocstatus)
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
    allocate( this%gamma_inv(n_particles), stat=allocstatus)
    this%allocated_gamma = .true.
    allocate( this%weight(n_particles), stat=allocstatus)
    this%allocated_weight = .true.
    allocate( this%part_index(n_particles), stat=allocstatus)
    this%allocated_index = .true.
   end select
  end subroutine

!=== Type bound procedures
  subroutine compute_gamma( this, pond_pot )
   class(base_species_T), intent(inout) :: this
   real(dp), intent(in), optional :: pond_pot(:)
   real(dp), allocatable :: temp(:)
   integer :: np

   np = this%how_many()
   allocate( temp(np), source=zero_dp )

   if ( this%allocated_px ) then
    temp(1:np) = temp(1:np) + this%call_component( PX_COMP )*this%call_component( PX_COMP )
   end if

   if ( this%allocated_py ) then
    temp(1:np) = temp(1:np) + this%call_component( PY_COMP )*this%call_component( PY_COMP )
   end if

   if ( this%allocated_pz ) then
    temp(1:np) = temp(1:np) + this%call_component( PZ_COMP )*this%call_component( PZ_COMP )
   end if

   if( present( pond_pot ) ) then
    temp(1:np) = temp(1:np) + pond_pot(1:np)
   end if

   temp(1:np) = one_dp/sqrt(one_dp + temp(1:np))
   call this%set_component(temp(1:np), INV_GAMMA_COMP, lb=1, ub=np)

  end subroutine

  pure function count_particles( this ) result( number )
   class(base_species_T), intent(in) :: this
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
    number = SIZE( this%gamma_inv(:) )
   else if ( this%allocated_index ) then
    number = SIZE( this%part_index(:) )
   else if ( this%allocated_weight ) then
    number = SIZE( this%weight(:) )
   else
    number = 0
   end if

  end function

  
  pure function how_many( this ) result(n_parts)
   !! Number of particles in the species
   class(base_species_T), intent(in) :: this
   integer :: n_parts

   n_parts = this%n_part
  
  end function

  subroutine set_charge_int( this, ch)
   class(base_species_T), intent(inout) :: this
   integer, intent(in) :: ch
   
   this%charge = real(ch, dp)
  end subroutine
  
  subroutine set_charge_real( this, ch)
   class(base_species_T), intent(inout) :: this
   real(dp), intent(in) :: ch
   
   this%charge = ch
  end subroutine

  
 subroutine set_part_number( this, n_parts)
  class(base_species_T), intent(inout) :: this
  integer, intent(in) :: n_parts

  this%n_part = n_parts
 end subroutine

!==== Procedures not bound to type ======

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