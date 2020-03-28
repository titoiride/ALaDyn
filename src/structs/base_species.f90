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
 use util, only: gasdev
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

  logical, allocatable :: initialized
  !! Flag that states if the species has been initialized
  logical :: empty
  !! Flag that states if there are particles
  real(dp) :: charge
  !! Particle charge
  integer, public :: n_part
  !! Number of particles
  integer :: dimensions
  !! Number of dimensions in which particles live
  real :: temperature
  !! Initial temperature given to the species

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
  
  real (dp), allocatable :: weight(:)
  !! Array containig the particle weights
  logical :: allocated_weight
  !! True if weight array is allocated
  
  integer, allocatable :: part_index(:)
  !! Array containig the particle index
  logical :: allocated_index
  !! True if index array is allocated

  contains
   procedure, public, pass :: add_data
   procedure, pass :: array_size
   procedure, private, pass :: call_particle_bounds
   procedure, private, pass :: call_particle_index_array
   procedure, private, pass :: call_particle_single
   procedure, pass :: compute_gamma
   procedure, pass, private :: copy_all
   procedure, pass, private :: copy_boundaries
   procedure, public, pass :: flatten
   procedure, pass :: how_many
   procedure, pass :: initialize_data
   procedure, pass :: new_species => new_species_abstract
   procedure, private :: pack_into_logical
   procedure, private :: pack_into_array
   procedure, public, pass :: reallocate
   procedure, public, pass :: redistribute
   procedure, pass :: set_charge_int
   procedure, pass :: set_charge_real
   procedure, pass :: set_part_number
   procedure, pass :: set_temperature
   procedure, public, pass :: total_size
   procedure(call_component_abstract), deferred, pass :: call_component
   procedure(copy_scalars_abstract), deferred, pass :: copy_scalars_from
   procedure(extend_abstract), deferred, pass :: extend
   procedure(sel_particles_bounds_abstract), deferred, pass :: sel_particles_bounds
   procedure(sel_particles_index_abstract), deferred, pass :: sel_particles_index
   procedure(set_component_abstract_real), deferred, pass :: set_component_real
   procedure(set_component_abstract_integer), deferred, pass :: set_component_integer
   procedure(sweep_abstract), deferred, pass :: sweep
   generic :: call_particle => call_particle_single, call_particle_bounds, call_particle_index_array
   generic :: copy => copy_all, copy_boundaries
   generic :: pack_into => pack_into_array, pack_into_logical
   generic :: sel_particles => sel_particles_bounds, sel_particles_index
   generic :: set_component => set_component_real, set_component_integer
   generic :: set_charge => set_charge_int, set_charge_real
  end type

  abstract interface
   pure function call_component_abstract( this, component, lb, ub ) result(comp)
    import :: base_species_T, dp
    implicit none
    class(base_species_T), intent(in) :: this
    integer, intent(in) :: component
    integer, intent(in), optional :: lb, ub
    real(dp), allocatable, dimension(:) :: comp
   end function
  end interface
  
  abstract interface
   subroutine copy_scalars_abstract( this, other )
    import base_species_T, dp
    implicit none
    class(base_species_T), intent(inout) :: this
    class(base_species_T),  intent(in)  :: other
   end subroutine
  end interface

  abstract interface
   subroutine extend_abstract( this, new_number )
    import base_species_T, dp
    class(base_species_T), intent(inout) :: this
    integer, intent(in) :: new_number
   end subroutine
  end interface

  abstract interface
   subroutine set_component_abstract_real( this, values, component, lb, ub )
    import :: base_species_T, dp
    implicit none
    class(base_species_T), intent(inout) :: this
    real(dp), intent(in), dimension(:) :: values
    integer, intent(in) :: component
    integer, intent(in), optional :: lb, ub
   end subroutine

   subroutine set_component_abstract_integer( this, values, component, lb, ub )
    import :: base_species_T, dp
    implicit none
    class(base_species_T), intent(inout) :: this
    integer, intent(in), dimension(:) :: values
    integer, intent(in) :: component
    integer, intent(in), optional :: lb, ub
   end subroutine
  end interface

  abstract interface
   subroutine sel_particles_index_abstract( this, out_sp, index_array )
    import base_species_T, dp
    class(base_species_T), intent(in) :: this
    class(base_species_T), intent(inout) :: out_sp
    integer, dimension(:), intent(in) :: index_array
   end subroutine
   
   subroutine sel_particles_bounds_abstract( this, out_sp, lower_bound, upper_bound )
    import base_species_T, dp
    class(base_species_T), intent(in) :: this
    class(base_species_T), intent(inout) :: out_sp
    integer, intent(in) :: lower_bound, upper_bound
   end subroutine

  end interface

  abstract interface

   subroutine sweep_abstract( this )
    import base_species_T, dp
    class(base_species_T), intent(inout) :: this
   end subroutine

  end interface

  interface assign
   module procedure :: assign_real_realdp
   module procedure :: assign_real_realsp
   module procedure :: assign_realsp_realsp
   module procedure :: assign_real_integer
   module procedure :: assign_integer_realdp
   module procedure :: assign_integer_realsp
   module procedure :: assign_integer_integer
  end interface
  contains
  
  !==== Constructor ===
  subroutine new_species_abstract( this, n_particles, curr_ndims )
   !! Constructor for the `species_new` type
   class(base_species_T), intent(inout) :: this
   integer, intent(in) :: n_particles, curr_ndims
   integer :: allocstatus
  
   if (n_particles < 0) then
    return
   end if
  
   if ( .not. allocated(this%initialized)) then
    allocate(this%initialized)
   end if
   this%initialized = .true.
   this%n_part = n_particles
   this%dimensions = curr_ndims
   if (n_particles == 0) then
    this%empty = .true.
    return
   end if
   this%empty = .false.
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
    !allocate( this%part_index(n_particles), stat=allocstatus)
    !this%allocated_index = .true.
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
    !allocate( this%part_index(n_particles), stat=allocstatus)
    !this%allocated_index = .true.
   
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
    !allocate( this%part_index(n_particles), stat=allocstatus)
    !this%allocated_index = .true.
   end select
  end subroutine

!=== Type bound procedures

  subroutine add_data( this, x_arr, y_arr, z_arr, &
   weightx_arr, weightyz_arr, loc_x0, loc_x, loc_y, loc_z, np_old)
   class( base_species_T), intent(inout) :: this
   real(dp), dimension(:), intent(in) :: x_arr, y_arr, z_arr
   real(dp), dimension(:), intent(in) :: weightx_arr
   real(dp), dimension(:, :), intent(in) :: weightyz_arr
   integer, intent(in) :: loc_x0, loc_x, loc_y, loc_z, np_old
   real(dp) :: u, t_x
   real(dp) :: whz
   integer :: p, dim, i, j, k

   t_x = this%temperature
   dim = this%dimensions
   p = np_old
   select case(dim)
   case(2)
    do k = 1, 1
     do j = 1, loc_y
      do i = loc_x0, loc_x
       p = p + 1
       this%x(p) = x_arr(i)
       this%y(p) = y_arr(j)
       call gasdev(u)
       this%px(p) = t_x*u
       call gasdev(u)
       this%py(p) = t_x*u
       this%weight(p) = weightx_arr(i)*weightyz_arr(j, k)
      end do
     end do
    end do
   case(3)
    do k = 1, loc_z
     do j = 1, loc_y
      do i = loc_x0, loc_x
       p = p + 1
       this%x(p) = x_arr(i)
       this%y(p) = y_arr(j)
       this%z(p) = z_arr(k)
       call gasdev(u)
       this%px(p) = t_x*u
       call gasdev(u)
       this%py(p) = t_x*u
       call gasdev(u)
       this%pz(p) = t_x*u
       whz = weightx_arr(i)*weightyz_arr(j, k)
       this%weight(p) = whz
      end do
     end do
    end do
   end select
  end subroutine

  subroutine call_particle_single( this, particles, index_in)
   class(base_species_T), intent(in) :: this
   real(dp), dimension(:), intent(inout) :: particles
   integer, intent(in) :: index_in

   select case(this%dimensions)
   case(1)
    particles(1) = this%x(index_in)
    particles(2) = this%px(index_in)
    particles(3) = this%weight(index_in)
    particles(4) = this%gamma_inv(index_in)
    particles(5) = this%part_index(index_in)
   case(2)
    particles(1) = this%x(index_in)
    particles(2) = this%y(index_in)
    particles(3) = this%px(index_in)
    particles(4) = this%py(index_in)
    particles(5) = this%weight(index_in)
    particles(6) = this%gamma_inv(index_in)
    particles(7) = this%part_index(index_in)
   case(3)
    particles(1) = this%x(index_in)
    particles(2) = this%y(index_in)
    particles(3) = this%z(index_in)
    particles(4) = this%px(index_in)
    particles(5) = this%py(index_in)
    particles(6) = this%pz(index_in)
    particles(7) = this%weight(index_in)
    particles(8) = this%gamma_inv(index_in)
    particles(9) = this%part_index(index_in)
   end select
  end subroutine

  subroutine call_particle_bounds( this, particles, lb, ub)
   class(base_species_T), intent(in) :: this
   real(dp), dimension(:, :), intent(inout) :: particles
   integer, intent(in) :: lb, ub

   select case(this%dimensions)
   case(1)
    particles(1:(ub-lb+1), 1) = this%x(lb:ub)
    particles(1:(ub-lb+1), 2) = this%px(lb:ub)
    particles(1:(ub-lb+1), 3) = this%weight(lb:ub)
    particles(1:(ub-lb+1), 4) = this%gamma_inv(lb:ub)
    particles(1:(ub-lb+1), 5) = this%part_index(lb:ub)
   case(2)
    particles(1:(ub-lb+1), 1) = this%x(lb:ub)
    particles(1:(ub-lb+1), 2) = this%y(lb:ub)
    particles(1:(ub-lb+1), 3) = this%px(lb:ub)
    particles(1:(ub-lb+1), 4) = this%py(lb:ub)
    particles(1:(ub-lb+1), 5) = this%weight(lb:ub)
    particles(1:(ub-lb+1), 6) = this%gamma_inv(lb:ub)
    particles(1:(ub-lb+1), 7) = this%part_index(lb:ub)
   case(3)
    particles(1:(ub-lb+1), 1) = this%x(lb:ub)
    particles(1:(ub-lb+1), 2) = this%y(lb:ub)
    particles(1:(ub-lb+1), 3) = this%z(lb:ub)
    particles(1:(ub-lb+1), 4) = this%px(lb:ub)
    particles(1:(ub-lb+1), 5) = this%py(lb:ub)
    particles(1:(ub-lb+1), 6) = this%pz(lb:ub)
    particles(1:(ub-lb+1), 7) = this%weight(lb:ub)
    particles(1:(ub-lb+1), 8) = this%gamma_inv(lb:ub)
    particles(1:(ub-lb+1), 9) = this%part_index(lb:ub)
   end select
  end subroutine

  subroutine call_particle_index_array( this, particles, index_in)
   class(base_species_T), intent(in) :: this
   real(dp), dimension(:, :), intent(inout) :: particles
   integer, dimension(:), intent(in) :: index_in
   integer :: n, size_ind, idx, k

   size_ind = SIZE(index_in, DIM=1)
   k = 1
   select case(this%dimensions)
   case(1)
    do n = 1, size_ind
     idx = index_in(n)
     particles(k, 1) = this%x(idx)
     particles(k, 2) = this%px(idx)
     particles(k, 3) = this%weight(idx)
     particles(k, 4) = this%gamma_inv(idx)
     particles(k, 5) = this%part_index(idx)
     k = k + 1
    end do
   case(2)
    do n = 1, size_ind
     idx = index_in(n)
     particles(k, 1) = this%x(idx)
     particles(k, 2) = this%y(idx)
     particles(k, 3) = this%px(idx)
     particles(k, 4) = this%py(idx)
     particles(k, 5) = this%weight(idx)
     particles(k, 6) = this%gamma_inv(idx)
     particles(k, 7) = this%part_index(idx)
     k = k + 1 
    end do
   case(3)
    do n = 1, size_ind
     idx = index_in(n)
     particles(k, 1) = this%x(idx)
     particles(k, 2) = this%y(idx)
     particles(k, 3) = this%z(idx)
     particles(k, 4) = this%px(idx)
     particles(k, 5) = this%py(idx)
     particles(k, 6) = this%pz(idx)
     particles(k, 7) = this%weight(idx)
     particles(k, 8) = this%gamma_inv(idx)
     particles(k, 9) = this%part_index(idx)
     k = k + 1
    end do
   end select
  end subroutine

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

  pure function array_size( this ) result( number )
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

  subroutine copy_all( this, other )
   class(base_species_T), intent(inout) :: this
   class(base_species_T), intent(in) :: other
   integer :: tot

   tot = other%n_part

   call this%reallocate(tot, other%dimensions)
   call this%set_charge(other%charge)

   if (other%allocated_x) then
    call assign(this%x, other%x(1:tot), 1, tot)
   end if
   if (other%allocated_y) then
    call assign(this%y, other%y(1:tot), 1, tot)
   end if
   if (other%allocated_z) then
    call assign(this%z, other%z(1:tot), 1, tot)
   end if
   if (other%allocated_px) then
    call assign(this%px, other%px(1:tot), 1, tot)
   end if
   if (other%allocated_py) then
    call assign(this%py, other%py(1:tot), 1, tot)
   end if
   if (other%allocated_pz) then
    call assign(this%pz, other%pz(1:tot), 1, tot)
   end if
   if (other%allocated_gamma) then
    call assign(this%gamma_inv, other%gamma_inv(1:tot), 1, tot)
   end if
   if (other%allocated_weight) then
    call assign(this%weight, other%weight(1:tot), 1, tot)
   end if
   if (other%allocated_index) then
    call assign(this%part_index, other%part_index(1:tot), 1, tot)
   end if
  end subroutine

  subroutine copy_boundaries( this, other, lower_bound, upper_bound )
   class(base_species_T), intent(inout) :: this
   class(base_species_T), intent(in) :: other
   integer, intent(in) :: lower_bound, upper_bound
   integer :: tot

   tot = upper_bound - lower_bound + 1

   call this%reallocate(tot, other%dimensions)
   call this%set_charge(other%charge)

   if (other%allocated_x) then
    call assign(this%x, other%x(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_y) then
    call assign(this%y, other%y(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_z) then
    call assign(this%z, other%z(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_px) then
    call assign(this%px, other%px(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_py) then
    call assign(this%py, other%py(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_pz) then
    call assign(this%pz, other%pz(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_gamma) then
    call assign(this%gamma_inv, other%gamma_inv(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_weight) then
    call assign(this%weight, other%weight(lower_bound:upper_bound), 1, tot)
   end if
   if (other%allocated_index) then
    call assign(this%part_index, other%part_index(lower_bound:upper_bound), 1, tot)
   end if
  end subroutine
  
  pure function how_many( this ) result(n_parts)
   !! Number of particles in the species
   class(base_species_T), intent(in) :: this
   integer :: n_parts

   n_parts = this%n_part
  
  end function

   pure function flatten( this ) result(flat_array)
   class(base_species_T), intent(in) :: this
   integer :: array_size, num_comps, i
   real(dp), allocatable :: temp(:, :), flat_array(:)

   array_size = this%how_many()
   num_comps = this%total_size()
   allocate(temp( array_size, num_comps ))

   i = 1
   if( this%allocated_x ) then
    temp( :, i ) = this%x(1:array_size)
    i = i + 1
   end if
   if( this%allocated_y ) then
    temp( :, i ) = this%y(1:array_size)
    i = i + 1
   end if
   if( this%allocated_z ) then
    temp( :, i ) = this%z(1:array_size)
    i = i + 1
   end if
   if( this%allocated_px ) then
    temp( :, i ) = this%px(1:array_size)
    i = i + 1
   end if
   if( this%allocated_py ) then
    temp( :, i ) = this%py(1:array_size)
    i = i + 1
   end if
   if( this%allocated_pz ) then
    temp( :, i ) = this%pz(1:array_size)
    i = i + 1
   end if
   if( this%allocated_gamma ) then
    temp( :, i ) = this%gamma_inv(1:array_size)
    i = i + 1
   end if
   if( this%allocated_weight ) then
    temp( :, i ) = this%weight(1:array_size)
    i = i + 1
   end if
   if( this%allocated_index ) then
    temp( :, i ) = this%part_index(1:array_size)
    i = i + 1
   end if

   flat_array = PACK( temp(:, :), .true. )

  end function

  subroutine initialize_data( this, x_arr, y_arr, z_arr, &
   weightx_arr, weightyz_arr, loc_x, loc_y, loc_z)
   class( base_species_T ), intent(inout) :: this
   real(dp), dimension(:), intent(in) :: x_arr, y_arr, z_arr
   real(dp), dimension(:), intent(in) :: weightx_arr
   real(dp), dimension(:, :), intent(in) :: weightyz_arr
   integer, intent(in) :: loc_x, loc_y, loc_z

   call this%add_data( x_arr, y_arr, z_arr, &
   weightx_arr, weightyz_arr, 1, loc_x, loc_y, loc_z, 0)

  end subroutine

  subroutine pack_into_logical( this, packed, mask )
   class(base_species_T), intent(in) :: this
   class(base_species_T), intent(inout) :: packed
   logical, intent(in) :: mask
   integer :: np
 
   np = this%how_many()
   call packed%sweep()
   call packed%new_species(np, this%dimensions)
   call packed%set_charge(this%charge)
 
   if( this%allocated_x ) then
    packed%x = PACK( this%x(1:np), mask)
   end if
   if( this%allocated_y ) then
    packed%y = PACK( this%y(1:np), mask)
   end if
   if( this%allocated_z ) then
    packed%z = PACK( this%z(1:np), mask)
   end if
   if( this%allocated_px ) then
    packed%px = PACK( this%px(1:np), mask)
   end if
   if( this%allocated_py ) then
    packed%py = PACK( this%py(1:np), mask)
   end if
   if( this%allocated_pz ) then
    packed%pz = PACK( this%pz(1:np), mask)
   end if
   if( this%allocated_gamma ) then
    packed%gamma_inv = PACK( this%gamma_inv(1:np), mask)
   end if
   if( this%allocated_weight ) then
    packed%weight = PACK( this%weight(1:np), mask)
   end if
   if( this%allocated_index ) then
    packed%part_index = PACK( this%part_index(1:np), mask)
   end if
 
   packed%n_part = packed%array_size()
 
  end subroutine
 
  subroutine pack_into_array( this, packed, mask )
   class(base_species_T), intent(in) :: this
   class(base_species_T), intent(inout) :: packed
   logical, intent(in) :: mask(:)
   integer :: tot_parts, np
 
   tot_parts = COUNT( mask )
   np = this%how_many()
   call packed%sweep()
   call packed%new_species(tot_parts, this%dimensions)
   call packed%set_charge(this%charge)
 
   if (tot_parts == 0) then
    return
   end if
   if( this%allocated_x ) then
    packed%x = PACK( this%x(1:np), mask(:) )
   end if
   if( this%allocated_y ) then
    packed%y = PACK( this%y(1:np), mask(:) )
   end if
   if( this%allocated_z ) then
    packed%z = PACK( this%z(1:np), mask(:) )
   end if
   if( this%allocated_px ) then
    packed%px = PACK( this%px(1:np), mask(:) )
   end if
   if( this%allocated_py ) then
    packed%py = PACK( this%py(1:np), mask(:) )
   end if
   if( this%allocated_pz ) then
    packed%pz = PACK( this%pz(1:np), mask(:) )
   end if
   if( this%allocated_gamma ) then
    packed%gamma_inv = PACK( this%gamma_inv(1:np), mask(:) )
   end if
   if( this%allocated_weight ) then
    packed%weight = PACK( this%weight(1:np), mask(:) )
   end if
   if( this%allocated_index ) then
    packed%part_index = PACK( this%part_index(1:np), mask(:) )
   end if
 
  end subroutine

  subroutine redistribute( this, flat_array, num_particles, dimensions )
   class(base_species_T), intent(inout) :: this
   real(dp), intent(in), dimension(:) :: flat_array
   integer, intent(in) :: num_particles, dimensions
   integer :: i

   i = 0
   call this%reallocate(num_particles, dimensions)
   if( this%allocated_x ) then
    this%x(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_y ) then
    this%y(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_z ) then
    this%z(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_px ) then
    this%px(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_py ) then
    this%py(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_pz ) then
    this%pz(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_gamma ) then
    this%gamma_inv(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_weight ) then
    this%weight(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if
   if( this%allocated_index ) then
    this%part_index(1:num_particles) = flat_array((i + 1): (i + num_particles))
    i = i + num_particles
   end if

  end subroutine

  subroutine reallocate(this, n_parts, n_dimensions)
   class(base_species_T), intent(inout) :: this
   integer :: n_parts, n_dimensions, sp_charge

   if ( allocated(this%initialized) ) then
    if (this%n_part < n_parts ) then
     sp_charge = this%charge
     call this%sweep()
     call this%new_species(n_parts, n_dimensions)
     call this%set_charge(sp_charge)
    end if
   else
    call this%new_species(n_parts, n_dimensions)
   end if

  end subroutine

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
  if ( n_parts > 0 ) then
   this%empty = .false.
  else if (n_parts == 0) then
   this%empty = .true. 
  else
   write(6, *) 'Error in part number'
  end if
 end subroutine

 subroutine set_temperature( this, temperature)
  class(base_species_T), intent(inout) :: this
  real(dp), intent(in) :: temperature

  this%temperature = temperature
 end subroutine

 pure function total_size( this ) result(size)
  class(base_species_T), intent(in) :: this
  integer :: size

  select case(this%dimensions)
  case(1)
   size = 5
  case(2)
   size = 7
  case(3)
   size = 9
  case default
   size = -1
  end select

 end function

!==== Procedures not bound to type ======

 subroutine assign_real_realdp( array, values, lb, ub, n_parts)
  real(dp), allocatable, dimension(:), intent(inout) :: array
  real(dp), dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = values(:)

 end subroutine

 subroutine assign_integer_realdp( array, values, lb, ub, n_parts)
  real(dp), allocatable, dimension(:), intent(inout) :: array
  integer, dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = real(values(:), dp)

 end subroutine

 subroutine assign_real_realsp( array, values, lb, ub, n_parts)
  real(sp), allocatable, dimension(:), intent(inout) :: array
  real(dp), dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = real(values(:), sp)

 end subroutine

 subroutine assign_realsp_realsp( array, values, lb, ub, n_parts)
  real(sp), allocatable, dimension(:), intent(inout) :: array
  real(sp), dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = values(:)

 end subroutine

 subroutine assign_integer_realsp( array, values, lb, ub, n_parts)
  real(sp), allocatable, dimension(:), intent(inout) :: array
  integer, dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = real(values(:), sp)

 end subroutine

 subroutine assign_real_integer( array, values, lb, ub, n_parts)
  integer, allocatable, dimension(:), intent(inout) :: array
  real(dp), dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = int(values(:))

 end subroutine

 subroutine assign_integer_integer( array, values, lb, ub, n_parts)
  integer, allocatable, dimension(:), intent(inout) :: array
  integer, dimension(:), intent(in)  :: values
  integer, intent(in)  :: lb, ub
  integer, intent(in), optional  :: n_parts
  integer :: size_value

  size_value = SIZE(values, DIM=1)
  if ( (ub - lb + 1) > size_value ) then
   write( 6, *) 'Assigning wrong value size'
  end if

  if ( .not. allocated(array) ) then
   allocate( array(n_parts) )
  end if

  array(lb:ub) = values(:)

 end subroutine

 subroutine check_array_1d( array, size_array )
  real(dp), allocatable, dimension(:), intent(inout) :: array
  integer, intent(in) :: size_array

  if ( allocated(array) ) then
   if ( SIZE(array) < size_array ) then
    deallocate(array)
    allocate(array(size_array))
   endif
  else
   allocate(array(size_array))
  end if
  
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
  case default
   cdir = -10
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
  case default
   pm_comp = -10
  end select
 end function

end module