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
!#define NEW_SP

 module struct_def
  use precision_def
  implicit none
  public

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
  end type

  interface species_new
   module procedure :: new_species_new
  end interface

  type grid
   integer :: ng
   !!Number of cells in a given direction of the grid
   integer :: p_ind(2)
   !!Minimum and maximum cell number of the grid
   real (dp) :: gmin
   !!Value of the corresponding axis at the minimum cell
   real (dp) :: gmax
   !!Value of the corresponding axis at the maximum cell
   integer :: min_cell
   !!Initial cell of the grid in absolute units (i.e. respect to the total grid)
   integer :: max_cell
   !!Final cell of the grid in absolute units (i.e. respect to the total grid)
  end type

  type sgrid
   integer :: sind(2)
   !!Initial and final stretched cell (sind(1) also coincides with the number of
   !!stretched cells)
   real (dp) :: smin
   !!Axis value on the boundary between stretched and unstretched grid (left side of the box)
   real (dp) :: smax
   !!Axis value on the boundary between stretched and unstretched grid (right side of the box)
  end type

  type index_array
   !! Type defining an array of consecutive integer numbers, useful as
   !! indices in arrays.
   integer, allocatable :: indices(:)
   contains
    procedure, public :: find_index
  end type

  interface index_array
   module procedure new_index_array
  end interface

  contains

   function new_index_array(length) result(this)
    !! Constructor for the index_array type
    integer, intent(in) :: length
    type(index_array) :: this
    integer :: i

    allocate( this%indices(length) )
    this%indices = [ (i, i = 1, length) ]
   end function

   subroutine find_index( index_in, mask )
    !! Type bound procedure that finds and pack all the array indices
    !! according to the given mask

    class(index_array), intent(inout) :: index_in
    logical, intent(in) :: mask(:)
    index_in%indices = PACK( index_in%indices, mask )
   end subroutine

   pure function how_many( this ) result(n_parts)
    class(species_new), intent(in) :: this
    integer :: n_parts
   
    n_parts = this%n_part
   
   end function
  
   function new_species_new( npt_max, curr_ndims ) result(this)
    integer, intent(in) :: npt_max, curr_ndims
    type(species_new) :: this
    integer :: allocstatus
   
    if (npt_max <= 0) then
     this%initialized = .false.
     this%n_part = 0
     return
    end if
   
    this%initialized = .true.
    this%n_part = npt_max
   
    select case(curr_ndims)
    
    case(1)
    
     allocate( this%x(npt_max), stat=allocstatus)
     allocate( this%px(npt_max), stat=allocstatus)
     allocate( this%weight(npt_max), stat=allocstatus)
     allocate( this%part_index(npt_max), stat=allocstatus)
    
    case(2)
    
     allocate( this%x(npt_max), stat=allocstatus)
     allocate( this%px(npt_max), stat=allocstatus)
     allocate( this%y(npt_max), stat=allocstatus)
     allocate( this%py(npt_max), stat=allocstatus)
     allocate( this%weight(npt_max), stat=allocstatus)
     allocate( this%part_index(npt_max), stat=allocstatus)
    
    case(3)
    
     allocate( this%x(npt_max), stat=allocstatus)
     allocate( this%px(npt_max), stat=allocstatus)
     allocate( this%y(npt_max), stat=allocstatus)
     allocate( this%py(npt_max), stat=allocstatus)
     allocate( this%z(npt_max), stat=allocstatus)
     allocate( this%pz(npt_max), stat=allocstatus)
     allocate( this%weight(npt_max), stat=allocstatus)
     allocate( this%part_index(npt_max), stat=allocstatus)
    end select
   end function
 end module