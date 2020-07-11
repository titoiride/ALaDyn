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

 module memory_pool
  !! Module that serves as container for all auxiliary arrays to avoid
  !! unnecessary allocations and multiple copies.

  use precision_def
  use base_species
  use particles_def
  use particles_aux_def

  implicit none
  public
  
  interface array_realloc_1d
   module procedure array_realloc_1d_dp
   module procedure array_realloc_1d_sp
   module procedure array_realloc_1d_logical
  end interface

  type memory_pool_t

  integer :: iter_count = 0
  integer :: used_memory = 0
  real(dp), allocatable, dimension(:) :: mp_xx_1d
  real(dp), allocatable, dimension(:, :) :: mp_xx_2d
  real(dp), allocatable, dimension(:, :, :) :: mp_xx_3d
  type(species_new), allocatable :: mp_species_new

  contains

  final :: final_clean

  end type

  contains

  subroutine final_clean( this )
   type(memory_pool_t), intent(inout) :: this

   write(6, *) 'Called memory pool destructor'
   if ( ALLOCATED(this%mp_xx_1d) ) then
    deallocate( this%mp_xx_1d )
   end if
   if ( ALLOCATED(this%mp_xx_2d) ) then
    deallocate( this%mp_xx_2d )
   end if
   if ( ALLOCATED(this%mp_xx_3d) ) then
    deallocate( this%mp_xx_3d )
   end if
   if ( ALLOCATED(this%mp_species_new) ) then
    deallocate( this%mp_species_new )
   end if

  end subroutine

  subroutine create_memory_pool( this )
   type(memory_pool_t), pointer, intent(inout) :: this

   if ( .not. ASSOCIATED(this) ) then
    write(6, *) 'Constructing the memory pool'
    ALLOCATE(this)
   end if

  end subroutine

  subroutine destroy_memory_pool( this )
   type(memory_pool_t), pointer, intent(inout) :: this

   if ( ASSOCIATED(this) ) then
    deallocate(this)
   end if

  end subroutine

  subroutine allocate_1d( this )
   type(memory_pool_t), pointer, intent(in) :: this

   allocate(this%mp_xx_1d(-5:5), source = 10.)

  end subroutine

  subroutine point_to_mp( this )
   type(memory_pool_t), pointer, intent(in) :: this
   real(dp), pointer, dimension(:) :: point_1d

   point_1d => this%mp_xx_1d

   write(6, *) "I'm pointing to ", point_1d(:)
   write(6, *) "My size", SIZE( point_1d, DIM=1)
   write(6, *) "Lower boundary", LBOUND( point_1d, DIM=1)
   write(6, *) "Upper boundary", UBOUND( point_1d, DIM=1)
   write(6, *) "I'm associated", ASSOCIATED(point_1d)
   write(6, *) "Now I'm ", point_1d(:)
   write(6, *) "My address", LOC(point_1d)
   write(6, *) "My address", LOC(this%mp_xx_1d)
   
  end subroutine
  !==============================================
  ! Subroutines for allocation
  !==============================================
  !DIR$ ATTRIBUTES INLINE :: array_realloc_1d_dp
  subroutine array_realloc_1d_dp(vdata, npt_new)
   real (dp), allocatable, intent (inout) :: vdata(:)
   integer, intent (in) :: npt_new
   integer :: allocstatus, deallocstatus

   if (allocated(vdata)) then
    if (size(vdata,1)<npt_new) then
     deallocate (vdata, stat=deallocstatus)
     allocate (vdata(1:npt_new), stat=allocstatus)
    end if
   else
    allocate (vdata(1:npt_new), stat=allocstatus)
   end if
  end subroutine
  !===========================
  !DIR$ ATTRIBUTES INLINE :: array_realloc_1d_sp
  subroutine array_realloc_1d_sp(vdata, npt_new)
   real (sp), allocatable, intent (inout) :: vdata(:)
   integer, intent (in) :: npt_new
   integer :: allocstatus, deallocstatus

   if (allocated(vdata)) then
    if (size(vdata,1)<npt_new) then
     deallocate (vdata, stat=deallocstatus)
     allocate (vdata(1:npt_new), stat=allocstatus)
    end if
   else
    allocate (vdata(1:npt_new), stat=allocstatus)
   end if
  end subroutine
  !===========================
  !DIR$ ATTRIBUTES INLINE :: array_realloc_1d_logical
  subroutine array_realloc_1d_logical(vdata, npt_new)
   logical, allocatable, intent (inout) :: vdata(:)
   integer, intent (in) :: npt_new
   integer :: allocstatus, deallocstatus

   if (allocated(vdata)) then
    if (size(vdata,1)<npt_new) then
     deallocate (vdata, stat=deallocstatus)
     allocate (vdata(1:npt_new), stat=allocstatus)
    end if
   else
    allocate (vdata(1:npt_new), stat=allocstatus)
   end if
  end subroutine
  !==========================================
  !DIR$ ATTRIBUTES INLINE :: xx_realloc
  subroutine xx_realloc( xx_in, npart_new, dimensions )
   real(dp), allocatable, dimension(:, :), intent(inout) :: xx_in
   integer, intent(in) :: npart_new, dimensions
   integer :: allocstatus, deallocstatus

   if ( ALLOCATED(xx_in) ) then   
    if (SIZE(xx_in, DIM=1) < npart_new .or. &
     SIZE(xx_in, DIM=2) < dimensions ) then
     deallocate(xx_in, stat=deallocstatus)
     allocate(xx_in(npart_new, dimensions), stat=allocstatus)
    end if
   else
    allocate(xx_in(npart_new, dimensions), stat=allocstatus)
   end if
  end subroutine
  !==========================================
 end module