!*****************************************************************************************************!
!                            Copyright 2008-2019  The ALaDyn Collaboration                            !
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

 module array_ralloc

  use pstruct_data
  use fstruct_data
  implicit none

  !--------------------------
  !DIR$ ATTRIBUTES INLINE :: p_realloc
 contains
  subroutine p_realloc(pdata, npt_new, ndv)
   type (species), intent (inout) :: pdata
   integer, intent (in) :: npt_new, ndv
   integer :: allocstatus, deallocstatus

   if (allocated(pdata%part)) then
    if (size(pdata%part,1)<npt_new) then
     deallocate (pdata%part, stat=deallocstatus)
     if (deallocstatus==0) allocate (pdata%part(1:npt_new,1:ndv), &
       stat=allocstatus)
    end if
   else
    allocate (pdata%part(1:npt_new,1:ndv), stat=allocstatus)
   end if
   pdata%part(:, :) = 0.0
  end subroutine
  !========================
  subroutine v_realloc(vdata, npt_new, ndv)
   real (dp), allocatable, intent (inout) :: vdata(:, :)
   integer, intent (in) :: npt_new, ndv
   integer :: allocstatus, deallocstatus

   if (allocated(vdata)) then
    if (size(vdata,1)<npt_new) then
     deallocate (vdata, stat=deallocstatus)
     allocate (vdata(1:npt_new,1:ndv), stat=allocstatus)
    end if
   else
    allocate (vdata(1:npt_new,1:ndv), stat=allocstatus)
   end if
   vdata(:, :) = 0.0
  end subroutine
  !===========================
 end module

