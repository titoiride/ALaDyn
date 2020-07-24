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

 module read_input

  use read_namelist
#if defined(JSON_FORTRAN)
  use read_json
#endif

  implicit none

 contains


  subroutine read_main_input
   logical :: exist_nml = .false.
   logical :: exist_json = .false.

   input_json_filename = 'input.json'
   input_namelist_filename = 'input.nml'
   inquire (file=input_namelist_filename, exist=exist_nml)

#if defined(JSON_FORTRAN)
   inquire (file=input_json_filename, exist=exist_json)   
   if (exist_json) then
    call read_input_json( input_json_filename, namelist)
   end if
#endif

   if (.not. exist_json .and. exist_nml) then
    write (6, *) '==========================================='
    write (6, *) 'No JSON input file has been found.&
    & Instead, the old nml format will be used.'
    write (6, *) '==========================================='
    call read_input_nml
   else
    write (6, *) 'No input file (.json or .nml) has been found'
    stop
   end if
  end subroutine

 end module
