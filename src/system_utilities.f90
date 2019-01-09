 !*****************************************************************************************************!
 !                            Copyright 2008-2018  The ALaDyn Collaboration                            !
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

 module system_utilities
 use mpi_var
 implicit none
 contains
 !--------------------------

 subroutine create_timestep_folder(iout)
 integer, intent(in) :: iout
 character(5) :: folderName

 write (folderName,'(i4.4)') iout
 folderName=TRIM(ADJUSTL(folderName))//char(0)
 if(pe0) call create_folder(folderName)

 end subroutine create_timestep_folder


 subroutine create_initial_folders
 character(12) :: folderDumpRestart
 character(12) :: folderDiagnostics
 folderDumpRestart="dumpRestart"//char(0)
 folderDiagnostics="diagnostics"//char(0)
 if(pe0) call create_folder(folderDumpRestart)
 if(pe0) call create_folder(folderDiagnostics)
 end subroutine create_initial_folders



 !---------------------------
 end module system_utilities
