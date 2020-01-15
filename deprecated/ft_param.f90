!*****************************************************************************************************!
!             Copyright 2008-2016 Pasquale Londrillo, Stefano Sinigardi, Andrea Sgattoni              !
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


  integer fftw_forward, fftw_backward
  parameter (fftw_forward=-1, fftw_backward=1)

  integer fftw_real_to_complex, fftw_complex_to_real
  parameter (fftw_real_to_complex=-1, fftw_complex_to_real=1)

  integer fftw_estimate, fftw_measure
  parameter (fftw_estimate=0, fftw_measure=1)

  integer fftw_out_of_place, fftw_in_place, fftw_use_wisdom
  parameter (fftw_out_of_place=0)
  parameter (fftw_in_place=8, fftw_use_wisdom=16)

  integer fftw_threadsafe
  parameter (fftw_threadsafe=128)

! Constants for the MPI wrappers:
  integer fftw_transposed_order, fftw_normal_order
  integer fftw_scrambled_input, fftw_scrambled_output
  parameter (fftw_transposed_order=1, fftw_normal_order=0)
  parameter (fftw_scrambled_input=8192)
  parameter (fftw_scrambled_output=16384)
  integer fftw_redft00

