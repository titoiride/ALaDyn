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

  subroutine read_bunch_namelist
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !C
   !C Reads bunch namelist for PWFA
   !C
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   !--- *** namelist *** ---!
   namelist /number_bunches/n_bunches, l_particles, &
     l_intdiagnostics_pwfa, l_intdiagnostics_classic, &
     l_embunchevolution, number_of_slices
   namelist /bunches/nb_tot, bunch_type, bunch_shape, particle_charge, &
     rhob, xc_bunch, yc_bunch, zc_bunch, gam, sxb, syb, epsy, epsz, dg, &
     charge_right, charge_left, sigma_cut_bunch, ppc_x_bunch, &
     ppc_y_bunch, ppc_z_bunch


   !--- reading number of bunches ---!
   open (nml_iounit, file=input_namelist_filename, status='old')
   l_particles = .false.
   l_intdiagnostics_pwfa = .false.
   l_intdiagnostics_classic = .true.
   l_embunchevolution = .true.
   number_of_slices = [ 10, 0, 0, 0 ]
   read (nml_iounit, number_bunches, iostat=nml_ierr)
   nml_error_message = 'NUMBER_BUNCHES'
   close (nml_iounit)
   if (nml_ierr>0) call print_at_screen_nml_error


   !--> initialization
   yc_bunch = 0.0
   zc_bunch = 0.0
   bunch_type = 1 !electron bunch
   bunch_shape = 1 !shape 1: bi-giassian
   !shape 2: trapezoidal (linear in Z, uniform with cutoff in R)
   !shape 3: trapezoidal-gaussian (linear in Z, gaussian in R)
   !shape 4: cylinder
   particle_charge = -1.
   rhob = 1.0 !relative density n_bunch/n_plasmabackground
   charge_right = -1.0
   charge_left = -1.0
   ppc_x_bunch = -1 !number of particle per cell 'x' direction :: this implies weighted option
   ppc_y_bunch = -1 !number of particle per cell 'x' direction :: this implies weighted option
   ppc_z_bunch = -1 !number of particle per cell 'x' direction :: this implies weighted option
   nb_tot = -1 !total number of bunch particles :: implies particle with same weight
   sigma_cut_bunch = 3. !standard cut at 3-rms

   open (nml_iounit, file=input_namelist_filename, status='old')
   read (nml_iounit, bunches, iostat=nml_ierr)
   nml_error_message = 'BUNCHES'
   close (nml_iounit)
   if (nml_ierr>0) call print_at_screen_nml_error

   call select_number_of_bunch_particles()
  end subroutine

  subroutine read_bunch_twiss_namelist
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !C
   !C Reads TWISS parameter for PWFA
   !C
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   namelist /twiss/l_twiss, alpha_twiss, beta_twiss

   open (nml_iounit, file=input_namelist_filename, status='old')
   l_twiss = .false. !bunch at waist
   read (nml_iounit, twiss, iostat=nml_ierr)
   nml_error_message = 'TWISS'
   close (nml_iounit)
   if (nml_ierr>0) call print_at_screen_nml_error
  end subroutine


  subroutine read_poloidal_field_namelist
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
   !C
   !C Reads Poloidal Field parameters for PWFA
   !C
   !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

   namelist /bpoloidal/l_bpoloidal, b_ex_poloidal, radius_poloidal

   open (nml_iounit, file=input_namelist_filename, status='old')
   l_bpoloidal = .false.
   b_ex_poloidal = 0.0
   radius_poloidal = 1.0
   read (nml_iounit, bpoloidal, iostat=nml_ierr)
   nml_error_message = 'BPOLOIDAL'
   close (nml_iounit)
   if (nml_ierr>0) call print_at_screen_nml_error
  end subroutine