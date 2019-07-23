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

 module rk_param

  use precision_def

  implicit none
  real (dp) :: a_rk(6), b_rk(6), c_rk(0:6)
  real (dp) :: a_lpf(4), b_lpf(4), c_lpf(4)
  integer :: iform, rk_form(5), lpf_form(4)

 contains

  subroutine lpf_sch_coeff

   a_lpf(1:4) = 1.
   b_lpf(1:4) = 1.
   lpf_form(1:4) = iform
   lpf_form(4) = 0
   b_lpf(1) = 1./(2.-2**(1./3.)) !alpha
   b_lpf(2) = 1. - 2.*b_lpf(1) !beta=1.-2*alpha
   b_lpf(3) = b_lpf(1) !alpha
   b_lpf(4) = 0.0
   a_lpf(1) = 0.5*b_lpf(1) !alpha/2
   a_lpf(2) = 0.5*(1.-b_lpf(1)) !1/2(1-alpha)
   a_lpf(3) = a_lpf(2)
   a_lpf(4) = a_lpf(1)
   c_lpf(1) = 0.0
   c_lpf(2) = a_lpf(1)
   c_lpf(3) = 0.5*b_lpf(2)
   c_lpf(4) = c_lpf(2)

  end subroutine

!--------------------------

  subroutine rk_tsch_coeff(rk)
   integer, intent (in) :: rk
   integer :: ip

   select case (rk)
   case (3)
    b_rk(1) = 1.
    b_rk(2) = 0.25
    b_rk(3) = 2./3.
    a_rk(1) = 0.0
    a_rk(2) = 1. - b_rk(2)
    a_rk(3) = 1. - b_rk(3)

   case (4)
    rk_form(1:4) = 2
    a_rk(1) = 1./6.
    a_rk(2) = 1./3.
    a_rk(3) = 1./3.
    a_rk(4) = a_rk(1)
    c_rk(0) = -4./3.
    c_rk(1) = 1./3.
    c_rk(2) = 2./3.
    c_rk(3) = c_rk(1)
    c_rk(4) = 1.
    b_rk(1) = 0.5
    b_rk(2) = 0.5
    b_rk(3) = 1.0
    b_rk(4) = 1.0/6.0

   case (5) !4th order five-stage optimized
    rk_form(1:5) = iform
    a_rk(1) = 0.4
    a_rk(2) = 0.0581597561
    a_rk(3) = 1.0746983979
    a_rk(4) = -0.5338616876
    a_rk(5) = 0.0

    b_rk(1) = -0.01
    b_rk(2) = 0.4818402439
    b_rk(3) = -0.8615395999
    b_rk(4) = 1.057293281
    b_rk(5) = 0.3324060746

    c_rk(1) = 0.0
    c_rk(2) = 0.39
    c_rk(3) = 0.53
    c_rk(4) = 0.6349990419
    c_rk(5) = 0.1337322374

    do ip = 2, 5
     c_rk(ip) = b_rk(ip-1) + c_rk(ip-1)
    end do

    do ip = 2, 5
     c_rk(ip) = a_rk(ip-1) + c_rk(ip)
    end do

   case (6) !4th-order six-stage optimized
    b_rk(1) = 0.10893125722541
    b_rk(2) = 0.13201701492152
    b_rk(3) = 0.38911623225517
    b_rk(4) = -0.59203884581148
    b_rk(5) = 0.47385028714844
    b_rk(6) = 0.48812405426094
    a_rk(1) = 0.17985400977138
    a_rk(2) = 0.14081893152111
    a_rk(3) = 0.08255631629428
    a_rk(4) = 0.65804425034331
    a_rk(5) = 0.31862993413251
    a_rk(6) = 0.
    c_rk(1) = 0.0
    do ip = 2, 6
     c_rk(ip) = b_rk(ip-1) + c_rk(ip-1)
    end do
    do ip = 2, 6
     c_rk(ip) = a_rk(ip-1) + c_rk(ip)
    end do

   end select

  end subroutine

 end module
