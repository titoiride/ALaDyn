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

! ===== File kept to store old unused functions (but that may be useful in some future) =====

  subroutine pwfa_density(wp, kp, i2, j2, k2, dir)
   !! WARNING: To be checked (Not used in the code)
   real (dp), intent (inout) :: wp(:, :, :)
   real (dp), intent (in) :: kp
   integer, intent (in) :: i2, j2, k2, dir
   integer :: ii, ik, i, j, k
   real (dp) :: kpx, sum0(2), temp

   ! in wp enters the beam charge density at staggered x-coordinate
   select case (dir)
   case (1)
    allocate (kern(i2))
    kpx = kp*dx
    kern = 0.0
    do i = 1, i2
     kern(i) = kpx*sin(kp*x(i))
    end do
    call ft_kern(kern, i2, -1)
    call ftw1d_sc(wp, i2, j2, k2, -1, 1, 1) !ft_sin along the x direction
    do k = 1, k2
     do j = 1, j2
      do i = 2, i2
       wp(i, j, k) = wp(i, j, k)*kern(i)
      end do
     end do
    end do
    call ftw1d_sc(wp, i2, j2, k2, 1, 1, 1)
    !=============== use n(x)=sin(kp*x)*sum_y<x[cos(kp*y)*nb(y)]
    !                         +cos(kp*x)*sum_y>x[sin(kp*y)*nb(y)]
   ! case (2)
   !  kpx = kp*dx
   !  temp = kp*sin(0.5*kpx)/(0.5*kpx)
   !  allocate (kern2(i2,2))
   !  do i = 1, i2
   !   kern2(i, 1) = sin(kp*x(i))
   !   kern2(i, 2) = cos(kp*x(i))
   !  end do
   !  do k = 1, k2
   !   do j = 1, j2
   !    w1_re(1:i2) = temp*wp(1:i2, j, k)
   !    wp(1:i2, j, k) = 0.0
   !    do i = 1, i2
   !     sum0 = 0.0
   !     ik = max(i, i2/2)
   !     do ii = ik, i2
   !      sum0(1) = sum0(1) + kern2(ii, 1)*w1_re(ii)
   !      sum0(2) = sum0(2) + kern2(ii, 2)*w1_re(ii)
   !     end do
   !     wp(i, j, k) = sum0(1)*kern2(i, 2) - sum0(2)*kern2(i, 1)
   !    end do
   !   end do
   !  end do
    end select
  end subroutine

  !================
  subroutine ft_kern(w, n1, is)
   !! WARNING: Still to be checked
   real (dp), intent (inout) :: w(:)
   integer, intent (in) :: n1, is
   integer :: ix, i2
   real (dp) :: sc
   integer :: ndb

   !!$PRAGMA C( DFFTW_EXECUTE )

   ndb = 2*n1
   sc = 1.
   if (is<0) then
    w1_re(1:n1) = w(1:n1)
    w1_re(n1+1) = 0.0
    do ix = n1 + 2, ndb
     w1_re(ix) = -w1_re(ndb+2-ix)
    end do
    call fftw_execute_dft_r2c(plan1, w1_re, w1_cplx)
    do ix = 1, n1
     i2 = 2*ix
     w(ix) = sc*w1_cplx(i2)
    end do
   else
    w1_re(1:n1) = w(1:n1)
    do ix = n1 + 2, ndb
     w1_re(ix) = w1_re(ndb+2-ix)
    end do
    w1_re(n1+1) = 0.5*(w1_re(n1)+w1_re(n1+2))
    call fftw_execute_dft_r2c(plan1, w1_re, w1_cplx)
    do ix = 1, n1
     i2 = 2*ix - 1
     w(ix) = sc*w1_cplx(i2)
    end do
   end if
  end subroutine