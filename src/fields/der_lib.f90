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

 module der_lib

  use grid_field_param
  use util

  implicit none

  real (dp) :: mat_env_coeff(5, 5)

  real (dp), allocatable :: mat_der_inv(:, :), fmat(:, :, :), &
    lpl_mat(:, :, :)
  real (dp), allocatable :: mat_env(:, :)
  logical :: derinv
 contains
!==========================================
  subroutine env_matrix_inv(jcr, lk0, dhx, b_diag, i1, n1p, j1, n2p, k1, &
    n3p)
   real (dp), intent (inout) :: jc(:, :, :, :)
   real (dp), intent (in) :: lk0, dhx, b_diag
   integer, intent (in) :: i1, n1p, j1, n2p, k1, n3p
   integer :: ii, i, j, k, ic
   real (dp) :: dx2_inv, aph2, b_1, c_1, b_n, a_n

   dx2_inv = dhx*dhx
   aph2 = 1./dx2_inv !dx*dx
   b_1 = b_diag + 2.
   c_1 = -b_1
   b_n = b_1
   a_n = c_1
   do k = k1, n3p
    do j = j1, n2p
     ii = 1
     i = i1
     ww0(ii, 1) = dhx*(jcr(i+1,j,k,1)-jcr(i,j,k,1)) + lk0*jc(i, j, k, 2)
     ww0(ii, 2) = dhx*(jcr(i+1,j,k,2)-jcr(i,j,k,2)) - lk0*jc(i, j, k, 1)
     do i = i1 + 1, n1p - 1
      ii = ii + 1
      ww0(ii, 1) = 0.5*dhx*(jcr(i+1,j,k,1)-jcr(i-1,j,k,1)) + &
        lk0*jcr(i, j, k, 2)
      ww0(ii, 2) = 0.5*dhx*(jcr(i+1,j,k,2)-jcr(i-1,j,k,2)) - &
        lk0*jcr(i, j, k, 1)
     end do
     ii = ii + 1
     i = n1p
     ww0(ii, 1) = dhx*(jcr(i,j,k,1)-jcr(i-1,j,k,1)) + lk0*jc(i, j, k, 2)
     ww0(ii, 2) = dhx*(jcr, j, k, 2) - jcr(i-1, j, k, 2)) - lk0*jcr(i, &
       j, k, 1)
    end do
   end do
   do ic = 1, 2
    do k = k1, n3p
     do j = j1, n2p
      do i = i1, n1p
       jc(i, j, k, ic) = aph2*jcr(i, j, k, ic)
      end do
     end do
    end do
   end do
  end subroutine
!===============================
  subroutine set_filt_mat(ng1, dir)
   integer, intent (in) :: ng1, dir
   integer :: i

   fmat(1: ng1, 1, dir) = filt_coeff(0)
   fmat(1: ng1, 2, dir) = 1.0
   fmat(1: ng1, 3, dir) = filt_coeff(0)
!        LU Factorize
   fmat(1, 2, dir) = 1.0/fmat(1, 2, dir)
   do i = 2, ng1
    fmat(i - 1, 3, dir) = fmat(i - 1, 3, dir)*fmat(i - 1, 2, dir)
    fmat(i, 2, dir) = fmat(i, 2, dir) - fmat(i - 1, 3, dir)*fmat(i, 1, &
      dir)
    fmat(i, 2, dir) = 1.0/fmat(i, 2, dir)
   end do
  end subroutine
!=============================
  subroutine trid_lpl_inv(n, dir)
   integer, intent (in) :: n, dir
   integer :: ix

   ix = 1
   dw(ix) = lpl_mat(ix, 2, dir)*dw(ix)
   do ix = 2, n
    dw(ix) = dw(ix) - lpl_mat(ix, 1, dir)*dw(ix - 1)
    dw(ix) = lpl_mat(ix, 2, dir)*dw(ix)
   end do
   do ix = n - 1, 1, -1
    dw(ix) = dw(ix) - lpl_mat(ix, 3, dir)*dw(ix + 1)
   end do
  end subroutine
!=========================
  subroutine trid_der_inv(n, ib)
   integer, intent (in) :: n, ib
   integer :: ix
   real (dp) :: alp0, bet0, gm, fact

   ix = 1
   dw(ix) = mat_der_inv(ix, 2)*dw(ix)
   do ix = 2, n
    dw(ix) = dw(ix) - mat_der_inv(ix, 1)*dw(ix - 1)
    dw(ix) = mat_der_inv(ix, 2)*dw(ix)
   end do
   do ix = n - 1, 1, -1
    dw(ix) = dw(ix) - mat_der_inv(ix, 3)*dw(ix + 1)
   end do
   if (ib /= 1) return
   alp0 = mat_der_inv(n, 3)
   bet0 = mat_der_inv(1, 1)
   gm = -1.0
   zp(2: n - 1) = 0.0
   zp(1) = gm
   zp(n) = alp0
!-----------call trid of z
   zp(1) = mat_der_inv(1, 2)*zp(1)
   do ix = 2, n
    zp(ix) = zp(ix) - mat_der_inv(ix, 1)*zp(ix - 1)
    zp(ix) = mat_der_inv(ix, 2)*zp(ix)
   end do
   do ix = n - 1, 1, -1
    zp(ix) = zp(ix) - mat_der_inv(ix, 3)*zp(ix + 1)
   end do
   fact = (dw(1) + bet0*dw(n)/gm)/(1.0 + zp(1) + bet0*zp(n)/gm)
   dw(1: n) = dw(1: n) - fact*zp(1: n)
   dw(n + 1) = dw(1)
!=======================
  end subroutine
!============================
  subroutine ftrid(n, ib, dir)
   integer, intent (in) :: n, ib, dir
   integer :: ix
   real (dp) :: alp0, bet0, gm, fact

   ix = 1
   dw(ix) = fmat(ix, 2, dir)*dw(ix)
   do ix = 2, n
    dw(ix) = dw(ix) - fmat(ix, 1, dir)*dw(ix - 1)
    dw(ix) = fmat(ix, 2, dir)*dw(ix)
   end do
   do ix = n - 1, 1, -1
    dw(ix) = dw(ix) - fmat(ix, 3, dir)*dw(ix + 1)
   end do
   if (ib /= 1) return
   alp0 = fmat(n, 3, dir)
   bet0 = fmat(1, 1, dir)
   gm = -1.0
   zp(2: n - 1) = 0.0
   zp(1) = gm
   zp(n) = alp0
!-----------call trid of z
   zp(1) = fmat(1, 2, dir)*zp(1)
   do ix = 2, n
    zp(ix) = zp(ix) - fmat(ix, 1, dir)*zp(ix - 1)
    zp(ix) = fmat(ix, 2, dir)*zp(ix)
   end do
   do ix = n - 1, 1, -1
    zp(ix) = zp(ix) - fmat(ix, 3, dir)*zp(ix + 1)
   end do
   fact = (dw(1) + bet0*dw(n)/gm)/(1.0 + zp(1) + bet0*zp(n)/gm)
   dw(1: n) = dw(1: n) - fact*zp(1: n)
   dw(n + 1) = dw(1)
  end subroutine
!=======================
  subroutine set_lpl_mat(ng1, dir)
   integer, intent (in) :: ng1, dir
   integer :: i

   lpl_mat(1: ng1, 1, dir) = 1.0
   lpl_mat(1: ng1, 2, dir) = -2.0
   lpl_mat(1: ng1, 3, dir) = 1.0
!        LU Factorize
   lpl_mat(1, 2, dir) = 1.0/lpl_mat(1, 2, dir)
   do i = 2, ng1
    lpl_mat(i - 1, 3, dir) = lpl_mat(i - 1, 3, dir)*lpl_mat(i - 1, 2, &
      dir)
    lpl_mat(i, 2, dir) = lpl_mat(i, 2, dir) - lpl_mat(i - 1, 3, dir)* &
      lpl_mat(i, 1, dir)
    lpl_mat(i, 2, dir) = 1.0/lpl_mat(i, 2, dir)
   end do
  end subroutine
!===================
!=======================
  subroutine set_radlpl_mat(rhg, rg, ng1)
   real (dp), intent (in) :: rhg(: ), rg(: )
   integer, intent (in) :: ng1
   integer :: i

   allocate (rmat(ng1, 3))

   i = 1
   rmat(i, 1) = 0.0
   rmat(i, 2) = -2.0
   rmat(i, 3) = 2.0
   do i = 2, ng1 - 1
    rmat(i, 1) = rhg(i - 1)/rg(i)
    rmat(i, 2) = -(rhg(i-1) + rhg(i))/rg(i)
    rmat(i, 3) = rhg(i)/rg(i)
   end do
   i = ng1 !pot(ng1+1)=0
   rmat(i, 1) = rhg(i - 1)/rg(i)
   rmat(i, 2) = -(rhg(i-1) + rhg(i))/rg(i)
!rmat(i,1)=-(rh(ng1)-rh(ng1-1))  linear extrapolation
!rmat(i,2)=(rh(ng1)-rh(ng1-1))
!        LU Factorize
   rmat(1, 2) = 1.0/rmat(2, 2)
   do i = 2, ng1
    rmat(i - 1, 3) = rmat(i - 1, 3)*rmat(i - 1, 2)
    rmat(i, 2) = rmat(i, 2) - rmat(i - 1, 3)*rmat(i, 1)
    rmat(i, 2) = 1.0/rmat(i, 2)
   end do
  end subroutine
!===================
!----------------------------------------------
!----------------------------------------------
!===============================================
  subroutine set_der_inv(ng, aph, aph1, ib)
   integer, intent (in) :: ng, ib
   real (dp), intent (in) :: aph, aph1
   integer :: i
   real (dp) :: ap, gm, bt

   mat_der_inv(1: ng, 1) = aph
   mat_der_inv(1: ng, 2) = 1.0
   mat_der_inv(1: ng, 3) = aph
   select case (ib)
   case (0)
    mat_der_inv(1, 1) = 0.0
    mat_der_inv(1, 3) = 0.0
    mat_der_inv(2, 1) = aph1
    mat_der_inv(2, 3) = aph1
    mat_der_inv(ng - 1, 1) = aph1
    mat_der_inv(ng - 1, 3) = aph1
    mat_der_inv(ng, 3) = 0.0
    mat_der_inv(ng, 1) = 0.0
   case (1)
!   the cyclic matrix
    ap = mat_der_inv(ng, 3)
    bt = mat_der_inv(1, 1)
    gm = -mat_der_inv(1, 2)
    mat_der_inv(1, 2) = mat_der_inv(1, 2) - gm
    mat_der_inv(ng, 2) = mat_der_inv(ng, 2) - ap*bt/gm
   end select
!        LU Factorize
   mat_der_inv(1, 2) = 1.0/mat_der_inv(1, 2)
   do i = 2, ng
    mat_der_inv(i - 1, 3) = mat_der_inv(i - 1, 3)*mat_der_inv(i - 1, 2)
    mat_der_inv(i, 2) = mat_der_inv(i, 2) - mat_der_inv(i - 1, 3)* &
      mat_der_inv(i, 1)
    mat_der_inv(i, 2) = 1.0/mat_der_inv(i, 2)
   end do

  end subroutine
  subroutine penta_diag_lufact(ng)

   integer, intent (in) :: ng
   integer :: i

   mat_env(1, 3) = 1.0/mat_env(1, 3)
   do i = 2, ng - 1
    if (i == 3) mat_env(i - 1, 4) = mat_env(i - 1, 4) - mat_env(i - 1, 2 &
      )*mat_env(i - 2, 5)
!A(i-1,i)=A(i-1,i)-A(i-1,i-2)*A(i-2,i)
    mat_env(i - 1, 4) = mat_env(i - 1, 4)*mat_env(i - 1, 3)
!A(i-1,i)=A(i-1,i-1)*A(i-1,i)
    mat_env(i, 3) = mat_env(i, 3) - mat_env(i - 1, 4)*mat_env(i, 2)
!A(i,i)=A(i,i)-A(i,i-1)*A(i-1,i)
    mat_env(i, 3) = 1.0/mat_env(i, 3)
   end do
   i = ng
   mat_env(i, 2) = mat_env(i, 2) - mat_env(i, 1)*mat_env(i - 2, 4)
!A(i,i-1)=A(i,i-1)-A(i,i-2)*A(i-2,i-1)
   mat_env(i, 3) = mat_env(i, 3) - mat_env(i - 1, 4)*mat_env(i, 2)
!A(i,i)=A(i,i)-A(i,i-1)*A(i-1,i)
   mat_env(i, 3) = 1.0/mat_env(i, 3)
  end subroutine
!==============================
  subroutine set_mat_env5(a, ng)
   integer, intent (in) :: ng
   real (dp), intent (in) :: a(5, 5)
   integer :: i, j
!==============
!          To invert a penta-diagonal matrix
!          coefficients  a(1,1:3)  first row
!          coefficients  a(2,1:4)  second row
!          coefficients  a(3,1:5)  interior rows
!
   allocate (amat(ng, ng))
   amat = 0.0

   amat(1, 1: 3) = a(1, 1: 3)
   amat(2, 1: 4) = a(2, 1: 4)
   do j = 3, ng - 2
    amat(j, j - 2: j + 2) = a(3, 1: 5)
   end do
   amat(ng - 1, ng - 3: ng) = a(4, 1: 4)
   amat(ng, ng - 2: ng) = a(5, 1: 3)
!        LU Factorize a penta-diagonal matrix
   call ludcmp(amat, ng)
   mat_env(1, 1: 2) = 0.0
   mat_env(1, 3: 5) = amat(1, 1: 3)
   mat_env(2, 2: 5) = amat(2, 1: 4)
   do j = 3, ng - 2
    mat_env(j, 1: 5) = amat(j, j - 2: j + 2)
   end do
   mat_env(ng - 1, 1: 4) = amat(ng - 1, ng - 3: ng)
   mat_env(ng, 1: 3) = amat(ng, ng - 2: ng)
   do i = 1, 5
    do j = 1, ng
     if (abs(mat_env(j,i)) < 1.e-08) mat_env(j, i) = 0.0
    end do
   end do
   deallocate (amat)
  end subroutine
!=======================
  subroutine set_mat_env2(bp, aph, ng)
   integer, intent (in) :: ng
   real (dp), intent (in) :: bp, aph
   integer :: i
   real (dp) :: ap2, a, b, c, b1, c1, d1, en, an, bn
!==============
!==============
   ap2 = bp*bp
   b1 = 1. + ap2*(1. - aph)*(1. - aph)
   c1 = 2.*ap2*aph*(1. - aph)
   d1 = ap2*aph*aph
!                      first row
   a = ap2*aph*(aph - 1.)
   b = 1. + ap2*(1. - 2.*aph*aph)
   c = ap2*aph*(aph + 1.)
!                    !interior rows
   en = d1
   bn = 1. + ap2*(1. + aph)*(1. + aph)
   an = -2.*ap2*aph*(1. + aph)
!                      last row

   mat_env(1, 1: 2) = 0.0
   mat_env(1, 3) = b1
   mat_env(1, 4) = c1
   mat_env(1, 5) = d1
   mat_env(2: ng - 1, 2) = a
   mat_env(2: ng - 1, 3) = b
   mat_env(2: ng - 1, 4) = c
   mat_env(ng, 1) = en
   mat_env(ng, 2) = an
   mat_env(ng, 3) = bn
   mat_env(ng, 4: 5) = 0.0
!        LU Factorize a tri-diagonal matrix

   mat_env(1, 3) = 1.0/mat_env(1, 3)
   do i = 2, ng - 1
    if (i == 3) mat_env(i - 1, 4) = mat_env(i - 1, 4) - mat_env(i - 1, 2 &
      )*mat_env(i - 2, 5)
!A(i-1,i)=A(i-1,i)-A(i-1,i-2)*A(i-2,i)
    mat_env(i - 1, 4) = mat_env(i - 1, 4)*mat_env(i - 1, 3)
!A(i-1,i)=A(i-1,i-1)*A(i-1,i)
    mat_env(i, 3) = mat_env(i, 3) - mat_env(i - 1, 4)*mat_env(i, 2)
!A(i,i)=A(i,i)-A(i,i-1)*A(i-1,i)
    mat_env(i, 3) = 1.0/mat_env(i, 3)
   end do
   i = ng
   mat_env(i, 2) = mat_env(i, 2) - mat_env(i, 1)*mat_env(i - 2, 4)
!A(i,i-1)=A(i,i-1)-A(i,i-2)*A(i-2,i-1)
   mat_env(i, 3) = mat_env(i, 3) - mat_env(i - 1, 4)*mat_env(i, 2)
!A(i,i)=A(i,i)-A(i,i-1)*A(i-1,i)
   mat_env(i, 3) = 1.0/mat_env(i, 3)
  end subroutine
!----------------------------------------------
  subroutine w_alloc(opt_der)
   real (dp), intent (in) :: opt_der
   real (dp) :: cfl_loc, as0, alp
   integer :: ndmx, nd1mx, ng
! sets logical bd flags
   xl_bd = .true.
   xr_bd = .true.
   yl_bd = .false.
   yr_bd = .false.
   zl_bd = .false.
   zr_bd = .false.
   if (pe0y) yl_bd = .true.
   if (pe1y) yr_bd = .true.
   if (pe0z) zl_bd = .true.
   if (pe1z) zr_bd = .true.

!============================
! allocate auxiliary arrays ww() dw() and
   ndmx = max(nx, ny, nz)
   nd1mx = max(nx_loc, nx1_loc)
   allocate (ww(ndmx+5), ww0(ndmx+5, 5+max(ny_loc,nz_loc)))
   allocate (wr(ndmx+6, 10), wl(ndmx+6, 10))
   allocate (var(ndmx+5, 10))
   var(:, : ) = 0.0
   wr(:, : ) = 0.0
   wl(:, : ) = 0.0
   ww(: ) = 0.0
   ww0(:, : ) = 0.0
   cfl_loc = 0.0
   ng = nx

   if (lpf_ord > 0) then
    cfl_loc = cfl
    if (ndim > 1) cfl_loc = cfl*yx_rat/sqrt(yx_rat*yx_rat + float(ndim) &
      - 1.)
   end if
   call set_mat_der(cfl_loc, opt_der, nx, ny, nz, ndim, ibx, der_ord, &
     ifilt, iform, comoving)
!============= k-space data for fft
   allocate (sf1(nx/2+1), sf2(nx/2+1))
   if (envelope) then
    allocate (mat_env(ng, 5))
    alp = dt*oml
    as0 = dt*dx_inv
    call set_env5_coeff(alp, as0)
    call set_mat_env5(mat_env_coeff, ng)
   end if
!============================
  end subroutine
!====================================
  subroutine set_env5_coeff(ap, a)
   real (dp), intent (in) :: ap, a
   real (dp) :: ap2, a2

   ap2 = ap*ap
   a2 = a*a
!======== The discrete operator 1+(dt*k0)*(dt*k0)+dt*dtD_xD_x-2*dt*D_x
!   First row
   mat_env_coeff(1, 1) = 1. + ap2 + 2*a + 0.5*a2
   mat_env_coeff(1, 2) = -a2 - 2.*a
   mat_env_coeff(1, 3) = 0.5*a2

!   second row
   mat_env_coeff(2, 1) = 0.5*a2 + a
   mat_env_coeff(2, 2) = 1. + ap2 + 0.25*a2 + a
   mat_env_coeff(2, 3) = -a
   mat_env_coeff(2, 4) = 0.25*a2

!   interior rows
   mat_env_coeff(3, 1) = 0.25*a2
   mat_env_coeff(3, 2) = a
   mat_env_coeff(3, 3) = 1. + ap2 - 0.5*a2
   mat_env_coeff(3, 4) = -a
   mat_env_coeff(3, 5) = 0.25*a2
!   n-1-th row
   mat_env_coeff(4, 1) = 0.25*a2
   mat_env_coeff(4, 2) = a
   mat_env_coeff(4, 3) = -0.5*a2
   mat_env_coeff(4, 4) = 1. + ap2 + 0.25*a2 - a
!   n-th row
   mat_env_coeff(5, 1) = 0.25*a2
   mat_env_coeff(5, 2) = a - 0.5*a2
   mat_env_coeff(5, 3) = 1 + ap2 - a
!=============================
  end subroutine
!====================
  subroutine matenv_inv(n, nc)
   integer, intent (in) :: n, nc
   integer :: ix, ic
! Invert a pentadiagonal matrix
! mat_env(n,5) LU factorized in set_mat_env5()

   do ic = 1, nc
    ix = 1
    ww0(ix, ic) = mat_env(ix, 3)*ww0(ix, ic)
    ix = 2
    ww0(ix, ic) = ww0(ix, ic) - mat_env(ix, 2)*ww0(ix - 1, ic)
    ww0(ix, ic) = mat_env(ix, 3)*ww0(ix, ic)
    do ix = 3, n
     ww0(ix, ic) = ww0(ix, ic) - mat_env(ix, 1)*ww0(ix - 2, ic) - &
       mat_env(ix, 2)*ww0(ix - 1, ic)
     ww0(ix, ic) = mat_env(ix, 3)*ww0(ix, ic)
    end do
    ix = n - 1
    ww0(ix, ic) = ww0(ix, ic) - mat_env(ix, 4)*ww0(ix + 1, ic)
    do ix = n - 2, 1, -1
     ww0(ix, ic) = ww0(ix, ic) - mat_env(ix, 4)*ww0(ix + 1, ic) - &
       mat_env(ix, 5)*ww0(ix + 2, ic)
    end do
   end do
  end subroutine
!=======================
  subroutine trid2_loc(bd, b1, c1, bn, an, n, nc)
   real (dp), intent (in) :: bd, b1, c1, bn, an
   integer, intent (in) :: n, nc
   integer :: k, ic
   real (dp) :: bet
!==========================
   do ic = 1, nc
    dw(1) = 0.0
    bet = b1
    ww0(1, ic) = ww0(1, ic) - ww0(2, ic)
    ww0(n, ic) = ww0(n, ic) - ww0(n - 1, ic)
    ww0(1, ic) = ww0(1, ic)/bet
    k = 2
    dw(k) = c1/bet
    bet = bd - dw(k)
    ww0(k, ic) = (ww0(k, ic) - ww0(k-1, ic))/bet
    do k = 3, n - 1
     dw(k) = 1./bet
     bet = bd - dw(k)
     ww0(k, ic) = (ww0(k, ic) - ww0(k-1, ic))/bet
    end do
    k = n
    dw(k) = 1./bet
    bet = bn - an*dw(k)
    ww0(k, ic) = (ww0(k, ic) - an*ww0(k-1, ic))/bet
    do k = n - 1, 1, -1
     ww0(k, ic) = ww0(k, ic) - dw(k + 1)*ww0(k + 1, ic)
    end do
   end do
  end subroutine
!===========================
  subroutine trid_odd_even_inv(a, b, c, b1, an, bn, n, ic1, ic2)
! subroutine trid_odd_even_inv(a,b,c,b1,c1,an,bn,n,ic1,ic2)
   real (dp), intent (in) :: a, b, c, b1, an, bn
! real(dp),intent(in) :: c1
   integer, intent (in) :: n, ic1, ic2
   integer :: k, ic, nh, i1, i2
   real (dp) :: bet
!==========================
! Solves
! a*ww(i-1)+b*ww(i)+c*ww(i+1)=u(i), i=2,3,..,n-1
! first order boundary clusure
!===============================
   nh = n/2
   do ic = ic1, ic2
    k = 1
    i2 = 2*k
    i1 = i2 - 1
    dw(k) = 0.0
    bet = b1
    ww0(i1, ic) = ww0(i1, ic)/bet
    ww0(i2, ic) = ww0(i2, ic)/bet
    do k = 2, nh - 1
     i2 = 2*k
     i1 = i2 - 1
     dw(k) = c/bet
     bet = b - a*dw(k)
     ww0(i1, ic) = (ww0(i1, ic) - a*ww0(i1-2, ic))/bet
     ww0(i2, ic) = (ww0(i2, ic) - a*ww0(i2-2, ic))/bet
    end do
    k = nh
    i2 = 2*k
    i1 = i2 - 1
    dw(k) = c/bet
    bet = bn - an*dw(k)
    ww0(i1, ic) = (ww0(i1, ic) - an*ww0(i1-2, ic))/bet
    ww0(i2, ic) = (ww0(i2, ic) - an*ww0(i2-2, ic))/bet
    do k = nh - 1, 1, -1
     i2 = 2*k
     i1 = i2 - 1
     ww0(i1, ic) = ww0(i1, ic) - dw(k + 1)*ww0(i1 + 2, ic)
     ww0(i2, ic) = ww0(i2, ic) - dw(k + 1)*ww0(i2 + 2, ic)
    end do
   end do
  end subroutine
!==================================
  subroutine cycle_trid_inv(a, b, c, n)
   real (dp), intent (in) :: a, b, c
   integer, intent (in) :: n
   integer :: k, ic, ic1
   real (dp) :: b1, bn
   real (dp) :: fact, alp, gm, bet
!==========================
! Solves
! a*ww(i-1)+b*ww(i)+c*ww(i+1)=u(i), i=1,2,..,n
! for periodic boundary clusure
!===============================
   alp = c
   bet = a
   gm = -b
   b1 = b - gm
   bn = b - alp*bet/gm
   ic = 1
   ic1 = ic + 1
!============= diagonal b1, b..,bn
   dw(1) = 0.0
   bet = b1
   ww0(1, ic) = ww0(1, ic)/bet
   do k = 2, n - 1
    dw(k) = c/bet
    bet = b - a*dw(k)
    ww0(k, ic) = (ww0(k, ic) - a*ww0(k-1, ic))/bet
   end do
   k = n
   dw(k) = c/bet
   bet = bn - a*dw(k)
   ww0(k, ic) = (ww0(k, ic) - a*ww0(k-1, ic))/bet
   do k = n - 1, 1, -1
    ww0(k, ic) = ww0(k, ic) - dw(k + 1)*ww0(k + 1, ic)
   end do
!----------------------
   ww0(1, ic1) = gm
   ww0(n, ic1) = alp
   ww0(2: n - 1, ic1) = 0.0
   bet = b1
   ww0(1, ic1) = ww0(1, ic1)/bet
   do k = 2, n - 1
    dw(k) = c/bet
    bet = b - a*dw(k)
    ww0(k, ic1) = (ww0(k, ic1) - a*ww0(k-1, ic1))/bet
   end do
   k = n
   dw(k) = c/bet
   bet = bn - a*dw(k)
   ww0(k, ic1) = (ww0(k, ic1) - a*ww0(k-1, ic1))/bet
   do k = n - 1, 1, -1
    ww0(k, ic1) = ww0(k, ic1) - dw(k + 1)*ww0(k + 1, ic1)
   end do
   fact = (ww0(1, ic) + bet*ww0(n, ic)/gm)/(1. + ww0(1, ic1) + bet*ww0(n &
     , ic1)/gm)
   ww0(1: n, ic) = ww0(1: n, ic) - fact*ww0(1: n, ic1)
  end subroutine
!=================================
  subroutine upper_trid(a, b, c, bn, cn, n, ic1, ic2)
   real (dp), intent (in) :: a, b, c, bn, cn
   integer, intent (in) :: n, ic1, ic2
   integer :: k, ic
!==========================
! Solves
! a*ww(i)+b*ww(i+1)+c*ww(i+2)=u(i), i=1,3,..,n-2
! at the n-1 row bn*ww(n-1)+cn*ww(n)=u(n-1)
! at the n row ww(n)=u(n)
! first order right boundary clusure
!===============================
   do ic = ic1, ic2
    k = n - 1
    ww0(k, ic) = ww0(k, ic) - cn*ww0(k + 1, ic)
    ww0(k, ic) = ww0(k, ic)/bn
    do k = n - 2, 1, -1
     ww0(k, ic) = ww0(k, ic) - b*ww0(k + 1, ic) - c*ww0(k + 2, ic)
     ww0(k, ic) = ww0(k, ic)/a
    end do
   end do
  end subroutine
!================
  subroutine trid_der1(a, b, c, b1, c1, an, bn, n, ic1, ic2, ord)
   real (dp), intent (in) :: a, b, c, b1, c1, an, bn
   integer, intent (in) :: n, ic1, ic2, ord
   integer :: k, ic
   real (dp) :: bet
!==========================
! Solves
! a*ww(i-1)+b*ww(i)+c*ww(i+1)=u(i), i=2,3,..,n-1
! at the first row b1*ww(1)+c1*ww(2)=u(1)
! at the n-last row an*ww(n-1)+bn*ww(n)=u(n)
! first order boundary clusure
!===============================
   if (ord > 0) then
    do ic = ic1, ic2
     ww0(1, ic) = ww0(1, ic) + ww0(2, ic)
     ww0(n, ic) = ww0(n, ic) + ww0(n - 1, ic)
    end do
   end if
!===================
   do ic = ic1, ic2
    dw(1) = 0.0
    bet = b1
    ww0(1, ic) = ww0(1, ic)/bet
    k = 2
    dw(k) = c1/bet
    bet = b - a*dw(k)
    ww0(k, ic) = (ww0(k, ic) - a*ww0(k-1, ic))/bet
    do k = 3, n - 1
     dw(k) = c/bet
     bet = b - a*dw(k)
     ww0(k, ic) = (ww0(k, ic) - a*ww0(k-1, ic))/bet
    end do
    k = n
    dw(k) = c/bet
    bet = bn - an*dw(k)
    ww0(k, ic) = (ww0(k, ic) - an*ww0(k-1, ic))/bet
    do k = n - 1, 1, -1
     ww0(k, ic) = ww0(k, ic) - dw(k + 1)*ww0(k + 1, ic)
    end do
   end do
  end subroutine
!================
  subroutine hder_inv(ngm, ib)
   integer, intent (in) :: ngm, ib
   integer :: ix
!============== ww collocated at half-integers; ww1 collocated at integers
! Enter ww(ngm+1) right boundary
   select case (ib)
   case (0)
    dw(1) = ww(1)
    do ix = 2, ngm - 1
     dw(ix) = avg_cmp*ww(ix)
    end do
    dw(ngm) = ww(ngm)
    if (derinv) call trid_der_inv(ngm, ib)
    ww(ngm) = -ww(ngm + 1)
    do ix = ngm - 1, 1, -1
     ww(ix) = ww(ix + 1) - dw(ix + 1)
    end do
!==========================
   case (1)
    dw(0) = ww(ngm + 1) !the current average <J>
!                    !periodicity conditions
    ww(0) = ww(ngm - 1)
    ww(ngm) = ww(1)
    do ix = 1, ngm - 1
     dw(ix) = avg_cmp*ww(ix)
    end do
    if (derinv) call trid_der_inv(ngm - 1, 1)
    ww(1) = dw(1)
    do ix = 2, ngm - 1
     ww(ix) = ww(ix - 1) + dw(ix)
    end do
    ww(0) = (dw(0) - sum(ww(1:ngm-1)))/real(ngm - 1, dp)
!ww(i)=DJ_i=J_{i+1/2)-J_{0+1/2}
!ww(0)=<J>-<DJ>= J_{0+1/2}
    do ix = 1, ngm - 1
     ww(ix) = ww(ix) + ww(0)
    end do
    ww(ngm) = ww(1)
   end select
!===========================
  end subroutine
  subroutine fluid_filter(fdata, i1, if2, jf1, jf2, kf1, kf2, ic1, ic2)
   real (dp), intent (inout) :: fdata(:, :, :, : )
   integer, intent (in) :: i1, if2, jf1, jf2, kf1, kf2, ic1, ic2
   real (dp) :: aph
   integer :: i, ii, j, k, ng, ic, j01, j02, k01, k02

!  filtered data overwrite input data
!===========================
!  filtering data along the x coordinte  i1=3  if2 assigned on input data
!  available at if2+1,if2+2 grid point
!===========================================
   j01 = jf1 - 2
   if (yl_bd) j01 = jf1
   j02 = jf2 + 2
   if (yr_bd) j02 = jf2
   if (ndim == 2) then
    k01 = kf1
    k02 = kf2
   end if
   if (ndim == 3) then
    k01 = kf1 - 2
    if (zl_bd) k01 = kf1
    k02 = kf2 + 2
    if (zr_bd) k02 = kf2
   end if
!===================
   aph = filt_coeff(0)
   ng = if2 - 4
   do ic = ic1, ic2
    do k = k01, k02
     do j = j01, j02
      do = i1 + 2, if2
      ii = i - 4
      ww0(ii, j) = filt_coeff(4)*fdata(i, j, k, ic) + filt_coeff(5)*( &
        fdata(i-1, j, k, ic) + fdata(i+1, j, k, ic)) + filt_coeff(6)*( &
        fdata(i-2, j, k, ic) + fdata(i+2, j, k, ic))
     end do
     i = i1 + 1
     fdata(i, j, k, ic) = 0.5*fdata(i, j, k, ic) + 0.25*(fdata(i+1, j, k &
       , ic) + fdata(i-1, j, k, ic))
!==============
    end do
    do j = j01, j02
     do i = i1 + 2, if2
      ii = i - 4
      fdata(i, j, k, ic) = ww0(ii, j)
     end do
    end do
   end do
  end do
!============== filtering along the y coordinate
  do ic = ic1, ic2
   do k = k01, k02
    do j = j01 + 2, j02 - 2 !jf1,jf2 for interior points
     do i = i1, if2
      ii = i - 2
      ww0(ii, j) = filt_coeff(4)*fdata(i, j, k, ic) + filt_coeff(5)*( &
        fdata(i, j-1, k, ic) + fdata(i, j+1, k, ic)) + filt_coeff(6)*( &
        fdata(i, j-2, k, ic) + fdata(i, j+2, k, ic))
     end do
    end do
!==============
    do j = j01 + 2, j02 - 2
     = i1, if2
     ii = i - 2
     fdata(i, j, k, ic) = ww0(ii, j)
    end do
   end do
  end do
 end do
 if (ndim < 3) return
!============== filtering along the z coordinate
 do ic = ic1, ic2
  do j = j01, j02
   do k = k01 + 2, k02 - 2
    do i = i1, if2
     ii = i - 2
     ww0(ii, k) = filt_coeff(4)*fdata(i, j, k, ic) + filt_coeff(5)*( &
       fdata(i, j, k-1, ic) + fdata(i, j, k+1, ic)) + filt_coeff(6)*( &
       fdata(i, j, k-2, ic) + fdata(i, j, k+2, ic))
    end do
   end do
!==============
   do k = k01 + 2, k02 - 2
    do i = i1, if2
     ii = i - 2
     fdata(i, j, k, ic) = ww0(ii, k)
    end do
   end do
  end do
 end do

end subroutine
!===================================
end module

