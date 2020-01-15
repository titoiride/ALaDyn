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
!===================================================
!     Local grid structure under mpi domain decomposition
!============== 
! grid [1:n]   np=n+2  extended domain [1:np+2]
! interior [3,np]   ghost [1:2], [np+1:np+2]  
! 
!             
!             overlapping grid points
!====================================================================
!                                     1-----2---- 3--- 4   |      pey+1
!                 1-----2----[3-------np-1--np]--np+1--np+2|    pey
!1-2--------------np-1--np---np+1                          |pey-1
!====================  Fill  extended grid data  ======================================
!      Right(pey+1)     [3:4] ==>      pey [np+1:np+2]  right ghost data
!      Left(pey-1)      [np+1:np+2]==> pey [1:2]        left  ghost data
!===================================


 module pic_rutil

  use array_alloc
  use grid_param
  use code_util
  use parallel

  implicit none
  real (dp), parameter :: tol0 = 0.02, tol_max = 0.1
  real (dp) :: loc_pstore(7)
  integer, parameter :: x_parity(6) = [ -1, 1, 1, -1, 1, 1 ]
  integer, parameter :: y_parity(6) = [ 1, -1, 1, 1, -1, 1 ]
  integer, parameter :: z_parity(6) = [ 1, 1, -1, 1, 1, -1 ]

 contains
!===========================
  subroutine field_xyzbd(ef, i01, i02, j01, j02, k01, k02, nc, stl, str)
   real (dp), intent (inout) :: ef(:, :, :, :)
   integer, intent (in) :: i01, i02, j01, j02, k01, k02, nc, stl, str
   integer :: ik, iy, ix, iz
   integer :: i1, i2, j1, j2, k1, k2
!==========================
! enter fields(1:n1,j01-1:j02+1,k01-1:k02+1,nc)
! Only for NON-PERIODIC boundaries
   j1 = j01
   j2 = j02
   k1 = k01
   k2 = k02
   i1 = i01
   i2 = i02
   if (prly) then
    j1 = j01 - stl
    j2 = j02 + str
   end if
   if (prlz) then
    k1 = k01 - stl
    k2 = k02 + str
   end if
   i1 = i01
   i2 = i02
   if (prlx) then
    i1 = i01 - stl
    i2 = i02 + str
   end if
   if (pex0) then
    if (ibx<2) then
     do ik = 1, nc
      do iz = k1, k2
       do iy = j1, j2
        ef(i01-1, iy, iz, ik) = 2.*ef(i01, iy, iz, ik) - &
          ef(i01+1, iy, iz, ik)
       end do
      end do
     end do
     i1 = i01 - 1
    end if
   end if
   if (pex1) then
    if (ibx==0) then
     do ik = 1, nc
      do iz = k1, k2
       do iy = j1, j2
        ef(i02+1, iy, iz, ik) = 2.*ef(i02, iy, iz, ik) - &
          ef(i02-1, iy, iz, ik)
       end do
      end do
     end do
    end if
    if (ibx==1) then
     do ik = 1, nc
      do iz = k1, k2
       do iy = j1, j2
        ef(i02+1, iy, iz, ik) = x_parity(ik)*ef(i02, iy, iz, ik)
       end do
      end do
     end do
     i2 = i02 + 1
    end if
   end if
   if (ndim<2) return
   if (pe0y) then
    if (iby==0) then
     do ik = 1, nc
      do iz = k1, k2
       do ix = i1, i2
        ef(ix, j01-1, iz, ik) = 2.*ef(ix, j01, iz, ik) - &
          ef(ix, j01+1, iz, ik)
       end do
      end do
     end do
    end if
    if (iby==1) then
     do ik = 1, nc
      do iz = k1, k2
       do ix = i1, i2
        ef(ix, j01-1, iz, ik) = y_parity(ik)*ef(ix, j01+1, iz, ik)
       end do
      end do
     end do
    end if
   end if
   if (pe1y) then
    if (iby==0) then
     do ik = 1, nc
      do iz = k1, k2
       do ix = i1, i2
        ef(ix, j02+1, iz, ik) = 2.*ef(ix, j02, iz, ik) - &
          ef(ix, j02-1, iz, ik)
       end do
      end do
     end do
    end if
    if (iby==1) then
     do ik = 1, nc
      do iz = k01, k02
       do ix = i1, i2
        ef(ix, j02+1, iz, ik) = y_parity(ik)*ef(ix, j02-1, iz, ik)
       end do
      end do
     end do
    end if
   end if
   if (ndim<3) return
   if (pe0z) then
    if (ibz==0) then
     do ik = 1, nc
      do iy = j1, j2
       do ix = i1, i2
        ef(ix, iy, k01-1, ik) = 2*ef(ix, iy, k01, ik) - &
          ef(ix, iy, k01+1, ik)
       end do
      end do
     end do
    end if
    if (ibz==1) then
     do ik = 1, nc
      do iy = j1, j2
       do ix = i1, i2
        ef(ix, iy, k01-1, ik) = z_parity(ik)*ef(ix, iy, k01+1, ik)
       end do
      end do
     end do
    end if
   end if
   if (pe1z) then
    if (ibz==0) then
     do ik = 1, nc
      do iy = j1, j2
       do ix = i1, i2
        ef(ix, iy, k02+1, ik) = 2*ef(ix, iy, k02, ik) - &
          ef(ix, iy, k02-1, ik)
       end do
      end do
     end do
    end if
    if (ibz==1) then
     do ik = 1, nc
      do iy = j1, j2
       do ix = i1, i2
        ef(ix, iy, k01+1, ik) = z_parity(ik)*ef(ix, iy, k01-1, ik)
       end do
      end do
     end do
    end if
   end if
  end subroutine
!===============================
  subroutine jc_xyzbd(curr, i1, n1, j1, n2, k1, n3, nc)
   real (dp), intent (inout) :: curr(:, :, :, :)
   integer, intent (in) :: i1, n1, j1, n2, k1, n3, nc
   integer :: ix, iy, iz, i0, j2, k2, ik
! Enter current data on extended ranges:
!========== Only for Periodic BDs period=n1-1
   j2 = n2
   k2 = n3
   if (ibx==0) then
    do ik = 1, nc
     curr(i1, j1:j2, k1:k2, ik) = curr(i1, j1:j2, k1:k2, ik) + &
       curr(i1-1, j1:j2, k1:k2, ik)
     curr(i1-1, j1:j2, k1:k2, ik) = 0.0

     curr(n1, j1:j2, k1:k2, ik) = curr(n1, j1:j2, k1:k2, ik) + &
       curr(n1+1, j1:j2, k1:k2, ik) + curr(n1+2, j1:j2, k1:k2, ik)
     curr(n1+1:n1+2, j1:j2, k1:k2, ik) = 0.0
    end do
    i0 = 1
   end if
   if (ndim<2) return
   if (pe0y) then
    if (iby==0) then
     do ik = 1, nc
      do iz = k1, k2
       do ix = i1, n1
        curr(ix, j1, iz, ik) = curr(ix, j1, iz, ik) + &
          curr(ix, j1-1, iz, ik) + curr(ix, j1-2, iz, ik)
        curr(ix, j1-2:j1-1, iz, ik) = 0.0
       end do
      end do
     end do
    end if
    if (iby==1) then
!Before norm  r*Jx(j)=rVxrho ==> rEx odd Jx(j-1)=-Jx(j+1)
!             r*Jr=rVrrho ==> rEr even   Jr(j-1/2)=jr(j+1/2)
     do iz = k1, k2
      do ix = i1, n1
       curr(ix, j1+1, iz, 1) = curr(ix, j1+1, iz, 1) + &
         curr(ix, j1-1, iz, 1)
       curr(ix, j1, iz, 2) = curr(ix, j1, iz, 2) - curr(ix, j1-1, iz, 2)
      end do
     end do
    end if
   end if
   if (pe1y) then
    if (iby==0) then
     do ik = 1, nc
      do iz = k1, k2
       do ix = i1, n1
        curr(ix, n2, iz, ik) = curr(ix, n2, iz, ik) + &
          curr(ix, n2+1, iz, ik) + curr(ix, n2+2, iz, ik)
        curr(ix, n2+1:n2+2, iz, ik) = 0.0
       end do
      end do
     end do
    end if
   end if
   if (ndim<3) return
   if (ibz==0) then
    if (pe0z) then
     do ik = 1, nc
      do iy = j1, j2
       do ix = i1, n1
        curr(ix, iy, k1, ik) = curr(ix, iy, k1, ik) + &
          curr(ix, iy, k1-1, ik)
        curr(ix, iy, k1-1, ik) = 0.0
       end do
      end do
     end do
    end if
    if (pe1z) then
     if (ibz==0) then
      do ik = 1, nc
       do iy = j1, j2
        do ix = i1, n1
         curr(ix, iy, n3, ik) = curr(ix, iy, n3, ik) + &
           curr(ix, iy, n3+1, ik) + curr(ix, iy, n3+2, ik)
         curr(ix, iy, n3+1:n3+2, ik) = 0.0
        end do
       end do
      end do
     end if
    end if
   end if
  end subroutine
!-----------------------------------------------
!
  subroutine den_zyxbd(rho, i1, i2, j1, j2, k1, k2, ik)
   real (dp), intent (inout) :: rho(:, :, :, :)
   integer, intent (in) :: i1, i2, j1, j2, k1, k2, ik
   integer :: ix, iy
! Enter current data on extended ranges:

!Enter data on the computational box [i1:n1p][j1:nyp][k1:nzp]
   if (ndim>2) then
    if (ibz==0) then
     if (pe0z) then
      do iy = j1, j2
       do ix = i1, i2
        rho(ix, iy, k1+1, ik) = rho(ix, iy, k1+1, ik) + &
          rho(ix, iy, k1-1, ik)
        rho(ix, iy, k1, ik) = rho(ix, iy, k1+1, ik)
       end do
      end do
     end if
     if (pe1z) then
      if (ibz<2) then
       do iy = j1, j2
        do ix = i1, i2
         rho(ix, iy, k2-1, ik) = rho(ix, iy, k2-1, ik) + &
           rho(ix, iy, k2+1, ik)
         rho(ix, iy, k2, ik) = rho(ix, iy, k2-1, ik)
        end do
       end do
      end if
     end if
    end if
   end if
!================
   if (ndim>1) then
    if (pe0y) then
     if (iby==0) then
      rho(i1:i2, j1+1, k1:k2, ik) = rho(i1:i2, j1+1, k1:k2, ik) + &
        rho(i1:i2, j1-1, k1:k2, ik)
      rho(i1:i2, j1, k1:k2, ik) = rho(i1:i2, j1+1, k1:k2, ik)
     end if
     if (iby==1) then
      rho(i1:i2, j1+1, k1:k2, ik) = rho(i1:i2, j1+1, k1:k2, ik) + &
        rho(i1:i2, j1-1, k1:k2, ik)
      rho(i1:i2, j1-1, k1:k2, ik) = 0.0
     end if
    end if
    if (pe1y) then
     if (iby<2) then
      rho(i1:i2, j2-1, k1:k2, ik) = rho(i1:i2, j2-1, k1:k2, ik) + &
        rho(i1:i2, j2+1, k1:k2, ik)
      rho(i1:i2, j2, k1:k2, ik) = rho(i1:i2, j2-1, k1:k2, ik)
     end if
    end if
   end if
!============== field data on [y_loc]
   if (ibx<2) then
    if (pex0) then
     rho(i1+1, j1:j2, k1:k2, ik) = rho(i1+1, j1:j2, k1:k2, ik) + &
       rho(i1-1, j1:j2, k1:k2, ik)
     rho(i1, j1:j2, k1:k2, ik) = rho(i1+1, j1:j2, k1:k2, ik)
    end if
    if (pex1) then
     rho(i2-1, j1:j2, k1:k2, ik) = rho(i2-1, j1:j2, k1:k2, ik) + &
       rho(i2+1, j1:j2, k1:k2, ik)
     rho(i2, j1:j2, k1:k2, ik) = rho(i2-1, j1:j2, k1:k2, ik)
    end if
   end if
  end subroutine
!====================
  subroutine fluid_left_xshift(fld, den_x, den_yz, i1, i2, j1, j2, k1, &
    k2, ic1, ic2, xsh)

   integer, intent (in) :: i1, i2, j1, j2, k1, k2, ic1, ic2, xsh
   real (dp), intent (inout) :: fld(:, :, :, :)
   real (dp), intent (in) :: den_x(:), den_yz(:, :)
   integer :: ic, ix, j, iy, iz, kk, lenws


   if (xsh==0) return
   lenws = (ic2+1-ic1)*(k2+1-k1)*(j2+1-j1)*xsh
   if (prlx) then
!Sends to x-left side xsh data (i1:i1+xsh-1)
!Recvs from x-right data(i2+1:i2+xsh)       i2=n1p-xsh
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 0, xsh - 1
        ix = i1 + j
        kk = kk + 1
        aux1(kk) = fld(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    call exchange_bdx_data(aux1, aux2, lenws, lenws, 3, right)
   end if
!
! shifts (i1+xsx:i2+xsh=n1p)=>(i1:i2)
   do ic = ic1, ic2 - 1
    do iz = k1, k2
     do iy = j1, j2
      fld(i1-1, iy, iz, ic) = fld(i1+xsh-1, iy, iz, ic)
      do ix = i1, i2
       fld(ix, iy, iz, ic) = fld(ix+xsh, iy, iz, ic)
      end do
      fld(i1, iy, iz, ic) = 0.5*fld(i1, iy, iz, ic) + &
        0.25*(fld(i1+1,iy,iz,ic)+fld(i1-1,iy,iz,ic))
      fld(i2+1:i2+xsh, iy, iz, ic) = 0.0
     end do
    end do
   end do
   ic = ic2
   do iz = k1, k2
    do iy = j1, j2
     fld(i1-1, iy, iz, ic) = fld(i1+xsh-1, iy, iz, ic)
     do ix = i1, i2
      fld(ix, iy, iz, ic) = fld(ix+xsh, iy, iz, ic)
     end do
     fld(i1, iy, iz, ic) = 0.5*fld(i1, iy, iz, ic) + &
       0.25*(fld(i1+1,iy,iz,ic)+fld(i1-1,iy,iz,ic))
     do ix = i2 + 1, i2 + xsh
      fld(ix, iy, iz, ic) = den_x(ix)*den_yz(iy, iz)
     end do
    end do
   end do
! now replaces (i2+1:i2+xsh=n1p)
   if (prlx) then
    if (pex1) aux2(1:lenws) = 0.0
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 1, xsh
        ix = i2 + j
        kk = kk + 1
        fld(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if

  end subroutine
!==================================
  subroutine fields_left_xshift(fld, i1, i2, j1, j2, k1, k2, ic1, ic2, &
    xsh)

   integer, intent (in) :: i1, i2, j1, j2, k1, k2, ic1, ic2, xsh
   real (dp) :: fld(:, :, :, :)
   integer :: ic, ix, j, iy, iz, kk, lenws


   if (xsh==0) return
   lenws = (ic2+1-ic1)*(k2+1-k1)*(j2+1-j1)*xsh
   if (prlx) then
!Sends to x-left side xsh data (i1:i1+xsh-1)
!Recvs from x-right data(i2+1:i2+xsh)       i2=n1p-sh
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 0, xsh - 1
        ix = i1 + j
        kk = kk + 1
        aux1(kk) = fld(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    call exchange_bdx_data(aux1, aux2, lenws, lenws, 3, right)
   end if
!
! shifts (i1+xsh:i2+xsh=n1p)=>(i1:i2)
   do ic = ic1, ic2
    do iz = k1, k2
     do iy = j1, j2
      do ix = i1, i2
       fld(ix, iy, iz, ic) = fld(ix+xsh, iy, iz, ic)
      end do
      fld(i2+1:i2+xsh, iy, iz, ic) = 0.0
     end do
    end do
   end do
! now replaces (i2+1:i2+xsh=n1p)
   if (prlx) then
    if (pex1) aux2(1:lenws) = 0.0
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 1, xsh
        ix = i2 + j
        kk = kk + 1
        fld(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if

  end subroutine
!==================================
  subroutine fill_flux_yzxbdsdata(flx, i1, i2, j1, j2, k1, k2, ic1, ic2, &
    str, stl)

   real (dp), intent (inout) :: flx(:, :, :, :)
   integer, intent (in) :: i1, i2, j1, j2, k1, k2, ic1, ic2, str, stl
   integer :: ic, ix, j, iy, iz, kk, iy1, iy2, iz1, iz2, lenws, lenwr

   iy1 = j1
   iy2 = j2
   iz1 = k1
   iz2 = k2
! WARNING str <4 , stl>3 allowed
!=======================
! Extends data to the y-left
! sends to the right iy=[j2:j2-str+1]
!recvs from the left iy=[j1-1:j1-str] sign=-1
!============================================
   if (prly) then
    if (str>0) then
     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 0, str - 1
        iy = j2 - j
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = flx(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 1, left)
     if (pe0y) then
      if (iby<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 1, str
        iy = j1 - j
        do ix = i1, i2
         kk = kk + 1
         flx(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iy1 = j1 - str
    end if
!---------------------
    if (stl>0) then
!=======================
! Extends data to the y-right, stl>2 allowed
!Sends to the left stl data[j1:j1+stl-1]
!Recvs from right [j2+1:j2+stl] f_data sign=+1

     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 0, stl - 1
        iy = j1 + j
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = flx(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws

!======================= next indx=1 cart dim=1 sign=+1
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 1, right)
     if (pe1y) then
      if (iby<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 1, stl
        iy = j + j2
        do ix = i1, i2
         kk = kk + 1
         flx(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iy2 = j2 + stl
    end if
   end if
!========================
   if (prlz) then
!=======================
! Extends data to the z-left
! sends to the right iz=[k2:k2-str+1]
!recvs from the left iz=[k1-1:k1-str] sign=-1
!============================================

    if (str>0) then
     kk = 0
     do ic = ic1, ic2
      do j = 0, str - 1
       iz = k2 - j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = flx(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 2, left)
     if (pe0z) then
      if (ibz<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do j = 1, str
       iz = k1 - j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         flx(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iz1 = k1 - str
    end if
    if (stl>0) then
! Extends data to the z-right
!Sends to the left stl data[k1:k1+stl-1]
!Recvs from right [k2+1:k2+stl] f_data sign=+1
     kk = 0
     do ic = ic1, ic2
      do j = 0, stl - 1
       iz = k1 + j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = flx(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 2, right)
     if (pe1z) then
      if (ibz<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do j = 1, stl
       iz = k2 + j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         flx(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iz2 = k2 + stl
    end if
   end if
!===============================
   if (.not. prlx) then
    if (ibx==2) then
     if (str>0) then
      do ic = ic1, ic2
       do iz = iz1, iz2
        do iy = iy1, iy2
         do j = 1, str
          flx(i1-j, iy, iz, ic) = flx(i2+1-j, iy, iz, ic)
         end do
        end do
       end do
      end do
     end if
     if (stl>0) then
      do ic = ic1, ic2
       do iz = iz1, iz2
        do iy = iy1, iy2
         do j = 1, stl
          flx(i2+j, iy, iz, ic) = flx(i1-1+j, iy, iz, ic)
         end do
        end do
       end do
      end do
     end if
    end if
    return
   end if
!====================================
! Extends data to the x-left
! sends to the right ix=[i2:i2-str+1]
!recvs from the left ix=[i1-1:i1-str] sign=-1
!============================================
   if (str>0) then
    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 0, str - 1
        ix = i2 - j
        kk = kk + 1
        aux1(kk) = flx(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, left)
    if (pex0) then
     if (ibx<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str
        ix = i1 - j
        kk = kk + 1
        flx(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if
!---------------------
   if (stl>0) then
!=======================
! Extends data to the x-right
!Sends to the left stl data[i1:i1+stl-1]
!Recvs from right [i2+1:i2+stl] f_data sign=+1

    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 0, stl - 1
        ix = i1 + j
        kk = kk + 1
        aux1(kk) = flx(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws

    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, right)
    if (pex1) then
     if (ibx<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, stl
        ix = j + i2
        kk = kk + 1
        flx(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if
! end data transfer
!===========================
  end subroutine
!=================================
  subroutine fill_ebfield_yzxbdsdata(fld, i1, i2, j1, j2, k1, k2, ic1, &
    ic2, str, stl)

   integer, intent (in) :: i1, i2, j1, j2, k1, k2, ic1, ic2, str, stl
   real (dp) :: fld(:, :, :, :)
   integer :: ic, ix, j, iy, iz, kk, iy1, iy2, iz1, iz2, lenws, lenwr

   iy1 = j1
   iy2 = j2
   iz1 = k1
   iz2 = k2
! WARNING str <3 , stl>2 allowed
!=======================
! Extends data to the y-left
! sends to the right iy=[j2:j2-str+1]
!recvs from the left iy=[j1-1:j1-str] sign=-1
!============================================
   if (prly) then
    if (str>0) then
     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 0, str - 1
        iy = j2 - j
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = fld(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 1, left)
     if (pe0y) then
      if (iby<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 1, str
        iy = j1 - j
        do ix = i1, i2
         kk = kk + 1
         fld(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iy1 = j1 - str
    end if
!---------------------
    if (stl>0) then
!=======================
! Extends data to the y-right, stl>2 allowed
!Sends to the left stl data[j1:j1+stl-1]
!Recvs from right [j2+1:j2+stl] f_data sign=+1

     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 0, stl - 1
        iy = j1 + j
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = fld(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws

!======================= next indx=1 cart dim=1 sign=+1
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 1, right)
     if (pe1y) then
      if (iby<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do iz = k1, k2
       do j = 1, stl
        iy = j + j2
        do ix = i1, i2
         kk = kk + 1
         fld(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iy2 = j2 + stl
    end if
   end if
!========================
   if (prlz) then
!=======================
! Extends data to the z-left
! sends to the right iz=[k2:k2-str+1]
!recvs from the left iz=[k1-1:k1-str] sign=-1
!============================================

    if (str>0) then
     kk = 0
     do ic = ic1, ic2
      do j = 0, str - 1
       iz = k2 - j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = fld(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 2, left)
     if (pe0z) then
      if (ibz<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do j = 1, str
       iz = k1 - j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         fld(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iz1 = k1 - str
    end if
    if (stl>0) then
! Extends data to the z-right
!Sends to the left stl data[k1:k1+stl-1]
!Recvs from right [k2+1:k2+stl] f_data sign=+1
     kk = 0
     do ic = ic1, ic2
      do j = 0, stl - 1
       iz = k1 + j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         aux1(kk) = fld(ix, iy, iz, ic)
        end do
       end do
      end do
     end do
     lenws = kk
     lenwr = lenws
     call exchange_bdx_data(aux1, aux2, lenws, lenwr, 2, right)
     if (pe1z) then
      if (ibz<2) then
       aux2(1:lenwr) = 0.0
      end if
     end if
     kk = 0
     do ic = ic1, ic2
      do j = 1, stl
       iz = k2 + j
       do iy = iy1, iy2
        do ix = i1, i2
         kk = kk + 1
         fld(ix, iy, iz, ic) = aux2(kk)
        end do
       end do
      end do
     end do
     iz2 = k2 + stl
    end if
   end if
!===============================


   if (.not. prlx) then
    if (ibx==2) then
     if (str>0) then
      do ic = ic1, ic2
       do iz = iz1, iz2
        do iy = iy1, iy2
         do j = 1, str
          fld(i1-j, iy, iz, ic) = fld(i2+1-j, iy, iz, ic)
         end do
        end do
       end do
      end do
     end if
     if (stl>0) then
      do ic = ic1, ic2
       do iz = iz1, iz2
        do iy = iy1, iy2
         do j = 1, stl
          fld(i2+j, iy, iz, ic) = fld(i1-1+j, iy, iz, ic)
         end do
        end do
       end do
      end do
     end if
    end if
    return
   end if
!====================================
! Extends data to the x-left
! sends to the right ix=[i2:i2-str+1]
!recvs from the left ix=[i1-1:i1-str] sign=-1
!============================================
   if (str>0) then
    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 0, str - 1
        ix = i2 - j
        kk = kk + 1
        aux1(kk) = fld(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, left)
    if (pex0) then
     if (ibx<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str
        ix = i1 - j
        kk = kk + 1
        fld(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if
!---------------------
   if (stl>0) then
!=======================
! Extends data to the x-right
!Sends to the left stl data[i1:i1+stl-1]
!Recvs from right [i2+1:i2+stl] f_data sign=+1

    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 0, stl - 1
        ix = i1 + j
        kk = kk + 1
        aux1(kk) = fld(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws

    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, right)
    if (pex1) then
     if (ibx<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    kk = 0
    do ic = ic1, ic2
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, stl
        ix = j + i2
        kk = kk + 1
        fld(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if
! end data transfer
!===========================
  end subroutine
!===============================
  subroutine fill_ebfield_xbdsdata(fld, i1, i2, j1, j2, k1, k2, ic1, &
    ic2, str, stl)
   integer, intent (in) :: i1, i2, j1, j2, k1, k2, ic1, ic2, str, stl
   real (dp) :: fld(:, :, :, :)
   integer :: ic, ix, j, iy, iz, kk, lenws, lenwr
!----------------------------
   if (.not. prlx) then
    if (ibx==2) then
     if (str>0) then
      do ic = ic1, ic2
       do iz = k1, k2
        do iy = j1, j2
         do j = 1, str
          fld(i1-j, iy, iz, ic) = fld(i2+1-j, iy, iz, ic)
         end do
        end do
       end do
      end do
     end if
     if (stl>0) then
      do ic = ic1, ic2
       do iz = k1, k2
        do iy = j1, j2
         do j = 1, stl
          fld(i2+j, iy, iz, ic) = fld(i1-1+j, iy, iz, ic)
         end do
        end do
       end do
      end do
     end if
    end if
    return
   end if
!====================================
   if (str>0) then
!Sends to x-right str data(i2+1-str)
!Recvs from the x-left str data(i1-str) sign=-1
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 0, str - 1
        ix = i2 - j
        kk = kk + 1
        aux1(kk) = fld(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, left)
    if (pex0) then
     if (ibx<2) aux2(1:lenws) = 0.0
    end if
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 1, str
        ix = i1 - j
        kk = kk + 1
        fld(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if
   if (stl>0) then
!Sends to left data(ii1:i1+stl-1)
!Recvs from the right str data(i2:i2+stl)
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 0, stl - 1
        ix = i1 + j
        kk = kk + 1
        aux1(kk) = fld(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, right)
    if (pex1) then
     if (ibx<2) aux2(1:lenws) = 0.0
    end if
    kk = 0
    do ic = ic1, ic2
     do iz = k1, k2
      do iy = j1, j2
       do j = 1, stl
        ix = i2 + j
        kk = kk + 1
        fld(ix, iy, iz, ic) = aux2(kk)
       end do
      end do
     end do
    end do
   end if
! end data transfer
!===========================
  end subroutine
!=====================
  subroutine fill_curr_yzxbdsdata(curr, i1, nxc, j1, nyc, k1, nzc, nc)
   integer, intent (in) :: i1, nxc, j1, nyc, k1, nzc, nc
   real (dp), intent (inout) :: curr(:, :, :, :)
   integer :: s1, s2, r1, r2, iy1, iy2, iz1, iz2, ix1, ix2
   integer :: ic, ix, j, iy, iz, kk, lenws, lenwr
   integer, parameter :: str = 3, str2 = 2
!================
! enter currents on a five-point extended stencil
!===========================
   iz1 = k1 - str2
   iz2 = nzc + str
   iy1 = j1 - str2
   iy2 = nyc + str
   if (ndim<3) then
    iz1 = k1
    iz2 = nzc
   end if
   if (ndim<2) then
    iy1 = k1
    iy2 = nyc
   end if
   ix1 = i1 - str2
   ix2 = nxc + str
   lenwr = str*nc*(ix2+1-ix1)*max(iz2+1-iz1, iy2+1-iy1)
   if (size(aux1)<lenwr) then
    deallocate (aux1, aux2)
    allocate (aux1(lenwr))
    allocate (aux2(lenwr))
   end if
   if (prly) then
! [j1-2:j1-1] left-y data
    s1 = j1 - str
    kk = 0
    do ic = 1, nc
     do iz = iz1, iz2
      do j = 1, str2
       iy = s1 + j
       do ix = ix1, ix2
        kk = kk + 1
        aux1(kk) = curr(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 1, right)
!=====================
! sends y=[j1-2;j1-1] str-1 data to left
!receives from right and adds data on y=[nyc-1:nyc] sign=+1
    r1 = nyc - str2
    kk = 0
    if (pe1y) then
     if (iby<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    do ic = 1, nc
     do iz = iz1, iz2
      do j = 1, str2
       iy = j + r1
       do ix = ix1, ix2
        kk = kk + 1
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + aux2(kk)
       end do
      end do
     end do
    end do
! sends y=[nyc+1:nyc+str] str=3 data to the right
!receives from left and adds data on y=[j1:j1+str2] sign=-1
    s2 = nyc
    kk = 0
    do ic = 1, nc
     do iz = iz1, iz2
      do j = 1, str
       iy = j + s2
       do ix = ix1, ix2
        kk = kk + 1
        aux1(kk) = curr(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 1, left)
!=====================
    r2 = j1 - 1
    kk = 0
    if (pe0y) then
     if (iby<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    do ic = 1, nc
     do iz = iz1, iz2
      do j = 1, str
       iy = j + r2
       do ix = ix1, ix2
        kk = kk + 1
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + aux2(kk)
       end do
      end do
     end do
    end do
   end if
! The reduced stencil of summed data
   iy1 = j1
   iy2 = nyc
!================
   if (prlz) then
!================
! sends z=[k1-2;k1-1] str-1 data to left
! receives from right and adds data on y=[nzc-1:nzc] sign=+1

    s1 = k1 - str
    kk = 0
    do ic = 1, nc
     do j = 1, str2
      iz = s1 + j
      do iy = iy1, iy2
       do ix = ix1, ix2
        kk = kk + 1
        aux1(kk) = curr(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 2, right)

    r1 = nzc - str2
    if (pe1z) then
     if (ibz<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    kk = 0
    do ic = 1, nc
     do j = 1, str2
      iz = j + r1
      do iy = iy1, iy2
       do ix = ix1, ix2
        kk = kk + 1
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + aux2(kk)
       end do
      end do
     end do
    end do
!
!================
! sends z=[nzc+1:nzc+str] str=3 data to the right
!receives from left and adds data on z=[k1:k1+str2] sign=-1
    s2 = nzc
    kk = 0
    do ic = 1, nc
     do j = 1, str
      iz = j + s2
      do iy = iy1, iy2
       do ix = ix1, ix2
        kk = kk + 1
        aux1(kk) = curr(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 2, left)
!================
! Recvs and adds on z[k1:k1:str2]
!================
    r2 = k1 - 1
    if (pe0z) then
     if (ibz<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    kk = 0
    do ic = 1, nc
     do j = 1, str
      iz = j + r2
      do iy = iy1, iy2
       do ix = ix1, ix2
        kk = kk + 1
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + aux2(kk)
       end do
      end do
     end do
    end do
! The reduced stencil of summed data
    iz1 = k1
    iz2 = nzc
   end if
!========================== prlx case
   if (prlx) then
! sends x=[i1-2;i1-1] str-1 data to left
! receives from right and adds data on x=[nxc-1:nxc] sign=+1
    s1 = i1 - str
    kk = 0
    do ic = 1, nc
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str2
        ix = s1 + j
        kk = kk + 1
        aux1(kk) = curr(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, right)
!=====================
! sends x=[i1-2;i1-1] str-1 data to left
!receives from right and adds data on x=[nxc-1:nxc] sign=+1
    r1 = nxc - str2
    kk = 0
    if (pex1) then
     if (ibx<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    do ic = 1, nc
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str2
        ix = j + r1
        kk = kk + 1
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + aux2(kk)
       end do
      end do
     end do
    end do
! sends x=[nxc+1:nxc+str] str=3 data to the right
!receives from left and adds data on x=[i1:i1+str2] sign=-1
    s2 = nxc
    kk = 0
    do ic = 1, nc
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str
        ix = j + s2
        kk = kk + 1
        aux1(kk) = curr(ix, iy, iz, ic)
       end do
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws
    call exchange_bdx_data(aux1, aux2, lenws, lenwr, 3, left)
!=====================
    r2 = i1 - 1
    kk = 0
    if (pex0) then
     if (ibx<2) then
      aux2(1:lenwr) = 0.0
     end if
    end if
    do ic = 1, nc
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str
        ix = j + r2
        kk = kk + 1
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + aux2(kk)
       end do
      end do
     end do
    end do
    return
   end if
!===================== No MPI x-decomposition
   if (ibx==2) then
! data nxc-1:nxc sums to i1-2:i1-1
! data i1:i1+2 sums to nxc+1:nxc+3
    s1 = i1 - str
    r1 = nxc - str2
    do ic = 1, nc
     do iz = iz1, iz2
      do iy = iy1, iy2
       do j = 1, str2
        curr(r1+j, iy, iz, ic) = curr(r1+j, iy, iz, ic) + &
          curr(s1+j, iy, iz, ic)
       end do
       do j = 1, str
        ix = i1 - 1 + j
        curr(ix, iy, iz, ic) = curr(ix, iy, iz, ic) + &
          curr(nxc+j, iy, iz, ic)
       end do
      end do
     end do
    end do
   end if
!=============================
  end subroutine
!=================
  subroutine traffic_size_eval(sp_loc, xl, xr, pel, per, ibd, ind, &
    npold, nsr, npnew)

   type (species), intent (in) :: sp_loc
   real (dp), intent (in) :: xl, xr
   logical, intent (in) :: pel, per
   integer, intent (in) :: ibd, ind, npold
   integer, intent (inout) :: nsr(4)
   integer, intent (inout) :: npnew
   integer :: n, p, q
   integer :: nl_send, nl_recv, nr_send, nr_recv, cdir
   real (dp) :: xp

   cdir = ind - 1
   if (ind==1) cdir = 3
   p = 0
   q = 0
   nr_recv = 0
   nl_recv = 0
   if (npold>0) then
    do n = 1, npold
     xp = sp_loc%part(n, ind)
     if (xp>xr) p = p + 1
     if (xp<xl) q = q + 1
    end do
   end if
   nr_send = p
   nl_send = q
   call sr_idata(nr_send, nl_recv, cdir, left)
!sends to right nr_send receives from left nl_recv
   call sr_idata(nl_send, nr_recv, cdir, right)
!sends to left nl_send  receives from right nr_recv

   if (ibd<2) then !NOT PERIODIC BD
    if (pel) nl_recv = 0
    if (per) then
     nr_recv = 0
     if (ibd==1) nr_send = 0
    end if
   end if
   npnew = npold + nl_recv + nr_recv - nl_send - nr_send
!================================
!if(pel)nl_send=0
!if(per)nr_send=0
   nsr(1) = nl_send
   nsr(2) = nr_send
   nsr(3) = nl_recv
   nsr(4) = nr_recv
!================================
  end subroutine
!======================================
  subroutine part_prl_exchange(sp_loc, pstore, xl, xr, xlmin, xrmax, &
    pel, per, ibd, dir, ndv, old_np, n_sr, npt)

   type (species), intent (inout) :: sp_loc
   real (dp), intent (inout) :: pstore(:, :)
   real (dp), intent (in) :: xl, xr, xlmin, xrmax
   logical, intent (in) :: pel, per
   integer, intent (in) :: ibd, dir, ndv, old_np, n_sr(4)
   integer, intent (out) :: npt
   integer, allocatable :: left_pind(:), right_pind(:)
   integer :: k, kk, n, p, q, ns, nr, cdir
   integer :: nl_send, nr_send, nl_recv, nr_recv, vxdir
   real (dp) :: xp
!================ dir are cartesian coordinate index (x,y,z)
!================ cdir are mpi-cartesian index (y,z,x)

   nl_send = n_sr(1)
   nr_send = n_sr(2)
   nl_recv = n_sr(3)
   nr_recv = n_sr(4)
   cdir = dir - 1
   if (dir==1) cdir = 3 !for x-direction
   vxdir = dir + ndim
!================== checks memory
   p = ndv*max(nl_send, nr_send)
   if (p>0) then
    if (size(aux1)<p) then
     deallocate (aux1)
     allocate (aux1(p))
     aux1 = 0.0
    end if
   end if
   q = ndv*max(nl_recv, nr_recv)
   if (q>0) then
    if (size(aux2)<q) then
     deallocate (aux2)
     allocate (aux2(q))
     aux2 = 0.0
    end if
   end if
!==================== copy remaining part => ebfp
   p = max(nr_send, 1)
   q = max(nl_send, 1)
   allocate (right_pind(p))
   allocate (left_pind(q))
   right_pind = 0
   left_pind = 0
   p = 0
   q = 0
   npt = 0
   if (ibd==1 .and. dir==1) then !reflecting on the right
    if (per) then
     do n = 1, old_np
      xp = sp_loc%part(n, dir)
      if (xp>xr) then
       sp_loc%part(n, dir) = xr - (xp-xr)
       sp_loc%part(n, vxdir) = -sp_loc%part(n, vxdir)
      end if
     end do
    end if
   end if
   do n = 1, old_np
    xp = sp_loc%part(n, dir)
    if (xp>xr) then
     p = p + 1
     right_pind(p) = n
    else if (xp<xl) then
     q = q + 1
     left_pind(q) = n
    else
     npt = npt + 1
     pstore(npt, 1:ndv) = sp_loc%part(n, 1:ndv)
    end if
   end do
!=======================
   ns = ndv*nr_send
   nr = ndv*nl_recv
   if (ibd<2) then !NON PERIODIC CASE
    if (per) ns = 0
    if (pel) nr = 0
   end if
   if (ns>0) then
    kk = 0
    if (per) then
!sends to the right only for Periodic boundary
     do k = 1, nr_send
      n = right_pind(k)
      loc_pstore(1:ndv) = sp_loc%part(n, 1:ndv)
      loc_pstore(dir) = loc_pstore(dir) + xlmin - xr
      do q = 1, ndv
       kk = kk + 1
       aux1(kk) = loc_pstore(q)
      end do
     end do
    else
     do k = 1, nr_send
      n = right_pind(k)
      loc_pstore(1:ndv) = sp_loc%part(n, 1:ndv)
      do q = 1, ndv
       kk = kk + 1
       aux1(kk) = loc_pstore(q)
      end do
     end do
    end if
   end if
   if (max(ns,nr)>0) call sr_pdata(aux1, aux2, ns, nr, cdir, left)
! sends ns data to the right
   if (nr>0) then !receives nr data from left
    kk = 0
    p = npt
    do n = 1, nl_recv
     p = p + 1
     do q = 1, ndv
      kk = kk + 1
      pstore(p, q) = aux2(kk)
     end do
    end do
    npt = p
   end if
   ns = ndv*nl_send
   nr = ndv*nr_recv
   if (ibd==0) then
    if (pel) ns = 0
    if (per) nr = 0
   end if
   if (ns>0) then
    kk = 0
    if (pel) then !only for periodic case
     do k = 1, nl_send
      n = left_pind(k)
      loc_pstore(1:ndv) = sp_loc%part(n, 1:ndv)
      loc_pstore(dir) = loc_pstore(dir) + xrmax - xl
      do q = 1, ndv
       kk = kk + 1
       aux1(kk) = loc_pstore(q)
      end do
     end do
    else
     do k = 1, nl_send
      n = left_pind(k)
      loc_pstore(1:ndv) = sp_loc%part(n, 1:ndv)
      do q = 1, ndv
       kk = kk + 1
       aux1(kk) = loc_pstore(q)
      end do
     end do
    end if
   end if
   if (max(ns,nr)>0) call sr_pdata(aux1, aux2, ns, nr, cdir, right)
! sends ns data to the left recieves nr data from right
   if (nr>0) then
    p = npt
    kk = 0
    do n = 1, nr_recv
     p = p + 1
     do q = 1, ndv
      kk = kk + 1
      pstore(p, q) = aux2(kk)
     end do
    end do
    npt = p
   end if
   if (allocated(left_pind)) deallocate (left_pind)
   if (allocated(right_pind)) deallocate (right_pind)
  end subroutine
!================
!=============================
  subroutine reset_all_part_dist(loc_sp, pstore, xl, xr, ib, np, ndv, &
    cin, ndm, np_new)
   type (species), intent (inout) :: loc_sp
   real (dp), intent (inout) :: pstore(:, :)
   real (dp), intent (in) :: xl, xr
   integer, intent (in) :: ib, np, ndv, cin, ndm
   integer, intent (out) :: np_new
   real (dp) :: xp, dxp
   integer :: n, p, pout
!-----------------------
   np_new = np
   p = 0
   pout = 0
   if (ib==2) then
    dxp = xr - xl
    do p = 1, np
     xp = loc_sp%part(p, cin)
     if (xp<xl) loc_sp%part(p, cin) = xp + dxp
     xp = loc_sp%part(p, cin)
     if (xp>xr) loc_sp%part(p, cin) = xp - dxp
    end do
    return
   end if
   if (ib==1) then
    do p = 1, np
     xp = loc_sp%part(p, cin)
     if (xp>xr) then
      loc_sp%part(p, cin) = xr - (xp-xr)
      loc_sp%part(p, cin+ndm) = -loc_sp%part(p, cin+ndm)
     end if
    end do
    return
   end if
   do n = 1, np
    xp = loc_sp%part(n, cin)
    if (xp<=xl) p = p + 1
    if (xp>xr) p = p + 1
   end do
   pout = p
   if (pout>0) then
    p = 0
    do n = 1, np
     xp = loc_sp%part(n, cin)
     if (xp>xl .and. xp<=xr) then
      p = p + 1
      pstore(p, 1:ndv) = loc_sp%part(n, 1:ndv)
     end if
    end do
    np_new = p
   end if
  end subroutine
!==============
  subroutine cell_part_dist(moving_wind)
   logical, intent (in) :: moving_wind
   integer :: ic, nspx, n, np, np_new, npout, np_new_allocate, ndv, &
     np_rs
   integer :: n_sr(4)
   real (dp) :: ymm, ymx, lbd_min, rbd_max
   real (dp) :: zmm, zmx
   real (dp) :: xmm, xmx

   ndv = nd2 + 1
!===================================
! In traffic_size_eval() Counts numbers of left-right exchanges
!nsr(1)=nl_send
!nsr(2)=nr_send
!nsr(3)=nl_recv
!nsr(4)=nr_recv
! ==> new particle number np_new= nl_recv-nl_send+ nr_recv-nr_send
!      In part_prl_exchange()    exchanges particle data by mpi_send_recv
!=====================================
   if (.not. moving_wind) then
    ymm = loc_ygrid(imody)%gmin
    ymx = loc_ygrid(imody)%gmax
    lbd_min = loc_ygrid(0)%gmin
    rbd_max = loc_ygrid(npe_yloc-1)%gmax
    if (prly) then
     do ic = 1, nsp_run
      n_sr = 0
      np = loc_npart(imody, imodz, imodx, ic)
      np_new = np
      call traffic_size_eval(spec(ic), ymm, ymx, pe0y, pe1y, iby, 2, np, &
        n_sr, np_new)
      np_new_allocate = max(1, np_new)
      np_rs = maxval(n_sr(1:4))
      if (np_rs>0) then
       call v_realloc(ebfp, np_new, ndv)
       call part_prl_exchange(spec(ic), ebfp, ymm, ymx, lbd_min, &
         rbd_max, pe0y, pe1y, iby, 2, ndv, np, n_sr, npout)
       if (npout/=np_new) then
        write (6, *) 'error in y-part count', mype, npout, np_new
        ier = 99
       end if
       call p_realloc(spec(ic), np_new, ndv)
       do n = 1, np_new
        spec(ic)%part(n, 1:ndv) = ebfp(n, 1:ndv)
       end do
       loc_npart(imody, imodz, imodx, ic) = np_new
      end if
     end do
    else
     do ic = 1, nsp_run
      np = loc_npart(imody, imodz, imodx, ic)
      if (np>0) then
       call reset_all_part_dist(spec(ic), ebfp, ymm, ymx, iby, np, ndv, &
         2, ndim, np_new)
       if (np_new<np) then
        loc_npart(imody, imodz, imodx, ic) = np_new
        do n = 1, np_new
         spec(ic)%part(n, 1:ndv) = ebfp(n, 1:ndv)
        end do
       end if
      end if
     end do
    end if
    if (ndim>2) then
     zmm = loc_zgrid(imodz)%gmin
     zmx = loc_zgrid(imodz)%gmax
     lbd_min = loc_zgrid(0)%gmin
     rbd_max = loc_zgrid(npe_zloc-1)%gmax
     if (prlz) then
      do ic = 1, nsp_run
       np = loc_npart(imody, imodz, imodx, ic)
       np_new = np
       n_sr = 0
       call traffic_size_eval(spec(ic), zmm, zmx, pe0z, pe1z, ibz, 3, &
         np, n_sr, np_new)
       np_new_allocate = max(1, np_new)
       np_rs = maxval(n_sr(1:4))
       if (np_rs>0) then
        call v_realloc(ebfp, np_new, ndv)
        call part_prl_exchange(spec(ic), ebfp, zmm, zmx, lbd_min, &
          rbd_max, pe0z, pe1z, ibz, 3, ndv, np, n_sr, npout)
        if (npout/=np_new) then
         write (6, *) 'error in x-part count', mype, npout, np_new
         ier = 99
        end if
        call p_realloc(spec(ic), np_new, ndv)
        do n = 1, np_new
         spec(ic)%part(n, 1:ndv) = ebfp(n, 1:ndv)
        end do
        loc_npart(imody, imodz, imodx, ic) = np_new
       end if
      end do
     else
      do ic = 1, nsp_run
       np = loc_npart(imody, imodz, imodx, ic)
       if (np>0) then
        call reset_all_part_dist(spec(ic), ebfp, zmm, zmx, ibz, np, ndv, &
          3, ndim, np_new)
        if (np_new<np) then
         loc_npart(imody, imodz, imodx, ic) = np_new
         do n = 1, np_new
          spec(ic)%part(n, 1:ndv) = ebfp(n, 1:ndv)
         end do
        end if
       end if
      end do
     end if
    end if
   end if !end of moving_window=false
!=====================
!In moving window all species leaving the computational box at the left
!x-boundary are removed
!==========================================
   nspx = nsp
   xmm = loc_xgrid(imodx)%gmin
   xmx = loc_xgrid(imodx)%gmax
   lbd_min = loc_xgrid(0)%gmin
   rbd_max = loc_xgrid(npe_xloc-1)%gmax
   if (prlx) then
    do ic = 1, nspx
     np = loc_npart(imody, imodz, imodx, ic)
     np_new = np
     n_sr = 0
     call traffic_size_eval(spec(ic), xmm, xmx, pex0, pex1, ibx, 1, np, &
       n_sr, np_new)
     np_new_allocate = max(1, np_new)
     np_rs = maxval(n_sr(1:4))
     if (np_rs>0) then
      call v_realloc(ebfp, np_new, ndv)
      call part_prl_exchange(spec(ic), ebfp, xmm, xmx, lbd_min, rbd_max, &
        pex0, pex1, ibx, 1, ndv, np, n_sr, npout)
      if (npout/=np_new) then
       write (6, *) 'error in x-part count', mype, npout, np_new
       ier = 99
      end if
      call p_realloc(spec(ic), np_new, ndv)
      do n = 1, np_new
       spec(ic)%part(n, 1:ndv) = ebfp(n, 1:ndv)
      end do
      loc_npart(imody, imodz, imodx, ic) = np_new
     end if
    end do
   else
    do ic = 1, nspx
     np = loc_npart(imody, imodz, imodx, ic)
     if (np>0) then
      call reset_all_part_dist(spec(ic), ebfp, xmm, xmx, ibx, np, ndv, &
        1, ndim, np_new)
      if (np_new<np) then
       loc_npart(imody, imodz, imodx, ic) = np_new
       do n = 1, np_new
        spec(ic)%part(n, 1:ndv) = ebfp(n, 1:ndv)
       end do
      end if
     end if
    end do
   end if
  end subroutine
!=========================
  subroutine cell_bpart_dist(moving_wind)
   logical, intent (in) :: moving_wind
   integer :: ic, n, np, np_new, np_new_allocate, npout, ndv, np_rs
   integer :: n_sr(4)
   real (dp) :: ymm, ymx
   real (dp) :: zmm, zmx
   real (dp) :: xmm, xmx
   real (dp) :: lbd_min, rbd_max

   ndv = nd2 + 1
!====================
   if (.not. moving_wind) then
    ymm = loc_ygrid(imody)%gmin
    ymx = loc_ygrid(imody)%gmax

    lbd_min = loc_ygrid(0)%gmin
    rbd_max = loc_ygrid(npe_yloc-1)%gmax

    if (prly) then
     do ic = 1, nsb
      n_sr = 0
      np = loc_nbpart(imody, imodz, imodx, ic)
      np_new = np
      call traffic_size_eval(bunch(ic), ymm, ymx, pe0y, pe1y, iby, 2, &
        np, n_sr, np_new)
      np_rs = maxval(n_sr(1:4))
      np_new_allocate = max(1, np_new)
      if (np_rs>0) then
       call v_realloc(ebfb, np_new, ndv)
       call part_prl_exchange(bunch(ic), ebfb, ymm, ymx, lbd_min, &
         rbd_max, pe0y, pe1y, iby, 2, ndv, np, n_sr, npout)
       if (npout/=np_new) then
        write (6, *) 'error in y-bpart count', mype, npout, np_new
        ier = 99
       end if
       call p_realloc(bunch(ic), np_new, ndv)
       do n = 1, np_new
        bunch(ic)%part(n, 1:ndv) = ebfb(n, 1:ndv)
       end do
       loc_nbpart(imody, imodz, imodx, ic) = np_new
      end if
     end do
    else
     do ic = 1, nsb
      np = loc_nbpart(imody, imodz, imodx, ic)
      if (np>0) then
       call reset_all_part_dist(bunch(ic), ebfb, ymm, ymx, iby, np, ndv, &
         2, ndim, np_new)
       if (np_new<np) then
        loc_nbpart(imody, imodz, imodx, ic) = np_new
        do n = 1, np_new
         bunch(ic)%part(n, 1:ndv) = ebfb(n, 1:ndv)
        end do
       end if
      end if
     end do
    end if
    zmm = loc_zgrid(imodz)%gmin
    zmx = loc_zgrid(imodz)%gmax
    lbd_min = loc_zgrid(0)%gmin
    rbd_max = loc_zgrid(npe_zloc-1)%gmax
    if (prlz) then
     do ic = 1, nsb
      n_sr = 0
      np = loc_nbpart(imody, imodz, imodx, ic)
      np_new = np
      n_sr = 0
      call traffic_size_eval(bunch(ic), zmm, zmx, pe0z, pe1z, ibz, 3, &
        np, n_sr, np_new)
      np_rs = maxval(n_sr(1:4))
      np_new_allocate = max(1, np_new)
      if (np_rs>0) then
       call v_realloc(ebfb, np_new, ndv)
       call part_prl_exchange(bunch(ic), ebfb, zmm, zmx, lbd_min, &
         rbd_max, pe0z, pe1z, ibz, 3, ndv, np, n_sr, npout)
       if (npout/=np_new) then
        ier = 99
       end if
       call p_realloc(bunch(ic), np_new, ndv)
       do n = 1, np_new
        bunch(ic)%part(n, 1:ndv) = ebfb(n, 1:ndv)
       end do
       loc_nbpart(imody, imodz, imodx, ic) = np_new
      end if
     end do
    else
     if (ndim==3) then
      do ic = 1, nsb
       np = loc_nbpart(imody, imodz, imodx, ic)
       if (np>0) then
        call reset_all_part_dist(bunch(ic), ebfb, zmm, zmx, ibz, np, &
          ndv, 3, ndim, np_new)
        if (np_new<np) then
         loc_nbpart(imody, imodz, imodx, ic) = np_new
         do n = 1, np_new
          bunch(ic)%part(n, 1:ndv) = ebfb(n, 1:ndv)
         end do
        end if
       end if
      end do
     end if
    end if
   end if !end of moving_window=false
!=====================
   xmm = loc_xgrid(imodx)%gmin
   xmx = loc_xgrid(imodx)%gmax
   lbd_min = loc_xgrid(0)%gmin
   rbd_max = loc_xgrid(npe_xloc-1)%gmax
   if (prlx) then
    do ic = 1, nsb
     n_sr = 0
     np = loc_nbpart(imody, imodz, imodx, ic)
     np_new = np
     n_sr = 0
     call traffic_size_eval(bunch(ic), xmm, xmx, pex0, pex1, ibx, 1, np, &
       n_sr, np_new)
     np_rs = maxval(n_sr(1:4))
     np_new_allocate = max(1, np_new)
     if (np_rs>0) then
      call v_realloc(ebfb, np_new, ndv)
      call part_prl_exchange(bunch(ic), ebfb, xmm, xmx, lbd_min, &
        rbd_max, pex0, pex1, ibx, 1, ndv, np, n_sr, npout)
      if (npout/=np_new) then
       ier = 99
      end if
      call p_realloc(bunch(ic), np_new, ndv)
      do n = 1, np_new
       bunch(ic)%part(n, 1:ndv) = ebfb(n, 1:ndv)
      end do
      loc_nbpart(imody, imodz, imodx, ic) = np_new
     end if
    end do
   else
    do ic = 1, nsb
     np = loc_nbpart(imody, imodz, imodx, ic)
     if (np>0) then
      call reset_all_part_dist(bunch(ic), ebfb, xmm, xmx, ibx, np, ndv, &
        1, ndim, np_new)
      if (np_new<np) then
       loc_nbpart(imody, imodz, imodx, ic) = np_new
       do n = 1, np_new
        bunch(ic)%part(n, 1:ndv) = ebfb(n, 1:ndv)
       end do
      end if
     end if
    end do
   end if
  end subroutine
!=====================
 end module
!==================================
