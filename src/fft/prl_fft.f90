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

 module prl_fft
  use parallel
  use fft_lib
  implicit none

  real (dp), allocatable :: fp1(:, :, :), fp2(:, :, :), faux1(:), &
    faux2(:)


 contains
  subroutine mpi_ftw_alloc(n1, n2, n2_loc, n3, n3_loc)
   integer, intent (in) :: n1, n2, n2_loc, n3, n3_loc
   integer :: n1_xloc, n2_xloc, lenw
   !=======================
   ! loc_grid allready define


   n1_xloc = n1/npe_yloc
   n2_xloc = n1/npe_zloc
   lenw = n1_xloc*n2_loc*n3_loc
   allocate (faux1(lenw), faux2(lenw))
   allocate (fp1(n1_xloc,n2,n3_loc))
   allocate (fp2(n2_xloc,n2_loc,n3))

  end subroutine
  !==========================
  subroutine mpi_ftw_dalloc

   if (allocated(faux1)) deallocate (faux1)
   if (allocated(faux2)) deallocate (faux2)
   if (allocated(fp1)) deallocate (fp1)
   if (allocated(fp2)) deallocate (fp2)

  end subroutine

  !===================
  subroutine swap_yx_3data(waux, wdata, n1_loc, n2, n3)

   integer, intent (in) :: n1_loc, n2, n3
   real (dp), intent (in) :: waux(:, :, :)
   real (dp), intent (out) :: wdata(:, :, :)

   integer :: pes, per, ip
   integer :: i1, j1, lenws, lenwr, tag, n2_xloc, iy, ix, iz
   integer :: kk
   !-----------------
   !From waux(1:n1,1:n2_xloc) to wdata(1:n1_loc,1:n2)
   n2_xloc = n2/npe_xloc

   do iz = 1, n3
    do iy = 1, n2_xloc
     j1 = iy + n2_xloc*imodx
     do ix = 1, n1_loc
      i1 = ix + n1_loc*imodx
      wdata(ix, j1, iz) = waux(i1, iy, iz)
     end do
    end do
   end do
   if (npe_yloc==1) return
   lenws = n1_loc*n2_xloc*n3
   lenwr = lenws
   if (size(faux1)<lenwr) then
    deallocate (faux1, faux2)
    allocate (faux1(lenws))
    allocate (faux2(lenwr))
   end if

   do ip = 1, npe_xloc - 1
    pes = xp_next(ip)
    i1 = n1_loc*pes
    per = xp_prev(ip)
    tag = 80 + ip
    kk = 0
    do iz = 1, n3
     do iy = 1, n2_xloc
      do ix = 1, n1_loc
       kk = kk + 1
       faux1(kk) = waux(ix+i1, iy, iz)
      end do
     end do
    end do
    call mpi_sendrecv(faux1(1), lenws, mpi_sd, pes, tag, faux2(1), &
      lenwr, mpi_sd, per, tag, comm_col(3), status, error)
    j1 = n2_xloc*per
    kk = 0
    do iz = 1, n3
     do iy = 1, n2_xloc
      do ix = 1, n1_loc
       kk = kk + 1
       wdata(ix, iy+j1, iz) = faux2(kk)
      end do
     end do
    end do
   end do

  end subroutine

  subroutine swap_xy_3data(wp1, wp2, n1_loc, n2_loc, n3_loc)

   integer, intent (in) :: n1_loc, n2_loc, n3_loc
   real (dp), intent (in) :: wp1(:, :, :)
   real (dp), intent (out) :: wp2(:, :, :)

   integer :: pes, per, ip
   integer :: i1, j1, lenws, lenwr, tag, iy, ix, iz
   integer :: kk
   !-----------------
   !From wp1(1:n1_loc,1:n2) to wp2(1:n1,1:n2_loc,ic)

   do iz = 1, n3_loc
    do iy = 1, n2_loc
     j1 = iy + n2_loc*imody
     do ix = 1, n1_loc
      i1 = ix + n1_loc*imody
      wp2(i1, iy, iz) = wp1(ix, j1, iz)
     end do
    end do
   end do
   if (npe_yloc==1) return
   lenws = n1_loc*n2_loc*n3_loc
   lenwr = lenws
   if (size(faux1)<lenwr) then
    deallocate (faux1, faux2)
    allocate (faux1(lenws))
    allocate (faux2(lenwr))
   end if

   do ip = 1, npe_yloc - 1
    pes = yp_next(ip)
    j1 = n2_loc*pes
    per = yp_prev(ip)
    tag = 80 + ip
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       faux1(kk) = wp1(ix, iy+j1, iz)
      end do
     end do
    end do
    call mpi_sendrecv(faux1(1), lenws, mpi_sd, pes, tag, faux2(1), &
      lenwr, mpi_sd, per, tag, comm_col(1), status, error)
    i1 = n1_loc*per
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       wp2(ix+i1, iy, iz) = faux2(kk)
      end do
     end do
    end do
   end do

  end subroutine

  subroutine swap_xz_3data(wp1, wp2, n1_loc, n2_loc, n3_loc)

   integer, intent (in) :: n1_loc, n2_loc, n3_loc
   real (dp), intent (in) :: wp1(:, :, :)
   real (dp), intent (out) :: wp2(:, :, :)

   integer :: pes, per, ip
   integer :: i1, lenw, tag, iy, ix, iz, iz1
   integer :: kk
   !-----------------
   !From wp1(1:n1_loc,1:n3) to wp2(1:n1,1:n3_loc)
   do iz = 1, n3_loc
    iz1 = iz + n3_loc*imodz
    do iy = 1, n2_loc
     do ix = 1, n1_loc
      i1 = ix + n1_loc*imodz
      wp2(i1, iy, iz) = wp1(ix, iy, iz1)
     end do
    end do
   end do
   if (npe_zloc==1) return
   lenw = n1_loc*n2_loc*n3_loc
   if (size(faux1)<lenw) then
    deallocate (faux1, faux2)
    allocate (faux1(lenw))
    allocate (faux2(lenw))
   end if

   do ip = 1, npe_zloc - 1
    pes = zp_next(ip)
    iz1 = n3_loc*pes
    per = zp_prev(ip)
    tag = 80 + ip
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       faux1(kk) = wp1(ix, iy, iz+iz1)
      end do
     end do
    end do
    call mpi_sendrecv(faux1(1), lenw, mpi_sd, pes, tag, faux2(1), lenw, &
      mpi_sd, per, tag, comm_col(2), status, error)
    i1 = n1_loc*per
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       wp2(ix+i1, iy, iz) = faux2(kk)
      end do
     end do
    end do
   end do

  end subroutine

  subroutine swap_yx_3data_inv(wdata, waux, n1_loc, n2, n3)

   integer, intent (in) :: n1_loc, n2, n3
   real (dp), intent (in) :: wdata(:, :, :)
   real (dp), intent (out) :: waux(:, :, :)

   integer :: pes, per, ip
   integer :: i1, j1, lenws, lenwr, tag, iy, ix, iz, n2_xloc
   integer :: kk
   !-----------------
   !enters wdata(n1_loc,n2,n3)
   !From wdata(1:n1_loc,1:n2) to waux(1:n1,1:n2_xloc)
   !n2_xloc=n2/npe_xloc

   n2_xloc = n2/npe_xloc
   kk = 0
   do iz = 1, n3
    do iy = 1, n2_xloc
     j1 = iy + n2_xloc*imodx
     do ix = 1, n1_loc
      i1 = ix + n1_loc*imodx
      kk = kk + 1
      waux(i1, iy, iz) = wdata(ix, j1, iz)
     end do
    end do
   end do
   if (npe_xloc==1) return
   lenwr = n1_loc*n2_xloc*n3
   lenws = lenwr
   if (size(faux1)<lenwr) then
    deallocate (faux1, faux2)
    allocate (faux1(lenws))
    allocate (faux2(lenwr))
   end if
   do ip = 1, npe_xloc - 1
    pes = xp_next(ip)
    j1 = n2_xloc*pes
    per = xp_prev(ip)
    tag = 40 + ip
    kk = 0
    do iz = 1, n3
     do iy = 1, n2_xloc
      do ix = 1, n1_loc
       kk = kk + 1
       faux1(kk) = wdata(ix, iy+j1, iz)
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws

    call mpi_sendrecv(faux1(1), lenws, mpi_sd, pes, tag, faux2(1), &
      lenwr, mpi_sd, per, tag, comm_col(3), status, error)
    i1 = n1_loc*per
    kk = 0
    do iz = 1, n3
     do iy = 1, n2_xloc
      do ix = 1, n1_loc
       kk = kk + 1
       waux(ix+i1, iy, iz) = faux2(kk)
      end do
     end do
    end do
   end do
  end subroutine

  subroutine swap_xy_3data_inv(wp2, wp1, n1_loc, n2_loc, n3_loc)

   integer, intent (in) :: n1_loc, n2_loc, n3_loc
   real (dp), intent (in) :: wp2(:, :, :)
   real (dp), intent (out) :: wp1(:, :, :)

   integer :: pes, per, ip
   integer :: i1, j1, lenws, lenwr, tag, iy, ix, iz
   integer :: kk
   !-----------------
   !From (1:n1,1:n2_loc) to (1:n1_loc,1:n2)

   kk = 0
   do iz = 1, n3_loc
    do iy = 1, n2_loc
     j1 = iy + n2_loc*imody
     do ix = 1, n1_loc
      i1 = ix + n1_loc*imody
      kk = kk + 1
      wp1(ix, j1, iz) = wp2(i1, iy, iz)
     end do
    end do
   end do
   if (npe_yloc==1) return
   lenwr = n1_loc*n2_loc*n3_loc
   lenws = lenwr
   if (size(faux1)<lenwr) then
    deallocate (faux1, faux2)
    allocate (faux1(lenws))
    allocate (faux2(lenwr))
   end if
   do ip = 1, npe_yloc - 1
    pes = yp_next(ip)
    i1 = n1_loc*pes
    per = yp_prev(ip)
    tag = 40 + ip
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       faux1(kk) = wp2(ix+i1, iy, iz)
      end do
     end do
    end do
    lenws = kk
    lenwr = lenws

    call mpi_sendrecv(faux1(1), lenws, mpi_sd, pes, tag, faux2(1), &
      lenwr, mpi_sd, per, tag, comm_col(1), status, error)
    j1 = n2_loc*per
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       wp1(ix, iy+j1, iz) = faux2(kk)
      end do
     end do
    end do
   end do
  end subroutine
  !=====================
  subroutine swap_xz_3data_inv(wp2, wp1, n1_loc, n2_loc, n3_loc)

   integer, intent (in) :: n1_loc, n2_loc, n3_loc
   real (dp), intent (in) :: wp2(:, :, :)
   real (dp), intent (out) :: wp1(:, :, :)

   integer :: pes, per, ip
   integer :: i1, lenw, tag, iy, ix, iz, iz1
   integer :: kk, k1
   !-----------------
   !From (1:n1,1:n3_loc) to (1:n1_loc,1:n3)

   do iz = 1, n3_loc
    k1 = iz + imodz*n3_loc
    do iy = 1, n2_loc
     do ix = 1, n1_loc
      i1 = ix + n1_loc*imodz
      wp1(ix, iy, k1) = wp2(i1, iy, iz)
     end do
    end do
   end do
   if (npe_zloc==1) return
   lenw = n1_loc*n2_loc*n3_loc
   if (size(faux1)<lenw) then
    deallocate (faux1, faux2)
    allocate (faux1(lenw))
    allocate (faux2(lenw))
   end if
   do ip = 1, npe_zloc - 1
    pes = zp_next(ip)
    i1 = n1_loc*pes
    per = zp_prev(ip)
    tag = 40 + ip
    kk = 0
    do iz = 1, n3_loc
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       faux1(kk) = wp2(ix+i1, iy, iz)
      end do
     end do
    end do

    call mpi_sendrecv(faux1(1), lenw, mpi_sd, pes, tag, faux2(1), lenw, &
      mpi_sd, per, tag, comm_col(2), status, error)
    k1 = n3_loc*per
    kk = 0
    do iz = 1, n3_loc
     iz1 = iz + k1
     do iy = 1, n2_loc
      do ix = 1, n1_loc
       kk = kk + 1
       wp1(ix, iy, iz1) = faux2(kk)
      end do
     end do
    end do
   end do
  end subroutine
  !====================================
  subroutine pftw2d_sc(w, n1, n2, n2_loc, n3, n3_loc, is, sym)
   real (dp), intent (inout) :: w(:, :, :)
   integer, intent (in) :: n1, n2, n2_loc, n3, n3_loc, is, sym
   integer :: n1_loc
    !performs a 2D FFT sin/cosine`on the (y,z) coordinates for a 3D data (x,y,z)
    !sym=1 for a sine transform
    !sym=2 for a cosine transform

   select case (is)
   case (-1) !
    if (n2<=2) return
    if (prly) then
     n1_loc = n1/npe_yloc
     call swap_xy_3data_inv(w, fp1, n1_loc, n2_loc, n3_loc)

     call ftw1d_sc(fp1, n1_loc, n2, n3_loc, is, 2, sym)
     call swap_xy_3data(fp1, w, n1_loc, n2_loc, n3_loc)
    else
     call ftw1d_sc(w, n1, n2, n3_loc, is, 2, sym)
    end if
    !=====================
    if (n3<=2) return
    if (npe_zloc>1) then
     n1_loc = n1/npe_zloc
     call swap_xz_3data_inv(w, fp2, n1_loc, n2_loc, n3_loc)
     call ftw1d_sc(fp2, n1_loc, n2_loc, n3, is, 3, sym)
     call swap_xz_3data(fp2, w, n1_loc, n2_loc, n3_loc)
    else
     call ftw1d_sc(w, n1, n2_loc, n3, is, 3, sym)
    end if
    !======== exit w(loc)
   case (1)
    ! enters w(loc)
    !========================
    if (n3>1) then
     if (npe_zloc>1) then
      n1_loc = n1/npe_zloc
      call swap_xz_3data_inv(w, fp2, n1_loc, n2_loc, n3_loc)
      call ftw1d_sc(fp2, n1_loc, n2_loc, n3, is, 3, sym)
      call swap_xz_3data(fp2, w, n1_loc, n2_loc, n3_loc)
     else
      call ftw1d_sc(w, n1, n2_loc, n3, is, 3, sym)
     end if
    end if
    !=================
    if (n2>1) then
     if (prly) then
      n1_loc = n1/npe_yloc
      call swap_xy_3data_inv(w, fp1, n1_loc, n2_loc, n3_loc)
      call ftw1d_sc(fp1, n1_loc, n2, n3_loc, is, 2, sym)
      call swap_xy_3data(fp1, w, n1_loc, n2_loc, n3_loc)
     else
      call ftw1d_sc(w, n1, n2, n3_loc, is, 2, sym)
     end if
    end if
   end select
   !===================
   !exit w(loc)
  end subroutine
  !========================
  subroutine pftw3d_sc(w, n1, n2, n2_loc, n3, n3_loc, is, sym)
   real (dp), intent (inout) :: w(:, :, :)
   integer, intent (in) :: n1, n2, n2_loc, n3, n3_loc, is, sym
   integer :: n1_loc

   select case (is)
   case (-1)
    call ftw1d_sc(w, n1, n2_loc, n3_loc, is, 1, sym)
    call pftw2d_sc(w, n1, n2, n2_loc, n3, n3_loc, is, sym)
   case (1)
    ! enters w(loc)
    call pftw2d_sc(w, n1, n2, n2_loc, n3, n3_loc, is, sym)
!   ========================
    call ftw1d_sc(w, n1, n2_loc, n3_loc, is, 1, sym)
   end select
   !===================
   !exit w(loc)
  end subroutine
  !============================
  subroutine pftw2d(w, n1, n2, n2_loc, n3, n3_loc, is)
   real (dp), intent (inout) :: w(:, :, :)
   integer, intent (in) :: n1, n2, n2_loc, n3, n3_loc, is
   integer :: n1_loc, dir

   select case (is)
   case (-1)
    if (n2<=2) return
    dir = 2
    if (prly) then
     n1_loc = n1/npe_yloc
     call swap_xy_3data_inv(w, fp1, n1_loc, n2_loc, n3_loc)
     call ftw1d(fp1, n1_loc, n2, n3_loc, is, dir)
     call swap_xy_3data(fp1, w, n1_loc, n2_loc, n3_loc)
    else
     call ftw1d(w, n1, n2, n3_loc, is, dir)
    end if
    !=====================
    if (n3<=2) return
    dir = 3
    if (npe_zloc>1) then
     n1_loc = n1/npe_zloc
     call swap_xz_3data_inv(w, fp2, n1_loc, n2_loc, n3_loc)
     call ftw1d(fp2, n1_loc, n2_loc, n3, is, dir)
     call swap_xz_3data(fp2, w, n1_loc, n2_loc, n3_loc)
    else
     call ftw1d(w, n1, n2_loc, n3, is, dir)
    end if
    !======== exit w(loc)
   case (1)
    ! enters w(loc)
    !========================
    if (n3>1) then
     dir = 3
     if (npe_zloc>1) then
      n1_loc = n1/npe_zloc
      call swap_xz_3data_inv(w, fp2, n1_loc, n2_loc, n3_loc)
      call ftw1d(fp2, n1_loc, n2_loc, n3, is, dir)
      call swap_xz_3data(fp2, w, n1_loc, n2_loc, n3_loc)
     else
      call ftw1d(w, n1, n2_loc, n3, is, dir)
     end if
    end if
    !=================
    if (n2>1) then
     dir = 2
     if (prly) then
      n1_loc = n1/npe_yloc
      call swap_xy_3data_inv(w, fp1, n1_loc, n2_loc, n3_loc)
      call ftw1d(fp1, n1_loc, n2, n3_loc, is, dir)
      call swap_xy_3data(fp1, w, n1_loc, n2_loc, n3_loc)
     else
      call ftw1d(w, n1, n2, n3_loc, is, dir)
     end if
    end if
   end select
   !===================
   !exit w(loc)
  end subroutine
  !====================================
  subroutine pftw3d(w, n1, n2, n2_loc, n3, n3_loc, is)
   real (dp), intent (inout) :: w(:, :, :)
   integer, intent (in) :: n1, n2, n2_loc, n3, n3_loc, is
   integer :: n1_loc

   select case (is)
   case (-1)

    call ftw1d(w, n1, n2_loc, n3_loc, is, 1)

    call pftw2d(w, n1, n2, n2_loc, n3, n3_loc, is)

   case (1)
    call pftw2d(w, n1, n2, n2_loc, n3, n3_loc, is)

    call ftw1d(w, n1, n2_loc, n3_loc, is, 1)
   end select
   !===================
   !exit w(loc)
  end subroutine
 end module