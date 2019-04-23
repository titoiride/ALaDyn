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
 module struct_def
  use precision_def
  implicit none

 type species
  real(dp),allocatable :: part(:,:)
 end type species


 type grid
  integer :: ng,p_ind(2)
  real(dp) :: gmin,gmax
 end type grid

 type sgrid
  integer :: sind(2)
  real(dp) :: smin,smax
 end type sgrid

 end module struct_def
!=======================
 module fstruct_data

 use precision_def

 implicit none

 real(dp),allocatable :: ebf(:,:,:,:),ebf_bunch(:,:,:,:),jc(:,:,:,:)
 real(dp),allocatable :: ebf0(:,:,:,:),ebf1(:,:,:,:)
 real(dp),allocatable :: ebf0_bunch(:,:,:,:),ebf1_bunch(:,:,:,:),jb(:,:,:,:)
 real(dp),allocatable :: env(:,:,:,:),env0(:,:,:,:),env1(:,:,:,:)
 real(dp),allocatable :: up(:,:,:,:),up0(:,:,:,:),up1(:,:,:,:),flux(:,:,:,:)
 real(dp),allocatable :: pot(:,:,:,:),fluid_x_profile(:),fluid_yz_profile(:,:)
 real(dp),allocatable :: aux1(:),aux2(:)
 end module fstruct_data
 !--------------------------
 module pstruct_data

 use precision_def

 use struct_def

 implicit none

 real(dp),allocatable :: ebfp(:,:),ebfb(:,:)
 real(dp),allocatable :: ebfp0(:,:),ebfp1(:,:)
 real(dp),allocatable :: pdata_tracking(:,:,:)
 real(dp),allocatable :: track_aux(:)
!
 real(dp),allocatable :: xpt(:,:),ypt(:,:),zpt(:,:),wghpt(:,:)
 real(dp),allocatable :: loc_ypt(:,:),loc_zpt(:,:),loc_wghyz(:,:,:)
 real(dp),allocatable :: loc_xpt(:,:),loc_wghx(:,:)
 type(species) :: spec(4),bunch(5)
 integer(hp_int),parameter :: ihx=3

 !=====================
 end module pstruct_data
!===================
 module array_wspace

  use pstruct_data
  use fstruct_data

 end module array_wspace

