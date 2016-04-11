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

 module pstruct_data
 use precision_def
 use struct_def

 implicit none

 type(species) :: spec(4),bunch(4)
 real(dp),allocatable :: ebfp(:,:),ebfb(:,:),ebfp0(:,:),ebfp1(:,:)
 real(dp),allocatable :: rho_0(:),xpt(:,:),ypt(:,:),zpt(:,:),wghpt(:,:)
 real(dp),allocatable :: loc_ypt(:,:),loc_zpt(:,:),loc_wghy(:,:),loc_wghz(:,:)
 real(dp),allocatable :: loc_xpt(:,:),loc_wghx(:,:)

 !--------------------------

 contains

 !--------------------------
!DIR$ ATTRIBUTES INLINE :: p_realloc
 subroutine p_realloc(pdata,ndv,npt_new)
  real(dp),intent(inout),allocatable  :: pdata(:,:)
  integer,intent(in) :: ndv,npt_new

  if(size(pdata,2) < npt_new) then
   deallocate(pdata)
   allocate(pdata(ndv,npt_new))
  endif
 end subroutine p_realloc
!===========================
!DIR$ ATTRIBUTES INLINE :: p_realloc_sp
 subroutine p_realloc_sp(pdata,npt_new)
  type(species),intent(inout) :: pdata
  integer,intent(in) :: npt_new
  integer :: i,AllocStatus, DeallocStatus

  if(allocated(pdata%part))then
   if(size(pdata%part) < npt_new)then
    deallocate(pdata%part,STAT=DeallocStatus)
    if(DeallocStatus==0)allocate(pdata%part(1:npt_new),STAT=AllocStatus)
   endif
  else
    allocate(pdata%part(1:npt_new),STAT=AllocStatus)
  endif
  do i=1,npt_new
   pdata%part(i)%cmp(:)=0.0
  end do
 end subroutine p_realloc_sp
!===========================
 subroutine p_alloc(npt_max,np_s,ns,lp,mid,ebfp_size,r_type,msize)
 integer,intent(in) :: npt_max,np_s(:),ns,lp,mid,ebfp_size,r_type
 integer,intent(inout) :: msize
 integer :: nsize,ic,npt,ndm_loc,AllocStatus
 integer :: i

 ndm_loc=ebfp_size-1
 select case(r_type)
 case(1)             !Plasma particles
  nsize=0
  do ic=1,ns
   npt=max(np_s(ic),1)
   allocate(spec(ic)%part(npt),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt
   do i=1,npt
    spec(ic)%part(i)%cmp(:)=0.0
   end do
  enddo
  if(mid>0)then
   allocate(ebfp(ebfp_size,npt_max),STAT=AllocStatus)
   nsize=nsize+ebfp_size*npt_max
   ebfp=0.0
  endif
  if(lp>2)then
   ic=1
   npt=max(np_s(ic),1)
   allocate(ebfp0(ndm_loc,npt),STAT=AllocStatus)
   nsize=nsize+ndm_loc*npt
   if(lp ==4)then
    allocate(ebfp1(ndm_loc,npt),STAT=AllocStatus)
    nsize=nsize+ndm_loc*npt
   endif
  endif
  msize=msize+nsize
 case(2)             !bunch particles
  nsize=0
  do ic=1,ns
   npt=np_s(ic)
   if(npt>0)allocate(bunch(ic)%part(npt),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt
  enddo
  if(mid>0)then
   allocate(ebfb(ebfp_size,npt_max),STAT=AllocStatus)
   nsize=nsize+ebfp_size*npt_max
   ebfb=0.0
  endif
  msize=msize+nsize
 end select
 end subroutine p_alloc
 !--------------------------
#if 0
  subroutine p_alloc_sp(npt_max,np_s,ns,lp,mid,r_type,msize)

  integer,intent(in) :: npt_max,np_s(:),ns,lp,mid,r_type
  integer,intent(inout) :: msize
  integer :: nsize,ic,npt,AllocStatus
  integer :: i

 select case(r_type)
 case(1)             !Plasma particles
  nsize=0
  do ic=1,ns
   npt=max(np_s(ic),1)
   allocate(spec(ic)%part(npt),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt
   do i=1,npt
    spec(ic)%part(i)%cmp(:)=0.0
   end do
  enddo
  if(mid>0)then
   allocate(ebfp%part(npt_max),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt_max
   do i=1,npt_max
    ebfp%part(i)%cmp(:)=0.0
   end do
  endif
  if(lp>2)then
   allocate(ebfp0%part(npt),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt
   do i=1,npt
    ebfp0%part(i)%cmp(1:P_ncmp)=0.0
   end do
   if(lp ==4)then
    allocate(ebfp1%part(npt),STAT=AllocStatus)
    nsize=nsize+P_ncmp*npt
    do i=1,npt
     ebfp1%part(i)%cmp(1:P_ncmp)=0.0
    end do
   endif
  endif
  msize=msize+nsize
 case(2)             !bunch particles
  nsize=0
  do ic=1,ns
   npt=np_s(ic)
   if(npt>0)allocate(bunch(ic)%part(npt),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt
   do i=1,npt
    bunch(ic)%part(i)%cmp(:)=0.0
   end do
  enddo
  if(mid>0)then
   allocate(ebfb%part(npt_max),STAT=AllocStatus)
   nsize=nsize+P_ncmp*npt_max
   do i=1,npt_max
    ebfb%part(i)%cmp(1:P_ncmp)=0.0
   end do
  endif
  msize=msize+nsize
 end select
 end subroutine p_alloc_sp
#endif
 !--------------------------
 end module pstruct_data

