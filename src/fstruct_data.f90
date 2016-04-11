 !*****************************************************************************************************!
 !             Copyright 2008-2016 Pasquale Londrillo, Stefano Sinigardi, Andrea Sgattoni              !
 !                                 Alberto Marocchino                                                  !
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

 module fstruct_data

 use precision_def

 implicit none

 real(dp),allocatable :: ebf_bunch_gammarange(:,:,:,:)
 real(dp),allocatable :: ebf1_bunch_gammarange(:,:,:,:)
 real(dp),allocatable :: ebf_gammarange(:,:,:,:)
 real(dp),allocatable :: jc_gammarange(:,:,:,:)

 real(dp),allocatable :: ebf(:,:,:,:),ebf_bunch(:,:,:,:),jc(:,:,:,:)
 real(dp),allocatable :: pot(:,:,:,:)
 real(dp),allocatable :: ebf0(:,:,:,:),ebf1(:,:,:,:)
 real(dp),allocatable :: ebf0_bunch(:,:,:,:),ebf1_bunch(:,:,:,:)
 integer,allocatable :: sp_count(:,:)
 real(dp),allocatable :: env(:,:,:,:),env0(:,:,:,:),env1(:,:,:,:),jb0(:,:,:,:)
 real(dp),allocatable :: up(:,:,:,:),up0(:,:,:,:),up1(:,:,:,:)
 real(dp),allocatable :: ub(:,:,:,:),ub0(:,:,:,:),ub1(:,:,:,:)
 !--------------------------

 contains

 !--------------------------

 subroutine v_alloc(n1,n2,n3,ncomp,njc,ndm,ns,ipot,lp,envlp,fsize)

 integer,intent(in) ::n1,n2,n3,ncomp,njc,ndm,ns,lp,ipot
 logical,intent(in) :: envlp
 integer,intent(inout) ::fsize
 integer :: njdim,ng,ng0,n1p,n2p,n3p,shx,AllocStatus,fsize_loc
 !==============
 ! ns =nsp for active ionization
 !======================
 shx=3
 n1p=n1+shx       !extended grid [1:n1p=n1+3=nx+5 ]
 n2p=n2+shx
 n3p=n3
 ng0=1+(n1-2)*(n2-2)
 if(ndm==3)then
  ng0=1+(n1-2)*(n2-2)*(n3-2)
  n3p=n3+shx
 endif
 ng=n1p*n2p*n3p
 fsize_loc=0

 ! allocates common arrays
 allocate(ebf(n1p,n2p,n3p,ncomp),STAT=AllocStatus)
 allocate(jc(n1p,n2p,n3p,njc),STAT=AllocStatus)

 if(ns>0)then
  allocate(sp_count(ng0,ns),STAT=AllocStatus)
  sp_count=0
 endif
 fsize_loc=fsize_loc+ng*ncomp+ng*njc+ns*ng0
 ebf=0.0
 jc=0.0
 if(ipot> 1)then
  allocate(pot(n1p,n2p,n3p,2),STAT=AllocStatus)
  pot=0.0
  fsize_loc=fsize_loc+2*ng
 endif
 if(lp>2)then
  allocate(ebf0(n1p,n2p,n3p,ncomp),STAT=AllocStatus)
  ebf0=0.0
  fsize_loc=fsize_loc+ng*ncomp
  if(lp >3)then
   allocate(ebf1(n1p,n2p,n3p,ncomp),STAT=AllocStatus)
   fsize_loc=fsize_loc+ng*ncomp
   ebf1=0.0
  endif
 endif
 if(envlp)then
  njdim=2
  allocate(env(n1p,n2p,n3p,njdim),STAT=AllocStatus)
  allocate(env0(n1p,n2p,n3p,njdim),STAT=AllocStatus)
  env=0.0
  env0=0.0
  fsize_loc=fsize_loc+2*njdim*ng
  if(lp >3)then
   allocate(env1(n1p,n2p,n3p,njdim),STAT=AllocStatus)
   env1=0.0
   fsize_loc=fsize_loc+njdim*ng
  endif
 endif
 fsize=fsize+fsize_loc  !sums all over memory alloc
 end subroutine v_alloc

 !--------------------------

 subroutine bv_alloc(n1,n2,n3,bcomp,ndm,ibch,fsize)

 integer,intent(in) ::n1,n2,n3,bcomp,ndm,ibch
 integer,intent(inout) ::fsize
 integer :: ng,n1p,n2p,n3p,shx,AllocStatus

 shx=3
 n1p=n1+shx       !x-grid ix=1,2 bd, 3:nx+2=n1p data n1p+1 bd
 n2p=n2+shx
 n3p=n3
 if(ndm==3)n3p=n3+shx
 ng=n1p*n2p*n3p
 allocate(ebf_bunch(n1p,n2p,n3p,bcomp),STAT=AllocStatus)
 ebf_bunch=0.0

 allocate(ebf_bunch_gammarange(n1p,n2p,n3p,bcomp),STAT=AllocStatus)
 allocate(ebf1_bunch_gammarange(n1p,n2p,n3p,bcomp),STAT=AllocStatus)
 allocate(ebf_gammarange(n1p,n2p,n3p,bcomp),STAT=AllocStatus)
 allocate(jc_gammarange(n1p,n2p,n3p,bcomp),STAT=AllocStatus)

 fsize=fsize+ng*bcomp
 if(ibch>0)then
  allocate(ebf1_bunch(n1p,n2p,n3p,bcomp),STAT=AllocStatus)
  ebf1_bunch=0.0
  fsize=fsize+ng*bcomp
 endif
 end subroutine bv_alloc
 !==================
 subroutine pbv_alloc(n1,n2,n3,bcomp,ndm,fsize)

 integer,intent(in) ::n1,n2,n3,bcomp,ndm
 integer,intent(inout) ::fsize
 integer :: ng,n1p,n2p,n3p,shx,AllocStatus

 shx=3
 n1p=n1+shx       !x-grid ix=1,2 bd, 3:nx+2=n1p data n1p+1 bd
 n2p=n2+shx
 n3p=n3
 if(ndm==3)n3p=n3+shx
 ng=n1p*n2p*n3p
 allocate(ebf_bunch(n1p,n2p,n3p,bcomp),STAT=AllocStatus)
 ebf_bunch=0.0
 fsize=fsize+ng*bcomp
!============ external B-field allocated
 allocate(ebf0_bunch(n1p,n2p,n3p,3),STAT=AllocStatus)
 ebf0_bunch=0.0
 fsize=fsize+3*ng
 end subroutine pbv_alloc

 !--------------------------

 subroutine v_dalloc

 if(allocated(ebf))deallocate(ebf)
 if(allocated(ebf0))deallocate(ebf0)
 if(allocated(ebf1))deallocate(ebf1)
 if(allocated(jc))deallocate(jc)

 if(allocated(ebf_bunch_gammarange)) deallocate(ebf_bunch_gammarange)
 if(allocated(ebf1_bunch_gammarange)) deallocate(ebf1_bunch_gammarange)
 if(allocated(ebf_gammarange)) deallocate(ebf_gammarange)
 if(allocated(jc_gammarange)) deallocate(jc_gammarange)

 end subroutine v_dalloc

 !--------------------------

 end module fstruct_data
