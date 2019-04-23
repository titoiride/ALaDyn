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


 module mpi_curr_interface

 use array_wspace
 use parallel

 implicit none

 integer(hp_int),parameter :: rt=1,lt=-1

 contains
!===============================
 subroutine jc_xyzbd(curr,i1,n1,j1,n2,k1,n3,nc)
 real(dp),intent(inout) :: curr(:,:,:,:)
 integer,intent(in) :: i1,n1,j1,n2,k1,n3,nc
 integer :: ix,iy,iz,i0,j2,k2,ik
 ! Enter current data on extended ranges:
 !========== Only for Periodic BDs period=n1-1
 j2=n2;k2=n3
 if(ibx==0)then
  do ik=1,nc
   curr(i1,j1:j2,k1:k2,ik)=curr(i1,j1:j2,k1:k2,ik)+curr(i1-1,j1:j2,k1:k2,ik)
   curr(i1-1,j1:j2,k1:k2,ik)=0.0

   curr(n1,j1:j2,k1:k2,ik)=curr(n1,j1:j2,k1:k2,ik)+ &
    curr(n1+1,j1:j2,k1:k2,ik)+ &
    curr(n1+2,j1:j2,k1:k2,ik)
   curr(n1+1:n1+2,j1:j2,k1:k2,ik)=0.0
  end do
  i0=1
 endif
 if(ndim < 2)return
 if(pe0y)then
  if(iby==0)then
   do ik=1,nc
    do iz=k1,k2
     do ix=i1,n1
      curr(ix,j1,iz,ik)=curr(ix,j1,iz,ik)+ &
       curr(ix,j1-1,iz,ik)+curr(ix,j1-2,iz,ik)
       curr(ix,j1-2:j1-1,iz,ik)=0.0
     end do
    end do
   end do
  endif
  if(iby==1)then
   !Before norm  r*Jx(j)=rVxrho ==> rEx odd Jx(j-1)=-Jx(j+1)
   !             r*Jr=rVrrho ==> rEr even   Jr(j-1/2)=jr(j+1/2)
   do iz=k1,k2
    do ix=i1,n1
     curr(ix,j1+1,iz,1)=curr(ix,j1+1,iz,1)+ &
      curr(ix,j1-1,iz,1)
     curr(ix,j1,iz,2)=curr(ix,j1,iz,2)- &
      curr(ix,j1-1,iz,2)
    end do
   end do
  endif
 endif
 if(pe1y)then
  if(iby ==0)then
   do ik=1,nc
    do iz=k1,k2
     do ix=i1,n1
      curr(ix,n2,iz,ik)=curr(ix,n2,iz,ik)+ &
       curr(ix,n2+1,iz,ik)+curr(ix,n2+2,iz,ik)
      curr(ix,n2+1:n2+2,iz,ik)=0.0
     end do
    end do
   end do
  endif
 endif
 if(ndim < 3)return
 if(ibz==0)then
  if(pe0z)then
   do ik=1,nc
    do iy=j1,j2
     do ix=i1,n1
      curr(ix,iy,k1,ik)=curr(ix,iy,k1,ik)+ &
       curr(ix,iy,k1-1,ik)
      curr(ix,iy,k1-1,ik)=0.0
     end do
    end do
   enddo
  endif
  if(pe1z)then
   if(ibz ==0)then
    do ik=1,nc
     do iy=j1,j2
      do ix=i1,n1
       curr(ix,iy,n3,ik)=curr(ix,iy,n3,ik)+ &
        curr(ix,iy,n3+1,ik)+curr(ix,iy,n3+2,ik)
        curr(ix,iy,n3+1:n3+2,ik)=0.0
      end do
     end do
    enddo
   endif
  endif
 endif
 end subroutine jc_xyzbd
 !-----------------------------------------------
 !
 subroutine den_zyxbd(rho,i1,i2,j1,j2,k1,k2,ik)
 real(dp),intent(inout) :: rho(:,:,:,:)
 integer,intent(in) :: i1,i2,j1,j2,k1,k2,ik
 integer :: ix,iy
 ! Enter current data on extended ranges:

 !Enter data on the computational box [i1:n1p][j1:nyp][k1:nzp]
 if(ndim>2)then
  if(ibz==0)then
   if(pe0z)then
    do iy=j1,j2
     do ix=i1,i2
      rho(ix,iy,k1+1,ik)=rho(ix,iy,k1+1,ik)+ &
       rho(ix,iy,k1-1,ik)
      rho(ix,iy,k1,ik)=rho(ix,iy,k1+1,ik)
     end do
    end do
   endif
   if(pe1z)then
    if(ibz <2)then
     do iy=j1,j2
      do ix=i1,i2
       rho(ix,iy,k2-1,ik)=rho(ix,iy,k2-1,ik)+ &
                          rho(ix,iy,k2+1,ik)
       rho(ix,iy,k2,ik)=rho(ix,iy,k2-1,ik)
      end do
     end do
    endif
   endif
  endif
 endif
 !================
 if(ndim>1)then
  if(pe0y)then
   if(iby==0)then
    rho(i1:i2,j1+1,k1:k2,ik)=rho(i1:i2,j1+1,k1:k2,ik)+ &
     rho(i1:i2,j1-1,k1:k2,ik)
    rho(i1:i2,j1,k1:k2,ik)=rho(i1:i2,j1+1,k1:k2,ik)
   endif
   if(iby==1)then
    rho(i1:i2,j1+1,k1:k2,ik)=rho(i1:i2,j1+1,k1:k2,ik)+ &
     rho(i1:i2,j1-1,k1:k2,ik)
    rho(i1:i2,j1-1,k1:k2,ik)=0.0
   endif
  endif
  if(pe1y)then
   if(iby <2)then
    rho(i1:i2,j2-1,k1:k2,ik)=rho(i1:i2,j2-1,k1:k2,ik)+ &
     rho(i1:i2,j2+1,k1:k2,ik)
    rho(i1:i2,j2,k1:k2,ik)=rho(i1:i2,j2-1,k1:k2,ik)
   endif
  endif
 endif
 !============== field data on [y_loc]
 if(ibx <2 )then
  if(pex0)then
   rho(i1+1,j1:j2,k1:k2,ik)=rho(i1+1,j1:j2,k1:k2,ik)+ &
    rho(i1-1,j1:j2,k1:k2,ik)
   rho(i1,j1:j2,k1:k2,ik)=rho(i1+1,j1:j2,k1:k2,ik)
  endif
  if(pex1)then
   rho(i2-1,j1:j2,k1:k2,ik)=rho(i2-1,j1:j2,k1:k2,ik)+ &
    rho(i2+1,j1:j2,k1:k2,ik)
   rho(i2,j1:j2,k1:k2,ik)=rho(i2-1,j1:j2,k1:k2,ik)
  endif
 endif
 end subroutine den_zyxbd
 !====================
 subroutine fill_curr_yzxbdsdata(curr,i1,nxc,j1,nyc,k1,nzc,nc)
 integer,intent(in) :: i1,nxc,j1,nyc,k1,nzc,nc
 real(dp),intent(inout) :: curr(:,:,:,:)
 integer :: s1,s2,r1,r2,iy1,iy2,iz1,iz2,ix1,ix2
 integer :: ic,ix,j,iy,iz,kk,lenws,lenwr
 integer,parameter :: str=3,str2=2
 !================
 ! enter currents on a five-point extended stencil
 !===========================
 iz1=k1-str2
 iz2=nzc+str
 iy1=j1-str2
 iy2=nyc+str
 if(ndim <3)then
  iz1=k1;iz2=nzc
 endif
 if(ndim <2)then
  iy1=k1;iy2=nyc
 endif
 ix1=i1-str2
 ix2=nxc+str
 lenwr=str*nc*(ix2+1-ix1)*max(iz2+1-iz1,iy2+1-iy1)
 if(size(aux1) <lenwr)then
  deallocate(aux1,aux2)
  allocate(aux1(lenwr))
  allocate(aux2(lenwr))
 endif
 if(prly)then
  ! [j1-2:j1-1] left-y data
  s1=j1-str
  kk=0
  do ic=1,nc
   do iz=iz1,iz2
    do j=1,str2
     iy=s1+j
     do ix=ix1,ix2
      kk=kk+1
      aux1(kk)=curr(ix,iy,iz,ic)
     end do
    end do
   end do
  end do
  lenws=kk
  lenwr=lenws
  call exchange_bdx_data(aux1,aux2,lenws,lenwr,1,rt)
  !=====================
  ! sends y=[j1-2;j1-1] str-1 data to left
  !receives from right and adds data on y=[nyc-1:nyc] sign=+1
  r1=nyc-str2
  kk=0
  if(pe1y)then
   if(iby < 2)then
    aux2(1:lenwr)=0.0
   endif
  endif
  do ic=1,nc
   do iz=iz1,iz2
    do j=1,str2
     iy=j+r1
     do ix=ix1,ix2
      kk=kk+1
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+aux2(kk)
     end do
    end do
   end do
  end do
  ! sends y=[nyc+1:nyc+str] str=3 data to the right
  !receives from left and adds data on y=[j1:j1+str2] sign=-1
  s2=nyc
  kk=0
  do ic=1,nc
   do iz=iz1,iz2
    do j=1,str
     iy=j+s2
     do ix=ix1,ix2
      kk=kk+1
      aux1(kk)=curr(ix,iy,iz,ic)
     end do
    end do
   end do
  end do
  lenws=kk
  lenwr=lenws
  call exchange_bdx_data(aux1,aux2,lenws,lenwr,1,lt)
  !=====================
  r2=j1-1
  kk=0
  if(pe0y)then
   if(iby < 2)then
    aux2(1:lenwr)=0.0
   endif
  endif
  do ic=1,nc
   do iz=iz1,iz2
    do j=1,str
     iy=j+r2
     do ix=ix1,ix2
      kk=kk+1
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+aux2(kk)
     end do
    end do
   end do
  end do
 endif
 ! The reduced stencil of summed data
 iy1=j1
 iy2=nyc
 !================
 if(prlz)then
  !================
  ! sends z=[k1-2;k1-1] str-1 data to left
  ! receives from right and adds data on y=[nzc-1:nzc] sign=+1

  s1=k1-str
  kk=0
  do ic=1,nc
   do j=1,str2
    iz=s1+j
    do iy=iy1,iy2
     do ix=ix1,ix2
      kk=kk+1
      aux1(kk)=curr(ix,iy,iz,ic)
     end do
    end do
   end do
  end do
  lenws=kk
  lenwr=lenws
  call exchange_bdx_data(aux1,aux2,lenws,lenwr,2,rt)

  r1=nzc-str2
  if(pe1z)then
   if(ibz < 2)then
    aux2(1:lenwr)=0.0
   endif
  endif
  kk=0
  do ic=1,nc
   do j=1,str2
    iz=j+r1
    do iy=iy1,iy2
     do ix=ix1,ix2
      kk=kk+1
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+aux2(kk)
     end do
    end do
   end do
  end do
  !
  !================
  ! sends z=[nzc+1:nzc+str] str=3 data to the right
  !receives from left and adds data on z=[k1:k1+str2] sign=-1
  s2=nzc
  kk=0
  do ic=1,nc
   do j=1,str
    iz=j+s2
    do iy=iy1,iy2
     do ix=ix1,ix2
      kk=kk+1
      aux1(kk)=curr(ix,iy,iz,ic)
     end do
    end do
   end do
  end do
  lenws=kk
  lenwr=lenws
  call exchange_bdx_data(aux1,aux2,lenws,lenwr,2,lt)
  !================
  ! Recvs and adds on z[k1:k1:str2]
  !================
  r2=k1-1
  if(pe0z)then
   if(ibz < 2)then
    aux2(1:lenwr)=0.0
   endif
  endif
  kk=0
  do ic=1,nc
   do j=1,str
    iz=j+r2
    do iy=iy1,iy2
     do ix=ix1,ix2
      kk=kk+1
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+aux2(kk)
     end do
    end do
   end do
  end do
  ! The reduced stencil of summed data
  iz1=k1
  iz2=nzc
 endif
 !========================== prlx case
 if(prlx)then
  ! sends x=[i1-2;i1-1] str-1 data to left
  ! receives from right and adds data on x=[nxc-1:nxc] sign=+1
  s1=i1-str
  kk=0
  do ic=1,nc
   do iz=iz1,iz2
    do iy=iy1,iy2
     do j=1,str2
      ix=s1+j
      kk=kk+1
      aux1(kk)=curr(ix,iy,iz,ic)
     end do
    end do
   end do
  end do
  lenws=kk
  lenwr=lenws
  call exchange_bdx_data(aux1,aux2,lenws,lenwr,3,rt)
  !=====================
  ! sends x=[i1-2;i1-1] str-1 data to left
  !receives from right and adds data on x=[nxc-1:nxc] sign=+1
  r1=nxc-str2
  kk=0
  if(pex1)then
   if(ibx < 2)then
    aux2(1:lenwr)=0.0
   endif
  endif
  do ic=1,nc
   do iz=iz1,iz2
    do iy=iy1,iy2
     do j=1,str2
      ix=j+r1
      kk=kk+1
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+aux2(kk)
     end do
    end do
   end do
  end do
  ! sends x=[nxc+1:nxc+str] str=3 data to the right
  !receives from left and adds data on x=[i1:i1+str2] sign=-1
  s2=nxc
  kk=0
  do ic=1,nc
   do iz=iz1,iz2
    do iy=iy1,iy2
     do j=1,str
      ix=j+s2
      kk=kk+1
      aux1(kk)=curr(ix,iy,iz,ic)
     end do
    end do
   end do
  end do
  lenws=kk
  lenwr=lenws
  call exchange_bdx_data(aux1,aux2,lenws,lenwr,3,lt)
  !=====================
  r2=i1-1
  kk=0
  if(pex0)then
   if(ibx < 2)then
    aux2(1:lenwr)=0.0
   endif
  endif
  do ic=1,nc
   do iz=iz1,iz2
    do iy=iy1,iy2
     do j=1,str
      ix=j+r2
      kk=kk+1
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+aux2(kk)
     end do
    end do
   end do
  end do
  return
 endif
 !===================== No MPI x-decomposition
 if(ibx==2)then
  ! data nxc-1:nxc sums to i1-2:i1-1
  ! data i1:i1+2 sums to nxc+1:nxc+3
  s1=i1-str
  r1=nxc-str2
  do ic=1,nc
   do iz=iz1,iz2
    do iy=iy1,iy2
     do j=1,str2
      curr(r1+j,iy,iz,ic)=curr(r1+j,iy,iz,ic)+curr(s1+j,iy,iz,ic)
     end do
     do j=1,str
      ix=i1-1+j
      curr(ix,iy,iz,ic)=curr(ix,iy,iz,ic)+curr(nxc+j,iy,iz,ic)
     end do
    end do
   end do
  end do
 endif
 !=============================
 end subroutine fill_curr_yzxbdsdata
 !=================
 end module mpi_curr_interface
