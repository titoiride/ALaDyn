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

 module init_beam_part_distrib

 use util
 use psolve
 use array_alloc
 use mpi_field_interface
 use mpi_curr_interface
 use init_grid_field
 use grid_part_connect

 implicit none

 real(dp),allocatable :: bpart(:,:)
 !--------------------------

 contains
 !=========================
 subroutine beam_data(ndm,np_tot)  !generates bpart(7,np_tot)
 integer,intent(in) :: ndm
 integer,intent(out) :: np_tot
 integer :: i,i1,i2,ip
 real(dp) :: cut,xh(5)
 integer :: nch
 logical :: sr
 !==========================================================
 ! bconf defines only the configuration of bunch charges
 !=======================
 ! default values
 !=======================
 np_tot=0
 do i=1,nsb
  np_tot=np_tot+nb_tot(i)
 end do
 cut=3.
 nch=2*ndm+1
 allocate(bpart(nch,np_tot))
  !---!
  i1=1
  charge=int(unit_charge(1),hp_int)      !nsb electron bunches
  part_ind=-1

  select case(ndm)
  case(2)
  do ip=1,nsb
   xh(ip)=xc_bunch(ip)
   wgh=real(j0_norm*jb_norm(ip),sp) !the bunch particle weight
   i2=i1+nb_tot(ip)-1
   !---Original Version---!
    call bunch_gen(ndm,i1,i2,sxb(ip),syb(ip),syb(ip),gam(ip),&
                   epsy(ip),epsz(ip),cut,dg(ip),bpart)
    do i=i1,i2
     bpart(1,i)=bpart(1,i)+xh(ip)       !x-shifting
     bpart(2,i)=bpart(2,i)+yc_bunch(ip) !y-shifting
     bpart(nch,i)=wgh_cmp
    end do
   i1=i2+1
  end do
 ! Pe0 p data are copied to all MPI tasks
  case(3)
  do ip=1,nsb
   xh(ip)=xc_bunch(ip)
   wgh=real(j0_norm*jb_norm(ip),sp) !the bunch particles weights
   i2=i1+nb_tot(ip)-1
   !---Original Version---!
    call bunch_gen(ndm,i1,i2,sxb(ip),syb(ip),syb(ip),gam(ip),&
    epsy(ip),epsz(ip),cut,dg(ip),bpart)
    do i=i1,i2
     bpart(1,i)=bpart(1,i)+xh(ip)       !x-shifting
     bpart(2,i)=bpart(2,i)+yc_bunch(ip) !y-shifting
     bpart(3,i)=bpart(3,i)+zc_bunch(ip) !z-shifting
     bpart(nch,i)=wgh_cmp
    end do
   i1=i2+1
  end do
  end select
 ! Pe0 p data are copied to all MPI tasks
 if(pe0)then
  sr=.true.
  do ip=1,npe-1
   call exchange_2d_grdata(sr,bpart,nch,np_tot,ip,ip+10)
  end do
 else
  sr=.false.
  call exchange_2d_grdata(sr,bpart,nch,np_tot,0,mype+10)
 endif
 !==============================
 end subroutine beam_data
 !===================
 subroutine MPI_beam_distribute(ndm)
 
 integer,intent(in) :: ndm

 integer :: i,ii,i1,j
 integer :: ic,p,ip,ipp,nb_loc
 real(dp) :: y1,y2,z1,z2,x1,x2
 integer :: nps_loc(nsb),npmax,np_tot
 !========= count particles on each (yz) MPI domain
 np_tot=sum(nb_tot(1:nsb))
 ! ALL MPI tasks do
 x1=loc_xgrid(imodx)%gmin
 x2=loc_xgrid(imodx)%gmax
 i1=0
 select case(ndm)
 case(2)
 ip=npe_zloc-1
 do ic=1,nsb
  do p=0,npe_xloc-1
   do ipp=0,npe_yloc-1
    y1=loc_ygrid(ipp)%gmin
    y2=loc_ygrid(ipp)%gmax
    loc_nbpart(ipp,ip,p,ic)=0
    do j=1,nb_tot(ic)
     i=i1+j
     if(bpart(2,i) >y1.and.bpart(2,i) <=y2)then
      if(bpart(1,i) >x1.and.bpart(1,i) <=x2)then
        loc_nbpart(ipp,ip,p,ic)=loc_nbpart(ipp,ip,p,ic)+1
      endif
     endif
    end do
   end do
  end do
  i1=i1+nb_tot(ic)
 end do
 nb_max=maxval(loc_nbpart(0:npe_yloc-1,0:npe_zloc-1,0:npe_xloc-1,1:nsb))
 nb_min=minval(loc_nbpart(0:npe_yloc-1,0:npe_zloc-1,0:npe_xloc-1,1:nsb))
 do ic=1,nsb
  do p=0,npe_xloc-1
   do ipp=0,npe_yloc-1
    i=ipp+npe_yloc*(ip+p*npe_zloc)
    if(loc_nbpart(ipp,ip,p,ic)==nb_max)pe_nbmax=i
    if(loc_nbpart(ipp,ip,p,ic)==nb_min)pe_nbmin=i
   end do
  end do
 end do
 !==================
 ! The local MPI task
 nps_loc(1:nsb)=loc_nbpart(imody,imodz,imodx,1:nsb)
 npmax=maxval(nps_loc(1:nsb))
 npmax=max(npmax,1)
 if(.not.allocated(bunch(1)%part))then
  allocate(bunch(1)%part(npmax,nd2+1))
 endif
 !=================================
 nb_loc=nps_loc(1)
 p=imodx
 ip=imodz
 ipp=imody
 y1=loc_ygrid(ipp)%gmin
 y2=loc_ygrid(ipp)%gmax
 !=========================
 ! Here 2D MPI decomp. allowed
 !===================================
 i1=0
 do ic=1,nsb
  ii=0
  do i=1,nb_tot(ic)
   j=i+i1
   if(bpart(2,j) >y1.and.bpart(2,j) <=y2)then
    if(bpart(1,j) >x1.and.bpart(1,j) <=x2)then
     ii=ii+1
     bunch(ic)%part(ii,1:nd2+1)=bpart(1:nd2+1,j)
    endif
   endif
  end do
  i1=i1+nb_tot(ic)
 enddo
 case(3)
 do ic=1,nsb
  do p=0,npe_xloc-1
   do ip=0,npe_zloc-1
    z1=loc_zgrid(ip)%gmin
    z2=loc_zgrid(ip)%gmax
    do ipp=0,npe_yloc-1
     y1=loc_ygrid(ipp)%gmin
     y2=loc_ygrid(ipp)%gmax
     loc_nbpart(ipp,ip,p,ic)=0
     do j=1,nb_tot(ic)
      i=i1+j
      if(bpart(2,i) >y1.and.bpart(2,i) <=y2)then
       if(bpart(3,i) >z1.and.bpart(3,i) <=z2)then
        if(bpart(1,i) >x1.and.bpart(1,i) <=x2)then
         loc_nbpart(ipp,ip,p,ic)=loc_nbpart(ipp,ip,p,ic)+1
        endif
       endif
      endif
     enddo
    end do
   end do
  end do
  i1=i1+nb_tot(ic)
 end do
 nb_max=maxval(loc_nbpart(0:npe_yloc-1,0:npe_zloc-1,0:npe_xloc-1,1:nsb))
 nb_min=minval(loc_nbpart(0:npe_yloc-1,0:npe_zloc-1,0:npe_xloc-1,1:nsb))
 do ic=1,nsb
  do p=0,npe_xloc-1
   do ip=0,npe_zloc-1
    do ipp=0,npe_yloc-1
     i=ipp+npe_yloc*(ip+p*npe_zloc)
     if(loc_nbpart(ipp,ip,p,ic)==nb_max)pe_nbmax=i
     if(loc_nbpart(ipp,ip,p,ic)==nb_min)pe_nbmin=i
    end do
   end do
  end do
 end do
 !==================
 ! The local MPI task
 nps_loc(1:nsb)=loc_nbpart(imody,imodz,imodx,1:nsb)
 npmax=maxval(nps_loc(1:nsb))
 !++++++++++++++++++++++++++++++++++++++
 npmax=max(npmax,1)
 if(.not.allocated(bunch(1)%part))then
  allocate(bunch(1)%part(npmax,nd2+1))
 endif
 !=================================
 nb_loc=nps_loc(1)
 p=imodx
 ip=imodz
 z1=loc_zgrid(ip)%gmin
 z2=loc_zgrid(ip)%gmax
 ipp=imody
 y1=loc_ygrid(ipp)%gmin
 y2=loc_ygrid(ipp)%gmax
 !=========================
 ! Here 3D MPI decomp. allowed
 !===================================
 i1=0
 do ic=1,nsb
  ii=0
  do i=1,nb_tot(ic)
   j=i+i1
   if(bpart(2,j) >y1.and.bpart(2,j) <=y2)then
    if(bpart(3,j) >z1.and.bpart(3,j) <=z2)then
     if(bpart(1,j) >x1.and.bpart(1,j) <=x2)then
      ii=ii+1
      bunch(ic)%part(ii,1:nd2+1)=bpart(1:nd2+1,j)
     endif
    endif
   endif
  end do
  i1=i1+nb_tot(ic)
 enddo
 end select
 !=================================
 end subroutine MPI_beam_distribute
 !========================
 subroutine beam_model_pot(pot,sx,sy,sz,b_am,i1,i2,j1,j2,k1,k2)

 real(dp),intent(inout) :: pot(:,:,:,:)
 real(dp),intent(in) :: sx,sy,sz,b_am
 integer,intent(in) :: i1,i2,j1,j2,k1,k2
 integer :: i,j,k,jj,kk
 real(dp) :: r2,brad2,sx2_inv,fact,r2max,pot0
 real(dp) :: xx,yy,zz
 !-----------------------
 brad2=sy*sy+sz*sz
 r2max=ymax*ymax+zmax*zmax
 pot0=log(r2max/brad2)
 sx2_inv=1./(2.*sx*sx)
 do k=k1,k2
  kk=k-2
  zz=loc_zg(kk,1,imodz)
  do j=j1,j2
   jj=j-2
   yy=loc_yg(jj,1,imody)
   r2=zz*zz+yy*yy
   if(r2 > brad2)then
    fact=log(r2/brad2)-pot0
   else
    fact=(r2/brad2-1.)-pot0
   endif
   do i=i1,i2
    xx=(x(i1)-xc_bunch(1))
    pot(i,j,k,1)=0.25*brad2*b_am*fact*exp(-xx*xx*sx2_inv)
   end do
  end do
 end do
 end subroutine beam_model_pot
 !=========================
 subroutine beam_inject

 integer :: ii,n,i1,i2,i2b,j1,k1,nyp,nzp
 integer :: id_ch,nptot,np,nb,n_st,ic,ft_mod,ft_sym
 integer :: nb_loc(1),ix,iy,iz
 real(dp) :: gam2,dmax(1),all_dmax(1)
 real(dp) :: xm,ym,zm

 gam2=gam0*gam0
 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 i2b=i2
 j1=loc_ygrid(imody)%p_ind(1)
 nyp=loc_ygrid(imody)%p_ind(2)
 k1=loc_zgrid(imodz)%p_ind(1)
 nzp=loc_zgrid(imodz)%p_ind(2)
 id_ch=nd2+1
 !=======================
 !!!!!!!!!!!!!!
 if(.not.allocated(ebf_bunch))then
  allocate(ebf_bunch(i2+3,nyp+3,nzp+3,nbfield))
  ebf_bunch(:,:,:,:)=0.0
 endif
 !=============================
 ! The fields of a moving e-bunch in vacuum
 ! E_x=-(DPhi/Dx)/gamma^2  , E_y=-DPhi/Dy   E_z=-DPhi/Dz
 ! Poisson eq.  D_xE_x+D_yE_y+D_zE_z=omp^2\rho (x-V_b*t_0,y,z)
 ! B_x=0   B_y=-V_b*E_z    B_z= V_b*E_y
 !=========================================
 call beam_data(ndim,nptot)
 ! Generates phase space coordinates for beam on bpart(7,np_tot)
 ! bpart() provisional storage in common to all MPI tasks
 xm=loc_xgrid(imodx)%gmin
 ym=loc_ygrid(imody)%gmin
 zm=loc_zgrid(imodz)%gmin
 !=======================
 call MPI_beam_distribute(ndim) !local bpart data are stored in bunch(1)%part struct for each MPI task
 nb=loc_nbpart(imody,imodz,imodx,1)
 if(nb >0)then
  if(allocated(spec(1)%part))then
   np=size(spec(1)%part,1)
   do n=1,np
    ebfp(n,1:id_ch)=spec(1)%part(n,1:id_ch)
   end do
   deallocate(spec(1)%part)
   allocate(spec(1)%part(np+nb,id_ch))
   do n=1,np
    spec(1)%part(n,1:id_ch)=ebfp(n,1:id_ch)
   end do
   deallocate(ebfp)
   do n=1,nb
    ii=n+np
    spec(1)%part(ii,1:id_ch)=bunch(1)%part(n,1:id_ch)
   end do
   allocate(ebfp(np+nb,id_ch))
   do ix=0,npe_xloc-1
    do iz=0,npe_zloc-1
     do iy=0,npe_yloc-1
      loc_npart(iy,iz,ix,1)=loc_npart(iy,iz,ix,1) +&
                           loc_nbpart(iy,iz,ix,1)
     end do
    end do
   end do
  else
   nb_loc(1)=nb
   call p_alloc(nb,nd2+1,nb_loc,nsb,LPf_ord,1,1,mem_psize)
   do n=1,nb
    spec(1)%part(n,1:id_ch)=bunch(1)%part(n,1:id_ch)
   end do
   do ix=0,npe_xloc-1
    do iz=0,npe_zloc-1
     do iy=0,npe_yloc-1
      loc_npart(iy,iz,ix,1)=loc_nbpart(iy,iz,ix,1) 
     end do
    end do
   end do
  endif
 endif
 if(allocated(bunch(1)%part))deallocate(bunch(1)%part)
 !==================== Computes the bunch inumber density
 jc(:,:,:,1)=0.0
 n_st=0
 do ic= 1,nsb
  np=loc_nbpart(imody,imodz,imodx,ic)
  if(np>0)call set_grid_charge(spec(ic),ebfp,jc,np,ndim,n_st,1,xm,ym,zm)
 end do
 !generates injc(1)=den(i,j,k)  jc(2)= Jx(i,j,k)
 if(prl)call fill_curr_yzxbdsdata(jc,i1,i2b,j1,nyp,k1,nzp,1)
 dmax(1)=0.0
 do iz=k1,nzp
  do iy=j1,nyp
   do ix=i1,i2
    dmax(1)=max(dmax(1),abs(jc(ix,iy,iz,1)))
   end do
  end do
 end do 
 all_dmax(1)=dmax(1)
 if(prl)call allreduce_dpreal(MAXV,dmax,all_dmax,1)
 !============================
 ! BUNCH grid DENSITY and current already normalized
 !In jc(1) beam density  in jc(2) beam Jx current
 jc(i1:i2,j1:nyp,k1:nzp,2)=bet0*jc(i1:i2,j1:nyp,k1:nzp,1)
 !=====================================================
 ! UNIFORM GRIDS ASSUMED
 !======================
 if(allocated(bpart))deallocate(bpart)
 ft_mod=2       !A sine transform along each coordinate
 ft_sym=1
 if(ndim==2)call FFT_2D_Psolv(jc,ompe,nx,nx_loc,ny,ny_loc,nz,nz_loc,&
                i1,i2b,j1,nyp,k1,nzp,ft_mod,ft_sym)
 if(ndim==3)call FFT_Psolv(jc,gam2,ompe,nx,nx_loc,ny,ny_loc,nz,nz_loc,&
                i1,i2b,j1,nyp,k1,nzp,ft_mod,ft_sym)
 !Solves Laplacian[pot]=ompe*rho
 !Beam potential in jc(1) 
 !===================
 !====================
 call fill_ebfield_yzxbdsdata(jc,i1,i2b,j1,nyp,k1,nzp,1,2,1,1)
 call initial_beam_fields(jc,ebf_bunch,i1,i2b,j1,nyp,k1,nzp,gam2,bet0)
 ! generates (Ex,Ey,Ez,By,Bz) bunch fields  Bx=ebf_bunc(4)=0
 ebf_bunch(:,:,:,4)=0.0
!========================================= Collect data
 ebf(i1:i2,j1:nyp,k1:nzp,1:nfield)=ebf(i1:i2,j1:nyp,k1:nzp,1:nfield)+ &
                                    ebf_bunch(i1:i2,j1:nyp,k1:nzp,1:nfield)
  deallocate(ebf_bunch)
 !========================================
 lp_end(1)=xc_bunch(1)+2.*sxb(1)
 !=====================================
 if(pe0)then
  !==================
  open(16,file='Initial_bunch_info.dat')
  write(16,'(a27,e11.4)')'Initial target x-position =',targ_in
  write(16,'(a20,e11.4)')'Plasma wave-length =',lambda_p
  write(16,*) '-------------------------------------'
  i1=1
   write(16,'(a13,i4)')'Bunch number ',i1
   write(16,'(a25,e11.4)')'Bunch particles per cell ',nb_per_cell
   write(16,'(a23,i8)')'Bunch particle number  ',nb_tot(i1)
   write(16,'(a17,3e11.4)')'Sizes and gamma  ',sxb(i1),syb(i1),gam(i1)
   write(16,'(a23,2e11.4)')'Transverse emittances= ',epsy(i1),epsz(i1)
   write(16,'(a21,e11.4)')'Initial xc-position= ',xc_bunch(i1)
   write(16,'(a20,e11.4)')'B charge    [pC] =  ',bunch_charge(i1)
   write(16,'(a21,e11.4)')'B_density/P_density  ', rhob(i1)
   write(16,'(a19,e11.4)')'Bunch max density  ', all_dmax(1)
  close(16)
 endif
 end subroutine beam_inject
!=====================================
 subroutine Bpulse(part_in)

 real(dp),intent(out) :: part_in
 integer :: i1,i2,i2b,j1,k1,nyp,nzp
 integer :: nptot,np,n_st,ic,ft_mod,ft_sym
 real(dp) :: gam2
 real(dp) :: xm,ym,zm
!======================================
 if(.not.Part)Part=.true.
 !!!!!!!!!!!!!!
 !=============================
 ! The fields of moving bunches in vacuum
 ! E_x=-(DPhi/Dx)/gamma^2  , E_y=-DPhi/Dy   E_z=-DPhi/Dz
 ! Poisson eq.  D_xE_x+D_yE_y+D_zE_z=omp^2\rho (x-V_b*t_0,y,z)
 ! B_x=0   B_y=-V_b*E_z    B_z= V_b*E_y
 !=========================================
 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 j1=loc_ygrid(imody)%p_ind(1)
 nyp=loc_ygrid(imody)%p_ind(2)
 k1=loc_zgrid(imodz)%p_ind(1)
 nzp=loc_zgrid(imodz)%p_ind(2)
 i2b=i2
 gam2=gam0*gam0
 !=======================
 call beam_data(ndim,nptot)
 ! Generates phase space coordinateis for all beams on bpart(7,np_tot)
 ! bpart in common to all MPI tasks
 xm=loc_xgrid(imodx)%gmin
 ym=loc_ygrid(imody)%gmin
 zm=loc_zgrid(imodz)%gmin
 !=======================
 call MPI_beam_distribute(ndim) !bpart are loaded on bunc() struct for each MPI task
 !==================== Computes the total bunches density
 jc(:,:,:,1)=0.0
 n_st=0
 do ic= 1,nsb
  np=loc_nbpart(imody,imodz,imodx,ic)
  call set_grid_charge(bunch(ic),ebfb,jc,np,ndim,n_st,1,xm,ym,zm)
 end do
 !generates ebf_bunc(1)=den(i,j,k)
 !on local MPI computational grid[i2,nyp,nzb]
 if(prl)call fill_curr_yzxbdsdata(jc,i1,i2b,j1,nyp,k1,nzp,1)
 !============================
 ! BUNCH grid DENSITY already normalized
 ! index=3 for bunch current: a same bet0 assumed for all bunches
 !===========================================
 jc(i1:i2b,j1:nyp,k1:nzp,2)=bet0*jc(i1:i2b,j1:nyp,k1:nzp,1)
 !Beam Ax potential in jc(2)    Ay=Az=0 assumed
 !In jc(1) beam density  in jc(2) beam Jx current
 !=====================================================
 ! UNIFORM GRIDS ASSUMED
 !======================
 if(allocated(bpart))deallocate(bpart)
 ft_mod=2       !A sine transform along each coordinate
 ft_sym=1
 if(ndim==2)call FFT_2D_Psolv(jc,ompe,nx,nx_loc,ny,ny_loc,nz,nz_loc,&
                i1,i2b,j1,nyp,k1,nzp,ft_mod,ft_sym)
 if(ndim==3)call FFT_Psolv(jc,gam2,ompe,nx,nx_loc,ny,ny_loc,nz,nz_loc,&
                i1,i2b,j1,nyp,k1,nzp,ft_mod,ft_sym)
 !Solves Laplacian[pot]=ompe*rho
 !Beam potential in jc(1) 
 !===================
 !====================
 call fill_ebfield_yzxbdsdata(jc,i1,i2b,j1,nyp,k1,nzp,1,2,1,1)
  call initial_beam_fields(jc,ebf_bunch,i1,i2b,j1,nyp,k1,nzp,gam2,bet0)
  ! generates (Ex,Ey,Ez,By,Bz) bunch fields
  if(ibeam>0)then
   ebf1_bunch(:,:,:,:)=0.0
  endif
 !========================================
 lp_end(1)=xc_bunch(1)+2.*sxb(1)
 !=====================================
 if(L_Bpoloidal) call set_poloidal_ex_fields( &
  ebf0_bunch,i1,i2b,j1,nyp,k1,nzp,B_ex_poloidal*T_unit,radius_poloidal)
 !=====================================
 part_in= lpx(7)+ lp_end(1)
 !============================
 !----------------------
 if(pe0)then
  !==================
  open(16,file='Initial_bunch_info.dat')
  write(16,'(a22,i4,e11.4)')' Plasma macro per cell',mp_per_cell(1),j0_norm
  write(16,'(a18,i4)')' Beam injected nb=',nsb
  write(16,'(a27,e11.4)')' Initial target x-position=',targ_in
  write(16,'(a20,e11.4)')' Plasma wave-length=',lambda_p
  write(16,*) '-------------------------------------'
  write(16,*)'L_part',L_particles
  do i1=1,nsb
   write(16,'(a13,i4)')' Bunch number',i1
   write(16,'(a12,i4)')' bunch type ',bunch_type(i1)
   write(16,'(a31,e11.4)')' relative bunch/particle weights',jb_norm(i1)
   write(16,'(a23,i8)')' bunch particle number ',nb_tot(i1)
   write(16,'(a17,3e11.4)')' sizes and gamma ',sxb(i1),syb(i1),gam(i1)
   write(16,'(a23,2e11.4)')' Transverse emittances=',epsy(i1),epsz(i1)
   write(16,'(a21,e11.4)')' Initial xc-position=',xc_bunch(i1)
   write(16,'(a20,2e11.4)')' b charge [pC],Qch =',bunch_charge(i1),reduced_charge(i1)
   write(16,'(a22,e11.4)')' bcharge_over_pcharge ', rhob(i1)
  end do
  close(16)
 endif
 end subroutine Bpulse
 !========================
 end module init_beam_part_distrib
