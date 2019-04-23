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

 module init_laser_field

 use array_wspace
 use init_grid_field

 implicit none

 contains
 !-------------------
 !+++++++++++++++++++++

 subroutine LP_pulse(lp_mod,part_in)

 integer,intent(in) :: lp_mod
 real(dp),intent(out) :: part_in
 integer :: ic,lp_ind,i1,i2
 real(dp) :: angle,shx_lp,sigm,eps,xm,tt,tau,tau1
 !===========================
 ! field grid index defined on set_pgrid

 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 xm=loc_xgrid(imodx)%gmin  ! => field index i1
 tt=0.0
 if(G_prof)then
  tau=w0_x*sqrt(2.)
  tau1=w1_x*sqrt(2.)
 else
  tau=0.5*w0_x   !cos^2() profile_
  tau1=0.5*w1_x   !cos^2() profile_
 endif
 lp_amp=oml*a0
 lp1_amp=om1*a1
!=====================
 lp_in(1)=xc_lp-tau
 lp_end(1)=xc_lp+tau
 eps=1./(oml*w0_y)
 sigm=lam0/w0_x
 angle=lpx(6)
 xf=xc_lp+t0_lp
  !=======================
 lp_ind=lp_mod
 shx_lp=0.0
 xc_loc(1)=xc_lp
 xf_loc(1)=xf
 if(Plane_wave)lp_ind=0
 if(angle >0.0)then
   shx_lp=lpx(7)
  if(lp_end(1) > xm)then
   call init_lp_fields(ebf,lp_amp,tt,t0_lp,w0_x,w0_y,xf,oml,&
                                       angle,shx_lp,lp_ind,i1,i2)
  endif
 else         !normal incidence
  if(lp_end(1) > xm)then
   call init_lp_inc0_fields(ebf,lp_amp,tt,t0_lp,w0_x,w0_y,xf,oml,lp_ind,i1,i2)
  endif
 endif
 if(nb_laser >1)then
  do ic=2,nb_laser
   lp_in(ic)=lp_in(ic-1)-lp_delay
   lp_end(ic)=lp_end(ic-1)-lp_delay
   xc_loc(ic)=xc_loc(ic-1)-lp_delay
   xf_loc(ic)=xc_loc(ic)+t0_lp
   if(lp_end(ic) > xm)call init_lp_inc0_fields(ebf,lp_amp,tt,t0_lp,w0_x,w0_y,xf_loc(ic),oml,&
                                              lp_ind,i1,i2)
  end do
 endif
!=================TWO-COLOR
 if(Two_color)then
  xc1_lp=xc_loc(nb_laser)-lp_offset
  xf1=xc1_lp+t1_lp
  
  lp_ionz_in=xc1_lp-tau1
  lp_ionz_end=xc1_lp+tau1
  call init_lp_inc0_fields(ebf,lp1_amp,tt,t1_lp,w1_x,w1_y,xf1,om1,&
                                              lp_ind,i1,i2)
  if(pe0)write(6,'(a30,e11.4)')'two-color activated at xc1_lp=',xc1_lp
 endif 
  
!==================================
 part_in=lp_end(1)+lpx(7)
 end subroutine LP_pulse
 !===========================
 subroutine CP_pulse(cp_mod,part_in)

 real(dp),intent(out) :: part_in
 integer,intent(in) :: cp_mod
 integer :: i1,i2,cp_ind
 real(dp) :: angle,shx_cp,eps,sigm,xm,tau,tt

 lp_amp=oml*a0/sqrt(2.)
 !-------------------------
 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 xm=loc_xgrid(imodx)%gmin  ! => field index i1
 tt=0.0
 tau=0.5*w0_x
 lp_in(1)=xc_lp-tau          !Pulse centroid xc_lp=xf-ts at time t=0 lp_in=xc-tau
 lp_end(1)=lp_in(1)+w0_x
 eps=1./(oml*w0_y)
 sigm=lam0/w0_x
 angle=lpx(6)
 xf=xc_lp+t0_lp
 shx_cp=0.0
 if(angle >0.0)shx_cp=lpx(7)
 if(lp_amp >0.0)then
  cp_ind=cp_mod
  if(Plane_wave)cp_ind=0
  !=======================
  call init_cp_fields(ebf,lp_amp,tt,t0_lp,w0_x,w0_y,xf,&
   angle,shx_cp,cp_ind,i1,i2)
  !=================def part distr points
  lp_end=xm
 endif
!======================
  part_in=lp_end(1)+lpx(7)
 end subroutine CP_pulse
 !-------------------------
 subroutine set_envelope(part_in)

 real(dp),intent(out) :: part_in
 integer :: ic,pw_ind,i1,i2,j1,nyp,k1,nzp
 real(dp) :: eps,sigm,xm,tt,tau,tau1,lp1_amp,loc_delay(3)

 lp_amp=a0
 lp1_amp=a1
 !===========================
 ! field grid index defined on set_pgrid

 i1=loc_xgrid(imodx)%p_ind(1)
 i2=loc_xgrid(imodx)%p_ind(2)
 j1=loc_ygrid(imody)%p_ind(1)
 nyp=loc_ygrid(imody)%p_ind(2)
 k1=loc_zgrid(imodz)%p_ind(1)
 nzp=loc_zgrid(imodz)%p_ind(2)
 xm=loc_xgrid(imodx)%gmin
 !=========================
 tt=0.0
 pw_ind=1
 if(Plane_wave)pw_ind=0
 if(G_prof)then
  tau=w0_x*sqrt(2.)
  tau1=w1_x*sqrt(2.)
 else
  tau=0.5*w0_x   !cos^2() profile_
  tau1=0.5*w1_x   !cos^2() profile_
 endif
 lp_in(1)=xc_lp-tau
 lp_end(1)=xc_lp+tau
 eps=1./(oml*w0_y)
 sigm=lam0/w0_x
 xf=xc_lp+t0_lp
  !=======================
 xc_loc(1)=xc_lp
 xf_loc(1)=xf
 env(:,:,:,:)=0.0
 if(lp_end(1) > xm)then
  if(G_prof)then
   call init_gprof_envelope_field(&
              env,a0,dt,tt,t0_lp,w0_x,w0_y,xf,oml,pw_ind,i1,i2)
  else
   call init_envelope_field(&
              env,a0,dt,tt,t0_lp,w0_x,w0_y,xf,oml,pw_ind,i1,i2)
   !call init_env_filtering(env,i1,i2,j1,nyp,k1,nzp)
  endif
 endif
 loc_delay(1)=lp_delay
 loc_delay(2)=lp_delay+lpy(1)
 loc_delay(3)=loc_delay(2)+lpy(2)
 if(nb_laser >1)then
  do ic=2,nb_laser
   lp_in(ic)=lp_in(ic-1)-loc_delay(ic-1)
   lp_end(ic)=lp_end(ic-1)-loc_delay(ic-1)
   xc_loc(ic)=xc_loc(ic-1)-loc_delay(ic-1)
   xf_loc(ic)=xc_loc(ic)+t0_lp
   if(lp_end(ic) > xm)then
    if(G_prof)then
     call init_gprof_envelope_field(&
              env,a0,dt,tt,t0_lp,w0_x,w0_y,xf_loc(ic),oml,pw_ind,i1,i2)
    else
     call init_envelope_field(&
              env,a0,dt,tt,t0_lp,w0_x,w0_y,xf_loc(ic),oml,pw_ind,i1,i2)
    endif
   endif
  end do
 endif
 if(Lpf_ord==4)then     !set initial first derivative
  env(i1:i2,j1:nyp,k1:nzp,3:4)=(env(i1:i2,j1:nyp,k1:nzp,1:2)-env(i1:i2,j1:nyp,k1:nzp,3:4))/dt
 endif 
!=================TWO-COLOR
 if(Two_color)then
  env1(:,:,:,:)=0.0
  xc1_lp=xc_loc(nb_laser)-lp_offset
  xf1=xc1_lp+t1_lp
   
  lp_ionz_in=xc1_lp-tau1
  lp_ionz_end=xc1_lp+tau1
  if(lp_ionz_end > xm)then
   if(G_prof)then
    call init_gprof_envelope_field(&
              env1,a1,dt,tt,t1_lp,w1_x,w1_y,xf1,om1,pw_ind,i1,i2)
   else
    call init_envelope_field(&
              env1,a1,dt,tt,t1_lp,w1_x,w1_y,xf1,om1,pw_ind,i1,i2)
   endif
  endif 
 endif
 !=======================
 !==========================
 ebf=0.0
 !=====================
 if(pe0)then
  open(26,file='Initial_env_info.dat')
  write(26,*)'number ',nb_laser,'LP envelope injected '
  write(26,*)' First pulse parameters'
  write(26,'(a21,e11.4)')' Focal x-position xf=',xf
  write(26,'(a21,e11.4)')' Centr x-position xc=',xc_lp
  write(26,'(a20,e11.4)')' Size (microns) lx=',w0_x
  write(26,'(a14,e11.4)')' FWHM (fs) lx=',tau_FWHM
  write(26,'(a19,e11.4)')' Waist(microns) ly=',w0_y
  write(26,'(a17,e11.4)')' FWHM(microns) = ',lp_rad
  if(nb_laser>1)write(26,'(a21,3e11.4)')' Delay among pulses= ',loc_delay(1:3)
  if(Two_color)then
   write(26,*)' Injection pulse parameters'
   write(26,'(a16,e11.4)')' amplitude  a1 =',a1
   write(26,'(a21,e11.4)')' Focal x-position xf=',xf1
   write(26,'(a21,e11.4)')' Centr x-position xc=',xc1_lp
   write(26,'(a20,e11.4)')' Size (microns) lx=',w1_x
   write(26,'(a19,e11.4)')' Waist(microns) ly=',w1_y
  endif
  close(26)
 endif
 part_in=lp_end(1)+lpx(7)
 !=================def part distr points
 end subroutine set_envelope

 end module init_laser_field
