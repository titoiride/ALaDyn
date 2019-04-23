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

 program ALaDyn

 use start_all
 use diag_part_and_fields
 use run_data_info
 use pic_out_util
 use pic_out
 use pic_evolve
 use env_evolve
 use bunch_evolve

 implicit none

 call start

!============= info related to initial data
 if(pe0)then
  write(6,*)'START OF RUN'
  write(6,*)'Running on ',npe,' cpu'
  write(6,'(a32,3i4)')' MPI decomposition along x-y-z: ',npe_xloc,npe_yloc,npe_zloc
 endif
 !====To activate/disactivate diagnostics in envar() en_data routines=========
 Diag=.true.
 if(iene==0)then
  Diag = .false.
  iene=1
 endif
 if(Part)then
  if(prl)then
   call part_numbers
   call max_pmemory_check()
  endif
 endif
 !if (pe0)then
 ! call initial_run_info(new_sim)
 !endif
!=============================
 call cpu_time(unix_time_now)
 unix_time_begin=unix_time_now
 unix_time_last_dump=unix_time_begin

 tnow=tstart
 ! iter=last_iter
 iter=0

 select case(mod_ord)
 case(1)
  call LP_cycle
 case(2)
  call ENV_cycle
 case(3)
  call BUNCH_cycle
 end select

 !call timing
 call mpi_barrier(comm,error)
 call final_run_info
 call end_parallel
 !--------------------------

 contains

 !--------------------------

 subroutine LP_cycle

 call data_out(jump)
 dt_loc=dt

 do while (tnow < tmax)
!=======================
  call LP_run(tnow,dt_loc,iter)

  call timing
  call data_out(jump)

  if (ier /= 0) then
   call error_message
   exit
  endif
  if (tnow+dt_loc >= tmax) dt_loc=tmax-tnow
 end do
 if (dump>0) call dump_data(iter,tnow)
 end subroutine LP_cycle

 !--------------------------

 subroutine ENV_cycle

  if(inject_ind==iter)then
   call beam_inject
   if(pe0)write(6,'(a24,e11.4)')' Injected beam at time =',tnow
  endif
 call data_out(jump)
!================
 tk_ind=0
 do while (tnow < tmax)
!=================================
  call ENV_run(tnow,dt_loc,iter)

  call timing           !iter=iter+1  tnow=tnow+dt_loc
  call data_out(jump)

  if (ier /= 0) then
   call error_message
   exit
  endif
  if (tnow+dt_loc >= tmax) dt_loc=tmax-tnow
 end do
 if (dump>0) call dump_data(iter,tnow)
 end subroutine ENV_cycle

 !--------------------------

 subroutine BUNCH_cycle

 !if(L_intdiagnostics_pwfa) then
  !call bunch_output_struct(tdia,dtdia,tout,dtout)
 !endif
 call bdata_out(jump)
 dt_loc=dt
 t_ind=0

 do while (tnow < tmax)
  call BUNCH_run(tnow,dt_loc,iter)

  call timing
  call bdata_out(jump)

  if (ier /= 0) then
   call error_message
   exit
  endif
  if (tnow+dt_loc >= tmax) dt_loc=tmax-tnow
 end do
 if (dump>0) call dump_data(iter,tnow)

 end subroutine BUNCH_cycle
 !======================
 subroutine data_out(jump)
 integer,intent(in) :: jump
 integer :: i,ic,idata


 idata=iout
 if(Diag)then
  if (tnow>=tdia) then
   ienout=ienout+1
   call envar(ienout,tnow)
   tdia=tdia+dtdia
   if (pe0) write(6,'(a10,i3,a10,e11.4)')' rms data ',ienout,' at time =',tnow
  endif
 endif

 if (tnow>=tout) then
  call create_timestep_folder(iout)
  tout=tout+dtout
  if(Diag)then
   if(pe0)call en_data(ienout,iter,idata)
  endif
!==================
  if(nvout>0)then
   if (mod_ord==2) then
    if(L_env_modulus)then
     i=0
     call env_fields_out(env,tnow,i,jump)
     if(Two_color)call env_fields_out(env1,tnow,-1,jump)    !EXIT |A|
    else
     if(Two_color)then
      do i=1,2
       call env_two_fields_out(env,env1,tnow,i,jump)
      end do
     else
      do i=1,2
       call env_fields_out(env,tnow,i,jump)     !EXIT [Ar,Ai]
      end do
     endif
    endif
   endif
   do i=1,nvout
    if(L_force_singlefile_output) then
     call fields_out(ebf,tnow,i,i,jump)      !i to label field name
    else
     call fields_out_new(ebf,tnow,i,i,jump)
    endif
   end do
  endif
  if(Part)then
   if (nden>0) then
    do i=1,nsp
     call prl_den_energy_interp(i)
     do ic=1,min(2,nden)
      call den_energy_out(tnow,i,ic,ic,jump)
     end do
    enddo
    if(nden > 2)then
     call set_wake_potential
     call den_energy_out(tnow,0,nden,1,jump)  !data on jc(1) for wake potential
    endif
   endif
  endif
  if(Hybrid)then
   do i=1,nfcomp
    call fluid_den_mom_out(up,tnow,i,nfcomp,jump)
   end do
  endif
  if(Ionization)call part_ionz_out(tnow)
  if(gam_min >1.)call part_high_gamma_out(gam_min,tnow)
  if (npout>0) then
   if (npout<=nsp) then
    call part_pdata_out(tnow,xp0_out,xp1_out,yp_out,npout,pjump)
   else
    do i=1,nsp
     call part_pdata_out(tnow,xp0_out,xp1_out,yp_out,i,pjump)
    end do
   endif
  endif
  call cpu_time(unix_time_now)
  if (pe0) then
   write(6,'(a10,i6,a10,e11.4,a10,e11.4)') 'iter = ',iter,' t = ',tnow,' dt = ',dt_loc
   write(6,*)' END DATA WRITE'
   write(6,'(a16,f12.3,a10,i15)')' Time elapsed = ',unix_time_now-unix_time_begin
  endif
  !if (dump>0 .and. time_interval_dumps < 0.0) then
  ! if (iter>0) call dump_data(iter,tnow)
  !endif
  iout=iout+1
 endif

 call cpu_time(unix_time_now)

 !if((unix_time_now - unix_time_last_dump) > time_interval_dumps .and. time_interval_dumps > 0.0) then
 ! call dump_data(iter,tnow)
 !endif

 end subroutine data_out


 !--------------------------
 subroutine bdata_out(jump)

 integer,intent(in) :: jump
 integer :: i,ic,idata

 idata=iout

 if(Diag)then
  if (tnow>=tdia) then
   ienout=ienout+1
   !uncomment for extracting lineout in enbvar
   !call collect_bunch_and_plasma_density(nsb,1) !data is now stored on jc(1)
   !---gather data for lineout---!
   !---> lineout total density
   if(L_intdiagnostics_pwfa) then
    call collect_bunch_and_plasma_density(0,1)
   endif
   !---!
   call enbvar(ienout,tnow)
   !---!
   call envar(ienout,tnow)
   tdia=tdia+dtdia
   if (pe0) then
    write(6,'(a10,i3,a10,e11.4)')' rms data ',ienout,' at time =',tnow
   endif
  endif
 endif

 if (tnow>=tout) then
  tout=tout+dtout

  call create_timestep_folder(iout)

  if(Diag)then
   if(pe0)then
    call en_data(ienout,iter,idata)
    call en_bdata(ienout,iter,idata)
   endif
  endif

  if(nvout>0)then
   do i=1,min(nvout,nfield)
    if(L_force_singlefile_output) then
     call fields_out(ebf,tnow,i,i,jump)    ! second index to label field
    else
     call fields_out_new(ebf,tnow,i,i,jump)
    endif
   end do

   if      (L_print_J_on_grid .AND. L_force_singlefile_output) then
     call fields_out(jc,tnow,1,0,jump)     ! 0 for Jx current
   else if  (L_print_J_on_grid .AND. L_force_singlefile_output) then
     call fields_out_new(jc,tnow,1,0,jump) ! 0  to label Jx field
   endif

   do i=1,nbfield
    if(L_force_singlefile_output) then
     call bfields_out(ebf_bunch,ebf1_bunch,tnow,i,jump)
    else
     call fields_out_new(ebf_bunch,tnow,i,i+6,jump)
    endif
   end do
  endif
  if(nden>0)then
   if(Hybrid)then
    do i=1,nfcomp
     call fluid_den_mom_out(up,tnow,i,nfcomp,jump)
    end do
   endif
   ic=0
   call prl_bden_energy_interp(ic)
   call bden_energy_out(tnow,1,jump)
   !============= bunch density
   do i=1,nsp
    call prl_den_energy_interp(i)
    do ic=1,min(2,nden)
     call den_energy_out(tnow,i,ic,ic,jump)
    end do
   enddo
   if(nden>2)then
    call set_wake_potential
    call den_energy_out(tnow,0,nden,1,jump)  !data on jc(1) for wake potential
   endif
  endif
  if(nbout> 0)then
   do i=1,nsb
    call part_bdata_out(tnow,i,pjump)
   end do
  endif
  if(Part)then
   if(Ionization)call part_ionz_out(tnow)
   if(gam_min >1.)call part_high_gamma_out(gam_min,tnow)
   if(npout> 0)then
    if(npout<=nsp)then
     call part_pdata_out(tnow,xp0_out,xp1_out,yp_out,npout,pjump)
    else
     do i=1,nsp
      call part_pdata_out(tnow,xp0_out,xp1_out,yp_out,i,pjump)
     end do
    endif
   endif
  endif
  if (pe0) then
   write(6,'(a10,i6,a10,e11.4,a10,e11.4)') 'iter = ',iter,' t = ',tnow,' dt = ',dt_loc
   write(6,*)' END B-DATA WRITE'
  endif
  if(dump>0 .and. time_interval_dumps < 0.)then
   if(iter>0)call dump_data(iter,tnow)
  endif
  iout=iout+1
 endif

 call cpu_time(unix_time_now)

 if((unix_time_now - unix_time_last_dump) > time_interval_dumps .and. time_interval_dumps > 0.0) then
  if(pe0) write(6,*) '3D dump data being written'
  call dump_data(iter,tnow)
 endif
 end subroutine bdata_out
 !--------------------------
 end program ALaDyn
