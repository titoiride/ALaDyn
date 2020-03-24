!*****************************************************************************************************!
!                            Copyright 2008-2020  The ALaDyn Collaboration                            !
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

module tracking

 contains
   !=============== Tracking particles============
  subroutine initial_tparticles_select(tx1, ty1)
  
   real (dp), intent (in) :: tx1, ty1
  
   integer :: np, p, ik, ndv, ik_max
   integer (hp_int) :: plab, last_ind
   integer, parameter :: tp_numb = 10
   real (dp) :: yy(tp_numb), ym, ymx
   !========================
   ! Define track_tot_nstep
   p = 0
   if (dt_loc>0.0) p = nint((t_out-t_in)/dt_loc)
   track_tot_nstep = nint(real(p,dp)/real(tkjump,dp))
   ndv = nd2 + 1
  
   ! Select particles on each mpi_task
   ! xp defined by tx1
   ! yp [-4,-3,-2.-1.,0,1,2,3,4], 9 test particles allocated
   !
   yy(1) = ty1
   do p = 2, tp_numb
    yy(p) = yy(p-1) + 1.
   end do
   ym = loc_ygrid(imody)%gmin
   ymx = loc_ygrid(imody)%gmax
   np = loc_npart(imody, imodz, imodx, 1)
   ik = 0
   if (allocated(spec(1)%part)) then
    deallocate (spec(1)%part)
    deallocate (ebfp)
    loc_npart(imody, imodz, imodx, 1) = 0
   end if
   select case (ndim)
   case (2)
    do p = 1, tp_numb
     if (ym<yy(p) .and. ymx>=yy(p)) then
      ik = ik + 1
     end if
    end do
    if (ik>0) then
     loc_npart(imody, imodz, imodx, 1) = ik
     allocate (spec(1)%part(ik,5))
     np = 0
     do p = 1, tp_numb
      if (ym<yy(p) .and. ymx>=yy(p)) then
       np = np + 1
       spec(1)%part(np, 1) = tx1
       spec(1)%part(np, 2) = yy(p)
       spec(1)%part(np, 3:4) = t0_pl(1)
       part_ind = int(np, hp_int)
       wgh = real(0.0, sp)
       charge = int(-1., hp_int)
       spec(1)%part(np, 5) = wgh_cmp
      end if
     end do
    end if
   case (3)
    return
   end select
   np = loc_npart(imody, imodz, imodx, 1)
   loc_tpart(mype+1) = np
   call intvec_distribute(np, loc_tpart, npe)
   track_tot_part = sum(loc_tpart(1:npe))
   ik_max = maxval(loc_tpart(1:npe))
   last_ind = 0
   if (mype==1) last_ind = loc_tpart(1)
   if (mype>1) last_ind = sum(loc_tpart(1:mype))
   !if(loc_tpart(mype+1)>0)write(6,*)'last particle index',mype,last_ind
   if (np>0) then
    do p = 1, np
     wgh_cmp = spec(1)%part(p, ndv)
     if (part_ind>0) part_ind = part_ind + last_ind
     spec(1)%part(p, ndv) = wgh_cmp
    end do
   end if
   allocate (track_aux(2*ndv*ik_max))
   if (.not. allocated(ebfp)) allocate (ebfp(ik_max,ndv))
   if (pe0) then
    allocate (pdata_tracking(ndv,track_tot_part,track_tot_nstep))
    write (6, *) '==== Initial track-Particle data==========='
    write (6, '(a19,i6)') '  tot_track_steps  ', track_tot_nstep
    write (6, '(a19,2i6)') '  tot_track_parts  ', track_tot_part, &
      tp_numb
    write (6, '(a18,i8)') '  short_int size  ', huge(plab)
    write (6, '(a20,e11.4)') '  track memory(MB) ', &
      1.e-06*real(4*ndv*track_tot_nstep*track_tot_part, dp)
   end if
  end subroutine

  subroutine t_particles_collect(time_ind)

   integer, intent (in) :: time_ind

   integer :: np, ik, ik_max, ip, p, ndv, ndvp, ipe, kk, ik1, ik2, nst
   logical :: sr
   real :: xm, ym, zm

   if (time_ind>track_tot_nstep) return
   xm = loc_xgrid(imodx)%gmin
   ym = loc_ygrid(imody)%gmin
   zm = loc_zgrid(imodz)%gmin
   np = loc_npart(imody, imodz, imodx, 1)
   nst = 0
   ndv = nd2 + 1
   ndvp = ndv
   if (stretch) nst = str_indx(imody, imodz)
   !===================
   np = loc_npart(imody, imodz, imodx, 1)
   ik = 0
   !=========== each pe counts track particle number
   do p = 1, np
    wgh_cmp = spec(1)%part(p, ndv)
    if (part_ind>0) ik = ik + 1
   end do
   loc_tpart(mype+1) = ik
   call intvec_distribute(ik, loc_tpart, npe)
   ik_max = maxval(loc_tpart(1:npe))
   if (ndvp*ik_max>size(track_aux)) then
    deallocate (track_aux)
    allocate (track_aux(ndvp*ik_max))
   end if
   !========= each pe stores tpart data in track_aux(ndv,loc_tpart)
   kk = 0
   do p = 1, np
    wgh_cmp = spec(1)%part(p, ndv)
    if (part_ind>0) then
     do ip = 1, ndv
      kk = kk + 1
      track_aux(kk) = spec(1)%part(p, ip)
     end do
    end if
   end do
   !=================
   !  all data are gathered onto pe0 
   if (pe0) then
    sr = .false.
    ik1 = 0
    do ipe = 1, npe - 1
     ik = loc_tpart(ipe+1)
     if (ik>0) then
      call exchange_1d_grdata(sr, track_aux, ik*ndvp, ipe, ipe+10)
      !pe0 receives from ipe ik sp_aux data and collects on track array
      kk = 0
      do p = 1, ik
       ik2 = ik1 + p
       do ip = 1, ndv - 1
        kk = kk + 1
        pdata_tracking(ip, ik2, time_ind) = track_aux(kk)
        !pdata_tracking(ip,ik2,time_ind)=track_aux(kk)
       end do
       kk = kk + 1
       wgh_cmp = track_aux(kk)
       pdata_tracking(ndv, ik2, time_ind) = real(part_ind, dp)
      end do
      ik1 = ik2
     end if
    end do
    loc_tpart(1) = ik2
   else
    ik = loc_tpart(mype+1)
    if (ik>0) then
     sr = .true.
     call exchange_1d_grdata(sr, track_aux, ik*ndvp, 0, mype+10) !sends ik data to pe0
    end if
   end if
  end subroutine
  !=================================================
 end module
