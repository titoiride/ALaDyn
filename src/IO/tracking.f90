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

 use parallel
 use pstruct_data
 use precision_def
 use system_utilities
 use array_alloc, only: array_realloc_1d
 use grid_part_lib, only: xx_realloc
 use util, only: endian
 use warnings, only: write_warning

 real(dp), allocatable, dimension(:, :), private, save :: track_aux
 real(sp), allocatable, dimension(:), private, save :: track_pdata
 integer, private, save :: iter_index
 logical, allocatable, dimension(:), private, save :: track_mask
 contains

  subroutine generate_track_index( spec_in, index_in, new_parts )
   type(species_new), intent(in) :: spec_in
   integer, allocatable, dimension(:), intent(inout) :: index_in
   integer, intent(in) :: new_parts
   integer :: max_index_loc(npe), new_parts_loc(npe)
   integer :: max_index, last_collective_index, tot_newparts, cum_sum
   type( track_data_t ) :: track_data
   type( index_array) :: new_inds

   track_data = spec_in%pick_track_params()
   max_index = track_data%highest_index
   max_index_loc( mype + 1 ) = max_index
   new_parts_loc( mype + 1 ) = new_parts

   call intvec_distribute( max_index, max_index_loc, npe)
   call intvec_distribute( new_parts, new_parts_loc, npe)

   last_collective_index = MAXVAL( max_index_loc(:) )

   tot_newparts = SUM( new_parts_loc(:) )

   if ( tot_newparts == 0 ) return

   new_inds = index_array(tot_newparts, lb=last_collective_index + 1)

   allocate ( index_in( new_parts ) )
   cum_sum = SUM( new_parts_loc(1:mype) )

   index_in(1:new_parts) = new_inds%indices((cum_sum + 1):( cum_sum + new_parts))

  end subroutine

  subroutine initialize_tracking( spec_in, spec_aux_in )
   type(species_new), dimension(:), intent(inout) :: spec_in
   type(species_aux), intent(inout) :: spec_aux_in
   logical :: tracking
   integer :: ic

   tracking = ANY(p_tracking, DIM=1)
   if (.not. tracking) return

   call create_tracking_folders

   iter_index = 0
   do ic = 1, nsp
    call spec_in(ic)%track( p_tracking(ic), allocate_now=.true. )
    call spec_in(ic)%set_tot_tracked_parts( 0 )
   end do

   call spec_aux_in%track( p_tracking(1), allocate_now=.true. )
   do ic = 1, nsp
    if ( spec_in(ic)%istracked() ) then
     call spec_in(ic)%set_track_params( txmin(ic), txmax(ic), &
                                        tymin(ic), tymax(ic), &
                                        tzmin(ic), tzmax(ic), &
                                        nkjump(ic))
    end if
   end do

   do ic = 1, nsp
    call initialize_particle_index( spec_in, ic)
   end do

  end subroutine

  !=============== Initialize Tracking particles============
  subroutine initialize_particle_index( spec_in, ic )
   type(species_new), dimension(:), intent(inout) :: spec_in
   integer, intent(in) :: ic
   integer, allocatable, dimension(:) :: parts, track_index
   integer :: np, part_jump
   type (track_data_t) :: track_data

   if ( .not. spec_in(ic)%istracked() ) return
   call xx_realloc(track_aux, np, spec_in(ic)%total_size())

   part_jump = 1
   np = spec_in(ic)%how_many()
   call array_realloc_1d( track_mask, np )
   track_mask(1:np) = .true.
   allocate ( parts(np), source=0 )

   track_data = spec_in(ic)%pick_track_params()
   part_jump = track_data%jump

   if (ndim>2) then

    associate( xx => spec_in(ic)%call_component(X_COMP, lb=1, ub=np), &
               yy => spec_in(ic)%call_component(Y_COMP, lb=1, ub=np), &
               zz => spec_in(ic)%call_component(Z_COMP, lb=1, ub=np))
               
               if ( allocated(track_data%xmin) ) then
                track_mask(1:np) = ( xx >= track_data%xmin .and. track_mask )
               end if
               if ( allocated(track_data%xmax) ) then
                track_mask(1:np) = ( xx <= track_data%xmax .and. track_mask )
               end if
               if ( allocated(track_data%ymin) ) then
                track_mask(1:np) = ( yy >= track_data%ymin .and. track_mask )
               end if
               if ( allocated(track_data%ymax) ) then
                track_mask(1:np) = ( yy <= track_data%ymax .and. track_mask )
               end if
               if ( allocated(track_data%zmin) ) then
                track_mask(1:np) = ( z >= track_data%zmin .and. track_mask )
               end if
               if ( allocated(track_data%zmax) ) then
                track_mask(1:np) = ( zz <= track_data%zmax .and. track_mask )
               end if
    end associate

    npt = COUNT(track_mask(1:np))

   else
    associate( xx => spec_in(ic)%call_component(X_COMP, lb=1, ub=np), &
               yy => spec_in(ic)%call_component(Y_COMP, lb=1, ub=np))

               if ( allocated(track_data%xmin) ) then
                track_mask(1:np) = ( xx >= track_data%xmin .and. track_mask )
               end if
               if ( allocated(track_data%xmax) ) then
                track_mask(1:np) = ( xx <= track_data%xmax .and. track_mask )
               end if
               if ( allocated(track_data%ymin) ) then
                track_mask(1:np) = ( yy >= track_data%ymin .and. track_mask )
               end if
               if ( allocated(track_data%ymax) ) then
                track_mask(1:np) = ( yy <= track_data%ymax .and. track_mask )
               end if

    end associate

    npt = COUNT(track_mask(1:np))

   end if

   call generate_track_index( spec_in(ic), track_index, npt)

   parts = UNPACK( track_index(1:npt), track_mask(1:np), parts(1:np) )

   write(6, *) '======================'
   write(6, *) 'Initial npt = ', npt
   write(6, *) '======================'
   call spec_in(ic)%set_component(parts, INDEX_COMP, lb=1, ub=np)
   call spec_in(ic)%set_tot_tracked_parts( npt )
   call spec_in(ic)%set_highest_track_index( MAXVAL(track_index) )

  end subroutine

  subroutine tracking_write_output(  spec_in, timenow, iter, pid )

  type(species_new), intent(in) :: spec_in
  real (dp), intent (in) :: timenow
  integer, intent (in) :: iter, pid

  character (12) :: fname
  character (8) :: foldername
  character (25) :: fname_out

  integer (dp) :: nptot_global_reduced
  integer :: npt, npt_loc(npe), cc
  integer :: ik, p, q, np, ip, ip_max, nptot, n_comp_out
  integer :: lenp, ip_loc(npe), ndv, i_end
  integer (offset_kind) :: disp
  type(index_array) :: out_parts
  type(track_data_t) :: track_data

  write (foldername, '(a8)') 'tracking'
  write( fname, '(a6, i1.1, a1, i4.4)') 'Track_', pid, '_', iter_index

  fname_out = foldername // '/' // fname // '.bin'

  ndv = spec_in%total_size()
  np = spec_in%how_many()
  n_comp_out = ndv

  track_data = spec_in%pick_track_params()

  npt = track_data%n_tracked

  npt_loc(mype + 1) = npt

  call intvec_distribute(npt, npt_loc, npe)

  if (ALL(npt_loc <= 0)) return
  ch = spec_in%pick_charge()
  out_parts = index_array(np)

  call xx_realloc(track_aux, npt, spec_in%total_size())
  call array_realloc_1d( track_mask, np )


  associate( inds => spec_in%call_component(INDEX_COMP, lb=1, ub=np) )

   track_mask(1:np) = (int(inds) > 0)

  end associate
  
  cc = COUNT(track_mask(1:np) )
  if (npt /= cc ) then
   call gdbattach
   call write_warning(text='Error in counting tracked particles', task=mype)
  end if

  call out_parts%find_index(track_mask(1:np))
  call spec_in%call_particle(track_aux, out_parts%indices(1:npt), tracking=.true.)

  ip_loc(mype+1) = npt
  ip = ip_loc(mype+1)

  call intvec_distribute(ip, ip_loc, npe)

  ! this differs from nptot_global since it represents just the reduced number of particles
  ! that will be present in the output (should be equal to nptot_global for p_jump=1)!
  nptot_global_reduced = 0
  do ik = 1, npe
   nptot_global_reduced = nptot_global_reduced + ip_loc(ik)
  end do
  if (nptot_global_reduced < 1e9) then
   nptot = int(nptot_global_reduced)
  else
   nptot = -1
  end if

  ip_max = ip
  if (pe0) ip_max = maxval(ip_loc(1:npe))
  lenp = n_comp_out*ip_loc(mype+1)
  call array_realloc_1d(track_pdata, lenp)
  ik = 0
  do p = 1, ip_loc(mype+1)
   do q = 1, n_comp_out
    ik = ik + 1
    track_pdata(ik) = real(track_aux(p,q), sp)
   end do
  end do
  if (ik /= lenp) write (6, '(a16,3i8)') 'wrong pdata size', mype, &
   lenp, ik

  call endian(i_end)

  disp = 0
  disp_col = 0
  if ( .not. pe0) then
   disp = mype + n_comp_out*sum(ip_loc(1:mype)) ! da usare con mpi_write_part
  end if

  disp = disp*4 ! sia gli int che i float sono di 4 bytes

  call mpi_write_part(track_pdata, lenp, ip, disp, 25, fname_out)
  iter_index = iter_index + 1
  if (pe0) then
   write(6, *)            '==========================================='
   write(6, '(a35,i2.2)') ' Tracking data written for species ', pid
   write(6, *)            '==========================================='
  end if
 end subroutine
 !================================
 subroutine track_out( spec_in, timenow, iter)
  type(species_new), intent(in), dimension(:) :: spec_in
  real(dp), intent(in) :: timenow
  integer, intent(in) :: iter
  integer :: ic

  do ic = 1, nsp
   if ( MOD(iter, every_track(ic) ) == 0 .and. spec_in(ic)%istracked()) then
    call tracking_write_output(spec_in(ic), timenow, iter, ic)
   end if
  end do
 end subroutine
   !=============== Tracking particles============
  ! subroutine initial_tparticles_select(tx1, ty1)
  
  !  real (dp), intent (in) :: tx1, ty1
  
  !  integer :: np, p, ik, ndv, ik_max
  !  integer (hp_int) :: plab, last_ind
  !  integer, parameter :: tp_numb = 10
  !  real (dp) :: yy(tp_numb), ym, ymx
  !  !========================
  !  ! Define track_tot_nstep
  !  p = 0
  !  if (dt_loc>0.0) p = nint((t_out-t_in)/dt_loc)
  !  track_tot_nstep = nint(real(p,dp)/real(tkjump,dp))
  !  ndv = nd2 + 1
  
  !  ! Select particles on each mpi_task
  !  ! xp defined by tx1
  !  ! yp [-4,-3,-2.-1.,0,1,2,3,4], 9 test particles allocated
  !  !
  !  yy(1) = ty1
  !  do p = 2, tp_numb
  !   yy(p) = yy(p-1) + 1.
  !  end do
  !  ym = loc_ygrid(imody)%gmin
  !  ymx = loc_ygrid(imody)%gmax
  !  np = loc_npart(imody, imodz, imodx, 1)
  !  ik = 0
  !  if (allocated(spec(1)%part)) then
  !   deallocate (spec(1)%part)
  !   deallocate (ebfp)
  !   loc_npart(imody, imodz, imodx, 1) = 0
  !  end if
  !  select case (ndim)
  !  case (2)
  !   do p = 1, tp_numb
  !    if (ym<yy(p) .and. ymx>=yy(p)) then
  !     ik = ik + 1
  !    end if
  !   end do
  !   if (ik>0) then
  !    loc_npart(imody, imodz, imodx, 1) = ik
  !    allocate (spec(1)%part(ik,5))
  !    np = 0
  !    do p = 1, tp_numb
  !     if (ym<yy(p) .and. ymx>=yy(p)) then
  !      np = np + 1
  !      spec(1)%part(np, 1) = tx1
  !      spec(1)%part(np, 2) = yy(p)
  !      spec(1)%part(np, 3:4) = t0_pl(1)
  !      part_ind = int(np, hp_int)
  !      wgh = real(0.0, sp)
  !      charge = int(-1., hp_int)
  !      spec(1)%part(np, 5) = wgh_cmp
  !     end if
  !    end do
  !   end if
  !  case (3)
  !   return
  !  end select
  !  np = loc_npart(imody, imodz, imodx, 1)
  !  loc_tpart(mype+1) = np
  !  call intvec_distribute(np, loc_tpart, npe)
  !  track_tot_part = sum(loc_tpart(1:npe))
  !  ik_max = maxval(loc_tpart(1:npe))
  !  last_ind = 0
  !  if (mype==1) last_ind = loc_tpart(1)
  !  if (mype>1) last_ind = sum(loc_tpart(1:mype))
  !  !if(loc_tpart(mype+1)>0)write(6,*)'last particle index',mype,last_ind
  !  if (np>0) then
  !   do p = 1, np
  !    wgh_cmp = spec(1)%part(p, ndv)
  !    if (part_ind>0) part_ind = part_ind + last_ind
  !    spec(1)%part(p, ndv) = wgh_cmp
  !   end do
  !  end if
  !  allocate (track_aux(2*ndv*ik_max))
  !  if (.not. allocated(ebfp)) allocate (ebfp(ik_max,ndv))
  !  if (pe0) then
  !   allocate (pdata_tracking(ndv,track_tot_part,track_tot_nstep))
  !   write (6, *) '==== Initial track-Particle data==========='
  !   write (6, '(a19,i6)') '  tot_track_steps  ', track_tot_nstep
  !   write (6, '(a19,2i6)') '  tot_track_parts  ', track_tot_part, &
  !     tp_numb
  !   write (6, '(a18,i8)') '  short_int size  ', huge(plab)
  !   write (6, '(a20,e11.4)') '  track memory(MB) ', &
  !     1.e-06*real(4*ndv*track_tot_nstep*track_tot_part, dp)
  !  end if
  ! end subroutine

  ! subroutine t_particles_collect(time_ind)

  !  integer, intent (in) :: time_ind

  !  integer :: np, ik, ik_max, ip, p, ndv, ndvp, ipe, kk, ik1, ik2, nst
  !  logical :: sr
  !  real :: xm, ym, zm

  !  if (time_ind>track_tot_nstep) return
  !  xm = loc_xgrid(imodx)%gmin
  !  ym = loc_ygrid(imody)%gmin
  !  zm = loc_zgrid(imodz)%gmin
  !  np = loc_npart(imody, imodz, imodx, 1)
  !  nst = 0
  !  ndv = nd2 + 1
  !  ndvp = ndv
  !  if (stretch) nst = str_indx(imody, imodz)
  !  !===================
  !  np = loc_npart(imody, imodz, imodx, 1)
  !  ik = 0
  !  !=========== each pe counts track particle number
  !  do p = 1, np
  !   wgh_cmp = spec(1)%part(p, ndv)
  !   if (part_ind>0) ik = ik + 1
  !  end do
  !  loc_tpart(mype+1) = ik
  !  call intvec_distribute(ik, loc_tpart, npe)
  !  ik_max = maxval(loc_tpart(1:npe))
  !  if (ndvp*ik_max>size(track_aux)) then
  !   deallocate (track_aux)
  !   allocate (track_aux(ndvp*ik_max))
  !  end if
  !  !========= each pe stores tpart data in track_aux(ndv,loc_tpart)
  !  kk = 0
  !  do p = 1, np
  !   wgh_cmp = spec(1)%part(p, ndv)
  !   if (part_ind>0) then
  !    do ip = 1, ndv
  !     kk = kk + 1
  !     track_aux(kk) = spec(1)%part(p, ip)
  !    end do
  !   end if
  !  end do
  !  !=================
  !  !  all data are gathered onto pe0 
  !  if (pe0) then
  !   sr = .false.
  !   ik1 = 0
  !   do ipe = 1, npe - 1
  !    ik = loc_tpart(ipe+1)
  !    if (ik>0) then
  !     call exchange_1d_grdata(sr, track_aux, ik*ndvp, ipe, ipe+10)
  !     !pe0 receives from ipe ik sp_aux data and collects on track array
  !     kk = 0
  !     do p = 1, ik
  !      ik2 = ik1 + p
  !      do ip = 1, ndv - 1
  !       kk = kk + 1
  !       pdata_tracking(ip, ik2, time_ind) = track_aux(kk)
  !       !pdata_tracking(ip,ik2,time_ind)=track_aux(kk)
  !      end do
  !      kk = kk + 1
  !      wgh_cmp = track_aux(kk)
  !      pdata_tracking(ndv, ik2, time_ind) = real(part_ind, dp)
  !     end do
  !     ik1 = ik2
  !    end if
  !   end do
  !   loc_tpart(1) = ik2
  !  else
  !   ik = loc_tpart(mype+1)
  !   if (ik>0) then
  !    sr = .true.
  !    call exchange_1d_grdata(sr, track_aux, ik*ndvp, 0, mype+10) !sends ik data to pe0
  !   end if
  !  end if
  ! end subroutine
  ! !=================================================

 end module
