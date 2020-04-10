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
 !! Module to manage the tracked particles

 use common_param, only: ref_nspec
 use parallel
 use pstruct_data
 use phys_param
 use system_utilities
! !#if defined(OPENPMD)
!  use hdf5io_class
! !#endif
 use array_alloc, only: array_realloc_1d
 use grid_part_lib, only: xx_realloc
 use util, only: endian
 use warnings, only: write_warning

 real(dp), allocatable, dimension(:, :), private, save :: track_aux
 real(sp), allocatable, dimension(:), private, save :: track_pdata
 integer, private, save, dimension(ref_nspec) :: iter_index
 logical, allocatable, dimension(:), private, save :: track_mask
 logical, save :: tracking_written
 character(len=8), public, parameter :: tracking_folder = 'tracking'
 character(len=25), private,&
  dimension(ref_nspec) :: track_dic
 integer, private, parameter, &
  dimension(ref_nspec) :: track_iounit = [(50 + i - 1, i = 1, ref_nspec)]

 contains

 !====== PARTICLE TRACKING UTILITIES ==========

  subroutine generate_track_index( spec_in, index_in, new_parts )
   !! This subroutine generate a shared array for particle index
   !! to avoid duplicated indexes.
   !! @warning
   !! This requires a sendrecieve, so it should be executed by all MPI tasks
   !! @endwarning

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
   
   allocate ( index_in( new_parts ), source=0 )
   if ( tot_newparts == 0 ) return

   new_inds = index_array(tot_newparts, lb=last_collective_index + 1)

   cum_sum = SUM( new_parts_loc(1:mype) )

   index_in(1:new_parts) = new_inds%indices((cum_sum + 1):( cum_sum + new_parts))

  end subroutine
  
  subroutine set_tracked_particle_index( spec_in, lowb, upb, ic)
   !! Assigns to every new tracked particle a unique index for identification

   type(species_new), dimension(:), intent(inout) :: spec_in
   integer, intent(in) :: lowb, upb, ic
   integer, allocatable, dimension(:) :: parts, track_index, temp_track_index
   integer :: npt, np, part_jump, npt_jump, i
   type(track_data_t) :: track_data

   if ( .not. spec_in(ic)%istracked() ) return

   np = upb - lowb + 1

   call array_realloc_1d( track_mask, np )
   track_mask(1:np) = .true.

   allocate ( parts(np), source=0 )

   track_data = spec_in(ic)%pick_track_params()
   part_jump = track_data%jump

   if (ndim>2) then

    associate( xx => spec_in(ic)%call_component(X_COMP, lb=lowb, ub=upb), &
               yy => spec_in(ic)%call_component(Y_COMP, lb=lowb, ub=upb), &
               zz => spec_in(ic)%call_component(Z_COMP, lb=lowb, ub=upb))
               
               if ( allocated(track_data%xmin) ) then
                track_mask(1:np) = ( xx(1:np) >= track_data%xmin .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%xmax) ) then
                track_mask(1:np) = ( xx(1:np) <= track_data%xmax .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%ymin) ) then
                track_mask(1:np) = ( yy(1:np) >= track_data%ymin .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%ymax) ) then
                track_mask(1:np) = ( yy(1:np) <= track_data%ymax .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%zmin) ) then
                track_mask(1:np) = ( zz(1:np) >= track_data%zmin .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%zmax) ) then
                track_mask(1:np) = ( zz(1:np) <= track_data%zmax .and. track_mask(1:np) )
               end if
    end associate

    npt = COUNT(track_mask(1:np))

   else
    associate( xx => spec_in(ic)%call_component(X_COMP, lb=lowb, ub=upb), &
               yy => spec_in(ic)%call_component(Y_COMP, lb=lowb, ub=upb))

               if ( allocated(track_data%xmin) ) then
                track_mask(1:np) = ( xx(1:np) >= track_data%xmin .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%xmax) ) then
                track_mask(1:np) = ( xx(1:np) <= track_data%xmax .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%ymin) ) then
                track_mask(1:np) = ( yy(1:np) >= track_data%ymin .and. track_mask(1:np) )
               end if
               if ( allocated(track_data%ymax) ) then
                track_mask(1:np) = ( yy(1:np) <= track_data%ymax .and. track_mask(1:np) )
               end if

    end associate

    npt = COUNT(track_mask(1:np))

   end if

   call generate_track_index( spec_in(ic), temp_track_index, npt)

   if (np <= 0 ) return

   allocate(track_index(npt), source=0)
   do i = 1, npt, part_jump
    track_index(i) = temp_track_index(i)
   end do
   npt_jump = CEILING(real(npt)/part_jump)
   parts = UNPACK( track_index(1:npt), track_mask(1:np), parts(1:np) )

   call spec_in(ic)%set_component(parts, INDEX_COMP, lb=lowb, ub=upb)
   call spec_in(ic)%set_tot_tracked_parts( npt_jump )
   call spec_in(ic)%set_highest_track_index( MAXVAL(track_index) )

  end subroutine

 !====== PARTICLE TRACKING INITIALIZATION ==========

  subroutine initialize_tracking( spec_in, spec_aux_in )
   !! Initializes the particle tracking at the beginning of
   !! the simulation

   type(species_new), dimension(:), intent(inout) :: spec_in
   type(species_aux), intent(inout) :: spec_aux_in
   logical :: tracking
   integer :: ic

   tracking_written = .true.
   tracking = ANY(p_tracking, DIM=1)
   if (.not. tracking) return

   call create_tracking_folders(tracking_folder)

   iter_index = 0
   do ic = 1, nsp
    call spec_in(ic)%track( p_tracking(ic), allocate_now=.true. )
    call spec_in(ic)%set_tot_tracked_parts( 0 )
    write(track_dic(ic), '(a20,i1.1,a4)') 'tracking_dictionary_',ic,'.dat'
    if (pe0) then
     open(unit=track_iounit(ic), file=tracking_folder//'/'//track_dic(ic), &
      form='formatted', status='new')
     write(track_iounit(ic), '(a9,a4,a4)') 'Iteration', '    ', 'Time'
     close(track_iounit(ic))
    end if
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

  subroutine initialize_particle_index( spec_in, ic )
   !! Initializes particle indexes at the beginning of the simulation

   type(species_new), dimension(:), intent(inout) :: spec_in
   integer, intent(in) :: ic
   integer :: np

   if ( .not. spec_in(ic)%istracked() ) return
   
   np = spec_in(ic)%how_many()
   call set_tracked_particle_index( spec_in, 1, np, ic)

  end subroutine

 !====== PARTICLE TRACKING I/O ==========

  subroutine tracking_write_output(  spec_in, timenow, pid )
   !! Tracking I/O with the same strategy as part_pdata_out_new.
   !! The printed file contains the particle phase space (2*n_dimension),
   !! the particle weight, Lorentz gamma and index

  type(species_new), intent(in) :: spec_in
  real (dp), intent (in) :: timenow
  integer, intent (in) :: pid

  character (12) :: fname
  character (8) :: foldername
  character (25) :: fname_out

  integer (dp) :: nptot_global_reduced
  integer :: npt, npt_loc(npe), cc
  integer :: ik, p, q, np, ip, ip_max, nptot, n_comp_out
  integer :: lenp, ip_loc(npe), ndv, i_end
  integer (offset_kind) :: disp, disp_col
  type(index_array) :: out_parts
  type(track_data_t) :: track_data
  real(dp) :: ch

  write (foldername, '(a8)') 'tracking'
  write( fname, '(a6, i1.1, a1, i4.4)') 'Track_', pid, '_', iter_index(pid)

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
  iter_index(pid) = iter_index(pid) + 1
  ! if (pe0) then
  !  write(6, *)            '==========================================='
  !  write(6, '(a35,i2.2)') ' Tracking data written for species ', pid
  !  write(6, *)            '==========================================='
  ! end if
  if (pe0) then
   open(unit=track_iounit(pid), file=tracking_folder//'/'//track_dic(pid), &
    form='formatted', status='old', position="append", action="write")
   write(track_iounit(pid), '(i6.6, e12.5)') iter_index(pid) - 1, timenow
   close(track_iounit(pid))
  end if
 end subroutine
 !================================
 subroutine track_out( spec_in, timenow, iter)
  !! Wrapper for the tracking I/O routine

  type(species_new), intent(in), dimension(:) :: spec_in
  real(dp), intent(in) :: timenow
  integer, intent(in) :: iter
  integer :: ic

  do ic = 1, nsp
   if (spec_in(ic)%istracked()) then
    if ( MOD(iter, every_track(ic) ) == 0 ) then
     call tracking_write_output(spec_in(ic), timenow, ic)
     tracking_written = .true.
    end if
   end if
  end do
  ! Warning: both dictionary writing and tracking_written flag
  ! must be switched to array if employing multiple species tracking
 end subroutine

! !#if defined(OPENPMD)
!  subroutine track_write_hdf5( spec_in, timenow, dt_loc_in, iter, ic )
!   type( species_new), intent(in) :: spec_in
!   real(dp), intent(in) :: timenow, dt_loc_in
!   integer, intent(in) :: iter, ic
!   integer :: np, npt, i
!   type(hdf5file) :: file
!   real(sp), allocatable, dimension(:) :: prova

!   !np = spec_in%how_many()
!   np = 100
!   npt = spec_in%pick_tot_tracked_parts()

!   allocate(prova(np))

!   do i=1, np
!    prova(i) = i
!   end do
!   ! Writing particle properties
!   call file%new(iter=iter,&
!    &particleName=spec_in%pick_name(),&
!    &unitDimension=unit_dimensions_mass,&
!    &records='mass',&
!    &unitSI=electron_mass_Kg,&
!    ! &iterationEncoding='groupBased',&
!    &component='')

!   call pwpart(communicator, file, electron_mass_norm, ierr)
  
!   call file%new(iter=iter,&
!    &particleName=spec_in%pick_name(),&
!    &unitDimension=unit_dimensions_mass,&
!    &records='charge',&
!    &unitSI=e_charge_C,&
!    ! &iterationEncoding='groupBased',&
!    &component='')
  
!   call pwpart(communicator, file, electron_charge_norm, ierr)

!   ! ! Writing particle phase space
!   ! ! Gamma
!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_gamma, &
!   !  &unitSI=one_dp,&
!   !  &records='gamma',&
!   !  &component='')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(INV_GAMMA_COMP, lb=1, ub=np), sp), npt, ierr)

!   ! ! Weight
!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_weight, &
!   !  &unitSI=one_dp,&
!   !  &records='weight',&
!   !  &component='')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(W_COMP, lb=1, ub=np), sp), npt, ierr)

!   ! ! Index
!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_index, &
!   !  &unitSI=one_dp,&
!   !  &records='index',&
!   !  &component='')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(INDEX_COMP, lb=1, ub=np), sp), npt, ierr)

!   ! X coord
!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_positions, &
!   !  &unitSI=unit_SI_position,&
!   !  &records='position',&
!   !  &component='x')

!   call file%new(iter=iter,&
!   &particleName='electron',&
!   &unitDimension=unit_dimensions_positions,&
!   &records='position',&
!   &component='x')

!   call pwpart(communicator, file,&
!    real(spec_in%call_component(X_COMP, lb=1, ub=np), sp), np, ierr)

!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_positions, &
!   !  &unitSI=unit_SI_position,&
!   !  &records='positionOffset',&
!   !  &component='x')

!   ! call pwpart(communicator, file, zero_dp, ierr)

!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_momentum, &
!   !  &records='momentum',&
!   !  &component='x')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(PX_COMP, lb=1, ub=np), sp), np, ierr)

!   ! if(ndim < 2) return 
  
!   ! ! Y coord
!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_positions, &
!   !  &unitSI=unit_SI_position,&
!   !  &records='position',&
!   !  &component='y')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_component(Y_COMP, lb=1, ub=np), sp), np, ierr)

!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_positions, &
!   !  &unitSI=unit_SI_position,&
!   !  &records='positionOffset',&
!   !  &component='y')

!   ! call pwpart(communicator, file, zero_dp, ierr)

!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_momentum, &
!   !  &records='momentum',&
!   !  &component='y')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(PY_COMP, lb=1, ub=np), sp), npt, ierr)

!   ! if(ndim < 3) return 

!   ! ! Y coord
!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_positions, &
!   !  &unitSI=unit_SI_position,&
!   !  &records='position',&
!   !  &component='z')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(Z_COMP, lb=1, ub=np), sp), npt, ierr)

!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_positions, &
!   !  &unitSI=unit_SI_position,&
!   !  &records='positionOffset',&
!   !  &component='z')

!   ! call pwpart(communicator, file, zero_dp, ierr)

!   ! call file%new(iter=iter,&
!   !  &time=real(timenow, sp),&
!   !  &dt=real(dt_loc_in, sp),&
!   !  &particleName=spec_in%pick_name(),&
!   !  ! &iterationEncoding='groupBased',&
!   !  &unitDimension=unit_dimensions_momentum, &
!   !  &records='momentum',&
!   !  &component='z')

!   ! call pwpart(communicator, file,&
!   !  real(spec_in%call_tracked_component(PZ_COMP, lb=1, ub=np), sp), npt, ierr)

!  end subroutine

!#endif
 end module
