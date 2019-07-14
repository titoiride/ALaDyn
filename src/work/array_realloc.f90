!=============================================
 module array_ralloc

  use array_wspace
  implicit none

!--------------------------
!DIR$ ATTRIBUTES INLINE :: p_realloc
 contains
  subroutine p_realloc(pdata, npt_new, ndv)
   type (species), intent (inout) :: pdata
   integer, intent (in) :: npt_new, ndv
   integer :: allocstatus, deallocstatus

   if (allocated(pdata%part)) then
    if (size(pdata%part,1)<npt_new) then
     deallocate (pdata%part, stat=deallocstatus)
     if (deallocstatus==0) allocate (pdata%part(1:npt_new,1:ndv), &
       stat=allocstatus)
    end if
   else
    allocate (pdata%part(1:npt_new,1:ndv), stat=allocstatus)
   end if
   pdata%part(:, :) = 0.0
  end subroutine
!========================
  subroutine v_realloc(vdata, npt_new, ndv)
   real (dp), allocatable, intent (inout) :: vdata(:, :)
   integer, intent (in) :: npt_new, ndv
   integer :: allocstatus, deallocstatus

   if (allocated(vdata)) then
    if (size(vdata,1)<npt_new) then
     deallocate (vdata, stat=deallocstatus)
     allocate (vdata(1:npt_new,1:ndv), stat=allocstatus)
    end if
   else
    allocate (vdata(1:npt_new,1:ndv), stat=allocstatus)
   end if
   vdata(:, :) = 0.0
  end subroutine
!===========================
 end module

