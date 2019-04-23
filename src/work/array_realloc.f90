!=============================================
 module array_ralloc

 use array_wspace
 implicit none

 !--------------------------
 !DIR$ ATTRIBUTES INLINE :: p_realloc
 contains
 subroutine p_realloc(pdata,npt_new,ndv)
  type(species),intent(inout)  :: pdata
  integer,intent(in) :: npt_new,ndv
  integer :: AllocStatus, DeallocStatus

 if(allocated(pdata%part))then
  if(size(pdata%part,1) < npt_new)then
   deallocate(pdata%part,STAT=DeallocStatus)
   if(DeallocStatus==0)allocate(pdata%part(1:npt_new,1:ndv),STAT=AllocStatus)
  endif
 else
  allocate(pdata%part(1:npt_new,1:ndv),STAT=AllocStatus)
 endif
 pdata%part(:,:)=0.0
 end subroutine p_realloc
 !========================
 subroutine v_realloc(vdata,npt_new,ndv)
 real(dp),allocatable,intent(inout) :: vdata(:,:)
 integer,intent(in) :: npt_new,ndv
 integer :: AllocStatus, DeallocStatus

 if(allocated(vdata))then
  if(size(vdata,1) < npt_new)then
   deallocate(vdata,STAT=DeallocStatus)
   allocate(vdata(1:npt_new,1:ndv),STAT=AllocStatus)
  endif
 else
  allocate(vdata(1:npt_new,1:ndv),STAT=AllocStatus)
 endif
 vdata(:,:)=0.0
 end subroutine v_realloc
 !===========================
 end module array_ralloc

